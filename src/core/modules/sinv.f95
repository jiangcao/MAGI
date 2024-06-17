module sinv
    ! selected inversion of structured sparse matrices
    ! The algorithm is based on the python library (SerinV)[https://github.com/vincent-maillou/serinv] 
    use parameters_mod,only:dp,c1i,czero,cone        
    use omp_lib
    !
    implicit none
    !
    public :: zbtasinv
    public :: zbtatrf
    public :: zbtatri
    !
    contains
    !
    ! computes an LU factorization of a complex block-tridiag-arrowhead (BTA) 
    ! matrix in place. It calls LAPACK `zgetrf` and `zgetri` functions on the
    ! dense blocks.
    subroutine zbtatrf( diag_blocksize, arrow_blocksize, n_diag_blocks, &
        A_diagonal_blocks,A_lower_diagonal_blocks,A_upper_diagonal_blocks, &
        A_arrow_bottom_blocks,A_arrow_right_blocks,A_arrow_tip_block, &
        ipiv_diagonal,ipiv_arrow_tip)
        ! in
        integer,intent(in) :: diag_blocksize, arrow_blocksize, n_diag_blocks
        ! in & out 
        complex(dp),intent(inout),dimension(:,:),target :: A_diagonal_blocks,A_lower_diagonal_blocks,A_upper_diagonal_blocks
        complex(dp),intent(inout),dimension(:,:),target :: A_arrow_bottom_blocks,A_arrow_right_blocks,A_arrow_tip_block
        ! out         
        integer,intent(out),dimension(diag_blocksize, n_diag_blocks) :: ipiv_diagonal
        integer,intent(out),dimension(arrow_blocksize)::ipiv_arrow_tip
        ! ---- local
        integer :: i,h,l,j
        complex(dp),allocatable,dimension(:,:) :: L_inv_temp, U_inv_temp
        integer :: info, nn
        integer, dimension(:), allocatable :: ipiv        
        complex(dp), dimension(:), allocatable :: work
        complex(dp), dimension(:,:), pointer :: A00,A01,A10,A11,A0N,AN0,A1N,AN1,ANN
        !        
        allocate( L_inv_temp(diag_blocksize, diag_blocksize) )
        allocate( U_inv_temp(diag_blocksize, diag_blocksize) )
        !
        nn = diag_blocksize
        allocate(work(nn*nn))
        allocate(ipiv(nn))
        !
        do i = 1, n_diag_blocks - 1
            h = i * diag_blocksize
            l = (i-1) * diag_blocksize + 1
            !
            A00 => A_diagonal_blocks(:, l : h)
            call zgetrf(nn, nn, A00, nn, ipiv, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetrf failed, info=', info
                call abort
            endif
            ipiv_diagonal(:,i) = ipiv
            !Compute lower factors
            U_inv_temp = czero
            call zlacpy('U',nn,nn,A00,nn,U_inv_temp,nn)            
            !
            call zgetri(nn, U_inv_temp, nn, ipiv, work, nn*nn, info)            
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info    
                call abort            
            end if
            ! L_{i+1, i} = A_{i+1, i} @ U{i, i}^{-1}
            A10 => A_lower_diagonal_blocks(:, l : h)
            A10 = matmul(A10 , U_inv_temp)
            ! L_{ndb+1, i} = A_{ndb+1, i} @ U{i, i}^{-1}
            AN0 => A_arrow_bottom_blocks(:, l : h)
            AN0 = matmul(AN0 , U_inv_temp)
            !             
            L_inv_temp = czero
            call zlacpy('L',nn,nn,A00,nn,L_inv_temp,nn)    
            do j=1,nn 
                L_inv_temp(j,j) = cone
            enddo
            !      
            call zgetri(nn, L_inv_temp, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info   
                call abort             
            end if
            ! U_{i, i+1} = L{i, i}^{-1} @ A_{i, i+1}
            A01 => A_upper_diagonal_blocks(:, l : h)
            A01 = matmul(L_inv_temp, A01)
            ! U_{i, ndb+1} = L{i, i}^{-1} @ A_{i, ndb+1}
            A0N => A_arrow_right_blocks(l : h, :) 
            A0N = matmul(L_inv_temp, A0N)            
            ! Update next diagonal block
            ! A_{i+1, i+1} = A_{i+1, i+1} - L_{i+1, i} @ U_{i, i+1}
            A11 => A_diagonal_blocks(:, (l+diag_blocksize) : (h+diag_blocksize)) 
            A11 = A11 - matmul( A10 , A01 )            
            ! Update next upper/lower blocks of the arrowhead
            ! A_{ndb+1, i+1} = A_{ndb+1, i+1} - L_{ndb+1, i} @ U_{i, i+1}
            AN1 => A_arrow_bottom_blocks(:, (l+diag_blocksize) : (h+diag_blocksize))
            AN1 = AN1 - matmul( AN0, A01 )
            ! A_{i+1, ndb+1} = A_{i+1, ndb+1} - L_{i+1, i} @ U_{i, ndb+1}
            A1N => A_arrow_right_blocks((l+diag_blocksize) : (h+diag_blocksize) , :)
            A1N = A1N - matmul( A10, A0N )
            ! Update the block at the tip of the arrowhead
            ! A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, i} @ U_{i, ndb+1}
            A_arrow_tip_block = A_arrow_tip_block - matmul( AN0, A0N )
        enddo 
        ! L_{ndb, ndb}, U_{ndb, ndb} = lu_dcmp(A_{ndb, ndb})
        h = n_diag_blocks * diag_blocksize
        l = (n_diag_blocks-1) * diag_blocksize + 1
        !
        A00 => A_diagonal_blocks(:, l : h)
        call zgetrf(nn, nn, A00, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            call abort
        endif
        ipiv_diagonal(:,n_diag_blocks) = ipiv
        !
        U_inv_temp = czero
        call zlacpy('U',nn,nn,A00,nn,U_inv_temp,nn)
        !
        call zgetri(nn, U_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info       
            call abort         
        end if
        ! L_{ndb+1, ndb} = A_{ndb+1, ndb} @ U_{ndb, ndb}^{-1}
        AN0 => A_arrow_bottom_blocks(:, l : h)
        AN0 = matmul(AN0 , U_inv_temp)        
        !
        L_inv_temp = czero
        call zlacpy('L',nn,nn,A00,nn,L_inv_temp,nn)
        do j=1,nn 
            L_inv_temp(j,j) = cone
        enddo
        !
        call zgetri(nn, L_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info     
            call abort           
        end if
        ! U_{ndb, ndb+1} = L_{ndb, ndb}^{-1} @ A_{ndb, ndb+1}
        A0N => A_arrow_right_blocks(l : h, :)
        A0N = matmul(L_inv_temp, A0N)
        ! A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, ndb} @ U_{ndb, ndb+1}
        ANN => A_arrow_tip_block
        ANN = ANN - matmul( AN0 , A0N )
        ! L_{ndb+1, ndb+1}, U_{ndb+1, ndb+1} = lu_dcmp(A_{ndb+1, ndb+1})
        nn = arrow_blocksize
        deallocate(work , ipiv)
        allocate(work(nn*nn))
        allocate(ipiv(nn))
        !
        call zgetrf(nn, nn, ANN, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            call abort
        endif
        ipiv_arrow_tip(:) = ipiv        
    end subroutine zbtatrf
    !
    ! computes the inverse of a complex block-tridiag-arrowhead (BTA) matrix 
    ! in place using the LU factorization computed by [[zbtatrf]]
    subroutine zbtatri( diag_blocksize, arrow_blocksize, n_diag_blocks, &
        A_diagonal_blocks,A_lower_diagonal_blocks,A_upper_diagonal_blocks, &
        A_arrow_bottom_blocks,A_arrow_right_blocks,A_arrow_tip_block, &
        ipiv_diagonal,ipiv_arrow_tip)
        ! in 
        integer,intent(in) :: diag_blocksize, arrow_blocksize, n_diag_blocks
        integer,intent(in),dimension(diag_blocksize, n_diag_blocks) :: ipiv_diagonal
        integer,intent(in),dimension(arrow_blocksize)::ipiv_arrow_tip
        ! in & out
        complex(dp),intent(inout),dimension(:,:),target :: A_diagonal_blocks,A_lower_diagonal_blocks,A_upper_diagonal_blocks
        complex(dp),intent(inout),dimension(:,:),target :: A_arrow_bottom_blocks,A_arrow_right_blocks,A_arrow_tip_block
        !        
        ! ---- local
        integer :: i,h,l,j
        complex(dp),allocatable,dimension(:,:) :: L_inv_temp, U_inv_temp, L10,LN0,U0N,U01,LN1,U1N
        integer :: info, nn      
        complex(dp), dimension(:), allocatable :: work
        integer, dimension(:), allocatable :: ipiv  
        complex(dp), dimension(:,:), pointer :: A10,A00,A01,A11,A1N,AN1,A0N,AN0,ANN    
        !
        nn = arrow_blocksize
        allocate(work(nn*nn))
        allocate(ipiv(nn))        
        ipiv = ipiv_arrow_tip
        ANN => A_arrow_tip_block
        call zgetri(nn, ANN, nn, ipiv, work, nn*nn, info)
        !
        nn = diag_blocksize
        deallocate(work,ipiv)
        allocate(work(nn*nn))
        allocate(ipiv(nn))        
        allocate(L_inv_temp(nn,nn))
        allocate(U_inv_temp(nn,nn))
        !
        l = diag_blocksize * (n_diag_blocks-1) + 1
        h = diag_blocksize * n_diag_blocks
        A11 => A_diagonal_blocks(:, l:h)
        ipiv = ipiv_diagonal(:,n_diag_blocks)
        U_inv_temp = czero
        L_inv_temp = czero
        call zlacpy('U',nn,nn,A11,nn, U_inv_temp ,nn)
        call zlacpy('L',nn,nn,A11,nn, L_inv_temp ,nn)
        do j=1,nn 
            L_inv_temp(j,j) = cone
        enddo
        call zgetri(nn, U_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info                
            call abort
        end if
        call zgetri(nn, L_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info                
            call abort
        end if
        allocate(LN1(diag_blocksize,diag_blocksize))
        allocate(U1N(diag_blocksize,diag_blocksize))
        AN1 => A_arrow_bottom_blocks(:, l:h)
        LN1 = AN1
        AN1 = - matmul(matmul(ANN , AN1), L_inv_temp)
        A1N => A_arrow_right_blocks(l:h, :)
        U1N = A1N
        A1N = - matmul(matmul(U_inv_temp, A1N), ANN)
        A11 => A_diagonal_blocks(:, l:h)
        A11 = matmul( U_inv_temp - &
                    matmul( A1N , LN1), &
                    L_inv_temp)
        allocate(U01(diag_blocksize,diag_blocksize))
        allocate(L10(diag_blocksize,diag_blocksize))
        allocate(LN0(arrow_blocksize,diag_blocksize))
        allocate(U0N(diag_blocksize,diag_blocksize))
        do i = n_diag_blocks-1, 1, -1
            l = diag_blocksize * (i-1) + 1
            h = diag_blocksize * i
            A00 => A_diagonal_blocks(:, l:h)
            ipiv = ipiv_diagonal(:,i)
            U_inv_temp = czero 
            L_inv_temp = czero
            call zlacpy('U',nn,nn,A00,nn, U_inv_temp ,nn)
            call zlacpy('L',nn,nn,A00,nn, L_inv_temp ,nn)
            do j=1,nn 
                L_inv_temp(j,j) = cone
            enddo
            call zgetri(nn, U_inv_temp, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info                
                call abort
            end if
            call zgetri(nn, L_inv_temp, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info                
                call abort
            end if
            ! off-diag part
            A10 => A_lower_diagonal_blocks(:, l:h)
            L10 = A10
            A11 => A_diagonal_blocks(:, (l+diag_blocksize):(h+diag_blocksize))
            A1N => A_arrow_right_blocks((l+diag_blocksize):(h+diag_blocksize) , :)
            AN0 => A_arrow_bottom_blocks(:, l:h)
            LN0 = AN0
            A10 = - matmul( A11 , L10)
            A10 = A10 - matmul( A1N , LN0 )
            A10 = matmul( A10 , L_inv_temp)
            !
            A01 => A_upper_diagonal_blocks(:, l:h)
            U01 = A01
            A0N => A_arrow_right_blocks(l:h, :)         
            U0N = A0N   
            AN1 => A_arrow_bottom_blocks(:, (l+diag_blocksize):(h+diag_blocksize))
            A01 = - matmul( U01 , A11 )            
            A01 = A01 - matmul( U0N, AN1 )
            A01 = matmul(U_inv_temp , A01)
            ! arrow part
            AN0 = - matmul( AN1 , L10 ) - matmul( ANN , LN0 )
            AN0 = matmul( AN0 , L_inv_temp)
            !
            A0N = - matmul( U01 , A1N ) - matmul( U0N , ANN)
            A0N = matmul( U_inv_temp , A0N )
            ! diag part
            A00 = U_inv_temp - matmul( A01,L10 ) - matmul( A0N, LN0 ) 
            A00 = matmul( A00,L_inv_temp )
        enddo                     
    !
    end subroutine zbtatri
    !
    !
    ! computes the selected inverse of a complex block-tridiag-arrowhead (BTA) matrix 
    ! in place without explicitly returning the LU and pivot.
    ! This function is not a simple wrapper of [[zbtatrf]] and [[zbtatri]], but a different
    ! implementation using the Schur complement formulation.
    subroutine zbtasinv( diag_blocksize, arrow_blocksize, n_diag_blocks, &
        A_diagonal_blocks,A_lower_diagonal_blocks,A_upper_diagonal_blocks, &
        A_arrow_bottom_blocks,A_arrow_right_blocks,A_arrow_tip_block )
        ! in
        integer,intent(in) :: diag_blocksize, arrow_blocksize, n_diag_blocks
        ! in & out 
        complex(dp),intent(inout),dimension(:,:),target :: A_diagonal_blocks,A_lower_diagonal_blocks,A_upper_diagonal_blocks
        complex(dp),intent(inout),dimension(:,:),target :: A_arrow_bottom_blocks,A_arrow_right_blocks,A_arrow_tip_block
        ! ---- local
        integer :: i,h,l,j
        complex(dp),allocatable,dimension(:,:) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
        integer :: info, nn
        complex(dp),dimension(:,:),pointer :: A00,A01,A10,A11,A0N,AN0,A1N,AN1,ANN
        !
        nn = diag_blocksize        
        !
        ! forward pass
        h = diag_blocksize
        l = 1
        A00 => A_diagonal_blocks(:, l : h)
        call invert_inplace( A00, nn )
        ANN => A_arrow_tip_block
        allocate(tmp1(nn,nn))
        allocate(tmp2(nn,nn))
        allocate(tmp3(arrow_blocksize,nn))
        allocate(tmp4(arrow_blocksize,arrow_blocksize))
        allocate(tmp5(arrow_blocksize,nn))
        allocate(tmp6(nn,arrow_blocksize))
        do i = 2, n_diag_blocks            
            h = h + diag_blocksize
            l = l + diag_blocksize
            A00 => A_diagonal_blocks(:, l-diag_blocksize : h-diag_blocksize)
            A11 => A_diagonal_blocks(:, l : h)
            A10 => A_lower_diagonal_blocks(:, l-diag_blocksize : h-diag_blocksize)
            A01 => A_upper_diagonal_blocks(:, l-diag_blocksize : h-diag_blocksize)
            AN1 => A_arrow_bottom_blocks(: , l : h)
            A1N => A_arrow_right_blocks(l : h , :)
            AN0 => A_arrow_bottom_blocks(: , l-diag_blocksize : h-diag_blocksize)
            A0N => A_arrow_right_blocks(l-diag_blocksize : h-diag_blocksize , :)
            tmp1 = matmul(A10, A00)
            tmp2 = matmul(tmp1, A01)
            A11 = A11 - tmp2
            call invert_inplace(A11, nn)
            tmp3 = matmul(AN0,A00)
            tmp5 = matmul(tmp3, A01)
            AN1 = AN1 - tmp5
            tmp6 = matmul(tmp1, A0N)
            A1N = A1N - tmp6
            tmp4 = matmul(tmp3, A0N)
            ANN = ANN - tmp4           
        enddo 
        tmp3 = matmul(AN1, A11)
        tmp4 = matmul(tmp3, A1N)
        ANN = ANN - tmp4
        call invert_inplace( ANN , arrow_blocksize )
        !
        ! backward pass        
        tmp6 = matmul(A11 , A1N)
        tmp3 = matmul(ANN , AN1)
        tmp2 = matmul(tmp6,tmp3)
        tmp1 = matmul(tmp2, A11) ! A11 A1N ANN AN1 A11        
        AN1 = - matmul(tmp3, A11) 
        A1N = - matmul(tmp6, ANN)
        A11 = A11 + tmp1   
        !           
        deallocate(tmp5)
        allocate(tmp5(nn,nn))
        !
        do i = n_diag_blocks - 1 , 1 , -1 
            h = h - diag_blocksize 
            l = l - diag_blocksize 
            A00 => A_diagonal_blocks(:, l : h)
            A11 => A_diagonal_blocks(:, l+diag_blocksize : h+diag_blocksize)
            A10 => A_lower_diagonal_blocks(:, l : h)
            A01 => A_upper_diagonal_blocks(:, l : h)
            AN0 => A_arrow_bottom_blocks(: , l : h)
            AN1 => A_arrow_bottom_blocks(: , l+diag_blocksize : h+diag_blocksize)
            A0N => A_arrow_right_blocks(l : h , :)
            A1N => A_arrow_right_blocks(l+diag_blocksize : h+diag_blocksize , :)
            !
            tmp1 = matmul(A01,A11)  + matmul(A0N,AN1)
            tmp6 = matmul(A01,A1N)  + matmul(A0N,ANN)
            tmp5 = matmul(tmp1,A10) + matmul(tmp6,AN0)
            tmp3 = matmul(AN1,A10) + matmul(ANN,AN0)
            tmp2 = matmul(A11,A10) + matmul(A1N,AN0)
            !
            A01 = - matmul(A00,tmp1)                 !!! <- A01 replaced by X01
            A0N = - matmul(A00,tmp6)                 !!! <- A0N replaced by X0N
            !            
            A10 = - matmul(tmp2,A00)                 !!! <- A10 replaced by X10
            AN0 = - matmul(tmp3,A00)                 !!! <- AN0 replaced by XN0
            !
            tmp1 = matmul(A00,tmp5)
            tmp2 = matmul(tmp1,A00)
            !
            A00 = A00 + tmp2
        enddo 
    end subroutine zbtasinv


    ! matrix inversion
    subroutine invert_inplace(A, nn)
        integer :: info, nn
        integer, dimension(:), allocatable :: ipiv
        complex(dp), dimension(nn, nn), intent(inout) :: A
        complex(dp), dimension(:), allocatable :: work
        allocate(work(nn*nn))
        allocate(ipiv(nn))
        call zgetrf(nn, nn, A, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            call abort
        endif
        call zgetri(nn, A, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info
            call abort
        end if
    end subroutine invert_inplace

end module sinv