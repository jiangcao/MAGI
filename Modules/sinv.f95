module sinv
    use parameters_mod,only:dp,c1i,czero,cone    
    use omp_lib
    !
    implicit none
    contains
    !
    ! computes an LU factorization of a complex block-tridiag-arrowhead (BTA) matrix in place
    ! it calls LAPACK `zgetrf` and `zgetri` functions
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
        integer :: i,h,l
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
            call zlacpy('L',nn,nn,A00,nn,L_inv_temp,nn)      
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
        ! L_{ndb+1, ndb} = A_{ndb+1, ndb} @ U_{ndb, ndb}^{-1}
        call zlacpy('U',nn,nn,A00,nn,U_inv_temp,nn)
        !
        call zgetri(nn, U_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info       
            call abort         
        end if
        AN0 => A_arrow_bottom_blocks(:, l : h)
        AN0 = matmul(AN0 , U_inv_temp)
        ! U_{ndb, ndb+1} = L_{ndb, ndb}^{-1} @ A_{ndb, ndb+1}
        !
        call zlacpy('L',nn,nn,A00,nn,L_inv_temp,nn)
        !
        call zgetri(nn, L_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info     
            call abort           
        end if
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
    ! computes the inverse of a complex block-tridiag-arrowhead (BTA) matrix in place
    ! using the LU factorization computed by [zbtatrf]
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
        ! out       
        ! ---- local
        integer :: i,h,l
        complex(dp),allocatable,dimension(:,:) :: L_inv_temp, U_inv_temp, L10,LN0,U01
        integer :: info, nn      
        complex(dp), dimension(:), allocatable :: work
        integer, dimension(:), allocatable :: ipiv  
        complex(dp), dimension(:,:), pointer :: A,A10,A00,A01,A11,A1N,AN1,A0N,AN0    
        !
        nn = arrow_blocksize
        allocate(work(nn*nn))
        allocate(ipiv(nn))        
        ipiv = ipiv_arrow_tip
        A => A_arrow_tip_block
        call zgetri(nn, A, nn, ipiv, work, nn*nn, info)
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
        A => A_diagonal_blocks(:, l:h)
        ipiv = ipiv_diagonal(:,n_diag_blocks)
        U_inv_temp = czero
        L_inv_temp = czero
        call zlacpy('U',nn,nn,A,nn, U_inv_temp ,nn)
        call zlacpy('L',nn,nn,A,nn, L_inv_temp ,nn)
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
        AN1 => A_arrow_bottom_blocks(:, l:h)
        AN1 = - matmul(matmul(A_arrow_tip_block , AN1), L_inv_temp)
        A1N => A_arrow_right_blocks(l:h, :)
        A1N = - matmul(matmul(U_inv_temp, A1N), A_arrow_tip_block)
        A11 => A_diagonal_blocks(:, l:h)
        A11 = matmul( U_inv_temp - &
                    matmul( A_arrow_right_blocks(l:h, :) , A_arrow_bottom_blocks(:, l:h)), &
                    L_inv_temp)
        allocate(U01(diag_blocksize,diag_blocksize))
        allocate(L10(diag_blocksize,diag_blocksize))
        allocate(LN0(arrow_blocksize,diag_blocksize))
        do i = n_diag_blocks-1, 1, -1
            l = diag_blocksize * (i-1) + 1
            h = diag_blocksize * i
            A00 => A_diagonal_blocks(:, l:h)
            ipiv = ipiv_diagonal(:,i)
            U_inv_temp = czero 
            L_inv_temp = czero
            call zlacpy('U',nn,nn,A00,nn, U_inv_temp ,nn)
            call zlacpy('L',nn,nn,A00,nn, L_inv_temp ,nn)
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
            A10 = matmul( - A11 , A10)
            A10 = A10 - matmul( A1N , AN0 )
            A10 = matmul( A10 , L_inv_temp)
            !
            A01 => A_upper_diagonal_blocks(:, l:h)
            U01 = A01
            A01 = matmul( - A01 , A11 )
            A0N => A_arrow_right_blocks(l:h, :)
            AN1 => A_arrow_bottom_blocks(:, (l+diag_blocksize):(h+diag_blocksize))
            A01 = A01 - matmul( A0N, AN1 )
            A01 = matmul(U_inv_temp , A01)
            ! arrow part
            AN0 = matmul( -AN1 , L10 ) - matmul( A_arrow_tip_block , LN0 )
            AN0 = matmul( AN0 , L_inv_temp)
            !
            A0N = matmul( -U01 , A1N ) - matmul( A0N , A_arrow_tip_block)
            A0N = matmul( U_inv_temp , A0N )
            ! diag part
            A00 = U_inv_temp - matmul( A01,L10 ) - matmul( A0N, LN0 ) 
            A00 = matmul(A00,L_inv_temp)
        enddo                     
    !
    end subroutine zbtatri
    
end module sinv