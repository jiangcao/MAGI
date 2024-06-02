module sinv
    use parameters_mod,only:dp,c1i,czero,cone    
    use omp_lib
    
    implicit none
    contains

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
        complex(dp), dimension(:,:), pointer :: A
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
            A => A_diagonal_blocks(:, l : h)
            ! L_{i, i}, U_{i, i} = lu_dcmp(A_{i, i})
            call zgetrf(nn, nn, A, nn, ipiv, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetrf failed, info=', info
                call abort
            endif
            ipiv_diagonal(:,i) = ipiv
            !Compute lower factors
            call zlacpy('U',nn,nn,A,nn,U_inv_temp,nn)            
            !
            call zgetri(nn, U_inv_temp, nn, ipiv, work, nn*nn, info)            
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info    
                call abort            
            end if
            ! L_{i+1, i} = A_{i+1, i} @ U{i, i}^{-1}
            A => A_lower_diagonal_blocks(:, l : h)
            A = matmul(A , U_inv_temp)
            ! L_{ndb+1, i} = A_{ndb+1, i} @ U{i, i}^{-1}
            A => A_arrow_bottom_blocks(:, l : h)
            A = matmul(A , U_inv_temp)
            ! Compute upper factors
            A => A_diagonal_blocks(:, l : h)
            call zlacpy('L',nn,nn,A,nn,L_inv_temp,nn)      
            !      
            call zgetri(nn, L_inv_temp, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info   
                call abort             
            end if
            ! U_{i, i+1} = L{i, i}^{-1} @ A_{i, i+1}
            A => A_upper_diagonal_blocks(:, l : h)
            A = matmul(L_inv_temp, A)
            ! U_{i, ndb+1} = L{i, i}^{-1} @ A_{i, ndb+1}
            A => A_arrow_right_blocks(l : h, :) 
            A = matmul(L_inv_temp, A)            
            ! Update next diagonal block
            ! A_{i+1, i+1} = A_{i+1, i+1} - L_{i+1, i} @ U_{i, i+1}
            A => A_diagonal_blocks(:, (l+diag_blocksize) : (h+diag_blocksize)) 
            A = A - matmul( A_lower_diagonal_blocks(:, l:h) , A_upper_diagonal_blocks(:,l:h) )            
            ! Update next upper/lower blocks of the arrowhead
            ! A_{ndb+1, i+1} = A_{ndb+1, i+1} - L_{ndb+1, i} @ U_{i, i+1}
            A => A_arrow_bottom_blocks(:, (l+diag_blocksize) : (h+diag_blocksize))
            A = A - matmul( A_arrow_bottom_blocks(:, l:h), A_upper_diagonal_blocks(:, l:h) )
            ! A_{i+1, ndb+1} = A_{i+1, ndb+1} - L_{i+1, i} @ U_{i, ndb+1}
            A => A_arrow_right_blocks((l+diag_blocksize) : (h+diag_blocksize) , :)
            A = A - matmul( A_lower_diagonal_blocks(:, l:h), A_arrow_right_blocks(l:h,:) )
            ! Update the block at the tip of the arrowhead
            ! A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, i} @ U_{i, ndb+1}
            A_arrow_tip_block = A_arrow_tip_block - matmul(A_arrow_bottom_blocks(:,l:h), A_arrow_right_blocks(l:h,:))
        enddo 
        ! L_{ndb, ndb}, U_{ndb, ndb} = lu_dcmp(A_{ndb, ndb})
        h = n_diag_blocks * diag_blocksize
        l = (n_diag_blocks-1) * diag_blocksize + 1
        !
        A => A_diagonal_blocks(:, l : h)
        call zgetrf(nn, nn, A, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            call abort
        endif
        ipiv_diagonal(:,n_diag_blocks) = ipiv
        ! L_{ndb+1, ndb} = A_{ndb+1, ndb} @ U_{ndb, ndb}^{-1}
        call zlacpy('U',nn,nn,A,nn,U_inv_temp,nn)
        !
        call zgetri(nn, U_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info       
            call abort         
        end if
        A => A_arrow_bottom_blocks(:, l : h)
        A = matmul(A , U_inv_temp)
        ! U_{ndb, ndb+1} = L_{ndb, ndb}^{-1} @ A_{ndb, ndb+1}
        A => A_diagonal_blocks(:, l : h)
        call zlacpy('L',nn,nn,A,nn,L_inv_temp,nn)
        !
        call zgetri(nn, L_inv_temp, nn, ipiv, work, nn*nn, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetri failed, info=', info     
            call abort           
        end if
        A => A_arrow_right_blocks(l : h, :)
        A = matmul(L_inv_temp, A)
        ! A_{ndb+1, ndb+1} = A_{ndb+1, ndb+1} - L_{ndb+1, ndb} @ U_{ndb, ndb+1}
        A => A_arrow_tip_block
        A = A - matmul( A_arrow_bottom_blocks(:, l:h) , A_arrow_right_blocks(l:h, :) )
        ! L_{ndb+1, ndb+1}, U_{ndb+1, ndb+1} = lu_dcmp(A_{ndb+1, ndb+1})
        nn = arrow_blocksize
        deallocate(work , ipiv)
        allocate(work(nn*nn))
        allocate(ipiv(nn))
        !
        call zgetrf(nn, nn, A, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            call abort
        endif
        ipiv_arrow_tip(:) = ipiv        
    end subroutine zbtatrf

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
        complex(dp),allocatable,dimension(:,:) :: L_inv_temp, U_inv_temp
        integer :: info, nn      
        complex(dp), dimension(:), allocatable :: work
        complex(dp), dimension(:,:), pointer :: A    
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
        A => A_arrow_bottom_blocks(:, l:h)
        A = - matmul(matmul(A_arrow_tip_block , A), L_inv_temp)
        A => A_arrow_right_blocks(l:h, :)
        A = - matmul(matmul(U_inv_temp, A), A_arrow_tip_block)
        A => A_diagonal_blocks(:, l;h)
        A = matmul( U_inv_temp - &
                    matmul( A_arrow_right_blocks(l:h, :) , A_arrow_bottom_blocks(:, l:h)), &
                    L_inv_temp)
        do i = n_diag_blocks-1, 1, -1
            l = diag_blocksize * (i-1) + 1
            h = diag_blocksize * i
            A => A_diagonal_blocks(:, l:h)
            ipiv = ipiv_diagonal(:,i)
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
            A => A_lower_diagonal_blocks(:, l:h)
            A = matmul( - A_diagonal_blocks(:, (l+diag_blocksize):(h+diag_blocksize)) , A)
            A = A - matmul( A_arrow_right_blocks((l+diag_blocksize):(h+diag_blocksize) , :) , A_arrow_bottom_blocks(:, l:h) )
            A = matmul( A , L_inv_temp)
            !
            A => A_upper_diagonal_blocks(:, l:h)
            

        enddo                     

    end subroutine zbtatri

end module sinv