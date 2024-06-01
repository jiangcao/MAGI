
module block_sparse 
    use parameters_mod
    implicit none 
    contains 
    !    
    ! Modified Blocked CSR format (bcsr) 
    ! ----------------------------------
    !
    ! in the data array v(:) the first N elements are the main diagonal of the
    ! matrix, followed by simply stacking CSR of each block (without the main
    ! diagonal) continuously, with an easy access to the i-th block on i-th
    ! diagonal. The sizes of diagonal blocks are defined by `block_sizes`. The
    ! diagonal blocks have to be square. The `ind_ptr(:,iblock,idiag)` is the 
    ! CSR row pointer of the corresponding block, pointing to the
    ! starting position of each row in the data array `v`. The `col_index` is
    ! the CSR column index (within each block, not the column of entire matrix). 
    !
    ! return a dense block matrix of size (block_size x block_size) filled with values from 
    !   array `v` , the `col_index` is like CSR index array but with column index within block matrix
    !   the `ind_ptr` is like a CSR `ind_ptr` but for each block
    !   the `iblock` is block index, `idiag` is off-diagonal index of the wanted block 
    !
    !   idiag = 0 : main diagonal blocks
    !   0 < idiag <= num_diag : upper off diagonal blocks
    !   num_diag < idiag <= 2*num_diag : lower off diagonal blocks
    subroutine get_block_from_bcsr(v,col_index,ind_ptr,block_size,block_start_index,num_blocks,num_diag,iblock,idiag,mat)
        complex(dp),intent(in) :: v(:)
        integer,intent(in) :: col_index(:),block_size,num_blocks,num_diag,block_start_index(:)
        integer,intent(in) :: iblock,idiag
        integer,intent(in) :: ind_ptr(:,:,:)
        complex(dp),intent(out) :: mat(block_size,block_size)
        ! ------
        integer :: i,j,ptr1,ptr2
        mat = czero
        do i=1,block_size
            ! get ind_ptr for the block row i
            ptr1 = ind_ptr(i,  iblock, idiag+1)
            ptr2 = ind_ptr(i+1,iblock, idiag+1)
            mat(i, col_index(ptr1:ptr2-1)) = v(ptr1:ptr2-1)
        enddo
        if (idiag == 0) then
            ! main diagonal block
            do i=1,block_size
                mat(i,i) = v( block_start_index(iblock) + i )
            enddo
        endif
    end subroutine get_block_from_bcsr
    !
    ! put a dense matrix values into the corresponding position of value array `v` , similar to `get_block_from_bcsr`    
    subroutine put_block_to_bcsr(v,col_index,ind_ptr,block_size,block_start_index,num_blocks,num_diag,iblock,idiag,mat)
        complex(dp),intent(inout) :: v(:)
        integer,intent(in) :: col_index(:),block_size,num_blocks,num_diag
        integer,intent(in) :: iblock,idiag,block_start_index(:)
        integer,intent(in) :: ind_ptr(:,:,:)
        complex(dp),intent(in) :: mat(block_size,block_size)
        ! ------
        integer :: i,j,ptr1,ptr2        
        do i=1,block_size
            ! get ind_ptr for the block row i
            ptr1 = ind_ptr(i,  iblock, idiag+1)
            ptr2 = ind_ptr(i+1,iblock, idiag+1)
            v(ptr1:ptr2-1) = mat(i, col_index(ptr1:ptr2-1))
        enddo
        if (idiag == 0) then
            ! main diagonal block
            do i=1,block_size
                v( block_start_index(iblock) + i ) = mat(i,i)  
            enddo
        endif
    end subroutine put_block_to_bcsr
    
    
end module block_sparse



