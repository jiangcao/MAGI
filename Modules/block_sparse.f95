
module block_sparse 
    use parameters_mod
    implicit none 
    contains 
    !    
    ! Modified Blocked CSR format (bcsr) 
    ! ----------------------------------
    !
    ! Note that the first N elements in the data array `v(:)` are the main diagonal 
    ! of matrix, followed by simply stacking CSR of each block continuously (exculding 
    ! the main diagonal elements! Otherwise, those elements are simply  
    ! overwritten by `v(1:N)`)
    ! The sizes of diagonal blocks are defined by `block_sizes`. The
    ! diagonal blocks have to be square. 
    ! The data pointer `ind_ptr` allows an easy and direct access to the i-th block 
    ! on i-th diagonal. 
    ! The `ind_ptr(row:row+1,iblock,idiag)` is the CSR row pointer of the corresponding block, 
    ! pointing to the starting position in the data array `v` of `row` and `row+1`. The `col_index` 
    ! is the CSR column index array (just within the block, NOT the global column index of the entire sparse matrix). 
    !
    ! Returns a dense block matrix of size (block_size x block_size) filled with values from 
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
    ! Puts a dense matrix values into the corresponding position of value array `v` , similar to `get_block_from_bcsr`    
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



