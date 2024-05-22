
module block_sparse 
    use parameters_mod
    implicit none 
    contains 
    !    
    ! return a dense block matrix of size (block_size x block_size) filled with values from 
    !   array `v` , the `col_index` is like CSR index array but with column index within block matrix
    !   the `ind_ptr` is like a CSR `ind_ptr` but for each block
    !   the `iblock` is block index, `idiag` is off-diagonal index of the wanted block 
    subroutine get_block_from_bcsr(v,col_index,ind_ptr,nnz,block_size,num_blocks,num_diag,iblock,idiag,mat)
        complex(dp),intent(in) :: v(nnz)
        integer,intent(in) :: col_index(nnz),nnz,block_size,num_blocks,num_diag,iblock,idiag
        integer,intent(in) :: ind_ptr(block_size+1,num_blocks,num_diag)
        complex(dp),intent(out) :: mat(block_size,block_size)
        ! ------
        integer :: i,j,ptr1,ptr2
        mat = czero
        do i=1,block_size
            ! get ind_ptr for the block row i
            ptr1 = ind_ptr(i,  iblock, idiag)
            ptr2 = ind_ptr(i+1,iblock, idiag)
            do j=ptr1,ptr2-1
                mat(i, col_index(j)) = v(j)
            enddo
        enddo
    end subroutine get_block_from_bcsr
    !
    ! put a dense matrix values into the corresponding position of value array `v` , similar to `get_block_from_bcsr`    
    subroutine put_block_to_bcsr(v,col_index,ind_ptr,nnz,block_size,num_blocks,num_diag,iblock,idiag,mat)
        complex(dp),intent(out) :: v(nnz)
        integer,intent(in) :: col_index(nnz),nnz,block_size,num_blocks,num_diag,iblock,idiag
        integer,intent(in) :: ind_ptr(block_size+1,num_blocks,num_diag)
        complex(dp),intent(in) :: mat(block_size,block_size)
        ! ------
        integer :: i,j,ptr1,ptr2        
        do i=1,block_size
            ! get ind_ptr for the block row i
            ptr1 = ind_ptr(i,  iblock, idiag)
            ptr2 = ind_ptr(i+1,iblock, idiag)
            do j=ptr1,ptr2-1
                v(j) = mat(i, col_index(j))
            enddo
        enddo
    end subroutine put_block_to_bcsr
    
    
end module block_sparse



