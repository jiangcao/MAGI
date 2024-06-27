!===============================================================================
! Copyright (C) 2023 Jiang Cao
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: Jiang Cao <jiacao@ethz.ch>
! Comment:
!  
! Maintenance:
!===============================================================================
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
    subroutine get_block_from_bcsr(v,col_index,ind_ptr,block_nrow,block_ncol,block_start_index,iblock,idiag,mat)
        complex(dp),intent(in) :: v(:)
        integer,intent(in) :: col_index(:),block_nrow,block_ncol
        integer,intent(in) :: iblock,idiag,block_start_index(:)
        integer,intent(in) :: ind_ptr(:,:,:)
        complex(dp),intent(out) :: mat(block_nrow,block_ncol)
        ! ------
        integer :: i,j,ptr1,ptr2
        mat = czero
        do i=1,block_nrow
            ! get ind_ptr for the block row i
            ptr1 = ind_ptr(i,  iblock, idiag+1)
            ptr2 = ind_ptr(i+1,iblock, idiag+1)
            mat(i, col_index(ptr1:ptr2-1)) = v(ptr1:ptr2-1)
        enddo
        if (idiag == 0) then
            ! main diagonal block
            do i=1,block_nrow
                mat(i,i) = v( block_start_index(iblock) + i )
            enddo
        endif
    end subroutine get_block_from_bcsr
    !
    ! Puts a dense matrix values into the corresponding position of value array `v` , similar to `get_block_from_bcsr`    
    subroutine put_block_to_bcsr(v,col_index,ind_ptr,block_size,block_start_index,iblock,idiag,mat)
        complex(dp),intent(inout) :: v(:)
        integer,intent(in) :: col_index(:),block_size
        integer,intent(in) :: iblock,idiag,block_start_index(:)
        integer,intent(in) :: ind_ptr(:,:,:)
        complex(dp),intent(in) :: mat(:,:)
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
    !    
    !
    !
    ! Modified Blocked Sparse Vector format (bmsv) 
    ! --------------------------------------------
    !
    ! The first N elements in the data array `v` are the diagonal elements in the matrix.
    ! Then, the remaining elements in the matrix is arranged by blocks. Each block is 
    ! stored in contiguous memory with the modified Sparse Vector format. The data is flattened
    ! into a vector. The integer array `IA` registers the number of zeros between consecutive 
    ! non-zero elements. The lookup table `block_ptr[1:2,i,j]` saves the pointer to the begining position of the
    ! block and the ending positions in `v` of the block [i,j].
    !
    subroutine get_block_from_bmsv(v,IA,block_ptr,block_nrow,block_ncol,block_start_index,block_ij,mat)
        complex(dp),intent(in) :: v(:)
        integer,intent(in) :: IA(:),block_nrow,block_ncol,block_start_index(:)
        integer,intent(in) :: block_ij(:)
        integer,intent(in) :: block_ptr(:,:,:)
        complex(dp),intent(out),target :: mat(block_nrow,block_ncol)
        ! ------
        integer :: i,ib,jb,row,col
        mat = czero
        ib=block_ij(1)
        jb=block_ij(2)
        row=1
        col=0
        do i=block_ptr(1,ib,jb), block_ptr(2,ib,jb)            
            col = col+ IA(i) + 1 
            do while (col > block_ncol)
                col = col - block_ncol
                row = row + 1
            enddo
            mat(row,col) = v(i)
        enddo
        if (ib == jb) then
            ! main diagonal block
            do i = 1, block_nrow
                mat(i,i) = v( block_start_index(ib) + i ) 
            enddo
        endif
    end subroutine get_block_from_bmsv
    !
    subroutine put_block_to_bmsv(v,IA,block_ptr,block_nrow,block_ncol,block_start_index,block_ij,mat)
        complex(dp),intent(inout) :: v(:)
        integer,intent(in) :: IA(:),block_nrow,block_ncol,block_start_index(:)
        integer,intent(in) :: block_ij(:)
        integer,intent(in) :: block_ptr(:,:,:)
        complex(dp),intent(in),target :: mat(block_nrow,block_ncol)
        ! ------
        integer :: i,ib,jb,row,col
        ib=block_ij(1)
        jb=block_ij(2)
        row=1
        col=0
        do i = block_ptr(1,ib,jb), block_ptr(2,ib,jb) 
            col = col+ IA(i) + 1 
            do while (col > block_ncol)
                col = col - block_ncol
                row = row + 1
            enddo           
            v(i) = mat(row,col)
        enddo
        if (ib == jb) then
            ! main diagonal block
            do i = 1, block_nrow
                v( block_start_index(ib) + i ) = mat(i,i) 
            enddo
        endif
    end subroutine put_block_to_bmsv
    
end module block_sparse



