
module gpu_polarization_m
#ifdef _CUDA
    ! custom CUDA kernel to compute the polarization integral
    use parameters_mod,only:dp,twopi 
    use cudafor
    implicit none
    contains
    !
    ! calculate summation
    subroutine polarization_sum(nen,M,dE,nop,x,y,res)
       ! in 
       ! Scalar arguments from the host should be passed by value
       integer,intent(in),value :: nen,nop,M
       ! Array arguments to kernels are automatically assumed to be on the device
       real(dp),intent(in) :: dE
       complex(dp),intent(in),dimension(nen,M) :: x,y
       ! out 
       complex(dp),intent(out),value :: res(M)
       ! --- local
       complex(dp)::mysum
       real(dp) :: weights
       integer :: tx, ty, ie, je
       ! the P4 tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                               
       weights = dE/twopi
       ! Threads in CUDA are organized in blocks with multiple blocks forming a larger grid.
       ! threadIdx is the index of a given thread within a block (x-dimension is leading 
       ! dimension). In CUDA Fortran, these indices are 1-based. 
       tx = threadIdx%x
       ty = threadIdx%y
       ! blockIdx is the index of a given block within the grid, and blockDim is the dimension 
       ! of the block (in number of threads). Here, we are assigning the 'i'th integral a 
       ! full x-row (a warp if blockDim%x == 32) of threads
       ie = (blockIdx%x - 1) * blockDim%y + ty
       ! We may have extra threads assigned beyond what we require. We can force a quick return
       ! for those threads. 
       if (ie > M) return
       ! In the following bit of code, each thread accumulates a local sum of several values of
       ! the 'i'th integral it is assigned to.  
       mysum = 0.d0
  
       do je = tx, nen, blockDim%x
          mysum = mysum + x(je,ie) * y(je,ie)
       end do

       ! Now, we need to reduce the values across threads assigned to the 'i'th integral.
       ! **Assuming** a full warp (32 threads) was assigned to each integral, we can
       ! complete this task using "shuffle" instructions. Alternative logic is required
       ! if multiple warps (more than 32 threads) within a block are assigned to a single 
       ! integral (using a shared memory buffer). 
       mysum = mysum + __shfl_down(mysum,1)
       mysum = mysum + __shfl_down(mysum,2)
       mysum = mysum + __shfl_down(mysum,4)
       mysum = mysum + __shfl_down(mysum,8)
       mysum = mysum + __shfl_down(mysum,16)
  
       ! First thread in each warp will contain final reduced result. It is tasked with
       ! assigning final result to res array.
       if (tx == 1) then
         res(ie) = weights * mysum 
       endif
    end subroutine polarization_sum
    !
#endif
end module gpu_polarization_m