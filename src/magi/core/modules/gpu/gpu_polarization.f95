
module gpu_polarization_m
    ! CUDA functions to compute the polarization integral
    use cublas
    implicit none
    integer, parameter :: dp=8
    contains
    !
    ! calculates independent-particle polarizability matrix at one frequency
    subroutine gpu_four_polarization_block(alpha,nm_dev,nen,dE,nop,block_size1,&
        block_size2,G_lesser,G_greater,G_retarded,jl,ki,L0)
        integer,intent(in) :: nm_dev,nen,nop,block_size1,block_size2
        integer,intent(in) :: jl(block_size1),ki(block_size2)
        real(dp),intent(in) :: dE, alpha 
        complex(dp),intent(in),dimension(nm_dev*nm_dev,nen) :: G_lesser,G_greater,G_retarded
        complex(dp),intent(out),dimension(block_size1,block_size2) :: L0
        ! ---
        real(dp) :: weights
        integer :: n
        complex(dp),dimension(:,:),allocatable :: a,b,c
        REAL(kind=dp), PARAMETER :: twopi = 6.2831853072_dp
        COMPLEX(kind=dp), PARAMETER :: czero  = dcmplx(0.0_dp,0.0_dp)
        COMPLEX(kind=dp), PARAMETER :: cone = dcmplx(1.0_dp,0.0_dp)
        complex(dp)::a1,a2
        ! the P4 tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                           
        weights = dE/twopi
        a1 = weights * alpha * 0.5_dp
        !
        n = nen-nop
        allocate(a(block_size1,n))
        allocate(b(n,block_size2))
        allocate(c(block_size1,block_size2))
        a = G_greater(jl, (nop+1):nen)         
        b = reshape( G_lesser(ki,1:n) , shape=[n,block_size2], order=[2,1] )
        call zgemm('n','n',block_size1,block_size2,n,a1,a,block_size1,b,n,czero,c,block_size1) 
        a = G_lesser(jl, (nop+1):nen)         
        b = reshape( G_greater(ki,1:n) , shape=[n,block_size2], order=[2,1] )
        call zgemm('n','n',block_size1,block_size2,n,-a1,a,block_size1,b,n,cone,c,block_size1) 
        !             
        ! L0 =    (1.0_dp - alpha) * ( sum( G_lesser_jl((nop+1):nen)   * conjg(G_retarded_ik(1:(nen-nop))) ) &
        !                             +  sum( G_retarded_jl((nop+1):nen) * G_lesser_ki(1:(nen-nop)) ) )  &
        !         + alpha * 0.5_dp * ( sum( G_greater_jl((nop+1):nen) * G_lesser_ki(1:(nen-nop)) )  & 
        !                             -  sum( G_lesser_jl((nop+1):nen)  * G_greater_ki(1:(nen-nop)) ) )  
        ! L0 = L0 * weights 
    end subroutine gpu_four_polarization_block
    !
    !
    ! calculates the summation of products of two complex arrays : 
    ! c(i,j) = sum( a(i, offset+1:n) + b(j, 1:n-offset) ) / dble(n)
    subroutine gpu_sum_complex(a,b,nnz,n,c)
        integer,intent(in) :: nnz, n
        complex(8),intent(in),dimension(nnz,n) :: a
        complex(8),intent(in),dimension(n,nnz) :: b
        complex(8),intent(out) :: c(nnz,nnz)
        complex(8),parameter::czero=dcmplx(0.0d0,0.0d0)
        complex(8),parameter::cone=dcmplx(1.0d0,0.0d0)
        ! ---       
        integer :: i,j,k
        !$omp target enter data map(to:a,b,c)  
        !$omp target data use_device_ptr(a,b,c)
        call zgemm('n','n',nnz,nnz,n,cone,a,nnz,b,n,czero,c,nnz) 
        !$omp end target data 
        !$omp target update from(c)     
        !$omp target exit data map(delete:a,b,c)
        !
    end subroutine gpu_sum_complex
    !
end module gpu_polarization_m