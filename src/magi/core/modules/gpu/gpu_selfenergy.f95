!===============================================================================
! Copyright (C) 2023 Jiang Cao
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jiacao <jiacao@ethz.ch>
! Comment:
!  
! Maintenance:
!===============================================================================
module gpu_selfenergy
    implicit none 
    
    integer,parameter :: size_of_complex=16 !8->single precision; 16->double precision

    contains

    subroutine gpu_selfenergy_GW(nen,nop,nm,num_ij,ij,copy_to_gpu,G_lesser,G_greater,G_retarded,&
                                W_lesser,W_greater,W_retarded,&
                                devPtrGG,devPtrGL,devPtrGR,&
                                sig_lesser,sig_greater,sig_retarded)
        integer,intent(in) :: nop,nen,nm,num_ij,ij(:)
        integer(8),intent(inout) :: devPtrGL, devPtrGG, devPtrGR
        complex(8), dimension(:,:), intent(in) :: G_greater, G_lesser, G_retarded
        complex(8), dimension(:), intent(in) :: W_greater, W_lesser, W_retarded
        complex(8), dimension(:,:), intent(inout) :: sig_lesser,sig_greater,sig_retarded
        complex(8), dimension(:), allocatable :: sig
        logical,intent(in) :: copy_to_gpu
        integer :: n,i,j        
        integer(8) :: devPtrA, devPtrB, devPtrC
        n = nen - nop 
        
        allocate(sig(nen))
        !
        if (copy_to_gpu) then 
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGG)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGL)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGR)        
            !copy data to GPU
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_lesser,nen,devPtrGL,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_greater,nen,devPtrGG,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_retarded,nen,devPtrGR,nen)
        endif
        !
        call cublas_alloc(nen, size_of_complex, devPtrA)
        !
        do i=1,num_ij                    
            ! sig^<[E] += G^<[E-dE] W^< - G^<[E+dE] conjg(W^>)    
            call cublas_zscal(nen, dcmplx(0.0,0.0), devPtrA, 1)                        
            call cublas_zcopy(n, devPtrGL + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrA+(nop*size_of_complex), 1)
            call cublas_zscal(n, W_lesser(ij(i)), devPtrA+(nop*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_greater(ij(i))), devPtrGL + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrA, 1)
            !copy data from GPU
            call cublas_get_matrix(nen,1,size_of_complex,devPtrA,nen,sig,nen)
            Sig_lesser(1:nen,ij(i)) = Sig_lesser(1:nen,ij(i)) + sig(1:nen)
            !
            ! sig^>[E] += G^>[E-dE] W^> - G^>[E+dE] conjg(W^<)    
            call cublas_zscal(nen, dcmplx(0.0,0.0), devPtrA, 1)                        
            call cublas_zcopy(n, devPtrGG + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrA+(nop*size_of_complex), 1)
            call cublas_zscal(n, W_greater(ij(i)), devPtrA+(nop*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_lesser(ij(i))), devPtrGG + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrA, 1)
            !copy data from GPU
            call cublas_get_matrix(nen,1,size_of_complex,devPtrA,nen,sig,nen)
            Sig_greater(1:nen,ij(i)) = Sig_greater(1:nen,ij(i)) + sig(1:nen)
            !
            ! sig^r[E] += G^<[E-dE] W^r + G^r[E-dE] W^< + G^r[E-dE] W^r - G^<[E+dE] conjg(W^r) - G^r[E+dE] conjg(W^>) - G^r[E-dE] conjg(W^r)   
            call cublas_zscal(nen, dcmplx(0.0,0.0), devPtrA, 1)                        
            call cublas_zcopy(n, devPtrGL + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrA+(nop*size_of_complex), 1)
            call cublas_zscal(n, W_retarded(ij(i)), devPtrA+(nop*size_of_complex), 1)
            call cublas_zaxpy(n, W_lesser(ij(i)),   devPtrGR + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrA+(nop*size_of_complex), 1)
            call cublas_zaxpy(n, W_retarded(ij(i)), devPtrGR + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrA+(nop*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_retarded(ij(i))), devPtrGL + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrA, 1)
            call cublas_zaxpy(n, -conjg(W_greater(ij(i))),  devPtrGR + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrA, 1)
            call cublas_zaxpy(n, -conjg(W_retarded(ij(i))), devPtrGR + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrA, 1)
            !copy data from GPU
            call cublas_get_matrix(nen,1,size_of_complex,devPtrA,nen,sig,nen)
            Sig_retarded(1:nen,ij(i)) = Sig_retarded(1:nen,ij(i)) + sig(1:nen)
        enddo 
        !
        call cublas_free(devPtrA)
    end subroutine gpu_selfenergy_GW 

    ! do i=1,m*m                                
    !     ie1 = max(nop,1) + 1
    !     ie2 = min(nen+nop,nen)
    !     Sig_lesser(ie1:ie2,ij(i))=Sig_lesser(ie1:ie2,ij(i)) + G_lesser((ie1-nop):(ie2-nop),ij(i)) * W_lesser(ij(i))                                
    !     Sig_greater(ie1:ie2,ij(i))=Sig_greater(ie1:ie2,ij(i)) + G_greater((ie1-nop):(ie2-nop),ij(i)) * W_greater(ij(i))   
    !     ref(ie1:ie2,ij(i))=ref(ie1:ie2,ij(i)) + &
    !                             G_lesser((ie1-nop):(ie2-nop),ij(i)) * W_retarded(ij(i)) + &                                      
    !                             G_retarded((ie1-nop):(ie2-nop),ij(i)) * W_lesser(ij(i)) + &
    !                             G_retarded((ie1-nop):(ie2-nop),ij(i)) * W_retarded(ij(i))                                                  
    !     !
    !     ie1 = max(-nop,1) + 1
    !     ie2 = min(nen-nop,nen)
    !     Sig_lesser(ie1:ie2,ij(i))=Sig_lesser(ie1:ie2,ij(i)) + G_lesser((ie1+nop):(ie2+nop),ij(i)) * W_greater(ij(i))   
    !     Sig_greater(ie1:ie2,ij(i))=Sig_greater(ie1:ie2,ij(i)) + G_greater((ie1+nop):(ie2+nop),ij(i)) * W_lesser(ij(i))   
    !     ref(ie1:ie2,ij(i))=ref(ie1:ie2,ij(i)) - &
    !                             G_lesser((ie1+nop):(ie2+nop),ij(i)) * conjg(W_retarded(ij(i))) - &                                      
    !                             G_retarded((ie1+nop):(ie2+nop),ij(i)) * conjg(W_greater(ij(i))) - &
    !                             G_retarded((ie1+nop):(ie2+nop),ij(i)) * conjg(W_retarded(ij(i)))     
    
    ! enddo

    subroutine gpu_selfenergy_GW_with_vertex()

    end subroutine gpu_selfenergy_GW_with_vertex 

    subroutine gpu_vector_scalar(n,devPtrV,scale)
        integer(8),intent(inout) :: devPtrV
        complex(8),intent(in) :: scale
        integer,intent(in) :: n
        call cublas_zscal(n, scale, devPtrV, 1)
    end subroutine gpu_vector_scalar


end module gpu_selfenergy