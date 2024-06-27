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
module gpu_selfenergy
    implicit none 
    
    integer,parameter :: size_of_complex=16 !8->single precision; 16->double precision

    contains

    subroutine gpu_add_selfenergy_GW(nen,nop,nm,num_ij,ij,copy_to_gpu,copy_to_cpu,G_lesser,G_greater,G_retarded,&
                                W_lesser,W_greater,W_retarded,&
                                devPtrGG,devPtrGL,devPtrGR,&
                                devPtrSigG,devPtrSigL,devPtrSigR,&
                                sig_lesser,sig_greater,sig_retarded)
        integer,intent(in) :: nop,nen,nm,num_ij,ij(:)
        integer(8),intent(inout) :: devPtrGL, devPtrGG, devPtrGR
        integer(8),intent(inout) :: devPtrSigL, devPtrSigG, devPtrSigR
        complex(8), dimension(:,:), intent(in) :: G_greater, G_lesser, G_retarded
        complex(8), dimension(:), intent(in) :: W_greater, W_lesser, W_retarded
        complex(8), dimension(:,:), intent(inout) :: sig_lesser,sig_greater,sig_retarded
        logical,intent(in) :: copy_to_gpu,copy_to_cpu
        integer :: n,i,j                
        n = nen - nop         
        !
        if (copy_to_gpu) then 
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGG)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGL)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGR)
            !   
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrSigG)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrSigL)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrSigR)        
            !copy data to GPU
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_lesser,nen,devPtrGL,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_greater,nen,devPtrGG,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_retarded,nen,devPtrGR,nen)
            !
            call cublas_set_matrix(nen,nm*nm,size_of_complex,Sig_lesser,nen,devPtrSigL,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,Sig_greater,nen,devPtrSigG,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,Sig_retarded,nen,devPtrSigR,nen)
        endif
        !
        do i=1,num_ij                   
            ! sig^<[E] += G^<[E-dE] W^< - G^<[E+dE] conjg(W^>)                                        
            call cublas_zaxpy(n, W_lesser(ij(i)), devPtrGL + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrSigL+(((ij(i)-1)*nen+nop)*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_greater(ij(i)) ), devPtrGL + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrSigL+(((ij(i)-1)*nen)*size_of_complex), 1)
            !
            ! sig^>[E] += G^>[E-dE] W^> - G^>[E+dE] conjg(W^<)    
            call cublas_zaxpy(n, W_greater(ij(i)), devPtrGG + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrSigG+(((ij(i)-1)*nen+nop)*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_lesser(ij(i))), devPtrGG + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrSigG+(((ij(i)-1)*nen)*size_of_complex), 1)
            !
            ! sig^r[E] += G^<[E-dE] W^r + G^r[E-dE] W^< + G^r[E-dE] W^r - G^<[E+dE] conjg(W^r) - G^r[E+dE] conjg(W^>) - G^r[E-dE] conjg(W^r)   
            call cublas_zaxpy(n, W_retarded(ij(i)), devPtrGL + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrSigR+(((ij(i)-1)*nen+nop)*size_of_complex), 1)
            call cublas_zaxpy(n, W_lesser(ij(i)),   devPtrGR + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrSigR+(((ij(i)-1)*nen+nop)*size_of_complex), 1)
            call cublas_zaxpy(n, W_retarded(ij(i)), devPtrGR + ((ij(i)-1) * nen)*size_of_complex, 1, devPtrSigR+(((ij(i)-1)*nen+nop)*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_retarded(ij(i))), devPtrGL + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrSigR+(((ij(i)-1)*nen)*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_greater(ij(i))),  devPtrGR + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrSigR+(((ij(i)-1)*nen)*size_of_complex), 1)
            call cublas_zaxpy(n, -conjg(W_retarded(ij(i))), devPtrGR + ((ij(i)-1) * nen + nop)*size_of_complex, 1, devPtrSigR+(((ij(i)-1)*nen)*size_of_complex), 1)
        enddo 
        !
        ! copy data from GPU
        if (copy_to_cpu) then 
            call cublas_get_matrix(nen,nm*nm,size_of_complex,devPtrSigG,nen,Sig_greater,nen)
            call cublas_get_matrix(nen,nm*nm,size_of_complex,devPtrSigL,nen,Sig_lesser,nen)
            call cublas_get_matrix(nen,nm*nm,size_of_complex,devPtrSigR,nen,sig_retarded,nen)  
        endif
    end subroutine gpu_add_selfenergy_GW 


    subroutine gpu_selfenergy_GW_with_vertex()

    end subroutine gpu_selfenergy_GW_with_vertex 

    subroutine gpu_vector_scalar(n,devPtrV,scale)
        integer(8),intent(inout) :: devPtrV
        complex(8),intent(in) :: scale
        integer,intent(in) :: n
        call cublas_zscal(n, scale, devPtrV, 1)
    end subroutine gpu_vector_scalar


end module gpu_selfenergy
