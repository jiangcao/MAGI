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
module gpu_polarization
    implicit none 
    integer,parameter :: size_of_complex=16 !8->single precision; 16->double precision

    contains 

    subroutine gpu_polarization(a1,a2,nop,nen,nm,num_jl,jl,num_ki,ki,copy_to_gpu,G_lesser,G_greater,G_retarded,G_advanced,&
        devPtrGL,devPtrGG,devPtrGR,devPtrGA,partial_P)
        complex(8),intent(in) :: a1,a2
        integer,intent(in) :: nop,nen,nm,jl(:),ki(:),num_jl,num_ki
        integer(8),intent(inout) :: devPtrGL, devPtrGG, devPtrGR, devPtrGA
        complex(8), dimension(:,:), intent(in) :: G_greater, G_lesser, G_retarded, G_advanced
        complex(8), dimension(:,:), intent(out) :: partial_P
        logical,intent(in) :: copy_to_gpu
        integer :: n, i        
        integer(8) :: devPtrA, devPtrB, devPtrC        
        !
        if (copy_to_gpu) then 
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGG)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGL)
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGR) 
            call cublas_alloc(nm*nm*nen, size_of_complex, devPtrGA)                
            !copy data to GPU
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_lesser,nen,devPtrGL,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_greater,nen,devPtrGG,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_retarded,nen,devPtrGR,nen)
            call cublas_set_matrix(nen,nm*nm,size_of_complex,G_advanced,nen,devPtrGA,nen)
        endif
        n = nen - nop
        call cublas_alloc(n*num_jl, size_of_complex, devPtrA)
        call cublas_alloc(n*num_ki, size_of_complex, devPtrB)
        call cublas_alloc(num_jl*num_ki, size_of_complex, devPtrC)
        ! G^>_jl G^<_ki
        do i = 1,num_jl 
            call cublas_zcopy(n, devPtrGG+((jl(i)-1) * nen + nop)*size_of_complex, 1, devPtrA+(i-1)*n*size_of_complex, 1)
        enddo
        do i = 1,num_ki                
            call cublas_zcopy(n, devPtrGL+((ki(i)-1) * nen)*size_of_complex      , 1, devPtrB+(i-1)*n*size_of_complex, 1)
        enddo
        call cublas_zgemm('t','n',num_jl,num_ki,n,a2,devPtrA,n,devPtrB,n,dcmplx(0.0,0.0),devPtrC,num_jl)
        ! G^<_jl G^>_ki
        do i = 1,num_jl 
            call cublas_zcopy(n, devPtrGL+((jl(i)-1) * nen + nop)*size_of_complex, 1, devPtrA+(i-1)*n*size_of_complex, 1)
        enddo
        do i = 1,num_ki                
            call cublas_zcopy(n, devPtrGG+((ki(i)-1) * nen)*size_of_complex      , 1, devPtrB+(i-1)*n*size_of_complex, 1)
        enddo
        call cublas_zgemm('t','n',num_jl,num_ki,n,-a2,devPtrA,n,devPtrB,n,dcmplx(1.0,0.0),devPtrC,num_jl)
        ! G^<_jl G^A_ki 
        do i = 1,num_jl 
            call cublas_zcopy(n, devPtrGL+((jl(i)-1) * nen + nop)*size_of_complex, 1, devPtrA+(i-1)*n*size_of_complex, 1)
        enddo
        do i = 1,num_ki                
            call cublas_zcopy(n, devPtrGA+((ki(i)-1) * nen)*size_of_complex      , 1, devPtrB+(i-1)*n*size_of_complex, 1)
        enddo
        call cublas_zgemm('t','n',num_jl,num_ki,n,a1,devPtrA,n,devPtrB,n,dcmplx(1.0,0.0),devPtrC,num_jl)
        ! G^R_jl G^<_ki
        do i = 1,num_jl 
            call cublas_zcopy(n, devPtrGR+((jl(i)-1) * nen + nop)*size_of_complex, 1, devPtrA+(i-1)*n*size_of_complex, 1)
        enddo
        do i = 1,num_ki                
            call cublas_zcopy(n, devPtrGL+((ki(i)-1) * nen)*size_of_complex      , 1, devPtrB+(i-1)*n*size_of_complex, 1)
        enddo
        call cublas_zgemm('t','n',num_jl,num_ki,n,a1,devPtrA,n,devPtrB,n,dcmplx(1.0,0.0),devPtrC,num_jl)

        !copy data from GPU
        call cublas_get_matrix(num_jl,num_ki,size_of_complex,devPtrC,num_jl,partial_P,num_ki)
    
        !Free GPU memory
        call cublas_free(devPtrA)
        call cublas_free(devPtrB)
        call cublas_free(devPtrC)
    end subroutine gpu_polarization

    subroutine test_cublaszgemm(m,n,A,B,C)
        integer,intent(in) :: m,n
        complex(8), dimension(:,:), intent(in) :: A, B
        complex(8), dimension(:,:), intent(out) :: C
        real(8) :: gflops
        real(8) :: tstart, tstop, elapsed_time
        integer,dimension(8) :: values
        integer(8) :: devPtrA, devPtrB, devPtrC

        call cublas_alloc(m*n, size_of_complex, devPtrA)
        call cublas_alloc(n*m, size_of_complex, devPtrB)
        call cublas_alloc(m*m, size_of_complex, devPtrC)
    
        call date_and_time(VALUES=values) !values(8) = milisecs of the second
        ! seed = values(8) !using value in milisecs as seeder
        ! call srand(seed) !not a std implementation, but i like it better.
                
        !copy data to GPU
        call cublas_set_matrix(m,n,size_of_complex,A,m,devPtrA,m)
        call cublas_set_matrix(n,m,size_of_complex,B,n,devPtrB,n)
        call cublas_set_matrix(m,m,size_of_complex,C,m,devPtrC,m)
    
        call cpu_time(tstart)
        !call SGEMM from CUBLAS
        call cublas_zgemm('n','n',m,m,n,dcmplx(1.0,0.0),devPtrA,m,devPtrB,n,dcmplx(0.0,0.0),devPtrC,m)
        call cpu_time(tstop)

        !copy data from GPU
        call cublas_get_matrix(m,m,size_of_complex,devPtrC,m,C,m)
        ! call cpu_time(tstop)
    
        elapsed_time = tstop - tstart !in seconds
    
        write(*,*) 'Matrix A ', m,'x',n
        write(*,*) 'Matrix B ', n,'x',m
        write(*,20) 'Elapsed time : ',elapsed_time, 'secs'
    
        gflops = 2*float(n)*float(n)*float(n)/(elapsed_time*1.0e9)
        write(*,10)  'Performance:', gflops,' GFLOPS'
        
    10 format(A12,2X,1F0.4,2X,A7)
    20 format(A15,2X,1F0.8,2X,A4)

        !Free GPU memory
        call cublas_free(devPtrA)
        call cublas_free(devPtrB)
        call cublas_free(devPtrC)
    end subroutine test_cublaszgemm

end module gpu_polarization
