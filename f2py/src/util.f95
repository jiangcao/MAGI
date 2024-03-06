! -*- f90 -*-

module parameters_mod   
    implicit none 
    !constants
    integer, parameter :: dp=8
    REAL(kind=dp), PARAMETER :: pi=3.14159265359d0
    REAL(kind=dp), PARAMETER :: twopi = 3.14159265359d0*2.0d0
    REAL(kind=dp), PARAMETER :: e_charge=1.6d-19            ! charge of an electron (C)
    REAL(kind=dp), PARAMETER :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)    
    REAL(kind=dp), PARAMETER :: light_speed=2.998d8           ! m/s
    REAL(kind=dp), PARAMETER :: m0_charge=5.6856D-16        ! eV s2 / cm2
    REAL(kind=dp), PARAMETER :: hbar=1.0546d-34     ! value of hbar=h/2pi (J s)
    REAL(kind=dp), PARAMETER :: hbar_eV=hbar/e_charge ! eV s    
    COMPLEX(kind=dp), PARAMETER :: cone = dcmplx(1.0d0,0.0d0)
    COMPLEX(kind=dp), PARAMETER :: czero  = dcmplx(0.0d0,0.0d0)
    COMPLEX(kind=dp), PARAMETER :: c1i  = dcmplx(0.0d0,1.0d0)     
    REAL(kind=dp), PARAMETER  :: BOLTZ = 8.61734d-05 !eV K-1 
end module parameters_mod


module linalg 
    implicit none 
    contains

    ! calculate eigen-values and eigen-vectors of a Hermitian matrix A
    !   upon return A will be modified and contains the eigen-vectors
    SUBROUTINE eigv(NN, A, V, W)
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(IN), DIMENSION(:, :) :: A
        REAL(8), INTENT(OUT), DIMENSION(NN) :: W
        COMPLEX(8), INTENT(OUT), DIMENSION(NN, NN) :: V
        ! -----                
        integer :: INFO, LWORK, lrwork
        complex(8), allocatable :: work(:)
        real(8), allocatable :: RWORK(:)
        !integer, allocatable :: iwork(:)
        lwork = max(1, 2*NN - 1)
        lrwork = max(1, 3*NN - 2)
        allocate (work(lwork))
        allocate (rwork(lrwork))
        V = A
        CALL zheev('V', 'U', NN, V, NN, W, WORK, LWORK, RWORK, INFO)
        !
        deallocate (work, rwork)
        if (INFO .ne. 0) then
            write (*, *) 'SEVERE WARNING: ZHEEV HAS FAILED. INFO=', INFO
            call abort()
        end if
    END SUBROUTINE eigv


    SUBROUTINE norm(vector, val)
        REAL(8),intent(in) :: vector(3)
        REAL(8),intent(out) :: val
        val = sqrt(dot_product(vector,vector))
    END SUBROUTINE norm

    ! vector cross-product
    SUBROUTINE cross(a, b, C)
        REAL(8), DIMENSION(3), INTENT(out) :: C
        REAL(8), DIMENSION(3), INTENT(IN) :: a, b
        C(1) = a(2)*b(3) - a(3)*b(2)
        C(2) = a(3)*b(1) - a(1)*b(3)
        C(3) = a(1)*b(2) - a(2)*b(1)
    END SUBROUTINE cross

    ! calculate eigen-values of a Hermitian matrix A
    subroutine eig(A, NN, W)
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(IN), DIMENSION(:, :) :: A
        REAL(8), INTENT(OUT), DIMENSION(NN) :: W
        ! -----        
        integer :: INFO, LWORK, lrwork
        complex(8), allocatable :: work(:)
        real(8), allocatable :: RWORK(:)        
        lwork = max(1, 2*NN - 1)
        lrwork = max(1, 3*NN - 2)
        allocate (work(lwork))
        allocate (rwork(lrwork))
        !
        CALL zheev('N', 'U', NN, A, NN, W, WORK, LWORK, RWORK, INFO)
        !
        deallocate (work, rwork)
        if (INFO .ne. 0) then
            write (*, *) 'SEVERE WARNING: ZHEEV HAS FAILED. INFO=', INFO
            call abort()
        end if        
    end subroutine eig

    ! matrix inversion
    subroutine invert(A, nn, invA)
        integer :: info
        integer, intent(in) :: nn
        integer, dimension(:), allocatable :: ipiv
        complex(8), dimension(nn, nn), intent(in) :: A
        complex(8), dimension(nn, nn), intent(out) :: invA
        complex(8), dimension(:), allocatable :: work
        COMPLEX(8), PARAMETER :: czero  = dcmplx(0.0d0,0.0d0)
        allocate (work(nn*nn))
        allocate (ipiv(nn))
        invA=A
        call zgetrf(nn, nn, invA, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            invA = czero
        else
            call zgetri(nn, invA, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info
                invA = czero
            end if
        end if
        deallocate (work)
        deallocate (ipiv)
    end subroutine invert

end module linalg
