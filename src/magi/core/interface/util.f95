! -*- f90 -*-
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
