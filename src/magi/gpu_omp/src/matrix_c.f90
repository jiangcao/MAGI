! Copyright (c) 2023 Jiang Cao, ETH Zurich
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
! 3. Neither the name of the copyright holder nor the names of its contributors
!    may be used to endorse or promote products derived from this software without
!    specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
module matrix_c
!! Complex Matrix Library
!! A 2D complex array, element to form a list/table of complex matrices
use cublas

    implicit none

    type type_matrix_complex
        complex(8), allocatable :: m(:, :)
    !! complex matrix
        integer :: size(2)
    !! matrix size
    end type type_matrix_complex

    interface operator(.m.)
        module procedure array_times_array_simple
    end interface

    interface operator(**)
        module procedure array_power, &
            array_transpose
    end interface

    interface operator(.md.)
        module procedure array_times_array_dagger
    end interface

    interface operator(.dm.)
        module procedure array_dagger_times_array
    end interface

    interface sizeof
        module procedure array_size, matrix_size, matrix_list_size, matrix_size_dim, matrix_list_size2
    end interface sizeof

    interface ReadTxt
        module procedure matrix_read
    end interface ReadTxt

    interface SaveTxt
        module procedure matrix_list_print, array_print
    end interface SaveTxt

    interface show
        module procedure array_print_on_screen
    end interface show

    interface malloc
        module procedure matrix_list_allocElem2, matrix_list_allocElem
        module procedure matrix_alloc, matrix_alloc2, array_alloc, array_alloc2
    end interface

    interface free
        module procedure matrix_free
    end interface

    interface eye
        module procedure array_eye
    end interface

    interface diag
        module procedure array_to_diag
    end interface

    interface trace
        module procedure array_trace, matrix_trace
    end interface

contains

!  =====  Allocation/Deallocation  =====

! allocate a matrix
    pure subroutine matrix_alloc(M, n, nn, source)
        implicit none
        type(type_matrix_complex), intent(out) :: M
        integer, intent(in):: n
        integer, intent(in), optional:: nn
        complex(8), intent(in), optional:: source(:, :)
        if (present(nn)) then
            call matrix_alloc2(M, (/n, nn/), source=source)
        else
            call matrix_alloc2(M, (/n, n/), source=source)
        end if
    end subroutine matrix_alloc

    ! allocate a matrix
    pure subroutine matrix_alloc2(M, n, source)
        implicit none
        type(type_matrix_complex), intent(inout) :: M
        integer, intent(in):: n(2)
        complex(8), intent(in), optional:: source(1:n(1), 1:n(2))
        if (.not. allocated(M%m)) then            
            allocate (M%m(n(1), n(2)))
        else            
            if ((M%size(1) == n(1)) .and. (M%size(2) == n(2))) then
            else                
                deallocate (M%m)
                allocate (M%m(n(1), n(2)))
            end if
        end if
        if (present(source)) then
            M%m(:, :) = source(:, :)
        else
            M%m(:, :) = dcmplx(0.0d0, 0.0d0)
        end if
        M%size = n
    end subroutine matrix_alloc2

    ! allocate an array
    pure subroutine array_alloc(M, n, nn)
        implicit none
        complex(8), intent(out), allocatable:: M(:, :)
        integer, intent(in):: n
        integer, intent(in), optional:: nn
        if (present(nn)) then
            call array_alloc2(M, (/n, nn/))
        else
            call array_alloc2(M, (/n, n/))
        end if
    end subroutine array_alloc

    pure subroutine array_alloc2(M, n)
!! This function allocates a 2D complex array. If the array is already allocated, this function
!! will resize the array to the new size. The allocated array is initiated to zero.
        implicit none
        complex(8), intent(out), allocatable:: M(:, :)
        integer, intent(in):: n(2)
        if (.not. allocated(M)) then
            allocate (M(n(1), n(2)))
        else
            if ((size(M, 1) == n(1)) .and. (size(M, 2) == n(2))) then
            else
                deallocate (M)
                allocate (M(n(1), n(2)))
            end if
        end if
        M = dcmplx(0.0d0, 0.0d0)
    end subroutine array_alloc2

    pure function array_eye(n) result(R)
        implicit none
        integer, intent(in) :: n
        complex(8)::R(n, n)
        INTEGER:: ii
        R = dcmplx(0.0d0, 0.0d0)
        forall (ii=1:n) R(ii, ii) = dcmplx(1.0d0, 0.0d0)
    end function array_eye

    pure subroutine matrix_list_allocElem(this, nx, nm, nn, source)
        implicit none
        integer, intent(in)  :: nx
        integer, intent(in)  :: nm(1:nx)
        integer, intent(in), optional:: nn(1:nx)
        type(type_matrix_complex), intent(out) :: this(1:nx)
        complex(8), intent(in), optional:: source(:, :, :) !! the source data to put into the matrices
        integer :: ii
        do ii = 1, nx
            if (present(nn)) then
                if (present(source)) then
                    call matrix_alloc2(this(ii), (/nm(ii), nn(ii)/), source=source(:, :, ii))
                else
                    call matrix_alloc2(this(ii), (/nm(ii), nn(ii)/))
                endif
            else
                if (present(source)) then
                    call matrix_alloc2(this(ii), (/nm(ii), nm(ii)/), source=source(:, :, ii))
                else
                    call matrix_alloc2(this(ii), (/nm(ii), nm(ii)/))
                endif
            end if
        end do
    end subroutine matrix_list_allocElem

    pure subroutine matrix_list_allocElem2(this, nx, n, source)
        implicit none
        integer, intent(in)  :: nx, n(2, 1:nx)
        type(type_matrix_complex), intent(inout) :: this(1:nx)
        complex(8), intent(in), optional:: source(:, :, :) !! the source data to put into the matrices
        integer :: ii
        do ii = 1, nx
            if (present(source)) then
                call matrix_alloc2(this(ii), n(1:2, ii), source=source(:, :, ii))
            else
                call matrix_alloc2(this(ii), n(1:2, ii))
            endif
        end do
    end subroutine matrix_list_allocElem2

    elemental subroutine matrix_free(this)
        implicit none
        type(type_matrix_complex), intent(out) :: this
        if (allocated(this%m)) deallocate (this%m)
    end subroutine matrix_free

!  =====  size, print  =====

    pure function matrix_list_size(list, dim) result(nm)
        implicit none
        type(type_matrix_complex), intent(in) :: list(:)
        integer, intent(in):: dim
        INTEGER  :: nm(size(list))
        integer :: ii
        forall (ii=1:size(list)) nm(ii) = list(ii)%size(dim)
    end function matrix_list_size

    pure function matrix_list_size2(list) result(nm)
        implicit none
        type(type_matrix_complex), intent(in) :: list(:)
        INTEGER  :: nm(2, size(list))
        integer :: ii
        forall (ii=1:size(list)) nm(:, ii) = list(ii)%size(:)
    end function matrix_list_size2

    pure function array_size(this) result(s)
        implicit none
        complex(8), intent(in) :: this(:, :)
        integer :: s(2), ii
        FORALL (ii=1:2) s(ii) = size(this, dim=ii)
    end function array_size

    pure function matrix_size(this) result(s)
        implicit none
        type(type_matrix_complex), intent(in) :: this
        integer :: s(2)
        s(:) = this%size
    end function matrix_size

    pure function matrix_size_dim(this, dim) result(s)
        implicit none
        type(type_matrix_complex), intent(in) :: this
        integer, intent(in) :: dim
        integer :: s
        s = this%size(dim)
    end function matrix_size_dim

    subroutine matrix_list_print(handle, this)
        implicit none
        type(type_matrix_complex), intent(in) :: this(:)
        integer, intent(in), optional:: handle
        integer :: ii, xx, yy
        if (present(handle)) then
            write (handle, '(1(i8))') size(this)
            do ii = 1, size(this)
                write (handle, '(2(i8))') this(ii)%size(:)
            end do
            write (handle, '(es15.4,es15.4)') (((this(ii)%m(xx, yy), &
                                                 xx=1, size(this(ii)%m, 1)), yy=1, size(this(ii)%m, 2)), ii=1, size(this))
        else
            print '(3(i8),es15.4,es15.4)', (((ii, xx, yy, this(ii)%m(xx, yy), &
                                              xx=1, size(this(ii)%m, 1)), yy=1, size(this(ii)%m, 2)), ii=1, size(this))
        end if
    end subroutine matrix_list_print

    subroutine array_print(handle, A)
        implicit none
        complex(8), intent(in) :: A(:, :)
        integer, intent(in) :: handle
        integer :: xx, yy
        write (handle, '(2(i8))') size(A, 1), size(A, 2)
        write (handle, '(es15.4,es15.4)') ((A(xx, yy), xx=1, size(A, 1)), yy=1, size(A, 2))
        write (handle, '(A)') "END"
    end subroutine array_print

    subroutine matrix_read(handle, A)
        implicit none
        type(type_matrix_complex), intent(out):: A
        integer, intent(in) :: handle
        integer :: xx, yy
        real(8) :: re, im
        character(len=100) :: s
        read (handle, *) xx, yy
        call matrix_alloc(A, xx, yy)
        read (handle, '(100A)') s
        do while (trim(s) /= "END")
            read (s, *) re, im
            A%m(xx, yy) = dcmplx(re, im)
            read (handle, '(100A)') s
        end do
    end subroutine matrix_read

    subroutine array_print_on_screen(A)
        implicit none
        complex(8), intent(in):: A(:, :)
        integer :: xx, yy
        do xx = 1, size(A, 1)
            print '(10(A,es8.1,",",es8.1,")"))', ("(", A(xx, yy), yy=1, size(A, 2))
        end do
    end subroutine array_print_on_screen

    pure function array_testHermitian(M) result(b)
        implicit none
        complex(8), intent(in) :: M(:, :)
        logical :: b
        integer :: i, j
        real(8), parameter :: TOL = 1.0D-10
        b = .true.
        do i = 1, size(M, 2)
            do j = 1, i
                if (abs(M(i, j) - conjg(M(j, i))) .gt. TOL) then
                    b = .false.
                    return
                end if
            end do
        end do
    end function array_testHermitian

!  =====  Multiplications and others  =====

    pure function array_times_array_dagger(A, B) result(C)
        implicit none
        complex(8), intent(in) :: A(:, :), B(:, :)
        complex(8) :: C(size(A, 1), size(B, 1))
        C = array_times_array(A, B, trA=.false., trB=.true., cjA=.false., cjB=.true.)
    end function array_times_array_dagger

    pure function array_dagger_times_array(A, B) result(C)
        implicit none
        complex(8), intent(in) :: A(:, :), B(:, :)
        complex(8) :: C(size(A, 2), size(B, 2))
        C = array_times_array(A, B, trA=.true., trB=.false., cjA=.true., cjB=.false.)
    end function array_dagger_times_array

    pure function array_times_array(A, B, trA, trB, cjA, cjB) result(C)
        implicit none
        complex(8), intent(in) :: A(:, :), B(:, :)
        LOGICAL(KIND=4), intent(in):: trA, trB
        complex(8) :: C(size(A, merge(2, 1, trA)), size(B, merge(1, 2, trB)))
        LOGICAL(KIND=4), intent(in), optional:: cjA, cjB
        integer :: lda, ldb, k, m, kb, n
        character :: ctrA, ctrB
        interface
            pure subroutine ZGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
                COMPLEX(8), intent(in):: ALPHA, BETA
                INTEGER, intent(in):: K, LDA, LDB, LDC, M, N
                CHARACTER, intent(in):: TRANSA, TRANSB
                COMPLEX(8), intent(in):: A(lda, *), B(ldb, *)
                COMPLEX(8), intent(inout):: C(ldc, *)
            end subroutine ZGEMM
        end interface
        lda = size(A, 1)
        ldb = size(B, 1)
        if (.not. trA) then
            k = size(A, 2)
            m = size(A, 1)
            ctrA = 'n'
        else
            k = size(A, 1)
            m = size(A, 2)
            if (present(cjA) .and. cjA) then
                ctrA = 'c'
            else
                ctrA = 't'
            end if
        end if
        if (.not. trB) then
            kb = size(B, 1)
            n = size(B, 2)
            ctrB = 'n'
        else
            kb = size(B, 2)
            n = size(B, 1)
            if (present(cjB) .and. cjB) then
                ctrB = 'c'
            else
                ctrB = 't'
            end if
        end if
        call zgemm(ctrA, ctrB, m, n, k, dcmplx(1.0d0, 0.0d0), A, lda, B, ldb, dcmplx(0.0d0, 0.0d0), C, m)
    end function array_times_array

    pure function array_times_array_simple(A, B) result(C)
        implicit none
        complex(8), intent(in) :: A(:, :), B(:, :)
        complex(8) :: C(size(A, 1), size(B, 2))
        C = array_times_array(A, B, .false., .false.)
    end function array_times_array_simple

    function array_power(A, n) result(C)
        implicit none
        complex(8), intent(in):: A(:, :)
        integer, intent(in) :: n
        complex(8) :: B(size(A, 1), size(A, 1))
        complex(8) :: C(size(A, 1), size(A, 1))
        integer :: ii
        if (n > 0) then
            B = A
            do ii = 2, n
                B = B.m.A
            end do
            C = B
        elseif (n == 0) then
            C = array_eye(size(A, dim=1))
        elseif (n == -1) then
            C = array_inverse(A)
        else
            C = array_inverse(A)
            B = C
            do ii = 2, -n
                B = B.m.C
            end do
            C = B
        end if
    end function array_power

    pure function array_transpose(A, t) result(C)
        implicit none
        complex(8), intent(in) :: A(:, :)
        character, intent(in):: t
        complex(8) :: C(size(A, 2), size(A, 1))
        if ((t == 't') .or. (t == 'T')) then
            C = Transpose(A)
        elseif ((t == 'c') .or. (t == 'C')) then
            C = Transpose(Conjg(A))
        end if
    end function array_transpose

    function array_eigen(A, B, eigvec, itype, uplo) result(eig)
        implicit none
        complex(8), intent(in) :: A(:, :)
        complex(8), intent(in), optional:: B(:, :)
        real(8) :: eig(size(A, 1))
        complex(8), intent(inout), optional :: eigvec(size(A, 1), size(A, 2))
        integer, intent(in), optional:: itype
        CHARACTER, intent(in), optional:: uplo
        integer:: LDA, N, LDB, lwork, INFO, itypeop
        CHARACTER:: jobz, uploop
        real(8):: RWORK(3*size(A, 2))
        complex(8):: work(1 + 4*size(A, 2) + size(A, 2)**2), C(size(A, 1), size(A, 2))
        C(:, :) = A(:, :)
        if (present(eigvec)) then
            jobz = 'V'
        else
            jobz = 'N'
        end if
        uploop = merge(uplo, 'U', present(uplo))
        itypeop = merge(itype, 1, present(itype))
        N = size(A, dim=2)
        LDA = size(A, dim=1)
        LWORK = size(WORK)
        if (present(B)) then
            LDB = size(B, dim=1)
            call zhegv(itypeop, jobz, uploop, N, C, LDA, B, LDB, eig, WORK, LWORK, RWORK, INFO)
            if (INFO .ne. 0) then
                print *, '@array_eigen ZHEGV fails with INFO=', INFO
                call abort()
            end if
        else
            LDB = LDA
            call zheev(jobz, uploop, N, C, LDA, eig, WORK, LWORK, RWORK, INFO)
            if (INFO .ne. 0) then
                print *, '@array_eigen ZHEEV fails with INFO=', INFO
                call abort()
            end if
        end if
        if (present(eigvec)) eigvec = C
    end function array_eigen

    pure function array_to_diag(A) result(diag)
        implicit none
        complex(8), intent(in):: A(:, :)
        complex(8) :: diag(size(A, 1))
        integer :: ii
        do concurrent(ii=1:size(A, 1))
            diag(ii) = A(ii, ii)
        end do
    end function array_to_diag

    function array_inverse2(A, UPLO)  ! for Hermitian matrix
        implicit none
        complex(8), intent(in) :: A(:, :)
        complex(8):: array_inverse2(size(A, dim=1), size(A, dim=1))
        CHARACTER, intent(in):: UPLO
        integer :: info, lda, lwork, n, nnz
        integer :: ipiv(size(A, 1))
        complex(8), allocatable :: work(:, :)
        n = size(A, 1)
        if (n /= size(A, 2)) then
            print *, '@array_inverse, size not square', n, size(A, 2)
            call abort()
        end if
        array_inverse2(:, :) = A(:, :)
        allocate (work(n*n, n*n))
        LDA = size(A, 2)
        call zhetrf(UPLO, n, array_inverse2, LDA, ipiv, WORK, size(WORK), info)
        if (info .ne. 0) then
            print *, '@array_inverse2 ZHETRF fails with INFO=', info
            call abort()
        end if
        call zhetri(UPLO, n, array_inverse2, LDA, ipiv, WORK, info)
        if (info .ne. 0) then
            print *, '@array_inverse2 ZHETRI fails with INFO=', info
            call abort()
        end if
    end function array_inverse2

    function array_inverse(A) ! for General matrix
        implicit none
        complex(8), intent(in) :: A(:, :)
        integer :: info, n
        integer :: ipiv(size(A, 1))
        complex(8), dimension(size(A, dim=1), size(A, dim=1)) :: array_inverse
        complex(8), allocatable :: work(:)
        n = size(A, 1)
        if (n /= size(A, 2)) then
            print *, '@array_inverse, size not square', n, size(A, 2)
            call abort()
        end if
        array_inverse(:, :) = A(:, :)
        allocate (work(n*n))
        call zgetrf(n, n, array_inverse, n, ipiv, info)
        if (info .ne. 0) then
            print *, '@array_inverse ZGETRF fails with INFO=', info
            call abort()
        end if
        call zgetri(n, array_inverse, n, ipiv, work, n*n, info)
        if (info .ne. 0) then
            print *, '@array_inverse ZGETRI fails with INFO=', info
            call abort()
        end if
    end function array_inverse

    !$omp declare target device_type(any)
    pure function array_trace(A) result(tr)
        implicit none
        complex(8), intent(in) :: A(:, :)
        complex(8):: tr
        integer :: ii
        tr = sum((/(A(ii, ii), ii=1, size(A, 1))/))
    end function array_trace

    elemental function matrix_trace(M) result(tr)
        implicit none
        type(type_matrix_complex), intent(in) :: M
        complex(8):: tr
        integer :: ii
        tr = sum((/(M%m(ii, ii), ii=1, M%size(1))/))
    end function matrix_trace

    subroutine matrix_copy(matrices, tab)
        implicit none
        type(type_matrix_complex), intent(in) :: matrices(:)
        complex(8), intent(out):: tab(:, :, :)
        integer :: i
        do concurrent(i=1:size(matrices))
            tab(1:matrices(i)%size(1), 1:matrices(i)%size(2), i) = matrices(i)%m(:, :)
        end do
    end subroutine matrix_copy

    subroutine triMUL_C(A, B, C, R, trA, trB, trC)
        complex(8), intent(in), dimension(:, :) :: A, B, C
        complex(8), intent(inout), allocatable :: R(:, :)
        character, intent(in) :: trA, trB, trC
        complex(8), allocatable, dimension(:, :) :: tmp
        integer :: n, m, k, kb
        if ((trA .ne. 'n') .and. (trA .ne. 'N') .and. (trA .ne. 't') .and. (trA .ne. 'T') &
            .and. (trA .ne. 'c') .and. (trA .ne. 'C')) then
            write (*, *) "ERROR in triMUL_C! trA is wrong: ", trA
            call abort()
        end if
        if ((trB .ne. 'n') .and. (trB .ne. 'N') .and. (trB .ne. 't') .and. (trB .ne. 'T') &
            .and. (trB .ne. 'c') .and. (trB .ne. 'C')) then
            write (*, *) "ERROR in triMUL_C! trB is wrong: ", trB
            call abort()
        end if
        if ((trA .eq. 'n') .or. (trA .eq. 'N')) then
            k = size(A, 2)
            m = size(A, 1)
        else
            k = size(A, 1)
            m = size(A, 2)
        end if
        if ((trB .eq. 'n') .or. (trB .eq. 'N')) then
            kb = size(B, 1)
            n = size(B, 2)
        else
            kb = size(B, 2)
            n = size(B, 1)
        end if
        if (k .ne. kb) then
            write (*, *) "ERROR in triMUL_C! Matrix dimension is wrong", k, kb
            call abort()
        end if
        call MUL_C(A, B, trA, trB, tmp)
        call MUL_C(tmp, C, 'n', trC, R)
        deallocate (tmp)
    end subroutine triMUL_C

    subroutine MUL_C(A, B, trA, trB, R)
        complex(8), intent(in) :: A(:, :), B(:, :)
        complex(8), intent(inout), allocatable :: R(:, :)
        CHARACTER, intent(in) :: trA, trB
        integer :: n, m, k, kb, lda, ldb
        if ((trA .ne. 'n') .and. (trA .ne. 'N') .and. (trA .ne. 't') .and. (trA .ne. 'T') &
            .and. (trA .ne. 'c') .and. (trA .ne. 'C')) then
            write (*, *) "ERROR in MUL_C! trA is wrong: ", trA
            call abort()
        end if
        if ((trB .ne. 'n') .and. (trB .ne. 'N') .and. (trB .ne. 't') .and. (trB .ne. 'T') &
            .and. (trB .ne. 'c') .and. (trB .ne. 'C')) then
            write (*, *) "ERROR in MUL_C! trB is wrong: ", trB
            call abort()
        end if
        lda = size(A, 1)
        ldb = size(B, 1)
        if ((trA .eq. 'n') .or. (trA .eq. 'N')) then
            k = size(A, 2)
            m = size(A, 1)
        else
            k = size(A, 1)
            m = size(A, 2)
        end if
        if ((trB .eq. 'n') .or. (trB .eq. 'N')) then
            kb = size(B, 1)
            n = size(B, 2)
        else
            kb = size(B, 2)
            n = size(B, 1)
        end if
        if (k .ne. kb) then
            write (*, *) "ERROR in MUL_C! Matrix dimension is wrong", k, kb
            call abort()
        end if
        if (allocated(R)) then
            if ((size(R, 1) .ne. m) .or. (size(R, 2) .ne. n)) then
                deallocate (R)
                Allocate (R(m, n))
            end if
        else
            Allocate (R(m, n))
        end if
        R = dcmplx(0.0d0, 0.0d0)
        call zgemm(trA, trB, m, n, k, dcmplx(1.0d0, 0.0d0), A, lda, B, ldb, dcmplx(0.0d0, 0.0d0), R, m)
    end subroutine MUL_C

end module matrix_c
