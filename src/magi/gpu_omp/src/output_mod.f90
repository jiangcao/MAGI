module Output
    use matrix_c, only: type_matrix_complex, sizeof
    implicit none

    private

    integer, parameter :: dp = 8

    public::write_spectrum_summed_over_k

contains

    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_summed_over_k(dataset, i, G, nen, en, nk, length, NB, NS, Lx, coeff, append)        
        character(len=*), intent(in) :: dataset
        type(type_matrix_complex), intent(in), dimension(length, nen, nk)::G
        integer, intent(in)::i, nen, length, NB, NS, nk
        real(dp), intent(in)::Lx, en(nen), coeff(2)
        logical, intent(in), optional::append
        ! ----
        integer:: ix,ie, j, ib, k, ik
        character(len=4) :: i_str
        character(len=8) :: fmt
        real(dp)::xx
        complex(dp)::tr
        logical append_
        append_ = .false.
        if (present(append)) append_ = append
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        if (append_) then
            open (unit=11, file=trim(dataset)//i_str//'.dat', status='unknown', position='append')
        else
            open (unit=11, file=trim(dataset)//i_str//'.dat', status='unknown')
        end if
        !        
        do ie = 1, nen
            xx = 0.0d0
            ix = 0
            do j = 1, length
                do k = 1, NS
                    ix = ix + 1
                    xx = xx + Lx
                    tr = 0.0d0
                    !!!$omp parallel default(shared) private(ik, ib) reduction(+:tr)
                    !!!$omp do
                    do ik = 1, nk
                        do ib = 1, nb
                            tr = tr + G(j, ie, ik)%m(ib + (k - 1)*NB, ib + (k - 1)*NB)
                        end do
                    end do
                    !!!$omp end do
                    !!!$omp end parallel
                    tr = tr/dble(nk)
                    write (11, '(4E18.4)') xx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)
                end do
            end do
            write (11, *)
        end do
        close (11)
    end subroutine write_spectrum_summed_over_k

end module Output
