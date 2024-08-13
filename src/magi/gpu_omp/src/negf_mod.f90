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
!
module negf_mod
    !! Non-equilibrium Green's function (NEGF) module, upper-level driver for solving the NEGF equations

    implicit none

    private

    integer, parameter :: dp = 8

    public :: negf_solve

contains

    subroutine negf_solve(nx, nen, nk, emin, emax, Hii, H1i, Sii, temp, mu, &
                          comm_size, comm_rank, local_NE, first_local_energy, nbnd, nslab, Lx)
        use matrix_c, only: type_matrix_complex, malloc, free, sizeof
        use cuda_rgf_mod, only: cuda_rgf_variableblock_forward, cuda_rgf_constblocksize,cuda_rgf_init,cuda_rgf_finish
        ! use rgf_mod, only: rgf_variableblock_forward
        use Output, only: write_spectrum_summed_over_k
        use omp_lib
        type(type_matrix_complex), intent(in), dimension(nx, nk)::Hii, Sii
        type(type_matrix_complex), intent(in), dimension(nx + 1, nk)::H1i
        integer(kind=4), intent(in) :: comm_size, comm_rank, local_NE, first_local_energy
        integer, intent(in)::nx  !! number of slabs
        integer, intent(in)::nen  !! number of energy points
        integer, intent(in)::nk  !! number of k points
        integer, intent(in)::nbnd  !! number of bands / orbitals per cell
        integer, intent(in)::nslab  !! number of cells in a slab
        real(dp), intent(in), dimension(2)::temp !! temperatures
        real(dp), intent(in), dimension(2):: mu !! chemical potentials
        real(dp), intent(in):: Lx !! Lx
        ! ----
        real(dp)::emin, emax !! min and max energy range
        integer::nm(2, nx), ik, ie, iter, NB, NS, i
        integer(kind=4)::ierr
        real(dp)::en(nen), dE, local_energies(local_NE)
        real(dp), dimension(:, :), allocatable::mul, mur, templ, tempr
        type(type_matrix_complex), dimension(nx, local_NE, nk)::sigma_lesser_ph, sigma_r_ph, G_r, G_lesser, G_greater
        type(type_matrix_complex), dimension(:),allocatable::Jdens, Gl, Gln
        real(dp), dimension(nen, nk)::tr, tre
        character(len=20) :: filename
        character(len=8) :: fmt
        character(len=4) :: rank_str
        logical::append
        !
    !    include "mpif.h"
        fmt = '(I4.4)'
        write (rank_str, fmt) comm_rank
        append = (comm_rank /= 0)
        nm = sizeof(Hii(:, 1))
        allocate (mul(nm(1, 1), nm(1, 1)))
        allocate (templ(nm(1, 1), nm(1, 1)))
        allocate (mur(nm(1, nx), nm(1, nx)))
        allocate (tempr(nm(1, nx), nm(1, nx)))
        mul = mu(1)
        mur = mu(2)
        templ = temp(1)
        tempr = temp(2)
        do ik = 1, nk
            do ie = 1, local_NE
                call malloc(sigma_lesser_ph(:, ie, ik), nx, nm)
                call malloc(sigma_r_ph(:, ie, ik), nx, nm)
                call malloc(G_r(:, ie, ik), nx, nm)
                call malloc(G_lesser(:, ie, ik), nx, nm)
                call malloc(G_greater(:, ie, ik), nx, nm)
            end do
        end do
        !
        if (comm_rank == 0) then
            print *, 'allocate memory DONE'
        end if
        !
        dE = (emax - emin)/dble(nen - 1)
        forall (ie=1:nen) En(ie) = emin + dble(ie - 1)*dE  ! global energy vector
        do ie = 1, local_NE
            local_energies(ie) = En(ie + first_local_energy - 1) ! local energy vector of rank-i
        end do
        !
        iter = 0
        !
        !$omp parallel default(shared) private(ie,ik,Jdens,Gl,Gln)
        allocate(Jdens(nx),Gl(nx),Gln(nx))        
        call malloc(Jdens, nx, nm)
        ! call malloc(Gl, nx, nm)
        ! call malloc(Gln, nx, nm)        
        call cuda_rgf_init(nm(1,1))
        !!$omp do        
        do ie = 1, local_NE
            do ik = 1, nk
                !  call rgf_variableblock_forward(nx, local_energies(ie), mul, mur, TEMPl, TEMPr, &
                !      Hii(:, ik), H1i(:, ik), Sii(:, ik), sigma_lesser_ph(:, ie, ik), &
                !      sigma_r_ph(:, ie, ik), G_r(:, ie, ik), G_lesser(:, ie, ik), G_greater(:, ie, ik), &
                !      Jdens, Gl, Gln, tr(ie, ik), tre(ie, ik))

               call cuda_rgf_constblocksize(nm(1,1), nx, local_energies(ie), mul, mur, TEMPl, TEMPr, &
                   Hii(:, ik), H1i(:, ik), Sii(:, ik), sigma_lesser_ph(:, ie, ik), &
                   sigma_r_ph(:, ie, ik), G_r(:, ie, ik), G_lesser(:, ie, ik), G_greater(:, ie, ik), &
                   Jdens, tr(ie, ik), tre(ie, ik))    
            end do
        end do        
        !!$omp end do
        call cuda_rgf_finish()
        call free(Jdens)
        ! call free(Gl)
        ! call free(Gln)
        deallocate(Jdens)
        ! deallocate(Gl,Gln)
        !$omp end parallel
        !
        NB=nbnd
        NS=nslab
        do i = 0, comm_size - 1
            if (i == comm_rank) then
                filename = 'ldos'
                call write_spectrum_summed_over_k(filename, iter, G_r, local_NE, local_energies, &
                                                  nk, nx, NB, NS, Lx, (/1.0d0, -2.0d0/), append)
                filename = 'ndos'
                call write_spectrum_summed_over_k(filename, iter, G_lesser, local_NE, local_energies, &
                                                  nk, nx, NB, NS, Lx, (/1.0d0, 1.0d0/), append)
                filename = 'pdos'
                call write_spectrum_summed_over_k(filename, iter, G_greater, local_NE, local_energies, &
                                                  nk, nx, NB, NS, Lx, (/1.0d0, -1.0d0/), append)
            end if
!            call MPI_Barrier(MPI_COMM_WORLD, ierr)
        end do
        !
        if (comm_rank == 0) then
            print *, 'free memory'
        end if
        deallocate (mul, mur, tempr, templ)
        !
        call free(sigma_lesser_ph)
        call free(sigma_r_ph)
        call free(G_r)
        call free(G_lesser)
        call free(G_greater)
        !
    end subroutine negf_solve

end module negf_mod
