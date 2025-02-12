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

module gw_sparse

use fft_mod, only: conv1d => conv1d2, corr1d => corr1d2 
use parameters_mod, only: c1i, cone, czero, pi, twopi, two
use rgf_mod, only: RGF_G

implicit none 

private

CONTAINS


subroutine green_rgf_solve_gw_ephoton_3d_mpi(alpha_mix,niter,NB,NS,nm,nx, &
   nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg,Pii,P1i,&
   polarization,intensity,hw,labs, comm_size, comm_rank, local_NE, first_local_energy)
  !
  integer ( kind = 4), intent(in) :: comm_size, comm_rank, local_NE, first_local_energy
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  !complex(8), intent(in):: V(nm*nx,nm*nx,nky*nkz)
  real(8), intent(in) :: polarization(3) ! light polarization vector 
  real(8), intent(in) :: intensity ! [W/m^2]
  logical, intent(in) :: labs ! whether to calculate Pi and absorption
  complex(8), intent(in):: Pii(nm,nm,3,nx,nky*nkz),P1i(nm,nm,3,nx,nky*nkz) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
  real(8), intent(in) :: hw ! hw is photon energy in eV
  ! -------- local variables
  real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*epsilon0*c0**3) 
  ! complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser,cur, g_r_i1
  ! complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  ! complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  ! complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser,P_greater,P_retarded
  ! complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:)::Sii
  complex(8),allocatable,dimension(:,:,:,:)::Mii,M1i
  real(8)::tr(local_NE,nky*nkz),tre(local_NE,nky*nkz)
  integer::ie,iter,i,ix,nopmax,nop,nopphot,iop,l,h,io,j, global_ij, row, col, ij
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk,ne
  complex(8)::dE, B(nm,nm)
  character ( len = 20 ) :: filename
  character ( len = 8 ) :: fmt
  character ( len = 4 ) :: rank_str
  character ( len = 8 ) :: cmethod
  logical append

  real(8) :: local_energies(local_NE)

  integer ( kind = 4 ) :: ierr, local_Nij, first_local_ij, local_NX, first_local_block
  real(8) :: local_sum_tr, global_sum_tr, local_sum_tre, global_sum_tre
  integer ( kind = 8 ) :: num_g, num_p, num_err, num, count

  integer reqs(8), stats(8)

  complex(8), pointer :: g_r_buf(:), g_r_by_energies(:, :, :, :, :), g_r_by_blocks(:, :, :, :)
  complex(8), pointer :: g_greater_extended_buf(:), g_greater_buf(:), g_greater_by_energies(:, :, :, :, :), g_greater_by_blocks(:, :, :, :)
  complex(8), pointer :: g_lesser_extended_buf(:), g_lesser_buf(:), g_lesser_by_energies(:, :, :, :, :), g_lesser_by_blocks(:, :, :, :)

  complex(8), pointer :: sigma_r_gw_buf(:), sigma_r_gw_by_energies(:, :, :, :, :), sigma_r_gw_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_greater_gw_buf(:), sigma_greater_gw_by_energies(:, :, :, :, :), sigma_greater_gw_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_lesser_gw_buf(:), sigma_lesser_gw_by_energies(:, :, :, :, :), sigma_lesser_gw_by_blocks(:, :, :, :)

  complex(8), pointer :: sigma_r_new_buf(:), sigma_r_new_by_energies(:, :, :, :, :), sigma_r_new_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_greater_new_buf(:), sigma_greater_new_by_energies(:, :, :, :, :), sigma_greater_new_by_blocks(:, :, :, :)
  complex(8), pointer :: sigma_lesser_new_buf(:), sigma_lesser_new_by_energies(:, :, :, :, :), sigma_lesser_new_by_blocks(:, :, :, :)

  complex(8), pointer :: P_retarded_buf(:), P_retarded_by_energies(:, :, :, :, :), P_retarded_by_blocks(:, :, :, :)
  complex(8), pointer :: P_greater_buf(:), P_greater_by_energies(:, :, :, :, :), P_greater_by_blocks(:, :, :, :)
  complex(8), pointer :: P_lesser_buf(:), P_lesser_by_energies(:, :, :, :, :), P_lesser_by_blocks(:, :, :, :)

  complex(8), pointer :: W_retarded_buf(:), W_retarded_by_energies(:, :, :, :, :), W_retarded_by_blocks(:, :, :, :)
  complex(8), pointer :: W_greater_buf(:), W_greater_by_energies(:, :, :, :, :), W_greater_by_blocks(:, :, :, :)
  complex(8), pointer :: W_lesser_buf(:), W_lesser_by_energies(:, :, :, :, :), W_lesser_by_blocks(:, :, :, :)

  real(8), allocatable, dimension(:, :, :, :, :) :: cur, jdens_local
  real(8), allocatable, dimension(:, :, :, :) :: tot_cur,tot_ecur,tot_cur_local,tot_ecur_local

  complex(8), pointer :: tmp0(:), tmp1(:), g_lesser_extended(:, :, :, :, :), g_greater_extended(:, :, :, :, :)

  complex(8), pointer :: g_greater_t_buf(:), g_greater_t_by_energies(:, :, :, :, :), g_greater_t_by_blocks(:, :, :, :)
  complex(8), pointer :: g_lesser_t_buf(:), g_lesser_t_by_energies(:, :, :, :, :), g_lesser_t_by_blocks(:, :, :, :)

  complex(8), pointer :: g_lesser_photon_left_send(:), g_lesser_photon_left_recv(:)
  complex(8), pointer :: g_lesser_photon_right_send(:), g_lesser_photon_right_recv(:)
  complex(8), pointer :: g_greater_photon_left_send(:), g_greater_photon_left_recv(:)
  complex(8), pointer :: g_greater_photon_right_send(:), g_greater_photon_right_recv(:)

  complex(8), pointer :: g_lesser_photon_buf(:), g_greater_photon_buf(:)
  complex(8), pointer :: g_lesser_photon(:, :, :, :), g_greater_photon(:, :, :, :)

  real(8) :: start, finish, it_start

  real(8), allocatable :: extended_local_energies(:)
  
  complex(8),allocatable::Ispec(:,:,:,:),Itot(:,:,:),Ispec_ik(:,:,:,:),Itot_ik(:,:,:) ! collision integral variables
  
  complex(8),allocatable,dimension(:,:,:)::Pi_retarded_ik,Pi_lesser_ik,Pi_greater_ik,Pi_retarded ! photon Pi self-energies

  local_Nij = (nm * nm) / comm_size
  first_local_ij = local_Nij * comm_rank + 1

  local_NX = nx / comm_size
  first_local_block = local_NX * comm_rank + 1

  fmt = '(I4.4)'
  write ( rank_str, fmt ) comm_rank
  append = (comm_rank /= 0)

  if (comm_rank == 0) then
    print *,'======== green_rgf_solve_gw_ephoton_3D ========'
  endif

  do ie = 1, local_NE
    local_energies(ie) = en(ie + first_local_energy - 1)
  enddo

  nk = nky * nkz
  ne = nen
  nopphot=floor(hw / (En(2)-En(1)))

  allocate(extended_local_energies(local_NE + 2 * nopphot))

  do ie = 1, local_NE
    local_energies(ie) = en(ie + first_local_energy - 1)
  enddo
  extended_local_energies = 0.0d0
  extended_local_energies(nopphot + 1:nopphot + local_NE) = local_energies(:)
  if (comm_rank - 1 >= 0) then
    do ie = 1, nopphot
      extended_local_energies(ie) = en(ie + first_local_energy - nopphot - 1)
    enddo
  endif
  if (comm_rank +1 < comm_size) then
    do ie = 1, nopphot
      extended_local_energies(ie + nopphot + local_NE) = en(local_NE + first_local_energy + ie - 1)
    enddo
  endif

  allocate(g_r_buf(nm * nm * local_NE * nk * nx))
  allocate(g_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(g_lesser_buf(nm * nm * nx * local_NE * nk))  

  ! For photon computation/communication
  allocate(g_lesser_photon_buf(nm * nm * nx * (local_NE + 2 * nopphot)))
  allocate(g_greater_photon_buf(nm * nm * nx * (local_NE + 2 * nopphot)))
  g_lesser_photon_buf = dcmplx(0.0d0, 0.0d0)
  g_greater_photon_buf = dcmplx(0.0d0, 0.0d0)
  g_lesser_photon(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot) => g_lesser_photon_buf
  g_greater_photon(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot) => g_greater_photon_buf
  g_lesser_photon_left_send(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(nm * nm * nx * nopphot + 1 : nm * nm * nx * 2 * nopphot)
  g_lesser_photon_left_recv(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(1 : nm * nm * nx * nopphot)
  g_lesser_photon_right_send(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(nm * nm * nx * local_NE + 1 : nm * nm * nx * (local_NE + nopphot))
  g_lesser_photon_right_recv(1:nm * nm * nx* nopphot) => g_lesser_photon_buf(nm * nm * nx * (local_NE + nopphot) + 1 : nm * nm * nx * (local_NE + 2 * nopphot))
  g_greater_photon_left_send(1:nm * nm * nx* nopphot) => g_greater_photon_buf(nm * nm * nx * nopphot + 1 : nm * nm * nx * 2 * nopphot)
  g_greater_photon_left_recv(1:nm * nm * nx* nopphot) => g_greater_photon_buf(1 : nm * nm * nx * nopphot)
  g_greater_photon_right_send(1:nm * nm * nx* nopphot) => g_greater_photon_buf(nm * nm * nx * local_NE + 1 : nm * nm * nx * (local_NE + nopphot))
  g_greater_photon_right_recv(1:nm * nm * nx* nopphot) => g_greater_photon_buf(nm * nm * nx * (local_NE + nopphot) + 1 : nm * nm * nx * (local_NE + 2 * nopphot))

  g_r_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_r_buf
  g_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_greater_buf
  g_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_lesser_buf

  g_r_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_r_buf
  g_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_greater_buf
  g_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_lesser_buf

  allocate(sigma_r_gw_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_greater_gw_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_lesser_gw_buf(nm * nm * nx * local_NE * nk))

  sigma_r_gw_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_r_gw_buf
  sigma_greater_gw_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_greater_gw_buf
  sigma_lesser_gw_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_lesser_gw_buf

  sigma_r_gw_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_r_gw_buf
  sigma_greater_gw_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_greater_gw_buf
  sigma_lesser_gw_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_lesser_gw_buf

  allocate(sigma_r_new_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_greater_new_buf(nm * nm * nx * local_NE * nk))
  allocate(sigma_lesser_new_buf(nm * nm * nx * local_NE * nk))

  sigma_r_new_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_r_new_buf
  sigma_greater_new_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_greater_new_buf
  sigma_lesser_new_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => sigma_lesser_new_buf

  sigma_r_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_r_new_buf
  sigma_greater_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_greater_new_buf
  sigma_lesser_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_lesser_new_buf

  allocate(P_retarded_buf(nm * nm * nx * local_NE * nk))
  allocate(P_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(P_lesser_buf(nm * nm * nx * local_NE * nk))

  P_retarded_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_retarded_buf
  P_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_greater_buf
  P_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_lesser_buf

  P_retarded_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_retarded_buf
  P_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_greater_buf
  P_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_lesser_buf

  allocate(W_retarded_buf(nm * nm * nx * local_NE * nk))  
  allocate(W_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(W_lesser_buf(nm * nm * nx * local_NE * nk))

  W_retarded_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_retarded_buf
  W_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_greater_buf
  W_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_lesser_buf

  W_retarded_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_retarded_buf
  W_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_greater_buf
  W_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_lesser_buf

  allocate(cur(nm, nm, nx, local_NE, nk))

  allocate(tmp0(nm * nm * nx * local_NE * nk))
  allocate(tmp1(nm * nm * nx * local_NE * nk))

  allocate(g_greater_t_buf(nm * nm * nx * local_NE * nk))
  allocate(g_lesser_t_buf(nm * nm * nx * local_NE * nk))

  g_greater_t_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_greater_t_buf
  g_lesser_t_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_lesser_t_buf

  g_greater_t_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_greater_t_buf
  g_lesser_t_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => g_lesser_t_buf

  sigma_greater_gw_buf = dcmplx(0.0d0,0.0d0)
  sigma_lesser_gw_buf = dcmplx(0.0d0,0.0d0)
  sigma_r_gw_buf = dcmplx(0.0d0,0.0d0)

  allocate(Mii(nm,nm,nx,nk))
  Mii(:,:,:,:)=czero
  do i=1,3
    Mii(:,:,:,:)=Mii(:,:,:,:)+ polarization(i) * Pii(:,:,i,:,:) 
  enddo
  if (comm_rank == 0) then
    print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity
  endif
  if (comm_rank == 0) then
    print *,'nop photon=',nopphot
  endif

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  do iter=0,niter
    if (comm_rank == 0) then
      print *,'+ iter=',iter
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    start = MPI_Wtime()
    it_start = start

    if (comm_rank == 0) then
      print *, 'Computing G ...'
    endif      
        
    allocate(jdens_local(nb,nb,nx*ns,local_NE,nk))
    allocate(tot_cur_local(nb,nb,nx*ns,nk))
    allocate(tot_ecur_local(nb,nb,nx*ns,nk))
    allocate(tot_cur(nb,nb,nx*ns,nk))
    allocate(tot_ecur(nb,nb,nx*ns,nk))
        
    do ik=1,nk
      if (comm_rank == 0) then
        print *, ' ik=', ik,'/',nk
      endif
      !$omp parallel default(shared) private(ie)
      !$omp do 
      do ie=1, local_NE     
        !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
        call RGF_G( &
          TEMP, nm, nx, local_energies(ie), mu, Hii(:,:,:,ik), H1i(:,:,:,ik), &
          sigma_lesser_gw_by_energies(:,:,:,ie,ik),sigma_greater_gw_by_energies(:,:,:,ie,ik), sigma_r_gw_by_energies(:,:,:,ie,ik), &
          g_lesser_by_energies(:,:,:,ie,ik), g_greater_by_energies(:,:,:,ie,ik), g_r_by_energies(:,:,:,ie,ik), &
          tre(ie,ik), tr(ie,ik), cur(:,:,:,ie,ik) ) 
      enddo
      !$omp end do 
      !$omp end parallel
      call calc_block_current(Hii(:,:,:,ik),g_lesser_by_energies(:,:,:,:,ik),cur(:,:,:,:,ik),local_NE,local_energies,spindeg,nb,ns,nm,nx,tot_cur_local(:,:,:,ik),tot_ecur_local(:,:,:,ik),jdens_local(:,:,:,:,ik))    
    enddo
    

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("G computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing G ...'
    endif
    start = finish

    do i = 0, comm_size - 1
      if (i == comm_rank) then
        filename = 'gw_ldos'
        call write_spectrum_summed_over_k(filename,iter,g_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/), append)
!        filename = 'gw_ndos'
!        call write_spectrum_summed_over_k(filename,iter,g_lesser_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_pdos'   
!        call write_spectrum_summed_over_k(filename,iter,g_greater_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
        filename = 'gw_Jdens'
        call write_current_spectrum_summed_over_kz(filename,iter,jdens_local,local_NE,local_energies,nx*NS,NB,Lx,nk, append)        
        filename = 'gw_trL'
        call write_transmission_spectrum_k(filename,iter,tr,local_NE,local_energies,nk, append)
        filename = 'gw_trR'
        call write_transmission_spectrum_k(filename,iter,tre,local_NE,local_energies,nk, append)
!        call write_dos_summed_over_k('gw_dos',iter,G_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx, append)
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo

    ! Reduce tr and tre
    local_sum_tr = sum(tr)
    local_sum_tre = sum(tre)
    call MPI_Reduce(local_sum_tr, global_sum_tr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(local_sum_tre, global_sum_tre, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! Reduce current and energy-current
    call MPI_Reduce(tot_cur_local, tot_cur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(tot_ecur_local, tot_ecur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (comm_rank == 0) then
      open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
      write(101,'(I4,2E16.6)') iter, global_sum_tr*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)/dble(nk), global_sum_tre*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)/dble(nk)
      close(101)
      call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
      call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)      
    endif
    
    deallocate(tot_cur,tot_ecur,jdens_local,tot_cur_local,tot_ecur_local)

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("G storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing P ...'
    endif
    start = finish

    g_r_buf = dcmplx( 0.0d0*dble(g_r_buf), aimag(g_r_buf))
    g_lesser_buf = dcmplx( 0.0d0*dble(g_lesser_buf), aimag(g_lesser_buf))
    g_greater_buf = dcmplx( 0.0d0*dble(g_greater_buf), aimag(g_greater_buf))

    g_greater_t_by_energies = reshape(g_greater_by_energies, shape(g_greater_t_by_energies), order=[2, 1, 3, 4, 5])
    g_lesser_t_by_energies = reshape(g_lesser_by_energies, shape(g_lesser_t_by_energies), order=[2, 1, 3, 4, 5])

    ! Redistribute G
    call energies_to_ijs_2(g_r_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_lesser_t_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(g_greater_t_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    P_lesser_buf = dcmplx(0.0d0,0.0d0)
    P_greater_buf = dcmplx(0.0d0,0.0d0)    
    P_retarded_buf = dcmplx(0.0d0,0.0d0)    

    nopmax=nen/2-1
    
    cmethod='fft'
    call MKL_SET_NUM_THREADS(4)

    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)             
    do i = 1, local_Nij
      global_ij = i + first_local_ij - 1
      col = (global_ij - 1) / nm + 1  ! convert to 0-based indexing, divide, and convert back to 1-based indexing
      row = mod(global_ij - 1, nm) + 1  ! convert to 0-based indexing, mod, and convert back to 1-based indexing
      l = max(row - ndiag, 1)
      h = min(nm, row + ndiag)
      if (col >= l .and. col <= h) then
      !$omp parallel default(shared) private(ix,iop,nop,ie,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
      !$omp do
      do ix = 1, nx
        do iqy=1,nky        
          do iqz=1,nkz
            iq=iqz + (iqy-1)*nkz
            do iky=1,nky
              do ikz=1,nkz              
                ik=ikz + (iky-1)*nkz
                ikzd=ikz-iqz + nkz/2
                ikyd=iky-iqy + nky/2
                if (ikzd<1)   ikzd=ikzd+nkz
                if (ikzd>nkz) ikzd=ikzd-nkz
                if (ikyd<1)   ikyd=ikyd+nky
                if (ikyd>nky) ikyd=ikyd-nky                
                if (nky==1)   ikyd=1
                if (nkz==1)   ikzd=1
                ikd=ikzd + (ikyd-1)*nkz
                
                P_lesser_by_blocks(:,iq,ix,i) = P_lesser_by_blocks(:,iq,ix,i) + &
                                                corr1d(nen,G_lesser_by_blocks(:,ik,ix,i),G_greater_t_by_blocks(:,ikd,ix,i),method=cmethod)
                
                P_greater_by_blocks(:,iq,ix,i) = P_greater_by_blocks(:,iq,ix,i) + & 
                                                corr1d(nen,G_greater_by_blocks(:,ik,ix,i),G_lesser_t_by_blocks(:,ikd,ix,i),method=cmethod)
                
                P_retarded_by_blocks(:,iq,ix,i) = P_retarded_by_blocks(:,iq,ix,i) + & 
                                                  corr1d(nen,G_lesser_by_blocks(:,ik,ix,i),conjg(G_r_by_blocks(:,ikd,ix,i)),method=cmethod) + &
                                                  corr1d(nen,G_r_by_blocks(:,ik,ix,i),G_lesser_t_by_blocks(:,ikd,ix,i),method=cmethod)
                              
                
              enddo
            enddo
          enddo
        enddo
      enddo
      !$omp end do 
      !$omp end parallel
    endif     
    enddo          

    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nk)
    P_lesser_buf = P_lesser_buf * dE
    P_greater_buf = P_greater_buf * dE
    P_retarded_buf = P_retarded_buf * dE

    ! Redistribute P
    call ijs_to_energies_2(P_retarded_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(P_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(P_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    ! Zero-out off diagonals
    !$omp parallel default(shared) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
    !$omp do
    do ik = 1, nk
      do ie = 1, local_NE
        do ix = 1, nx
          do i=1,nm      
            l = max(i-ndiag,1)
            h = min(nm,i+ndiag)
            if (l > 1) then
              P_lesser_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_greater_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_retarded_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
            endif
            if (h < nm) then
              P_lesser_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_greater_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              P_retarded_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0) 
            endif    
          enddo      
        enddo
      enddo
    enddo    
    !$omp end do 
    !$omp end parallel

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print *, 'ndiag=', ndiag
      print '("P computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing P ...'
    endif
    start = finish

!    do i = 0, comm_size - 1
!      if (i == comm_rank) then
!        filename = 'gw_PR'
!        call write_spectrum_summed_over_k(filename,iter,P_retarded_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_PL'
!        call write_spectrum_summed_over_k(filename,iter,P_lesser_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_PG'
!        call write_spectrum_summed_over_k(filename,iter,P_greater_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!      endif
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("P storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing W ...'
    endif
    start = finish

    do iq=1,nky*nkz

      if (comm_rank == 0) then
        print *, ' iq=', iq,'/',nk
      endif

      !$omp parallel default(shared) private(nop, ie)
      !$omp do 
      do nop = max(-nopmax + nen/2, first_local_energy), min(nopmax + nen/2, first_local_energy + local_NE - 1)
        ie = nop - first_local_energy + 1

        call green_calc_w_full( &
          0, nm, nx, Vii(:,:,:,iq), V1i(:,:,:,iq), &
          p_lesser_by_energies(:,:,:,ie,iq), p_greater_by_energies(:,:,:,ie,iq), p_retarded_by_energies(:,:,:,ie,iq), &
          w_lesser_by_energies(:,:,:,ie,iq), w_greater_by_energies(:,:,:,ie,iq), w_retarded_by_energies(:,:,:,ie,iq))

!        call green_rgf_calc_w( &
!          0, nm, nx, Vii(:,:,:,iq), V1i(:,:,:,iq), &
!          p_lesser_by_energies(:,:,:,ie,iq), p_greater_by_energies(:,:,:,ie,iq), p_retarded_by_energies(:,:,:,ie,iq), &
!          w_lesser_by_energies(:,:,:,ie,iq), w_greater_by_energies(:,:,:,ie,iq), w_retarded_by_energies(:,:,:,ie,iq))  
      enddo
      !$omp end do 
      !$omp end parallel
    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("W computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing W ...'
    endif
    start = finish

!    do i = 0, comm_size - 1
!      if (i == comm_rank) then
!        filename = 'gw_WR'
!        call write_spectrum_summed_over_k(filename,iter,W_retarded_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_WL'
!        call write_spectrum_summed_over_k(filename,iter,W_lesser_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_WG'
!        call write_spectrum_summed_over_k(filename,iter,W_greater_by_energies,local_NE,local_energies-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!      endif
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("W storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing SigGW ...'
    endif
    start = finish

    call energies_to_ijs_2(W_retarded_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(W_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call energies_to_ijs_2(W_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)    

    Sigma_greater_new_buf = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new_buf = dcmplx(0.0d0,0.0d0)
    Sigma_r_new_buf = dcmplx(0.0d0,0.0d0)
  
    cmethod='fft'
    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)          
    do i = 1, local_Nij
      global_ij = i + first_local_ij - 1
      col = (global_ij - 1) / nm + 1  ! convert to 0-based indexing, divide, and convert back to 1-based indexing
      row = mod(global_ij - 1, nm) + 1  ! convert to 0-based indexing, mod, and convert back to 1-based indexing
      l = max(row - ndiag, 1)
      h = min(nm, row + ndiag)
      if (col >= l .and. col <= h) then
      !$omp parallel default(shared) private(ix,iop,nop,ie,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
      !$omp do  
      do ix=1,nx     
        do iky=1,nky
          do ikz=1,nkz   
          
            ik=ikz+(iky-1)*nkz
            
            do iqy=1,nky
              do iqz=1,nkz              
              
                iq=iqz + (iqy-1)*nkz
                
                ikzd=ikz-iqz + nkz/2
                ikyd=iky-iqy + nky/2            
                if (ikzd<1)   ikzd=ikzd+nkz
                if (ikzd>nkz) ikzd=ikzd-nkz
                if (ikyd<1)   ikyd=ikyd+nky
                if (ikyd>nky) ikyd=ikyd-nky        
                if (nky==1)   ikyd=1
                if (nkz==1)   ikzd=1        
                
                ikd=ikzd + (ikyd-1)*nkz
                
                sigma_lesser_new_by_blocks(1:nen, ik, ix, i) = sigma_lesser_new_by_blocks(1:nen, ik, ix, i) + &
                                                    conv1d(nen,G_lesser_by_blocks(:, ikd, ix, i),W_lesser_by_blocks(:, iq, ix, i),method=cmethod)
                
                Sigma_greater_new_by_blocks(1:nen, ik, ix, i) = Sigma_greater_new_by_blocks(1:nen, ik, ix, i) + &
                                                    conv1d(nen,G_greater_by_blocks(:, ikd, ix, i),W_greater_by_blocks(:, iq, ix, i),method=cmethod)
                
                Sigma_r_new_by_blocks(1:nen, ik, ix, i) = Sigma_r_new_by_blocks(1:nen, ik, ix, i) + & 
                                                    conv1d(nen,G_lesser_by_blocks(:, ikd, ix, i), W_retarded_by_blocks(:, iq, ix, i),method=cmethod) + &
                                                    conv1d(nen,G_r_by_blocks(:, ikd, ix, i), W_lesser_by_blocks(:, iq, ix, i),method=cmethod) + &
                                                    conv1d(nen,G_r_by_blocks(:, ikd, ix, i), W_retarded_by_blocks(:, iq, ix, i) ,method=cmethod)
                
              
                
              enddo
            enddo            
          enddo
        enddo
      enddo
      !$omp end do
      !$omp end parallel
    endif      
    enddo

    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nk) 
    Sigma_lesser_new_buf = Sigma_lesser_new_buf  * dE
    Sigma_greater_new_buf = Sigma_greater_new_buf * dE
    Sigma_r_new_buf = Sigma_r_new_buf * dE
    Sigma_r_new_buf = dcmplx( dble(Sigma_r_new_buf), aimag(Sigma_greater_new_buf-Sigma_lesser_new_buf)/2.0d0 )

    call ijs_to_energies_2(sigma_r_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(sigma_lesser_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(sigma_greater_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    ! Zero-out off diagonals
    !$omp parallel default(shared) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
    !$omp do
    do ik = 1, nk
      do ie = 1, local_NE
        do ix = 1, nx
          do i=1,nm      
            l = max(i-ndiag,1)
            h = min(nm,i+ndiag)
            if (l > 1) then
              Sigma_lesser_new_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_greater_new_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_r_new_by_energies(i,1:l-1,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
            endif
            if (h < nm) then
              Sigma_lesser_new_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_greater_new_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0)
              Sigma_r_new_by_energies(i,h+1:nm,ix,ie,ik) = dcmplx(0.0d0, 0.0d0) 
            endif    
          enddo      
        enddo
      enddo
    enddo
    !$omp end do
    !$omp end parallel   
 
    ! symmetrize the selfenergies
    ! do ie=1,nen
    do ie = 1, local_NE
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new_by_energies(:,:,ix,ie,ik))
          Sigma_r_new_by_energies(:,:,ix,ie,ik) = (Sigma_r_new_by_energies(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new_by_energies(:,:,ix,ie,ik))
          Sigma_lesser_new_by_energies(:,:,ix,ie,ik) = (Sigma_lesser_new_by_energies(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new_by_energies(:,:,ix,ie,ik))
          Sigma_greater_new_by_energies(:,:,ix,ie,ik) = (Sigma_greater_new_by_energies(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo

    ! mixing with the previous one

    Sigma_r_gw_buf = Sigma_r_gw_buf + alpha_mix * (Sigma_r_new_buf - Sigma_r_gw_buf)
    Sigma_lesser_gw_buf  = Sigma_lesser_gw_buf + alpha_mix * (Sigma_lesser_new_buf - Sigma_lesser_gw_buf)
    Sigma_greater_gw_buf = Sigma_greater_gw_buf + alpha_mix * (Sigma_greater_new_buf -Sigma_greater_gw_buf)    
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("SigGW computation time = ", F0.3 ," seconds.")', finish-start
      print *, 'Storing SigGW ...'
    endif
    start = finish

!    do i = 0, comm_size - 1
!      if (i == comm_rank) then
!        filename = 'gw_SigR'
!        call write_spectrum_summed_over_k(filename,iter,Sigma_r_gw_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_SigL'
!        call write_spectrum_summed_over_k(filename,iter,Sigma_lesser_gw_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        filename = 'gw_SigG'
!        call write_spectrum_summed_over_k(filename,iter,Sigma_greater_gw_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!      endif
!      call MPI_Barrier(MPI_COMM_WORLD, ierr)
!    enddo
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("SigGW storage time = ", F0.3 ," seconds.")', finish-start
      print *, 'Computing GW scattering rates ...'
    endif
    start = finish
      
    ! We need to convert G back to energies     
    call ijs_to_energies_2(g_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
    call ijs_to_energies_2(g_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

    !!!!! calculate collision integral
    
    allocate(Ispec(nm,nm,nx,local_NE))
    allocate(Itot(nm,nm,nx))
    allocate(Ispec_ik(nm,nm,nx,local_NE))
    allocate(Itot_ik(nm,nm,nx))
    !
    Ispec=czero
    Itot=czero
    do ik=1,nk
      call calc_block_collision(sigma_lesser_gw_by_energies(:,:,:,:,ik),sigma_greater_gw_by_energies(:,:,:,:,ik),G_lesser_by_energies(:,:,:,:,ik),G_greater_by_energies(:,:,:,:,ik),local_NE,local_energies,spindeg,nm,nx,Itot_ik,Ispec_ik)
      Ispec=Ispec+Ispec_ik
      Itot=Itot+Itot_ik
    enddo
    
    do i = 0, comm_size - 1
      if (i == comm_rank) then
        filename = 'gw_Scat'
        call write_spectrum(filename,iter,Ispec,local_NE,local_energies,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/),append)    
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo    
        
    finish = MPI_Wtime()
    if (comm_rank == 0) then
      print '("Scattering computation time = ", F0.3 ," seconds.")', finish-start      
    endif
    
    if (iter>=5) then

      if (comm_rank == 0) then
        print *, 'Computing SigEPhoton ...'
      endif
      start = finish

      ! call energies_to_ijs(sigma_r_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call energies_to_ijs(sigma_lesser_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call energies_to_ijs(sigma_greater_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

      ! call ijs_to_energies_2(g_r_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call ijs_to_energies_2(g_lesser_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call ijs_to_energies_2(g_greater_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

      do ik=1, nk

        write (rank_str, fmt) ik

        if (comm_rank == 0) then
          print *, ' ik=', ik,'/',nk
        endif

        ! do i = 0, comm_size - 1
        !   if (i == comm_rank) then
        !     filename = 'gw_GL_ik' // rank_str // '_'
        !     call write_spectrum_summed_over_k(filename,iter,g_lesser_by_energies(:,:,:,:,ik:ik),local_NE,local_energies,1,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
        !     filename = 'gw_GG_ik' // rank_str // '_'  
        !     call write_spectrum_summed_over_k(filename,iter,g_greater_by_energies(:,:,:,:,ik:ik),local_NE,local_energies,1,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
        !   endif
        !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
        ! enddo

        ! Copy local_NE energies to the buffers
        g_lesser_photon_buf = dcmplx(0.0d0,0.0d0)
        g_greater_photon_buf = dcmplx(0.0d0,0.0d0)
        g_lesser_photon(:, :, :, nopphot + 1 : nopphot + local_NE) = g_lesser_by_energies(:, :, :, :, ik)
        g_greater_photon(:, :, :, nopphot + 1 : nopphot + local_NE) = g_greater_by_energies(:, :, :, :, ik)

        ! Communicate up to 2*nopphot energies to the left and right
        if (comm_rank - 1 >= 0) then
          call MPI_Isend(g_lesser_photon_left_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 0, MPI_COMM_WORLD, reqs(1), ierr)
          call MPI_Irecv(g_lesser_photon_left_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 1, MPI_COMM_WORLD, reqs(2), ierr)
          call MPI_Isend(g_greater_photon_left_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 2, MPI_COMM_WORLD, reqs(3), ierr)
          call MPI_Irecv(g_greater_photon_left_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank-1, 3, MPI_COMM_WORLD, reqs(4), ierr)
        endif
        if (comm_rank + 1 < comm_size) then
          call MPI_Isend(g_lesser_photon_right_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 1, MPI_COMM_WORLD, reqs(5), ierr)
          call MPI_Irecv(g_lesser_photon_right_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 0, MPI_COMM_WORLD, reqs(6), ierr)
          call MPI_Isend(g_greater_photon_right_send, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 3, MPI_COMM_WORLD, reqs(7), ierr)
          call MPI_Irecv(g_greater_photon_right_recv, nm*nm*nx*nopphot, MPI_DOUBLE_COMPLEX, comm_rank+1, 2, MPI_COMM_WORLD, reqs(8), ierr)
        endif
  
        if (comm_rank - 1 >= 0) then
          call MPI_Waitall(4, reqs(1:4), MPI_STATUSES_IGNORE, ierr)
        endif
        if (comm_rank +1 < comm_size) then
          call MPI_Waitall(4, reqs(5:8), MPI_STATUSES_IGNORE, ierr)
        endif

        ! do i = 0, comm_size - 1
        !   if (i == comm_rank) then
        !     filename = 'gw_GL_ext_ik' // rank_str // '_'
        !     call write_spectrum_summed_over_k(filename,iter,g_lesser_photon,local_NE + 2 * nopphot,extended_local_energies,1,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
        !     filename = 'gw_GG_ext_ik' // rank_str // '_'   
        !     call write_spectrum_summed_over_k(filename,iter,g_greater_photon,local_NE + 2 * nopphot,extended_local_energies,1,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
        !   endif
        !   call MPI_Barrier(MPI_COMM_WORLD, ierr)
        ! enddo

        call calc_sigma_ephoton_monochromatic_ext( &
          nm, NX, local_NE, En, nopphot, Mii(:,:,:,ik), &
          g_lesser_photon(:,:,:,:), g_greater_photon(:,:,:,:), &
          Sigma_r_new_by_energies(:,:,:,:,ik), Sigma_lesser_new_by_energies(:,:,:,:,ik), Sigma_greater_new_by_energies(:,:,:,:,ik), &
          pre_fact, intensity, hw)
      enddo

      ! call ijs_to_energies(sigma_r_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call ijs_to_energies(sigma_lesser_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)
      ! call ijs_to_energies(sigma_greater_new_buf, tmp0, tmp1, nm, nx, nen, nk, local_Nij, local_NE, first_local_ij, first_local_energy, comm_rank, comm_size)

      Sigma_r_gw_buf       = Sigma_r_gw_buf + Sigma_r_new_buf 
      Sigma_lesser_gw_buf  = Sigma_lesser_gw_buf + Sigma_lesser_new_buf 
      Sigma_greater_gw_buf = Sigma_greater_gw_buf + Sigma_greater_new_buf

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      finish = MPI_Wtime()

      if (comm_rank == 0) then
        print '("SigEPhoton computation time = ", F0.3 ," seconds.")', finish-start
        print *, 'Storing SigEphoton ...'
      endif
      start = finish

!      do i = 0, comm_size - 1
!        if (i == comm_rank) then
!          filename = 'eph_SigR_'
!          call write_spectrum_summed_over_k(filename,iter,sigma_r_new_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!          filename = 'eph_SigL'
!          call write_spectrum_summed_over_k(filename,iter,sigma_lesser_new_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!          filename = 'eph_SigG'
!          call write_spectrum_summed_over_k(filename,iter,sigma_greater_new_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
!        endif
!        call MPI_Barrier(MPI_COMM_WORLD, ierr)
!      enddo

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      finish = MPI_Wtime()

      if (comm_rank == 0) then
        print '("SigEPhoton storage time = ", F0.3 ," seconds.")', finish-start
        print *, 'Computing Ephoton scattering rates ...'
      endif
      
      !!!! calculate e-photon scattering rates
            
      start = finish
      
      Ispec=czero
      Itot=czero
      do ik=1,nk
        call calc_block_collision(sigma_lesser_new_by_energies(:,:,:,:,ik),sigma_greater_new_by_energies(:,:,:,:,ik),G_lesser_by_energies(:,:,:,:,ik),G_greater_by_energies(:,:,:,:,ik),local_NE,local_energies,spindeg,nm,nx,Itot_ik,Ispec_ik)
        Ispec=Ispec+Ispec_ik
        Itot=Itot+Itot_ik
      enddo
      do i = 0, comm_size - 1
        if (i == comm_rank) then
          filename = 'eph_Scat'
          call write_spectrum(filename,iter,Ispec,local_NE,local_energies,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/),append)    
        endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
      enddo                 
      
      finish = MPI_Wtime()

      if (comm_rank == 0) then
        print '("EPhoton scattering time = ", F0.3 ," seconds.")', finish-start        
      endif 

    endif
    
    deallocate(Ispec,Itot,Ispec_ik,Itot_ik) 
    
    ! make sure self-energy is continuous near leads (by copying edge block)
    ! do ix=1,2
    !   sigma_r_gw_by_energies(:,:,ix,:,:)=Sigma_r_gw_by_energies(:,:,3,:,:)
    !   sigma_lesser_gw_by_energies(:,:,ix,:,:)=Sigma_lesser_gw_by_energies(:,:,3,:,:)
    !   sigma_greater_gw_by_energies(:,:,ix,:,:)=Sigma_greater_gw_by_energies(:,:,3,:,:)
    ! enddo
    ! do ix=1,2
    !   Sigma_r_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_r_gw_by_energies(:,:,nx-2,:,:)
    !   Sigma_lesser_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_lesser_gw_by_energies(:,:,nx-2,:,:)
    !   Sigma_greater_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_greater_gw_by_energies(:,:,nx-2,:,:)
    ! enddo

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    finish = MPI_Wtime()

    if (comm_rank == 0) then
      print '("Iteration ", I3.3, " time = ", F0.3 ," seconds.")', iter, finish-it_start
      print *, ''
    endif

  enddo
  
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  start = MPI_Wtime()
  it_start = start

  if (comm_rank == 0) then
    print *, 'Computing G ...'
  endif      
      
  allocate(jdens_local(nb,nb,nx*ns,local_NE,nk))
  allocate(tot_cur_local(nb,nb,nx*ns,nk))
  allocate(tot_ecur_local(nb,nb,nx*ns,nk))
  allocate(tot_cur(nb,nb,nx*ns,nk))
  allocate(tot_ecur(nb,nb,nx*ns,nk))
      
  do ik=1,nk
    if (comm_rank == 0) then
      print *, ' ik=', ik,'/',nk
    endif
    !$omp parallel default(shared) private(ie)
    !$omp do 
    do ie=1, local_NE     
      !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
      call green_RGF_RS( &
        TEMP, nm, nx, local_energies(ie), mu, Hii(:,:,:,ik), H1i(:,:,:,ik), &
        sigma_lesser_gw_by_energies(:,:,:,ie,ik),sigma_greater_gw_by_energies(:,:,:,ie,ik), sigma_r_gw_by_energies(:,:,:,ie,ik), &
        g_lesser_by_energies(:,:,:,ie,ik), g_greater_by_energies(:,:,:,ie,ik), g_r_by_energies(:,:,:,ie,ik), &
        tre(ie,ik), tr(ie,ik), cur(:,:,:,ie,ik) ) 
    enddo
    !$omp end do 
    !$omp end parallel
    call calc_block_current(Hii(:,:,:,ik),g_lesser_by_energies(:,:,:,:,ik),cur(:,:,:,:,ik),local_NE,local_energies,spindeg,nb,ns,nm,nx,tot_cur_local(:,:,:,ik),tot_ecur_local(:,:,:,ik),jdens_local(:,:,:,:,ik))    
  enddo
  

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  finish = MPI_Wtime()

  if (comm_rank == 0) then
    print '("G computation time = ", F0.3 ," seconds.")', finish-start
    print *, 'Storing G ...'
  endif
  start = finish

  do i = 0, comm_size - 1
    if (i == comm_rank) then
      filename = 'gw_ldos'
      call write_spectrum_summed_over_k(filename,iter,g_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/), append)
      filename = 'gw_ndos'
      call write_spectrum_summed_over_k(filename,iter,g_lesser_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/), append)
      filename = 'gw_pdos'   
      call write_spectrum_summed_over_k(filename,iter,g_greater_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/), append)
      filename = 'gw_Jdens'
      call write_current_spectrum_summed_over_kz(filename,iter,jdens_local,local_NE,local_energies,nx*NS,NB,Lx,nk, append)        
      filename = 'gw_trL'
      call write_transmission_spectrum_k(filename,iter,tr*spindeg,local_NE,local_energies,nk, append)
      filename = 'gw_trR'
      call write_transmission_spectrum_k(filename,iter,tre*spindeg,local_NE,local_energies,nk, append)
      call write_dos_summed_over_k('gw_dos',iter,G_r_by_energies,local_NE,local_energies,nk,nx,NB,NS,Lx, append)
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  enddo

  ! Reduce tr and tre
  local_sum_tr = sum(tr)
  local_sum_tre = sum(tre)
  call MPI_Reduce(local_sum_tr, global_sum_tr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(local_sum_tre, global_sum_tre, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  ! Reduce current and energy-current
  call MPI_Reduce(tot_cur_local, tot_cur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call MPI_Reduce(tot_ecur_local, tot_ecur, nb*nb*nx*ns*nk, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (comm_rank == 0) then
    open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, global_sum_tr*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)/dble(nk), global_sum_tre*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)/dble(nk)
    close(101)
    call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
    call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)      
  endif
  
  deallocate(tot_cur,tot_ecur,jdens_local,tot_cur_local,tot_ecur_local)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  finish = MPI_Wtime()

  if (comm_rank == 0) then
    print '("G storage time = ", F0.3 ," seconds.")', finish-start    
  endif

  deallocate( cur, tmp0, tmp1)
  ! deallocate(g_r_buf, g_lesser_extended_buf, g_greater_extended_buf)
  deallocate(g_r_buf, g_lesser_buf, g_greater_buf)
  deallocate(g_lesser_photon_buf, g_greater_photon_buf)
  deallocate(sigma_lesser_gw_buf, sigma_greater_gw_buf, sigma_r_gw_buf)
  deallocate(sigma_lesser_new_buf, sigma_greater_new_buf, sigma_r_new_buf)
  deallocate(P_retarded_buf, P_lesser_buf, P_greater_buf)
  deallocate(W_retarded_buf, W_lesser_buf, W_greater_buf)
  deallocate(Mii)

  deallocate(extended_local_energies)
end subroutine green_rgf_solve_gw_ephoton_3d_mpi




end module gw_sparse
