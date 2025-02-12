!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao, Alexandros Nikolaos Ziogas
!!!!!!!!!!!!!!!! DATE: 02/2023

module green_rgf

use green,only:get_OBC_blocks_for_W,get_dL_OBC_for_W

implicit none 

private

public :: green_rgf_cms, green_rgf_solve_gw_1d !, green_rgf_calc_w
public :: green_rgf_solve_gw_3d
public :: green_rgf_solve_gw_ephoton_3d, green_rgf_solve_gw_ephoton_3d_mpi_kpar, green_rgf_solve_gw_ephoton_3d_ijs

complex(8), parameter :: cone = dcmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = dcmplx(0.0d0,0.0d0)
real(8), parameter :: hbar=1.0546d-34 ! m^2 kg / s
real(8), parameter :: m0=9.109d-31 ! kg
real(8), parameter :: eps0=8.854d-12 ! C/V/m 
real(8), parameter :: c0=2.998d8 ! m/s
real(8), parameter :: e0=1.6022d-19 ! C
REAL(8), PARAMETER :: pi = 3.14159265359d0
REAL(8), PARAMETER :: tpi = 3.14159265359d0*2.0d0

include "mpif.h"

CONTAINS



! calculate e-photon polarization self-energies in the monochromatic assumption, for independent electron-hole pair
subroutine calc_pi_ephoton_monochromatic(nm,length,nen,En,nop,M,G_lesser,G_greater,Pi_retarded,Pi_lesser,Pi_greater)
integer,intent(in)::nm,length,nen,nop
real(8),intent(in)::en(nen)
complex(8),intent(in),dimension(nm,nm,length)::M ! e-photon interaction matrix
complex(8),intent(in),dimension(nm,nm,length,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm,nm,length)::Pi_retarded,Pi_lesser,Pi_greater
!---------
integer::ie,ix
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
complex(8)::dE
  Pi_lesser=0.0d0
  Pi_greater=0.0d0
  Pi_retarded=0.0d0  
  dE= dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi ) 
  ! Pi^<>(hw) = \Sum_E M G^<>(E) M G^><(E - hw) M
  !$omp parallel default(shared) private(ix,ie,A,B) 
  allocate(B(nm,nm))
  allocate(A(nm,nm))  
  !$omp do
  do ix=1,length
    do ie=1,nen
      if ((ie-nop>=1).and.(ie-nop<=nen)) then
        ! Pi^<(hw) = \sum_E M G<(E) M G>(E-hw)
        call zgemm('n','n',nm,nm,nm,cone,M,nm,G_lesser(:,:,ix,ie),nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,M(:,:,ix),nm,czero,A,nm) 
        call zgemm('n','n',nm,nm,nm,cone,A,nm,G_greater(:,:,ix,ie-nop),nm,cone,Pi_lesser(:,:,ix),nm)         
        ! Pi^>(hw) = \sum_E M G>(E) M G<(E-hw)   
        call zgemm('n','n',nm,nm,nm,cone,M,nm,G_greater(:,:,ix,ie),nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,M(:,:,ix),nm,czero,A,nm) 
        call zgemm('n','n',nm,nm,nm,cone,A,nm,G_lesser(:,:,ix,ie-nop),nm,cone,Pi_greater(:,:,ix),nm)             
      endif
    enddo  
  enddo
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Pi_lesser=Pi_lesser*dE
  Pi_greater=Pi_greater*dE
  Pi_retarded = dcmplx(0.0d0*dble(Pi_retarded),aimag(Pi_greater-Pi_lesser)/2.0d0)
end subroutine calc_pi_ephoton_monochromatic


subroutine energies_to_blocks(buf, tmp0, tmp1, NM, NX, NE, NK, local_NX, local_NE, first_local_NX, first_local_NE, comm_rank, comm_size)
  complex(8), target, intent (inout) :: buf(:), tmp0(:), tmp1(:)
  integer(kind = 4), intent ( in ) :: NM, NX, NE, NK, local_NX, local_NE, first_local_NX, first_local_NE, comm_rank, comm_size

  complex(8), pointer :: p0(:, :, :, :, :), p1(:, :, :, :, :), p2(:, :, :, :, :, :), p3(:, :, :, :, :, :)
  integer(kind = 4) :: count, ierr, i, j, ix, ie, ik, r, src_idx, dst_idx

  character (len = 8) :: fmt
  character (len = 20) :: filename
  character (len = 4) :: rank_str

  fmt = '(I4.4)'
  write ( rank_str, fmt ) comm_rank

  ! Assume (global) g(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))

  ! 1. Tranpose to g(1:NM, 1:NM, 1:NK, 1:local_NE, 1:NX, 1:comm_size (NE))
  ! !$omp parallel do default(shared) private(i, j, ix, ie, ik, src_idx, dst_idx)
  ! do ie=1,local_ne
  !   do ix=1,nx
  !     do ik=1,nk
  !       do j=1,nm
  !         do i=1,nm
  !           src_idx = 1 + (i-1) + (j-1) * nm + (ix-1) * nm * nm + (ie-1) * nm * nm * nx + (ik-1) * nm * nm * nx * local_ne
  !           dst_idx = 1 + (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-1) * nm * nm * nk + (ix-1) * nm * nm * nk * local_ne
  !           tmp0(dst_idx) = buf(src_idx)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! !$omp end parallel do
  p0(1:nm, 1:nm, 1:nx, 1:local_ne, 1:nk) => buf
  p1(1:nm, 1:nm, 1:nk, 1:local_ne, 1:nx) => tmp0
  p1 = reshape(p0, shape(p1), order = [1, 2, 5, 4, 3])

  ! filename = 'G_reshaped_buf_r' // rank_str // '_'
  ! open(unit=11, file=trim(filename)//'.dat', status='unknown')
  ! do ix=1,nx
  !   do ie=1,local_ne
  !     do ik=1,nk
  !       do j=1,nm
  !         do i=1,nm
  !           src_idx = 1 + (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-1) * nm * nm * nk + (ix-1) * nm * nm * nk * local_ne
  !           write(11, *) ik, ie-1+first_local_ne, ix, j, i, tmp0(src_idx)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(11)

  ! 2a. Interpret as g(1:NM, 1:NM, 1:NK, 1:local_NE, 1:local_NX, 1:comm_size (NX), 1:comm_size (NE))
  ! 2b. Redistribute to g(1:NM, 1:NM, 1:NK, 1:local_NE, 1:local_NX, 1:comm_size (NE), 1:comm_size (NX))
  count = NM * NM * NK * local_NE * local_NX
  call MPI_Alltoall(tmp0, count, MPI_DOUBLE_COMPLEX, tmp1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! filename = 'G_alltoall_buf_r' // rank_str // '_'
  ! open(unit=11, file=trim(filename)//'.dat', status='unknown')
  ! do ix=1,local_nx
  !   do ie=1,ne
  !     do ik=1,nk
  !       do j=1,nm
  !         do i=1,nm
  !           src_idx = 1 + (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-1) * nm * nm * nk + (ix-1) * nm * nm * nk * ne
  !           write(11, *) ik, ie, ix-1+first_local_nx, j, i, tmp1(src_idx)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(11)

  ! 3. Tranpose to g(1:NM, 1:NM, 1:local_NX, 1:local_NE, 1:comm_size (NE), 1:NK, 1:comm_size (NX))
  ! buf = dcmplx(0.0d0, 0.0d0)
  ! !$omp parallel do default(shared) private(i, j, ix, ie, ik, r, src_idx, dst_idx)
  ! do ie=1,local_ne
  !   do r=1,comm_size
  !     do ix=1,local_nx
  !       do ik=1,nk
  !         do j=1,nm
  !           do i=1,nm
  !             dst_idx = 1 + (i-1) + (j-1) * nm + (ix-1) * nm * nm + (ie-1) * nm * nm * local_nx + (r-1) * nm * nm * local_nx * local_ne + (ik-1) * nm * nm * local_nx * local_ne * comm_size
  !             src_idx = 1 + (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-1) * nm * nm * nk + (ix-1) * nm * nm * nk * local_ne + (r-1) * nm * nm * nk * local_ne * local_nx
  !             buf(dst_idx) = tmp1(src_idx)
  !           enddo
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! !$omp end parallel do
  p2(1:nm, 1:nm, 1:nk, 1:local_ne, 1:local_nx, 1:comm_size) => tmp1
  p3(1:nm, 1:nm, 1:local_nx, 1:local_ne, 1:comm_size, 1:nk) => buf
  p3 = reshape(p2, shape(p3), order = [1, 2, 6, 4, 3, 5])

end subroutine energies_to_blocks


subroutine blocks_to_energies(buf, tmp0, tmp1, NM, NX, NE, NK, local_NX, local_NE, first_local_NX, first_local_NE, comm_rank, comm_size)
  complex(8), target, intent (inout) :: buf(:), tmp0(:), tmp1(:)
  integer(kind = 4), intent ( in ) :: NM, NX, NE, NK, local_NX, local_NE, first_local_NX, first_local_NE, comm_rank, comm_size

  complex(8), pointer :: p0(:, :, :, :, :), p1(:, :, :, :, :), p2(:, :, :, :, :, :), p3(:, :, :, :, :, :)
  integer(kind = 4) :: count, ierr, i, j, ix, ie, ik, r, src_idx, dst_idx

  character (len = 8) :: fmt
  character (len = 20) :: filename
  character (len = 4) :: rank_str

  fmt = '(I4.4)'
  write ( rank_str, fmt ) comm_rank

  ! Assume (global) g(1:NM, 1:NM, 1:local_NX, 1:NE, 1:NK, 1:comm_size (NX))

  ! 1. Tranpose to g(1:NM, 1:NM, 1:local_NX, 1:NK, 1:NE, 1:comm_size (NX))
  p0(1:nm, 1:nm, 1:local_nx, 1:ne, 1:nk) => buf
  p1(1:nm, 1:nm, 1:local_nx, 1:nk, 1:ne) => tmp0
  p1 = reshape(p0, shape(p1), order = [1, 2, 3, 5, 4])

  ! 2a. Interpret as g(1:NM, 1:NM, 1:local_NX, 1:NK, 1:local_NE, 1:comm_size (NE), 1:comm_size (NX))
  ! 2b. Redistribute to g(1:NM, 1:NM, 1:local_NX, 1:NK, 1:local_NE, 1:comm_size (NX), 1:comm_size (NE))
  count = NM * NM * local_NX * NK * local_NE
  call MPI_Alltoall(tmp0, count, MPI_DOUBLE_COMPLEX, tmp1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! 3. Tranpose to g(1:NM, 1:NM, 1:local_NX, 1:comm_size(NX), 1:local_NE, 1:NK, 1:comm_size (NE))
  p2(1:nm, 1:nm, 1:local_nx, 1:nk, 1:local_ne, 1:comm_size) => tmp1
  p3(1:nm, 1:nm, 1:local_nx, 1:comm_size, 1:local_ne, 1:nk) => buf
  p3 = reshape(p2, shape(p3), order = [1, 2, 3, 6, 5, 4])

end subroutine blocks_to_energies


subroutine energies_to_ijs(buf, tmp0, tmp1, NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size)
  complex(8), target, intent (inout) :: buf(:), tmp0(:), tmp1(:)
  integer(kind = 4), intent ( in ) :: NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size

  complex(8), pointer :: p0(:, :, :, :, :), p1(:, :, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  ! complex(8), pointer :: p0(:, :, :, :), p1(:, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  integer(kind = 4) :: count, ierr, i, j, ix, ie, ik, r, src_idx, dst_idx

  ! Assume (global) g(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))

  ! 1a. Interpret as g(1:NM*NM, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  ! 1b. Interpret as g(1:local_Nij, 1:comm_size (Nij), 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  ! 1c. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (Nij), 1:comm_size (NE))
  p0(1:local_nij, 1:comm_size, 1:nx, 1:local_ne, 1:nk) => buf
  p1(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp0
  p1 = reshape(p0, shape(p1), order = [1, 5, 2, 3, 4])
  ! p0(1:nm*nm, 1:nx, 1:local_ne, 1:nk) => buf
  ! p1(1:nx, 1:local_ne, 1:nk, 1:nm*nm) => tmp0
  ! p1 = reshape(p0, shape(p1), order = [4, 1, 2, 3])

  ! 2. Redistribute to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE), 1:comm_size (Nij))
  count = local_nij * nx * local_ne * nk
  call MPI_Alltoall(tmp0, count, MPI_DOUBLE_COMPLEX, tmp1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! 3. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:comm_size (NE), 1:NK, 1:comm_size (Nij))
  p2(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp1
  p3(1:local_nij, 1:nx, 1:local_ne, 1:comm_size, 1:nk) => buf
  p3 = reshape(p2, shape(p3), order = [1, 2, 3, 5, 4])
  ! p2(1:nx, 1:local_ne, 1:nk, 1:local_nij, 1:comm_size) => tmp1
  ! p3(1:local_nij, 1:nx, 1:local_ne, 1:comm_size, 1:nk) => buf
  ! p3 = reshape(p2, shape(p3), order = [2, 3, 5, 1, 4])

end subroutine energies_to_ijs


subroutine energies_to_ijs_2(buf, tmp0, tmp1, NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size)
  complex(8), target, intent (inout) :: buf(:), tmp0(:), tmp1(:)
  integer(kind = 4), intent ( in ) :: NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size

  complex(8), pointer :: p0(:, :, :, :, :), p1(:, :, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  ! complex(8), pointer :: p0(:, :, :, :), p1(:, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  integer(kind = 4) :: count, ierr, i, j, ix, ie, ik, r, src_idx, dst_idx

  ! Assume (global) g(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))

  ! 1a. Interpret as g(1:NM*NM, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  ! 1b. Interpret as g(1:local_Nij, 1:comm_size (Nij), 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  ! 1c. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (Nij), 1:comm_size (NE))
  p0(1:local_nij, 1:comm_size, 1:nx, 1:local_ne, 1:nk) => buf
  p1(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp0
  p1 = reshape(p0, shape(p1), order = [1, 5, 2, 3, 4])
  ! p0(1:nm*nm, 1:nx, 1:local_ne, 1:nk) => buf
  ! p1(1:nx, 1:local_ne, 1:nk, 1:nm*nm) => tmp0
  ! p1 = reshape(p0, shape(p1), order = [4, 1, 2, 3])

  ! 2. Redistribute to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE), 1:comm_size (Nij))
  count = local_nij * nx * local_ne * nk
  call MPI_Alltoall(tmp0, count, MPI_DOUBLE_COMPLEX, tmp1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! 3. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:comm_size (NE), 1:NK, 1:comm_size (Nij))
  p2(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp1
  p3(1:local_ne, 1:comm_size, 1:nk, 1:nx, 1:local_nij) => buf
  p3 = reshape(p2, shape(p3), order = [5, 4, 1, 3, 2])
  ! p2(1:nx, 1:local_ne, 1:nk, 1:local_nij, 1:comm_size) => tmp1
  ! p3(1:local_nij, 1:nx, 1:local_ne, 1:comm_size, 1:nk) => buf
  ! p3 = reshape(p2, shape(p3), order = [2, 3, 5, 1, 4])

end subroutine energies_to_ijs_2


subroutine ijs_to_energies(buf, tmp0, tmp1, NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size)
  complex(8), target, intent (inout) :: buf(:), tmp0(:), tmp1(:)
  integer(kind = 4), intent ( in ) :: NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size

  complex(8), pointer :: p0(:, :, :, :, :), p1(:, :, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  integer(kind = 4) :: count, ierr, i, j, ix, ie, ik, r, src_idx, dst_idx

  ! Assume (global) g(1:local_Nij, 1:NX, 1:NE, 1:NK, 1:comm_size (Nij))

  ! 1a. Interpret as g(1:local_Nij, 1:NX, 1:local_NE, 1:comm_size (NE), 1:NK, 1:comm_size (Nij))
  ! 1b. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE), 1:comm_size (Nij))
  p0(1:local_nij, 1:nx, 1:local_ne, 1:comm_size, 1:nk) => buf
  p1(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp0
  p1 = reshape(p0, shape(p1), order = [1, 2, 3, 5, 4])

  ! 2. Redistribute to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (Nij), 1:comm_size (NE))
  count = local_Nij * nx * local_NE * nk
  call MPI_Alltoall(tmp0, count, MPI_DOUBLE_COMPLEX, tmp1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! 3. Transpose to g(1:local_Nij, 1:comm_size (Nij), 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  p2(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp1
  p3(1:local_nij, 1:comm_size, 1:nx, 1:local_ne, 1:nk) => buf
  p3 = reshape(p2, shape(p3), order = [1, 3, 4, 5, 2])

end subroutine ijs_to_energies


subroutine ijs_to_energies_2(buf, tmp0, tmp1, NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size)
  complex(8), target, intent (inout) :: buf(:), tmp0(:), tmp1(:)
  integer(kind = 4), intent ( in ) :: NM, NX, NE, NK, local_Nij, local_NE, first_local_Nij, first_local_NE, comm_rank, comm_size

  complex(8), pointer :: p0(:, :, :, :, :), p1(:, :, :, :, :), p2(:, :, :, :, :), p3(:, :, :, :, :)
  integer(kind = 4) :: count, ierr, i, j, ix, ie, ik, r, src_idx, dst_idx

  ! Assume (global) g(1:local_Nij, 1:NX, 1:NE, 1:NK, 1:comm_size (Nij))

  ! 1a. Interpret as g(1:local_Nij, 1:NX, 1:local_NE, 1:comm_size (NE), 1:NK, 1:comm_size (Nij))
  ! 1b. Transpose to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE), 1:comm_size (Nij))
  p0(1:local_ne, 1:comm_size, 1:nk, 1:nx, 1:local_nij) => buf
  p1(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp0
  p1 = reshape(p0, shape(p1), order = [3, 5, 4, 2, 1])

  ! 2. Redistribute to g(1:local_Nij, 1:NX, 1:local_NE, 1:NK, 1:comm_size (Nij), 1:comm_size (NE))
  count = local_Nij * nx * local_NE * nk
  call MPI_Alltoall(tmp0, count, MPI_DOUBLE_COMPLEX, tmp1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! 3. Transpose to g(1:local_Nij, 1:comm_size (Nij), 1:NX, 1:local_NE, 1:NK, 1:comm_size (NE))
  p2(1:local_nij, 1:nx, 1:local_ne, 1:nk, 1:comm_size) => tmp1
  p3(1:local_nij, 1:comm_size, 1:nx, 1:local_ne, 1:nk) => buf
  p3 = reshape(p2, shape(p3), order = [1, 3, 4, 5, 2])

end subroutine ijs_to_energies_2


! Redistribute g(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK) to g(1:NM, 1:NM, 1:NK, 1:NE, 1:local_NX)
subroutine g2p(g, p, NM, NX, NE, NK, local_NX, local_NE, first_local_NX, first_local_NE, comm_rank)
  complex(8), intent (in) :: g(:, :, :, :, :)
  complex(8), intent (inout) :: p(:)
  complex(8), allocatable :: buf(:)
  integer ( kind = 4 ), intent ( in ) :: NM, NX, NE, NK, local_NX, local_NE, first_local_NX, first_local_NE, comm_rank
  ! complex ( 8 ), allocatable :: buf(:, :, :, :, :)
  complex ( 8 ), allocatable :: buf0(:), buf1(:)
  integer ( kind = 4 ) :: count, ierr, i, j, ix, ie, ik, idx
  character ( len = 8 ) :: fmt
  character ( len = 20 ) :: filename
  character ( len = 4 ) :: rank_str

  ! Assume g(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK)

  ! 1. Tranpose to g(1:NM, 1:NM, 1:NK, 1:local_NE, 1:NX)
  allocate(buf(NM * NM * NK * local_NE * NX))
  call transpose_g(g, buf, NM, NX, NE, NK, local_NE)

  ! 2. Redistribute to g(1:NM, 1:NM, 1:NK, 1:NE, 1:local_NX)
  count = NM * NM * NK * local_NE * local_NX
  call MPI_Alltoall(buf, count, MPI_DOUBLE_COMPLEX, p, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  deallocate(buf)


  ! fmt = '(I4.4)'
  ! write ( rank_str, fmt ) comm_rank

  ! ! Assume g(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK)

  ! ! 1. Reshape to g(1:NM, 1:NM, 1:NK, 1:local_NE, 1:NX)
  ! ! allocate(buf(NM, NM, NK, local_NE, NX))
  ! allocate(buf0(NM * NM * NK * local_NE * NX))
  ! ! buf = reshape(g, shape(buf), order = [1, 2, 5, 4, 3])
  ! do ix=1,nx
  !   do ie=first_local_ne,first_local_ne + local_ne - 1
  !     do ik=1,nk
  !       do j=1,nm
  !         do i=1,nm
  !           idx = (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-first_local_NE) * nm * nm * nk + (ix-1) * nm * nm * nk * local_NE
  !           buf0(idx) = g(i,j,ix,ie-first_local_NE,ik)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo

  ! filename = 'G_reshaped_buf_r' // rank_str // '_'
  ! open(unit=11, file=trim(filename)//'.dat', status='unknown')
  ! ! do ik=1,nk
  ! !   do ie=first_local_ne,first_local_ne + local_ne - 1
  ! !     do ix=1,nx
  ! !       do j=1,nm
  ! !         do i=1,nm
  ! !           write(11, *) ik, ie, ix, j, i, buf(i,j,ik,ie-first_local_ne,ix)
  ! !         enddo
  ! !       enddo
  ! !     enddo
  ! !   enddo
  ! ! enddo
  ! do ix=1,nx
  !   do ie=first_local_ne,first_local_ne + local_ne - 1
  !     do ik=1,nk
  !       do j=1,nm
  !         do i=1,nm
  !           idx = (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-first_local_NE) * nm * nm * nk + (ix-1) * nm * nm * nk * local_NE
  !           write(11, *) ik, ie, ix, j, i, buf0(idx)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(11)

  ! ! 2. Redistribute to g(1:NM, 1:NM, 1:NK, 1:NE, 1:local_NX)
  ! allocate(buf1(NM * NM * NK * NE * local_NX))
  ! count = NM * NM * NK * local_NE * local_NX
  ! call MPI_Alltoall(buf0, count, MPI_DOUBLE_COMPLEX, buf1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)


  ! filename = 'P_alltoall_buf_r' // rank_str // '_'
  ! open(unit=11, file=trim(filename)//'.dat', status='unknown')
  ! ! do ik=1,nk
  ! !   do ie=first_local_ne,first_local_ne + local_ne - 1
  ! !     do ix=1,nx
  ! !       do j=1,nm
  ! !         do i=1,nm
  ! !           write(11, *) ik, ie, ix, j, i, buf(i,j,ik,ie-first_local_ne,ix)
  ! !         enddo
  ! !       enddo
  ! !     enddo
  ! !   enddo
  ! ! enddo
  ! do ix=first_local_NX,first_local_NX + local_NX - 1
  !   do ie=1,ne
  !     do ik=1,nk
  !       do j=1,nm
  !         do i=1,nm
  !           idx = (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-1) * nm * nm * nk + (ix-first_local_NX) * nm * nm * nk * ne
  !           write(11, *) ik, ie, ix, j, i, buf1(idx)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo
  ! close(11)

  ! deallocate(g)
  ! allocate(g(NM, NM, NK, NE, local_NX))
  ! do ix=first_local_NX,first_local_NX + local_NX - 1
  !   do ie=1,ne
  !     do ik=1,nk
  !       do j=1,nm
  !         do i=1,nm
  !           idx = (i-1) + (j-1) * nm + (ik-1) * nm * nm + (ie-1) * nm * nm * nk + (ix-first_local_NX) * nm * nm * nk * ne
  !           g(i,j,ik,ie,ix) = buf1(idx)
  !         enddo
  !       enddo
  !     enddo
  !   enddo
  ! enddo

  ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! print *, 'rank = ', comm_rank, ', ierr = ', ierr

  ! deallocate(buf0, buf1)

endsubroutine g2p


! Redistribute p(1:NM, 1:NM, 1:NK, 1:NE, 1:local_NX) to p(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK)
subroutine p2g(p, g, NM, NX, NE, NK, local_NX, local_NE)
  complex ( 8 ), allocatable, intent ( inout ) :: p(:, :, :, :, :), g(:, :, :, :, :)
  integer ( kind = 4 ), intent ( in ) :: NM, NX, NE, NK, local_NX, local_NE
  complex ( 8 ), allocatable :: buf0(:, :, :, :, :), buf1(:, :, :, :, :)
  integer ( kind = 4 ) :: count, ierr

  ! Assume p(1:NM, 1:NM, 1:NK, 1:NE, 1:local_NX)

  ! 1. Reshape to p(1:NM, 1:NM, 1:NK, 1:local_NX, 1:NE)
  allocate(buf0(NM, NM, NK, local_NX, NE))
  buf0 = reshape(p, shape(buf0), order = [1, 2, 3, 5, 4])

  ! 2. Redistribute to p(1:NM, 1:NM, 1:NK, 1:NX, 1:local_NE)
  allocate(buf1(NM, NM, NK, NX, local_NE))
  count = NM * NM * NK * local_NX * local_NE
  call MPI_Alltoall(buf0, count, MPI_DOUBLE_COMPLEX, buf1, count, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierr)

  ! 3. Reshape to p(1:NM, 1:NM, 1:NX, 1:local_NE, 1:NK)
  ! deallocate(p)
  ! allocate(p(NM, NM, NX, local_NE, NK))
  g = reshape(buf1, shape(g), order = [1, 2, 5, 3, 4])

  deallocate(buf0, buf1)

endsubroutine p2g

! calculate e-phonon self-energies in the deformation potential approximation
subroutine calc_sigma_ephonon_DPA(nm,nx,nen,nmode,Dac,Dop,hw,en,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater)
integer,intent(in)::nm,nx,nen
integer,intent(in)::nmode ! number of optical phonon modes (Einstein model)
real(8),intent(in)::Dac
real(8),intent(in)::Dop(nmode)
real(8),intent(in)::hw(nmode)
real(8),intent(in)::en(nen)
complex(8),intent(in),dimension(nm,nm,nx,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm,nm,nx,nen)::Sig_retarded,Sig_lesser,Sig_greater
!---------
integer::ie,ix
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
real(8)::nbose

end subroutine calc_sigma_ephonon_DPA

! calculate e-photon self-energies in the monochromatic assumption
subroutine calc_sigma_ephoton_monochromatic(nm,length,nen,En,nop,Mii,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,pre_fact,intensity,hw)
integer,intent(in)::nm,length,nen,nop
real(8),intent(in)::en(nen),pre_fact,intensity,hw
complex(8),intent(in),dimension(nm,nm,length)::Mii ! e-photon interaction matrix blocks
complex(8),intent(in),dimension(nm,nm,length,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm,nm,length,nen)::Sig_retarded,Sig_lesser,Sig_greater
!---------
integer::ie,ix
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
  Sig_lesser=czero
  Sig_greater=czero
  Sig_retarded=czero
  ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M
  !           ~ M [ G^<>(E -+ hw) + G^<>(E +- hw)] M * N
  !$omp parallel default(none) private(ie,A,B,ix) shared(length,nop,nen,nm,G_lesser,G_greater,Sig_lesser,Sig_greater,Mii)
  allocate(B(nm,nm))
  allocate(A(nm,nm))  
  !$omp do
  do ie=1,nen
    do ix=1,length
      ! Sig^<(E)
      A = czero
      if (ie-nop>=1) A =A+ G_lesser(:,:,ix,ie-nop)
      if (ie+nop<=nen) A =A+ G_lesser(:,:,ix,ie+nop)
      call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
      call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
      Sig_lesser(:,:,ix,ie) = A 
      ! Sig^>(E)
      A = czero
      if (ie-nop>=1) A =A+ G_greater(:,:,ix,ie-nop)
      if (ie+nop<=nen) A =A+ G_greater(:,:,ix,ie+nop)
      call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
      call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
      Sig_greater(:,:,ix,ie) = A
    enddo
  enddo  
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Sig_greater=Sig_greater*pre_fact*intensity/hw**2
  Sig_lesser=Sig_lesser*pre_fact*intensity/hw**2
  Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
end subroutine calc_sigma_ephoton_monochromatic


! calculate e-photon self-energies in the monochromatic assumption
subroutine calc_sigma_ephoton_monochromatic_ext(nm,length,nen,En,nop,Mii,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,pre_fact,intensity,hw)
  integer,intent(in)::nm,length,nen,nop
  real(8),intent(in)::en(nen),pre_fact,intensity,hw
  complex(8),intent(in),dimension(nm,nm,length)::Mii ! e-photon interaction matrix blocks
  complex(8),intent(in),dimension(nm,nm,length,nen+2*nop)::G_lesser,G_greater
  complex(8),intent(inout),dimension(nm,nm,length,nen)::Sig_retarded,Sig_lesser,Sig_greater
  !---------
  integer::ie,ix,ie2
  complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
    Sig_lesser=czero
    Sig_greater=czero
    Sig_retarded=czero
    ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M
    !           ~ M [ G^<>(E -+ hw) + G^<>(E +- hw)] M * N
    !$omp parallel default(none) private(ie,A,B,ix,ie2) shared(length,nop,nen,nm,G_lesser,G_greater,Sig_lesser,Sig_greater,Mii)
    allocate(B(nm,nm))
    allocate(A(nm,nm))  
    !$omp do
    do ie=1,nen
      ie2 = ie+nop
      do ix=1,length
        ! Sig^<(E)
        A = czero
        ! if (ie-nop>=1) A =A+ G_lesser(:,:,ix,ie2-nop)
        ! if (ie+nop<=nen) A =A+ G_lesser(:,:,ix,ie2+nop)
        A =A+ G_lesser(:,:,ix,ie2-nop)
        A =A+ G_lesser(:,:,ix,ie2+nop)
        call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
        Sig_lesser(:,:,ix,ie) = A 
        ! Sig^>(E)
        A = czero
        ! if (ie-nop>=1) A =A+ G_greater(:,:,ix,ie2-nop)
        ! if (ie+nop<=nen) A =A+ G_greater(:,:,ix,ie2+nop)
        A =A+ G_greater(:,:,ix,ie2-nop)
        A =A+ G_greater(:,:,ix,ie2+nop)
        call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
        Sig_greater(:,:,ix,ie) = A
      enddo
    enddo  
    !$omp end do
    deallocate(A,B)
    !$omp end parallel
    Sig_greater=Sig_greater*pre_fact*intensity/hw**2
    Sig_lesser=Sig_lesser*pre_fact*intensity/hw**2
    Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
  end subroutine calc_sigma_ephoton_monochromatic_ext


subroutine calc_sigma_ephoton_monochromatic_nx(nm,length,nen,En,nop,Mii,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,pre_fact,intensity,hw)
  integer,intent(in)::nm,length,nen,nop
  real(8),intent(in)::en(nen),pre_fact,intensity,hw
  complex(8),intent(in),dimension(nm,nm,length)::Mii ! e-photon interaction matrix blocks
  complex(8),intent(in),dimension(nm,nm,nen,length)::G_lesser,G_greater
  complex(8),intent(inout),dimension(nm,nm,nen,length)::Sig_retarded,Sig_lesser,Sig_greater
  !---------
  integer::ie,ix
  complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
    Sig_lesser=czero
    Sig_greater=czero
    Sig_retarded=czero
    ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M
    !           ~ M [ G^<>(E -+ hw) + G^<>(E +- hw)] M * N
    !$omp parallel default(none) private(ie,A,B,ix) shared(length,nop,nen,nm,G_lesser,G_greater,Sig_lesser,Sig_greater,Mii)
    allocate(B(nm,nm))
    allocate(A(nm,nm))  
    !$omp do
    do ix=1,length
      do ie=1,nen
        ! Sig^<(E)
        A = czero
        if (ie-nop>=1) A =A+ G_lesser(:,:,ie-nop,ix)
        if (ie+nop<=nen) A =A+ G_lesser(:,:,ie+nop,ix)
        call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
        Sig_lesser(:,:,ie,ix) = A 
        ! Sig^>(E)
        A = czero
        if (ie-nop>=1) A =A+ G_greater(:,:,ie-nop,ix)
        if (ie+nop<=nen) A =A+ G_greater(:,:,ie+nop,ix)
        call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
        call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
        Sig_greater(:,:,ie,ix) = A
      enddo
    enddo  
    !$omp end do
    deallocate(A,B)
    !$omp end parallel
    Sig_greater=Sig_greater*pre_fact*intensity/hw**2
    Sig_lesser=Sig_lesser*pre_fact*intensity/hw**2
    Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
  end subroutine calc_sigma_ephoton_monochromatic_nx



subroutine green_rgf_solve_gw_ephoton_3d(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg,Pii,P1i,polarization,intensity,hw,labs)
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
  real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser, g_r_i1
  real(8),allocatable,dimension(:,:,:,:,:)::cur,jdens
  real(8),allocatable,dimension(:,:,:,:)::tot_cur,tot_ecur
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser,P_greater,P_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_i1,W_greater_i1,W_retarded_i1
  complex(8),allocatable,dimension(:,:,:)::Sii
  complex(8),allocatable,dimension(:,:,:,:)::Mii,M1i
  complex(8),allocatable,dimension(:,:,:)::Pi_retarded_ik,Pi_lesser_ik,Pi_greater_ik,Pi_retarded
  complex(8),allocatable::Ispec(:,:,:,:),Itot(:,:,:),Ispec_ik(:,:,:,:),Itot_ik(:,:,:)
  real(8)::tr(nen,nky*nkz),tre(nen,nky*nkz)
  integer::ie,iter,i,ix,nopmax,nop,nopphot,iop,l,h,io
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk
  complex(8)::dE, B(nm,nm)
  !
  print *,'======== green_rgf_solve_gw_ephoton_3D ========'
  nk=nky*nkz
  allocate(g_r(nm,nm,nx,nen,nk))
  allocate(g_r_i1(nm,nm,nx,nen,nk))
  allocate(g_greater(nm,nm,nx,nen,nk))
  allocate(g_lesser(nm,nm,nx,nen,nk))
  allocate(cur(nm,nm,nx,nen,nk))
  allocate(jdens(nb,nb,nx*ns,nen,nk))
  allocate(tot_cur(nb,nb,nx*ns,nk))
  allocate(tot_ecur(nb,nb,nx*ns,nk))
  allocate(sigma_lesser_gw(nm,nm,nx,nen,nk))
  allocate(sigma_greater_gw(nm,nm,nx,nen,nk))
  allocate(sigma_r_gw(nm,nm,nx,nen,nk))
  allocate(sigma_lesser_new(nm,nm,nx,nen,nk))
  allocate(sigma_greater_new(nm,nm,nx,nen,nk))
  allocate(sigma_r_new(nm,nm,nx,nen,nk))
  allocate(P_lesser(nm,nm,nx,nen,nk))
  allocate(P_greater(nm,nm,nx,nen,nk))
  allocate(P_retarded(nm,nm,nx,nen,nk))
!  allocate(P_lesser_1i(nm,nm,nx,nen,nk))
!  allocate(P_greater_1i(nm,nm,nx,nen,nk))
!  allocate(P_retarded_1i(nm,nm,nx,nen,nk))
  allocate(W_lesser(nm,nm,nx,nen,nk))
  allocate(W_greater(nm,nm,nx,nen,nk))
  allocate(W_retarded(nm,nm,nx,nen,nk))
!  allocate(W_lesser_i1(nm,nm,nx,nen,nk))
!  allocate(W_greater_i1(nm,nm,nx,nen,nk))
!  allocate(W_retarded_i1(nm,nm,nx,nen,nk))  
  allocate(Mii(nm,nm,nx,nk))
  allocate(Ispec(nm,nm,nx,nen))
  allocate(Itot(nm,nm,nx))
  allocate(Ispec_ik(nm,nm,nx,nen))
  allocate(Itot_ik(nm,nm,nx))
  Mii(:,:,:,:)=czero
  do i=1,3
    Mii(:,:,:,:)=Mii(:,:,:,:)+ polarization(i) * Pii(:,:,i,:,:) 
  enddo   
  print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity  
  nopphot=floor(hw / (En(2)-En(1)))
  print *,'nop photon=',nopphot
  do iter=0,niter
    print *,'+ iter=',iter
    !
    print *, 'calc G'          
    do ik=1,nk
      print *, ' ik=', ik,'/',nk
      !$omp parallel default(none) private(ie) shared(ik,nen,temp,nm,nx,en,mu,Hii,H1i,sigma_lesser_gw,sigma_greater_gw,sigma_r_gw,g_lesser,g_greater,g_r,tre,tr,cur,g_r_i1)    
      !$omp do 
      do ie=1,nen     
        !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
        call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii(:,:,:,ik),H1i(:,:,:,ik),sigma_lesser_gw(:,:,:,ie,ik),sigma_greater_gw(:,:,:,ie,ik),&
                        & sigma_r_gw(:,:,:,ie,ik),g_lesser(:,:,:,ie,ik),g_greater(:,:,:,ie,ik),g_r(:,:,:,ie,ik),tr(ie,ik),&
                        & tre(ie,ik),cur(:,:,:,ie,ik),g_r_i1(:,:,:,ie,ik)) 
      enddo
      !$omp end do 
      !$omp end parallel
      call calc_block_current(Hii(:,:,:,ik),G_lesser(:,:,:,:,ik),cur(:,:,:,:,ik),nen,en,spindeg,nb,ns,nm,nx,tot_cur(:,:,:,ik),tot_ecur(:,:,:,ik),jdens(:,:,:,:,ik))      
    enddo
    !
    call write_spectrum_summed_over_k('gw_ldos',iter,g_r,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
    call write_spectrum_summed_over_k('gw_ndos',iter,g_lesser,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
    call write_spectrum_summed_over_k('gw_pdos',iter,g_greater,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/))     
    call write_dos_summed_over_k('gw_dos',iter,G_r,nen,En,nk,nx,NB,NS,Lx)
    call write_current_spectrum_summed_over_kz('gw_Jdens',iter,jdens,nen,En,nx*ns,NB,Lx,nk) 
    call write_current_spectrum_summed_over_kz('gw_Jx',iter,cur,nen,En,nx,NB*ns,Lx*ns,nk)  
    call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
    call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)
    call write_transmission_spectrum_k('gw_trR',iter,tr,nen,En,nk)
    call write_transmission_spectrum_k('gw_trL',iter,tre,nen,En,nk)
    open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
    close(101)
    g_r = dcmplx( 0.0d0*dble(g_r), aimag(g_r))
    g_lesser = dcmplx( 0.0d0*dble(g_lesser), aimag(g_lesser))
    g_greater = dcmplx( 0.0d0*dble(g_greater), aimag(g_greater))
    !        
    print *, 'calc P'    
    nopmax=nen/2-10      
    P_lesser = dcmplx(0.0d0,0.0d0)
    P_greater = dcmplx(0.0d0,0.0d0)    
    P_retarded = dcmplx(0.0d0,0.0d0)    
!    P_lesser_1i = dcmplx(0.0d0,0.0d0)
!    P_greater_1i = dcmplx(0.0d0,0.0d0)    
!    P_retarded_1i = dcmplx(0.0d0,0.0d0)  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do         
    do nop=-nopmax,nopmax
      iop=nop+nen/2  
      do iqy=1,nky        
        do iqz=1,nkz
          iq=iqz + (iqy-1)*nkz
          do ix=1,nx
            do i=1,nm      
              l=max(i-ndiag,1)
              h=min(nm,i+ndiag)   
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
                  do ie = max(nop+1,1),min(nen,nen+nop)           
                    P_lesser(i,l:h,ix,iop,iq) = P_lesser(i,l:h,ix,iop,iq) + G_lesser(i,l:h,ix,ie,ik) * G_greater(l:h,i,ix,ie-nop,ikd)
                    P_greater(i,l:h,ix,iop,iq) = P_greater(i,l:h,ix,iop,iq) + G_greater(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd)        
                    P_retarded(i,l:h,ix,iop,iq) = P_retarded(i,l:h,ix,iop,iq) + (G_lesser(i,l:h,ix,ie,ik) * conjg(G_r(i,l:h,ix,ie-nop,ikd)) & 
                                              & + G_r(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd))        
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo      
    enddo          
    !$omp end do 
    !$omp end parallel
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nk)
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=P_retarded*dE
    call write_spectrum_summed_over_k('gw_PR',iter,P_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PG',iter,P_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc W'                
    do iq=1,nky*nkz
      print *, ' iq=', iq,'/',nk
      !$omp parallel default(none) private(nop) shared(iq,nopmax,nen,nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,w_lesser,w_greater,w_retarded)    
      !$omp do 
      do nop=-nopmax+nen/2,nopmax+nen/2 
        !call green_calc_w_full(0,nm,nx,Vii(:,:,:,iq),V1i(:,:,:,iq),p_lesser(:,:,:,nop,iq),p_greater(:,:,:,nop,iq),p_retarded(:,:,:,nop,iq),w_lesser(:,:,:,nop,iq),w_greater(:,:,:,nop,iq),w_retarded(:,:,:,nop,iq))
        
        call green_rgf_calc_w(0,nm,nx,Vii(:,:,:,iq),V1i(:,:,:,iq),p_lesser(:,:,:,nop,iq),p_greater(:,:,:,nop,iq),p_retarded(:,:,:,nop,iq),w_lesser(:,:,:,nop,iq),w_greater(:,:,:,nop,iq),w_retarded(:,:,:,nop,iq))
        
      enddo
      !$omp end do 
      !$omp end parallel
    enddo
    call write_spectrum_summed_over_k('gw_WR',iter,W_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WG',iter,W_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc SigGW'
    nopmax=nen/2-10
    Sigma_greater_new = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new = dcmplx(0.0d0,0.0d0)
    Sigma_r_new = dcmplx(0.0d0,0.0d0)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,w_lesser,w_greater,w_retarded,sigma_lesser_new,sigma_greater_new,sigma_r_new,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do    
    do ie=1,nen      
      do iky=1,nky
        do ikz=1,nkz
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then   
              ik=ikz+(iky-1)*nkz
              do iqy=1,nky
                do iqz=1,nkz              
                  iq=iqz + (iqy-1)*nkz
                  iop=nop+nen/2                
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2            
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky        
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1        
                  ikd=ikzd + (ikyd-1)*nkz
                  do ix=1,nx
                    do i=1,nm
                      l=max(i-ndiag,1)
                      h=min(nm,i+ndiag)                                                
                      Sigma_lesser_new(i,l:h,ix,ie,ik)=Sigma_lesser_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq)
                      Sigma_greater_new(i,l:h,ix,ie,ik)=Sigma_greater_new(i,l:h,ix,ie,ik)+G_greater(i,l:h,ix,ie-nop,ikd)*W_greater(i,l:h,ix,iop,iq)
                      Sigma_r_new(i,l:h,ix,ie,ik)=Sigma_r_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq)    
                     enddo
                  enddo
                enddo
              enddo
            endif            
          enddo
        enddo
      enddo      
    enddo
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nk) 
    Sigma_lesser_new = Sigma_lesser_new  * dE
    Sigma_greater_new= Sigma_greater_new * dE
    Sigma_r_new=Sigma_r_new* dE
    Sigma_r_new = dcmplx( dble(Sigma_r_new), aimag(Sigma_greater_new-Sigma_lesser_new)/2.0d0 )    
    ! symmetrize the selfenergies
    do ie=1,nen
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new(:,:,ix,ie,ik))
          Sigma_r_new(:,:,ix,ie,ik) = (Sigma_r_new(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new(:,:,ix,ie,ik))
          Sigma_lesser_new(:,:,ix,ie,ik) = (Sigma_lesser_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new(:,:,ix,ie,ik))
          Sigma_greater_new(:,:,ix,ie,ik) = (Sigma_greater_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo
    ! mixing with the previous one
    Sigma_r_gw= Sigma_r_gw+ alpha_mix * (Sigma_r_new -Sigma_r_gw)
    Sigma_lesser_gw  = Sigma_lesser_gw+ alpha_mix * (Sigma_lesser_new -Sigma_lesser_gw)
    Sigma_greater_gw = Sigma_greater_gw+ alpha_mix * (Sigma_greater_new -Sigma_greater_gw)  
    ! 
    call write_spectrum_summed_over_k('gw_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !!!!! calculate collision integral
    Ispec=czero
    Itot=czero
    do ik=1,nk
      call calc_block_collision(sigma_lesser_gw(:,:,:,:,ik),sigma_greater_gw(:,:,:,:,ik),G_lesser(:,:,:,:,ik),G_greater(:,:,:,:,ik),nen,en,spindeg,nm,nx,Itot_ik,Ispec_ik)
      Ispec=Ispec+Ispec_ik
      Itot=Itot+Itot_ik
    enddo
    call write_spectrum('gw_Scat',iter,Ispec,nen,En,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/))    
    !
    if (iter>=(niter-15)) then
      print *, 'calc SigEPhoton'
      do ik=1,nk
        print *, ' ik=', ik,'/',nk
        call calc_sigma_ephoton_monochromatic(nm,nx,nen,En,nopphot,Mii(:,:,:,ik),G_lesser(:,:,:,:,ik),G_greater(:,:,:,:,ik),Sigma_r_new(:,:,:,:,ik),Sigma_lesser_new(:,:,:,:,ik),Sigma_greater_new(:,:,:,:,ik),pre_fact,intensity,hw)
      enddo
      Sigma_r_gw       = Sigma_r_gw + Sigma_r_new 
      Sigma_lesser_gw  = Sigma_lesser_gw + Sigma_lesser_new 
      Sigma_greater_gw = Sigma_greater_gw + Sigma_greater_new  
      call write_spectrum_summed_over_k('eph_SigR',iter,Sigma_r_new,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('eph_SigL',iter,Sigma_lesser_new,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('eph_SigG',iter,Sigma_greater_new,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    endif
    ! make sure self-energy is continuous near leads (by copying edge block)
    do ix=1,2
      Sigma_r_gw(:,:,ix,:,:)=Sigma_r_gw(:,:,3,:,:)
      Sigma_lesser_gw(:,:,ix,:,:)=Sigma_lesser_gw(:,:,3,:,:)
      Sigma_greater_gw(:,:,ix,:,:)=Sigma_greater_gw(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw(:,:,nx-ix+1,:,:)=Sigma_r_gw(:,:,nx-2,:,:)
      Sigma_lesser_gw(:,:,nx-ix+1,:,:)=Sigma_lesser_gw(:,:,nx-2,:,:)
      Sigma_greater_gw(:,:,nx-ix+1,:,:)=Sigma_greater_gw(:,:,nx-2,:,:)
    enddo    
  enddo
  !!! calculate G for the last time
  print *, 'calc G'          
  do ik=1,nk
    print *, ' ik=', ik,'/',nk
    !$omp parallel default(none) private(ie) shared(ik,nen,temp,nm,nx,en,mu,Hii,H1i,sigma_lesser_gw,sigma_greater_gw,sigma_r_gw,g_lesser,g_greater,g_r,tre,tr,cur,g_r_i1)    
    !$omp do 
    do ie=1,nen     
      !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
      call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii(:,:,:,ik),H1i(:,:,:,ik),sigma_lesser_gw(:,:,:,ie,ik),sigma_greater_gw(:,:,:,ie,ik),&
                      & sigma_r_gw(:,:,:,ie,ik),g_lesser(:,:,:,ie,ik),g_greater(:,:,:,ie,ik),g_r(:,:,:,ie,ik),tr(ie,ik),&
                      & tre(ie,ik),cur(:,:,:,ie,ik),g_r_i1(:,:,:,ie,ik)) 
    enddo
    !$omp end do 
    !$omp end parallel
    call calc_block_current(Hii(:,:,:,ik),G_lesser(:,:,:,:,ik),cur(:,:,:,:,ik),nen,en,spindeg,nb,ns,nm,nx,tot_cur(:,:,:,ik),tot_ecur(:,:,:,ik),jdens(:,:,:,:,ik))      
  enddo
  !
  call write_spectrum_summed_over_k('gw_ldos',iter,g_r,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
  call write_spectrum_summed_over_k('gw_ndos',iter,g_lesser,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
  call write_spectrum_summed_over_k('gw_pdos',iter,g_greater,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/))     
  call write_dos_summed_over_k('gw_dos',iter,G_r,nen,En,nk,nx,NB,NS,Lx)
  call write_current_spectrum_summed_over_kz('gw_Jdens',iter,jdens,nen,En,nx*ns,NB,Lx,nk)  
  call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
  call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)
  call write_transmission_spectrum_k('gw_trR',iter,tr*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
  call write_transmission_spectrum_k('gw_trL',iter,tre*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
  open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
  write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
  close(101)
  !!!!!!!!!
  deallocate(g_r,g_lesser,g_greater,cur)
  deallocate(sigma_lesser_gw,sigma_greater_gw,sigma_r_gw)   
  deallocate(sigma_lesser_new,sigma_greater_new,sigma_r_new)   
  deallocate(P_retarded,P_lesser,P_greater)
  !deallocate(P_retarded_1i,P_lesser_1i,P_greater_1i)
  deallocate(W_retarded,W_lesser,W_greater)
  !deallocate(W_retarded_i1,W_lesser_i1,W_greater_i1)  
  deallocate(Mii)
  deallocate(Ispec,Itot)
  deallocate(Ispec_ik,Itot_ik)
end subroutine green_rgf_solve_gw_ephoton_3d



subroutine green_rgf_solve_gw_ephoton_3d_ijs(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg,Pii,P1i,polarization,intensity,hw,labs,lnogw, comm_size, comm_rank, local_NE, first_local_energy)
  use fft_mod, only : conv1d => conv1d2, corr1d => corr1d2
  integer ( kind = 4), intent(in) :: comm_size, comm_rank, local_NE, first_local_energy
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  !complex(8), intent(in):: V(nm*nx,nm*nx,nky*nkz)
  real(8), intent(in) :: polarization(3) ! light polarization vector 
  real(8), intent(in) :: intensity ! [W/m^2]
  logical, intent(in) :: labs ! whether to calculate Pi and absorption
  logical, intent(in) :: lnogw ! whether to activate gw
  complex(8), intent(in):: Pii(nm,nm,3,nx,nky*nkz),P1i(nm,nm,3,nx,nky*nkz) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
  real(8), intent(in) :: hw ! hw is photon energy in eV
  ! -------- local variables
  real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
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

  ! For debugging/validation
  ! complex(8), pointer :: g_r_buf2(:), g_r_by_energies2(:, :, :, :, :), g_r_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: g_greater_buf2(:), g_greater_by_energies2(:, :, :, :, :), g_greater_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: g_lesser_buf2(:), g_lesser_by_energies2(:, :, :, :, :), g_lesser_by_blocks2(:, :, :, :, :)

  ! complex(8), pointer :: P_retarded_buf2(:), P_retarded_by_energies2(:, :, :, :, :), P_retarded_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: P_greater_buf2(:), P_greater_by_energies2(:, :, :, :, :), P_greater_by_blocks2(:, :, :, :, :)
  ! complex(8), pointer :: P_lesser_buf2(:), P_lesser_by_energies2(:, :, :, :, :), P_lesser_by_blocks2(:, :, :, :, :)  

  local_Nij = (nm * nm) / comm_size
  first_local_ij = local_Nij * comm_rank + 1

  local_NX = nx / comm_size
  first_local_block = local_NX * comm_rank + 1

  ! num_g = nm * nm * local_NE * nkz * nky * local_NX
  ! num_p = nm * nm * local_NE * nkz * nky * local_NX

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

  ! real(8) :: extended_local_energies(local_NE + 2 * nopphot)
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
    ! do ie = nopphot + local_NE + 1, local_NE + 2 * nopphot
    do ie = 1, nopphot
      extended_local_energies(ie + nopphot + local_NE) = en(local_NE + first_local_energy + ie - 1)
    enddo
  endif

  allocate(g_r_buf(nm * nm * local_NE * nk * nx))
  allocate(g_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(g_lesser_buf(nm * nm * nx * local_NE * nk))
  ! allocate(g_greater_extended_buf(nm * nm * nx * (local_NE + 2 * nopphot) * nk))
  ! allocate(g_lesser_extended_buf(nm * nm * nx * (local_NE + 2 * nopphot) * nk))
  ! count = nm * nm * nx * nopphot * nk
  ! g_greater_buf(1:nm * nm * local_NE * nk * nx) => g_greater_extended_buf(count:count + nm * nm * local_NE * nk * nx - 1)
  ! g_lesser_buf(1:nm * nm * local_NE * nk * nx) => g_lesser_extended_buf(count:count + nm * nm * local_NE * nk * nx - 1)

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

  ! g_lesser_extended_buf = dcmplx(0.0d0,0.0d0)
  ! g_greater_extended_buf = dcmplx(0.0d0,0.0d0)
  ! g_lesser_extended(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot, 1:nk) => g_lesser_extended_buf
  ! g_greater_extended(1:nm, 1:nm, 1:nx, 1:local_NE + 2 * nopphot, 1:nk) => g_greater_extended_buf

  g_r_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_r_buf
  g_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_greater_buf
  g_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_lesser_buf

  ! g_r_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_r_buf
  ! g_greater_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_greater_buf
  ! g_lesser_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_lesser_buf
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

  ! sigma_r_new_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_r_new_buf
  ! sigma_greater_new_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_greater_new_buf
  ! sigma_lesser_new_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => sigma_lesser_new_buf
  sigma_r_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_r_new_buf
  sigma_greater_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_greater_new_buf
  sigma_lesser_new_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => sigma_lesser_new_buf
if(.not. lnogw) then
  allocate(P_retarded_buf(nm * nm * nx * local_NE * nk))
  allocate(P_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(P_lesser_buf(nm * nm * nx * local_NE * nk))

  P_retarded_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_retarded_buf
  P_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_greater_buf
  P_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => P_lesser_buf

  ! P_retarded_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => P_retarded_buf
  ! P_greater_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => P_greater_buf
  ! P_lesser_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => P_lesser_buf
  P_retarded_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_retarded_buf
  P_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_greater_buf
  P_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => P_lesser_buf

  allocate(W_retarded_buf(nm * nm * nx * local_NE * nk))  
  allocate(W_greater_buf(nm * nm * nx * local_NE * nk))
  allocate(W_lesser_buf(nm * nm * nx * local_NE * nk))

  W_retarded_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_retarded_buf
  W_greater_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_greater_buf
  W_lesser_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => W_lesser_buf

  ! W_retarded_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => W_retarded_buf
  ! W_greater_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => W_greater_buf
  ! W_lesser_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => W_lesser_buf
  W_retarded_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_retarded_buf
  W_greater_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_greater_buf
  W_lesser_by_blocks(1:NE, 1:nk, 1:NX, 1:local_Nij) => W_lesser_buf
endif

  allocate(cur(nm, nm, nx, local_NE, nk))

  allocate(tmp0(nm * nm * nx * local_NE * nk))
  allocate(tmp1(nm * nm * nx * local_NE * nk))

  allocate(g_greater_t_buf(nm * nm * nx * local_NE * nk))
  allocate(g_lesser_t_buf(nm * nm * nx * local_NE * nk))

  g_greater_t_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_greater_t_buf
  g_lesser_t_by_energies(1:nm, 1:nm, 1:nx, 1:local_NE, 1:nk) => g_lesser_t_buf

  ! g_greater_t_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_greater_t_buf
  ! g_lesser_t_by_blocks(1:local_Nij, 1:NX, 1:NE, 1:nk) => g_lesser_t_buf
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
  ! nopphot=floor(hw / (En(2)-En(1)))
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
      write(101,'(I4,2E16.6)') iter, global_sum_tr*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), global_sum_tre*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
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
if (.not. lnogw) then
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
                
                
!                do nop=-nopmax,nopmax
!                  iop=nop+nen/2  
!                  do ie = max(nop+1,1),min(nen,nen+nop)
!                    ! P_lesser_by_blocks(i,l:h,ix,iop,iq) = P_lesser_by_blocks(i,l:h,ix,iop,iq) + G_lesser_by_blocks(i,l:h,ix,ie,ik) * G_greater_by_blocks(l:h,i,ix,ie-nop,ikd)
!                    ! P_greater_by_blocks(i,l:h,ix,iop,iq) = P_greater_by_blocks(i,l:h,ix,iop,iq) + G_greater_by_blocks(i,l:h,ix,ie,ik) * G_lesser_by_blocks(l:h,i,ix,ie-nop,ikd)        
!                    ! P_retarded_by_blocks(i,l:h,ix,iop,iq) = P_retarded_by_blocks(i,l:h,ix,iop,iq) + (G_lesser_by_blocks(i,l:h,ix,ie,ik) * conjg(G_r_by_blocks(i,l:h,ix,ie-nop,ikd)) & 
!                    !                           & + G_r_by_blocks(i,l:h,ix,ie,ik) * G_lesser_by_blocks(l:h,i,ix,ie-nop,ikd))
                    
!                    ! P_lesser_by_blocks(i,ix,iop,iq) = P_lesser_by_blocks(i,ix,iop,iq) + G_lesser_by_blocks(i,ix,ie,ik) * G_greater_t_by_blocks(i,ix,ie-nop,ikd)
!                    ! P_greater_by_blocks(i,ix,iop,iq) = P_greater_by_blocks(i,ix,iop,iq) + G_greater_by_blocks(i,ix,ie,ik) * G_lesser_t_by_blocks(i,ix,ie-nop,ikd)        
!                    ! P_retarded_by_blocks(i,ix,iop,iq) = P_retarded_by_blocks(i,ix,iop,iq) + (G_lesser_by_blocks(i,ix,ie,ik) * conjg(G_r_by_blocks(i, ix,ie-nop,ikd)) & 
!                    !                           & + G_r_by_blocks(i,ix,ie,ik) * G_lesser_t_by_blocks(i,ix,ie-nop,ikd))

!                    P_lesser_by_blocks(iop, iq, ix, i) = P_lesser_by_blocks(iop, iq, ix, i) + G_lesser_by_blocks(ie, ik, ix, i) * G_greater_t_by_blocks(ie-nop, ikd, ix, i)
!                    P_greater_by_blocks(iop, iq, ix, i) = P_greater_by_blocks(iop, iq, ix, i) + G_greater_by_blocks(ie, ik, ix, i) * G_lesser_t_by_blocks(ie-nop, ikd, ix, i)        
!                    P_retarded_by_blocks(iop, iq, ix, i) = P_retarded_by_blocks(iop, iq, ix, i) + (G_lesser_by_blocks(ie, ik, ix, i) * conjg(G_r_by_blocks(ie-nop, ikd, ix, i)) & 
!                                              & + G_r_by_blocks(ie, ik, ix, i) * G_lesser_t_by_blocks(ie-nop, ikd, ix, i))  
                    
!                  enddo
!                enddo
                
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
!        call green_calc_w_full( &
!          0, nm, nx, Vii(:,:,:,iq), V1i(:,:,:,iq), &
!          p_lesser_by_energies(:,:,:,ie,iq), p_greater_by_energies(:,:,:,ie,iq), p_retarded_by_energies(:,:,:,ie,iq), &
!          w_lesser_by_energies(:,:,:,ie,iq), w_greater_by_energies(:,:,:,ie,iq), w_retarded_by_energies(:,:,:,ie,iq))

        call green_rgf_calc_w( &
          0, nm, nx, Vii(:,:,:,iq), V1i(:,:,:,iq), &
          p_lesser_by_energies(:,:,:,ie,iq), p_greater_by_energies(:,:,:,ie,iq), p_retarded_by_energies(:,:,:,ie,iq), &
          w_lesser_by_energies(:,:,:,ie,iq), w_greater_by_energies(:,:,:,ie,iq), w_retarded_by_energies(:,:,:,ie,iq))  
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
                
                sigma_lesser_new_by_blocks(:, ik, ix, i) = sigma_lesser_new_by_blocks(:, ik, ix, i) + &
                                                    conv1d(nen,G_lesser_by_blocks(:, ikd, ix, i),W_lesser_by_blocks(:, iq, ix, i),method=cmethod)
                
                Sigma_greater_new_by_blocks(:, ik, ix, i) = Sigma_greater_new_by_blocks(:, ik, ix, i) + &
                                                    conv1d(nen,G_greater_by_blocks(:, ikd, ix, i),W_greater_by_blocks(:, iq, ix, i),method=cmethod)
                
                Sigma_r_new_by_blocks(:, ik, ix, i) = Sigma_r_new_by_blocks(:, ik, ix, i) + & 
                                                    conv1d(nen,G_lesser_by_blocks(:, ikd, ix, i), W_retarded_by_blocks(:, iq, ix, i),method=cmethod) + &
                                                    conv1d(nen,G_r_by_blocks(:, ikd, ix, i), W_lesser_by_blocks(:, iq, ix, i),method=cmethod) + &
                                                    conv1d(nen,G_r_by_blocks(:, ikd, ix, i), W_retarded_by_blocks(:, iq, ix, i) ,method=cmethod)

                !Sigma_r_new_by_blocks(:, ik, ix, i) = Sigma_r_new_by_blocks(:, ik, ix, i) + sum(G_r_by_blocks(:, ikd, ix, i))* Vii(row,col,ix,iq)
                
                
!                do nop= -nopmax,nopmax
!                  iop=nop+nen/2
!                  do ie=1,nen     
!                    if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then
!                    ! do i=1,nm
!                    !   l=max(i-ndiag,1)
!                    !   h=min(nm,i+ndiag)                                                
!                      ! sigma_lesser_new_by_blocks(i,ix,ie,ik)=Sigma_lesser_new_by_blocks(i,ix,ie,ik)+G_lesser_by_blocks(i,ix,ie-nop,ikd)*W_lesser_by_blocks(i,ix,iop,iq)
!                      ! Sigma_greater_new_by_blocks(i,ix,ie,ik)=Sigma_greater_new_by_blocks(i,ix,ie,ik)+G_greater_by_blocks(i,ix,ie-nop,ikd)*W_greater_by_blocks(i,ix,iop,iq)
!                      ! Sigma_r_new_by_blocks(i,ix,ie,ik)=Sigma_r_new_by_blocks(i,ix,ie,ik)+G_lesser_by_blocks(i,ix,ie-nop,ikd)*W_retarded_by_blocks(i,ix,iop,iq) &
!                      !                                     &  +G_r_by_blocks(i,ix,ie-nop,ikd)*W_lesser_by_blocks(i,ix,iop,iq) &
!                      !                                     &  +G_r_by_blocks(i,ix,ie-nop,ikd)*W_retarded_by_blocks(i,ix,iop,iq)  
                      
!                      sigma_lesser_new_by_blocks(ie, ik, ix, i)=Sigma_lesser_new_by_blocks(ie, ik, ix, i)+G_lesser_by_blocks(ie-nop, ikd, ix, i)*W_lesser_by_blocks(iop, iq, ix, i)
                      
!                      Sigma_greater_new_by_blocks(ie, ik, ix, i)=Sigma_greater_new_by_blocks(ie, ik, ix, i)+G_greater_by_blocks(ie-nop, ikd, ix, i)*W_greater_by_blocks(iop, iq, ix, i)
                      
!                      Sigma_r_new_by_blocks(ie, ik, ix, i)=Sigma_r_new_by_blocks(ie, ik, ix, i)+G_lesser_by_blocks(ie-nop, ikd, ix, i)*W_retarded_by_blocks(iop, iq, ix, i) &
!                                                          &  +G_r_by_blocks(ie-nop, ikd, ix, i)*W_lesser_by_blocks(iop, iq, ix, i) &
!                                                          &  +G_r_by_blocks(ie-nop, ikd, ix, i)*W_retarded_by_blocks(iop, iq, ix, i)  
!                    endif
!                  enddo
!                enddo
                
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
  endif ! end computing GW
    if (iter>=(niter-5)) then

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
    do ix=1,2
      sigma_r_gw_by_energies(:,:,ix,:,:)=Sigma_r_gw_by_energies(:,:,3,:,:)
      sigma_lesser_gw_by_energies(:,:,ix,:,:)=Sigma_lesser_gw_by_energies(:,:,3,:,:)
      sigma_greater_gw_by_energies(:,:,ix,:,:)=Sigma_greater_gw_by_energies(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_r_gw_by_energies(:,:,nx-2,:,:)
      Sigma_lesser_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_lesser_gw_by_energies(:,:,nx-2,:,:)
      Sigma_greater_gw_by_energies(:,:,nx-ix+1,:,:)=Sigma_greater_gw_by_energies(:,:,nx-2,:,:)
    enddo

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
    write(101,'(I4,2E16.6)') iter, global_sum_tr*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), global_sum_tre*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
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
  if (.not. lnogw) then
    deallocate(P_retarded_buf, P_lesser_buf, P_greater_buf)
    deallocate(W_retarded_buf, W_lesser_buf, W_greater_buf)
  endif
  deallocate(Mii)

  deallocate(extended_local_energies)
end subroutine green_rgf_solve_gw_ephoton_3d_ijs




subroutine green_rgf_solve_gw_ephoton_3d_mpi_kpar(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg,Pii,P1i,polarization,intensity,hw,labs)
use, intrinsic :: iso_fortran_env
include "mpif.h"
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  real(8), intent(in) :: polarization(3) ! light polarization vector 
  real(8), intent(in) :: intensity ! [W/m^2]
  logical, intent(in) :: labs ! whether to calculate Pi and absorption
  complex(8), intent(in):: Pii(nm,nm,3,nx,nky*nkz),P1i(nm,nm,3,nx,nky*nkz) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
  real(8), intent(in) :: hw ! hw is photon energy in eV
  ! -------- local variables
  real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
  real(8),allocatable,dimension(:,:,:,:,:)::cur,cur_local_k,jdens,jdens_local_k
  real(8),allocatable,dimension(:,:,:,:)::tot_cur,tot_ecur,tot_cur_local_k,tot_ecur_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser
  complex(8),allocatable,dimension(:,:,:,:,:)::P_r,P_greater,P_lesser
  complex(8),allocatable,dimension(:,:,:,:,:)::W_r,W_greater,W_lesser
  complex(8),allocatable,dimension(:,:,:,:,:)::Sigma_r_gw,Sigma_lesser_gw,Sigma_greater_gw
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r_local_k,g_greater_local_k,g_lesser_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r_local_x,g_greater_local_x,g_lesser_local_x
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw_local_k,sigma_greater_gw_local_k,sigma_r_gw_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw_local_x,sigma_greater_gw_local_x,sigma_r_gw_local_x
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new_local_k,sigma_greater_new_local_k,sigma_r_new_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new_local_x,sigma_greater_new_local_x,sigma_r_new_local_x
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_local_k,P_greater_local_k,P_retarded_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_local_x,P_greater_local_x,P_retarded_local_x  
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_local_k,W_greater_local_k,W_retarded_local_k
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_local_x,W_greater_local_x,W_retarded_local_x  
  complex(8),allocatable,dimension(:,:,:,:,:)::sbuf,rbuf
  complex(8),allocatable::Ispec(:,:,:,:),Itot(:,:,:),Ispec_ik(:,:,:,:),Itot_ik(:,:,:)
  complex(8),allocatable,dimension(:,:,:)::Pi_retarded_ik,Pi_lesser_ik,Pi_greater_ik,Pi_retarded
  complex(8),dimension(nm,nm,nx,nky*nkz)::Mii!,M1i
  real(8),allocatable,dimension(:,:)::tr,tre,tr_local_k,tre_local_k
  integer::ie,iter,i,ix,nopmax,nop,nopphot,iop,l,h,io
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk,ndata
  complex(8)::dE, B(nm,nm)
  ! --- mpi
  integer(kind=int32) :: rank,num_proc
  integer(kind=int32) :: ierror
  integer::scount,rcount
  integer::nnode,inode,nnx
  ! Get the individual process (rank)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, num_proc, ierror)  
  !
  call MPI_Barrier(MPI_COMM_WORLD,ierror)
  if (rank == 0) then
    print *,'======== green_rgf_solve_gw_ephoton_3D_MPI ========'
  endif
  !
  nnode=num_proc
  nk=nky*nkz/nnode  ! each node does nk k point
  nnx=nx/nnode ! each node does nnx x points
  !  
  !
  allocate(g_r_local_k(nm,nm,nx,nen,nk))  
  allocate(g_greater_local_k(nm,nm,nx,nen,nk))
  allocate(g_lesser_local_k(nm,nm,nx,nen,nk))        
  !
  allocate(sigma_lesser_gw_local_k(nm,nm,nx,nen,nk))
  allocate(sigma_greater_gw_local_k(nm,nm,nx,nen,nk))
  allocate(sigma_r_gw_local_k(nm,nm,nx,nen,nk)) 
  sigma_lesser_gw_local_k = czero
  sigma_greater_gw_local_k = czero
  sigma_r_gw_local_k = czero
  !
  Mii(:,:,:,:)=czero  ! e-photon coupling matrix
  do i=1,3
    Mii(:,:,:,:)=Mii(:,:,:,:)+ polarization(i) * Pii(:,:,i,:,:) 
  enddo   
  call MPI_Barrier(MPI_COMM_WORLD,ierror)
  if (rank == 0) then
    print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity  
  endif
  nopphot=floor(hw / (En(2)-En(1)))
  if (rank == 0) then
    print *,'nop photon=',nopphot
  endif
  do iter=0,niter
    if (rank == 0) then
      print *,'+ iter=',iter
      !
      print *, 'calc G'      
    endif
    !
    allocate(cur_local_k(nm,nm,nx,nen,nk))      
    allocate(jdens_local_k(nb,nb,nx*ns,nen,nk))
    allocate(tot_cur_local_k(nb,nb,nx*ns,nk))
    allocate(tot_ecur_local_k(nb,nb,nx*ns,nk))
    allocate(tr_local_k(nen,nk))  
    allocate(tre_local_k(nen,nk))  
    !
    do ik=1,nk
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      !print *, ' ik=', ik+(rank)*nk,'/',nky*nkz      
      !$omp parallel default(none) private(ie) shared(ik,nen,temp,nm,nk,nx,rank,en,mu,Hii,H1i,sigma_lesser_gw_local_k,sigma_greater_gw_local_k,sigma_r_gw_local_k,g_lesser_local_k,g_greater_local_k,g_r_local_k,tre_local_k,tr_local_k,cur_local_k)    
      !$omp do 
      do ie=1,nen             
        call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii(:,:,:,ik+(rank)*nk),H1i(:,:,:,ik+(rank)*nk),&
                        & sigma_lesser_gw_local_k(:,:,:,ie,ik),sigma_greater_gw_local_k(:,:,:,ie,ik),&
                        & sigma_r_gw_local_k(:,:,:,ie,ik),g_lesser_local_k(:,:,:,ie,ik),g_greater_local_k(:,:,:,ie,ik),&
                        & g_r_local_k(:,:,:,ie,ik),tr_local_k(ie,ik),&
                        & tre_local_k(ie,ik),cur_local_k(:,:,:,ie,ik)) 
      enddo
      !$omp end do 
      !$omp end parallel  
      call calc_block_current(Hii(:,:,:,ik+(rank)*nk),G_lesser_local_k(:,:,:,:,ik),cur_local_k(:,:,:,:,ik),nen,en,spindeg,nb,ns,nm,nx,tot_cur_local_k(:,:,:,ik),tot_ecur_local_k(:,:,:,ik),jdens_local_k(:,:,:,:,ik))    
    enddo          
!    call write_spectrum_summed_over_k('gw_ldos_r'//string(rank),iter,g_r_local_k,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
!    call write_spectrum_summed_over_k('gw_ndos_r'//string(rank),iter,g_lesser_local_k,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
!    call write_spectrum_summed_over_k('gw_pdos_r'//string(rank),iter,g_greater_local_k,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
!    call write_current_spectrum_summed_over_kz('gw_Jdens_r'//string(rank),iter,jdens_local_k,nen,En,nx*ns,NB,Lx,nk)  
!    call write_transmission_spectrum_k('gw_trR_r'//string(rank),iter,tr_local_k*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
!    call write_transmission_spectrum_k('gw_trL_r'//string(rank),iter,tre_local_k*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
    !!! reduce to sum up k points
    allocate(g_r(nm,nm,nx,nen,nk))  
    allocate(g_greater(nm,nm,nx,nen,nk))
    allocate(g_lesser(nm,nm,nx,nen,nk))  
    allocate(jdens(nb,nb,nx*ns,nen,nk))
    allocate(tot_cur(nb,nb,nx*ns,nk))
    allocate(tot_ecur(nb,nb,nx*ns,nk))
    allocate(tr(nen,nk))  
    allocate(tre(nen,nk))  
    g_r=czero
    g_lesser=czero
    g_greater=czero
    jdens=0.0d0
    tot_cur=0.0d0
    tot_ecur=0.0d0
    tr=0.0d0
    tre=0.0d0
    if (rank == 0) then
      print *,' reduce'
    endif    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(g_r_local_k, g_r, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    g_r=g_r/dble(nnode)        
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(g_lesser_local_k, g_lesser, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    g_lesser=g_lesser/dble(nnode)               
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(g_greater_local_k, g_greater, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    g_greater=g_greater/dble(nnode)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(tr_local_k, tr, nen*nk, MPI_DOUBLE, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    tr=tr/dble(nnode)       
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(tre_local_k, tre, nen*nk, MPI_DOUBLE, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    tre=tre/dble(nnode)       
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(jdens_local_k, jdens, nb*nb*nx*ns*nen*nk, MPI_DOUBLE, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    jdens=jdens/dble(nnode)       
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(tot_cur_local_k, tot_cur, nb*nb*nx*ns*nk, MPI_DOUBLE, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    tot_cur=tot_cur/dble(nnode)       
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(tot_ecur_local_k, tot_ecur, nb*nb*nx*ns*nk, MPI_DOUBLE, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    tot_ecur=tot_ecur/dble(nnode)       
    if (rank == 0) then
      call write_spectrum_summed_over_k('gw_ldos',iter,g_r,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
      call write_spectrum_summed_over_k('gw_ndos',iter,g_lesser,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
      call write_spectrum_summed_over_k('gw_pdos',iter,g_greater,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
      call write_dos_summed_over_k('gw_dos',iter,G_r,nen,En,nk,nx,NB,NS,Lx)
      call write_current_spectrum_summed_over_kz('gw_Jdens',iter,jdens,nen,En,nx*ns,NB,Lx,nk)  
      call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
      call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)      
      call write_transmission_spectrum_k('gw_trR',iter,tr*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
      call write_transmission_spectrum_k('gw_trL',iter,tre*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
      open(unit=101,file='Id_iteration.dat',status='unknown',position='append')
      write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
      close(101)
    endif
    deallocate(g_r,g_greater,g_lesser,tot_cur,tot_ecur,jdens,tr,tre)
    deallocate(jdens_local_k,tot_cur_local_k,tot_ecur_local_k,tr_local_k,tre_local_k,cur_local_k)
    !!!
    g_r_local_k = dcmplx( 0.0d0*dble(g_r_local_k), aimag(g_r_local_k))
    g_lesser_local_k = dcmplx( 0.0d0*dble(g_lesser_local_k), aimag(g_lesser_local_k))
    g_greater_local_k = dcmplx( 0.0d0*dble(g_greater_local_k), aimag(g_greater_local_k))
    !        
    allocate(g_r_local_x(nm,nm,nen,nky*nkz,nnx),stat=ierror)  
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    allocate(g_greater_local_x(nm,nm,nen,nky*nkz,nnx),stat=ierror)
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    allocate(g_lesser_local_x(nm,nm,nen,nky*nkz,nnx),stat=ierror)
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    !!! all-to-all communication to make G local in x
    !!! G(:,:,ix,:,ik=rank) --> buff(:,:,:) --> G(:,:,ix=rank,:,ik)
    if (rank == 0) then
      print *,' all-to-all'
    endif
    allocate(sbuf(nm,nm,nen,nx,nk),stat=ierror)
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    allocate(rbuf(nm,nm,nen,nky*nkz,nnx),stat=ierror)    
    if (ierror .ne. 0) then
      print*,rank,'allocate fail'
      stop
    endif
    !!!
    scount=nm*nm*nen*nnx*nk
    rcount=nm*nm*nen*nnx*nk
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_r_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    if (ierror.ne.0) print *,rank,'ierr=',ierror
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    g_r_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !!
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_lesser_local_k(:,:,ix,:,ik)
      enddo    
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    if (ierror.ne.0) print *,rank,'ierr=',ierror
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    g_lesser_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)        
    !!
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_greater_local_k(:,:,ix,:,ik)
      enddo    
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    if (ierror.ne.0) print *,rank,'ierr=',ierror
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    g_greater_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)       
    call MPI_Barrier(MPI_COMM_WORLD,ierror) 
    deallocate(sbuf)
    deallocate(rbuf)
    !!!
    if (rank == 0) then     
      print *, 'calc P'    
    endif
    nopmax=nen/2-10              
    allocate(P_lesser_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(P_greater_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(P_retarded_local_x(nm,nm,nen,nky*nkz,nnx))
    P_lesser_local_x = dcmplx(0.0d0,0.0d0)
    P_greater_local_x = dcmplx(0.0d0,0.0d0)    
    P_retarded_local_x = dcmplx(0.0d0,0.0d0)        
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)           
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iky,ikyd,ik,ikd,iqz,iqy,iq) shared(ndiag,nopmax,P_lesser_local_x,P_greater_local_x,P_retarded_local_x,nen,En,nm,G_lesser_local_x,G_greater_local_x,G_r_local_x,nx,nkz,nky,nk,nnx)    
    !$omp do         
    do nop=-nopmax,nopmax
      iop=nop+nen/2  
      do iqy=1,nky        
        do iqz=1,nkz
          iq=iqz + (iqy-1)*nkz
          do ix=1,nnx  ! x parallelized by MPI
            do i=1,nm      
              l=max(i-ndiag,1)
              h=min(nm,i+ndiag)   
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
                  do ie = max(nop+1,1),min(nen,nen+nop)           
                    P_lesser_local_x(i,l:h,iop,iq,ix) = P_lesser_local_x(i,l:h,iop,iq,ix)  & 
                          & + G_lesser_local_x(i,l:h,ie,ik,ix) * G_greater_local_x(l:h,i,ie-nop,ikd,ix)
                    P_greater_local_x(i,l:h,iop,iq,ix) = P_greater_local_x(i,l:h,iop,iq,ix) &
                          & + G_greater_local_x(i,l:h,ie,ik,ix) * G_lesser_local_x(l:h,i,ie-nop,ikd,ix)        
                    P_retarded_local_x(i,l:h,iop,iq,ix) = P_retarded_local_x(i,l:h,iop,iq,ix) &
                          & + (G_lesser_local_x(i,l:h,ie,ik,ix) * conjg(G_r_local_x(i,l:h,ie-nop,ikd,ix)) & 
                          & + G_r_local_x(i,l:h,ie,ik,ix) * G_lesser_local_x(l:h,i,ie-nop,ikd,ix))        
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo                    
    !$omp end do 
    !$omp end parallel             
    deallocate(G_r_local_x,G_lesser_local_x,G_greater_local_x)
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nky*nkz)
    P_lesser_local_x  =P_lesser_local_x*dE
    P_greater_local_x =P_greater_local_x*dE
    P_retarded_local_x=P_retarded_local_x*dE
    !
    allocate(P_lesser_local_k(nm,nm,nx,nen,nk))
    allocate(P_greater_local_k(nm,nm,nx,nen,nk))
    allocate(P_retarded_local_k(nm,nm,nx,nen,nk))
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make P local in k
    !!! P(:,:,ix=rank,:,ik) --> buff(:,:,:) --> P(:,:,ix,:,ik=rank)
    allocate(sbuf(nm,nm,nen,nky*nkz,nnx))        
    allocate(rbuf(nm,nm,nen,nx,nk))  
    !!!    
    sbuf(:,:,:,:,:)=P_retarded_local_x(:,:,:,:,:)      
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    do ik=1,nk
      do ix=1,nx    
        P_retarded_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !    
    sbuf(:,:,:,:,:)=P_lesser_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    do ik=1,nk
      do ix=1,nx
        P_lesser_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !
    sbuf(:,:,:,:,:)=P_greater_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    do ik=1,nk
      do ix=1,nx
        P_greater_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    deallocate(sbuf,rbuf)
    !!!    
    deallocate(P_retarded_local_x,P_lesser_local_x,P_greater_local_x)  
    allocate(P_r(nm,nm,nx,nen,nk))  
    allocate(P_greater(nm,nm,nx,nen,nk))
    allocate(P_lesser(nm,nm,nx,nen,nk))
    P_r=czero
    P_lesser=czero
    P_greater=czero    
    if (rank == 0) then
      print *,' reduce'
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(P_retarded_local_k, P_r, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    P_r=P_r/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(P_lesser_local_k, P_lesser, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    P_lesser=P_lesser/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(P_greater_local_k, P_greater, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    P_greater=P_greater/dble(nnode)   
    !
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    if (rank == 0) then  
      call write_spectrum_summed_over_k('gw_PR',iter,P_r,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_PG',iter,P_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    endif
    !
    deallocate(P_r,P_greater,P_lesser)
    if (rank == 0) then  
      print *, 'calc W'                
    endif
    allocate(W_lesser_local_k(nm,nm,nx,nen,nk))
    allocate(W_greater_local_k(nm,nm,nx,nen,nk))
    allocate(W_retarded_local_k(nm,nm,nx,nen,nk))
    do iq=1,nk            
      !$omp parallel default(none) private(nop) shared(rank,nk,iq,nopmax,nen,nm,nx,Vii,V1i,p_lesser_local_k,p_greater_local_k,p_retarded_local_k,w_lesser_local_k,w_greater_local_k,w_retarded_local_k)    
      !$omp do 
      do nop=-nopmax+nen/2,nopmax+nen/2 
        call green_calc_w_full(0,nm,nx,Vii(:,:,:,iq+(rank)*nk),V1i(:,:,:,iq+(rank)*nk),p_lesser_local_k(:,:,:,nop,iq),p_greater_local_k(:,:,:,nop,iq),p_retarded_local_k(:,:,:,nop,iq),w_lesser_local_k(:,:,:,nop,iq),w_greater_local_k(:,:,:,nop,iq),w_retarded_local_k(:,:,:,nop,iq))
      enddo
      !$omp end do 
      !$omp end parallel      
    enddo
    deallocate(P_retarded_local_k,P_lesser_local_k,P_greater_local_k)
    allocate(W_r(nm,nm,nx,nen,nk))  
    allocate(W_greater(nm,nm,nx,nen,nk))
    allocate(W_lesser(nm,nm,nx,nen,nk))
    W_r=czero
    W_lesser=czero
    W_greater=czero   
    if (rank == 0) then
      print *,' reduce'
    endif 
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(W_retarded_local_k, W_r, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    W_r=W_r/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(W_lesser_local_k, W_lesser, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    W_lesser=W_lesser/dble(nnode)   
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(W_greater_local_k, W_greater, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    W_greater=W_greater/dble(nnode)       
    if (rank == 0) then  
        call write_spectrum_summed_over_k('gw_WR',iter,W_r       ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
        call write_spectrum_summed_over_k('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
        call write_spectrum_summed_over_k('gw_WG',iter,W_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))    
    endif
    deallocate(W_r,W_greater,W_lesser)    
    !
    allocate(W_lesser_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(W_greater_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(W_retarded_local_x(nm,nm,nen,nky*nkz,nnx))
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make W local in x
    !!! W(:,:,ix,:,ik=rank) --> buff(:,:,:) --> W(:,:,ix=rank,:,ik)
    allocate(sbuf(nm,nm,nen,nx,nk))
    allocate(rbuf(nm,nm,nen,nky*nkz,nnx))        
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=W_retarded_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    W_retarded_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=W_lesser_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    W_lesser_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !   
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=W_greater_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    W_greater_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    deallocate(sbuf,rbuf)
    deallocate(W_retarded_local_k,W_lesser_local_k,W_greater_local_k)
    !!!
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    if (rank == 0) then 
      print *, 'calc SigGW'
    endif
    !        
    allocate(g_r_local_x(nm,nm,nen,nky*nkz,nnx))  
    allocate(g_greater_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(g_lesser_local_x(nm,nm,nen,nky*nkz,nnx))
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make G local in x
    !!! G(:,:,ix,:,ik=rank) --> buff(:,:,:) --> G(:,:,ix=rank,:,ik)
    allocate(sbuf(nm,nm,nen,nx,nk))
    allocate(rbuf(nm,nm,nen,nky*nkz,nnx))        
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_r_local_k(:,:,ix,:,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    g_r_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_lesser_local_k(:,:,ix,:,ik)
      enddo
    enddo    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    g_lesser_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    do ik=1,nk
      do ix=1,nx
        sbuf(:,:,:,ix,ik)=g_greater_local_k(:,:,ix,:,ik)
      enddo
    enddo    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    g_greater_local_x(:,:,:,:,:)=rbuf(:,:,:,:,:)    
    !
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    deallocate(sbuf,rbuf)
    !!!
    allocate(Sigma_lesser_new_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(Sigma_greater_new_local_x(nm,nm,nen,nky*nkz,nnx))
    allocate(Sigma_r_new_local_x(nm,nm,nen,nky*nkz,nnx))
    nopmax=nen/2-10
    Sigma_greater_new_local_x = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new_local_x = dcmplx(0.0d0,0.0d0)
    Sigma_r_new_local_x = dcmplx(0.0d0,0.0d0)    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,w_lesser_local_x,w_greater_local_x,w_retarded_local_x,sigma_lesser_new_local_x,sigma_greater_new_local_x,sigma_r_new_local_x,nen,En,nm,G_lesser_local_x,G_greater_local_x,G_r_local_x,nx,nkz,nky,nk,nnx)    
    !$omp do    
    do ie=1,nen      
      do iky=1,nky
        do ikz=1,nkz
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then   
              ik=ikz+(iky-1)*nkz
              do iqy=1,nky
                do iqz=1,nkz              
                  iq=iqz + (iqy-1)*nkz
                  iop=nop+nen/2                
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2            
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky        
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1        
                  ikd=ikzd + (ikyd-1)*nkz
                  do ix=1,nnx ! x via MPI parallel
                    do i=1,nm
                      l=max(i-ndiag,1)
                      h=min(nm,i+ndiag)                                                
                      Sigma_lesser_new_local_x(i,l:h,ie,ik,ix)=Sigma_lesser_new_local_x(i,l:h,ie,ik,ix)&
                          & + G_lesser_local_x(i,l:h,ie-nop,ikd,ix) * W_lesser_local_x(i,l:h,iop,iq,ix)
                      Sigma_greater_new_local_x(i,l:h,ie,ik,ix)=Sigma_greater_new_local_x(i,l:h,ie,ik,ix)&
                          & + G_greater_local_x(i,l:h,ie-nop,ikd,ix) * W_greater_local_x(i,l:h,iop,iq,ix)
                      Sigma_r_new_local_x(i,l:h,ie,ik,ix)=Sigma_r_new_local_x(i,l:h,ie,ik,ix)&
                          & + G_lesser_local_x(i,l:h,ie-nop,ikd,ix) * W_retarded_local_x(i,l:h,iop,iq,ix) &
                          & + G_r_local_x(i,l:h,ie-nop,ikd,ix) * W_lesser_local_x(i,l:h,iop,iq,ix) &
                          & + G_r_local_x(i,l:h,ie-nop,ikd,ix) * W_retarded_local_x(i,l:h,iop,iq,ix)    
                     enddo
                  enddo
                enddo
              enddo
            endif            
          enddo
        enddo
      enddo      
    enddo
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nky*nkz) 
    Sigma_lesser_new_local_x = Sigma_lesser_new_local_x  * dE
    Sigma_greater_new_local_x= Sigma_greater_new_local_x * dE
    Sigma_r_new_local_x=Sigma_r_new_local_x* dE
    Sigma_r_new_local_x = dcmplx( dble(Sigma_r_new_local_x), aimag(Sigma_greater_new_local_x-Sigma_lesser_new_local_x)/2.0d0 )    
    deallocate(W_retarded_local_x,W_lesser_local_x,W_greater_local_x)
    deallocate(G_r_local_x,G_lesser_local_x,G_greater_local_x)
    if (rank == 0) then
      print *,' all-to-all'
    endif
    !!! all-to-all communication to make Sigma local in k
    !!! Sig(:,:,ix=rank,:,ik) --> buff(:,:,:) --> Sig(:,:,ix,:,ik=rank)
    allocate(sigma_lesser_new_local_k(nm,nm,nx,nen,nk))
    allocate(sigma_greater_new_local_k(nm,nm,nx,nen,nk))
    allocate(sigma_r_new_local_k(nm,nm,nx,nen,nk)) 
    allocate(rbuf(nm,nm,nen,nx,nk))
    allocate(sbuf(nm,nm,nen,nky*nkz,nnx))        
    !    
    sbuf(:,:,:,:,:)=Sigma_r_new_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    do ik=1,nk
      do ix=1,nx      
        Sigma_r_new_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !
    sbuf(:,:,:,:,:)=Sigma_lesser_new_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    do ik=1,nk
      do ix=1,nx
        Sigma_lesser_new_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    !
    sbuf(:,:,:,:,:)=Sigma_greater_new_local_x(:,:,:,:,:)    
    call MPI_Barrier(MPI_COMM_WORLD,ierror)    
    call MPI_ALLTOALL(sbuf, scount, MPI_C_DOUBLE_COMPLEX, rbuf, rcount, MPI_C_DOUBLE_COMPLEX, MPI_COMM_WORLD, ierror)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    do ik=1,nk    
      do ix=1,nx
        Sigma_greater_new_local_k(:,:,ix,:,ik)=rbuf(:,:,:,ix,ik)
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    deallocate(rbuf,sbuf)
    deallocate(Sigma_r_new_local_x,Sigma_lesser_new_local_x,Sigma_greater_new_local_x)
    !!!
    ! symmetrize the selfenergies
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    do ie=1,nen
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new_local_k(:,:,ix,ie,ik))
          Sigma_r_new_local_k(:,:,ix,ie,ik) = (Sigma_r_new_local_k(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new_local_k(:,:,ix,ie,ik))
          Sigma_lesser_new_local_k(:,:,ix,ie,ik) = (Sigma_lesser_new_local_k(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new_local_k(:,:,ix,ie,ik))
          Sigma_greater_new_local_k(:,:,ix,ie,ik) = (Sigma_greater_new_local_k(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo
    ! mixing with the previous one
    Sigma_r_gw_local_k= Sigma_r_gw_local_k+ alpha_mix * (Sigma_r_new_local_k -Sigma_r_gw_local_k)
    Sigma_lesser_gw_local_k  = Sigma_lesser_gw_local_k+ alpha_mix * (Sigma_lesser_new_local_k -Sigma_lesser_gw_local_k)
    Sigma_greater_gw_local_k = Sigma_greater_gw_local_k+ alpha_mix * (Sigma_greater_new_local_k -Sigma_greater_gw_local_k)  
    !
    allocate(Ispec(nm,nm,nx,nen))
    allocate(Itot(nm,nm,nx))
    allocate(Ispec_ik(nm,nm,nx,nen))
    allocate(Itot_ik(nm,nm,nx))
    !!!!! calculate collision integral
    Ispec=czero
    Itot=czero
    do ik=1,nk
      call calc_block_collision(sigma_lesser_gw_local_k(:,:,:,:,ik),sigma_greater_gw_local_k(:,:,:,:,ik),G_lesser_local_k(:,:,:,:,ik),G_greater_local_k(:,:,:,:,ik),nen,en,spindeg,nm,nx,Itot_ik,Ispec_ik)
      Ispec=Ispec+Ispec_ik
      Itot=Itot+Itot_ik
    enddo        
    ! 
    allocate(sigma_r_gw(nm,nm,nx,nen,nk))  
    allocate(sigma_greater_gw(nm,nm,nx,nen,nk))
    allocate(sigma_lesser_gw(nm,nm,nx,nen,nk))      
    sigma_r_gw=czero
    sigma_lesser_gw=czero
    sigma_greater_gw=czero    
    Ispec_ik=czero
    if (rank == 0) then
      print *,' reduce'
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(sigma_r_gw_local_k, sigma_r_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    sigma_r_gw=sigma_r_gw/dble(nnode)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(sigma_lesser_gw_local_k, sigma_lesser_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    sigma_lesser_gw=sigma_lesser_gw/dble(nnode)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(sigma_greater_gw_local_k, sigma_greater_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    sigma_greater_gw=sigma_greater_gw/dble(nnode)
    call MPI_Barrier(MPI_COMM_WORLD,ierror)
    call MPI_Reduce(Ispec, Ispec_ik, nm*nm*nx*nen, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    Ispec_ik=Ispec_ik/dble(nnode)
    !
    if ( rank==0 ) then
      call write_spectrum_summed_over_k('gw_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum_summed_over_k('gw_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
      call write_spectrum('gw_Scat',iter,Ispec_ik,nen,En,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/)) 
    endif    
    !
    if (iter>=(niter-5)) then
      if ( rank==0 ) then
        print *, 'calc SigEPhoton'
      endif      
      do ik=1,nk
        call MPI_Barrier(MPI_COMM_WORLD,ierror)
        !print *, ' ik=', ik+(rank)*nk,'/',nky*nkz
        call calc_sigma_ephoton_monochromatic(nm,nx,nen,En,nopphot,Mii(:,:,:,ik+(rank)*nk),G_lesser_local_k(:,:,:,:,ik),G_greater_local_k(:,:,:,:,ik),Sigma_r_new_local_k(:,:,:,:,ik),Sigma_lesser_new_local_k(:,:,:,:,ik),Sigma_greater_new_local_k(:,:,:,:,ik),pre_fact,intensity,hw)
      enddo      
      Sigma_r_gw_local_k       = Sigma_r_gw_local_k + Sigma_r_new_local_k 
      Sigma_lesser_gw_local_k  = Sigma_lesser_gw_local_k + Sigma_lesser_new_local_k 
      Sigma_greater_gw_local_k = Sigma_greater_gw_local_k + Sigma_greater_new_local_k 
      !!!!! calculate collision integral      
      Ispec=czero
      Itot=czero
      do ik=1,nk
        call calc_block_collision(sigma_lesser_new_local_k(:,:,:,:,ik),sigma_greater_new_local_k(:,:,:,:,ik),G_lesser_local_k(:,:,:,:,ik),G_greater_local_k(:,:,:,:,ik),nen,en,spindeg,nm,nx,Itot_ik,Ispec_ik)
        Ispec=Ispec+Ispec_ik
        Itot=Itot+Itot_ik
      enddo
      !!!!! calculate light absorption 
      if (labs) then        
        allocate(Pi_retarded(nm,nm,nx))        
        allocate(Pi_retarded_ik(nm,nm,nx))
        allocate(Pi_lesser_ik(nm,nm,nx))
        allocate(Pi_greater_ik(nm,nm,nx))                
        do i=1,floor(4.0d0/(En(2)-En(1)))
          Pi_retarded=czero        
          do ik=1,nk
            call calc_pi_ephoton_monochromatic(nm,nx,nen,En,i,Mii(:,:,:,ik+(rank)*nk),G_lesser_local_k(:,:,:,:,ik),G_greater_local_k(:,:,:,:,ik),Pi_retarded_ik,Pi_lesser_ik,Pi_greater_ik)  
            Pi_retarded=Pi_retarded+Pi_retarded_ik
          enddo
          Pi_retarded_ik=czero
          call MPI_Reduce(Pi_retarded, Pi_retarded_ik, nm*nm*nx, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
          Pi_retarded_ik=Pi_retarded_ik/dble(nnode)
          if ( rank==0 ) then
            call write_trace('eph_absorp',iter,Pi_retarded_ik,nx,NB,Lx,(/1.0d0/dble(nk),-1.0d0/dble(nk)/),E=dble(i)*(En(2)-En(1)))          
          endif
          call MPI_Barrier(MPI_COMM_WORLD,ierror)
        enddo        
        deallocate(Pi_retarded)
        deallocate(Pi_lesser_ik,Pi_greater_ik,Pi_retarded_ik)
      endif
      ! 
      sigma_r_gw=czero
      sigma_lesser_gw=czero
      sigma_greater_gw=czero   
      Ispec_ik=czero 
      if (rank == 0) then
        print *,' reduce'
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      call MPI_Reduce(sigma_r_new_local_k, sigma_r_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                  MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      sigma_r_gw=sigma_r_gw/dble(nnode)
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      call MPI_Reduce(sigma_lesser_new_local_k, sigma_lesser_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                  MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      sigma_lesser_gw=sigma_lesser_gw/dble(nnode)
      call MPI_Barrier(MPI_COMM_WORLD,ierror)
      call MPI_Reduce(sigma_greater_new_local_k, sigma_greater_gw, nm*nm*nx*nen*nk, MPI_C_DOUBLE_COMPLEX, &
                  MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      sigma_greater_gw=sigma_greater_gw/dble(nnode) 
      call MPI_Reduce(Ispec, Ispec_ik, nm*nm*nx*nen, MPI_C_DOUBLE_COMPLEX, &
                MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      Ispec_ik=Ispec_ik/dble(nnode)
      if ( rank==0 ) then
          call write_spectrum_summed_over_k('eph_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
          call write_spectrum_summed_over_k('eph_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
          call write_spectrum_summed_over_k('eph_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
          call write_spectrum('eph_Scat',iter,Ispec_ik,nen,En,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/)) 
      endif      
    endif  
    deallocate(Ispec,Itot,Ispec_ik,Itot_ik)  
    !!! make sure self-energy is continuous near leads by copying edge block
    do ix=1,2
      Sigma_r_gw_local_k(:,:,ix,:,:)=Sigma_r_gw_local_k(:,:,3,:,:)
      Sigma_lesser_gw_local_k(:,:,ix,:,:)=Sigma_lesser_gw_local_k(:,:,3,:,:)
      Sigma_greater_gw_local_k(:,:,ix,:,:)=Sigma_greater_gw_local_k(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw_local_k(:,:,nx-ix+1,:,:)=Sigma_r_gw_local_k(:,:,nx-2,:,:)
      Sigma_lesser_gw_local_k(:,:,nx-ix+1,:,:)=Sigma_lesser_gw_local_k(:,:,nx-2,:,:)
      Sigma_greater_gw_local_k(:,:,nx-ix+1,:,:)=Sigma_greater_gw_local_k(:,:,nx-2,:,:)
    enddo    
    deallocate(sigma_r_new_local_k,sigma_lesser_new_local_k,sigma_greater_new_local_k)     
    deallocate(sigma_r_gw,sigma_lesser_gw,sigma_greater_gw)
  enddo  
  deallocate(g_r_local_k,g_lesser_local_k,g_greater_local_k)  
  deallocate(sigma_lesser_gw_local_k,sigma_greater_gw_local_k,sigma_r_gw_local_k)         
end subroutine green_rgf_solve_gw_ephoton_3d_mpi_kpar




subroutine green_rgf_solve_gw_3d(alpha_mix,niter,NB,NS,nm,nx,nky,nkz,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg)
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag,nky,nkz
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx,nky*nkz)::Hii,H1i,Vii,V1i
  !complex(8), intent(in):: V(nm*nx,nm*nx,nky*nkz)
  ! -------- local variables
  complex(8),allocatable,dimension(:,:,:,:,:)::g_r,g_greater,g_lesser, g_r_i1
  real(8),allocatable,dimension(:,:,:,:,:)::cur,jdens
  real(8),allocatable,dimension(:,:,:,:)::tot_cur,tot_ecur
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  complex(8),allocatable,dimension(:,:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser,P_greater,P_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:,:,:)::W_lesser_i1,W_greater_i1,W_retarded_i1  
  complex(8),allocatable::Ispec(:,:,:,:),Itot(:,:,:),Ispec_ik(:,:,:,:),Itot_ik(:,:,:)
  real(8)::tr(nen,nky*nkz),tre(nen,nky*nkz)
  integer::ie,iter,i,ix,nopmax,nop,iop,l,h,io
  integer::ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,nk
  complex(8)::dE, B(nm,nm)
  !
  print *,'======== green_rgf_solve_gw_3D ========'
  nk=nky*nkz
  allocate(g_r(nm,nm,nx,nen,nk))
  allocate(g_r_i1(nm,nm,nx,nen,nk))
  allocate(g_greater(nm,nm,nx,nen,nk))
  allocate(g_lesser(nm,nm,nx,nen,nk))
  allocate(cur(nm,nm,nx,nen,nk))
  allocate(jdens(nb,nb,nx*ns,nen,nk))
  allocate(tot_cur(nb,nb,nx*ns,nk))
  allocate(tot_ecur(nb,nb,nx*ns,nk))
  allocate(sigma_lesser_gw(nm,nm,nx,nen,nk))
  allocate(sigma_greater_gw(nm,nm,nx,nen,nk))
  allocate(sigma_r_gw(nm,nm,nx,nen,nk))
  allocate(sigma_lesser_new(nm,nm,nx,nen,nk))
  allocate(sigma_greater_new(nm,nm,nx,nen,nk))
  allocate(sigma_r_new(nm,nm,nx,nen,nk))
  allocate(P_lesser(nm,nm,nx,nen,nk))
  allocate(P_greater(nm,nm,nx,nen,nk))
  allocate(P_retarded(nm,nm,nx,nen,nk))
!  allocate(P_lesser_1i(nm,nm,nx,nen,nk))
!  allocate(P_greater_1i(nm,nm,nx,nen,nk))
!  allocate(P_retarded_1i(nm,nm,nx,nen,nk))
  allocate(W_lesser(nm,nm,nx,nen,nk))
  allocate(W_greater(nm,nm,nx,nen,nk))
  allocate(W_retarded(nm,nm,nx,nen,nk))
!  allocate(W_lesser_i1(nm,nm,nx,nen,nk))
!  allocate(W_greater_i1(nm,nm,nx,nen,nk))
!  allocate(W_retarded_i1(nm,nm,nx,nen,nk))  
  allocate(Ispec(nm,nm,nx,nen))
  allocate(Itot(nm,nm,nx))
  allocate(Ispec_ik(nm,nm,nx,nen))
  allocate(Itot_ik(nm,nm,nx))
  do iter=0,niter
    print *,'+ iter=',iter
    !
    print *, 'calc G'      
    do ik=1,nk
      print *, ' ik=', ik,'/',nk
      !$omp parallel default(none) private(ie) shared(ik,nen,temp,nm,nx,en,mu,Hii,H1i,sigma_lesser_gw,sigma_greater_gw,sigma_r_gw,g_lesser,g_greater,g_r,tre,tr,cur,g_r_i1)    
      !$omp do 
      do ie=1,nen     
        !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
        call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii(:,:,:,ik),H1i(:,:,:,ik),sigma_lesser_gw(:,:,:,ie,ik),sigma_greater_gw(:,:,:,ie,ik),&
                        & sigma_r_gw(:,:,:,ie,ik),g_lesser(:,:,:,ie,ik),g_greater(:,:,:,ie,ik),g_r(:,:,:,ie,ik),tr(ie,ik),&
                        & tre(ie,ik),cur(:,:,:,ie,ik),g_r_i1(:,:,:,ie,ik)) 
      enddo
      !$omp end do 
      !$omp end parallel
      call calc_block_current(Hii(:,:,:,ik),G_lesser(:,:,:,:,ik),cur(:,:,:,:,ik),nen,en,spindeg,nb,ns,nm,nx,tot_cur(:,:,:,ik),tot_ecur(:,:,:,ik),jdens(:,:,:,:,ik))
    enddo
    !
    call write_spectrum_summed_over_k('gw_ldos',iter,g_r,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
    call write_spectrum_summed_over_k('gw_ndos',iter,g_lesser,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
    call write_spectrum_summed_over_k('gw_pdos',iter,g_greater,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
    call write_current_spectrum_summed_over_kz('gw_Jdens',iter,jdens,nen,En,nx*ns,NB,Lx,nk)        
    call write_current_summed_over_k('gw_I',iter,tot_cur,nx*ns,NB,Lx,nk)
    call write_current_summed_over_k('gw_EI',iter,tot_ecur,nx*ns,NB,Lx,nk)
    call write_transmission_spectrum_k('gw_trR',iter,tr*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
    call write_transmission_spectrum_k('gw_trL',iter,tre*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk),nen,En,nk)
    open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)/dble(nk)
    close(101)
    g_r = dcmplx( 0.0d0*dble(g_r), aimag(g_r))
    g_lesser = dcmplx( 0.0d0*dble(g_lesser), aimag(g_lesser))
    g_greater = dcmplx( 0.0d0*dble(g_greater), aimag(g_greater))
    !        
    print *, 'calc P'    
    nopmax=nen/2-10      
    P_lesser = dcmplx(0.0d0,0.0d0)
    P_greater = dcmplx(0.0d0,0.0d0)    
    P_retarded = dcmplx(0.0d0,0.0d0)    
!    P_lesser_1i = dcmplx(0.0d0,0.0d0)
!    P_greater_1i = dcmplx(0.0d0,0.0d0)    
!    P_retarded_1i = dcmplx(0.0d0,0.0d0)  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do         
    do nop=-nopmax,nopmax
      iop=nop+nen/2  
      do iqy=1,nky        
        do iqz=1,nkz
          iq=iqz + (iqy-1)*nkz
          do ix=1,nx
            do i=1,nm      
              l=max(i-ndiag,1)
              h=min(nm,i+ndiag)   
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
                  do ie = max(nop+1,1),min(nen,nen+nop)           
                    P_lesser(i,l:h,ix,iop,iq) = P_lesser(i,l:h,ix,iop,iq) + G_lesser(i,l:h,ix,ie,ik) * G_greater(l:h,i,ix,ie-nop,ikd)
                    P_greater(i,l:h,ix,iop,iq) = P_greater(i,l:h,ix,iop,iq) + G_greater(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd)        
                    P_retarded(i,l:h,ix,iop,iq) = P_retarded(i,l:h,ix,iop,iq) + (G_lesser(i,l:h,ix,ie,ik) * conjg(G_r(i,l:h,ix,ie-nop,ikd)) & 
                                              & + G_r(i,l:h,ix,ie,ik) * G_lesser(l:h,i,ix,ie-nop,ikd))        
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo      
    enddo          
    !$omp end do 
    !$omp end parallel
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg /dble(nk)
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=P_retarded*dE
    call write_spectrum_summed_over_k('gw_PR',iter,P_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_PG',iter,P_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc W'                
    do iq=1,nky*nkz
      print *, ' iq=', iq,'/',nk
      !$omp parallel default(none) private(nop) shared(iq,nopmax,nen,nm,nx,Vii,V1i,p_lesser,p_greater,p_retarded,w_lesser,w_greater,w_retarded)    
      !$omp do 
      do nop=-nopmax+nen/2,nopmax+nen/2 
        !call green_calc_w_full(0,nm,nx,Vii(:,:,:,iq),V1i(:,:,:,iq),p_lesser(:,:,:,nop,iq),p_greater(:,:,:,nop,iq),p_retarded(:,:,:,nop,iq),w_lesser(:,:,:,nop,iq),w_greater(:,:,:,nop,iq),w_retarded(:,:,:,nop,iq))
        
        call green_rgf_calc_w(0,nm,nx,Vii(:,:,:,iq),V1i(:,:,:,iq),p_lesser(:,:,:,nop,iq),p_greater(:,:,:,nop,iq),p_retarded(:,:,:,nop,iq),w_lesser(:,:,:,nop,iq),w_greater(:,:,:,nop,iq),w_retarded(:,:,:,nop,iq))
      enddo
      !$omp end do 
      !$omp end parallel
    enddo
    call write_spectrum_summed_over_k('gw_WR',iter,W_retarded,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_WG',iter,W_greater ,nen,En-en(nen/2),nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc SigGW'
    nopmax=nen/2-10
    Sigma_greater_new = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new = dcmplx(0.0d0,0.0d0)
    Sigma_r_new = dcmplx(0.0d0,0.0d0)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd) shared(ndiag,nopmax,w_lesser,w_greater,w_retarded,sigma_lesser_new,sigma_greater_new,sigma_r_new,nen,En,nm,G_lesser,G_greater,G_r,nx,nkz,nky,nk)    
    !$omp do    
    do ie=1,nen      
      do iky=1,nky
        do ikz=1,nkz
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then   
              ik=ikz+(iky-1)*nkz
              do iqy=1,nky
                do iqz=1,nkz              
                  iq=iqz + (iqy-1)*nkz
                  iop=nop+nen/2                
                  ikzd=ikz-iqz + nkz/2
                  ikyd=iky-iqy + nky/2            
                  if (ikzd<1)   ikzd=ikzd+nkz
                  if (ikzd>nkz) ikzd=ikzd-nkz
                  if (ikyd<1)   ikyd=ikyd+nky
                  if (ikyd>nky) ikyd=ikyd-nky        
                  if (nky==1)   ikyd=1
                  if (nkz==1)   ikzd=1        
                  ikd=ikzd + (ikyd-1)*nkz
                  do ix=1,nx
                    do i=1,nm
                      l=max(i-ndiag,1)
                      h=min(nm,i+ndiag)                                                
                      Sigma_lesser_new(i,l:h,ix,ie,ik)=Sigma_lesser_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq)
                      Sigma_greater_new(i,l:h,ix,ie,ik)=Sigma_greater_new(i,l:h,ix,ie,ik)+G_greater(i,l:h,ix,ie-nop,ikd)*W_greater(i,l:h,ix,iop,iq)
                      Sigma_r_new(i,l:h,ix,ie,ik)=Sigma_r_new(i,l:h,ix,ie,ik)+G_lesser(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_lesser(i,l:h,ix,iop,iq) &
                                                          &  +G_r(i,l:h,ix,ie-nop,ikd)*W_retarded(i,l:h,ix,iop,iq)    
                     enddo
                  enddo
                enddo
              enddo
            endif            
          enddo
        enddo
      enddo      
    enddo
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nk) 
    Sigma_lesser_new = Sigma_lesser_new  * dE
    Sigma_greater_new= Sigma_greater_new * dE
    Sigma_r_new=Sigma_r_new* dE
    Sigma_r_new = dcmplx( dble(Sigma_r_new), aimag(Sigma_greater_new-Sigma_lesser_new)/2.0d0 )    
    ! symmetrize the selfenergies
    do ie=1,nen
      do ik=1,nk
        do ix=1,nx
          B(:,:)=transpose(Sigma_r_new(:,:,ix,ie,ik))
          Sigma_r_new(:,:,ix,ie,ik) = (Sigma_r_new(:,:,ix,ie,ik) + B(:,:))/2.0d0    
          B(:,:)=transpose(Sigma_lesser_new(:,:,ix,ie,ik))
          Sigma_lesser_new(:,:,ix,ie,ik) = (Sigma_lesser_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
          B(:,:)=transpose(Sigma_greater_new(:,:,ix,ie,ik))
          Sigma_greater_new(:,:,ix,ie,ik) = (Sigma_greater_new(:,:,ix,ie,ik) + B(:,:))/2.0d0
        enddo
      enddo
    enddo
    ! mixing with the previous one
    Sigma_r_gw= Sigma_r_gw+ alpha_mix * (Sigma_r_new -Sigma_r_gw)
    Sigma_lesser_gw  = Sigma_lesser_gw+ alpha_mix * (Sigma_lesser_new -Sigma_lesser_gw)
    Sigma_greater_gw = Sigma_greater_gw+ alpha_mix * (Sigma_greater_new -Sigma_greater_gw)  
    ! make sure self-energy is continuous near leads (by copying edge block)
    do ix=1,2
      Sigma_r_gw(:,:,ix,:,:)=Sigma_r_gw(:,:,3,:,:)
      Sigma_lesser_gw(:,:,ix,:,:)=Sigma_lesser_gw(:,:,3,:,:)
      Sigma_greater_gw(:,:,ix,:,:)=Sigma_greater_gw(:,:,3,:,:)
    enddo
    do ix=1,2
      Sigma_r_gw(:,:,nx-ix+1,:,:)=Sigma_r_gw(:,:,nx-2,:,:)
      Sigma_lesser_gw(:,:,nx-ix+1,:,:)=Sigma_lesser_gw(:,:,nx-2,:,:)
      Sigma_greater_gw(:,:,nx-ix+1,:,:)=Sigma_greater_gw(:,:,nx-2,:,:)
    enddo
    call write_spectrum_summed_over_k('gw_SigR',iter,Sigma_r_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigL',iter,Sigma_lesser_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum_summed_over_k('gw_SigG',iter,Sigma_greater_gw,nen,En,nk,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !!!!! calculate collision integral
    Ispec=czero
    Itot=czero
    do ik=1,nk
      call calc_block_collision(sigma_lesser_gw(:,:,:,:,ik),sigma_greater_gw(:,:,:,:,ik),G_lesser(:,:,:,:,ik),G_greater(:,:,:,:,ik),nen,en,spindeg,nm,nx,Itot_ik,Ispec_ik)
      Ispec=Ispec+Ispec_ik
      Itot=Itot+Itot_ik
    enddo
    call write_spectrum('gw_Scat',iter,Ispec,nen,En,nx,NB,NS,Lx,(/1.0d0/dble(nk),1.0d0/dble(nk)/))    
    !
  enddo
  deallocate(g_r,g_lesser,g_greater,cur)
  deallocate(sigma_lesser_gw,sigma_greater_gw,sigma_r_gw)   
  deallocate(sigma_lesser_new,sigma_greater_new,sigma_r_new)   
  deallocate(P_retarded,P_lesser,P_greater)
  !deallocate(P_retarded_1i,P_lesser_1i,P_greater_1i)
  deallocate(W_retarded,W_lesser,W_greater)
  !deallocate(W_retarded_i1,W_lesser_i1,W_greater_i1)  
  deallocate(tot_cur,tot_ecur,jdens)
  deallocate(Ispec,Itot)
  deallocate(Ispec_ik,Itot_ik)
end subroutine green_rgf_solve_gw_3d


subroutine green_rgf_solve_gw_1d(alpha_mix,niter,NB,NS,nm,nx,ndiag,Lx,nen,en,temp,mu,Hii,H1i,Vii,V1i,spindeg)
  integer,intent(in)::nm,nx,nen,niter,NB,NS,ndiag
  real(8),intent(in)::en(nen),temp(2),mu(2),Lx,alpha_mix,spindeg
  complex(8),intent(in),dimension(nm,nm,nx)::Hii,H1i,Vii,V1i
  ! -------- local variables
  complex(8),allocatable,dimension(:,:,:,:)::g_r,g_greater,g_lesser, g_r_i1
  real(8),allocatable,dimension(:,:,:,:)::cur,jdens
  real(8),allocatable,dimension(:,:,:)::tot_cur,tot_ecur
  complex(8),allocatable,dimension(:,:,:,:)::sigma_lesser_gw,sigma_greater_gw,sigma_r_gw
  complex(8),allocatable,dimension(:,:,:,:)::sigma_lesser_new,sigma_greater_new,sigma_r_new
  complex(8),allocatable,dimension(:,:,:,:)::P_lesser,P_greater,P_retarded
  !complex(8),allocatable,dimension(:,:,:,:)::P_lesser_1i,P_greater_1i,P_retarded_1i
  complex(8),allocatable,dimension(:,:,:,:)::W_lesser,W_greater,W_retarded
  !complex(8),allocatable,dimension(:,:,:,:)::W_lesser_i1,W_greater_i1,W_retarded_i1
  complex(8),allocatable,dimension(:,:,:)::Sii
  real(8)::tr(nen),tre(nen)
  integer::ie,iter,i,ix,nopmax,nop,iop,l,h,io
  complex(8)::dE
  !
  print *,'====== green_rgf_solve_gw_1D ======'
  allocate(g_r(nm,nm,nx,nen))
  allocate(g_r_i1(nm,nm,nx,nen))
  allocate(g_greater(nm,nm,nx,nen))
  allocate(g_lesser(nm,nm,nx,nen))
  allocate(cur(nm,nm,nx,nen))
  allocate(jdens(nb,nb,nx*ns,nen))
  allocate(tot_cur(nb,nb,nx*ns))
  allocate(tot_ecur(nb,nb,nx*ns))
  allocate(sigma_lesser_gw(nm,nm,nx,nen))
  allocate(sigma_greater_gw(nm,nm,nx,nen))
  allocate(sigma_r_gw(nm,nm,nx,nen))
  allocate(sigma_lesser_new(nm,nm,nx,nen))
  allocate(sigma_greater_new(nm,nm,nx,nen))
  allocate(sigma_r_new(nm,nm,nx,nen))
  allocate(P_lesser(nm,nm,nx,nen))
  allocate(P_greater(nm,nm,nx,nen))
  allocate(P_retarded(nm,nm,nx,nen))
  !allocate(P_lesser_1i(nm,nm,nx,nen))
  !allocate(P_greater_1i(nm,nm,nx,nen))
  !allocate(P_retarded_1i(nm,nm,nx,nen))
  allocate(W_lesser(nm,nm,nx,nen))
  allocate(W_greater(nm,nm,nx,nen))
  allocate(W_retarded(nm,nm,nx,nen))
  !allocate(W_lesser_i1(nm,nm,nx,nen))
  !allocate(W_greater_i1(nm,nm,nx,nen))
  !allocate(W_retarded_i1(nm,nm,nx,nen))  
  do iter=0,niter
    print *,'+ iter=',iter
    !
    print *, 'calc G'      
    do ie=1,nen     
      if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
      call green_RGF_RS(TEMP,nm,nx,En(ie),mu,Hii,H1i,sigma_lesser_gw(:,:,:,ie),sigma_greater_gw(:,:,:,ie),&
                      & sigma_r_gw(:,:,:,ie),g_lesser(:,:,:,ie),g_greater(:,:,:,ie),g_r(:,:,:,ie),tr(ie),&
                      & tre(ie),cur(:,:,:,ie),g_r_i1(:,:,:,ie))
      call calc_block_current(Hii(:,:,:),G_lesser(:,:,:,:),cur(:,:,:,:),nen,en,spindeg,nb,ns,nm,nx,tot_cur(:,:,:),tot_ecur(:,:,:),jdens(:,:,:,:))                      
    enddo
    !
    call write_spectrum('gw_ldos',iter,g_r,nen,En,nx,NB,NS,Lx,(/1.0d0,-2.0d0/))  
    call write_spectrum('gw_ndos',iter,g_lesser,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))       
    call write_spectrum('gw_pdos',iter,g_greater,nen,En,nx,NB,NS,Lx,(/1.0d0,-1.0d0/)) 
    call write_current_spectrum('gw_Jdens',iter,jdens,nen,En,nx*ns,NB,Lx)              
    call write_transmission_spectrum('gw_trR',iter,tr(:)*spindeg,nen,En)
    call write_transmission_spectrum('gw_trL',iter,tre(:)*spindeg,nen,En)
    open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, sum(tre)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg), sum(tr)*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)
    close(101)
    g_r = dcmplx( 0.0d0*dble(g_r), aimag(g_r))
    g_lesser = dcmplx( 0.0d0*dble(g_lesser), aimag(g_lesser))
    g_greater = dcmplx( 0.0d0*dble(g_greater), aimag(g_greater))
    !        
    print *, 'calc P'    
    nopmax=nen/2-10  
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	 * spindeg   
    P_lesser(:,:,:,:) = dcmplx(0.0d0,0.0d0)
    P_greater(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
    P_retarded(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
!    P_lesser_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)
!    P_greater_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)    
!    P_retarded_1i(:,:,:,:) = dcmplx(0.0d0,0.0d0)  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    !$omp parallel default(none) private(ix,l,h,iop,nop,ie,i) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm,G_lesser,G_greater,G_r,nx)    
    !$omp do 
    do ix=1,nx
      do i=1,nm      
        l=max(i-ndiag,1)
        h=min(nm,i+ndiag)   
        do nop=-nopmax,nopmax
          iop=nop+nen/2  
          do ie = max(nop+1,1),min(nen,nen+nop)           
              P_lesser(i,l:h,ix,iop) = P_lesser(i,l:h,ix,iop) + G_lesser(i,l:h,ix,ie) * G_greater(l:h,i,ix,ie-nop)
              P_greater(i,l:h,ix,iop) = P_greater(i,l:h,ix,iop) + G_greater(i,l:h,ix,ie) * G_lesser(l:h,i,ix,ie-nop)        
              P_retarded(i,l:h,ix,iop) = P_retarded(i,l:h,ix,iop) + (G_lesser(i,l:h,ix,ie) * conjg(G_r(i,l:h,ix,ie-nop)) & 
                                       & + G_r(i,l:h,ix,ie) * G_lesser(l:h,i,ix,ie-nop))        
          enddo
        enddo
      enddo
    enddo
    !$omp end do 
    !$omp end parallel
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=P_retarded*dE
    call write_spectrum('gw_PR',iter,P_retarded,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_PL',iter,P_lesser  ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_PG',iter,P_greater ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc W'
    do ie=1,nen
      !if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',nen
      call green_calc_w_full(0,nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie))
      !call green_rgf_calc_w(nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),p_lesser_1i(:,:,:,ie),p_greater_1i(:,:,:,ie),p_retarded_1i(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie),w_lesser_i1(:,:,:,ie),w_greater_i1(:,:,:,ie),w_retarded_i1(:,:,:,ie))
      !call green_rgf_w(nm,nx,Vii,V1i,p_lesser(:,:,:,ie),p_greater(:,:,:,ie),p_retarded(:,:,:,ie),p_lesser_1i(:,:,:,ie),p_greater_1i(:,:,:,ie),p_retarded_1i(:,:,:,ie),w_lesser(:,:,:,ie),w_greater(:,:,:,ie),w_retarded(:,:,:,ie))
    enddo
    call write_spectrum('gw_WR',iter,W_retarded,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_WL',iter,W_lesser  ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_WG',iter,W_greater ,nen,En-en(nen/2),nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    !
    print *, 'calc SigGW'
    nopmax=nen/2-10
    Sigma_greater_new = dcmplx(0.0d0,0.0d0)
    Sigma_lesser_new = dcmplx(0.0d0,0.0d0)
    Sigma_r_new = dcmplx(0.0d0,0.0d0)
    dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)    
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)    
    !$omp parallel default(none) private(io,ix,l,h,iop,nop,ie,i) shared(ndiag,nopmax,w_lesser,w_greater,w_retarded,sigma_lesser_new,sigma_greater_new,sigma_r_new,nen,En,nm,G_lesser,G_greater,G_r,nx)    
    !$omp do
    do ix=1,nx
      do i=1,nm
        l=max(i-ndiag,1)
        h=min(nm,i+ndiag)  
        do ie=1,nen      
          do nop= -nopmax,nopmax    
            if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then  
              iop=nop+nen/2                
              Sigma_lesser_new(i,l:h,ix,ie)=Sigma_lesser_new(i,l:h,ix,ie)+G_lesser(i,l:h,ix,ie-nop)*W_lesser(i,l:h,ix,iop)
              Sigma_greater_new(i,l:h,ix,ie)=Sigma_greater_new(i,l:h,ix,ie)+G_greater(i,l:h,ix,ie-nop)*W_greater(i,l:h,ix,iop)
              Sigma_r_new(i,l:h,ix,ie)=Sigma_r_new(i,l:h,ix,ie)+G_lesser(i,l:h,ix,ie-nop)*W_retarded(i,l:h,ix,iop) &
                                                            &  +G_r(i,l:h,ix,ie-nop)*W_lesser(i,l:h,ix,iop) &
                                                            &  +G_r(i,l:h,ix,ie-nop)*W_retarded(i,l:h,ix,iop)    
            endif
          enddo
        enddo
      enddo    
    enddo
    !$omp end do
    !$omp end parallel
    Sigma_lesser_new = Sigma_lesser_new  * dE
    Sigma_greater_new= Sigma_greater_new * dE
    Sigma_r_new=Sigma_r_new* dE
    Sigma_r_new = dcmplx( dble(Sigma_r_new), aimag(Sigma_greater_new-Sigma_lesser_new)/2.0d0 )    
    ! mixing with the previous one
    Sigma_r_gw= Sigma_r_gw+ alpha_mix * (Sigma_r_new -Sigma_r_gw)
    Sigma_lesser_gw  = Sigma_lesser_gw+ alpha_mix * (Sigma_lesser_new -Sigma_lesser_gw)
    Sigma_greater_gw = Sigma_greater_gw+ alpha_mix * (Sigma_greater_new -Sigma_greater_gw)  
    ! make sure self-energy is continuous near leads (by copying edge block)
    do ix=1,2
      Sigma_r_gw(:,:,ix,:)=Sigma_r_gw(:,:,3,:)
      Sigma_lesser_gw(:,:,ix,:)=Sigma_lesser_gw(:,:,3,:)
      Sigma_greater_gw(:,:,ix,:)=Sigma_greater_gw(:,:,3,:)
    enddo
    do ix=1,2
      Sigma_r_gw(:,:,nx-ix+1,:)=Sigma_r_gw(:,:,nx-2,:)
      Sigma_lesser_gw(:,:,nx-ix+1,:)=Sigma_lesser_gw(:,:,nx-2,:)
      Sigma_greater_gw(:,:,nx-ix+1,:)=Sigma_greater_gw(:,:,nx-2,:)
    enddo
    call write_spectrum('gw_SigR',iter,Sigma_r_gw,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_SigL',iter,Sigma_lesser_gw,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('gw_SigG',iter,Sigma_greater_gw,nen,En,nx,NB,NS,Lx,(/1.0d0,1.0d0/))
  enddo
  deallocate(g_r,g_lesser,g_greater,cur)
  deallocate(sigma_lesser_gw,sigma_greater_gw,sigma_r_gw)   
  deallocate(sigma_lesser_new,sigma_greater_new,sigma_r_new)   
  deallocate(P_retarded,P_lesser,P_greater)
  !deallocate(P_retarded_1i,P_lesser_1i,P_greater_1i)
  deallocate(W_retarded,W_lesser,W_greater)
  !deallocate(W_retarded_i1,W_lesser_i1,W_greater_i1)  
  deallocate(tot_cur,tot_ecur,jdens)  
end subroutine green_rgf_solve_gw_1d



subroutine full2block(fullA,A,NB,length,offdiag)
complex(8),intent(in)::fullA(NB*length,NB*length)
complex(8),intent(out)::A(NB,NB,length)
integer,intent(in)::NB,length
integer,intent(in),optional::offdiag
integer::i,l,h
A=0.0d0
do i=1,length
  l=(i-1)*NB+1
  h=i*NB
  if ((present(offdiag)).and.((i+offdiag) <= length)) then ! n-th off-diagonal blocks
    A(:,:,i)=fullA(l:h,l+nb*offdiag:h+nb*offdiag)
  else ! diagonal blocks
    A(:,:,i)=fullA(l:h,l:h)
  endif
enddo
end subroutine full2block


subroutine block2full(A,fullA,NB,length,offdiag)
complex(8),intent(in)::A(NB,NB,length)
complex(8),intent(out)::fullA(NB*length,NB*length)
integer,intent(in)::NB,length
integer,intent(in),optional::offdiag
integer::i,l,h
do i=1,length
  l=(i-1)*NB+1
  h=i*NB  
  if (present(offdiag)) then
    if ((i+offdiag) <= length) then ! n-th off-diagonal blocks
      fullA(l:h,l+nb*offdiag:h+nb*offdiag) = transpose(conjg(A(:,:,i)))
      fullA(l+nb*offdiag:h+nb*offdiag,l:h) = A(:,:,i)
    endif
  else ! diagonal blocks
    fullA(l:h,l:h)=A(:,:,i)
  endif
enddo
end subroutine block2full


subroutine green_calc_w_full(NBC,nm,nx,Vii,V1i,PLii,PGii,PRii,WLii,WGii,WRii)
integer,intent(in)::nm,nx,NBC
!complex(8),intent(in),dimension(nm*nx,nm*nx) :: V
complex(8),intent(in),dimension(nm,nm,nx) :: PLii,PGii,PRii
complex(8),intent(in),dimension(nm,nm,nx) :: Vii, V1i
!complex(8),intent(in),dimension(nm,nm,nx) :: PL1i,PG1i,PR1i
complex(8),intent(inout),dimension(nm,nm,nx) :: WLii,WGii,WRii!,WLi1,WGi1,WRi1
! --------- local
complex(8),allocatable,dimension(:,:)::B,M,V
complex(8),allocatable,dimension(:,:)::PL,PR,PG, WL,WR,WG
integer::i,NT
NT=nm*nx! total size
allocate(V(NT,NT))
allocate(PR(NT,NT))
allocate(PL(NT,NT))
allocate(PG(NT,NT))
allocate(WR(NT,NT))
allocate(WL(NT,NT))
allocate(WG(NT,NT))
allocate(B(NT,NT))
allocate(M(NT,NT))  
PR=czero
PL=czero
PG=czero
V=czero
call block2full(Vii,V,nm,nx)  
call block2full(V1i,V,nm,nx,1)  
call block2full(PRii,PR,nm,nx)  
call block2full(PLii,PL,nm,nx)  
call block2full(PGii,PG,nm,nx)  
!! M = I - V P^r
call zgemm('n','n',NT,NT,NT,-cone,V,NT,PR,NT,czero,M,NT)
!
do i=1,NT
   M(i,i) = 1.0d0 + M(i,i)
enddo    
!!!! calculate W^r = (I - V P^r)^-1 V    
call invert(M,NT) ! M -> xR
call zgemm('n','n',NT,NT,NT,cone,M,NT,V,NT,czero,WR,NT)           
! calculate W^< and W^> = W^r P^<> W^r dagger
call zgemm('n','n',NT,NT,NT,cone,WR,NT,PL,NT,czero,B,NT) 
call zgemm('n','c',NT,NT,NT,cone,B,NT,WR,NT,czero,WL,NT) 
call zgemm('n','n',NT,NT,NT,cone,WR,NT,PG,NT,czero,B,NT) 
call zgemm('n','c',NT,NT,NT,cone,B,NT,WR,NT,czero,WG,NT)  
call full2block(WR,WRii,nm,nx)
call full2block(WL,WLii,nm,nx)
call full2block(WG,WGii,nm,nx)
deallocate(B,M,V)  
deallocate(PR,PL,PG,WR,WL,WG)
end subroutine green_calc_w_full


subroutine green_rgf_calc_w(NBC,nm,nx,Vii,V1i,PL,PG,PR,WL,WG,WR,PL1i,PG1i,PR1i,WLi1,WGi1,WRi1)
implicit none
integer,intent(in)::nm,nx,NBC
complex(8),intent(in),dimension(nm,nm,nx) :: Vii,V1i,PL,PG,PR ! diagonal blocks of V and P matrices
complex(8),intent(in),dimension(nm,nm,nx),optional :: PL1i,PG1i,PR1i ! 1st off-diagonal blocks of P
complex(8),intent(inout),dimension(nm,nm,nx) :: WL,WG,WR ! diagonal blocks of W 
complex(8),intent(inout),dimension(nm,nm,nx),optional :: WLi1,WGi1,WRi1 ! 1st off-diagonal blocks of W
! --------
complex(8)::A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),sig(nm,nm),AL(nm,nm),AG(nm,nm),BL(nm,nm),BG(nm,nm)
complex(8),allocatable,dimension(:,:,:)::S,M,LL,LL1i,LG,LG1i,VV,M1i,Mi1,S1i,Si1,WRi1_,WLi1_,WGi1_,PR1i_,PL1i_,PG1i_
complex(8),allocatable,dimension(:,:,:)::xlr,wlr,wln,wlp,xR
complex(8),dimension(:,:),allocatable::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,&
    PL00,PL01,PL10,PG00,PG01,PG10,LL00,LL01,LL10,LG00,LG01,LG10
complex(8),dimension(:,:),allocatable::VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,&
    MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N
complex(8),dimension(:,:),allocatable::dM11,xR11,dLL11,dLG11,dV11
complex(8),dimension(:,:),allocatable::dMnn,xRnn,dLLnn,dLGnn,dVnn
integer::i,NL,NR,NT,LBsize,RBsize
integer::ix
real(8)::condL,condR
NL=nm ! left contact block size
NR=nm ! right contact block size
NT=nm*nx! total size
!print *,' RGF W'
allocate(VV(nm,nm,nx))
VV = Vii
allocate(M(nm,nm,nx))
allocate(M1i(nm,nm,nx))
allocate(Mi1(nm,nm,nx))
allocate(S(nm,nm,nx))
allocate(S1i(nm,nm,nx))
allocate(Si1(nm,nm,nx))
allocate(LL(nm,nm,nx))
allocate(LG(nm,nm,nx))
allocate(LL1i(nm,nm,nx))
allocate(LG1i(nm,nm,nx))
allocate(PR1i_(nm,nm,nx))
allocate(PL1i_(nm,nm,nx))
allocate(PG1i_(nm,nm,nx))
if (present(PR1i)) then
  PR1i_=PR1i
  PL1i_=PL1i
  PG1i_=PG1i
else
  PR1i_=czero
  PL1i_=czero
  PG1i_=czero
endif
!! S = V P^r
! Si,i = Vi,i Pi,i + Vi,i+1 Pi+1,i + Vi,i-1 Pi-1,i
do ix=1,nx
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,ix),nm,PR(:,:,ix),nm,czero,S(:,:,ix),nm)
  if (present(PR1i)) then
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,ix),nm,PR1i(:,:,ix),nm,cone,S(:,:,ix),nm) 
    if (ix==1) then
      call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,ix),nm,PR1i(:,:,ix),nm,cone,S(:,:,ix),nm) 
    else
      call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,ix-1),nm,PR1i(:,:,ix-1),nm,cone,S(:,:,ix),nm) 
    endif
  endif
enddo
!
do i=1,nx-1
  ! Si+1,i = Vi+1,i+1 Pi+1,i + Vi+1,i Pi,i
  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,PR(:,:,i),nm,czero,S1i(:,:,i),nm) 
  if (present(PR1i)) then
     call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,PR1i(:,:,i),nm,cone,S1i(:,:,i),nm) 
  endif
  !!
  ! Si,i+1 = Vi,i Pi,i+1 + Vi,i+1 Pi+1,i+1
  call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PR(:,:,i+1),nm,czero,Si1(:,:,i),nm)   
  if (present(PR1i)) then
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PR1i(:,:,i),nm,cone,Si1(:,:,i),nm)
  endif
  !!  
enddo
Si1(:,:,nx)=Si1(:,:,nx-1)
S1i(:,:,nx)=S1i(:,:,nx-1)
M   = -S
M1i = -S1i
Mi1 = -Si1
do ix=1,nx
  do i=1,nm
      M(i,i,ix) = 1.0d0 + M(i,i,ix)
  enddo
enddo
deallocate(S,S1i,Si1)
!! LL=V P^l V
!! LLi,i
do i=1,nx
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PL(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LL(:,:,i),nm)
  !
  if (i<nx) then
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL(:,:,i+1),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL(:,:,i),nm)
  endif
  if (i>1) then
    call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i-1),nm,PL(:,:,i-1),nm,czero,A,nm)
    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,i-1),nm,cone,LL(:,:,i),nm)
  endif
  !
!  if (present(PL1i)) then
!    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL1i(:,:,i),nm,czero,A,nm)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LL(:,:,i),nm)
!    !
!    call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,PL1i(:,:,i),nm,czero,A,nm)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL(:,:,i),nm)
!    !
!    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PL1i(:,:,max(i-1,1)),nm,czero,A,nm)
!    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LL(:,:,i),nm)      
!    !
!    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PL1i(:,:,max(i-1,1)),nm,czero,A,nm)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LL(:,:,i),nm)
!  endif
enddo
!! LLi+1,i = Vi+1,i P^li,i Vi,i + Vi+1,i+1 P^li+1,i+1 Vi+1,i
do i=1,nx
  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LL1i(:,:,i),nm)
  if (i<nx) then
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,PL(:,:,i+1),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL1i(:,:,i),nm)
  endif
enddo
!! LG=V P^g V'    
!! LGi,i
do i=1,nx
  call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PG(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LG(:,:,i),nm)
  !
  if (i<nx) then
    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG(:,:,i+1),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG(:,:,i),nm)
  endif
  if (i>1) then    
    call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i-1),nm,PG(:,:,i-1),nm,czero,A,nm)
    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,i-1),nm,cone,LG(:,:,i),nm)
  endif
  !
!  if (present(PG1i)) then
!    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG1i(:,:,i),nm,czero,A,nm)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LG(:,:,i),nm)
!    !
!    call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,PG1i(:,:,i),nm,czero,A,nm)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG(:,:,i),nm)
!    !
!    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PG1i(:,:,max(i-1,1)),nm,czero,A,nm)
!    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LG(:,:,i),nm)
!    !
!    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PG1i(:,:,max(i-1,1)),nm,czero,A,nm)
!    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LG(:,:,i),nm)
!  endif
enddo
!! LGi+1,i = Vi+1,i P^gi,i Vi,i + Vi+1,i+1 P^gi+1,i+1 Vi+1,i 
do i=1,nx
  call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG(:,:,i),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LG1i(:,:,i),nm)
  if (i<nx) then
    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,PG(:,:,i+1),nm,czero,A,nm)
    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG1i(:,:,i),nm)
  endif
enddo
!
if (NBC>0) then ! apply open-boundary correction to all matrices
  !
  LBsize=NL*NBC
  RBsize=NR*NBC
  !
  allocate(V00 (LBsize,LBsize))
  allocate(V01 (LBsize,LBsize))
  allocate(V10 (LBsize,LBsize))
  allocate(M00 (LBsize,LBsize))
  allocate(M01 (LBsize,LBsize))
  allocate(M10 (LBsize,LBsize))
  allocate(PR00(LBsize,LBsize))
  allocate(PR01(LBsize,LBsize))
  allocate(PR10(LBsize,LBsize))
  allocate(PG00(LBsize,LBsize))
  allocate(PG01(LBsize,LBsize))
  allocate(PG10(LBsize,LBsize))
  allocate(PL00(LBsize,LBsize))
  allocate(PL01(LBsize,LBsize))
  allocate(PL10(LBsize,LBsize))
  allocate(LG00(LBsize,LBsize))
  allocate(LG01(LBsize,LBsize))
  allocate(LG10(LBsize,LBsize))
  allocate(LL00(LBsize,LBsize))
  allocate(LL01(LBsize,LBsize))
  allocate(LL10(LBsize,LBsize))
  allocate(dM11(LBsize,LBsize))
  allocate(xR11(LBsize,LBsize))
  allocate(dV11(LBsize,LBsize))
  allocate(dLL11(LBsize,LBsize))
  allocate(dLG11(LBsize,LBsize))
  !  
  allocate(VNN  (RBsize,RBsize))
  allocate(VNN1 (RBsize,RBsize))
  allocate(Vn1n (RBsize,RBsize))
  allocate(Mnn  (RBsize,RBsize))
  allocate(Mnn1 (RBsize,RBsize))
  allocate(Mn1n (RBsize,RBsize))
  allocate(PRnn (RBsize,RBsize))
  allocate(PRnn1(RBsize,RBsize))
  allocate(PRn1n(RBsize,RBsize))
  allocate(PGnn (RBsize,RBsize))
  allocate(PGnn1(RBsize,RBsize))
  allocate(PGn1n(RBsize,RBsize))
  allocate(PLnn (RBsize,RBsize))
  allocate(PLnn1(RBsize,RBsize))
  allocate(PLn1n(RBsize,RBsize))
  allocate(LGnn (RBsize,RBsize))
  allocate(LGnn1(RBsize,RBsize))
  allocate(LGn1n(RBsize,RBsize))
  allocate(LLnn (RBsize,RBsize))
  allocate(LLnn1(RBsize,RBsize))
  allocate(LLn1n(RBsize,RBsize))
  allocate(dMnn (RBsize,RBsize))
  allocate(xRnn (RBsize,RBsize))
  allocate(dLLnn(RBsize,RBsize))
  allocate(dLGnn(RBsize,RBsize))
  allocate(dVnn(RBsize,RBsize))
  ! left OBC
  call get_OBC_blocks_for_W(NL,Vii(:,:,1),V1i(:,:,1),PR(:,:,1),PR1i_(:,:,1),&
      PL(:,:,1),PL1i_(:,:,1),PG(:,:,1),PG1i_(:,:,1),NBC,&
      V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
      LL00,LL01,LL10,LG00,LG01,LG10)
  ! right OBC  
  call get_OBC_blocks_for_W(NR,Vii(:,:,nx),transpose(conjg(V1i(:,:,nx))),PR(:,:,nx),&
      transpose(PR1i_(:,:,nx)),PL(:,:,nx),-transpose(conjg(PL1i_(:,:,nx))),&
      PG(:,:,nx),-transpose(conjg(PG1i_(:,:,nx))),NBC,&
      VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,&
      LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N)
  !
  ! Correct first and last block of M to account for elements in the contacts
  M(:,:,1)=M(:,:,1) - matmul(V10,PR01)
  M(:,:,nx)=M(:,:,nx) - matmul(VNN1,PRN1N)
  ! 
  ! Correct first and last block of LL to account for elements in the contacts
  LL(:,:,1)=LL(:,:,1) + matmul(matmul(V10,PL00),V01) + &
    matmul(matmul(V10,PL01),V00) + matmul(matmul(V00,PL10),V01)
  !  
  LL(:,:,nx)=LL(:,:,nx) + &
    matmul(matmul(VNN,PLNN1),VN1N) + matmul(matmul(VNN1,PLN1N),VNN) + &
    matmul(matmul(VNN1,PLNN),VN1N)
  !  
  ! Correct first and last block of LG to account for elements in the contacts
  LG(:,:,1)=LG(:,:,1) + matmul(matmul(V10,PG00),V01) + &
    matmul(matmul(V10,PG01),V00) + matmul(matmul(V00,PG10),V01)
  LG(:,:,nx)=LG(:,:,nx) + &
    matmul(matmul(VNN,PGNN1),VN1N) + matmul(matmul(VNN1,PGN1N),VNN) + matmul(matmul(VNN1,PGNN),VN1N)
  !  
  ! WR/WL/WG OBC Left
  call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
  ! WR/WL/WG OBC right
  call open_boundary_conditions(NR,MNN,MNN1,MN1N,VN1N,xRNN,dMNN,dVNN,condR)
  !
  if (condL<1.0d-6) then   
    !
    call get_dL_OBC_for_W(NL,xR11,LL00,LL01,LG00,LG01,M10,'L', dLL11,dLG11)
    !
    M(:,:,1)=M(:,:,1) - dM11
    VV(:,:,1)=VV(:,:,1) - dV11    
    LL(:,:,1)=LL(:,:,1) + dLL11
    LG(:,:,1)=LG(:,:,1) + dLG11    
  endif
  if (condR<1.0d-6) then    
    !
    call get_dL_OBC_for_W(NR,xRNN,LLNN,LLN1N,LGNN,LGN1N,MNN1,'R', dLLNN,dLGNN)
    !
    M(:,:,nx)=M(:,:,nx) - dMNN
    VV(:,:,nx)=VV(:,:,nx)- dVNN
    LL(:,:,nx)=LL(:,:,nx) + dLLNN
    LG(:,:,nx)=LG(:,:,nx) + dLGNN    
  endif
  !
  deallocate(V00,V01,V10)
  deallocate(M00,M01,M10)
  deallocate(PR00,PR01,PR10)
  deallocate(PG00,PG01,PG10)
  deallocate(PL00,PL01,PL10)
  deallocate(LG00,LG01,LG10)
  deallocate(LL00,LL01,LL10)
  deallocate(VNN,VNN1,Vn1n)
  deallocate(Mnn,Mnn1,Mn1n)
  deallocate(PRnn,PRnn1,PRn1n)
  deallocate(PGnn,PGnn1,PGn1n)
  deallocate(PLnn,PLnn1,PLn1n)
  deallocate(LGnn,LGnn1,LGn1n)
  deallocate(LLnn,LLnn1,LLn1n)
  deallocate(dM11,xR11,dLL11,dLG11,dV11)
  deallocate(dMnn,xRnn,dLLnn,dLGnn,dVnn)
  !
end if 
!
allocate(xlr(nm,nm,nx)) ! left-connected xR
allocate(wlr(nm,nm,nx)) ! left-connected Wr
allocate(wln(nm,nm,nx)) ! left-connected W<
allocate(wlp(nm,nm,nx)) ! left-connected W>
allocate(xR(nm,nm,nx))  ! full xR
allocate(WRi1_(nm,nm,nx)) ! off-diagonal block of WR
!print*,' first pass, from right to left ... '
A=M(:,:,nx)
call invert(A,nm)
xlr(:,:,nx)=A
call zgemm('n','n',nm,nm,nm,cone,A,nm,VV(:,:,nx),nm,czero,Wlr(:,:,nx),nm) 
!
call zgemm('n','n',nm,nm,nm,cone,A,nm,LL(:,:,nx),nm,czero,B,nm) 
call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wln(:,:,nx),nm) 
!
call zgemm('n','n',nm,nm,nm,cone,A,nm,LG(:,:,nx),nm,czero,B,nm) 
call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wlp(:,:,nx),nm) 
!
do ix=nx-1,1,-1
  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,xlr(:,:,ix+1),nm,czero,B,nm) ! B -> Mn,n+1 xlr -> MxR
  call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,C,nm)             ! C -> Mn,n+1 xlr Mn+1,n
  A = M(:,:,ix) - C
  call invert(A,nm) ! A -> xlr
  xlr(:,:,ix)=A     
  call zgemm('n','n',nm,nm,nm,cone,B,nm,V1i(:,:,ix),nm,czero,C,nm) ! C -> MxR Vn+1,n  
  D=VV(:,:,ix) - C                                                 ! D -> Vn,n - Mn,n+1 xlr Vn+1,n
  call zgemm('n','n',nm,nm,nm,cone,A,nm,D,nm,czero,Wlr(:,:,ix),nm) ! xlr ( Vn,n - Mn,n+1 xlr Vn+1,n )
  !
!  !AL=MxR*LL10  
!  !AG=MxR*LG10          
!  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,xlr(:,:,ix+1),nm,czero,B,nm) ! B -> Mn,n+1 xlr -> MxR
!  AL=matmul(B, LL1i(:,:,ix))
!  AG=matmul(B, LG1i(:,:,ix))
!  !  
!  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,wln(:,:,ix+1),nm,czero,B,nm) 
!  call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,sig,nm)   ! sig -> Mn,n+1 w<n+1 Mn+1,n
!  C = LL(:,:,ix) + sig - (AL-transpose(conjg(AL)))
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
!  call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,wln(:,:,ix),nm)       
!  !
!  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,wlp(:,:,ix+1),nm,czero,B,nm) 
!  call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,sig,nm)   ! sig -> Mn,n+1 w>n+1 Mn+1,n    
!  C = LG(:,:,ix) + sig - (AG-transpose(conjg(AG)))
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
!  call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,wlp(:,:,ix),nm)     
  !  
enddo
!print*,' second pass, from left to right ... '
WR(:,:,1)=Wlr(:,:,1)
call zgemm('n','t',nm,nm,nm,cone,WR(:,:,1),nm,M1i(:,:,1),nm,czero,B,nm)
C = transpose((V1i(:,:,1))) - B
call zgemm('n','t',nm,nm,nm,cone,C,nm,xlr(:,:,2),nm,czero,B,nm)  ! B -> (V12 - WR1 M12) xlr2
WRi1_(:,:,1) = B
XR(:,:,1)=xlr(:,:,1)
!WL(:,:,1)=Wln(:,:,1)
!WG(:,:,1)=Wlp(:,:,1)
do ix=2,nx
  call zgemm('n','n',nm,nm,nm,cone,xlr(:,:,ix),nm,M1i(:,:,ix-1),nm,czero,B,nm) ! B -> xlr1 * M10
  call zgemm('n','n',nm,nm,nm,cone,B,nm,WRi1_(:,:,ix-1),nm,czero,C,nm)         ! C -> xlr1 * M10 WR0 M01 xlr1
  WR(:,:,ix)=Wlr(:,:,ix)  - C  
  call zgemm('n','n',nm,nm,nm,cone,B,nm,XR(:,:,ix-1),nm,czero,A,nm)            ! A -> xlr1 * M10 * XR0 
  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix-1),nm,xlr(:,:,ix),nm,czero,C,nm) ! C -> M01 * xlr1
  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)                       ! D -> xlr1 * M10 * XR0 * M01 * xlr1
  XR(:,:,ix)=xlr(:,:,ix)  + D
  !   
!  ! AL=xlr1 LL10 XR0' M01 xlr1'
!  ! BL=xlr1 M10 XR0 M01 wln1  
!  call zgemm('c','n',nm,nm,nm,cone,xR(:,:,ix-1),nm,Mi1(:,:,ix-1),nm,czero,D,nm)   ! D -> XR0' M01
!  call zgemm('n','c',nm,nm,nm,cone,D,nm,xlr(:,:,ix),nm,czero,A,nm)    ! A -> XR0' M01 xlr1'
!  !
!  AL = matmul( matmul(xlr(:,:,ix), LL1i(:,:,ix-1)), A )
!  AG = matmul( matmul(xlr(:,:,ix), LG1i(:,:,ix-1)), A )
!  !
!  call zgemm('n','n',nm,nm,nm,cone,B,nm,xR(:,:,ix-1),nm,czero,D,nm)   ! D -> xlr1 * M10 * XR0
!  call zgemm('n','n',nm,nm,nm,cone,D,nm,Mi1(:,:,ix-1),nm,czero,A,nm)  ! A -> xlr1 * M10 * XR0 * M01 
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,wln(:,:,ix),nm,czero,BL,nm)   ! BL-> xlr1 * M10 * XR0 * M01 * wln1
!  !
!  call zgemm('n','n',nm,nm,nm,cone,B,nm,WL(:,:,ix-1),nm,czero,A,nm)            ! A -> xlr1 M10 WL0
!  call zgemm('n','c',nm,nm,nm,cone,Mi1(:,:,ix-1),nm,xlr(:,:,ix),nm,czero,C,nm) ! C -> M01 xlr1'
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)            ! D -> xlr1 M10 WL0 M01 xlr1'
!  WL(:,:,ix) = wln(:,:,ix) !+ D + (BL-transpose(conjg(BL))) - (AL-transpose(conjg(AL)))
!  !
!  ! AG=xlr1 LG10 XR0' M01 xlr1'
!  ! BG=xlr1 M10 XR0 M01 wlp1    
!  !
!  call zgemm('n','n',nm,nm,nm,cone,B,nm,xR(:,:,ix-1),nm,czero,D,nm)    ! D -> xlr1 * M10 * XR0 
!  call zgemm('n','n',nm,nm,nm,cone,D,nm,Mi1(:,:,ix-1),nm,czero,A,nm)   ! A -> xlr1 * M10 * XR0 * M01
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,wlp(:,:,ix),nm,czero,BG,nm)    ! BG-> xlr1 * M10 * XR0 * M01 * wlp1
!  !
!  call zgemm('n','n',nm,nm,nm,cone,B,nm,WG(:,:,ix-1),nm,czero,A,nm)    ! A -> xlr1 * M10 * WG0 
!  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)               ! D -> xlr1 * M10 * WG0 * M01 xlr1' 
!  WG(:,:,ix) = wlp(:,:,ix) !+ D + (BG-transpose(conjg(BG))) - (AG-transpose(conjg(AG)))
  !
  if (ix<nx) then
    call zgemm('n','t',nm,nm,nm,cone,WR(:,:,ix),nm,M1i(:,:,ix),nm,czero,B,nm)
    C = transpose((V1i(:,:,ix))) - B
    call zgemm('n','t',nm,nm,nm,cone,C,nm,xlr(:,:,ix+1),nm,czero,B,nm)  ! B -> (V12 - WR1 M12) xlr2
    WRi1_(:,:,ix) = B        
  endif
  !     
enddo
!
if (present(WRi1)) then
  WRi1= WRi1_
end if
!
WRi1_ = dcmplx(0.0d0*dble(WRi1_),aimag(WRi1_))
do ix=1,nx
  call zgemm('n','n',nm,nm,nm,cone,WR(:,:,ix),nm,PL(:,:,ix),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,WR(:,:,ix),nm,czero,B,nm)
  if (ix<nx) then
    call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix),nm,PL(:,:,ix+1),nm,czero,A,nm)
    call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix),nm,cone,B,nm)
  else
    call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix-1),nm,PL(:,:,nx),nm,czero,A,nm)
    call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix-1),nm,cone,B,nm)  
  endif
  call zgemm('t','n',nm,nm,nm,cone,WRi1_(:,:,max(ix-1,1)),nm,PL(:,:,max(ix-1,1)),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,-WRi1_(:,:,max(ix-1,1)),nm,cone,B,nm)
  WL(:,:,ix)=B
  call zgemm('n','n',nm,nm,nm,cone,WR(:,:,ix),nm,PG(:,:,ix),nm,czero,A,nm)
  call zgemm('n','c',nm,nm,nm,cone,A,nm,WR(:,:,ix),nm,czero,B,nm)
  if (ix<nx) then
    call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix),nm,PG(:,:,ix+1),nm,czero,A,nm)
    call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix),nm,cone,B,nm)
  else
    call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix-1),nm,PG(:,:,nx),nm,czero,A,nm)
    call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix-1),nm,cone,B,nm)  
  endif
  call zgemm('t','n',nm,nm,nm,cone,WRi1_(:,:,max(ix-1,1)),nm,PG(:,:,max(ix-1,1)),nm,czero,A,nm)
  call zgemm('n','n',nm,nm,nm,cone,A,nm,-WRi1_(:,:,max(ix-1,1)),nm,cone,B,nm)
  WG(:,:,ix)=B
enddo
!
deallocate(M,M1i,Mi1,LL,LG,VV,LL1i,LG1i)
deallocate(wln,wlp,wlr,xlr,Xr,WRi1_)
deallocate(PR1i_,PL1i_,PG1i_)
end subroutine green_rgf_calc_w




!subroutine green_calc_w_banded(NBC,nm,nx,Vii,V1i,PLii,PGii,PRii,WLii,WGii,WRii)
!integer,intent(in)::nm,nx,NBC
!!complex(8),intent(in),dimension(nm*nx,nm*nx) :: V
!complex(8),intent(in),dimension(nm,nm,nx) :: PLii,PGii,PRii
!complex(8),intent(in),dimension(nm,nm,nx) :: Vii, V1i
!!complex(8),intent(in),dimension(nm,nm,nx) :: PL1i,PG1i,PR1i
!complex(8),intent(inout),dimension(nm,nm,nx) :: WLii,WGii,WRii!,WLi1,WGi1,WRi1
!! --------- local
!complex(8),allocatable,dimension(:,:)::B,M,V, A
!complex(8),allocatable,dimension(:,:)::PL,PR,PG, WL,WR,WG
!integer,allocatable::IPIV(:)
!integer::i,NT, LDAB, info, LDB,NRHS
!NT=nm*nx! total size
!LDAB= 2*nm*2+nm*2+1
!LDB=NT
!NRHS=NT
!allocate(IPIV(NT))
!allocate(V(NT,NT))
!allocate(PR(NT,NT))
!allocate(PL(NT,NT))
!allocate(PG(NT,NT))
!allocate(WR(NT,NT))
!allocate(WL(NT,NT))
!allocate(WG(NT,NT))
!allocate(B(NT,NT))
!allocate(M(NT,NT))  
!PR=czero
!PL=czero
!PG=czero
!V=czero
!call block2full(Vii,V,nm,nx)  
!call block2full(V1i,V,nm,nx,1)  
!call block2full(PRii,PR,nm,nx)  
!call block2full(PLii,PL,nm,nx)  
!call block2full(PGii,PG,nm,nx)  
!!! M = I - V P^r
!call zgemm('n','n',NT,NT,NT,-cone,V,NT,PR,NT,czero,M,NT)
!!
!do i=1,NT
!   M(i,i) = 1.0d0 + M(i,i)
!enddo  
!! put M into band storage 
!allocate(A(LDAB,NT))
!do i=1,NT
!  A(max(i-nm*2,1)-i+nm*4:min(i+nm*2,NT)-i+nm*4,i) = M(max(i-nm*2,1):min(i+nm*2,NT),i)
!enddo
!M=czero
!do i=1,NT
!   M(i,i) = 1.0d0 
!enddo    
!call ZGBSV( NT, nm*2, nm*2, NT, A, LDAB, IPIV, M, LDB, info )
!if (info .ne. 0) then
!  print*,'SEVERE warning: ZGBSV failed, info=',info
!  M=0.0d0
!endif  
!deallocate(A,IPIV)
!! M contains (I - V P^r)^-1
!!!!! calculate W^r = (I - V P^r)^-1 V      
!call zgemm('n','n',NT,NT,NT,cone,M,NT,V,NT,czero,WR,NT)           
!! calculate W^< and W^> = W^r P^<> W^r dagger
!call zgemm('n','n',NT,NT,NT,cone,WR,NT,PL,NT,czero,B,NT) 
!call zgemm('n','c',NT,NT,NT,cone,B,NT,WR,NT,czero,WL,NT) 
!call zgemm('n','n',NT,NT,NT,cone,WR,NT,PG,NT,czero,B,NT) 
!call zgemm('n','c',NT,NT,NT,cone,B,NT,WR,NT,czero,WG,NT)  
!call full2block(WR,WRii,nm,nx)
!call full2block(WL,WLii,nm,nx)
!call full2block(WG,WGii,nm,nx)
!deallocate(B,M,V)  
!deallocate(PR,PL,PG,WR,WL,WG)
!end subroutine green_calc_w_banded


subroutine identity(A,n)
  integer, intent(in) :: n        
  complex(8), dimension(n,n), intent(inout) :: A
  integer :: i
  A = dcmplx(0.0d0,0.0d0)
  do i = 1,n
    A(i,i) = dcmplx(1.0d0,0.0d0)
  end do
end subroutine identity



!!!! RGF for diagonal blocks of G^r,<,>
subroutine green_RGF_RS(TEMP,nm,nx,E,mu,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,GRi1,GLi1,GGi1)     
    integer,intent(in)::nm,nx
    real(8),intent(in) :: E,mu(2),temp(2)
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph ! diag blocks of scattering SE
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Hii,H1i ! diag blocks of Overlap and H and 1st off-diag blocks of H
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: ldos,ndens,pdens ! diag blocks of GFs    
    real(8),intent(inout),dimension(nm,nm,Nx) :: cur ! current density
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx),optional :: GRi1,GLi1,GGi1 ! off-diag blocks (i,i+1) of GFs
    real(8),intent(inout) :: tr,tre ! current spectrum on the Right and Left contacts
    ! ------- 
    COMPLEX(8) :: H00(nm,nm),H10(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),S00(nm,nm),G00(nm,nm),GBB(nm,nm),GN0(nm,nm),Gn(nm,nm),Gp(nm,nm)
    COMPLEX(8) :: sig(nm,nm),sigmal(nm,nm),sigmar(nm,nm)
    COMPLEX(8) :: z
    integer::i,j,k,l
    real(8)::tim,mul,mur,templ,tempr,tmp1
    COMPLEX(8), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:) ! left-connected green function
    complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
    complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
    REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
    mul=mu(1)
    mur=mu(2)
    templ=temp(1)
    tempr=temp(2)
    !
    z=dcmplx(E,0.0d-6)
    !
    allocate(Gl(nm,nm,nx))
    allocate(Gln(nm,nm,nx))
    allocate(Glp(nm,nm,nx))
    !
    Gln=0.0D0
    Glp=0.0D0
    Gl=0.0D0
    ldos=0.0d0
    ndens=0.0d0
    pdens=0.0d0
    cur=0.0d0
    S00=0.0d0
    do i=1,nm
        S00(i,i)=1.0d0
    enddo
    ! self energy on the left contact        
    H00(:,:)=Hii(:,:,1)+sigma_r_ph(:,:,1)
    H10(:,:)=H1i(:,:,1)
    !
    call sancho(NM,E,S00,H00,transpose(conjg(H10)),G00,GBB)
    !
    tmp1=maxval(abs(transpose(G00) - G00)) / maxval(abs(G00))
    if (tmp1 > 2e-2) then
        print *, 'E=',E, '|G00^t-G00|=', tmp1
    endif
    !
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmal,nm)      
    sig(:,:)=-(sigmal(:,:)-transpose(conjg(sigmal(:,:))))*ferm((E-mul)/(BOLTZ*TEMPl))+sigma_lesser_ph(:,:,1)
    A=z*S00-H00-sigmal
    !                
    call invert(A,nm)
    Gl(:,:,1)=A(:,:)
    !
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
    Gln(:,:,1)=C(:,:)
    Do l=2,nx-1                
        H00(:,:)=Hii(:,:,l)+sigma_r_ph(:,:,l)
        H10(:,:)=H1i(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,l-1),nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
        A=z*S00-H00-C   
        !
        call invert(A,nm)
        Gl(:,:,l)=A(:,:)
        !
        sig=Gln(:,:,l-1)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)             
        C(:,:)=C(:,:)+sigma_lesser_ph(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
        Gln(:,:,l)=Gn(:,:)
    enddo
    ! self energy on the right contact        
    H00(:,:)=Hii(:,:,nx)+sigma_r_ph(:,:,nx)
    H10(:,:)=H1i(:,:,nx)
    !
    call sancho(NM,E,S00,H00,H10,G00,GBB)
    !
    call zgemm('c','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmar,nm)  
    H10(:,:)=H1i(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,nx-1),nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
    G00=z*S00-H00-sigmar-C   
    !
    call invert(G00,nm)
    !
    ldos(:,:,nx)=G00(:,:)!cmplx(0.0d0,1.0d0)*(G00(:,:)-transpose(conjg(G00(:,:))))
    sig=Gln(:,:,nx-1)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)  ! C=H10 Gl< H01    
    ! B=Sig< S00
    sig(:,:)=-(sigmar(:,:)-transpose(conjg(sigmar(:,:))))*ferm((E-mur)/(BOLTZ*TEMPl))+C(:,:)+sigma_lesser_ph(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
    ! G<00 = G00 sig< G00'
    ndens(:,:,nx)=Gn(:,:)    
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    pdens(:,:,nx)=Gp(:,:)
    A=-(sigmar-transpose(conjg(sigmar)))*ferm((E-mur)/(BOLTZ*TEMPl))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmar-transpose(conjg(sigmar)))*(ferm((E-mur)/(BOLTZ*TEMPl))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm        
      tim=tim-dble(B(i,i)-C(i,i))        
    enddo
    tr=tim
    ! transmission
    !-------------------------
    do l=nx-1,1,-1
        H10(:,:)=H1i(:,:,l)
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,H10,nm,Gl(:,:,l),nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,B,nm,beta,C,nm) 
        A=Gln(:,:,l)          
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,A,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,G00,nm,B,nm,beta,A,nm)
        B=C+A
        call zgemm('c','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i+1,i     
        cur(:,:,l)=dble(A(:,:))     
        !-------------------------
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G<
        A=Gln(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g< H10 G'
        B=C+A
        if (present(GLi1)) then 
          GLi1(:,:,l)=B
        endif
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i,i+1
        cur(:,:,l)=cur(:,:,l)-dble(A(:,:))        
        !-------------------------
        if (present(GGi1)) then
          A=Gp
          call zgemm('c','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G>
          A=Glp(:,:,l)
          call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g> H10 G
          B=C+A         
          GGi1(:,:,l)=B
        endif        
        !-------------------------
        D(:,:)= Gl(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
        if (present(GRi1)) then 
          GRi1(:,:,l)=GN0
        endif
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
        G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i
        ldos(:,:,l)=G00(:,:)
        !-------------------------
        A(:,:)=Gn(:,:)     
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        Gn(:,:)= Gln(:,:,l) + C(:,:)
        A(:,:)=Gln(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)         
        Gn(:,:)= Gn(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)          
        Gn(:,:)= Gn(:,:)+C(:,:)!     					 !!! G<_i,i
        !-------------------------
!        A(:,:)=Gp(:,:)
!        call zgemm('c','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
!        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
!        !
!        Gp(:,:)= Glp(:,:,l) + C(:,:)
!        A(:,:)=Glp(:,:,l)
!        call zgemm('c','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        Gp(:,:)= Gp(:,:)+C(:,:)!                     			 
!        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)     
!        Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G>_i,i
        !-------------------------    
        ndens(:,:,l)=Gn(:,:)
        pdens(:,:,l)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    enddo
    !
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    A=-(sigmal-transpose(conjg(sigmal)))*ferm((E-mul)/(BOLTZ*TEMPr))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmal-transpose(conjg(sigmal)))*(ferm((E-mul)/(BOLTZ*TEMPr))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
      tim=tim+dble(B(i,i)-C(i,i))
    enddo
    tre=tim
    deallocate(Gl)
    deallocate(Gln)
    deallocate(Glp)
end subroutine green_RGF_RS

!!!! Coupled-Mode-Space RGF for diagonal blocks of G^r,<,>
subroutine green_RGF_CMS(TEMP,nm,nx,E,mu,Sii,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,ndens,pdens,ldos,tr,tre,cur,GRi1,GLi1,GGi1)         
    integer,intent(in)::nm,nx
    real(8),intent(in) :: E,mu(2),temp(2)
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph ! diag blocks of scattering SE
    COMPLEX(8),intent(in),dimension(nm,nm,Nx) :: Sii,Hii,H1i ! diag blocks of Overlap and H and 1st off-diag blocks of H
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx) :: ldos,ndens,pdens ! diag blocks of GFs
    real(8),intent(inout),dimension(nm,nm,Nx) :: cur
    COMPLEX(8),intent(inout),dimension(nm,nm,Nx),optional :: GRi1,GLi1,GGi1 ! off-diag blocks (i,i+1) of GFs
    real(8),intent(inout) :: tr,tre ! current spectrum on the Right and Left contacts
    ! ------- 
    COMPLEX(8) :: H00(nm,nm),H10(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),S00(nm,nm),G00(nm,nm),GBB(nm,nm),GN0(nm,nm),Gn(nm,nm),Gp(nm,nm)
    COMPLEX(8) :: sig(nm,nm),sigmal(nm,nm),sigmar(nm,nm)
    COMPLEX(8) :: z
    integer::i,j,k,l
    real(8)::tim,mul,mur,templ,tempr
    COMPLEX(8), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:) ! left-connected green function
    complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
    complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
    REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
    mul=mu(1)
    mur=mu(2)
    templ=temp(1)
    tempr=temp(2)
    !
    z=dcmplx(E,0.0d-6)
    !
    allocate(Gl(nm,nm,nx))
    allocate(Gln(nm,nm,nx))
    allocate(Glp(nm,nm,nx))
    !
    Gln=0.0D0
    Glp=0.0D0
    Gl=0.0D0
    ldos=0.0d0
    ndens=0.0d0
    pdens=0.0d0
    cur=0.0d0
    ! self energy on the left contact
    S00(:,:)=Sii(:,:,1)
    call zgemm('n','n',nm,nm,nm,alpha,sigma_r_ph(:,:,1),nm,S00,nm,beta,B,nm)
    H00(:,:)=Hii(:,:,1)+B(:,:)!sigma_r_ph(:,:,1)
    H10(:,:)=H1i(:,:,1)
    !
    call sancho(NM,E,S00,H00,transpose(conjg(H10)),G00,GBB)
    !
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmal,nm)  
    call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser_ph(:,:,1),nm,S00,nm,beta,B,nm)
    sig(:,:)=-(sigmal(:,:)-transpose(conjg(sigmal(:,:))))*ferm((E-mul)/(BOLTZ*TEMPl))+B(:,:)
    A=z*S00-H00-sigmal
    !                
    call invert(A,nm)
    Gl(:,:,1)=A(:,:)
    !
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
    Gln(:,:,1)=C(:,:)
    Do l=2,nx-1
        S00(:,:)=Sii(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,sigma_r_ph(:,:,l),nm,S00,nm,beta,B,nm)
        H00(:,:)=Hii(:,:,l)+B(:,:)
        H10(:,:)=H1i(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,l-1),nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
        A=z*S00-H00-C   
        !
        call invert(A,nm)
        Gl(:,:,l)=A(:,:)
        !
        sig=Gln(:,:,l-1)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser_ph(:,:,l),nm,S00,nm,beta,B,nm)
        C(:,:)=C(:,:)+B(:,:)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
        Gln(:,:,l)=Gn(:,:)
    enddo
    ! self energy on the right contact
    S00(:,:)=Sii(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,sigma_r_ph(:,:,nx),nm,S00,nm,beta,B,nm)
    H00(:,:)=Hii(:,:,nx)+B(:,:)
    H10(:,:)=H1i(:,:,nx)
    !
    call sancho(NM,E,S00,H00,H10,G00,GBB)
    !
    call zgemm('c','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmar,nm)  
    H10(:,:)=H1i(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,nx-1),nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
    G00=z*S00-H00-sigmar-C   
    !
    call invert(G00,nm)
    !
    ldos(:,:,nx)=G00(:,:)!dcmplx(0.0d0,1.0d0)*(G00(:,:)-transpose(conjg(G00(:,:))))
    sig=Gln(:,:,nx-1)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)  ! C=H10 Gl< H01
    call zgemm('n','n',nm,nm,nm,alpha,sigma_lesser_ph(:,:,nx),nm,S00,nm,beta,B,nm)
    ! B=Sig< S00
    sig(:,:)=-(sigmar(:,:)-transpose(conjg(sigmar(:,:))))*ferm((E-mur)/(BOLTZ*TEMPl))+C(:,:)+B(:,:)
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
    ! G<00 = G00 sig< G00'
    ndens(:,:,nx)=Gn(:,:)    
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    pdens(:,:,nx)=Gp(:,:)
    A=-(sigmar-transpose(conjg(sigmar)))*ferm((E-mur)/(BOLTZ*TEMPl))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmar-transpose(conjg(sigmar)))*(ferm((E-mur)/(BOLTZ*TEMPl))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
        do j=i,i!1,nm
            tim=tim-dble(B(i,j)-C(i,j))
        enddo
    enddo
    tr=tim
    ! transmission
    !-------------------------
    do l=nx-1,1,-1
        H10(:,:)=H1i(:,:,l)
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,H10,nm,Gl(:,:,l),nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,B,nm,beta,C,nm) 
        A=Gln(:,:,l)          
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,A,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,G00,nm,B,nm,beta,A,nm)
        B=C+A
        call zgemm('c','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i+1,i     
        cur(:,:,l)=dble(A(:,:))     
        !-------------------------
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G<
        A=Gln(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g< H10 G'
        B=C+A
        if (present(GLi1)) then 
          GLi1(:,:,l)=B
        endif
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i,i+1
        cur(:,:,l)=cur(:,:,l)-dble(A(:,:))        
        !-------------------------
        if (present(GGi1)) then
          A=Gp
          call zgemm('c','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G>
          A=Glp(:,:,l)
          call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g> H10 G
          B=C+A         
          GGi1(:,:,l)=B
        endif        
        !-------------------------
        D(:,:)= Gl(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1
        if (present(GRi1)) then 
          GRi1(:,:,l)=GN0
        endif
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
        G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i
        ldos(:,:,l)=G00(:,:)
        !-------------------------
        A(:,:)=Gn(:,:)     
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        Gn(:,:)= Gln(:,:,l) + C(:,:)
        A(:,:)=Gln(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)         
        Gn(:,:)= Gn(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)          
        Gn(:,:)= Gn(:,:)+C(:,:)!     					 !!! G<_i,i
        !-------------------------
        A(:,:)=Gp(:,:)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        !
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        !
        Gp(:,:)= Glp(:,:,l) + C(:,:)
        A(:,:)=Glp(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        !
        Gp(:,:)= Gp(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)     
        Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G>_i,i
        !-------------------------    
        ndens(:,:,l)=Gn(:,:)
        pdens(:,:,l)=Gp(:,:)
    enddo
    !
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    A=-(sigmal-transpose(conjg(sigmal)))*ferm((E-mul)/(BOLTZ*TEMPr))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmal-transpose(conjg(sigmal)))*(ferm((E-mul)/(BOLTZ*TEMPr))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
        tim=tim+dble(B(i,i)-C(i,i))
    enddo
    tre=tim
    deallocate(Gl)
    deallocate(Gln)
    deallocate(Glp)
end subroutine green_RGF_CMS


! find the inverse of a band matrix A by solving a system of linear equations
! on exit, A contains the band matrix of inv(A)
subroutine invert_banded(A,nn,nb)
integer,intent(in)::nn,nb
complex(8),intent(inout)::A(3*nb+1,nn)
complex(8),allocatable::work(:),B(:,:),X(:,:)
integer,allocatable::ipiv(:)
integer::info,lda,lwork,ldb,i,nrhs
allocate(ipiv(nn))
allocate(work(nn*nn))
lda=3*nb+1
call zgbtrf(nn,nn,nb,nb,A,lda,ipiv,info)
if (info.ne.0) then
  print*,'SEVERE warning: zgbtrf failed, info=',info
  call abort
endif
ldb=1
allocate(B(ldb,nn))
allocate(X(lda,nn))
nrhs=ldb
do i=1,nn
  B=0.0d0
  B(1,i)=1.0d0
  call zgbtrs('N',nn,nb,nb,nrhs,A,lda,ipiv,B,ldb,info)
  if (info.ne.0) then
    print*,'SEVERE warning: zgbtrs failed, info=',info
    call abort
  endif
  X(1:nb*2+1,i)=B(1,i-nb:i+nb)
enddo
A=X
deallocate(B,work,ipiv,X)
end subroutine invert_banded


subroutine invert(A,nn)  
  integer :: info,lda,lwork,nn      
  integer, dimension(:), allocatable :: ipiv
  complex(8), dimension(nn,nn),intent(inout) :: A
  complex(8), dimension(:), allocatable :: work
  allocate(work(nn*nn))
  allocate(ipiv(nn))
  call zgetrf(nn,nn,A,nn,ipiv,info)
  if (info.ne.0) then
    print*,'SEVERE warning: zgetrf failed, info=',info
    A=czero
  else
    call zgetri(nn,A,nn,ipiv,work,nn*nn,info)
    if (info.ne.0) then
      print*,'SEVERE warning: zgetri failed, info=',info
      A=czero
    endif
  endif
  deallocate(work)
  deallocate(ipiv)
end subroutine invert


Function ferm(a)
    Real (8) a,ferm
    ferm=1.0d0/(1.0d0+Exp(a))
End Function ferm


! Sancho-Rubio 
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)    
    complex(8), parameter :: alpha = dcmplx(1.0d0,0.0d0)
    complex(8), parameter :: beta  = dcmplx(0.0d0,0.0d0)
    integer i,j,k,nm,nmax
    COMPLEX(8) :: z
    real(8) :: E,error
    REAL(8) :: TOL=1.0D-100  ! [eV]
    COMPLEX(8), INTENT(IN) ::  S00(nm,nm), H00(nm,nm), H10(nm,nm)
    COMPLEX(8), INTENT(OUT) :: G00(nm,nm), GBB(nm,nm)
    COMPLEX(8), ALLOCATABLE :: A(:,:), B(:,:), C(:,:), tmp(:,:)
    COMPLEX(8), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)
    !COMPLEX(8), ALLOCATABLE :: WORK(:)
    !COMPLEX(8), EXTERNAL :: ZLANGE
    Allocate( H_BB(nm,nm) )
    Allocate( H_SS(nm,nm) )
    Allocate( H_01(nm,nm) )
    Allocate( H_10(nm,nm) )
    Allocate( Id(nm,nm) )
    Allocate( A(nm,nm) )
    Allocate( B(nm,nm) )
    Allocate( C(nm,nm) )
    Allocate( tmp(nm,nm) )
    nmax=100
    z = dcmplx(E,1.0d-3)
    Id=0.0d0
    tmp=0.0d0
    do i=1,nm
        Id(i,i)=1.0d0
        tmp(i,i)=dcmplx(0.0d0,1.0d0)
    enddo
    H_BB = H00
    H_10 = H10
    H_01 = TRANSPOSE( CONJG( H_10 ) )
    H_SS = H00
    do i = 1, nmax
        A = z*S00 - H_BB
        !
        call invert(A,nm)
        !
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 
        H_SS = H_SS + C
        H_BB = H_BB + C
        call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,C,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_01,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,A,nm)  
        H_10 = C    
        H_BB = H_BB + A
        call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 
        H_01 = C 
        ! NORM --> inspect the diagonal of A
        error=0.0d0
        DO k=1,nm
            DO j=1,nm
                error=error+sqrt(aimag(C(k,j))**2+Dble(C(k,j))**2)
            END DO
        END DO	
        tmp=H_SS
        IF ( abs(error) < TOL ) THEN	
            EXIT
        ELSE
        END IF
        IF (i .EQ. nmax) THEN
            write(*,*) 'SEVERE warning: nmax reached in sancho!!!',error
            !call abort
            H_SS=H00
            H_BB=H00
        END IF
    enddo
    G00 = z*S00 - H_SS
    !
    call invert(G00,nm)
    !
    GBB = z*S00 - H_BB
    !
    call invert(GBB,nm)
    !
    Deallocate( tmp )
    Deallocate( A )
    Deallocate( B )
    Deallocate( C )
    Deallocate( H_BB )
    Deallocate( H_SS )
    Deallocate( H_01 )
    Deallocate( H_10 )
    Deallocate( Id )
end subroutine sancho


! write spectrum into file (pm3d map)
subroutine write_spectrum(dataset,i,G,nen,en,length,NB,NS,Lx,coeff,append)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(NB*NS,NB*NS,length,nen)
integer, intent(in)::i,nen,length,NB,NS
real(8), intent(in)::Lx,en(nen),coeff(2)
logical, intent(in), optional::append
integer:: ie,j,ib,k
complex(8)::tr
logical append_
append_ = .false.
if (present(append)) append_ = append
if (append_) then
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown',position='append')
else
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
endif

do ie = 1,nen
    do j = 1,length
      do k=1,NS         
        tr=0.0d0         
        do ib=1,nb
            tr = tr+ G(ib+(k-1)*NB,ib+(k-1)*NB,j,ie)            
        end do
        write(11,'(4E18.4)') (j-1)*Lx*NS+(k-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
      enddo
    end do
    write(11,*)    
end do
close(11)
end subroutine write_spectrum


! calculate scattering collision integral from the self-energy, for diagonal blocks
! I = Sig> G^< - Sig< G^>
subroutine calc_block_collision(Sig_lesser,Sig_greater,G_lesser,G_greater,nen,en,spindeg,nm,nx,I,Ispec)
complex(8),intent(in),dimension(nm,nm,nx,nen)::G_greater,G_lesser,Sig_lesser,Sig_greater
real(8),intent(in)::en(nen),spindeg
integer,intent(in)::nen,nm,nx
complex(8),intent(out)::I(nm,nm,nx) ! collision integral
complex(8),intent(out),optional::Ispec(nm,nm,nx,nen) ! collision integral spectrum
!----
complex(8),allocatable::B(:,:)
integer::ie,ix
real(8),parameter::tpi=6.28318530718  
  allocate(B(nm,nm))
  I=dcmplx(0.0d0,0.0d0)
  do ix=1,nx
    do ie=1,nen    
      call zgemm('n','n',nm,nm,nm,cone,Sig_greater(:,:,ix,ie),nm,G_lesser(:,:,ix,ie),nm,czero,B,nm)
      call zgemm('n','n',nm,nm,nm,-cone,Sig_lesser(:,:,ix,ie),nm,G_greater(:,:,ix,ie),nm,cone,B,nm) 
      I(:,:,ix)=I(:,:,ix)+B(:,:)
      if (present(Ispec)) Ispec(:,:,ix,ie)=B(:,:)*spindeg    
    enddo
  enddo
  I(:,:,:)=I(:,:,:)*dble(en(2)-en(1))/tpi*spindeg
  deallocate(B)
end subroutine calc_block_collision


! calculate bond current within each block using I_ij = H_ij G<_ji - H_ji G^<_ij
subroutine calc_block_current(H,G_lesser,cur,nen,en,spindeg,nb,ns,nm,nx,tot_cur,tot_ecur,jdens)
complex(8),intent(in)::H(nm,nm,nx),G_lesser(nm,nm,nx,nen)
real(8),intent(in)::cur(nm,nm,nx,nen)
real(8),intent(in)::en(nen),spindeg
integer,intent(in)::nen,nm,nx,nb,ns ! number of E and device dimension
real(8),intent(out)::tot_cur(nb,nb,nx*ns) ! total bond current density
real(8),intent(out),optional::tot_ecur(nb,nb,nx*ns) ! total bond energy current density
real(8),intent(out),optional::jdens(nb,nb,nx*ns,nen) ! energy resolved bond current density
!----
complex(8),allocatable::B(:,:)
integer::ie,io,jo,ix,ib,i,j,l
real(8),parameter::tpi=6.28318530718  
  allocate(B(nb,nb))
  tot_cur=0.0d0  
  tot_ecur=0.0d0
  do ie=1,nen
    do ix=1,nx
      do ib=1,(ns-1)
        do io=1,nb
          do jo=1,nb
            i=io+(ib-1)*nb
            j=jo+(ib)*nb
            B(io,jo)=H(i,j,ix)*G_lesser(j,i,ix,ie) - H(j,i,ix)*G_lesser(i,j,ix,ie)
          enddo
        enddo    
        B=B*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)
        l=(ix-1)*ns+ib
        if (present(jdens)) jdens(:,:,l,ie) = dble(B)
        if (present(tot_ecur)) tot_ecur(:,:,l)=tot_ecur(:,:,l) + en(ie)*dble(B)
        tot_cur(:,:,l)=tot_cur(:,:,l) + dble(B)          
      enddo
      l=ix*ns            
      if (present(jdens)) jdens(:,:,l,ie) = dble(B)
      if (present(tot_ecur)) tot_ecur(:,:,l)=tot_ecur(:,:,l) + en(ie)*dble(B)
      tot_cur(:,:,l)=tot_cur(:,:,l) + dble(B)      
    enddo
  enddo
  deallocate(B)
end subroutine calc_block_current



! write current spectrum into file (pm3d map)
subroutine write_current_spectrum(dataset,i,cur,nen,en,length,NB,Lx)
  character(len=*), intent(in) :: dataset
  real(8), intent(in),dimension(:,:,:,:) :: cur  
  integer, intent(in)::i,nen,length,NB
  real(8), intent(in)::Lx,en(nen)  
  integer:: ie,j,ib,jb
  real(8)::tr  
  open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
  do ie = 1,nen
    do j = 1,length      
      tr=0.0d0                
      do ib=1,nb  
        do jb=1,nb        
          tr = tr+ cur(ib,jb,j,ie)
        enddo                        
      enddo      
      write(11,'(3E18.4)') dble(j-1)*Lx, en(ie), dble(tr)
    enddo
    write(11,*)    
  enddo
  close(11)
end subroutine write_current_spectrum



! write current spectrum into file (pm3d map)
subroutine write_current_spectrum_summed_over_kz(dataset,i,cur,nen,en,length,NB,Lx,nk, append)
  character(len=*), intent(in) :: dataset
  real(8), intent(in) :: cur(:,:,:,:,:)
  integer, intent(in)::i,nen,length,NB,nk
  real(8), intent(in)::Lx,en(nen)
  integer:: ie,j,ib,jb,ik
  real(8)::tr
  logical, intent(in), optional::append
  logical append_
  append_ = .false.
  if (present(append)) append_ = append

  if (append_) then
    open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown',position='append')
  else
    open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  endif  

  do ie = 1,nen
    do j = 1,length-1
      tr=0.0d0
      !$omp parallel default(shared) private(ik, ib, jb) reduction(+:tr)
      !$omp do         
      do ik=1,nk
        do ib=1,nb  
          do jb=1,nb        
            tr = tr+ cur(ib,jb,j,ie,ik)
          enddo                        
        enddo
      enddo
      !$omp end do 
      !$omp end parallel
      write(11,'(3E18.4)') dble(j-1)*Lx, en(ie), dble(tr)/dble(nk)
    enddo
    write(11,*)    
  enddo
  close(11)
end subroutine write_current_spectrum_summed_over_kz


! write current into file 
subroutine write_current_summed_over_k(dataset,i,cur,nx,NB,Lx,nk)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: cur(:,:,:,:)
integer, intent(in)::i,nx,NB,nk
real(8), intent(in)::Lx
integer:: j,ib,jb,ii,ik
real(8)::tr
  open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
  do ii = 1,nx-1
    tr=0.0d0          
    do ib=1,nb  
      do jb=1,nb       
        do ik=1,nk
          tr = tr + cur(ib,jb,ii,ik)
        enddo
      enddo                        
    end do
    write(11,'(2E18.4)') dble(ii-1)*Lx, tr/dble(nk)
  end do
end subroutine write_current_summed_over_k


! write spectrum into file (pm3d map)
subroutine write_spectrum_summed_over_k(dataset,i,G,nen,en,nk,length,NB,NS,Lx,coeff,append)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(NB*NS,NB*NS,length,nen,nk)
integer, intent(in)::i,nen,length,NB,NS,nk
real(8), intent(in)::Lx,en(nen),coeff(2)
logical, intent(in), optional::append
integer:: ie,j,ib,k,ik
complex(8)::tr
logical append_
append_ = .false.
if (present(append)) append_ = append

if (append_) then
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown',position='append')
else
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
endif
! open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
  do j = 1,length
    do k=1,NS         
      tr=0.0d0
      !$omp parallel default(shared) private(ik, ib) reduction(+:tr)
      !$omp do       
      do ik=1,nk
        do ib=1,nb
          tr = tr+ G(ib+(k-1)*NB,ib+(k-1)*NB,j,ie,ik)            
        enddo
      enddo
      !$omp end do 
      !$omp end parallel
      tr=tr/dble(nk)
      write(11,'(4E18.4)') (j-1)*Lx*NS+(k-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
    enddo
  end do
  write(11,*)    
end do
close(11)
end subroutine write_spectrum_summed_over_k

! write dos into file
subroutine write_dos_summed_over_k(dataset,i,G,nen,en,nk,length,NB,NS,Lx,append)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(NB*NS,NB*NS,length,nen,nk)
integer, intent(in)::i,nen,length,NB,NS,nk
real(8), intent(in)::en(nen),Lx
logical, intent(in), optional::append
integer:: ie,j,ib,k,ik
complex(8)::tr
logical append_
append_ = .false.
if (present(append)) append_ = append
if (append_) then
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown',position='append')
else
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
endif

do ie = 1,nen
  tr=0.0d0         
  do j=1,length
    do k=1,NS             
      do ik=1,nk
        do ib=1,nb
          tr = tr+ G(ib+(k-1)*NB,ib+(k-1)*NB,j,ie,ik)            
        end do
      enddo            
    enddo
  end do  
  tr=tr/dble(nk)/dble(length)/dble(NS)/Lx
  write(11,'(2E18.4)') en(ie), aimag(tr)*(-2.0d0)        
end do
close(11)
end subroutine write_dos_summed_over_k

! write transmission spectrum into file
subroutine write_transmission_spectrum(dataset,i,tr,nen,en)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: tr(:)
integer, intent(in)::i,nen
real(8), intent(in)::en(nen)
integer:: ie,j,ib
open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen    
  write(11,'(2E18.4)') en(ie), dble(tr(ie))      
end do
close(11)
end subroutine write_transmission_spectrum


! write transmission spectrum into file
subroutine write_transmission_spectrum_k(dataset,i,tr,nen,en,nk, append)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: tr(:,:)
integer, intent(in)::i,nen,nk
real(8), intent(in)::en(nen)
integer:: ie,j,ib,ik
real(8):: sumtr(nen)

logical, intent(in), optional::append
logical append_
append_ = .false.
if (present(append)) append_ = append

if (append_) then
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown',position='append')
else
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
endif
! open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
sumtr=0.0d0
do ik=1,nk
  !open(unit=11,file=trim(dataset)//'kz'//TRIM(STRING(ik))//'_'//TRIM(STRING(i))//'.dat',status='unknown')
  !do ie = 1,nen    
  !  write(11,'(2E18.4)') en(ie), dble(tr(ie,ik))      
  !end do
  !close(11)
  sumtr=sumtr+tr(:,ik)
enddo
sumtr=sumtr/dble(nk)
! open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen    
  write(11,'(2E18.4)') en(ie), dble(sumtr(ie))      
end do
close(11)
end subroutine write_transmission_spectrum_k



! write trace of diagonal blocks
subroutine write_trace(dataset,i,G,length,NB,Lx,coeff,E)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:,:)
integer, intent(in)::i,length,NB
real(8), intent(in)::Lx,coeff(2)
real(8), intent(in),optional::E
integer:: ie,j,ib,ix
complex(8)::tr
open(unit=11,file=trim(dataset)//'_'//TRIM(STRING(i))//'.dat',status='unknown', position="append", action="write")
do j = 1,length
    tr=0.0d0          
    do ib=1,nb
        tr = tr+ G(ib,ib,j)            
    end do
    if (.not.(present(E))) then
     write(11,'(3E18.4)') (j-1)*Lx, dble(tr)*coeff(1), aimag(tr)*coeff(2)        
    else
     write(11,'(4E18.4)') (j-1)*Lx, E, dble(tr)*coeff(1), aimag(tr)*coeff(2)         
    endif
end do
write(11,*)
close(11)
end subroutine write_trace



FUNCTION STRING(inn)
  IMPLICIT NONE
  INTEGER, PARAMETER :: POS= 4
  INTEGER, INTENT(IN) :: inn
  CHARACTER(LEN=POS) :: STRING
  !............................................................
  INTEGER :: cifra, np, mm, num  
  IF (inn > (10**POS)-1) stop "ERRORE: (inn > (10**3)-1)  in STRING"
  num= inn
  DO np= 1, POS
     mm= pos-np
     cifra= num/(10**mm)            
     STRING(np:np)= ACHAR(48+cifra)
     num= num - cifra*(10**mm)
  END DO
END FUNCTION STRING

subroutine open_boundary_conditions(nm,M00,M01,M10,V10,xR,dM,dV,cond)
integer,intent(in)::nm
complex(8),intent(in),dimension(nm,nm)::M00,M01,M10,V10
complex(8),intent(out),dimension(nm,nm)::xR,dM,dV
real(8),intent(out)::cond
complex(8),dimension(nm,nm)::tmp1
call surface_function(nm,M00,M01,M10,xR,cond);
!dM=M01*xR*M10
call zgemm('n','n',nm,nm,nm,cone,M01,nm,xR,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dM,nm)
!dV=M01*xR*V10
call zgemm('n','n',nm,nm,nm,cone,M01,nm,xR,nm,czero,tmp1,nm)
call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,V10,nm,czero,dV,nm)
end subroutine open_boundary_conditions


! a slightly modified version of sancho
subroutine surface_function(nm,M00,M01,M10,SF,cond)
integer,intent(in)::nm
complex(8),intent(in),dimension(nm,nm)::M00,M01,M10
complex(8),intent(out),dimension(nm,nm)::SF
real(8),intent(out)::cond
real(8)::cond_limit
integer::max_iteration,IC
complex(8),dimension(:,:),allocatable::alpha,beta,Eps,Eps_surf,inv_element,a_i_b,b_i_a,i_alpha,i_beta
allocate(alpha(nm,nm))
allocate(beta(nm,nm))
allocate(Eps(nm,nm))
allocate(Eps_surf(nm,nm))
allocate(inv_element(nm,nm))
allocate(i_alpha(nm,nm))
allocate(i_beta(nm,nm))
allocate(a_i_b(nm,nm))
allocate(b_i_a(nm,nm))
cond=1.0d10;
cond_limit=1.0d-10;
max_iteration=5000;
IC=1;
alpha=M01
beta=M10
Eps=M00
Eps_surf=M00
do while ((cond>cond_limit).and.(IC<max_iteration))      
    inv_element=Eps
    call invert(inv_element,nm)
    i_alpha=matmul(inv_element,alpha)
    i_beta=matmul(inv_element,beta)
    a_i_b=matmul(alpha,i_beta)
    b_i_a=matmul(beta,i_alpha)
    Eps=Eps-a_i_b-b_i_a
    Eps_surf=Eps_surf-a_i_b
    alpha=matmul(alpha,i_alpha)
    beta=matmul(beta,i_beta)
    !
    cond=sum(abs(alpha)+abs(beta))/2.0d0;
    !
    IC=IC+1;
end do
if (cond>cond_limit) then 
  write(*,*) 'SEVERE warning: nmax reached in surface function!!!',cond
endif
call invert(Eps_surf,nm)
SF=Eps_surf
deallocate(alpha,beta,Eps,Eps_surf,inv_element,a_i_b,b_i_a,i_alpha,i_beta)
end subroutine surface_function


FUNCTION eigv(NN, A)
implicit none
INTEGER, INTENT(IN) :: NN
COMPLEX(8), INTENT(INOUT), DIMENSION(:,:) :: A
REAL(8) :: eigv(NN)
real(8) :: W(1:NN)
integer :: INFO,LWORK,liwork, lrwork
complex(8), allocatable :: work(:)
real(8), allocatable :: RWORK(:)
!integer, allocatable :: iwork(:) 
lwork= max(1,2*NN-1)
lrwork= max(1,3*NN-2)
allocate(work(lwork))
allocate(rwork(lrwork))
!
CALL zheev( 'V','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )
!
deallocate(work,rwork)
if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
   call abort
endif
eigv(:)=W(:)
END FUNCTION eigv

end module green_rgf
