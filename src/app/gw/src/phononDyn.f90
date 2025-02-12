!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 01/2023
! Phonon dynamics
MODULE phononDyn

IMPLICIT NONE

private
!!real(8),allocatable::C(:,:,:,:,:,:,:)  ! 2nd order IFC, (dir,dir,atom1,atom2,cell_x,cell_y,cell_z) [eV/ang^2]
real(8),allocatable::W(:,:,:,:,:,:,:)  ! mass scaled dynamic matrix (dir,dir,atom1,atom2,cell_x,cell_y,cell_z) [1/sec^2]
real(8),allocatable::atom_pos(:,:)  ! atom positions (unit cell) [ang]
real(8),allocatable::atom_mass(:)  ! atom mass (a.u.)
integer,allocatable::atom_type(:)  ! atom type
character(4),allocatable::atom_symbol(:)  ! element symbol
real(8),allocatable::Zeu(:,:,:)  ! Effective charge tensor
real(8),allocatable::eps_inf(:,:)  ! HF dielectric tensor
REAL(8), DIMENSION(3) :: alpha,beta,gamm,xhat,yhat,b1,b2
REAL(8) :: Lx, Ly,cell(3,3) ! in Ang
INTEGER :: zmin,zmax,ymin,ymax,xmin,xmax,nx,ny,nz,NA,NB
! neighbor index range 
! nb is the number of phonon bands = number of atoms NA*3
COMPLEX(8), PARAMETER :: zzero = dcmplx(0.0d0,0.0d0)
COMPLEX(8), PARAMETER :: z1j = dcmplx(0.0d0,1.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0


CONTAINS

SUBROUTINE phon_free_memory()
!!  if (allocated(C)) deallocate(C)
  if (allocated(W)) deallocate(W)
  if (allocated(atom_pos)) deallocate(atom_pos)
  if (allocated(atom_mass)) deallocate(atom_mass)
  if (allocated(atom_type)) deallocate(atom_type)
  if (allocated(atom_symbol)) deallocate(atom_symbol)
  if (allocated(Zeu)) deallocate(Zeu)
  if (allocated(eps_inf)) deallocate(eps_inf)
END SUBROUTINE phon_free_memory

! load from a force-constant file
SUBROUTINE phon_load_from_file(fid,fileformat)
implicit none
integer, intent(in) :: fid
character(len=*), intent(in)::fileformat
  
  select case (fileformat)
  
    case ('qe')
      call phon_load_from_file_qe(fid)
    
    case default
      call phon_load_from_file_qe(fid)
  end select
END SUBROUTINE phon_load_from_file


! load from a QuantumEspresso ifc file produced by q2r.x
SUBROUTINE phon_load_from_file_qe(fid)
implicit none
integer, intent(in) :: fid
  
END SUBROUTINE phon_load_from_file_qe

END MODULE phononDyn

