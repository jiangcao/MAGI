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
module setup_mod
    USE utilities
    USE parameters_mod,only : dp
    USE wannierHam3d, only :  w90_load_from_file

    implicit none
    character(8)  :: date
    character(10) :: time
    integer:: fu,rc
    
    logical :: lread_input_bse=.false.
    logical :: lread_input_post=.false.
    logical :: lread_input_gw=.false.
    logical :: lread_input_ham=.false.
    
    !system variables
    integer :: nkx=1
    integer :: nky=1
    integer :: nkz=1
    logical :: reorder_axis=.false.
    integer :: ncpu=1
    integer ::xyz(3)=(/1,2,3/)

    !bse_input variables        
    integer :: ncb=1
    integer :: nvb=1
    integer :: nex=1 ! number of exciton states we want to plot
    integer :: nx=50
    integer :: ny=50
    integer :: nz=10
    integer :: exciton_q(3) = 0
    real(kind=dp) :: potscale=1.0_dp
    real(kind=dp) :: F=0.0_dp
    real(kind=dp) :: xi=6.5_dp
    real(kind=dp) :: epsilon=1.0_dp
    real(kind=dp) :: k0=1.0_dp
    real(kind=dp) :: E_cutoff=4.0_dp
    real(kind=dp) :: sig
    real(kind=dp) :: e1(3)=(/1_dp,0_dp,0_dp/)
    real(kind=dp) :: dx
    real(kind=dp) :: dy
    real(kind=dp) :: dz
    real(kind=dp) :: E_scissor=0.0_dp
    real(kind=dp) :: rsmear=3.0_dp ! Angstrom
    logical :: lwcenter=.true.
    logical :: lbse=.true.
    logical :: lplotexciton=.false.
    logical :: lread_screened_coulomb=.true.
    logical :: lexchange=.true.
    logical :: ldirect=.true.

    !input variables from input_post file
    character(len=100) :: dataset
    integer  ::number=0
    
    
    !input variables from input file
    integer :: ns=1
    integer ::  pottype=1
    integer :: niter=1
    integer :: length=1
    integer :: nen=1    
    integer :: ndiag=1
    real(kind=dp) :: emax
    real(kind=dp) :: emin
    real(kind=dp) :: eps_screen
    real(kind=dp) :: r0
    real(kind=dp) :: mus
    real(kind=dp) :: mud
    real(kind=dp) :: temps
    real(kind=dp) :: tempd
    real(kind=dp) :: alpha_mix
    real(kind=dp) :: ky_shift
    real(kind=dp) :: kz_shift
    real(kind=dp) :: hw
    real(kind=dp) :: intensity
    real(kind=dp) :: polarization(3)
    real(kind=dp) :: ita
    logical :: ltrans=.false.
    logical :: lreadpot=.false.
    logical :: lrcoulomb=.false.
    logical :: ldiag=.true.
    logical :: lrgf=.false.
    logical :: lkz=.false.
    logical :: lephot=.false.
    logical :: lnogw=.false.
    logical :: labs=.false.
    logical :: lvertex=.false.    
    logical :: lflatband=.false.

        
    ! MPI variables
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) comm_size
    integer ( kind = 4 ) comm_rank
    integer ( kind = 4 ) local_Nkz
    integer ( kind = 4 ) local_Nky
    integer ( kind = 4 ) local_NE
    integer ( kind = 4 ) first_local_energy

    namelist /system/ NKX,NKY,NKZ,ncpu,reorder_axis,xyz
    namelist /input_bse/ ncb,nvb,F,epsilon,k0,E_cutoff,sig,lwcenter,lbse,e1,lplotexciton,nex,nx,ny,nz,dx,dy,dz,E_scissor,nen,xi,lread_screened_coulomb,lexchange,ldirect,rsmear,exciton_q
    namelist /input_post/ dataset,number
    namelist /input_gw/ ns,niter,length,emin,emax,nen,eps_screen,r0,mus,mud,temps,tempd,alpha_mix,ky_shift,kz_shift,hw,intensity,polarization,ita,ltrans,lreadpot,lrcoulomb,ldiag,lrgf,lkz,lephot,lnogw,labs,lvertex,lflatband,lbse,pottype,ndiag


    CONTAINS
    
    SUBROUTINE begin_setup()
            
        include "mpif.h"
        call MPI_Init(ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierr)
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        
        if (comm_rank == 0) then
          print *, 'Comm Size =', comm_size
        else
          print *, 'Comm Rank =', comm_rank
        endif
        
        call MPI_Barrier(MPI_COMM_WORLD, ierr)

        call date_and_time(DATE=date,TIME=time)
        print '(a,2x,a)', date, time

        open(action='read', file='input', iostat=rc, newunit=fu)
        read(nml=system, iostat=rc, unit=fu)
        close(fu)                              

        if (lread_input_bse) then
            open(action='read', file='input', iostat=rc, newunit=fu)
            read(nml=input_bse, iostat=rc, unit=fu)
            close(fu)                              
        end if

        if (lread_input_post) then
            open(action='read', file='input', iostat=rc, newunit=fu)
            read(nml=input_post, iostat=rc, unit=fu)
            close(fu)
        endif

        if (lread_input_gw) then
            open(action='read', file='input', iostat=rc, newunit=fu)
            read(nml=input_gw, iostat=rc, unit=fu)
            close(fu)
        endif

        if (lread_input_ham) then
            open(unit=10,file='ham_dat',status='unknown')
            call w90_load_from_file(10, reorder_axis,xyz)
            close(10)   
        endif
        
        
        call omp_set_num_threads(ncpu)
        call mkl_set_num_threads(ncpu)

    END SUBROUTINE begin_setup

end module setup_mod
