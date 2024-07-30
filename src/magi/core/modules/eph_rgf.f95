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
module eph_rgf 
    use parameters_mod
    use output
    use legendre
    use observ
    use open_boundary
    use omp_lib
    use rgf
    implicit none 

    contains


    subroutine solve_eph_rgf_3d(nx,mm,nm, nen, energies, nphiy, nphiz, niter, scba_tol, spindeg, mus, mud, temps, tempd, &
        Hii, H1i, Sii, G_r, G_lesser, G_greater, Jdens, tr, tre, output_files,nb,Lx)
        integer, intent(in) :: mm !! max size of blocks
        integer, intent(in) :: nx !! lenght of the device    
        integer, intent(in) :: nen !! number of energies  
        integer, intent(in) :: niter !! max number of SCBA iterations
        real(dp), intent(in) :: scba_tol !! SCBA error tolerance 
        integer, intent(in) :: nb, nphiy, nphiz 
        real(dp), intent(in) :: Lx 
        real(dp), intent(in) :: spindeg !! spin degeneracy
        complex(dp), intent(in) :: Hii(mm,mm,nx,nphiy*nphiz), H1i(mm,mm,nx+1,nphiy*nphiz), Sii(mm,mm,nx,nphiy*nphiz) !! H and overlap
        real(dp), intent(in):: energies(nen)
        real(dp), intent(in):: mus,mud,temps,tempd     
        integer, intent(in) :: nm(nx) !! size of each block
        logical, intent(in) :: output_files
        complex(dp), intent(out), dimension(mm,mm,nx,nen,nphiy*nphiz) :: G_greater, G_lesser, G_r, Jdens
        real(dp), intent(out) :: tr(nen,nphiy*nphiz), tre(nen,nphiy*nphiz)   
        ! ----
        complex(dp),allocatable,dimension(:,:,:,:) :: sigma_lesser_ph, sigma_r_ph
        complex(dp),allocatable,dimension(:,:,:,:) :: sigma_lesser_ph_new, sigma_r_ph_new
        real(dp), dimension(mm,mm) :: mul, mur, TEMPr, TEMPl
        character(len=50) :: dataset_name
        real(dp) :: scba_error
        integer::iter
        integer::ik
        !
        allocate(sigma_lesser_ph(mm,mm,nx,nen), source=czero)
        allocate(sigma_r_ph(mm,mm,nx,nen), source=czero)
        allocate(sigma_lesser_ph_new(mm,mm,nx,nen), source=czero)
        allocate(sigma_r_ph_new(mm,mm,nx,nen), source=czero)
        mul(:,:)=mus
        mur(:,:)=mud
        TEMPl(:,:)=temps
        TEMPr(:,:)=tempd
        !
        scba_error=1.0_dp
        iter=0
        !
        print '(a8,f15.4,a8,f15.4)', 'mus=',mus,'mud=',mud
        !
        do while ( (scba_error>=scba_tol).and.(iter<=niter)  )
            !
            do ik=1,nphiy*nphiz
                call rgf_energies(nx,mm,nm, nen, energies, mul, mur, TEMPl, TEMPr, &
                    Hii(:,:,:,ik), H1i(:,:,:,ik), Sii(:,:,:,ik), &
                    sigma_lesser_ph(:,:,:,:), sigma_r_ph(:,:,:,:), &
                    G_r(:,:,:,:,ik), G_lesser(:,:,:,:,ik), G_greater(:,:,:,:,ik), &
                    Jdens(:,:,:,:,ik), tr(:,ik), tre(:,ik), verbose=.false.)
            enddo
            !
            iter=iter+1
        enddo
        !
        tr = tr *e_charge/twopi/hbar*e_charge*dble(spindeg)/dble(nphiz)/dble(nphiy)    
        tre = tre *e_charge/twopi/hbar*e_charge*dble(spindeg)/dble(nphiz)/dble(nphiy)    
        !
        if (output_files) then 
            dataset_name = 'eph_ldos'
            call write_rgf_spectrum_summed_over_kz(dataset_name,iter,G_r(:,:,:,:,:),nen,energies,nphiy*nphiz,nx,nm,nb,Lx,[1.0d0,-2.0d0])
        endif
    end subroutine solve_eph_rgf_3d


end module eph_rgf