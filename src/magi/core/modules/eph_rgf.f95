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
    use gw_dense, only: map_kq_2d
    implicit none 

    contains

    ! solve SCBA with electron-phonon interaction under simple deformation potential approximation
    subroutine solve_eph_rgf_3d(nx,mm,nm, nen, energies, nphiy, nphiz, niter, scba_tol, &
        spindeg, num_phmode, Dop, Nop, alpha_mix, mus, mud, temps, tempd, &
        Hii, H1i, Sii, G_r, G_lesser, G_greater, Jdens, tr, tre, output_files,nb,Lx,midgap,nelec,pelec)
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
        real(dp), intent(in):: alpha_mix !! SCBA mixing parameter
        real(dp), intent(in):: mus,mud,temps,tempd     
        integer, intent(in) :: nm(nx) !! size of each block
        integer,intent(in) :: num_phmode  !! number of optical phonon modes
        real(dp),intent(in) :: Dop(num_phmode)  !! optical deformation potentials
        integer,intent(in)  :: Nop(num_phmode)  !! optical phonon freq. in unit of energy step        
        logical, intent(in) :: output_files
        real(dp),intent(in) :: midgap(nx)  
        complex(dp), intent(out), dimension(mm,mm,nx,nen) :: G_greater, G_lesser, G_r, Jdens
        real(dp), intent(out) :: tr(nen,nphiy*nphiz), tre(nen,nphiy*nphiz)   
        real(dp),intent(out) ::nelec(nx*mm),pelec(nx*mm)
        ! ----
        complex(dp),allocatable,dimension(:,:,:,:) :: sigma_lesser_ph, sigma_r_ph
        complex(dp), dimension(mm,mm,nx,nen) :: G_greater_k, G_lesser_k, G_r_k, Jdens_k
        complex(dp),allocatable,dimension(:,:,:,:) :: sigma_lesser_ph_new, sigma_r_ph_new, sigma_greater_ph_new
        real(dp), dimension(mm,mm) :: mul, mur, TEMPr, TEMPl
        character(len=50) :: dataset_name
        real(dp) :: scba_error, n_bose, dE, dkt, tr_old(nen,nphiy*nphiz)
        integer::iter
        integer::ik, ph_mode
        !
        allocate(sigma_lesser_ph(mm,mm,nx,nen), source=czero)
        allocate(sigma_r_ph(mm,mm,nx,nen), source=czero)
        allocate(sigma_lesser_ph_new(mm,mm,nx,nen), source=czero)
        allocate(sigma_greater_ph_new(mm,mm,nx,nen), source=czero)
        allocate(sigma_r_ph_new(mm,mm,nx,nen), source=czero)
        mul(:,:)=mus
        mur(:,:)=mud
        TEMPl(:,:)=temps
        TEMPr(:,:)=tempd
        !
        scba_error=1.0_dp
        tr_old=0.0_dp
        iter=0
        dE = energies(2) - energies(1)
        dkt= 1.0_dp/dble(nphiy)/dble(nphiz)
        !
        print '(a8,f15.4,a8,f15.4)', 'mus=',mus,'mud=',mud
        !
        do while ( (scba_error>=scba_tol).and.(iter<=niter) )
            !
            print *,''
            print *,'  calc G ...'  
            G_r=czero
            G_greater=czero
            G_lesser=czero
            tr_old=tr
            do ik=1,nphiy*nphiz
                call rgf_energies(nx,mm,nm, nen, energies, mul, mur, TEMPl, TEMPr, &
                    Hii(:,:,:,ik), H1i(:,:,:,ik), Sii(:,:,:,ik), &
                    sigma_lesser_ph(:,:,:,:), sigma_r_ph(:,:,:,:), &
                    G_r_k(:,:,:,:), G_lesser_k(:,:,:,:), G_greater_k(:,:,:,:), &
                    Jdens_k(:,:,:,:), tr(:,ik), tre(:,ik), verbose=.false.)
                G_r = G_r + G_r_k * dkt
                G_lesser = G_lesser + G_lesser_k * dkt
                G_greater= G_greater + G_greater_k * dkt
                Jdens = Jdens + Jdens_k * dkt    
            enddo
            tr = tr *dE/twopi*e_charge/twopi/hbar*e_charge*dble(spindeg) * dkt  
            tre = tre *dE/twopi*e_charge/twopi/hbar*e_charge*dble(spindeg) * dkt
            !
            open(unit=101,file='eph_Id_iteration.dat',status='unknown',position='append')
            write(101,'(I4,2E16.6)') iter, -sum(tr(:,:)), sum(tre(:,:))
            close(101)
            !
            ! compute the e-phonon self-energies
            print *,'  calc Sig ...'  
            sigma_lesser_ph_new = czero
            sigma_greater_ph_new = czero
            do ph_mode=1,num_phmode
                ! Bose-Einstein
                n_bose=1.0_dp/(EXP((dble(Nop(ph_mode))*dE)/(BOLTZ*((temps+tempd)/2.0_dp)))-1.0_dp)            
                call selfenergy_eph_rgf_simple(nm=mm,nx=nx,nen=nen,en=energies,nop=Nop(ph_mode),Dop=Dop(ph_mode),&
                                G_lesser=G_lesser,G_greater=G_greater,&
                                Sig_lesser=sigma_lesser_ph_new,Sig_greater=sigma_greater_ph_new,&
                                n_bose=n_bose,init_selfenergy=.false.)
            enddo
            ! scba error
            sigma_r_ph_new = dcmplx( 0.0_dp, aimag(sigma_greater_ph_new - sigma_lesser_ph_new)/2.0d0 )
            scba_error = sqrt( sum( abs(sigma_r_ph_new - sigma_r_ph)**2 ) / sum( abs(sigma_r_ph_new)**2 ) )
            scba_error = ( scba_error + sum( abs(tr-tr_old)**2 ) / sum( abs(tr)**2 ) ) / 2.0_dp
            open(unit=101,file='eph_scba_error.dat',status='unknown',position='append')
            write(101,'(I4,E16.6)') iter, scba_error
            close(101)            
            write(*,'("+ iter=",I8,"  error=",2E16.6)') iter, scba_error
            write(*,'("   IDS=",2E16.6)') -sum(tr(:,:)), sum(tre(:,:))
            iter=iter+1
            ! mixing self-energies with the previous ones
            sigma_r_ph = sigma_r_ph+ alpha_mix * (sigma_r_ph_new -sigma_r_ph)
            sigma_lesser_ph = sigma_lesser_ph+ alpha_mix * (sigma_lesser_ph_new -sigma_lesser_ph)      
        enddo
        !
        print *,'  calc G last time'  
        G_r=czero
        G_greater=czero
        G_lesser=czero
        do ik=1,nphiy*nphiz
            call rgf_energies(nx,mm,nm, nen, energies, mul, mur, TEMPl, TEMPr, &
                Hii(:,:,:,ik), H1i(:,:,:,ik), Sii(:,:,:,ik), &
                sigma_lesser_ph(:,:,:,:), sigma_r_ph(:,:,:,:), &
                G_r_k(:,:,:,:), G_lesser_k(:,:,:,:), G_greater_k(:,:,:,:), &
                Jdens_k(:,:,:,:), tr(:,ik), tre(:,ik), verbose=.false.)
            G_r = G_r + G_r_k * dkt
            G_lesser = G_lesser + G_lesser_k * dkt
            G_greater = G_greater + G_greater_k * dkt
            Jdens = Jdens + Jdens_k * dkt    
        enddo
        ! current spectra from the leads
        tr = tr *dE/twopi*e_charge/twopi/hbar*e_charge*dble(spindeg) * dkt    
        tre = tre *dE/twopi*e_charge/twopi/hbar*e_charge*dble(spindeg) * dkt    
        ! output files
        if (output_files) then 
            dataset_name = 'eph_ldos'
            call write_rgf_spectrum(dataset_name,iter,G_r,nen,energies,nx,nm,nb,Lx,[1.0d0,-2.0d0])
            dataset_name = 'eph_cur'
            call write_rgf_spectrum(dataset_name,iter,Jdens,nen,energies,nx,nm,nb,Lx,[1.0d0,0.0d0])
        endif
        ! compute the charges
        call calc_charge_rgf(G_lesser,G_greater,nen,energies,nm,mm,nx,nelec,pelec,midgap)
        nelec = nelec * dble(spindeg)
        pelec = pelec * dble(spindeg)
        write(*,'( "  total charges in device =", E16.6,"(n)", E16.6,"(p)" )') sum(nelec),sum(pelec)
        !
    end subroutine solve_eph_rgf_3d

    ! calculate simple (deformation potential approximation) e-phonon self-energies 
    subroutine selfenergy_eph_rgf_simple(nm,nx,nen,En,nop,Dop,G_lesser,G_greater,&
        Sig_lesser,Sig_greater,n_bose,init_selfenergy)
        integer,intent(in)::nm,nen,nx
        integer,intent(in)::nop !! phonon freq. in unit of energy discretization step
        real(dp),intent(in)::en(nen),n_bose 
        real(dp),intent(in)::Dop !! deformation potential 
        logical,intent(in)::init_selfenergy !! initialize the self-energies to zero
        complex(dp),intent(in),dimension(nm,nm,nx,nen)::G_lesser,G_greater !! Green's functions
        complex(dp),intent(inout),dimension(nm,nm,nx,nen)::Sig_lesser,Sig_greater !! accumulate the e-ph self-energies         
        !---------
        integer::ie,i,ix
        real(8)::dE 
        if (init_selfenergy) then
            Sig_lesser = czero
            Sig_greater = czero
        endif
        dE = (en(2)-en(1)) / twopi                     
        ! Sig^<>(E,k) = Dop^2 [ N G^<>(E -+ hw,k-+q) + (N+1) G^<>(E +- hw,k+-q)]        
        !$omp parallel default(shared) private(ie,ix,i)         
        !$omp do
        do ie=1,nen
            ! Sig^<(E,k)
            if (ie-nop>=1) then 
                do ix=1,nx
                    do i=1,nm
                        Sig_lesser(i,i,ix,ie) = Sig_lesser(i,i,ix,ie) + G_lesser(i,i,ix,ie-nop) * n_bose * Dop * dE 
                    enddo
                enddo
            endif                
            if (ie+nop<=nen) then 
                do ix=1,nx
                    do i=1,nm
                        Sig_lesser(i,i,ix,ie) = Sig_lesser(i,i,ix,ie) + G_lesser(i,i,ix,ie+nop) * (n_bose+1.0_dp) * Dop * dE            
                    enddo
                enddo
            endif
            !
            ! Sig^>(E,k)            
            if (ie-nop>=1) then 
                do ix=1,nx
                    do i=1,nm
                        Sig_greater(i,i,ix,ie) = Sig_greater(i,i,ix,ie) + G_greater(i,i,ix,ie-nop) * (n_bose+1.0_dp) * Dop * dE   
                    enddo
                enddo
            endif
            if (ie+nop<=nen) then 
                do ix=1,nx
                    do i=1,nm
                        Sig_greater(i,i,ix,ie) = Sig_greater(i,i,ix,ie) + G_greater(i,i,ix,ie+nop) * n_bose * Dop * dE   
                    enddo
                enddo
            endif			  
        enddo
        !$omp end do   
        !$omp end parallel
    end subroutine selfenergy_eph_rgf_simple

end module eph_rgf
