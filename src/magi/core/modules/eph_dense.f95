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
module eph_dense 
    use parameters_mod
    use output
    use legendre
    use observ
    use open_boundary
    use omp_lib
    use gw_dense
    implicit none 

    contains

    subroutine solve_eph_3d(niter,scba_tol,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
        alpha_mix,nen,nsub,En,nb,ns,nphiy,nphiz,Ham,H00lead,H10lead,T,&
        ndiag,num_lead,Dop,nop,midgap,&
        flatband,output_files,G_retarded,G_lesser,G_greater,tr,nelec,pelec)
        ! 
        integer, intent(in) :: nen, nsub, nb, ns,niter,nm_dev,length, nphiz, nphiy, num_lead
        real(dp), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg,scba_tol
        complex(dp),intent(in) :: Ham(nm_dev,nm_dev,nphiy*nphiz)
        complex(dp),intent(in) :: H00lead(NB*NS,NB*NS,num_lead,nphiy*nphiz)
        complex(dp),intent(in) :: H10lead(NB*NS,NB*NS,num_lead,nphiy*nphiz)
        complex(dp),intent(in) :: T(NB*NS,nm_dev,num_lead,nphiy*nphiz)        
        integer,intent(in)::ndiag
        logical,intent(in)::flatband
        logical,intent(in) :: output_files
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nen,nsub,nphiy*nphiz) ::  G_retarded,G_lesser,G_greater        
        real(dp),intent(in) :: Dop  
        real(dp),intent(in) :: midgap(nm_dev)  
        integer,intent(in) :: Nop  
        real(dp),intent(out) ::Tr(nen,num_lead) ! current spectrum on leads    
        real(dp),intent(out) ::nelec(nm_dev),pelec(nm_dev)
        ! ------        
        complex(dp),dimension(:,:,:,:),allocatable ::  Sig_retarded,Sig_lesser,Sig_greater
        complex(dp),dimension(:,:,:,:),allocatable ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
        complex(dp),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
        complex(dp),allocatable,dimension(:,:):: B ! tmp matrix
        real(dp),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:),sumcur(:,:,:),sumtot_cur(:,:),sumtot_ecur(:,:)
        complex(dp),allocatable::Ispec(:,:,:),Itot(:,:)    
        real(dp),allocatable::Te(:,:,:) ! transmission matrix spectrum
        real(dp),allocatable::sumTr(:,:) ! current spectrum on leads summed over k
        real(dp),allocatable::sumTe(:,:,:) ! transmission matrix spectrum summed over k
        integer :: iter,ie
        integer :: i,j,nm,l,h,iop,ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,isub        
        complex(dp) :: dE
        real(dp)::mu(2),temp(2)
        real(dp)::weights(nsub),xen(nsub)
        real(dp)::scba_error
        complex(dp),allocatable::Scat_spec(:,:,:,:) ! collision integral spectrum
        complex(dp),allocatable::Scat(:,:) ! collision integral
        integer :: nqy,nqz,ik_start,ik_end
        real(dp) :: n_bose
        logical :: gamma_q

        allocate(Sig_retarded(nm_dev,nm_dev,nen,nphiy*nphiz), source=czero)
        allocate(Sig_lesser(nm_dev,nm_dev,nen,nphiy*nphiz), source=czero)
        allocate(Sig_greater(nm_dev,nm_dev,nen,nphiy*nphiz), source=czero)
        
        allocate(Sig_retarded_new(nm_dev,nm_dev,nen,nphiy*nphiz), source=czero)
        allocate(Sig_lesser_new(nm_dev,nm_dev,nen,nphiy*nphiz), source=czero)
        allocate(Sig_greater_new(nm_dev,nm_dev,nen,nphiy*nphiz), source=czero)                 
        
        scba_error=1.0d0
        !
        print *,'============ green_solve_eph_3D ============'
        dE = En(2) - En(1)
        call gaulegf(0.0d0, dble(dE), xen, weights, nsub) ! obtain the Legendre ordinates and weights    
        !
        allocate(siglead(NB*NS,NB*NS,nen,num_lead,nphiy*nphiz))
        ! get leads sigma
        do ikz=1, nphiy*nphiz
            siglead(:,:,:,1,ikz) = Sig_retarded(1:NB*NS,1:NB*NS,:,ikz)
            siglead(:,:,:,2,ikz) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:,ikz)  
        enddo
        !
        allocate(tot_cur(nm_dev,nm_dev))
        allocate(tot_ecur(nm_dev,nm_dev))
        allocate(sumtot_cur(nm_dev,nm_dev))
        allocate(sumtot_ecur(nm_dev,nm_dev))
        allocate(cur(nm_dev,nm_dev,nen))
        allocate(sumcur(nm_dev,nm_dev,nen))
        allocate(Ispec(nm_dev,nm_dev,nen))
        allocate(Itot(nm_dev,nm_dev))    
        allocate(te(nen,num_lead,num_lead))
        allocate(sumtr(nen,num_lead))
        allocate(sumte(nen,num_lead,num_lead))
        if (flatband) then
            mu=(mus+mud)/2.0d0
            temp=(temps+tempd)/2.0d0
        else
            mu=(/ mus, mud /)
            temp=(/ temps, tempd /)
        endif
        iter=0
        print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
        do while ( (scba_error>=scba_tol).and.(iter<=niter)  )
            print *,'+ iter=',iter,'error=',scba_error
            print *,'  calc G'  
            sumtot_cur=0.0d0
            sumtot_ecur=0.0d0
            sumcur=0.0d0
            sumTr=0.0d0
            sumTe=0.0d0
            do ikz=1,nphiy*nphiz
                do isub=1,nsub
                    call calc_gf(nen,En+xen(isub),num_lead,nm_dev,(/nb*ns,nb*ns/),nb*ns,&
                        Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),&
                        T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),&
                        G_retarded(:,:,:,isub,ikz),G_lesser(:,:,:,isub,ikz),&
                        G_greater(:,:,:,isub,ikz),Tr,Te,mu,temp,flatband)
                    !call write_spectrum('ldos_kz'//string(ikz)//'_',iter,G_retarded(:,:,:,ikz),nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
                    call calc_bond_current(Ham(:,:,ikz),G_lesser(:,:,:,isub,ikz),nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)    
                    call write_current_spectrum('last_Jdens',ikz,cur,nen,en,length,NB,Lx)    
                    call write_spectrum('last_ldos',ikz,G_retarded(:,:,:,:,ikz),nen,nsub,En,xen,length,NB,Lx,(/1.0d0,-2.0d0/))
                    sumcur=sumcur + cur * weights(isub)
                    sumtot_cur=sumtot_cur + tot_cur * weights(isub)
                    sumtot_ecur=sumtot_ecur + tot_ecur * weights(isub)
                    sumTr=sumTr + Tr * weights(isub)
                    sumTe=sumTe + Te * weights(isub)
                enddo
            enddo
            call calc_charge(G_lesser,G_greater,nen,nsub,nphiy*nphiz,En,NS,NB,nm_dev,nelec,pelec,midgap)
            nelec = nelec * dble(spindeg)
            pelec = pelec * dble(spindeg)
            write(*,'( "  total charges in device =", E16.6,"(n)", E16.6,"(p)" )') sum(nelec),sum(pelec)
            sumcur=sumcur/dble(nphiy)/dble(nphiz)
            sumtot_cur=sumtot_cur/dble(nphiy)/dble(nphiz)
            sumtot_ecur=sumtot_ecur/dble(nphiy)/dble(nphiz)
            sumTr=sumTr/dble(nphiz)/dble(nphiy)
            sumTe=sumTe/dble(nphiz)/dble(nphiy)
            sumTr = sumTr *e_charge/twopi/hbar*e_charge*dble(spindeg)
            if (output_files) then
                if (flatband) then
                    print *,'flatband'
                    ! call write_spectrum_per_kz('eph_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
                    ! call write_spectrum_per_kz('eph_gamma-centered_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=0.0_dp,at_kz=0.0_dp)
                    call write_spectrum_summed_over_kz('eph_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
            
                else
                    call write_spectrum_summed_over_kz('eph_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
                    call write_spectrum_summed_over_kz('eph_ndos',iter,G_lesser(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
                    call write_spectrum_summed_over_kz('eph_pdos',iter,G_greater(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
                endif
                call write_current_spectrum('eph_Jdens',iter,sumcur,nen,en,length,NB,Lx)
                call write_current('eph_I',iter,sumtot_cur,length,NB,NS,Lx)
                call write_current('eph_EI',iter,sumtot_ecur,length,NB,NS,Lx)
                call write_transmission_spectrum('eph_trL',iter,sumTr(:,1),nen,En)
                call write_transmission_spectrum('eph_trR',iter,sumTr(:,2),nen,En)
                ! call write_transmission_spectrum('eph_TE_LR',iter,sumTe(:,1,2)*spindeg,nen,En)
                ! call write_transmission_spectrum('eph_TE_RL',iter,sumTe(:,2,1)*spindeg,nen,En)            
            endif            
            open(unit=101,file='eph_Id_iteration.dat',status='unknown',position='append')
            write(101,'(I4,2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
            close(101)
            write(*,'(I4,"  IDS=",2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
            !
            G_retarded=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
            G_lesser=dcmplx(0.0d0*dble(G_lesser),aimag(G_lesser))
            G_greater=dcmplx(0.0d0*dble(G_greater),aimag(G_greater))            
            !                      
            gamma_q=.false.
            n_bose=1.0_dp/(EXP((dble(Nop)*dE)/(BOLTZ*(sum(TEMP)/2.0_dp)))-1.0_dp)
            call selfenergy_eph_simple(nm=nm_dev,nen=nen,en=en,nop=Nop,nky=nphiy,nkz=nphiz,Dop=Dop,&
                                G_lesser=G_lesser,G_greater=G_greater,&
                                Sig_lesser=sig_lesser_new,Sig_greater=sig_greater_new,&
                                n_bose=n_bose,gamma_q=gamma_q)
            !
            Sig_retarded_new = dcmplx( 0.0_dp, aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
            ! 
            if (output_files) then
                print *,'  calc collision integral'
                allocate(Scat_spec(nm_dev,nm_dev,nen,nsub),source=czero)
                allocate(Scat(nm_dev,nm_dev),source=czero)
                call calc_collision(Sig_lesser_new,Sig_greater_new,G_lesser,G_greater,nen,en,nsub,nphiy*nphiz,spindeg,nm_dev,Scat,Scat_spec)
                call write_spectrum('eph_Scat',iter,Scat_spec,nen,nsub,En,xen,length,NB,Lx,(/1.0d0,1.0d0/))                                
                deallocate(Scat,Scat_spec)
            endif
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! symmetrize the selfenergies
            !$omp parallel default(shared) private(ie,ik,B)
            allocate(B(nm_dev,nm_dev))
            !$omp do
            do ie=1,nen
                do ik=1,nphiy*nphiz
                    B(:,:)=transpose(Sig_retarded_new(:,:,ie,ik))
                    Sig_retarded_new(:,:,ie,ik) = (Sig_retarded_new(:,:,ie,ik) + B(:,:))/2.0d0    
                    B(:,:)=transpose(Sig_lesser_new(:,:,ie,ik))
                    Sig_lesser_new(:,:,ie,ik) = (Sig_lesser_new(:,:,ie,ik) + B(:,:))/2.0d0
                    B(:,:)=transpose(Sig_greater_new(:,:,ie,ik))
                    Sig_greater_new(:,:,ie,ik) = (Sig_greater_new(:,:,ie,ik) + B(:,:))/2.0d0
                enddo
            enddo
            !$omp end do
            deallocate(B)
            !$omp end parallel
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            scba_error = sum( abs(Sig_retarded_new - Sig_retarded)**2 ) / sum( abs(Sig_retarded_new)**2 )
            open(unit=101,file='eph_scba_error.dat',status='unknown',position='append')
            write(101,'(I4,E16.6)') iter, scba_error
            close(101)
            iter=iter+1
            ! mixing with previous ones
            Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
            Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
            Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)              
            !
            if (.not. flatband) then
                ! get leads sigma
                do iqz=1,nphiy*nphiz
                    siglead(:,:,:,1,iqz) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:,iqz)
                    siglead(:,:,:,2,iqz) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:,iqz)    
                enddo            
            endif
            if (flatband) then            
                ! call write_spectrum_per_kz('eph_SigR',iter,Sig_retarded,nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
            endif 
            ! call write_spectrum_summed_over_kz('eph_SigR',iter,Sig_retarded,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            !  call write_spectrum_summed_over_kz('SigL',iter,Sig_lesser,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
            !  call write_spectrum_summed_over_kz('SigG',iter,Sig_greater,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
        end do  
        !
        deallocate(Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
        if (iter == niter) then 
            print *, "warning: max number of iterations reached!"        
        endif
        ! calculate GF for the last time    
        print *, 'calc G for the last time'  
        sumtot_cur=0.0d0
        sumtot_ecur=0.0d0
        sumcur=0.0d0
        sumTr=0.0
        sumTe=0.0
        do ikz=1,nphiy*nphiz
            do isub=1,nsub
                call calc_gf(nen,En+xen(isub),num_lead,nm_dev,(/nb*ns,nb*ns/),nb*ns,&
                    Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),&
                    T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),&
                    G_retarded(:,:,:,isub,ikz),G_lesser(:,:,:,isub,ikz),&
                    G_greater(:,:,:,isub,ikz),Tr,Te,mu,temp,flatband)
                !call write_spectrum('ldos_kz'//string(ikz)//'_',iter,G_retarded(:,:,:,ikz),nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
                call calc_bond_current(Ham(:,:,ikz),G_lesser(:,:,:,isub,ikz),nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)    
                !call write_current_spectrum('Jdens_kz'//string(ikz)//'_',iter,cur,nen,en,length,NB,Lx)    
                sumcur=sumcur + cur * weights(isub)
                sumtot_cur=sumtot_cur + tot_cur * weights(isub)
                sumtot_ecur=sumtot_ecur + tot_ecur * weights(isub)
                sumTr=sumTr + Tr * weights(isub)
                sumTe=sumTe + Te * weights(isub)
            enddo
        enddo
        sumcur=sumcur/dble(nphiy)/dble(nphiz)
        sumtot_cur=sumtot_cur/dble(nphiy)/dble(nphiz)
        sumtot_ecur=sumtot_ecur/dble(nphiy)/dble(nphiz)
        sumTr=sumTr/dble(nphiz)/dble(nphiy)
        sumTe=sumTe/dble(nphiz)/dble(nphiy)
        sumTr=sumTr *e_charge/twopi/hbar*e_charge*dble(spindeg)
        if (flatband) then
            print *,'flatband'
            ! call write_spectrum_per_kz('eph_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
            ! call write_spectrum_per_kz('eph_gamma-centered_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=0.0_dp,at_kz=0.0_dp)
            call write_spectrum_summed_over_kz('eph_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))

        else
            call write_spectrum_summed_over_kz('eph_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
            call write_spectrum_summed_over_kz('eph_ndos',iter,G_lesser(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            call write_spectrum_summed_over_kz('eph_pdos',iter,G_greater(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
        endif
        call write_current_spectrum('eph_Jdens',iter,sumcur,nen,en,length,NB,Lx)
        call write_current('eph_I',iter,sumtot_cur,length,NB,NS,Lx)
        call write_current('eph_EI',iter,sumtot_ecur,length,NB,NS,Lx)
        call write_transmission_spectrum('eph_trL',iter,sumTr(:,1),nen,En)
        call write_transmission_spectrum('eph_trR',iter,sumTr(:,2),nen,En)
        ! call write_transmission_spectrum('eph_TE_LR',iter,sumTe(:,1,2)*spindeg,nen,En)
        ! call write_transmission_spectrum('eph_TE_RL',iter,sumTe(:,2,1)*spindeg,nen,En)        
        open(unit=101,file='eph_Id_iteration.dat',status='unknown',position='append')
        write(101,'(I4,2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
        close(101)
        write(*,'(I4,"  IDS=",2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
        !
        tr=sumTr
        deallocate(siglead)
        deallocate(sumcur,cur,tot_cur,tot_ecur,sumtot_cur,sumtot_ecur)
        deallocate(Ispec,Itot)
        deallocate(Te,sumTr,sumTe)    
        deallocate(Sig_greater,Sig_lesser,Sig_retarded)
    end subroutine solve_eph_3d

    ! calculate simple (deformation potential approximation) e-phonon self-energies 
    subroutine selfenergy_eph_simple(nm,nen,En,nop,Dop,nky,nkz,G_lesser,G_greater,&
        Sig_lesser,Sig_greater,n_bose,gamma_q)
        integer,intent(in)::nm,nen,nky,nkz
        integer,intent(in)::nop ! phonon freq. in unit of energy discretization step
        real(dp),intent(in)::en(nen),n_bose 
        logical,intent(in)::gamma_q ! whether q=0
        real(dp),intent(in)::Dop ! deformation potential at q for ith-mode
        complex(dp),intent(in),dimension(nm,nm,nen,nky*nkz)::G_lesser,G_greater ! Green's functions
        complex(dp),intent(inout),dimension(nm,nm,nen,nky*nkz)::Sig_lesser,Sig_greater ! accumulate the e-ph self-energies 
        !---------
        integer::ie,ikd,ik,iq,nq,nk,i
        real(8)::dE 
        Sig_lesser = czero
        Sig_greater = czero
        dE = (en(2)-en(1)) / twopi             
        nk=nky*nkz
        if (gamma_q) then 			
            nq=1
		else
			nq=nky*nkz
		endif
        ! Sig^<>(E,k) = Dop^2 [ N G^<>(E -+ hw,k-+q) + (N+1) G^<>(E +- hw,k+-q)]        
        !$omp parallel default(shared) private(ie,ik,ikd,i)         
        !$omp do
        do ie=1,nen
			do ik=1,nk
                do iq=1,nq
                    ! Sig^<(E,k) 
                    if (gamma_q) then 
                        ikd = ik					
                    else
                        ikd = map_kq_2d(-1,ik,iq,nky,nkz)					
                    endif
                    if (ie-nop>=1) then 
                        do i=1,nm
                            Sig_lesser(i,i,ie,ik) = Sig_lesser(i,i,ie,ik) + G_lesser(i,i,ie-nop,ikd) * n_bose * Dop**2 * dE 
                        enddo
                    endif
                    if (gamma_q) then 
                        ikd = ik
                    else
                        ikd = map_kq_2d(+1,ik,iq,nky,nkz)
                    endif
                    if (ie+nop<=nen) then 
                        do i=1,nm
                            Sig_lesser(i,i,ie,ik) = Sig_lesser(i,i,ie,ik) + G_lesser(i,i,ie+nop,ikd) * (n_bose+1.0_dp) * Dop**2 * dE            
                        enddo
                    endif
                    !
                    ! Sig^>(E,k)
                    if (gamma_q) then 
                        ikd = ik
                    else
                        ikd = map_kq_2d(-1,ik,iq,nky,nkz)
                    endif
                    if (ie-nop>=1) then 
                        do i=1,nm
                            Sig_greater(i,i,ie,ik) = Sig_greater(i,i,ie,ik) + G_greater(i,i,ie-nop,ikd) * (n_bose+1.0_dp) * Dop**2 * dE   
                        enddo
                    endif
                    if (gamma_q) then 
                        ikd = ik 
                    else                
                        ikd = map_kq_2d(+1,ik,iq,nky,nkz)
                    endif
                    if (ie+nop<=nen) then 
                        do i=1,nm
                            Sig_greater(i,i,ie,ik) = Sig_greater(i,i,ie,ik) + G_greater(i,i,ie+nop,ikd) * n_bose * Dop**2 * dE   
                        enddo
                    endif
                enddo
			enddo  
        enddo
        !$omp end do   
        !$omp end parallel
    end subroutine selfenergy_eph_simple

    ! calculate e-photon/phonon self-energies for single mode in thermal equilibrium 
    subroutine selfenergy_eph_mono(nm,nen,En,nop,nky,nkz,nqy,nqz,ik_start,ik_end,iq_in,M,G_lesser,G_greater,&
        Sig_lesser,Sig_greater,n_bose,gamma_q)
        integer,intent(in)::nm,nen,nky,nkz,nqy,nqz,iq_in,ik_start,ik_end
        integer,intent(in)::nop ! phonon freq. in unit of energy discretization step
        real(8),intent(in)::en(nen),n_bose 
        logical,intent(in)::gamma_q ! whether q=0
        complex(8),intent(in),dimension(nm,nm,nky*nkz,nqy*nqz)::M ! list of the interaction matrix at all k q pairs
        complex(8),intent(in),dimension(nm,nm,nen,nky*nkz)::G_lesser,G_greater ! Green's functions
        complex(8),intent(inout),dimension(nm,nm,nen,nky*nkz)::Sig_lesser,Sig_greater ! accumulate the e-ph self-energies 
        !---------
        integer::ie,ikd,ik ,iq,iqd
        complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix  
        real(8)::dE 
        Sig_lesser = czero
        Sig_greater = czero
        dE = (en(2)-en(1)) / twopi             
        if (gamma_q) then 
			iq=1
			iqd=1
		else
			iq=iq_in
			iqd=nqy*nqz+1-iq
		endif
        ! Sig^<>(E,k) = M_{k,-+q} [ N G^<>(E -+ hw,k-+q) + (N+1) G^<>(E +- hw,k+-q)] M_{k-+q,+-q}       
        !$omp parallel default(shared) private(ie,A,B,ik,ikd) 
        allocate(B(nm,nm))
        allocate(A(nm,nm))                                
        !$omp do
        do ie=1,nen
			do ik=ik_start,ik_end
				! Sig^<(E,k)
				A = czero            
				if (gamma_q) then 
					ikd = ik					
				else
					ikd = map_kq_2d(-1,ik,iq,nky,nkz)					
				endif
				if (ie-nop>=1) A =A+ G_lesser(:,:,ie-nop,ikd) * n_bose
				if (gamma_q) then 
					ikd = ik
				else
					ikd = map_kq_2d(+1,ik,iq,nky,nkz)
				endif
				if (ie+nop<=nen) A =A+ G_lesser(:,:,ie+nop,ikd) * (n_bose+1.0_dp)
				call zgemm('n','n',nm,nm,nm,cone,M(:,:,ik,iq),nm,A,nm,czero,B,nm) 
				call zgemm('n','n',nm,nm,nm,cone,B,nm,M(:,:,ikd,iqd),nm,czero,A,nm)     
				Sig_lesser(:,:,ie,ik) = Sig_lesser(:,:,ie,ik) + A * dE            
				!
				! Sig^>(E,k)
				A = czero
				if (gamma_q) then 
					ikd = ik
				else
					ikd = map_kq_2d(-1,ik,iq,nky,nkz)
				endif
				if (ie-nop>=1) A =A+ G_greater(:,:,ie-nop,ikd) * (n_bose+1.0_dp)
				if (gamma_q) then 
					ikd = ik 
				else                
					ikd = map_kq_2d(+1,ik,iq,nky,nkz)
				endif
				if (ie+nop<=nen) A =A+ G_greater(:,:,ie+nop,ikd) * n_bose
				call zgemm('n','n',nm,nm,nm,cone,M(:,:,ik,iq),nm,A,nm,czero,B,nm) 
				call zgemm('n','n',nm,nm,nm,cone,B,nm,M(:,:,ikd,iqd),nm,czero,A,nm)     
				Sig_greater(:,:,ie,ik) = Sig_greater(:,:,ie,ik) + A * dE               
			enddo  
        enddo
        !$omp end do        
        deallocate(A,B)
        !$omp end parallel
    end subroutine selfenergy_eph_mono

end module eph_dense
