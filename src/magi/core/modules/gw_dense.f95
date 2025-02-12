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
module gw_dense 
    use parameters_mod
    use output
    use legendre
    use observ
    use open_boundary
    use omp_lib
    implicit none 

    contains

    ! driver for iterating G -> P -> W -> Sig 
    ! memory saving version of green_solve_gw_1D 
    !   the full matrix P and W are needed for only one energy in this implementation,
    !   they are computed per energy point, and the contribution to selfenergy is 
    !   immediately added to sigma_x_new matrices.
    subroutine solve_gw_1D_memsaving(niter,nm_dev,Lx,length,spindeg,temp,mu,&
        alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,V,&
        ndiag,encut,Egap,vertex,bse,flatband,output_files,&
        G_retarded,G_lesser,G_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new,&
        current,transmission,W0_retarded,W0_lesser,W0_greater)  
        integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
        integer, intent(in) :: ndiag
        real(dp), intent(in) :: En(nen), temp(2), mu(2), alpha_mix, Lx, spindeg, Egap
        complex(dp),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
        complex(dp), intent(in):: V(nm_dev,nm_dev)    
        real(dp),intent(in)::encut(2) ! intraband and interband cutoff for P
        real(dp),intent(out)::current(nen,2) ! current spectrum on leads
        real(dp),intent(out)::transmission(nen,2) ! transmission matrix
        logical, intent(in) :: vertex, bse, flatband, output_files
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nen) ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
        complex(dp),intent(out),dimension(nm_dev,nm_dev) ::  W0_retarded,W0_lesser,W0_greater        
        !----
        complex(dp),dimension(nm_dev,nm_dev,nen) ::  Sig_retarded,Sig_lesser,Sig_greater
        complex(dp),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
        complex(dp),allocatable,dimension(:,:):: B ! tmp matrix
        real(dp),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:)
        real(dp),allocatable::wen(:) ! energy vector for P and W
        integer,allocatable::nops(:) ! discretized energy for P and W        
        real(dp),allocatable::Te(:,:,:),tr(:,:) ! transmission matrix spectrum 
        integer :: iter,ie,iop,nnop,nnop1,nnop2
        integer :: i,j,nm,l,h,nop
        complex(dp),allocatable::Ispec(:,:,:),Itot(:,:)
        complex(dp),allocatable,dimension(:,:) ::  P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater        
        real(dp) :: start,finish, time_P, time_W, time_sigma
        complex(dp) :: dE, epsilon
        real(dp)::nelec(2),pelec(2)
        !
        print *,'================= green_solve_gw_1D_memsaving ================='        
        print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
        print '(a8,f15.4,a8,f15.4)', 'T_s=',temp(1),'T_d=',temp(2)
        ! build the energy vector for P and W
        dE= En(2)-En(1) 
        nnop1=floor(min(encut(1),Egap)/dble(dE)) ! intraband exclude encut(1), include 0 
        nnop2=floor((min(encut(2),(maxval(En)-minval(En))) - Egap)/dble(dE))  ! interband , include Egap
        nnop=nnop1*2-1+nnop2*2 ! + and - freq.
        allocate(nops(nnop))
        allocate(wen(nnop))
        do iop=1,nnop1*2-1
            nops(iop+nnop2) = iop - nnop1    
        enddo
        do iop=1,nnop2
            nops(nnop2+1-iop) = -iop+1 - floor(Egap/dble(dE))
            nops(nnop2+nnop1*2-1+iop) = -nops(nnop2+1-iop)
        enddo
        wen(:) = dble(nops(:))*dble(dE)
        print *,'---------------------------------------------------------------'
        print *, ' Energy cutoff: intra-band    inter-band    Eg (eV)' 
        print '(A12,3F14.4)',' ',encut,egap
        ! print *, ' Nop='
        ! print '(10I5)',nops
        ! print *, ' Eop= (eV)'
        ! print '(6F8.3)',wen
        print *,'---------------------------------------------------------------'
        !
        allocate(siglead(NB*NS,NB*NS,nen,2))
        ! get leads sigma
        siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
        siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)  
        allocate(B(nm_dev,nm_dev))
        allocate(tot_cur(nm_dev,nm_dev))
        allocate(tot_ecur(nm_dev,nm_dev))
        allocate(cur(nm_dev,nm_dev,nen))
        allocate(Ispec(nm_dev,nm_dev,nen))
        allocate(Itot(nm_dev,nm_dev))
        allocate(te(nen,2,2))
        allocate(tr(nen,2))
        !
        allocate(P_lesser(nm_dev,nm_dev))
        allocate(P_greater(nm_dev,nm_dev))
        allocate(P_retarded(nm_dev,nm_dev)) 
        allocate(W_lesser(nm_dev,nm_dev))
        allocate(W_greater(nm_dev,nm_dev))
        allocate(W_retarded(nm_dev,nm_dev)) 
        !
        do iter=0,niter
            print *,'+ iter=',iter  
            print *, 'calc G'  
            start = omp_get_wtime()
            call calc_gf(nen,En,2,nm_dev,[nb*ns,nb*ns],nb*ns,&
                            Ham(:,:),H00lead(:,:,:),H10lead(:,:,:),Siglead(:,:,:,:),&
                            T(:,:,:),Sig_retarded(:,:,:),Sig_lesser(:,:,:),Sig_greater(:,:,:),&
                            G_retarded(:,:,:),G_lesser(:,:,:),&
                            G_greater(:,:,:),Tr,Te,mu,temp,flatband)                        
            finish = omp_get_wtime()
            print '("  G computation time = ", F0.3 ," seconds.")', finish-start
            start = finish                    
            if (output_files) then
                call write_current_spectrum('gw_Jdens',iter,cur,nen,en,length,NB,Lx)
                call write_current('gw_I',iter,tot_cur,length,NB,NS,Lx)
                call write_current('gw_EI',iter,tot_ecur,length,NB,NS,Lx)
                call write_spectrum_nosub('gw_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
                call write_spectrum_nosub('gw_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
                call write_spectrum_nosub('gw_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))
                call write_transmission_spectrum('gw_trL',iter,Tr(:,1)*spindeg,nen,En)
                call write_transmission_spectrum('gw_trR',iter,Tr(:,2)*spindeg,nen,En)
                call write_transmission_spectrum('gw_TE_LR',iter,Te(:,1,2)*spindeg,nen,En)
                call write_transmission_spectrum('gw_TE_RL',iter,Te(:,2,1)*spindeg,nen,En)    
            endif
            !        
            ! empty sigma_x_new matrices for accumulation
            sig_retarded_new=czero
            sig_lesser_new=czero
            sig_greater_new=czero
            print *, 'calc P, solve W, add to Sigma_new'     
            print *,'ndiag=',min(ndiag,nm_dev)
            !
            !print *,'   i / n :  Nop   Eop (eV)'
            time_P = 0.0_dp
            time_W = 0.0_dp
            time_sigma = 0.0_dp
            !
            do iop=1,nnop              
            !print '(I5,A,I5,A,I5,F8.3)',iop,'/',nnop,':',nops(iop),wen(iop)    
            nop=nops(iop)
            P_lesser = czero
            P_greater = czero
            P_retarded = czero
            start = omp_get_wtime()
            if ((vertex).and.(iter>0)) then
                if (iop==1) print *, '  vertex on '
                call calc_P_vertex_correction(vertex,nm_dev,nen,nop,ndiag,(en(2)-en(1)),&
                        G_retarded,G_lesser,G_greater,W0_retarded,W0_lesser,W0_greater,&
                        P_retarded,P_lesser,P_greater)
            else
                call calc_P_vertex_correction(.false.,nm_dev,nen,nop,ndiag,(en(2)-en(1)),&
                        G_retarded,G_lesser,G_greater,W0_retarded,W0_lesser,W0_greater,&
                        P_retarded,P_lesser,P_greater)
            endif
            finish = omp_get_wtime()
            time_P = time_P + finish-start
            !             
            dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )* spindeg    
            P_lesser=P_lesser*dE
            P_greater=P_greater*dE  
            ! P_retarded=P_retarded*dE      
            P_retarded=dcmplx(0.0_dp*dble(P_retarded), 0.5_dp*aimag(P_greater-P_lesser))
            !         
            ! calculate W
            start = omp_get_wtime()
            !
            call calc_w(1,NB,NS,nm_dev,P_retarded,P_lesser,P_greater,V,W_retarded,W_lesser,W_greater)
            !
            finish = omp_get_wtime()
            time_W = time_W + finish-start               
            !        
            !
            if (iop == (nnop1+nnop2)) then
                ! store static W for the vertex correction in the next iteration
                W0_retarded = W_retarded
                W0_lesser = W_lesser
                W0_greater = W_greater                
            endif          
            !
            start = omp_get_wtime()
            ! Accumulate the GW to Sigma
            ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)  
            !$omp parallel default(shared) private(l,h,i,ie) 
            !$omp do
            do i=1,nm_dev
                l=max(i-ndiag,1)
                h=min(nm_dev,i+ndiag)           
                do ie=1,nen
                if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then 
                    Sig_lesser_new(i,l:h,ie)=Sig_lesser_new(i,l:h,ie)+G_lesser(i,l:h,ie-nop)*W_lesser(i,l:h)
                    Sig_greater_new(i,l:h,ie)=Sig_greater_new(i,l:h,ie)+G_greater(i,l:h,ie-nop)*W_greater(i,l:h)
                    Sig_retarded_new(i,l:h,ie)=Sig_retarded_new(i,l:h,ie)+G_lesser(i,l:h,ie-nop)*W_retarded(i,l:h) + &                                      
                                            G_retarded(i,l:h,ie-nop)*W_lesser(i,l:h) + &
                                            G_retarded(i,l:h,ie-nop)*W_retarded(i,l:h)                                               
                endif     
                enddo   
            enddo
            !$omp end do
            !$omp end parallel    
            finish = omp_get_wtime()
            time_sigma = time_sigma + finish - start      
            enddo  
            close(199)    
            print '("  P computation time = ", F0.3 ," seconds.")', time_P
            print '("  W computation time = ", F0.3 ," seconds.")', time_W
            print '("  Sigma computation time = ", F0.3 ," seconds.")', time_sigma
            !
            dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)                
            Sig_lesser_new = Sig_lesser_new  * dE
            Sig_greater_new= Sig_greater_new * dE
            Sig_retarded_new=Sig_retarded_new* dE
            !
            Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
            ! symmetrize the selfenergies
            do ie=1,nen
              B(:,:)=transpose(Sig_retarded_new(:,:,ie))
              Sig_retarded_new(:,:,ie) = (Sig_retarded_new(:,:,ie) + B(:,:))/2.0d0
              B(:,:)=transpose(Sig_lesser_new(:,:,ie))
              Sig_lesser_new(:,:,ie) = dcmplx((dble(Sig_lesser_new(:,:,ie)) - dble(B(:,:)))/2.0d0, (dimag(Sig_lesser_new(:,:,ie)) + dimag(B(:,:)))/2.0d0)
              B(:,:)=transpose(Sig_greater_new(:,:,ie))
              Sig_greater_new(:,:,ie) = dcmplx((dble(Sig_greater_new(:,:,ie)) - dble(B(:,:)))/2.0d0, (dimag(Sig_greater_new(:,:,ie)) + dimag(B(:,:)))/2.0d0)
            enddo
            !        
            ! mixing with the previous one
            Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
            Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
            Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)            
            ! get leads sigma
            siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
            siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)    
            !
            if (output_files) then
                call write_spectrum_nosub('gw_SigR',iter,Sig_retarded,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
                call write_spectrum_nosub('gw_SigL',iter,Sig_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
                call write_spectrum_nosub('gw_SigG',iter,Sig_greater,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
                !call write_matrix_summed_overE('Sigma_r',iter,Sig_retarded,nen,en,length,NB,(/1.0,1.0/))
                !!!! calculate collision integral
                ! call calc_collision(Sig_lesser_new,Sig_greater_new,G_lesser,G_greater,nen,en,spindeg,nm_dev,Itot,Ispec)
                ! call write_spectrum_nosub('gw_Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
            endif
        enddo             
        if (output_files) then 
            open(unit=11,file='WR0.dat',status='unknown')
            do i=1, nm_dev
                do j=1, nm_dev
                    write(11,'(2I6,2E15.4)') i,j, dble(W0_retarded(i,j)), aimag(W0_retarded(i,j))
                end do
                write(11,*)
            end do
            close(11)
        endif
        !! last step  
        print *, 'calc G last time ...'  
        !
        call calc_gf(nen,En,2,nm_dev,[nb*ns,nb*ns],nb*ns,&
                            Ham(:,:),H00lead(:,:,:),H10lead(:,:,:),Siglead(:,:,:,:),&
                            T(:,:,:),Sig_retarded(:,:,:),Sig_lesser(:,:,:),Sig_greater(:,:,:),&
                            G_retarded(:,:,:),G_lesser(:,:,:),&
                            G_greater(:,:,:),Tr,Te,mu,temp,flatband)                        
        !   
        current=tr
        transmission(:,1)=te(:,1,2)
        transmission(:,2)=te(:,2,1)     
        if (output_files) then
            call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
            call write_current_spectrum('gw_Jdens',iter,cur,nen,en,length,NB,Lx)
            call write_current('gw_I',iter,tot_cur,length,NB,NS,Lx)
            call write_current('gw_EI',iter,tot_ecur,length,NB,NS,Lx)
            call write_spectrum_nosub('gw_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
            call write_spectrum_nosub('gw_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
            call write_spectrum_nosub('gw_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))                    
        endif
        !
    end subroutine solve_gw_1D_memsaving
    !
    !
    ! calculate polarization with 1st order vertex correction 
    subroutine calc_P_vertex_correction(lvertex,nm_dev,nen,nop,ndiag,dE,G_retarded,G_lesser,G_greater,W_retarded,W_lesser,W_greater,P_retarded,P_lesser,P_greater)
        integer, intent(in) :: nm_dev, nen, nop, ndiag
        real(dp), intent(in) :: dE
        logical, intent(in) :: lvertex
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater ! electron GF
        complex(dp),intent(in),dimension(nm_dev,nm_dev) ::  W_retarded,W_lesser,W_greater ! W_0 static screened Coulomb interaction
        complex(dp),intent(out),dimension(nm_dev,nm_dev) ::  P_retarded,P_lesser,P_greater ! 1st order vertex-corrected polarization
        ! ---------- local variables
        integer :: i,j,k, ie, n, ie1,ie2
        real(dp)::alpha
        alpha = 0.0_dp
        ie1=max(nop+1,1)
        ie2=min(nen,nen+nop) 
        ! P0_ijk = -i G_ik * G_kj
        ! P1_mnk = P0_mnk + i \sum_ij P0_ijk W0_ij G_mi * G_jn
        ! if we only need after-all P1_nnk, we just plug P0 into P1 and get
        ! P1_nnk = -i G_nk * G_kn + i \sum_ij (-i G_ik * G_kj) W0_ij G_ni * G_jn
        ! the range of n,k,i,j should be smaller than ndiag
        !$omp parallel default(shared) private(k,n,i,j,ie) 
        !$omp do
        do n = 1, nm_dev        
            do k = 1, nm_dev         
                if ((abs(n-k) <= ndiag)) then
                ! RPA polarization P0      
                do ie = max(nop+1,1),min(nen,nen+nop)                        
                    P_lesser(n,k) = P_lesser(n,k) + G_lesser(n,k,ie) * G_greater(k,n,ie-nop)
                    P_greater(n,k) = P_greater(n,k) + G_greater(n,k,ie) * G_lesser(k,n,ie-nop) 
            !        P_retarded(n,k) = P_retarded(n,k) + G_lesser(n,k,ie) * conjg(G_retarded(n,k,ie-nop)) + &
            !                          G_retarded(n,k,ie) * G_lesser(k,n,ie-nop)        
                enddo
                !
                if (lvertex) then
                    ! 1st order vertex correction    
                    do i = max(1,k-ndiag), min(nm_dev,k+ndiag)        
                        do j = max(1,k-ndiag), min(nm_dev,k+ndiag)          
                            if ((abs(i-k) <= ndiag).and.(abs(j-k) <= ndiag)) then                                                   
                                P_lesser(n,k) = P_lesser(n,k) + c1i* sum(G_lesser(i,k,ie1:ie2) * G_greater(k,j,ie1-nop:ie2-nop)) &
                                                                * W_lesser(i,j) * sum(G_lesser(n,i,ie1:ie2) * G_greater(j,n,ie1-nop:ie2-nop)) 
                                P_greater(n,k) = P_greater(n,k) + c1i* sum(G_greater(i,k,ie1:ie2) * G_lesser(k,j,ie1-nop:ie2-nop)) &
                                                                * W_greater(i,j) * sum(G_greater(n,i,ie1:ie2) * G_lesser(j,n,ie1-nop:ie2-nop)) 
                                ! P_retarded(n,k) = P_retarded(n,k) + &
                                !   c1i* sum(G_lesser(i,k,ie1:ie2) * conjg(G_retarded(j,k,ie1-nop:ie2-nop))) &
                                !   * W_retarded(i,j) * sum(G_lesser(n,i,ie1:ie2) * conjg(G_retarded(n,j,ie1-nop:ie2-nop))) +&
                                !   c1i* sum(G_retarded(i,k,ie1:ie2) * G_lesser(k,j,ie1-nop:ie2-nop)) &
                                !   * W_retarded(i,j) * sum(G_retarded(n,i,ie1:ie2) * G_lesser(n,j,ie1-nop:ie2-nop))                                  
                            endif
                        enddo
                    enddo
                endif
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        ! P_retarded(:,:) = alpha * P_retarded(:,:) + (1.0_dp - alpha) * 0.5_dp * ( P_greater(:,:) - P_lesser(:,:) )
    end subroutine calc_P_vertex_correction
    !
    !
    subroutine selfenergy_gw(nm_dev,nen,nsub,nphiy,nphiz,nb,ns,ndiag,length, en, V,flatband,spindeg,sum_method,&
                             G_retarded,G_lesser,G_greater, &
                             Sig_retarded_new,Sig_lesser_new,Sig_greater_new,W0_retarded)
        !
        use fft_mod, only : conv1d_fock, corr1d => corr1d2  
        use parameters_mod
        !
        integer,intent(in):: nm_dev,nen,nsub,nphiy,nphiz,nb,length,ndiag,ns
        real(dp),dimension(nen),intent(in)::en
        real(dp),intent(in)::spindeg
        logical,intent(in) :: flatband
        character(len=*),intent(in)::sum_method
        complex(dp), intent(in):: V(nm_dev,nm_dev,nphiy*nphiz)
        complex(dp),dimension(nm_dev,nm_dev,nen,nsub,nphiy*nphiz),intent(in) ::  G_retarded,G_lesser,G_greater
        complex(dp),dimension(nm_dev,nm_dev,nen,nphiy*nphiz),intent(out) ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nphiy*nphiz) ::  W0_retarded
        ! -------
        complex(dp),dimension(:,:,:,:),allocatable ::  P_retarded,P_lesser,P_greater
        complex(dp),dimension(:,:,:,:),allocatable ::  W_retarded,W_lesser,W_greater        
        complex(dp),dimension(nen) :: tmp
        complex(dp) :: dE        
        real(dp)::weights(nsub),xen(nsub)        
        integer::nw ! number of >=0 energy points in P and W
        integer::ie,ik,ikd,iq,i,j,l,h,isub,nop,nopmax
        !
        dE = En(2) - En(1)
        call gaulegf(0.0d0, dble(dE), xen, weights, nsub) ! obtain the Legendre ordinates and weights    
        nopmax = nen/2 - 1 ! max nop in P and W
        nw = nopmax + 1 ! to include 0 and positive frequencies
        ! G -> P
        print *, '  calc P'
        allocate(P_retarded(nm_dev,nm_dev,nw,nphiy*nphiz),source=czero)
        allocate(P_lesser(nm_dev,nm_dev,nw,nphiy*nphiz),source=czero)
        allocate(P_greater(nm_dev,nm_dev,nw,nphiy*nphiz),source=czero)
        !                        
        ! Pij^<>(hw,kz') = \int_dE Gij^<>(E,kz) * Gji^><(E-hw,kz-kz')
        ! Pij^r(hw,kz')  = \int_dE Gij^<(E,kz) * Gji^a(E-hw,kz-kz') + Gij^r(E,kz) * Gji^<(E-hw,kz-kz')        
        do iq=1,nphiy*nphiz     
            do ik=1,nphiy*nphiz                                                
                ikd = map_kq_2d(-1,ik,iq,nphiy,nphiz)                    
                !$omp parallel default(shared) private(l,h,i,j,isub,tmp) 
                !$omp do        
                do i = 1, nm_dev        
                    l=max(i-ndiag,1)
                    h=min(nm_dev,i+ndiag)                        
                    do j = l,h
                        do isub = 1,nsub
                            tmp=corr1d(nen,G_lesser(i,j,:,isub,ik),G_greater(j,i,:,isub,ikd),method=sum_method) * weights(isub)
                            P_lesser(i,j,:,iq) = P_lesser(i,j,:,iq)   + tmp(nen/2:nen/2+nopmax)
                            tmp=corr1d(nen,G_greater(i,j,:,isub,ik),G_lesser(j,i,:,isub,ikd),method=sum_method) * weights(isub)         
                            P_greater(i,j,:,iq) = P_greater(i,j,:,iq) + tmp(nen/2:nen/2+nopmax)
                            tmp= corr1d(nen,G_lesser(i,j,:,isub,ik),conjg(G_retarded(i,j,:,isub,ikd)),method=sum_method) * weights(isub) &
                               + corr1d(nen,G_retarded(i,j,:,isub,ik),G_lesser(j,i,:,isub,ikd),method=sum_method) * weights(isub) 
                            P_retarded(i,j,:,iq) = P_retarded(i,j,:,iq) + tmp(nen/2:nen/2+nopmax)
                        enddo
                    enddo
                enddo
                !$omp end do
                !$omp end parallel                
            enddo                    
        enddo                
        dE = dcmplx(0.0d0 , -1.0d0 / 2.0d0 / pi ) * spindeg /dble(nphiy)/dble(nphiz) 
        P_lesser=dE*P_lesser
        P_greater=dE*P_greater
        P_retarded=dE*P_retarded
        ! P -> W
        print *, '  calc W'
        allocate(W_retarded(nm_dev,nm_dev,nw,nphiy*nphiz),source = czero)
        allocate(W_lesser(nm_dev,nm_dev,nw,nphiy*nphiz),source = czero)
        allocate(W_greater(nm_dev,nm_dev,nw,nphiy*nphiz),source = czero)        
        do iq=1,nphiy*nphiz        
            !$omp parallel default(shared) private(nop)
            !$omp do
            do nop=1,nw
                if (flatband) then
                    call calc_w(0,NB,NS,nm_dev,P_retarded(:,:,nop,iq),P_lesser(:,:,nop,iq),P_greater(:,:,nop,iq),&
                                V(:,:,iq),W_retarded(:,:,nop,iq),W_lesser(:,:,nop,iq),W_greater(:,:,nop,iq))
                else      
                    call calc_w(1,NB,NS,nm_dev,P_retarded(:,:,nop,iq),P_lesser(:,:,nop,iq),P_greater(:,:,nop,iq),&
                                V(:,:,iq),W_retarded(:,:,nop,iq),W_lesser(:,:,nop,iq),W_greater(:,:,nop,iq))
                endif
            enddo
            !$omp end do
            !$omp end parallel
        enddo  
        deallocate(P_retarded,P_lesser,P_greater)
        ! W -> Sigma
        print *, '  calc Sig'
        !   hw should go from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)   
        !   but W<>(-hw) = - conjg( W><(hw) )    
        !   so Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw) + G^<>_ij(E+hw) W^><_ij(hw)  
        Sig_greater_new = dcmplx(0.0d0,0.0d0)
        Sig_lesser_new = dcmplx(0.0d0,0.0d0)
        Sig_retarded_new = dcmplx(0.0d0,0.0d0)     
        do ik=1,nphiy*nphiz    
            do iq=1,nphiy*nphiz    
                ikd = map_kq_2d(-1,ik,iq,nphiy,nphiz)                 
                !$omp parallel default(shared) private(l,h,i,j,isub)
                !$omp do  
                do i = 1,nm_dev   
                    l=max(i-ndiag,1)
                    h=min(nm_dev,i+ndiag)       
                    do j = l,h
                        do isub = 1,nsub
                            Sig_lesser_new(i,j,:,ik) = Sig_lesser_new(i,j,:,ik) &
                                + conv1d_fock(nen,nw,G_lesser(i,j,:,isub,ikd),W_lesser(i,j,:,iq),W_greater(i,j,:,iq),method=sum_method) * weights(isub)                                                             
                            Sig_greater_new(i,j,:,ik) = Sig_greater_new(i,j,:,ik) &
                                + conv1d_fock(nen,nw,G_greater(i,j,:,isub,ikd),W_greater(i,j,:,iq),W_lesser(i,j,:,iq),method=sum_method) * weights(isub) 
                                
                           Sig_retarded_new(i,j,:,ik) = Sig_retarded_new(i,j,:,ik) &
                               + conv1d_fock(nen,nw,G_lesser(i,j,:,isub,ikd),W_retarded(i,j,:,iq),W_retarded(i,j,:,iq),method=sum_method) * weights(isub) &
                               + conv1d_fock(nen,nw,G_retarded(i,j,:,isub,ikd),W_lesser(i,j,:,iq),W_greater(i,j,:,iq),method=sum_method) * weights(isub) &
                               + conv1d_fock(nen,nw,G_retarded(i,j,:,isub,ikd),W_retarded(i,j,:,iq),W_retarded(i,j,:,iq),method=sum_method) * weights(isub)                                               
                        enddo 
                    enddo
                enddo      
                !$omp end do
                !$omp end parallel
            enddo            
        enddo        
        dE = dcmplx(0.0d0, 1.0d0/twopi) /dble(nphiy)/dble(nphiz)
        Sig_lesser_new = Sig_lesser_new  * dE
        Sig_greater_new= Sig_greater_new * dE
        Sig_retarded_new=Sig_retarded_new* dE
        Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
        !    
        W0_retarded = W_retarded(:,:,1,:)
        deallocate(W_lesser,W_greater,W_retarded)
        !
    end subroutine selfenergy_gw


    ! 3D GW solver with two periodic directions (y,z)
    ! iterating G -> P -> W -> Sig 
    subroutine solve_gw_3D(niter,scba_tol,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
        alpha_mix,nen,nsub,En,nb,ns,nphiy,nphiz,Ham,H00lead,H10lead,T,V,&
        ndiag,num_lead,flatband,output_files,sum_method,G_retarded,G_lesser,G_greater,W0_retarded,tr)
        !
        ! use fft_mod, only : conv1d => conv1d2, corr1d => corr1d2  
        use parameters_mod
        !  
        integer, intent(in) :: nen, nsub, nb, ns,niter,nm_dev,length, nphiz, nphiy, num_lead
        real(dp), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg,scba_tol
        complex(dp),intent(in) :: Ham(nm_dev,nm_dev,nphiy*nphiz),H00lead(NB*NS,NB*NS,num_lead,nphiy*nphiz),H10lead(NB*NS,NB*NS,num_lead,nphiy*nphiz),T(NB*NS,nm_dev,num_lead,nphiy*nphiz)
        complex(dp), intent(in):: V(nm_dev,nm_dev,nphiy*nphiz)
        integer,intent(in)::ndiag
        logical,intent(in)::flatband
        logical,intent(in) :: output_files
        character(len=*),intent(in)::sum_method
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nen,nsub,nphiy*nphiz) ::  G_retarded,G_lesser,G_greater
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nphiy*nphiz) ::  W0_retarded
        real(dp),intent(out) ::Tr(nen,num_lead) ! current spectrum on leads    
        !------
        complex(dp),dimension(:,:,:,:),allocatable ::  P_retarded,P_lesser,P_greater
        complex(dp),dimension(:,:,:,:),allocatable ::  W_retarded,W_lesser,W_greater
        complex(dp),dimension(:,:,:,:),allocatable ::  Sig_retarded,Sig_lesser,Sig_greater
        complex(dp),dimension(:,:,:,:),allocatable ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
        complex(dp),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
        complex(dp),allocatable,dimension(:,:):: B ! tmp matrix
        real(dp),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:),sumcur(:,:,:),sumtot_cur(:,:),sumtot_ecur(:,:)
        complex(dp),allocatable::Ispec(:,:,:),Itot(:,:)    
        real(dp),allocatable::Te(:,:,:) ! transmission matrix spectrum
        real(dp),allocatable::sumTr(:,:) ! current spectrum on leads summed over k
        real(dp),allocatable::sumTe(:,:,:) ! transmission matrix spectrum summed over k
        integer :: iter,ie,nopmax
        integer :: i,j,nm,nop,l,h,iop,ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,isub        
        complex(dp) :: dE
        real(dp)::nelec(2),mu(2),pelec(2),temp(2)
        real(dp)::weights(nsub),xen(nsub)
        real(dp)::scba_error
        complex(dp),allocatable::Scat_spec(:,:,:,:) ! collision integral spectrum
        complex(dp),allocatable::Scat(:,:) ! collision integral
        
        allocate(Sig_retarded(nm_dev,nm_dev,nen,nphiy*nphiz),Sig_lesser(nm_dev,nm_dev,nen,nphiy*nphiz),Sig_greater(nm_dev,nm_dev,nen,nphiy*nphiz))
        
        scba_error=1.0d0
        Sig_retarded = czero
        sig_lesser = czero
        sig_greater = czero
        !
        print *,'============ green_solve_gw_3D ============'
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
            temp=(/temps,tempd/)
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
            sumcur=sumcur/dble(nphiy)/dble(nphiz)
            sumtot_cur=sumtot_cur/dble(nphiy)/dble(nphiz)
            sumtot_ecur=sumtot_ecur/dble(nphiy)/dble(nphiz)
            sumTr=sumTr/dble(nphiz)/dble(nphiy)
            sumTe=sumTe/dble(nphiz)/dble(nphiy)
            sumTr = sumTr *e_charge/twopi/hbar*e_charge*dble(spindeg)
            if (output_files) then
                if (flatband) then
                    print *,'flatband'
                    ! call write_spectrum_per_kz('gw_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
                    ! call write_spectrum_per_kz('gw_gamma-centered_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=0.0_dp,at_kz=0.0_dp)
                    call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
            
                else
                    call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
                    call write_spectrum_summed_over_kz('gw_ndos',iter,G_lesser(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
                    call write_spectrum_summed_over_kz('gw_pdos',iter,G_greater(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
                endif
                call write_current_spectrum('gw_Jdens',iter,sumcur,nen,en,length,NB,Lx)
                call write_current('gw_I',iter,sumtot_cur,length,NB,NS,Lx)
                call write_current('gw_EI',iter,sumtot_ecur,length,NB,NS,Lx)
                call write_transmission_spectrum('gw_trL',iter,sumTr(:,1),nen,En)
                call write_transmission_spectrum('gw_trR',iter,sumTr(:,2),nen,En)
                ! call write_transmission_spectrum('gw_TE_LR',iter,sumTe(:,1,2)*spindeg,nen,En)
                ! call write_transmission_spectrum('gw_TE_RL',iter,sumTe(:,2,1)*spindeg,nen,En)            
            endif            
            open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
            write(101,'(I4,2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
            close(101)
            write(*,'(I4,"  IDS=",2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
            !
            G_retarded=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
            G_lesser=dcmplx(0.0d0*dble(G_lesser),aimag(G_lesser))
            G_greater=dcmplx(0.0d0*dble(G_greater),aimag(G_greater))
            !                    
            allocate(Sig_retarded_new(nm_dev,nm_dev,nen,nphiy*nphiz),source=czero)
            allocate(Sig_lesser_new(nm_dev,nm_dev,nen,nphiy*nphiz),source=czero)
            allocate(Sig_greater_new(nm_dev,nm_dev,nen,nphiy*nphiz),source=czero)
            !          
            call selfenergy_gw(nm_dev,nen,nsub,nphiy,nphiz,nb,ns,ndiag,length, en, V,flatband,spindeg,sum_method,&
                             G_retarded,G_lesser,G_greater, &
                             Sig_retarded_new,Sig_lesser_new,Sig_greater_new,W0_retarded)
            !
            if (output_files) then
                print *,'  calc collision integral'
                allocate(Scat_spec(nm_dev,nm_dev,nen,nsub),source=czero)
                allocate(Scat(nm_dev,nm_dev),source=czero)
                call calc_collision(Sig_lesser_new,Sig_greater_new,G_lesser,G_greater,nen,en,nsub,nphiy*nphiz,spindeg,nm_dev,Scat,Scat_spec)
                call write_spectrum('gw_Scat',iter,Scat_spec,nen,nsub,En,xen,length,NB,Lx,(/1.0d0,1.0d0/))                                
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
            open(unit=101,file='gw_scba_error.dat',status='unknown',position='append')
            write(101,'(I4,E16.6)') iter, scba_error
            close(101)
            iter=iter+1
            ! mixing with previous ones
            Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
            Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
            Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
            !
            deallocate(Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
            !
            if (.not. flatband) then
                ! get leads sigma
                do iqz=1,nphiy*nphiz
                    siglead(:,:,:,1,iqz) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:,iqz)
                    siglead(:,:,:,2,iqz) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:,iqz)    
                enddo            
            endif
            if (flatband) then            
                ! call write_spectrum_per_kz('gw_SigR',iter,Sig_retarded,nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
            endif 
            ! call write_spectrum_summed_over_kz('gw_SigR',iter,Sig_retarded,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            !  call write_spectrum_summed_over_kz('SigL',iter,Sig_lesser,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
            !  call write_spectrum_summed_over_kz('SigG',iter,Sig_greater,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
        end do  
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
        if (flatband) then
            print *,'flatband'
            ! call write_spectrum_per_kz('gw_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
            ! call write_spectrum_per_kz('gw_gamma-centered_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=0.0_dp,at_kz=0.0_dp)
            call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))

        else
            call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
            call write_spectrum_summed_over_kz('gw_ndos',iter,G_lesser(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            call write_spectrum_summed_over_kz('gw_pdos',iter,G_greater(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
        endif
        call write_current_spectrum('gw_Jdens',iter,sumcur,nen,en,length,NB,Lx)
        call write_current('gw_I',iter,sumtot_cur,length,NB,NS,Lx)
        call write_current('gw_EI',iter,sumtot_ecur,length,NB,NS,Lx)
        call write_transmission_spectrum('gw_trL',iter,sumTr(:,1)*spindeg,nen,En)
        call write_transmission_spectrum('gw_trR',iter,sumTr(:,2)*spindeg,nen,En)
        ! call write_transmission_spectrum('gw_TE_LR',iter,sumTe(:,1,2)*spindeg,nen,En)
        ! call write_transmission_spectrum('gw_TE_RL',iter,sumTe(:,2,1)*spindeg,nen,En)
        sumTr = sumTr *e_charge/twopi/hbar*e_charge*dble(spindeg)
        open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
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
    end subroutine solve_gw_3D
  
    function map_kq(sgn,ik,iq,nk)
        integer,intent(in)::ik,iq,nk,sgn 
        integer::map_kq 
        integer::ikd
        ikd = ik + sgn * (iq - nk/2)            
        if (ikd<1) ikd=ikd+nk
        if (ikd>nk) ikd=ikd-nk
        if (nk==1) ikd=1
        map_kq=ikd 
    end function map_kq

    function map_kq_2d(sgn,ik,iq,nky,nkz)
        integer,intent(in)::ik,iq,nky,nkz,sgn 
        integer::map_kq_2d 
        integer::ikd,ikyd,ikzd,iqy,iqz,iky,ikz 
        !
        iqz = mod(iq-1,nkz)+1
        iqy = (iq-1) / nkz +1
        ikz = mod(ik-1,nkz)+1
        iky = (ik-1) / nkz +1
        ikzd= map_kq(sgn,ikz,iqz,nkz)
        ikyd= map_kq(sgn,iky,iqy,nky)        
        ikd = ikzd + (ikyd-1) * nkz
        map_kq_2d = ikd
    end function map_kq_2d

    Function ferm(a)
        Real(dp),intent(in):: a
        Real(dp):: ferm
        ferm=1.0d0/(1.0d0+Exp(a))
    End Function ferm

    ! calculate Gr and G<>
    subroutine calc_gf(ne,E,num_lead,nm_dev,nm_lead,max_nm_lead,Ham,lead_H00,lead_H10,Siglead,T,&
        Scat_Sig_retarded,Scat_Sig_lesser,Scat_Sig_greater,G_retarded,G_lesser,G_greater,&
        cur,te,mu,temp,flatband)
        use parameters_mod
        implicit none
        integer, intent(in) :: num_lead ! number of leads/contacts
        integer, intent(in) :: nm_dev   ! size of device Hamiltonian
        integer, intent(in) :: nm_lead(:) ! size of lead Hamiltonians
        integer, intent(in) :: max_nm_lead ! max size of lead Hamiltonians
        real(kind=dp), intent(in) :: E(:)  ! energy vector
        real(kind=dp), intent(out) :: cur(ne,num_lead)  ! current spectrum on leads
        real(kind=dp), intent(out) :: te(ne,num_lead,num_lead)  ! transmission matrix
        integer, intent(in) :: ne ! number of energy points
        complex(kind=dp), intent(in) :: Ham(nm_dev,nm_dev) ! system Hamiltonian
        complex(kind=dp), intent(in) :: lead_H00(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian diagonal blocks
        complex(kind=dp), intent(in) :: lead_H10(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian off-diagonal blocks
        complex(kind=dp), intent(in) :: Siglead(max_nm_lead,max_nm_lead,ne,num_lead) ! lead sigma_r scattering
        complex(kind=dp), intent(in) :: T(max_nm_lead,nm_dev,num_lead)  ! coupling matrix between leads and device
        complex(kind=dp), intent(in) :: Scat_Sig_retarded(nm_dev,nm_dev,ne) ! scattering Selfenergy
        complex(kind=dp), intent(in) :: Scat_Sig_lesser(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(in) :: Scat_Sig_greater(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(out) :: G_retarded(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(out) :: G_lesser(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(out) :: G_greater(nm_dev,nm_dev,ne)
        logical,intent(in)::flatband
        real(kind=dp), intent(in) :: mu(num_lead), temp(num_lead)
        ! ----
        integer :: i,j,nm,ie,io
        complex(kind=dp), allocatable, dimension(:,:) :: S00,G00,GBB,A,sig,sig_lesser,sig_greater,B,C,Hii,H1i
        complex(kind=dp), allocatable, dimension(:,:,:) :: gamma_lead
        real(kind=dp) :: fd
        complex(kind=dp):: z     
        !
        cur=0.0d0
        te=0.0d0
        !
        !$omp parallel default(shared) private(z,sig,ie,sig_lesser,sig_greater,B,C,i,nm,S00,G00,GBB,A,fd,gamma_lead,Hii,H1i) 
        allocate(sig(nm_dev,nm_dev))  
        allocate(gamma_lead(nm_dev,nm_dev,num_lead))  
        allocate(sig_lesser(nm_dev,nm_dev))
        allocate(sig_greater(nm_dev,nm_dev))          
        allocate(B(nm_dev,nm_dev))
        allocate(C(nm_dev,nm_dev))
        !$omp do        
        do ie = 1, ne
            z=dcmplx(E(ie),0.0d0)
            G_retarded(:,:,ie) = - Ham(:,:) - Scat_Sig_retarded(:,:,ie)             
            sig_lesser(:,:) = dcmplx(0.0d0,0.0d0)      
            sig_greater(:,:) = dcmplx(0.0d0,0.0d0)      
            ! compute and add contact self-energies    
            !  open(unit=101,file='sancho_gbb.dat',status='unknown',position='append')
            !  open(unit=102,file='sancho_g00.dat',status='unknown',position='append')
            !  open(unit=103,file='sancho_sig.dat',status='unknown',position='append')
            do i = 1,num_lead
                nm = nm_lead(i)    
                allocate(Hii(nm,nm))
                allocate(H1i(nm,nm))
                allocate(S00(nm,nm))
                allocate(G00(nm,nm))
                allocate(GBB(nm,nm))
                allocate(A(nm_dev,nm))    
                call identity(S00,nm)            
                !
                if (flatband) then
                    G00 = -S00*c1i
                else        
                    Hii = lead_H00(1:nm,1:nm,i) + siglead(1:nm,1:nm,ie,i)
                    H1i = lead_H10(1:nm,1:nm,i)
                    call sancho(NM,E(ie),S00,Hii,H1i,G00,GBB)
                endif
                !
                call zgemm('c','n',nm_dev,nm,nm,cone,T(1:nm,1:nm_dev,i),nm,G00,nm,czero,A,nm_dev) 
                call zgemm('n','n',nm_dev,nm_dev,nm,cone,A,nm_dev,T(1:nm,1:nm_dev,i),nm,czero,sig,nm_dev)  
                !    write(101,'(i4,2E15.4)') i, E(ie), -aimag(trace(GBB,nm))*2.0d0
                !    write(102,'(i4,2E15.4)') i, E(ie), -aimag(trace(G00,nm))*2.0d0
                !    write(103,'(i4,2E15.4)') i, E(ie), -aimag(trace(sig,nm_dev))*2.0d0
                G_retarded(:,:,ie) = G_retarded(:,:,ie) - sig(:,:)
                !
                fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))
                B(:,:) = conjg(sig(:,:))
                C(:,:) = transpose(B(:,:))
                B(:,:) = sig(:,:) - C(:,:)
                sig_lesser(:,:) = sig_lesser(:,:) - B(:,:)*fd
                sig_greater(:,:) = sig_greater(:,:) + B(:,:)*(1.0d0-fd)            
                gamma_lead(:,:,i)= B(:,:)             
                deallocate(S00,G00,GBB,A,Hii,H1i)
            end do  
            !  close(101)
            !  close(102)
            !  close(103)
            do i = 1,nm_dev
                G_retarded(i,i,ie) = G_retarded(i,i,ie) + z 
            end do
            !
            call invert_inplace(G_retarded(:,:,ie),nm_dev) 
            !
            sig_lesser = sig_lesser + Scat_Sig_lesser(:,:,ie)
            sig_greater = sig_greater + Scat_Sig_greater(:,:,ie)     
            !
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_lesser,nm_dev,czero,B,nm_dev) 
            call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
            G_lesser(:,:,ie) = C   
            !
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_greater,nm_dev,czero,B,nm_dev) 
            call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
            G_greater(:,:,ie) = C         
            !
            ! calculate current spec and/or transmission at each lead/contact
            do i=1,num_lead                                  
                fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))
                call zgemm('n','n',nm_dev,nm_dev,nm_dev,dcmplx(1.0d0-fd,0.0d0),&
                    gamma_lead(:,:,i),nm_dev,G_lesser(:,:,ie),nm_dev,czero,B,nm_dev)
                call zgemm('n','n',nm_dev,nm_dev,nm_dev,dcmplx(fd,0.0d0),&
                    gamma_lead(:,:,i),nm_dev,G_greater(:,:,ie),nm_dev,cone,B,nm_dev)
                do io=1,nm_dev
                    cur(ie,i)=cur(ie,i)+ dble(B(io,io))
                enddo            
                do j=1,num_lead                      
                    if (j.ne.i) then
                        call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,gamma_lead(:,:,i),nm_dev,G_retarded(:,:,ie),nm_dev,czero,B,nm_dev)
                        call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,gamma_lead(:,:,j),nm_dev,czero,C,nm_dev)
                        call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,C,nm_dev,G_retarded(:,:,ie),nm_dev,czero,B,nm_dev)
                        do io=1,nm_dev
                            te(ie,i,j)=te(ie,i,j)- dble(B(io,io)) ! Gamma = i[Sig^r - Sig^r \dagger] , hence the -1
                        enddo
                    endif
                enddo
            enddo
        end do  
        !$omp end do
        deallocate(sig,B,C,sig_lesser,sig_greater,gamma_lead)
        !$omp end parallel
    end subroutine calc_gf


    subroutine calc_w(NBC,NB,NS,nm_dev,PR,PL,PG,V,WR,WL,WG)
        integer,intent(in)::nm_dev,NB,NS,NBC
        complex(dp),intent(in),dimension(nm_dev,nm_dev)::PR,PL,PG,V
        complex(dp),intent(out),dimension(nm_dev,nm_dev)::WR,WL,WG
        ! --------- local
        complex(dp),allocatable,dimension(:,:)::B,S,M,LL,LG,VV
        complex(dp),dimension(:,:),allocatable::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,&
            PL00,PL01,PL10,PG00,PG01,PG10,LL00,LL01,LL10,LG00,LG01,LG10
        complex(dp),dimension(:,:),allocatable::VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,&
            MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N
        complex(dp),dimension(:,:),allocatable::dM11,xR11,dLL11,dLG11,dV11
        complex(dp),dimension(:,:),allocatable::dMnn,xRnn,dLLnn,dLGnn,dVnn
        integer::i,NL,NR,NT,LBsize,RBsize
        real(dp)::condL,condR
        NL=NB*NS ! left contact block size
        NR=NB*NS ! right contact block size
        NT=nm_dev! total size
        LBsize=NL*NBC
        RBsize=NR*NBC
        if (NBC>0) then
          allocate(B(NT,NT))
          allocate(M(NT,NT))
          allocate(S(NT,NT))
          allocate(LL(NT,NT))
          allocate(LG(NT,NT))
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
          !
          call get_OBC_blocks_for_W(NL,V(1:NL,1:NL),V(1:NL,NL+1:2*NL),PR(1:NL,1:NL),PR(1:NL,NL+1:2*NL),&
              PL(1:NL,1:NL),PL(1:NL,NL+1:2*NL),PG(1:NL,1:NL),PG(1:NL,NL+1:2*NL),NBC,&
              V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
              LL00,LL01,LL10,LG00,LG01,LG10)
          !    
          call get_OBC_blocks_for_W(NR,V(NT-NR+1:NT,NT-NR+1:NT),transpose(conjg(V(NT-NR+1:NT,NT-2*NR+1:NT-NR))),PR(NT-NR+1:NT,NT-NR+1:NT),&
              transpose(PR(NT-NR+1:NT,NT-2*NR+1:NT-NR)),PL(NT-NR+1:NT,NT-NR+1:NT),-transpose(conjg(PL(NT-NR+1:NT,NT-2*NR+1:NT-NR))),&
              PG(NT-NR+1:NT,NT-NR+1:NT),-transpose(conjg(PG(NT-NR+1:NT,NT-2*NR+1:NT-NR))),NBC,&
              VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,&
              LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N)
          !
          !! S = V P^r
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,V,nm_dev,PR,nm_dev,czero,S,nm_dev)
          ! Correct first and last block to account for elements in the contacts
          S(1:LBsize,1:LBsize)=S(1:LBsize,1:LBsize) + matmul(V10,PR01)
          S(NT-RBsize+1:NT,NT-RBsize+1:NT)=S(NT-RBsize+1:NT,NT-RBsize+1:NT) + matmul(VNN1,PRN1N)
          !
          M = -S
          do i=1,nm_dev
             M(i,i) = 1.0d0 + M(i,i)
          enddo
          deallocate(S)
          !! LL=V P^l V'
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,V,nm_dev,PL,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V,nm_dev,czero,LL,nm_dev) 
          !Correct first and last block to account for elements in the contacts
          LL(1:LBsize,1:LBsize)=LL(1:LBsize,1:LBsize) + matmul(matmul(V10,PL00),V01) + &
            matmul(matmul(V10,PL01),V00) + matmul(matmul(V00,PL10),V01)
          !  
          LL(NT-RBsize+1:NT,NT-RBsize+1:NT)=LL(NT-RBsize+1:NT,NT-RBsize+1:NT) + &
            matmul(matmul(VNN,PLNN1),VN1N) + matmul(matmul(VNN1,PLN1N),VNN) + &
            matmul(matmul(VNN1,PLNN),VN1N)
          !
          !! LG=V P^g V'    
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,V,nm_dev,PG,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V,nm_dev,czero,LG,nm_dev) 
          !Correct first and last block to account for elements in the contacts
          LG(1:LBsize,1:LBsize)=LG(1:LBsize,1:LBsize) + matmul(matmul(V10,PG00),V01) + &
            matmul(matmul(V10,PG01),V00) + matmul(matmul(V00,PG10),V01)
          LG(NT-RBsize+1:NT,NT-RBsize+1:NT)=LG(NT-RBsize+1:NT,NT-RBsize+1:NT) + &
            matmul(matmul(VNN,PGNN1),VN1N) + matmul(matmul(VNN1,PGN1N),VNN) + &
            matmul(matmul(VNN1,PGNN),VN1N)
            
          ! WR/WL/WG OBC Left
          call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
          ! WR/WL/WG OBC right
          call open_boundary_conditions(NR,MNN,MNN1,MN1N,VN1N,xRNN,dMNN,dVNN,condR)
          allocate(VV(nm_dev,nm_dev))
          VV = V
          if (condL<1.0d-6) then   
              !
            !   call get_dL_OBC_for_W(NL,xR11,LL00,LL01,LG00,LG01,M10,'L', dLL11,dLG11)
              !
              M(1:LBsize,1:LBsize)=M(1:LBsize,1:LBsize) - dM11
              VV(1:LBsize,1:LBsize)=V(1:LBsize,1:LBsize) - dV11    
            !   LL(1:LBsize,1:LBsize)=LL(1:LBsize,1:LBsize) + dLL11
            !   LG(1:LBsize,1:LBsize)=LG(1:LBsize,1:LBsize) + dLG11    
          endif
          if (condR<1.0d-6) then    
              !
            !   call get_dL_OBC_for_W(NR,xRNN,LLNN,LLN1N,LGNN,LGN1N,MNN1,'R', dLLNN,dLGNN)
              !
              M(NT-RBsize+1:NT,NT-RBsize+1:NT)=M(NT-RBsize+1:NT,NT-RBsize+1:NT) - dMNN
              VV(NT-RBsize+1:NT,NT-RBsize+1:NT)=V(NT-RBsize+1:NT,NT-RBsize+1:NT)- dVNN
            !   LL(NT-RBsize+1:NT,NT-RBsize+1:NT)=LL(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLLNN
            !   LG(NT-RBsize+1:NT,NT-RBsize+1:NT)=LG(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLGNN    
          endif
          !!!! calculate W^r = (I - V P^r)^-1 V    
          call invert_inplace(M,nm_dev) ! M -> xR
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,VV,nm_dev,czero,WR,nm_dev)           
          ! calculate W^< and W^> = W^r P^<> W^r dagger
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,LL,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,WL,nm_dev) 
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,LG,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,WG,nm_dev)  
          deallocate(M,LL,LG,B,VV)
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
        else ! no OBC correction
          allocate(B(NT,NT))
          allocate(M(NT,NT))    
          !! M = I - V P^r
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,PR,nm_dev,czero,M,nm_dev)
          !
          do i=1,nm_dev
             M(i,i) = 1.0d0 + M(i,i)
          enddo    
          !!!! calculate W^r = (I - V P^r)^-1 V    
          call invert_inplace(M,nm_dev) ! M -> xR
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,V,nm_dev,czero,WR,nm_dev)           
          ! calculate W^< and W^> = W^r P^<> W^r dagger
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PL,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WL,nm_dev) 
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PG,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WG,nm_dev)  
          deallocate(B,M)  
        endif 
    end subroutine calc_w

    
    Function trace(A,nn) 
        integer :: nn,i        
        complex(8), dimension(nn,nn),intent(in) :: A
        complex(8) :: trace, tr
        tr=dcmplx(0.0d0,0.0d0)
        do i=1,nn
            tr=tr+A(i,i)
        enddo
        trace=tr
    end function trace

end module gw_dense


