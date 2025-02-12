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

use fft_mod, only: conv1d => conv1d2, corr1d => corr1d2 
use parameters_mod, only: dp, c1i, cone, czero, pi, twopi, two, hbar, e_charge
use green, only: green_calc_g

implicit none 

private

public:: green_solve_gw_1D

CONTAINS

! driver subroutine for iterating the GW self-energy SCBA loop: G -> P -> W -> Sig 
subroutine green_solve_gw_1D(scba_tol,niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,conv_method,mid_bandgap,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,lead_coupling,V,&
  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
  W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,charge)
  !
  integer, intent(in) :: nen,nb,ns,niter,nm_dev,length
  real(dp), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg,scba_tol, mid_bandgap(nm_dev)
  complex(dp),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),lead_coupling(NB*NS,nm_dev,2)
  complex(dp), intent(in):: V(nm_dev,nm_dev)
  logical,intent(in)::ldiag
  character(len=*),intent(in)::conv_method
  complex(dp),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater
  complex(dp),intent(inout),dimension(nm_dev,nm_dev,nen) ::  Sig_retarded,Sig_lesser,Sig_greater
  complex(dp),intent(inout),dimension(nm_dev,nm_dev,nen) ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
  complex(dp),intent(inout),dimension(nm_dev,nm_dev,nen) ::  P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater
  real(dp),dimension(nm_dev),intent(inout)::charge
  !---- local variables
  complex(dp),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
  complex(dp),allocatable,dimension(:,:):: B ! tmp matrix
  real(dp),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:)
  real(dp),allocatable::Tr(:,:) ! current spectrum on leads
  real(dp),allocatable::Te(:,:,:) ! transmission matrix spectrum
  integer :: iter,ie,nopmax,iep
  integer :: i,j,nm,nop,l,h,iop,ndiag
  complex(dp),allocatable::Ispec(:,:,:),Itot(:,:)  
  real(dp)::scba_error
  real(dp)::E_hartree
  complex(dp) :: dE
  real(dp)::nelec(2),mu(2),pelec(2), tmp, G_lesser_sum
  real(dp),dimension(nm_dev)::ndens,ndens_intrinsic,pdens,pdens_intrinsic
  complex(dp)::tmp2(nen)
  integer::ie_gap(nm_dev),ie_gap_old(nm_dev) ! energy index of mid-bandgap energy
  !
  print *,'============ green_solve_gw_1D ============'  
  allocate(siglead(NB*NS,NB*NS,nen,2))
  allocate(wen(nen))
  wen(:)=en(:)-en(nen/2)
  ! get leads sigma
  siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
  siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)  
  allocate(B(nm_dev,nm_dev))
  allocate(tot_cur(nm_dev,nm_dev))
  allocate(tot_ecur(nm_dev,nm_dev))
  allocate(cur(nm_dev,nm_dev,nen))
  allocate(Ispec(nm_dev,nm_dev,nen))
  allocate(Itot(nm_dev,nm_dev))
  allocate(tr(nen,2))
  mu=(/ mus, mud /)
  print '(a5,f15.4,a5,f15.6)', 'mus=',mu(1),'mud=',mu(2)
  !$omp parallel default(shared) private(i)  
  !$omp do
  do i=1,nm_dev          
    ie_gap(i) = (mid_bandgap(i)-en(1)) / (en(2) - en(1)) + 1
  enddo        
  !$omp end do
  !$omp end parallel
  !
  do iter=0,niter
    print '("+ iter=",i8)',iter    
    print *, '  calc G'  
    call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,&
      Ham,H00lead,H10lead,Siglead,lead_coupling,&
      Sig_retarded,Sig_lesser,Sig_greater,&
      G_retarded,G_lesser,G_greater,mu=mu,temp=(/temps,tempd/),cur=Tr)
    !
    call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
    !
    call write_current_spectrum('gw_Jdens',iter,cur,nen,en,length,NB,Lx)
    call write_current('gw_I',iter,tot_cur,length,NB,NS,Lx)
    call write_current('gw_EI',iter,tot_ecur,length,NB,NS,Lx)
    call write_spectrum('gw_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
    call write_spectrum('gw_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('gw_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))
    call write_transmission_spectrum('gw_trL',iter,Tr(:,1)*spindeg,nen,En)
    call write_transmission_spectrum('gw_trR',iter,Tr(:,2)*spindeg,nen,En)
    !call write_matrix_summed_overE('Gr',iter,G_retarded,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('G_r',iter,G_retarded,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('G_l',iter,G_lesser,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('G_g',iter,G_greater,nen,en,length,NB,(/1.0,1.0/))
    !
    open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
    write(101,'(I4,2E16.6)') iter, - sum(Tr(:,1))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg), &
      sum(Tr(:,2))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)
    close(101)
    write(6,'(I4," ID=",2E16.6)') iter, - sum(Tr(:,1))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg), &
      sum(Tr(:,2))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)
    !        
    call write_potprofile('gw_midgap',iter,(ie_gap(:)-1)*(en(2)-en(1))+en(1),length,NB,Lx)
    !$omp parallel default(shared) private(i)  
    !$omp do
    do i=1,nm_dev      
      ndens(i) = dimag(sum(G_lesser(i,i,ie_gap(i):))) * (en(2)-en(1)) / twopi * spindeg 
      pdens(i) =-dimag(sum(G_greater(i,i,:ie_gap(i)))) * (en(2)-en(1)) / twopi * spindeg 
    enddo
    !$omp end do
    !$omp end parallel
    call write_charge('gw_n',iter,ndens,length,NB,Lx)
    call write_charge('gw_p',iter,pdens,length,NB,Lx)
    if (iter < 1) then
      ndens_intrinsic(:) = ndens(:)
      pdens_intrinsic(:) = pdens(:)
    endif
    !        
    print *, '  calc P'  
    ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
    ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
    P_lesser=czero
    P_greater=czero
    P_retarded=czero
    !
    ndiag=NB*NS
    if (ldiag) ndiag=0  
    !
    nopmax=nen/2-1
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / twopi )	* spindeg    
    !$omp parallel default(shared) private(l,h,i,j)  
    !$omp do
    do i=1,nm_dev
      l=max(i-ndiag,1)
      h=min(nm_dev,i+ndiag)
      do j=l,h
        P_lesser(i,j,:) =corr1d(nen,G_lesser(i,j,:),G_greater(j,i,:),method=conv_method)
        P_greater(i,j,:)=corr1d(nen,G_greater(i,j,:),G_lesser(j,i,:),method=conv_method)        
      enddo
    enddo
    !$omp end do
    !$omp end parallel
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE
    P_retarded=(P_greater - P_lesser) / two
    !
!    call write_spectrum('PR',0,P_retarded,nen,wen,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('PL',0,P_lesser,nen,wen,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('PG',0,P_greater,nen,wen,length,NB,Lx,(/1.0d0,1.0d0/))    
    !
    print *, '  calc W'  
    W_lesser=czero
    W_greater=czero
    W_retarded=czero
    !
    !$omp parallel default(shared) private(nop)
    !$omp do
    do nop=1,nen
      call green_calc_w(1,NB,NS,nm_dev,P_retarded(:,:,nop),P_lesser(:,:,nop),P_greater(:,:,nop),V,W_retarded(:,:,nop),W_lesser(:,:,nop),W_greater(:,:,nop))
    enddo
    !$omp end do
    !$omp end parallel
!    call write_spectrum('WR',0,W_retarded,nen,wen,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('WL',0,W_lesser,  nen,wen,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('WG',0,W_greater, nen,wen,length,NB,Lx,(/1.0d0,1.0d0/))    
    !
    print *, '  calc SigGW'
    Sig_greater_new = dcmplx(0.0d0,0.0d0)
    Sig_lesser_new = dcmplx(0.0d0,0.0d0)
    Sig_retarded_new = dcmplx(0.0d0,0.0d0)  
    !
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)
    !$omp parallel default(shared) private(l,h,i,j)
    !$omp do
    do i=1,nm_dev
      l=max(i-ndiag,1)
      h=min(nm_dev,i+ndiag)
      do j=l,h
        Sig_lesser_new(i,j,:)  =conv1d(nen,G_lesser(i,j,:),W_lesser(i,j,:),method=conv_method)
        Sig_greater_new(i,j,:) =conv1d(nen,G_greater(i,j,:),W_greater(i,j,:),method=conv_method)
        Sig_retarded_new(i,j,:)=conv1d(nen,G_lesser(i,j,:),W_retarded(i,j,:),method=conv_method) +&
                                conv1d(nen,G_retarded(i,j,:),W_lesser(i,j,:),method=conv_method) +&
                                conv1d(nen,G_retarded(i,j,:),W_retarded(i,j,:),method=conv_method)                                
                                
      enddo
    enddo                            
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0, (En(2)-En(1))/twopi)        
    Sig_lesser_new = Sig_lesser_new  * dE
    Sig_greater_new= Sig_greater_new * dE    
    Sig_retarded_new=Sig_retarded_new* dE
    !
    Sig_retarded_new = dcmplx( dble(Sig_retarded_new), dimag(Sig_greater_new-Sig_lesser_new)/two )
    Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), dimag(Sig_lesser_new) )
    Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), dimag(Sig_greater_new) )
    !    
    ! Hartree potential    
    ie_gap_old(:) = ie_gap(:)
    if (iter>=1) then
      !$omp parallel default(shared) private(i,E_hartree)  
      !$omp do
      do i=1,nm_dev      
        E_hartree = dble(sum(V(i,:) * ( ndens(:) - ndens_intrinsic(:) + (pdens(:) - pdens_intrinsic(:)))))        
        ! left boundary correction
        if (i <= NB*NS) then 
          E_hartree = E_hartree + dble(sum(V(i+NB*NS,:NB*NS) * ( ndens(:NB*NS) - ndens_intrinsic(:NB*NS) + (pdens(:NB*NS) - pdens_intrinsic(:NB*NS)))))
        endif 
        ! right boundary correction
        if (i > (nm_dev-NB)) then 
          E_hartree = E_hartree + dble(sum(V(i-NB,nm_dev-NB+1:) * ( ndens(nm_dev-NB+1:) - ndens_intrinsic(nm_dev-NB+1:) + (pdens(nm_dev-NB+1:) - pdens_intrinsic(nm_dev-NB+1:)))))
        endif 
        !
        Sig_retarded_new(i,i,:) = Sig_retarded_new(i,i,:) + E_hartree
        ie_gap(i) = ( (mid_bandgap(i)+E_hartree-en(1)) / (en(2) - en(1)) + 1 - ie_gap_old(i) ) * alpha_mix + ie_gap_old(i)        
        ie_gap(i) = max(1, ie_gap(i))
        ie_gap(i) = min(nen, ie_gap(i))
      enddo        
      !$omp end do
      !$omp end parallel      
    endif
    !
    ! symmetrize the self-energies
    do ie=1,nen
      B(:,:)=transpose(Sig_retarded_new(:,:,ie))
      Sig_retarded_new(:,:,ie) = (Sig_retarded_new(:,:,ie) + B(:,:))/2.0d0
      B(:,:)=transpose(Sig_lesser_new(:,:,ie))
      Sig_lesser_new(:,:,ie) = dcmplx((dble(Sig_lesser_new(:,:,ie)) - dble(B(:,:)))/2.0d0, (dimag(Sig_lesser_new(:,:,ie)) + dimag(B(:,:)))/2.0d0)
      B(:,:)=transpose(Sig_greater_new(:,:,ie))
      Sig_greater_new(:,:,ie) = dcmplx((dble(Sig_greater_new(:,:,ie)) - dble(B(:,:)))/2.0d0, (dimag(Sig_greater_new(:,:,ie)) + dimag(B(:,:)))/2.0d0)
    enddo
    !
    scba_error = sqrt(sum( abs(Sig_retarded_new - Sig_retarded)**2 )) / sqrt(sum( abs(Sig_retarded_new)**2 ))
    open(unit=101,file='gw_scba_error.dat',status='unknown',position='append')
    write(101,'(I4,E16.6)') iter, scba_error
    close(101)
    write(6,'(I4," error=",F16.10)') iter, scba_error
    ! mixing with the previous one and update
    Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
    Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
    Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)    
    !
    ! get leads sigma
    siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
    siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)      
    !
!    call write_spectrum('SigR',0,Sig_retarded,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('SigL',0,Sig_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('SigG',0,Sig_greater,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
    !
    ! calculate collision integral
    call calc_collision(Sig_lesser_new,Sig_greater_new,G_lesser,G_greater,nen,en,spindeg,nm_dev,Itot,Ispec)
    call write_spectrum('gw_Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0d0,0.0d0/))
    !
    if (scba_error < scba_tol) then
      exit
    endif
    if (iter == niter) then 
      print *, "warning: max number of iterations reached!"
    endif
  enddo ! iter loop                 
  !  
  ! calculate GF for the last time      
  print *, 'calc G for the last time'  
  call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,&
      Ham,H00lead,H10lead,Siglead,lead_coupling,&
      Sig_retarded,Sig_lesser,Sig_greater,&
      G_retarded,G_lesser,G_greater,mu=mu,temp=(/temps,tempd/),cur=Tr)
  !
  call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
  !
  call write_current_spectrum('gw_Jdens',iter,cur,nen,en,length,NB,Lx)
  call write_current('gw_I',iter,tot_cur,length,NB,NS,Lx)
  call write_current('gw_EI',iter,tot_ecur,length,NB,NS,Lx)
  call write_spectrum('gw_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
  call write_spectrum('gw_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
!    call write_spectrum('gw_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))
  call write_transmission_spectrum('gw_trL',iter,Tr(:,1)*spindeg,nen,En)
  call write_transmission_spectrum('gw_trR',iter,Tr(:,2)*spindeg,nen,En)
  !
  open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
  write(101,'(I4,2E16.6)') iter, - sum(Tr(:,1))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg), &
    sum(Tr(:,2))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)
  close(101)
  write(6,'(I4,2E16.6)') iter, - sum(Tr(:,1))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg), &
    sum(Tr(:,2))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)
  !
  charge(:) = ndens(:) - pdens(:)
end subroutine green_solve_gw_1D


subroutine principle_integral()
!    do ie=1,nen
!      if (ie .ne. (nopmax+1)) then
!        tmp2(ie) = 1.0d0 / (ie-nopmax-1) / twopi      
!      endif
!    enddo
!    tmp2(nopmax+1)=0.0d0
!    !$omp parallel default(shared) private(l,h,i,j,G_lesser_sum)
!    !$omp do
!    do i=1,nm_dev
!      l=max(i-ndiag,1)
!      h=min(nm_dev,i+ndiag)
!      do j=l,h
!        G_lesser_sum = sum( aimag(G_lesser(i,j,:)) ) * dE
        
!        Sig_retarded_new(i,j,:) = dcmplx( dble(conv1d(nen, c1i * (Sig_greater_new(i,j,:) - Sig_lesser_new(i,j,:)), tmp2, method='sum')) &
!                                        + dble(V(i,j) * G_lesser_sum), aimag(Sig_greater_new(i,j,:) - Sig_lesser_new(i,j,:))/two ) 
        
!      enddo
!    enddo
!    !$omp end do
!    !$omp end parallel
end subroutine principle_integral


! calculate matrix blocks for the Open Boundary Condition of W
subroutine get_OBC_blocks_for_W(n,v_00,v_01,pR_00,pR_01,pR_10,pL_00,pL_01,pG_00,pG_01,NBC,&
    V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
    LL00,LL01,LL10,LG00,LG01,LG10)
  integer,intent(in)::n,NBC
  complex(8),intent(in),dimension(n,n)::v_00,v_01,pR_00,pR_01,pR_10,pL_00,pL_01,pG_00,pG_01
  complex(8),intent(out),dimension(n*NBC,n*NBC)::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
      LL00,LL01,LL10,LG00,LG01,LG10
  complex(8),dimension(n*NBC,n*NBC)::II
  !
  select case (NBC)
    !
    case(1)
      !
      V00=v_00
      V01=v_01
      V10=transpose(conjg(V01))
      !
      PR00=pR_00
      PR01=pR_01
      PR10=pR_10
      !
      PL00=pL_00;
      PL01=pL_01;
      PL10= - transpose(conjg(PL01))
      !
      PG00=pG_00
      PG01=pG_01
      PG10= - transpose(conjg(PG01))
      !    
    case(2)
      !
      V00(1:n,1:n)=v_00 
      V00(1:n,n+1:2*n)=v_01
      V00(n+1:2*n,1:n)=transpose(conjg(v_01))
      V00(n+1:2*n,n+1:2*n)= v_00
      V01=czero
      V01(n+1:2*n,1:n)=v_01
      V10=transpose(conjg(V01))
      !
      PR00(1:n,1:n)=pR_00 
      PR00(1:n,n+1:2*n)=pR_01
      PR00(n+1:2*n,1:n)=pR_10
      PR00(n+1:2*n,n+1:2*n)= pR_00
      PR01=czero
      PR01(n+1:2*n,1:n)=pR_01
      PR10=-transpose(conjg(PR01))
      !
      PG00(1:n,1:n)=pG_00 
      PG00(1:n,n+1:2*n)=pG_01
      PG00(n+1:2*n,1:n)=-transpose(conjg(pG_01))
      PG00(n+1:2*n,n+1:2*n)= pG_00
      PG01=czero
      PG01(n+1:2*n,1:n)=pG_01
      PG10=-transpose(conjg(PG01))
      !
      PL00(1:n,1:n)=pL_00 
      PL00(1:n,n+1:2*n)=pL_01
      PL00(n+1:2*n,1:n)=-transpose(conjg(pL_01))
      PL00(n+1:2*n,n+1:2*n)= pL_00
      PL01=czero
      PL01(n+1:2*n,1:n)=pL_01
      PL10=-transpose(conjg(PL01))    
    !
  end select
  !
  call identity(II,NBC*N)
  M00=II*dcmplx(1.0d0,1d-10)-matmul(V10,PR01)
  M00=M00-matmul(V00,PR00)-matmul(V01,PR10)
  M01=-matmul(V00,PR01)-matmul(V01,PR00)
  M10=-matmul(V10,PR00)-matmul(V00,PR10)
  !
  LL00=matmul(matmul(V10,PL00),V01)+matmul(matmul(V10,PL01),V00)
  LL00=LL00+matmul(matmul(V00,PL10),V01)+matmul(matmul(V00,PL00),V00)
  LL00=LL00+matmul(matmul(V00,PL01),V10)
  LL00=LL00+matmul(matmul(V01,PL10),V00)+matmul(matmul(V01,PL00),V10)
  !
  LL01=matmul(matmul(V10,PL01),V01)+matmul(matmul(V00,PL00),V01)
  LL01=LL01+matmul(matmul(V00,PL01),V00)+matmul(matmul(V01,PL10),V01)
  LL01=LL01+matmul(matmul(V01,PL00),V00)
  !
  LL10=-transpose(conjg(LL01))
  !
  LG00=matmul(matmul(V10,PG00),V01)+matmul(matmul(V10,PG01),V00)
  LG00=LG00+matmul(matmul(V00,PG10),V01)+matmul(matmul(V00,PG00),V00)
  LG00=LG00+matmul(matmul(V00,PG01),V10)
  LG00=LG00+matmul(matmul(V01,PG10),V00)+matmul(matmul(V01,PG00),V10)
  !
  LG01=matmul(matmul(V10,PG01),V01)
  LG01=LG01+matmul(matmul(V00,PG00),V01)+matmul(matmul(V00,PG01),V00)
  LG01=LG01+matmul(matmul(V01,PG10),V01)+matmul(matmul(V01,PG00),V00)
  !
  LG10=-transpose(conjg(LG01))  
end subroutine get_OBC_blocks_for_W


! calculate corrections to the L matrix blocks for the Open Boundary Condition
subroutine get_dL_OBC_for_W(nm,xR,LL00,LL01,LG00,LG01,M10,typ, dLL11,dLG11)
  integer,intent(in)::nm
  character(len=*),intent(in)::typ
  complex(8),intent(in),dimension(nm,nm)::xR,LL00,LL01,LG00,LG01,M10
  complex(8),intent(out),dimension(nm,nm)::dLL11,dLG11
  ! -----
  complex(8),dimension(nm,nm)::AL,AG,FL,FG,A,V,iV,yL_NN,wL_NN,yG_NN,wG_NN,tmp1,tmp2
  complex(8),dimension(nm)::E
  integer::i,j
  !!!! AL=M10*xR*LL01;
  !!!! AG=M10*xR*LG01;
  call zgemm('n','n',nm,nm,nm,cone,M10,nm,xR,nm,czero,tmp1,nm)
  call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,LL01,nm,czero,AL,nm)
  call zgemm('n','n',nm,nm,nm,cone,M10,nm,xR,nm,czero,tmp1,nm)
  call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,LG01,nm,czero,AG,nm)
  !!!! FL=xR*(LL00-(AL-AL'))*xR';
  !!!! FG=xR*(LG00-(AG-AG'))*xR';
  call zgemm('n','n',nm,nm,nm,cone,xR,nm,(LL00-(AL-transpose(conjg(AL)))),nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,FL,nm)
  call zgemm('n','n',nm,nm,nm,cone,xR,nm,(LG00-(AG-transpose(conjg(AG)))),nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,FG,nm)
  !
  call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,V,nm)  
  do i=1,nm
    V(i,i)=V(i,i)+dcmplx(0.0d0,1.0d-4)  ! 1i*1e-4 added to stabilize matrix
  enddo
  E=eigv(nm,V)
  iV=V
  call invert(iV,nm)
  !lesser component
  call zgemm('n','n',nm,nm,nm,cone,iV,nm,FL,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,iV,nm,czero,yL_NN,nm)
  yL_NN=yL_NN/(1.0d0 - sum(E*conjg(E)))
  call zgemm('n','n',nm,nm,nm,cone,V,nm,yL_NN,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,V,nm,czero,wL_NN,nm)
  !refinement iteration
  call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,tmp1,nm)
  call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,wL_NN,nm,czero,tmp2,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp2,nm,M10,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,tmp2,nm)
  wL_NN=FL+tmp2
  !
  call zgemm('n','n',nm,nm,nm,cone,M10,nm,wL_NN,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dLL11,nm)
  dLL11=dLL11-(AL-transpose(conjg(AL)))
  !greater component
  call zgemm('n','n',nm,nm,nm,cone,iV,nm,FG,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,iV,nm,czero,yG_NN,nm)
  yG_NN=yG_NN/(1.0d0 - sum(E*conjg(E)))
  call zgemm('n','n',nm,nm,nm,cone,V,nm,yG_NN,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,V,nm,czero,wG_NN,nm)
  !refinement iteration
  call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,tmp1,nm)
  call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,wG_NN,nm,czero,tmp2,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp2,nm,M10,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,tmp2,nm)
  wG_NN=FG+tmp2
  !
  call zgemm('n','n',nm,nm,nm,cone,M10,nm,wG_NN,nm,czero,tmp1,nm)
  call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dLG11,nm)
  dLG11=dLG11-(AG-transpose(conjg(AG)))
end subroutine get_dL_OBC_for_W


subroutine green_calc_w(NBC,NB,NS,nm_dev,PR,PL,PG,V,WR,WL,WG)
  integer,intent(in)::nm_dev,NB,NS,NBC
  complex(8),intent(in),dimension(nm_dev,nm_dev)::PR,PL,PG,V
  complex(8),intent(out),dimension(nm_dev,nm_dev)::WR,WL,WG
  ! --------- local
  complex(8),allocatable,dimension(:,:)::B,S,M,LL,LG,VV
  complex(8),dimension(:,:),allocatable::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,&
      PL00,PL01,PL10,PG00,PG01,PG10,LL00,LL01,LL10,LG00,LG01,LG10
  complex(8),dimension(:,:),allocatable::VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,&
      MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N
  complex(8),dimension(:,:),allocatable::dM11,xR11,dLL11,dLG11,dV11
  complex(8),dimension(:,:),allocatable::dMnn,xRnn,dLLnn,dLGnn,dVnn
  integer::i,NL,NR,NT,LBsize,RBsize
  real(8)::condL,condR
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
    call get_OBC_blocks_for_W(NL,V(1:NL,1:NL),V(1:NL,NL+1:2*NL),PR(1:NL,1:NL),PR(1:NL,NL+1:2*NL),PR(NL+1:2*NL,1:NL),&
        PL(1:NL,1:NL),PL(1:NL,NL+1:2*NL),PG(1:NL,1:NL),PG(1:NL,NL+1:2*NL),NBC,&
        V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
        LL00,LL01,LL10,LG00,LG01,LG10)
    !    
    call get_OBC_blocks_for_W(NR,V(NT-NR+1:NT,NT-NR+1:NT),V(NT-2*NR+1:NT-NR,NT-NR+1:NT),PR(NT-NR+1:NT,NT-NR+1:NT),PR(NT-NR+1:NT,NT-NR+1:NT),&
        PR(NT-2*NR+1:NT-NR,NT-NR+1:NT),PL(NT-NR+1:NT,NT-NR+1:NT),PL(NT-2*NR+1:NT-NR,NT-NR+1:NT),&
        PG(NT-NR+1:NT,NT-NR+1:NT),PG(NT-2*NR+1:NT-NR,NT-NR+1:NT),NBC,&
        VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,&
        LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N)
    !
    !! S = V P^r
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,V,nm_dev,PR,nm_dev,czero,S,nm_dev)
    !
    !! Correct first and last block to account for elements in the contacts
    S(1:LBsize,1:LBsize)=S(1:LBsize,1:LBsize) + matmul(V10,PR01)
    S(NT-RBsize+1:NT,NT-RBsize+1:NT)=S(NT-RBsize+1:NT,NT-RBsize+1:NT) + matmul(VNN1,PRN1N)
    !
    M = -S
    do i=1,nm_dev
       M(i,i) = 1.0d0 + M(i,i)
    enddo
    deallocate(S)
    !
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
    !  
    ! WR/WL/WG OBC Left
    call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
    !
    ! WR/WL/WG OBC right
    call open_boundary_conditions(NR,MNN,MNN1,MN1N,VN1N,xRNN,dMNN,dVNN,condR)
    allocate(VV(nm_dev,nm_dev))
    VV = V
    if (condL<1.0d-6) then   
        !
        !call get_dL_OBC_for_W(NL,xR11,LL00,LL01,LG00,LG01,M10,'L', dLL11,dLG11)
        !
        M(1:LBsize,1:LBsize)=M(1:LBsize,1:LBsize) - dM11
        VV(1:LBsize,1:LBsize)=V(1:LBsize,1:LBsize) - dV11    
!        LL(1:LBsize,1:LBsize)=LL(1:LBsize,1:LBsize) + dLL11
!        LG(1:LBsize,1:LBsize)=LG(1:LBsize,1:LBsize) + dLG11    
    endif
    if (condR<1.0d-6) then    
        !
        !call get_dL_OBC_for_W(NR,xRNN,LLNN,LLN1N,LGNN,LGN1N,MNN1,'R', dLLNN,dLGNN)
        !
        M(NT-RBsize+1:NT,NT-RBsize+1:NT)=M(NT-RBsize+1:NT,NT-RBsize+1:NT) - dMNN
        VV(NT-RBsize+1:NT,NT-RBsize+1:NT)=V(NT-RBsize+1:NT,NT-RBsize+1:NT)- dVNN
!        LL(NT-RBsize+1:NT,NT-RBsize+1:NT)=LL(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLLNN
!        LG(NT-RBsize+1:NT,NT-RBsize+1:NT)=LG(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLGNN    
    endif
    !
    !! calculate W^r = (I - V P^r)^-1 V    
    call invert(M,nm_dev) ! M -> xR
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,VV,nm_dev,czero,WR,nm_dev)           
    !
    !! calculate W^< and W^> = W^r P^<> W^r dagger
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
    call invert(M,nm_dev) ! M -> xR
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,V,nm_dev,czero,WR,nm_dev)           
    ! calculate W^< and W^> = W^r P^<> W^r dagger
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PL,nm_dev,czero,B,nm_dev) 
    call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WL,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PG,nm_dev,czero,B,nm_dev) 
    call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WG,nm_dev)  
    !
    deallocate(B,M)  
  endif 
end subroutine green_calc_w

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

FUNCTION eigv(NN, A)
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

  CALL zheev( 'V','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )

  deallocate(work,rwork)
  if (INFO.ne.0)then
     write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
     call abort
  endif
  eigv(:)=W(:)
END FUNCTION eigv


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

subroutine identity(A,n)
  integer, intent(in) :: n        
  complex(8), dimension(n,n), intent(inout) :: A
  integer :: i
  A = dcmplx(0.0d0,0.0d0)
  do i = 1,n
    A(i,i) = dcmplx(1.0d0,0.0d0)
  end do
end subroutine identity

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


subroutine expand_size_bycopy(A,nm,nb,add)
  complex(8),intent(inout)::A(nm,nm)
  integer, intent(in)::nm,add,nb
  integer::i,nm0,l,l2
  nm0=nm-nb*add*2
  A(1:add*nb,:)=0.0d0
  A(:,1:add*nb)=0.0d0
  A(add*nb+nm0+1:nm,:)=0.0d0
  A(:,add*nb+nm0+1:nm)=0.0d0
  do i=0,add-1
    A(i*nb+1:i*nb+nb,i*nb+1:i*nb+nm0)=A(add*nb+1:add*nb+nb,add*nb+1:add*nb+nm0)
    A(i*nb+1:i*nb+nm0,i*nb+1:i*nb+nb)=A(add*nb+1:add*nb+nm0,add*nb+1:add*nb+nb)
    l=add*nb+nm0+i*nb
    l2=add*nb+i*nb+nb
    A(l+1:l+nb,l2+1:l2+nm0)=A(add*nb+nm0-nb+1:add*nb+nm0,add*nb+1:add*nb+nm0)
    A(l2+1:l2+nm0,l+1:l+nb)=A(add*nb+1:add*nb+nm0,add*nb+nm0-nb+1:add*nb+nm0)  
  enddo
end subroutine expand_size_bycopy


subroutine homogenize_selfenergy(Sig_retarded,Sig_lesser,Sig_greater,nm_dev,nen,NB,NS)
  integer,intent(in) :: nm_dev, nen, NB, NS
  complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  Sig_retarded,Sig_lesser,Sig_greater
  integer::ie
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    call expand_size_bycopy(Sig_retarded(:,:,ie),nm_dev,NB,NS+1)
    call expand_size_bycopy(Sig_lesser(:,:,ie),nm_dev,NB,NS+1)
    call expand_size_bycopy(Sig_greater(:,:,ie),nm_dev,NB,NS+1)
  enddo
end subroutine homogenize_selfenergy


!-------------------------  OBSERVABLES  -------------------------------

! calculate bond current using I_ij = H_ij G<_ji - H_ji G^<_ij
subroutine calc_bond_current(H,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
  complex(8),intent(in)::H(nm_dev,nm_dev),G_lesser(nm_dev,nm_dev,nen)
  real(8),intent(in)::en(nen),spindeg
  integer,intent(in)::nen,nm_dev ! number of E and device dimension
  real(8),intent(out)::tot_cur(nm_dev,nm_dev) ! total bond current density
  real(8),intent(out),optional::tot_ecur(nm_dev,nm_dev) ! total bond energy current density
  real(8),intent(out),optional::cur(nm_dev,nm_dev,nen) ! energy resolved bond current density
  !----
  complex(8),allocatable::B(:,:)
  integer::ie,io,jo
  real(8),parameter::tpi=6.28318530718  
  allocate(B(nm_dev,nm_dev))
  tot_cur=0.0d0  
  tot_ecur=0.0d0
  do ie=1,nen
    do io=1,nm_dev
      !$omp parallel default(shared) private(jo)  
      !$omp do
      do jo=1,nm_dev
        B(io,jo)=H(io,jo)*G_lesser(jo,io,ie) - H(jo,io)*G_lesser(io,jo,ie)
      enddo
      !$omp end do
      !$omp end parallel
    enddo    
    B=B*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)
    if (present(cur)) cur(:,:,ie) = dble(B)
    if (present(tot_ecur)) tot_ecur=tot_ecur+ en(ie)*dble(B)
    tot_cur=tot_cur+ dble(B)          
  enddo
  deallocate(B)
end subroutine calc_bond_current

! calculate scattering collision integral from the self-energy
! I = Sig> G^< - Sig< G^>
subroutine calc_collision(Sig_lesser,Sig_greater,G_lesser,G_greater,nen,en,spindeg,nm_dev,I,Ispec)
  complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::G_greater,G_lesser,Sig_lesser,Sig_greater
  real(8),intent(in)::en(nen),spindeg
  integer,intent(in)::nen,nm_dev
  complex(8),intent(out)::I(nm_dev,nm_dev) ! collision integral
  complex(8),intent(out),optional::Ispec(nm_dev,nm_dev,nen) ! collision integral spectrum
  !----
  complex(8),allocatable::B(:,:)
  integer::ie
  real(8),parameter::tpi=6.28318530718  
  allocate(B(nm_dev,nm_dev))
  I=dcmplx(0.0d0,0.0d0)
  do ie=1,nen
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,Sig_greater(:,:,ie),nm_dev,G_lesser(:,:,ie),nm_dev,czero,B,nm_dev)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Sig_lesser(:,:,ie),nm_dev,G_greater(:,:,ie),nm_dev,cone,B,nm_dev) 
    I(:,:)=I(:,:)+B(:,:)
    if (present(Ispec)) Ispec(:,:,ie)=B(:,:)*spindeg
  enddo
  I(:,:)=I(:,:)*dble(en(2)-en(1))/tpi*spindeg
  deallocate(B)
end subroutine calc_collision


!----------------------------  OUTPUTS  --------------------------------


! write current into file 
subroutine write_current(dataset,i,cur,length,NB,NS,Lx)
  character(len=*), intent(in) :: dataset
  real(8), intent(in) :: cur(:,:)
  integer, intent(in)::i,length,NB,NS
  real(8), intent(in)::Lx
  integer:: j,ib,jb,ii
  real(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ii = 1,length-1
    tr=0.0d0          
    do ib=1,nb        
      do j=ii,min(ii+NS-1,length-1)
        do jb=1,nb       
          tr = tr+ cur((ii-1)*nb+ib,j*nb+jb)
        enddo
      enddo                        
    end do
    write(11,'(2E18.4)') dble(ii)*Lx, tr
  end do
end subroutine write_current


! write charge into file 
subroutine write_charge(dataset,i,charge,length,NB,Lx)
  character(len=*), intent(in) :: dataset
  real(8), intent(in) :: charge(:)
  integer, intent(in)::i,length,NB
  real(8), intent(in)::Lx
  integer:: j,ib,jb,ii
  real(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ii = 1,length
    tr=0.0d0          
    do ib=1,nb        
      tr = tr+ charge((ii-1)*nb+ib)        
    end do
    write(11,'(2E18.4)') dble(ii-1)*Lx, tr
  end do
  close(11)
end subroutine write_charge

! write potential profile into file 
subroutine write_potprofile(dataset,i,pot,length,NB,Lx)
  character(len=*), intent(in) :: dataset
  real(8), intent(in) :: pot(:)
  integer, intent(in)::i,length,NB
  real(8), intent(in)::Lx
  integer:: j,ib,jb,ii
  real(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ii = 1,length
    tr=0.0d0          
    do ib=1,nb        
      tr = tr+ pot((ii-1)*nb+ib)        
    end do
    write(11,'(2E18.4)') dble(ii-1)*Lx, tr/dble(nb)
  end do
  close(11)
end subroutine write_potprofile


! write current spectrum into file (pm3d map)
subroutine write_current_spectrum(dataset,i,cur,nen,en,length,NB,Lx)
  character(len=*), intent(in) :: dataset
  real(8), intent(in) :: cur(:,:,:)
  integer, intent(in)::i,nen,length,NB
  real(8), intent(in)::Lx,en(nen)
  integer:: ie,j,ib,jb
  real(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ie = 1,nen
      do j = 1,length-1
          tr=0.0d0          
          do ib=1,nb  
            do jb=1,nb        
              tr = tr+ cur((j-1)*nb+ib,j*nb+jb,ie)
            enddo                        
          end do
          write(11,'(3E18.4)') dble(j)*Lx, en(ie), tr
      end do
      write(11,*)    
  end do
  close(11)
end subroutine write_current_spectrum


! write spectrum into file (pm3d map)
subroutine write_spectrum(dataset,i,G,nen,en,length,NB,Lx,coeff)
  character(len=*), intent(in) :: dataset
  complex(8), intent(in) :: G(:,:,:)
  integer, intent(in)::i,nen,length,NB
  real(8), intent(in)::Lx,en(nen),coeff(2)
  integer:: ie,j,ib
  complex(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ie = 1,nen
      do j = 1,length
          tr=0.0d0          
          do ib=1,nb
              tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie)            
          end do
          write(11,'(4E18.4)') (j-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
      end do
      write(11,*)    
  end do
  close(11)
end subroutine write_spectrum


! write transmission spectrum into file
subroutine write_transmission_spectrum(dataset,i,tr,nen,en)
  character(len=*), intent(in) :: dataset
  real(8), intent(in) :: tr(:)
  integer, intent(in)::i,nen
  real(8), intent(in)::en(nen)
  integer:: ie,j,ib
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ie = 1,nen    
    write(11,'(2E18.4)') en(ie), dble(tr(ie))      
  end do
  close(11)
end subroutine write_transmission_spectrum

FUNCTION STRING(inn)
  INTEGER, PARAMETER :: POS= 4
  INTEGER, INTENT(IN) :: inn
  CHARACTER(LEN=POS) :: STRING  
  INTEGER :: cifra, np, mm, num  
  IF (inn > (10**POS)-1) stop "ERROR: (inn > (10**3)-1)  in STRING"
  num= inn
  DO np= 1, POS
      mm= pos-np
      cifra= num/(10**mm)            
      STRING(np:np)= ACHAR(48+cifra)
      num= num - cifra*(10**mm)
  END DO
END FUNCTION STRING


end module gw_dense
