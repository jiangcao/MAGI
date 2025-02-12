!!!!!!!!!!!!!!!! AUTHOR: Jiang Cao
!!!!!!!!!!!!!!!! DATE: 11/2022

module green

implicit none 

private

public :: green_calc_g,green_calc_w
public :: green_solve_gw_1D,green_solve_gw_2D
public :: green_solve_ephoton_freespace_1D
public :: green_solve_gw_ephoton_1D
public :: green_solve_gw_1D_memsaving,green_solve_gw_2D_memsaving
public :: get_OBC_blocks_for_W,get_dL_OBC_for_W
public :: green_solve_gw_1D_supermemsaving

complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), parameter :: hbar=1.0546d-34 ! m^2 kg / s
real(8), parameter :: m0=9.109d-31 ! kg
real(8), parameter :: eps0=8.854d-12 ! C/V/m 
real(8), parameter :: c0=2.998d8 ! m/s
real(8), parameter :: e0=1.6022d-19 ! C
REAL(8), PARAMETER :: pi = 3.14159265359d0

CONTAINS

! driver for solving the GW and e-photon SCBA together   
subroutine green_solve_gw_ephoton_1D(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,midgap,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,V,&
  Pmn,polarization,intensity,hw,labs,&
  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
  W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,encut,Egap)
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg,Egap,midgap(2)
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: V(nm_dev,nm_dev)
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
real(8), intent(in) :: polarization(3) ! light polarization vector 
real(8), intent(in) :: intensity ! [W/m^2]
real(8), intent(in) :: hw ! hw is photon energy in eV
complex(8), intent(in):: Pmn(nm_dev,nm_dev,3) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
logical,intent(in)::ldiag
real(8),intent(in)::encut(2) ! intraband and interband cutoff for P
logical, intent(in) :: labs ! whether to calculate Pi and absorption
!----
complex(8),dimension(:,:),allocatable ::  Pi_retarded,Pi_lesser,Pi_greater,M
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:)
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
integer::iter,ndiagmin,nop,i
  print *,'======================================='
  print *,'====== green_solve_gw_ephoton_1D ======'
  print *,'======================================='
  print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity  
  nop=floor(hw / (En(2)-En(1)))
  print *,'nop=',nop
  allocate(tot_cur(nm_dev,nm_dev))
  allocate(tot_ecur(nm_dev,nm_dev))
  allocate(cur(nm_dev,nm_dev,nen))
  allocate(Ispec(nm_dev,nm_dev,nen))
  allocate(Itot(nm_dev,nm_dev))    
  do iter=0,niter
    ndiagmin=NB*(min(NS,iter))
    if (ldiag) ndiagmin = 0
    call green_solve_gw_1D_memsaving(5,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,midgap,&
            alpha_mix,nen,En,nb,ns,Ham(:,:),H00lead(:,:,:),H10lead(:,:,:),T(:,:,:),V(:,:),&
            G_retarded(:,:,:),G_lesser(:,:,:),G_greater(:,:,:),&
            Sig_retarded(:,:,:),Sig_lesser(:,:,:),Sig_greater(:,:,:),&
            Sig_retarded_new(:,:,:),Sig_lesser_new(:,:,:),Sig_greater_new(:,:,:),ldiag,encut,Egap,ndiagmin=ndiagmin)
    call green_solve_ephoton_freespace_1D(0,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
            0.0d0,nen,En,nb,ns,Ham(:,:),H00lead(:,:,:),H10lead(:,:,:),T(:,:,:),&
            Pmn(:,:,:),polarization,intensity,hw,.false.,&
            G_retarded(:,:,:),G_lesser(:,:,:),G_greater(:,:,:),Sig_retarded(:,:,:),Sig_lesser(:,:,:),Sig_greater(:,:,:),&
            Sig_retarded_new(:,:,:),Sig_lesser_new(:,:,:),Sig_greater_new(:,:,:))       
    ! combine e-photon Sig to GW Sig
    Sig_retarded = Sig_retarded+ Sig_retarded_new 
    Sig_lesser  = Sig_lesser+ Sig_lesser_new 
    Sig_greater = Sig_greater+ Sig_greater_new 
    call write_spectrum('gw_eph_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
    !call write_spectrum('gw_eph_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
    !call write_spectrum('gw_eph_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))
    call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
    call write_current_spectrum('gw_eph_Jdens',iter,cur,nen,en,length,NB,Lx)
    call write_current('gw_eph_I',iter,tot_cur,length,NB,1,Lx)
    call write_current('gw_eph_EI',iter,tot_ecur,length,NB,1,Lx)    
  enddo  
  if (labs) then
    print *, 'calc Pi'
    allocate(Pi_retarded(nm_dev,nm_dev))
    allocate(Pi_lesser(nm_dev,nm_dev))
    allocate(Pi_greater(nm_dev,nm_dev))
    allocate(M(nm_dev,nm_dev))
    M=dcmplx(0.0d0,0.0d0)  
    do i=1,3
      M=M+ polarization(i) * Pmn(:,:,i) 
    enddo  
    open(unit=99,file='gw_eph_totabs'//TRIM(STRING(iter))//'.dat',status='unknown')
    do i=1,floor(4.0d0/(En(2)-En(1)))
        call calc_pi_gw_ephoton(spindeg,nm_dev,length,nen,En,i,M,G_lesser,G_greater,G_retarded,V,Pi_retarded,Pi_lesser,Pi_greater)  
        call write_trace('gw_eph_absorp',iter,Pi_retarded,length,NB,Lx,(/1.0d0,-1.0d0/),E=dble(i)*(En(2)-En(1)))
        write(99,*) dble(i)*(En(2)-En(1)) , -aimag(trace(Pi_retarded,nm_dev))
    enddo
    close(99) 
    deallocate(Pi_lesser,Pi_greater,Pi_retarded,M)
  endif
  deallocate(cur,tot_cur,tot_ecur)
  deallocate(Ispec,Itot)  
end subroutine green_solve_gw_ephoton_1D

! driver for solving the e-photon SCBA
subroutine green_solve_ephoton_freespace_1D(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,&
  Pmn,polarization,intensity,hw,labs,&
  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, Lx,spindeg,alpha_mix 
real(8), intent(in) :: hw ! hw is photon energy in eV
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: Pmn(nm_dev,nm_dev,3) ! momentum matrix [eV] (multiplied by light-speed, Pmn=c0*p)
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
real(8), intent(in) :: polarization(3) ! light polarization vector 
real(8), intent(in) :: intensity ! [W/m^2]
logical, intent(in) :: labs ! whether to calculate Pi and absorption
real(8), parameter :: pre_fact=((hbar/m0)**2)/(2.0d0*eps0*c0**3) 
!---------
complex(8),dimension(:,:),allocatable ::  Pi_retarded,Pi_lesser,Pi_greater
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:)
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
integer::ie,nop,i,iter,j
complex(8),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable::M(:,:) ! e-photon coupling matrix
real(8)::Nphot,mu(2)
  print *,'====== green_solve_ephoton_freespace_1D ======'
  allocate(tot_cur(nm_dev,nm_dev))
  allocate(tot_ecur(nm_dev,nm_dev))
  allocate(cur(nm_dev,nm_dev,nen))
  allocate(Ispec(nm_dev,nm_dev,nen))
  allocate(Itot(nm_dev,nm_dev))  
  allocate(Pi_retarded(nm_dev,nm_dev))
  allocate(Pi_lesser(nm_dev,nm_dev))
  allocate(Pi_greater(nm_dev,nm_dev))
  mu=(/ mus, mud /)
  allocate(siglead(NB*NS,NB*NS,nen,2))  
  ! get leads sigma
  siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
  siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)  
  do i=1,NB*NS
    do j=1,NB*NS
       if(i.ne.j) then
         siglead(i,j,:,:)=dcmplx(dble(siglead(i,j,:,:)),0.d0*aimag(siglead(i,j,:,:)))
       endif
    enddo
  enddo
  allocate(M(nm_dev,nm_dev))
  M=dcmplx(0.0d0,0.0d0)
  print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
  do i=1,3
    M=M+ polarization(i) * Pmn(:,:,i) 
  enddo    
  !print *, 'pre_fact=', pre_fact
  print '(a8,f15.4,a8,e15.4)','hw=',hw,'I=',intensity  
  nop=floor(hw / (En(2)-En(1)))
  print *,'nop=',nop
  print '(a15,3f8.2)','polarization=',polarization
  print '(a8,f15.4)','dE(meV)=',(En(2)-En(1))*1.0d3
  do iter=0,niter
    ! empty files for sancho 
!    open(unit=101,file='sancho_gbb.dat',status='unknown')
!    close(101)
!    open(unit=101,file='sancho_g00.dat',status='unknown')
!    close(101)
!    open(unit=101,file='sancho_sig.dat',status='unknown')
!    close(101)
    print *, 'calc G'  
    call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham,H00lead,H10lead,Siglead,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,mu=mu,temp=(/temps,tempd/))    
    call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
    call write_current_spectrum('eph_Jdens',iter,cur,nen,en,length,NB,Lx)
    call write_current('eph_I',iter,tot_cur,length,NB,1,Lx)
    call write_current('eph_EI',iter,tot_ecur,length,NB,1,Lx)
    call write_spectrum('eph_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
    call write_spectrum('eph_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
    call write_spectrum('eph_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))
    !
    print *, 'calc Sig'
    call calc_sigma_ephoton_monochromatic(nm_dev,length,nen,En,nop,M,G_lesser,G_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
    Sig_retarded_new=Sig_retarded_new*pre_fact*intensity/hw**2
    Sig_greater_new=Sig_greater_new*pre_fact*intensity/hw**2
    Sig_lesser_new=Sig_lesser_new*pre_fact*intensity/hw**2
    ! mixing with the previous one
    Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
    Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
    Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
    ! get leads sigma
    siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
    siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)  
!    call write_spectrum('SigR',iter,Sig_retarded,nen,En,length,NB,Lx,(/1.0,1.0/))
!    call write_spectrum('SigL',iter,Sig_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
!    call write_spectrum('SigG',iter,Sig_greater,nen,En,length,NB,Lx,(/1.0,1.0/))
    ! calculate collision integral
    call calc_collision(Sig_lesser_new,Sig_greater_new,G_lesser,G_greater,nen,en,spindeg,nm_dev,Itot,Ispec)
    call write_spectrum('eph_Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
    ! calculate absorption
    if (labs) then
      print *, 'calc Pi'
      open(unit=99,file='eph_totabs'//TRIM(STRING(iter))//'.dat',status='unknown')
      do i=1,floor(4.0d0/(En(2)-En(1)))
          call calc_pi_ephoton_monochromatic(nm_dev,length,nen,En,i,M,G_lesser,G_greater,Pi_retarded,Pi_lesser,Pi_greater)  
          call write_trace('eph_absorp',iter,Pi_retarded,length,NB,Lx,(/1.0d0,-1.0d0/),E=dble(i)*(En(2)-En(1)))
          write(99,*) dble(i)*(En(2)-En(1)) , -aimag(trace(Pi_retarded,nm_dev))
      enddo
      close(99)
    endif
  enddo      
  deallocate(Pi_lesser,Pi_greater,Pi_retarded)
  deallocate(M,siglead)
  deallocate(cur,tot_cur,tot_ecur)
  deallocate(Ispec,Itot)
end subroutine green_solve_ephoton_freespace_1D

! calculate e-photon self-energies in the monochromatic assumption
subroutine calc_sigma_ephoton_monochromatic(nm_dev,length,nen,En,nop,M,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater)
integer,intent(in)::nm_dev,length,nen,nop
real(8),intent(in)::en(nen)
complex(8),intent(in),dimension(nm_dev,nm_dev)::M ! e-photon interaction matrix
complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen)::Sig_retarded,Sig_lesser,Sig_greater
!---------
integer::ie
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
  Sig_lesser=0.0d0
  Sig_greater=0.0d0
  Sig_retarded=0.0d0  
  ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M
  !           ~ M [ G^<>(E -+ hw) + G^<>(E +- hw)] M * N
  !$omp parallel default(none) private(ie,A,B) shared(nop,nen,nm_dev,G_lesser,G_greater,Sig_lesser,Sig_greater,M)
  allocate(B(nm_dev,nm_dev))
  allocate(A(nm_dev,nm_dev))  
  !$omp do
  do ie=1,nen
    ! Sig^<(E)
    A = 0.0d0
    if (ie-nop>=1) A =A+ G_lesser(:,:,ie-nop)
    if (ie+nop<=nen) A =A+ G_lesser(:,:,ie+nop)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,A,nm_dev,czero,B,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev)     
    Sig_lesser(:,:,ie) = A    
    ! Sig^>(E)
    A = 0.0d0
    if (ie-nop>=1) A =A+ G_greater(:,:,ie-nop)
    if (ie+nop<=nen) A =A+ G_greater(:,:,ie+nop)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,A,nm_dev,czero,B,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev)     
    Sig_greater(:,:,ie) = A    
  enddo  
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
end subroutine calc_sigma_ephoton_monochromatic


! calculate e-photon polarization of interacting electron-hole pair
subroutine calc_pi_gw_ephoton(spindeg,nm_dev,length,nen,En,nop,M,G_lesser,G_greater,G_retarded,V,Pi_retarded,Pi_lesser,Pi_greater)
integer,intent(in)::nm_dev,length,nen,nop
real(8),intent(in)::en(nen),spindeg
complex(8),intent(in),dimension(nm_dev,nm_dev)::M ! e-photon interaction matrix
complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::G_lesser,G_greater,G_retarded
complex(8),intent(in),dimension(nm_dev,nm_dev)::V ! bare Coulomb matrix
complex(8),intent(inout),dimension(nm_dev,nm_dev)::Pi_retarded,Pi_lesser,Pi_greater
!---------
integer::ie,i
complex(8),allocatable::B(:,:),A(:,:),PR(:,:),WR(:,:) ! tmp matrix
complex(8)::dE
  Pi_lesser=0.0d0
  Pi_greater=0.0d0
  Pi_retarded=0.0d0  
  dE= dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi ) * spindeg
  ! independent quasi-particle polarization
  ! Pi^<>(hw) = \Sum_E M G^<>(E) M G^><(E - hw) M
  !$omp parallel default(none) private(ie,A,B) shared(nop,nen,nm_dev,G_lesser,G_greater,Pi_lesser,Pi_greater,M)
  allocate(B(nm_dev,nm_dev))
  allocate(A(nm_dev,nm_dev))  
  !$omp do
  do ie=1,nen
    if ((ie-nop>=1).and.(ie-nop<=nen)) then
      ! Pi^<(hw) = \sum_E M G<(E) M G>(E-hw)
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,G_lesser(:,:,ie),nm_dev,czero,B,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,A,nm_dev,G_greater(:,:,ie-nop),nm_dev,cone,Pi_lesser,nm_dev)         
      ! Pi^>(hw) = \sum_E M G>(E) M G<(E-hw)   
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,G_greater(:,:,ie),nm_dev,czero,B,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,A,nm_dev,G_lesser(:,:,ie-nop),nm_dev,cone,Pi_greater,nm_dev)             
    endif
  enddo  
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Pi_lesser=Pi_lesser*dE
  Pi_greater=Pi_greater*dE
  !!
  ! calculate PR
  allocate(PR(nm_dev,nm_dev))
  allocate(WR(nm_dev,nm_dev)) 
  allocate(B(nm_dev,nm_dev))
  allocate(A(nm_dev,nm_dev)) 
  PR = czero ! PR
  !$omp parallel default(none) private(i,ie) shared(nop,nen,PR,nm_dev,G_lesser,G_retarded)  
  !$omp do
  do i = 1, nm_dev        
    do ie = max(nop+1,1),min(nen,nen+nop)                                                 
      PR(i,i) = PR(i,i) + G_lesser(i,i,ie) * conjg(G_retarded(i,i,ie-nop)) +&
                        & G_retarded(i,i,ie) * G_lesser(i,i,ie-nop)        
    enddo
  enddo    
  !$omp end do
  !$omp end parallel
  dE = dcmplx(0.0d0, -( En(2) - En(1) ) / 2.0d0 / pi * spindeg )
  PR=-PR*dE 
  !
  A = czero ! PR static
  !$omp parallel default(none) private(i,ie) shared(nop,nen,A,nm_dev,G_lesser,G_retarded)  
  !$omp do
  do i = 1, nm_dev        
    do ie = 1,nen                                                
      A(i,i) = A(i,i) + G_lesser(i,i,ie) * conjg(G_retarded(i,i,ie)) +&
                      & G_retarded(i,i,ie) * G_lesser(i,i,ie)        
    enddo
  enddo    
  !$omp end do
  !$omp end parallel
  A=A*dE
  ! calculate WR static
  !! B = I - V P^r
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,A,nm_dev,czero,B,nm_dev)
  !
  do i=1,nm_dev
     B(i,i) = 1.0d0 + B(i,i)
  enddo    
  !!!! calculate W^r = (I - V P^r)^-1 V    
  call invert(B,nm_dev) 
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V,nm_dev,czero,WR,nm_dev) 
  ! diagonal approximation BSE
  do i=1,nm_dev
    Pi_lesser(i,i)=Pi_lesser(i,i) / (1.0d0 - WR(i,i)*PR(i,i)* dcmplx(0.0d0,1.0d0)) 
    Pi_greater(i,i)=Pi_greater(i,i) / (1.0d0 - conjg(WR(i,i))*conjg(PR(i,i))* dcmplx(0.0d0,1.0d0)) 
  enddo
  !
  Pi_retarded = dcmplx(0.0d0*dble(Pi_retarded),aimag(Pi_greater-Pi_lesser)/2.0d0)
  deallocate(PR)
  deallocate(A,B,WR)
end subroutine calc_pi_gw_ephoton


! calculate e-photon polarization self-energies in the monochromatic assumption, for independent electron-hole pair
subroutine calc_pi_ephoton_monochromatic(nm_dev,length,nen,En,nop,M,G_lesser,G_greater,Pi_retarded,Pi_lesser,Pi_greater)
integer,intent(in)::nm_dev,length,nen,nop
real(8),intent(in)::en(nen)
complex(8),intent(in),dimension(nm_dev,nm_dev)::M ! e-photon interaction matrix
complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm_dev,nm_dev)::Pi_retarded,Pi_lesser,Pi_greater
!---------
integer::ie
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix
complex(8)::dE
  Pi_lesser=0.0d0
  Pi_greater=0.0d0
  Pi_retarded=0.0d0  
  dE= dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi ) 
  ! Pi^<>(hw) = \Sum_E M G^<>(E) M G^><(E - hw) M
  !$omp parallel default(none) private(ie,A,B) shared(nop,nen,nm_dev,G_lesser,G_greater,Pi_lesser,Pi_greater,M)
  allocate(B(nm_dev,nm_dev))
  allocate(A(nm_dev,nm_dev))  
  !$omp do
  do ie=1,nen
    if ((ie-nop>=1).and.(ie-nop<=nen)) then
      ! Pi^<(hw) = \sum_E M G<(E) M G>(E-hw)
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,G_lesser(:,:,ie),nm_dev,czero,B,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,A,nm_dev,G_greater(:,:,ie-nop),nm_dev,cone,Pi_lesser,nm_dev)         
      ! Pi^>(hw) = \sum_E M G>(E) M G<(E-hw)   
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,G_greater(:,:,ie),nm_dev,czero,B,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,A,nm_dev,G_lesser(:,:,ie-nop),nm_dev,cone,Pi_greater,nm_dev)             
    endif
  enddo  
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Pi_lesser=Pi_lesser*dE
  Pi_greater=Pi_greater*dE
  Pi_retarded = dcmplx(0.0d0*dble(Pi_retarded),aimag(Pi_greater-Pi_lesser)/2.0d0)
end subroutine calc_pi_ephoton_monochromatic


! calculate e-photon self-energies in spontaneous emission and ADD onto the Sig_r<>
subroutine calc_sigma_ephoton_monochromatic_spontaneous_emission(nm_dev,length,nen,En,nop,M,prefactor,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater)
integer,intent(in)::nm_dev,length,nen,nop
real(8),intent(in)::en(nen)
complex(8),intent(in)::prefactor
complex(8),intent(in),dimension(nm_dev,nm_dev)::M ! e-photon interaction matrix
complex(8),intent(in),dimension(nm_dev,nm_dev,nen)::G_lesser,G_greater
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen)::Sig_retarded,Sig_lesser,Sig_greater
!---------
integer::ie
complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix  
  ! Sig^<>(E) = M G^<>(E +- hw) M  
  !$omp parallel default(none) private(ie,A,B) shared(nop,nen,nm_dev,G_lesser,G_greater,Sig_lesser,Sig_greater,M,prefactor)
  allocate(B(nm_dev,nm_dev))
  allocate(A(nm_dev,nm_dev))  
  !$omp do
  do ie=1,nen
    ! Sig^<(E)    
    if (ie+nop<=nen) then 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,G_lesser(:,:,ie+nop),nm_dev,czero,B,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev)     
      Sig_lesser(:,:,ie) =Sig_lesser(:,:,ie)+ A(:,:) * prefactor
    endif    
    ! Sig^>(E)    
    if (ie-nop>=1) then 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,G_greater(:,:,ie-nop),nm_dev,czero,B,nm_dev) 
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,A,nm_dev)     
      Sig_greater(:,:,ie) =Sig_greater(:,:,ie)+ A(:,:) * prefactor   
    endif
  enddo  
  !$omp end do
  deallocate(A,B)
  !$omp end parallel
  Sig_retarded = dcmplx(0.0d0*dble(Sig_retarded),aimag(Sig_greater-Sig_lesser)/2.0d0)
end subroutine calc_sigma_ephoton_monochromatic_spontaneous_emission


! memory saving version of 2D GW solver with one periodic direction (z)
! iterating G -> P -> W -> Sig 
subroutine green_solve_gw_2D_memsaving(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,nphiz,Ham,H00lead,H10lead,T,V,&
  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,encut,Egap,ndiagmin,writeGF)
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length, nphiz
integer, intent(in), optional :: ndiagmin
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg,Egap
complex(8),intent(in) :: Ham(nm_dev,nm_dev,nphiz),H00lead(NB*NS,NB*NS,2,nphiz),H10lead(NB*NS,NB*NS,2,nphiz),T(NB*NS,nm_dev,2,nphiz)
complex(8), intent(in):: V(nm_dev,nm_dev,nphiz)
logical,intent(in)::ldiag
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen,nphiz) ::  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
real(8),intent(in)::encut(2) ! intraband and interband cutoff for P
logical, intent(in), optional :: writeGF 
!------
complex(8),allocatable,dimension(:,:) ::  P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater
complex(8),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),cur_k(:,:,:,:),sumtot_cur(:,:),sumtot_ecur(:,:)
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
logical :: lwriteGF
real(8),allocatable::wen(:) ! energy vector for P and W
integer,allocatable::nops(:) ! discretized energy for P and W
real(8),allocatable::Tr(:,:) ! current spectrum on leads
real(8),allocatable::Te(:,:,:) ! transmission matrix spectrum
integer :: iter,ie,iop,nnop,nnop1,nnop2
integer :: i,j,nm,nop,l,h,ndiag,ikz,iqz,ikzd
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
if (present(writeGF)) then
  lwriteGF=writeGF
else
  lwriteGF=.false.
endif
print *,'================= green_solve_gw_2D_memsaving ================='
allocate(B(nm_dev,nm_dev))
allocate(tot_cur(nm_dev,nm_dev))
allocate(tot_ecur(nm_dev,nm_dev))
allocate(sumtot_cur(nm_dev,nm_dev))
allocate(sumtot_ecur(nm_dev,nm_dev))
allocate(cur(nm_dev,nm_dev,nen))
allocate(cur_k(nm_dev,nm_dev,nen,nphiz))
allocate(Ispec(nm_dev,nm_dev,nen))
allocate(Itot(nm_dev,nm_dev))
!
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
print *,'Nkz=',nphiz
print *,'NEn=',nen
print *,'ND=',nm_dev
print *,'dE=',(en(2)-en(1))
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
print *, ' Encut: intra    inter    Eg (eV)' 
print '(A6,3F8.3)',' ',encut,egap
!print *, ' Nop='
!print '(10I5)',nops
print *, ' w= (eV)'
print '(6F8.3)',wen
print *, '--------------------------------------------------------------'
!
allocate(siglead(NB*NS,NB*NS,nen,2,nphiz))
! get leads sigma
do ikz=1,nphiz
  siglead(:,:,:,1,ikz) = Sig_retarded(1:NB*NS,1:NB*NS,:,ikz)
  siglead(:,:,:,2,ikz) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:,ikz)  
enddo
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
  sumtot_cur=0.0d0
  sumtot_ecur=0.0d0
  do ikz=1,nphiz
    print *, ' ikz=', ikz
    call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),G_retarded(:,:,:,ikz),G_lesser(:,:,:,ikz),G_greater(:,:,:,ikz),mu=mu,temp=(/temps,tempd/))
    call write_spectrum('gw_ldos_kz'//string(ikz)//'_',iter,G_retarded(:,:,:,ikz),nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
    call calc_bond_current(Ham(:,:,ikz),G_lesser(:,:,:,ikz),nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
    cur_k(:,:,:,ikz)=cur
    call write_current_spectrum('gw_Jdens_kz'//string(ikz)//'_',iter,cur,nen,en,length,NB,Lx)    
    sumtot_cur=sumtot_cur+tot_cur
    sumtot_ecur=sumtot_ecur+tot_ecur
  enddo
  call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded,nen,En,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
  call write_spectrum_summed_over_kz('gw_ndos',iter,G_lesser,nen,En,nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum_summed_over_kz('gw_pdos',iter,G_greater,nen,En,nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
  call write_current_spectrum_summed_over_kz('gw_Jdens',iter,cur_k,nen,En,nphiz,length,NB,Lx)
  call write_current('gw_I',iter,sumtot_cur,length,NB,NS,Lx)
  call write_current('gw_EI',iter,sumtot_ecur,length,NB,NS,Lx)
  !
  ! empty sigma_x_new matrices for accumulation
  sig_retarded_new=czero
  sig_lesser_new=czero
  sig_greater_new=czero
  print *, 'calc P, solve W, add to Sigma_new'     
  ndiag=NB*(min(NS,iter))
  if (ldiag) ndiag=0  
  if (present(ndiagmin)) ndiag=max(ndiagmin,ndiag)
  if (lwriteGF) ndiag=nm_dev
  print *,'ndiag=',min(ndiag,nm_dev)
  !
  !print *,'   i / n :  Nop   Eop (eV)'
  do iop=1,nnop        
    !print '(I5,A,I5,A,I5,F8.3)',iop,'/',nnop,':',nops(iop),wen(iop)  
    nop=nops(iop)
    do iqz=1,nphiz
      !print *, ' iqz=', iqz
      P_lesser=czero
      P_greater=czero
      P_retarded=czero
      do ie = max(nop+1,1),min(nen,nen+nop) 
        do ikz=1,nphiz          
          ikzd=ikz-iqz + nphiz/2
          if (ikzd<1) ikzd=ikzd+nphiz
          if (ikzd>nphiz) ikzd=ikzd-nphiz
          !$omp parallel default(none) private(l,h,i) shared(ndiag,P_lesser,P_greater,P_retarded,nm_dev,G_lesser,G_greater,G_retarded,nop,ie,ikz,ikzd)
          !$omp do
          do i = 1, nm_dev        
              l=max(i-ndiag,1)
              h=min(nm_dev,i+ndiag)              
              P_lesser(i,l:h) = P_lesser(i,l:h) + G_lesser(i,l:h,ie,ikz) * G_greater(l:h,i,ie-nop,ikzd)
              P_greater(i,l:h) = P_greater(i,l:h) + G_greater(i,l:h,ie,ikz) * G_lesser(l:h,i,ie-nop,ikzd)        
              P_retarded(i,l:h) = P_retarded(i,l:h) + &
                  & (G_lesser(i,l:h,ie,ikz) * conjg(G_retarded(i,l:h,ie-nop,ikzd)) + G_retarded(i,l:h,ie,ikz) * G_lesser(l:h,i,ie-nop,ikzd))        
          enddo
          !$omp end do
          !$omp end parallel
        enddo
      enddo
      dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi ) * spindeg  
      P_lesser=dE*P_lesser
      P_greater=dE*P_greater
      P_retarded=dE*P_retarded
      ! calculate W
      call green_calc_w(0,NB,NS,nm_dev,P_retarded,P_lesser,P_greater,V,W_retarded,W_lesser,W_greater)
      !      
      if (lwriteGF) then
        call write_matrix('W_r',0,W_retarded(:,:),wen(iop),length,NB,(/1.0d0,1.0d0/))
      endif
      !
      ! Accumulate the GW to Sigma
      ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)  
      do ie=1,nen
        do ikz=1,nphiz          
          if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then                
            ikzd=ikz-iqz + nphiz/2            
            if (ikzd<1) ikzd=ikzd+nphiz
            if (ikzd>nphiz) ikzd=ikzd-nphiz
            !$omp parallel default(none) private(l,h,i) shared(ndiag,nop,Sig_lesser_new,Sig_greater_new,Sig_retarded_new,W_lesser,W_greater,W_retarded,nm_dev,G_lesser,G_greater,G_retarded,ie,ikz,ikzd)  
            !$omp do
            do i = 1,nm_dev   
              l=max(i-ndiag,1)
              h=min(nm_dev,i+ndiag)       
              Sig_lesser_new(i,l:h,ie,ikz)=Sig_lesser_new(i,l:h,ie,ikz)+G_lesser(i,l:h,ie-nop,ikzd)*W_lesser(i,l:h)
              Sig_greater_new(i,l:h,ie,ikz)=Sig_greater_new(i,l:h,ie,ikz)+G_greater(i,l:h,ie-nop,ikzd)*W_greater(i,l:h)
              Sig_retarded_new(i,l:h,ie,ikz)=Sig_retarded_new(i,l:h,ie,ikz)+G_lesser(i,l:h,ie-nop,ikzd)*W_retarded(i,l:h) 
              Sig_retarded_new(i,l:h,ie,ikz)=Sig_retarded_new(i,l:h,ie,ikz)+G_retarded(i,l:h,ie-nop,ikzd)*W_lesser(i,l:h) 
              Sig_retarded_new(i,l:h,ie,ikz)=Sig_retarded_new(i,l:h,ie,ikz)+G_retarded(i,l:h,ie-nop,ikzd)*W_retarded(i,l:h)          
            enddo      
            !$omp end do
            !$omp end parallel           
          endif
        enddo
      enddo
    enddo
  enddo  
  dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)  
  Sig_lesser_new = Sig_lesser_new  * dE
  Sig_greater_new= Sig_greater_new * dE
  Sig_retarded_new=Sig_retarded_new* dE
  Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
  !!! Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
  !!! Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
  !
  ! symmetrize the selfenergies
  do ie=1,nen
    do ikz=1,nphiz
      B(:,:)=transpose(Sig_retarded_new(:,:,ie,ikz))
      Sig_retarded_new(:,:,ie,ikz) = (Sig_retarded_new(:,:,ie,ikz) + B(:,:))/2.0d0    
      B(:,:)=transpose(Sig_lesser_new(:,:,ie,ikz))
      Sig_lesser_new(:,:,ie,ikz) = (Sig_lesser_new(:,:,ie,ikz) + B(:,:))/2.0d0
      B(:,:)=transpose(Sig_greater_new(:,:,ie,ikz))
      Sig_greater_new(:,:,ie,ikz) = (Sig_greater_new(:,:,ie,ikz) + B(:,:))/2.0d0
    enddo
  enddo
  ! mixing with previous ones
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
  Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
  ! get leads sigma
  do iqz=1,nphiz
    siglead(:,:,:,1,iqz) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:,iqz)
    siglead(:,:,:,2,iqz) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:,iqz)    
  enddo
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    do iqz=1,nphiz
      call expand_size_bycopy(Sig_retarded(:,:,ie,iqz),nm_dev,NB,2)
      call expand_size_bycopy(Sig_lesser(:,:,ie,iqz),nm_dev,NB,2)
      call expand_size_bycopy(Sig_greater(:,:,ie,iqz),nm_dev,NB,2)
    enddo
  enddo
  call write_spectrum_summed_over_kz('gw_SigR',iter,Sig_retarded,nen,En,nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
!  call write_spectrum_summed_over_kz('gw_SigL',iter,Sig_lesser,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('gw_SigG',iter,Sig_greater,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
end do  
deallocate(siglead,B)
deallocate(cur_k,cur,tot_cur,tot_ecur,sumtot_cur,sumtot_ecur)
deallocate(Ispec,Itot)
deallocate(P_retarded,P_lesser,P_greater)
deallocate(W_retarded,W_lesser,W_greater)
deallocate(wen,nops)
end subroutine green_solve_gw_2D_memsaving


! 2D GW solver with one periodic direction (z)
! iterating G -> P -> W -> Sig 
subroutine green_solve_gw_2D(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,nphiz,Ham,H00lead,H10lead,T,V,&
  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
  W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag)
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length, nphiz
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg
complex(8),intent(in) :: Ham(nm_dev,nm_dev,nphiz),H00lead(NB*NS,NB*NS,2,nphiz),H10lead(NB*NS,NB*NS,2,nphiz),T(NB*NS,nm_dev,2,nphiz)
complex(8), intent(in):: V(nm_dev,nm_dev,nphiz)
logical,intent(in)::ldiag
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen,nphiz) ::  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
!------
complex(8),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:),cur_k(:,:,:,:),sumtot_cur(:,:),sumtot_ecur(:,:)
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
integer :: iter,ie,nopmax
integer :: i,j,nm,nop,l,h,iop,ndiag,ikz,iqz,ikzd
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
print *,'====== green_solve_gw_2D ======'
allocate(siglead(NB*NS,NB*NS,nen,2,nphiz))
siglead=dcmplx(0.0d0,0.0d0)
allocate(B(nm_dev,nm_dev))
allocate(tot_cur(nm_dev,nm_dev))
allocate(tot_ecur(nm_dev,nm_dev))
allocate(sumtot_cur(nm_dev,nm_dev))
allocate(sumtot_ecur(nm_dev,nm_dev))
allocate(cur(nm_dev,nm_dev,nen))
allocate(cur_k(nm_dev,nm_dev,nen,nphiz))
allocate(Ispec(nm_dev,nm_dev,nen))
allocate(Itot(nm_dev,nm_dev))
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
do iter=0,niter
  ! empty files for sancho 
!  open(unit=101,file='sancho_gbb.dat',status='unknown')
!  close(101)
!  open(unit=101,file='sancho_g00.dat',status='unknown')
!  close(101)
!  open(unit=101,file='sancho_sig.dat',status='unknown')
!  close(101)
  print *,'+ iter=',iter
  print *, 'calc G'  
  sumtot_cur=0.0d0
  sumtot_ecur=0.0d0
  do ikz=1,nphiz
    print *, ' ikz=', ikz
    call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),G_retarded(:,:,:,ikz),G_lesser(:,:,:,ikz),G_greater(:,:,:,ikz),mu=mu,temp=(/temps,tempd/))
    call write_spectrum('ldos_kz'//string(ikz)//'_',iter,G_retarded(:,:,:,ikz),nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
    call calc_bond_current(Ham(:,:,ikz),G_lesser(:,:,:,ikz),nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
    cur_k(:,:,:,ikz)=cur
    call write_current_spectrum('Jdens_kz'//string(ikz)//'_',iter,cur,nen,en,length,NB,Lx)    
    sumtot_cur=sumtot_cur+tot_cur
    sumtot_ecur=sumtot_ecur+tot_ecur
  enddo
  call write_spectrum_summed_over_kz('ldos',iter,G_retarded,nen,En,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
  call write_spectrum_summed_over_kz('ndos',iter,G_lesser,nen,En,nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum_summed_over_kz('pdos',iter,G_greater,nen,En,nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
  call write_current_spectrum_summed_over_kz('Jdens_',iter,cur_k,nen,En,nphiz,length,NB,Lx)
  call write_current('I',iter,sumtot_cur,length,NB,NS,Lx)
  call write_current('EI',iter,sumtot_ecur,length,NB,NS,Lx)
  G_retarded=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
  G_lesser=dcmplx(0.0d0*dble(G_lesser),aimag(G_lesser))
  G_greater=dcmplx(0.0d0*dble(G_greater),aimag(G_greater))
  !        
  print *, 'calc P'
  !
  nopmax=nen/2-10  
  ndiag=NB*(min(NS*2,iter))
  if (ldiag) ndiag=0  
  print *,'ndiag=',min(ndiag,nm_dev)
  ! Pij^<>(hw,kz') = \int_dE Gij^<>(E,kz) * Gji^><(E-hw,kz-kz')
  ! Pij^r(hw,kz')  = \int_dE Gij^<(E,kz) * Gji^a(E-hw,kz-kz') + Gij^r(E,kz) * Gji^<(E-hw,kz-kz')
  !$omp parallel default(none) private(l,h,iop,nop,ie,ikz,ikzd,iqz,dE,i,j) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded,nphiz)
  !$omp do
  do nop=-nopmax,nopmax
    do iqz=1,nphiz
      iop=nop+nen/2
      P_lesser(:,:,iop,iqz) = dcmplx(0.0d0,0.0d0)
      P_greater(:,:,iop,iqz) = dcmplx(0.0d0,0.0d0)    
      P_retarded(:,:,iop,iqz) = dcmplx(0.0d0,0.0d0)    
      do ie = max(nop+1,1),min(nen,nen+nop) 
        do ikz=1,nphiz
          dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )   
          do i = 1, nm_dev        
              l=max(i-ndiag,1)
              h=min(nm_dev,i+ndiag)
              ikzd=ikz-iqz + nphiz/2
              if (ikzd<1) ikzd=ikzd+nphiz
              if (ikzd>nphiz) ikzd=ikzd-nphiz
              P_lesser(i,l:h,iop,iqz) = P_lesser(i,l:h,iop,iqz) + G_lesser(i,l:h,ie,ikz) * G_greater(l:h,i,ie-nop,ikzd)
              P_greater(i,l:h,iop,iqz) = P_greater(i,l:h,iop,iqz) + G_greater(i,l:h,ie,ikz) * G_lesser(l:h,i,ie-nop,ikzd)        
              P_retarded(i,l:h,iop,iqz) = P_retarded(i,l:h,iop,iqz) + (G_lesser(i,l:h,ie,ikz) * conjg(G_retarded(i,l:h,ie-nop,ikzd)) + G_retarded(i,l:h,ie,ikz) * G_lesser(l:h,i,ie-nop,ikzd))        
          enddo
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi ) * spindeg  
  P_lesser=dE*P_lesser
  P_greater=dE*P_greater
  P_retarded=dE*P_retarded
  call write_spectrum_summed_over_kz('PR',iter,P_retarded,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
!  call write_spectrum_summed_over_kz('PL',iter,P_lesser  ,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('PG',iter,P_greater ,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
  !
  print *, 'calc W'
  !
  do nop=-nopmax+nen/2,nopmax+nen/2   
    do iqz=1,nphiz
      call green_calc_w(2,NB,NS,nm_dev,P_retarded(:,:,nop,iqz),P_lesser(:,:,nop,iqz),P_greater(:,:,nop,iqz),V,W_retarded(:,:,nop,iqz),W_lesser(:,:,nop,iqz),W_greater(:,:,nop,iqz))
    enddo
  enddo
  call write_spectrum_summed_over_kz('WR',iter,W_retarded,nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
!  call write_spectrum_summed_over_kz('WL',iter,W_lesser,  nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('WG',iter,W_greater, nen,En-en(nen/2),nphiz,length,NB,Lx,(/1.0,1.0/))
  !
  print *, 'calc SigGW'
  !  
  ndiag=NB*(min(NS*2,iter))
  if (ldiag) ndiag=0  
  print *,'ndiag=',min(ndiag,nm_dev)
  nopmax=nen/2-10
  Sig_greater_new = dcmplx(0.0d0,0.0d0)
  Sig_lesser_new = dcmplx(0.0d0,0.0d0)
  Sig_retarded_new = dcmplx(0.0d0,0.0d0)      
  ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)
  !$omp parallel default(none) private(l,h,nop,ie,i,j,iop,ikz,ikzd,iqz) shared(ndiag,nopmax,Sig_lesser_new,Sig_greater_new,Sig_retarded_new,W_lesser,W_greater,W_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded,dE,nphiz)
  !$omp do
  do ie=1,nen
    do ikz=1,nphiz
      do nop= -nopmax,nopmax    
        if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then      
          do iqz=1,nphiz
            iop=nop+nen/2
            ikzd=ikz-iqz + nphiz/2            
            if (ikzd<1) ikzd=ikzd+nphiz
            if (ikzd>nphiz) ikzd=ikzd-nphiz
            do i = 1,nm_dev   
              l=max(i-ndiag,1)
              h=min(nm_dev,i+ndiag)       
              Sig_lesser_new(i,l:h,ie,ikz)=Sig_lesser_new(i,l:h,ie,ikz)+G_lesser(i,l:h,ie-nop,ikzd)*W_lesser(i,l:h,iop,iqz)
              Sig_greater_new(i,l:h,ie,ikz)=Sig_greater_new(i,l:h,ie,ikz)+G_greater(i,l:h,ie-nop,ikzd)*W_greater(i,l:h,iop,iqz)
              Sig_retarded_new(i,l:h,ie,ikz)=Sig_retarded_new(i,l:h,ie,ikz)+G_lesser(i,l:h,ie-nop,ikzd)*W_retarded(i,l:h,iop,iqz) &
                                             +G_retarded(i,l:h,ie-nop,ikzd)*W_lesser(i,l:h,iop,iqz) &
                                             +G_retarded(i,l:h,ie-nop,ikzd)*W_retarded(i,l:h,iop,iqz)          
            enddo      
          enddo
        endif
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)
  Sig_lesser_new = Sig_lesser_new  * dE
  Sig_greater_new= Sig_greater_new * dE
  Sig_retarded_new=Sig_retarded_new* dE
  Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
  !!! Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
  !!! Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
  !
  ! symmetrize the selfenergies
  do ie=1,nen
    do ikz=1,nphiz
      B(:,:)=transpose(Sig_retarded_new(:,:,ie,ikz))
      Sig_retarded_new(:,:,ie,ikz) = (Sig_retarded_new(:,:,ie,ikz) + B(:,:))/2.0d0    
      B(:,:)=transpose(Sig_lesser_new(:,:,ie,ikz))
      Sig_lesser_new(:,:,ie,ikz) = (Sig_lesser_new(:,:,ie,ikz) + B(:,:))/2.0d0
      B(:,:)=transpose(Sig_greater_new(:,:,ie,ikz))
      Sig_greater_new(:,:,ie,ikz) = (Sig_greater_new(:,:,ie,ikz) + B(:,:))/2.0d0
    enddo
  enddo
  ! mixing with previous ones
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
  Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
  ! get leads sigma
  do iqz=1,nphiz
    siglead(:,:,:,1,iqz) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:,iqz)
    siglead(:,:,:,2,iqz) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:,iqz)    
  enddo
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    do iqz=1,nphiz
      call expand_size_bycopy(Sig_retarded(:,:,ie,iqz),nm_dev,NB,2)
      call expand_size_bycopy(Sig_lesser(:,:,ie,iqz),nm_dev,NB,2)
      call expand_size_bycopy(Sig_greater(:,:,ie,iqz),nm_dev,NB,2)
    enddo
  enddo
  call write_spectrum_summed_over_kz('SigR',iter,Sig_retarded,nen,En,nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
!  call write_spectrum_summed_over_kz('SigL',iter,Sig_lesser,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum_summed_over_kz('SigG',iter,Sig_greater,nen,En,nphiz,length,NB,Lx,(/1.0,1.0/))
end do  
deallocate(siglead,B)
deallocate(cur_k,cur,tot_cur,tot_ecur,sumtot_cur,sumtot_ecur)
deallocate(Ispec,Itot)
end subroutine green_solve_gw_2D


! driver for iterating G -> P -> W -> Sig 
! super memory saving version of green_solve_gw_1D , only keeping the diagonal
! blocks of size (NBxNB) of G and Sigma in memory
subroutine green_solve_gw_1D_supermemsaving(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,midgap,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,V,&
  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,encut,Egap,ndiagmin,writeGF)
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
integer, intent(in), optional :: ndiagmin
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg, Egap,midgap(2)
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: V(nm_dev,nm_dev)
logical,intent(in)::ldiag
real(8),intent(in)::encut(2) ! intraband and interband cutoff for P
complex(8),intent(inout),dimension(NB,NB,length,nen) ::  G_retarded,G_lesser,G_greater,&
     Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
logical, intent(in), optional :: writeGF 
!----
complex(8),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
complex(8),allocatable::cur(:,:,:,:)
real(8),allocatable::tot_cur(:,:),tot_ecur(:,:)
real(8),allocatable::wen(:)  ! energy vector for P and W
integer,allocatable::nops(:) ! discretized energy for P and W
real(8),allocatable::Tr(:,:) ! current spectrum on leads
real(8),allocatable::Te(:,:,:) ! transmission matrix spectrum 
integer :: iter,ie,iop,nnop,nnop1,nnop2
integer :: i,j,nm,l,h,ndiag,nop,ix
logical :: lwriteGF
complex(8),allocatable::Ispec(:,:,:,:),Itot(:,:)
complex(8),allocatable,dimension(:,:) ::  P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
if (present(writeGF)) then
  lwriteGF=writeGF
else
  lwriteGF=.false.
endif
print *,'================= green_solve_gw_1D_supermemsaving ================='
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
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
print *, ' Encut: intra    inter    Eg (eV)' 
print '(A6,3F8.3)',' ',encut,egap
!print *, ' Nop='
!print '(10I5)',nops
print *, ' Eop= (eV)'
print '(6F8.3)',wen
print *, '--------------------------------------------------------------'
!
allocate(siglead(NB*NS,NB*NS,nen,2))
siglead=czero
! get leads sigma
do ie=1,nen
  do i=1,NS
    l=(i-1)*NB+1
    h=i*NB
    siglead(l:h,l:h,ie,1) = Sig_retarded(:,:,i,ie)
    siglead(l:h,l:h,ie,2) = Sig_retarded(:,:,length-NS+i,ie)    
  enddo
enddo
!
allocate(B(nm_dev,nm_dev))
allocate(tot_cur(nm_dev,nm_dev))
allocate(tot_ecur(nm_dev,nm_dev))
allocate(cur(nb,nb,length,nen))
allocate(Ispec(nb,nb,length,nen))
allocate(Itot(nm_dev,nm_dev))
allocate(tr(nen,2))
allocate(te(nen,2,2))
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
  call green_calc_g_block(nen,En,2,nm_dev,nb,length,(/nb*ns,nb*ns/),nb*ns,spindeg,Ham,H00lead,H10lead,Siglead,T,&
    Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,cur,tot_cur,tot_ecur,Itot,Ispec,&
    cur=Tr,te=Te,mu=mu,temp=(/temps,tempd/)) 
  ! 
  call write_current_spectrum_block('gw_Jdens',iter,cur,nen,en,length,NB,Lx)
  call write_current('gw_I',iter,tot_cur,length,NB,NS,Lx)
  call write_current('gw_EI',iter,tot_ecur,length,NB,NS,Lx)
  call write_spectrum_block('gw_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0,-2.0/))
  call write_spectrum_block('gw_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
  call write_spectrum_block('gw_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0,-1.0/))
  call write_transmission_spectrum('gw_trL',iter,Tr(:,1)*spindeg,nen,En)
  call write_transmission_spectrum('gw_trR',iter,Tr(:,2)*spindeg,nen,En)
  call write_transmission_spectrum('gw_TE_LR',iter,Te(:,1,2)*spindeg,nen,En)
  call write_transmission_spectrum('gw_TE_RL',iter,Te(:,2,1)*spindeg,nen,En)
  call write_spectrum_block('gw_Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0,1.0/))
  G_retarded=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
  G_lesser=dcmplx(0.0d0*dble(G_lesser),aimag(G_lesser))
  G_greater=dcmplx(0.0d0*dble(G_greater),aimag(G_greater))
  if (lwriteGF) then
 !   call write_matrix_E('G_r',0,G_retarded,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('G_l',0,G_lesser,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('G_g',0,G_greater,nen,en,length,NB,(/1.0,1.0/))
  endif
  !        
  ! empty sigma_x_new matrices for accumulation
  sig_retarded_new=czero
  sig_lesser_new=czero
  sig_greater_new=czero
  print *, 'calc P, solve W, add to Sigma_new'     
  ndiag=NB*(min(NS,iter))
  if (ldiag) ndiag=0  
  if (present(ndiagmin)) ndiag=max(ndiagmin,ndiag)
  if (lwriteGF) ndiag=nm_dev
  print *,'ndiag=',min(ndiag,nm_dev)
  !
  !print *,'   i / n :  Nop   Eop (eV)'
  do iop=1,nnop        
    !print '(I5,A,I5,A,I5,F8.3)',iop,'/',nnop,':',nops(iop),wen(iop)    
    nop=nops(iop)
    P_lesser = czero
    P_greater = czero
    P_retarded = czero
    !$omp parallel default(none) private(ix,l,h,i,ie) shared(length,NB,NS,ndiag,nop,nen,P_lesser,P_greater,P_retarded,nm_dev,G_lesser,G_greater,G_retarded)  
    !$omp do
    do ix=1,length
      do i = 1, NB        
        do ie = max(nop+1,1),min(nen,nen+nop)                   
          l=max(i-ndiag,1)
          h=min(NB,i+ndiag)                            
          P_lesser(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) = P_lesser(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) &
            + G_lesser(i,l:h,ix,ie) * G_greater(l:h,i,ix,ie-nop)
          P_greater(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) = P_greater(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) &
            + G_greater(i,l:h,ix,ie) * G_lesser(l:h,i,ix,ie-nop) 
          P_retarded(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) = P_retarded(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) &
            + G_lesser(i,l:h,ix,ie) * conjg(G_retarded(i,l:h,ix,ie-nop)) &
            + G_retarded(i,l:h,ix,ie) * G_lesser(l:h,i,ix,ie-nop)        
        enddo
      enddo    
    enddo
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )* spindeg    
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE  
    P_retarded=P_retarded*dE
    if (lwriteGF) then
      !call write_matrix('P_r',0,P_retarded(:,:),wen(iop),length,NB,(/1.0,1.0/))
    endif
    !
    ! calculate W
    call green_calc_w(0,NB,NS,nm_dev,P_retarded,P_lesser,P_greater,V,W_retarded,W_lesser,W_greater)
    !
    if (lwriteGF) then
    !  call write_matrix('W_r',0,W_retarded(:,:),wen(iop),length,NB,(/1.0,1.0/))
    endif
    !
    ! Accumulate the GW to Sigma
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)  
    !$omp parallel default(none) private(ix,l,h,i,ie) shared(NB,NS,length,ndiag,nop,nen,Sig_lesser_new,Sig_greater_new,Sig_retarded_new,W_lesser,W_greater,W_retarded,nm_dev,G_lesser,G_greater,G_retarded)  
    !$omp do
    do ix=1,length
      do i=1,NB
        l=max(i-ndiag,1)
        h=min(NB,i+ndiag)           
        do ie=1,nen
          if ((ie .gt. max(nop,1)).and.(ie .lt. (nen+nop))) then 
            Sig_lesser_new(i,l:h,ix,ie)=Sig_lesser_new(i,l:h,ix,ie)+G_lesser(i,l:h,ix,ie-nop)*W_lesser(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB)
            Sig_greater_new(i,l:h,ix,ie)=Sig_greater_new(i,l:h,ix,ie)+G_greater(i,l:h,ix,ie-nop)*W_greater(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB)
            Sig_retarded_new(i,l:h,ix,ie)=Sig_retarded_new(i,l:h,ix,ie)+G_lesser(i,l:h,ix,ie-nop)*W_retarded(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) + &
                                       G_retarded(i,l:h,ix,ie-nop)*W_lesser(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB) + &
                                       G_retarded(i,l:h,ix,ie-nop)*W_retarded(i+(ix-1)*NB,l+(ix-1)*NB:h+(ix-1)*NB)          
          endif     
        enddo   
      enddo ! i
    enddo ! ix
    !$omp end do
    !$omp end parallel    
  enddo ! nop                               
  !
  dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)                
  Sig_lesser_new = Sig_lesser_new  * dE
  Sig_greater_new= Sig_greater_new * dE
  Sig_retarded_new=Sig_retarded_new* dE
  !
  Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
  ! symmetrize the selfenergies
  do ie=1,nen
    do ix=1,length      
      Sig_retarded_new(:,:,ix,ie) = (Sig_retarded_new(:,:,ix,ie) + transpose(Sig_retarded_new(:,:,ix,ie)))/2.0d0          
      Sig_lesser_new(:,:,ix,ie) = (Sig_lesser_new(:,:,ix,ie) + transpose(Sig_lesser_new(:,:,ix,ie)))/2.0d0      
      Sig_greater_new(:,:,ix,ie) = (Sig_greater_new(:,:,ix,ie) + transpose(Sig_greater_new(:,:,ix,ie)))/2.0d0
    enddo
  enddo  
  !
  if (lwriteGF) then
!    call write_matrix_E('Sigma_r',0,Sig_retarded_new,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('Sigma_l',0,Sig_lesser_new,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('Sigma_g',0,Sig_greater_new,nen,en,length,NB,(/1.0,1.0/))
  endif
  ! mixing with the previous one
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
  Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)    
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    do ix=1,NS
      Sig_retarded(:,:,ix,ie)=Sig_retarded(:,:,NS+1,ie)
      Sig_lesser(:,:,ix,ie)=Sig_lesser(:,:,NS+1,ie)
      Sig_greater(:,:,ix,ie)=Sig_greater(:,:,NS+1,ie)
    enddo
    do ix=length-NS+1,length
      Sig_retarded(:,:,ix,ie)=Sig_retarded(:,:,length-NS,ie)
      Sig_lesser(:,:,ix,ie)=Sig_lesser(:,:,length-NS,ie)
      Sig_greater(:,:,ix,ie)=Sig_greater(:,:,length-NS,ie)
    enddo
  enddo
  ! get leads sigma
  do i=1,NS
    l=(i-1)*NB+1
    h=i*NB
    siglead(l:h,l:h,:,1) = Sig_retarded(:,:,i,:)
    siglead(l:h,l:h,:,2) = Sig_retarded(:,:,length-NS+i,:)    
  enddo
  !
!  call write_spectrum('gw_SigR',iter,Sig_retarded,nen,En,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum('gw_SigL',iter,Sig_lesser,nen,En,length,NB,Lx,(/1.0,1.0/))
!  call write_spectrum('gw_SigG',iter,Sig_greater,nen,En,length,NB,Lx,(/1.0,1.0/))
enddo                
deallocate(siglead)
deallocate(B,cur,tot_cur,tot_ecur)
deallocate(Ispec,Itot,Tr,Te)
deallocate(P_retarded,P_lesser,P_greater)
deallocate(W_retarded,W_lesser,W_greater)
deallocate(wen,nops)
end subroutine green_solve_gw_1D_supermemsaving


subroutine block2full(A,fullA,NB,length)
complex(8),intent(in)::A(NB,NB,length)
complex(8),intent(out)::fullA(NB*length,NB*length)
integer,intent(in)::NB,length
integer::i,l,h
fullA=0.0d0
do i=1,length
  l=(i-1)*NB+1
  h=i*NB
  fullA(l:h,l:h)=A(:,:,i)
enddo
end subroutine block2full


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


! calculate Gr and optionally G<>
! calculate the full GFs but only save the diagonal blocks of GF and current to save memory
! the array structure is different 
subroutine green_calc_g_block(ne,E,num_lead,nm_dev,NB,NX,nm_lead,max_nm_lead,spindeg,Ham,H00,H10,Siglead,T,&
  Scat_Sig_retarded_diag,Scat_Sig_lesser_diag,Scat_Sig_greater_diag,G_retarded_diag,G_lesser_diag,G_greater_diag,&
  jdens,tot_cur,tot_ecur,Itot,Ispec,cur,te,mu,temp)
integer, intent(in) :: num_lead ! number of leads/contacts
integer, intent(in) :: nm_dev   ! size of device Hamiltonian
integer, intent(in) :: nm_lead(num_lead) ! size of lead Hamiltonians
integer, intent(in) :: max_nm_lead ! max size of lead Hamiltonians
integer, intent(in) :: NB, NX
real(8), intent(in) :: E(ne)  ! energy vector
real(8), intent(in) :: spindeg 
real(8), intent(out),optional :: cur(ne,num_lead)  ! current spectrum on leads
real(8), intent(out),optional :: te(ne,num_lead,num_lead)  ! transmission matrix
integer, intent(in) :: ne
complex(8), intent(in) :: Ham(nm_dev,nm_dev)
complex(8), intent(in) :: H00(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian diagonal blocks
complex(8), intent(in) :: H10(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian off-diagonal blocks
complex(8), intent(in) :: Siglead(max_nm_lead,max_nm_lead,ne,num_lead) ! lead sigma_r scattering
complex(8), intent(in) :: T(max_nm_lead,nm_dev,num_lead)  ! coupling matrix between leads and device
complex(8), intent(in), dimension(NB,NB,NX,ne) :: Scat_Sig_retarded_diag,Scat_Sig_lesser_diag,Scat_Sig_greater_diag ! scattering Selfenergy
complex(8), intent(inout), dimension(NB,NB,NX,ne) :: G_retarded_diag,G_lesser_diag,G_greater_diag ! Green's functions
complex(8),intent(out), dimension(NB,NB,NX,ne)::jdens     ! current spec ( only between neighboring slice )
real(8),intent(out), dimension(nm_dev,nm_dev)::tot_cur,tot_ecur ! currents integrated over E
complex(8),intent(out), dimension(NB,NB,NX,ne)::Ispec  ! collision spec ( only diag block )
complex(8),intent(out), dimension(nm_dev,nm_dev)::Itot ! collision integrated over E
real(8), intent(in), optional :: mu(num_lead), temp(num_lead)
integer :: i,j,nm,ie,io,jo
complex(8), allocatable, dimension(:,:) :: S00,G00,GBB,A,sig,sig_lesser,sig_greater,B,C
complex(8), allocatable, dimension(:,:,:) :: gamma_lead
complex(8), allocatable, dimension(:,:) :: Scat_Sig_retarded,Scat_Sig_lesser,Scat_Sig_greater
complex(8), allocatable, dimension(:,:) :: G_retarded,G_lesser,G_greater
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
real(8),parameter::tpi=6.28318530718  
real(8) :: fd, dE
logical :: solve_Gr
dE=E(2)-E(1)
allocate(Scat_Sig_lesser(nm_dev,nm_dev))  ! full matrix of the whole device 
allocate(Scat_Sig_greater(nm_dev,nm_dev))
allocate(Scat_Sig_retarded(nm_dev,nm_dev)) 
allocate(G_lesser(nm_dev,nm_dev))
allocate(G_greater(nm_dev,nm_dev))
allocate(G_retarded(nm_dev,nm_dev)) 
!
allocate(sig(nm_dev,nm_dev))  
jdens=czero
tot_cur=0.0d0
tot_ecur=0.0d0
Itot=czero
solve_Gr = .true.
if (present(cur)) then
  cur=0.0d0
endif
if (present(te)) then
  te=0.0d0
endif
if ((present(cur)).or.(present(te))) then
 allocate(gamma_lead(nm_dev,nm_dev,num_lead))  
endif
do ie = 1, ne
  if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',ne
  ! convert diagonal blocks to full matrix
  call block2full(Scat_Sig_lesser_diag(:,:,:,ie),Scat_Sig_lesser,NB,NX)
  call block2full(Scat_Sig_greater_diag(:,:,:,ie),Scat_Sig_greater,NB,NX)
  call block2full(Scat_Sig_retarded_diag(:,:,:,ie),Scat_Sig_retarded,NB,NX)
  if (.not. solve_Gr) then 
    call block2full(G_retarded_diag(:,:,:,ie),G_retarded,NB,NX)
  else
    G_retarded(:,:) = - Ham(:,:) - Scat_Sig_retarded(:,:) 
  endif
  if (ie .eq. 1) then
    allocate(sig_lesser(nm_dev,nm_dev))
    allocate(sig_greater(nm_dev,nm_dev))          
    allocate(B(nm_dev,nm_dev))
    allocate(C(nm_dev,nm_dev))
  end if
  sig_lesser(:,:) = dcmplx(0.0d0,0.0d0)      
  sig_greater(:,:) = dcmplx(0.0d0,0.0d0)      
  ! compute and add contact self-energies    
  do i = 1,num_lead
    NM = nm_lead(i)    
    allocate(S00(nm,nm))
    allocate(G00(nm,nm))
    allocate(GBB(nm,nm))
    allocate(A(nm_dev,nm))    
    call identity(S00,nm)
    call sancho(NM,E(ie),S00,H00(1:nm,1:nm,i)+siglead(1:nm,1:nm,ie,i),H10(1:nm,1:nm,i),G00,GBB)
    call zgemm('c','n',nm_dev,nm,nm,cone,T(1:nm,1:nm_dev,i),nm,G00,nm,czero,A,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm,cone,A,nm_dev,T(1:nm,1:nm_dev,i),nm,czero,sig,nm_dev)  
    if (solve_Gr) G_retarded(:,:) = G_retarded(:,:) - sig(:,:)
    fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))		
    B(:,:) = conjg(sig(:,:))
    C(:,:) = transpose(B(:,:))
    B(:,:) = sig(:,:) - C(:,:)
    sig_lesser(:,:) = sig_lesser(:,:) - B(:,:)*fd	        
    sig_greater(:,:) = sig_greater(:,:) + B(:,:)*(1.0d0-fd)	 
    if ((present(te)).or.(present(cur))) then
      gamma_lead(:,:,i)= B(:,:) 
    endif       
    deallocate(S00,G00,GBB,A)
  end do  
  if (solve_Gr) then
    do i = 1,nm_dev
      G_retarded(i,i) = G_retarded(i,i) + dcmplx(E(ie),0.0d0)
    end do
    !
    call invert(G_retarded,nm_dev) 
    !
    call full2block(G_retarded,G_retarded_diag(:,:,:,ie),NB,NX)
  endif
  sig_lesser = sig_lesser + Scat_Sig_lesser
  sig_greater = sig_greater + Scat_Sig_greater     
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded,nm_dev,sig_lesser,nm_dev,czero,B,nm_dev) 
  call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded,nm_dev,czero,C,nm_dev)
  G_lesser = C
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded,nm_dev,sig_greater,nm_dev,czero,B,nm_dev) 
  call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded,nm_dev,czero,C,nm_dev)
  G_greater = C      
  call full2block(G_lesser,G_lesser_diag(:,:,:,ie),NB,NX)
  call full2block(G_greater,G_greater_diag(:,:,:,ie),NB,NX)
  !
  ! compute the bond current inside device by using I_ij = H_ij G<_ji - H_ji G^<_ij
  do io=1,nm_dev
    do jo=1,nm_dev
      B(io,jo)=Ham(io,jo)*G_lesser(jo,io) - Ham(jo,io)*G_lesser(io,jo)
    enddo
  enddo    
  B=B*(E(2)-E(1))*e0/tpi/hbar*e0*dble(spindeg)
  call full2block(B, jdens(:,:,:,ie), NB, NX, offdiag=1)
  tot_ecur=tot_ecur+ E(ie)*dble(B)
  tot_cur=tot_cur+ dble(B)  
  ! compute the collision integral
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,Sig_greater(:,:),nm_dev,G_lesser(:,:),nm_dev,czero,B,nm_dev)
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Sig_lesser(:,:),nm_dev,G_greater(:,:),nm_dev,cone,B,nm_dev) 
  Itot(:,:)=Itot(:,:)+B(:,:)
  call full2block(B*spindeg, Ispec(:,:,:,ie), NB, Nx) 
  !
  if ((present(cur)).or.(present(te))) then
    ! calculate current spec and/or transmission at each lead/contact
    do i=1,num_lead                      
      if (present(cur)) then
        fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))		
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,(1.0d0-fd),gamma_lead(:,:,i),nm_dev,G_lesser,nm_dev,czero,B,nm_dev)
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,fd,gamma_lead(:,:,i),nm_dev,G_greater,nm_dev,cone,B,nm_dev)
        do io=1,nm_dev
          cur(ie,i)=cur(ie,i)+ dble(B(io,io))
        enddo
      endif        
      if (present(te)) then
        do j=1,num_lead                      
          if (j.ne.i) then
            call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,gamma_lead(:,:,i),nm_dev,G_retarded,nm_dev,czero,B,nm_dev)
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,gamma_lead(:,:,j),nm_dev,czero,C,nm_dev)
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,C,nm_dev,G_retarded,nm_dev,czero,B,nm_dev)
            do io=1,nm_dev
              te(ie,i,j)=te(ie,i,j)- dble(B(io,io)) ! Gamma = i[Sig^r - Sig^r \dagger] , hence the -1
            enddo
          endif
        enddo ! j
      endif
    enddo ! i
  endif
end do  ! ie
Itot=Itot*dble(E(2)-E(1))/tpi*spindeg
deallocate(sig)
deallocate(Scat_Sig_retarded,Scat_Sig_greater,Scat_Sig_lesser)
deallocate(G_retarded,G_lesser,G_greater)
deallocate(B,C,sig_lesser,sig_greater)
if ((present(cur)).or.(present(te))) then
 deallocate(gamma_lead)
endif
end subroutine green_calc_g_block



! driver for iterating G -> P -> W -> Sig 
! memory saving version of green_solve_gw_1D 
!   the full matrix P and W over energy are not needed in this implementation
!   they are computed per energy point, and the contribution to selfenergy is 
!   immediately added to sigma_x_new matrices, can be extended more easily to 
!   shared memory parallelization and GPU.
subroutine green_solve_gw_1D_memsaving(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,midgap,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,V,&
  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,encut,Egap,ndiagmin,writeGF)
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
integer, intent(in), optional :: ndiagmin
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg, Egap,midgap(2)
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: V(nm_dev,nm_dev)
logical,intent(in)::ldiag
real(8),intent(in)::encut(2) ! intraband and interband cutoff for P
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
logical, intent(in), optional :: writeGF 
!----
complex(8),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:)
real(8),allocatable::wen(:) ! energy vector for P and W
integer,allocatable::nops(:) ! discretized energy for P and W
real(8),allocatable::Tr(:,:) ! current spectrum on leads
real(8),allocatable::Te(:,:,:) ! transmission matrix spectrum 
integer :: iter,ie,iop,nnop,nnop1,nnop2
integer :: i,j,nm,l,h,ndiag,nop
logical :: lwriteGF
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
complex(8),allocatable,dimension(:,:) ::  P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
if (present(writeGF)) then
  lwriteGF=writeGF
else
  lwriteGF=.false.
endif
print *,'================= green_solve_gw_1D_memsaving ================='
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
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
print *, ' Encut: intra    inter    Eg (eV)' 
print '(A6,3F8.3)',' ',encut,egap
print *, ' Nop='
print '(10I5)',nops
print *, ' Eop= (eV)'
print '(6F8.3)',wen
print *, '--------------------------------------------------------------'
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
allocate(tr(nen,2))
allocate(te(nen,2,2))
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
  call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham,H00lead,H10lead,Siglead,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,cur=Tr,te=Te,mu=mu,temp=(/temps,tempd/))
!  if (iter == 0) then     
!    call calc_n_electron(G_lesser,G_greater,nen,En,NS,NB,nm_dev,nelec,pelec,midgap)  ! calculate N and P at contacts  
!    print '(a8,f15.4,a8,f15.4)', 'Ns=',nelec(1),'Nd=',nelec(2)
!    print '(a8,f15.4,a8,f15.4)', 'Ps=',pelec(1),'Pd=',pelec(2)
!  else    
!    call calc_fermi_level(G_retarded,nelec,pelec,nen,En,NS,NB,nm_dev,(/temps,tempd/),mu,midgap)    
!    mu=(/ mus, mud /) + 0.5*(mu-(/ mus, mud /)) ! move Fermi level at contacts
!    print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)    
!  end if  
  ! 
  call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
  call write_current_spectrum('gw_Jdens',iter,cur,nen,en,length,NB,Lx)
  call write_current('gw_I',iter,tot_cur,length,NB,NS,Lx)
  call write_current('gw_EI',iter,tot_ecur,length,NB,NS,Lx)
  call write_spectrum('gw_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
  call write_spectrum('gw_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('gw_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))
  call write_transmission_spectrum('gw_trL',iter,Tr(:,1)*spindeg,nen,En)
  call write_transmission_spectrum('gw_trR',iter,Tr(:,2)*spindeg,nen,En)
  call write_transmission_spectrum('gw_TE_LR',iter,Te(:,1,2)*spindeg,nen,En)
  call write_transmission_spectrum('gw_TE_RL',iter,Te(:,2,1)*spindeg,nen,En)
  G_retarded(:,:,:)=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
  G_lesser(:,:,:)=dcmplx(0.0d0*dble(G_lesser),aimag(G_lesser))
  G_greater(:,:,:)=dcmplx(0.0d0*dble(G_greater),aimag(G_greater))
  !call write_matrix_summed_overE('Gr',iter,G_retarded,nen,en,length,NB,(/1.0,1.0/))
  if (lwriteGF) then
    call write_matrix_E('G_r',0,G_retarded,nen,en,length,NB,(/1.0d0,1.0d0/))
    !call write_matrix_E('G_l',0,G_lesser,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('G_g',0,G_greater,nen,en,length,NB,(/1.0,1.0/))
  endif
  !        
  ! empty sigma_x_new matrices for accumulation
  sig_retarded_new=czero
  sig_lesser_new=czero
  sig_greater_new=czero
  print *, 'calc P, solve W, add to Sigma_new'     
  ndiag=NB*(min(NS,iter))
  if (ldiag) ndiag=0  
  if (present(ndiagmin)) ndiag=max(ndiagmin,ndiag)
  if (lwriteGF) ndiag=nm_dev
  print *,'ndiag=',min(ndiag,nm_dev)
  !
  !print *,'   i / n :  Nop   Eop (eV)'
  do iop=1,nnop        
    !print '(I5,A,I5,A,I5,F8.3)',iop,'/',nnop,':',nops(iop),wen(iop)    
    nop=nops(iop)
    P_lesser = czero
    P_greater = czero
    P_retarded = czero
    !$omp parallel default(none) private(l,h,i,ie) shared(ndiag,nop,nen,P_lesser,P_greater,P_retarded,nm_dev,G_lesser,G_greater,G_retarded)  
    !$omp do
    do i = 1, nm_dev        
      do ie = max(nop+1,1),min(nen,nen+nop)                   
        l=max(i-ndiag,1)
        h=min(nm_dev,i+ndiag)                            
        P_lesser(i,l:h) = P_lesser(i,l:h) + G_lesser(i,l:h,ie) * G_greater(l:h,i,ie-nop)
        P_greater(i,l:h) = P_greater(i,l:h) + G_greater(i,l:h,ie) * G_lesser(l:h,i,ie-nop) 
        P_retarded(i,l:h) = P_retarded(i,l:h) + G_lesser(i,l:h,ie) * conjg(G_retarded(i,l:h,ie-nop)) +&
                  & G_retarded(i,l:h,ie) * G_lesser(l:h,i,ie-nop)        
      enddo
    enddo    
    !$omp end do
    !$omp end parallel
    dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )* spindeg    
    P_lesser=P_lesser*dE
    P_greater=P_greater*dE  
    P_retarded=P_retarded*dE
    if (lwriteGF) then
      call write_matrix('P_r',0,P_retarded(:,:),wen(iop),length,NB,(/1.0d0,1.0d0/))
    endif
    !
    ! calculate W
    call green_calc_w(0,NB,NS,nm_dev,P_retarded,P_lesser,P_greater,V,W_retarded,W_lesser,W_greater)
    !
    if (lwriteGF) then
      call write_matrix('W_r',0,W_retarded(:,:),wen(iop),length,NB,(/1.0d0,1.0d0/))
    endif
    !
    ! Accumulate the GW to Sigma
    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)  
    !$omp parallel default(none) private(l,h,i,ie) shared(ndiag,nop,nen,Sig_lesser_new,Sig_greater_new,Sig_retarded_new,W_lesser,W_greater,W_retarded,nm_dev,G_lesser,G_greater,G_retarded)  
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
  enddo                                
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
    Sig_lesser_new(:,:,ie) = (Sig_lesser_new(:,:,ie) + B(:,:))/2.0d0
    B(:,:)=transpose(Sig_greater_new(:,:,ie))
    Sig_greater_new(:,:,ie) = (Sig_greater_new(:,:,ie) + B(:,:))/2.0d0
  enddo
  !!!Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
  !!!Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
  !
  if (lwriteGF) then
    call write_matrix_E('Sigma_r',0,Sig_retarded_new,nen,en,length,NB,(/1.0d0,1.0d0/))
    !call write_matrix_E('Sigma_l',0,Sig_lesser_new,nen,en,length,NB,(/1.0,1.0/))
    !call write_matrix_E('Sigma_g',0,Sig_greater_new,nen,en,length,NB,(/1.0,1.0/))
  endif
  ! mixing with the previous one
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
  Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)    
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    call expand_size_bycopy(Sig_retarded(:,:,ie),nm_dev,NB,2)
    call expand_size_bycopy(Sig_lesser(:,:,ie),nm_dev,NB,2)
    call expand_size_bycopy(Sig_greater(:,:,ie),nm_dev,NB,2)
  enddo
  ! get leads sigma
  siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
  siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)    
  !
  call write_spectrum('gw_SigR',iter,Sig_retarded,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('gw_SigL',iter,Sig_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('gw_SigG',iter,Sig_greater,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  !call write_matrix_summed_overE('Sigma_r',iter,Sig_retarded,nen,en,length,NB,(/1.0,1.0/))
  !!!! calculate collision integral
  call calc_collision(Sig_lesser_new,Sig_greater_new,G_lesser,G_greater,nen,en,spindeg,nm_dev,Itot,Ispec)
  call write_spectrum('gw_Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
enddo                
deallocate(siglead)
deallocate(B,cur,tot_cur,tot_ecur)
deallocate(Ispec,Itot,Tr,Te)
deallocate(P_retarded,P_lesser,P_greater)
deallocate(W_retarded,W_lesser,W_greater)
deallocate(wen,nops)
end subroutine green_solve_gw_1D_memsaving

! driver for iterating G -> P -> W -> Sig 
subroutine green_solve_gw_1D(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
  alpha_mix,nen,En,nb,ns,Ham,H00lead,H10lead,T,V,invV,&
  G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
  W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
  Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag)
integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg
complex(8),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
complex(8), intent(in):: V(nm_dev,nm_dev), invV(nm_dev,nm_dev)
logical,intent(in)::ldiag
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,Sig_retarded_new,Sig_lesser_new,Sig_greater_new
complex(8),intent(inout),dimension(nm_dev,nm_dev,nen*2+1) ::  P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater
!----
complex(8),allocatable::siglead(:,:,:,:) ! lead scattering sigma_retarded
complex(8),allocatable,dimension(:,:):: B ! tmp matrix
real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:)
integer :: iter,ie,nopmax
integer :: i,j,nm,nop,l,h,iop,ndiag
complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
complex(8) :: dE
real(8)::nelec(2),mu(2),pelec(2)
print *,'====== green_solve_gw_1D ======'
allocate(siglead(NB*NS,NB*NS,nen,2))
allocate(wen(nen*2-1))
wen(nen:1:-1)=-(en-minval(en))
wen(nen:nen*2-1)=(en-minval(en))
! get leads sigma
siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)  
allocate(B(nm_dev,nm_dev))
allocate(tot_cur(nm_dev,nm_dev))
allocate(tot_ecur(nm_dev,nm_dev))
allocate(cur(nm_dev,nm_dev,nen))
allocate(Ispec(nm_dev,nm_dev,nen))
allocate(Itot(nm_dev,nm_dev))
mu=(/ mus, mud /)
print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
do iter=0,niter
    print *,'+ iter=',iter
  ! empty files for sancho 
!  open(unit=101,file='sancho_gbb.dat',status='unknown')
!  close(101)
!  open(unit=101,file='sancho_g00.dat',status='unknown')
!  close(101)
!  open(unit=101,file='sancho_sig.dat',status='unknown')
!  close(101)
  print *, 'calc G'  
  call green_calc_g(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham,H00lead,H10lead,Siglead,T,Sig_retarded,Sig_lesser,Sig_greater,G_retarded,G_lesser,G_greater,mu=mu,temp=(/temps,tempd/))
 ! if (iter == 0) then     
 !   call calc_n_electron(G_lesser,G_greater,nen,En,NS,NB,nm_dev,nelec,pelec)    
 ! else    
 !   call calc_fermi_level(G_retarded,nelec,pelec,nen,En,NS,NB,nm_dev,(/temps,tempd/),mu)    
 !   mu=(/ mus, mud /) - 0.2*sum(mu-(/ mus, mud /))/2.0d0 ! move Fermi level because Sig_GW shifts slightly the energies
 !   print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)    
 ! end if  
  call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
  call write_current_spectrum('gw_Jdens',iter,cur,nen,en,length,NB,Lx)
  call write_current('gw_I',iter,tot_cur,length,NB,NS,Lx)
  call write_current('gw_EI',iter,tot_ecur,length,NB,NS,Lx)
  call write_spectrum('gw_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
  call write_spectrum('gw_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('gw_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))
  !call write_matrix_summed_overE('Gr',iter,G_retarded,nen,en,length,NB,(/1.0,1.0/))
  !call write_matrix_E('G_r',iter,G_retarded,nen,en,length,NB,(/1.0,1.0/))
  !call write_matrix_E('G_l',iter,G_lesser,nen,en,length,NB,(/1.0,1.0/))
  !call write_matrix_E('G_g',iter,G_greater,nen,en,length,NB,(/1.0,1.0/))
  !        
  print *, 'calc P'  
  ! Pij^<>(hw) = \int_dE Gij^<>(E) * Gji^><(E-hw)
  ! Pij^r(hw)  = \int_dE Gij^<(E) * Gji^a(E-hw) + Gij^r(E) * Gji^<(E-hw)
  P_lesser=czero
  P_greater=czero
  P_retarded=czero
  ndiag=NB*(min(NS,iter))
  if (ldiag) ndiag=0  
  nopmax=nen/2-1
  dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi )	* spindeg
  print *,'ndiag=',min(ndiag,nm_dev)
  !$omp parallel default(none) private(l,h,i,j) shared(ndiag,nopmax,P_lesser,P_greater,P_retarded,nen,En,nm_dev,G_lesser,G_greater,G_retarded)  
  !$omp do
  do i=1,nm_dev
    l=max(i-ndiag,1)
    h=min(nm_dev,i+ndiag)
    do j=l,h
      P_lesser(i,j,:) =corr1d(nen,G_lesser(i,j,:),G_greater(j,i,:),method='simple')
      P_greater(i,j,:)=corr1d(nen,G_greater(i,j,:),G_lesser(j,i,:),method='simple')
      P_retarded(i,j,:)=corr1d(nen,G_lesser(i,j,:),conjg(G_retarded(i,j,:)),method='simple')+&
                        corr1d(nen,G_retarded(i,j,:),G_lesser(j,i,:),method='simple')
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  P_lesser=P_lesser*dE
  P_greater=P_greater*dE
  P_retarded=P_retarded*dE
  call write_spectrum('PR',iter,P_retarded,nen*2-1,wen,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('PL',iter,P_lesser,nen*2-1,wen,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('PG',iter,P_greater,nen*2-1,wen,length,NB,Lx,(/1.0d0,1.0d0/))
  !call write_matrix_summed_overE('P_r',iter,P_retarded(:,:,nen/2+1:nen/2+nen),nen,en-en(nen/2),length,NB,(/1.0,1.0/))
!  call write_matrix_E('P_r',iter,P_retarded(:,:,nen/2+1:nen/2+nen),nen,en-en(nen/2),length,NB,(/1.0,1.0/))
  !call write_matrix_E('P_l',iter,P_lesser(:,:,nen/2+1:nen/2+nen),nen,en-en(nen/2),length,NB,(/1.0,1.0/))
  !call write_matrix_E('P_g',iter,P_greater(:,:,nen/2+1:nen/2+nen),nen,en-en(nen/2),length,NB,(/1.0,1.0/))
  !
  print *, 'calc W'  
  W_lesser=czero
  W_greater=czero
  W_retarded=czero
  do nop=1,nen*2-1
    print '(I5,A,I5)',nop,'/',nen*2-1    
    !! B = -V P
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,P_retarded(:,:,nop),nm_dev,czero,B,nm_dev)
    !! B = I-VP    
    do i=1,nm_dev
      B(i,i) = 1.0d0 + B(i,i)
    enddo        
    !!!! calculate W^r = (I - V P^r)^-1 V       
    call invert(B,nm_dev)
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V,nm_dev,czero,W_retarded(:,:,nop),nm_dev)           
    ! calculate W^< and W^> = W^r P^<> W^r dagger
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,W_retarded(:,:,nop),nm_dev,P_lesser(:,:,nop),nm_dev,czero,B,nm_dev) 
    call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,W_retarded(:,:,nop),nm_dev,czero,W_lesser(:,:,nop),nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,W_retarded(:,:,nop),nm_dev,P_greater(:,:,nop),nm_dev,czero,B,nm_dev) 
    call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,W_retarded(:,:,nop),nm_dev,czero,W_greater(:,:,nop),nm_dev)   
  enddo
  call write_spectrum('WR',iter,W_retarded,nen*2-1,wen,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('WL',iter,W_lesser,  nen*2-1,wen,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('WG',iter,W_greater, nen*2-1,wen,length,NB,Lx,(/1.0d0,1.0d0/))
  !call write_matrix_summed_overE('W_r',iter,W_retarded(:,:,nen/2+1:nen/2+nen),nen,en,length,NB,(/1.0,1.0/))
!  call write_matrix_E('W_r',iter,W_retarded(:,:,nen/2+1:nen/2+nen),nen,en-en(nen/2),length,NB,(/1.0,1.0/))
  !call write_matrix_E('W_g',iter,W_greater(:,:,nen/2+1:nen/2+nen),nen,en-en(nen/2),length,NB,(/1.0,1.0/))
  !call write_matrix_E('W_l',iter,W_lesser(:,:,nen/2+1:nen/2+nen),nen,en-en(nen/2),length,NB,(/1.0,1.0/))
  !
  print *, 'calc SigGW'
  Sig_greater_new = dcmplx(0.0d0,0.0d0)
  Sig_lesser_new = dcmplx(0.0d0,0.0d0)
  Sig_retarded_new = dcmplx(0.0d0,0.0d0)
  dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi)      
  !ndiag=NB*NS*2  
  ndiag=NB*(min(NS,iter))
  if (ldiag) ndiag=0  
  print *,'ndiag=',min(ndiag,nm_dev)
  ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)
  !$omp parallel default(none) private(l,h,i,j) shared(ndiag,nen,Sig_lesser_new,Sig_greater_new,Sig_retarded_new,W_lesser,W_greater,W_retarded,nm_dev,G_lesser,G_greater,G_retarded)  
  !$omp do
  do i=1,nm_dev
    l=max(i-ndiag,1)
    h=min(nm_dev,i+ndiag)
    do j=l,h
      Sig_lesser_new(i,j,:)  =conv1d(nen,G_lesser(i,j,:),W_lesser(i,j,:),method='simple')
      Sig_greater_new(i,j,:) =conv1d(nen,G_greater(i,j,:),W_greater(i,j,:),method='simple')
      Sig_retarded_new(i,j,:)=conv1d(nen,G_lesser(i,j,:),W_retarded(i,j,:),method='simple') +&
                              conv1d(nen,G_retarded(i,j,:),W_lesser(i,j,:),method='simple') +&
                              conv1d(nen,G_retarded(i,j,:),W_retarded(i,j,:),method='simple')
    enddo
  enddo                            
  !$omp end do
  !$omp end parallel
  Sig_lesser_new = Sig_lesser_new  * dE
  Sig_greater_new= Sig_greater_new * dE
  Sig_retarded_new=Sig_retarded_new* dE
  Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
  ! symmetrize the selfenergies
  do ie=1,nen
    B(:,:)=transpose(Sig_retarded_new(:,:,ie))
    Sig_retarded_new(:,:,ie) = (Sig_retarded_new(:,:,ie) + B(:,:))/2.0d0    
    B(:,:)=transpose(Sig_lesser_new(:,:,ie))
    Sig_lesser_new(:,:,ie) = (Sig_lesser_new(:,:,ie) + B(:,:))/2.0d0
    B(:,:)=transpose(Sig_greater_new(:,:,ie))
    Sig_greater_new(:,:,ie) = (Sig_greater_new(:,:,ie) + B(:,:))/2.0d0
  enddo
  !!!Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
  !!!Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
!  call write_matrix_E('Sigma_r',iter,Sig_retarded_new,nen,en,length,NB,(/1.0,1.0/))
  !call write_matrix_E('Sigma_l',iter,Sig_lesser_new,nen,en,length,NB,(/1.0,1.0/))
  !call write_matrix_E('Sigma_g',iter,Sig_greater_new,nen,en,length,NB,(/1.0,1.0/))
  ! mixing with the previous one
  Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
  Sig_lesser  = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
  Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)    
  ! make sure self-energy is continuous near leads (by copying edge block)
  do ie=1,nen
    call expand_size_bycopy(Sig_retarded(:,:,ie),nm_dev,NB,3)
    call expand_size_bycopy(Sig_lesser(:,:,ie),nm_dev,NB,3)
    call expand_size_bycopy(Sig_greater(:,:,ie),nm_dev,NB,3)
  enddo
  ! get leads sigma
  siglead(:,:,:,1) = Sig_retarded(1:NB*NS,1:NB*NS,:)
  siglead(:,:,:,2) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)    
!  do i=1,NB*NS
!    do j=1,NB*NS
!       if(i.ne.j) then
!         siglead(i,j,:,:)=dcmplx(0.0*dble(siglead(i,j,:,:)),0.d0*aimag(siglead(i,j,:,:)))
!       endif
!    enddo
!  enddo
  !
  call write_spectrum('gw_SigR',iter,Sig_retarded,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('gw_SigL',iter,Sig_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  call write_spectrum('gw_SigG',iter,Sig_greater,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
  !call write_matrix_summed_overE('Sigma_r',iter,Sig_retarded,nen,en,length,NB,(/1.0,1.0/))
  !!!! calculate collision integral
  !call calc_collision(Sig_lesser_new,Sig_greater_new,G_lesser,G_greater,nen,en,spindeg,nm_dev,Itot,Ispec)
  !call write_spectrum('gw_Scat',iter,Ispec,nen,En,length,NB,Lx,(/1.0,1.0/))
enddo                
deallocate(siglead)
deallocate(B,cur,tot_cur,tot_ecur)
deallocate(Ispec,Itot)
deallocate(wen)
end subroutine green_solve_gw_1D


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

! calculate number of electrons and holes from G< and G> 
subroutine calc_n_electron(G_lesser,G_greater,nen,E,NS,NB,nm_dev,nelec,pelec,midgap)
complex(8), intent(in) :: G_lesser(nm_dev,nm_dev,nen)
complex(8), intent(in) :: G_greater(nm_dev,nm_dev,nen)
real(8), intent(in)    :: E(nen),midgap(2)
integer, intent(in)    :: NS,NB,nm_dev,nen
real(8), intent(out)   :: nelec(2),pelec(2)
real(8)::dE
integer::i,j
nelec=0.0d0
pelec=0.0d0
dE=E(2)-E(1)
do i=1,nen
  do j=1,NS*NB
    if (E(i)>midgap(1))then
      nelec(1)=nelec(1)+aimag(G_lesser(j,j,i))*dE
    else
      pelec(1)=pelec(1)-aimag(G_greater(j,j,i))*dE
    endif
  enddo
enddo
do i=1,nen
  do j=nm_dev-NS*NB+1,nm_dev
    if (E(i)>midgap(2))then
      nelec(2)=nelec(2)+aimag(G_lesser(j,j,i))*dE
    else
      pelec(2)=pelec(2)-aimag(G_greater(j,j,i))*dE
    endif
  enddo
enddo
end subroutine calc_n_electron

! determine the quasi-fermi level from the Gr and electron/hole number
subroutine calc_fermi_level(G_retarded,nelec,pelec,nen,En,NS,NB,nm_dev,Temp,mu,midgap)
real(8),intent(in)::Temp(2)
real(8),intent(out)::mu(2)
complex(8), intent(in) :: G_retarded(nm_dev,nm_dev,nen)
real(8), intent(in)    :: En(nen)
integer, intent(in)    :: NS,NB,nm_dev,nen
real(8), intent(in)    :: nelec(2),pelec(2),midgap(2)
real(8)::dE,n(nen),p(nen),fd,mun,mup
real(8),allocatable::dos(:),K(:),Q(:),fermi_derivative(:,:)
REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
integer::nefermi,i,j,l
allocate(dos(nen))
allocate(K(nen))
allocate(Q(nen))
!nefermi=1+2*floor(10.0*boltz*maxval(temp)/(En(2)-En(1))) ! number of energy points for the fermi derivative function
nefermi=nen
dE=En(2)-En(1)
do i=1,2
  dos=0.0d0
  K=0.0d0
  Q=0.0d0
  n=0.0d0
  p=0.0d0
  do j=1,nen
    if (i==1) then
      do l=1,NS*NB
        dos(j) = dos(j)+aimag(G_retarded(l,l,j))
      enddo
    else
      do l=nm_dev-NS*NB+1,nm_dev
        dos(j) = dos(j)+aimag(G_retarded(l,l,j))
      enddo
    endif    
    dos(j)=-2.0d0*dos(j)
    if ((j>1).and.(En(j)>midgap(i))) K(j)=K(j-1)+dos(j)*dE    
    if ((j>1).and.(En(nen-j+1)<midgap(i))) Q(nen-j+1)=Q(nen-j+2)+dos(j)*dE    
  enddo
  ! search for the Fermi level
!  if (dE<(BOLTZ*TEMP(i))) then
!    allocate(fermi_derivative(nefermi,2))
!    do j=1,nefermi
!      fd=ferm(dble(dE)*dble(j-nefermi/2-1)/(BOLTZ*TEMP(i)))
!      fermi_derivative(j,i) = -fd*(1.0d0-fd)/(BOLTZ*TEMP(i))    
!    enddo
!    do j=1,nen
!      n(j)=-sum(K(max(j-nefermi/2+1,1):min(j+nefermi/2,nen))*fermi_derivative(max(nefermi/2-j,1):min(nen-j,nefermi),i))*dE
!    enddo
!    deallocate(fermi_derivative)
!  else ! energy grid too coarse
!    ! approximate F-D to step-function
    n = K
    p = Q
!  endif
  n=n-nelec(i)  
  p=p-pelec(i)
  do j=2,nen
    if ((n(j)>=0.0).and.(n(j-1)<=0.0)) then
      mun=En(j)
      exit
    endif
  enddo
  do j=nen-1,1,-1
    if ((p(j)>=0.0).and.(p(j+1)<=0.0)) then
      mup=En(j)
      exit
    endif
  enddo
  if (nelec(i)>pelec(i)) then
    mu(i)=mun
  else
    mu(i)=mup
  endif
enddo
deallocate(dos,K)
end subroutine calc_fermi_level


! calculate Gr and optionally G<>
subroutine green_calc_g(ne,E,num_lead,nm_dev,nm_lead,max_nm_lead,Ham,H00,H10,Siglead,T,Scat_Sig_retarded,Scat_Sig_lesser,Scat_Sig_greater,G_retarded,G_lesser,G_greater,cur,te,mu,temp,mode)
integer, intent(in) :: num_lead ! number of leads/contacts
integer, intent(in) :: nm_dev   ! size of device Hamiltonian
integer, intent(in) :: nm_lead(num_lead) ! size of lead Hamiltonians
integer, intent(in) :: max_nm_lead ! max size of lead Hamiltonians
real(8), intent(in) :: E(ne)  ! energy vector
real(8), intent(out),optional :: cur(ne,num_lead)  ! current spectrum on leads
real(8), intent(out),optional :: te(ne,num_lead,num_lead)  ! transmission matrix
integer, intent(in) :: ne
complex(8), intent(in) :: Ham(nm_dev,nm_dev)
complex(8), intent(in) :: H00(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian diagonal blocks
complex(8), intent(in) :: H10(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian off-diagonal blocks
complex(8), intent(in) :: Siglead(max_nm_lead,max_nm_lead,ne,num_lead) ! lead sigma_r scattering
complex(8), intent(in) :: T(max_nm_lead,nm_dev,num_lead)  ! coupling matrix between leads and device
complex(8), intent(in) :: Scat_Sig_retarded(nm_dev,nm_dev,ne) ! scattering Selfenergy
complex(8), intent(in) :: Scat_Sig_lesser(nm_dev,nm_dev,ne)
complex(8), intent(in) :: Scat_Sig_greater(nm_dev,nm_dev,ne)
complex(8), intent(inout) :: G_retarded(nm_dev,nm_dev,ne)
complex(8), intent(inout), optional :: G_lesser(nm_dev,nm_dev,ne)
complex(8), intent(inout), optional :: G_greater(nm_dev,nm_dev,ne)
real(8), intent(in), optional :: mu(num_lead), temp(num_lead)
character(len=*),intent(in), optional :: mode
integer :: i,j,nm,ie,io
complex(8), allocatable, dimension(:,:) :: S00,G00,GBB,A,sig,sig_lesser,sig_greater,B,C
complex(8), allocatable, dimension(:,:,:) :: gamma_lead
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
REAL(8), PARAMETER  :: BOLTZ=8.61734d-05 !eV K-1
real(8) :: fd
logical :: solve_Gr
allocate(sig(nm_dev,nm_dev))  
solve_Gr = .true.
if (present(mode).and.(mode=='use_gr')) then
  solve_Gr = .false.
endif
if (present(cur)) then
  cur=0.0d0
endif
if (present(te)) then
  te=0.0d0
endif
if ((present(cur)).or.(present(te))) then
 allocate(gamma_lead(nm_dev,nm_dev,num_lead))  
endif
do ie = 1, ne
  if (mod(ie,100)==0) print '(I5,A,I5)',ie,'/',ne
  if (solve_Gr) G_retarded(:,:,ie) = - Ham(:,:) - Scat_Sig_retarded(:,:,ie) 
  if ((present(G_lesser)).or.(present(G_greater))) then    
    if (ie .eq. 1) then
      allocate(sig_lesser(nm_dev,nm_dev))
      allocate(sig_greater(nm_dev,nm_dev))          
      allocate(B(nm_dev,nm_dev))
      allocate(C(nm_dev,nm_dev))
    end if
    sig_lesser(:,:) = dcmplx(0.0d0,0.0d0)      
    sig_greater(:,:) = dcmplx(0.0d0,0.0d0)      
  end if    
  ! compute and add contact self-energies    
!  open(unit=101,file='sancho_gbb.dat',status='unknown',position='append')
!  open(unit=102,file='sancho_g00.dat',status='unknown',position='append')
!  open(unit=103,file='sancho_sig.dat',status='unknown',position='append')
  do i = 1,num_lead
    NM = nm_lead(i)    
    allocate(S00(nm,nm))
    allocate(G00(nm,nm))
    allocate(GBB(nm,nm))
    allocate(A(nm_dev,nm))    
    call identity(S00,nm)
    call sancho(NM,E(ie),S00,H00(1:nm,1:nm,i)+siglead(1:nm,1:nm,ie,i),H10(1:nm,1:nm,i),G00,GBB)
    call zgemm('c','n',nm_dev,nm,nm,cone,T(1:nm,1:nm_dev,i),nm,G00,nm,czero,A,nm_dev) 
    call zgemm('n','n',nm_dev,nm_dev,nm,cone,A,nm_dev,T(1:nm,1:nm_dev,i),nm,czero,sig,nm_dev)  
!    write(101,'(i4,2E15.4)') i, E(ie), -aimag(trace(GBB,nm))*2.0d0
!    write(102,'(i4,2E15.4)') i, E(ie), -aimag(trace(G00,nm))*2.0d0
!    write(103,'(i4,2E15.4)') i, E(ie), -aimag(trace(sig,nm_dev))*2.0d0
    if (solve_Gr) G_retarded(:,:,ie) = G_retarded(:,:,ie) - sig(:,:)
    if ((present(G_lesser)).or.(present(G_greater))) then      
      fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))		
      B(:,:) = conjg(sig(:,:))
      C(:,:) = transpose(B(:,:))
      B(:,:) = sig(:,:) - C(:,:)
      sig_lesser(:,:) = sig_lesser(:,:) - B(:,:)*fd	        
      sig_greater(:,:) = sig_greater(:,:) + B(:,:)*(1.0d0-fd)	 
      if ((present(te)).or.(present(cur))) then
        gamma_lead(:,:,i)= B(:,:) 
      endif       
    end if
    deallocate(S00,G00,GBB,A)
  end do  
!  close(101)
!  close(102)
!  close(103)
  if (solve_Gr) then
    do i = 1,nm_dev
      G_retarded(i,i,ie) = G_retarded(i,i,ie) + dcmplx(E(ie),0.0d0)
    end do
    !
    call invert(G_retarded(:,:,ie),nm_dev) 
    !
  endif
  if ((present(G_lesser)).or.(present(G_greater))) then    
    sig_lesser = sig_lesser + Scat_Sig_lesser(:,:,ie)
    sig_greater = sig_greater + Scat_Sig_greater(:,:,ie)     
    if (present(G_lesser)) then
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_lesser,nm_dev,czero,B,nm_dev) 
      call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
      G_lesser(:,:,ie) = C
    end if
    if (present(G_greater)) then
      call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_greater,nm_dev,czero,B,nm_dev) 
      call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
      G_greater(:,:,ie) = C      
    end if      
    if ((present(cur)).or.(present(te))) then
      ! calculate current spec and/or transmission at each lead/contact
      do i=1,num_lead                      
        if (present(cur)) then
          fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))		
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,(1.0d0-fd),gamma_lead(:,:,i),nm_dev,G_lesser(:,:,ie),nm_dev,czero,B,nm_dev)
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,fd,gamma_lead(:,:,i),nm_dev,G_greater(:,:,ie),nm_dev,cone,B,nm_dev)
          do io=1,nm_dev
            cur(ie,i)=cur(ie,i)+ dble(B(io,io))
          enddo
        endif        
        if (present(te)) then
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
        endif
      enddo
    endif
  end if 
end do  
deallocate(sig)
if ((present(G_lesser)).or.(present(G_greater))) then      
  deallocate(B,C,sig_lesser,sig_greater)
end if
if ((present(cur)).or.(present(te))) then
 deallocate(gamma_lead)
endif
end subroutine green_calc_g

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
      do jo=1,nm_dev
        B(io,jo)=H(io,jo)*G_lesser(jo,io,ie) - H(jo,io)*G_lesser(io,jo,ie)
      enddo
    enddo    
    B=B*(En(2)-En(1))*e0/tpi/hbar*e0*dble(spindeg)
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
      do jb=1,nb       
        do j=ii,min(ii+NS-1,length-1)
          tr = tr+ cur((ii-1)*nb+ib,j*nb+jb)
        enddo
      enddo                        
    end do
    write(11,'(2E18.4)') dble(ii)*Lx, tr
  end do
end subroutine write_current

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

! write current spectrum into file (pm3d map)
subroutine write_current_spectrum_block(dataset,i,cur,nen,en,length,NB,Lx)
  character(len=*), intent(in) :: dataset
  complex(8), intent(in) :: cur(:,:,:,:)
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
              tr = tr+ cur(ib,jb,j,ie)
            enddo                        
          enddo
          write(11,'(3E18.4)') dble(j)*Lx, en(ie), dble(tr)
      enddo
      write(11,*)    
  enddo
  close(11)
end subroutine write_current_spectrum_block

! write trace of diagonal blocks
subroutine write_trace(dataset,i,G,length,NB,Lx,coeff,E)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:)
integer, intent(in)::i,length,NB
real(8), intent(in)::Lx,coeff(2)
real(8), intent(in),optional::E
integer:: ie,j,ib
complex(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown', position="append", action="write")
do j = 1,length
    tr=0.0d0          
    do ib=1,nb
        tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib)            
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

! write spectrum into file (pm3d map)
subroutine write_spectrum_block(dataset,i,G,nen,en,length,NB,Lx,coeff)
  character(len=*), intent(in) :: dataset
  complex(8), intent(in) :: G(:,:,:,:)
  integer, intent(in)::i,nen,length,NB
  real(8), intent(in)::Lx,en(nen),coeff(2)
  integer:: ie,j,ib
  complex(8)::tr
  open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
  do ie = 1,nen
      do j = 1,length
          tr=0.0d0          
          do ib=1,nb
              tr = tr+ G(ib,ib,j,ie)            
          end do
          write(11,'(4E18.4)') (j-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
      end do
      write(11,*)    
  end do
  close(11)
end subroutine write_spectrum_block

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

! write a matrix summed over energy index into a file
subroutine write_matrix_summed_overE(dataset,i,G,nen,en,length,NB,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:,:)
integer, intent(in)::i,nen,length,NB
real(8), intent(in)::en(nen),coeff(2)
integer:: ie,j,ib,l
complex(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do l=1,length*NB
  do j = 1,length*NB
    tr=0.0d0          
    do ie=1,nen  
        tr = tr+ G(l,j,ie)            
    end do
    tr=tr/dble(nen)
    write(11,'(2I8,2E18.4)') l,j, dble(tr)*coeff(1), aimag(tr)*coeff(2)        
  end do
  write(11,*)    
end do
close(11)
end subroutine write_matrix_summed_overE

! write a matrix for one energy index into a file
subroutine write_matrix(dataset,i,G,en,length,NB,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:)
integer, intent(in)::i,length,NB
real(8), intent(in)::en,coeff(2)
integer:: ie,j,ib,l
complex(8)::tr
logical :: lexist
inquire(file=trim(dataset)//TRIM(STRING(i))//'.dat', exist=lexist)
if (lexist) then
    open(11, file=trim(dataset)//TRIM(STRING(i))//'.dat', status="old", position="append", action="write")
else
    open(11, file=trim(dataset)//TRIM(STRING(i))//'.dat', status="new", action="write")
end if
do l=1,length*NB
    do j = 1,length*NB
        tr = G(l,j)            
        write(11,'(E18.6,2I8,2E18.6)') en,l,j, dble(tr)*coeff(1), aimag(tr)*coeff(2)        
    end do
end do
write(11,*)    
close(11)
end subroutine write_matrix

! write a matrix for all energy index into a file
subroutine write_matrix_E(dataset,i,G,nen,en,length,NB,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:,:)
integer, intent(in)::i,nen,length,NB
real(8), intent(in)::en(nen),coeff(2)
integer:: ie,j,ib,l
complex(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie=1,nen  
    do l=1,length*NB
        do j = 1,length*NB
            tr = G(l,j,ie)            
            write(11,'(E18.6,2I8,2E18.6)') en(ie),l,j, dble(tr)*coeff(1), aimag(tr)*coeff(2)        
        end do
    end do
    write(11,*)    
end do
close(11)
end subroutine write_matrix_E


! write current spectrum into file (pm3d map)
subroutine write_current_spectrum_summed_over_kz(dataset,i,cur,nen,en,nphiz,length,NB,Lx)
character(len=*), intent(in) :: dataset
real(8), intent(in) :: cur(:,:,:,:)
integer, intent(in)::i,nen,length,NB,nphiz
real(8), intent(in)::Lx,en(nen)
integer:: ie,j,ib,jb,ikz
real(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
    do j = 1,length-1
        tr=0.0d0          
        do ib=1,nb  
          do jb=1,nb        
            do ikz=1,nphiz
              tr = tr+ cur((j-1)*nb+ib,j*nb+jb,ie,ikz)
            enddo
          enddo                        
        end do
        write(11,'(3E18.4)') dble(j-1)*Lx, en(ie), tr
    end do
    write(11,*)    
end do
close(11)
end subroutine write_current_spectrum_summed_over_kz

! write spectrum into file (pm3d map)
subroutine write_spectrum_summed_over_kz(dataset,i,G,nen,en,nkz,length,NB,Lx,coeff)
character(len=*), intent(in) :: dataset
complex(8), intent(in) :: G(:,:,:,:)
integer, intent(in)::i,nen,length,NB,nkz
real(8), intent(in)::Lx,en(nen),coeff(2)
integer:: ie,j,ib,ikz
complex(8)::tr
open(unit=11,file=trim(dataset)//TRIM(STRING(i))//'.dat',status='unknown')
do ie = 1,nen
    do j = 1,length
        tr=0.0d0          
        do ib=1,nb
          do ikz=1,nkz
            tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie,ikz)            
          enddo
        enddo
        tr=tr/dble(nkz)
        write(11,'(4E18.4)') dble(j-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
    end do
    write(11,*)    
enddo
close(11)
end subroutine write_spectrum_summed_over_kz


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


! Sancho-Rubio 
subroutine sancho(nm,E,S00,H00,H10,G00,GBB)
  complex(8), parameter :: alpha = cmplx(1.0d0,0.0d0)
  complex(8), parameter :: beta  = cmplx(0.0d0,0.0d0)
  integer i,j,k,nm,nmax
  COMPLEX(8) :: z
  real(8) :: E,error
  REAL(8) :: TOL=1.0D-10  ! [eV]
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
  nmax=200
  z = cmplx(E,1.0d-5)
  Id=0.0d0
  tmp=0.0d0
  do i=1,nm
     Id(i,i)=1.0d0
     tmp(i,i)=cmplx(0.0d0,1.0d0)
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


subroutine subspace_invert(nm,A,add_nm,in_method)
integer, intent(in) :: nm ! dimension of subspace 
integer, intent(in) :: add_nm ! additional dimension 
complex(8), intent(inout) :: A(nm,nm)
character(len=*), intent(in),optional :: in_method
character(len=10) :: method
complex(8), dimension(:,:), allocatable :: invA,V00,V10,sigmal,sigmar,S00,G00,Gbb,B
integer::nm_lead,nm_dev
if (present(in_method)) then
  method=in_method
else
  method='direct'
endif

select case (trim(method))
  case('direct')    
    allocate(invA(nm+2*add_nm,nm+2*add_nm))
    invA=cmplx(0.0d0,0.0d0)
    invA(add_nm+1:add_nm+nm,add_nm+1:add_nm+nm)=A    
    invA(1:add_nm,add_nm+1:2*add_nm)=A(1:add_nm,add_nm+1:2*add_nm)
    invA(add_nm+1:2*add_nm,1:add_nm)=A(add_nm+1:2*add_nm,1:add_nm)    
    invA(add_nm+nm+1:2*add_nm+nm,nm+1:add_nm+nm)=A(nm-add_nm+1:nm,nm-2*add_nm+1:nm-add_nm)
    invA(nm+1:add_nm+nm,add_nm+nm+1:2*add_nm+nm)=A(nm-2*add_nm+1:nm-add_nm,nm-add_nm+1:nm)
    call invert(invA,nm+2*add_nm)
    A=invA(add_nm+1:add_nm+nm,add_nm+1:add_nm+nm)
    deallocate(invA)
  case('sancho')       
    allocate(V00(add_nm,add_nm))
    allocate(V10(add_nm,add_nm))
    allocate(S00(add_nm,add_nm))
    allocate(G00(add_nm,add_nm))
    allocate(GBB(add_nm,add_nm))
    allocate(B(add_nm,add_nm))
    allocate(sigmal(add_nm,add_nm))
    allocate(sigmar(add_nm,add_nm))
    nm_lead=add_nm
    nm_dev=nm
    ! get OBC on left  
    V00 = - A(1:nm_lead,1:nm_lead) 
    V10 = - A(nm_lead+1:2*nm_lead,1:nm_lead)
    call identity(S00,nm_lead)
    call sancho(nm_lead,0.0d0,S00,V00,V10,G00,GBB)
    call zgemm('n','n',nm_lead,nm_lead,nm_lead,cone,V10,nm_lead,G00,nm_lead,czero,B,nm_lead) 
    call zgemm('n','c',nm_lead,nm_lead,nm_lead,cone,B,nm_lead,V10,nm_lead,czero,sigmal,nm_lead)  
    ! get OBC on right
    call sancho(nm_lead,0.0d0,S00,V00,transpose(conjg(V10)),G00,GBB)
    call zgemm('c','n',nm_lead,nm_lead,nm_lead,cone,V10,nm_lead,G00,nm_lead,czero,B,nm_lead) 
    call zgemm('n','n',nm_lead,nm_lead,nm_lead,cone,B,nm_lead,V10,nm_lead,czero,sigmar,nm_lead)  
    !    
    A(1:nm_lead,1:nm_lead) = A(1:nm_lead,1:nm_lead) + sigmal
    A(nm_dev-nm_lead+1:nm_dev,nm_dev-nm_lead+1:nm_dev) = A(nm_dev-nm_lead+1:nm_dev,nm_dev-nm_lead+1:nm_dev) + sigmar
    !
    call invert(A,nm_dev)
    deallocate(V00,V10,sigmar,sigmal,G00,Gbb,B,S00)
  
end select  

end subroutine subspace_invert


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


FUNCTION STRING(inn)
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

! Z(nop) = sum_ie X(ie) Y(ie-nop)
function corr1d(n,X,Y,method) result(Z)
integer, intent(in)::n
character(len=*),intent(in)::method
complex(8),intent(in)::X(n),Y(n)
complex(8)::Z(2*n-1)
integer::i,ie
select case (trim(method))
  case('index')
    Z=czero
    do i=-n+1,n-1
      do ie = max(i+1,1),min(n,n+i)
        Z(i+n)=Z(i+n) + X(ie)*Y(ie-i)
      enddo
    enddo
  case('simple')
    do i=-n+1,n-1
      Z(i+n)=sum(X(max(1,1+i):min(n+i,n))*Y(max(1-i,1):min(n,n-i)))
    enddo
  case('fft')  
  
end select
end function corr1d

! Z(ie) = sum_nop X(ie-nop) Y(nop)
function conv1d(n,X,Y,method) result(Z)
integer, intent(in)::n
character(len=*),intent(in)::method
complex(8),intent(in)::X(n),Y(n*2-1)
complex(8)::Z(n)
complex(8),allocatable,dimension(:)::x_in
integer::i,ie
select case (trim(method))
  case('index')
    Z=czero
    do ie=1,n
      do i= -n+1,n-1
        if ((ie .gt. max(i,1)).and.(ie .lt. min(n,(n+i)))) then
          Z(ie)=Z(ie) + X(ie-i)*Y(i+n)
        endif
      enddo
    enddo
  case('simple')
    do i=1,n
      Z(i)=sum(X(n:1:-1)*Y(i:i+n-1))
    enddo
  case('fft')
    allocate(X_in(n*2-1))
    X_in(1:n)=X
    X_in(n+1:n*2-1)=czero
    call do_mkl_dfti_conv(n,X_in,Y,Z)
    deallocate(X_in)
end select
end function conv1d

subroutine do_mkl_dfti_conv(n,X_in,Y_in,Z_out)
! 1D complex to complex
Use MKL_DFTI
integer :: n
Complex(8) :: X_in(n),Y_in(n),Z_out(n)
Complex(8) :: X_out(n),Y_out(n),Z_in(n)
type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
Integer :: Status
! Perform a complex to complex transform
Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n )
Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
Status = DftiCommitDescriptor( My_Desc1_Handle )
Status = DftiComputeForward( My_Desc1_Handle, X_in, X_out )
Status = DftiComputeForward( My_Desc1_Handle, Y_in, Y_out )
!
Z_in(:) = X_out(:) * Y_out(:)
!
Status = DftiComputeBackward( My_Desc1_Handle, Z_in, Z_out )
Status = DftiFreeDescriptor(My_Desc1_Handle)
end subroutine do_mkl_dfti_conv


! calculate matrix blocks for the Open Boundary Condition of W
subroutine get_OBC_blocks_for_W(n,v_00,v_01,pR_00,pR_01,pL_00,pL_01,pG_00,pG_01,NBC,&
    V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
    LL00,LL01,LL10,LG00,LG01,LG10)
integer,intent(in)::n,NBC
complex(8),intent(in),dimension(n,n)::v_00,v_01,pR_00,pR_01,pL_00,pL_01,pG_00,pG_01
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
    PR10=transpose(PR01)
    !
    PL00=pL_00;
    PL01=pL_01;
    PL10=-transpose(conjg(PL01))
    !
    PG00=pG_00
    PG01=pG_01
    PG10=-transpose(conjg(PG01))
        
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
    PR00(n+1:2*n,1:n)=transpose(pR_01)
    PR00(n+1:2*n,n+1:2*n)= pR_00
    PR01=czero
    PR01(n+1:2*n,1:n)=pR_01
    PR10=transpose(PR01)
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
! ---------
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
    matmul(matmul(VNN,PGNN1),VN1N) + matmul(matmul(VNN1,PGN1N),VNN) + matmul(matmul(VNN1,PGNN),VN1N)
    
  ! WR/WL/WG OBC Left
  call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
  ! WR/WL/WG OBC right
  call open_boundary_conditions(NR,MNN,MNN1,MN1N,VN1N,xRNN,dMNN,dVNN,condR)
  allocate(VV(nm_dev,nm_dev))
  VV = V
  if (condL<1.0d-6) then   
      !
      call get_dL_OBC_for_W(NL,xR11,LL00,LL01,LG00,LG01,M10,'L', dLL11,dLG11)
      !
      M(1:LBsize,1:LBsize)=M(1:LBsize,1:LBsize) - dM11
      VV(1:LBsize,1:LBsize)=V(1:LBsize,1:LBsize) - dV11    
      LL(1:LBsize,1:LBsize)=LL(1:LBsize,1:LBsize) + dLL11
      LG(1:LBsize,1:LBsize)=LG(1:LBsize,1:LBsize) + dLG11    
  endif
  if (condR<1.0d-6) then    
      !
      call get_dL_OBC_for_W(NR,xRNN,LLNN,LLN1N,LGNN,LGN1N,MNN1,'R', dLLNN,dLGNN)
      !
      M(NT-RBsize+1:NT,NT-RBsize+1:NT)=M(NT-RBsize+1:NT,NT-RBsize+1:NT) - dMNN
      VV(NT-RBsize+1:NT,NT-RBsize+1:NT)=V(NT-RBsize+1:NT,NT-RBsize+1:NT)- dVNN
      LL(NT-RBsize+1:NT,NT-RBsize+1:NT)=LL(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLLNN
      LG(NT-RBsize+1:NT,NT-RBsize+1:NT)=LG(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLGNN    
  endif
  !!!! calculate W^r = (I - V P^r)^-1 V    
  call invert(M,nm_dev) ! M -> xR
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
  call invert(M,nm_dev) ! M -> xR
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,V,nm_dev,czero,WR,nm_dev)           
  ! calculate W^< and W^> = W^r P^<> W^r dagger
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PL,nm_dev,czero,B,nm_dev) 
  call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WL,nm_dev) 
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PG,nm_dev,czero,B,nm_dev) 
  call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WG,nm_dev)  
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


end module green
