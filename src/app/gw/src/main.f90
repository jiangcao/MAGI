PROGRAM main
USE wannierHam, only : NB, w90_load_from_file, w90_free_memory,Ly, w90_MAT_DEF, CBM,VBM,eig,w90_MAT_DEF_ribbon_simple, w90_ribbon_add_peierls, w90_MAT_DEF_full_device, invert, Lx,w90_MAT_DEF_dot,w90_dot_add_peierls, w90_bare_coulomb_full_device,kt_CBM,spin_deg,w90_momentum_full_device,w90_bare_coulomb_blocks,Eg
use green, only : green_calc_g, green_solve_gw_1D,green_solve_gw_2D,green_solve_ephoton_freespace_1D,green_solve_gw_ephoton_1D,green_solve_gw_1D_memsaving,green_solve_gw_2D_memsaving,green_solve_gw_1D_supermemsaving
use green_rgf, only : green_rgf_solve_gw_1d,green_RGF_CMS
implicit none
real(8), parameter :: pi=3.14159265359d0
integer :: NS, nm, ie, ne, width,nkx,i,j,k,axis,num_B,ib,ncpu,xyz(3),length,num_vac
integer,allocatable:: orb_vac(:)
real(8) :: ky, emax, emin
real(8), allocatable::phix(:),ek(:,:),B(:),en(:)
complex(8), allocatable :: H00(:,:),H10(:,:),H01(:,:),tmpV(:,:)
complex(8), allocatable,dimension(:,:,:) :: Hii,H1i,Vii,V1i
complex(8), allocatable :: H00ld(:,:,:,:),H10ld(:,:,:,:),T(:,:,:,:),V(:,:,:),Ham(:,:,:),invV(:,:,:)
complex(8), allocatable,dimension(:,:,:,:) :: G_retarded,G_lesser,G_greater
complex(8), allocatable,dimension(:,:,:,:) :: P_retarded,P_lesser,P_greater
complex(8), allocatable,dimension(:,:,:,:) :: W_retarded,W_lesser,W_greater
complex(8), allocatable,dimension(:,:,:,:) :: Sig_retarded,Sig_lesser,Sig_greater
complex(8), allocatable,dimension(:,:,:,:) :: Sig_retarded_new,Sig_lesser_new,Sig_greater_new
complex(8), allocatable :: Pmn(:,:,:,:)

logical :: reorder_axis, ltrans, lreadpot, lqdot, lkz, lephot, lnogw, lrgf,ldiag,lrcoulomb,labs, lvertex
integer :: nen
complex(8), parameter :: cone = cmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = cmplx(0.0d0,0.0d0)
real(8), allocatable :: pot(:)
integer, allocatable :: cell_index(:,:)
integer :: nm_dev, iter, niter, nkz,ikz,ndiag
real(8) :: eps_screen, mud,mus,temps,tempd, alpha_mix, dkz,kz, r0,potscale,encut(2)
real(8) :: intensity,hw,midgap(2),polaris(3)
num_vac=0 ! number of vacancies
open(unit=10,file='input',status='unknown')
read(10,*) ns
read(10,*) width
read(10,*) nkx
read(10,*) num_B
allocate(B(num_B))
read(10,*) B
if (B(1)>0) then
  num_vac=floor(B(1))
  allocate(orb_vac(num_vac))
  read(10,*) orb_vac ! orbit index 
endif
read(10,*) axis
read(10,*) ncpu
read(10,*) reorder_axis
if (reorder_axis) then
    read(10,*) xyz
end if    
read(10,*) lqdot
read(10,*) ltrans
if (ltrans) then
    read(10,*) length
    read(10,*) emin, emax, nen
    read(10,*) lreadpot
    if (lreadpot) then
      read(10,*) potscale
    endif
    read(10,*) niter
    read(10,*) eps_screen
    read(10,*) r0
    read(10,*) mus,mud
    read(10,*) temps, tempd
    read(10,*) alpha_mix
    read(10,*) encut(1:2) ! intraband interband cutoff energies for P and W
    read(10,*) lrcoulomb
    read(10,*) ldiag
    read(10,*) lrgf
    read(10,*) lvertex
    read(10,*) lkz
    if (lkz) then
      read(10,*) nkz
    endif
    read(10,*) lephot
    if (lephot) then
      read(10,*) hw,intensity
      read(10,*) polaris(:)
      read(10,*) lnogw
      read(10,*) labs
    endif
end if
close(10)

open(unit=10,file='ham_dat',status='unknown')
call w90_load_from_file(10, reorder_axis,xyz)
close(10)

call omp_set_num_threads(ncpu)

if (ltrans) then    
  if (.not.lrgf) then    
    if (lkz) then ! 2d case
    print *, 'Build the full device H'
      print *, 'length=',length
      allocate(Ham(nb*length,nb*length,nkz))
      allocate(   V(nb*length,nb*length,nkz))      
      allocate(pot(length))
      allocate(H00ld(nb*NS,nb*NS,2,nkz))
      allocate(H10ld(nb*NS,nb*NS,2,nkz))
      allocate(T(nb*ns,nb*length,2,nkz))
      T = dcmplx(0.0d0,0.0d0)
      dkz=2.0d0*pi/Ly / dble(nkz)
      ! contact Ham blocks
      open(unit=11,file='ek.dat',status='unknown')
      do ikz=1,nkz
        kz=-pi/Ly + dble(ikz)*dkz
        call w90_MAT_DEF(H00ld(:,:,1,ikz),H10ld(:,:,1,ikz),0.0d0, kz,NS)
        !
        H10ld(:,:,2,ikz) = H10ld(:,:,1,ikz)
        H10ld(:,:,1,ikz) = transpose(conjg(H10ld(:,:,1,ikz)))
        H00ld(:,:,2,ikz) = H00ld(:,:,1,ikz)              
        T(1:nb*ns,1:nb*ns,1,ikz) = H10ld(:,:,1,ikz)
        T(1:nb*ns,nb*(length-ns)+1:nb*length,2,ikz) = H10ld(:,:,2,ikz)    
        ! write bands into ek.dat        
        allocate(phix(nkx))
        allocate(ek(nb*ns,nkx))
        allocate(H01(nb*ns,nb*ns))
        allocate(H00(nb*ns,nb*ns))
        phix=(/(i, i=1,nkx, 1)/) / dble(nkx-1) * pi * 2.0d0 - pi        
        do i=1,nkx
            H00(:,:) = H00ld(:,:,1,ikz)            
            H01(:,:) = transpose(H10ld(:,:,1,ikz))
            H01(:,:) = conjg(H01(:,:))
            H00 = H00 + exp(+dcmplx(0.0d0,1.0d0)*phix(i))*H10ld(:,:,1,ikz)
            H00 = H00 + exp(-dcmplx(0.0d0,1.0d0)*phix(i))*H01               
            ek(:,i) = eig(NB*NS,H00)
        end do    
        deallocate(H00,H01)    
        do j=1,nb*NS
            do i=1,nkx
                write(11,'(3E18.8)') kz*Ly, phix(i), ek(j,i)
            end do
            write(11,*)
        end do            
        deallocate(ek,phix)        
      enddo
      close(11)
      pot(:) = 0.0d0
      if (lreadpot) then
          open(unit=10,file='pot_dat',status='unknown')
          do i = 1,length
              read(10,*) pot(i)
          end do
          close(10)
      end if         
      allocate(en(nen))
      allocate(G_retarded(nb*length,nb*length,nen,nkz))
      allocate(G_lesser(nb*length,nb*length,nen,nkz))
      allocate(G_greater(nb*length,nb*length,nen,nkz))      
      
!      allocate(P_retarded(nb*length,nb*length,nen,nkz))
!      allocate(P_lesser(nb*length,nb*length,nen,nkz))
!      allocate(P_greater(nb*length,nb*length,nen,nkz))    
      
!      allocate(W_retarded(nb*length,nb*length,nen,nkz))
!      allocate(W_lesser(nb*length,nb*length,nen,nkz))
!      allocate(W_greater(nb*length,nb*length,nen,nkz))    
      
      allocate(Sig_retarded(nb*length,nb*length,nen,nkz))
      allocate(Sig_lesser(nb*length,nb*length,nen,nkz))
      allocate(Sig_greater(nb*length,nb*length,nen,nkz))
      
      allocate(Sig_retarded_new(nb*length,nb*length,nen,nkz))
      allocate(Sig_lesser_new(nb*length,nb*length,nen,nkz))
      allocate(Sig_greater_new(nb*length,nb*length,nen,nkz))
      
      Sig_retarded = dcmplx(0.0d0,0.0d0)
      Sig_lesser = dcmplx(0.0d0,0.0d0)
      Sig_greater = dcmplx(0.0d0,0.0d0)
      
      en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin
      
      do ikz=1,nkz
        kz=-pi/Ly + dble(ikz)*dkz
        ! device Ham matrix
        call w90_MAT_DEF_full_device(Ham(:,:,ikz),kz,length)      
        ! Coulomb operator
        call w90_bare_coulomb_full_device(V(:,:,ikz),kz,length,eps_screen,r0,ldiag)      
      enddo
      open(unit=11,file='V.dat',status='unknown')
      do i=1, size(V,1)
          do j=1, size(V,2)
              write(11,'(2I6,2F15.4)') i,j, dble(V(i,j,nkz/2+1)), aimag(V(i,j,nkz/2+1))
          end do
          write(11,*)
      end do
      close(11)
      
      open(unit=11,file='Ham.dat',status='unknown')
      do i=1, size(Ham,1)
          do j=1, size(Ham,2)
              write(11,'(2I6,2F15.4)') i,j, dble(Ham(i,j,1)), aimag(Ham(i,j,1))
          end do
          write(11,*)
      end do
      close(11)
      
      ! add on potential    
      do j = 1,length
          do ib = 1,nb
              Ham((j-1)*nb+ib,(j-1)*nb+ib,:)=Ham((j-1)*nb+ib,(j-1)*nb+ib,:) + pot(j)*potscale
          end do
      end do      
      do ib = 1,nb*NS
          H00ld(ib,ib,1,:)=H00ld(ib,ib,1,:)+pot(1)
          H00ld(ib,ib,2,:)=H00ld(ib,ib,2,:)+pot(length)
      end do
      nm_dev=nb*length    
      
      midgap=(/ (CBM+VBM)/2.0d0,(CBM+VBM)/2.0d0 /)
      midgap= midgap + (/ pot(1), pot(length)/)
      !
!      call green_solve_gw_2D(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,&
!        alpha_mix,nen,En,nb,ns,nkz,Ham,H00ld,H10ld,T,V,&
!        G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
!        W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
!        Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag)
      call green_solve_gw_2D_memsaving(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,&
          alpha_mix,nen,En,nb,ns,nkz,Ham,H00ld,H10ld,T,V,&
          G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,&
          Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,encut,Eg)
      deallocate(pot)
      deallocate(H00ld)
      deallocate(H10ld) 
      deallocate(Ham)    
      deallocate(T)
      deallocate(G_retarded)
      deallocate(G_lesser)
      deallocate(G_greater)
!      deallocate(P_retarded)
!      deallocate(P_lesser)
!      deallocate(P_greater)
!      deallocate(W_retarded)
!      deallocate(W_lesser)
!      deallocate(W_greater)
      deallocate(Sig_retarded)
      deallocate(Sig_lesser)
      deallocate(Sig_greater)
      deallocate(Sig_retarded_new)
      deallocate(Sig_lesser_new)
      deallocate(Sig_greater_new)
      deallocate(en)    
      deallocate(V)    
    else ! 1d case
      print *, 'Build the full device H'
      print *, 'length=',length
      allocate(Ham(nb*length,nb*length,1))
      allocate(   V(nb*length,nb*length,1))      
      allocate(pot(length))
      allocate(H00ld(nb*NS,nb*NS,2,1))
      allocate(H10ld(nb*NS,nb*NS,2,1))
      allocate(T(nb*ns,nb*length,2,1))
      ! contact Ham blocks
      call w90_MAT_DEF(H00ld(:,:,1,1),H10ld(:,:,1,1),0.0d0, kt_CBM,NS)
      !
      H10ld(:,:,2,1) = H10ld(:,:,1,1)
      H10ld(:,:,1,1) = transpose(conjg(H10ld(:,:,1,1)))
      H00ld(:,:,2,1) = H00ld(:,:,1,1)      
      T = dcmplx(0.0d0,0.0d0)
      T(1:nb*ns,1:nb*ns,1,1) = H10ld(:,:,1,1)
      T(1:nb*ns,nb*(length-ns)+1:nb*length,2,1) = H10ld(:,:,2,1)      
      ! write bands into ek.dat        
      open(unit=11,file='ek.dat',status='unknown')
      allocate(phix(nkx))
      allocate(ek(nb*ns,nkx))
      allocate(H01(nb*ns,nb*ns))
      allocate(H00(nb*ns,nb*ns))
      phix=(/(i, i=1,nkx, 1)/) / dble(nkx-1) * pi * 2.0d0 - pi        
      ikz=1
      do i=1,nkx
          H00(:,:) = H00ld(:,:,1,ikz)            
          H01(:,:) = transpose(H10ld(:,:,1,ikz))
          H01(:,:) = conjg(H01(:,:))
          H00 = H00 + exp(+dcmplx(0.0d0,1.0d0)*phix(i))*H10ld(:,:,1,ikz)
          H00 = H00 + exp(-dcmplx(0.0d0,1.0d0)*phix(i))*H01               
          ek(:,i) = eig(NB*NS,H00)
      end do    
      deallocate(H00,H01)    
      do j=1,nb*NS
          do i=1,nkx
              write(11,'(3E18.8)') phix(i), ek(j,i)
          end do
          write(11,*)
      end do            
      deallocate(ek,phix)
      close(11)
      pot(:) = 0.0d0
      if (lreadpot) then
          open(unit=10,file='pot_dat',status='unknown')
          do i = 1,length
              read(10,*) pot(i)
          end do
          close(10)
      end if         
      pot=pot*potscale

      allocate(en(nen))
      
      en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin
      
      ! device Ham matrix
      call w90_MAT_DEF_full_device(Ham(:,:,1),kt_CBM,length,NS)      
      !
      ! Coulomb operator
      if (lrcoulomb) then
        call w90_bare_coulomb_full_device(V(:,:,1),0.0d0,length,eps_screen,r0,ldiag,method='fromfile')      
      else
        call w90_bare_coulomb_full_device(V(:,:,1),0.0d0,length,eps_screen,r0,ldiag,method='pointlike')      
      endif
      !
      ! apply vacancies to Ham and V
      if (num_vac>0) then
        do i=1,num_vac
          Ham(orb_vac(i),:,1)=0.0d0
          Ham(:,orb_vac(i),1)=0.0d0
          Ham(orb_vac(i),orb_vac(i),1)=1e4
          V(orb_vac(i),:,1)=0.0d0
          V(:,orb_vac(i),1)=0.0d0
        enddo
      endif
      !
      open(unit=11,file='V.dat',status='unknown')
      do i=1, size(V,1)
          do j=1, size(V,2)
              write(11,'(2I6,2E15.4)') i,j, dble(V(i,j,1)), aimag(V(i,j,1))
          end do
          write(11,*)
      end do
      close(11)
      !
      open(unit=11,file='Ham.dat',status='unknown')
      do i=1, size(Ham,1)
          do j=1, size(Ham,2)
              write(11,'(2I6,2E15.4)') i,j, dble(Ham(i,j,1)), aimag(Ham(i,j,1))
          end do
          write(11,*)
      end do
      close(11)
      !
      ! add on potential    
      do j = 1,length
          do ib = 1,nb
              Ham((j-1)*nb+ib,(j-1)*nb+ib,1)=Ham((j-1)*nb+ib,(j-1)*nb+ib,1) + pot(j)
          end do
      end do      
      do ib = 1,nb*NS
          H00ld(ib,ib,1,1)=H00ld(ib,ib,1,1)+pot(1)
          H00ld(ib,ib,2,1)=H00ld(ib,ib,2,1)+pot(length)
      end do
      nm_dev=nb*length  
      
      midgap=(/ (CBM+VBM)/2.0d0,(CBM+VBM)/2.0d0 /)
      midgap= midgap + (/ pot(1), pot(length)/)
      !  
      ! momentum operator   
      allocate(Pmn(nb*length,nb*length,3,1))
      call w90_momentum_full_device(Pmn(:,:,:,1),0.0d0,length,NS,'approx')
      !
      ! apply vacancies to P
      if (num_vac>0) then
        do i=1,num_vac
          Pmn(orb_vac(i),:,:,1)=0.0d0
          Pmn(:,orb_vac(i),:,1)=0.0d0          
        enddo
      endif
      !
      open(unit=11,file='Px.dat',status='unknown')      
      do i=1, size(Pmn,1)
          do j=1, size(Pmn,2)
              write(11,'(2I6,2E15.4)') i,j, dble(Pmn(i,j,1,1)), aimag(Pmn(i,j,1,1))
          end do
          write(11,*)
      end do
      close(11)    
      !
      open(unit=11,file='Py.dat',status='unknown')      
      do i=1, size(Pmn,1)
          do j=1, size(Pmn,2)
              write(11,'(2I6,2E15.4)') i,j, dble(Pmn(i,j,2,1)), aimag(Pmn(i,j,2,1))
          end do
          write(11,*)
      end do
      close(11)
      !
      ! e photon
      if (lephot) then
        allocate(G_retarded(nb*length,nb*length,nen,1))
        allocate(G_lesser(nb*length,nb*length,nen,1))
        allocate(G_greater(nb*length,nb*length,nen,1))
         
        allocate(Sig_retarded(nb*length,nb*length,nen,1))
        allocate(Sig_lesser(nb*length,nb*length,nen,1))
        allocate(Sig_greater(nb*length,nb*length,nen,1))
        
        allocate(Sig_retarded_new(nb*length,nb*length,nen,1))
        allocate(Sig_lesser_new(nb*length,nb*length,nen,1))
        allocate(Sig_greater_new(nb*length,nb*length,nen,1))
        
        Sig_retarded = dcmplx(0.0d0,0.0d0)
        Sig_lesser = dcmplx(0.0d0,0.0d0)
        Sig_greater = dcmplx(0.0d0,0.0d0)

        if (.not. lnogw) then
          call green_solve_gw_ephoton_1D(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,midgap,&
              alpha_mix,nen,En,nb,ns,Ham(:,:,1),H00ld(:,:,:,1),H10ld(:,:,:,1),T(:,:,:,1),V(:,:,1),&
              Pmn(:,:,:,1),polaris,intensity,hw,labs,&
              G_retarded(:,:,:,1),G_lesser(:,:,:,1),G_greater(:,:,:,1),P_retarded(:,:,:,1),P_lesser(:,:,:,1),P_greater(:,:,:,1),&
              W_retarded(:,:,:,1),W_lesser(:,:,:,1),W_greater(:,:,:,1),Sig_retarded(:,:,:,1),Sig_lesser(:,:,:,1),Sig_greater(:,:,:,1),&
              Sig_retarded_new(:,:,:,1),Sig_lesser_new(:,:,:,1),Sig_greater_new(:,:,:,1),ldiag,encut,Eg,lvertex=lvertex)
         
         else
         
           call green_solve_ephoton_freespace_1D(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,&
              alpha_mix,nen,En,nb,ns,Ham(:,:,1),H00ld(:,:,:,1),H10ld(:,:,:,1),T(:,:,:,1),&
              Pmn(:,:,:,1),polaris,intensity,hw,labs,&
              G_retarded(:,:,:,1),G_lesser(:,:,:,1),G_greater(:,:,:,1),Sig_retarded(:,:,:,1),Sig_lesser(:,:,:,1),Sig_greater(:,:,:,1),&
              Sig_retarded_new(:,:,:,1),Sig_lesser_new(:,:,:,1),Sig_greater_new(:,:,:,1))
         endif
         
      else
        !  call green_solve_gw_1D(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,&
        !    alpha_mix,nen,En,nb,ns,Ham(:,:,1),H00ld(:,:,:,1),H10ld(:,:,:,1),T(:,:,:,1),V(:,:,1),invV(:,:,1),&
        !    G_retarded(:,:,:,1),G_lesser(:,:,:,1),G_greater(:,:,:,1),P_retarded(:,:,:,1),P_lesser(:,:,:,1),P_greater(:,:,:,1),&
        !    W_retarded(:,:,:,1),W_lesser(:,:,:,1),W_greater(:,:,:,1),Sig_retarded(:,:,:,1),Sig_lesser(:,:,:,1),Sig_greater(:,:,:,1),&
        !    Sig_retarded_new(:,:,:,1),Sig_lesser_new(:,:,:,1),Sig_greater_new(:,:,:,1),ldiag)
        
        allocate(G_retarded(nb*length,nb*length,nen,1))
        allocate(G_lesser(nb*length,nb*length,nen,1))
        allocate(G_greater(nb*length,nb*length,nen,1))
         
        allocate(Sig_retarded(nb*length,nb*length,nen,1))
        allocate(Sig_lesser(nb*length,nb*length,nen,1))
        allocate(Sig_greater(nb*length,nb*length,nen,1))
        
        allocate(Sig_retarded_new(nb*length,nb*length,nen,1))
        allocate(Sig_lesser_new(nb*length,nb*length,nen,1))
        allocate(Sig_greater_new(nb*length,nb*length,nen,1))
        
        Sig_retarded = dcmplx(0.0d0,0.0d0)
        Sig_lesser = dcmplx(0.0d0,0.0d0)
        Sig_greater = dcmplx(0.0d0,0.0d0)
        
        call green_solve_gw_1D_memsaving(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,midgap,&
            alpha_mix,nen,En,nb,ns,Ham(:,:,1),H00ld(:,:,:,1),H10ld(:,:,:,1),T(:,:,:,1),V(:,:,1),&
            G_retarded(:,:,:,1),G_lesser(:,:,:,1),G_greater(:,:,:,1),Sig_retarded(:,:,:,1),Sig_lesser(:,:,:,1),Sig_greater(:,:,:,1),&
            Sig_retarded_new(:,:,:,1),Sig_lesser_new(:,:,:,1),Sig_greater_new(:,:,:,1),ldiag,encut,Eg,writeGF=.false.,lvertex=lvertex)

!        allocate(G_retarded(nb,nb,length,nen))
!        allocate(G_lesser(nb,nb,length,nen))
!        allocate(G_greater(nb,nb,length,nen))
         
!        allocate(Sig_retarded(nb,nb,length,nen))
!        allocate(Sig_lesser(nb,nb,length,nen))
!        allocate(Sig_greater(nb,nb,length,nen))
        
!        allocate(Sig_retarded_new(nb,nb,length,nen))
!        allocate(Sig_lesser_new(nb,nb,length,nen))
!        allocate(Sig_greater_new(nb,nb,length,nen))
        
!        Sig_retarded = dcmplx(0.0d0,0.0d0)
!        Sig_lesser = dcmplx(0.0d0,0.0d0)
!        Sig_greater = dcmplx(0.0d0,0.0d0)

!        call green_solve_gw_1D_supermemsaving(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,midgap,&
!          alpha_mix,nen,En,nb,ns,Ham(:,:,1),H00ld(:,:,:,1),H10ld(:,:,:,1),T(:,:,:,1),V(:,:,1),&
!          G_retarded,G_lesser,G_greater,Sig_retarded,Sig_lesser,Sig_greater,&
!          Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,encut,Eg,writeGF=.false.)
      endif
      
      
      deallocate(pot)
      deallocate(H00ld)
      deallocate(H10ld) 
      deallocate(Ham)    
      deallocate(T)
      deallocate(G_retarded)
      deallocate(G_lesser)
      deallocate(G_greater)
 !     deallocate(P_retarded)
 !     deallocate(P_lesser)
 !     deallocate(P_greater)
 !     deallocate(W_retarded)
 !     deallocate(W_lesser)
 !     deallocate(W_greater)
      deallocate(Sig_retarded)
      deallocate(Sig_lesser)
      deallocate(Sig_greater)
      deallocate(Sig_retarded_new)
      deallocate(Sig_lesser_new)
      deallocate(Sig_greater_new)
      deallocate(en)    
      deallocate(V)    
      if (lephot) deallocate(Pmn)
    
    endif 
        
  else
    ! Long device, use RGF
    print *, '~~~~~~~~~~~~~~~~~ RGF ~~~~~~~~~~~~~~~~~'
    
    print *, 'Build the full device H'
    print *, 'length=',length*NS
    nm=NB*NS
    allocate(Hii(nm,nm,length))
    allocate(H1i(nm,nm,length))
    allocate(Vii(nm,nm,length))
    allocate(V1i(nm,nm,length))    
    allocate(pot(length))    
    ! get Ham blocks
    call w90_MAT_DEF(Hii(:,:,1),H1i(:,:,1),0.0d0, kt_CBM,NS)
    ! build device Ham 
    do i=2,length
      Hii(:,:,i)=Hii(:,:,1)
      H1i(:,:,i)=H1i(:,:,1)
    enddo
    pot(:) = 0.0d0
    if (lreadpot) then
        open(unit=10,file='pot_dat',status='unknown')
        do i = 1,length
            read(10,*) pot(i)
        end do
        close(10)
    end if         
    allocate(en(nen))
    en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin
    ! add on potential    
    do j = 1,length
        do ib = 1,nb
            Hii(ib,ib,j)=Hii(ib,ib,j) + pot(j)*potscale
        end do
    end do      
    ! coulomb operator blocks
    call w90_bare_coulomb_blocks(Vii(:,:,1),V1i(:,:,1),0.0d0, 0.0d0,eps_screen,r0,ns,ldiag)
    !
    do i=2,length
      Vii(:,:,i)=Vii(:,:,1)
      V1i(:,:,i)=V1i(:,:,1)
    enddo    
    if (ldiag) then
      ndiag=1
    else
      ndiag=NB*NS
    endif
    call green_rgf_solve_gw_1d(alpha_mix,niter,NB,NS,nm,length,ndiag,Lx,nen,en,(/temps,tempd/),(/mus,mud/),Hii,H1i,Vii,V1i,dble(spin_deg))
    deallocate(Hii,H1i,Vii,V1i)
    deallocate(pot)
    deallocate(en)
  endif        
end if
if (allocated(B)) deallocate(B)
if (allocated(orb_vac)) deallocate(orb_vac)
call w90_free_memory
print *, 'End of program'
END PROGRAM main
