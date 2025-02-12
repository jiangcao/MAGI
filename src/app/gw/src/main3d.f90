! Copyright (c) 2023 Jiang Cao, ETH Zurich 
! All rights reserved.
!
PROGRAM main

USE wannierHam3d, only : NB, w90_load_from_file, w90_free_memory,Ly,Lz,CBM,VBM,eig,w90_MAT_DEF_full_device, invert, Lx, w90_bare_coulomb_full_device,kt_CBM,spin_deg,Eg,w90_MAT_DEF,w90_bare_coulomb_blocks,w90_momentum_blocks
use green, only : green_calc_g, green_solve_gw_3D
use gw_dense, only : green_solve_gw_1D
use green_rgf, only : green_rgf_solve_gw_1d,green_RGF_CMS,green_rgf_solve_gw_3d,green_rgf_solve_gw_ephoton_3d, green_rgf_solve_gw_ephoton_3d_mpi_kpar, green_rgf_solve_gw_ephoton_3d_ijs
implicit none

real(8), parameter :: pi=3.14159265359d0
integer :: NS, nm, ie, ne, width,nkx,nky,i,j,k,axis,num_B,ib,ncpu,xyz(3),length,num_vac,iky,ik
integer,allocatable:: orb_vac(:)
real(8) :: ky, emax, emin
real(8), allocatable::phix(:),ek(:,:),B(:),en(:)
complex(8), allocatable :: H00(:,:),H10(:,:),H01(:,:),tmpV(:,:)
complex(8), allocatable,dimension(:,:,:,:) :: Hii,H1i,Vii,V1i
complex(8), allocatable,dimension(:,:,:,:,:) :: Pii,P1i
complex(8), allocatable :: H00ld(:,:,:,:),H10ld(:,:,:,:),T(:,:,:,:),V(:,:,:),Ham(:,:,:),invV(:,:,:)
complex(8), allocatable,dimension(:,:,:,:) :: G_retarded,G_lesser,G_greater
complex(8), allocatable,dimension(:,:,:,:) :: P_retarded,P_lesser,P_greater
complex(8), allocatable,dimension(:,:,:,:) :: W_retarded,W_lesser,W_greater
complex(8), allocatable,dimension(:,:,:,:) :: Sig_retarded,Sig_lesser,Sig_greater
complex(8), allocatable,dimension(:,:,:,:) :: Sig_retarded_new,Sig_lesser_new,Sig_greater_new
complex(8), allocatable :: Pmn(:,:,:,:)

logical :: reorder_axis, ltrans, lreadpot, lqdot, lkz, lephot, lnogw, lrgf,ldiag,lrcoulomb,labs
character(len=5)::conv_method
integer :: nen , pottype
complex(8), parameter :: cone = dcmplx(1.0d0,0.0d0)
complex(8), parameter :: czero  = dcmplx(0.0d0,0.0d0)
real(8), allocatable :: pot(:)
integer, allocatable :: cell_index(:,:)
integer :: nm_dev, iter, niter, nkz,ikz,ndiag,nk
real(8) :: eps_screen, mud,mus,temps,tempd, alpha_mix, dkz,kz, r0,potscale,encut(2),dky
real(8) :: intensity,hw,midgap(2),polaris(3), ky_shift, kz_shift
real(8) :: scba_tol
real(8), allocatable :: mid_bandgap(:)
real(8), allocatable,dimension(:) ::charge

! MPI variables
integer ( kind = 4 ) ierr
integer ( kind = 4 ) comm_size
integer ( kind = 4 ) comm_rank
integer ( kind = 4 ) local_Nkz
integer ( kind = 4 ) local_Nky
integer ( kind = 4 ) local_NE
integer ( kind = 4 ) first_local_energy

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
      read(10,*) potscale, pottype
    endif
    read(10,*) niter
    read(10,*) scba_tol
    read(10,*) eps_screen
    read(10,*) r0
    read(10,*) mus,mud
    read(10,*) temps, tempd
    read(10,*) alpha_mix
    read(10,*) encut(1:2) ! intraband interband cutoff energies for P and W
    read(10,*) lrcoulomb
    read(10,*) ldiag
    read(10,*) lrgf
    read(10,*) conv_method
    read(10,*) lkz
    if (lkz) then
      read(10,*) nky,nkz
      read(10,*) ky_shift, kz_shift
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

call MPI_Barrier(MPI_COMM_WORLD, ierr)

if (mod(nen, comm_size) /= 0) then
  if (comm_rank == 0) then
    print *, 'Number of energies (', nen, ') must be divisible by comm_size (', comm_size, ')'
  endif
  call MPI_Finalize(ierr)
  stop
endif

local_Nkz = nkz
local_Nky = nky
local_NE = nen / comm_size
first_local_energy = local_NE * comm_rank + 1

open(unit=10,file='ham_dat',status='unknown')
call w90_load_from_file(10, reorder_axis,xyz)
close(10)

call omp_set_num_threads(ncpu)

if (ltrans) then    
  if (.not.lrgf) then    

      if (.not. lkz) then ! 1d case
        nky=1
        nkz=1
      endif
      print *, 'Build the full device H'
      print *, 'length=',length
      allocate(Ham(nb*length,nb*length,nky*nkz))
      allocate(  V(nb*length,nb*length,nky*nkz))      
      allocate(H00ld(nb*NS,nb*NS,2,nky*nkz))
      allocate(H10ld(nb*NS,nb*NS,2,nky*nkz))
      allocate(T(nb*ns,nb*length,2,nky*nkz))
      T = dcmplx(0.0d0,0.0d0)
      if (nkz>1) then
        dkz=2.0d0*pi/Lz / dble(nkz)
      else
        dkz=pi/Lz
      endif
      if (nky>1) then
        dky=2.0d0*pi/Ly / dble(nky)
      else
        dky=pi/Ly
      endif
      ! contact Ham blocks
      open(unit=11,file='ek.dat',status='unknown')
      do iky=1,nky
        ky=-pi/Ly + dble(iky)*dky
        do ikz=1,nkz
          kz=-pi/Lz + dble(ikz)*dkz
          ik=ikz+(iky-1)*nkz
          call w90_MAT_DEF(H00ld(:,:,1,ik),H10ld(:,:,1,ik),0.0d0, ky,kz,NS)
          !
          H10ld(:,:,2,ik) = H10ld(:,:,1,ik)
          H10ld(:,:,1,ik) = transpose(conjg(H10ld(:,:,1,ik)))
          H00ld(:,:,2,ik) = H00ld(:,:,1,ik)              
          T(1:nb*ns,1:nb*ns,1,ik) = H10ld(:,:,1,ik)
          T(1:nb*ns,nb*(length-ns)+1:nb*length,2,ik) = H10ld(:,:,2,ik)    
          ! write bands into ek.dat        
          allocate(phix(nkx))
          allocate(ek(nb*ns,nkx))
          allocate(H01(nb*ns,nb*ns))
          allocate(H00(nb*ns,nb*ns))
          phix=(/(i, i=1,nkx, 1)/) / dble(nkx-1) * pi * 2.0d0 - pi        
          do i=1,nkx
              H00(:,:) = H00ld(:,:,1,ik)            
              H01(:,:) = transpose(H10ld(:,:,1,ik))
              H01(:,:) = conjg(H01(:,:))
              H00 = H00 + exp(+dcmplx(0.0d0,1.0d0)*phix(i))*H10ld(:,:,1,ik)
              H00 = H00 + exp(-dcmplx(0.0d0,1.0d0)*phix(i))*H01               
              ek(:,i) = eig(NB*NS,H00)
          enddo    
          deallocate(H00,H01)    
          do j=1,nb*NS
              do i=1,nkx
                  write(11,'(4E18.8)') ky*Ly, kz*Lz, phix(i), ek(j,i)
              enddo
              write(11,*)
          enddo            
          deallocate(ek,phix)        
        enddo
      enddo
      close(11)

      allocate(en(nen))
      allocate(G_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(G_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(G_greater(nb*length,nb*length,nen,nky*nkz))      
      
      allocate(P_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(P_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(P_greater(nb*length,nb*length,nen,nky*nkz))    
      
      allocate(W_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(W_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(W_greater(nb*length,nb*length,nen,nky*nkz))    
      
      allocate(Sig_retarded(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_lesser(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_greater(nb*length,nb*length,nen,nky*nkz))
      
      allocate(Sig_retarded_new(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_lesser_new(nb*length,nb*length,nen,nky*nkz))
      allocate(Sig_greater_new(nb*length,nb*length,nen,nky*nkz))
      
      Sig_retarded = dcmplx(0.0d0,0.0d0)
      Sig_lesser = dcmplx(0.0d0,0.0d0)
      Sig_greater = dcmplx(0.0d0,0.0d0)
      
      en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin

      do iky=1,nky
        ky=-pi/Ly + dble(iky)*dky
        do ikz=1,nkz
          kz=-pi/Lz + dble(ikz)*dkz
          ik=ikz+(iky-1)*nkz
          ! device Ham matrix
          call w90_MAT_DEF_full_device(Ham(:,:,ik),ky,kz,length)      
          ! Coulomb operator
          call w90_bare_coulomb_full_device(V(:,:,ik),ky,kz,length,eps_screen,r0,ldiag)      
        enddo
      enddo
      open(unit=11,file='V.dat',status='unknown')
      do i=1, size(V,1)
          do j=1, size(V,2)
              write(11,'(2I6,2F15.4)') i,j, dble(V(i,j,nky*nkz/2+1)), aimag(V(i,j,nky*nkz/2+1))
          end do
          write(11,*)
      end do
      close(11)
      
      open(unit=11,file='Ham.dat',status='unknown')
      do i=1, size(Ham,1)
          do j=1, size(Ham,2)
              write(11,'(2I6,2F15.4)') i,j, dble(Ham(i,j,nky*nkz/2+1)), aimag(Ham(i,j,nky*nkz/2+1))
          end do
          write(11,*)
      end do
      close(11)
      nm_dev=nb*length    
      allocate(mid_bandgap(nm_dev))
      mid_bandgap(:) = (CBM+VBM)/2.0d0      
      
      allocate(charge(nm_dev))
      
      ! add on potential    
      if (lreadpot) then
        open(unit=10,file='pot_dat',status='unknown')
        if (pottype == 1) then
          allocate(pot(length))
          pot(:) = 0.0d0
          do i = 1,length
              read(10,*) pot(i)
          enddo
        else 
          allocate(pot(length*NB))
          pot(:) = 0.0d0
          do i = 1,length*NB
              read(10,*) pot(i)
          enddo
        endif
        close(10)
        pot = pot*potscale
        if (pottype == 1) then
           do j = 1,length
               do ib = 1,nb
                   Ham((j-1)*nb+ib,(j-1)*nb+ib,:)=Ham((j-1)*nb+ib,(j-1)*nb+ib,:) + pot(j)
                   mid_bandgap((j-1)*nb+ib) = mid_bandgap((j-1)*nb+ib) + pot(j)
               end do               
           end do      
        else  
           do j = 1,length
               do ib = 1,nb
                   Ham((j-1)*nb+ib,(j-1)*nb+ib,:)=Ham((j-1)*nb+ib,(j-1)*nb+ib,:) + pot((j-1)*nb+ib)
                   mid_bandgap((j-1)*nb+ib) = mid_bandgap((j-1)*nb+ib) + pot((j-1)*nb+ib)
               end do
           end do
        endif
        do ib = 1,nb*NS
            H00ld(ib,ib,1,:)=H00ld(ib,ib,1,:)+pot(1)
            H00ld(ib,ib,2,:)=H00ld(ib,ib,2,:)+pot(length)
        end do
        deallocate(pot)
      endif
      write(20,*) mid_bandgap      
      !
      if ((nkz == 1) .and. (nky==1)) then
      
        call green_solve_gw_1D(scba_tol,niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,conv_method,mid_bandgap,&
          alpha_mix,nen,En,nb,ns,Ham,H00ld,H10ld,T,V,&
          G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
          W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
          Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,charge)
          
      else
      
        call green_solve_gw_3D(niter,nm_dev,Lx,length,dble(spin_deg),temps,tempd,mus,mud,&
          alpha_mix,nen,En,nb,ns,nky,nkz,Ham,H00ld,H10ld,T,V,&
          G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
          W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
          Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag)
          
      endif

   else
    ! Long device, use RGF


    if (comm_rank == 0) then

      print *, '~~~~~~~~~~~~~~~~~ RGF ~~~~~~~~~~~~~~~~~'    
      print *, 'Build the full device H'
      print *, 'length=',length*NS,'uc =',length*NS*Lx/1.0d1,'(nm)'
    endif
    if (.not. lkz) then
      nky=1
      nkz=1
      ky_shift=0.0d0
      kz_shift=0.0d0
    endif
    
  
    nm=NB*NS
    nk=nky*nkz
    allocate(Hii(nm,nm,length,nk))
    allocate(H1i(nm,nm,length,nk))
!    allocate(V(nm*length,nm*length,nky*nkz))   
    allocate(Vii(nm,nm,length,nk))
    allocate(V1i(nm,nm,length,nk))    
    allocate(Pii(nm,nm,3,length,nk))
    allocate(P1i(nm,nm,3,length,nk))    
    if (nkz>1) then
      dkz=2.0d0*pi/Lz / dble(nkz)
    else
      dkz=pi/Lz
    endif
    if (nky>1) then
      dky=2.0d0*pi/Ly / dble(nky)
    else
      dky=pi/Ly
    endif

    if (comm_rank == 0) then

      open(unit=10,file='Hii.dat',status='unknown')
      open(unit=11,file='H1i.dat',status='unknown')
      open(unit=20,file='Vii.dat',status='unknown')
      open(unit=21,file='V1i.dat',status='unknown')
    endif

    do iky=1,nky
      ky=-pi/Ly + dble(iky)*dky + ky_shift*2.0d0*pi/Ly
      do ikz=1,nkz
        kz=-pi/Lz + dble(ikz)*dkz + kz_shift*2.0d0*pi/Ly
        ik=ikz+(iky-1)*nkz

        if (comm_rank == 0) then
          write(10,'(A,I6,2F15.4)') '# ik',ik,ky*Ly,kz*Lz
          write(11,'(A,I6,2F15.4)') '# ik',ik,ky*Ly,kz*Lz
          write(20,'(A,I6,2F15.4)') '# ik',ik,ky*Ly,kz*Lz
          write(21,'(A,I6,2F15.4)') '# ik',ik,ky*Ly,kz*Lz
        endif

        ! get Ham blocks
        call w90_MAT_DEF(Hii(:,:,1,ik),H1i(:,:,1,ik),0.0d0, ky,kz,NS)
        !
        call w90_bare_coulomb_blocks(Vii(:,:,1,ik),V1i(:,:,1,ik),0.0d0,ky,kz,eps_screen,r0,ns,ldiag)
        !
        call w90_momentum_blocks(Pii(:,:,:,1,ik),P1i(:,:,:,1,ik),0.0d0,ky,kz,NS,'approx')

        ! write Ham blocks

        if (comm_rank == 0) then

          do i=1,nm
            do j=1,nm
              write(10,'(2I6,2F15.6)') i,j,dble(Hii(i,j,1,ik)),aimag(Hii(i,j,1,ik))
              write(11,'(2I6,2F15.6)') i,j,dble(H1i(i,j,1,ik)),aimag(H1i(i,j,1,ik))
              write(20,'(2I6,2F15.6)') i,j,dble(Vii(i,j,1,ik)),aimag(Vii(i,j,1,ik))
              write(21,'(2I6,2F15.6)') i,j,dble(V1i(i,j,1,ik)),aimag(V1i(i,j,1,ik))
            enddo
          enddo
          write(10,*)
          write(11,*)
          write(20,*)
          write(21,*)
        endif

        ! build device Ham 
        do i=2,length
          Hii(:,:,i,ik)=Hii(:,:,1,ik)
          H1i(:,:,i,ik)=H1i(:,:,1,ik)
          !
          Vii(:,:,i,ik)=Vii(:,:,1,ik)
          V1i(:,:,i,ik)=V1i(:,:,1,ik)
          !
          Pii(:,:,:,i,ik)=Pii(:,:,:,1,ik)
          P1i(:,:,:,i,ik)=P1i(:,:,:,1,ik)
        enddo
      enddo
    enddo


    if (comm_rank == 0) then

      close(10)
      close(11)
      close(20)
      close(21)

      open(unit=11,file='ek.dat',status='unknown')
    endif

    !

    allocate(H00(nb*ns,nb*ns))
    allocate(phix(nkx))
    allocate(ek(nb*ns,nkx))
    do iky=1,nky
      ky=-pi/Ly + dble(iky)*dky +ky_shift*2.0d0*pi/Ly
      do ikz=1,nkz
        kz=-pi/Lz + dble(ikz)*dkz +kz_shift*2.0d0*pi/Ly
        ik=ikz+(iky-1)*nkz          
        ! write bands into ek.dat                                      
        phix=(/(i, i=1,nkx, 1)/) / dble(nkx-1) * pi * 2.0d0 - pi        
        do i=1,nkx
            H00(:,:) = Hii(:,:,1,ik)                          
            H00 = H00 + exp(+dcmplx(0.0d0,1.0d0)*phix(i))*H1i(:,:,1,ik)
            H00 = H00 + exp(-dcmplx(0.0d0,1.0d0)*phix(i))*conjg(transpose(H1i(:,:,1,ik)))
            ek(:,i) = eig(NB*NS,H00)
        enddo  
        
        if (comm_rank == 0) then
          do j=1,nb*NS
              do i=1,nkx
                  write(11,'(4E18.8)') ky*Ly, kz*Lz, phix(i), ek(j,i)
              enddo
              write(11,*)
          enddo
        endif  
                    
      enddo
    enddo
    deallocate(H00)    
    deallocate(ek,phix)
    
    if (comm_rank == 0) then
      close(11)
    endif


    !
    if (lreadpot) then
      open(unit=10,file='pot_dat',status='unknown')
      if (pottype == 1) then
        allocate(pot(length*NS))
        pot(:) = 0.0d0
        do i = 1,length*NS
            read(10,*) pot(i)
        enddo
      else 
        allocate(pot(length*NS*NB))
        pot(:) = 0.0d0
        do i = 1,length*NS*NB
            read(10,*) pot(i)
        enddo
      endif
      close(10)
      pot=pot*potscale
      ! add on potential    
      if (pottype == 1) then
        do ik=1,nk
          do j = 1,length
            do k = 1,NS
              do ib = 1,nb
                Hii(ib+(k-1)*NB,ib+(k-1)*NB,j,ik)=Hii(ib+(k-1)*NB,ib+(k-1)*NB,j,ik) + pot((j-1)*NS+k)
              enddo
            enddo
          enddo      
        enddo
      else
        do ik=1,nk
          do j = 1,length
            do k = 1,NS
              do ib = 1,nb
                Hii(ib+(k-1)*NB,ib+(k-1)*NB,j,ik)=Hii(ib+(k-1)*NB,ib+(k-1)*NB,j,ik) + pot((j-1)*NS*nb+(k-1)*nb+ib)
              enddo
            enddo
          enddo      
        enddo
      endif
      deallocate(pot)
    endif         

    allocate(en(nen))
    en=(/(i, i=1,nen, 1)/) / dble(nen) * (emax-emin) + emin

    if (ldiag) then
      ndiag=0
    else
      ndiag=NB
    endif

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    if (comm_rank == 0) then
      print *, 'Starting simulation ...'
    endif

    if (lephot) then

      if (comm_size == 1) then
        call green_rgf_solve_gw_ephoton_3d(alpha_mix,niter,NB,NS,nm,length,nky,nkz,ndiag,Lx,nen,en,(/temps,tempd/),(/mus,mud/),Hii,H1i,Vii,V1i,dble(spin_deg),Pii,P1i,polaris,intensity,hw,labs)
      else
        call green_rgf_solve_gw_ephoton_3d_ijs(alpha_mix,niter,NB,NS,nm,length,nky,nkz,ndiag,Lx,nen,en,(/temps,tempd/),(/mus,mud/),Hii,H1i,Vii,V1i,dble(spin_deg),Pii,P1i,polaris,intensity,hw,labs,lnogw, comm_size, comm_rank, local_NE, first_local_energy)

      endif
    else
      if (comm_size == 1) then
        call green_rgf_solve_gw_3d(alpha_mix,niter,NB,NS,nm,length,nky,nkz,ndiag,Lx,nen,en,(/temps,tempd/),(/mus,mud/),Hii,H1i,Vii,V1i,dble(spin_deg))
      else
        call green_rgf_solve_gw_ephoton_3d_ijs(alpha_mix,niter,NB,NS,nm,length,nky,nkz,ndiag,Lx,nen,en,(/temps,tempd/),(/mus,mud/),Hii,H1i,Vii,V1i,dble(spin_deg),Pii,P1i,polaris,0.0d0,1.0d0,.false.,lnogw, comm_size, comm_rank, local_NE, first_local_energy)        
      endif
    endif
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
  endif        
end if

call w90_free_memory

if (comm_rank == 0) then
  print *, 'End of program'
endif
call MPI_Finalize( ierr )

END PROGRAM main
