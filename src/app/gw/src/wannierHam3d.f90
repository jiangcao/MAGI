! Copyright (c) 2023 Jiang Cao, ETH Zurich 
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE. 
!

MODULE wannierHam3d

IMPLICIT NONE 

private 
COMPLEX(8), ALLOCATABLE :: Hr(:,:,:,:,:)
COMPLEX(8), ALLOCATABLE :: rmn(:,:,:,:,:,:) ! position operator
COMPLEX(8), ALLOCATABLE :: pmn(:,:,:,:,:,:) ! momentum operator
COMPLEX(8), ALLOCATABLE :: Vmn(:,:,:,:,:) ! Coulomb operator
real(8), allocatable :: wannier_center(:,:)
REAL(8), DIMENSION(3) :: alpha,beta,gamm,xhat,yhat,zhat,b1,b2
REAL(8) :: Lx, Ly,Lz,cell(3,3) ! in Ang
REAL(8) :: CBM, VBM, Eg ! CB minimum and VB maximum
REAL(8) :: kt_CBM, kt_VBM ! transverse kt value of CBM and VBM location
INTEGER :: ymin,ymax,xmin,xmax,zmin,zmax,nb,nx,ny,nz,nvb, spin_deg
! neighbor index range 
! nb is the number of Wannier orbitals in a unit cell
COMPLEX(8), PARAMETER :: zzero = dcmplx(0.0d0,0.0d0)
COMPLEX(8), PARAMETER :: z1j = dcmplx(0.0d0,1.0d0)
REAL(8), PARAMETER :: pi = 3.14159265359d0
REAL(8), PARAMETER :: m0=5.6856D-16 ! eV s2 / cm2
REAL(8), PARAMETER :: hbar=6.58211899D-16 ! eV s
REAL(8), PARAMETER :: c0=2.998d8 ! m/s

public :: w90_free_memory, w90_load_from_file,w90_MAT_DEF
public :: NB, Lx, Ly,Lz, Nvb, VBM, CBM, kt_CBM, kt_VBM, Eg, spin_deg
public :: eig,cross,eigv,b1,b2,norm,wannier_center,alpha,beta, invert
public :: w90_MAT_DEF_full_device, w90_bare_coulomb_full_device
public :: w90_bare_coulomb_blocks
public :: w90_momentum_full_device,w90_momentum_blocks

CONTAINS

SUBROUTINE w90_free_memory()
    deallocate(Hr)
    if (allocated(rmn)) deallocate(rmn)
    if (allocated(pmn)) deallocate(pmn)
    if (allocated(Vmn)) deallocate(Vmn)
    deallocate(wannier_center)
END SUBROUTINE w90_free_memory

SUBROUTINE w90_load_from_file(fid,lreorder_axis,axis)
use, intrinsic :: iso_fortran_env
include "mpif.h"
integer, intent(in) :: fid
logical, intent(in), optional :: lreorder_axis
integer, intent(in), optional :: axis(3)
integer :: n,i,nkx,nky,nkz
character(len=40) :: line, comment
REAL(8) :: dky,dkz,aux1(3),aux2(3,3)
integer, allocatable :: ind(:)
REAL(8), allocatable :: ham(:,:), energ(:,:), aux3(:,:)

integer(kind=int32) :: rank
integer(kind=int32) :: ierror
! Get the individual process (rank)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

    read(fid,*) nvb, spin_deg ! number of VBs, spin-degeneracy
    read(fid,*) comment
    read(fid,*) alpha
    read(fid,*) beta
    read(fid,*) gamm
    cell(:,1) = alpha
    cell(:,2) = beta
    cell(:,3) = gamm
    read(fid,*) comment
    if (present(lreorder_axis).and.(lreorder_axis)) then
        aux2(:,:) = cell(:,axis)
        cell = aux2
        alpha = cell(:,1)
        beta  = cell(:,2)
        gamm  = cell(:,3)   
    end if
    read(fid,*) n
    allocate(ham(n,7))    
    do i=1,n
        read(fid,*) ham(i,:)
    end do    
    if (present(lreorder_axis).and.(lreorder_axis)) then
        allocate(aux3(n,7))
        aux3 = ham
        aux3(:,1:3) = ham(:,axis)
        ham = aux3        
        deallocate(aux3)
    end if
    xmin=minval(ham(:,1))
    xmax=maxval(ham(:,1))
    ymin=minval(ham(:,2))
    ymax=maxval(ham(:,2))
    zmin=minval(ham(:,3))
    zmax=maxval(ham(:,3))
    nb = maxval(ham(:,4))
    nx = xmax-xmin+1
    ny = ymax-ymin+1
    nz = zmax-zmin+1
    xhat = alpha/norm(alpha)
    yhat = - cross(xhat,gamm)
    yhat = yhat/norm(yhat)
    zhat = gamm/norm(gamm)
    Ly=abs(dot_product(beta,yhat)); ! L is in unit of A
    Lx=abs(dot_product(alpha,xhat));
    Lz=norm(gamm)
if (rank == 0) then    
    print '(A40)', 'reading Wannier H from file, info:'
    print '(8a5)', 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax', 'nb', 'nvb'
    print '(8i5)', xmin, xmax, ymin, ymax,zmin,zmax, nb, nvb
    print '(3a5)', 'Lx', 'Ly', 'Lz'
    print '(3f5.1)', Lx, Ly, Lz   
endif    
    b1=cross(beta,gamm)/dot_product(alpha,cross(beta,gamm))
    b2 = cross(gamm,alpha)/dot_product(beta,cross(gamm,alpha))    
    allocate(Hr(nb,nb,nx,ny,nz))
    Hr(:,:,:,:,:) = zzero
    do i = 1,n 
        Hr(floor(ham(i,4)),floor(ham(i,5)),floor(ham(i,1))-xmin+1,&
        &floor(ham(i,2))-ymin+1,floor(ham(i,3))-zmin+1) = ham(i,6) + z1j*ham(i,7);
    end do    
    !print *, 'Find CBM and VBM'
    nkx=floor((1.0d0/Lx)/0.01d0)+1
    nky=floor((1.0d0/Ly)/0.01d0)+1    
    nkz=floor((1.0d0/Lz)/0.01d0)+1
    if (nx==1) nkx=1
    if (ny==1) nky=1
    if (nz==1) nkz=1
    if (ny .eq. 1) then
        nky = 1
    end if
    if (nx .eq. 1) then
        nkx = 1
    end if
    if (nz .eq. 1) then
        nkz = 1
    end if
    if (nky > 1) then 
        dky = 1.0d0 / dble(nky-1) * 2 * pi / Ly
    else 
        dky = 2 * pi / Ly
    end if
    if (nkz > 1) then 
        dkz = 1.0d0 / dble(nkz-1) * 2 * pi / Lz
    else 
        dkz = 2 * pi / Lz
    end if
    allocate(energ(nb,nkx*nky*nkz))
    allocate(ind(nb))
    call w90_PLOT_BZ(nkx,nky,nkz,energ)
        
    ind = maxloc(energ,2)
    kt_vbm = mod(ind(nvb),nkx)*dky - pi / Ly ! k transverse corresponding to VBMax      
    ind = minloc(energ,2)
    kt_cbm = mod(ind(nvb+1),nkx)*dky - pi / Ly ! k transverse corresponding to CBM
    
    if (ny .eq. 1) then
        kt_cbm=0.0d0
        kt_vbm=0.0d0
    end if
    ind = maxloc(energ,2)
    VBM = energ(nvb,ind(nvb))
    ind = minloc(energ,2)
    CBM = energ(nvb+1,ind(nvb+1))
    Eg = CBM-VBM
if (rank == 0) then    
    print '(3A8)','CBM','VBM','Eg'
    print '(3f8.3)', CBM, VBM, CBM-VBM
    print '(3A8)','kt_CBM','kt_VBM','2pi/Ly'
    print '(2f8.3)', kt_cbm/(2.0*pi/Ly), kt_vbm/(2.0*pi/Ly)
    print '(A40)', 'reading Wannier centers from file'
endif    
    read(fid,*) comment    
    allocate(wannier_center(3,nb))
    do i=1,nb
        read(fid,*) wannier_center(:,i)
    end do
    if (present(lreorder_axis).and.(lreorder_axis)) then
        allocate(aux3(3,nb))
        aux3 = wannier_center
        aux3(1:3,:) = wannier_center(axis,:)
        wannier_center = aux3        
        deallocate(aux3)
    end if
    ! bring wannier_center into the first unit-cell, only works for orth.
    ! to-do: make this general
    do i=1,NB
        if (wannier_center(1,i)<0) wannier_center(:,i)=wannier_center(:,i)+alpha
        if (wannier_center(1,i)>Lx) wannier_center(:,i)=wannier_center(:,i)-alpha
        if (wannier_center(2,i)<0) wannier_center(:,i)=wannier_center(:,i)+beta
        if (wannier_center(2,i)>Ly) wannier_center(:,i)=wannier_center(:,i)-beta
        if (wannier_center(3,i)<0) wannier_center(:,i)=wannier_center(:,i)+gamm
        if (wannier_center(3,i)>Lz) wannier_center(:,i)=wannier_center(:,i)-gamm
    end do
    ! center z around 0
    wannier_center(3,:)=wannier_center(3,:)-(sum(wannier_center(3,:)/dble(NB)))
    deallocate(ham)
    deallocate(energ)
    deallocate(ind)
END SUBROUTINE w90_load_from_file


SUBROUTINE w90_PLOT_BZ(nkx,nky,nkz,EN)
implicit none
integer, intent(in) :: nkx,nky,nkz
real(8), dimension(nb,nkx*nky*nkz), intent(out) :: EN
integer :: i,j,l
real(8) :: dkx, dky, dkz,kx, ky, kz
complex(8), dimension(NB,NB) :: Hii    
    if (nkx > 1) then
        dkx = 1.0d0 / dble(nkx-1) * 2 * pi / Lx
    else 
        dkx = pi / Lx
    end if
    if (nky > 1) then
        dky = 1.0d0 / dble(nky-1) * 2 * pi / Ly
    else
        dky = pi / Ly
    end if
    if (nkz > 1) then
        dkz = 1.0d0 / dble(nkz-1) * 2 * pi / Lz
    else
        dkz = pi / Lz
    end if
    do i = 1,nkx
      do j = 1,nky
        do l = 1,nkz
          kx = i*dkx - pi / Lx
          ky = j*dky - pi / Ly
          kz = l*dkz - pi / Lz
        
          call w90_MAT_DEF_3D(Hii, kx,ky,kz)    
        
          EN(1:NB,i+(j-1)*nkx+(l-1)*nkx*nky) = eig(NB,Hii)
        enddo
      enddo
    enddo
END SUBROUTINE w90_PLOT_BZ

!!! construct the diagonal and off-diagonal blocks H(I,I), H(I+1,I)
SUBROUTINE w90_MAT_DEF(Hii,H1i,kx, ky,kz,ns)
! ky in [2pi/Ang]
implicit none
integer, intent(in) :: ns
COMPLEX(8), INTENT(OUT), DIMENSION(NB*ns,NB*ns) :: Hii, H1i
real(8), intent(in) :: ky,kx,kz
integer :: i,j,k,l
real(8), dimension(3) :: kv, r
Hii(:,:) = zzero
H1i(:,:) = zzero
do i = 1,ns
    do k = 1,ns    
        do j = ymin,ymax
          do l = zmin,zmax
            kv = kx*xhat + ky*yhat + kz*zhat
            r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
            if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) )           
            end if                 
            if (((i-k+ns) <= xmax) .and. ((i-k+ns) >= xmin)) then                   
                H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + & 
                & Hr(:,:,i-k-xmin+ns+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) )            
            end if
          enddo
        end do                 
    end do
end do  
END SUBROUTINE w90_MAT_DEF

!!! construct the fully periodic Hamiltonian matrix
SUBROUTINE w90_MAT_DEF_3D(Hii,kx,ky,kz)
! kx and ky is in unit of [2pi/Ang]
implicit none
REAL(8), INTENT(IN) :: kx, ky, kz
COMPLEX(8), INTENT(OUT), DIMENSION(NB,NB) :: Hii
real(8), dimension(3) :: kv, r
integer :: i,j,l
    Hii(:,:) = zzero    
    kv = kx*xhat + ky*yhat + kz*zhat
    do i = xmin,xmax           
      do j = ymin,ymax
        do l = zmin,zmax
            r = dble(i)*alpha + dble(j)*beta + dble(l)*gamm                         
            Hii(:,:) = Hii(:,:) + Hr(:,:,i-xmin+1,j-ymin+1,l-zmin+1) &
              &* exp(-z1j * dot_product(r,kv))             
        enddo
      enddo      
    enddo  
    !!! Hii = (Hii + transpose(conjg(Hii))) / 2.0d0     
END SUBROUTINE w90_MAT_DEF_3D


!!! construct the full-device Hamiltonian Matrix
SUBROUTINE w90_MAT_DEF_full_device(Ham,ky,kz,length,NS)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky,kz
complex(8), intent(inout), dimension(NB*length,NB*length) :: Ham
integer :: i,j, k,l
real(8), dimension(3) :: kv, r
complex(8) :: phi
Ham = dcmplx(0.0d0,0.0d0)
do i = 1, length
  do k = 1, length
    do j = ymin,ymax
      do l = zmin,zmax
        kv = ky*yhat + kz*zhat
        r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
        phi = dcmplx( 0.0d0, - dot_product(r,kv) )
        if (present(NS)) then
          if ((i-k <= min(NS,xmax) ) .and. (i-k >= max(-NS,xmin) )) then                
            Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
          endif                 
        else
          if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
            Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
          endif                 
        endif
      enddo
    enddo
  enddo
enddo
END SUBROUTINE w90_MAT_DEF_full_device



!!! construct the diagonal and off-diagonal blocks V(I,I), V(I+1,I)
SUBROUTINE w90_bare_coulomb_blocks(Hii,H1i,kx, ky,kz,eps,r0,NS,ldiag)
! ky in [2pi/Ang]
implicit none
integer, intent(in) :: ns
COMPLEX(8), INTENT(OUT), DIMENSION(NB*ns,NB*ns) :: Hii, H1i
real(8), intent(in) :: ky,kx,kz
real(8), intent(in) :: eps ! dielectric constant
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
logical, intent(in) :: ldiag ! include diagonal 
integer :: i,j,k,l
real(8), dimension(3) :: kv, r
Hii(:,:) = zzero
H1i(:,:) = zzero
kv = kx*xhat + ky*yhat + kz*zhat
do i = 1,ns
    do k = 1,ns    
        do j = ymin,ymax
          do l = zmin,zmax            
            r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                               
            Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                    & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )                                   
            H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + & 
                    & bare_coulomb(i-k+ns,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )                        
          enddo
        end do        
    end do
end do  
END SUBROUTINE w90_bare_coulomb_blocks

!!! construct the bare Coulomb Matrix for the full-device
SUBROUTINE w90_bare_coulomb_full_device(V,ky,kz,length,eps,r0,ldiag,NS,method)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky,kz, eps ! dielectric constant / to reduce V
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
logical, intent(in) :: ldiag ! include diagonal 
character(len=*),intent(in),optional :: method
character(len=20) :: cmethod
complex(8), intent(out), dimension(NB*length,NB*length) :: V
integer :: i,j, k,l
real(8), dimension(3) :: kv, r
V = dcmplx(0.0d0,0.0d0)
if (present(method)) then
    cmethod=method
else
    cmethod='pointlike'
endif
select case(trim(cmethod))
  case('pointlike')
    do i = 1, length
        do k = 1, length
            do j = ymin,ymax
              do l = zmin,zmax
                kv = ky*yhat + kz*zhat
                r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
                if (present(NS)) then
                    if ((i-k <= NS ) .and. (i-k >= -NS )) then                
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )           
                    end if                 
                else
                    V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                       & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )                           
                end if              
              enddo
            end do
        end do
    end do
  case('fromfile')
    call w90_load_coulomb_blocks
    do i = 1, length
        do k = 1, length
            do j = ymin,ymax
              do l = zmin,zmax
                kv = ky*yhat + kz*zhat
                r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
                if (present(NS)) then
                    if ((i-k <= NS ) .and. (i-k >= -NS )) then                
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & Vmn(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) ) / eps          
                    endif                 
                else
                    if (abs(i-k)<4)then
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & Vmn(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp(-z1j* dot_product(r,kv) ) / eps                           
                    else
                      V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & bare_coulomb(i-k,j,l,eps,r0,ldiag) * exp(-z1j* dot_product(r,kv) )    
                    endif
                endif
              enddo
            enddo
        enddo
    enddo
end select
END SUBROUTINE w90_bare_coulomb_full_device
! function to calculate the bare coulomb potential for wannier orbitals between the (0,0) and (a1,a2) cells
FUNCTION bare_coulomb(a1,a2,a3,eps,r0,ldiag)
implicit none
integer, intent(in) :: a1, a2, a3
real(8), dimension(NB,NB) :: bare_coulomb
real(8), intent(in) :: eps ! dielectric constant
real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
logical, intent(in) :: ldiag ! include diagonal 
real(8), parameter :: pi=3.14159265359d0
real(8), parameter :: e=1.6d-19            ! charge of an electron (C)
real(8), parameter :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)
real(8) :: r(3),normr
real(8) :: maxV
integer :: i,j  
do i=1,NB
    do j=1,NB
        r = dble(a1)*alpha + dble(a2)*beta + dble(a3)*gamm + wannier_center(:,i) - wannier_center(:,j)
        normr = norm(r)
        if (normr >0.0d0) then
          bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*normr*1.0d-10) * tanh(normr/r0)  ! in eV
        else
          if (ldiag) then
            bare_coulomb(i,j) = (e)/(4.0d0*pi*epsilon0*eps*1.0d-10) * (1.0d0/r0) ! self-interaction 
          else
            bare_coulomb(i,j) = 0.0d0
          endif
        endif
    end do
end do 
END FUNCTION bare_coulomb


!!! construct the diagonal and off-diagonal blocks P(I,I), P(I+1,I)
SUBROUTINE w90_momentum_blocks(Hii,H1i,kx,ky,kz,NS,method)
! ky in [2pi/Ang]
implicit none
integer, intent(in) :: ns
COMPLEX(8), INTENT(OUT), DIMENSION(NB*ns,NB*ns,3) :: Hii, H1i ! momentum matrix block [eV]
real(8), intent(in) :: ky,kx,kz
character(len=*),intent(in)::method
integer :: i,j,k,l,v
real(8), dimension(3) :: kv, r
complex(8) :: phi
Hii(:,:,:) = zzero
H1i(:,:,:) = zzero
call calc_momentum_operator(method)
kv = kx*xhat + ky*yhat + kz*zhat
do v=1,3 ! cart direction
  do i = 1,ns
    do k = 1,ns    
      do j = ymin,ymax
        do l = zmin,zmax            
          r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm        
          phi = dcmplx( 0.0d0, - dot_product(r,kv) ) 
          if ((i-k <= xmax ) .and. (i-k >= xmin )) then                      
            Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) = Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) + &
                  & pmn(v,:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
          endif
          if (((i-k+ns) <= xmax) .and. ((i-k+ns) >= xmin)) then                                          
            H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) = H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) + & 
                  & pmn(v,:,:,i-k-xmin+1+NS,j-ymin+1,l-zmin+1) * exp( phi )                        
          endif
        enddo
      enddo        
    enddo
  enddo  
enddo
END SUBROUTINE w90_momentum_blocks


SUBROUTINE w90_momentum_full_device(Ham,ky,kz,length,NS,method)
implicit none
integer, intent(in) :: length
integer, intent(in), optional :: NS
real(8), intent(in) :: ky,kz
complex(8), intent(inout), dimension(NB*length,NB*length,3) :: Ham ! momentum matrix [eV]
character(len=*),intent(in)::method
integer :: i,j, k,v,l
real(8), dimension(3) :: kv, r
complex(8) :: phi
Ham = dcmplx(0.0d0,0.0d0)
call calc_momentum_operator(method)
do v=1,3 ! cart direction
  do i = 1, length
    do k = 1, length
      do j = ymin,ymax
        do l = zmin,zmax
          kv = ky*yhat + kz*zhat
          r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
          phi = dcmplx( 0.0d0, - dot_product(r,kv) )
          if (present(NS)) then
              if ((i-k <= min(NS,xmax) ) .and. (i-k >= max(-NS,xmin) )) then                
                  Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) + &
                  & pmn(v,:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
              end if                 
          else
              if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                  Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) = Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb,v) + &
                  & pmn(v,:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
              end if                 
          end if
        enddo
      end do
    end do
  end do
enddo
END SUBROUTINE w90_momentum_full_device
!
SUBROUTINE calc_momentum_operator(method)
implicit none
character(len=*),intent(in)::method
integer::io,jo,ix,iy,iz,mx,my,mz,mo
real(8)::r(3)
complex(8)::pre_fact
if (.not.(allocated(pmn))) allocate(pmn(3,NB,NB,nx,ny,nz)) ! [eV]
pre_fact=1.0d-8 * dcmplx(0.0d0,1.0d0) *m0/hbar*1.0d2 * c0  ! multiply light-speed so c0*pmn in energy eV 
select case (method)
  case ('approx')
  ! use wannier centers, point-like orbitals
    pmn = 0.0d0
    do io=1,NB
      do jo=1,NB
        do ix=xmin,xmax
          do iy=ymin,ymax
            do iz=zmin,zmax
              r = dble(ix)*alpha + dble(iy)*beta + dble(iz)*gamm + wannier_center(:,io) - wannier_center(:,jo)
              pmn(:,io,jo,ix-xmin+1,iy-ymin+1,iz-zmin+1)=Hr(io,jo,ix-xmin+1,iy-ymin+1,iz-zmin+1)*r
            enddo
          enddo
        enddo
      enddo
    enddo
    pmn=pmn*pre_fact
  case ('exact')
  ! use position operator : im_0/hbar sum_{R'l} H_{nl}(R-R') r_{lm}(R') - r_{nl}(R-R') H_{lm}(R')
    pmn = 0.0d0
    do io=1,NB
      do jo=1,NB
        do ix=xmin,xmax
          do iy=ymin,ymax
            do iz=zmin,zmax
              do mo=1,NB
                do mx=xmin,xmax
                  do my=ymin,ymax
                    do mz=zmin,zmax
                      if (((ix-mx)>=xmin).and.((ix-mx)<=xmax).and.((iy-my)>=ymin).and.((iy-my)<=ymax).and.((iz-mz)>=zmin).and.((iz-mz)<=zmax)) then
                        pmn(:,io,jo,ix-xmin+1,iy-ymin+1,iz-zmin+1)=pmn(:,io,jo,ix-xmin+1,iy-ymin+1,iz-zmin+1)&
                        &+Hr(io,mo,ix-mx-xmin+1,iy-my-ymin+1,iz-zmin+1)*rmn(:,mo,jo,mx-xmin+1,my-ymin+1,mz-zmin+1)&
                        &-rmn(:,io,mo,ix-mx-xmin+1,iy-my-ymin+1,iz-mz-zmin+1)*Hr(mo,jo,mx-xmin+1,my-ymin+1,mz-zmin+1)
                      endif
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
    pmn=pmn*pre_fact
  case default
    print *, 'Unknown method!!'
    call abort
end select
END SUBROUTINE

FUNCTION norm(vector)
implicit none
REAL(8) :: vector(3),norm
norm = sqrt(dot_product(vector,vector))
END FUNCTION

FUNCTION cross(a, b)
implicit none
REAL(8), DIMENSION(3) :: cross
REAL(8), DIMENSION(3), INTENT(IN) :: a, b
cross(1) = a(2) * b(3) - a(3) * b(2)
cross(2) = a(3) * b(1) - a(1) * b(3)
cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

FUNCTION eig(NN, A)
implicit none
INTEGER, INTENT(IN) :: NN
COMPLEX(8), INTENT(INOUT), DIMENSION(:,:) :: A
REAL(8) :: eig(NN)
real(8) :: W(1:NN)
integer :: INFO,LWORK,liwork, lrwork
complex(8), allocatable :: work(:)
real(8), allocatable :: RWORK(:)
!integer, allocatable :: iwork(:) 
lwork= max(1,2*NN-1)
lrwork= max(1,3*NN-2)
allocate(work(lwork))
allocate(rwork(lrwork))

CALL zheev( 'N','U', NN, A, NN, W, WORK, LWORK, RWORK, INFO )

deallocate(work,rwork)
if (INFO.ne.0)then
   write(*,*)'SEVERE WARNING: ZHEEV HAS FAILED. INFO=',INFO
   call abort
endif
eig(:)=W(:)
END FUNCTION eig

FUNCTION eigv(NN, A)
implicit none
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
  implicit none      
  integer :: info,lda,lwork,nn      
  integer, dimension(:), allocatable :: ipiv
  complex(8), dimension(nn,nn),intent(inout) :: A
  complex(8), dimension(:), allocatable :: work
  allocate(work(nn*nn))
  allocate(ipiv(nn))
  call zgetrf(nn,nn,A,nn,ipiv,info)
  call zgetri(nn,A,nn,ipiv,work,nn*nn,info)
  deallocate(work)
  deallocate(ipiv)
end subroutine invert

subroutine w90_load_coulomb_blocks()
implicit none
integer::i
real(8)::reV(NB)
real(8),parameter::ry2ev=13.6d0
real(8)::tmp(NB,NB)
allocate(Vmn(NB,NB,nx,ny,nz))
Vmn=dcmplx(0.0d0,0.0d0)

open(unit=10,file='V_CNT_0_0_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+1,1,1)=dcmplx(reV,0.0d0)*ry2ev
    Vmn(i,i,-xmin+1,1,1)=dcmplx(max(reV(i)*ry2ev,0.01),0.0d0)
enddo
close(10)
!! symmetrize V
tmp=transpose(conjg(Vmn(:,:,-xmin+1,1,1)))
Vmn(:,:,-xmin+1,1,1)=(Vmn(:,:,-xmin+1,1,1)+tmp(:,:))/2.0d0

open(unit=10,file='V_CNT_0_1_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+2,1,1)=dcmplx(reV,0.0d0)*ry2ev    
enddo
close(10)

open(unit=10,file='V_CNT_0_2_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+3,1,1)=dcmplx(reV,0.0d0)*ry2ev    
enddo
close(10)

open(unit=10,file='V_CNT_0_3_dat',status='unknown')
do i=1,NB
    read(10,*) reV(1:NB)
    Vmn(:,i,-xmin+4,1,1)=dcmplx(reV,0.0d0)*ry2ev    
enddo
close(10)

do i=1,3
    Vmn(:,:,-xmin+1-i,1,1) = transpose(conjg(Vmn(:,:,-xmin+1+i,1,1)))
enddo
end subroutine w90_load_coulomb_blocks

END MODULE wannierHam3d
