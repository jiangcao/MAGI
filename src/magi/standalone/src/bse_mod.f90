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

!
! a part of this module is an implementation of PRB 94, 245434 (2016) B. Scharf
! et al.
! 

module bse_mod
use linalg,only:norm
use parameters_mod,only:dp,twopi,pi
implicit none 
private
public::gaussian,screenedpot,barepot,efieldpotential,rijvaluesfunc,dHdk,exciton_wavefunction_simple,exciton_wavefunction_grid

contains

real(dp) function gaussian(x,mu,sig) 
implicit none
! Gaussian function
real(dp), intent(in)::x,mu,sig
    gaussian = 1.0d0/sig/sqrt(twopi)*exp(-1.0d0/2.0d0*(x-mu)**2/sig**2)
end function gaussian

real(dp) function screenedpot(r,r0,epsilon0,epsilon,e)
!Function to get the screened potential W(Rijvvd) in eV given Rijvvd, r0, epsilon
!epsilon0, e
implicit none
real(dp),intent(in) :: r(1:3),r0,epsilon0,epsilon,e
real(dp)::aux3,aux4,aux0,wrijvvd
if (norm(r) .lt. 1.0e-5) then ! Ang
    screenedpot=0.0d0
else
    aux3=norm(r)/r0
    aux4=aux3/(1.0d0+aux3);
    aux0=0.577216d0-log(2.0d0);
    wrijvvd=-(log(aux4)+aux0*exp(-aux3))*((e)/(4*pi*epsilon0*r0*1d-10));
   !wrijvvd=((e)/(8.*epsilon0.*r0.*1e-10)).*(StruveH0((epsilon.*norm(r))/r0)-bessely(0,((epsilon.*norm(r))/r0))); % in eV, equation 4
   screenedpot = wrijvvd   
end if
end function screenedpot



real(dp) function barepot(r,epsilon0,epsilon,e)
!Function to get the bare potential V(Rijvvd) in eV given rihvvd, epsilon0,e and
!epsilon
implicit none
real(dp),intent(in) :: r(3),epsilon0,epsilon,e
real(dp)::vrijvvd
if (norm(r) .lt. 1.0e-5) then ! Ang
    vrijvvd=0.0d0
else
    vrijvvd=(e)/(4.0d0*pi*epsilon0*epsilon*norm(r)*1.0d-10);  ! in eV
end if
barepot=vrijvvd
end function barepot

! real(dp) function efieldpotential(r,F,L,e,k0,xhat)
! function to get the electric field contribution to the potential
! real(dp),intent(in)::r(3),F,L,e,k0,xhat(3)
! efieldpotential=e*F*dot_product(r,xhat)*(tanh(k0*(0.25d0-(dot_product(r,xhat)/L)**2))) ! implementation of equation 7
! end function efieldpotential


real(dp) function efieldpotential(rvvd,e1,F,k0,i,j,alpha,beta,NKX,NKY)
!Function to get the electric field contribution to the potential
real(dp),intent(in)::alpha(3),beta(3),k0,rvvd(3)
real(dp),intent(in)::F ! E field strength
real(dp),intent(in)::e1(3)  ! unit vector along the E field
integer,intent(in)::i,j,NKX,NKY
real(dp)::minRv,minStark,Rvect(3),R,Stark,sintheta
real(dp)::normR
integer::sci,scj
normR=norm(rvvd)
if (normR.eq.0) then
    efieldpotential=0.0d0
else
    minRv = dot_product(e1,rvvd)/norm(e1)
    sintheta=sin(acos(minRv/normR))
    minStark=-F*minRv*(tanh(k0*(0.25d0-(minRv/dble(normR*NKX*sintheta))**2)))
    do sci=-1,1           ! super-cell index (loop over nearest-neighbors) 
        do scj=-1,1       ! to find the minimum R_ij under the periodic boundary
            Rvect(:)=dble(i+sci*NKX)*alpha(:)+dble(j+scj*NKY)*beta(:)
            R=dot_product(Rvect,e1)/norm(e1)
            Stark=-F*R*(tanh(k0*(0.25d0-(R/dble(normR*NKX*sintheta))**2)))
            if (abs(Stark) .lt. abs(minStark) ) then
                minStark=Stark
            end if
        end do
    end do
    efieldpotential=minStark
end if
end function efieldpotential

subroutine rijvaluesfunc(r,i,j,rij,iv,ivd,nb,wanniercenter,lwcenter)
    !Funtion to get the Rijvvd values given indices i j k l iv and ivd and the
    !locations of the wannier centers
    implicit none
    real(dp),intent(in)::rij(3),wanniercenter(1:3,1:nb)
    integer,intent(in)::i,j,iv,ivd,nb   
    real(dp),intent(out)::r(3)
    logical,intent(in)::lwcenter
    if (lwcenter) then
        if ((i.eq.0).and.(j.eq.0)) then
            r(:)=rij(:)
        else
            r(:)=rij(:) +(wanniercenter(:,iv)-wanniercenter(:,ivd))
        end if
    else
        r(:)=rij(:) 
    endif
end subroutine rijvaluesfunc

! compute the dipole-matrix element by finite-difference method dH/dkx at
! kx=(H(kx+1)-H(kx))/(k(x+1)-k(x))
complex(dp) function dHdk(nb,nkx,nky,ivb,icb,kx,ky,ak,Hiimat,nvb,kv)
implicit none
integer,intent(in) :: nb,ivb,icb,kx,ky,nvb,nkx,nky
complex(dp),intent(in) :: ak(nb,nb,nkx,nky),Hiimat(nkx,nky,nb,nb)
real(dp),intent(in) :: kv(3,nkx,nky)
real(dp)::dk
integer::v,vd
complex(dp),dimension(nb,nb)::Hvvd,Hvvdprev,dkHvvd
Hvvd=Hiimat(kx,ky,:,:)
if(kx .eq. 1) then
    Hvvdprev=Hvvd
    Hvvd=Hiimat(2,ky,:,:)
    dk=norm(kv(:,kx+1,ky)-kv(:,kx,ky))
else
    Hvvdprev=Hiimat(kx-1,ky,:,:)    
    dk=norm(kv(:,kx,ky)-kv(:,kx-1,ky))
end if
dkHvvd=(Hvvd-Hvvdprev)/dk; ! H is in eV, k is in 1/Ang, so dH/dk is in eV.Ang
dHdk=0.0d0
do v=1,nb
  do vd=1,nb
    dHdk = dHdk + conjg(ak(v,ivb,kx,ky))*dkHvvd(v,vd)*ak(vd,nvb+icb,kx,ky)  ! dvck is in eV.Ang
  end do
end do
end function dHdk

subroutine exciton_wavefunction_simple(asvckmat,s,NKX,NKY,ncb,nvb,nx,nb,nvbtot,ak,chi,kv,rij)
    implicit none
    complex(dp),intent(in)::asvckmat(nvb,ncb,NKX,NKY,nx)
    complex(dp),intent(in)::ak(nb,nb,NKX,NKY)
    real(dp),intent(in)::kv(3,NKX,NKY),rij(3,NKX,NKY)
    integer,intent(in)::s,NKX,NKY,ncb,nvb,nx,nb,nvbtot
    complex(dp),intent(out)::chi(NB,NB,NKX,NKY) ! exciton wavefunction
    integer::v,vd,icb,ivb,ikx,iky,ix,iy,vb0
    real(dp)::Ri(3)
    complex(dp)::phi
    vb0=nvbtot-nvb
    chi=dcmplx(0.0d0,0.0d0)
    !$omp parallel default(none) private(iy,ix,Ri,ikx,iky,phi,v,vd,ivb,icb) shared(NKX,NKY,nb,nvb,ncb,kv,rij,chi,ak,asvckmat,s,nvbtot,vb0)
    !$omp do
    do iy=1,NKY
        do ix=1,NKX
            Ri(:)=rij(:,ix,iy)
            do ikx=1,NKX
                do iky=1,NKY
                    phi=exp(dcmplx(0.0d0, dot_product(kv(:,ikx,iky),Ri)))
                    do v=1,nb
                        do vd=1,nb
                            do ivb=1,nvb
                                do icb=1,ncb
                                    chi(v,vd,ix,iy)=chi(v,vd,ix,iy)+asvckmat(ivb,icb,ikx,iky,s)*conjg(ak(v,ivb+vb0,ikx,iky))*ak(vd,icb+nvbtot,ikx,iky)*phi
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
end subroutine exciton_wavefunction_simple

subroutine exciton_wavefunction_grid(NKX,NKY,nb,chi,chixyz,rij,grid,npt,rsmear,wannier_center)
    implicit none    
    real(dp),intent(in)::rij(3,NKX,NKY)
    integer,intent(in)::NKX,NKY,nb,npt
    real(dp),intent(in)::grid(3,npt) ! real-space grid for generating exciton wavefunction
    real(dp),intent(in)::rsmear ! smearing radius
    real(dp),intent(in)::wannier_center(3,nb) ! wannier function centers
    complex(dp),intent(in)::chi(NB,NB,NKX,NKY) ! exciton wavefunction
    complex(dp),intent(out)::chixyz(npt) ! exciton wavefunction on the grid    
    integer::ix,iy,v,vd,ikx,iky,ipt,ngauss,ig
    real(dp)::r,g    
    chixyz=dcmplx(0.0d0,0.0d0)
    !$omp parallel default(none) private(ipt,ix,iy,v,vd,r,g,ikx,iky) shared(NKX,NKY,nb,chi,chixyz,rij,rsmear,grid,wannier_center,npt)
    !$omp do
    do ipt=1,npt 
        do ix=1,NKX
            do iy=1,NKY                            
                do v=1,nb
                    do vd=1,nb
                        r=norm(grid(:,ipt) - (rij(:,ix,iy)+wannier_center(:,vd)-wannier_center(:,v))) 
                        if (r .lt. 10.0d0*rsmear) then                    
                            chixyz(ipt)=chixyz(ipt)+chi(v,vd,ix,iy)*gaussian(r,0.0d0,rsmear)                            
                        end if
                    end do
                end do
            end do
        end do    
    end do
    !$omp end do
    !$omp end parallel    
end subroutine exciton_wavefunction_grid

end module bse_mod
