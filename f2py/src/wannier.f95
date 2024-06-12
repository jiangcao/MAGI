! -*- f90 -*-

! include "../../Modules/parameters.f95"

module linalg
    implicit none 
    contains

    ! calculate eigen-values of a Hermitian matrix A
    FUNCTION eig(NN, A)
        INTEGER, INTENT(IN) :: NN
        COMPLEX(8), INTENT(INOUT), DIMENSION(:, :) :: A
        ! -----
        REAL(8) :: eig(NN)
        real(8) :: W(1:NN)
        integer :: INFO, LWORK, liwork, lrwork
        complex(8), allocatable :: work(:)
        real(8), allocatable :: RWORK(:)
        !integer, allocatable :: iwork(:)
        lwork = max(1, 2*NN - 1)
        lrwork = max(1, 3*NN - 2)
        allocate (work(lwork))
        allocate (rwork(lrwork))
        !
        CALL zheev('N', 'U', NN, A, NN, W, WORK, LWORK, RWORK, INFO)
        !
        deallocate (work, rwork)
        if (INFO .ne. 0) then
            write (*, *) 'SEVERE WARNING: ZHEEV HAS FAILED. INFO=', INFO
            call abort()
        end if
        eig(:) = W(:)
    END FUNCTION eig
end module linalg


MODULE wannierHam
use linalg
use parameters_mod

IMPLICIT NONE 

CONTAINS

    SUBROUTINE load_from_file(fname,lreorder_axis,axis,nb,nx,ny,nz,hr,wannier_center,n_range,cell,length)    
        character(len=*), intent(in) :: fname
        logical, intent(in) :: lreorder_axis
        integer, intent(in) :: axis(3),nb,nx,ny,nz
        complex(8), intent(out) :: hr(nb,nb,nx,ny,nz)
        real(8), intent(out) :: wannier_center(3,nb)
        real(8), intent(out) :: cell(3,3)  
        real(8), intent(out) :: length(3)    
        integer, intent(out) :: n_range(9) ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        ! ----
        integer :: fid
        integer :: n,i,rc
        character(len=40) :: comment
        REAL(8) :: aux2(3,3)
        real(8), dimension(3):: alpha,beta,gamm,xhat,yhat,zhat    
        REAL(8), allocatable :: ham(:,:), energ(:,:), aux3(:,:)
        open(action='read', file=trim(fname), iostat=rc, newunit=fid)
        read(fid,*) n_range(8), n_range(9) ! number of VBs, spin-degeneracy
        read(fid,*) comment
        read(fid,*) alpha
        read(fid,*) beta
        read(fid,*) gamm
        cell(:,1) = alpha
        cell(:,2) = beta
        cell(:,3) = gamm
        read(fid,*) comment
        if (lreorder_axis) then        
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
        if (lreorder_axis) then
            allocate(aux3(n,7))
            aux3 = ham
            aux3(:,1:3) = ham(:,axis)
            ham = aux3        
            deallocate(aux3)
        end if
        n_range(1)=nint(maxval(ham(:,4))) ! number of WFs
        n_range(2)=nint(minval(ham(:,1))) ! xmin
        n_range(3)=nint(maxval(ham(:,1))) ! xmax
        n_range(4)=nint(minval(ham(:,2))) ! ymin
        n_range(5)=nint(maxval(ham(:,2))) ! ymax
        n_range(6)=nint(minval(ham(:,3))) ! zmin
        n_range(7)=nint(maxval(ham(:,3))) ! zmax        
        !
        Hr(:,:,:,:,:) = czero
        do i = 1,n 
            Hr(nint(ham(i,4)),nint(ham(i,5)),nint(ham(i,1))-n_range(2)+1,&
            &nint(ham(i,2))-n_range(4)+1,nint(ham(i,3))-n_range(6)+1) = &
            dcmplx( ham(i,6) , ham(i,7) )
        end do      
        read(fid,*) comment        
        do i=1,n_range(1)
            read(fid,*) wannier_center(:,i)
        end do    
        if (lreorder_axis) then
            allocate(aux3(3,n_range(1)))
            aux3 = wannier_center
            aux3(1:3,:) = wannier_center(axis,:)
            wannier_center = aux3        
            deallocate(aux3)
        end if    
        deallocate(ham)
        close(fid)
        xhat = alpha/norm(alpha)
        yhat = - cross(xhat,gamm)
        yhat = yhat/norm(yhat)
        zhat = gamm/norm(gamm)
        length(2)=abs(dot_product(beta,yhat)); ! L is in unit of A
        length(1)=abs(dot_product(alpha,xhat));
        length(3)=norm(gamm)
    END SUBROUTINE load_from_file


    !!! construct the diagonal and off-diagonal blocks H(I,I), H(I+1,I)
    SUBROUTINE block_mat_def(Hii,H1i,kx, ky,kz,nb,ns,n_range,hr,cell)
        ! ky in [2pi/Ang]
        implicit none
        integer, intent(in) :: ns,nb
        integer, intent(in) :: n_range(9) ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        COMPLEX(8), INTENT(OUT), DIMENSION(nb*ns,nb*ns) :: Hii, H1i
        COMPLEX(8), INTENT(in):: hr(:,:,:,:,:)
        real(8), intent(in) :: cell(3,3)
        real(8), intent(in) :: ky,kx,kz
        integer :: i,j,k,l
        real(8), dimension(3) :: kv, r,xhat, yhat, zhat
        Hii(:,:) = czero
        H1i(:,:) = czero
        associate(nb => n_range(1))
        associate(xmin => n_range(2))
        associate(xmax => n_range(3))
        associate(ymin => n_range(4))
        associate(ymax => n_range(5))
        associate(zmin => n_range(6))
        associate(zmax => n_range(7))
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        xhat = alpha/norm(alpha)
        yhat = - cross(xhat,gamm)
        yhat = yhat/norm(yhat)
        zhat = gamm/norm(gamm)
        do i = 1,ns
            do k = 1,ns    
                do j = ymin,ymax
                do l = zmin,zmax
                    kv = kx*xhat + ky*yhat + kz*zhat
                    r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
                    if ((i-k <= xmax ) .and. (i-k >= xmin )) then                
                        Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = &
                        Hii(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                        & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp(-c1i* dot_product(r,kv) )           
                    end if                 
                    if (((i-k+ns) <= xmax) .and. ((i-k+ns) >= xmin)) then                   
                        H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = &
                        H1i(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + & 
                        & Hr(:,:,i-k-xmin+ns+1,j-ymin+1,l-zmin+1) * exp(-c1i* dot_product(r,kv) )            
                    end if
                enddo
                end do                 
            end do
        end do  
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
    END SUBROUTINE block_mat_def



    !!! construct the full-device Hamiltonian Matrix
    SUBROUTINE full_device_mat_def(Ham,ky,kz,nb,ns,length,hr,cell,n_range)
        implicit none
        integer, intent(in) :: length
        integer, intent(in) :: NS
        integer, intent(in) :: NB
        real(8), intent(in) :: ky,kz
        integer, intent(in) :: n_range(:) ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        complex(8), intent(out), dimension(nb*length,nb*length) :: Ham
        complex(8), intent(in) :: hr(:,:,:,:,:)
        real(8), intent(in) :: cell(:,:)    
        !---
        integer :: i,j, k,l
        real(8), dimension(3) :: kv, r,xhat, yhat, zhat
        complex(8) :: phi
        Ham = dcmplx(0.0d0,0.0d0)
        associate(nb => n_range(1))
        associate(xmin => n_range(2))
        associate(xmax => n_range(3))
        associate(ymin => n_range(4))
        associate(ymax => n_range(5))
        associate(zmin => n_range(6))
        associate(zmax => n_range(7))
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        xhat = alpha/norm(alpha)
        yhat = - cross(xhat,gamm)
        yhat = yhat/norm(yhat)
        zhat = gamm/norm(gamm)
        do i = 1, length
            do k = 1, length
                do j = ymin,ymax
                    do l = zmin,zmax
                        kv = ky*yhat + kz*zhat
                        r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                   
                        phi = dcmplx( 0.0d0, - dot_product(r,kv) )        
                        if ((i-k <= min(NS,xmax) ) .and. (i-k >= max(-NS,xmin) )) then                
                            Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = &
                            Ham(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                            & Hr(:,:,i-k-xmin+1,j-ymin+1,l-zmin+1) * exp( phi )           
                        endif                         
                    enddo
                enddo
            enddo
        enddo
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
    END SUBROUTINE full_device_mat_def


    !!! construct the bare Coulomb Matrix for the full-device
    SUBROUTINE full_device_bare_coulomb(V,ky,kz,length,eps,r0,ldiag,ns,nb,method,n_range,wannier_center,cell)
        implicit none
        integer, intent(in) :: length
        integer, intent(in) :: NS,NB
        real(8), intent(in) :: ky,kz, eps ! dielectric constant / to reduce V
        real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
        logical, intent(in) :: ldiag ! include diagonal         
        integer, intent(in) :: n_range(:) ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        real(8), intent(in) :: wannier_center(:,:)
        real(8), intent(in) :: cell(:,:)
        character(len=*),intent(in) :: method        
        complex(8), intent(out), dimension(NB*length,NB*length) :: V
        integer :: i,j, k,l
        real(8), dimension(3) :: kv, r,xhat, yhat, zhat
        V = dcmplx(0.0d0,0.0d0)                
        associate(xmin => n_range(2))
        associate(xmax => n_range(3))
        associate(ymin => n_range(4))
        associate(ymax => n_range(5))
        associate(zmin => n_range(6))
        associate(zmax => n_range(7))
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        !
        xhat = alpha/norm(alpha)
        yhat = - cross(xhat,gamm)
        yhat = yhat/norm(yhat)
        zhat = gamm/norm(gamm)            
        !
        select case(trim(method))
        case('pointlike')
            do i = 1, length
                do k = 1, length
                    do j = ymin,ymax
                    do l = zmin,zmax
                        kv = ky*yhat + kz*zhat
                        r =  dble(i-k)*alpha + dble(j)*beta + dble(l)*gamm                                           
                        if ((i-k <= NS ) .and. (i-k >= -NS )) then                
                            V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) = V(((i-1)*nb+1):i*nb,((k-1)*nb+1):k*nb) + &
                                & bare_coulomb(i-k,j,l,eps,r0,nb,ldiag,wannier_center,cell) * exp(-c1i* dot_product(r,kv) )           
                        end if                                         
                    enddo
                    end do
                end do
            end do
        case('fromfile')
            
        end select
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate        
    END SUBROUTINE full_device_bare_coulomb

    ! function to calculate the bare coulomb potential for wannier orbitals between the (0,0) and (a1,a2) cells
    FUNCTION bare_coulomb(a1,a2,a3,eps,r0,nb,ldiag,wannier_center,cell)
        implicit none
        integer, intent(in) :: a1, a2, a3, nb
        real(8), dimension(NB,NB) :: bare_coulomb
        real(8), intent(in) :: eps ! dielectric constant
        real(8), intent(in) :: r0 ! length [ang] to remove singularity of 1/r
        logical, intent(in) :: ldiag ! include diagonal 
        real(8), intent(in) :: wannier_center(:,:)
        real(8), intent(in) :: cell(:,:)
        real(8), parameter :: pi=3.14159265359d0
        real(8), parameter :: e=1.6d-19            ! charge of an electron (C)
        real(8), parameter :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)
        real(8) :: r(3),normr
        real(8) :: maxV
        integer :: i,j  
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        do i=1,NB
            do j=1,NB
                r = dble(a1)*alpha + dble(a2)*beta + dble(a3)*gamm + wannier_center(:,i) - wannier_center(:,j)
                normr = norm(r)
                if (normr >0.0d0) then
                bare_coulomb(i,j) = (e_charge)/(4.0d0*pi*epsilon0*eps*normr*1.0d-10) * tanh(normr/r0)  ! in eV
                else
                if (ldiag) then
                    bare_coulomb(i,j) = (e_charge)/(4.0d0*pi*epsilon0*eps*1.0d-10) * (1.0d0/r0) ! self-interaction 
                else
                    bare_coulomb(i,j) = 0.0d0
                endif
                endif
            end do
        end do 
        end associate
        end associate
        end associate
    END FUNCTION bare_coulomb

    
    SUBROUTINE plot_bands_bz(nkx,nky,nkz,EN,hr,nb,length,cell,n_range)
        implicit none
        integer, intent(in) :: nkx,nky,nkz,nb
        real(8), intent(in) :: length(3),cell(3,3)
        complex(8), intent(in) :: hr(:,:,:,:,:)
        integer, intent(in) :: n_range(:) ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        real(8), dimension(nb,nkx*nky*nkz), intent(out) :: EN
        integer :: i,j,l
        real(8) :: dkx, dky, dkz,kx, ky, kz
        complex(8), dimension(NB,NB) :: Hii    
        associate(Lx => length(1))
        associate(Ly => length(2))
        associate(Lz => length(3))
        if (nkx > 1) then
            dkx = 3.0_dp / dble(nkx) * pi / Lx
        else 
            dkx = 1.5_dp * pi / Lx
        end if
        if (nky > 1) then
            dky = 3.0_dp / dble(nky) * pi / Ly
        else
            dky = 1.5_dp * pi / Ly
        end if
        if (nkz > 1) then
            dkz = 3.0_dp / dble(nkz) * pi / Lz
        else
            dkz = 1.5_dp * pi / Lz
        end if
        do i = 1,nkx
            do j = 1,nky
                do l = 1,nkz
                    kx = dble(i-1)*dkx - 1.5_dp * pi / Lx
                    ky = dble(j-1)*dky - 1.5_dp * pi / Ly
                    kz = dble(l-1)*dkz - 1.5_dp * pi / Lz
                    !
                    call mat_def_periodic(Hii, kx,ky,kz, nb,hr,cell,n_range)    
                    !
                    EN(1:NB,l + (j-1)*nkz + (i-1)*nkz*nky) = eig(NB,Hii)
                enddo
            enddo
        enddo
        end associate
        end associate
        end associate
    END SUBROUTINE plot_bands_bz


    !!! construct the fully periodic Hamiltonian matrix
    SUBROUTINE mat_def_periodic(Hii,kx,ky,kz,nb,hr,cell,n_range)
        ! kx, ky, and kz are in unit of [1/Ang] and in cartesian coordinate 
        integer, intent(in) :: nb
        REAL(8), INTENT(IN) :: kx, ky, kz, cell(3,3)
        complex(8), intent(in) :: hr(:,:,:,:,:)
        integer, intent(in) :: n_range(:) ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        COMPLEX(8), INTENT(OUT), DIMENSION(NB,NB) :: Hii
        real(8), dimension(3) :: kv, r
        real(8)::xhat(3),yhat(3),zhat(3)
        integer :: i,j,l
        Hii(:,:) = czero    
        call get_cart(cell,xhat,yhat,zhat)
        kv = kx*xhat + ky*yhat + kz*zhat
        associate(xmin => n_range(2))
        associate(xmax => n_range(3))
        associate(ymin => n_range(4))
        associate(ymax => n_range(5))
        associate(zmin => n_range(6))
        associate(zmax => n_range(7))
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        do i = xmin,xmax           
            do j = ymin,ymax
                do l = zmin,zmax
                    r = dble(i)*alpha + dble(j)*beta + dble(l)*gamm                         
                    Hii(:,:) = Hii(:,:) + Hr(:,:,i-xmin+1,j-ymin+1,l-zmin+1) &
                              * exp(-c1i * dot_product(r,kv))             
                enddo
            enddo      
        enddo           
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate   
        end associate
        end associate
        end associate   
    END SUBROUTINE mat_def_periodic
    

    SUBROUTINE calc_momentum_operator(pmn,method,nb,nx,ny,nz,hr,cell,n_range,wannier_center,rmn)
        implicit none
        character(len=*),intent(in)::method
        integer, intent(in) :: n_range(:),nb,nx,ny,nz ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        REAL(8), INTENT(IN) :: cell(3,3)
        complex(8), intent(in) :: hr(:,:,:,:,:)
        real(8),intent(in),optional::wannier_center(:,:)
        complex(8),dimension(:,:,:,:,:,:),intent(in),optional::rmn 
        integer::io,jo,ix,iy,iz,mx,my,mz,mo
        real(8)::r(3)
        complex(8)::pre_fact
        complex(8),dimension(3,NB,NB,nx,ny,nz),intent(out)::pmn ! [eV]
        pre_fact=1.0d-8 * dcmplx(0.0d0,1.0d0) *m0_charge/hbar_eV*1.0d2 * light_speed  ! multiply light-speed so c0*pmn in energy eV 
        associate(xmin => n_range(2))
        associate(xmax => n_range(3))
        associate(ymin => n_range(4))
        associate(ymax => n_range(5))
        associate(zmin => n_range(6))
        associate(zmax => n_range(7))
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
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
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate   
        end associate
        end associate
        end associate   
    END SUBROUTINE calc_momentum_operator
    

    SUBROUTINE w90_momentum_full_device(Ham,ky,kz,length,NS,n_range,nb,cell,pmn)
        implicit none
        integer, intent(in) :: n_range(:),nb ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        REAL(8), INTENT(IN) :: cell(3,3)        
        integer, intent(in) :: length
        integer, intent(in), optional :: NS
        real(8), intent(in) :: ky,kz
        complex(8), intent(out), dimension(NB*length,NB*length,3) :: Ham ! momentum matrix [eV]        
        integer :: i,j, k,v,l
        real(8), dimension(3) :: kv, r
        complex(8) :: phi
        complex(8),intent(in),dimension(:,:,:,:,:,:)::pmn ! [eV] generated by [[calc_momentum_operator]]
        real(8)::xhat(3),yhat(3),zhat(3)
        associate(xmin => n_range(2))
        associate(xmax => n_range(3))
        associate(ymin => n_range(4))
        associate(ymax => n_range(5))
        associate(zmin => n_range(6))
        associate(zmax => n_range(7))
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        call get_cart(cell,xhat,yhat,zhat)
        Ham = dcmplx(0.0d0,0.0d0)        
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
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate   
        end associate
        end associate
        end associate   
    END SUBROUTINE w90_momentum_full_device



    !!! construct the diagonal and off-diagonal blocks P(I,I), P(I+1,I)
    SUBROUTINE w90_momentum_blocks(Hii,H1i,kx,ky,kz,NS,n_range,nb,cell,pmn)
        ! ky in [2pi/Ang]
        implicit none
        integer, intent(in) :: n_range(:),nb ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
        REAL(8), INTENT(IN) :: cell(3,3)
        integer, intent(in) :: ns
        COMPLEX(8), INTENT(OUT), DIMENSION(NB*ns,NB*ns,3) :: Hii, H1i ! momentum matrix block [eV]
        real(8), intent(in) :: ky,kx,kz
        complex(8),intent(in),dimension(:,:,:,:,:,:)::pmn ! [eV] generated by [[calc_momentum_operator]]
        integer :: i,j,k,l,v
        real(8), dimension(3) :: kv, r
        complex(8) :: phi
        real(8)::xhat(3),yhat(3),zhat(3)
        associate(xmin => n_range(2))
        associate(xmax => n_range(3))
        associate(ymin => n_range(4))
        associate(ymax => n_range(5))
        associate(zmin => n_range(6))
        associate(zmax => n_range(7))
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        call get_cart(cell,xhat,yhat,zhat)
        Hii(:,:,:) = dcmplx(0.0d0,0.0d0)  
        H1i(:,:,:) = dcmplx(0.0d0,0.0d0)  
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
        end associate
        end associate
        end associate
        end associate
        end associate
        end associate   
        end associate
        end associate
        end associate   
    END SUBROUTINE w90_momentum_blocks
    
    
    FUNCTION norm(vector)
        REAL(8) :: vector(3),norm
        norm = sqrt(dot_product(vector,vector))
    END FUNCTION

    ! vector cross-product
    FUNCTION cross(a, b)
        REAL(8), DIMENSION(3) :: cross
        REAL(8), DIMENSION(3), INTENT(IN) :: a, b
        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)
    END FUNCTION cross


    subroutine get_cart(cell,xhat,yhat,zhat)
        real(8),intent(in)::cell(3,3)
        real(8),intent(out)::xhat(3),yhat(3),zhat(3)
        associate(alpha => cell(:,1))
        associate(beta => cell(:,2))
        associate(gamm => cell(:,3))
        !
        xhat = alpha/norm(alpha)
        yhat = - cross(xhat,gamm)
        yhat = yhat/norm(yhat)
        zhat = gamm/norm(gamm)   
        end associate
        end associate
        end associate
    end subroutine get_cart


END MODULE wannierHam





