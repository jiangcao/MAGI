! -*- f90 -*-

module parameters_mod
    
    implicit none 

    !constants
    integer, parameter :: dp=8
    REAL(kind=dp), PARAMETER :: pi=3.14159265359d0
    REAL(kind=dp), PARAMETER :: twopi = 3.14159265359d0*2.0d0
    REAL(kind=dp), PARAMETER :: e_charge=1.6d-19            ! charge of an electron (C)
    REAL(kind=dp), PARAMETER :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)    
    REAL(kind=dp), PARAMETER :: light_speed=2.998d8           ! m/s
    REAL(kind=dp), PARAMETER :: m0_charge=5.6856D-16        ! eV s2 / cm2
    REAL(kind=dp), PARAMETER :: hbar=1.0546d-34     ! value of hbar=h/2pi (J s)
    REAL(kind=dp), PARAMETER :: hbar_eV=hbar/e_charge ! eV s    
    COMPLEX(kind=dp), PARAMETER :: cone = dcmplx(1.0d0,0.0d0)
    COMPLEX(kind=dp), PARAMETER :: czero  = dcmplx(0.0d0,0.0d0)
    COMPLEX(kind=dp), PARAMETER :: c1i  = dcmplx(0.0d0,1.0d0)     
    REAL(kind=dp), PARAMETER  :: BOLTZ = 8.61734d-05 !eV K-1
   
end module parameters_mod

MODULE wannierHam

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


END MODULE wannierHam





