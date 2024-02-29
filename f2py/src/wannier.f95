! -*- f90 -*-

MODULE wannierHam

use parameters_mod

IMPLICIT NONE 

CONTAINS

SUBROUTINE load_from_file(fname,lreorder_axis,axis,nb,nx,ny,nz,hr,wannier_center,n_range,cell)    
    character(len=*), intent(in) :: fname
    logical, intent(in) :: lreorder_axis
    integer, intent(in) :: axis(3),nb,nx,ny,nz
    complex(8), intent(out) :: hr(nb,nb,nx,ny,nz)
    real(8), intent(out) :: wannier_center(3,nb)
    real(8), intent(out) :: cell(3,3)    
    integer, intent(out) :: n_range(9) ! [nb,xmin,xmax,ymin,ymax,zmin,zmax,nvb,nspin]
    ! ----
    integer :: fid
    integer :: n,i,rc
    character(len=40) :: comment
    REAL(8) :: aux2(3,3),alpha(3),beta(3),gamm(3)    
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
    n_range(1) = nint(maxval(ham(:,4))) ! number of WFs
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
    real(8), intent(in) :: cell(3,3)    
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





