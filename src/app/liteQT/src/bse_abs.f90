PROGRAM BSE_ABS
USE setup_mod
use bse_mod
use sort, only: dmerge_sort
use parameters_mod
USE linalg, only : eigv, eigv_feast,cross,norm
USE utilities
use xsf_mod
USE wannierHam3d, only : NB, w90_load_from_file, w90_free_memory,Ly,w90_mat_def_kv,b1,b2,&
                         wannier_center,alpha,beta,totnvb=>nvb,spin_deg,cell,n_atom,atoms_pos_cart,name_atoms

implicit none
integer, parameter :: maxnatom = 6
integer :: nk,i,j,l,iv,ivd,sci,scj,a1,a2,a3,a4,a3p,a4p,vb0,v,vd,newsize,npt,N
integer :: ivb,icb,ikx,iky,jvb,jcb,jkx,jky,iqx,iqy,s,ib,jb,ikz,jkz,m
real(kind=dp) :: ky,kx,minR,realpart,imagpart,a
complex(kind=dp)::Kernel_d,Kernel_x,tmp1
real(kind=dp), dimension(3) :: R, minRv, rvvd, qv,dk
integer, allocatable :: indice(:,:), ind_keep(:), index2(:,:)
real(kind=dp), allocatable :: kv(:,:,:),ek(:,:,:),D(:),rij(:,:,:),rijvvd(:,:,:,:,:),wrijvvd(:,:,:,:),vrijvvd(:,:,:,:),q(:,:,:),ediff(:,:,:,:),omega(:),efieldpot(:,:,:,:),hbarw(:),absp(:),grid(:,:)
complex(kind=dp), allocatable ::kernel(:,:),Hiimat(:,:,:,:),ak(:,:,:,:),Hii(:,:),w_rspace(:,:,:,:)
complex(kind=dp), allocatable::wsum(:,:,:,:),vsum(:,:),Knew(:,:),asvckmat (:,:,:,:,:),sumas(:),chixyz(:),chi(:,:,:,:)
integer, allocatable:: asvck_reordered_index(:)
real(kind=dp):: partsum
real(kind=dp),allocatable:: asvck(:),atoms_pos_cart_shifted(:,:)
! 
real(dp):: origin(3)        ! origin
character(len=2):: nat(3)
integer:: iat
real(dp):: primvec(3,3)        ! primitive lattice vectors
real(dp):: convvec(3,3)        ! conventional lattice vectors
real(dp):: coor(3,maxnatom)    ! atomic coordinates
real(dp),allocatable:: values(:,:,:)

lread_input_bse = .true.
lread_input_post = .true.
lread_input_gw = .false.
lread_input_ham = .false.

call begin_setup()

open(unit=10,file='ham_dat',status='unknown')
call w90_load_from_file(10)
close(10)

N=Nkx*Nky
print '(3A15)','epsilon=','k0=','xi='
print '(3F15.4)', epsilon,k0,xi
print '(3A15)','E_cutoff=','nkx=','nky='
print '(F15.4,2I15)', E_cutoff,nkx,nky
print '(3A15)','E_scissor=','ncb=','nvb='
print '(F15.4,2I15)', E_scissor,ncb,nvb

allocate(Hii(nb,nb))
allocate(D(nb))
allocate(Hiimat(NKX,NKY,nb,nb))    ! store all hamiltonians
Hiimat = dcmplx(0.0d0,0.0)
allocate(ek(nb,NKX,NKY))           ! store all eigen energies
ek=0.0d0
allocate(ak(nb,nb,NKX,NKY))        ! store all eigen vectors
ak = dcmplx(0.0d0,0.0d0)
allocate(kv(3,NKX,NKY))            ! store all the k vectors
kv=0.0d0
do i = 1,NKX
    do j = 1,NKY
        if (nky==1) then 
            kv(:,i,j) = 1.0d0/dble(NKX)*dble(i)*twopi*b1      ! frac. coord. to cartesian           
        else
            if (nkx==1) then 
                kv(:,i,j) = 1.0d0/dble(NKY)*dble(j)*twopi*b2      ! frac. coord. to cartesian       
            else
                kv(:,i,j) = 1.0d0/dble(NKX)*dble(i)*twopi*b1 + 1.0d0/dble(NKY)*dble(j)*twopi*b2      ! frac. coord. to cartesian       
            endif
        endif
        call w90_mat_def_kv(Hii,kv(:,i,j))   ! calculate Hamiltonian for a single kx and ky
        Hiimat(i,j,:,:)=Hii(:,:)
        D = eigv(nb, Hii)                ! Eigen value problem  
        ek(:,i,j) = D      
        ak(:,:,i,j) = Hii(:,:);           ! obtain eigen vectors corresponding to the eigen energies
    end do
end do ! ek(band,kx,ky)

open(unit=11,file='kpoints.dat',status='unknown')
open(unit=12,file='energies.dat',status='unknown')
do i = 1,NKX
    do j = 1,NKY
        write(11,'(3F15.4)') kv(:,i,j)
        do l = 1,NB
            write(12,'(2F15.4,I6,F15.4)') kv(1:2,i,j),l,ek(l,i,j)
        end do
    end do
end do    
close(11)
close(12)

allocate(rij (3,NKX,NKY));             ! Rij=Ri-Rj
allocate(rijvvd (3,nb,nb,NKX,NKY) );   ! Rijvvd= Ri-Rj+wannier center(iv)-wanniercenter(ivd)
allocate(vrijvvd (nb,nb,NKX,NKY) );    ! bare potential 
allocate(wrijvvd (nb,nb,NKX,NKY) );    ! screened potential
allocate(efieldpot(nb,nb,NKX,NKY) );   ! e field potential


print *,'calculate the W and V matrices in real space'
do i=0,NKX-1
  do j=0,NKY-1                
    minRv = dble(i)*alpha + dble(j)*beta    
    minR = norm(minRv)
    do sci=-1,1           ! super-cell index (loop over nearest-neighbors) 
        do scj=-1,1       ! to find the minimum R_ij under the periodic boundary condition
            R=dble(i+ sci*NKX )*alpha + dble(j+ scj*NKY )*beta
            if (norm(R) < minR ) then
                minR=norm(R)
                minRv=R                                
            end if
        end do
    end do
    R=minRv
    rij(:,i+1,j+1)=R                        
    do iv=1,nb
        do ivd=1,nb                                   
            call rijvaluesfunc(rvvd,i,j,R,iv,ivd,nb,wannier_center,lwcenter)            
            rijvvd(1:3,iv,ivd,i+1,j+1)=rvvd                            
            efieldpot(iv,ivd,i+1,j+1)=efieldpotential(R,e1,F,k0,i,j,alpha,beta,NKX,NKY)  ! add an electric field potential
            wrijvvd(iv,ivd,i+1,j+1)=screenedpot(rvvd,twopi*xi,epsilon0,epsilon,e_charge)+efieldpot(iv,ivd,i+1,j+1) ! in eV
            vrijvvd(iv,ivd,i+1,j+1)=barepot(rvvd,epsilon0,epsilon,e_charge)                  ! in eV      
        end do
    end do                                    
  end do
end do

open(unit=11,file='wr.dat',status='unknown')
open(unit=12,file='vr.dat',status='unknown')
open(unit=13,file='efr.dat',status='unknown')
do i = 1,NKX
    do j = 1,NKY
      do iv=1,nb
        do ivd=1,nb
          write(11,'(2I8,2F15.4)') iv,ivd,norm(rijvvd(:,iv,ivd,i,j)),wrijvvd(iv,ivd,i,j)
          write(12,'(2I8,2F15.4)') iv,ivd,norm(rijvvd(:,iv,ivd,i,j)),vrijvvd(iv,ivd,i,j)
          write(13,'(2I8,2F15.4)') iv,ivd,norm(rijvvd(:,iv,ivd,i,j)),efieldpot(iv,ivd,i,j)
        end do
       end do
    end do
end do    
close(11)
close(12)
close(13)

if (lread_screened_coulomb) then
    print *,'read in the W in real space from WR file'
    allocate(w_rspace(nb,nb,NKX,NKY));  ! realspace screened potentials
    open(unit=11,file='rspace_'//trim(dataset)//TRIM(STRING(number))//'_grid.dat',status='unknown')
    R=0.0_dp
    do ikx=1,nkx
        do iky=1,nky
            do ib=1,nb
                do jb=1,nb
                    read(11,*) R(1:2), realpart, imagpart
                    w_rspace(ib,jb,ikx,iky) = dcmplx(realpart, imagpart) 
                    !write (101,*) ikx, iky, norm(rij(:,ikx,iky)-R)
                enddo            
            enddo
        enddo
    enddo  
    close(11)

    ! open(unit=11,file='diff_'//trim(dataset)//TRIM(STRING(number))//'_grid.dat',status='unknown')
    ! do ikx=1,nkx
    !     do iky=1,nky
    !         do ib=1,nb
    !             do jb=1,nb
    !             write(11,*) rij(1:2,ikx,iky), real(w_rspace(ib,jb,ikx,iky)-wrijvvd(ib,jb,ikx,iky)), &
    !                                          aimag(w_rspace(ib,jb,ikx,iky)-wrijvvd(ib,jb,ikx,iky))
    !             enddo            
    !         enddo
    !     enddo
    ! enddo  
    ! close(11)
    do i = 1,NKX
        do j = 1,NKY
          do iv=1,nb
            do ivd=1,nb
              wrijvvd(iv,ivd,i,j) = real(w_rspace(iv,ivd,i,j))+efieldpot(iv,ivd,i,j) 
            end do
           end do
        end do
    end do 
endif

if (lbse) then
    print *,'calculate fourier transform of W and V, wsum and vsum'

    allocate(wsum(nb,nb,NKX*2,NKY*2));  ! sum of screened potentials
    allocate(vsum(nb,nb));              ! sum of bare potentials
    allocate(q(3,NKX*2,NKY*2));         ! grid having k-k'
    do ikx=1-NKX,NKX
        do iky=1-NKY,NKY
            qv=1.0d0/dble(NKX)*dble(ikx)*twopi*b1 + 1.0d0/dble(NKY)*dble(iky)*twopi*b2; ! frac. coord. to cartesian       
            q(:,ikx+NKX,iky+NKY)=qv;
        end do
    end do
    wsum = dcmplx(0.0d0,0.0d0)
    vsum = dcmplx(0.0d0,0.0d0)
    !$omp parallel default(none) private(ikx,iky,i,j)shared(N,wrijvvd,rij,q,wsum,nkx,nky)
    !$omp do
    do ikx=1,NKX*2
        do iky=1,NKY*2         
            do i=1,NKX
                do j=1,NKY
                    wsum(:,:,ikx,iky)=wsum(:,:,ikx,iky) + exp(-dcmplx(0.0d0,1.0d0)*dot_product(q(:,ikx,iky),rij(1:3,i,j)))*wrijvvd(:,:,i,j)/dble(N)   ! part of eqn 5                 
            end do
            end do        
        end do
    end do
    !$omp end do
    !$omp end parallel

    do i=1,NKX
        do j=1,NKY
            vsum=vsum+vrijvvd(:,:,i,j)/dble(N)   ! part of eqn 6
        end do
    end do

    open(unit=11,file='wq.dat',status='unknown')
    do i = 1,NKX*2
        do j = 1,NKY*2
            write(11,'(4F15.4)') q(1,i,j),q(2,i,j),real(wsum(2,2,i,j)),imag(wsum(2,2,i,j))
        end do
    end do    
    close(11)
    
end if

print *,'build the BSE'

a=norm(cross(alpha,beta))          ! area of unit cell in Ang^2

vb0=totnvb-nvb
allocate(ediff (nvb,ncb,NKX,NKY) )
allocate(index2(4,nvb*ncb*NKX*NKY))
newsize=0
do a1=1,nvb
    do a2=1,ncb
        do a3=1,NKX
            do a4=1,NKY                
                a3p=a3+exciton_q(1)
                if (a3p>NKX) a3p = a3p - NKX
                if (a3p<1) a3p = a3p + NKX
                a4p=a4+exciton_q(2)
                if (a4p>NKX) a4p = a4p - NKY
                if (a4p<1) a4p = a4p + NKY
                ediff(a1,a2,a3,a4) = ek(a2+totnvb,a3p,a4p) - ek(vb0+a1,a3,a4)   ! ek=ecb-evb                            
                if (ediff(a1,a2,a3,a4) < E_cutoff) then
                    newsize=newsize+1
                    index2(:,newsize) = (/a1,a2,a3,a4/)
                endif
            end do
        end do
    end do
end do
print *, 'nnz=', newsize
allocate(Knew(newsize,newsize) )
knew = dcmplx(0.0d0,0.0d0)

if (lbse) then       
    print *, '   calculate Kernel'
    !$omp parallel default(shared) private(i,ivb,icb,ikx,iky,j,jvb,jcb,jkx,jky,dk,iqx,iqy,Kernel_d,Kernel_x,v,vd) 
    !$omp do
    do i = 1,newsize
        
        ivb=index2(1,i)
        icb=index2(2,i)
        ikx=index2(3,i)
        iky=index2(4,i)
        
        do j = 1,newsize        
            
            jvb=index2(1,j)
            jcb=index2(2,j)
            jkx=index2(3,j)
            jky=index2(4,j)
                    
            dk=kv(:,ikx,iky)-kv(:,jkx,jky)
            iqx = (ikx - jkx)
            iqy = (iky - jky)

            iqx = iqx+NKX  ! Obtain the indices iqx and iqy where q=k-k'
            iqy = iqy+NKY
            if (norm(dk(:) - q(:,iqx,iqy)) .gt. 1e-10)  then ! terminate if q#k-k'
                print *, "k-k' ~= q !!"
                print *,dk
                print *,q(:,iqx,iqy)
                stop
            end if
            !Implementations of equations 5 and 6 from the paper PRB 94, 245434
            Kernel_d=dcmplx(0.0d0,0.0d0) ! direct term
            Kernel_x=dcmplx(0.0d0,0.0d0) ! exchange term        
            do v=1,nb
                do vd=1,nb
                    if (ldirect) then
                        Kernel_d = Kernel_d + (conjg(ak(v,icb+totnvb,ikx,iky))*ak(v,jcb+totnvb,jkx,jky))*wsum(v,vd,iqx,iqy)*(ak(vd,ivb+vb0,ikx,iky)*conjg(ak(vd,jvb+vb0,jkx,jky))) ! use the indices iqx and iqy to get correct wsum values
                    endif
                    
                !    print *, (conjg(ak(v,icb+totnvb,ikx,iky))*ak(v,jcb+totnvb,jkx,jky))*wsum(v,vd,iqx,iqy)*(ak(vd,ivb+vb0,ikx,iky)*conjg(ak(vd,jvb+vb0,jkx,jky)))
                    if (lexchange) then
                        Kernel_x = Kernel_x + (conjg(ak(v,icb+totnvb,ikx,iky))*ak(v,ivb+vb0,ikx,iky))*vsum(v,vd)*(ak(vd,jcb+totnvb,jkx,jky)*conjg(ak(vd,jvb+vb0,jkx,jky)))
                    endif
                end do
            end do
                            
            Knew(i,j) = -Kernel_d + Kernel_x ! obtain the kernel for the BSE eigen value problem   

        end do
        Knew(i,i) = Knew(i,i) + dcmplx(ediff(ivb,icb,ikx,iky) + E_scissor,0.0d0)   ! add the ek value to the diagnol elements of the kernel            
    end do
    !$omp end do
    !$omp end parallel
else
    do i = 1,newsize
        
        ivb=index2(1,i)
        icb=index2(2,i)
        ikx=index2(3,i)
        iky=index2(4,i)
                
        Knew(i,i) = dcmplx( ediff(ivb,icb,ikx,iky) + E_scissor,0.0d0 )   ! add the ek value to the diagnol elements of the kernel            
    end do
end if

if (any(abs(Knew - conjg(transpose(Knew))) .gt. 1e-7)) then
    print *,'kernel not Hermitian!!'
    print *, maxval(abs(Knew - conjg(transpose(Knew))))
    stop
end if

if (lbse) then
deallocate(wsum)
deallocate(vsum)
deallocate(q)
end if

print *,'solve the BSE'

allocate(omega(newsize))

call mkl_set_num_threads(ncpu)

! if (newsize < 20000) then 
    omega=eigv(newsize,Knew)  ! solve the bse eigen value problem
    m=newsize ! m is number of eigenvalues
! else
    ! omega=eigv_feast(newsize, Knew, 0.0_dp, E_cutoff * 0.6, m, newsize/2)    ! m is number of eigenvalues
! endif

allocate( asvckmat (nvb,ncb,NKX,NKY,m) )
asvckmat=czero
do s=1,m
    do i=1,newsize

        ivb=index2(1,i)
        icb=index2(2,i)
        ikx=index2(3,i)
        iky=index2(4,i)

        asvckmat(ivb,icb,ikx,iky,s) = Knew(i,s)

    end do
end do

open(unit=11,file='omega.dat',status='unknown')
do s=1,m
    write(11,'(F15.4)') omega(s)    
end do    
close(11)

open(unit=11,file='asvck.dat',status='unknown')
write(11,*) '# exciton state, ivb, icb, kvec(ikx,iky), re(asvck), im(asvck)' 
 do s=1,nex
    do ivb=1,nvb
        do icb=1,ncb
            do ikx=1,nkx
                do iky=1,nky
                    if (abs(asvckmat(ivb,icb,ikx,iky,s)) > 0.0_dp) then
                        write(11,'(3I10, 4F15.4)') s, ivb, icb, kv(1:2,ikx, iky), dble(asvckmat(ivb,icb,ikx,iky,s)),aimag(asvckmat(ivb,icb,ikx,iky,s))
                    endif
                enddo
            enddo
        enddo
    enddo
 enddo
close(11)


!sort and reorder asvckmat entries
allocate( asvck(newsize) )
allocate( asvck_reordered_index(newsize) )

open(unit=11,file='asvck_reordered.dat',status='unknown')
write(11,*) '# exciton state , ivb, icb, kvec(ikx,iky), amplitude' 
do s=1,nex
    asvck = abs(Knew(:,s)) 
    !insert efficient sorting algorithm here
    call dmerge_sort(asvck, asvck_reordered_index)    
    !
    partsum = 0.0_dp
    i=newsize+1
    do while ( partsum < 0.9 )  ! 
        i = i-1
        j = asvck_reordered_index(i)
        partsum = partsum + asvck(j)**2
        !
        ivb=index2(1,j)
        icb=index2(2,j)
        ikx=index2(3,j)
        iky=index2(4,j)
        write(11,'(3I10, 3F15.4)') s, ivb, icb, kv(1:2, ikx,iky), asvck(j)
    enddo    
enddo
close(11)

print *, 'compute the absorbance spectrum'

allocate(sumas(newsize))
sumas(:)=0.0d0
allocate(hbarw(nen))
hbarw=(/(i, i=1,nen, 1)/) / dble(nen) * E_cutoff * 0.6 ! photon energies in eV
allocate(absp(nen))

do ivb=1,nvb
    do icb=1,ncb
        do ikx=1,NKX
            do iky=1,NKY                                                                          
                tmp1 = dHdk(nb,nkx,nky,ivb+vb0,icb,ikx,iky,ak,Hiimat,totnvb,kv)
                do s=1,m
                    sumas(s)=sumas(s) + asvckmat(ivb,icb,ikx,iky,s) * tmp1
                    ! multiplication of asvck and dvck summed over v,c,k                    
                end do                
            end do
        end do
    end do
end do

sumas = abs(sumas)**2
absp=0.0d0

do s=1,m    
  do i=1,nen
    ! summed over exciton states with the delta function    
    absp(i)=absp(i)+(4.0d0*pi*pi)*7.297d-3*spin_deg/A/dble(N)*dble(sumas(s))/hbarw(i) * gaussian(hbarw(i),omega(s),sig) 
    ! eps_0 = e^2/(2 alpha h c) with alpha the fine-structure constant = 7.297*10^-3. Implementation of eqn 8
  end do
end do

open(unit=11,file='absp.dat',status='unknown')
do i = 1,size(hbarw)   
    write(11,'(2F15.4)') hbarw(i), absp(i)    
end do    
close(11)
open(unit=11,file='sumas.dat',status='unknown')
do s=1,newsize
    write(11,'(2F15.4)') omega(s) ,   dble(sumas(s))
end do    
close(11)

if (lplotexciton) then
    
    npt=nx*ny*nz
    allocate(grid(3,npt))
    allocate(chi(nb,nb,nkx,nky))
    allocate(chixyz(size(grid,2)))
    s=0
    do i=1,nz
        do j=1,ny
            do l=1,nx
                s=s+1
                grid(:,s) = (/ dble(l-nx/2)*dx, dble(j-ny/2)*dy, dble(i-nz/2)*dz /) 
            enddo
        enddo
    enddo

    allocate(values(nx,ny,nz))
    allocate(atoms_pos_cart_shifted(3,n_atom))
    do i=1,n_atom
        atoms_pos_cart_shifted(:,i) = atoms_pos_cart(:,i) !- sum(atoms_pos_cart,2)/dble(n_atom) !+ sum(cell,2)/2.0_dp 
    enddo    

    do s=1,nex ! exciton state we want to plot
    
        print *, 'compute the exciton wavefunction'
        call exciton_wavefunction_simple(asvckmat,s,NKX,NKY,ncb,nvb,newsize,nb,totnvb,ak,chi,kv,rij)
        print *, 'compute the exciton wavefunction on the grid'    
        call exciton_wavefunction_grid(NKX,NKY,nb,chi,chixyz,rij,grid,npt,rsmear,wannier_center)                        

        open(unit=11,file='chi_3d_'//string(s)//'.xsf',status='unknown')

        primvec = cell 
        
        call xsf_write_header (primvec,primvec,atoms_pos_cart,n_atom, name_atoms, 11)
 
        origin = -(/dble(nx)*dx, dble(ny)*dy, dble(nz)*dz/) / 2.0_dp 
        
        values=reshape(abs(chixyz), (/nx,ny,nz/) )

        primvec = 0.0_dp
        primvec(1,1) = dble(nx)*dx
        primvec(2,2) = dble(ny)*dy
        primvec(3,3) = dble(nz)*dz

        if (nky==1) then 
            primvec(:,2)=primvec(:,2)/2.0_dp
            origin(2)=origin(2)/2.0_dp
        endif

        origin = origin + sum(atoms_pos_cart,2)/dble(n_atom) 
        
        call xsf_write_3ddatablock(values, nx, ny, nz, origin, primvec, 11)
        close(11)

        open(unit=11,file='chi_xy_'//string(s)//'.dat',status='unknown')
        do i=1,nx
            do j=1,ny
                write(11,'(3E25.8)') grid(1:2,i+(j-1)*nx), sum(values(i,j,:)**2)    
            end do
        end do
        close(11)
        
        open(unit=11,file='chi_xz_'//string(s)//'.dat',status='unknown')
        do i=1,nx
            do j=1,nz
                write(11,'(3E25.8)') grid(1,i+(j-1)*ny*nx),grid(3,i+(j-1)*ny*nx), sum(values(i,:,j)**2)    
            end do
        end do
        close(11)
    enddo

    deallocate(values,atoms_pos_cart_shifted)
    deallocate(grid)
    deallocate(chi)
    deallocate(chixyz)
end if

print *, 'Free memory'
deallocate(Hii)
deallocate(D)
deallocate(Hiimat) 
deallocate(ek)
deallocate(ak)
deallocate(kv)
deallocate(rij)
deallocate(rijvvd)
deallocate(wrijvvd)
deallocate(vrijvvd)
deallocate(efieldpot)
deallocate(ediff)
deallocate(index2)
deallocate(Knew)
deallocate(omega)
deallocate(asvckmat)
deallocate(asvck)
deallocate(asvck_reordered_index)
deallocate(sumas)
deallocate(hbarw)
deallocate(absp)
if (lread_screened_coulomb) then
deallocate(w_rspace)

endif

call w90_free_memory()

call MPI_Finalize( ierr )

END PROGRAM BSE_ABS


