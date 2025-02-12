module post_mod
  use parameters_mod, only: pi,dp
  use wannierHam3d,only : alpha,beta,gamm,b1,b2,Lz,Ly,xhat,yhat,zhat
  use utilities
  use linalg
  implicit none 
  private
  public:: W_to_real_space
  
  contains


  subroutine W_to_real_space(dataset, number, nky, nkz, nen, nb ,w_rspace)
      implicit none
      character(len=*), intent(in) :: dataset
      complex(dp), allocatable, intent(out) :: w_rspace(:,:,:,:)
      integer, intent(in):: number,nky,nkz,nen,nb
      integer:: N,i,j,sci,scj,iky,ikz,ib,jb
      real(dp), dimension(3) :: qv,R,minRv
      real(dp)::minR,ky,kz,dky,dkz
      real(dp), allocatable::W_real(:,:,:,:),W_imag(:,:,:,:),q(:,:,:), rij(:,:,:)
      complex(dp), allocatable :: W(:,:,:,:)
      complex(dp)::tr    

      allocate(W_real(nb,nb,nky,nkz))
      allocate(W_imag(nb,nb,nky,nkz))
      open(unit=11,file=trim(dataset)//TRIM(STRING(number))//'_grid.dat',status='unknown')
      do iky=1,nky
          do ikz=1,nkz
              do ib=1,nb
                  do jb=1,nb
                    read(11,*) W_real(ib,jb,iky,ikz), W_imag(ib,jb,iky,ikz) 
                  enddo
              enddo
          enddo
      enddo    
      close(11)
      allocate(q(3,nky,nkz))        ! k/q grid 
      allocate(rij (3,nky,nkz))     ! Rij=Ri-Rj
      allocate(w(NB,NB,nky,nkz))
      if (.not.(allocated(w_rspace))) then
          allocate(w_rspace(NB,NB,nky,nkz))
      else
          deallocate(w_rspace)
          allocate(w_rspace(NB,NB,nky,nkz))
      endif
      W = cmplx(W_real,W_imag) 
      deallocate(W_real,W_imag)
      w_rspace=cmplx(0.0d0,0.0d0)
      !
      do i=0,NKy-1
          do j=0,NKz-1                
              minRv = dble(i)*beta + dble(j)*gamm
              minR = norm(minRv)
              do sci=-1,1           ! super-cell index (loop over nearest-neighbors) 
                  do scj=-1,1       ! to find the minimum R_ij under the periodic boundary condition
                      R=dble(i+ sci*NKy )*beta + dble(j+ scj*NKz )*gamm
                      if (norm(R) .lt. minR ) then
                          minR=norm(R)
                          minRv=R                                
                      end if
                  end do
              end do
              R=minRv
              rij(:,i+1,j+1)=R                        
          end do
      end do
    
      N=nky*nkz

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

      do iky=1,nky
        ky=-pi/Ly + dble(iky)*dky
        do ikz=1,nkz
          kz=-pi/Lz + dble(ikz)*dkz
          q(:,iky,ikz)=ky*yhat + kz*zhat
        enddo
      enddo
    
      print *,'calculate ',dataset,' in real space'
      !$omp parallel default(none) private(iky,ikz,i,j)shared(N,w,rij,nky,nkz,w_rspace,q)
      !$omp do
      do iky=1,nky
        do ikz=1,nkz         
            do i=1,nky
                do j=1,nkz
                    w_rspace(:,:,i,j)=w_rspace(:,:,i,j) + exp(dcmplx(0.0d0,1.0d0)*dot_product(q(:,iky,ikz),rij(1:3,i,j)))*w(:,:,iky,ikz)/dble(N)                   
                end do
              end do        
        end do
      end do
      !$omp end do
      !$omp end parallel

      open(unit=11,file='rspace_'//trim(dataset)//TRIM(STRING(number))//'_grid.dat',status='unknown')
      do iky=1,nky
          do ikz=1,nkz            
              do ib=1,nb
                do jb=1,nb
                    tr= w_rspace(ib,jb,iky,ikz)
                    write(11,'(4E18.4)') rij(1:2,iky,ikz), real(tr),aimag(tr)
                enddo
              enddo      
          enddo
      enddo    
      close(11)


      open(unit=11,file='rspace_trace_'//trim(dataset)//TRIM(STRING(number))//'_grid.dat',status='unknown')
      do iky=1,nky
          do ikz=1,nkz 
              tr=0.0_dp           
              do ib=1,nb
                do jb=1,nb
                    tr=tr+ w_rspace(ib,jb,iky,ikz)                  
                enddo
              enddo      
              tr=tr / dble(nb)/dble(nb)
              write(11,'(4E18.4)') rij(1:2,iky,ikz), real(tr),aimag(tr)
          enddo
      enddo    
      close(11)

      deallocate(q,w,rij)

  end subroutine W_to_real_space


end module post_mod
