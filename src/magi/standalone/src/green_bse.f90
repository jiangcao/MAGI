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
module green_bse
use linalg
use parameters_mod,only:dp,twopi,pi,e_charge,c0,epsilon0,hbar,c1i,czero,cone

implicit none 

private

public :: green_bse_solve, green_bse_fullsolve, green_bse_fullsolve_opt

include "mpif.h"

CONTAINS

    subroutine four_polarization(alpha,nm_dev,nen,en,nop,ndiag,G_lesser,G_greater,G_retarded,i,j,k,l,L0)
       integer,intent(in) :: nm_dev,nen,nop,ndiag, i,j,k,l
       real(dp),intent(in) :: en(nen), alpha 
       complex(dp),intent(in),dimension(nm_dev,nm_dev,nen) :: G_lesser,G_greater,G_retarded
       complex(dp),intent(out) :: L0
       ! ---
       real(dp) :: dE, weights, xen
       integer :: ie, isub, ik, ikd
       ! the P4 IPA tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
       dE = ( En(2) - En(1) )          
       weights=dE/twopi
       !                 
       ! calculate P4_IPA from GG
       L0 =    (1.0_dp - alpha) * ( sum( G_lesser(j,l,(nop+1):nen) * conjg(G_retarded(i,k,1:(nen-nop))) ) + &
                                    sum( G_retarded(j,l,(nop+1):nen) * G_lesser(k,i,1:(nen-nop)) ) )   + &
                 alpha * 0.5_dp * ( sum( G_greater(j,l,(nop+1):nen) * G_lesser(k,i,1:(nen-nop)) ) - & 
                                    sum( G_lesser(j,l,(nop+1):nen)  * G_greater(k,i,1:(nen-nop)) ) )  
       L0 = L0 * weights 
    end subroutine four_polarization


    subroutine four_polarization_old(alpha,nm_dev,nen,en,nop,ndiag,G_lesser,G_greater,G_retarded,i,j,k,l,L0)
       integer,intent(in) :: nm_dev,nen,nop,ndiag, i,j,k,l
       real(dp),intent(in) :: en(nen), alpha 
       complex(dp),intent(in),dimension(nm_dev,nm_dev,nen) :: G_lesser,G_greater,G_retarded
       complex(dp),intent(out) :: L0
       ! ---
       real(dp) :: dE, weights, xen
       integer :: ie, isub, ik, ikd
       ! the P4 IPA tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
       dE = ( En(2) - En(1) )          
       weights=dE/twopi
       !                 
       ! calculate P4_IPA from GG
       L0=czero                   
       ! calculate P4_IPA from GG       
       do ie=nop+1,nen                            
           L0 = L0 + &
                   (1.0_dp - alpha) * ( G_lesser(j,l,ie) * conjg(G_retarded(i,k,ie-nop)) + &
                                       G_retarded(j,l,ie) * G_lesser(k,i,ie-nop) )   + &
                       alpha * 0.5_dp * ( G_greater(j,l,ie) * G_lesser(k,i,ie-nop) - & 
                                       G_lesser(j,l,ie)  * G_greater(k,i,ie-nop) )  
       enddo 
       L0 = L0 * weights 
    end subroutine four_polarization_old


! solve the full Bethe-Salpeter Equation optimized
    subroutine green_bse_fullsolve_opt(alpha,spindeg,nm_dev,ndiag,nen,En,nop,G_lesser,G_greater,G_retarded,W,V,solve,P_retarded,epsilon_M)        
        integer,intent(in)::nm_dev,nen,nop,ndiag
        real(dp),intent(in)::en(nen),spindeg,alpha
        logical,intent(in),optional::solve
        !integer, intent(out)::nn ! size of the system
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen):: G_lesser,G_greater,G_retarded ! electron GFs
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: P_retarded ! 2-point polarization function with interacting electron-hole at frequency [[nop]]        
        
        !complex(dp),intent(out),dimension(nm_dev*(ndiag*2+1) ,nm_dev*(ndiag*2+1)):: system ! system matrix
        !complex(dp),intent(out),dimension(nm_dev*(ndiag*2+1) ,nm_dev*(ndiag*2+1)),optional:: L0,M ! L0 matrix
        complex(dp),intent(out),dimension(nm_dev,nm_dev) :: epsilon_M ! macroscopic dielectric function
        !complex(dp),intent(out),dimension(nm_dev,nm_dev,nm_dev,nm_dev),optional:: P4_retarded ! 4-point polarization function with interacting electron-hole 
        !---------
        complex(dp),dimension(:,:),allocatable :: Lmat ! two-particle Green's function 
        complex(dp),dimension(:,:),allocatable :: Mmat ! 4-point Kernel
        complex(dp),dimension(:,:),allocatable :: Amat ! system matrix        
        complex(dp) :: epsM, L0ijkl
        logical::lsolve
        real(dp) :: start, finish
        integer :: N,i,j,k,l,p,q,ie,row,col, it, ii,jj,nn
        integer,allocatable::table(:,:)
        !          
        lsolve=.true.
        if(present(solve)) then 
            lsolve = solve
        endif
        N = nm_dev*nm_dev
        allocate(table(2,N))
        ! construct the table of reordered indices   
        it=0
        ! first put the i=j        
        do i=1,nm_dev
            it=it+1
            table(:,it) = [i,i]          
        enddo
        ! then put the others, but first within the ndiag
        do i=1,nm_dev
            do j=1,nm_dev               
                if (i/=j) then                     
                    if (abs(i-j)<=ndiag) then
                        it=it+1
                        table(:,it) = [i,j]                    
                    endif
                endif                    
            enddo
        enddo
        nn=it
        ! then put the others, but outside ndiag
        do i=1,nm_dev
            do j=1,nm_dev
                if (i/=j) then 
                    if (abs(i-j)>ndiag) then
                        it=it+1
                        table(:,it) = [i,j]
                    endif
                endif                    
            enddo
        enddo
        if (it/=N) then 
            print *, 'ERROR!'
            call abort
        endif
        N = nn ! resize the problem
        print *, 'nm_dev=',nm_dev
        print *, 'resized system size=',N 
        ! start computation
        allocate(Mmat(N,N), source=czero)        
        allocate(Lmat(N,N), source=czero)     
        allocate(Amat(N,N), source=czero)                
        !
        start = MPI_Wtime()
        print *,'  start computation L0_ijkl = G_jl G_ki ...'        
        !$omp parallel default(shared) private(row,col,L0ijkl,i,j,k,l)
        !$omp do
        do row=1,N 
            do col=1,N
                i=table(1,row)
                j=table(2,row)
                k=table(1,col)
                l=table(2,col)
                if ((abs(i-k)<=ndiag).and.(abs(j-l)<=ndiag).and.(abs(j-k)<=ndiag).and.&
                    (abs(i-l)<=ndiag).and.(abs(i-j)<=ndiag).and.(abs(k-l)<=ndiag).and.&
                    (j>0).and.(j<=nm_dev).and.(l>0).and.(l<=nm_dev)) then 
                    call four_polarization(alpha,nm_dev,nen,en,nop,ndiag,&
                        G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                    Lmat(row,col) = L0ijkl * spindeg
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel 
        !
        !$omp parallel default(shared) private(row,col,i,j,k,l)
        !$omp do        
        do row=1,N                        
            do col=1,N
                i=table(1,row)
                j=table(2,row)
                k=table(1,col)
                l=table(2,col)           
                if ((i==j).and.(k==l)) then                        
                    Mmat(row,col) = Mmat(row,col) - c1i *  V(i,k) * spindeg                        
                endif 
                if ((i==k).and.(j==l).and.(j>0).and.(j<=nm_dev).and.(l>0).and.(l<=nm_dev)) then                        
                    Mmat(row,col) = Mmat(row,col) + c1i *  W(i,j)
                endif 
            enddo
        enddo    
        !$omp end do
        !$omp end parallel 
        !            
        finish = MPI_Wtime()
        print '("  computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
        print *,'  start computation -L0 K'
        !  
        call zgemm('n','n',N,N,N,-cone,Lmat,N,Mmat,N,czero,Amat,N)         
        !
        ! (I - L0 K) -> A
        !$omp parallel default(shared) private(i)
        !$omp do  
        do i=1,N 
            Amat(i,i) = Amat(i,i) + dcmplx(1.0_dp, 0.0_dp)
        enddo  
        !$omp end do
        !$omp end parallel 
        !        
        finish = MPI_Wtime()
        print '("  computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
        !        
        if (lsolve) then             
            N=nn
            print *,'  start invert (I - L0 K)'
            !
            call invert(Amat(1:N,1:N),N)
            !
            finish = MPI_Wtime()
            print '("  computation time = ", F0.3 ," seconds.")', finish-start
            start = finish
            print *,'  start computation L = (I - L0 K) \ L0  '
            !
            call zgemm('n','n',N,N,N,cone,Amat(1:N,1:N),N,Lmat(1:N,1:N),N,czero,Mmat(1:N,1:N),N)                 
            !
            finish = MPI_Wtime()
            print '("  computation time = ", F0.3 ," seconds.")', finish-start
            start = finish
            !
            !$omp parallel default(shared) private(row,col,i,k)
            !$omp do
            do row=1,nm_dev
                do col=1,nm_dev
                    i=table(1,row)
                    k=table(1,col)
                    P_retarded(i,k) =  - c1i * Mmat(row,col)                
                enddo
            enddo                          
            !$omp end do
            !$omp end parallel      
            !        
           ! ! ! calculate RPA epsilon and output to file
           ! call zgemm('n','n',nm_dev,nm_dev,nm_dev,c1i,V,nm_dev,Lmat(1:nm_dev,1:nm_dev),nm_dev,czero,epsilon_M,nm_dev) 
           ! do j=1,nm_dev 
           !     epsilon_M(j,j) = epsilon_M(j,j) + cone
           ! enddo      
           ! call invert(epsilon_M,nm_dev)        
           ! open(unit=99,file='rpa_epsilonM.dat',status='unknown', position="append", action="write")    
           ! epsM = sum( epsilon_M(nm_dev/2,1:nm_dev) )
           ! write(99,*) dble(nop)*(En(2)-En(1)) , - aimag(epsM), dble(epsM) ! - Im \epsilon^{-1}
           ! close(99)
           ! !
           ! ! ! calculate BSE epsilon_M and output to file        
           ! call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,P_retarded,nm_dev,czero,epsilon_M,nm_dev) 
           ! do j=1,nm_dev 
           !     epsilon_M(j,j) = epsilon_M(j,j) + cone
           ! enddo      
           ! open(unit=99,file='bse_epsilonM.dat',status='unknown', position="append", action="write")    
           ! epsM = sum( epsilon_M(nm_dev/2, :) )
           ! write(99,*) dble(nop)*(En(2)-En(1)) , aimag(epsM), dble(epsM) ! Im \epsilon_M -> absorption
           ! close(99)
        endif
        !                
        deallocate(Mmat,Lmat,Amat)
    end subroutine green_bse_fullsolve_opt
  

! solve the full Bethe-Salpeter Equation
subroutine green_bse_fullsolve(spindeg,nm_dev,ndiag,nen,En,nop,G_lesser,G_greater,G_retarded,W_retarded,V,P_retarded,luse_pr,P4_retarded)
  integer,intent(in)::nm_dev,nen,nop,ndiag
  real(dp),intent(in)::en(nen),spindeg
  complex(dp),intent(in),dimension(nm_dev,nm_dev,nen),optional:: G_lesser,G_greater,G_retarded ! electron GFs
  complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W_retarded ! W_0 static screened Coulomb interaction
  complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction
  complex(dp),intent(inout),dimension(nm_dev,nm_dev):: P_retarded ! 2-point polarization function with interacting electron-hole at frequency [[nop]]
  logical,intent(in)::luse_pr ! whether to use the P_retarded as input for P 
  complex(dp),intent(out),dimension(nm_dev,nm_dev,nm_dev,nm_dev),optional:: P4_retarded ! 4-point polarization function with interacting electron-hole 
  !---------
  complex(dp),dimension(:,:),allocatable :: Lmat ! two-particle Green's function 
  complex(dp),dimension(:,:),allocatable :: Mmat ! 4-point Kernel
  complex(dp),dimension(:,:),allocatable :: Amat ! 
  complex(dp) :: dE
  real(dp) :: start, finish, alpha
  integer :: N,i,j,k,l,p,q,ie,row,col, ne_margin
  logical :: lexchange
  !  
  alpha = 0.99_dp
  ne_margin = nen/20 ! margin of energy window
  N = nm_dev*nm_dev
  dE = ( En(2) - En(1) ) / twopi * spindeg 
  !
  allocate(Lmat(N,N))
  allocate(Mmat(N,N))
  allocate(Amat(N,N))
  print *,'  start computation L0_ijkl = G_jl G_ki ...'
  start = MPI_Wtime()
  Lmat=czero      
  !
  !$omp parallel default(shared) private(i,j,k,l,row,col,ie,lexchange)
  !$omp do
  do i=1,nm_dev
    do j=max(1,i-ndiag),min(nm_dev,i+ndiag)
      do k=max(1,i-ndiag),min(nm_dev,i+ndiag)
        do l=max(1,i-ndiag),min(nm_dev,i+ndiag)           
          if ((abs(j-l)<=ndiag).and.(abs(j-k)<=ndiag).and.(abs(k-l)<=ndiag)) then
            row= (i-1)*nm_dev + j                
            col= (k-1)*nm_dev + l
            lexchange = ((i==j) .and. (k==l))
            if ((present(G_lesser) .and. present(G_retarded)) .and. ((.not. luse_pr) .or. (.not. lexchange))) then
              ! calculate P4_IPA from -iGG
              do ie=nop+1,nen
                Lmat(row,col) = Lmat(row,col) + &
                        (1.0_dp - alpha) * ( G_lesser(j,l,ie) * conjg(G_retarded(i,k,ie-nop)) + &
                                             G_retarded(j,l,ie) * G_lesser(k,i,ie-nop) ) + &
                          alpha * 0.5_dp * ( G_greater(j,l,ie) * G_lesser(k,i,ie-nop) - & 
                                             G_lesser(j,l,ie)  * G_greater(k,i,ie-nop) ) 
              enddo          
              Lmat(row,col) = Lmat(row,col) * dE 
            else
              ! take P2_IPA entered through P_retarded 
              if (luse_pr .and. lexchange) then
                Lmat(row,col) = P_retarded(j,l)
              endif
            endif            
          endif
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  !
  finish = MPI_Wtime()
  print '("  computation time = ", F0.3 ," seconds.")', finish-start
  start = finish
  !
  Mmat=czero  
  print *,'  start computation -L0 K'
  !
  !$omp parallel default(shared) private(i,j,row,col)
  !$omp do
  do i=1,nm_dev
    do j=1,nm_dev
      ! exchange part
      row= (i-1)*nm_dev + i                
      col= (j-1)*nm_dev + j
      Mmat(row,col) = Mmat(row,col) - c1i *  V(i,j) * spindeg
      ! direct part
      row= (i-1)*nm_dev + j
      col= row 
      Mmat(row,col) = Mmat(row,col) + c1i *  W_retarded(i,j) 
    enddo
  enddo    
  !$omp end do
  !$omp end parallel  
  !call save_matrix('bse_M.dat',N, Mmat)
  !call save_matrix('bse_L0.dat',N, Lmat)
  !  
  call zgemm('n','n',N,N,N,-cone,Lmat,N,Mmat,N,czero,Amat,N) 
  !
  finish = MPI_Wtime()
  print '("  computation time = ", F0.3 ," seconds.")', finish-start
  start = finish
  ! (I - L0 K) -> A
  do i=1,N 
    Amat(i,i) = Amat(i,i) + dcmplx(1.0_dp, 0.0_dp)
  enddo  
  print *,'  start invert (I - L0 K)'
  !
  call invert(Amat,N)
  !
  finish = MPI_Wtime()
  print '("  computation time = ", F0.3 ," seconds.")', finish-start
  start = finish
  !
  print *,'  start computation L = (I - L0 K) \ L0  '
  !
  call zgemm('n','n',N,N,N,cone,Amat,N,Lmat,N,czero,Mmat,N) 
  !
  !call save_matrix('bse_L.dat',N, Mmat)
  !
  P_retarded = czero
  !$omp parallel default(shared) private(i,j,row,col)
  !$omp do
  do i=1,nm_dev
    do j=1,nm_dev      
      row= (i-1)*nm_dev + i                
      col= (j-1)*nm_dev + j
      P_retarded(i,j) = - c1i * Mmat(row,col)                
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  if (present(P4_retarded)) then 
    !$omp parallel default(shared) private(i,j,k,l,row,col)
    !$omp do
    do i=1,nm_dev
      do j=1,nm_dev      
        do k=1,nm_dev
          do l=1,nm_dev
            row= (i-1)*nm_dev + j                
            col= (k-1)*nm_dev + l
            P4_retarded(i,j,k,l) = - c1i * Mmat(row,col) 
          enddo
        enddo
      enddo
    enddo
    !$omp end do
    !$omp end parallel
  endif 
  finish = MPI_Wtime()
  print '("  computation time = ", F0.3 ," seconds.")', finish-start
  !
  deallocate(Lmat,Mmat,Amat)
end subroutine green_bse_fullsolve


! solve the Bethe-Salpeter Equation under approximation
subroutine green_bse_solve(spindeg,nm_dev,nen,En,nop,G_lesser,G_greater,G_retarded,W_retarded,V,P_retarded)
integer,intent(in)::nm_dev,nen,nop
real(dp),intent(in)::en(nen),spindeg
complex(dp),intent(in),dimension(nm_dev,nm_dev,nen):: G_lesser,G_greater,G_retarded ! electron GF
complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W_retarded ! W_0 static screened Coulomb interaction
complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction
complex(dp),intent(out),dimension(nm_dev,nm_dev):: P_retarded ! 2-point polarization function with interacting electron-hole
!---------
integer :: ie,i,j,k,l,m,n,p
integer :: nn,nm,pp,pq
complex(dp),allocatable :: L0xx(:,:),A(:,:),B(:,:),Kxx(:,:), Mxx(:,:), Sxx(:,:)
complex(dp) :: Qijkl, L0Kdd
real(dp) :: dE ,alpha 
  allocate( L0xx(nm_dev,nm_dev) )
  allocate( Mxx(nm_dev,nm_dev) )
  allocate( Sxx(nm_dev,nm_dev) )
  allocate( Kxx(nm_dev,nm_dev) )
  allocate( A(nm_dev,nm_dev) )
  L0xx = czero
  alpha = 0.99_dp
  !$omp parallel default(shared) private(i,j,ie)
  !$omp do
  do i=1,nm_dev
    do j=1,nm_dev
      do ie=nop+1,nen
        L0xx(i,j) = L0xx(i,j) + &
                    (1.0_dp - alpha) * ( G_lesser(j,i,ie) * conjg(G_retarded(i,j,ie-nop)) + &
                                           G_retarded(j,i,ie) * G_lesser(j,i,ie-nop) ) + &
                        alpha * 0.5_dp * ( G_greater(j,i,ie) * G_lesser(j,i,ie-nop) - & 
                                           G_lesser(j,i,ie)  * G_greater(j,i,ie-nop) ) 
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  dE = (en(2) - en(1))
  L0xx = L0xx * dE * spindeg
  Kxx(:,:) = - c1i*V(:,:)
  do i=1,nm_dev
    Kxx(i,i) = Kxx(i,i) + c1i*W_retarded(i,i)
  enddo
  !
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,L0xx,nm_dev,Kxx,nm_dev,czero,A,nm_dev) 
  do i=1,nm_dev
    A(i,i) = A(i,i) + cone
  enddo
  !
  Sxx = czero
  Mxx = czero  
  !$omp parallel default(shared) private(i,j,k,l,Qijkl,L0Kdd,p,ie)
  !$omp do
  do i=1,nm_dev
    do j=1,nm_dev
      do k=1,nm_dev
        do l=1,nm_dev
          if (k/=l) then
            Qijkl=czero
            do ie=nop+1,nen
                ! L0xd * L0dx
              Qijkl = Qijkl + (1.0_dp - alpha) * ( G_lesser(i,l,ie) * conjg(G_retarded(i,k,ie-nop))*G_lesser(l,j,ie) * conjg(G_retarded(k,j,ie-nop)) + &
                                                   G_retarded(i,l,ie) * G_lesser(k,i,ie-nop) * G_retarded(l,j,ie) * G_lesser(j,k,ie-nop) ) + &
                      alpha * 0.5d0*( G_greater(i,l,ie) * G_lesser(k,i,ie-nop) * G_greater(l,j,ie) * G_lesser(j,k,ie-nop) &
                                    - G_lesser(i,l,ie) * G_greater(k,i,ie-nop) * G_lesser(l,j,ie) * G_greater(j,k,ie-nop) )
                                                                      
            enddo
            Qijkl = - c1i * Qijkl * dE
            !
            Qijkl = Qijkl * W_retarded(k,l) 
            L0Kdd=czero
            do ie=nop+1,nen
              L0Kdd = L0Kdd + &
                    (1.0_dp - alpha) * ( G_lesser(k,k,ie) * conjg(G_retarded(l,l,ie-nop)) + &
                                           G_retarded(k,k,ie) * G_lesser(l,l,ie-nop) ) + &
                          alpha *0.5d0*( G_greater(k,k,ie) * G_lesser(l,l,ie-nop)   &
                                      - G_lesser(k,k,ie) * G_greater(l,l,ie-nop)   ) 
                                                                      
            enddo
            L0Kdd = L0Kdd * dE
            !
            L0Kdd = cone - c1i * L0Kdd * W_retarded(k,l) 
            Qijkl = Qijkl / L0Kdd
            Sxx(i,j) = Sxx(i,j) + Qijkl
            do p=1,nm_dev
              Mxx(i,p) = Mxx(i,p) - c1i * Qijkl * ( W_retarded(j,j) - V(j,p) ) 
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel
  !
  ! call save_matrix('Sxx.dat',nm_dev, Sxx)
  ! call save_matrix('Mxx.dat',nm_dev, Mxx)
  ! call save_matrix('A.dat',nm_dev, A)
  ! call save_matrix('Kxx.dat',nm_dev, Kxx)
  ! call save_matrix('L0xx.dat',nm_dev, L0xx)  
  !
  A(:,:) = A(:,:) - Mxx(:,:)
  call invert(A,nm_dev)
  Sxx(:,:) = L0xx(:,:) - Sxx(:,:)
  call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,A,nm_dev,Sxx,nm_dev,czero,P_retarded,nm_dev) 
  !
  ! call save_matrix('pr.dat',nm_dev, P_retarded)  
  !
  deallocate(L0xx,Mxx,Sxx,Kxx,A)
end subroutine green_bse_solve


! save a complex matrix to a file in row-column-value format for non-zero entries
subroutine save_matrix(filename, nm, Mat)
    character(len=*),intent(in)::filename ! file name
    integer,intent(in)::nm ! number of bands
    complex(dp), intent(in) :: Mat(nm,nm) ! matrix        
    ! ----
    integer::i,j        
    open(unit=11, file=filename, status='unknown')                
    do i=1,nm
        do j=1,nm
            if ( abs(Mat(i,j)) .gt. 0.0d0 ) then
                write(11,'(2I10,2E18.6)') i,j,dble(Mat(i,j)),aimag(Mat(i,j))
            endif
        enddo
        write(11,*)
    enddo        
    close(11)
end subroutine save_matrix

end module green_bse
