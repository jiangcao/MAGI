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

module rgf_mod

use parameters_mod, only: dp, c1i, cone, czero, pi, twopi, two, BOLTZ

implicit none 

private

CONTAINS


!! RGF with fixed blocksize
subroutine RGF_G(temperature,nm,nx,E,mu,Hii,H1i,sigma_lesser_ph,sigma_greater_ph,sigma_r_ph,&
   tr,tre,ndens,pdens,ldos,cur,GRi1,GLi1,GGi1)     
    integer,intent(in)::nm,nx
    real(8),intent(in) :: E,mu(2),temperature(2)
    complex(8),intent(in),dimension(nm,nm,Nx) :: sigma_lesser_ph,sigma_greater_ph,sigma_r_ph ! diag blocks of scattering SE
    complex(8),intent(in),dimension(nm,nm,Nx) :: Hii,H1i ! diag blocks of Overlap and H and 1st off-diag blocks of H
    real(8),intent(inout) :: tr,tre ! current spectrum on the Right and Left contacts
    complex(8),intent(inout),dimension(nm,nm,Nx),optional :: ldos,ndens,pdens ! diag blocks of GFs    
    real(8),intent(inout),dimension(nm,nm,Nx),optional :: cur ! current density
    complex(8),intent(inout),dimension(nm,nm,Nx),optional :: GRi1,GLi1,GGi1 ! off-diag blocks (i,i+1) of GFs    
    ! ------- local variables
    complex(8) :: H00(nm,nm),H10(nm,nm),A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),S00(nm,nm),G00(nm,nm),GBB(nm,nm),GN0(nm,nm),Gn(nm,nm),Gp(nm,nm)
    complex(8) :: sig(nm,nm),sigmal(nm,nm),sigmar(nm,nm)
    complex(8) :: z
    integer::i,j,k,l
    real(8)::tim,mul,mur,templ,tempr
    complex(8), allocatable :: Gl(:,:,:),Gln(:,:,:),Glp(:,:,:) ! left-connected green function
    complex(8), parameter :: alpha = cone
    complex(8), parameter :: beta  = czero
    !
    mul=mu(1)
    mur=mu(2)
    templ=temperature(1)
    tempr=temperature(2)
    !
    z=dcmplx(E,0.0d-6)
    !
    allocate(Gl(nm,nm,nx))
    allocate(Gln(nm,nm,nx))
    allocate(Glp(nm,nm,nx))
    !
    Gln=0.0D0
    Glp=0.0D0
    Gl=0.0D0
    ldos=0.0d0
    ndens=0.0d0
    pdens=0.0d0
    cur=0.0d0
    S00=0.0d0
    do i=1,nm
        S00(i,i)=1.0d0
    enddo
    ! self energy on the left contact        
    H00(:,:)=Hii(:,:,1)+sigma_r_ph(:,:,1)
    H10(:,:)=H1i(:,:,1)
    !
    call sancho(NM,E,S00,H00,transpose(conjg(H10)),G00,GBB)
    !
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmal,nm)      
    sig(:,:)=-(sigmal(:,:)-transpose(conjg(sigmal(:,:))))*ferm((E-mul)/(BOLTZ*TEMPl))+sigma_lesser_ph(:,:,1)
    A=z*S00-H00-sigmal
    !                
    call invert(A,nm)
    Gl(:,:,1)=A(:,:)
    !
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm) 
    Gln(:,:,1)=C(:,:)
    !-------------------------
    ! forward pass
    Do l=2,nx-1                
        H00(:,:)=Hii(:,:,l)+sigma_r_ph(:,:,l)
        H10(:,:)=H1i(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,l-1),nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
        A=z*S00-H00-C   
        !
        call invert(A,nm)
        Gl(:,:,l)=A(:,:)
        !
        sig=Gln(:,:,l-1)
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)             
        C(:,:)=C(:,:)+sigma_lesser_ph(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,C,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,A,nm,beta,Gn,nm)
        Gln(:,:,l)=Gn(:,:)
    enddo
    ! self energy on the right contact        
    H00(:,:)=Hii(:,:,nx)+sigma_r_ph(:,:,nx)
    H10(:,:)=H1i(:,:,nx)
    !
    call sancho(NM,E,S00,H00,H10,G00,GBB)
    !
    call zgemm('c','n',nm,nm,nm,alpha,H10,nm,G00,nm,beta,A,nm) 
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,H10,nm,beta,sigmar,nm)  
    H10(:,:)=H1i(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,Gl(:,:,nx-1),nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)
    G00=z*S00-H00-sigmar-C   
    !
    call invert(G00,nm)
    !
    ldos(:,:,nx)=G00(:,:)
    sig=Gln(:,:,nx-1)
    call zgemm('n','n',nm,nm,nm,alpha,H10,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,H10,nm,beta,C,nm)  ! C=H10 Gl< H01    
    ! B=Sig< S00
    sig(:,:)=-(sigmar(:,:)-transpose(conjg(sigmar(:,:))))*ferm((E-mur)/(BOLTZ*TEMPl))+C(:,:)+sigma_lesser_ph(:,:,nx)
    call zgemm('n','n',nm,nm,nm,alpha,G00,nm,sig,nm,beta,B,nm) 
    call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,Gn,nm) 
    ! G<00 = G00 sig< G00'
    ndens(:,:,nx)=Gn(:,:)    
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    pdens(:,:,nx)=Gp(:,:)
    A=-(sigmar-transpose(conjg(sigmar)))*ferm((E-mur)/(BOLTZ*TEMPl))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmar-transpose(conjg(sigmar)))*(ferm((E-mur)/(BOLTZ*TEMPl))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm        
      tim=tim-dble(B(i,i)-C(i,i))        
    enddo
    tr=tim
    ! transmission
    !-------------------------
    ! backward pass
    do l=nx-1,1,-1
        H10(:,:)=H1i(:,:,l)
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,H10,nm,Gl(:,:,l),nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,B,nm,beta,C,nm) 
        A=Gln(:,:,l)          
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,A,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,G00,nm,B,nm,beta,A,nm)
        B=C+A
        call zgemm('c','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i+1,i =     
        cur(:,:,l)=dble(A(:,:))     
        !-------------------------
        A=Gn
        call zgemm('n','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G<
        A=Gln(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g< H10 G'
        B=C+A
        if (present(GLi1)) then 
          GLi1(:,:,l)=B
        endif
        call zgemm('n','n',nm,nm,nm,alpha,H10,nm,B,nm,beta,A,nm)      !!! G<_i,i+1 = 
        cur(:,:,l)=cur(:,:,l)-dble(A(:,:))        
        !-------------------------
        if (present(GGi1)) then
          A=Gp
          call zgemm('c','c',nm,nm,nm,alpha,Gl(:,:,l),nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)   ! g H10 G>
          A=Glp(:,:,l)
          call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
          call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,A,nm) ! g> H10 G
          B=C+A         
          GGi1(:,:,l)=B
        endif        
        !-------------------------
        D(:,:)= Gl(:,:,l)
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,G00,nm,beta,GN0,nm)      !!! G_i,i+1 = 
        if (present(GRi1)) then 
          GRi1(:,:,l)=GN0
        endif
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,A,nm)
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)     
        G00(:,:)=Gl(:,:,l)+C(:,:)                                       !!! G_i,i = 
        ldos(:,:,l)=G00(:,:)
        !-------------------------
        A(:,:)=Gn(:,:)     
        call zgemm('n','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
        Gn(:,:)= Gln(:,:,l) + C(:,:)
        A(:,:)=Gln(:,:,l)
        call zgemm('n','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)         
        Gn(:,:)= Gn(:,:)+C(:,:)!                     			 
        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
        call zgemm('n','c',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)          
        Gn(:,:)= Gn(:,:)+C(:,:)!     					 !!! G<_i,i = 
        !-------------------------
!        A(:,:)=Gp(:,:)
!        call zgemm('c','c',nm,nm,nm,alpha,D,nm,H10,nm,beta,B,nm)  
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        call zgemm('n','n',nm,nm,nm,alpha,C,nm,H10,nm,beta,A,nm)
!        call zgemm('n','n',nm,nm,nm,alpha,A,nm,D,nm,beta,C,nm)
!        !
!        Gp(:,:)= Glp(:,:,l) + C(:,:)
!        A(:,:)=Glp(:,:,l)
!        call zgemm('c','n',nm,nm,nm,alpha,GN0,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,A,nm,beta,C,nm)     
!        !
!        Gp(:,:)= Gp(:,:)+C(:,:)!                     			 
!        call zgemm('n','c',nm,nm,nm,alpha,A,nm,H10,nm,beta,B,nm) 
!        call zgemm('n','n',nm,nm,nm,alpha,B,nm,GN0,nm,beta,C,nm)     
!        Gp(:,:)= Gp(:,:)+C(:,:)!     					 !!! G>_i,i
        !-------------------------    
        ndens(:,:,l)=Gn(:,:)
!        pdens(:,:,l)=Gp(:,:)
        pdens(:,:,l)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    enddo
    !
    Gp(:,:)=Gn(:,:)+(G00(:,:)-transpose(conjg(G00(:,:))))
    A=-(sigmal-transpose(conjg(sigmal)))*ferm((E-mul)/(BOLTZ*TEMPr))
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gp,nm,beta,B,nm)
    A=-(sigmal-transpose(conjg(sigmal)))*(ferm((E-mul)/(BOLTZ*TEMPr))-1.0d0)
    call zgemm('n','n',nm,nm,nm,alpha,A,nm,Gn,nm,beta,C,nm)
    tim=0.0d0
    do i=1,nm
      tim=tim+dble(B(i,i)-C(i,i))
    enddo
    tre=tim
end subroutine RGF_G




subroutine RGF_W(NBC,nm,nx,Vii,V1i,PL,PG,PR,WL,WG,WR,PL1i,PG1i,PR1i,WLi1,WGi1,WRi1)
   integer,intent(in)::nm,nx,NBC
   complex(8),intent(in),dimension(nm,nm,nx) :: Vii,V1i ! diagonal and 1st off-diagonal blocks of V  
   complex(8),intent(in),dimension(nm,nm,nx) :: PL,PG,PR ! diagonal blocks of P matrices
   complex(8),intent(in),dimension(nm,nm,nx),optional :: PL1i,PG1i,PR1i ! 1st off-diagonal blocks of P
   complex(8),intent(inout),dimension(nm,nm,nx) :: WL,WG,WR ! diagonal blocks of W 
   complex(8),intent(inout),dimension(nm,nm,nx),optional :: WLi1,WGi1,WRi1 ! 1st off-diagonal blocks of W
   ! ------ local variables
   complex(8)::A(nm,nm),B(nm,nm),C(nm,nm),D(nm,nm),sig(nm,nm),AL(nm,nm),AG(nm,nm),BL(nm,nm),BG(nm,nm)
   complex(8),allocatable,dimension(:,:,:)::S,M,LL,LL1i,LG,LG1i,VV,M1i,Mi1,S1i,Si1,WRi1_,WLi1_,WGi1_,PR1i_,PL1i_,PG1i_
   complex(8),allocatable,dimension(:,:,:)::xlr,wlr,wln,wlp,xR
   complex(8),dimension(:,:),allocatable::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,&
       PL00,PL01,PL10,PG00,PG01,PG10,LL00,LL01,LL10,LG00,LG01,LG10
   complex(8),dimension(:,:),allocatable::VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,&
       MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N
   complex(8),dimension(:,:),allocatable::dM11,xR11,dLL11,dLG11,dV11
   complex(8),dimension(:,:),allocatable::dMnn,xRnn,dLLnn,dLGnn,dVnn
   integer::i,NL,NR,NT,LBsize,RBsize
   integer::ix
   real(8)::condL,condR
   NL=nm ! left contact block size
   NR=nm ! right contact block size
   NT=nm*nx! total size   
   allocate(VV(nm,nm,nx))
   VV = Vii
   allocate(M(nm,nm,nx))
   allocate(M1i(nm,nm,nx))
   allocate(Mi1(nm,nm,nx))
   allocate(S(nm,nm,nx))
   allocate(S1i(nm,nm,nx))
   allocate(Si1(nm,nm,nx))
   allocate(LL(nm,nm,nx))
   allocate(LG(nm,nm,nx))
   allocate(LL1i(nm,nm,nx))
   allocate(LG1i(nm,nm,nx))
   allocate(PR1i_(nm,nm,nx))
   allocate(PL1i_(nm,nm,nx))
   allocate(PG1i_(nm,nm,nx))
   if (present(PR1i)) then
     PR1i_=PR1i
     PL1i_=PL1i
     PG1i_=PG1i
   else
     PR1i_=czero
     PL1i_=czero
     PG1i_=czero
   endif
   ! S = V P^r
   ! Si,i = Vi,i Pi,i + Vi,i+1 Pi+1,i + Vi,i-1 Pi-1,i
   do ix=1,nx
     call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,ix),nm,PR(:,:,ix),nm,czero,S(:,:,ix),nm)
     if (present(PR1i)) then
       call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,ix),nm,PR1i(:,:,ix),nm,cone,S(:,:,ix),nm) 
       if (ix==1) then
         call zgemm('n','t',nm,nm,nm,cone,V1i(:,:,ix),nm,PR1i(:,:,ix),nm,cone,S(:,:,ix),nm) 
       else
         call zgemm('n','t',nm,nm,nm,cone,V1i(:,:,ix-1),nm,PR1i(:,:,ix-1),nm,cone,S(:,:,ix),nm) 
       endif
     endif
   enddo
   !
   do i=1,nx-1
     ! Si+1,i = Vi+1,i+1 Pi+1,i + Vi+1,i Pi,i
     call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,PR(:,:,i),nm,czero,S1i(:,:,i),nm) 
     if (present(PR1i)) then
        call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,PR1i(:,:,i),nm,cone,S1i(:,:,i),nm) 
     endif
     !
     ! Si,i+1 = Vi,i Pi,i+1 + Vi,i+1 Pi+1,i+1
     call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PR(:,:,i+1),nm,czero,Si1(:,:,i),nm)   
     if (present(PR1i)) then
       call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PR1i(:,:,i),nm,cone,Si1(:,:,i),nm)
     endif
     !  
   enddo
   Si1(:,:,nx)=Si1(:,:,nx-1)
   S1i(:,:,nx)=S1i(:,:,nx-1)
   M   = -S
   M1i = -S1i
   Mi1 = -Si1
   do ix=1,nx
     do i=1,nm
         M(i,i,ix) = cone + M(i,i,ix)
     enddo
   enddo
   deallocate(S,S1i,Si1)
   !
   ! LL=V P^< V
   ! LLi,i
   do i=1,nx
     call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PL(:,:,i),nm,czero,A,nm)
     call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LL(:,:,i),nm)
     !
     if (i<nx) then
       call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL(:,:,i+1),nm,czero,A,nm)
       call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL(:,:,i),nm)
     endif
     if (i>1) then
       call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i-1),nm,PL(:,:,i-1),nm,czero,A,nm)
       call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,i-1),nm,cone,LL(:,:,i),nm)
     endif
     !
   !  if (present(PL1i)) then
   !    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL1i(:,:,i),nm,czero,A,nm)
   !    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LL(:,:,i),nm)
   !    !
   !    call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,PL1i(:,:,i),nm,czero,A,nm)
   !    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL(:,:,i),nm)
   !    !
   !    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PL1i(:,:,max(i-1,1)),nm,czero,A,nm)
   !    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LL(:,:,i),nm)      
   !    !
   !    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PL1i(:,:,max(i-1,1)),nm,czero,A,nm)
   !    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LL(:,:,i),nm)
   !  endif
   enddo
   ! LLi+1,i = Vi+1,i P^li,i Vi,i + Vi+1,i+1 P^li+1,i+1 Vi+1,i
   do i=1,nx
     call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,PL(:,:,i),nm,czero,A,nm)
     call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LL1i(:,:,i),nm)
     if (i<nx) then
       call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,PL(:,:,i+1),nm,czero,A,nm)
       call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LL1i(:,:,i),nm)
     endif
   enddo
   ! LG=V P^> V'    
   ! LGi,i
   do i=1,nx
     call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PG(:,:,i),nm,czero,A,nm)
     call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LG(:,:,i),nm)
     !
     if (i<nx) then
       call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG(:,:,i+1),nm,czero,A,nm)
       call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG(:,:,i),nm)
     endif
     if (i>1) then    
       call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i-1),nm,PG(:,:,i-1),nm,czero,A,nm)
       call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,i-1),nm,cone,LG(:,:,i),nm)
     endif
     !
   !  if (present(PG1i)) then
   !    call zgemm('c','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG1i(:,:,i),nm,czero,A,nm)
   !    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LG(:,:,i),nm)
   !    !
   !    call zgemm('n','c',nm,nm,nm,cone,Vii(:,:,i),nm,PG1i(:,:,i),nm,czero,A,nm)
   !    call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG(:,:,i),nm)
   !    !
   !    call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i),nm,PG1i(:,:,max(i-1,1)),nm,czero,A,nm)
   !    call zgemm('n','c',nm,nm,nm,cone,A,nm,V1i(:,:,max(i-1,1)),nm,cone,LG(:,:,i),nm)
   !    !
   !    call zgemm('n','c',nm,nm,nm,cone,V1i(:,:,max(i-1,1)),nm,PG1i(:,:,max(i-1,1)),nm,czero,A,nm)
   !    call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,cone,LG(:,:,i),nm)
   !  endif
   enddo
   !! LGi+1,i = Vi+1,i P^>i,i Vi,i + Vi+1,i+1 P^>i+1,i+1 Vi+1,i 
   do i=1,nx
     call zgemm('n','n',nm,nm,nm,cone,V1i(:,:,i),nm,PG(:,:,i),nm,czero,A,nm)
     call zgemm('n','n',nm,nm,nm,cone,A,nm,Vii(:,:,i),nm,czero,LG1i(:,:,i),nm)
     if (i<nx) then
       call zgemm('n','n',nm,nm,nm,cone,Vii(:,:,i+1),nm,PG(:,:,i+1),nm,czero,A,nm)
       call zgemm('n','n',nm,nm,nm,cone,A,nm,V1i(:,:,i),nm,cone,LG1i(:,:,i),nm)
     endif
   enddo
   !
   if (NBC>0) then ! apply open-boundary correction to all matrices
     !
     LBsize=NL*NBC
     RBsize=NR*NBC
     !
     allocate(V00 (LBsize,LBsize))
     allocate(V01 (LBsize,LBsize))
     allocate(V10 (LBsize,LBsize))
     allocate(M00 (LBsize,LBsize))
     allocate(M01 (LBsize,LBsize))
     allocate(M10 (LBsize,LBsize))
     allocate(PR00(LBsize,LBsize))
     allocate(PR01(LBsize,LBsize))
     allocate(PR10(LBsize,LBsize))
     allocate(PG00(LBsize,LBsize))
     allocate(PG01(LBsize,LBsize))
     allocate(PG10(LBsize,LBsize))
     allocate(PL00(LBsize,LBsize))
     allocate(PL01(LBsize,LBsize))
     allocate(PL10(LBsize,LBsize))
     allocate(LG00(LBsize,LBsize))
     allocate(LG01(LBsize,LBsize))
     allocate(LG10(LBsize,LBsize))
     allocate(LL00(LBsize,LBsize))
     allocate(LL01(LBsize,LBsize))
     allocate(LL10(LBsize,LBsize))
     allocate(dM11(LBsize,LBsize))
     allocate(xR11(LBsize,LBsize))
     allocate(dV11(LBsize,LBsize))
     allocate(dLL11(LBsize,LBsize))
     allocate(dLG11(LBsize,LBsize))
     !  
     allocate(VNN  (RBsize,RBsize))
     allocate(VNN1 (RBsize,RBsize))
     allocate(Vn1n (RBsize,RBsize))
     allocate(Mnn  (RBsize,RBsize))
     allocate(Mnn1 (RBsize,RBsize))
     allocate(Mn1n (RBsize,RBsize))
     allocate(PRnn (RBsize,RBsize))
     allocate(PRnn1(RBsize,RBsize))
     allocate(PRn1n(RBsize,RBsize))
     allocate(PGnn (RBsize,RBsize))
     allocate(PGnn1(RBsize,RBsize))
     allocate(PGn1n(RBsize,RBsize))
     allocate(PLnn (RBsize,RBsize))
     allocate(PLnn1(RBsize,RBsize))
     allocate(PLn1n(RBsize,RBsize))
     allocate(LGnn (RBsize,RBsize))
     allocate(LGnn1(RBsize,RBsize))
     allocate(LGn1n(RBsize,RBsize))
     allocate(LLnn (RBsize,RBsize))
     allocate(LLnn1(RBsize,RBsize))
     allocate(LLn1n(RBsize,RBsize))
     allocate(dMnn (RBsize,RBsize))
     allocate(xRnn (RBsize,RBsize))
     allocate(dLLnn(RBsize,RBsize))
     allocate(dLGnn(RBsize,RBsize))
     allocate(dVnn(RBsize,RBsize))
     ! left OBC
     call get_OBC_blocks_for_W(NL,Vii(:,:,1),V1i(:,:,1),PR(:,:,1),PR1i_(:,:,1),&
         PL(:,:,1),PL1i_(:,:,1),PG(:,:,1),PG1i_(:,:,1),NBC,&
         V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
         LL00,LL01,LL10,LG00,LG01,LG10)
     ! right OBC  
     call get_OBC_blocks_for_W(NR,Vii(:,:,nx),transpose(conjg(V1i(:,:,nx))),PR(:,:,nx),&
         transpose(PR1i_(:,:,nx)),PL(:,:,nx),-transpose(conjg(PL1i_(:,:,nx))),&
         PG(:,:,nx),-transpose(conjg(PG1i_(:,:,nx))),NBC,&
         VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,&
         LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N)
     !
     ! Correct first and last blocks of M to account for elements in the contacts
     M(:,:,1)=M(:,:,1) - matmul(V10,PR01)
     M(:,:,nx)=M(:,:,nx) - matmul(VNN1,PRN1N)
     ! 
     ! Correct first and last blocks of LL to account for elements in the contacts
     LL(:,:,1)=LL(:,:,1) + matmul(matmul(V10,PL00),V01) + &
       matmul(matmul(V10,PL01),V00) + matmul(matmul(V00,PL10),V01)
     !  
     LL(:,:,nx)=LL(:,:,nx) + &
       matmul(matmul(VNN,PLNN1),VN1N) + matmul(matmul(VNN1,PLN1N),VNN) + &
       matmul(matmul(VNN1,PLNN),VN1N)
     !  
     ! Correct first and last blocks of LG to account for elements in the contacts
     LG(:,:,1)=LG(:,:,1) + matmul(matmul(V10,PG00),V01) + &
       matmul(matmul(V10,PG01),V00) + matmul(matmul(V00,PG10),V01)
     LG(:,:,nx)=LG(:,:,nx) + &
       matmul(matmul(VNN,PGNN1),VN1N) + matmul(matmul(VNN1,PGN1N),VNN) + matmul(matmul(VNN1,PGNN),VN1N)
     !  
     ! WR/WL/WG OBC Left
     call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
     ! WR/WL/WG OBC right
     call open_boundary_conditions(NR,MNN,MNN1,MN1N,VN1N,xRNN,dMNN,dVNN,condR)
     !
     if (condL<1.0d-6) then   
       !
       ! call get_dL_OBC_for_W(NL,xR11,LL00,LL01,LG00,LG01,M10,'L', dLL11,dLG11)
       !
       M(:,:,1)=M(:,:,1) - dM11
       VV(:,:,1)=VV(:,:,1) - dV11    
!       LL(:,:,1)=LL(:,:,1) + dLL11
!       LG(:,:,1)=LG(:,:,1) + dLG11    
     endif
     if (condR<1.0d-6) then    
       !
       ! call get_dL_OBC_for_W(NR,xRNN,LLNN,LLN1N,LGNN,LGN1N,MNN1,'R', dLLNN,dLGNN)
       !
       M(:,:,nx)=M(:,:,nx) - dMNN
       VV(:,:,nx)=VV(:,:,nx)- dVNN
!       LL(:,:,nx)=LL(:,:,nx) + dLLNN
!       LG(:,:,nx)=LG(:,:,nx) + dLGNN    
     endif
     !
     deallocate(V00,V01,V10)
     deallocate(M00,M01,M10)
     deallocate(PR00,PR01,PR10)
     deallocate(PG00,PG01,PG10)
     deallocate(PL00,PL01,PL10)
     deallocate(LG00,LG01,LG10)
     deallocate(LL00,LL01,LL10)
     deallocate(VNN,VNN1,Vn1n)
     deallocate(Mnn,Mnn1,Mn1n)
     deallocate(PRnn,PRnn1,PRn1n)
     deallocate(PGnn,PGnn1,PGn1n)
     deallocate(PLnn,PLnn1,PLn1n)
     deallocate(LGnn,LGnn1,LGn1n)
     deallocate(LLnn,LLnn1,LLn1n)
     deallocate(dM11,xR11,dLL11,dLG11,dV11)
     deallocate(dMnn,xRnn,dLLnn,dLGnn,dVnn)
     !
   end if ! NBC>0
   !
   allocate(xlr(nm,nm,nx)) ! left-connected xR
   allocate(wlr(nm,nm,nx)) ! left-connected Wr
   allocate(wln(nm,nm,nx)) ! left-connected W<
   allocate(wlp(nm,nm,nx)) ! left-connected W>
   allocate(xR(nm,nm,nx))  ! full xR
   allocate(WRi1_(nm,nm,nx)) ! off-diagonal block of WR
   !print*,' first pass, from right to left ... '
   A=M(:,:,nx)
   call invert(A,nm)
   xlr(:,:,nx)=A
   call zgemm('n','n',nm,nm,nm,cone,A,nm,VV(:,:,nx),nm,czero,Wlr(:,:,nx),nm) 
   !
   call zgemm('n','n',nm,nm,nm,cone,A,nm,LL(:,:,nx),nm,czero,B,nm) 
   call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wln(:,:,nx),nm) 
   !
   call zgemm('n','n',nm,nm,nm,cone,A,nm,LG(:,:,nx),nm,czero,B,nm) 
   call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,Wlp(:,:,nx),nm) 
   !
   do ix=nx-1,1,-1
     call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,xlr(:,:,ix+1),nm,czero,B,nm) ! B -> Mn,n+1 xlr -> MxR
     call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,C,nm)             ! C -> Mn,n+1 xlr Mn+1,n
     A = M(:,:,ix) - C
     call invert(A,nm) ! A -> xlr
     xlr(:,:,ix)=A     
     call zgemm('n','n',nm,nm,nm,cone,B,nm,V1i(:,:,ix),nm,czero,C,nm) ! C -> MxR Vn+1,n  
     D=VV(:,:,ix) - C                                                 ! D -> Vn,n - Mn,n+1 xlr Vn+1,n
     call zgemm('n','n',nm,nm,nm,cone,A,nm,D,nm,czero,Wlr(:,:,ix),nm) ! xlr ( Vn,n - Mn,n+1 xlr Vn+1,n )
     !
   !  !AL=MxR*LL10  
   !  !AG=MxR*LG10          
   !  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,xlr(:,:,ix+1),nm,czero,B,nm) ! B -> Mn,n+1 xlr -> MxR
   !  AL=matmul(B, LL1i(:,:,ix))
   !  AG=matmul(B, LG1i(:,:,ix))
   !  !  
   !  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,wln(:,:,ix+1),nm,czero,B,nm) 
   !  call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,sig,nm)   ! sig -> Mn,n+1 w<n+1 Mn+1,n
   !  C = LL(:,:,ix) + sig - (AL-transpose(conjg(AL)))
   !  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
   !  call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,wln(:,:,ix),nm)       
   !  !
   !  call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix),nm,wlp(:,:,ix+1),nm,czero,B,nm) 
   !  call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix),nm,czero,sig,nm)   ! sig -> Mn,n+1 w>n+1 Mn+1,n    
   !  C = LG(:,:,ix) + sig - (AG-transpose(conjg(AG)))
   !  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,B,nm) 
   !  call zgemm('n','c',nm,nm,nm,cone,B,nm,A,nm,czero,wlp(:,:,ix),nm)     
     !  
   enddo
   !print*,' second pass, from left to right ... '
   WR(:,:,1)=Wlr(:,:,1)
   call zgemm('n','t',nm,nm,nm,cone,WR(:,:,1),nm,M1i(:,:,1),nm,czero,B,nm)
   C = transpose((V1i(:,:,1))) - B
   call zgemm('n','t',nm,nm,nm,cone,C,nm,xlr(:,:,2),nm,czero,B,nm)  ! B -> (V12 - WR1 M12) xlr2
   WRi1_(:,:,1) = B
   XR(:,:,1)=xlr(:,:,1)
   !WL(:,:,1)=Wln(:,:,1)
   !WG(:,:,1)=Wlp(:,:,1)
   do ix=2,nx
     call zgemm('n','n',nm,nm,nm,cone,xlr(:,:,ix),nm,M1i(:,:,ix-1),nm,czero,B,nm) ! B -> xlr1 * M10
     call zgemm('n','n',nm,nm,nm,cone,B,nm,WRi1_(:,:,ix-1),nm,czero,C,nm)         ! C -> xlr1 * M10 WR0 M01 xlr1
     WR(:,:,ix)=Wlr(:,:,ix)  - C  
     call zgemm('n','n',nm,nm,nm,cone,B,nm,XR(:,:,ix-1),nm,czero,A,nm)            ! A -> xlr1 * M10 * XR0 
     call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix-1),nm,xlr(:,:,ix),nm,czero,C,nm) ! C -> M01 * xlr1
     call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)                       ! D -> xlr1 * M10 * XR0 * M01 * xlr1
     XR(:,:,ix)=xlr(:,:,ix)  + D
     !   
   !  ! AL=xlr1 LL10 XR0' M01 xlr1'
   !  ! BL=xlr1 M10 XR0 M01 wln1  
   !  call zgemm('c','n',nm,nm,nm,cone,xR(:,:,ix-1),nm,Mi1(:,:,ix-1),nm,czero,D,nm)   ! D -> XR0' M01
   !  call zgemm('n','c',nm,nm,nm,cone,D,nm,xlr(:,:,ix),nm,czero,A,nm)    ! A -> XR0' M01 xlr1'
   !  !
   !  AL = matmul( matmul(xlr(:,:,ix), LL1i(:,:,ix-1)), A )
   !  AG = matmul( matmul(xlr(:,:,ix), LG1i(:,:,ix-1)), A )
   !  !
   !  call zgemm('n','n',nm,nm,nm,cone,B,nm,xR(:,:,ix-1),nm,czero,D,nm)   ! D -> xlr1 * M10 * XR0
   !  call zgemm('n','n',nm,nm,nm,cone,D,nm,Mi1(:,:,ix-1),nm,czero,A,nm)  ! A -> xlr1 * M10 * XR0 * M01 
   !  call zgemm('n','n',nm,nm,nm,cone,A,nm,wln(:,:,ix),nm,czero,BL,nm)   ! BL-> xlr1 * M10 * XR0 * M01 * wln1
   !  !
   !  call zgemm('n','n',nm,nm,nm,cone,B,nm,WL(:,:,ix-1),nm,czero,A,nm)            ! A -> xlr1 M10 WL0
   !  call zgemm('n','c',nm,nm,nm,cone,Mi1(:,:,ix-1),nm,xlr(:,:,ix),nm,czero,C,nm) ! C -> M01 xlr1'
   !  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)            ! D -> xlr1 M10 WL0 M01 xlr1'
   !  WL(:,:,ix) = wln(:,:,ix) !+ D + (BL-transpose(conjg(BL))) - (AL-transpose(conjg(AL)))
   !  !
   !  ! AG=xlr1 LG10 XR0' M01 xlr1'
   !  ! BG=xlr1 M10 XR0 M01 wlp1    
   !  !
   !  call zgemm('n','n',nm,nm,nm,cone,B,nm,xR(:,:,ix-1),nm,czero,D,nm)    ! D -> xlr1 * M10 * XR0 
   !  call zgemm('n','n',nm,nm,nm,cone,D,nm,Mi1(:,:,ix-1),nm,czero,A,nm)   ! A -> xlr1 * M10 * XR0 * M01
   !  call zgemm('n','n',nm,nm,nm,cone,A,nm,wlp(:,:,ix),nm,czero,BG,nm)    ! BG-> xlr1 * M10 * XR0 * M01 * wlp1
   !  !
   !  call zgemm('n','n',nm,nm,nm,cone,B,nm,WG(:,:,ix-1),nm,czero,A,nm)    ! A -> xlr1 * M10 * WG0 
   !  call zgemm('n','n',nm,nm,nm,cone,A,nm,C,nm,czero,D,nm)               ! D -> xlr1 * M10 * WG0 * M01 xlr1' 
   !  WG(:,:,ix) = wlp(:,:,ix) !+ D + (BG-transpose(conjg(BG))) - (AG-transpose(conjg(AG)))
     !
     if (ix<nx) then
       call zgemm('n','t',nm,nm,nm,cone,WR(:,:,ix),nm,M1i(:,:,ix),nm,czero,B,nm)
       C = transpose((V1i(:,:,ix))) - B
       call zgemm('n','t',nm,nm,nm,cone,C,nm,xlr(:,:,ix+1),nm,czero,B,nm)  ! B -> (V12 - WR1 M12) xlr2
       WRi1_(:,:,ix) = B        
     endif
     !     
   enddo
   !
   if (present(WRi1)) then
     WRi1= WRi1_
   end if
   !
   WRi1_ = dcmplx(0.0d0*dble(WRi1_),aimag(WRi1_))
   do ix=1,nx
     call zgemm('n','n',nm,nm,nm,cone,WR(:,:,ix),nm,PL(:,:,ix),nm,czero,A,nm)
     call zgemm('n','c',nm,nm,nm,cone,A,nm,WR(:,:,ix),nm,czero,B,nm)
     if (ix<nx) then
       call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix),nm,PL(:,:,ix+1),nm,czero,A,nm)
       call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix),nm,cone,B,nm)
     else
       call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix-1),nm,PL(:,:,nx),nm,czero,A,nm)
       call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix-1),nm,cone,B,nm)  
     endif
     call zgemm('t','n',nm,nm,nm,cone,WRi1_(:,:,max(ix-1,1)),nm,PL(:,:,max(ix-1,1)),nm,czero,A,nm)
     call zgemm('n','n',nm,nm,nm,cone,A,nm,-WRi1_(:,:,max(ix-1,1)),nm,cone,B,nm)
     WL(:,:,ix)=B
     call zgemm('n','n',nm,nm,nm,cone,WR(:,:,ix),nm,PG(:,:,ix),nm,czero,A,nm)
     call zgemm('n','c',nm,nm,nm,cone,A,nm,WR(:,:,ix),nm,czero,B,nm)
     if (ix<nx) then
       call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix),nm,PG(:,:,ix+1),nm,czero,A,nm)
       call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix),nm,cone,B,nm)
     else
       call zgemm('n','n',nm,nm,nm,cone,WRi1_(:,:,ix-1),nm,PG(:,:,nx),nm,czero,A,nm)
       call zgemm('n','t',nm,nm,nm,cone,A,nm,-WRi1_(:,:,ix-1),nm,cone,B,nm)  
     endif
     call zgemm('t','n',nm,nm,nm,cone,WRi1_(:,:,max(ix-1,1)),nm,PG(:,:,max(ix-1,1)),nm,czero,A,nm)
     call zgemm('n','n',nm,nm,nm,cone,A,nm,-WRi1_(:,:,max(ix-1,1)),nm,cone,B,nm)
     WG(:,:,ix)=B
   enddo
   !
   deallocate(M,M1i,Mi1,LL,LG,VV,LL1i,LG1i)
   deallocate(wln,wlp,wlr,xlr,Xr,WRi1_)
   deallocate(PR1i_,PL1i_,PG1i_)
end subroutine RGF_W



end module rgf_mod
