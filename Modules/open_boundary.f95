
module open_boundary
    use parameters_mod
    implicit none 
    contains

    ! Sancho-Rubio 
    subroutine sancho(nm,E,S00,H00,H10,G00,GBB)
        use parameters_mod
        complex(8), parameter :: alpha = dcmplx(1.0d0,0.0d0)
        complex(8), parameter :: beta  = dcmplx(0.0d0,0.0d0)
        integer i,j,k,nmax
        COMPLEX(8) :: z
        real(8),intent(in) :: E
        real(8) :: error
        REAL(8) :: TOL=1.0D-10  ! [eV]
        integer, intent(in) :: nm
        COMPLEX(8), INTENT(IN) ::  S00(nm,nm), H00(nm,nm), H10(nm,nm)
        COMPLEX(8), INTENT(OUT) :: G00(nm,nm), GBB(nm,nm)
        COMPLEX(8), ALLOCATABLE :: A(:,:), B(:,:), C(:,:), tmp(:,:), G(:,:)
        COMPLEX(8), ALLOCATABLE :: H_BB(:,:), H_SS(:,:), H_01(:,:), H_10(:,:), Id(:,:)
        !COMPLEX(8), ALLOCATABLE :: WORK(:)
        !COMPLEX(8), EXTERNAL :: ZLANGE
        Allocate( H_BB(nm,nm) )
        Allocate( H_SS(nm,nm) )
        Allocate( H_01(nm,nm) )
        Allocate( H_10(nm,nm) )
        Allocate( Id(nm,nm) )
        Allocate( A(nm,nm) )
        Allocate( B(nm,nm) )
        Allocate( C(nm,nm) )
        Allocate( G(nm,nm) )
        Allocate( tmp(nm,nm) )
        nmax=200
        z = dcmplx(E,1.0d-5)
        Id=0.0d0
        tmp=0.0d0
        do i=1,nm
        Id(i,i)=1.0d0
        tmp(i,i)=dcmplx(0.0d0,1.0d0)
        enddo
        H_BB = H00
        H_10 = H10
        H_01 = TRANSPOSE( CONJG( H_10 ) )
        H_SS = H00
        do i = 1, nmax
        A = z*S00 - H_BB
        !
        call invert_inplace(A,nm)      
        !
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_10,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 
        H_SS = H_SS + C
        H_BB = H_BB + C
        call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,C,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,A,nm,H_01,nm,beta,B,nm) 
        call zgemm('n','n',nm,nm,nm,alpha,H_10,nm,B,nm,beta,A,nm)  
        H_10 = C    
        H_BB = H_BB + A
        call zgemm('n','n',nm,nm,nm,alpha,H_01,nm,B,nm,beta,C,nm) 
        H_01 = C 
        ! NORM --> inspect the diagonal of A
        error=0.0d0
        DO k=1,nm
        DO j=1,nm
            error=error+sqrt(aimag(C(k,j))**2+Dble(C(k,j))**2)
        END DO
        END DO
        tmp=H_SS
        IF ( abs(error) < TOL ) THEN
            EXIT
        ELSE
        END IF
        IF (i .EQ. nmax) THEN
            write(*,*) 'SEVERE warning: nmax reached in sancho!!!',error        
            H_SS=H00
            H_BB=H00
        END IF
        enddo
        G00 = z*S00 - H_SS
        !
        call invert_inplace(G00,nm)    
        !
        GBB = z*S00 - H_BB
        !
        call invert_inplace(GBB,nm)    
        !
        Deallocate( tmp )
        Deallocate( A )
        Deallocate( B )
        Deallocate( C )
        Deallocate( H_BB )
        Deallocate( H_SS )
        Deallocate( H_01 )
        Deallocate( H_10 )
        Deallocate( Id )
        Deallocate( G )
    end subroutine sancho



    ! matrix inversion
    subroutine invert_inplace(A, nn)
        integer :: info, nn
        integer, dimension(:), allocatable :: ipiv
        complex(8), dimension(nn, nn), intent(inout) :: A
        complex(8), dimension(:), allocatable :: work
        COMPLEX(8), PARAMETER :: czero  = dcmplx(0.0d0,0.0d0)
        allocate (work(nn*nn))
        allocate (ipiv(nn))
        call zgetrf(nn, nn, A, nn, ipiv, info)
        if (info .ne. 0) then
            print *, 'SEVERE warning: zgetrf failed, info=', info
            A = czero
        else
            call zgetri(nn, A, nn, ipiv, work, nn*nn, info)
            if (info .ne. 0) then
                print *, 'SEVERE warning: zgetri failed, info=', info
                A = czero
            end if
        end if
        deallocate (work)
        deallocate (ipiv)
    end subroutine invert_inplace


    subroutine identity(A,n)
        integer, intent(in) :: n        
        complex(8), dimension(n,n), intent(out) :: A
        integer :: i
        A = dcmplx(0.0d0,0.0d0)
        do i = 1,n
          A(i,i) = dcmplx(1.0d0,0.0d0)
        end do
    end subroutine identity



    subroutine expand_size_bycopy(A,nm,nb,add)
        complex(8),intent(inout)::A(nm,nm)
        integer, intent(in)::nm,add,nb
        integer::i,nm0,l,l2
        nm0=nm-nb*add*2
        A(1:add*nb,:)=0.0d0
        A(:,1:add*nb)=0.0d0
        A(add*nb+nm0+1:nm,:)=0.0d0
        A(:,add*nb+nm0+1:nm)=0.0d0
        do i=0,add-1
          A(i*nb+1:i*nb+nb,i*nb+1:i*nb+nm0)=A(add*nb+1:add*nb+nb,add*nb+1:add*nb+nm0)
          A(i*nb+1:i*nb+nm0,i*nb+1:i*nb+nb)=A(add*nb+1:add*nb+nm0,add*nb+1:add*nb+nb)
          l=add*nb+nm0+i*nb
          l2=add*nb+i*nb+nb
          A(l+1:l+nb,l2+1:l2+nm0)=A(add*nb+nm0-nb+1:add*nb+nm0,add*nb+1:add*nb+nm0)
          A(l2+1:l2+nm0,l+1:l+nb)=A(add*nb+1:add*nb+nm0,add*nb+nm0-nb+1:add*nb+nm0)  
        enddo
    end subroutine expand_size_bycopy


    ! calculate matrix blocks for the Open Boundary Condition of W
    subroutine get_OBC_blocks_for_W(n,v_00,v_01,pR_00,pR_01,pL_00,pL_01,pG_00,pG_01,NBC,&
        V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
        LL00,LL01,LL10,LG00,LG01,LG10)
        integer,intent(in)::n,NBC
        complex(8),intent(in),dimension(n,n)::v_00,v_01,pR_00,pR_01,pL_00,pL_01,pG_00,pG_01
        complex(8),intent(out),dimension(n*NBC,n*NBC)::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
            LL00,LL01,LL10,LG00,LG01,LG10
        complex(8),dimension(n*NBC,n*NBC)::II
    !
        select case (NBC)
        !
        case(1)
            !
            V00=v_00
            V01=v_01
            V10=transpose(conjg(V01))
            !
            PR00=pR_00
            PR01=pR_01
            PR10=transpose(PR01)
            !
            PL00=pL_00;
            PL01=pL_01;
            PL10= - transpose(conjg(PL01))
            !
            PG00=pG_00
            PG01=pG_01
            PG10= - transpose(conjg(PG01))
                
        case(2)
            !
            V00(1:n,1:n)=v_00 
            V00(1:n,n+1:2*n)=v_01
            V00(n+1:2*n,1:n)=transpose(conjg(v_01))
            V00(n+1:2*n,n+1:2*n)= v_00
            V01=czero
            V01(n+1:2*n,1:n)=v_01
            V10=transpose(conjg(V01))
            !
            PR00(1:n,1:n)=pR_00 
            PR00(1:n,n+1:2*n)=pR_01
            PR00(n+1:2*n,1:n)=transpose(pR_01)
            PR00(n+1:2*n,n+1:2*n)= pR_00
            PR01=czero
            PR01(n+1:2*n,1:n)=pR_01
            PR10=transpose(PR01)
            !
            PG00(1:n,1:n)=pG_00 
            PG00(1:n,n+1:2*n)=pG_01
            PG00(n+1:2*n,1:n)=-transpose(conjg(pG_01))
            PG00(n+1:2*n,n+1:2*n)= pG_00
            PG01=czero
            PG01(n+1:2*n,1:n)=pG_01
            PG10=-transpose(conjg(PG01))
            !
            PL00(1:n,1:n)=pL_00 
            PL00(1:n,n+1:2*n)=pL_01
            PL00(n+1:2*n,1:n)=-transpose(conjg(pL_01))
            PL00(n+1:2*n,n+1:2*n)= pL_00
            PL01=czero
            PL01(n+1:2*n,1:n)=pL_01
            PL10=-transpose(conjg(PL01))    
        !
        end select
        !
        call identity(II,NBC*N)
        M00=II*dcmplx(1.0d0,1d-10)-matmul(V10,PR01)
        M00=M00-matmul(V00,PR00)-matmul(V01,PR10)
        M01=-matmul(V00,PR01)-matmul(V01,PR00)
        M10=-matmul(V10,PR00)-matmul(V00,PR10)
        !
        LL00=matmul(matmul(V10,PL00),V01)+matmul(matmul(V10,PL01),V00)
        LL00=LL00+matmul(matmul(V00,PL10),V01)+matmul(matmul(V00,PL00),V00)
        LL00=LL00+matmul(matmul(V00,PL01),V10)
        LL00=LL00+matmul(matmul(V01,PL10),V00)+matmul(matmul(V01,PL00),V10)
        !
        LL01=matmul(matmul(V10,PL01),V01)+matmul(matmul(V00,PL00),V01)
        LL01=LL01+matmul(matmul(V00,PL01),V00)+matmul(matmul(V01,PL10),V01)
        LL01=LL01+matmul(matmul(V01,PL00),V00)
        !
        LL10=-transpose(conjg(LL01))
        !
        LG00=matmul(matmul(V10,PG00),V01)+matmul(matmul(V10,PG01),V00)
        LG00=LG00+matmul(matmul(V00,PG10),V01)+matmul(matmul(V00,PG00),V00)
        LG00=LG00+matmul(matmul(V00,PG01),V10)
        LG00=LG00+matmul(matmul(V01,PG10),V00)+matmul(matmul(V01,PG00),V10)
        !
        LG01=matmul(matmul(V10,PG01),V01)
        LG01=LG01+matmul(matmul(V00,PG00),V01)+matmul(matmul(V00,PG01),V00)
        LG01=LG01+matmul(matmul(V01,PG10),V01)+matmul(matmul(V01,PG00),V00)
        !
        LG10=-transpose(conjg(LG01))  
    end subroutine get_OBC_blocks_for_W


    ! calculate corrections to the L matrix blocks for the Open Boundary Condition
    subroutine get_dL_OBC_for_W(nm,xR,LL00,LL01,LG00,LG01,M10,typ, dLL11,dLG11)
        integer,intent(in)::nm
        character(len=*),intent(in)::typ
        complex(8),intent(in),dimension(nm,nm)::xR,LL00,LL01,LG00,LG01,M10
        complex(8),intent(out),dimension(nm,nm)::dLL11,dLG11
        ! -----
        complex(8),dimension(nm,nm)::AL,AG,FL,FG,A,V,iV,yL_NN,wL_NN,yG_NN,wG_NN,tmp1,tmp2
        complex(8),dimension(nm)::E
        integer::i,j
        !!!! AL=M10*xR*LL01;
        !!!! AG=M10*xR*LG01;
        call zgemm('n','n',nm,nm,nm,cone,M10,nm,xR,nm,czero,tmp1,nm)
        call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,LL01,nm,czero,AL,nm)
        call zgemm('n','n',nm,nm,nm,cone,M10,nm,xR,nm,czero,tmp1,nm)
        call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,LG01,nm,czero,AG,nm)
        !!!! FL=xR*(LL00-(AL-AL'))*xR';
        !!!! FG=xR*(LG00-(AG-AG'))*xR';
        call zgemm('n','n',nm,nm,nm,cone,xR,nm,(LL00-(AL-transpose(conjg(AL)))),nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,FL,nm)
        call zgemm('n','n',nm,nm,nm,cone,xR,nm,(LG00-(AG-transpose(conjg(AG)))),nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,FG,nm)
        !
        call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,V,nm)  
        do i=1,nm
        V(i,i)=V(i,i)+dcmplx(0.0d0,1.0d-4)  ! 1i*1e-4 added to stabilize matrix
        enddo
        E=eigv(nm,V)
        iV=V
        call invert_inplace(iV,nm)
        !lesser component
        call zgemm('n','n',nm,nm,nm,cone,iV,nm,FL,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,iV,nm,czero,yL_NN,nm)
        yL_NN=yL_NN/(1.0d0 - sum(E*conjg(E)))
        call zgemm('n','n',nm,nm,nm,cone,V,nm,yL_NN,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,V,nm,czero,wL_NN,nm)
        !refinement iteration
        call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,tmp1,nm)
        call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,wL_NN,nm,czero,tmp2,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp2,nm,M10,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,tmp2,nm)
        wL_NN=FL+tmp2
        !
        call zgemm('n','n',nm,nm,nm,cone,M10,nm,wL_NN,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dLL11,nm)
        dLL11=dLL11-(AL-transpose(conjg(AL)))
        !greater component
        call zgemm('n','n',nm,nm,nm,cone,iV,nm,FG,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,iV,nm,czero,yG_NN,nm)
        yG_NN=yG_NN/(1.0d0 - sum(E*conjg(E)))
        call zgemm('n','n',nm,nm,nm,cone,V,nm,yG_NN,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,V,nm,czero,wG_NN,nm)
        !refinement iteration
        call zgemm('n','n',nm,nm,nm,cone,xR,nm,M10,nm,czero,tmp1,nm)
        call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,wG_NN,nm,czero,tmp2,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp2,nm,M10,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,xR,nm,czero,tmp2,nm)
        wG_NN=FG+tmp2
        !
        call zgemm('n','n',nm,nm,nm,cone,M10,nm,wG_NN,nm,czero,tmp1,nm)
        call zgemm('n','c',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dLG11,nm)
        dLG11=dLG11-(AG-transpose(conjg(AG)))
    end subroutine get_dL_OBC_for_W


    FUNCTION eigv(NN, A)
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

    subroutine open_boundary_conditions(nm,M00,M01,M10,V10,xR,dM,dV,cond)
        integer,intent(in)::nm
        complex(8),intent(in),dimension(nm,nm)::M00,M01,M10,V10
        complex(8),intent(out),dimension(nm,nm)::xR,dM,dV
        real(8),intent(out)::cond
        complex(8),dimension(nm,nm)::tmp1
        call surface_function(nm,M00,M01,M10,xR,cond);
        !dM=M01*xR*M10
        call zgemm('n','n',nm,nm,nm,cone,M01,nm,xR,nm,czero,tmp1,nm)
        call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,M10,nm,czero,dM,nm)
        !dV=M01*xR*V10
        call zgemm('n','n',nm,nm,nm,cone,M01,nm,xR,nm,czero,tmp1,nm)
        call zgemm('n','n',nm,nm,nm,cone,tmp1,nm,V10,nm,czero,dV,nm)
    end subroutine open_boundary_conditions


    ! a slightly modified version of sancho
    subroutine surface_function(nm,M00,M01,M10,SF,cond)
        integer,intent(in)::nm
        complex(8),intent(in),dimension(nm,nm)::M00,M01,M10
        complex(8),intent(out),dimension(nm,nm)::SF
        real(8),intent(out)::cond
        real(8)::cond_limit
        integer::max_iteration,IC
        complex(8),dimension(:,:),allocatable::alpha,beta,Eps,Eps_surf,inv_element,a_i_b,b_i_a,i_alpha,i_beta
        allocate(alpha(nm,nm))
        allocate(beta(nm,nm))
        allocate(Eps(nm,nm))
        allocate(Eps_surf(nm,nm))
        allocate(inv_element(nm,nm))
        allocate(i_alpha(nm,nm))
        allocate(i_beta(nm,nm))
        allocate(a_i_b(nm,nm))
        allocate(b_i_a(nm,nm))
        cond=1.0d10;
        cond_limit=1.0d-10;
        max_iteration=5000;
        IC=1;
        alpha=M01
        beta=M10
        Eps=M00
        Eps_surf=M00
        do while ((cond>cond_limit).and.(IC<max_iteration))      
            inv_element=Eps
            call invert_inplace(inv_element,nm)
            i_alpha=matmul(inv_element,alpha)
            i_beta=matmul(inv_element,beta)
            a_i_b=matmul(alpha,i_beta)
            b_i_a=matmul(beta,i_alpha)
            Eps=Eps-a_i_b-b_i_a
            Eps_surf=Eps_surf-a_i_b
            alpha=matmul(alpha,i_alpha)
            beta=matmul(beta,i_beta)
            !
            cond=sum(abs(alpha)+abs(beta))/2.0d0;
            !
            IC=IC+1;
        end do
        if (cond>cond_limit) then 
        write(*,*) 'SEVERE warning: nmax reached in surface function!!!',cond
        endif
        call invert_inplace(Eps_surf,nm)
        SF=Eps_surf
        deallocate(alpha,beta,Eps,Eps_surf,inv_element,a_i_b,b_i_a,i_alpha,i_beta)
    end subroutine surface_function
    

end module open_boundary


