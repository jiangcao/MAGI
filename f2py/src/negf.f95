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


module fft_mod
    use parameters_mod
    implicit none 
    
    CONTAINS
    
    ! Z(nop) = sum_ie X(ie) Y(ie-nop)
    function corr1d(n,X,Y,method) result(Z)
    integer, intent(in)::n
    character(len=*),intent(in)::method
    complex(8),intent(in)::X(n),Y(n)
    complex(8)::Z(2*n-1)
    complex(8),allocatable,dimension(:) :: X_, Y_
    integer::i,ie
    select case (trim(method))
      case default ! explicit index
        Z=czero
        do i=-n+1,n-1
          do ie = max(i+1,1),min(n,n+i)
            Z(i+n)=Z(i+n) + X(ie)*Y(ie-i)
          enddo
        enddo
      case('simple')
        do i=-n+1,n-1
          Z(i+n)=sum(X(max(1,1+i):min(n+i,n))*Y(max(1-i,1):min(n,n-i)))
        enddo
      case('fft')  
        allocate(X_(n*2-1))
        allocate(Y_(n*2-1))
        X_=czero ! pad by zero
        Y_=czero
        X_(1:n) = X
        Y_(1:n) = Y(n:1:-1)
        
        call do_mkl_dfti_conv(n*2-1,X_,Y_,Z)
        
        deallocate(X_,Y_)
    end select
    end function corr1d
    
    
    ! Z(nop) = sum_ie X(ie) Y(ie-nop)
    function corr1d2(n,X,Y,method) result(Z)
    integer, intent(in)::n
    character(len=*),intent(in)::method
    complex(8),intent(in)::X(n),Y(n)
    complex(8)::Z(n)
    complex(8),allocatable,dimension(:) :: X_, Y_, Z_
    integer::i,ie
    select case (trim(method))
      case default ! explicit index
        Z=czero
        do i=-n/2+1,n/2-1
          do ie = max(i+1,1),min(n,n+i)
            Z(i+n/2)=Z(i+n/2) + X(ie)*Y(ie-i)
          enddo
        enddo
      case('simple')
        do i=-n/2+1,n/2-1
          Z(i+n/2)=sum(X(max(1,1+i):min(n+i,n))*Y(max(1-i,1):min(n,n-i)))
        enddo
      case('fft')  
        allocate(X_(n*2-1))
        allocate(Y_(n*2-1))
        allocate(Z_(n*2-1))
        X_=czero ! pad by zero
        Y_=czero
        X_(1:n) = X
        Y_(1:n) = Y(n:1:-1)
        
        call do_mkl_dfti_conv(n*2-1,X_,Y_,Z_)
        
        Z=Z_(n-n/2:n+n/2-1)
        deallocate(X_,Y_,Z_)
    end select
    end function corr1d2
    
    ! Z(ie) = sum_nop X(ie-nop) Y(nop)
    function conv1d(n,X,Y,method) result(Z)
    integer, intent(in)::n
    character(len=*),intent(in)::method
    complex(8),intent(in)::X(n),Y(n*2-1)
    complex(8)::Z(n)
    complex(8),allocatable,dimension(:)::x_in, y_in, z_in
    integer::i,ie
    select case (trim(method))
      case default ! explicit index
        Z=czero
        do ie=1,n
          do i= -n+1,n-1
            if ((ie .gt. max(i,1)).and.(ie .lt. min(n,(n+i)))) then
              Z(ie)=Z(ie) + X(ie-i)*Y(i+n)
            endif
          enddo
        enddo
      case('simple')
        do i=1,n
          Z(i)=sum(X(n:1:-1)*Y(i:i+n-1))
        enddo
      case('fft')
        allocate(X_in(n*2-1))
        allocate(Y_in(n*2-1))
        allocate(Z_in(n*2-1))
        X_in=czero
        X_in(1:n)=X    
        Y_in=cshift(Y,-n)    
        call do_mkl_dfti_conv(n*2-1,Y_in,X_in,Z_in)
        Z=Z_in(1:n)
        deallocate(X_in,Y_in,Z_in)
    end select
    end function conv1d
    
    
    ! Z(ie) = sum_nop X(ie-nop) Y(nop)
    function conv1d2(n,X,Y,method) result(Z)
    integer, intent(in)::n
    character(len=*),intent(in)::method
    complex(8),intent(in)::X(n),Y(n)
    complex(8)::Z(n)
    complex(8),allocatable,dimension(:)::x_in, y_in, z_in,tmp
    integer::i,ie,iop
    select case (trim(method))
      case default ! explicit index
        Z = czero  
        do i= -n/2+1,n/2-1  
          iop = i+n/2
          Z((max(i,1)+1):min(n,(n+i)))=Z(max(i,1)+1:min(n,(n+i))) + X((max(i,1)+1-i):(min(n,(n+i))-i))*Y(iop)
    !      do ie=max(i,1)+1,min(n,(n+i))      
    !        !if ((ie .gt. max(i,1)).and.(ie .lt. (n+i))) then
    !          Z(ie)=Z(ie) + X(ie-i)*Y(iop)
    !        !endif
    !      enddo
        enddo    
      case('simple')
        allocate(Y_in(n*2-1))    
        Y_in=czero
        Y_in(n/2+1:n/2+n-1)=Y(1:n-1)    
        do i=1,n
          Z(i) = sum( Y_in(i:i+n-1)*X(n:1:-1) )
        enddo
        deallocate(Y_in)
      case('fft')
        allocate(X_in(n*2-1))
        allocate(Y_in(n*2-1))
        allocate(Z_in(n*2-1))    
        X_in=czero
        Y_in=czero
        X_in(1:n)=X    
        Y_in(1:n/2)=Y(n/2:n-1)
        Y_in(n*2-n/2+1:n*2-1)=Y(1:n/2-1)
         
        call do_mkl_dfti_conv(n*2-1,Y_in,X_in,Z_in)
        
        Z=Z_in(1:n)
        deallocate(X_in,Y_in,Z_in)
    end select
    end function conv1d2
    
    subroutine do_mkl_dfti_conv(n,X_in,Y_in,Z_out)
    ! 1D complex to complex
    Use MKL_DFTI
    integer :: n
    Complex(8) :: X_in(n),Y_in(n),Z_out(n)
    Complex(8) :: X_out(n),Y_out(n),Z_in(n)
    type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle, My_Desc2_Handle
    Integer :: Status
    ! Perform a complex to complex transform
    Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n )
    Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
    Status = DftiCommitDescriptor( My_Desc1_Handle )
    Status = DftiComputeForward( My_Desc1_Handle, X_in, X_out )
    Status = DftiComputeForward( My_Desc1_Handle, Y_in, Y_out )
    !
    Z_in(:) = X_out(:) * Y_out(:)
    !
    Status = DftiComputeBackward( My_Desc1_Handle, Z_in, Z_out )
    Status = DftiFreeDescriptor(My_Desc1_Handle)
    Z_out(:) = Z_out(:) / dble(n)
    end subroutine do_mkl_dfti_conv
    
end module fft_mod
    





module gf_dense 
    use parameters_mod
    implicit none 

    contains


    ! 3D GW solver with two periodic directions (y,z)
    ! iterating G -> P -> W -> Sig 
    subroutine solve_gw_3D(niter,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
        alpha_mix,nen,En,nb,ns,nphiy,nphiz,Ham,H00lead,H10lead,T,V,&
        G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,&
        W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater,&
        Sig_retarded_new,Sig_lesser_new,Sig_greater_new,ldiag,flatband)
    !
    use fft_mod, only : conv1d => conv1d2, corr1d => corr1d2  
    use parameters_mod
    !  
    integer, intent(in) :: nen, nb, ns,niter,nm_dev,length, nphiz, nphiy
    real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg
    complex(8),intent(in) :: Ham(nm_dev,nm_dev,nphiy*nphiz),H00lead(NB*NS,NB*NS,2,nphiy*nphiz),H10lead(NB*NS,NB*NS,2,nphiy*nphiz),T(NB*NS,nm_dev,2,nphiy*nphiz)
    complex(8), intent(in):: V(nm_dev,nm_dev,nphiy*nphiz)
    logical,intent(in)::ldiag
    logical,intent(in)::flatband
    complex(8),intent(out),dimension(nm_dev,nm_dev,nen,nphiy*nphiz) ::  G_retarded,G_lesser,G_greater
    complex(8),intent(out),dimension(nm_dev,nm_dev,nen,nphiy*nphiz) ::  P_retarded,P_lesser,P_greater
    complex(8),intent(out),dimension(nm_dev,nm_dev,nen,nphiy*nphiz) ::  W_retarded,W_lesser,W_greater
    complex(8),intent(out),dimension(nm_dev,nm_dev,nen,nphiy*nphiz) ::  Sig_retarded,Sig_lesser,Sig_greater
    complex(8),intent(out),dimension(nm_dev,nm_dev,nen,nphiy*nphiz) ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
    !------
    complex(8),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
    complex(8),allocatable,dimension(:,:):: B ! tmp matrix
    real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:),sumcur(:,:,:),sumtot_cur(:,:),sumtot_ecur(:,:)
    complex(8),allocatable::Ispec(:,:,:),Itot(:,:)
    real(8),allocatable::Tr(:,:) ! current spectrum on leads
    real(8),allocatable::Te(:,:,:) ! transmission matrix spectrum
    real(8),allocatable::sumTr(:,:) ! current spectrum on leads summed over k
    real(8),allocatable::sumTe(:,:,:) ! transmission matrix spectrum summed over k
    integer :: iter,ie,nopmax
    integer :: i,j,nm,nop,l,h,iop,ndiag,ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd        
    complex(8) :: dE
    real(8)::nelec(2),mu(2),pelec(2),temp(2)
    !
    print *,'============ green_solve_gw_3D ============'
    allocate(siglead(NB*NS,NB*NS,nen,2,nphiy*nphiz))
    ! get leads sigma
    do ikz=1, nphiy*nphiz
        siglead(:,:,:,1,ikz) = Sig_retarded(1:NB*NS,1:NB*NS,:,ikz)
        siglead(:,:,:,2,ikz) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:,ikz)  
    enddo
    allocate(B(nm_dev,nm_dev))
    allocate(tot_cur(nm_dev,nm_dev))
    allocate(tot_ecur(nm_dev,nm_dev))
    allocate(sumtot_cur(nm_dev,nm_dev))
    allocate(sumtot_ecur(nm_dev,nm_dev))
    allocate(cur(nm_dev,nm_dev,nen))
    allocate(sumcur(nm_dev,nm_dev,nen))
    allocate(Ispec(nm_dev,nm_dev,nen))
    allocate(Itot(nm_dev,nm_dev))
    allocate(tr(nen,2))
    allocate(te(nen,2,2))
    allocate(sumtr(nen,2))
    allocate(sumte(nen,2,2))
    if (flatband) then
        mu=(mus+mud)/2.0d0
        temp=(temps+tempd)/2.0d0
    else
        mu=(/ mus, mud /)
        temp=(/temps,tempd/)
    endif
    
    print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
    do iter=0,niter  
        print *,'+ iter=',iter
        print *, 'calc G'  
        sumtot_cur=0.0d0
        sumtot_ecur=0.0d0
        sumcur=0.0d0
        do ikz=1,nphiy*nphiz
        !  print *, ' ik=', ikz,'/',nphiy*nphiz
        call calc_gf(nen,En,2,nm_dev,(/nb*ns,nb*ns/),nb*ns,Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),&
            T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),G_retarded(:,:,:,ikz),G_lesser(:,:,:,ikz),&
            G_greater(:,:,:,ikz),Tr,Te,mu,temp,flatband)
        !call write_spectrum('ldos_kz'//string(ikz)//'_',iter,G_retarded(:,:,:,ikz),nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
        call calc_bond_current(Ham(:,:,ikz),G_lesser(:,:,:,ikz),nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)    
        !call write_current_spectrum('Jdens_kz'//string(ikz)//'_',iter,cur,nen,en,length,NB,Lx)    
        sumcur=sumcur+cur
        sumtot_cur=sumtot_cur+tot_cur
        sumtot_ecur=sumtot_ecur+tot_ecur
        sumTr=sumTr+Tr
        sumTe=sumTe+Te
        enddo
        sumcur=sumcur/dble(nphiy)/dble(nphiz)
        sumtot_cur=sumtot_cur/dble(nphiy)/dble(nphiz)
        sumtot_ecur=sumtot_ecur/dble(nphiy)/dble(nphiz)
        sumTr=sumTr/dble(nphiz)/dble(nphiy)
        sumTe=sumTe/dble(nphiz)/dble(nphiy)
        if (flatband) then
            print *,'flatband'
            ! call write_spectrum_per_kz('gw_ldos',iter,G_retarded,nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
            ! call write_spectrum_per_kz('gw_gamma-centered_ldos',iter,G_retarded,nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=0.0_dp,at_kz=0.0_dp)
            call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
    
        else
            call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
            call write_spectrum_summed_over_kz('gw_ndos',iter,G_lesser,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            call write_spectrum_summed_over_kz('gw_pdos',iter,G_greater,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
        endif
        ! call write_current_spectrum('gw_Jdens',iter,sumcur,nen,en,length,NB,Lx)
        ! call write_current('gw_I',iter,sumtot_cur,length,NB,NS,Lx)
        ! call write_current('gw_EI',iter,sumtot_ecur,length,NB,NS,Lx)
        ! call write_transmission_spectrum('gw_trL',iter,sumTr(:,1)*spindeg,nen,En)
        ! call write_transmission_spectrum('gw_trR',iter,sumTr(:,2)*spindeg,nen,En)
        ! call write_transmission_spectrum('gw_TE_LR',iter,sumTe(:,1,2)*spindeg,nen,En)
        ! call write_transmission_spectrum('gw_TE_RL',iter,sumTe(:,2,1)*spindeg,nen,En)
        open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
        write(101,'(I4,2E16.6)') iter, sum(sumTr(:,1))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg), &
                                    sum(sumTr(:,2))*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)
        close(101)
        !
        G_retarded=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
        G_lesser=dcmplx(0.0d0*dble(G_lesser),aimag(G_lesser))
        G_greater=dcmplx(0.0d0*dble(G_greater),aimag(G_greater))
        !        
        print *, 'calc P'
        !
        nopmax=nen/2-1  
        ndiag=nm_dev !NB*(min(NS*2,iter))
        if (ldiag) ndiag=0  
        print *,'ndiag=',ndiag
        ! Pij^<>(hw,kz') = \int_dE Gij^<>(E,kz) * Gji^><(E-hw,kz-kz')
        ! Pij^r(hw,kz')  = \int_dE Gij^<(E,kz) * Gji^a(E-hw,kz-kz') + Gij^r(E,kz) * Gji^<(E-hw,kz-kz')
        !$omp parallel default(shared) private(l,h,ikz,ikzd,iqz,dE,i,j,iky,ikyd,iqy,ik,iq,ikd) 
        !$omp do
        !do nop=-nopmax,nopmax
        do iq=1,nphiy*nphiz
    !    do iqy=1,nphiy
    !      do iqz=1,nphiz
            !iop=nop+nen/2
    !        iq=iqz+(iqy-1)*nphiz
            iqz = mod(iq-1,nphiz)+1
            iqy = (iq-1) / nphiz +1
            P_lesser(:,:,:,iq) = dcmplx(0.0d0,0.0d0)
            P_greater(:,:,:,iq) = dcmplx(0.0d0,0.0d0)    
            P_retarded(:,:,:,iq) = dcmplx(0.0d0,0.0d0)    
            ! do ie = max(nop+1,1),min(nen,nen+nop) 
                do iky=1,nphiy
                do ikz=1,nphiz              
                    ik=ikz + (iky-1)*nphiz
                    do i = 1, nm_dev        
                    l=max(i-ndiag,1)
                    h=min(nm_dev,i+ndiag)
                    ikzd=ikz-iqz + nphiz/2
                    ikyd=iky-iqy + nphiy/2
                    if (ikzd<1) ikzd=ikzd+nphiz
                    if (ikzd>nphiz) ikzd=ikzd-nphiz
                    if (ikyd<1) ikyd=ikyd+nphiy
                    if (ikyd>nphiy) ikyd=ikyd-nphiy   
                    if (nphiy==1)   ikyd=1
                    if (nphiz==1)   ikzd=1             
                    ikd=ikzd + (ikyd-1)*nphiz
                    do j = l,h
                        P_lesser(i,j,:,iq) = P_lesser(i,j,:,iq) + corr1d(nen,G_lesser(i,j,:,ik),G_greater(j,i,:,ikd),method='fft') 
                        P_greater(i,j,:,iq) = P_greater(i,j,:,iq) + corr1d(nen,G_greater(i,j,:,ik),G_lesser(j,i,:,ikd),method='fft')         
                        P_retarded(i,j,:,iq) = P_retarded(i,j,:,iq) + corr1d(nen,G_lesser(i,j,:,ik),conjg(G_retarded(i,j,:,ikd)),method='fft') &
                                                                    + corr1d(nen,G_retarded(i,j,:,ik),G_lesser(j,i,:,ikd),method='fft') 
                    enddo
                    enddo
                enddo
                enddo
            ! enddo
        ! enddo
        enddo
        !enddo
        !$omp end do
        !$omp end parallel
        dE = dcmplx(0.0d0 , -1.0d0*( En(2) - En(1) ) / 2.0d0 / pi ) * spindeg /dble(nphiy)/dble(nphiz) 
        P_lesser=dE*P_lesser
        P_greater=dE*P_greater
        P_retarded=dE*P_retarded
        if (flatband) then
            print *,'flatband'
            ! call write_spectrum_per_kz('PR',iter,P_retarded,nen,En-en(nen/2),nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
        endif 
        ! call write_spectrum_summed_over_kz('PR',iter,P_retarded,nen,En-en(nen/2),nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
    !  call write_spectrum_summed_over_kz('PL',iter,P_lesser  ,nen,En-en(nen/2),nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
    !  call write_spectrum_summed_over_kz('PG',iter,P_greater ,nen,En-en(nen/2),nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
        !
        print *, 'calc W'
        !
        do iq=1,nphiy*nphiz
        !  print *, ' iq=', iq,'/',nphiy*nphiz
        !$omp parallel default(shared) private(nop)
        !$omp do
        do nop=-nopmax+nen/2,nopmax+nen/2       
            if (flatband) then
                ! call green_calc_w(0,NB,NS,nm_dev,P_retarded(:,:,nop,iq),P_lesser(:,:,nop,iq),P_greater(:,:,nop,iq),V(:,:,iq),W_retarded(:,:,nop,iq),W_lesser(:,:,nop,iq),W_greater(:,:,nop,iq))
            else      
                ! call green_calc_w(1,NB,NS,nm_dev,P_retarded(:,:,nop,iq),P_lesser(:,:,nop,iq),P_greater(:,:,nop,iq),V(:,:,iq),W_retarded(:,:,nop,iq),W_lesser(:,:,nop,iq),W_greater(:,:,nop,iq))
            endif
        enddo
        !$omp end do
        !$omp end parallel
        enddo  
        if (flatband) then
            print *,'flatband'
            ! call write_spectrum_per_kz('WR',iter,W_retarded,nen,En-en(nen/2),nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            if (iter==0) then
                ! call write_W_per_kz('WR',iter,W_retarded,nen,En-en(nen/2),nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/),V)
            else
                ! call write_W_per_kz('WR',iter,W_retarded,nen,En-en(nen/2),nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            endif
        endif 
        ! call write_spectrum_summed_over_kz('WR',iter,W_retarded,nen,En-en(nen/2),nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
    !  call write_spectrum_summed_over_kz('WL',iter,W_lesser,  nen,En-en(nen/2),nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
    !  call write_spectrum_summed_over_kz('WG',iter,W_greater, nen,En-en(nen/2),nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
        !
        print *, 'calc SigGW'
        !  
        ndiag=nm_dev
        if (ldiag) ndiag=0  
        print *,'ndiag=',ndiag
        nopmax=nen/2-1
        Sig_greater_new = dcmplx(0.0d0,0.0d0)
        Sig_lesser_new = dcmplx(0.0d0,0.0d0)
        Sig_retarded_new = dcmplx(0.0d0,0.0d0)      
        ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)
        !$omp parallel default(shared) private(l,h,i,j,ikz,ikzd,iqz,iky,ikyd,iqy,ik,iq,ikd)
        !$omp do  
        do ik=1,nphiy*nphiz
    !  do iky=1,nphiy
    !    do ikz=1,nphiz      
    !      ik=ikz+(iky-1)*nphiz
            ikz = mod(ik-1,nphiz)+1
            iky = (ik-1) / nphiz +1
            do iqy=1,nphiy
            do iqz=1,nphiz
                iq=iqz+(iqy-1)*nphiz         
                ikzd=ikz-iqz + nphiz/2            
                if (ikzd<1) ikzd=ikzd+nphiz
                if (ikzd>nphiz) ikzd=ikzd-nphiz
                ikyd=iky-iqy + nphiy/2            
                if (ikyd<1) ikyd=ikyd+nphiy
                if (ikyd>nphiy) ikyd=ikyd-nphiy
                if (nphiy==1)   ikyd=1
                if (nphiz==1)   ikzd=1
                ikd=ikzd+(ikyd-1)*nphiz
                do i = 1,nm_dev   
                l=max(i-ndiag,1)
                h=min(nm_dev,i+ndiag)       
                do j = l,h
                    Sig_lesser_new(i,j,:,ik)=Sig_lesser_new(i,j,:,ik) + conv1d(nen,G_lesser(i,j,:,ikd),W_lesser(i,j,:,iq),method='fft') 
                    Sig_greater_new(i,j,:,ik)=Sig_greater_new(i,j,:,ik) + conv1d(nen,G_greater(i,j,:,ikd),W_greater(i,j,:,iq),method='fft') 
                    Sig_retarded_new(i,j,:,ik)=Sig_retarded_new(i,j,:,ik) &
                                                + conv1d(nen,G_lesser(i,j,:,ikd),W_retarded(i,j,:,iq),method='fft') &
                                                + conv1d(nen,G_retarded(i,j,:,ikd),W_lesser(i,j,:,iq),method='fft') &
                                                + conv1d(nen,G_retarded(i,j,:,ikd),W_retarded(i,j,:,iq),method='fft')                                               
                enddo
                enddo      
            enddo
            enddo        
    !    enddo
        enddo
        !$omp end do
        !$omp end parallel
        dE = dcmplx(0.0d0, (En(2)-En(1))/2.0d0/pi) /dble(nphiy)/dble(nphiz)
        Sig_lesser_new = Sig_lesser_new  * dE
        Sig_greater_new= Sig_greater_new * dE
        Sig_retarded_new=Sig_retarded_new* dE
        Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
        !!! Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
        !!! Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
        !
        ! symmetrize the selfenergies
        do ie=1,nen
        do ikz=1,nphiy*nphiz
            B(:,:)=transpose(Sig_retarded_new(:,:,ie,ikz))
            Sig_retarded_new(:,:,ie,ikz) = (Sig_retarded_new(:,:,ie,ikz) + B(:,:))/2.0d0    
            B(:,:)=transpose(Sig_lesser_new(:,:,ie,ikz))
            Sig_lesser_new(:,:,ie,ikz) = (Sig_lesser_new(:,:,ie,ikz) + B(:,:))/2.0d0
            B(:,:)=transpose(Sig_greater_new(:,:,ie,ikz))
            Sig_greater_new(:,:,ie,ikz) = (Sig_greater_new(:,:,ie,ikz) + B(:,:))/2.0d0
        enddo
        enddo
        ! mixing with previous ones
        Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
        Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
        Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
        if (.not. flatband) then
        ! get leads sigma
        do iqz=1,nphiy*nphiz
            siglead(:,:,:,1,iqz) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:,iqz)
            siglead(:,:,:,2,iqz) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:,iqz)    
        enddo
        if (length>3) then
            ! make sure self-energy is continuous near leads (by copying edge block)
            do ie=1,nen
                do iqz=1,nphiy*nphiz
                    call expand_size_bycopy(Sig_retarded(:,:,ie,iqz),nm_dev,NB,2)
                    call expand_size_bycopy(Sig_lesser(:,:,ie,iqz),nm_dev,NB,2)
                    call expand_size_bycopy(Sig_greater(:,:,ie,iqz),nm_dev,NB,2)
                enddo
            enddo
        endif
        endif
        if (flatband) then
            print *,'flatband'
            ! call write_spectrum_per_kz('gw_SigR',iter,Sig_retarded,nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
        endif 
        ! call write_spectrum_summed_over_kz('gw_SigR',iter,Sig_retarded,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
    !  call write_spectrum_summed_over_kz('SigL',iter,Sig_lesser,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
    !  call write_spectrum_summed_over_kz('SigG',iter,Sig_greater,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
    end do  
    deallocate(siglead,B)
    deallocate(sumcur,cur,tot_cur,tot_ecur,sumtot_cur,sumtot_ecur)
    deallocate(Ispec,Itot)
    deallocate(Tr,Te,sumTr,sumTe)
    end subroutine solve_gw_3D
  
  

    Function ferm(a)
        Real(8),intent(in):: a
        Real(8):: ferm
        ferm=1.0d0/(1.0d0+Exp(a))
    End Function ferm

    ! calculate Gr and G<>
    subroutine calc_gf(ne,E,num_lead,nm_dev,nm_lead,max_nm_lead,Ham,lead_H00,lead_H10,Siglead,T,&
        Scat_Sig_retarded,Scat_Sig_lesser,Scat_Sig_greater,G_retarded,G_lesser,G_greater,&
        cur,te,mu,temp,flatband)
        use parameters_mod
        implicit none
        integer, intent(in) :: num_lead ! number of leads/contacts
        integer, intent(in) :: nm_dev   ! size of device Hamiltonian
        integer, intent(in) :: nm_lead(:) ! size of lead Hamiltonians
        integer, intent(in) :: max_nm_lead ! max size of lead Hamiltonians
        real(kind=dp), intent(in) :: E(:)  ! energy vector
        real(kind=dp), intent(out) :: cur(ne,num_lead)  ! current spectrum on leads
        real(kind=dp), intent(out) :: te(ne,num_lead,num_lead)  ! transmission matrix
        integer, intent(in) :: ne ! number of energy points
        complex(kind=dp), intent(in) :: Ham(nm_dev,nm_dev) ! system Hamiltonian
        complex(kind=dp), intent(in) :: lead_H00(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian diagonal blocks
        complex(kind=dp), intent(in) :: lead_H10(max_nm_lead,max_nm_lead,num_lead) ! lead Hamiltonian off-diagonal blocks
        complex(kind=dp), intent(in) :: Siglead(max_nm_lead,max_nm_lead,ne,num_lead) ! lead sigma_r scattering
        complex(kind=dp), intent(in) :: T(max_nm_lead,nm_dev,num_lead)  ! coupling matrix between leads and device
        complex(kind=dp), intent(in) :: Scat_Sig_retarded(nm_dev,nm_dev,ne) ! scattering Selfenergy
        complex(kind=dp), intent(in) :: Scat_Sig_lesser(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(in) :: Scat_Sig_greater(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(out) :: G_retarded(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(out) :: G_lesser(nm_dev,nm_dev,ne)
        complex(kind=dp), intent(out) :: G_greater(nm_dev,nm_dev,ne)
        logical,intent(in)::flatband
        real(kind=dp), intent(in) :: mu(num_lead), temp(num_lead)
        ! ----
        integer :: i,j,nm,ie,io
        complex(kind=dp), allocatable, dimension(:,:) :: S00,G00,GBB,A,sig,sig_lesser,sig_greater,B,C,Hii,H1i
        complex(kind=dp), allocatable, dimension(:,:,:) :: gamma_lead
        real(kind=dp) :: fd
        complex(kind=dp):: z     
        !
        cur=0.0d0
        te=0.0d0
        !
        !$omp parallel default(shared) private(z,sig,ie,sig_lesser,sig_greater,B,C,i,nm,S00,G00,GBB,A,fd,gamma_lead,Hii,H1i) 
        allocate(sig(nm_dev,nm_dev))  
        allocate(gamma_lead(nm_dev,nm_dev,num_lead))  
        allocate(sig_lesser(nm_dev,nm_dev))
        allocate(sig_greater(nm_dev,nm_dev))          
        allocate(B(nm_dev,nm_dev))
        allocate(C(nm_dev,nm_dev))
        !$omp do        
        do ie = 1, ne
            z=dcmplx(E(ie),0.0d0)
            G_retarded(:,:,ie) = - Ham(:,:) - Scat_Sig_retarded(:,:,ie)             
            sig_lesser(:,:) = dcmplx(0.0d0,0.0d0)      
            sig_greater(:,:) = dcmplx(0.0d0,0.0d0)      
            ! compute and add contact self-energies    
            !  open(unit=101,file='sancho_gbb.dat',status='unknown',position='append')
            !  open(unit=102,file='sancho_g00.dat',status='unknown',position='append')
            !  open(unit=103,file='sancho_sig.dat',status='unknown',position='append')
            do i = 1,num_lead
                nm = nm_lead(i)    
                allocate(Hii(nm,nm))
                allocate(H1i(nm,nm))
                allocate(S00(nm,nm))
                allocate(G00(nm,nm))
                allocate(GBB(nm,nm))
                allocate(A(nm_dev,nm))    
                call identity(S00,nm)            
                !
                if (flatband) then
                    G00 = -S00*c1i
                else        
                    Hii = lead_H00(1:nm,1:nm,i) + siglead(1:nm,1:nm,ie,i)
                    H1i = lead_H10(1:nm,1:nm,i)
                    call sancho(NM,E(ie),S00,Hii,H1i,G00,GBB)
                endif
                !
                call zgemm('c','n',nm_dev,nm,nm,cone,T(1:nm,1:nm_dev,i),nm,G00,nm,czero,A,nm_dev) 
                call zgemm('n','n',nm_dev,nm_dev,nm,cone,A,nm_dev,T(1:nm,1:nm_dev,i),nm,czero,sig,nm_dev)  
                !    write(101,'(i4,2E15.4)') i, E(ie), -aimag(trace(GBB,nm))*2.0d0
                !    write(102,'(i4,2E15.4)') i, E(ie), -aimag(trace(G00,nm))*2.0d0
                !    write(103,'(i4,2E15.4)') i, E(ie), -aimag(trace(sig,nm_dev))*2.0d0
                G_retarded(:,:,ie) = G_retarded(:,:,ie) - sig(:,:)
                !
                fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))
                B(:,:) = conjg(sig(:,:))
                C(:,:) = transpose(B(:,:))
                B(:,:) = sig(:,:) - C(:,:)
                sig_lesser(:,:) = sig_lesser(:,:) - B(:,:)*fd
                sig_greater(:,:) = sig_greater(:,:) + B(:,:)*(1.0d0-fd)            
                gamma_lead(:,:,i)= B(:,:)             
                deallocate(S00,G00,GBB,A,Hii,H1i)
            end do  
            !  close(101)
            !  close(102)
            !  close(103)
            do i = 1,nm_dev
                G_retarded(i,i,ie) = G_retarded(i,i,ie) + z 
            end do
            !
            call invert_inplace(G_retarded(:,:,ie),nm_dev) 
            !
            sig_lesser = sig_lesser + Scat_Sig_lesser(:,:,ie)
            sig_greater = sig_greater + Scat_Sig_greater(:,:,ie)     
            !
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_lesser,nm_dev,czero,B,nm_dev) 
            call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
            G_lesser(:,:,ie) = C   
            !
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,G_retarded(:,:,ie),nm_dev,sig_greater,nm_dev,czero,B,nm_dev) 
            call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,G_retarded(:,:,ie),nm_dev,czero,C,nm_dev)
            G_greater(:,:,ie) = C         
            !
            ! calculate current spec and/or transmission at each lead/contact
            do i=1,num_lead                                  
                fd = ferm((E(ie)-mu(i))/(BOLTZ*TEMP(i)))
                call zgemm('n','n',nm_dev,nm_dev,nm_dev,dcmplx(1.0d0-fd,0.0d0),&
                    gamma_lead(:,:,i),nm_dev,G_lesser(:,:,ie),nm_dev,czero,B,nm_dev)
                call zgemm('n','n',nm_dev,nm_dev,nm_dev,dcmplx(fd,0.0d0),&
                    gamma_lead(:,:,i),nm_dev,G_greater(:,:,ie),nm_dev,cone,B,nm_dev)
                do io=1,nm_dev
                    cur(ie,i)=cur(ie,i)+ dble(B(io,io))
                enddo            
                do j=1,num_lead                      
                    if (j.ne.i) then
                        call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,gamma_lead(:,:,i),nm_dev,G_retarded(:,:,ie),nm_dev,czero,B,nm_dev)
                        call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,gamma_lead(:,:,j),nm_dev,czero,C,nm_dev)
                        call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,C,nm_dev,G_retarded(:,:,ie),nm_dev,czero,B,nm_dev)
                        do io=1,nm_dev
                            te(ie,i,j)=te(ie,i,j)- dble(B(io,io)) ! Gamma = i[Sig^r - Sig^r \dagger] , hence the -1
                        enddo
                    endif
                enddo
            enddo
        end do  
        !$omp end do
        deallocate(sig,B,C,sig_lesser,sig_greater,gamma_lead)
        !$omp end parallel
    end subroutine calc_gf


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



    ! calculate bond current using I_ij = H_ij G<_ji - H_ji G^<_ij
    subroutine calc_bond_current(H,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
        complex(8),intent(in)::H(nm_dev,nm_dev),G_lesser(nm_dev,nm_dev,nen)
        real(8),intent(in)::en(nen),spindeg
        integer,intent(in)::nen,nm_dev ! number of E and device dimension
        real(8),intent(out)::tot_cur(nm_dev,nm_dev) ! total bond current density
        real(8),intent(out)::tot_ecur(nm_dev,nm_dev) ! total bond energy current density
        real(8),intent(out)::cur(nm_dev,nm_dev,nen) ! energy resolved bond current density
        !----
        complex(8),allocatable::B(:,:)
        integer::ie,io,jo
        real(8),parameter::tpi=6.28318530718  
          allocate(B(nm_dev,nm_dev))
          tot_cur=0.0d0  
          tot_ecur=0.0d0
          do ie=1,nen
            do io=1,nm_dev
              do jo=1,nm_dev
                B(io,jo)=H(io,jo)*G_lesser(jo,io,ie) - H(jo,io)*G_lesser(io,jo,ie)
              enddo
            enddo    
            B=B*(En(2)-En(1))*e_charge/twopi/hbar*e_charge*dble(spindeg)
            cur(:,:,ie) = dble(B)
            tot_ecur=tot_ecur+ en(ie)*dble(B)
            tot_cur=tot_cur+ dble(B)          
          enddo
          deallocate(B)
    end subroutine calc_bond_current

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


    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_per_kz(dataset,i,G,nen,en,nky,nkz,length,NB,Lx,coeff,at_ky,at_kz)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:,:) ! (m,m,e,k) kz is the fast-running index in k
        integer, intent(in)::i,nen,length,NB,nky,nkz
        real(8), intent(in)::Lx,en(nen),coeff(2)
        real(8), intent(in), optional::at_ky,at_kz
        integer:: ie,j,ib,ikz,iky,k,kb
        real(8):: kcenter(2),dky,dkz
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        kcenter=0.0d0
        if (present(at_ky)) then
            kcenter(1)=at_ky
        endif
        if (present(at_kz)) then
            kcenter(2)=at_kz
        endif
        dky=1.0d0/dble(nky)
        dkz=1.0d0/dble(nkz)
        open(unit=11,file=trim(dataset)//i_str//'_ky.dat',status='unknown')
        do ie = 1,nen
            ikz=max(min(nkz/2+ floor(kcenter(2)/dkz)+1 , nkz) , 1)            
            do iky=1,nky
                tr=0.0d0          
                do j = 1,length
                    do ib=1,nb
                        tr = tr+G((j-1)*nb+ib,(j-1)*nb+ib,ie,ikz+nkz*(iky-1))            
                    enddo
                end do
                write(11,'(4E18.4)') dble(iky)/dble(nky), en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            enddo
            write(11,*)    
        enddo
        close(11)

        open(unit=11,file=trim(dataset)//i_str//'_kz.dat',status='unknown')
        do ie = 1,nen
            iky=max(min(nky/2+ floor(kcenter(1)/dky)+1 , nkz) , 1)           
            do ikz=1,nkz
                tr=0.0d0          
                do j = 1,length
                    do ib=1,nb
                        tr = tr+G((j-1)*nb+ib,(j-1)*nb+ib,ie,ikz+nkz*(iky-1))            
                    enddo
                end do
                write(11,'(4E18.4)') dble(ikz)/dble(nkz), en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            enddo
            write(11,*)    
        enddo
        close(11)
    end subroutine write_spectrum_per_kz


    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_summed_over_kz(dataset,i,G,nen,en,nkz,length,NB,Lx,coeff)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:,:)
        integer, intent(in)::i,nen,length,NB,nkz
        real(8), intent(in)::Lx,en(nen),coeff(2)
        integer:: ie,j,ib,ikz
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do j = 1,length
                tr=0.0d0          
                do ib=1,nb
                do ikz=1,nkz
                    tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie,ikz)            
                enddo
                enddo
                tr=tr/dble(nkz)
                write(11,'(4E18.4)') dble(j-1)*Lx, en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            end do
            write(11,*)    
        enddo
        close(11)
    end subroutine write_spectrum_summed_over_kz



end module gf_dense

