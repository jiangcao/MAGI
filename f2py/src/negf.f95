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
    


module output 
    implicit none
    contains 

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
                write(11,'(4E20.6)') dble(iky)/dble(nky), en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
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
                write(11,'(4E20.6)') dble(ikz)/dble(nkz), en(ie), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
            enddo
            write(11,*)    
        enddo
        close(11)
    end subroutine write_spectrum_per_kz


    ! write spectrum into file (pm3d map)
    subroutine write_spectrum_summed_over_kz(dataset,i,G,nen,nsub,en,ensub,nkz,length,NB,Lx,coeff)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:,:,:)
        integer, intent(in)::i,nen,length,NB,nkz,nsub
        real(8), intent(in)::Lx,en(nen),coeff(2),ensub(nsub)
        integer:: ie,j,ib,ikz,isub
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do isub = 1,nsub
                do j = 1,length
                    tr=0.0d0          
                    do ib=1,nb
                    do ikz=1,nkz
                        tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie,isub,ikz)            
                    enddo
                    enddo
                    tr=tr/dble(nkz)
                    write(11,'(4E20.6)') dble(j-1)*Lx, en(ie)+ensub(isub), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
                end do
                write(11,*)   
            enddo 
        enddo
        close(11)
    end subroutine write_spectrum_summed_over_kz


    ! write spectrum into file (pm3d map)
    subroutine write_spectrum(dataset,i,G,nen,nsub,en,ensub,length,NB,Lx,coeff)
        character(len=*), intent(in) :: dataset
        complex(8), intent(in) :: G(:,:,:,:)
        integer, intent(in)::i,nen,length,NB,nsub
        real(8), intent(in)::Lx,en(nen),coeff(2),ensub(nsub)
        integer:: ie,j,ib,isub
        complex(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do isub = 1,nsub
                do j = 1,length
                    tr=0.0d0          
                    do ib=1,nb                    
                        tr = tr+ G((j-1)*nb+ib,(j-1)*nb+ib,ie,isub)                                
                    enddo                    
                    write(11,'(4E20.6)') dble(j-1)*Lx, en(ie)+ensub(isub), dble(tr)*coeff(1), aimag(tr)*coeff(2)        
                end do
                write(11,*)   
            enddo 
        enddo
        close(11)
    end subroutine write_spectrum


    ! write current into file 
    subroutine write_current(dataset,i,cur,length,NB,NS,Lx)
        character(len=*), intent(in) :: dataset
        real(8), intent(in) :: cur(:,:)
        integer, intent(in)::i,length,NB,NS
        real(8), intent(in)::Lx
        integer:: j,ib,jb,ii
        real(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ii = 1,length-1
            tr=0.0d0          
            do ib=1,nb  
            do jb=1,nb       
                do j=ii,min(ii+NS-1,length-1)
                tr = tr+ cur((ii-1)*nb+ib,j*nb+jb)
                enddo
            enddo                        
            end do
            write(11,'(2E20.6)') dble(ii)*Lx, tr
        end do
    end subroutine write_current
    
    ! write current spectrum into file (pm3d map)
    subroutine write_current_spectrum(dataset,i,cur,nen,en,length,NB,Lx)
        character(len=*), intent(in) :: dataset
        real(8), intent(in) :: cur(:,:,:)
        integer, intent(in)::i,nen,length,NB
        real(8), intent(in)::Lx,en(nen)
        integer:: ie,j,ib,jb
        real(8)::tr
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen
            do j = 1,length-1
                tr=0.0d0          
                do ib=1,nb  
                do jb=1,nb        
                    tr = tr+ cur((j-1)*nb+ib,j*nb+jb,ie)
                enddo                        
                end do
                write(11,'(3E20.6)') dble(j)*Lx, en(ie), tr
            end do
            write(11,*)    
        end do
        close(11)
    end subroutine write_current_spectrum

    ! write transmission spectrum into file
    subroutine write_transmission_spectrum(dataset,i,tr,nen,en)
        character(len=*), intent(in) :: dataset
        real(8), intent(in) :: tr(:)
        integer, intent(in)::i,nen
        real(8), intent(in)::en(nen)
        integer:: ie,j,ib
        character(len=4) :: i_str
        character(len=8) :: fmt
        fmt = '(I4.4)'
        write (i_str, fmt) i 
        open(unit=11,file=trim(dataset)//i_str//'.dat',status='unknown')
        do ie = 1,nen    
        write(11,'(2E20.6)') en(ie), dble(tr(ie))      
        end do
        close(11)
    end subroutine write_transmission_spectrum

end module output


module legendre
    contains 
    ! gauleg.f90     P145 Numerical Recipes in Fortran
    ! compute x(i) and w(i)  i=1,n  Legendre ordinates and weights
    ! on interval -1.0 to 1.0 (length is 2.0)
    ! use ordinates and weights for Gauss Legendre integration
    !
    subroutine gaulegf(x1, x2, x, w, n)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: x1, x2
        double precision, dimension(n), intent(out) :: x, w
        integer :: i, j, m
        double precision :: p1, p2, p3, pp, xl, xm, z, z1
        double precision, parameter :: eps=3.d-14
            
        m = (n+1)/2
        xm = 0.5d0*(x2+x1)
        xl = 0.5d0*(x2-x1)
        do i=1,m
        z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
        z1 = 0.0
        do while(abs(z-z1) .gt. eps)
            p1 = 1.0d0
            p2 = 0.0d0
            do j=1,n
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
            end do
            pp = n*(z*p1-p2)/(z*z-1.0d0)
            z1 = z
            z = z1 - p1/pp
        end do
        x(i) = xm - xl*z
        x(n+1-i) = xm + xl*z
        w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
        w(n+1-i) = w(i)
        end do
    end subroutine gaulegf
end module legendre



module observ
    use parameters_mod
    use legendre
    implicit none
    contains

    ! calculate number of electrons and holes from G< and G> 
    subroutine calc_charge(G_lesser,G_greater,nen,nsub,nk,E,NS,NB,nm_dev,nelec,pelec,midgap)
        complex(dp), intent(in) :: G_lesser(nm_dev,nm_dev,nen,nsub,nk)
        complex(dp), intent(in) :: G_greater(nm_dev,nm_dev,nen,nsub,nk)
        real(dp), intent(in)    :: E(nen),midgap(nm_dev)
        integer, intent(in)    :: NS,NB,nm_dev,nen,nsub,nk
        real(dp), intent(out)   :: nelec(nm_dev),pelec(nm_dev)
        real(dp)::dE, weights(nsub), xen(nsub)
        integer::ie,j,isub,ik
        dE=E(2)-E(1)
        call gaulegf(0.0d0, dble(dE), xen, weights, nsub) ! obtain the Legendre ordinates and weights    
        weights=weights/dble(nk)
        nelec=0.0d0
        pelec=0.0d0        
        do ik=1,nk
            !$omp parallel default(shared) private(j,isub,ie) 
            !$omp do 
            do j=1,nm_dev
                do isub=1,nsub
                    do ie=1,nen
                        if (E(ie)>midgap(j))then
                            nelec(j)=nelec(j)+aimag(G_lesser(j,j,ie,isub,ik))*weights(isub)
                        else
                            pelec(j)=pelec(j)-aimag(G_greater(j,j,ie,isub,ik))*weights(isub)
                        endif
                    enddo
                enddo
            enddo
            !$omp end do
            !$omp end parallel
        enddo
    end subroutine calc_charge

    ! calculate scattering collision integral from the self-energy
    ! I = sum_E Sig> G^< - Sig< G^>
    subroutine calc_collision(Sig_lesser,Sig_greater,G_lesser,G_greater,nen,en,nk,spindeg,nm_dev,I,Ispec)
        complex(8),intent(in),dimension(nm_dev,nm_dev,nen,nk)::G_greater,G_lesser,Sig_lesser,Sig_greater
        real(8),intent(in)::en(nen),spindeg
        integer,intent(in)::nen,nm_dev,nk
        complex(8),intent(out)::I(nm_dev,nm_dev) ! collision integral
        complex(8),intent(out),optional::Ispec(nm_dev,nm_dev,nen) ! collision integral spectrum
        !----
        complex(8),allocatable::B(:,:)
        real(dp)::dE
        integer::ie,ik        
        I=dcmplx(0.0d0,0.0d0)
        if (present(Ispec)) then 
            Ispec=dcmplx(0.0d0,0.0d0)
        endif
        dE=en(2)-en(1)
        do ik=1,nk
            !$omp parallel default(shared) private(B,ie) 
            allocate(B(nm_dev,nm_dev))
            !$omp do 
            do ie=1,nen
                call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,Sig_greater(:,:,ie,ik),nm_dev,G_lesser(:,:,ie,ik),nm_dev,czero,B,nm_dev)
                call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Sig_lesser(:,:,ie,ik),nm_dev,G_greater(:,:,ie,ik),nm_dev,cone,B,nm_dev) 
                I(:,:)=I(:,:)+B(:,:)
                if (present(Ispec)) then 
                    Ispec(:,:,ie)=B(:,:)*spindeg+Ispec(:,:,ie)
                endif
            enddo            
            !$omp end do
            deallocate(B)
            !$omp end parallel
        enddo
        I(:,:)=I(:,:)*dE/twopi*spindeg/dble(nk)        
    end subroutine calc_collision


    ! calculate bond current using I_ij = H_ij G<_ji - H_ji G^<_ij
    subroutine calc_bond_current(H,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
        complex(8),intent(in)::H(nm_dev,nm_dev),G_lesser(nm_dev,nm_dev,nen)
        real(8),intent(in)::en(nen),spindeg
        integer,intent(in)::nen,nm_dev ! number of E and device dimension
        real(8),intent(out)::tot_cur(nm_dev,nm_dev) ! total bond current density
        real(8),intent(out)::tot_ecur(nm_dev,nm_dev) ! total bond energy current density
        real(8),intent(out),optional::cur(nm_dev,nm_dev,nen) ! energy resolved bond current density
        !----
        complex(8),allocatable::B(:,:)
        integer::ie,io,jo        
        allocate(B(nm_dev,nm_dev))
        tot_cur=0.0d0  
        tot_ecur=0.0d0
        do ie=1,nen
            !$omp parallel default(shared) private(jo,io) 
            !$omp do 
            do io=1,nm_dev
                do jo=1,nm_dev
                B(io,jo)=H(io,jo)*G_lesser(jo,io,ie) - H(jo,io)*G_lesser(io,jo,ie)
                enddo
            enddo  
            !$omp end do
            !$omp end parallel  
            B=B*e_charge/twopi/hbar*e_charge*dble(spindeg)
            if (present(cur)) cur(:,:,ie) = dble(B)
            tot_ecur=tot_ecur+ en(ie)*dble(B)
            tot_cur=tot_cur+ dble(B)          
        enddo
        deallocate(B)
    end subroutine calc_bond_current

end module observ


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


module gf_dense 
    use parameters_mod
    use output
    use legendre
    use observ
    use open_boundary
    implicit none 

    contains

    ! 3D GW solver with two periodic directions (y,z)
    ! iterating G -> P -> W -> Sig 
    subroutine solve_gw_3D(niter,scba_tol,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
        alpha_mix,nen,nsub,En,nb,ns,nphiy,nphiz,Ham,H00lead,H10lead,T,V,&
        ndiag,num_lead,flatband,output_files,G_retarded,G_lesser,G_greater,W0_retarded,tr)
    !
        use fft_mod, only : conv1d => conv1d2, corr1d => corr1d2  
        use parameters_mod
        !  
        integer, intent(in) :: nen, nsub, nb, ns,niter,nm_dev,length, nphiz, nphiy, num_lead
        real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg,scba_tol
        complex(8),intent(in) :: Ham(nm_dev,nm_dev,nphiy*nphiz),H00lead(NB*NS,NB*NS,num_lead,nphiy*nphiz),H10lead(NB*NS,NB*NS,num_lead,nphiy*nphiz),T(NB*NS,nm_dev,num_lead,nphiy*nphiz)
        complex(8), intent(in):: V(nm_dev,nm_dev,nphiy*nphiz)
        integer,intent(in)::ndiag
        logical,intent(in)::flatband
        logical,intent(in) :: output_files
        complex(8),intent(out),dimension(nm_dev,nm_dev,nen,nsub,nphiy*nphiz) ::  G_retarded,G_lesser,G_greater
        complex(8),intent(out),dimension(nm_dev,nm_dev,nphiy*nphiz) ::  W0_retarded
        real(8),intent(out) ::Tr(nen,num_lead) ! current spectrum on leads    
        !------
        complex(8),dimension(:,:,:,:),allocatable ::  P_retarded,P_lesser,P_greater
        complex(8),dimension(:,:,:,:),allocatable ::  W_retarded,W_lesser,W_greater
        complex(8),dimension(:,:,:,:),allocatable ::  Sig_retarded,Sig_lesser,Sig_greater
        complex(8),dimension(:,:,:,:),allocatable ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
        complex(8),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
        complex(8),allocatable,dimension(:,:):: B ! tmp matrix
        real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:),sumcur(:,:,:),sumtot_cur(:,:),sumtot_ecur(:,:)
        complex(8),allocatable::Ispec(:,:,:),Itot(:,:)    
        real(8),allocatable::Te(:,:,:) ! transmission matrix spectrum
        real(8),allocatable::sumTr(:,:) ! current spectrum on leads summed over k
        real(8),allocatable::sumTe(:,:,:) ! transmission matrix spectrum summed over k
        integer :: iter,ie,nopmax
        integer :: i,j,nm,nop,l,h,iop,ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,isub        
        complex(8) :: dE
        real(8)::nelec(2),mu(2),pelec(2),temp(2)
        real(8)::weights(nsub),xen(nsub)
        real(8)::scba_error
        
        allocate(Sig_retarded(nm_dev,nm_dev,nen,nphiy*nphiz),Sig_lesser(nm_dev,nm_dev,nen,nphiy*nphiz),Sig_greater(nm_dev,nm_dev,nen,nphiy*nphiz))
        
        scba_error=1.0d0
        Sig_retarded = czero
        sig_lesser = czero
        sig_greater = czero
        !
        print *,'============ green_solve_gw_3D ============'
        dE = En(2) - En(1)
        call gaulegf(0.0d0, dble(dE), xen, weights, nsub) ! obtain the Legendre ordinates and weights    
        !
        allocate(siglead(NB*NS,NB*NS,nen,num_lead,nphiy*nphiz))
        ! get leads sigma
        do ikz=1, nphiy*nphiz
            siglead(:,:,:,1,ikz) = Sig_retarded(1:NB*NS,1:NB*NS,:,ikz)
            siglead(:,:,:,2,ikz) = Sig_retarded(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:,ikz)  
        enddo
        !
        allocate(tot_cur(nm_dev,nm_dev))
        allocate(tot_ecur(nm_dev,nm_dev))
        allocate(sumtot_cur(nm_dev,nm_dev))
        allocate(sumtot_ecur(nm_dev,nm_dev))
        allocate(cur(nm_dev,nm_dev,nen))
        allocate(sumcur(nm_dev,nm_dev,nen))
        allocate(Ispec(nm_dev,nm_dev,nen))
        allocate(Itot(nm_dev,nm_dev))    
        allocate(te(nen,num_lead,num_lead))
        allocate(sumtr(nen,num_lead))
        allocate(sumte(nen,num_lead,num_lead))
        if (flatband) then
            mu=(mus+mud)/2.0d0
            temp=(temps+tempd)/2.0d0
        else
            mu=(/ mus, mud /)
            temp=(/temps,tempd/)
        endif
        iter=0
        print '(a8,f15.4,a8,f15.4)', 'mus=',mu(1),'mud=',mu(2)
        do while ( (scba_error>=scba_tol).and.(iter<=niter)  )
            print *,'+ iter=',iter,'error=',scba_error
            print *, 'calc G'  
            sumtot_cur=0.0d0
            sumtot_ecur=0.0d0
            sumcur=0.0d0
            sumTr=0.0d0
            sumTe=0.0d0
            do ikz=1,nphiy*nphiz
                do isub=1,nsub
                    call calc_gf(nen,En+xen(isub),num_lead,nm_dev,(/nb*ns,nb*ns/),nb*ns,&
                        Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),&
                        T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),&
                        G_retarded(:,:,:,isub,ikz),G_lesser(:,:,:,isub,ikz),&
                        G_greater(:,:,:,isub,ikz),Tr,Te,mu,temp,flatband)
                    !call write_spectrum('ldos_kz'//string(ikz)//'_',iter,G_retarded(:,:,:,ikz),nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
                    call calc_bond_current(Ham(:,:,ikz),G_lesser(:,:,:,isub,ikz),nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)    
                    call write_current_spectrum('last_Jdens',ikz,cur,nen,en,length,NB,Lx)    
                    call write_spectrum('last_ldos',ikz,G_retarded(:,:,:,:,ikz),nen,nsub,En,xen,length,NB,Lx,(/1.0d0,-2.0d0/))
                    sumcur=sumcur + cur * weights(isub)
                    sumtot_cur=sumtot_cur + tot_cur * weights(isub)
                    sumtot_ecur=sumtot_ecur + tot_ecur * weights(isub)
                    sumTr=sumTr + Tr * weights(isub)
                    sumTe=sumTe + Te * weights(isub)
                enddo
            enddo
            sumcur=sumcur/dble(nphiy)/dble(nphiz)
            sumtot_cur=sumtot_cur/dble(nphiy)/dble(nphiz)
            sumtot_ecur=sumtot_ecur/dble(nphiy)/dble(nphiz)
            sumTr=sumTr/dble(nphiz)/dble(nphiy)
            sumTe=sumTe/dble(nphiz)/dble(nphiy)
            sumTr = sumTr *e_charge/twopi/hbar*e_charge*dble(spindeg)
            if (output_files) then
                if (flatband) then
                    print *,'flatband'
                    ! call write_spectrum_per_kz('gw_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
                    ! call write_spectrum_per_kz('gw_gamma-centered_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=0.0_dp,at_kz=0.0_dp)
                    call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
            
                else
                    call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
                    call write_spectrum_summed_over_kz('gw_ndos',iter,G_lesser(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
                    call write_spectrum_summed_over_kz('gw_pdos',iter,G_greater(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
                endif
                call write_current_spectrum('gw_Jdens',iter,sumcur,nen,en,length,NB,Lx)
                call write_current('gw_I',iter,sumtot_cur,length,NB,NS,Lx)
                call write_current('gw_EI',iter,sumtot_ecur,length,NB,NS,Lx)
                call write_transmission_spectrum('gw_trL',iter,sumTr(:,1),nen,En)
                call write_transmission_spectrum('gw_trR',iter,sumTr(:,2),nen,En)
                ! call write_transmission_spectrum('gw_TE_LR',iter,sumTe(:,1,2)*spindeg,nen,En)
                ! call write_transmission_spectrum('gw_TE_RL',iter,sumTe(:,2,1)*spindeg,nen,En)            
            endif            
            open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
            write(101,'(I4,2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
            close(101)
            write(*,'(I4,"  IDS=",2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
            !
            G_retarded=dcmplx(0.0d0*dble(G_retarded),aimag(G_retarded))
            G_lesser=dcmplx(0.0d0*dble(G_lesser),aimag(G_lesser))
            G_greater=dcmplx(0.0d0*dble(G_greater),aimag(G_greater))
            !        
            print *, 'calc P'
            allocate(P_retarded(nm_dev,nm_dev,nen,nphiy*nphiz),P_lesser(nm_dev,nm_dev,nen,nphiy*nphiz),P_greater(nm_dev,nm_dev,nen,nphiy*nphiz))
            !
            nopmax=nen/2-1           
            ! print *,'ndiag=',ndiag
            ! Pij^<>(hw,kz') = \int_dE Gij^<>(E,kz) * Gji^><(E-hw,kz-kz')
            ! Pij^r(hw,kz')  = \int_dE Gij^<(E,kz) * Gji^a(E-hw,kz-kz') + Gij^r(E,kz) * Gji^<(E-hw,kz-kz')        
            do iq=1,nphiy*nphiz    
                iqz = mod(iq-1,nphiz)+1
                iqy = (iq-1) / nphiz +1
                P_lesser(:,:,:,iq) = dcmplx(0.0d0,0.0d0)
                P_greater(:,:,:,iq) = dcmplx(0.0d0,0.0d0)    
                P_retarded(:,:,:,iq) = dcmplx(0.0d0,0.0d0)                
                do iky=1,nphiy
                    do ikz=1,nphiz              
                        ik=ikz + (iky-1)*nphiz
                        ikzd=ikz-iqz + nphiz/2
                        ikyd=iky-iqy + nphiy/2
                        if (ikzd<1) ikzd=ikzd+nphiz
                        if (ikzd>nphiz) ikzd=ikzd-nphiz
                        if (ikyd<1) ikyd=ikyd+nphiy
                        if (ikyd>nphiy) ikyd=ikyd-nphiy   
                        if (nphiy==1)   ikyd=1
                        if (nphiz==1)   ikzd=1             
                        ikd=ikzd + (ikyd-1)*nphiz
                        !$omp parallel default(shared) private(l,h,i,j,isub) 
                        !$omp do        
                        do i = 1, nm_dev        
                            l=max(i-ndiag,1)
                            h=min(nm_dev,i+ndiag)                        
                            do j = l,h
                                do isub = 1,nsub
                                    P_lesser(i,j,:,iq) = P_lesser(i,j,:,iq) + corr1d(nen,G_lesser(i,j,:,isub,ik),G_greater(j,i,:,isub,ikd),method='fft') * weights(isub)
                                    P_greater(i,j,:,iq) = P_greater(i,j,:,iq) + corr1d(nen,G_greater(i,j,:,isub,ik),G_lesser(j,i,:,isub,ikd),method='fft') * weights(isub)         
                                    P_retarded(i,j,:,iq) = P_retarded(i,j,:,iq) + corr1d(nen,G_lesser(i,j,:,isub,ik),conjg(G_retarded(i,j,:,isub,ikd)),method='fft') * weights(isub) &
                                                                                + corr1d(nen,G_retarded(i,j,:,isub,ik),G_lesser(j,i,:,isub,ikd),method='fft') * weights(isub) 
                                enddo
                            enddo
                        enddo
                        !$omp end do
                        !$omp end parallel
                    enddo
                enddo                    
            enddo                
            dE = dcmplx(0.0d0 , -1.0d0 / 2.0d0 / pi ) * spindeg /dble(nphiy)/dble(nphiz) 
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
            allocate(W_retarded(nm_dev,nm_dev,nen,nphiy*nphiz),W_lesser(nm_dev,nm_dev,nen,nphiy*nphiz),W_greater(nm_dev,nm_dev,nen,nphiy*nphiz))
            !
            do iq=1,nphiy*nphiz        
                !$omp parallel default(shared) private(nop)
                !$omp do
                do nop=-nopmax+nen/2,nopmax+nen/2       
                    if (flatband) then
                        call calc_w(0,NB,NS,nm_dev,P_retarded(:,:,nop,iq),P_lesser(:,:,nop,iq),P_greater(:,:,nop,iq),&
                                    V(:,:,iq),W_retarded(:,:,nop,iq),W_lesser(:,:,nop,iq),W_greater(:,:,nop,iq))
                    else      
                        call calc_w(1,NB,NS,nm_dev,P_retarded(:,:,nop,iq),P_lesser(:,:,nop,iq),P_greater(:,:,nop,iq),&
                                    V(:,:,iq),W_retarded(:,:,nop,iq),W_lesser(:,:,nop,iq),W_greater(:,:,nop,iq))
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
            deallocate(P_greater,P_retarded,P_lesser)
            print *, 'calc SigGW'
            allocate(Sig_retarded_new(nm_dev,nm_dev,nen,nphiy*nphiz),Sig_lesser_new(nm_dev,nm_dev,nen,nphiy*nphiz),Sig_greater_new(nm_dev,nm_dev,nen,nphiy*nphiz))
            !          
            ! print *,'ndiag=',ndiag
            nopmax=nen/2-1
            Sig_greater_new = dcmplx(0.0d0,0.0d0)
            Sig_lesser_new = dcmplx(0.0d0,0.0d0)
            Sig_retarded_new = dcmplx(0.0d0,0.0d0)      
            ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)        
            do ik=1,nphiy*nphiz    
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
                        !$omp parallel default(shared) private(l,h,i,j,isub)
                        !$omp do  
                        do i = 1,nm_dev   
                            l=max(i-ndiag,1)
                            h=min(nm_dev,i+ndiag)       
                            do j = l,h
                                do isub = 1,nsub
                                    Sig_lesser_new(i,j,:,ik)=Sig_lesser_new(i,j,:,ik) + conv1d(nen,G_lesser(i,j,:,isub,ikd),W_lesser(i,j,:,iq),method='fft') * weights(isub)
                                    Sig_greater_new(i,j,:,ik)=Sig_greater_new(i,j,:,ik) + conv1d(nen,G_greater(i,j,:,isub,ikd),W_greater(i,j,:,iq),method='fft') * weights(isub) 
                                    Sig_retarded_new(i,j,:,ik)=Sig_retarded_new(i,j,:,ik) &
                                                                + conv1d(nen,G_lesser(i,j,:,isub,ikd),W_retarded(i,j,:,iq),method='fft') * weights(isub) &
                                                                + conv1d(nen,G_retarded(i,j,:,isub,ikd),W_lesser(i,j,:,iq),method='fft') * weights(isub) &
                                                                + conv1d(nen,G_retarded(i,j,:,isub,ikd),W_retarded(i,j,:,iq),method='fft') * weights(isub)                                               
                                enddo 
                            enddo
                        enddo      
                        !$omp end do
                        !$omp end parallel
                    enddo
                enddo            
            enddo        
            dE = dcmplx(0.0d0, 1.0d0/twopi) /dble(nphiy)/dble(nphiz)
            Sig_lesser_new = Sig_lesser_new  * dE
            Sig_greater_new= Sig_greater_new * dE
            Sig_retarded_new=Sig_retarded_new* dE
            Sig_retarded_new = dcmplx( dble(Sig_retarded_new), aimag(Sig_greater_new-Sig_lesser_new)/2.0d0 )
            !!! Sig_lesser_new = dcmplx( 0.0d0*dble(Sig_lesser_new), aimag(Sig_lesser_new) )
            !!! Sig_greater_new = dcmplx( 0.0d0*dble(Sig_greater_new), aimag(Sig_greater_new) )
            W0_retarded = W_retarded(:,:,nen/2,:)
            deallocate(W_lesser,W_greater,W_retarded)
            !
            ! symmetrize the selfenergies
            !$omp parallel default(shared) private(ie,ikz,B)
            allocate(B(nm_dev,nm_dev))
            !$omp do
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
            !$omp end do
            deallocate(B)
            !$omp end parallel
            if (length>3) then
                ! make sure self-energy is continuous near leads (by copying edge block)
                !$omp parallel default(shared) private(ie,iqz)
                !$omp do
                do ie=1,nen
                    do iqz=1,nphiy*nphiz
                        call expand_size_bycopy(Sig_retarded_new(:,:,ie,iqz),nm_dev,NB,2)
                        call expand_size_bycopy(Sig_lesser_new(:,:,ie,iqz),nm_dev,NB,2)
                        call expand_size_bycopy(Sig_greater_new(:,:,ie,iqz),nm_dev,NB,2)
                    enddo
                enddo
                !$omp end do        
                !$omp end parallel
            endif
            scba_error = sum( abs(Sig_retarded_new - Sig_retarded)**2 ) / sum( abs(Sig_retarded_new)**2 )
            open(unit=101,file='gw_scba_error.dat',status='unknown',position='append')
            write(101,'(I4,E16.6)') iter, scba_error
            close(101)
            iter=iter+1
            ! mixing with previous ones
            Sig_retarded = Sig_retarded+ alpha_mix * (Sig_retarded_new -Sig_retarded)
            Sig_lesser = Sig_lesser+ alpha_mix * (Sig_lesser_new -Sig_lesser)
            Sig_greater = Sig_greater+ alpha_mix * (Sig_greater_new -Sig_greater)  
            !
            deallocate(Sig_retarded_new,Sig_lesser_new,Sig_greater_new)
            !
            if (.not. flatband) then
                ! get leads sigma
                do iqz=1,nphiy*nphiz
                    siglead(:,:,:,1,iqz) = Sig_retarded(2*NB*NS+1:3*NB*NS,2*NB*NS+1:3*NB*NS,:,iqz)
                    siglead(:,:,:,2,iqz) = Sig_retarded(nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,nm_dev-3*NB*NS+1:nm_dev-2*NB*NS,:,iqz)    
                enddo            
            endif
            if (flatband) then            
                ! call write_spectrum_per_kz('gw_SigR',iter,Sig_retarded,nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,1.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
            endif 
            ! call write_spectrum_summed_over_kz('gw_SigR',iter,Sig_retarded,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            !  call write_spectrum_summed_over_kz('SigL',iter,Sig_lesser,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
            !  call write_spectrum_summed_over_kz('SigG',iter,Sig_greater,nen,En,nphiy*nphiz,length,NB,Lx,(/1.0,1.0/))
        end do  
        if (iter == niter) then 
            print *, "warning: max number of iterations reached!"
        
        endif
        ! calculate GF for the last time    
        print *, 'calc G for the last time'  
        sumtot_cur=0.0d0
        sumtot_ecur=0.0d0
        sumcur=0.0d0
        sumTr=0.0
        sumTe=0.0
        do ikz=1,nphiy*nphiz
            do isub=1,nsub
                call calc_gf(nen,En+xen(isub),num_lead,nm_dev,(/nb*ns,nb*ns/),nb*ns,&
                    Ham(:,:,ikz),H00lead(:,:,:,ikz),H10lead(:,:,:,ikz),Siglead(:,:,:,:,ikz),&
                    T(:,:,:,ikz),Sig_retarded(:,:,:,ikz),Sig_lesser(:,:,:,ikz),Sig_greater(:,:,:,ikz),&
                    G_retarded(:,:,:,isub,ikz),G_lesser(:,:,:,isub,ikz),&
                    G_greater(:,:,:,isub,ikz),Tr,Te,mu,temp,flatband)
                !call write_spectrum('ldos_kz'//string(ikz)//'_',iter,G_retarded(:,:,:,ikz),nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
                call calc_bond_current(Ham(:,:,ikz),G_lesser(:,:,:,isub,ikz),nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)    
                !call write_current_spectrum('Jdens_kz'//string(ikz)//'_',iter,cur,nen,en,length,NB,Lx)    
                sumcur=sumcur + cur * weights(isub)
                sumtot_cur=sumtot_cur + tot_cur * weights(isub)
                sumtot_ecur=sumtot_ecur + tot_ecur * weights(isub)
                sumTr=sumTr + Tr * weights(isub)
                sumTe=sumTe + Te * weights(isub)
            enddo
        enddo
        sumcur=sumcur/dble(nphiy)/dble(nphiz)
        sumtot_cur=sumtot_cur/dble(nphiy)/dble(nphiz)
        sumtot_ecur=sumtot_ecur/dble(nphiy)/dble(nphiz)
        sumTr=sumTr/dble(nphiz)/dble(nphiy)
        sumTe=sumTe/dble(nphiz)/dble(nphiy)
        if (flatband) then
            print *,'flatband'
            ! call write_spectrum_per_kz('gw_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=kt_cbm/(twopi/Ly),at_kz=ktz_cbm/(twopi/Lz))
            ! call write_spectrum_per_kz('gw_gamma-centered_ldos',iter,G_retarded(:,:,:,1,:),nen,En,nphiy,nphiz,length,NB,Lx,(/1.0d0,-2.0d0/),at_ky=0.0_dp,at_kz=0.0_dp)
            call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))

        else
            call write_spectrum_summed_over_kz('gw_ldos',iter,G_retarded(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-2.0d0/))
            call write_spectrum_summed_over_kz('gw_ndos',iter,G_lesser(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,1.0d0/))
            call write_spectrum_summed_over_kz('gw_pdos',iter,G_greater(:,:,:,:,:),nen,nsub,En,xen,nphiy*nphiz,length,NB,Lx,(/1.0d0,-1.0d0/))
        endif
        call write_current_spectrum('gw_Jdens',iter,sumcur,nen,en,length,NB,Lx)
        call write_current('gw_I',iter,sumtot_cur,length,NB,NS,Lx)
        call write_current('gw_EI',iter,sumtot_ecur,length,NB,NS,Lx)
        call write_transmission_spectrum('gw_trL',iter,sumTr(:,1)*spindeg,nen,En)
        call write_transmission_spectrum('gw_trR',iter,sumTr(:,2)*spindeg,nen,En)
        ! call write_transmission_spectrum('gw_TE_LR',iter,sumTe(:,1,2)*spindeg,nen,En)
        ! call write_transmission_spectrum('gw_TE_RL',iter,sumTe(:,2,1)*spindeg,nen,En)
        sumTr = sumTr *e_charge/twopi/hbar*e_charge*dble(spindeg)
        open(unit=101,file='gw_Id_iteration.dat',status='unknown',position='append')
        write(101,'(I4,2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
        close(101)
        write(*,'(I4,"  IDS=",2E16.6)') iter, -sum(sumTr(:,1)), sum(sumTr(:,2))
        !
        tr=sumTr
        deallocate(siglead)
        deallocate(sumcur,cur,tot_cur,tot_ecur,sumtot_cur,sumtot_ecur)
        deallocate(Ispec,Itot)
        deallocate(Te,sumTr,sumTe)    
        deallocate(Sig_greater,Sig_lesser,Sig_retarded)
    end subroutine solve_gw_3D
  
    function map_kq(sgn,ik,iq,nk)
        integer,intent(in)::ik,iq,nk,sgn 
        integer::map_kq 
        integer::ikd
        ikd = ik + sgn * iq + nk/2            
        if (ikd<1) ikd=ikd+nk
        if (ikd>nk) ikd=ikd-nk
        if (nk==1) ikd=1
        map_kq=ikd 
    end function map_kq

    function map_kq_2d(sgn,ik,iq,nky,nkz)
        integer,intent(in)::ik,iq,nky,nkz,sgn 
        integer::map_kq_2d 
        integer::ikd,ikyd,ikzd,iqy,iqz,iky,ikz 
        !
        iqz = mod(iq-1,nkz)+1
        iqy = (iq-1) / nkz +1
        ikz = mod(ik-1,nkz)+1
        iky = (ik-1) / nkz +1
        ikzd= map_kq(sgn,ikz,iqz,nkz)
        ikyd= map_kq(sgn,iky,iqy,nky)        
        ikd = ikzd + (ikyd-1) * nkz
        map_kq_2d = ikd
    end function map_kq_2d


    subroutine solve_eph(niter,scba_tol,nm_dev,Lx,length,spindeg,temps,tempd,mus,mud,&
        alpha_mix,nen,nsub,En,nb,ns,nphiy,nphiz,Ham,H00lead,H10lead,T,&
        ndiag,num_lead,flatband,output_files,G_retarded,G_lesser,G_greater,tr)
    ! 
        integer, intent(in) :: nen, nsub, nb, ns,niter,nm_dev,length, nphiz, nphiy, num_lead
        real(8), intent(in) :: En(nen), temps,tempd, mus, mud, alpha_mix,Lx,spindeg,scba_tol
        complex(8),intent(in) :: Ham(nm_dev,nm_dev,nphiy*nphiz),H00lead(NB*NS,NB*NS,num_lead,nphiy*nphiz),H10lead(NB*NS,NB*NS,num_lead,nphiy*nphiz),T(NB*NS,nm_dev,num_lead,nphiy*nphiz)        
        integer,intent(in)::ndiag
        logical,intent(in)::flatband
        logical,intent(in) :: output_files
        complex(8),intent(out),dimension(nm_dev,nm_dev,nen,nsub,nphiy*nphiz) ::  G_retarded,G_lesser,G_greater        
        real(8),intent(out) ::Tr(nen,num_lead) ! current spectrum on leads    
        !------        
        complex(8),dimension(:,:,:,:),allocatable ::  Sig_retarded,Sig_lesser,Sig_greater
        complex(8),dimension(:,:,:,:),allocatable ::  Sig_retarded_new,Sig_lesser_new,Sig_greater_new
        complex(8),allocatable::siglead(:,:,:,:,:) ! lead scattering sigma_retarded
        complex(8),allocatable,dimension(:,:):: B ! tmp matrix
        real(8),allocatable::cur(:,:,:),tot_cur(:,:),tot_ecur(:,:),wen(:),sumcur(:,:,:),sumtot_cur(:,:),sumtot_ecur(:,:)
        complex(8),allocatable::Ispec(:,:,:),Itot(:,:)    
        real(8),allocatable::Te(:,:,:) ! transmission matrix spectrum
        real(8),allocatable::sumTr(:,:) ! current spectrum on leads summed over k
        real(8),allocatable::sumTe(:,:,:) ! transmission matrix spectrum summed over k
        integer :: iter,ie,nopmax
        integer :: i,j,nm,nop,l,h,iop,ikz,iqz,ikzd,iky,iqy,ikyd,ik,iq,ikd,isub        
        complex(8) :: dE
        real(8)::mu(2),temp(2)
        real(8)::weights(nsub),xen(nsub)
        real(8)::scba_error


    end subroutine solve_eph

    ! calculate e-photon/phonon self-energies for single mode in thermal equilibrium 
    subroutine selfenergy_eph_mono(nm,nen,En,nop,nphiy,nphiz,ik,iq,M,G_lesser,G_greater,&
        Sig_lesser,Sig_greater,n_bose,gamma_q)
    ! 
        integer,intent(in)::nm,nen,nop,nphiy,nphiz,iq,ik
        real(8),intent(in)::en(nen),n_bose
        logical,intent(in)::gamma_q
        complex(8),intent(in),dimension(nm,nm)::M ! interaction matrix at q
        complex(8),intent(in),dimension(nm,nm,nen,nphiy*nphiz)::G_lesser,G_greater
        complex(8),intent(out),dimension(nm,nm,nen,nphiy*nphiz)::Sig_lesser,Sig_greater
        !---------
        integer::ie,ikd 
        complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix       
        ! Sig^<>(E,k) = M_{-q} [ N G^<>(E -+ hw,k-+q) + (N+1) G^<>(E +- hw,k+-q)] M_q       
        !$omp parallel default(shared) private(ie,A,B) 
        allocate(B(nm,nm))
        allocate(A(nm,nm))                                
        !$omp do
        do ie=1,nen
            ! Sig^<(E,k)
            A = czero            
            if (gamma_q) then 
                ikd = ik
            else
                ikd = map_kq_2d(-1,ik,iq,nphiy,nphiz)
            endif
            if (ie-nop>=1) A =A+ G_lesser(:,:,ie-nop,ikd) * n_bose
            if (gamma_q) then 
                ikd = ik
            else
                ikd = map_kq_2d(+1,ik,iq,nphiy,nphiz)
            endif
            if (ie+nop<=nen) A =A+ G_lesser(:,:,ie+nop,ikd) * (n_bose+1.0_dp)
            call zgemm('n','n',nm,nm,nm,cone,M,nm,A,nm,czero,B,nm) 
            call zgemm('n','n',nm,nm,nm,cone,B,nm,M,nm,czero,A,nm)     
            Sig_lesser(:,:,ie,ik) = Sig_lesser(:,:,ie,ik) + A             
            !
            ! Sig^>(E,k)
            A = czero
            if (gamma_q) then 
                ikd = ik
            else
                ikd = map_kq_2d(-1,ik,iq,nphiy,nphiz)
            endif
            if (ie-nop>=1) A =A+ G_greater(:,:,ie-nop,ikd) * (n_bose+1.0_dp)
            if (gamma_q) then 
                ikd = ik 
            else                
                ikd = map_kq_2d(+1,ik,iq,nphiy,nphiz)
            endif
            if (ie+nop<=nen) A =A+ G_greater(:,:,ie+nop,ikd) * n_bose
            call zgemm('n','n',nm,nm,nm,cone,M,nm,A,nm,czero,B,nm) 
            call zgemm('n','n',nm,nm,nm,cone,B,nm,M,nm,czero,A,nm)     
            Sig_greater(:,:,ie,ik) = Sig_greater(:,:,ie,ik) + A                
        enddo  
        !$omp end do        
        deallocate(A,B)
        !$omp end parallel
    end subroutine selfenergy_eph_mono
    


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



    subroutine calc_w(NBC,NB,NS,nm_dev,PR,PL,PG,V,WR,WL,WG)
        integer,intent(in)::nm_dev,NB,NS,NBC
        complex(8),intent(in),dimension(nm_dev,nm_dev)::PR,PL,PG,V
        complex(8),intent(out),dimension(nm_dev,nm_dev)::WR,WL,WG
        ! --------- local
        complex(8),allocatable,dimension(:,:)::B,S,M,LL,LG,VV
        complex(8),dimension(:,:),allocatable::V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,&
            PL00,PL01,PL10,PG00,PG01,PG10,LL00,LL01,LL10,LG00,LG01,LG10
        complex(8),dimension(:,:),allocatable::VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,&
            MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N
        complex(8),dimension(:,:),allocatable::dM11,xR11,dLL11,dLG11,dV11
        complex(8),dimension(:,:),allocatable::dMnn,xRnn,dLLnn,dLGnn,dVnn
        integer::i,NL,NR,NT,LBsize,RBsize
        real(8)::condL,condR
        NL=NB*NS ! left contact block size
        NR=NB*NS ! right contact block size
        NT=nm_dev! total size
        LBsize=NL*NBC
        RBsize=NR*NBC
        if (NBC>0) then
          allocate(B(NT,NT))
          allocate(M(NT,NT))
          allocate(S(NT,NT))
          allocate(LL(NT,NT))
          allocate(LG(NT,NT))
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
          !
          call get_OBC_blocks_for_W(NL,V(1:NL,1:NL),V(1:NL,NL+1:2*NL),PR(1:NL,1:NL),PR(1:NL,NL+1:2*NL),&
              PL(1:NL,1:NL),PL(1:NL,NL+1:2*NL),PG(1:NL,1:NL),PG(1:NL,NL+1:2*NL),NBC,&
              V00,V01,V10,PR00,PR01,PR10,M00,M01,M10,PL00,PL01,PL10,PG00,PG01,PG10,&
              LL00,LL01,LL10,LG00,LG01,LG10)
          !    
          call get_OBC_blocks_for_W(NR,V(NT-NR+1:NT,NT-NR+1:NT),transpose(conjg(V(NT-NR+1:NT,NT-2*NR+1:NT-NR))),PR(NT-NR+1:NT,NT-NR+1:NT),&
              transpose(PR(NT-NR+1:NT,NT-2*NR+1:NT-NR)),PL(NT-NR+1:NT,NT-NR+1:NT),-transpose(conjg(PL(NT-NR+1:NT,NT-2*NR+1:NT-NR))),&
              PG(NT-NR+1:NT,NT-NR+1:NT),-transpose(conjg(PG(NT-NR+1:NT,NT-2*NR+1:NT-NR))),NBC,&
              VNN,VNN1,VN1N,PRNN,PRNN1,PRN1N,MNN,MNN1,MN1N,PLNN,PLNN1,PLN1N,PGNN,PGNN1,PGN1N,&
              LLNN,LLNN1,LLN1N,LGNN,LGNN1,LGN1N)
          !
          !! S = V P^r
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,V,nm_dev,PR,nm_dev,czero,S,nm_dev)
          ! Correct first and last block to account for elements in the contacts
          S(1:LBsize,1:LBsize)=S(1:LBsize,1:LBsize) + matmul(V10,PR01)
          S(NT-RBsize+1:NT,NT-RBsize+1:NT)=S(NT-RBsize+1:NT,NT-RBsize+1:NT) + matmul(VNN1,PRN1N)
          !
          M = -S
          do i=1,nm_dev
             M(i,i) = 1.0d0 + M(i,i)
          enddo
          deallocate(S)
          !! LL=V P^l V'
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,V,nm_dev,PL,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V,nm_dev,czero,LL,nm_dev) 
          !Correct first and last block to account for elements in the contacts
          LL(1:LBsize,1:LBsize)=LL(1:LBsize,1:LBsize) + matmul(matmul(V10,PL00),V01) + &
            matmul(matmul(V10,PL01),V00) + matmul(matmul(V00,PL10),V01)
          !  
          LL(NT-RBsize+1:NT,NT-RBsize+1:NT)=LL(NT-RBsize+1:NT,NT-RBsize+1:NT) + &
            matmul(matmul(VNN,PLNN1),VN1N) + matmul(matmul(VNN1,PLN1N),VNN) + &
            matmul(matmul(VNN1,PLNN),VN1N)
          !
          !! LG=V P^g V'    
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,V,nm_dev,PG,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,V,nm_dev,czero,LG,nm_dev) 
          !Correct first and last block to account for elements in the contacts
          LG(1:LBsize,1:LBsize)=LG(1:LBsize,1:LBsize) + matmul(matmul(V10,PG00),V01) + &
            matmul(matmul(V10,PG01),V00) + matmul(matmul(V00,PG10),V01)
          LG(NT-RBsize+1:NT,NT-RBsize+1:NT)=LG(NT-RBsize+1:NT,NT-RBsize+1:NT) + &
            matmul(matmul(VNN,PGNN1),VN1N) + matmul(matmul(VNN1,PGN1N),VNN) + &
            matmul(matmul(VNN1,PGNN),VN1N)
            
          ! WR/WL/WG OBC Left
          call open_boundary_conditions(NL,M00,M10,M01,V01,xR11,dM11,dV11,condL)
          ! WR/WL/WG OBC right
          call open_boundary_conditions(NR,MNN,MNN1,MN1N,VN1N,xRNN,dMNN,dVNN,condR)
          allocate(VV(nm_dev,nm_dev))
          VV = V
          if (condL<1.0d-6) then   
              !
            !   call get_dL_OBC_for_W(NL,xR11,LL00,LL01,LG00,LG01,M10,'L', dLL11,dLG11)
              !
              M(1:LBsize,1:LBsize)=M(1:LBsize,1:LBsize) - dM11
              VV(1:LBsize,1:LBsize)=V(1:LBsize,1:LBsize) - dV11    
            !   LL(1:LBsize,1:LBsize)=LL(1:LBsize,1:LBsize) + dLL11
            !   LG(1:LBsize,1:LBsize)=LG(1:LBsize,1:LBsize) + dLG11    
          endif
          if (condR<1.0d-6) then    
              !
            !   call get_dL_OBC_for_W(NR,xRNN,LLNN,LLN1N,LGNN,LGN1N,MNN1,'R', dLLNN,dLGNN)
              !
              M(NT-RBsize+1:NT,NT-RBsize+1:NT)=M(NT-RBsize+1:NT,NT-RBsize+1:NT) - dMNN
              VV(NT-RBsize+1:NT,NT-RBsize+1:NT)=V(NT-RBsize+1:NT,NT-RBsize+1:NT)- dVNN
            !   LL(NT-RBsize+1:NT,NT-RBsize+1:NT)=LL(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLLNN
            !   LG(NT-RBsize+1:NT,NT-RBsize+1:NT)=LG(NT-RBsize+1:NT,NT-RBsize+1:NT) + dLGNN    
          endif
          !!!! calculate W^r = (I - V P^r)^-1 V    
          call invert_inplace(M,nm_dev) ! M -> xR
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,VV,nm_dev,czero,WR,nm_dev)           
          ! calculate W^< and W^> = W^r P^<> W^r dagger
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,LL,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,WL,nm_dev) 
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,LG,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,M,nm_dev,czero,WG,nm_dev)  
          deallocate(M,LL,LG,B,VV)
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
        else ! no OBC correction
          allocate(B(NT,NT))
          allocate(M(NT,NT))    
          !! M = I - V P^r
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,PR,nm_dev,czero,M,nm_dev)
          !
          do i=1,nm_dev
             M(i,i) = 1.0d0 + M(i,i)
          enddo    
          !!!! calculate W^r = (I - V P^r)^-1 V    
          call invert_inplace(M,nm_dev) ! M -> xR
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,M,nm_dev,V,nm_dev,czero,WR,nm_dev)           
          ! calculate W^< and W^> = W^r P^<> W^r dagger
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PL,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WL,nm_dev) 
          call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,WR,nm_dev,PG,nm_dev,czero,B,nm_dev) 
          call zgemm('n','c',nm_dev,nm_dev,nm_dev,cone,B,nm_dev,WR,nm_dev,czero,WG,nm_dev)  
          deallocate(B,M)  
        endif 
    end subroutine calc_w



    
end module gf_dense


module bse_dense    
    use parameters_mod,only:dp,twopi,pi,e_charge,epsilon0,m0_charge,hbar,c1i,czero,cone
    use legendre
    implicit none
    contains

    subroutine four_polarization(alpha,nm_dev,nen,nsub,en,nop,nk,ndiag,G_lesser,G_greater,G_retarded,L0,i,j,k,l)
        integer,intent(in) :: nm_dev,nen,nsub,nop,ndiag ,nk, i,j,k,l
        real(dp),intent(in) :: en(nen), alpha 
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen,nsub,nk) :: G_lesser,G_greater,G_retarded
        complex(dp),intent(out) :: L0
        ! ---
        real(dp) :: dE, weights(nsub), xen(nsub)
        integer :: ie, isub, ik, ikd
        ! the P4 IPA tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
        dE = ( En(2) - En(1) )  
        call gaulegf(0.0d0, dble(dE), xen, weights, nsub) ! obtain the Legendre ordinates and weights    
        weights=weights/twopi
        !                 
        L0=czero                   
        ! calculate P4_IPA from GG
        if (nk == 1) then 
            ik=1
            ikd=1
        endif
        do isub=1,nsub
            do ie=nop+1,nen                            
                L0 = L0 + &
                        (1.0_dp - alpha) * ( G_lesser(j,l,ie,isub,ik) * conjg(G_retarded(i,k,ie-nop,isub,ikd)) + &
                                            G_retarded(j,l,ie,isub,ik) * G_lesser(k,i,ie-nop,isub,ikd) )  * weights(isub) + &
                            alpha * 0.5_dp * ( G_greater(j,l,ie,isub,ik) * G_lesser(k,i,ie-nop,isub,ikd) - & 
                                            G_lesser(j,l,ie,isub,ik)  * G_greater(k,i,ie-nop,isub,ikd) )  * weights(isub) 
            enddo 
        enddo                                                                            
    end subroutine four_polarization

    ! solve the full Bethe-Salpeter Equation
    subroutine bse_fullsolve(alpha,spindeg,nm_dev,ndiag,nen,nsub,En,nop,nk,G_lesser,G_greater,G_retarded,W,V,solve,P_retarded,system,epsilon_M,L0,M,nn)
        use gf_dense, only: invert_inplace
        integer,intent(in)::nm_dev,nen,nop,ndiag,nsub,nk
        real(dp),intent(in)::en(nen),spindeg,alpha
        logical,intent(in),optional::solve
        integer, intent(out)::nn ! size of the system
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen,nsub,nk):: G_lesser,G_greater,G_retarded ! electron GFs
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: P_retarded ! 2-point polarization function with interacting electron-hole at frequency [[nop]]        
        complex(dp),intent(out),dimension(nm_dev*nm_dev ,nm_dev*nm_dev):: system ! system matrix
        complex(dp),intent(out),dimension(nm_dev*nm_dev ,nm_dev*nm_dev),optional:: L0,M ! L0 matrix
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
        integer :: N,i,j,k,l,p,q,ie,row,col, it, ii,jj
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
            table(:,it) = (/i,i/)            
        enddo
        ! then put the others, but first within the ndiag
        do i=1,nm_dev
            do j=1,nm_dev               
                if (i/=j) then                     
                    if (abs(i-j)<=ndiag) then
                        it=it+1
                        table(:,it) = (/i,j/)                    
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
                        table(:,it) = (/i,j/)
                    endif
                endif                    
            enddo
        enddo
        if (it/=N) then 
            print *, 'ERROR!'
            call abort
        endif
        N = it
        ! start computation
        allocate(Mmat(N,N), source=czero)        
        allocate(Lmat(N,N), source=czero)     
        allocate(Amat(N,N), source=czero)                
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
                    call four_polarization(alpha,nm_dev,nen,nsub,en,nop,nk,ndiag,G_lesser,G_greater,G_retarded,L0ijkl,i,j,k,l)
                    Lmat(row,col) = L0ijkl
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
        !$omp parallel default(shared) private(i,j)
        !$omp do
        do i=1,N
            do j=1,N
                system(i,j) = Amat(N-i+1,N-j+1) ! flip the matrix 
                L0(i,j) = Lmat(N-i+1,N-j+1) ! flip the matrix 
                M(i,j) = Mmat(N-i+1,N-j+1) ! flip the matrix 
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        !        
        if (lsolve) then 
            N=nn
            print *,'  start invert (I - L0 K)'
            !
            call invert_inplace(Amat(1:N,1:N),N)
            !
            print *,'  start computation L = (I - L0 K) \ L0  '
            !
            call zgemm('n','n',N,N,N,cone,Amat(1:N,1:N),N,Lmat(1:N,1:N),N,czero,Mmat(1:N,1:N),N)                 
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
            ! ! calculate RPA epsilon and output to file
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,c1i,V,nm_dev,Lmat(1:nm_dev,1:nm_dev),nm_dev,czero,epsilon_M,nm_dev) 
            do j=1,nm_dev 
                epsilon_M(j,j) = epsilon_M(j,j) + cone
            enddo      
            call invert_inplace(epsilon_M,nm_dev)        
            open(unit=99,file='rpa_epsilonM.dat',status='unknown', position="append", action="write")    
            epsM = sum( epsilon_M(nm_dev/2,1:nm_dev) )
            write(99,*) dble(nop)*(En(2)-En(1)) , - aimag(epsM), dble(epsM) ! - Im \epsilon^{-1}
            close(99)
            !
            ! ! calculate BSE epsilon_M and output to file        
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,P_retarded,nm_dev,czero,epsilon_M,nm_dev) 
            do j=1,nm_dev 
                epsilon_M(j,j) = epsilon_M(j,j) + cone
            enddo      
            open(unit=99,file='bse_epsilonM.dat',status='unknown', position="append", action="write")    
            epsM = sum( epsilon_M(nm_dev/2,1:nm_dev) )
            write(99,*) dble(nop)*(En(2)-En(1)) , aimag(epsM), dble(epsM) ! Im \epsilon_M -> absorption
            close(99)
        endif
        !                
        deallocate(Mmat,Lmat,Amat)
    end subroutine bse_fullsolve
  


    ! solve the full Bethe-Salpeter Equation
    subroutine bse_fullsolve_orig(alpha,spindeg,nm_dev,ndiag,nen,nsub,En,nop,nk,G_lesser,G_greater,G_retarded,W,V,P_retarded,system,epsilon_M)
        use gf_dense, only: invert_inplace
        integer,intent(in)::nm_dev,nen,nop,ndiag,nsub,nk
        real(dp),intent(in)::en(nen),spindeg,alpha
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen,nsub,nk):: G_lesser,G_greater,G_retarded ! electron GFs
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: P_retarded ! 2-point polarization function with interacting electron-hole at frequency [[nop]]    
        complex(dp),intent(out),dimension(nm_dev*nm_dev ,nm_dev*nm_dev ):: system ! system matrix
        complex(dp),intent(out),dimension(nm_dev,nm_dev) :: epsilon_M ! macroscopic dielectric function
        !---------
        complex(dp),dimension(:,:),allocatable :: Lmat ! two-particle Green's function 
        complex(dp),dimension(:,:),allocatable :: Mmat ! 4-point Kernel
        complex(dp),dimension(:,:),allocatable :: Amat ! 
        complex(dp) :: epsM, L0ijkl
        real(dp) :: start, finish
        integer :: N,i,j,k,l,p,q,ie,row,col, ne_margin
        logical :: lexchange
        !                  
        N = nm_dev*nm_dev        
        !
        allocate(Lmat(N,N), source=czero)
        allocate(Mmat(N,N), source=czero)
        allocate(Amat(N,N), source=czero)
        print *,'  start computation L0_ijkl = G_jl G_ki ...'
        !$omp parallel default(shared) private(i,j,k,l,row,col,L0ijkl)
        !$omp do
        do i=1,nm_dev
            do j=max(1,i-ndiag),min(nm_dev,i+ndiag)
                do k=max(1,i-ndiag),min(nm_dev,i+ndiag)
                    do l=max(1,i-ndiag),min(nm_dev,i+ndiag)           
                        if ((abs(j-k)<=ndiag).and.(abs(k-l)<=ndiag).and.(abs(j-l)<=ndiag)) then 
                            row= (j-1)*nm_dev + i               
                            col= (l-1)*nm_dev + k
                            call four_polarization(alpha,nm_dev,nen,nsub,en,nop,nk,ndiag,G_lesser,G_greater,G_retarded,L0ijkl,i,j,k,l)
                            Lmat(row,col) = L0ijkl                
                        endif
                    enddo
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        !
        !        
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
                row= (j-1)*nm_dev + i
                col= row 
                Mmat(row,col) = Mmat(row,col) + c1i *  W(i,j) 
            enddo
        enddo    
        !$omp end do
        !$omp end parallel  
        !call save_matrix('bse_M.dat',N, Mmat)
        !call save_matrix('bse_L0.dat',N, Lmat)
        !  
        call zgemm('n','n',N,N,N,-cone,Lmat,N,Mmat,N,czero,Amat,N) 
        !
        ! (I - L0 K) -> A
        do i=1,N 
        Amat(i,i) = Amat(i,i) + dcmplx(1.0_dp, 0.0_dp)
        enddo  
        !
        !$omp parallel default(shared) private(i,j)
        !$omp do
        do i=1,N
            do j=1,N
                system(i,j) = Amat(N-i+1,N-j+1) ! flip the matrix 
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        !
        print *,'  start invert (I - L0 K)'
        !
        call invert_inplace(Amat,N)
        !
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
                Amat(i,j) = Lmat(row,col)   
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        !
        ! ! calculate RPA epsilon and output to file
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,c1i,V,nm_dev,Amat(1:nm_dev,1:nm_dev),nm_dev,czero,epsilon_M,nm_dev) 
        do j=1,nm_dev 
            epsilon_M(j,j) = epsilon_M(j,j) + cone
        enddo      
        call invert_inplace(epsilon_M,nm_dev)        
        open(unit=99,file='rpa_epsilonM.dat',status='unknown', position="append", action="write")    
        epsM = sum( epsilon_M(nm_dev/2,1:nm_dev) )
        write(99,*) dble(nop)*(En(2)-En(1)) , - aimag(epsM), dble(epsM) ! - Im \epsilon^{-1}
        close(99)
        !
        ! ! calculate BSE epsilon_M and output to file        
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,P_retarded,nm_dev,czero,epsilon_M,nm_dev) 
        do j=1,nm_dev 
            epsilon_M(j,j) = epsilon_M(j,j) + cone
        enddo      
        open(unit=99,file='bse_epsilonM.dat',status='unknown', position="append", action="write")    
        epsM = sum( epsilon_M(nm_dev/2,1:nm_dev) )
        write(99,*) dble(nop)*(En(2)-En(1)) , aimag(epsM), dble(epsM) ! Im \epsilon_M -> absorption
        close(99)
        !
        deallocate(Lmat,Mmat,Amat)
    end subroutine bse_fullsolve_orig
  

    ! solve the Bethe-Salpeter Equation under approximation
    subroutine bse_solve(alpha,spindeg,nm_dev,ndiag,nen,nsub,En,nop,nk,G_lesser,G_greater,G_retarded,W,V,P_retarded)
        use gf_dense, only: invert_inplace
        integer,intent(in)::nm_dev,nen,nop,nsub,nk,ndiag
        real(dp),intent(in)::en(nen),spindeg,alpha
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen,nsub,nk):: G_lesser,G_greater,G_retarded ! electron GF
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: P_retarded ! 2-point polarization function with interacting electron-hole
        !---------
        integer :: ie,i,j,k,l,m,n,p
        integer :: nn,nm,pp,pq
        complex(dp),allocatable :: L0xx(:,:),A(:,:),B(:,:),Kxx(:,:), Mxx(:,:), Sxx(:,:)
        complex(dp) :: Qijkl, L0Kdd, epsM
        real(dp) :: dE  
        real(dp) :: weights(nsub), xen(nsub)
        integer :: isub, ik, ikd
        ! the P4 IPA tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
        dE = ( En(2) - En(1) )  
        call gaulegf(0.0d0, dble(dE), xen, weights, nsub) ! obtain the Legendre ordinates and weights    
        weights=weights/twopi
        !
        ik=1
        ikd=1
        !
        allocate( L0xx(nm_dev,nm_dev) , source=czero)
        allocate( Mxx(nm_dev,nm_dev) , source=czero)
        allocate( Sxx(nm_dev,nm_dev) , source=czero)
        allocate( Kxx(nm_dev,nm_dev) , source=czero)
        allocate( A(nm_dev,nm_dev) , source=czero)        
        print *,'  start computation L0_xx'
        !$omp parallel default(shared) private(i,j,ie,isub)
        !$omp do
        do i=1,nm_dev
            do j=1,nm_dev
                if (abs(i-j)<=ndiag) then
                    do isub=1,nsub
                        do ie=nop+1,nen
                            L0xx(i,j) = L0xx(i,j) + &
                                        (1.0_dp - alpha) * ( G_lesser(j,i,ie,isub,ik) * conjg(G_retarded(i,j,ie-nop,isub,ikd)) + &
                                                            G_retarded(j,i,ie,isub,ik) * G_lesser(j,i,ie-nop,isub,ikd) )* weights(isub) + &
                                            alpha * 0.5_dp * ( G_greater(j,i,ie,isub,ik) * G_lesser(j,i,ie-nop,isub,ikd) - & 
                                                            G_lesser(j,i,ie,isub,ik)  * G_greater(j,i,ie-nop,isub,ikd) )* weights(isub) 
                        enddo
                    enddo
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel         
        L0xx = L0xx 
        Kxx(:,:) = - c1i*V(:,:)*spindeg
        do i=1,nm_dev
            Kxx(i,i) = Kxx(i,i) + c1i*W(i,i)
        enddo
        !
        print *,'  start computation L0_xx K_xx'
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,L0xx,nm_dev,Kxx,nm_dev,czero,A,nm_dev) 
        do i=1,nm_dev
            A(i,i) = A(i,i) + cone
        enddo
        !
        print *,'  start computation S and M'
        !$omp parallel default(shared) private(i,j,k,l,Qijkl,L0Kdd,p,ie,isub)
        !$omp do
        do i=1,nm_dev
            do j=1,nm_dev
                do k=1,nm_dev
                    do l=1,nm_dev
                        if (k/=l) then
                            if ((abs(i-k)<=ndiag).and.(abs(j-l)<=ndiag).and.(abs(j-k)<=ndiag).and.&
                                (abs(i-l)<=ndiag).and.(abs(i-j)<=ndiag).and.(abs(k-l)<=ndiag)) then 
                                Qijkl=czero
                                do isub=1,nsub
                                    do ie=nop+1,nen
                                        ! L0xd * L0dx
                                        Qijkl = Qijkl + (1.0_dp - alpha) * ( G_lesser(i,l,ie,isub,ik) * conjg(G_retarded(i,k,ie-nop,isub,ikd))*G_lesser(l,j,ie,isub,ik) * conjg(G_retarded(k,j,ie-nop,isub,ikd)) + &
                                                                            G_retarded(i,l,ie,isub,ik) * G_lesser(k,i,ie-nop,isub,ikd) * G_retarded(l,j,ie,isub,ik) * G_lesser(j,k,ie-nop,isub,ikd) )* weights(isub) + &
                                                alpha * 0.5d0*( G_greater(i,l,ie,isub,ik) * G_lesser(k,i,ie-nop,isub,ikd) * G_greater(l,j,ie,isub,ik) * G_lesser(j,k,ie-nop,isub,ikd) &
                                                                - G_lesser(i,l,ie,isub,ik) * G_greater(k,i,ie-nop,isub,ikd) * G_lesser(l,j,ie,isub,ik) * G_greater(j,k,ie-nop,isub,ikd) )* weights(isub)
                                                                                            
                                    enddo
                                enddo                                
                                !
                                Qijkl = - c1i * Qijkl * W(k,l) 
                                L0Kdd=czero
                                do isub=1,nsub
                                    do ie=nop+1,nen
                                        L0Kdd = L0Kdd + &
                                                (1.0_dp - alpha) * ( G_lesser(k,k,ie,isub,ik) * conjg(G_retarded(l,l,ie-nop,isub,ikd)) + &
                                                                    G_retarded(k,k,ie,isub,ik) * G_lesser(l,l,ie-nop,isub,ikd) )* weights(isub) + &
                                                    alpha *0.5d0*( G_greater(k,k,ie,isub,ik) * G_lesser(l,l,ie-nop,isub,ikd)   &
                                                                - G_lesser(k,k,ie,isub,ik) * G_greater(l,l,ie-nop,isub,ikd)   )* weights(isub) 
                                                                                                
                                    enddo
                                enddo                                
                                !
                                L0Kdd = cone - c1i * L0Kdd * W(k,l) 
                                Qijkl = Qijkl / L0Kdd
                                Sxx(i,j) = Sxx(i,j) + Qijkl
                                do p=1,nm_dev                                    
                                    Mxx(i,p) = Mxx(i,p) - c1i * Qijkl * ( W(j,j) - V(j,p)*spindeg )                                     
                                enddo
                            endif
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
        !
        print *,'  start computation Axx^{-1}'
        call invert_inplace(A,nm_dev)
        Sxx(:,:) = L0xx(:,:) - Sxx(:,:)
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,A,nm_dev,Sxx,nm_dev,czero,Mxx,nm_dev) 
        P_retarded = - c1i * Mxx
        !
        ! call save_matrix('pr.dat',nm_dev, P_retarded)  
        !
        ! ! calculate RPA epsilon and output to file
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,c1i,V,nm_dev,L0xx(1:nm_dev,1:nm_dev),nm_dev,czero,A,nm_dev) 
        do j=1,nm_dev 
            A(j,j) = A(j,j) + cone
        enddo      
        call invert_inplace(A,nm_dev)        
        open(unit=99,file='rpa_epsilonM.dat',status='unknown', position="append", action="write")    
        epsM = sum( A(nm_dev/2,1:nm_dev) )
        write(99,*) dble(nop)*(En(2)-En(1)) , - aimag(epsM), dble(epsM) ! - Im \epsilon^{-1}
        close(99)
        !
        ! ! calculate BSE epsilon_M and output to file        
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,V,nm_dev,P_retarded,nm_dev,czero,A,nm_dev) 
        do j=1,nm_dev 
            A(j,j) = A(j,j) + cone
        enddo      
        open(unit=99,file='bse_epsilonM.dat',status='unknown', position="append", action="write")    
        epsM = sum( A(nm_dev/2,1:nm_dev) )
        write(99,*) dble(nop)*(En(2)-En(1)) , aimag(epsM), dble(epsM) ! Im \epsilon_M -> absorption
        close(99)
        deallocate(L0xx,Mxx,Sxx,Kxx,A)
    end subroutine bse_solve

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
    
end module bse_dense



module matrix_c 
    use parameters_mod
    implicit none 
    contains 


    pure function trace(A) result(tr)
        implicit none
        complex(dp), intent(in) :: A(:, :)
        complex(dp):: tr
        integer :: ii
        tr = sum((/(A(ii, ii), ii=1, size(A, 1))/))
    end function trace


    subroutine triMUL_C(A, B, C, R, trA, trB, trC)
        complex(8), intent(in), dimension(:, :) :: A, B, C
        complex(8), intent(out), allocatable :: R(:, :)
        character, intent(in) :: trA, trB, trC
        complex(8), allocatable, dimension(:, :) :: tmp
        integer :: n, m, k, kb
        if ((trA .ne. 'n') .and. (trA .ne. 'N') .and. (trA .ne. 't') .and. (trA .ne. 'T') &
            .and. (trA .ne. 'c') .and. (trA .ne. 'C')) then
            write (*, *) "ERROR in triMUL_C! trA is wrong: ", trA
            call abort()
        end if
        if ((trB .ne. 'n') .and. (trB .ne. 'N') .and. (trB .ne. 't') .and. (trB .ne. 'T') &
            .and. (trB .ne. 'c') .and. (trB .ne. 'C')) then
            write (*, *) "ERROR in triMUL_C! trB is wrong: ", trB
            call abort()
        end if
        if ((trA .eq. 'n') .or. (trA .eq. 'N')) then
            k = size(A, 2)
            m = size(A, 1)
        else
            k = size(A, 1)
            m = size(A, 2)
        end if
        if ((trB .eq. 'n') .or. (trB .eq. 'N')) then
            kb = size(B, 1)
            n = size(B, 2)
        else
            kb = size(B, 2)
            n = size(B, 1)
        end if
        if (k .ne. kb) then
            write (*, *) "ERROR in triMUL_C! Matrix dimension is wrong", k, kb
            call abort()
        end if
        call MUL_C(A, B, trA, trB, tmp)
        call MUL_C(tmp, C, 'n', trC, R)
        deallocate (tmp)
    end subroutine triMUL_C

    subroutine MUL_C(A, B, trA, trB, R)
        complex(8), intent(in) :: A(:, :), B(:, :)
        complex(8), intent(out), allocatable :: R(:, :)
        CHARACTER, intent(in) :: trA, trB
        integer :: n, m, k, kb, lda, ldb
        if ((trA .ne. 'n') .and. (trA .ne. 'N') .and. (trA .ne. 't') .and. (trA .ne. 'T') &
            .and. (trA .ne. 'c') .and. (trA .ne. 'C')) then
            write (*, *) "ERROR in MUL_C! trA is wrong: ", trA
            call abort()
        end if
        if ((trB .ne. 'n') .and. (trB .ne. 'N') .and. (trB .ne. 't') .and. (trB .ne. 'T') &
            .and. (trB .ne. 'c') .and. (trB .ne. 'C')) then
            write (*, *) "ERROR in MUL_C! trB is wrong: ", trB
            call abort()
        end if
        lda = size(A, 1)
        ldb = size(B, 1)
        if ((trA .eq. 'n') .or. (trA .eq. 'N')) then
            k = size(A, 2)
            m = size(A, 1)
        else
            k = size(A, 1)
            m = size(A, 2)
        end if
        if ((trB .eq. 'n') .or. (trB .eq. 'N')) then
            kb = size(B, 1)
            n = size(B, 2)
        else
            kb = size(B, 2)
            n = size(B, 1)
        end if
        if (k .ne. kb) then
            write (*, *) "ERROR in MUL_C! Matrix dimension is wrong", k, kb
            call abort()
        end if
        if (allocated(R)) then
            if ((size(R, 1) .ne. m) .or. (size(R, 2) .ne. n)) then
                deallocate (R)
                Allocate (R(m, n))
            end if
        else
            Allocate (R(m, n))
        end if        
        call zgemm(trA, trB, m, n, k, dcmplx(1.0d0, 0.0d0), A, lda, B, ldb, dcmplx(0.0d0, 0.0d0), R, m)
    end subroutine MUL_C

end module matrix_c



module rgf
    use parameters_mod
    use open_boundary 
    use gf_dense,only: invert => invert_inplace
    use matrix_c, only: MUL_C, triMUL_C, trace
    use omp_lib
    implicit none 
    contains

    !!  Fermi distribution function
    elemental Function ferm(a)
        Real(dp), intent(in) ::  a
        real(dp) :: ferm
        ferm = 1.0d0/(1.0d0 + Exp(a))
    End Function ferm

    !! RGF for a batch of energies
    subroutine rgf_energies(nx,mm,nm, nen, energies, mul, mur, TEMPl, TEMPr, Hii, H1i, Sii, sigma_lesser_ph, &
        sigma_r_ph, G_r, G_lesser, G_greater, Jdens, tr, tre, verbose)
        !!  Recursive Green's solver, solves these two equations together and compute the current
        !!  $$[zI-H-\Sigma^r] G^r = I$$
        !!  $$G^{<>} = G^r \Sigma^{<>} (G^r)^\dagger$$
        !!  $$J = [H,G^<]$$         
        integer, intent(in) :: mm !! max size of blocks
        integer, intent(in) :: nx !! lenght of the device    
        integer, intent(in) :: nen !! number of energies  
        complex(dp), intent(in) :: Hii(mm,mm,nx), H1i(mm,mm,nx + 1), Sii(mm,mm,nx), sigma_lesser_ph(mm,mm,nx,nen), sigma_r_ph(mm,mm,nx,nen)
        real(dp), intent(in)       :: energies(nen), mul(:, :), mur(:, :), TEMPr(:, :), TEMPl(:, :)    
        integer, intent(in) :: nm(nx) !! size of each block
        logical, intent(in) :: verbose
        complex(dp), intent(out) :: G_greater(mm,mm,nx,nen), G_lesser(mm,mm,nx,nen), G_r(mm,mm,nx,nen), Jdens(mm,mm,nx,nen)            
        real(dp), intent(out)      :: tr(nen), tre(nen)    
        integer :: ie 
        print *, 'calc G'
        !$omp parallel default(shared) private(ie)
        !$omp do
        do ie = 1,nen 
            call rgf_std(nx,mm,nm, energies(ie), mul, mur, TEMPl, TEMPr, Hii, H1i, Sii, sigma_lesser_ph(:,:,:,ie), &
                sigma_r_ph(:,:,:,ie), G_r(:,:,:,ie), G_lesser(:,:,:,ie), G_greater(:,:,:,ie), Jdens(:,:,:,ie), tr(ie), tre(ie), verbose)
        enddo
        !$omp end do
        !$omp end parallel 
    end subroutine rgf_energies

    subroutine rgf_std(nx,mm,nm, En, mul, mur, TEMPl, TEMPr, Hii, H1i, Sii, sigma_lesser_ph, &
            sigma_r_ph, G_r, G_lesser, G_greater, Jdens, tr, tre, verbose)
        !!  Recursive Green's solver, solves these two equations together and compute the current
        !!  $$[zI-H-\Sigma^r] G^r = I$$
        !!  $$G^{<>} = G^r \Sigma^{<>} (G^r)^\dagger$$
        !!  $$J = [H,G^<]$$         
        integer, intent(in) :: mm !! max size of blocks
        integer, intent(in) :: nx !! lenght of the device    
        complex(dp), intent(in) :: Hii(mm,mm,nx), H1i(mm,mm,nx + 1), Sii(mm,mm,nx), sigma_lesser_ph(mm,mm,nx), sigma_r_ph(mm,mm,nx)
        real(dp), intent(in) :: En, mul(:, :), mur(:, :), TEMPr(:, :), TEMPl(:, :)    
        integer, intent(in) :: nm(nx) !! size of each block
        logical, intent(in) :: verbose
        complex(dp), intent(out) :: G_greater(mm,mm,nx), G_lesser(mm,mm,nx), G_r(mm,mm,nx), Jdens(mm,mm,nx)            
        real(dp), intent(out) :: tr, tre    
        !---- local variables
        complex(dp) :: Gl(mm,mm,nx), Gln(mm,mm,nx)    
        integer    :: M, M1, ii, jj
        complex(dp) :: z
        real(dp)    :: tim, start, finish, start_0
        complex(dp), allocatable :: sig(:, :), H00(:, :), H10(:, :)
        complex(dp), allocatable :: A(:, :), B(:, :), C(:, :), G00(:, :), GBB(:, :), sigmar(:, :), sigmal(:, :), GN0(:, :)
        !                
        z = dcmplx(En, 0.0d0)
        ! on the left contact
        ii = 1
        M = nm(ii)
        allocate (H00(M, M))
        allocate (H10(M, M))
        allocate (G00(M, M))
        allocate (GBB(M, M))
        allocate (sigmal(M, M))
        allocate (sig(M, M))
        Gl = czero
        Gln= czero
        !
        start = omp_get_wtime()
        start_0=start
        !
        !! $$H00 = H(i,i) + \Sigma_{ph}(i) * S(i,i)$$
        call MUL_c(sigma_r_ph(1:M,1:M,ii), Sii(1:M,1:M,ii), 'n', 'n', B)
        !
        H00 = Hii(1:M,1:M,ii) + B
        H10 = H1i(1:M,1:M,ii)
        call sancho(M, En, Sii(1:M,1:M,ii), H00, transpose(conjg(H10)), G00, GBB)
        !
        if (verbose) then 
        !$omp critical
            open (unit=10, file='sancho_g00.dat', position='append')
            write (10, *) En, 2, -aimag(trace(G00))
            close (10)
            open (unit=10, file='sancho_gbb.dat', position='append')
            write (10, *) En, 2, -aimag(trace(Gbb))
            close (10)
        !$omp end critical
        endif
        !
        !! $$\Sigma^R = H_{i,i+1} * G_{00} * H_{i+1,i}$$
        !! $$Gl(i) = [E*S_{i,i} - H00 - \Sigma_R]^{-1}$$
        call triMUL_c(H10, G00, H10, sigmal, 'n', 'n', 'c')
        B = z* Sii(1:M,1:M,ii) - H00 - sigmal
        call invert(B, M)
        Gl(1:M,1:M,ii) = B
        !
        !! $$Gln(i) = Gl(i) * [\Sigma_{ph}^<(i)*S(i,i) + (-(\Sigma^R - \Sigma_R^\dagger)*ferm(..))] * Gl(i)^\dagger$$
        call MUL_c(sigma_lesser_ph(1:M,1:M,ii), Sii(1:M,1:M,ii), 'n', 'n', B)
        sig = -(sigmal - transpose(conjg(sigmal)))*ferm((En - mur)/(BOLTZ*TEMPr))
        !
        sig = sig + B
        call triMUL_c(Gl(1:M,1:M,ii), sig, Gl(1:M,1:M,ii), B, 'n', 'n', 'c')
        Gln(1:M,1:M,ii) = B
        deallocate (G00, GBB, sig, H10)
        !
        finish = omp_get_wtime()
        if (verbose) print *, "--- left contact took seconds", finish - start
        start = finish
        !
        allocate (A(M, M))
        ! inside device l -> r
        do ii = 2, nx - 1
            M1= M
            M = nm(ii)
            if (size(H00, 1) .ne. M) then
                deallocate (H00, A)
                allocate (H00(M, M))
                allocate (A(M, M))
            end if
            call MUL_c(sigma_r_ph(1:M,1:M,ii), Sii(1:M,1:M,ii), 'n', 'n', B)
            H00 = Hii(1:M,1:M,ii) + B
            !
            !! $$H00 = H(i,i) + \Sigma_{ph}(i) * S(i,i)$$
            !! $$Gl(i) = [E*S(i,i) - H00 - H(i,i-1) * Gl(i-1) * H(i-1,i)]^{-1}$$
            call triMUL_c(H1i(1:M,1:M1,ii), Gl(1:M1,1:M1,ii - 1), H1i(1:M,1:M1,ii), B, 'n', 'n', 'c')
            A = z*Sii(1:M,1:M,ii) - H00 - B
            call invert(A, M)
            Gl(1:M,1:M,ii) = A
            !
            !! $$Gln(i) = Gl(i) * [\Sigma_{ph}^<(i)*S(i,i) + H(i,i+1)*Gln(i+1)*H(i+1,i)] * Gl(i)^\dagger$$
            call triMUL_c(H1i(1:M,1:M1,ii), Gln(1:M1,1:M1,ii - 1), H1i(1:M,1:M1,ii), B, 'n', 'n', 'c')
            call MUL_c(sigma_lesser_ph(1:M,1:M,ii), Sii(1:M,1:M,ii), 'n', 'n', A)
            B = B + A
            call triMUL_c(Gl(1:M,1:M,ii), B, Gl(1:M,1:M,ii), A, 'n', 'n', 'c')
            Gln(1:M,1:M,ii) = A
        end do
        !
        finish = omp_get_wtime()
        if (verbose) print *, "--- first pass took seconds", finish - start
        start = finish
        !
        ! on the right contact
        ii = nx
        M1= M
        M = nm(ii)
        allocate (H10(M, M))
        allocate (G00(M, M))
        allocate (GBB(M, M))
        allocate (sig(M, M))
        allocate (sigmar(M, M))
        if (size(H00, 1) .ne. M) then
            deallocate (H00)
            allocate (H00(M, M))
        end if
        !
        call MUL_c(sigma_r_ph(1:M,1:M,ii), Sii(1:M,1:M,ii), 'n', 'n', B)
        H00 = Hii(1:M,1:M,ii) + B
        H10 = H1i(1:M,1:M,nx + 1)
        !
        call sancho(M, En, Sii(1:M,1:M,ii), H00, H10, G00, GBB)
        !
        call triMUL_c(H10, G00, H10, sigmar, 'c', 'n', 'n')
        !
        if (verbose) then 
        !$omp critical
            open (unit=10, file='sancho_g00.dat', position='append')
            write (10, *) En, 1, -aimag(trace(G00))
            close (10)
            open (unit=10, file='sancho_gbb.dat', position='append')
            write (10, *) En, 1, -aimag(trace(Gbb))
            close (10)
        !$omp end critical
        endif
        !
        call triMUL_c(H1i(1:M1,1:M,nx), Gl(1:M1,1:M1,nx - 1), H1i(1:M1,1:M,nx), B, 'n', 'n', 'c')
        A = z*Sii(1:M,1:M,ii) - H00 - B - sigmar
        !
        call invert(A, M)
        G_r(1:M,1:M,ii) = A
        Gl(1:M,1:M,ii) = A
        !
        !! $$\Sigma^< = \Sigma_11^< + \Sigma_{ph}^< + \Sigma_s^<$$
        call triMUL_c(H1i(1:M1,1:M,nx), Gln(1:M1,1:M1,nx - 1), H1i(1:M1,1:M,nx), B, 'n', 'n', 'c')
        call MUL_c(sigma_lesser_ph(1:M,1:M,nx), Sii(1:M,1:M,nx), 'n', 'n', A)
        sig = -(sigmar - transpose(conjg(sigmar)))*ferm((En - mul)/(BOLTZ*TEMPl))
        sig = sig + A + B
        !
        !! $$G^< = G * \Sigma^< * G^\dagger$$
        call triMUL_c(G_r(1:M,1:M,ii), sig, G_r(1:M,1:M,ii), B, 'n', 'n', 'c')
        !
        G_lesser(1:M,1:M,ii) = B
        G_greater(1:M,1:M,ii) = G_lesser(1:M,1:M,ii) + (G_r(1:M,1:M,ii) - transpose(conjg(G_r(1:M,1:M,ii))))
        !
        A = -(sigmar - transpose(conjg(sigmar)))*ferm((En - mul)/(BOLTZ*TEMPl))
        call MUL_c(A, G_greater(1:M,1:M,ii), 'n', 'n', B)
        A = -(sigmar - transpose(conjg(sigmar)))*(ferm((En - mul)/(BOLTZ*TEMPl)) - 1.0d0)
        call MUL_c(A, G_lesser(1:M,1:M,ii), 'n', 'n', C)
        !
        Jdens(1:M,1:M,ii) = B - C
        !
        tim = 0.0d0
        do jj = 1, M
            tim = tim + dble(Jdens(jj,jj,ii))
        end do
        tr = tim ! transmission
        deallocate (sigmar, sig, G00, GBB, H10)
        allocate (GN0(M, M))
        !
        !
        finish = omp_get_wtime()
        if (verbose) print *, "--- right contact took seconds", finish - start
        start = finish
        !
        ! inside device r -> l
        do ii = nx - 1, 1, -1
            M1= M
            M = nm(ii)
            !! $$A = G^<(i+1) * H(i+1,i) * Gl(i)^\dagger + G(i+1) * H(i+1,i) * Gln(i)$$
            call triMUL_c(G_lesser(1:M1,1:M1,ii + 1), H1i(1:M1,1:M,ii), Gl(1:M,1:M,ii), A, 'n', 'n', 'c')
            call triMUL_c(G_r(1:M1,1:M1,ii + 1), H1i(1:M1,1:M,ii), Gln(1:M,1:M,ii), B, 'n', 'n', 'n')
            A = A + B
            !! $$B = H(i,i+1) * A$$
            !! $$Jdens(i) = -2 * B$$
            call MUL_c(H1i(1:M1,1:M,ii), A, 'c', 'n', B)
            Jdens(1:M,1:M,ii) = -2.0d0*B(:, :)
            !
            !! $$GN0 = Gl(i) * H(i,i+1) * G(i+1)$$
            !! $$G(i) = Gl(i) + GN0 * H(i+1,i) * Gl(i)$$
            call MUL_c(Gl(1:M,1:M,ii), H1i(1:M1,1:M,ii), 'n', 'c', B)
            call MUL_c(B, G_r(1:M1,1:M1,ii + 1), 'n', 'n', GN0)
            call MUL_c(GN0, H1i(1:M1,1:M,ii), 'n', 'n', C)
            call MUL_c(C, Gl(1:M,1:M,ii), 'n', 'n', A)
            G_r(1:M,1:M,ii) = Gl(1:M,1:M,ii) + A
            !
            !! $$G^<(i) = Gln(i) + Gl(i) * H(i,i+1) * G^<(i+1) * H(i+1,i) *Gl(i)^\dagger$$
            call MUL_c(Gl(1:M,1:M,ii), H1i(1:M1,1:M,ii), 'n', 'c', B)
            call MUL_c(B, G_lesser(1:M1,1:M1,ii + 1), 'n', 'n', C)
            call MUL_c(C, H1i(1:M1,1:M,ii), 'n', 'n', A)
            call MUL_c(A, Gl(1:M,1:M,ii), 'n', 'c', C)
            G_lesser(1:M,1:M,ii) = Gln(1:M,1:M,ii) + C
            !
            !! $$G^<(i) = G^<(i) + GN0 * H(i+1,i) * Gln(i)$$
            call MUL_c(GN0, H1i(1:M1,1:M,ii), 'n', 'n', B)
            call MUL_c(B, Gln(1:M,1:M,ii), 'n', 'n', C)
            G_lesser(1:M,1:M,ii) = G_lesser(1:M,1:M,ii) + C
            !
            !! $$G^<(i) = G^<(i) + Gln(i) * H(i,i+1) * GN0$$
            call MUL_c(Gln(1:M,1:M,ii), H1i(1:M1,1:M,ii), 'n', 'c', B)
            call MUL_c(B, GN0, 'n', 'c', C)
            G_lesser(1:M,1:M,ii) = G_lesser(1:M,1:M,ii) + C
            !
            !! $$G^>(i) = G^<(i) + [G(i) - G(i)^\dagger]$$
            G_greater(1:M,1:M,ii) = G_lesser(1:M,1:M,ii) + (G_r(1:M,1:M,ii) - transpose(conjg(G_r(1:M,1:M,ii))))
        end do
        !
        finish = omp_get_wtime()
        if (verbose) print *, "--- second pass took seconds", finish - start
        start = finish
        !
        ii = 1
        M = nm(ii)
        ! on the left contact
        A = -(sigmal - transpose(conjg(sigmal)))*ferm((En - mur)/(BOLTZ*TEMPr))
        call MUL_c(A, G_greater(1:M,1:M,ii), 'n', 'n', B)
        A = -(sigmal - transpose(conjg(sigmal)))*(ferm((En - mur)/(BOLTZ*TEMPr)) - 1.0d0)
        call MUL_c(A, G_lesser(1:M,1:M,ii), 'n', 'n', C)
        tim = 0.0d0
        do jj = 1, M
            tim = tim + dble(B(jj, jj) - C(jj, jj))
        end do
        tre = tim
        deallocate (B, A, C, GN0, sigmal)
        !           
    end subroutine rgf_std


    ! calculate e-photon/phonon self-energies in the monochromatic assumption
    subroutine selfenergy_eph_mono(nm,nx,nen,En,nop,Mii,M1i,Mi1,G_lesser,G_greater,&
        Sig_lesser,Sig_greater,n_bose)
        integer,intent(in)::nm,nx,nen,nop
        real(8),intent(in)::en(nen),n_bose
        complex(8),intent(in),dimension(nm,nm,nx)::Mii ! interaction matrix diag blocks
        complex(8),intent(in),dimension(nm,nm,nx+1)::M1i,Mi1 ! interaction matrix 1st offdiag blocks
        complex(8),intent(in),dimension(nm,nm,nx,nen)::G_lesser,G_greater
        complex(8),intent(out),dimension(nm,nm,nx,nen)::Sig_lesser,Sig_greater
        !---------
        integer::ie,ix
        complex(8),allocatable::B(:,:),A(:,:) ! tmp matrix        
        ! Sig^<>(E) = M [ N G^<>(E -+ hw) + (N+1) G^<>(E +- hw)] M        
        !$omp parallel default(shared) private(ie,A,B,ix) 
        allocate(B(nm,nm))
        allocate(A(nm,nm))  
        !$omp do
        do ie=1,nen
            do ix=1,nx
                ! Sig^<(E)
                ! i,i = i,i @ i,i @ i,i
                A = czero
                if (ie-nop>=1) A =A+ G_lesser(:,:,ix,ie-nop) * n_bose
                if (ie+nop<=nen) A =A+ G_lesser(:,:,ix,ie+nop) * (n_bose+1.0_dp)
                call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
                call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
                Sig_lesser(:,:,ix,ie) = Sig_lesser(:,:,ix,ie) + A 
                ! i,i = i,i-1 @ i-1,i-1 @ i-1,i
                A = czero
                if (ie-nop>=1) A =A+ G_lesser(:,:,max(1,ix-1),ie-nop) * n_bose
                if (ie+nop<=nen) A =A+ G_lesser(:,:,max(1,ix-1),ie+nop) * (n_bose+1.0_dp)
                call zgemm('n','n',nm,nm,nm,cone,M1i(:,:,ix),nm,A,nm,czero,B,nm) 
                call zgemm('n','n',nm,nm,nm,cone,B,nm,Mi1(:,:,ix),nm,czero,A,nm)
                Sig_lesser(:,:,ix,ie) = Sig_lesser(:,:,ix,ie) + A 
                ! i,i = i,i+1 @ i+1,i+1 @ i+1,i
                A = czero
                if (ie-nop>=1) A =A+ G_lesser(:,:,min(nx,ix+1),ie-nop) * n_bose
                if (ie+nop<=nen) A =A+ G_lesser(:,:,min(nx,ix+1),ie+nop) * (n_bose+1.0_dp)
                call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix+1),nm,A,nm,czero,B,nm) 
                call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix+1),nm,czero,A,nm)     
                Sig_lesser(:,:,ix,ie) = Sig_lesser(:,:,ix,ie) + A 
                !
                ! Sig^>(E)
                ! i,i = i,i @ i,i @ i,i
                A = czero
                if (ie-nop>=1) A =A+ G_greater(:,:,ix,ie-nop) * (n_bose+1.0_dp)
                if (ie+nop<=nen) A =A+ G_greater(:,:,ix,ie+nop) * n_bose
                call zgemm('n','n',nm,nm,nm,cone,Mii(:,:,ix),nm,A,nm,czero,B,nm) 
                call zgemm('n','n',nm,nm,nm,cone,B,nm,Mii(:,:,ix),nm,czero,A,nm)     
                Sig_greater(:,:,ix,ie) = Sig_greater(:,:,ix,ie) + A
                ! i,i = i,i-1 @ i-1,i-1 @ i-1,i
                A = czero
                if (ie-nop>=1) A =A+ G_greater(:,:,max(1,ix-1),ie-nop) * (n_bose+1.0_dp)
                if (ie+nop<=nen) A =A+ G_greater(:,:,max(1,ix-1),ie+nop) * n_bose
                call zgemm('n','n',nm,nm,nm,cone,M1i(:,:,ix),nm,A,nm,czero,B,nm) 
                call zgemm('n','n',nm,nm,nm,cone,B,nm,Mi1(:,:,ix),nm,czero,A,nm)     
                Sig_greater(:,:,ix,ie) = Sig_greater(:,:,ix,ie) + A
                ! i,i = i,i-1 @ i-1,i-1 @ i-1,i
                A = czero
                if (ie-nop>=1) A =A+ G_greater(:,:,min(nx,ix+1),ie-nop) * (n_bose+1.0_dp)
                if (ie+nop<=nen) A =A+ G_greater(:,:,min(nx,ix+1),ie+nop) * n_bose
                call zgemm('n','n',nm,nm,nm,cone,Mi1(:,:,ix+1),nm,A,nm,czero,B,nm) 
                call zgemm('n','n',nm,nm,nm,cone,B,nm,M1i(:,:,ix+1),nm,czero,A,nm)     
                Sig_greater(:,:,ix,ie) = Sig_greater(:,:,ix,ie) + A
            enddo
        enddo  
        !$omp end do
        deallocate(A,B)
        !$omp end parallel
    end subroutine selfenergy_eph_mono
    

     
    
end module rgf
