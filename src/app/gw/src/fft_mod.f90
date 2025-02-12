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

module fft_mod
    use parameters_mod
    implicit none 
    
    CONTAINS
    
    ! Z(nop) = sum_ie X(ie) Y(ie-nop)
    function corr1d(n,X,Y,method) result(Z)
      integer, intent(in)::n
      character(len=*),intent(in)::method
      complex(dp),intent(in)::X(n),Y(n)
      complex(dp)::Z(2*n-1)
      complex(dp),allocatable,dimension(:) :: X_, Y_
      integer::i,ie
      select case (trim(method))
        case default ! explicit index
          Z=czero
          do i=-n+1,n-1
            do ie = max(i+1,1),min(n,n+i)
              Z(i+n)=Z(i+n) + X(ie)*Y(ie-i)
            enddo
          enddo
        case('sum')
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
      complex(dp),intent(in)::X(n),Y(n)
      complex(dp)::Z(n)
      complex(dp),allocatable,dimension(:) :: X_, Y_, Z_
      integer::i,ie
      select case (trim(method))
        case default ! explicit index
          Z=czero
          do i=-n/2+1,n/2-1
            do ie = max(i+1,1),min(n,n+i)
              Z(i+n/2)=Z(i+n/2) + X(ie)*Y(ie-i)
            enddo
          enddo
        case('sum')
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
      complex(dp),intent(in)::X(n),Y(n*2-1)
      complex(dp)::Z(n)
      complex(dp),allocatable,dimension(:)::x_in, y_in, z_in
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
        case('sum')
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
    ! using symmetry of Y, Y(-w) = - conjg( Y(w) )
    function conv1d_fock(n,nop,X,Y1,Y2,method) result(Z)
        integer, intent(in)::n,nop
        character(len=*),intent(in)::method
        complex(dp),intent(in)::X(n),Y1(nop),Y2(nop)
        !
        complex(dp)::Z(n),Y(n)
        complex(dp),allocatable,dimension(:)::x_in, y_in, z_in
        integer::i,iop
        !        
        Y=czero
        Y(n/2:n/2+nop-1) = Y1(1:nop)
        do i=2,nop
            Y(n/2+1-i) = - conjg(Y2(i))    
        enddo                
        select case (trim(method))
          case default ! explicit index
            Z = czero  
            do i= -n/2+1,n/2-1  
              iop = i+n/2
              Z((max(i,1)+1):min(n,(n+i)))=Z(max(i,1)+1:min(n,(n+i))) + X((max(i,1)+1-i):(min(n,(n+i))-i))*Y(iop)
            enddo    
          case('sum')
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
            ! 
            call do_mkl_dfti_conv(n*2-1,Y_in,X_in,Z_in)
            !
            Z=Z_in(1:n)
            deallocate(X_in,Y_in,Z_in)
        end select
    end function conv1d_fock
    
    ! Z(ie) = sum_nop X(ie-nop) Y(nop)
    function conv1d2(n,X,Y,method) result(Z)
      integer, intent(in)::n
      character(len=*),intent(in)::method
      complex(dp),intent(in)::X(n),Y(n)
      complex(dp)::Z(n)
      complex(dp),allocatable,dimension(:)::x_in, y_in, z_in
      integer::i,iop
      select case (trim(method))
        case default ! explicit index
          Z = czero  
          do i= -n/2+1,n/2-1  
            iop = i+n/2
            Z((max(i,1)+1):min(n,(n+i)))=Z(max(i,1)+1:min(n,(n+i))) + X((max(i,1)+1-i):(min(n,(n+i))-i))*Y(iop)
          enddo    
        case('sum')
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
      Complex(dp) :: X_in(n),Y_in(n),Z_out(n)
      Complex(dp) :: X_out(n),Y_out(n),Z_in(n)
     type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
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
!
    subroutine do_mkl_dfti_fft(n,x_in,x_out,forward)
   ! 1D complex to complex
     Use MKL_DFTI
      integer :: n
      Complex(dp),intent(in) :: X_in(n)
      Complex(dp),intent(out) :: x_out(n)
      logical,intent(in) :: forward
     type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle
     Integer :: Status
     ! Perform a complex to complex transform
     Status = DftiCreateDescriptor( My_Desc1_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 1, n )
     Status = DftiSetValue( My_Desc1_Handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
     Status = DftiCommitDescriptor( My_Desc1_Handle )
     if ( forward ) then 
       Status = DftiComputeForward( My_Desc1_Handle, x_in, x_out )
     else 
       Status = DftiComputeBackward( My_Desc1_Handle, x_in, x_out )
     endif 
     Status = DftiFreeDescriptor(My_Desc1_Handle)
     x_out(:) = x_out(:) / dble(n)
    end subroutine do_mkl_dfti_fft
    
end module fft_mod
    


