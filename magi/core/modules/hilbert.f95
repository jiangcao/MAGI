
module hilbert
    use parameters_mod
    use fft_mod,only:do_mkl_dfti_fft

    implicit none 
    
    CONTAINS

    !!  Compute the hilbert transform with FFT
    subroutine hilbert_transform_fft(nse,nomega,nomegasf,sf_chi0,chi0)
      integer,intent(in) :: nomega,nomegasf,nse
      complex(dp),intent(inout) :: sf_chi0(nse,nomegasf)
      complex(dp),intent(inout) :: chi0(nse,nomega)
      ! using method proposed by Lucas [AIP Advances 2, 032144 (2012)]


    end subroutine hilbert_transform_fft


    !!  Compute the Hilbert Transform with matrix-matrix multiplication
    !!  the computational cost scales with O(nomegasf*nomega)
    !!  $\chi^0( r r' ; \omega_j) = \sum_i t_{ji} \chi^s( r r' ; \omega_i)
    !!  $\chi$ represents a real-space response function between r and r' at frequency $\omega$
    !!  $t_{ji}$ is obtained from the [[calc_kkweight]] function
    !!  `nse` is the number of specified elements (r r') in the matrix $\chi$ to compute HT
    subroutine hilbert_transform_mmm(nse,nomega,nomegasf,kkweight,sf_chi0,chi0)
      integer,intent(in) :: nomega,nomegasf,nse
      complex(dp) :: kkweight(nomegasf,nomega)
      complex(dp),intent(inout) :: sf_chi0(nse,nomegasf)
      complex(dp),intent(inout) :: chi0(nse,nomega)
      ! using method proposed by Miyake [PRB 61, 7172] 
      ! and Shishkin & Kresse [PRB 74, 035101]
      ! Compute chi0 = MATMUL(chi0_sf,kkweight)
      call ZGEMM('N','N',nse,nomega,nomegasf,cone,sf_chi0,nse,kkweight,nomegasf,czero,chi0,nse)
    end subroutine hilbert_transform_mmm

    !   This subroutine is copied and adapted from Abinit 'm_chi0tk.F90'
    !!  Calculate frequency dependent weights needed to perform the Hilbert transform
    !!  Subroutine needed to implement the calculation
    !!  of the polarizability using the spectral representation as proposed in:
    !!  PRB 74, 035101 (2006) [[cite:Shishkin2006]]
    !!  and PRB 61, 7172 (2000) [[cite:Miyake2000]]
    !!
    !! INPUTS
    !! nsp = number of frequencies where the imaginary part of the polarizability is evaluated
    !! ne = number of frequencies for the polarizability
    !! omegasp(nsp) = real frequencies for the imaginary part of the polarizability
    !! omegae(ne) =  imaginary frequencies for the polarizability
    !! delta = small imaginary part used to avoid poles, input variables
    !!
    !! OUTPUT
    !! kkweight(nsp,ne) = frequency dependent weights Eq A1 PRB 74, 035101 (2006) [[cite:Shishkin2006]]
    !!
    subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)
      integer,intent(in) :: ne,nsp
      real(dp),intent(in) :: delta,omegamax
      real(dp),intent(in) :: omegasp(nsp)
      complex(dp),intent(in) :: omegae(ne)
      complex(dp),intent(out) :: kkw(nsp,ne)
      ! --- local variables
      integer :: isp,je
      real(dp) :: eta,xx1,xx2,den1,den2
      complex(dp) :: c1,c2,wt
      !
      kkw(:,:)=czero
      !
      do je=1,ne
        eta=delta
        wt=omegae(je)
        ! Not include shift at omega==0, what about metallic systems?
        if (abs(real(omegae(je)))<tol6 .and. abs(aimag(wt))<tol6) eta=tol12
        !  Not include shift along the imaginary axis
        if (abs(aimag(wt))>tol6) eta=zero
        do isp=1,nsp
          if (isp==1) then
            ! Skip negative point, should check that this would not lead to spurious effects
            c1=czero
            den1=one
          else
            xx1=omegasp(isp-1)
            xx2=omegasp(isp)
            den1= xx2-xx1
            c1= -(wt-xx1+c1i*eta)*log( (wt-xx2+c1i*eta)/(wt-xx1+c1i*eta) )&
      &          +(wt+xx1-c1i*eta)*log( (wt+xx2-c1i*eta)/(wt+xx1-c1i*eta) )
            c1= c1/den1
          end if
          xx1=omegasp(isp)
          if (isp==nsp) then
            ! Skip last point should check that this would not lead to spurious effects
            xx2=omegamax
          else
            xx2=omegasp(isp+1)
          end if
          den2=xx2-xx1
          c2=  (wt-xx2+c1i*eta)*log( (wt-xx2+c1i*eta)/(wt-xx1+c1i*eta) )&
      &        -(wt+xx2-c1i*eta)*log( (wt+xx2-c1i*eta)/(wt+xx1-c1i*eta) )
          c2= c2/den2
          kkw(isp,je)=  c1/den1 + c2/den2
        end do
      end do
      end subroutine calc_kkweight

end module hilbert      