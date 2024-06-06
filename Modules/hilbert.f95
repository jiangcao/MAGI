
module hilbert
    use parameters_mod
    use fft_mod
    implicit none 
    
    CONTAINS
    

    !!  Compute the hilbert transform.
    !!
    !! INPUTS
    !! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$
    !! nomega=number of frequencies in $\chi0(q,\omega)$.
    !! max_rest,min_res=max and min resonant transition energy (for this q-point)
    !! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
    subroutine hilbert_transformation(npwe,nomega,nomegasf,my_wl,my_wr,kkweight,sf_chi0,chi0,spmeth)
      integer,intent(in) :: spmeth,nomega,nomegasf,my_wl,my_wr,npwe
      complex(dp),intent(in) :: kkweight(nomegasf,nomega)
      complex(dp),intent(inout) :: sf_chi0(npwe,npwe,my_wl:my_wr)
      complex(dp),intent(inout) :: chi0(npwe,npwe,nomega)
      !Local variables-------------------------------
      !scalars
      integer :: ig2,my_nwp
      character(len=500) :: msg
      !arrays
      complex(dp),allocatable :: A_g1wp(:,:),H_int(:,:),my_kkweight(:,:)
      ! using method proposed by Miyake [PRB 61, 7172] 
      ! and Shishkin & Kresse [PRB 74, 035101]
      my_nwp = my_wr - my_wl +1
      !$omp parallel private(my_kkweight, A_g1wp, H_int, ig2)
      allocate(my_kkweight(my_wl:my_wr,nomega))
      my_kkweight = kkweight(my_wl:my_wr,:)
      allocate(A_g1wp(npwe, my_nwp))
      allocate(H_int(npwe, nomega))
      !$omp do
      do ig2=1,npwe
        A_g1wp = sf_chi0(:,ig2,:)
        ! Compute H_int = MATMUL(A_g1wp,my_kkweight)
        call XGEMM('N','N',npwe,nomega,my_nwp,cone,A_g1wp,npwe,my_kkweight,my_nwp,czero,H_int,npwe)
        chi0(:,ig2,:) = H_int
      end do
      deallocate(my_kkweight)
      deallocate(A_g1wp)
      deallocate(H_int)
      !$omp end parallel
    end subroutine hilbert_transformation

    !
    !!  Calculate frequency dependent weights needed to perform the Hilbert transform
    !!  Subroutine needed to implement the calculation
    !!  of the polarizability using the spectral representation as proposed in:
    !!  PRB 74, 035101 (2006) [[cite:Shishkin2006]]
    !!  and PRB 61, 7172 (2000) [[cite:Miyake2000]]
    !!
    !! INPUTS
    !! nsp=number of frequencies where the imaginary part of the polarizability is evaluated
    !! ne=number of frequencies for the polarizability (same as in epsilon^-1)
    !! omegasp(nsp)=real frequencies for the imaginary part of the polarizability
    !! omegae(ne)= imaginary frequencies for the polarizability
    !! delta=small imaginary part used to avoid poles, input variables
    !!
    !! OUTPUT
    !! kkweight(nsp,ne)=frequency dependent weights Eq A1 PRB 74, 035101 (2006) [[cite:Shishkin2006]]
    !!
    subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)
      integer,intent(in) :: ne,nsp
      real(dp),intent(in) :: delta,omegamax
      real(dp),intent(in) :: omegasp(nsp)
      complex(dp),intent(in) :: omegae(ne)
      complex(dp),intent(out) :: kkw(nsp,ne)
      !Local variables-------------------------------
      !scalars
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