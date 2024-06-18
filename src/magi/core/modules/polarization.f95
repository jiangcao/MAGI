!===============================================================================
! Copyright (C) 2023 Jiang Cao
!
! This program is distributed under the terms of the GNU General Public License.
! See the file `LICENSE' in the root directory of this distribution, or obtain 
! a copy of the License at <https://www.gnu.org/licenses/gpl-3.0.txt>.
!
! Author: jiacao <jiacao@ethz.ch>
! Comment:
!  
! Maintenance:
!===============================================================================
module polarization 
    use parameters_mod,only:dp,twopi,pi,e_charge,epsilon0,m0_charge,hbar,c1i,czero,cone    
    use omp_lib
    use fft_mod, only : corr1d => corr1d2  
    implicit none
    contains
    !
    ! calculates independent-particle polarizability matrix at one frequency
    pure subroutine four_polarization(alpha,nm_dev,nen,en,nop,ndiag,&
       G_lesser,G_greater,G_retarded,i,j,k,l,L0)
       integer,intent(in) :: nm_dev,nen,nop,ndiag, i,j,k,l
       real(dp),intent(in) :: en(nen), alpha 
       complex(dp),intent(in),dimension(nm_dev,nm_dev,nen) :: G_lesser,G_greater,G_retarded
       complex(dp),intent(out) :: L0
       ! ---
       real(dp) :: dE, weights, xen
       integer :: ie, isub, ik, ikd
       ! the P4 tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
       dE = ( En(2) - En(1) )          
       weights = dE/twopi
       !                        
       L0 =    (1.0_dp - alpha) * ( sum( G_lesser(j,l,(nop+1):nen)   * conjg(G_retarded(i,k,1:(nen-nop))) ) &
                                 +  sum( G_retarded(j,l,(nop+1):nen) * G_lesser(k,i,1:(nen-nop)) ) )  &
               + alpha * 0.5_dp * ( sum( G_greater(j,l,(nop+1):nen) * G_lesser(k,i,1:(nen-nop)) )  & 
                                 -  sum( G_lesser(j,l,(nop+1):nen)  * G_greater(k,i,1:(nen-nop)) ) )  
       L0 = L0 * weights 
    end subroutine four_polarization
    !
    ! calculate independent-particle polarizability matrix at multiply frequencies
    subroutine four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
        G_lesser,G_greater,G_retarded,i,j,k,l,L0)        
        integer,intent(in) :: nm_dev,nen,nnop,nop(nnop),ndiag, i,j,k,l
        real(dp),intent(in) :: en(nen), alpha 
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen),target :: G_lesser,G_greater,G_retarded
        complex(dp),intent(out) :: L0(nnop)
        ! ---
        complex(dp),dimension(nen) :: Gl,Gg,Gr
        complex(dp),dimension(nen) :: Gl_down,Gg_down,Ga_down
        real(dp) :: dE, weights, xen, a1,a2
        integer :: ie, isub, ik, ikd
        complex(dp),dimension(nen) :: tmp
        ! the P4 tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
        dE = ( En(2) - En(1) )                               
        weights = dE/twopi
        a1=(1.0_dp - alpha)*weights
        a2=(alpha * 0.5_dp)*weights
        !                        
        Gl(1:nen) = G_lesser(j,l,1:nen)
        Gg(1:nen) = G_greater(j,l,1:nen)
        Gr(1:nen) = G_retarded(j,l,1:nen)
        !
        Gl_down(1:nen) = G_lesser(k,i,1:nen)                    
        Gg_down(1:nen) = G_greater(k,i,1:nen)
        Ga_down(1:nen) = conjg(G_retarded(i,k,1:nen))
        ! calculate P4_IPA from GG
        tmp = corr1d(nen,Gl,Ga_down,method='fft') * a1
        tmp = tmp  + corr1d(nen,Gr,Gl_down,method='fft') * a1
        tmp = tmp  + corr1d(nen,Gg,Gl_down,method='fft') * a2
        tmp = tmp  - corr1d(nen,Gl,Gg_down,method='fft') * a2
        L0(1:nnop) = tmp(nop(1:nnop)+nen/2)        
        !
    end subroutine four_polarization_fft
    !
    

end module polarization