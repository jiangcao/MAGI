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
    subroutine calc_collision(Sig_lesser,Sig_greater,G_lesser,G_greater,nen,en,nsub,nk,spindeg,nm_dev,I,Ispec)
        complex(8),intent(in),dimension(nm_dev,nm_dev,nen,nsub,nk)::G_greater,G_lesser
        complex(8),intent(in),dimension(nm_dev,nm_dev,nen,nk)::Sig_lesser,Sig_greater
        real(8),intent(in)::en(nen),spindeg
        integer,intent(in)::nen,nm_dev,nk,nsub
        complex(8),intent(out)::I(nm_dev,nm_dev) ! collision integral
        complex(8),intent(out),optional::Ispec(nm_dev,nm_dev,nen,nsub) ! collision integral spectrum
        !----
        complex(8),allocatable::B(:,:)
        real(dp)::dE
        integer::ie,ik,isub
        real(dp)::weights(nsub), xen(nsub)        
        dE=En(2)-En(1)
        call gaulegf(0.0d0, dble(dE), xen, weights, nsub) ! obtain the Legendre ordinates and weights    
        weights=weights/dble(nk)/twopi*spindeg
        !
        I=dcmplx(0.0d0,0.0d0)
        if (present(Ispec)) then 
            Ispec=dcmplx(0.0d0,0.0d0)
        endif        
        do ik=1,nk
            do isub=1,nsub
                !$omp parallel default(shared) private(B,ie) 
                allocate(B(nm_dev,nm_dev))
                !$omp do 
                do ie=1,nen
                    call zgemm('n','n',nm_dev,nm_dev,nm_dev,cone,Sig_greater(:,:,ie,ik),nm_dev,G_lesser(:,:,ie,isub,ik),nm_dev,czero,B,nm_dev)
                    call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Sig_lesser(:,:,ie,ik),nm_dev,G_greater(:,:,ie,isub,ik),nm_dev,cone,B,nm_dev) 
                    I(:,:) = I(:,:) + B(:,:)*weights(isub)
                    if (present(Ispec)) then 
                        Ispec(:,:,ie,isub) = B(:,:)*weights(isub) + Ispec(:,:,ie,isub)
                    endif
                enddo            
                !$omp end do
                deallocate(B)
                !$omp end parallel
            enddo
        enddo            
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


