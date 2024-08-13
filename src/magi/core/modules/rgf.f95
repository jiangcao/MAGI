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
module rgf
    use parameters_mod
    use open_boundary 
    use gw_dense,only: invert => invert_inplace
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
        ! print *, 'calc G'
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
        sig = -(sigmal - transpose(conjg(sigmal)))*ferm((En - mul)/(BOLTZ*TEMPl))
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
        sig = -(sigmar - transpose(conjg(sigmar)))*ferm((En - mur)/(BOLTZ*TEMPr))
        sig = sig + A + B
        !
        !! $$G^< = G * \Sigma^< * G^\dagger$$
        call triMUL_c(G_r(1:M,1:M,ii), sig, G_r(1:M,1:M,ii), B, 'n', 'n', 'c')
        !
        G_lesser(1:M,1:M,ii) = B
        G_greater(1:M,1:M,ii) = G_lesser(1:M,1:M,ii) + (G_r(1:M,1:M,ii) - transpose(conjg(G_r(1:M,1:M,ii))))
        !
        A = -(sigmar - transpose(conjg(sigmar)))*ferm((En - mur)/(BOLTZ*TEMPr))
        call MUL_c(A, G_greater(1:M,1:M,ii), 'n', 'n', B)
        A = -(sigmar - transpose(conjg(sigmar)))*(ferm((En - mur)/(BOLTZ*TEMPr)) - 1.0d0)
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
        A = -(sigmal - transpose(conjg(sigmal)))*ferm((En - mul)/(BOLTZ*TEMPl))
        call MUL_c(A, G_greater(1:M,1:M,ii), 'n', 'n', B)
        A = -(sigmal - transpose(conjg(sigmal)))*(ferm((En - mul)/(BOLTZ*TEMPl)) - 1.0d0)
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
        Sig_lesser = czero
        Sig_greater= czero
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
