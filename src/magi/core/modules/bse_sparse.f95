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
module bse_sparse    
    use parameters_mod,only:dp,twopi,pi,e_charge,epsilon0,m0_charge,hbar,c1i,czero,cone,light_speed    
    use omp_lib
    use polarization
    use sinv
    use observ
    use output
    use gw_dense, only : calc_w,solve_gw => solve_gw_1D_memsaving,calc_gf
    use eph_dense, only : selfenergy_eph_mono
    implicit none
    contains

    subroutine reshape_BTA_block2stack(num_blocks,blocksize,nm_dev,&
        Adiag, Aupper, Alower, Alowerarrow, Aupperarrow,&
        out_Adiag, out_Aupper, out_Alower, out_Alowerarrow, out_Aupperarrow)
        ! input
        integer,intent(in)::num_blocks,blocksize,nm_dev
        complex(dp),intent(in),dimension(:,:):: Adiag
        complex(dp),intent(in),dimension(:,:):: Aupper,Alower 
        complex(dp),intent(in),dimension(:,:):: Alowerarrow
        complex(dp),intent(in),dimension(:,:):: Aupperarrow 
        ! output
        complex(dp),intent(out),dimension(num_blocks,blocksize,blocksize):: out_Adiag
        complex(dp),intent(out),dimension(num_blocks-1,blocksize,blocksize):: out_Aupper,out_Alower 
        complex(dp),intent(out),dimension(num_blocks,nm_dev,blocksize):: out_Alowerarrow
        complex(dp),intent(out),dimension(num_blocks,blocksize,nm_dev):: out_Aupperarrow 
        !
        integer::ib,i,j 
        ! !$omp parallel default(shared) private(ib,i,j)
        ! !$omp do
        ! do ib=1,num_blocks            
        !     do i=1,blocksize
        !         do j=1,blocksize
        !             out_Adiag(ib,i,j) = Adiag(i,j+(ib-1)*blocksize)
        !             if (ib<num_blocks) then 
        !                 out_Aupper(ib,i,j) = Aupper(i,j+(ib-1)*blocksize)
        !                 out_Alower(ib,i,j) = Alower(i,j+(ib-1)*blocksize)
        !             endif                     
        !         enddo 
        !         do j=1,nm_dev
        !             out_Alowerarrow(ib,j,i) = Alowerarrow(j,i+(ib-1)*blocksize)
        !             out_Aupperarrow(ib,i,j) = Aupperarrow(i+(ib-1)*blocksize,j)
        !         enddo 
        !     enddo 
        ! enddo            
        ! !$omp end do
        ! !$omp end parallel  

        out_Adiag = reshape(Adiag, [num_blocks,blocksize,blocksize], order=[2,3,1])
        out_Aupper = reshape(Aupper, [num_blocks-1,blocksize,blocksize], order=[2,3,1])
        out_Alower = reshape(Alower, [num_blocks-1,blocksize,blocksize], order=[2,3,1])
        out_Aupperarrow = reshape(Aupperarrow, [num_blocks,blocksize,nm_dev], order=[2,1,3])
        out_Alowerarrow = reshape(Alowerarrow, [num_blocks,nm_dev,blocksize], order=[2,3,1])
    end subroutine reshape_BTA_block2stack


    ! preprocessing the sparsity pattern and decide the block_size and num_blocks in the BTA matrix
    subroutine bse_sparse_pre(nm_dev,ndiag,N,nnz,table,blocksize,num_blocks)
        ! input
        integer,intent(in)::nm_dev
        integer,intent(in)::ndiag 
        ! output
        integer,intent(out)::table(2,nm_dev*nm_dev) ! lookup table from 2-body to 1-body            
        integer,intent(out)::blocksize,num_blocks
        integer,intent(out)::N! size of the reduced 2-body system
        integer(8),intent(out)::nnz! nonzero       
        ! ---- local
        integer::i,j,k,l,it,bandwidth,col,row,NT              
        table = 0
        N = nm_dev**2 - (nm_dev-ndiag-1)*(nm_dev-ndiag) ! compressed system size ~ 2*nm_dev*ndiag-ndiag*ndiag                        
        ! construct a lookup table of reordered indices 
        ! tip， first put the i=j        
        do i=1,nm_dev            
            table(:,i) = [i,i]
        enddo
        ! then put the others, but within the ndiag
        it=nm_dev+1
        do i=1,nm_dev
            l = max(1,i-ndiag)
            k = min(nm_dev,i+ndiag)
            do  j= l , (i-1)
                table(:,it + j-l) = [i,j]                                
            enddo
            it=it + i-l
            do j= (i+1) , k                               
                table(:,it + j-i-1) = [i,j]                                    
            enddo
            it=it + k-i
        enddo
        if ((it-1)/=N) then 
            print *, '  ERROR!'
            call abort
        endif        
        print '("  nm_dev=                ", I10)', nm_dev
        print '("  ndiag =                ", I10)', ndiag
        print '("  resized system size=   ", I10)', N 
        ! determine coordinates of nnz
        nnz=0
        bandwidth=0
        do row = 1,N 
            do col = 1,N         
                i=table(1,row)
                j=table(2,row)
                k=table(1,col)
                l=table(2,col)    
                if ((abs(i-k)<=ndiag).and.(abs(j-l)<=ndiag).and.(abs(j-k)<=ndiag).and.&
                    (abs(i-l)<=ndiag).and.(abs(i-j)<=ndiag).and.(abs(k-l)<=ndiag)) then              
                    nnz=nnz+1 
                    if ((col>nm_dev).and.(row>nm_dev).and.(abs(col-row)>bandwidth)) then  
                        bandwidth = abs(col-row) 
                    endif
                endif
            enddo 
        enddo
        blocksize = bandwidth 
        num_blocks = ceiling( dble(N - nm_dev) / blocksize )  
        NT = blocksize * num_blocks         
        print '("  total arrow size=      ", I10)', NT
        print '("  arrow block size=      ", I10)', blocksize
        print '("  arrow number of blocks=", I10)', num_blocks
        print '("  nonzero elements=      ", F0.3, " Million")', dble(nnz)/1e6
        print '("  nonzero ratio =        ", F0.3 ," %")', dble(nnz)/(dble(NT+nm_dev)**2)*100
    end subroutine bse_sparse_pre

    ! solve the Bethe-Salpeter Equation with selected inversion 
    subroutine bse_sparse_solve(method,alpha,spindeg,nm_dev,ndiag,nen,nsub,En,nops,nnop,nk,&
        G_lesser,G_greater,G_retarded,W,V,solve_sigma,with_vertex,nb,ns,&
        P_retarded,P0_retarded,sig_retarded,sig_lesser,sig_greater)        
        ! in
        character(len=*),intent(in)::method
        integer,intent(in)::nm_dev,nen ! device dimension, number of energies        
        integer,intent(in)::nnop,nops(nnop) ! number of optical energies, optical energies in unit of energy interval
        integer,intent(in)::nsub,nk,ndiag ! number of legendre sub-energy nodes, number of k points, number of offdiagonals
        real(dp),intent(in)::en(nen),spindeg,alpha ! energy grid, spin degeneracy, empirical parameter
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen,nsub,nk):: G_lesser,G_greater,G_retarded ! electron Green Functions
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb operator
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb operator
        logical,intent(in),optional::solve_sigma , with_vertex
        integer,intent(in),optional::ns,nb ! NEEDED if solve_sigma, number of cells inside lead supercell, number of WF basis
        ! out 
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nnop):: P_retarded ! 2-point polarization function with interacting electron-hole
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nnop):: P0_retarded ! RPA polarization function 
        complex(dp),dimension(nm_dev,nm_dev,nen),intent(out) ::  Sig_retarded,Sig_lesser,Sig_greater
        !---- local
        complex(dp),allocatable,dimension(:,:):: Adiag,Aupper,Alower,Alowerarrow,Aupperarrow,Atip 
        complex(dp),allocatable,dimension(:,:):: Ktip
        complex(dp),allocatable,dimension(:,:,:):: vertex 
        complex(dp),allocatable,dimension(:)  :: Kdiag
        complex(dp),allocatable,dimension(:,:,:):: Ldiag,Lupper,Llower,Llowerarrow,Lupperarrow,Ltip 
        integer::N,blocksize,num_blocks,NT,local_nnop,iop,row,col,fliped_col,fliped_row,i,j,k,i1,i2,local_nops(1)
        integer(8):: nnz
        integer,allocatable,dimension(:,:)::table
        integer,allocatable,dimension(:,:)::ipiv_diagonal
        integer,allocatable,dimension(:)::ipiv_arrow_tip        
        real(dp) :: start, finish
        integer::l,h,ie,isub,ik,nop,ib,nepoch,ie1,ie2
        complex(dp),allocatable,dimension(:,:) ::  P_lesser,P_greater,W_retarded,W_lesser,W_greater,tmp 
        logical::lsolve_sigma, lwith_vertex
        complex::dE
        !
        if (present(solve_sigma)) then 
            lsolve_sigma = solve_sigma 
        else
            lsolve_sigma = .false. 
        endif
        if (present(with_vertex)) then 
            lwith_vertex = with_vertex 
        else
            lwith_vertex = .false. 
        endif
        print *,'====================== bse_sparse_solve ======================='                 
        ! pre-process the sparsity pattern of system
        print *, " pre-process ... "
        allocate(table(2,nm_dev*nm_dev))
        call bse_sparse_pre(nm_dev,ndiag,N,nnz,table,blocksize,num_blocks)
        ! prepare the memory
        NT = blocksize * num_blocks
        if (trim(method) == 'fft') then 
            local_nnop = nnop
        else 
            local_nnop = 1 ! compute one optical energy at a time
        endif 
        print *, " init memory ... "
        allocate(Ldiag(blocksize,NT,local_nnop), source=czero) 
        allocate(Lupper(blocksize,NT-blocksize,local_nnop), source=czero)
        allocate(Llower(blocksize,NT-blocksize,local_nnop), source=czero) 
        allocate(Llowerarrow(nm_dev,NT,local_nnop), source=czero) 
        allocate(Lupperarrow(NT,nm_dev,local_nnop), source=czero)              
        allocate(Ltip(nm_dev,nm_dev,local_nnop), source=czero)  
        allocate(Kdiag(NT), source=czero)  ! diagonal of Kernel
        allocate(Ktip(nm_dev,nm_dev), source=czero)  ! dense tip block of Kernel
        allocate(Adiag(blocksize,NT), source=czero) 
        allocate(Aupper(blocksize,NT-blocksize), source=czero)
        allocate(Alower(blocksize,NT-blocksize), source=czero) 
        allocate(Alowerarrow(nm_dev,NT), source=czero) 
        allocate(Aupperarrow(NT,nm_dev), source=czero)              
        allocate(Atip(nm_dev,nm_dev), source=czero)          
        allocate(ipiv_diagonal(blocksize,num_blocks), source=0)
        allocate(ipiv_arrow_tip(nm_dev), source=0)
        if (lsolve_sigma) then 
            sig_retarded = czero 
            sig_lesser = czero 
            sig_greater = czero            
            allocate(P_lesser(nm_dev,nm_dev), source=czero)
            allocate(P_greater(nm_dev,nm_dev), source=czero)
            allocate(W_lesser(nm_dev,nm_dev), source=czero)
            allocate(W_greater(nm_dev,nm_dev), source=czero)
            allocate(W_retarded(nm_dev,nm_dev), source=czero) 
            allocate(tmp(nm_dev,nm_dev))
            allocate(vertex(nm_dev,nm_dev,nm_dev), source=czero)
        endif 
        !
        if (local_nnop > 1) then  
            ! build BTA blocks of RPA polarization L0 and 2-body interaction kernal K                     
            call bse_sparse_build(method,alpha,spindeg,nm_dev,ndiag,nen,En,nops,local_nnop,blocksize,num_blocks,N,table,&
                                    G_lesser,G_greater,G_retarded,W,V,&
                                    Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag)       
            print *, " start selected inversion " 
            start = omp_get_wtime()                                     
            do iop = 1,local_nnop
                write(*, '(A)', advance="no") '.'
                ! build system matrix blocks (I - L0 @ K)
                ! print *, "  build system "
                call bse_sparse_build_system(blocksize,num_blocks,nm_dev,local_nnop,iop,&
                                Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip,&
                                Kdiag, Ktip,&
                                Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip)
                ! selected inversion of the system matrix
                ! call zbtatrf( blocksize, nm_dev, num_blocks, &
                !             Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                !             ipiv_diagonal,ipiv_arrow_tip)
                ! call zbtatri( blocksize, nm_dev, num_blocks, &
                !             Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                !             ipiv_diagonal,ipiv_arrow_tip)
                call zbtasinv(blocksize, nm_dev, num_blocks, &
                            Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip)
                ! compute P_retarded
                ! 
                Atip = matmul( Atip , Ltip(:,:,iop) )
                Atip = Atip + matmul( Alowerarrow, Lupperarrow(:,:,iop) ) 
                !$omp parallel default(shared) private(row,col,fliped_row,fliped_col,i,k)
                !$omp do
                do row=1,nm_dev
                    do col=1,nm_dev
                        fliped_row = nm_dev - col + 1
                        fliped_col = nm_dev - row + 1 
                        i=table(1,row)
                        k=table(1,col)
                        P_retarded(i,k,iop) =  - c1i * Atip(fliped_row,fliped_col)                
                        P0_retarded(i,k,iop) =  - c1i * Ltip(fliped_row,fliped_col,iop)    
                    enddo
                enddo                          
                !$omp end do
                !$omp end parallel  
            enddo
            finish = omp_get_wtime()
            print *
            print '("  computation time = ", F0.3 ," seconds.")', finish-start
        else 
            do iop = 1,nnop    
                print *," "   
                print '("  iop = ", I10, "/", I10 )', iop, nnop         
                ! build BTA blocks of RPA polarization L0 and 2-body interaction kernal K  
                local_nops = nops(iop)       
                call bse_sparse_build(method,alpha,spindeg,nm_dev,ndiag,nen,En,local_nops,1,blocksize,num_blocks,N,table,&
                                    G_lesser,G_greater,G_retarded,W,V,&
                                    Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag)                   
                ! build system matrix blocks (I - L0 @ K)
                start = omp_get_wtime()
                call bse_sparse_build_system(blocksize,num_blocks,nm_dev,1,1,&
                                Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip, &
                                Kdiag, Ktip,&
                                Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip)
                ! selected inversion of the system matrix
                ! call zbtatrf( blocksize, nm_dev, num_blocks, &
                !             Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                !             ipiv_diagonal,ipiv_arrow_tip)
                ! call zbtatri( blocksize, nm_dev, num_blocks, &
                !             Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                !             ipiv_diagonal,ipiv_arrow_tip)
                call zbtasinv(blocksize, nm_dev, num_blocks, &
                            Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip)
                finish = omp_get_wtime()
                print '("  selected inversion time = ", F0.3 ," seconds.")', finish-start            
                !                
                if (lsolve_sigma .and. lwith_vertex) then 
                    start = omp_get_wtime()                    
                    print '("  Start vertex computation ... ")'
                    ! compute vertex Gamma_ijk = L_ijkk , hence the upper-arrow block of L
                    Llowerarrow(:,:,1) = matmul( Atip , Llowerarrow(:,:,1) )
                    do ib = 1,num_blocks 
                        l = (ib-1)*blocksize + 1
                        h =  ib * blocksize 
                        if (ib > 1) then 
                            Llowerarrow(:,l:h,1) = Llowerarrow(:,l:h,1) + &
                                matmul( Alowerarrow(:, (l-blocksize):(h-blocksize)) , Lupper(:, (l-blocksize):(h-blocksize),1) )
                        end if 
                        Llowerarrow(:,l:h,1) = Llowerarrow(:,l:h,1) + &
                                matmul( Alowerarrow(:, l:h) , Ldiag(:, l:h,1) )
                        if (ib < num_blocks) then 
                            Llowerarrow(:,l:h,1) = Llowerarrow(:,l:h,1) + &
                                matmul( Alowerarrow(:, (l+blocksize):(h+blocksize)) , Llower(:, l:h,1) )
                        endif 
                        !
                    enddo 
                    !$omp parallel default(shared) private(row,col,fliped_row,fliped_col,i,j,k)
                    !$omp do
                    do row=1,NT
                        do col=1,nm_dev
                            fliped_row = nm_dev - col + 1
                            fliped_col = NT - row + 1 
                            if ( (row+nm_dev) <= N ) then
                                i=table(1,row+nm_dev)
                                j=table(2,row+nm_dev)
                                k=table(1,col)
                                if ((i>0).and.(j>0).and.(k>0)) then 
                                    vertex(i,j,k) =  Llowerarrow(fliped_row,fliped_col,1)                
                                endif
                            endif
                        enddo
                    enddo                          
                    !$omp end do
                    !$omp end parallel  
                    finish = omp_get_wtime()
                    print '("  vertex computation time = ", F0.3 ," seconds.")', finish-start 
                    start = finish
                endif 
                ! compute P_retarded P_ij = L_iijj , hence the arrowtip block of L                
                Atip = matmul( Atip , Ltip(:,:,1) )
                Atip = Atip + matmul( Alowerarrow, Lupperarrow(:,:,1) ) 
                !$omp parallel default(shared) private(row,col,fliped_row,fliped_col,i,k)
                !$omp do
                do row=1,nm_dev
                    do col=1,nm_dev
                        fliped_row = nm_dev - col + 1
                        fliped_col = nm_dev - row + 1 
                        i=table(1,row)
                        k=table(1,col)
                        P_retarded(i,k,iop) =  - c1i * Atip(fliped_row,fliped_col)    
                        P0_retarded(i,k,iop) =  - c1i * Ltip(fliped_row,fliped_col,1)    
                        if (lsolve_sigma .and. lwith_vertex) then  
                            vertex(i,i,k) = Atip(fliped_row,fliped_col)           
                        endif
                    enddo
                enddo                          
                !$omp end do
                !$omp end parallel  
                !
                if (lsolve_sigma) then 
                    !                
                    ! solve W
                    isub=1
                    ik=1
                    nop=nops(iop)
                    P_lesser = czero                     
                    ! P_retarded(:,:,iop) = dcmplx(0.0_dp, -abs(aimag(P_retarded(:,:,iop))))
                    P_greater = 2.0_dp * c1i * aimag(P_retarded(:,:,iop)) 
                    !
                    print '("  Start W computation ... ")'
                    start = finish
                    !
                    call calc_w(1,NB,NS,nm_dev,P_retarded(:,:,iop),P_lesser,P_greater,V,W_retarded,W_lesser,W_greater)
                    !
                    finish = omp_get_wtime()
                    print '("  W computation time = ", F0.3 ," seconds.")', finish-start 
                    !
                    start = finish
                    write(*, '(A)', advance="no") "  Start Sigma computation ... "                    
                    !
                    ! Accumulate the contribution from this optical frequency to
                    ! the self-energy
                    !
                    ! hw from -inf to +inf: Sig^<>_ij(E) = (i/2pi) \int_dhw G^<>_ij(E-hw) W^<>_ij(hw)  
                    if (lwith_vertex) then 
                        !
                        ! include vertex, $\Sigma = i G*W*\Gamma$
                        !
                        !$omp parallel default(shared) private(i1,i2,l,i,j,k,ie,ie1,ie2) 
                        !$omp do
                        do i=1,nm_dev                        
                            i1=max(i-ndiag,1)
                            i2=min(nm_dev,i+ndiag)   
                            do concurrent (j=i1:i2, l=i1:i2, k=i1:i2)    
                                if ((abs(l-j)<ndiag).and.(abs(l-k)<ndiag).and.(abs(j-k)<ndiag)) then   
                                    ie1 = max(nop,1) + 1
                                    ie2 = min(nen+nop,nen)
                                    Sig_lesser(i,j,ie1:ie2)=Sig_lesser(i,j,ie1:ie2) + G_lesser(i,l,(ie1-nop):(ie2-nop),isub,ik) * W_lesser(i,k) * vertex(l,j,k)                                
                                    Sig_greater(i,j,ie1:ie2)=Sig_greater(i,j,ie1:ie2) + G_greater(i,l,(ie1-nop):(ie2-nop),isub,ik) * W_greater(i,k) * vertex(l,j,k)   
                                    Sig_retarded(i,j,ie1:ie2)=Sig_retarded(i,j,ie1:ie2) + &
                                                            G_lesser(i,l,(ie1-nop):(ie2-nop),isub,ik) * W_retarded(i,k) * vertex(l,j,k) + &                                      
                                                            G_retarded(i,l,(ie1-nop):(ie2-nop),isub,ik) * W_lesser(i,k) * vertex(l,j,k) + &
                                                            G_retarded(i,l,(ie1-nop):(ie2-nop),isub,ik) * W_retarded(i,k) * vertex(l,j,k)                                                  
                                    !
                                    ie1 = max(-nop,1) + 1
                                    ie2 = min(nen-nop,nen)
                                    Sig_lesser(i,j,ie1:ie2)=Sig_lesser(i,j,ie1:ie2) + G_lesser(i,l,(ie1+nop):(ie2+nop),isub,ik) * W_greater(i,k) * vertex(l,j,k)   
                                    Sig_greater(i,j,ie1:ie2)=Sig_greater(i,j,ie1:ie2) + G_greater(i,l,(ie1+nop):(ie2+nop),isub,ik) * W_lesser(i,k) * vertex(l,j,k)   
                                    Sig_retarded(i,j,ie1:ie2)=Sig_retarded(i,j,ie1:ie2) - &
                                                            G_lesser(i,l,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_retarded(i,k)) * vertex(l,j,k) - &                                      
                                                            G_retarded(i,l,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_greater(i,k)) * vertex(l,j,k) - &
                                                            G_retarded(i,l,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_retarded(i,k)) * vertex(l,j,k)     
                                endif
                            enddo
                        enddo
                        !$omp end do
                        !$omp end parallel       
                    else
                        !
                        ! approx vertex as delta, $\Sigma = i G*W$
                        !
                        !$omp parallel default(shared) private(i1,i2,i,j,ie,ie1,ie2) 
                        !$omp do
                        do i=1,nm_dev                        
                            i1=max(i-ndiag,1)
                            i2=min(nm_dev,i+ndiag)   
                            do concurrent (j=i1:i2)    
                                ie1 = max(nop,0) + 1
                                ie2 = min(nen+nop,nen)
                                Sig_lesser(i,j,ie1:ie2)=Sig_lesser(i,j,ie1:ie2) + G_lesser(i,j,(ie1-nop):(ie2-nop),isub,ik) * W_lesser(i,j)                                
                                Sig_greater(i,j,ie1:ie2)=Sig_greater(i,j,ie1:ie2) + G_greater(i,j,(ie1-nop):(ie2-nop),isub,ik) * W_greater(i,j)   
                                Sig_retarded(i,j,ie1:ie2)=Sig_retarded(i,j,ie1:ie2) + &
                                                        G_lesser(i,j,(ie1-nop):(ie2-nop),isub,ik) * W_retarded(i,j) + &                                      
                                                        G_retarded(i,j,(ie1-nop):(ie2-nop),isub,ik) * W_lesser(i,j) + &
                                                        G_retarded(i,j,(ie1-nop):(ie2-nop),isub,ik) * W_retarded(i,j)                                                  
                                ! negative frequency part by symmetry
                                if (nop /= 0) then
                                    ie1 = max(-nop,0) + 1
                                    ie2 = min(nen-nop,nen)
                                    Sig_lesser(i,j,ie1:ie2)=Sig_lesser(i,j,ie1:ie2) - G_lesser(i,j,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_greater(i,j))   
                                    Sig_greater(i,j,ie1:ie2)=Sig_greater(i,j,ie1:ie2) - G_greater(i,j,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_lesser(i,j))   
                                    Sig_retarded(i,j,ie1:ie2)=Sig_retarded(i,j,ie1:ie2) - &
                                                            G_lesser(i,j,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_retarded(i,j)) - &                                      
                                                            G_retarded(i,j,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_greater(i,j)) - &
                                                            G_retarded(i,j,(ie1+nop):(ie2+nop),isub,ik) * conjg(W_retarded(i,j))     
                                endif
                            enddo
                        enddo
                        !$omp end do
                        !$omp end parallel 
                    endif 
                    finish = omp_get_wtime()
                    print *,""
                    print '("  Sigma computation time = ", F0.3 ," seconds.")', finish-start   
                endif
            enddo
            if (lsolve_sigma) then 
                dE = dcmplx( 0.0_dp, (En(2)-En(1))/twopi )               
                Sig_lesser  = Sig_lesser  * dE
                Sig_greater = Sig_greater * dE
                Sig_retarded=Sig_retarded * dE
                !
                Sig_retarded = dcmplx( dble(Sig_retarded), aimag(Sig_greater-Sig_lesser)/2.0_dp )
                !
                ! symmetrize the selfenergies
                do ie=1,nen
                    tmp(:,:)=transpose(Sig_retarded(:,:,ie))
                    Sig_retarded(:,:,ie) = (Sig_retarded(:,:,ie) + tmp(:,:))/2.0_dp  
                    tmp(:,:)=transpose(Sig_lesser(:,:,ie))
                    Sig_lesser(:,:,ie) = (Sig_lesser(:,:,ie) + tmp(:,:))/2.0_dp
                    tmp(:,:)=transpose(Sig_greater(:,:,ie))
                    Sig_greater(:,:,ie) = (Sig_greater(:,:,ie) + tmp(:,:))/2.0_dp
                enddo                
            endif
        endif          
    end subroutine bse_sparse_solve

    ! driver function for solving BSE with SCBA iteration
    subroutine bse_sparse_solve_scba(method,niter,nm_dev,Lx,length,spindeg,temp,mu,&
        alpha_mix,nen,En,nops,nnop,nb,ns,Ham,H00lead,H10lead,T,V,&
        ndiag,encut,egap,vertex,bse_sigma,flatband,output_files,inj_photon,nphot,m_phot,n_bose_phot,&
        G_retarded,G_lesser,G_greater,&
        current,transmission,P_retarded,P0_retarded) 
        ! in 
        character(len=*),intent(in)::method
        integer, intent(in) :: nen, nb, ns,niter,nm_dev,length
        integer, intent(in) :: ndiag
        integer,intent(in) :: nnop,nops(nnop) !! number of optical energies, optical energies in unit of energy interval        
        real(dp), intent(in) :: En(nen), temp(2), mu(2), alpha_mix, Lx, spindeg, egap
        complex(dp),intent(in) :: Ham(nm_dev,nm_dev),H00lead(NB*NS,NB*NS,2),H10lead(NB*NS,NB*NS,2),T(NB*NS,nm_dev,2)
        complex(dp), intent(in):: V(nm_dev,nm_dev)    
        real(dp),intent(in) :: encut(2) !! intraband and interband cutoff 
        logical, intent(in) :: vertex, bse_sigma, flatband, output_files
        !
        logical, intent(in) :: inj_photon !! photon injection
        complex(dp), intent(in),optional :: m_phot(nm_dev,nm_dev,1,1) !! e-photon interaction H of size (N,N,nk,nq)
        integer,intent(in),optional :: nphot !! photon energy in unit of energy interval 
        real(dp), intent(in),optional :: n_bose_phot !! Bose number of photon mode
        ! out 
        real(dp),intent(out)::current(nen,2) !! current spectrum on leads
        real(dp),intent(out)::transmission(nen,2) !! transmission matrix
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nen) ::  G_retarded,G_lesser,G_greater
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nnop) ::  P_retarded  ! bse 
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nnop) ::  P0_retarded ! rpa 
        ! --- local 
        real(dp),allocatable :: cur(:,:,:),tot_cur(:,:),tot_ecur(:,:)
        complex(dp),allocatable,dimension(:,:) :: W0_retarded,W0_lesser,W0_greater
        complex(dp),allocatable,dimension(:,:,:,:)::siglead
        complex(dp),allocatable,dimension(:,:,:)::Sig_r,Sig_l,Sig_g
        complex(dp),allocatable,dimension(:,:,:)::Sig_r_bse,Sig_l_bse,Sig_g_bse
        real(dp),dimension(nen,2,2)::te
        integer::iter
        real(dp)::E_phot
        real(dp), parameter :: pre_fact=((hbar/m0_kg)**2)*e_charge/(2.0d0*epsilon0*light_speed**3)
        !
        print *, " init memory ... "
        allocate(W0_greater(nm_dev,nm_dev))
        allocate(W0_lesser(nm_dev,nm_dev))
        allocate(W0_retarded(nm_dev,nm_dev))
        allocate(Sig_r(nm_dev,nm_dev,nen))
        allocate(Sig_l(nm_dev,nm_dev,nen))
        allocate(Sig_g(nm_dev,nm_dev,nen))
        allocate(Sig_r_bse(nm_dev,nm_dev,nen))
        allocate(Sig_l_bse(nm_dev,nm_dev,nen))
        allocate(Sig_g_bse(nm_dev,nm_dev,nen))
        allocate(siglead(NB*NS,NB*NS,nen,2))
        allocate(tot_cur(nm_dev,nm_dev))
        allocate(tot_ecur(nm_dev,nm_dev))
        allocate(cur(nm_dev,nm_dev,nen))
        ! solve G0W0 
        call solve_gw(&
                niter=0,nm_dev=nm_dev,lx=Lx,length=length,spindeg=spindeg,&
                temp=temp,mu=mu,alpha_mix=alpha_mix,&
                nen=nen,en=en,nb=nb,ns=ns,&
                ham=ham,h00lead=h00lead,h10lead=h10lead,t=t,v=v,&
                ndiag=ndiag,encut=encut,egap=egap,flatband=.False.,vertex=.False.,bse=.False.,output_files=output_files,&
                G_retarded=G_retarded,G_lesser=G_lesser,G_greater=G_greater,&
                Sig_retarded_new=Sig_r,sig_lesser_new=Sig_l,sig_greater_new=Sig_g,&
                current=current,transmission=transmission,&
                W0_retarded=W0_retarded,W0_lesser=W0_lesser,W0_greater=W0_greater )
        ! solve BSE
        call bse_sparse_solve(&
                method=method,alpha=0.99_dp,spindeg=spindeg,&
                nm_dev=nm_dev,ndiag=ndiag,nen=nen,nsub=1,&
                en=en,nops=nops,nnop=nnop,nk=1,&
                g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,&
                w=W0_retarded,v=v,solve_sigma=bse_sigma,with_vertex=vertex,nb=nb,ns=ns,&
                P_retarded=P_retarded,P0_retarded=P0_retarded,sig_retarded=Sig_r_bse,sig_lesser=sig_l_bse,sig_greater=sig_g_bse)
        if (bse_sigma) then   
            ! update self-energy 
            Sig_l = aimag( Sig_l_bse ) * c1i
            Sig_g = aimag( Sig_g_bse ) * c1i
            Sig_r = dble(Sig_r) + aimag(Sig_g - Sig_l) / 2.0_dp * c1i
            ! get leads sigma        
            siglead(:,:,:,1) = Sig_r(1:NB*NS,1:NB*NS,:)
            siglead(:,:,:,2) = Sig_r(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)    
            !
            print *, 'calc G ...'  
            !
            call calc_gf(nen,En,2,nm_dev,[nb*ns,nb*ns],nb*ns,&
                        Ham(:,:),H00lead(:,:,:),H10lead(:,:,:),Siglead(:,:,:,:),&
                        T(:,:,:),Sig_r(:,:,:),Sig_l(:,:,:),Sig_g(:,:,:),&
                        G_retarded(:,:,:),G_lesser(:,:,:),&
                        G_greater(:,:,:),current,Te,mu,temp,flatband)   
        endif 
        ! 
        if (inj_photon) then 
            !
            print *, 'calc Photon self-energy ...'      
            sig_l_bse = czero
            sig_g_bse = czero            
            E_phot = nphot * (En(2) - En(1))
            !
            print *,'n_bose_phot=',n_bose_phot
            print *,'E_phot (eV)=',E_phot
            print *,'max |M|^2 * n_bose_phot * pre-factor = ', maxval( abs( M_phot ) )**2 * n_bose_phot * pre_fact / E_phot 
            !
            call selfenergy_eph_mono(nm_dev,nen,En,nphot,1,1,1,1,1,1,1,M_phot,G_lesser,G_greater,&
                        Sig_l_bse,Sig_g_bse,n_bose=n_bose_phot,gamma_q=.True.)  
            !
            sig_l_bse = sig_l_bse * pre_fact / E_phot
            sig_g_bse = sig_g_bse * pre_fact / E_phot
            Sig_r_bse = aimag(Sig_g_bse - Sig_l_bse) / 2.0_dp * c1i
            ! accumulate 
            Sig_r = Sig_r + aimag( Sig_r_bse ) * c1i
            Sig_l = Sig_l + aimag( Sig_l_bse ) * c1i
            Sig_g = Sig_g + aimag( Sig_g_bse ) * c1i
            ! get leads sigma        
            siglead(:,:,:,1) = Sig_r(1:NB*NS,1:NB*NS,:)
            siglead(:,:,:,2) = Sig_r(nm_dev-NB*NS+1:nm_dev,nm_dev-NB*NS+1:nm_dev,:)    
            !
        endif
        print *, 'calc G last time ...'  
        !
        call calc_gf(nen,En,2,nm_dev,[nb*ns,nb*ns],nb*ns,&
                    Ham(:,:),H00lead(:,:,:),H10lead(:,:,:),Siglead(:,:,:,:),&
                    T(:,:,:),Sig_r(:,:,:),Sig_l(:,:,:),Sig_g(:,:,:),&
                    G_retarded(:,:,:),G_lesser(:,:,:),&
                    G_greater(:,:,:),current,Te,mu,temp,flatband)    
        !                                      
        transmission(:,1)=te(:,1,2)
        transmission(:,2)=te(:,2,1)     
        iter=niter
        if (output_files) then
            call calc_bond_current(Ham,G_lesser,nen,en,spindeg,nm_dev,tot_cur,tot_ecur,cur)
            call write_current_spectrum('bse_Jdens',iter,cur,nen,en,length,NB,Lx)
            call write_current('bse_I',iter,tot_cur,length,NB,NS,Lx)
            call write_current('bse_EI',iter,tot_ecur,length,NB,NS,Lx)
            call write_spectrum_nosub('bse_ldos',iter,G_retarded,nen,En,length,NB,Lx,(/1.0d0,-2.0d0/))
            call write_spectrum_nosub('bse_ndos',iter,G_lesser,nen,En,length,NB,Lx,(/1.0d0,1.0d0/))
            call write_spectrum_nosub('bse_pdos',iter,G_greater,nen,En,length,NB,Lx,(/1.0d0,-1.0d0/))                    
            call write_transmission_spectrum('bse_trL',iter,current(:,1)*spindeg,nen,En)
            call write_transmission_spectrum('bse_trR',iter,current(:,2)*spindeg,nen,En)
            call write_transmission_spectrum('bse_TE_LR',iter,Te(:,1,2)*spindeg,nen,En)
            call write_transmission_spectrum('bse_TE_RL',iter,Te(:,2,1)*spindeg,nen,En)    
        endif
        !
    end subroutine bse_sparse_solve_scba
  
    ! build the Bethe-Salpeter Equation system matrix to invert from L0 and Kernel matrices
    !   A = ( I - L0 @ K )
    subroutine bse_sparse_build_system(blocksize,num_blocks,nm_dev,nnop,iop,&
        Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip, Kdiag, Ktip,&
        Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip)
        ! input
        integer,intent(in)::blocksize,num_blocks,nm_dev,nnop,iop
        complex(dp),intent(in),dimension(:,:,:):: Ldiag
        complex(dp),intent(in),dimension(:,:,:):: Lupper,Llower ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(in),dimension(:,:,:):: Llowerarrow
        complex(dp),intent(in),dimension(:,:,:):: Lupperarrow ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(in),dimension(:,:,:):: Ltip ! dense tip block of 2-point polarization function with interacting electron-hole at frequency [[nop]]                        
        complex(dp),intent(in),dimension(:):: Kdiag ! diagonal of Kernel
        complex(dp),intent(in),dimension(:,:):: Ktip ! dense tip block of Kernel
        ! output
        complex(dp),intent(out),dimension(blocksize,blocksize*num_blocks):: Adiag
        complex(dp),intent(out),dimension(blocksize,blocksize*(num_blocks-1)):: Aupper,Alower 
        complex(dp),intent(out),dimension(nm_dev,blocksize*num_blocks):: Alowerarrow
        complex(dp),intent(out),dimension(blocksize*num_blocks,nm_dev):: Aupperarrow 
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: Atip 
        ! --- local
        integer:: ib,i,j,N 
        if ((iop>0).and.(iop<=nnop)) then
            N = blocksize*num_blocks
            ! A_xx
            call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Ltip(:,:,iop),nm_dev,Ktip,nm_dev,czero,Atip,nm_dev)  
            !$omp parallel default(shared) private(i)
            !$omp do
            do i=1,nm_dev
                Atip(i,i) = Atip(i,i) + cone 
            enddo 
            !$omp end do
            !$omp end parallel  
            ! A_xd = - L_xd * K_dd
            !$omp parallel default(shared) private(i,j)
            !$omp do   
            do i = 1,nm_dev
                do concurrent (j = 1:N) 
                    Alowerarrow(i,j) = - Llowerarrow(i,j,iop) * Kdiag(j)
                enddo
            enddo                   
            !$omp end do
            !$omp end parallel  
            ! A_dx = - L_dx * K_xx            
            call zgemm('n','n',N,nm_dev,nm_dev,-cone,Lupperarrow(:,:,iop),N,Ktip,nm_dev,czero,Aupperarrow,N)
            ! A_dd
            ! diagonal blocks         
            !$omp parallel default(shared) private(i,j)
            !$omp do  
            do i=1,blocksize
                do concurrent (j=1:blocksize*num_blocks)
                    Adiag(i,j) = - Ldiag(i,j,iop) * Kdiag(j)
                enddo               
            enddo                        
            !$omp end do
            !$omp end parallel  
            !$omp parallel default(shared) private(ib,i)
            !$omp do  
            do ib=1,num_blocks
                do concurrent (i=1:blocksize)
                    Adiag(i,i+(ib-1)*blocksize) = Adiag(i,i+(ib-1)*blocksize) + cone
                enddo 
            enddo                      
            !$omp end do
            !$omp end parallel 
            ! upper and lower diagonal blocks
            !$omp parallel default(shared) private(i,j)
            !$omp do  
            do i=1,blocksize
                do concurrent (j=1:blocksize*(num_blocks-1))
                    Aupper(i,j) = - Lupper(i,j,iop) * Kdiag(j+blocksize)
                    Alower(i,j) = - Llower(i,j,iop) * Kdiag(j)
                enddo 
            enddo                       
            !$omp end do
            !$omp end parallel 
            !
        else 
            print *,"iop not correct!"
            call abort
        endif
    end subroutine bse_sparse_build_system

    ! check the sparse BTA system matrix blocks against dense matrix
    subroutine bse_sparse_check_system( & 
        tol,alpha,spindeg,nm_dev,ndiag,nen,En,nop,nnop,iop,blocksize,num_blocks,N,table,&
        G_lesser,G_greater,G_retarded,W,V, &
        Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip, &
        Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip &
        )
        ! input 
        integer,intent(in)::nm_dev,nen,nnop,nop(nnop),ndiag,N,table(2,N)
        integer,intent(in)::blocksize, num_blocks, iop ! arrow block size and number of blocks (excluding tip block)
        real(dp),intent(in)::en(nen),spindeg,alpha                
        real(dp),intent(in)::tol
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen):: G_lesser,G_greater,G_retarded ! electron GFs    
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction    
        complex(dp),intent(in),dimension(:,:):: Adiag,Aupper,Alower,Alowerarrow,Aupperarrow,Atip 
        complex(dp),intent(in),dimension(:,:,:):: Ldiag,Lupper,Llower,Llowerarrow,Lupperarrow,Ltip  
        ! ---- local 
        complex(dp),dimension(:,:),allocatable :: Lmat,Amat ! two-particle Green's function 
        complex(dp),dimension(:,:),allocatable :: Mmat ! 4-point Kernel
        integer :: i,j,k,l,row,col,NT,ib,p,q 
        complex(dp) :: L0ijkl,tmp
        real(dp) :: error        
        !
        NT = blocksize * num_blocks + nm_dev
        print *, " start checking the system matrix blocks"
        allocate(Mmat(NT,NT), source=czero)        
        allocate(Lmat(NT,NT), source=czero)     
        allocate(Amat(NT,NT), source=czero)  
        !           
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
                    (abs(i-l)<=ndiag).and.(abs(i-j)<=ndiag).and.(abs(k-l)<=ndiag)) then 
                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                        G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                    Lmat(NT-col+1,NT-row+1) = L0ijkl * spindeg
                endif
            enddo
        enddo
        !$omp end do
        !$omp end parallel 
        !
        ! set L to zero for elements outside of arrowhead structure
        do row = 1, blocksize*num_blocks
            do col = 1, blocksize*num_blocks
                ib = (row-1) / blocksize
                p = ib * blocksize 
                if (col > (p + blocksize * 2).and.(col <= (p - blocksize))) then 
                    Lmat(row,col) = czero
                endif
            enddo 
        enddo 
        !$omp parallel default(shared) private(row,col,i,j,k,l)
        !$omp do        
        do row=1,N                        
            do col=1,N
                i=table(1,row)
                j=table(2,row)
                k=table(1,col)
                l=table(2,col)           
                if ((i==j).and.(k==l)) then                        
                    Mmat(NT-col+1,NT-row+1) = Mmat(NT-col+1,NT-row+1) - c1i *  V(i,k) * spindeg                        
                endif 
                if ((i==k).and.(j==l)) then                        
                    Mmat(NT-col+1,NT-row+1) = Mmat(NT-col+1,NT-row+1) + c1i *  W(i,j)
                endif 
            enddo
        enddo    
        !$omp end do
        !$omp end parallel 
        !        
        call zgemm('n','n',NT,NT,NT,-cone,Lmat,NT,Mmat,NT,czero,Amat,NT)         
        !
        ! (I - L0 K) -> A
        !$omp parallel default(shared) private(i)
        !$omp do  
        do i=1,NT 
            Amat(i,i) = Amat(i,i) + cone
        enddo  
        !$omp end do
        !$omp end parallel
        !
        ! check tip block   
        do row=NT-nm_dev+1,NT
            do col=NT-nm_dev+1,NT
                error = abs(Lmat(row,col) - Ltip(row - (NT-nm_dev),col - (NT-nm_dev),iop)) 
                if (error > tol) then 
                    print *, "ERROR in L tip block", row, col, error  
                    call abort
                endif                 
            enddo
        enddo
        do row=NT-nm_dev+1,NT
            do col=NT-nm_dev+1,NT
                error = abs(Amat(row,col) - Atip(row - (NT-nm_dev),col - (NT-nm_dev))) 
                if (error > tol) then 
                    print *, "ERROR in A tip block", row, col, error  
                    call abort
                endif 
            enddo
        enddo
        ! check diagonal blocks
        do ib=1,num_blocks 
            do i=1,blocksize
                do j=1,blocksize
                    row = (ib-1) * blocksize + i 
                    col = (ib-1) * blocksize + j 
                    error = abs(Lmat(row,col) - Ldiag(i,j+(ib-1)*blocksize,iop))
                    if (error > tol) then 
                        print *, "ERROR in L diagonal block", ib, i, j, row, col, error 
                        call abort
                    endif
                    error = abs(Amat(row,col) - Adiag(i,j+(ib-1)*blocksize))
                    if (error > tol) then 
                        print *, "ERROR in A diagonal block", ib, i, j, row, col, error 
                        call abort
                    endif
                enddo 
            enddo 
        enddo 
        ! check upper/lower diagonal blocks
        do ib=1, num_blocks-1
            do i=1, blocksize
                do j=1, blocksize
                    row = (ib-1) * blocksize + i 
                    col = ib * blocksize + j 
                    error = abs(Lmat(row,col) - Lupper(i,col-blocksize,iop))
                    if (error > tol) then 
                        print *, "ERROR in L upper diagonal block", ib, i, j, row, col, error       
                        call abort                      
                    endif
                    error = abs(Amat(row,col) - Aupper(i,col-blocksize))
                    if (error > tol) then 
                        print *, "ERROR in A upper diagonal block", ib, i, j, row, col, error    
                        call abort                         
                    endif           
                    row = ib * blocksize + i         
                    col = (ib-1) * blocksize + j
                    error = abs(Lmat(row,col) - Llower(i,col,iop))
                    if (error > tol) then 
                        print *, "ERROR in L lower diagonal block", ib, i, j, row, col, error 
                        call abort                            
                    endif 
                    error = abs(Amat(row,col) - Alower(i,col))
                    if (error > tol) then 
                        print *, "ERROR in A lower diagonal block", ib, i, j, row, col, error         
                        call abort                    
                    endif 
                enddo 
            enddo 
        enddo 
        ! check arrow upper/lower blocks 
        do ib=1,num_blocks 
            do i=1,blocksize
                do j=1,nm_dev
                    col = NT - nm_dev + j 
                    row = (ib-1)*blocksize + i 
                    !
                    error = abs(Lmat(row,col) - Lupperarrow(row,j,iop))
                    if (error > tol) then 
                        print *, "ERROR in L arrow upper blocks", ib, i, j, row, col , error                        
                        call abort
                    endif
                    error = abs(Amat(row,col) - Aupperarrow(i+(ib-1)*blocksize,j))
                    if (error > tol) then 
                        print *, "ERROR in A arrow upper blocks", ib, i, j, row, col , error    
                        print *, Amat(row,col) , Aupperarrow(i+(ib-1)*blocksize,j)
                        call abort                    
                    endif
                    !
                    col = (ib-1)*blocksize + i 
                    row = NT - nm_dev + j
                    !
                    error = abs(Lmat(row,col) - Llowerarrow(j, i+(ib-1)*blocksize,iop))
                    if (error > tol) then 
                        print *, "ERROR in L arrow lower blocks", ib, i, j, row, col , error
                        print *, Lmat(row,col) , Llowerarrow(j, i+(ib-1)*blocksize,iop)
                        call abort
                    endif
                    error = abs(Amat(row,col) - Alowerarrow(j, i+(ib-1)*blocksize))
                    if (error > tol) then 
                        print *, "ERROR in A arrow lower blocks", ib, i, j, row, col , error
                        print *, Amat(row,col) , Alowerarrow(j, i+(ib-1)*blocksize)
                        call abort
                    endif
                    !
                enddo 
            enddo 
        enddo 
        print *, " DONE CHECK, all checks pass "
    end subroutine bse_sparse_check_system

    ! build the Bethe-Salpeter Equation L0 and Kernel matrices
    subroutine bse_sparse_build(method,alpha,spindeg,nm_dev,ndiag,nen,En,nop,nnop,blocksize,num_blocks,N,table,&
        G_lesser,G_greater,G_retarded,W,V,&
        Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag)      
        !
        ! use gpu_polarization  
        ! input
        character(len=*),intent(in)::method
        integer,intent(in)::nm_dev,nen,nnop,nop(nnop),ndiag,N,table(2,N)
        integer,intent(in)::blocksize, num_blocks ! arrow block size and number of blocks (excluding tip block)
        real(dp),intent(in)::en(nen),spindeg,alpha                
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen):: G_lesser,G_greater,G_retarded ! electron GFs    
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction    
        ! output
        complex(dp),intent(out),dimension(blocksize,blocksize*num_blocks,nnop):: Ldiag
        complex(dp),intent(out),dimension(blocksize,blocksize*(num_blocks-1),nnop):: Lupper,Llower ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(out),dimension(nm_dev,blocksize*num_blocks,nnop):: Llowerarrow
        complex(dp),intent(out),dimension(blocksize*num_blocks,nm_dev,nnop):: Lupperarrow ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nnop):: Ltip ! dense tip block of 2-point polarization function with interacting electron-hole at frequency [[nop]]                        
        complex(dp),intent(out),dimension(blocksize*num_blocks):: Kdiag ! diagonal of Kernel
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: Ktip ! dense tip block of Kernel
        ! ---- local
        complex(dp) :: L0ijkl(nnop), P_ijkl              
        real(dp) :: start, finish
        complex(dp) :: a1, a2 
        integer :: i,j,k,l,p,q,ie,row,col,it,iop,ib,NT,nepoch,fliped_row,fliped_col     
        integer(8) :: devPtrA, devPtrB, devPtrC, devPtrGL,devPtrGG,devPtrGR,devPtrGA           
        !
        a1 = 1.0 - alpha
        a2 = alpha * 0.5
        !
        Ltip = czero
        Ldiag = czero
        Lupper = czero
        Llower = czero
        Lupperarrow = czero
        Llowerarrow = czero
        NT = blocksize * num_blocks   
        !
        if (trim(method)=='gpu_sum') then
            ! print *, "gpu polarization prepare"
            ! call gpu_polarization_prepare(&
            !     nen=nen,nm=nm_dev,&
            !     G_lesser=g_lesser,G_greater=g_greater,G_retarded=g_retarded,&
            !     devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA)
        endif
        !
        nepoch = N / 50
        start = omp_get_wtime()              
        write(*, '(A)', advance="no") '  start computation L0_ijkl = G_jl G_ki '                 
        !$omp parallel default(shared) private(row,col,i,j,k,l,L0ijkl,ib,p,q,fliped_row,fliped_col,devPtrA, devPtrB, devPtrC)
        ! call gpu_polarization_calc_ijkl(&
        !                                 a1=a1,a2=a2,nop=nop(1),nen=nen,nm=nm_dev,i=1,j=1,k=1,l=1,&
        !                                 devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,&
        !                                 devPtrA=devPtrA,devPtrB=devPtrB,devPtrC=devPtrC,&
        !                                 P_ijkl = P_ijkl, malloc_gpu=.true.)
        !$omp do        
        do row = 1,N 
            if (mod(row-1,nepoch)==0) write(*, '(A)', advance="no") '.'                 
            do col = 1,N         
                i=table(1,row)
                j=table(2,row)
                k=table(1,col)
                l=table(2,col)    
                if ((abs(i-k)<=ndiag).and.(abs(j-l)<=ndiag).and.(abs(j-k)<=ndiag).and.&
                    (abs(i-l)<=ndiag).and.(abs(i-j)<=ndiag).and.(abs(k-l)<=ndiag)) then                   
                    ! need to flip the row and col when putting into arrowhead structure               
                    fliped_row = NT + nm_dev - col + 1
                    fliped_col = NT + nm_dev - row + 1
                    if (fliped_col > NT) then 
                        if (fliped_row > NT) then 
                            ! tip block
                            select case (trim(method))
                            case('fft')
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)

                            case('gpu_sum') 
                                ! do iop=1,nnop
                                !     call gpu_polarization_calc_ijkl(&
                                !         a1=a1,a2=a2,nop=nop(iop),nen=nen,nm=nm_dev,i=i,j=j,k=k,l=l,&
                                !         devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,&
                                !         devPtrA=devPtrA,devPtrB=devPtrB,devPtrC=devPtrC,&
                                !         P_ijkl = P_ijkl, malloc_gpu=.false.)
                                !     L0ijkl(iop) = P_ijkl * ( En(2) - En(1) ) / twopi
                                ! enddo
                                
                            case default    
                                do iop=1,nnop
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo

                            end select         
                            Ltip(fliped_row - NT,fliped_col - NT,1:nnop) = L0ijkl * spindeg
                        else
                            ! upper arrow block 
                            select case (trim(method))
                            case('fft')
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)

                            case('gpu_sum') 
                                ! do iop=1,nnop
                                !     call gpu_polarization_calc_ijkl(&
                                !         a1=a1,a2=a2,nop=nop(iop),nen=nen,nm=nm_dev,i=i,j=j,k=k,l=l,&
                                !         devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,&
                                !         devPtrA=devPtrA,devPtrB=devPtrB,devPtrC=devPtrC,&
                                !         P_ijkl = P_ijkl, malloc_gpu=.false.)
                                !     L0ijkl(iop) = P_ijkl * ( En(2) - En(1) ) / twopi
                                ! enddo

                            case default    
                                do iop=1,nnop
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo

                            end select 
                            Lupperarrow(fliped_row,fliped_col - NT,1:nnop) = L0ijkl * spindeg          
                        endif 
                    else 
                        if (fliped_row > NT) then 
                            ! lower arrow block 
                            select case (trim(method))
                            case('fft')
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)

                            case('gpu_sum') 
                                ! do iop=1,nnop
                                !     call gpu_polarization_calc_ijkl(&
                                !         a1=a1,a2=a2,nop=nop(iop),nen=nen,nm=nm_dev,i=i,j=j,k=k,l=l,&
                                !         devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,&
                                !         devPtrA=devPtrA,devPtrB=devPtrB,devPtrC=devPtrC,&
                                !         P_ijkl = P_ijkl, malloc_gpu=.false.)
                                !     L0ijkl(iop) = P_ijkl * ( En(2) - En(1) ) / twopi
                                ! enddo
                                               
                            case default    
                                do iop=1,nnop
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo

                            end select
                            Llowerarrow(fliped_row - NT,fliped_col,1:nnop) = L0ijkl * spindeg   
                        else 
                            ib = (fliped_row-1) / blocksize
                            p = ib * blocksize 
                            q = p + blocksize
                            if ((fliped_col > p).and.(fliped_col <= q)) then 
                                ! diag block 
                                select case (trim(method))
                                case('fft')
                                    call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)

                                case('gpu_sum') 
                                    ! do iop=1,nnop
                                    !     call gpu_polarization_calc_ijkl(&
                                    !         a1=a1,a2=a2,nop=nop(iop),nen=nen,nm=nm_dev,i=i,j=j,k=k,l=l,&
                                    !         devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,&
                                    !         devPtrA=devPtrA,devPtrB=devPtrB,devPtrC=devPtrC,&
                                    !         P_ijkl = P_ijkl, malloc_gpu=.false.)
                                    !     L0ijkl(iop) = P_ijkl * ( En(2) - En(1) ) / twopi
                                    ! enddo
                                                    
                                case default    
                                    do iop=1,nnop
                                        call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                    enddo

                                end select
                                Ldiag(fliped_row - p, fliped_col, 1:nnop) = L0ijkl * spindeg   
                            else
                                if ((fliped_col > q).and.(fliped_col <= (q+blocksize))) then                             
                                    ! upper diag block 
                                    select case (trim(method))
                                    case('fft')
                                        call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                                    
                                    case('gpu_sum') 
                                        ! do iop=1,nnop
                                        !     call gpu_polarization_calc_ijkl(&
                                        !         a1=a1,a2=a2,nop=nop(iop),nen=nen,nm=nm_dev,i=i,j=j,k=k,l=l,&
                                        !         devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,&
                                        !         devPtrA=devPtrA,devPtrB=devPtrB,devPtrC=devPtrC,&
                                        !         P_ijkl = P_ijkl, malloc_gpu=.false.)
                                        !     L0ijkl(iop) = P_ijkl * ( En(2) - En(1) ) / twopi
                                        ! enddo

                                    case default    
                                        do iop=1,nnop
                                            call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                        enddo

                                    end select
                                    Lupper(fliped_row - p, fliped_col - blocksize,1:nnop) = L0ijkl * spindeg                               
                                endif
                                if ((fliped_col > (p-blocksize)).and.(fliped_col <= p)) then                             
                                    ! lower diag block
                                    select case (trim(method))
                                    case('fft')
                                        call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)

                                    case('gpu_sum') 
                                        ! do iop=1,nnop
                                        !     call gpu_polarization_calc_ijkl(&
                                        !         a1=a1,a2=a2,nop=nop(iop),nen=nen,nm=nm_dev,i=i,j=j,k=k,l=l,&
                                        !         devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA,&
                                        !         devPtrA=devPtrA,devPtrB=devPtrB,devPtrC=devPtrC,&
                                        !         P_ijkl = P_ijkl, malloc_gpu=.false.)
                                        !     L0ijkl(iop) = P_ijkl * ( En(2) - En(1) ) / twopi
                                        ! enddo

                                    case default    
                                        do iop=1,nnop
                                            call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                        enddo

                                    end select
                                    Llower(fliped_row - p, fliped_col,1:nnop ) = L0ijkl * spindeg                            
                                endif 
                            endif 
                        endif 
                    endif
                endif 
            enddo
        enddo         
        !$omp end do        
        ! call cublas_free(devPtrA)
        ! call cublas_free(devPtrB)
        ! call cublas_free(devPtrC)
        !$omp end parallel        
        if (trim(method)=='gpu_sum') then
            print *, "gpu polarization finish, free gpu memory"
            ! call gpu_polarization_finalize(&
            !     devPtrGL=devPtrGL,devPtrGG=devPtrGG,devPtrGR=devPtrGR,devPtrGA=devPtrGA)
        endif
        !
        Ktip=czero
        Kdiag=czero
        !$omp parallel default(shared) private(row,col,i,j,k,l,fliped_row,fliped_col)
        !$omp do        
        do row=1,N                        
            do col=1,N
                i=table(1,row)
                j=table(2,row)
                k=table(1,col)
                l=table(2,col)           
                fliped_row = NT + nm_dev - col + 1
                fliped_col = NT + nm_dev - row + 1           
                if ((i==j).and.(k==l)) then                               
                    Ktip(fliped_row - NT,fliped_col - NT) = Ktip(fliped_row - NT,fliped_col - NT) - c1i *  V(i,k) * spindeg        
                    if ((i==k).and.(j==l).and.(row<=nm_dev)) then                
                        Ktip(fliped_row - NT,fliped_row - NT) = Ktip(fliped_row - NT,fliped_row - NT) + c1i *  W(i,j)
                    endif
                endif 
                if ((i==k).and.(j==l).and.(row>nm_dev)) then                                                             
                    Kdiag(fliped_row) = Kdiag(fliped_row) + c1i *  W(i,j)
                endif 
            enddo
        enddo    
        !$omp end do
        !$omp end parallel 
        finish = omp_get_wtime()
        print *
        print '("  L0 and K computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
    end subroutine bse_sparse_build

end module bse_sparse
