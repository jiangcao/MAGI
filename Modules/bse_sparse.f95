
module bse_sparse    
    use parameters_mod,only:dp,twopi,pi,e_charge,epsilon0,m0_charge,hbar,c1i,czero,cone    
    use omp_lib
    use polarization
    use sinv
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
        out_Adiag = reshape(Adiag, [num_blocks,blocksize,blocksize], order=[3,1,2])
        out_Aupper = reshape(Aupper, [num_blocks-1,blocksize,blocksize], order=[3,1,2])
        out_Alower = reshape(Alower, [num_blocks-1,blocksize,blocksize], order=[3,1,2])
        out_Aupperarrow = reshape(Aupperarrow, [num_blocks,blocksize,nm_dev], order=[2,1,3])
        out_Alowerarrow = reshape(Alowerarrow, [num_blocks,nm_dev,blocksize], order=[3,1,2])
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
    subroutine bse_sparse_solve(method,alpha,spindeg,nm_dev,ndiag,nen,nsub,En,nops,nnop,nk,G_lesser,G_greater,G_retarded,W,V,P_retarded)        
        ! in
        character(len=*),intent(in)::method
        integer,intent(in)::nm_dev,nen ! device dimension, number of energies
        integer,intent(in)::nnop,nops(nnop) ! number of optical energies, optical energies in unit of energy interval
        integer,intent(in)::nsub,nk,ndiag ! number of legendre sub-energy nodes, number of k points, number of offdiagonals
        real(dp),intent(in)::en(nen),spindeg,alpha ! energy grid, spin degeneracy, empirical parameter
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen,nsub,nk):: G_lesser,G_greater,G_retarded ! electron Green Functions
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb operator
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb operator
        ! out 
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nnop):: P_retarded ! 2-point polarization function with interacting electron-hole
        !---- local
        complex(dp),allocatable,dimension(:,:):: Adiag,Aupper,Alower,Alowerarrow,Aupperarrow,Atip 
        complex(dp),allocatable,dimension(:,:):: Ktip 
        complex(dp),allocatable,dimension(:)  :: Kdiag
        complex(dp),allocatable,dimension(:,:,:):: Ldiag,Lupper,Llower,Llowerarrow,Lupperarrow,Ltip 
        integer::N,blocksize,num_blocks, NT, local_nnop, iop, row, col, fliped_col, fliped_row, i, k
        integer(8):: nnz
        integer,allocatable,dimension(:,:)::table
        integer,allocatable,dimension(:,:)::ipiv_diagonal
        integer,allocatable,dimension(:)::ipiv_arrow_tip
        real(dp) :: start, finish
        !
        ! pre-process the sparsity pattern of system
        print *, " pre-process ... "
        allocate(table(2,nm_dev*nm_dev))
        call bse_sparse_pre(nm_dev,ndiag,N,nnz,table,blocksize,num_blocks)
        ! prepare the memory
        NT = blocksize * num_blocks
        if (trim(method) == 'batched') then 
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
        !
        if (local_nnop > 1) then  
            ! build BTA blocks of RPA polarization L0 and 2-body interaction kernal K                     
            call bse_sparse_build(alpha,spindeg,nm_dev,ndiag,nen,En,nops,local_nnop,blocksize,num_blocks,N,table,&
                                    G_lesser,G_greater,G_retarded,W,V,&
                                    Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag)       
            print *, " start selected inversion " 
            start = omp_get_wtime()                                     
            do iop = 1,local_nnop
                write(*, '(A)', advance="no") '.'
                ! build system matrix blocks (I - L0 @ K)
                ! print *, "  build system "
                call bse_sparse_build_system(blocksize,num_blocks,nm_dev,local_nnop,iop,&
                                Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip, &
                                Kdiag, Ktip,&
                                Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip)
                ! selected inversion of the system matrix
                ! print *, "  factorize "
                call zbtatrf( blocksize, nm_dev, num_blocks, &
                            Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                            ipiv_diagonal,ipiv_arrow_tip)
                ! print *, "  invert "
                call zbtatri( blocksize, nm_dev, num_blocks, &
                            Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                            ipiv_diagonal,ipiv_arrow_tip)
                ! call zbtasinv(blocksize, nm_dev, num_blocks, &
                !             Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip)
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
                ! build BTA blocks of RPA polarization L0 and 2-body interaction kernal K         
                call bse_sparse_build(alpha,spindeg,nm_dev,ndiag,nen,En,nops(iop),1,blocksize,num_blocks,N,table,&
                                    G_lesser,G_greater,G_retarded,W,V,&
                                    Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag)                   
                ! build system matrix blocks (I - L0 @ K)
                call bse_sparse_build_system(blocksize,num_blocks,nm_dev,1,1,&
                                Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip, &
                                Kdiag, Ktip,&
                                Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip)
                ! selected inversion of the system matrix
                call zbtatrf( blocksize, nm_dev, num_blocks, &
                            Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                            ipiv_diagonal,ipiv_arrow_tip)
                call zbtatri( blocksize, nm_dev, num_blocks, &
                            Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip, &
                            ipiv_diagonal,ipiv_arrow_tip)
                ! call zbtasinv(blocksize, nm_dev, num_blocks, &
                !             Adiag,Alower,Aupper,Alowerarrow,Aupperarrow,Atip)
                ! compute P_retarded
                ! 
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
                    enddo
                enddo                          
                !$omp end do
                !$omp end parallel  
            enddo
        endif 
    end subroutine bse_sparse_solve

  
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
        print *, " check the system matrix"
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
        print *, "DONE CHECK"
    end subroutine bse_sparse_check_system


    ! build the Bethe-Salpeter Equation L0 and Kernel matrices
    subroutine bse_sparse_build(alpha,spindeg,nm_dev,ndiag,nen,En,nop,nnop,blocksize,num_blocks,N,table,&
        G_lesser,G_greater,G_retarded,W,V,&
        Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag)        
        ! input
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
        complex(dp) :: L0ijkl(nnop)              
        real(dp) :: start, finish
        integer :: i,j,k,l,p,q,ie,row,col,it,iop,ib,NT,nepoch,fliped_row,fliped_col                
        !
        Ltip = czero
        Ldiag = czero
        Lupper = czero
        Llower = czero
        Lupperarrow = czero
        Llowerarrow = czero
        NT = blocksize * num_blocks   
        !
        nepoch = N / 400
        start = omp_get_wtime()              
        print *,'  start computation L0_ijkl = G_jl G_ki ...'                 
        !$omp parallel default(shared) private(row,col,i,j,k,l,L0ijkl,ib,p,q,fliped_row,fliped_col)
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
                            if (nnop>10) then
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                            else    
                                do iop=1,nnop
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo
                            endif         
                            Ltip(fliped_row - NT,fliped_col - NT,1:nnop) = L0ijkl * spindeg
                        else
                            ! upper arrow block 
                            if (nnop>10) then
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                            else    
                                do iop=1,nnop
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo
                            endif 
                            Lupperarrow(fliped_row,fliped_col - NT,1:nnop) = L0ijkl * spindeg          
                        endif 
                    else 
                        if (fliped_row > NT) then 
                            ! lower arrow block 
                            if (nnop>10) then
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                            else    
                                do iop=1,nnop
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo
                            endif 
                            Llowerarrow(fliped_row - NT,fliped_col,1:nnop) = L0ijkl * spindeg   
                        else 
                            ib = (fliped_row-1) / blocksize
                            p = ib * blocksize 
                            q = p + blocksize
                            if ((fliped_col > p).and.(fliped_col <= q)) then 
                                ! diag block 
                                if (nnop>10) then
                                    call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                                else    
                                    do iop=1,nnop
                                        call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                    enddo
                                endif 
                                Ldiag(fliped_row - p, fliped_col, 1:nnop) = L0ijkl * spindeg   
                            else
                                if ((fliped_col > q).and.(fliped_col <= (q+blocksize))) then                             
                                    ! upper diag block 
                                    if (nnop>10) then
                                        call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                                    else    
                                        do iop=1,nnop
                                            call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                        enddo
                                    endif 
                                    Lupper(fliped_row - p, fliped_col - blocksize,1:nnop) = L0ijkl * spindeg                               
                                endif
                                if ((fliped_col > (p-blocksize)).and.(fliped_col <= p)) then                             
                                    ! lower diag block
                                    if (nnop>10) then
                                        call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                                    else    
                                        do iop=1,nnop
                                            call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                        enddo
                                    endif 
                                    Llower(fliped_row - p, fliped_col,1:nnop ) = L0ijkl * spindeg                            
                                endif 
                            endif 
                        endif 
                    endif
                endif 
            enddo
        enddo         
        !$omp end do
        !$omp end parallel 
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
        print '("  computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
    end subroutine bse_sparse_build

end module bse_sparse
