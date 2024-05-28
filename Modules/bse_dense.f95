
module bse_dense    
    use parameters_mod,only:dp,twopi,pi,e_charge,epsilon0,m0_charge,hbar,c1i,czero,cone
    use legendre
    use omp_lib
    implicit none
    contains



    subroutine four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,G_lesser,G_greater,G_retarded,i,j,k,l,L0)
        use fft_mod, only : corr1d => corr1d2  
        integer,intent(in) :: nm_dev,nen,nnop,nop(nnop),ndiag, i,j,k,l
        real(dp),intent(in) :: en(nen), alpha 
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen),target :: G_lesser,G_greater,G_retarded
        complex(dp),intent(out) :: L0(nnop)
        ! ---
        complex(dp),dimension(nen) :: Gl,Gg,Gr
        complex(dp),dimension(nen) :: Gl_down,Gg_down,Gr_down
        real(dp) :: dE, weights, xen, a1,a2
        integer :: ie, isub, ik, ikd
        complex(dp),dimension(nen) :: tmp
        ! the P4 IPA tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
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
        Gr_down(1:nen) = G_retarded(k,i,1:nen)
        ! calculate P4_IPA from GG
        tmp = corr1d(nen,Gl,conjg(Gr_down),method='fft') * a1
        tmp = tmp  + corr1d(nen,Gr,Gl_down,method='fft') * a1
        tmp = tmp  + corr1d(nen,Gg,Gl_down,method='fft') * a2
        tmp = tmp  - corr1d(nen,Gl,Gg_down,method='fft') * a2
        L0(1:nnop) = tmp(nop(1:nnop)+nen/2)        
        !
    end subroutine four_polarization_fft


    subroutine bse_sparse_pre(nm_dev,ndiag,N,nnz,table,blocksize,num_blocks)
        integer,intent(in)::nm_dev
        integer,intent(in)::ndiag 
        integer,intent(out)::table(2,nm_dev*nm_dev) ! lookup table from 2-body to 1-body            
        integer,intent(out)::blocksize,num_blocks
        integer,intent(out)::N! size of the reduced 2-body system
        integer(8),intent(out)::nnz! nonzero       
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
            print *, 'ERROR!'
            call abort
        endif        
        print *, 'nm_dev=', nm_dev
        print *, 'resized system size=', N 
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
        blocksize = bandwidth/2 
        num_blocks = ceiling( dble(N - nm_dev) / blocksize )  
        NT = blocksize * num_blocks         
        print '("  total arrow size=", I20)', NT
        print '("  arrow block size=", I20)', blocksize
        print '("  arrow number of blocks=", I20)', num_blocks
        print '("  nonzero elements=", F0.3, " Million")', dble(nnz)/1e6
        print '("  nonzero ratio = ", F0.3 ," %")', dble(nnz)/(dble(NT+nm_dev)**2)*100
    end subroutine bse_sparse_pre

    ! build the Bethe-Salpeter Equation system matrix to invert from L0 and Kernel matrices
    !   A = ( I - L0 @ K )
    subroutine bse_sparse_build_system(blocksize,num_blocks,nm_dev,&
        Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip, Kdiag, Ktip,&
        Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip)
        integer,intent(in)::blocksize,num_blocks,nm_dev
        complex(dp),intent(in),dimension(blocksize,blocksize*num_blocks):: Ldiag
        complex(dp),intent(in),dimension(blocksize,blocksize*(num_blocks-1)):: Lupper,Llower ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(in),dimension(nm_dev,blocksize*num_blocks):: Llowerarrow
        complex(dp),intent(in),dimension(blocksize*num_blocks,nm_dev):: Lupperarrow ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(in),dimension(nm_dev,nm_dev):: Ltip ! dense tip block of 2-point polarization function with interacting electron-hole at frequency [[nop]]                        
        complex(dp),intent(in),dimension(blocksize*num_blocks+nm_dev):: Kdiag ! diagonal of Kernel
        complex(dp),intent(in),dimension(nm_dev,nm_dev):: Ktip ! dense tip block of Kernel
        !
        complex(dp),intent(out),dimension(blocksize,blocksize*num_blocks):: Adiag
        complex(dp),intent(out),dimension(blocksize,blocksize*(num_blocks-1)):: Aupper,Alower 
        complex(dp),intent(out),dimension(nm_dev,blocksize*num_blocks):: Alowerarrow
        complex(dp),intent(out),dimension(blocksize*num_blocks,nm_dev):: Aupperarrow 
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: Atip 
        ! --- 
        integer:: ib,i,j,N 
        N = blocksize*num_blocks
        ! A_xx
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Ltip,nm_dev,Ktip,nm_dev,czero,Atip,nm_dev)  
        do concurrent (i=1:nm_dev)
            Atip(i,i) = Atip(i,i) + cone 
        enddo 
        ! A_xd = - L_xd * K_dd
        do concurrent (i = 1:nm_dev)
            do concurrent (j = 1:N) 
                Alowerarrow(i,j) = - Llowerarrow(i,j) * Kdiag(j)
            enddo
        enddo        
        ! A_dx = - L_dx * K_xx
        call zgemm('n','n',N,nm_dev,nm_dev,-cone,Lupperarrow,N,Ktip,nm_dev,czero,Aupperarrow,N)
        ! A_dd
        ! diagonal blocks         
        do concurrent (i=1:blocksize)
            do concurrent (j=1:blocksize*num_blocks)
                Adiag(i,j) = - Ldiag(i,j) * Kdiag(j)
            enddo               
        enddo 
        do concurrent (ib=1:num_blocks)
            do concurrent (i=1:blocksize)
                Adiag(i,i+(ib-1)*blocksize) = Adiag(i,i+(ib-1)*blocksize) + cone
            enddo 
        enddo
        ! upper and lower diagonal blocks
        do concurrent (i=1:blocksize)
            do concurrent (j=1:blocksize*(num_blocks-1))
                Aupper(i,j) = - Lupper(i,j) * Kdiag(j+blocksize)
                Alower(i,j) = - Llower(i,j) * Kdiag(j)
            enddo 
        enddo 
    end subroutine bse_sparse_build_system

    ! build the Bethe-Salpeter Equation L0 and Kernel matrices
    subroutine bse_sparse_build(alpha,spindeg,nm_dev,ndiag,nen,En,nop,nnop,blocksize,num_blocks,N,table,&
        G_lesser,G_greater,G_retarded,W,V,&
        Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag)        
        integer,intent(in)::nm_dev,nen,nnop,nop(nnop),ndiag,N,table(2,N)
        integer,intent(in)::blocksize, num_blocks ! arrow block size and number of blocks (excluding tip block)
        real(dp),intent(in)::en(nen),spindeg,alpha                
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen):: G_lesser,G_greater,G_retarded ! electron GFs    
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction    
        complex(dp),intent(out),dimension(blocksize,blocksize*num_blocks,nnop):: Ldiag
        complex(dp),intent(out),dimension(blocksize,blocksize*(num_blocks-1),nnop):: Lupper,Llower ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(out),dimension(nm_dev,blocksize*num_blocks,nnop):: Llowerarrow
        complex(dp),intent(out),dimension(blocksize*num_blocks,nm_dev,nnop):: Lupperarrow ! dense blocks of 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        complex(dp),intent(out),dimension(nm_dev,nm_dev,nnop):: Ltip ! dense tip block of 2-point polarization function with interacting electron-hole at frequency [[nop]]                        
        complex(dp),intent(out),dimension(blocksize*num_blocks):: Kdiag ! diagonal of Kernel
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: Ktip ! dense tip block of Kernel
        !---------
        complex(dp) :: L0ijkl(nnop)              
        real(dp) :: start, finish
        integer :: i,j,k,l,p,q,ie,row,col,it,iop,ib,NT,nepoch,fliped_row,fliped_col        
        !
        print *,'  init memory ...'
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
        !$omp parallel default(shared) private(row,col,i,j,k,l,L0ijkl,ib,p,q)
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
                    if (row<=nm_dev) then 
                        if (col<=nm_dev) then 
                            ! tip block
                            if (nnop>10) then
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                            else    
                                do concurrent (iop=1:nnop)
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo
                            endif         
                            Ltip(nm_dev-row+1,nm_dev-col+1,1:nnop) = L0ijkl * spindeg
                        else
                            ! upper arrow block 
                            if (nnop>10) then
                                call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                            else    
                                do concurrent (iop=1:nnop)
                                    call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                enddo
                            endif 
                            Lupperarrow(NT-(col-nm_dev)+1,nm_dev-row+1,1:nnop) = L0ijkl * spindeg          
                        endif 
                    else 
                        if (col<=nm_dev) then 
                            ! lower arrow block 
                            call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                            G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                            Llowerarrow(nm_dev-col+1,NT-(row-nm_dev)+1,1:nnop) = L0ijkl * spindeg   
                        else 
                            ib = (row-nm_dev-1) / blocksize
                            p = ib * blocksize + nm_dev
                            q = p + blocksize
                            if ((col>p).and.(col<=q)) then 
                                ! diag block 
                                if (nnop>10) then
                                    call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                                else    
                                    do concurrent (iop=1:nnop)
                                        call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                    enddo
                                endif 
                                Ldiag(blocksize - mod(row-nm_dev-1,blocksize), NT-(col-nm_dev)+1, 1:nnop) = L0ijkl * spindeg   
                            else
                                if ((col>q).and.(col<=(q+blocksize))) then                             
                                    ! upper diag block 
                                    if (nnop>10) then
                                        call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                                    else    
                                        do concurrent (iop=1:nnop)
                                            call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                        enddo
                                    endif 
                                    Lupper(blocksize - (row-p-1), NT-(col-nm_dev)+1,1:nnop) = L0ijkl * spindeg                               
                                endif
                                if ((col>(p-blocksize)).and.(col<=p)) then                             
                                    ! lower diag block
                                    if (nnop>10) then
                                        call four_polarization_fft(alpha,nm_dev,nen,en,nop,nnop,ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                                    else    
                                        do concurrent (iop=1:nnop)
                                            call four_polarization(alpha,nm_dev,nen,en,nop(iop),ndiag,&
                                                    G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl(iop))
                                        enddo
                                    endif 
                                    Llower(blocksize - (row-p-1), NT-(col-nm_dev)+1-blocksize,1:nnop ) = L0ijkl * spindeg                            
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
        !$omp parallel default(shared) private(row,col,i,j,k,l)
        !$omp do        
        do row=1,N                        
            do col=1,N
                i=table(1,row)
                j=table(2,row)
                k=table(1,col)
                l=table(2,col)           
                if ((i==j).and.(k==l)) then     
                    fliped_row = nm_dev-row+1
                    fliped_col = nm_dev-col+1                   
                    Ktip(fliped_row,fliped_col) = Ktip(fliped_row,fliped_col) - c1i *  V(i,k) * spindeg        
                    if ((i==k).and.(j==l).and.(row<=nm_dev)) then                
                        Ktip(fliped_row,fliped_row) = Ktip(fliped_row,fliped_row) + c1i *  W(i,j)
                    endif
                endif 
                if ((i==k).and.(j==l).and.(row>nm_dev)) then                        
                    fliped_row = N-row+1                    
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
  
    pure subroutine four_polarization(alpha,nm_dev,nen,en,nop,ndiag,G_lesser,G_greater,G_retarded,i,j,k,l,L0)
       integer,intent(in) :: nm_dev,nen,nop,ndiag, i,j,k,l
       real(dp),intent(in) :: en(nen), alpha 
       complex(dp),intent(in),dimension(nm_dev,nm_dev,nen) :: G_lesser,G_greater,G_retarded
       complex(dp),intent(out) :: L0
       ! ---
       real(dp) :: dE, weights, xen
       integer :: ie, isub, ik, ikd
       ! the P4 IPA tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
       dE = ( En(2) - En(1) )          
       weights = dE/twopi
       !                        
       ! calculate P4_IPA from GG
       L0 =    (1.0_dp - alpha) * ( sum( G_lesser(j,l,(nop+1):nen)   * conjg(G_retarded(i,k,1:(nen-nop))) ) &
                                 +  sum( G_retarded(j,l,(nop+1):nen) * G_lesser(k,i,1:(nen-nop)) ) )  &
               + alpha * 0.5_dp * ( sum( G_greater(j,l,(nop+1):nen) * G_lesser(k,i,1:(nen-nop)) )  & 
                                 -  sum( G_lesser(j,l,(nop+1):nen)  * G_greater(k,i,1:(nen-nop)) ) )  
       L0 = L0 * weights 
    end subroutine four_polarization

    ! solve the full Bethe-Salpeter Equation
    subroutine bse_fullsolve(alpha,spindeg,nm_dev,ndiag,nen,En,nop,G_lesser,G_greater,G_retarded,W,V,P_retarded,nn)
        use gw_dense, only: invert_inplace
        integer,intent(in)::nm_dev,nen,nop,ndiag
        real(dp),intent(in)::en(nen),spindeg,alpha        
        integer, intent(out)::nn ! size of the reduced 2-body system
        complex(dp),intent(in),dimension(nm_dev,nm_dev,nen):: G_lesser,G_greater,G_retarded ! electron GFs
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: W ! W_0 static screened Coulomb interaction
        complex(dp),intent(in),dimension(nm_dev,nm_dev) :: V ! bare Coulomb interaction
        complex(dp),intent(out),dimension(nm_dev,nm_dev):: P_retarded ! 2-point polarization function with interacting electron-hole at frequency [[nop]]                
        !complex(dp),intent(out),dimension(nm_dev,nm_dev,nm_dev,nm_dev),optional:: P4_retarded ! 4-point polarization function with interacting electron-hole 
        !---------
        complex(dp),dimension(:,:),allocatable :: Lmat ! two-particle Green's function 
        complex(dp),dimension(:,:),allocatable :: Mmat ! 4-point Kernel
        complex(dp),dimension(:,:),allocatable :: Amat ! system matrix        
        complex(dp) :: epsM, L0ijkl        
        real(dp) :: start, finish
        integer :: N,i,j,k,l,p,q,ie,row,col, it, ii,jj
        integer,allocatable::table(:,:)        
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
        N = nn ! resize the problem
        print *, 'nm_dev=',nm_dev
        print *, 'resized system size=',N 
        ! start computation
        allocate(Mmat(N,N), source=czero)        
        allocate(Lmat(N,N), source=czero)     
        allocate(Amat(N,N), source=czero)  
        !
        start = omp_get_wtime()              
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
                    call four_polarization(alpha,nm_dev,nen,en,nop,ndiag,&
                        G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                    Lmat(row,col) = L0ijkl * spindeg
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
                if ((i==k).and.(j==l)) then                        
                    Mmat(row,col) = Mmat(row,col) + c1i *  W(i,j)
                endif 
            enddo
        enddo    
        !$omp end do
        !$omp end parallel 
        !            
        finish = omp_get_wtime()
        print '("  computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
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
        finish = omp_get_wtime()
        print '("  computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
        print *,'  start invert (I - L0 K)'
        !
        call invert_inplace(Amat(1:N,1:N),N)
        !
        finish = omp_get_wtime()
        print '("  computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
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
        finish = omp_get_wtime()
        print '("  computation time = ", F0.3 ," seconds.")', finish-start
        start = finish
        !                
        deallocate(Mmat,Lmat,Amat)
    end subroutine bse_fullsolve
  


    ! solve the full Bethe-Salpeter Equation
    subroutine bse_fullsolve_orig(alpha,spindeg,nm_dev,ndiag,nen,nsub,En,nop,nk,G_lesser,G_greater,G_retarded,W,V,P_retarded,system,epsilon_M)
        use gw_dense, only: invert_inplace
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
                            call four_polarization(alpha,nm_dev,nen,en,nop,ndiag,&
                                G_lesser,G_greater,G_retarded,i,j,k,l,L0ijkl)
                            Lmat(row,col) = L0ijkl * spindeg               
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
        use gw_dense, only: invert_inplace
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
                if ( abs(Mat(i,j)) > 0.0d0 ) then
                    write(11,'(2I10,2E18.6)') i,j,dble(Mat(i,j)),aimag(Mat(i,j))
                endif
            enddo
            write(11,*)
        enddo        
        close(11)
    end subroutine save_matrix
    
end module bse_dense



