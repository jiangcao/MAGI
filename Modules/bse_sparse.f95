
module bse_sparse    
    use parameters_mod,only:dp,twopi,pi,e_charge,epsilon0,m0_charge,hbar,c1i,czero,cone    
    use omp_lib
    use polarization
    implicit none
    contains

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
            print *, 'ERROR!'
            call abort
        endif        
        print *, 'nm_dev=', nm_dev
        print *, 'ndiag =', ndiag
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
        N = blocksize*num_blocks
        ! A_xx
        call zgemm('n','n',nm_dev,nm_dev,nm_dev,-cone,Ltip(:,:,iop),nm_dev,Ktip,nm_dev,czero,Atip,nm_dev)  
        do concurrent (i=1:nm_dev)
            Atip(i,i) = Atip(i,i) + cone 
        enddo 
        ! A_xd = - L_xd * K_dd
        do concurrent (i = 1:nm_dev)
            do concurrent (j = 1:N) 
                Alowerarrow(i,j) = - Llowerarrow(i,j,iop) * Kdiag(j)
            enddo
        enddo        
        ! A_dx = - L_dx * K_xx
        call zgemm('n','n',N,nm_dev,nm_dev,-cone,Lupperarrow(:,:,iop),N,Ktip,nm_dev,czero,Aupperarrow,N)
        ! A_dd
        ! diagonal blocks         
        do concurrent (i=1:blocksize)
            do concurrent (j=1:blocksize*num_blocks)
                Adiag(i,j) = - Ldiag(i,j,iop) * Kdiag(j)
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
                Aupper(i,j) = - Lupper(i,j,iop) * Kdiag(j+blocksize)
                Alower(i,j) = - Llower(i,j,iop) * Kdiag(j)
            enddo 
        enddo 
    end subroutine bse_sparse_build_system

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
                                do iop=1,nnop
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
                                do iop=1,nnop
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
                                    do iop=1,nnop
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
                                        do iop=1,nnop
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
                                        do iop=1,nnop
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
        !$omp parallel default(shared) private(row,col,i,j,k,l,fliped_row,fliped_col)
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
                    fliped_row = NT-(row-nm_dev)+1                    
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
