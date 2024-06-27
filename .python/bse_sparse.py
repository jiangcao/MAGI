import numpy as np 
from polarization import four_polarization

# preprocessing the sparsity pattern and decide the block_size and num_blocks in the BTA matrix
def bse_sparse_pre(nm_dev,ndiag):
    N = nm_dev**2 - (nm_dev-ndiag-1)*(nm_dev-ndiag) # compressed system size ~ 2*nm_dev*ndiag-ndiag*ndiag                        
    table=np.zeros((2,N),dtype=int) 
    inverse_table=np.zeros((nm_dev,nm_dev),dtype=int)    
    # construct a lookup table of reordered indices 
    # tip for the "exchange" space， where we put the i=j
    for i in range(nm_dev):            
        table[:,i] = [i,i]    
        inverse_table[i,i] = i
    # then put the "direct" space, but within the ndiag (interaction range)
    it=nm_dev
    for i in range(nm_dev):
        l = max(0,i-ndiag)
        k = min(nm_dev-1,i+ndiag)
        for j in range( l , k+1 ):
            if (i != j):
                table[:,it] = [i,j]
                inverse_table[i,j] = it
                it += 1                                
        
    if (it!=N): 
        print(f'ERROR!, it={it}, N={N}', flush=True)
                    
    # determine coordinates of nnz
    nnz=0
    bandwidth=0
    for row in range(N):
        for col in range(N):         
            i=table[1,row]
            j=table[2,row]
            k=table[1,col]
            l=table[2,col]    
            if ((abs(i-k)<=ndiag) and (abs(j-l)<=ndiag) and (abs(j-k)<=ndiag) and 
                (abs(i-l)<=ndiag) and (abs(i-j)<=ndiag) and (abs(k-l)<=ndiag)):              
                nnz+=1 
                if ((col>nm_dev) and (row>nm_dev) and (abs(col-row)>bandwidth)):  
                    bandwidth = abs(col-row) 
                
    blocksize = bandwidth 
    num_blocks = np.ceil( (N - nm_dev) / blocksize )  
    NT = blocksize * num_blocks         
    print ("  total arrow size=", NT)
    print ("  arrow block size=", blocksize)
    print ("  arrow number of blocks=", num_blocks)
    print ("  nonzero elements=", nnz/1e6," Million")
    print ("  nonzero ratio = ", nnz/(NT+nm_dev)**2*100," %")
    
    return (N,nnz,table,inverse_table,blocksize,num_blocks)

def bse_sparse_build(alpha,spindeg,nm_dev,ndiag,nen,en,nop,blocksize,num_blocks,N,table,
        G_lesser,G_greater,G_retarded,W,V):

        (N,nnz,table,inverse_table,blocksize,num_blocks) = bse_sparse_pre(nm_dev,ndiag)
           
        
        print('  init memory ...')
        Ltip = np.zeros((nm_dev,nm_dev),dtype="complex")
        Ldiag = np.zeros((blocksize,blocksize*num_blocks),dtype="complex")
        Lupper = np.zeros((blocksize,blocksize*(num_blocks-1)),dtype="complex")
        Llower = np.zeros((blocksize,blocksize*(num_blocks-1)),dtype="complex")
        Lupperarrow = np.zeros((blocksize*num_blocks,nm_dev),dtype="complex")
        Llowerarrow = np.zeros((nm_dev,blocksize*num_blocks),dtype="complex")
        
        NT = blocksize * num_blocks   
        
        print('  start computation L0_ijkl = G_jl G_ki ...')                 
        for row in range(0,N): 
            
            for col in range(0,N):         

                i=table[0,row]
                j=table[1,row]
                k=table[0,col]
                l=table[1,col]    
                if ((abs(i-k)<=ndiag) and (abs(j-l)<=ndiag) and (abs(j-k)<=ndiag) and 
                    (abs(i-l)<=ndiag) and (abs(i-j)<=ndiag) and (abs(k-l)<=ndiag)):                    
                    # need to flip the row and col when putting into arrowhead structure               
                    fliped_row = NT + nm_dev - col - 1
                    fliped_col = NT + nm_dev - row - 1
                    if (fliped_col >= NT) : 
                        if (fliped_row >= NT) : 
                            L0ijkl = four_polarization(alpha,nm_dev,nen,en,nop,ndiag,
                                            G_lesser,G_greater,G_retarded,i,j,k,l)
                                   
                            Ltip[fliped_row - NT,fliped_col - NT] = L0ijkl * spindeg
                        else:
                            # upper arrow block 
                            L0ijkl = four_polarization(alpha,nm_dev,nen,en,nop,ndiag,
                                            G_lesser,G_greater,G_retarded,i,j,k,l)
                            Lupperarrow[fliped_row,fliped_col - NT] = L0ijkl * spindeg          
                         
                    else: 
                        if (fliped_row >= NT) : 
                            # lower arrow block 
                            L0ijkl = four_polarization(alpha,nm_dev,nen,en,nop,ndiag,
                                            G_lesser,G_greater,G_retarded,i,j,k,l)
                            Llowerarrow[fliped_row - NT,fliped_col] = L0ijkl * spindeg   
                        else: 
                            ib = fliped_row // blocksize
                            p = ib * blocksize 
                            q = p + blocksize
                            if ((fliped_col >= p)and(fliped_col < q)):  
                                # diag block 
                                L0ijkl = four_polarization(alpha,nm_dev,nen,en,nop,ndiag,
                                            G_lesser,G_greater,G_retarded,i,j,k,l)
                                Ldiag[fliped_row - p, fliped_col] = L0ijkl * spindeg   
                            else:
                                if ((fliped_col >= q)and(fliped_col < (q+blocksize))) :
                                    # upper diag block 
                                    L0ijkl = four_polarization(alpha,nm_dev,nen,en,nop,ndiag,
                                            G_lesser,G_greater,G_retarded,i,j,k,l)
                                    Lupper[fliped_row - p, fliped_col - blocksize] = L0ijkl * spindeg                               
                                
                                if ((fliped_col >= (p-blocksize))and(fliped_col < p)) :
                                    # lower diag block
                                    L0ijkl = four_polarization(alpha,nm_dev,nen,en,nop,ndiag,
                                            G_lesser,G_greater,G_retarded,i,j,k,l) 
                                    Llower[fliped_row - p, fliped_col] = L0ijkl * spindeg                            
                                
        Ktip=np.zeros((nm_dev,nm_dev),dtype="complex")
        Kdiag=np.zeros((blocksize*num_blocks),dtype="complex")
             
        for row in range(0,N):
            for col in range(0,N):
                i=table[0,row]
                j=table[1,row]
                k=table[0,col]
                l=table[1,col]           
                fliped_row = NT + nm_dev - col - 1
                fliped_col = NT + nm_dev - row - 1           
                if ((i==j) and (k==l)) :                               
                    Ktip[fliped_row - NT, fliped_col - NT] += - 1j *  V[i,k] * spindeg        
                    if ((i==k) and (j==l) and (row<nm_dev)):
                        Ktip[fliped_row - NT,fliped_row - NT] += 1j *  W[i,j]
                                     
                if ((i==k) and (j==l) and (row>=nm_dev)): 
                    Kdiag[fliped_row] += 1j *  W[i,j]
                
        return ( Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag )


def bse_sparse_build_system(blocksize,num_blocks,nm_dev,
        Ldiag, Lupper, Llower, Llowerarrow, Lupperarrow, Ltip, Kdiag, Ktip):
        
        Adiag = np.zeros((blocksize,blocksize*num_blocks),dtype="complex")
        Aupper = np.zeros((blocksize,blocksize*(num_blocks-1)),dtype="complex")
        Alower = np.zeros((blocksize,blocksize*(num_blocks-1)),dtype="complex") 
        Alowerarrow = np.zeros((nm_dev,blocksize*num_blocks),dtype="complex")
        Aupperarrow = np.zeros((blocksize*num_blocks,nm_dev),dtype="complex")
        Atip = np.zeros((nm_dev,nm_dev),dtype="complex") 
    
        N = blocksize*num_blocks
        # A_xx
        Atip = - Ltip @ Ktip + np.eye(nm_dev)
        
        # A_xd = - L_xd * K_dd
        for i in range(0,nm_dev):
            for j in range(0,N): 
                Alowerarrow[i,j] = - Llowerarrow[i,j] * Kdiag[j]
            
        # A_dx = - L_dx * K_xx
        Aupperarrow = - Lupperarrow @ Ktip        
        
        # A_dd
        # diagonal blocks         
        for i in range (0,blocksize):
            for j in range (0, blocksize*num_blocks):
                Adiag[i,j] = - Ldiag[i,j] * Kdiag[j]

        for ib in range(0,num_blocks):
            for i in range(0,blocksize):
                Adiag[i,i+ib*blocksize] += 1.0
            
        # upper and lower diagonal blocks
        for i in range(0,blocksize):
            for j in range(0,blocksize*(num_blocks-1)):
                Aupper[i,j] = - Lupper[i,j] * Kdiag[j+blocksize]
                Alower[i,j] = - Llower[i,j] * Kdiag[j]
            
        return (Adiag, Aupper, Alower, Alowerarrow, Aupperarrow, Atip)