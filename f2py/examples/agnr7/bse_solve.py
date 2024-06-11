import numpy as np
import time


def four_polarization_loop(
    alpha: float,
    nm_dev: int,
    nen: int,
    en: np.ndarray,
    nop: int,
    ndiag: int,
    G_lesser: np.ndarray,
    G_greater: np.ndarray,
    G_retarded: np.ndarray,
    i: int,
    j: int,
    k: int,
    l: int,
) -> complex:
    """
    Computes the P4 IPA tensor using Green's functions G_lesser, G_greater, and G_retarded.

    Args:
        alpha (float): Mixing parameter.
        nm_dev (int): Dimension of the matrix.
        nen (int): Number of energy points.
        en (np.ndarray): Energy values.
        nop (int): Offset for energy index.
        ndiag (int): Not used in this function.
        G_lesser (np.ndarray): Lesser Green's function.
        G_greater (np.ndarray): Greater Green's function.
        G_retarded (np.ndarray): Retarded Green's function.
        i, j, k, l (int): Indices.

    Returns:
        complex: P4 IPA tensor value.
    """
    dE = en[1] - en[0]
    weights = dE / (2 * np.pi)
    L0 = 0.0j

    for ie in range(nop, nen):
        L0 += (1.0 - alpha) * (
            G_lesser[j, l, ie] * np.conj(G_retarded[i, k, ie - nop])
            + G_retarded[j, l, ie] * G_lesser[k, i, ie - nop]
        ) + alpha * 0.5 * (
            G_greater[j, l, ie] * G_lesser[k, i, ie - nop]
            - G_lesser[j, l, ie] * G_greater[k, i, ie - nop]
        )

    L0 *= weights
    return L0

def four_polarization_dot(alpha,nm_dev,nen,En,nop,ndiag,G_lesser,G_greater,G_retarded,i,j,k,l):
    # nop: number of discretizqtino of the electrons grid
    # nen: number of electrons energies
    # the P4 Independent-Particle tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
       dE = En[1] - En[0]           
       weights = dE/2.0/np.pi                               
       # calculate P4_IPA from GG
       L0 = (  (1.0 - alpha) * ( np.dot( G_lesser[j,l,nop:nen]   , np.conj(G_retarded[i,k,0:(nen-nop)]) )
                               + np.dot( G_retarded[j,l,nop:nen] , G_lesser[k,i,0:(nen-nop)] ) ) 
               + alpha * 0.5 * ( np.dot( G_greater[j,l,nop:nen]  , G_lesser[k,i,0:(nen-nop)]  )    
                               - np.dot( G_lesser[j,l,nop:nen]  , G_greater[k,i,0:(nen-nop)] ) )   )
       L0 = L0 * weights
       return L0


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
    # then put the others, but within the ndiag
    it=nm_dev
    for i in range(nm_dev):
        l = max(0,i-ndiag)
        k = min(nm_dev-1,i+ndiag)
        for j in range( l , k+1 ):
            if (i != j):
                table[:,it] = [i,j]
                inverse_table[i,j] = it
                it += 1                                
        
    if ((it)!=N): 
        print(f'ERROR!, it={it}, N={N}', flush=True)
                    
    # determine coordinates of nnz
    nnz=0
    bandwidth=0
    for row in range(N):
        for col in range(N):         
            i=table[0,row]
            j=table[1,row]
            k=table[0,col]
            l=table[1,col]    
            if ((abs(i-k)<=ndiag) and (abs(j-l)<=ndiag) and (abs(j-k)<=ndiag) and 
                (abs(i-l)<=ndiag) and (abs(i-j)<=ndiag) and (abs(k-l)<=ndiag)):              
                nnz+=1 
                if ((col>nm_dev) and (row>nm_dev) and (abs(col-row)>bandwidth)):  
                    bandwidth = abs(col-row) 
                
    blocksize = bandwidth//2 
    num_blocks = np.ceil( (N - nm_dev) / blocksize )  
    NT = blocksize * num_blocks         
    print ("  total arrow size=", NT)
    print ("  arrow block size=", blocksize)
    print ("  arrow number of blocks=", num_blocks)
    print ("  nonzero elements=", nnz/1e6," Million")
    print ("  nonzero ratio = ", nnz/(NT+nm_dev)**2*100," %")
    
    return (N,nnz,table,inverse_table,blocksize,num_blocks)


def bse_fullsolve(
    alpha: float,
    spindeg: float,
    nm_dev: int,
    ndiag: int,
    nen: int,
    en: np.ndarray,
    nop: int,
    G_lesser: np.ndarray,
    G_greater: np.ndarray,
    G_retarded: np.ndarray,
    W: np.ndarray,
    V: np.ndarray,
) -> np.ndarray:
    """
    Solves the full Bethe-Salpeter Equation.

    Args:
        alpha (float): Mixing parameter.
        spindeg (float): Spin degeneracy.
        nm_dev (int): Dimension of the matrix.
        ndiag (int): Not used in this function.
        nen (int): Number of energy points.
        en (np.ndarray): Energy values.
        nop (int): Offset for energy index.
        G_lesser (np.ndarray): Lesser Green's function.
        G_greater (np.ndarray): Greater Green's function.
        G_retarded (np.ndarray): Retarded Green's function.
        W (np.ndarray): Static screened Coulomb interaction.
        V (np.ndarray): Bare Coulomb interaction.

    Returns:
        np.ndarray: 2-point polarization function with interacting electron-hole at frequency [[nop]].
    """

    (N,nnz,table,inverse_table,blocksize,num_blocks) = bse_sparse_pre(nm_dev,ndiag)

    # Initialize matrices
    Lmat = np.zeros((N, N), dtype=complex)
    Mmat = np.zeros((N, N), dtype=complex)
    Amat = np.zeros((N, N), dtype=complex)

    # Put G on GPU

    # Create mem for Lmat based on nnz

    # Compute L0_ijkl = G_jl * G_ki
    tic = time.perf_counter()
    for row in range(N):
        for col in range(N):
            i, j = table[:, row]
            k, l = table[:, col]
            if (abs(i - k) <= ndiag) and (abs(j - l) <= ndiag):
                # print(f"Start polarization at: i={i}, j={j}, k={k}, l={l}, N={N}", flush=True)
                L0ijkl = four_polarization_dot(
                    alpha,
                    nm_dev,
                    nen,
                    en,
                    nop,
                    ndiag,
                    G_lesser,
                    G_greater,
                    G_retarded,
                    i,
                    j,
                    k,
                    l,
                )
                Lmat[row, col] = L0ijkl * spindeg
                # Store as COO index

    # Get COO to LIL

    # Slice LIL into blocks

    print(f"Finished polarization", flush=True)

    # Construct Mmat
    # -> bse_sparse_build()
    # -- Construct BTA array representation
    for row in range(N):
        for col in range(N):
            i, j = table[:, row]
            k, l = table[:, col]
            if i == j and k == l:
                Mmat[row, col] -= 1j * V[i, k] * spindeg
            if i == k and j == l:
                Mmat[row, col] += 1j * W[i, j]
    toc = time.perf_counter()
    print(f"L0_ijkl = G_jl * G_ki: {toc - tic} seconds", flush=True)

    # -> bse_sparse_build_system()
    # Compute (I - L0 K) -> A
    # -- Modified to block-operations
    tic = time.perf_counter()
    np.fill_diagonal(Amat, 1.0 + 0.0j)
    Amat -= Lmat @ Mmat
    toc = time.perf_counter()
    print(f"-L0 K: {toc - tic} seconds", flush=True)

    # Invert (I - L0 K)
    # Call SINV
    tic = time.perf_counter()
    Amat_inv = np.linalg.inv(Amat)
    toc = time.perf_counter()
    print(f"(I - L0 K)^-1: {toc - tic} seconds", flush=True)

    # Compute L = (I - L0 K) \ L0
    # Solve done with block-matrix mult
    tic = time.perf_counter()
    Lmat = Amat_inv @ Lmat
    toc = time.perf_counter()
    print(f"L = (I - L0 K) \\ L0: {toc - tic} seconds", flush=True)

    # Compute P_retarded
    # (Extract the tip of L)
    P_retarded = np.zeros((nm_dev, nm_dev), dtype=complex)
    for row in range(nm_dev):
        for col in range(nm_dev):
            i, k = table[:, row]
            P_retarded[i, k] = -1j * Mmat[row, col]

    return P_retarded, N
