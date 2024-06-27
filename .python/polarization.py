import numpy as np

def four_polarization(alpha,nm_dev,nen,En,nop,ndiag,G_lesser,G_greater,G_retarded,i,j,k,l,L0):
    # the P4 Independent-Particle tensor is computed from $P4(q,E') = \sum_{k} \int dE G(E,k) G(E-E',k-q)                
       dE = En[1] - En[0]           
       weights = dE/2.0/np.pi                               
       # calculate P4_IPA from GG
       L0 = (  (1.0 - alpha) * ( np.dot( G_lesser[j,l,nop:nen]   , np.conj(G_retarded[i,k,0:(nen-nop)]) ) &
                               + np.dot( G_retarded[j,l,nop:nen] , G_lesser[k,i,0:(nen-nop)] ) ) 
               + alpha * 0.5 * ( np.dot( G_greater[j,l,nop:nen]  , G_lesser[k,i,0:(nen-nop)]  )    
                               - np.dot( G_lesser[j,l,nop:nen]  , G_greater[k,i,0:(nen-nop)] ) )   )
       L0 = L0 * weights
       return L0