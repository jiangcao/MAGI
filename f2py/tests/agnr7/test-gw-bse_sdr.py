
import numpy as np
# from util import linalg
# from negf import gf_dense, fft_mod, bse_dense
# from wannier import wannierham
import matplotlib.pyplot as plt

from sdr.lu.lu_factorize import lu_factorize_tridiag_arrowhead
from sdr.lu.lu_selected_inversion import lu_sinv_tridiag_arrowhead
from sdr.lu.lu_solve import lu_slv_tridiag_arrowhead
from sdr.utils import matrix_generation
from sdr.utils.matrix_transform import (cut_to_blocktridiag_arrowhead,
                                        from_arrowhead_arrays_to_dense,
                                        from_dense_to_arrowhead_arrays)

if __name__=='__main__':   

   nb=14
   ns = 2
   length = 6 
  

   nm_dev = nb*length
   ndiag = nm_dev

   npzfile = np.load('system_L.npz')
   print(npzfile.files)
   A = npzfile['A'] 
   L = npzfile['L'] 

   print(A.shape[0],A.shape[1])
   # plt.matshow(np.real(A))
   # plt.show()
   # plt.matshow(np.imag(A))
   # plt.show()
   #A = np.real(A)   


   plt.spy(L)
   plt.show()

   import scipy.linalg as la

   diag_blocksize = npzfile['blocksize']
   arrow_blocksize = npzfile['blocksize']
   ndiag = npzfile['ndiag']


   print('diag blocksize=',diag_blocksize)
   print('tip blocksize=',arrow_blocksize)
   print('number of blocks=',(A.shape[0]-arrow_blocksize)//diag_blocksize)
   print('ndiag=',ndiag)

   B_arrow_tip_block = np.array(L[-arrow_blocksize:,-arrow_blocksize:])

   B = np.zeros((A.shape[0],arrow_blocksize),dtype=A.dtype)
   B[-arrow_blocksize:,:] = B_arrow_tip_block

   # plt.spy(B)
   # plt.show()

   X_ref = la.inv(A)

   P_ref = la.solve(A, B)
   P_ref = - 1j* P_ref[-arrow_blocksize:,:]

   # plt.spy(P_ref)
   # plt.show()


   (
        A_diagonal_blocks,
        A_lower_diagonal_blocks,
        A_upper_diagonal_blocks,
        A_arrow_bottom_blocks,
        A_arrow_right_blocks,
        A_arrow_tip_block,
    ) = from_dense_to_arrowhead_arrays(A, diag_blocksize, arrow_blocksize)

   (
        L_diagonal_blocks,
        L_lower_diagonal_blocks,
        L_arrow_bottom_blocks,
        U_diagonal_blocks,
        U_upper_diagonal_blocks,
        U_arrow_right_blocks,
        P_arrow_tip_blocks
    ) = lu_factorize_tridiag_arrowhead(
        A_diagonal_blocks,
        A_lower_diagonal_blocks,
        A_upper_diagonal_blocks,
        A_arrow_bottom_blocks,
        A_arrow_right_blocks,
        A_arrow_tip_block,
    )
 

   (
        X_sdr_diagonal_blocks,
        X_sdr_lower_diagonal_blocks,
        X_sdr_upper_diagonal_blocks,
        X_sdr_arrow_bottom_blocks,
        X_sdr_arrow_right_blocks,
        X_sdr_arrow_tip_block,
    ) = lu_sinv_tridiag_arrowhead(
        L_diagonal_blocks,
        L_lower_diagonal_blocks,
        L_arrow_bottom_blocks,
        U_diagonal_blocks,
        U_upper_diagonal_blocks,
        U_arrow_right_blocks,
    )

   X_sdr_arrow_tip_block = X_sdr_arrow_tip_block @ np.transpose(P_arrow_tip_blocks)

   X_sdr = from_arrowhead_arrays_to_dense(
        X_sdr_diagonal_blocks,
        X_sdr_lower_diagonal_blocks,
        X_sdr_upper_diagonal_blocks,
        X_sdr_arrow_bottom_blocks,
        X_sdr_arrow_right_blocks,
        X_sdr_arrow_tip_block,
    )


   X_diff = np.array(X_ref[-arrow_blocksize:,-arrow_blocksize:]) - X_sdr_arrow_tip_block


   # plt.spy(X_sdr)
   # plt.show()

   print('tip max abs diff=', np.max(np.abs(X_diff)) , '/ tip max abs Xref=', np.max(np.abs(X_ref[-arrow_blocksize:,-arrow_blocksize:])) )

   for i in range((A.shape[0]-arrow_blocksize)//diag_blocksize):
      X_diff = np.array(X_ref[diag_blocksize*i:diag_blocksize*(i+1),diag_blocksize*i:diag_blocksize*(i+1)]) - X_sdr_diagonal_blocks[:,diag_blocksize*i:diag_blocksize*(i+1)]
      print(i,'diag block max abs diff=', np.max(np.abs(X_diff)), '/ max abs Xref=', np.max(np.abs(X_sdr_diagonal_blocks[:,diag_blocksize*i:diag_blocksize*(i+1)])) ) 

   P_sdr = - 1j* X_sdr @ B

   P_diff = P_ref - P_sdr[-arrow_blocksize:,:]

   # plt.matshow(np.imag(P_diff))
   # plt.show()

   print('max abs diff=', np.max(np.abs(P_diff)) )
   print('max abs Pref=', np.max(np.abs(P_ref)) )
   print('max abs Psdr=', np.max(np.abs(P_sdr)) )

   
