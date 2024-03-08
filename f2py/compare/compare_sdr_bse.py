
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
#import pytest

from util import linalg
from negf import gf_dense, fft_mod, bse_dense
from wannier import wannierham
import matplotlib.pyplot as plt

#from lu_decompose import lu_dcmp_ndiags_arrowhead
#from sdr.lu.lu_selected_inversion import lu_sinv_ndiags_arrowhead

if __name__=='__main__':   

   nb=14
   nx=21
   ny=1
   nz=1

   hr,wannier_center,n_range,cell = wannierham.load_from_file(fname='ham_dat',lreorder_axis=False,axis=[1,2,3],nb=nb,nx=nx,ny=ny,nz=nz)

   ns = 2
   length = 10 
   nen = 6400
   nky=1
   nkz=1
   nk=nky*nkz
   niter=1
   eps_screen=2.5
   r0=3.0
   emin=-10.0
   emax=4.0

   v = np.zeros((nb*length,nb*length,nk), dtype='complex')  
   ham = np.zeros((nb*length,nb*length,nk), dtype='complex')  

   h00,h10 = wannierham.block_mat_def(kx=0.0, ky=0.0, kz=0.0,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
   ham[:,:,0] = wannierham.full_device_mat_def(ky=0.0,kz=0.0,nb=nb,ns=2,length=length,hr=hr,cell=cell,n_range=n_range)
   v[:,:,0] = wannierham.full_device_bare_coulomb(ky=0.0,kz=0.0,length=length,eps=eps_screen,r0=r0,ldiag=True,nb=nb,ns=ns,method='pointlike',n_range=n_range,wannier_center=wannier_center,cell=cell)

   energies = np.linspace(emin,emax,nen)

   dim_lead = np.ones(2)* nb*ns
   temp =  np.ones(2)* 300.0
   mu = np.array( [-1, -1.2] )

   siglead = np.zeros((nb*ns,nb*ns,nen,2,nk), dtype='complex')

   lead_h10 = np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
   lead_h10[:,:,0,0] = np.transpose( np.conjugate( h10 ) ) 
   lead_h10[:,:,1,0] = h10

   lead_h00 = np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
   lead_h00[:,:,0,0] = h00
   lead_h00[:,:,1,0] = h00

   lead_coupling = np.zeros((nb*ns,nb*length,2,nk), dtype='complex')
   lead_coupling[0:nb*ns,0:nb*ns,0,0] = lead_h10[:,:,0,0]
   lead_coupling[0:nb*ns,nb*(length-ns):nb*length,1,0] = lead_h10[:,:,1,0]

   sig_r = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   sig_l = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   sig_g = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')

   # G_retarded,G_lesser,G_greater,cur,te = gf_dense.calc_gf(ne=nen,e=energies,num_lead=2,nm_dev=nb*length,nm_lead=dim_lead,max_nm_lead=nb*ns,
   #                                                          ham=ham,lead_h00=lead_h00,lead_h10=lead_h10,
   #                                                          siglead=siglead,t=lead_coupling,
   #                                                          scat_sig_retarded=sig_r,scat_sig_lesser=sig_l,scat_sig_greater=sig_g,
   #                                                          mu=mu,temp=temp,flatband=False)


   # plt.plot(energies, cur[:,0] )
   # plt.show()

   G_retarded,G_lesser,G_greater,W0 = gf_dense.solve_gw_3d(niter=niter,nm_dev=nb*length,lx=4.26,length=length,spindeg=2.0,
                                                        temps=300.0,tempd=300.0,mus=mu[0],mud=mu[1],alpha_mix=0.5,
                                                        nen=nen,en=energies,nb=nb,ns=ns,nphiy=nky,nphiz=nkz,
                                                        ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
                                                        ldiag=True,flatband=False)

  #for nop in [-3.05, -1.3]:
        #print( nop )
        #P_retarded = bse_dense.bse_fullsolve(spindeg=2.0,nm_dev=nb*length,ndiag=2,nen=nen,en=energies,nop=nop,
                                          #  g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,
                                          #  w_retarded=W0[:,:,0],v=v[:,:,0])  
        #P_retarded,system = bse_dense.bse_solve(spindeg=2.0,nm_dev=nb*length,nen=nen,en=energies,nop=nop,g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,w_retarded=v[:,:,0],v=v[:,:,0]) 
        
        #M1 = bse_dense.bse_fullsolve.system #M1=system matrix=(I - L0 K) 
        #M2 = bse_dense.bse_fullsolve.Amat   #M2=after inversion
        #M3 = P_retarded                     #M3=after multiplication



   #plt.spy(system)
   #plt.savefig('pattern.png')

   #norm_M1 = np.linalg.norm(M1, 2) #2-norm of M1
   #norm_M2 = np.linalg.norm(M2, 2) #2-norm of M2 
   #norm_M3 = np.linalg.norm(M3, 2) #2-norm of M3
   
   #comp_time_m1 = bse_fullsolve.comp_time_m1
   #comp_time_m2 = bse_fullsolve.comp_time_m2d
   #comp_time_m3 = bse_fullsolve.comp_time_m3
   
   #start = omp_get_wtime()
   #finish = omp_get_wtime()
   #comp_time1 = finish - start
   #print *, "M1 Calculation duration in seconds", comp_time1


#compare with SDR - LU decomp
#L_sdr, U_sdr = lu_dcmp_ndiags_arrowhead(A, nblocks, diag_blocksize, arrow_blocksize)
#X_sdr = lu_sinv_ndiags_arrowhead(L_sdr, U_sdr, ndiags, diag_blocksize, arrow_blocksize)
