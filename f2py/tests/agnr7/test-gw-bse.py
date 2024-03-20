
import numpy as np
# from util import linalg
from negf import gf_dense, fft_mod, bse_dense
from wannier import wannierham
import matplotlib.pyplot as plt

if __name__=='__main__':   

   nb=14
   nx=21
   ny=1
   nz=1

   hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',lreorder_axis=False,axis=[1,2,3],nb=nb,nx=nx,ny=ny,nz=nz)

   Lz=L[2]
   Ly=L[1]
   Lx=L[0]

   ns = 2
   length = 6 
   nen = 500
   nsub = 5
   nky=1
   nkz=1
   nk=nky*nkz
   niter=0
   eps_screen=1.0
   r0=3.0
   emin=-10.0
   emax=4.0

   v = np.zeros((nb*length,nb*length,nk), dtype='complex')  
   ham = np.zeros((nb*length,nb*length,nk), dtype='complex')  

   h00,h10 = wannierham.block_mat_def(kx=0.0, ky=0.0, kz=0.0,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
   ham[:,:,0] = wannierham.full_device_mat_def(ky=0.0,kz=0.0,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)
   v[:,:,0] = wannierham.full_device_bare_coulomb(ky=0.0,kz=0.0,length=length,eps=eps_screen,r0=r0,ldiag=True,nb=nb,ns=ns,method='pointlike',n_range=n_range,wannier_center=wannier_center,cell=cell)

   energies = np.linspace(emin,emax,nen)

   dim_lead = np.ones(2)* nb*ns
   temp =  np.ones(2)* 7.0
   mu = np.array( [-2.85, -2.85] )

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

   # sig_r = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   # sig_l = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   # sig_g = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   
   ndiag=nb

   G_retarded,G_lesser,G_greater,W0,tr = gf_dense.solve_gw_3d(scba_tol=1e-3,niter=niter,nm_dev=nb*length,lx=Lx,length=length,spindeg=2.0,
                                                       temps=temp[0],tempd=temp[1],mus=mu[0],mud=mu[1],alpha_mix=0.5,
                                                       nen=nen,nsub=nsub,en=energies,nb=nb,ns=ns,nphiy=nky,nphiz=nkz,
                                                       ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
                                                       ndiag=ndiag,num_lead=2,flatband=False,output_files=True)

   nop=200         
   P_retarded1,system1,epsilon1,L1,M1,nn = bse_dense.bse_fullsolve(alpha=0.99,spindeg=2.0,ndiag=ndiag,nm_dev=nb*length,nen=nen,nsub=nsub,en=energies,nop=nop,nk=nk,g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,w=W0[:,:,0],v=v[:,:,0])                                                        
   plt.spy(system1)
   plt.savefig('pattern1.png')   
   # plt.show()
   plt.spy(L1)
   plt.savefig('L-pattern1.png')  
   plt.spy(M1)
   plt.savefig('M-pattern1.png')     
   # plt.matshow(np.real(M1))
   # plt.show()
   # plt.matshow(np.imag(M1))
   # plt.show()

   P_retarded2,system2,epsilon2 = bse_dense.bse_fullsolve_orig(alpha=0.99,spindeg=2.0,ndiag=ndiag,nm_dev=nb*length,nen=nen,nsub=nsub,en=energies,nop=nop,nk=nk,g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,w=W0[:,:,0],v=v[:,:,0])                                                        
   plt.spy(system2)
   plt.savefig('pattern2.png')   
   

   print('Max error=', np.max(np.abs(P_retarded1-P_retarded2)))
   print('Max element in 1=', np.max(np.abs(P_retarded1)))
   print('Max element in 2=', np.max(np.abs(P_retarded2)))

   blocksize = nb*length
   ndiag = nb

   # A = np.zeros((blocksize*ndiag*2,blocksize*ndiag*2),dtype='complex')
   # L = np.zeros((blocksize*ndiag*2,blocksize*ndiag*2),dtype='complex')
   # A[-nn:,-nn:] = system1[:nn,:nn]
   # for i in range(blocksize*ndiag*2-nn):
   #    A[i,i]=1.0
   # L[-nn:,-nn:] = L1[:nn,:nn]

   np.savez('system_L.npz', A=system1, L=L1, ndiag=ndiag,blocksize=blocksize)




   # P_retarded3 = bse_dense.bse_solve(alpha=0.99,spindeg=2.0,ndiag=ndiag,nm_dev=nb*length,nen=nen,nsub=nsub,en=energies,nop=nop,nk=nk,g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,w=W0[:,:,0],v=v[:,:,0])                                                        
   # print('Max error=', np.max(np.abs(P_retarded3-P_retarded2)))
   # print('Max element in 3=', np.max(np.abs(P_retarded3)))

#    dE = energies[1] - energies[0]
#    for nop in range(10,int(4.0/dE)):
#        print( nop, dE*nop )
#        P_retarded3 = bse_dense.bse_fullsolve(alpha=0.99,spindeg=2.0,ndiag=ndiag,nm_dev=nb*length,
#                                          nen=nen,nsub=nsub,en=energies,nop=nop,nk=nk,
#                                          g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,
#                                          w=W0[:,:,0],v=v[:,:,0])                                                        
   
  
  
   
