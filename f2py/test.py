
import numpy as np
from util import linalg
from negf import gf_dense
from wannier import wannierham
import matplotlib.pyplot as plt

if __name__=='__main__':   

   nb=14
   nx=21
   ny=1
   nz=1

   hr,wannier_center,n_range,cell = wannierham.load_from_file(fname='ham_dat',lreorder_axis=False,axis=[1,2,3],nb=nb,nx=nx,ny=ny,nz=nz)

   ns = 2

   h00,h10 = wannierham.block_mat_def(kx=0.0, ky=0.0, kz=0.0,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)

   length = 20 

   # ham = wannierham.full_device_mat_def(ky=0.0,kz=0.0,nb=nb,ns=2,length=length,hr=hr,cell=cell,n_range=n_range)

   nen = 2000
   nky=1
   nkz=1
   nk=nky*nkz

   v = np.zeros((nb*length,nb*length,nk), dtype='complex')  
   ham = np.zeros((nb*length,nb*length,nk), dtype='complex')  
   ham[:,:,0] = wannierham.full_device_mat_def(ky=0.0,kz=0.0,nb=nb,ns=2,length=length,hr=hr,cell=cell,n_range=n_range)

   energies = np.linspace(-4,2,nen)

   dim_lead = np.ones(2)* nb*ns
   temp =  np.ones(2)* 300.0
   mu = np.array( [-1, -1.2] )

   siglead = np.zeros((nb*ns,nb*ns,nen,2,nk), dtype='complex')

   lead_h10 = np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
   lead_h10[:,:,0,0] = np.transpose( np.conjugate(h10) ) 
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

   G_retarded,G_lesser,G_greater,P_retarded,P_lesser,P_greater,W_retarded,W_lesser,W_greater,Sig_retarded,Sig_lesser,Sig_greater, Sig_retarded_new,Sig_lesser_new,Sig_greater_new = gf_dense.solve_gw_3d(niter=10,nm_dev=nb*length,lx=4.26,length=length,spindeg=2.0,temps=300.0,tempd=300.0,mus=mu[0],mud=mu[1],alpha_mix=0.5,nen=nen,en=energies,nb=nb,ns=ns,nphiy=nky,nphiz=nkz,ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,ldiag=True,flatband=False)