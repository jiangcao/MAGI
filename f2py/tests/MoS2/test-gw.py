
import numpy as np
from util import linalg
from negf import gf_dense, fft_mod
from wannier import wannierham
import matplotlib.pyplot as plt

if __name__=='__main__':   

   nb=11
   nx=29
   ny=29
   nz=1

   hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',lreorder_axis=False,axis=[1,2,3],nb=nb,nx=nx,ny=ny,nz=nz)

   ns = 2
   length = 10 
   nen = 1600 # number of energy points
   nsub = 4 # number of Legendre nodes in each interval
   nky=9
   nkz=1
   nk=nky*nkz
   niter=10
   eps_screen=2.5
   r0=3.0
   emin=-10.0
   emax=4.0

   Lz=L[2]
   Ly=L[1]
   dkz=2.0*np.pi/Lz / nkz
   dky=2.0*np.pi/Ly / nky

   energies = np.linspace(emin,emax,nen)

   dim_lead = np.ones(2)* nb*ns
   temp =  np.ones(2)* 300.0
   mu = np.array( [-1, -1.2] )

   v = np.zeros((nb*length,nb*length,nk), dtype='complex')  
   ham = np.zeros((nb*length,nb*length,nk), dtype='complex')  
   # siglead = np.zeros((nb*ns,nb*ns,nen,2,nk), dtype='complex')
   lead_h00 = np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
   lead_h10 = np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
   lead_coupling = np.zeros((nb*ns,nb*length,2,nk), dtype='complex')

   for iky in range(nky):
      for ikz in range(nkz):
         ik = ikz + iky*nkz
         ky=-np.pi/Ly + dky*iky
         kz=-np.pi/Lz + dkz*ikz
         ham[:,:,ik] = wannierham.full_device_mat_def(ky=ky,kz=kz,nb=nb,ns=2,length=length,hr=hr,cell=cell,n_range=n_range)
         v[:,:,ik] = wannierham.full_device_bare_coulomb(ky=kz,kz=kz,length=length,eps=eps_screen,r0=r0,ldiag=True,nb=nb,ns=ns,
                                                      method='pointlike',n_range=n_range,wannier_center=wannier_center,cell=cell)
   
         h00,h10 = wannierham.block_mat_def(kx=0.0, ky=ky, kz=kz,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
   
         lead_h10[:,:,0,ik] = np.transpose( np.conjugate(h10) ) 
         lead_h10[:,:,1,ik] = h10
         lead_h00[:,:,0,ik] = h00
         lead_h00[:,:,1,ik] = h00
         lead_coupling[0:nb*ns,0:nb*ns,0,ik] = lead_h10[:,:,0,ik]
         lead_coupling[0:nb*ns,nb*(length-ns):nb*length,1,ik] = lead_h10[:,:,1,ik]

   # sig_r = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   # sig_l = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   # sig_g = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
   # G_retarded,G_lesser,G_greater,cur,te = gf_dense.calc_gf(ne=nen,e=energies,num_lead=2,nm_dev=nb*length,nm_lead=dim_lead,max_nm_lead=nb*ns,
   #                                                          ham=ham,lead_h00=lead_h00,lead_h10=lead_h10,
   #                                                          siglead=siglead,t=lead_coupling,
   #                                                          scat_sig_retarded=sig_r,scat_sig_lesser=sig_l,scat_sig_greater=sig_g,
   #                                                          mu=mu,temp=temp,flatband=False)
   # plt.plot(energies, cur[:,0] )
   # plt.show()

   G_retarded,G_lesser,G_greater,W0 = gf_dense.solve_gw_3d(niter=niter,nm_dev=nb*length,lx=4.26,length=length,spindeg=2.0,
                                                        temps=300.0,tempd=300.0,mus=mu[0],mud=mu[1],alpha_mix=0.5,
                                                        nen=nen,nsub=nsub,en=energies,nb=nb,ns=ns,nphiy=nky,nphiz=nkz,
                                                        ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
                                                        ldiag=False,flatband=False)