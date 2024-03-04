
import numpy as np
from util import linalg
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
   length = 10 
   nen = 1200
   nsub = 3
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
   mu = np.array( [-2.0, -2.0] )

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

   G_retarded,G_lesser,G_greater,W0 = gf_dense.solve_gw_3d(niter=niter,nm_dev=nb*length,lx=Lx,length=length,spindeg=2.0,
                                                        temps=temp[0],tempd=temp[1],mus=mu[0],mud=mu[1],alpha_mix=0.5,
                                                        nen=nen,nsub=nsub,en=energies,nb=nb,ns=ns,nphiy=nky,nphiz=nkz,
                                                        ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
                                                        ldiag=True,flatband=False)

   for nop in range(10,int(5.0/(energies[1]-energies[0])),4):
       print( nop,nop*(energies[1]-energies[0]) )
    #    P4 = bse_dense.four_polarization_dense(nm_dev=nb*length,nen=nen,nsub=nsub,nk=nk,en=energies,nop=nop,ndiag=0,
    #                                           g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded)
       
    #    plt.spy(np.abs(P4))
    #    plt.show()

       P_retarded, system = bse_dense.bse_fullsolve(spindeg=2.0,nm_dev=nb*length,ndiag=4,nen=nen,nsub=nsub,en=energies,nop=nop,nk=nk,
                                            g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,
                                            w_retarded=W0[:,:,0],v=v[:,:,0])                                                        
       plt.spy(system)
       plt.savefig('pattern.png')
#    P_retarded = bse_dense.bse_solve(spindeg=2.0,nm_dev=nb*length,nen=nen,en=energies,nop=nop,g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,w_retarded=v[:,:,0],v=v[:,:,0])                                                        
   