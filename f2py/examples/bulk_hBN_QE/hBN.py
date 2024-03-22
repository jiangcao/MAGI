# %%

import numpy as np
# from util import linalg
from negf import gf_dense
from wannier import wannierham
import matplotlib.pyplot as plt
import os
os.environ["OMP_NUM_THREADS"] = "128"
# %%
nb=16
nx=5
ny=15
nz=15

hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',lreorder_axis=True,axis=[3,2,1],nb=nb,nx=nx,ny=ny,nz=nz)

ns = 2
length = 8
nen = 500 # number of energy points
nsub = 3 # number of Legendre nodes in each interval

niter=50
eps_screen=4.0
r0=3.0
emin=-15.0
emax=35.0

Lz=L[2]
Ly=L[1]
Lx=L[0]

dim_lead = np.ones(2)* nb*ns
temp =  np.ones(2)* 300.0
mu = np.array( [11.5, 11.2] )
   

# %%
ndiag_list = np.array([0,nb,nb*2],dtype='i')
nk_list = np.array([1,4,6,8],dtype='i')
for nky in nk_list:    
    nkz = nky
    nk=nky*nkz
    dkz=2.0*np.pi/Lz / nkz
    dky=2.0*np.pi/Ly / nky
    for ndiag in ndiag_list:   
        

        if (ndiag==0):
            ldiag=True
        else:
            ldiag=False

        v = np.zeros((nb*length,nb*length,nk), dtype=np.complex128)  
        ham = np.zeros((nb*length,nb*length,nk), dtype=np.complex128)  
        
        lead_h00 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
        lead_h10 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
        lead_coupling = np.zeros((nb*ns,nb*length,2,nk), dtype=np.complex128)

        for iky in range(nky):
            for ikz in range(nkz):
                ik = ikz + iky*nkz
                ky=-np.pi/Ly + dky*iky
                kz=-np.pi/Lz + dkz*ikz
                ham[:,:,ik] = wannierham.full_device_mat_def(ky=ky,kz=kz,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)
                v[:,:,ik] = wannierham.full_device_bare_coulomb(ky=kz,kz=kz,length=length,eps=eps_screen,r0=r0,ldiag=ldiag,nb=nb,ns=ns,
                                                            method='pointlike',n_range=n_range,wannier_center=wannier_center,cell=cell)
            
                h00,h10 = wannierham.block_mat_def(kx=0.0, ky=ky, kz=kz,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
            
                lead_h10[:,:,0,ik] = np.transpose( np.conjugate(h10) ) 
                lead_h10[:,:,1,ik] = h10
                lead_h00[:,:,0,ik] = h00
                lead_h00[:,:,1,ik] = h00
                lead_coupling[0:nb*ns,0:nb*ns,0,ik] = lead_h10[:,:,0,ik]
                lead_coupling[0:nb*ns,nb*(length-ns):nb*length,1,ik] = lead_h10[:,:,1,ik]


        nen_list = np.array([5000,8000,10000],dtype='i')
        nsub_list = np.array([1,2],dtype='i')

        ID_list = np.zeros((2,nsub_list.shape[0],nen_list.shape[0]))

        for i, toten in enumerate(nen_list):
            print(toten)
            for j, nsub in enumerate(nsub_list):
                nen = toten
                energies = np.linspace(emin,emax,nen)
                print('ndiag=',ndiag)
                print('Nky=',nky,'Nkz=',nkz)
                print('nen=',nen,'nsub=',nsub)
                G_retarded,G_lesser,G_greater,W0,tr = gf_dense.solve_gw_3d(scba_tol=1e-4,niter=niter,nm_dev=nb*length,lx=Lx,length=length,spindeg=2.0,
                                                            temps=temp[0],tempd=temp[1],mus=mu[0],mud=mu[1],alpha_mix=0.5,
                                                            nen=nen,nsub=nsub,en=energies,nb=nb,ns=ns,nphiy=nky,nphiz=nkz,
                                                            ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
                                                            ndiag=ndiag,flatband=False,output_files=False)
                print('NEGF done')                                                       
                print('current=', -np.sum(tr[:,0]) , np.sum(tr[:,1]))                                                       
                ID_list[:,j,i]= [-np.sum(tr[:,0]) , np.sum(tr[:,1]) ]    
                
                Gr_diag = G_retarded.diagonal()
                Gn_diag = G_lesser.diagonal()
                
                np.savez('run_nk'+str(nk)+'_ndiag'+str(ndiag)+'.npz', 
                                    ID_list=ID_list,
                                    nen_list=nen_list,
                                    nsub_list=nsub_list,
                                    tr=tr,
                                    energies=energies,
                                    Gr_diag=Gr_diag,
                                    Gn_diag=Gn_diag)

# %%



