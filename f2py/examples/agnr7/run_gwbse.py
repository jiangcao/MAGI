#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:38:35 2024

@author: jiacao
"""


import numpy as np
from negf import gw_dense,bse_dense
from wannier import wannierham
import matplotlib.pyplot as plt
import os
os.environ["OMP_NUM_THREADS"] = "28"


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
   nen = 3200
   nsub = 1
   nky=1
   nkz=1
   nk=nky*nkz
   niter=0
   eps_screen=1.0
   r0=3.0
   emin=-10.0
   emax= 4.0
   temp =  np.ones(2)* 300.0
   mu = np.array( [-2.25,-2.25 ] )

   ndiag=nb*ns

   if (ndiag==0):
       ldiag=True
   else:
       ldiag=False

   v = np.zeros((nb*length,nb*length,nk), dtype=np.complex128)  
   ham = np.zeros((nb*length,nb*length,nk), dtype=np.complex128)  

   h00,h10 = wannierham.block_mat_def(kx=0.0, ky=0.0, kz=0.0,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
   ham[:,:,0] = wannierham.full_device_mat_def(ky=0.0,kz=0.0,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)
   v[:,:,0] = wannierham.full_device_bare_coulomb(ky=0.0,kz=0.0,length=length,eps=eps_screen,r0=r0,ldiag=ldiag,nb=nb,ns=ns,method='pointlike',n_range=n_range,wannier_center=wannier_center,cell=cell)

   dim_lead = np.ones(2)* nb*ns
   siglead = np.zeros((nb*ns,nb*ns,nen,2,nk), dtype=np.complex128)

   lead_h10 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
   lead_h10[:,:,0,0] = np.transpose( np.conjugate( h10 ) ) 
   lead_h10[:,:,1,0] = h10

   lead_h00 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
   lead_h00[:,:,0,0] = h00
   lead_h00[:,:,1,0] = h00

   lead_coupling = np.zeros((nb*ns,nb*length,2,nk), dtype=np.complex128)
   lead_coupling[0:nb*ns,0:nb*ns,0,0] = lead_h10[:,:,0,0]
   lead_coupling[0:nb*ns,nb*(length-ns):nb*length,1,0] = lead_h10[:,:,1,0]

   egap=1.67
   encut=[7.0,7.0]
   eps_M=np.zeros(nen//4,dtype='complex')

   energies = np.linspace(emin,emax,nen)
   print(nen,nsub)
   G_retarded,G_lesser,G_greater,Sig_r,Sig_l,Sig_g,tr,te,W0_r,W0_l,W0_g = gw_dense.solve_gw_1d_memsaving(
                 niter=niter,nm_dev=nb*length,lx=Lx,length=length,spindeg=2.0,
                 temp=temp,mu=mu,alpha_mix=0.5,
                 nen=nen,en=energies,nb=nb,ns=ns,
                 ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
                 ndiag=ndiag,encut=encut,egap=egap,flatband=False,vertex=False,bse=False,output_files=True)
   print('done')                                                       
   print('current=', -np.sum(tr[:,0]) , np.sum(tr[:,1]))                                                       

   ID_list = [-np.sum(tr[:,0]) , np.sum(tr[:,1]) ]    

   Gr_diag = G_retarded.diagonal()
   Gn_diag = G_lesser.diagonal()
   ldos = np.imag(Gr_diag)
   ndos = np.imag(Gn_diag)

            
   dE=energies[1]-energies[0]
   for iop in range(nen//4):

       print(iop*4, "Ephot=", iop*dE)

       P_r,nn = bse_dense.bse_fullsolve(
               alpha=0.99,spindeg=2.0,nm_dev=nb*length,ndiag=ndiag,nen=nen,en=energies,nop=iop*4,
               g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,
               w=W0_r,v=v)

       epsilon_M = v[:,:,0] @ P_r
       eps_M[iop] = np.sum( epsilon_M[ nb*length//2, nb*ns:(nb*length-nb*ns) ] )
       print('epsilon_2=', np.imag(eps_M[iop]) )


   np.savez('run_ndiag'+str(ndiag)+'.npz', 
                        ID_list=ID_list,
                        tr=tr,
                        te=te,
                        energies=energies,
                        ldos=ldos,
                        ndos=ndos,
                        eps_M=eps_M)
