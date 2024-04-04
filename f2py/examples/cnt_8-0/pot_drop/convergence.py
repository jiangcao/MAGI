#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:38:35 2024

@author: jiacao
"""


import numpy as np
# from util import linalg
from negf import gf_dense, fft_mod, bse_dense
from wannier import wannierham
import matplotlib.pyplot as plt
import os
os.environ["OMP_NUM_THREADS"] = "128"


if __name__=='__main__':   

   nb=32
   nx=21
   ny=1
   nz=1

   hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',lreorder_axis=False,axis=[1,2,3],nb=nb,nx=nx,ny=ny,nz=nz)

   Lz=L[2]
   Ly=L[1]
   Lx=L[0]

   ns = 2
   length = 8 
   nen = 2000
   nsub = 3
   nky=1
   nkz=1
   nk=nky*nkz
   niter=50
   eps_screen=4.0
   r0=3.0
   emin=-15.0
   emax= 4.0
   temp =  np.ones(2)* 300.0

   pot_drop = 0.4
   pot = np.zeros(length) 
   pot[:ns*2]  = 0.0
   pot[-ns*2:] = -pot_drop
   pot[ns*2:-ns*2] = np.linspace(0,-pot_drop,length-ns*4) 

   mu = np.array( [-3.0, -3.0-pot_drop] )

   ndiag_list = np.array([nb*3,nb*4,nb*5,nb*6,nb*7,nb*8,nb,nb*2],dtype='i')
   for ndiag in ndiag_list:   
      print('ndiag=',ndiag)
      if (ndiag==0):
         ldiag=True
      else:
         ldiag=False

      v = np.zeros((nb*length,nb*length,nk), dtype=np.complex128)  
      ham = np.zeros((nb*length,nb*length,nk), dtype=np.complex128)  

      h00,h10 = wannierham.block_mat_def(kx=0.0, ky=0.0, kz=0.0,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
      ham[:,:,0] = wannierham.full_device_mat_def(ky=0.0,kz=0.0,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)
      v[:,:,0] = wannierham.full_device_bare_coulomb(ky=0.0,kz=0.0,length=length,eps=eps_screen,r0=r0,ldiag=ldiag,nb=nb,ns=ns,method='pointlike',n_range=n_range,wannier_center=wannier_center,cell=cell)

      for ix in range(length):
          ham[ix*nb:(ix+1)*nb,ix*nb:(ix+1)*nb,0] += np.diag(np.ones(nb)*pot[ix])

      dim_lead = np.ones(2)* nb*ns
      siglead = np.zeros((nb*ns,nb*ns,nen,2,nk), dtype=np.complex128)

      lead_h10 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
      lead_h10[:,:,0,0] = np.transpose( np.conjugate( h10 ) ) 
      lead_h10[:,:,1,0] = h10

      lead_h00 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
      lead_h00[:,:,0,0] = h00 + np.diag(np.ones(nb*ns)* pot[0])
      lead_h00[:,:,1,0] = h00 + np.diag(np.ones(nb*ns)* pot[-1])

      lead_coupling = np.zeros((nb*ns,nb*length,2,nk), dtype=np.complex128)
      lead_coupling[0:nb*ns,0:nb*ns,0,0] = lead_h10[:,:,0,0]
      lead_coupling[0:nb*ns,nb*(length-ns):nb*length,1,0] = lead_h10[:,:,1,0]

      nen_list = np.array([3000,5000,10000,15000],dtype='i')
      nsub_list = np.array([2,3],dtype='i')
      

      ID_list = np.zeros((2,nsub_list.shape[0],nen_list.shape[0]))


      for i, toten in enumerate(nen_list):
         print(toten)
         for j, nsub in enumerate(nsub_list):
            nen = toten
            if (nen*nsub > 30000):
               break
            energies = np.linspace(emin,emax,nen)
            print(nen,nsub)
            G_retarded,G_lesser,G_greater,W0,tr = gf_dense.solve_gw_3d(scba_tol=1e-4,niter=niter,nm_dev=nb*length,lx=Lx,length=length,spindeg=2.0,
                                                            temps=temp[0],tempd=temp[1],mus=mu[0],mud=mu[1],alpha_mix=0.5,
                                                            nen=nen,nsub=nsub,en=energies,nb=nb,ns=ns,nphiy=nky,nphiz=nkz,
                                                            ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
                                                            ndiag=ndiag,num_lead=2,flatband=False,output_files=True)
            print('NEGF done')                                                       
            print('current=', -np.sum(tr[:,0]) , np.sum(tr[:,1]))                                                       
            ID_list[:,j,i]= [-np.sum(tr[:,0]) , np.sum(tr[:,1]) ]    

      np.savez('run_ndiag'+str(ndiag)+'.npz', ID_list=ID_list,
                           nen_list=nen_list,
                           nsub_list=nsub_list)
