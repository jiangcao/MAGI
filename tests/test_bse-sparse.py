#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 2024

@author: Jiang Cao
"""


import numpy as np
from negf import gw_dense,bse_dense,bse_sparse,parameters_mod
from wannier import wannierham
import matplotlib.pyplot as plt
import os
os.environ["OMP_NUM_THREADS"] = "128"
# from sinv import sinv_tridiagonal_arrowhead
# from serinv.sequential import ddbtasinv
import time

if __name__=='__main__':   

   nb=14
   nx=21
   ny=1
   nz=1

   hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',
                                                                lreorder_axis=False,axis=[1,2,3],
                                                                nb=nb,nx=nx,ny=ny,nz=nz)

   Lz=L[2]
   Ly=L[1]
   Lx=L[0]

   ns = 2
   length = 10
   nen = 800
   nsub = 1
   nky=1
   nkz=1
   nk=nky*nkz
   niter=0
   eps_screen=1.0
   r0=3.0
   emin=-6.0
   emax= 2.0
   temp =  np.ones(2)* 300.0
   mu = np.array( [-2.25,-2.25 ] )
   light_polar = np.array( [1.0, 1.0, 0.0] )
   light_polar = light_polar / np.linalg.norm(light_polar)
   hw_phot = 2.5 # eV
   intensity = 1e10 # W/m2

   n_bose_phot = intensity / hw_phot / parameters_mod.e_charge
   print('n_bose_phot=',n_bose_phot)

   ndiag=nb*2

   if (ndiag==0):
       ldiag=True
   else:
       ldiag=False

   nm_dev=nb*length

   v = np.zeros((nm_dev,nm_dev,nk), dtype=np.complex128)  
   ham = np.zeros((nm_dev,nm_dev,nk), dtype=np.complex128)  
   M_phot = np.zeros((nm_dev,nm_dev,nk,nk), dtype=np.complex128)  

   h00,h10 = wannierham.block_mat_def(kx=0.0, ky=0.0, kz=0.0,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
   ham[:,:,0] = wannierham.full_device_mat_def(ky=0.0,kz=0.0,nb=nb,ns=ns,length=length,hr=hr,cell=cell,
                                               n_range=n_range)
   v[:,:,0] = wannierham.full_device_bare_coulomb(ky=0.0,kz=0.0,length=length,eps=eps_screen,r0=r0,
                                                  ldiag=ldiag,nb=nb,ns=ns,method='pointlike',n_range=n_range,
                                                  wannier_center=wannier_center,cell=cell)
   pmn = wannierham.calc_momentum_operator(method='approx',nb=nb,nx=nx,ny=ny,nz=nz,hr=hr,cell=cell,
                                           n_range=n_range,wannier_center=wannier_center,
                                           rmn=np.zeros((nm_dev,nm_dev,3)))
   pij = wannierham.w90_momentum_full_device(ky=0,kz=0,length=length,ns=ns,n_range=n_range,nb=nb,cell=cell,
                                             pmn=pmn)
   
   for i in range(3):
      M_phot[:,:,0,0] += pij[:,:,i] * light_polar[i] 
   
   print('max |M_phot|=',np.max(np.abs(M_phot)))

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

   energies = np.linspace(emin,emax,nen)
   print('num energies=',nen,', num sub-energies=',nsub)
   
   dE=energies[1]-energies[0]
   n_phot = int(hw_phot / dE)
   nstep=4
   eps_M=np.zeros(nen//nstep,dtype='complex')
   eps_M_sinv=np.zeros(nen//nstep,dtype='complex')
   eps_M_sinv2=np.zeros(nen//nstep,dtype='complex')
   eps_en=np.zeros(nen//nstep)

   Ephmin = 0.5
   nnop= int( (4.0-Ephmin)/dE/nstep )
   nops=np.arange(nnop)*nstep + int(Ephmin / dE)
   print('num E_opt = ' , nnop)
   print('E_optical = ' , nops * dE)

   nm_dev=nb*length

   start_t = time.time()
   ( G_retarded,G_lesser,G_greater,
        tr,te,P_r )= bse_sparse.bse_sparse_solve_scba(method='sum',niter=0,nm_dev=nm_dev,lx=Lx,length=length,spindeg=2.0,temp=temp,mu=mu,
        alpha_mix=0.5,nen=nen,en=energies,nops=nops,nnop=nnop,nb=nb,ns=ns,ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
        ndiag=ndiag,encut=encut,egap=egap,vertex=False,bse=False,flatband=False,output_files=True,
        inj_photon=True,nphot=n_phot,m_phot=M_phot,n_bose_phot=n_bose_phot) 
   print('done')                                                       
   print('current=', -np.sum(tr[:,0]) , np.sum(tr[:,1]))                                                       

   ID_list = [-np.sum(tr[:,0]) , np.sum(tr[:,1]) ]    

   Gr_diag = G_retarded.diagonal()
   Gn_diag = G_lesser.diagonal()
   ldos = np.imag(Gr_diag)
   ndos = np.imag(Gn_diag)

   for iop in range(nnop):
      epsilon_M = np.eye(nm_dev) -  v[:,:,0] @ P_r[:,:,iop]
      eps_M_sinv[iop] = np.sum( epsilon_M[ nm_dev//2, nb*ns:(nm_dev-nb*ns) ] )
      eps_en[iop]=nops[iop]*dE
      print('- E=',eps_en[iop],np.abs( np.imag(eps_M_sinv[iop]) ) )   

   finish_t = time.time()
   print("! BSE solver takes %s second for all optical energies " % (finish_t - start_t))
   bse_time = (finish_t - start_t)

   plt.plot(eps_en,np.abs(np.imag(eps_M_sinv)))
   plt.savefig('fig_epsilon2.png')

   np.savez('data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'.npz', 
                        ID_list=ID_list,
                        tr=tr,
                        te=te,
                        energies=energies,
                        ldos=ldos,
                        ndos=ndos,
                        eps_M=eps_M,
                        eps_M_sinv=eps_M_sinv,
                        eps_en=eps_en,
                        bse_time=bse_time,
                        intensity=intensity,
                        n_phot=n_phot,
                        Lx=Lx,
                        cell=cell,
                        egap=egap)

