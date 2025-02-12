#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 2024

@author: Jiang Cao
"""


import numpy as np
import numpy.linalg as npla
from negf import gw_dense,bse_sparse
from wannier import wannierham
import os
os.environ["OMP_NUM_THREADS"] = "28"
os.environ["MKL_NUM_THREADS"] = "28"
os.environ["OPENBLAS_NUM_THREADS"] = "28"
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
   length = 30
   nen = 1800
   nstep= 1
   nsub = 1
   nky=1
   nkz=1
   nk=nky*nkz
   niter=0
   eps_screen=12.0
   r0=3.0
   emin=-12.0
   emax= 6.0
   temp =  np.ones(2)* 300.0
   # mu = np.array( [-2.25,-2.25 ] )
   mu = np.array( [-1.9,-1.9 ] )

   pot_drop = 0.1 # V
   pot= np.concatenate( [[0], -pot_drop*np.arange(length-2)/(length-3), [-pot_drop]] )
#    print(pot)s
   mu[1]=mu[1]+pot[-1]

   ndiag= nb * 1
   nm_dev=nb * length

   if (ndiag==0):
       ldiag=True
   else:
       ldiag=False

   # inject light 
   light_polar = np.array( [1.0, 1.0, 0.0] )
   light_polar = light_polar / np.linalg.norm(light_polar)
   hw_phots = np.arange(10,40) * 0.1 #2.5 # eV
   intensity = 1e16 # W/m2

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
   
   for i in range(length):
       for ib in range(nb):
           ham[i*nb+ib,i*nb+ib,0] += pot[i]

   for i in range(3):
      M_phot[:,:,0,0] += pij[:,:,i] * light_polar[i] 
   
   # print('max |M_phot|=',np.max(np.abs(M_phot)))

   dim_lead = np.ones(2)* nb*ns
   siglead = np.zeros((nb*ns,nb*ns,nen,2,nk), dtype=np.complex128)

   lead_h10 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
   lead_h10[:,:,0,0] = np.transpose( np.conjugate( h10 ) ) 
   lead_h10[:,:,1,0] = h10

   lead_h00 = np.zeros((nb*ns,nb*ns,2,nk), dtype=np.complex128)
   lead_h00[:,:,0,0] = h00
   lead_h00[:,:,1,0] = h00

   for ib in range( nb * ns ): 
       lead_h00[ib,ib,1,0] += pot[-1]

   lead_coupling = np.zeros((nb*ns,nb*length,2,nk), dtype=np.complex128)
   lead_coupling[0:nb*ns,0:nb*ns,0,0] = lead_h10[:,:,0,0]
   lead_coupling[0:nb*ns,nb*(length-ns):nb*length,1,0] = lead_h10[:,:,1,0]

   egap=1.67
   encut=[7.0,7.0]

   energies = np.linspace(emin,emax,nen)
   print('num energies=',nen,', num sub-energies=',nsub, flush=True)
   

   for hw_phot in hw_phots:

      n_bose_phot = intensity / hw_phot / 1.6e-19 
      print('hw phot=',hw_phot, flush=True)

      # define optical energy grid
      dE=energies[1]-energies[0]
      n_phot = int(hw_phot / dE)
      
      Ephmin = 0.0
      nnop= int( (4.0-Ephmin)/dE/nstep )
      nops=np.arange(nnop)*nstep + int(Ephmin / dE)
      print('num E_opt = ' , nnop, flush=True)
      # print('E_optical = ' , nops * dE)

      eps_M=np.zeros(nnop,dtype='complex')
      eps_M_sinv=np.zeros(nnop,dtype='complex')
      eps_2=np.zeros(nnop,dtype='complex')
      eps_rpa=np.zeros(nnop,dtype='complex')
      eps_M_sinv2=np.zeros(nnop,dtype='complex')
      eps_en=np.zeros(nnop)

      start_t = time.time()

      ( G_retarded,G_lesser,G_greater,current,transmission,P_r,P0_r ) = bse_sparse.bse_sparse_solve_scba(
         method='fft',niter=0,nm_dev=nm_dev,lx=Lx,length=length,spindeg=2.0,temp=temp,mu=mu,
         alpha_mix=0.5,nen=nen,en=energies,nops=nops,nnop=nnop,nb=nb,ns=ns,
         ham=ham,h00lead=lead_h00,h10lead=lead_h10,t=lead_coupling,v=v,
         ndiag=ndiag,encut=encut,egap=egap,
         vertex=True,bse_sigma=True,flatband=False,output_files=True,
         inj_photon=True,nphot=n_phot,m_phot=M_phot,n_bose_phot=n_bose_phot ) 

      # (G_retarded,G_lesser,G_greater,W0_r,tr) = gw_dense.solve_gw_3d(
      #                      niter=0,scba_tol=1.0,nm_dev=nm_dev,lx=Lx,length=length,spindeg=2.0,
      #                      temps=temp[0],tempd=temp[1],mus=mu[0],mud=mu[1],
      #                      alpha_mix=0.5,nen=nen,nsub=1,en=energies,nb=nb,ns=ns,
      #                      nphiy=nky,nphiz=nkz,ham=ham,h00lead=lead_h00,h10lead=lead_h10,
      #                      t=lead_coupling,v=v,
      #                      ndiag=ndiag,num_lead=2,flatband=False,output_files=True)
      
      # (G_retarded,G_lesser,G_greater,
      #    Sig_retarded_new,Sig_lesser_new,Sig_greater_new,
      #    current,transmission,W0_r,W0_lesser,W0_greater) = gw_dense.solve_gw_1d_memsaving(
      #          niter=0,nm_dev=nm_dev,lx=Lx,length=length,spindeg=2.0,temp=temp,mu=mu,
      #          alpha_mix=0.5,nen=nen,en=energies,nb=nb,ns=ns,ham=ham,h00lead=lead_h00,h10lead=lead_h10,
      #          t=lead_coupling,v=v,
      #          ndiag=ndiag,encut=encut,egap=egap,vertex=False,bse=False,flatband=False,output_files=True)
      
      ID_list = [ -np.sum(current[:,0]) , np.sum(current[:,1]) ]    
      print('done', flush=True)                                                       
      print('current=', ID_list , flush=True)                                                      

      Gr_diag = G_retarded.diagonal()
      Gn_diag = G_lesser.diagonal()
      ldos = np.imag(Gr_diag)
      ndos = np.imag(Gn_diag)
      # # fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'.npz'
      fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'_potdrop'+str(pot_drop)+'_eps'+str(eps_screen)+'.npz'
      # print("save data in ",fname)
      # np.savez(fname, 
      #          current=current,
      #          transmission=transmission,
      #          energies=energies, 
      #          ldos=ldos,
      #          ndos=ndos,
      #          Lx=Lx,
      #          cell=cell,
      #          egap=egap)
      
      # (P_r,P0_r,Sig_r,Sig_l,Sig_g) = bse_sparse.bse_sparse_solve(
      #                method='fft',alpha=0.99,spindeg=2.0,
      #                nm_dev=nm_dev,ndiag=ndiag,nen=nen,nsub=1,
      #                en=energies,nops=nops,nnop=nnop,nk=1,
      #                g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,
      #                w=W0_r,v=v,solve_sigma=True,with_vertex=True,nb=nb,ns=ns)    

      finish_t = time.time()
      print("! BSE solver takes %s second for all optical energies " % (finish_t - start_t), flush=True)
      bse_time = (finish_t - start_t)

      for iop in range(nnop):
         epsilon_2 = np.eye(nm_dev) -  v[:,:,0] @ P0_r[:,:,iop]
         epsilon_M = np.eye(nm_dev) -  v[:,:,0] @ P_r[:,:,iop]
         epsilon_rpa = npla.inv(np.eye(nm_dev) -  v[:,:,0] @ P0_r[:,:,iop])
         eps_M_sinv[iop] = np.sum( epsilon_M[ nm_dev//2, nb*ns:(nm_dev-nb*ns) ] )
         eps_rpa[iop] = np.sum( epsilon_rpa[ nm_dev//2, nb*ns:(nm_dev-nb*ns) ] )
         eps_2[iop] = np.sum( epsilon_2[ nm_dev//2, nb*ns:(nm_dev-nb*ns) ] )
         eps_en[iop]=nops[iop]*dE
         # print('- E=',eps_en[iop],np.abs( np.imag(eps_M_sinv[iop]) ),np.abs( np.imag(eps_rpa[iop]) ) )   

      # fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'.npz'
      print("save data in ",fname, flush=True)
      np.savez(fname, 
               ID_list=ID_list,
               current=current,
               transmission=transmission,
               energies=energies,
               ldos=ldos,
               ndos=ndos,
               Lx=Lx,
               cell=cell,
               egap=egap,
               # eps_M=eps_M,
               eps_2=eps_2,
               eps_M_sinv=eps_M_sinv,
               eps_rpa=eps_rpa,
               eps_en=eps_en,
               bse_time=bse_time,
               intensity=intensity,
               n_phot=n_phot,            
               )

