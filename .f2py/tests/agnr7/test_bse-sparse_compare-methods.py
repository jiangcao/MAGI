#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 2024

@author: Jiang Cao
"""


import numpy as np
from negf import gw_dense,bse_dense,bse_sparse,parameters_mod
from wannier import wannierham
import matplotlib.pyplot as plt
import os
os.environ["OMP_NUM_THREADS"] = "128"
# from sinv import sinv_tridiagonal_arrowhead
from serinv.sequential import ddbtasinv
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
   nen = 1800
   nsub = 1
   nky=1
   nkz=1
   nk=nky*nkz
   niter=0
   eps_screen=1.0
   r0=3.0
   emin=-12.0
   emax= 6.0
   temp =  np.ones(2)* 300.0
   mu = np.array( [-2.25,-2.25 ] )
   light_polar = np.array( [1.0, 0.0, 0.0] )
   hw_phot = 2.0 # eV
   intensity = 1e13 # W/m2
   n_bose_phot = intensity / hw_phot / parameters_mod.e_charge
   pre_factor = ((parameters_mod.hbar/parameters_mod.m0_kg)**2 * parameters_mod.e_charge**2)/(2.0*parameters_mod.epsilon0*parameters_mod.light_speed**3*hw_phot*parameters_mod.e_charge)
   print('pre_factor=',pre_factor)
   print('n_bose_phot=',n_bose_phot)

   ndiag=nb*2

   if (ndiag==0):
       ldiag=True
   else:
       ldiag=False

   nm_dev=nb*length

   v = np.zeros((nm_dev,nm_dev,nk), dtype=np.complex128)  
   ham = np.zeros((nm_dev,nm_dev,nk), dtype=np.complex128)  
   M_phot = np.zeros((nm_dev,nm_dev,nk), dtype=np.complex128)  

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
      M_phot[:,:,0] = pij[:,:,i] * light_polar[i] * pre_factor

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
   print(nen,nsub)
   
   dE=energies[1]-energies[0]
   n_phot = int(hw_phot / dE)
   nstep=4
   eps_M=np.zeros(nen//nstep,dtype='complex')
   eps_M_sinv=np.zeros(nen//nstep,dtype='complex')
   eps_M_sinv2=np.zeros(nen//nstep,dtype='complex')
   eps_en=np.zeros(nen//nstep)

   Ephmin = 0.8
   nnop= int( (4.0-Ephmin)/dE/nstep )
   nops=np.arange(nnop)*nstep + int(Ephmin / dE)
   print('Num E_opt = ' , nnop)
   print('E_optical = ' , nops * dE)


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


   print("=============================")
   print("--- Start full BSE solver ---")
   print("=============================")
   start_t = time.time()

   for iop in range(nnop):
      print("iop=",iop+1,"/",nnop,"E_optical=",nops[iop]*dE)      
      
      P_r,nn = bse_dense.bse_fullsolve(
               alpha=0.99,spindeg=2.0,nm_dev=nm_dev,ndiag=ndiag,nen=nen,
               en=energies,nop=nops[iop],
               g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,
               w=W0_r,v=v )      

      epsilon_M = np.eye(nm_dev) - v[:,:,0] @ ( - 1j* np.abs(np.imag(P_r)) )

      eps_M[iop] = np.sum( epsilon_M[ nm_dev//2, nb*ns:(nm_dev-nb*ns) ] )
      eps_en[iop]=nops[iop]*dE
      print('- Eopt=',eps_en[iop],'epsilon_2=', np.abs( np.imag(eps_M[iop]) ) )
      print(" ")      

   finish_t = time.time()
   print("! dense BSE solver takes %s second for all optical energies " % (finish_t - start_t))
   


   
   print("===============================")
   print("--- Start sparse BSE solver ---")
   print("===============================")   
      
   start_t = time.time()

   P_r,Sig_r,Sig_l,Sig_g = bse_sparse.bse_sparse_solve(
                  method='sum',alpha=0.99,spindeg=2.0,
                  nm_dev=nm_dev,ndiag=ndiag,nen=nen,nsub=1,
                  en=energies,nops=nops,nnop=nnop,nk=1,
                  g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,
                  w=W0_r,v=v,
                  solve_sigma=True,with_vertex=True,nb=nb,ns=ns)        
   
   finish_t = time.time()   

   for iop in range(nnop):
      epsilon_M = np.eye(nm_dev) -  v[:,:,0] @ P_r[:,:,iop]
      eps_M_sinv[iop] = np.sum( epsilon_M[ nm_dev//2, nb*ns:(nm_dev-nb*ns) ] )
      eps_en[iop]=nops[iop]*dE
      print('- Eopt=',eps_en[iop],np.abs( np.imag(eps_M_sinv[iop]) ),
            ', ref=', np.abs( np.imag(eps_M[iop]) ),
            ', err=', (np.abs( np.imag(eps_M_sinv[iop]) ) - np.abs( np.imag(eps_M[iop]) ))/np.abs( np.imag(eps_M[iop]) ) )     

   print("")
   print("! sparse BSE solver takes %s second for %s optical energies " % ((finish_t - start_t) , nnop))            
   



   
   # print("--- Start pre-processing for solver ---")

   # N,nnz,table,blocksize,num_blocks = bse_sparse.bse_sparse_pre(nm_dev=nm_dev,ndiag=ndiag)
   # # resize
   # table=table[:,:N]

   # Ldiag,Lupper,Llower,Lupperarrow,Llowerarrow,Ltip,Ktip,Kdiag = bse_sparse.bse_sparse_build(
   #             method='sum',
   #             alpha=0.99,spindeg=2.0,nm_dev=nm_dev,ndiag=ndiag,nen=nen,en=energies,nop=nops,nnop=nnop,
   #             g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,w=W0_r,v=v,
   #             blocksize=blocksize,num_blocks=num_blocks,n=N,table=table)



   # print("================================")
   # print("--- Start solver with SerinV ---")
   # print("================================")
   
   # start_t = time.time()

   # for iop in range(nnop): 
   #    Adiag,Aupper,Alower,Alowerarrow,Aupperarrow,Atip = bse_sparse.bse_sparse_build_system(
   #        blocksize=blocksize, num_blocks=num_blocks, nnop=nnop, iop=iop+1, nm_dev=nm_dev,
   #        ldiag=Ldiag, lupper=Lupper, llower=Llower, llowerarrow=Llowerarrow, 
   #        lupperarrow=Lupperarrow, ltip=Ltip, kdiag=Kdiag, ktip=Ktip )
      
   #    newAdiag,newAupper,newAlower,newAlowerarrow,newAupperarrow = bse_sparse.reshape_bta_block2stack(
   #        num_blocks=num_blocks,blocksize=blocksize,nm_dev=nm_dev,
   #        adiag=Adiag, aupper=Aupper, alower=Alower, alowerarrow=Alowerarrow, aupperarrow=Aupperarrow)      

   #    bse_sparse.bse_sparse_check_system( tol=1e-6, 
   #        alpha=0.99,spindeg=2.0,nm_dev=nm_dev,ndiag=ndiag,nen=nen,en=energies,
   #        nop=nops,nnop=nnop,iop=iop+1,
   #        blocksize=blocksize,num_blocks=num_blocks,n=N,table=table,
   #        g_lesser=G_lesser,g_greater=G_greater,g_retarded=G_retarded,w=W0_r,v=v, 
   #        adiag=Adiag, aupper=Aupper, alower=Alower, alowerarrow=Alowerarrow, aupperarrow=Aupperarrow, atip=Atip,
   #        ldiag=Ldiag, lupper=Lupper, llower=Llower, llowerarrow=Llowerarrow, lupperarrow=Lupperarrow, ltip=Ltip )


   #    try:
   #       start_t2 = time.time()
   #       (
   #          X_diagonal_blocks,
   #          X_lower_diagonal_blocks,
   #          X_upper_diagonal_blocks,
   #          X_arrow_bottom_blocks,
   #          X_arrow_right_blocks,
   #          X_arrow_tip_block,
   #       ) = ddbtasinv(
   #          A_diagonal_blocks=newAdiag,
   #          A_lower_diagonal_blocks=newAlower,
   #          A_arrow_bottom_blocks=newAlowerarrow,
   #          A_arrow_tip_block=Atip,
   #          A_upper_diagonal_blocks=newAupper,
   #          A_arrow_right_blocks=newAupperarrow,
   #       )
   #       finish_t2 = time.time()
   #       print(" SerinV takes %s second for one SINV " % (finish_t2 - start_t2))

   #    except Exception as error:

   #       print("! Fails with error: ", error ) 

   #    tmp = - 1j* X_arrow_tip_block @ Ltip[:,:,iop]
   #    for i in range(num_blocks):
   #       tmp += - 1j* X_arrow_bottom_blocks[i,:,:] @ Lupperarrow[i*blocksize:(i+1)*blocksize,:,iop]

   #    P_r = np.zeros((nm_dev,nm_dev),dtype='complex')
   #    for row in range(nm_dev):
   #       for col in range(nm_dev):
   #          fliped_row = nm_dev - col - 1
   #          fliped_col = nm_dev - row - 1 
   #          i=table[0,row] - 1
   #          k=table[0,col] - 1
   #          P_r[i,k] =  tmp[fliped_row,fliped_col]                
      
   #    epsilon_M = np.eye(nm_dev) -  v[:,:,0] @ P_r

   #    eps_M_sinv2[iop] = np.sum( epsilon_M[ nm_dev//2, nb*ns:(nm_dev-nb*ns) ] )
      
   #    print('- Eopt=',eps_en[iop],'epsilon_2=', np.abs( np.imag(eps_M_sinv2[iop]) ),', ref=', np.abs( np.imag(eps_M[iop]) ) )
                  
   # finish_t = time.time()

   # print("! SerinV solver takes %s second in total " % (finish_t - start_t))




   print("save checkpoint.")
   np.savez('chkpt_len'+str(length)+'_ndiag'+str(ndiag)+'.npz', 
                        ID_list=ID_list,
                        tr=tr,
                        te=te,
                        energies=energies,
                        ldos=ldos,
                        ndos=ndos,
                        eps_M=eps_M,
                        eps_en=eps_en,
                        eps_M_sinv=eps_M_sinv,
                        eps_M_sinv2=eps_M_sinv2
            )

