
import numpy as np
from negf import rgf
from wannier import wannierham
import matplotlib.pyplot as plt

if __name__=='__main__':   

   nb=104
   nx=15
   ny=1
   nz=1

   hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',lreorder_axis=False,axis=[1,2,3],nb=nb,nx=nx,ny=ny,nz=nz)

   ns = 2
   length = 5
   nen = 100 # number of energy points
   nsub = 3 # number of Legendre nodes in each interval
   nky=1
   nkz=1
   nk=nky*nkz
   niter=10
   eps_screen=2.5
   r0=3.0
   emin=-10.0
   emax=4.0

   Lz=L[2]
   Ly=L[1]
   Lx=L[0]
   dkz=2.0*np.pi/Lz / nkz
   dky=2.0*np.pi/Ly / nky

   energies = np.linspace(emin,emax,nen)


   dev_h00 = np.zeros((nb*ns,nb*ns,length), dtype='complex')
   dev_s00 = np.zeros((nb*ns,nb*ns,length), dtype='complex')
   dev_h10 = np.zeros((nb*ns,nb*ns,length+1), dtype='complex')

   sig_l = np.zeros((nb*ns,nb*ns,length), dtype='complex')
   sig_r = np.zeros((nb*ns,nb*ns,length), dtype='complex')

   mul = np.ones((nb*ns,nb*ns), dtype='double') * -1.5
   mur = np.ones((nb*ns,nb*ns), dtype='double') * -1.6

   templ  = np.ones((nb*ns,nb*ns), dtype='double') * 300.0
   tempr  = np.ones((nb*ns,nb*ns), dtype='double') * 300.0

   tr  = np.zeros(nen, dtype='double')
   tre  = np.zeros(nen, dtype='double')

   nm = np.ones(length, dtype='int') * nb*ns

   ky=0.0
   kz=0.0
         
   h00,h10 = wannierham.block_mat_def(kx=0.0, ky=ky, kz=kz,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)
   
   for ix in range(length):
      dev_h00[:,:,ix] = h00
      dev_h10[:,:,ix] = h10
      dev_s00[:,:,ix] = np.eye(nb*ns)
   dev_h10[:,:,length] = h10
   
   for ie in range(nen):
      print(ie, energies[ie])
      G_r, G_lesser, G_greater, Jdens, tr[ie], tre[ie] = rgf.rgf_std(nx=length,mm=nb*ns,nm=nm, en=energies[ie], mul=mul, mur=mur, templ=templ, tempr=tempr, 
                                                            hii=dev_h00, h1i=dev_h10, sii=dev_s00, 
                                                            sigma_lesser_ph=sig_l, sigma_r_ph=sig_r, verbose=False)
      
   plt.plot(energies,tr)
   plt.plot(energies,tre)
   plt.show()
      