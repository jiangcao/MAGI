
import numpy as np
from util import linalg
from negf import gf_dense, fft_mod
from wannier import wannierham
import matplotlib.pyplot as plt
import time

if __name__=='__main__':   

   nb=11
   nx=1
   ny=29
   nz=29

   hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',lreorder_axis=True,axis=[3,2,1],nb=nb,nx=nx,ny=ny,nz=nz)

   ns = 1
   length = 1 

   nky=301
   nkz=301
   nk=nky*nkz

   Lz=L[2]
   Ly=L[1]
   Lx=L[0]
   
   kz_min = -1.5*np.pi/Lz
   kz_max =  1.5*np.pi/Lz
   ky_min = -1.5*np.pi/Ly
   ky_max =  1.5*np.pi/Ly

   dkz=(kz_max-kz_min) / nkz
   dky=(ky_max-ky_min) / nky

   ham = np.zeros((nb*length,nb*length), dtype='complex')  
   
   ek = np.zeros((nb*length,nk), dtype='double')     

   start = time.time()
   for iky in range(nky):
      for ikz in range(nkz):
         ik = ikz + iky*nkz
         ky=ky_min + dky*iky
         kz=kz_min + dkz*ikz
         if (nkz==1):
            kz=0.0
         if (nky==1):
            ky=0.0

         ham = wannierham.full_device_mat_def(ky=ky,kz=kz,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)

         ek[:,ik] = np.real(np.linalg.eigvalsh(ham))
         ek[:,ik] = np.sort(ek[:,ik])
   finish = time.time()      
   print('python done in seconds: ', finish-start)

   x=np.linspace(kz_min,kz_max, nkz)
   y=np.linspace(ky_min,ky_max, nky)

   nvb = n_range[-2]

   for ib in range(nvb-1,nvb+1):
      val=np.reshape( ek[ib,:], (nky,nkz))
      plt.contourf(x,y,val)      
      plt.title('band #'+str(ib+1))
      plt.savefig('p'+str(ib)+'.png')

   start = time.time()
   ek2 = wannierham.plot_bands_bz(nkx=1,nky=nky,nkz=nkz,hr=hr,nb=nb,length=L,cell=cell,n_range=n_range)
   finish = time.time()      
   print('fortran done in seconds: ', finish-start)

   for ib in range(nvb-1,nvb+1):
      val=np.reshape( ek2[ib,:], (nky,nkz))
      plt.contourf(x,y,val)
      plt.title('band #'+str(ib+1))
      plt.savefig('f'+str(ib)+'.png')

   print('max difference in Ek: ', np.max(abs( ek2 - ek )) )
         

   