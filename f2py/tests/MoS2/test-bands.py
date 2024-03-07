
import numpy as np
from util import linalg
from negf import gf_dense, fft_mod
from wannier import wannierham
import matplotlib.pyplot as plt

if __name__=='__main__':   

   nb=11
   nx=1
   ny=29
   nz=29

   hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname='ham_dat',lreorder_axis=True,axis=[3,2,1],nb=nb,nx=nx,ny=ny,nz=nz)

   ns = 1
   length = 1 

   nky=61
   nkz=61
   nk=nky*nkz

   Lz=L[2]
   Ly=L[1]
   Lx=L[0]
   
   kz_min = -1.5*np.pi/Lz
   kz_max =  1.5*np.pi/Lz
   ky_min = -1.5*np.pi/Lz
   ky_max =  1.5*np.pi/Lz

   dkz=(kz_max-kz_min) / nkz
   dky=(ky_max-ky_min) / nky

   ham = np.zeros((nb*length,nb*length,nk), dtype='complex')  
   
   ek = np.zeros((nb*length,nk), dtype='double')     

   for iky in range(nky):
      for ikz in range(nkz):
         ik = ikz + iky*nkz
         ky=ky_min + dky*iky
         kz=kz_min + dkz*ikz
         if (nkz==1):
            kz=0.0
         if (nky==1):
            ky=0.0

         ham[:,:,ik] = wannierham.full_device_mat_def(ky=ky,kz=kz,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)

         ek[:,ik] = np.real(np.linalg.eigvalsh(ham[:,:,ik]))
         ek[:,ik] = np.sort(ek[:,ik])

   x=np.linspace(kz_min,kz_max, nkz)
   y=np.linspace(ky_min,ky_max, nky)

   for ib in range(nb):
      val=np.reshape( ek[ib,:], (nky,nkz))
      plt.contourf(x,y,val)
      plt.colorbar()
      plt.title('band #'+str(ib+1))
      plt.show()



         

   