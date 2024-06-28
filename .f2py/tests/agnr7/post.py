import matplotlib.pyplot as plt
import numpy as np

ndiag= 28
length=10
nen = 800
n_phot = 613

path='./'
fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'.npz'

f=np.load(path+fname)

eps_en = f['eps_en']
eps_M_sinv = f['eps_M_sinv']

plt.plot(eps_en,np.imag(eps_M_sinv),'.-')
plt.show()
plt.savefig('fig_epsilon2.png')
