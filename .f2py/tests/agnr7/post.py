import matplotlib.pyplot as plt
import numpy as np

ndiag= 28
length=10
path='./'
fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'.npz'
f=np.load(path+fname)

eps_en = f['eps_en']
eps_M_sinv = f['eps_M_sinv']

plt.plot(eps_en,np.imag(eps_M_sinv),'.-')
plt.show()
