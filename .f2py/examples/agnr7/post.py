import matplotlib.pyplot as plt
import numpy as np

ndiag= 14
path='/usr/scratch2/tortin16/jiacao/BSE_calc/agnr7/python/'
path='./'
f=np.load(path+'/checkpoint_ndiag'+str(ndiag)+'.npz')

#plt.plot(f['eps_en'],np.abs(np.imag(f['eps_M'])))
#plt.show()

plt.plot(np.imag(f['Kdiag']))
plt.plot(np.real(f['Kdiag']))
plt.show()
#
plt.spy(f['Ktip'])
plt.show()
Ltip = f['Ltip']


plt.spy(f['Ldiag'])
plt.show()

plt.spy(f['Lupper'])
plt.show()

plt.spy(f['Llower'])
plt.show()

plt.spy(f['Lupperarrow'])
plt.show()

plt.spy(f['Llowerarrow'])
plt.show()

plt.spy(f['Ltip'])
plt.show()

#
#
#
plt.spy(f['Adiag'])
plt.show()

plt.spy(f['Aupper'])
plt.show()

plt.spy(f['Alower'])
plt.show()

plt.spy(f['Aupperarrow'])
plt.show()

plt.spy(f['Alowerarrow'])
plt.show()

plt.spy(f['Atip'])
plt.show()


