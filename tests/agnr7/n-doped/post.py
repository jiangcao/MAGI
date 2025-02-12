import matplotlib.pyplot as plt
import numpy as np
ndiag= 14
length=30
nen = 1800
n_phot = 99
mu = np.array( [-1.9,-1.9 ] )
eps_screen=4.0

pot_drop = 0.1 # V

pot= np.concatenate( [[0], -pot_drop*np.arange(length-2)/(length-3), [-pot_drop]] )
mu[1]=mu[1]+pot[-1]

hw_phots = np.arange(10,35) * 0.1 #2.5 # eV

# fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'.npz'
fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'_potdrop'+str(pot_drop)+'.npz'
# fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'.npz'
# fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'_potdrop'+str(pot_drop)+'_eps'+str(eps_screen)+'.npz'
f=np.load('./'+fname)
energies = f['energies']
Egap = 2.9
IDs1 = []
IDs2 = []

for hw_phot in hw_phots:
    dE=energies[1]-energies[0]
    n_phot = int(hw_phot / dE)

    # fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'.npz'
    fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'_potdrop'+str(pot_drop)+'.npz'
    # fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'.npz'
    # fname='data_len'+str(length)+'_ndiag'+str(ndiag)+'_nen'+str(nen)+'_nphot'+str(n_phot)+'_mul'+str(mu[0])+'mur'+str(mu[1])+'_potdrop'+str(pot_drop)+'_eps'+str(eps_screen)+'.npz'

    f=np.load('./'+fname)

    ID = f['ID_list']
    print(hw_phot, ID)
    if (hw_phot < Egap):
        factor = 1e-3
    else:        
        factor = 1.0

    IDs1.append(ID[0]/factor)
    IDs2.append(ID[1]/factor)

    # plot epsilon spec
    eps_en = f['eps_en']
    eps_M_sinv = f['eps_M_sinv']
    eps_rpa = f['eps_rpa']            
    plt.subplot(2,2,1)
    # plt.plot(eps_en,np.abs(np.imag(eps_rpa)),'--',label='RPA')
    plt.plot(eps_en,np.abs(np.imag(eps_M_sinv)),'-',label='BSE')    
    plt.ylim(1e-3,2)
    plt.xlim(np.min(hw_phots), np.max(hw_phots))
    plt.arrow(hw_phot,2,0,-1,shape='full',width=0.05,color='r')
    plt.arrow(Egap,-0.5,0,2.5)
    plt.xlabel('Photon Energy')
    plt.yscale('log')
    # plt.legend()
    # fig.savefig('fig_epsilon2'+fname+'.png')
    # plt.show()

    # plot contact current spec
    tr = f['current'] / factor
    plt.subplot(2,2,2)    
    plt.fill_between(energies, -tr[:,0],label='Left')
    plt.fill_between(energies, tr[:,1],label='Right')
    plt.legend()
    plt.xlabel('Electron Energy')
    plt.xlim(-6,2)
    

    # plot LDOS colormap
    ldos = f['ldos']*2
    plt.subplot(2,2,4)
    plt.plot(energies,-ldos[:,0],label='Left')
    plt.plot(energies,-ldos[:,-4],label='Right')
    plt.legend()
    plt.xlabel('Electron Energy')
    plt.xlim(-6,2)
    plt.ylim(0,3)
    plt.arrow(-1.8-hw_phot/2.0,1.0,hw_phot,0,width=0.1,color='r',length_includes_head=True)
    
    plt.subplot(2,2,3)
    plt.plot(hw_phots[:len(IDs1)],(np.abs(np.array(IDs1)+np.array(IDs2)))/2.0,'o-')
    plt.yscale('log')
    plt.xlim(np.min(hw_phots), np.max(hw_phots))
    plt.xlabel('Photon Energy')

    plt.tight_layout()
    plt.show()