

# %%
import numpy as np
import pickle
from scipy.linalg import eig
from sklearn.neighbors import NearestNeighbors
from negf import gw_dense, fft_mod
from wannier import wannierham
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import scipy.constants as sc
from numba import njit, prange
import math
from scipy import interpolate



nb=16
nb_L=8
nx=41
ny=41
nz=10
ns_L=2


folder='/usr/scratch/mont-fort1/sem24f6/MAGI/f2py/QE_local-tests/'
hamiltonian='gr_hbn_hr.dat'
hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname=folder+hamiltonian,lreorder_axis=False,axis=[3,2,1],nb=nb,nx=nx,ny=ny,nz=nz)
folder_save='/usr/scratch/mont-fort1/sem24f6/MAGI/f2py/QE_local-tests/Save_heteroside2/'
label='hetero0'


# %%
#parameters Wannier centers
alt=[0,3] #Altitude of different layers
delta_alt=1. #delta around those layer to find the Wannier centers

interaction_photon=False #Activate or not the scattering 
interaction_phonon=False
# Global variables for the GF calculation
ns = ns_L
length = 2
nen = 2000 # number of energy points
nky=12
nkz=12
nkx=12
nk=nky*nkz


emin=-15.0
emax=30.0
energies = np.linspace(emin,emax,nen)
dE=energies[1]-energies[0]
#definition of k-space vectors, real space vectores, ...
Lz=L[2]
Ly=L[1]
Lx=L[0]
dkz=2.0*np.pi/Lz / (nkz)
dky=2.0*np.pi/Ly / (nky)

alpha=cell[:,0]
beta=cell[:,1]
gamma=cell[:,2]
V=np.dot(np.cross(beta,gamma),alpha)
b_1=np.cross(beta,gamma)/V*2*np.pi
b_2=np.cross(gamma,alpha)/V*2*np.pi
b_3=np.cross(alpha,beta)/V*2*np.pi
xhat=alpha/np.linalg.norm(alpha)
yhat=-np.cross(xhat,gamma)
yhat=yhat/np.linalg.norm(yhat)
zhat=gamma/np.linalg.norm(gamma)


K=((b_2+b_3)/3)
Gamma=np.zeros(3)
Z=b_1/2
Q=(b_2+b_3)/3+b_1/2
M=(b_3/2)



#parameters of the leads
dim_lead = np.ones(2)* nb*ns
temp =  np.ones(2)* 300.0
EF=8
DEF=0.0
mu = np.array( [EF+DEF, EF-DEF] )


# Variables of the phonon interaction calculations
Dop2=1e-4 #interaction electron - phonon optic
Dac2=np.array([1e-4,1e-4,1e-4]) #interaction electron - phonon acoustic
Ac_fact=np.array([84.9385484086987e10,1120.28324417723e10,649.5300760665194e10]) #Proportionnality frequencies/norm of q for acoustic modes


frequences=np.array([47696980067800.9,47696980067800.9,47696980067800,47696980067800,26051964600200,26051964600200,6210977229936.221,1514179755168.0798,1514179755168.0798])
Nphop=1/(np.exp(sc.h*frequences/(sc.k*temp[0]))-1)
l_freq=len(frequences)
nop_2=np.maximum(np.rint(frequences * sc.h / dE / sc.e),1)
# Variables of the photon interaction calculations

polarization_direction=np.array([1,0,0])
ph_freq=1e15
en_photon=sc.h*ph_freq/sc.e

nop=int(np.rint(en_photon/(dE)))
J=1e12
V=1e-3
n_photon=J/sc.e*V/sc.c/en_photon
factor=(sc.hbar/sc.m_e)**2/(2*V*sc.epsilon_0*sc.c**2)/(en_photon)*sc.e
epsilon=1e-2

# %% [markdown]
# Some useful functions:

# %%
def ind2coord(N,ik,starty,startz,stepy,stepz,yhat,zhat):
    iky=ik//N
    ikz=ik%N
    return (starty+stepy*iky)*yhat+(startz+stepz*ikz)*zhat

ind2coord=np.vectorize(ind2coord)


def bedistribution(energy,temperature):
    return 1/(np.exp(energy/sc.k/temperature)-1)

bedistribution_vect=np.vectorize(bedistribution)
bedistribution_vect.excluded.add(1)

def k_dist(kpoints_list):

    k_dist=np.zeros(kpoints_list.shape[0])
    for ik in range(kpoints_list.shape[0]):
        k=ind2coord(N=nkz,ik=kpoints_list[ik],starty=-np.pi/Ly,startz=-np.pi/Lz,stepy=dky,stepz=dkz,yhat=yhat,zhat=zhat)
        k_dist[ik]=np.linalg.norm(k)
    return k_dist

def layer(list_coord,delta,L, nbr_of_layers,ns,alt):
    intermed=[]
    for j in range(nbr_of_layers):
        tresh=alt[j]
        
        intermedbis=[]
        for k in range(len(list_coord[0])):
            if abs(tresh-np.dot(zhat,list_coord[:,k]))<delta:
                intermedbis.append(int(k))
        intermed.append(np.array(intermedbis))
    return np.array([np.array(intermed)+i*len(list_coord[0]) for i in range(ns)])

# # Initialization of the GF calculation over the whole grid










def energies_top(kx,ky,layer_top):
    h00,h10 = wannierham.block_mat_def(kx=kx, ky=ky, kz=0.0,nb=nb,ns=ns_L,n_range=n_range,hr=hr,cell=cell)
    h00_top=h00[layer_top,:][:,layer_top]
    h10_top=h10[layer_top,:][:,layer_top]

    h=h00_top+np.exp(-1j*kx*Lx*ns_L)*h10_top+np.exp(1j*kx*Lx*ns_L)*h10_top.conj().T
    e=np.real(eig(h)[0])
    e=np.sort(e)
    return e


def energies_bot(kx,ky,layer_bot):
    h00,h10 = wannierham.block_mat_def(kx=kx, ky=ky, kz=0.0,nb=nb,ns=ns_L,n_range=n_range,hr=hr,cell=cell)
    
    h00_bot=h00[layer_bot,:][:,layer_bot]

    h10_bot=h10[layer_bot,:][:,layer_bot]

    h=h00_bot+np.exp(-1j*kx*Lx*ns_L)*h10_bot+np.exp(1j*kx*Lx*ns_L)*h10_bot.conj().T
    
    e=np.real(eig(h)[0])
    e=np.sort(e)
    return e

def energies_full(kx,ky):
    h00,h10 = wannierham.block_mat_def(kx=kx, ky=ky, kz=0.0,nb=nb,ns=ns_L,n_range=n_range,hr=hr,cell=cell)

    h00_top=h00
    h10_top=h10

    h=h00_top+np.exp(-1j*kx*Lx*ns_L)*h10_top+np.exp(1j*kx*Lx*ns_L)*h10_top.conj().T
    e=np.real(eig(h)[0])
    e=np.sort(e)
    return e




def energies_fct(path,wannier_center,delta,L, nbr_of_layers,ns_,nb_L,alt):
    layer_top=np.zeros(ns_L*nb_L,dtype='i')
    layer_bot=np.zeros(ns_L*nb_L,dtype='i')
    layer_list=layer(list_coord=wannier_center, delta=delta,L=L, nbr_of_layers=nbr_of_layers,ns=ns_,alt=alt)

    
    for k in range(ns_):
        layer_top[nb_L*k:nb_L*(k+1)]=layer_list[k][0]
        layer_bot[nb_L*k:nb_L*(k+1)]=layer_list[k][1]
    energy_top_list=np.zeros((len(path),ns_*nb_L))
    energy_bot_list=np.zeros((len(path),ns_*nb_L))
    energy_full_list=np.zeros((len(path),nb_L*ns_*nbr_of_layers))
    for k in range(len(path)):
        energy_top_list[k]=energies_top(path[k][0],path[k][1],layer_top)
        energy_bot_list[k]=energies_bot(path[k][0],path[k][1],layer_bot)
        energy_full_list[k]=energies_full(path[k][0],path[k][1])

    return energy_top_list,energy_bot_list,energy_full_list


def ham_block(path,wannier_center,delta,L, nbr_of_layers,ns_,nb_L,alt):
    layer_top=np.zeros(ns_L*nb_L,dtype='i')
    layer_bot=np.zeros(ns_L*nb_L,dtype='i')
    layer_list=layer(list_coord=wannier_center, delta=delta,L=L, nbr_of_layers=nbr_of_layers,ns=ns_,alt=alt)

    
    for k in range(ns_):
        layer_top[nb_L*k:nb_L*(k+1)]=layer_list[k][0]
        layer_bot[nb_L*k:nb_L*(k+1)]=layer_list[k][1]
    h00_top=np.zeros((len(path),len(layer_top),len(layer_top)),dtype='complex')
    h10_top=np.zeros((len(path),len(layer_top),len(layer_top)),dtype='complex')
    h00_bot=np.zeros((len(path),len(layer_bot),len(layer_bot)),dtype='complex')
    h10_bot=np.zeros((len(path),len(layer_bot),len(layer_bot)),dtype='complex')
    h00_full=np.zeros((len(path),nb_L*nbr_of_layers*ns_,nb_L*nbr_of_layers*ns_),dtype='complex')
    h10_full=np.zeros((len(path),nb_L*nbr_of_layers*ns_,nb_L*nbr_of_layers*ns_),dtype='complex')
    for j in range(len(path)):
        kx=path[j][0]
        ky=path[j][1]
        h00,h10 = wannierham.block_mat_def(kx=kx, ky=ky, kz=0.,nb=nb,ns=ns_L,n_range=n_range,hr=hr,cell=cell)
        h00_top[j]=h00[layer_top,:][:,layer_top]
        h10_top[j]=h10[layer_top,:][:,layer_top]
        h00_bot[j]=h00[layer_bot,:][:,layer_bot]
        h10_bot[j]=h10[layer_bot,:][:,layer_bot]
        h00_full[j]=h00
        h10_full[j]=h10
    return h00_top,h10_top,h00_bot,h10_bot,h00_full,h10_full





band_struct=np.zeros((nb_L*ns_L,nkx))



dkx=2*np.pi/nkx/Lx

path=np.zeros((nkx*nky,2))
dky=2*np.pi/Ly/nky

for iky in range(nky):
    for ikx in range(nkx):
        ik = ikx + iky*nkx
        kx=-np.pi/Lx + dkx*ikx
        ky=-np.pi/Ly + dky*iky
        path[ik] = [kx,ky]

h00_top,h10_top,h00_bot,h10_bot,h00_full,h10_full=ham_block(path=path,wannier_center=wannier_center,delta=delta_alt,L=L,nbr_of_layers= 2,ns_=ns_L,nb_L=nb_L,alt=alt)

# for iky in range(nky):
#     for ikx in range(nkx):
#         ik = ikx + iky*nkx
        
#         kx=path[ik][0]
        

#         h=h00_top[ik,:,:]+np.exp(-1j*kx*Lx*ns_L)*h10_top[ik,:,:]+np.exp(1j*kx*Lx*ns_L)*h10_top[ik,:,:].conj().T
#         e=np.real(eig(h)[0])
#         e=np.sort(e)
#         band_struct[:,ikx]=e
#     fig=plt.figure()
#     for ind in range(nb_L*ns_L):
#         plt.plot(band_struct[ind])
#     plt.savefig(folder_save+'test_bs.pdf')
#     plt.show()
##Calculation out of plan
nb=16
nb_L=8
nx=10
ny=41
nz=41
ns_L=2

hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname=folder+hamiltonian,lreorder_axis=True,axis=[3,2,1],nb=nb,nx=nx,ny=ny,nz=nz)

nkz=nkx
nk=nky*nkz


ns=1
Lz=L[2]
Ly=L[1]
Lx=L[0]
dkz=2.0*np.pi/Lz / (nkz)
dky=2.0*np.pi/Ly / (nky)

alpha=cell[:,0]
beta=cell[:,1]
gamma=cell[:,2]
V=np.dot(np.cross(beta,gamma),alpha)
b_1=np.cross(beta,gamma)/V*2*np.pi
b_2=np.cross(gamma,alpha)/V*2*np.pi
b_3=np.cross(alpha,beta)/V*2*np.pi
xhat=alpha/np.linalg.norm(alpha)
yhat=-np.cross(xhat,gamma)
yhat=yhat/np.linalg.norm(yhat)
zhat=gamma/np.linalg.norm(gamma)


K=((b_2+b_3)/3)
Gamma=np.zeros(3)
Z=b_1/2
Q=(b_2+b_3)/3+b_1/2
# L=(b_1+b_3)/2
M=(b_3/2)

def make_dk_list(k_list,nstep):
    list=[]
    for k in range(len(k_list)-1):
        dk=(k_list[k+1]-k_list[k])/(nstep-1)
        list.append(dk)
    dk=(k_list[0]-k_list[-1])/(nstep-1)
    list.append(dk)
    return list

k_list=[K,Gamma,M]

dk_list=make_dk_list(k_list,nstep=nkz)
l=len(dk_list)







# %% [markdown]
# # Non coherent transport :

# %%
#Calculation of Pmn in real space for photon interaction :

pmn_r = wannierham.calc_momentum_operator(method='approx',nb=nb,nx=nx,ny=ny,nz=nz,hr=hr,cell=cell,n_range=n_range,wannier_center=wannier_center,rmn=np.zeros((3,nb,nb,nx,ny,nz)))

# %%

#Initialization of the initial array without putting to 0 the self-energies


ham = np.zeros((nb*length,nb*length,nk), dtype='complex')  
lead_h00 = np.zeros((nb*ns_L,nb*ns_L,2,nk), dtype='complex')
lead_h10 = np.zeros((nb*ns_L,nb*ns_L,2,nk), dtype='complex')
lead_coupling = np.zeros((nb*ns_L,nb*length,2,nk), dtype='complex')
M_mat= np.zeros((nb*length,nb*length,nk,1), dtype='complex')

cur=np.zeros((nk,nen,2),dtype='complex')
G_retarded=np.zeros((nk,nb*length,nb*length,nen),dtype='complex')
G_lesser=np.zeros((nk,nb*length,nb*length,nen),dtype='complex')
G_greater=np.zeros((nk,nb*length,nb*length,nen),dtype='complex')
te=np.zeros((nk,nen,2,2),dtype='complex')

#Initialization of the initial array of the self-energies if restart = True or not


# restart=True

# if restart == False:
#     sig_g=sig_g_tempo
#     sig_l=sig_l_tempo
#     siglead=siglead_tempo
#     sig_r=sig_r_tempo

# if restart==True:

sig_r = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
sig_l = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
sig_r_photon = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
sig_l_photon = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
sig_g_photon = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')

sig_g = np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
siglead=np.zeros((nb*ns_L,nb*ns_L,nen,2),dtype='complex')



nac=np.zeros((len(Ac_fact),nk),dtype='i')
nac = np.maximum(np.rint(Ac_fact[:, np.newaxis] * sc.h * k_dist(np.arange(nk)) / dE / sc.e),1)




for iky in range(nky):
        for ikz in range(nkz):
            ik = ikz + iky*nkz
            ky=-np.pi/Ly + dky*iky
            kz=-np.pi/Lz + dkz*ikz
            
            ham[:,:,ik] = wannierham.full_device_mat_def(ky=ky,kz=kz,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)
            # v[:,:,ik] = wannierham.full_device_bare_coulomb(ky=kz,kz=kz,length=length,eps=eps_screen,r0=r0,ldiag=True,nb=nb,ns=ns,
            #                                             method='pointlike',n_range=n_range,wannier_center=wannier_center,cell=cell)

            
            Pmn = wannierham.w90_momentum_full_device(ky=ky,kz=kz,length=length,ns=ns,n_range=n_range,nb=nb,cell=cell,pmn=pmn_r)

            M_mat[:,:,ik,0]=np.tensordot(polarization_direction, Pmn, axes=(0, -1))


            # h00,h10 = wannierham.block_mat_def(kx=0.0, ky=ky, kz=kz,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)

            lead_h10[0:nb*ns_L,0:nb*ns_L,0,ik] = np.transpose( np.conjugate(h10_full[ik,:,:]) ) 
            lead_h10[nb*ns_L-(nb*ns_L):nb*ns_L,nb*ns_L-(nb*ns_L):nb*ns_L,1,ik] = h10_full[ik,:,:]
            lead_h00[0:(nb*ns_L),0:(nb*ns_L),0,ik] = h00_full[ik,:,:]
            lead_h00[nb*ns_L-(nb*ns_L):nb*ns_L,nb*ns_L-(nb*ns_L):nb*ns_L,1,ik] = h00_full[ik,:,:]
            lead_coupling[0:nb*ns_L,0:nb,0,ik] = lead_h10[:,0:nb,0,ik]
            lead_coupling[0:nb*ns_L,nb*(length-ns):nb*length,1,ik] = lead_h10[:,0:nb,1,ik]




delta=1
niter=0


# %%

print('end Ham loop')
while delta> epsilon and niter<30:
    for iky in range(nky):
        for ikz in range(nkz):
            ik = ikz + iky*nkz
            ky=-np.pi/Ly + dky*iky
            kz=-np.pi/Lz + dkz*ikz

            
            G_retarded[ik],G_lesser[ik],G_greater[ik],cur[ik],te[ik] = gw_dense.calc_gf(ne=nen,e=energies,num_lead=2,nm_dev=nb*length,nm_lead=dim_lead,max_nm_lead=nb*ns_L,
                                                            ham=ham[:,:,ik],lead_h00=lead_h00[:,:,:,ik],lead_h10=lead_h10[:,:,:,ik],
                                                            siglead=siglead,t=lead_coupling[:,:,:,ik],
                                                            scat_sig_retarded=sig_r[:,:,:,ik],scat_sig_lesser=sig_l[:,:,:,ik],scat_sig_greater=sig_g[:,:,:,ik],
                                                            mu=mu,temp=temp,flatband=False)
    
    print('end_GF')
    #calculation of the photon SE:

    if interaction_photon==False and interaction_phonon==False:
        break
    if interaction_photon:
        sig_l_photon,sig_g_photon=gw_dense.selfenergy_eph_mono(nm=nb*length,nen=nen,en=energies,nop=nop,nky=nky,nkz=nkz,nqy=1,nqz=1,ik_start=1,ik_end=nky*nkz,iq_in=1,m=M_mat,g_lesser=np.transpose(G_lesser,(1,2,3,0)),g_greater=np.transpose(G_greater,(1,2,3,0)),
                                            n_bose=n_photon,gamma_q=True)
    
        sig_l_photon,sig_g_photon=factor*sig_l_photon,factor*sig_g_photon
        sig_l=sig_l_photon
        sig_g=sig_g_photon
        print(np.max(np.abs(sig_l_photon.flat)))

                                            
        print('end photon SE calculation')
    
    

    sig_r_old=sig_r


    sum_inf=np.zeros((nb*ns,nb*ns,nen,l_freq),dtype='complex')
    sum_sup=np.zeros((nb*ns,nb*ns,nen,l_freq),dtype='complex')


    sig_l_int=np.zeros((nb*length,nb*length,nen,nk),dtype='complex')
    sig_g_int=np.zeros((nb*length,nb*length,nen,nk),dtype='complex')



    sum_g_less=np.zeros((nb*length,nb*length,nen,nk),dtype='complex')
    sum_g_great=np.zeros((nb*length,nb*length,nen,nk),dtype='complex')


    Number_acoust_mode=len(Dac2)

    if interaction_phonon:
    
        for iky in range(nky):
            for ikz in range(nkz):
                ik = iky * nkz + ikz
                for iqy in range(-nky//8, nky//8):
                    for iqz in range(-nkz//8, nkz//8):
                        new_iky = (iky + iqy) % nky
                        new_ikz = (ikz + iqz) % nkz
                        idx = new_iky * nkz + new_ikz
                        sum_g_less[...,ik] += G_lesser[idx]
                        sum_g_great[...,ik] += G_greater[idx]
                        for ind in range(Number_acoust_mode):
                            nE=int(nac[ind, (nky//2+iqy)*nkz+(nkz//2+iqz)])
                            n = bedistribution(sc.e * dE * nE, temp[0])   
                            sig_l_int[...,nE:,ik]+=Dac2[ind] *n*G_lesser[idx, ...][...,:-nE]
                            sig_l_int[...,:-nE,ik]+=Dac2[ind] *(n+1)*G_lesser[idx, ...][...,nE:]
                            sig_g_int[...,:-nE,ik]+=Dac2[ind] *n*G_greater[idx, ...][...,nE: ]
                            sig_g_int[...,nE:,ik]+=Dac2[ind] *(n+1)*G_greater[idx, ...][...,:-nE]
        print('end big loop')
        #calculation of optical phonons SE

        for f in range(l_freq):
            nE=int(nop_2[f])
            n_phop=Nphop[f]

            sig_l_int[...,nE:,:]+=Dop2*(n_phop*sum_g_less[...,:-nE,:])
            sig_l_int[...,:-nE,:]+=Dop2*(n_phop+1)*sum_g_less[...,nE:,:]

            sig_g_int[...,:-nE,:]+=Dop2*(n_phop*sum_g_great[...,nE:,:])
            sig_g_int[...,nE:,:]+=Dop2*((n_phop+1)*sum_g_great[...,:-nE,:])




    sig_l+=sig_l_int
    sig_g+=sig_g_int


    
    sig_r=-1j/2 *np.imag(sig_l-sig_g)
    siglead[...,0]=sig_r[:nb*ns,:nb*ns,:,0]
    siglead[...,1]=sig_r[(length-ns)*nb:,(length-ns)*nb:,:,-1]
    if niter>0:  
        delta=np.max((np.linalg.norm(sig_r-sig_r_old, axis=(0,1))/np.linalg.norm(sig_r_old, axis=(0,1))).flat)

        print('Iteration : ', niter)
        print('Delta : ', delta)
    else:
        print('Iteration : ', niter)
        print('Delta : ', 'first iteration')
    niter+=1
print('End Green-function loop')
# sig_g_tempo=sig_g
# sig_l_tempo=sig_l
# siglead_tempo=siglead
# sig_r_tempo=sig_r




######################################################################################################
#####################################################################################################
band_struct=np.zeros((nb*ns_L,nk,ns_L))

# for ibz in range(ns_L):
for iky in range(nky):
    for ikx in range(nkx):
        ik = ikx + iky*nkx
        
        kx=path[ik][0]        
        phiz=kx*Lz*ns_L
        h=lead_h00[:,:,0,ik]+np.exp(+1j*phiz)*lead_h10[:,:,0,ik] + np.exp(-1j*phiz)*lead_h10[:,:,0,ik].conj().T
        e=np.real(eig(h)[0])
        e=np.sort(e)
        band_struct[:,ik,0]=e


######################################################################################################
#####################################################################################################

fig=plt.figure()
x=np.arange(nk)
A=-np.sum(np.imag(np.diagonal(G_retarded,axis1=1,axis2=2)),2)
pcm = plt.pcolormesh(x,energies,A[:,:].T, 
                    vmax=10,
                   cmap='hot')
for ind in range(nb*ns_L):
    # for ibz in range(ns_L):
    plt.plot(band_struct[ind,:,0],'r-',alpha=0.5)
plt.colorbar(pcm)
plt.title('Density of state over the whole grid')
plt.xlabel('k points indices')
plt.ylabel('Energies')

plt.savefig(folder_save+'DOS-BZ.pdf')


######################################################################################################
#CUT BY CLOSEST NEIGHBOUR METHOD
#####################################################################################################

# %% [markdown]
# # Cut with closest neighbour method


nk_cut=101

k_list=[K,Gamma,M]
dk_list=make_dk_list(k_list,nstep=nk_cut)

k_cut_list=np.zeros((len(dk_list)*nk_cut,3))


for ind in range(len(dk_list)):
    start_k=k_list[ind]
    dk=dk_list[ind]
    k_cut_list[ind*nk_cut:(ind+1)*nk_cut]=(np.array([ start_k+i*dk for i in range(nk_cut)]))



k_points_coords=np.array([ind2coord(N=nkz,ik=ik,starty=-np.pi/Ly,startz=-np.pi/Lz,stepy=dky,stepz=dkz,yhat=yhat,zhat=zhat) for ik in np.arange(nk)])

nn = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(k_points_coords)


replength=1 #nbr of cells around the unit cell taken into account

intermediate_dist=np.zeros(2*(2*replength+1)**2)
intermediate_ind=np.zeros(2*(2*replength+1)**2)
ind_final=np.zeros(len(k_cut_list),dtype='i')
distances_final=np.zeros(len(k_cut_list))
for k_pt in range(len(k_cut_list)):
    a=0
    for i in np.arange(-replength,replength+1):
        for j in np.arange(-replength,replength+1):
            intermediate_dist[a], intermediate_ind[a] = nn.kneighbors([k_cut_list[k_pt]+i*b_2+j*b_3])
            intermediate_dist[a+1], intermediate_ind[a+1] = nn.kneighbors([-k_cut_list[k_pt]+i*b_2+j*b_3])
            a+=2
    f=np.argmin(intermediate_dist)
    ind_final[k_pt]=intermediate_ind.flat[f]
    distances_final[k_pt]=intermediate_dist.flat[f]


# %%
nk_cut=len(k_cut_list) #redefinition of nk_cut
cur_nn=cur[ind_final,...]

G_retarded_diag_nn=(np.diagonal(G_retarded,axis1=2,axis2=1))[ind_final,...]
G_lesser_diag_nn=(np.diagonal(G_lesser,axis1=2,axis2=1))[ind_final,...]
G_greater_diag_nn=(np.diagonal(G_greater,axis1=2,axis2=1))[ind_final,...]
R_in_nn=np.zeros((nk_cut,nb*length,nb*length,nen),dtype='complex')
R_out_nn=np.zeros((nk_cut,nb*length,nen),dtype='complex')
te_nn=te[ind_final,...]


# for j in range(len(k_cut_list)):
#     ik = int(ind_final[j])
#     k=ind2coord(N=nkz,ik=ik,starty=-np.pi/Ly,startz=-np.pi/Lz,stepy=dky,stepz=dkz,yhat=yhat,zhat=zhat)
    
#     ky=k[0]
#     kz=k[1]
    

    

#     G_retarded_nn,G_lesser_nn,G_greater_nn,cur_nn[j],te_nn[j] = gw_dense.calc_gf(ne=nen,e=energies,num_lead=2,nm_dev=nb*length,nm_lead=dim_lead,max_nm_lead=nb*ns,
#                                                     ham=ham[:,:,ik],lead_h00=lead_h00[:,:,:,ik],lead_h10=lead_h10[:,:,:,ik],
#                                                     siglead=siglead,t=lead_coupling[:,:,:,ik],
#                                                     scat_sig_retarded=sig_r[:,:,:,ik],scat_sig_lesser=sig_l[:,:,:,ik],scat_sig_greater=sig_g[:,:,:,ik],
#                                                     mu=mu,temp=temp,flatband=False)
#     G_retarded_diag_nn[j]=np.diagonal(G_retarded_nn,axis1=1,axis2=0).T
#     G_lesser_diag_nn[j]=np.diagonal(G_lesser_nn,axis1=1,axis2=0).T
#     G_greater_diag_nn[j]=np.diagonal(G_greater_nn,axis1=1,axis2=0).T

#     G_retarded_balist_nn,G_lesser_balist_nn,G_greater_balist_nn,cur_balist_nn[j],te_balist_nn[j] = gw_dense.calc_gf(ne=nen,e=energies,num_lead=2,nm_dev=nb*length,nm_lead=dim_lead,max_nm_lead=nb*ns,
#                                                     ham=ham[:,:,ik],lead_h00=lead_h00[:,:,:,ik],lead_h10=lead_h10[:,:,:,ik],
#                                                     siglead=siglead,t=lead_coupling[:,:,:,ik],
#                                                     scat_sig_retarded=np.zeros((nb*length,nb*length,nen)),scat_sig_lesser=np.zeros((nb*length,nb*length,nen)),scat_sig_greater=np.zeros((nb*length,nb*length,nen)),
#                                                     mu=mu,temp=temp,flatband=False)
#     G_retarded_diag_balist_nn[j]=np.diagonal(G_retarded_balist_nn,axis1=1,axis2=0).T
#     G_lesser_diag_balist_nn[j]=np.diagonal(G_lesser_balist_nn,axis1=1,axis2=0).T
#     G_greater_diag_balist_nn[j]=np.diagonal(G_greater_balist_nn,axis1=1,axis2=0).T

## can also be obtain by taking the indices of final_ind (much faster)

# %%
sig_l_nn=sig_l[...,ind_final]


# %%


######################################################################################################
#####################################################################################################

fig=plt.figure()

x=np.arange(nk_cut)
A=np.log(np.sum(np.imag(np.diagonal(sig_l_nn,axis1=0,axis2=1)),-1).T)

for i in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][i],'--',color='white')
pcm = plt.pcolormesh(x,energies,A[:,:].T, 
                    # vmin=-12, vmax=1,
                   cmap='hot')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.savefig(folder_save+'sig_l_nn.pdf')



# %%
######################################################################################################
#CUT BY INTERPOLATION
#####################################################################################################

nk_cut=101



# dk_list=[(Gamma-K)/(nkz-1),(Z-Gamma)/(nkz-1),(Q-Z)/(nkz-1),(L-Q)/(nkz-1),(Gamma-L)/(nkz-1),(M-gamma)/(nkz-1),(K-M)/(nkz-1)]
dk_list=make_dk_list(k_list,nstep=nk_cut)

k_cut_list=np.zeros((len(dk_list)*nk_cut,3))


for ind in range(len(dk_list)):
    start_k=k_list[ind]
    dk=dk_list[ind]
    k_cut_list[ind*nk_cut:(ind+1)*nk_cut]=(np.array([ start_k+i*dk for i in range(nk_cut)]))



nk_cut=len(k_cut_list)

k_points_coords=np.array([ind2coord(N=nkz,ik=ik,starty=-np.pi/Ly,startz=-np.pi/Lz,stepy=dky,stepz=dkz,yhat=yhat,zhat=zhat) for ik in np.arange(nk)])

K_Y=-3*np.pi/Ly + dky*np.arange(3*nky)
K_Z=-3*np.pi/Lz + dkz*np.arange(3*nkz)
K_list_large=K_Y[:,np.newaxis]*yhat+K_Z[:,np.newaxis]*zhat
G_lesser_large=np.zeros((3*nky,3*nkz,nen))

G_retarded_diag=np.sum(np.imag(np.diagonal(G_retarded,axis1=1,axis2=2)),2).reshape(nky,nkz,nen)
G_lesser_diag=np.sum(np.imag(np.diagonal(G_lesser,axis1=1,axis2=2)),2).reshape(nky,nkz,nen)
sig_l_calc=(np.sum(np.imag(np.diagonal(sig_l,axis1=0,axis2=1)),-1).T).reshape(nky,nkz,nen)

sig_l_photon_calc=(np.sum(np.imag(np.diagonal(sig_l_photon,axis1=0,axis2=1)),-1).T).reshape(nky,nkz,nen)
cur_calc=np.real(cur[:,:,0].reshape(nky,nkz,nen))


def larger(f,n_rep,nky,nkz):
    if len(f.shape)>2:
        large=np.zeros((3*nky,3*nkz,*f.shape[2:]))
    else:
        large=np.zeros((3*nky,3*nkz))
    for i in range(n_rep):
        for j in range(n_rep):
            large[i*nky:(i+1)*nky,j*nkz:(j+1)*nkz]=f
    return large


G_retarded_large=larger(G_retarded_diag,3,nky,nkz)

G_lesser_large=larger(G_lesser_diag,3,nky,nkz)
cur_large=larger(cur_calc,3,nky,nkz)
sig_l_calc_large=larger(sig_l_calc,3,nky,nkz)
sig_l_photon_calc_large=larger(sig_l_photon_calc,3,nky,nkz)
G_retarded_calc=np.zeros((nk_cut,nen),dtype='complex')
G_lesser_calc=np.zeros((nk_cut,nen),dtype='complex')
sig_l_interp_calc=np.zeros((nk_cut,nen),dtype='complex')
sig_l_photon_interp_calc=np.zeros((nk_cut,nen),dtype='complex')
cur_calc_2=np.zeros((nk_cut,nen),dtype='complex')

for ie in range(nen):

    G_lesser_ie=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],G_lesser_large[:,:,ie])
    G_retarded_ie=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],G_retarded_large[:,:,ie])
    sig_l_ie=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],sig_l_calc_large[:,:,ie])
    sig_l_photon_ie=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],sig_l_photon_calc_large[:,:,ie])
    cur_ie=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],cur_large[:,:,ie])
    for k in range(nk_cut):
        G_lesser_calc[k,ie]=G_lesser_ie(k_cut_list[k,0],k_cut_list[k,1])
        cur_calc_2[k,ie]=cur_ie(k_cut_list[k,0],k_cut_list[k,1])
        G_retarded_calc[k,ie]=G_retarded_ie(k_cut_list[k,0],k_cut_list[k,1])
        sig_l_interp_calc[k,ie]=sig_l_ie(k_cut_list[k,0],k_cut_list[k,1])
        sig_l_photon_interp_calc[k,ie]=sig_l_photon_ie(k_cut_list[k,0],k_cut_list[k,1])









# %%
def interp(f,cut):
    large=larger(f=f.reshape(nky,nkz,nen),n_rep=3,nky=nky,nkz=nkz)
    int=np.zeros((len(cut),nen),dtype='complex')

    for ie in range(nen):
        int_ie=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],large[:,:,ie])
        for k in range(len(cut)):
            int[k,ie]=int_ie(cut[k,0],cut[k,1])
    return int

# %%

band_struct_interpolate=np.zeros((nb*ns_L,len(k_cut_list)))
for k in range(nb*ns_L):
    band_struct_larger_k=larger(np.reshape(band_struct[k],(nky,nkz)),3,nky,nkz)
    band_struct_k=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],band_struct_larger_k)
    for ik in range(len(k_cut_list)):
        
        band_struct_interpolate[k,ik]=band_struct_k(k_cut_list[ik,0],k_cut_list[ik,1])

# %%
######################################################################################################
#####################################################################################################
fig=plt.figure()
for k in range(nb*ns_L):
    plt.plot(band_struct_interpolate[k])

#plt.show()

# %%
######################################################################################################
#####################################################################################################
fig=plt.figure()
pcm = plt.pcolormesh(np.arange(nk_cut),energies,np.log(-np.real(G_retarded_calc).T), 
                    # vmin=-10, vmax=4,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.ylim(0, 14)
plt.savefig(folder_save+'Greaterdedinterp014.pdf')


######################################################################################################
#####################################################################################################

fig=plt.figure()
pcm = plt.pcolormesh(np.arange(nk_cut),energies,np.real(cur_calc_2.T), 
                    # vmin=-1e-4, vmax=4,
                   cmap=cm.seismic)

for k in range(nb*ns_L):
    plt.plot(band_struct_interpolate[k],'--',color='black',alpha=0.1)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.ylim(0, 20)
plt.savefig(folder_save+'current.pdf')


######################################################################################################
#####################################################################################################

fig=plt.figure()
pcm = plt.pcolormesh(np.arange(nk_cut),energies,np.log(-np.real(G_retarded_calc).T), 
                    vmin=-10, vmax=4,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.ylim(0, 20)
plt.savefig(folder_save+'Gretardeddedinterp.pdf')
# %%

######################################################################################################
#####################################################################################################

fig=plt.figure()
pcm = plt.pcolormesh(np.arange(nk_cut),energies,np.log(np.real(G_lesser_calc).T), 
                    vmin=-10,vmax=4,
                   cmap='hot')

plt.axhline(7.7,color='green',linestyle='--')
plt.axhline(8.2,color='green',linestyle='--')
for k in range(nb*ns_L):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.axhline(EF+en_photon/2, color='white', linestyle='--')
plt.ylim(0,14)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')

plt.colorbar(pcm)
plt.savefig(folder_save+'glesserinterp.pdf')
# %%


######################################################################################################
#####################################################################################################

fig=plt.figure()
pcm = plt.pcolormesh(np.arange(nk_cut),energies,np.log(np.real(sig_l_interp_calc).T), 
                    vmin=-12,vmax=3,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.savefig(folder_save+'sigl-interp.pdf')
# %%

######################################################################################################
#####################################################################################################


fig=plt.figure()
pcm = plt.pcolormesh(np.arange(nk_cut),energies,np.log(np.real(sig_l_photon_interp_calc).T), 
                    vmin=-12,vmax=-4,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.savefig(folder_save+'sigl_photon-interp.pdf')
plt.colorbar(pcm)


######################################################################################################
#calculation of the current
#####################################################################################################


nEF=int((EF-emin)/dE)
I_elec=np.zeros((nb*length,nb*length))
Jhole=np.zeros((nb*length,nb*length))

def current_int(ham,i,j,G_lesser,nk,nEF):
    I=0
    J=0
    for ik in range(nk):
        for ie in range(nen):
            if ie>=nEF:
                I+=(ham[i,j,ik]*G_lesser[ik,j,i,ie]-ham[j,i,ik]*G_lesser[ik,i,j,ie])*(energies[ie]-EF)
            else:
                J+=np.sum(ham[i,j,ik]*G_lesser[ik,j,i,ie]-ham[j,i,ik]*G_lesser[ik,i,j,ie])*(energies[ie]-EF)
        
    return I,J

for i in range (length*nb):
    for j in range(length*nb):
    
        I_elec[i,j],Jhole[i,j]=current_int(ham,i,j,G_lesser,nky*nkz,nEF)

ind=layer(list_coord=wannier_center,delta=delta_alt,L=L, nbr_of_layers=2,ns=length,alt=alt)

energy_current_electron=np.zeros((length*2))
energy_current_hole=np.zeros((length*2))
for i in range(length):
    for j in range(2):
        if j==0:
            energy_current_electron[j+i*2]+=np.sum(I_elec[ind[i,j],:][:,ind[i,j+1]])
            energy_current_hole[j+i*2]+=np.sum(Jhole[ind[i,j],:][:,ind[i,j+1]])
        else:
            for l in range(i+1,length):
                energy_current_electron[j+i*2]+=np.sum(I_elec[ind[i,j],:][:,ind[l,j]])
                energy_current_hole[j+i*2]+=np.sum(Jhole[ind[i,j],:][:,ind[l,j]])

fig=plt.figure()
x=np.arange(2*length)
plt.plot(x,energy_current_electron)
plt.title('energy current electron')
plt.xlabel('layer_number')
plt.ylabel('Energy current')
plt.savefig(folder_save+'energy current electron.pdf')
fig=plt.figure()
plt.plot(x,energy_current_hole)
plt.title('energy current hole')

plt.xlabel('layer_number')
plt.ylabel('Energy current')
plt.savefig(folder_save+'energy current electron.pdf')


######################################################################################################
#CALCULATION OF THE SCATTERING
#####################################################################################################

# %% [markdown]
# ### R_in and R_out

# %%

R_out=-np.trace(np.matmul(np.transpose(np.imag(sig_g),(3,2,0,1)),np.transpose(np.imag(G_lesser),(0,3,1,2))),axis1=-2,axis2=-1)
R_in=-np.trace(np.matmul(np.transpose(np.imag(sig_l),(3,2,0,1)),np.transpose(np.imag(G_greater),(0,3,1,2))),axis1=-2,axis2=-1)

# %%

R_out_phonon=-np.trace(np.matmul(np.transpose(np.imag(sig_g-sig_g_photon),(3,2,0,1)),np.transpose(np.imag(G_lesser),(0,3,1,2))),axis1=-2,axis2=-1)
R_in_phonon=-np.trace(np.matmul(np.transpose(np.imag(sig_l-sig_l_photon),(3,2,0,1)),np.transpose(np.imag(G_greater),(0,3,1,2))),axis1=-2,axis2=-1)

# %%

R_out_photon=-np.trace(np.matmul(np.transpose(np.imag(sig_g_photon),(3,2,0,1)),np.transpose(np.imag(G_lesser),(0,3,1,2))),axis1=-2,axis2=-1)
R_in_photon=-np.trace(np.matmul(np.transpose(np.imag(sig_l_photon),(3,2,0,1)),np.transpose(np.imag(G_greater),(0,3,1,2))),axis1=-2,axis2=-1)

# %%
def e_current_dens(R_in,R_out,k,E,EF):

    a=(R_in[k,E]-R_out[k,E])
    b=(energies[E]-EF)
    
    return (R_in[k,E]-R_out[k,E])*(energies[E]-EF)
e_current_dens=np.vectorize(e_current_dens, excluded=['R_in'])
e_current_dens.excluded.add(0)
e_current_dens.excluded.add(1)
e_current_dens.excluded.add(2)


# %%
e_current_dens_list=np.zeros((nk,nen))

for k in range(nk):
    e_current_dens_list[k]=((e_current_dens(R_in,R_out,k,np.arange(nen),EF=8.1)))


e_current_dens_list_photon=np.zeros((nk,nen))

for k in range(nk):
    e_current_dens_list_photon[k]=((e_current_dens(R_in_photon,R_out_photon,k,np.arange(nen),EF=8.1)))


e_current_dens_list_phonon=np.zeros((nk,nen))

for k in range(nk):
    e_current_dens_list_phonon[k]=((e_current_dens(R_in_phonon,R_out_phonon,k,np.arange(nen),EF=8.1)))

# %%
e_current_dens_list_cut_interp=interp(e_current_dens_list,k_cut_list)
e_current_dens_list_cut_interp_photon=interp(e_current_dens_list_photon,k_cut_list)
e_current_dens_list_cut_interp_phonon=interp(e_current_dens_list_phonon,k_cut_list)




# %%
######################################################################################################
#####################################################################################################
fig=plt.figure()

x=np.arange(nk_cut)




pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(e_current_dens_list_cut_interp).T,1e-20)), 
                     vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')

plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.ylim(4,13)


plt.savefig(folder_save+'e-current-interp.pdf')
######################################################################################################
#####################################################################################################
fig=plt.figure()
x=np.arange(nk_cut)



pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(e_current_dens_list_cut_interp_photon).T,1e-20)), 
                     vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')

plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.ylim(4,13)

plt.savefig(folder_save+'e-current-interp-photon.pdf')


x=np.arange(nk_cut)
######################################################################################################
#####################################################################################################

fig=plt.figure()
pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(e_current_dens_list_cut_interp_phonon).T,1e-20)), 
                     vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')

plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.ylim(4,13)

plt.savefig(folder_save+'e-current-interp-phonon.pdf')

######################################################################################################
#####################################################################################################

# %%
R_in_cut=interp(R_in,k_cut_list)
R_out_cut=interp(R_out,k_cut_list)
R_in_phonon_cut=interp(R_in_phonon,k_cut_list)
R_out_phonon_cut=interp(R_out_phonon,k_cut_list)
R_in_photon_cut=interp(R_in_photon,k_cut_list)
R_out_photon_cut=interp(R_out_photon,k_cut_list)

# %%
######################################################################################################
#####################################################################################################
fig=plt.figure()

x=np.arange(nk_cut)
pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_in_cut),1e-20)).T, 
                    vmin=-15,vmax=-3,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)

plt.colorbar(pcm)
plt.ylim(0,14)
plt.axhline(EF,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.savefig(folder_save+'Rin-cut.pdf')
#plt.show()
######################################################################################################
#####################################################################################################
fig=plt.figure()

x=np.arange(nk_cut)
pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_out_cut),1e-20)).T, 
                    vmin=-15,vmax=-3,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)

plt.colorbar(pcm)
plt.ylim(0,14)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'Routcut.pdf')
#plt.show()


# %%
######################################################################################################
#####################################################################################################
fig=plt.figure()

x=np.arange(nk_cut)
pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_in_photon_cut),1e-20)).T, 
                    vmin=-15,vmax=-4,
                   cmap='hot')
plt.colorbar(pcm)
plt.ylim(0,14)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
#plt.gca().set_aspect("equal")
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)

plt.savefig(folder_save+'rinphoton.pdf')
#plt.show()
######################################################################################################
#####################################################################################################
fig=plt.figure()

x=np.arange(nk_cut)
pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_out_photon_cut),1e-20)).T, 
                    vmin=-15,vmax=-4,
                   cmap='hot')
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.colorbar(pcm)
plt.ylim(0,14)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')#plt.gca().set_aspect("equal")
for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)

plt.savefig(folder_save+'routphoton.pdf')
#plt.show()

######################################################################################################
#####################################################################################################
# %%

fig=plt.figure()

x=np.arange(nk_cut)

pcm = plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_out_cut),1e-20)).T, 
                    vmin=-10,vmax=0,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,14)

plt.savefig(folder_save+'routcut2.pdf')
#plt.show()

######################################################################################################
#####################################################################################################
# %%
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut=(e_current_dens_list)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)

plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,20)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'ecurrentcutnn.pdf')

######################################################################################################
#####################################################################################################
# %%
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut=(e_current_dens_list)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)

plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,14)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'ecurrentcutnn014.pdf')
#plt.show()

######################################################################################################
#####################################################################################################
# %%
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_phonon=(e_current_dens_list_phonon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_phonon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)

plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,20)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'ecurrentphononnn.pdf')
#plt.show()
# %%
######################################################################################################
#####################################################################################################
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_phonon=(e_current_dens_list_phonon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_phonon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)

plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,14)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'ecurrentphononnn014.pdf')

# %%

######################################################################################################
#####################################################################################################

fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_photon=(e_current_dens_list_photon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_photon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)
plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,20)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'ecurrentphotonnn.pdf')
#plt.show()

######################################################################################################
#####################################################################################################

fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_photon=(e_current_dens_list_photon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_photon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)
plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,14)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'ecurrentphotonnn014.pdf')

# %%
######################################################################################################
#####################################################################################################
fig=plt.figure()
x=np.arange(nk_cut)
A=np.log(np.maximum(np.abs(np.real(R_in_phonon)[ind_final]),1e-20))
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-12,vmax=-4,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)
plt.colorbar(pcm)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.ylim(0,14)
plt.savefig(folder_save+'rinphonon-nn.pdf')
#plt.gca().set_aspect("equal")
#plt.show()
fig=plt.figure()
A=np.log(np.maximum(np.abs(np.real(R_out_phonon[ind_final])),1e-20))
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-12,vmax=-4,
                   cmap='hot')
for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)
plt.axhline(7.7,color='red',linestyle='--')
plt.axhline(8.2,color='red',linestyle='--')
plt.axhline(EF-en_photon/2,color='white',linestyle='--')
plt.axhline(EF+en_photon/2,color='white',linestyle='--')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)

plt.ylim(0,14)

#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'routphononnn.pdf')

#plt.show()

# %%
######################################################################################################
#####################################################################################################
fig =plt.figure()
x=np.arange(nk_cut)
A=np.log(np.maximum(np.abs(np.real(R_in_photon)[ind_final]),1e-20))
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-12,vmax=-4,
                   cmap='hot')

for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)
plt.colorbar(pcm)

plt.ylim(0,14)
plt.savefig(folder_save+'rinphotonnn.pdf')
#plt.gca().set_aspect("equal")
#plt.show()
fig=plt.figure()
x=np.arange(nk_cut)
A=np.log(np.maximum(np.abs(np.real(R_out_photon)[ind_final]),1e-20))
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm = plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-15,vmax=-4,
                   cmap='hot')

plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
for k in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)
plt.ylim(0,14)
plt.savefig(folder_save+'routphotonnn.pdf')
#plt.gca().set_aspect("equal")
#plt.show()



######################################################################################################
#####################################################################################################


fig, ax = plt.subplots(1, 2, figsize=(15, 7))

x = np.arange(nk_cut)
A = -np.sum(np.imag(G_retarded_diag_nn), axis=-1)

for i in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][i],'--',color='white')

pcm = plt.pcolormesh(x, energies, A.T, vmin=0, vmax=10, cmap='hot')
fig.colorbar(pcm)
plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies')
plt.ylim(0, 20)

# plt.ylim(0,20)
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'Band-structure-cut-closest.pdf')
#plt.show()
######################################################################################################
#####################################################################################################


# %%

fig=plt.figure()

x = np.arange(nk_cut)
A = np.log(np.sum(np.imag(G_lesser_diag_nn), axis=-1))

pcm = plt.pcolormesh(x, energies, A.T, vmin=-10, vmax=4, cmap='hot')
fig.colorbar(pcm)
for i in range(nb*ns_L):
    plt.plot(x,band_struct[:,ind_final][i],'--',color='white')
# ax[0].set_title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies')
plt.ylim(0, 14)

# plt.ylim(0,20)
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'Delec-cut.pdf')

######################################################################################################
#####################################################################################################


fig=plt.figure()
x = np.arange(nk_cut)
A = np.sum(np.imag(G_lesser_diag_nn), axis=-1)

pcm = plt.pcolormesh(x, energies, A.T, vmin=0, vmax=2, cmap='hot')
fig.colorbar(pcm)
# plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies (eV)')
plt.ylim(0, 14)
plt.savefig(folder_save+'Delec-nn.pdf')
#plt.show()
######################################################################################################
#####################################################################################################
fig=plt.figure()
x = np.arange(nk_cut)
A = np.log(np.sum(np.imag(G_greater_diag_nn), axis=-1))

pcm = plt.pcolormesh(x, energies, A.T, vmin=-12, vmax=2, cmap='hot')
fig.colorbar(pcm)
# plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies (eV)')
plt.ylim(0, 20)
plt.savefig(folder_save+'Dos-nn.pdf')
######################################################################################################
#####################################################################################################

fig= plt.figure()
fig=plt.figure()
x = np.arange(nk_cut)
A = np.sum(np.imag(G_greater_diag_nn), axis=-1)

pcm = plt.pcolormesh(x, energies, np.log(A.T), vmin=-12, vmax=2, cmap='hot')
fig.colorbar(pcm)
# plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies (eV)')
plt.ylim(0, 14)
plt.savefig(folder_save+'Dos014-nn.pdf')


######################################################################################################
#####################################################################################################










######################################################################################################
#CURRENT CALCULATION
#####################################################################################################


I=np.zeros((nb*length,nb*length))
def current(ham,i,j,G_lesser,nk):
    I=0
    for ik in range(nk):
        I+=np.sum(ham[i,j,ik]*G_lesser[ik,j,i,:]-ham[j,i,ik]*G_lesser[ik,i,j,:])
        
    return I

for i in range (nb*length):
    for j in range(nb*length):
        I[i,j]=current(ham,i,j,G_lesser,nky*nkz)







# %%
nEF=int((EF-emin)/dE)
I_=np.zeros((nb*length,nb*length))
J=np.zeros((nb*length,nb*length))
def current_int(ham,i,j,G_lesser,nk,nEF):
    I=0
    J=0
    for ik in range(nk):
        I+=np.sum(ham[i,j,ik]*G_lesser[ik,j,i,nEF:]-ham[j,i,ik]*G_lesser[ik,i,j,nEF:])
        J+=np.sum(ham[i,j,ik]*G_lesser[ik,j,i,:nEF]-ham[j,i,ik]*G_lesser[ik,i,j,:nEF])
        
    return I,J

for i in range (length*nb):
    for j in range(length*nb):
        I_[i,j],J[i,j]=current_int(ham,i,j,G_lesser,nky*nkz,nEF)

######################################################################################################
#SAVE OF VARIABLES
#####################################################################################################

filename = folder_save + 'G_lesser'
fileObject = open(filename, 'wb')


pickle.dump(G_lesser,fileObject)
fileObject.close()


filename = folder_save + 'G_retarded'
fileObject = open(filename, 'wb')

pickle.dump(G_retarded,fileObject)
fileObject.close()


filename = folder_save + 'R_in'
fileObject = open(filename, 'wb')

pickle.dump(R_in,fileObject)
fileObject.close()

filename = folder_save + 'R_out'
fileObject = open(filename, 'wb')

pickle.dump(R_out,fileObject)
fileObject.close()

filename = folder_save + 'R_in_photon'
fileObject = open(filename, 'wb')

pickle.dump(R_in_photon,fileObject)
fileObject.close()

filename = folder_save + 'R_out_photon'
fileObject = open(filename, 'wb')

pickle.dump(R_out_photon,fileObject)
fileObject.close()


filename = folder_save + 'I'
fileObject = open(filename, 'wb')

pickle.dump(I,fileObject)
fileObject.close()

filename = folder_save + 'I_'
fileObject = open(filename, 'wb')

pickle.dump(I_,fileObject)
fileObject.close()

filename = folder_save + 'J'
fileObject = open(filename, 'wb')

pickle.dump(J,fileObject)
fileObject.close()
