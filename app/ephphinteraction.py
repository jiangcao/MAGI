# Script realized by Loris CROS during its semester project with Dr. Jiang CAO, during spring 2024. End the 3/06/2024. If questions contact the locros@student.ethz.ch
# This script permits to calculate quantum transport observable taking into account electron-phonon and electron-photon interactions.

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





#############################################################################
#Parameters to modify depending on the studied system (Graphene, hBn, ...)
#############################################################################
nb=16
nx=10
ny=41
nz=41
folder='/usr/scratch/mont-fort1/sem24f6/MAGI/f2py/QE_local-tests/'




material='graphen'

if material=='graphene':
    hamiltonian='ham_dat'
    alt=[0,1]#Altitude of different layers
    delta_alt=0.4#delta around those layer to find the Wannier centers

if material=='hetero':
    hamiltonian='gr_hbn_hr.dat'
    alt=[0,3]
    delta_alt=1.0   


hr,wannier_center,n_range,cell,L= wannierham.load_from_file(fname=folder+hamiltonian,lreorder_axis=True,axis=[3,2,1],nb=nb,nx=nx,ny=ny,nz=nz)
folder_save='/usr/scratch/mont-fort1/sem24f6/MAGI/f2py/QE_local-tests/Save_hetero4/'
label='hetero0'
# %%

interaction_photon=True
interaction_phonon=False


#############################################################################
#############################################################################

# Global variables for the GF calculation
ns= 2
length= 2
nen= 3000 # number of energy points
nky=12
nkz=12
nk=nky*nkz


emin=-15.0
emax=30.0
energies= np.linspace(emin,emax,nen)
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
# L=(b_1+b_3)/2
M=(b_3/2)



#parameters of the leads
dim_lead= np.ones(2)* nb*ns
temp=  np.ones(2)* 300.0
EF=8
DEF=0.0
mu= np.array( [EF+DEF, EF-DEF] )


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





## Calculation for the balistic case a long the path k_list

ham_balist= np.zeros((nb*length,nb*length,nk), dtype='complex')  
lead_h00_balist= np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
lead_h10_balist= np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
lead_coupling_balist= np.zeros((nb*ns,nb*length,2,nk), dtype='complex')

sig_r_balist= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
sig_l_balist= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
sig_g_balist= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
siglead_balist=np.zeros((nb*ns,nb*ns,nen,2),dtype='complex')


for j in range(l):
    for ikz in range(nkz):
        ik= ikz+j*nkz
        k=dk_list[j]*ikz+k_list[j]
        ky=k[1]
        kz=k[0]
        
        ham_balist[:,:,ik]= wannierham.full_device_mat_def(ky=ky,kz=kz,nb=nb,ns=1,length=length,hr=hr,cell=cell,n_range=n_range)

        h00_balist,h10_balist= wannierham.block_mat_def(kx=0.0, ky=ky, kz=kz,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)

        lead_h10_balist[:,:,0,ik]= np.transpose( np.conjugate(h10_balist) ) 
        lead_h10_balist[:,:,1,ik]= h10_balist
        lead_h00_balist[:,:,0,ik]= h00_balist
        lead_h00_balist[:,:,1,ik]= h00_balist
        lead_coupling_balist[0:nb*ns,0:nb*ns,0,ik]= lead_h10_balist[:,:,0,ik]
        lead_coupling_balist[0:nb*ns,nb*(length-ns):nb*length,1,ik]= lead_h10_balist[:,:,1,ik]
        


G_retarded_diag_balist=np.zeros((l*nkz,nb*length,nen),dtype='complex')
G_lesser_diag_balist=np.zeros((l*nkz,nb*length,nen),dtype='complex')
te_balist=np.zeros((l*nkz,nen,2,2),dtype='complex')
for j in range(l):
    for ikz in range(nkz):
        ik= ikz+j*nkz
        k=dk_list[j]*ikz+k_list[j]
        ky=k[1]
        kz=k[0]
        

        

        G_retarded_balist,G_lesser_balist,G_greater_balist,cur_balist,te_balist[ik]= gw_dense.calc_gf(ne=nen,e=energies,num_lead=2,nm_dev=nb*length,nm_lead=dim_lead,max_nm_lead=nb*ns,
                                                        ham=ham_balist[:,:,ik],lead_h00=lead_h00_balist[:,:,:,ik],lead_h10=lead_h10_balist[:,:,:,ik],
                                                        siglead=siglead_balist,t=lead_coupling_balist[:,:,:,ik],
                                                        scat_sig_retarded=sig_r_balist[:,:,:,ik],scat_sig_lesser=sig_l_balist[:,:,:,ik],scat_sig_greater=sig_g_balist[:,:,:,ik],
                                                        mu=mu,temp=temp,flatband=False)
        G_retarded_diag_balist[ik]=np.diagonal(G_retarded_balist,axis1=1,axis2=0).T
        G_lesser_diag_balist[ik]=np.diagonal(G_lesser_balist,axis1=1,axis2=0).T


##Calculation of teh band structure in the lead

fig=plt.figure()

band_struct=np.zeros((nb*ns,nkz*l))
kz=np.linspace(0,2*np.pi/Lz/3,nkz*l)
ky=np.linspace(1,np.pi/Lz,nky)
for i in range(nkz*l):
    h=lead_h00_balist[:,:,0,i]+lead_h10_balist[:,:,0,i]+lead_h10_balist[:,:,0,i].conj().T
    e=np.real(eig(h)[0])
    e=np.sort(e)
    band_struct[:,i]=e

for j in range(0,nb*ns):
    plt.plot(kz,band_struct[j,:])
plt.savefig(folder_save+'Band-structure-lead.pdf')



# %% [markdown]
# # Non coherent transport :

# %%
#Calculation of Pmn in real space for photon interaction :

pmn_r= wannierham.calc_momentum_operator(method='approx',nb=nb,nx=nx,ny=ny,nz=nz,hr=hr,cell=cell,n_range=n_range,wannier_center=wannier_center,rmn=np.zeros((3,nb,nb,nx,ny,nz)))

# %%

#Initialization of the initial array without putting to 0 the self-energies


ham= np.zeros((nb*length,nb*length,nk), dtype='complex')  
lead_h00= np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
lead_h10= np.zeros((nb*ns,nb*ns,2,nk), dtype='complex')
lead_coupling= np.zeros((nb*ns,nb*length,2,nk), dtype='complex')
M_mat= np.zeros((nb*length,nb*length,nk,1), dtype='complex')

cur=np.zeros((nk,nen,2),dtype='complex')
G_retarded=np.zeros((nk,nb*length,nb*length,nen),dtype='complex')
G_lesser=np.zeros((nk,nb*length,nb*length,nen),dtype='complex')
G_greater=np.zeros((nk,nb*length,nb*length,nen),dtype='complex')
te=np.zeros((nk,nen,2,2),dtype='complex')

#Initialization of the initial array of the self-energies if restart= True or not


restart=True

if restart== False:
    sig_g=sig_g_tempo
    sig_l=sig_l_tempo
    siglead=siglead_tempo
    sig_r=sig_r_tempo

if restart==True:

    sig_r= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
    sig_l= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
    sig_r_photon= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
    sig_l_photon= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
    sig_g_photon= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')

    sig_g= np.zeros((nb*length,nb*length,nen,nk), dtype='complex')
    siglead=np.zeros((nb*ns,nb*ns,nen,2),dtype='complex')



    nac=np.zeros((len(Ac_fact),nk),dtype='i')
    nac= np.maximum(np.rint(Ac_fact[:, np.newaxis] * sc.h * k_dist(np.arange(nk)) / dE / sc.e),1)


#initialization of the values in the hamiltonian

for iky in range(nky):
        for ikz in range(nkz):
            ik= ikz + iky*nkz
            ky=-np.pi/Ly + dky*iky
            kz=-np.pi/Lz + dkz*ikz

            
            
            ham[:,:,ik]= wannierham.full_device_mat_def(ky=ky,kz=kz,nb=nb,ns=ns,length=length,hr=hr,cell=cell,n_range=n_range)

            Pmn= wannierham.w90_momentum_full_device(ky=ky,kz=kz,length=length,ns=ns,n_range=n_range,nb=nb,cell=cell,pmn=pmn_r)

            M_mat[:,:,ik,0]=np.tensordot(polarization_direction, Pmn, axes=(0, -1))


            h00,h10= wannierham.block_mat_def(kx=0.0, ky=ky, kz=kz,nb=nb,ns=ns,n_range=n_range,hr=hr,cell=cell)

            lead_h10[:,:,0,ik]= np.transpose( np.conjugate(h10) ) 
            lead_h10[:,:,1,ik]= h10
            lead_h00[:,:,0,ik]= h00
            lead_h00[:,:,1,ik]= h00
            lead_coupling[0:nb*ns,0:nb*ns,0,ik]= lead_h10[:,:,0,ik]
            lead_coupling[0:nb*ns,nb*(length-ns):nb*length,1,ik]= lead_h10[:,:,1,ik]


delta=1
niter=0


# %%

print('end Ham loop')
while delta> epsilon and niter<30:
    for iky in range(nky):
        for ikz in range(nkz):
            ik= ikz + iky*nkz
            ky=-np.pi/Ly + dky*iky
            kz=-np.pi/Lz + dkz*ikz

            
            G_retarded[ik],G_lesser[ik],G_greater[ik],cur[ik],te[ik]= gw_dense.calc_gf(ne=nen,e=energies,num_lead=2,nm_dev=nb*length,nm_lead=dim_lead,max_nm_lead=nb*ns,
                                                            ham=ham[:,:,ik],lead_h00=lead_h00[:,:,:,ik],lead_h10=lead_h10[:,:,:,ik],
                                                            siglead=siglead,t=lead_coupling[:,:,:,ik],
                                                            scat_sig_retarded=sig_r[:,:,:,ik],scat_sig_lesser=sig_l[:,:,:,ik],scat_sig_greater=sig_g[:,:,:,ik],
                                                            mu=mu,temp=temp,flatband=False)
    
    print('end_GF')
    #calculation of the photon SE:
    if (interaction_photon==False and interaction_phonon==False):
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

sig_g_tempo=sig_g
sig_l_tempo=sig_l
siglead_tempo=siglead
sig_r_tempo=sig_r


# %%
band_struct=np.zeros((nb*ns,nk))

for i in range(nk):
    h=lead_h00[:,:,0,i]+lead_h10[:,:,0,i]+lead_h10[:,:,0,i].conj().T
    e=np.real(eig(h)[0])
    e=np.sort(e)
    band_struct[:,i]=e




fig=plt.figure()
x=np.arange(nk)
A=-np.sum(np.imag(np.diagonal(G_retarded,axis1=1,axis2=2)),2)
pcm= plt.pcolormesh(x,energies,A[:,:].T, 
                    vmax=10,
                   cmap='hot')

plt.colorbar(pcm)
plt.title('Density of state over the whole grid')
plt.xlabel('k points indices')
plt.ylabel('Energies')

plt.savefig(folder_save+'DOS-BZ.pdf')




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

nn= NearestNeighbors(n_neighbors=1, algorithm='auto').fit(k_points_coords)


replength=0 #nbr of cells around the unit cell taken into account

intermediate_dist=np.zeros(2*(2*replength+1)**2)
intermediate_ind=np.zeros(2*(2*replength+1)**2)
ind_final=np.zeros(len(k_cut_list),dtype='i')
distances_final=np.zeros(len(k_cut_list))
for k_pt in range(len(k_cut_list)):
    a=0
    for i in np.arange(-replength,replength+1):
        for j in np.arange(-replength,replength+1):
            intermediate_dist[a], intermediate_ind[a]= nn.kneighbors([k_cut_list[k_pt]+i*b_2+j*b_3])
            intermediate_dist[a+1], intermediate_ind[a+1]= nn.kneighbors([-k_cut_list[k_pt]+i*b_2+j*b_3])
            a+=2
    f=np.argmin(intermediate_dist)
    ind_final[k_pt]=intermediate_ind.flat[f]
    distances_final[k_pt]=intermediate_dist.flat[f]




nk_cut=len(k_cut_list)
cur_nn=np.zeros((nk_cut,nen,2),dtype='complex')

G_retarded_diag_nn=(np.diagonal(G_retarded,axis1=2,axis2=1))[ind_final,...]
G_lesser_diag_nn=(np.diagonal(G_lesser,axis1=2,axis2=1))[ind_final,...]
G_greater_diag_nn=(np.diagonal(G_greater,axis1=2,axis2=1))[ind_final,...]
R_in_nn=np.zeros((nk_cut,nb*length,nb*length,nen),dtype='complex')
R_out_nn=np.zeros((nk_cut,nb*length,nen),dtype='complex')
te_nn=te[ind_final,...]




# %%
sig_l_nn=sig_l[...,ind_final]


# %%




fig=plt.figure()

x=np.arange(nk_cut)
A=np.log(np.sum(np.imag(np.diagonal(sig_l_nn,axis1=0,axis2=1)),-1).T)

for i in range(nb*ns):
    plt.plot(x,band_struct[:,ind_final][i],'--',color='white')
pcm= plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-12, vmax=1,
                   cmap='hot')
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.savefig(folder_save+'sig_l_nn.pdf')



# %%


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

band_struct_interpolate=np.zeros((nb*ns,len(k_cut_list)))
for k in range(nb*ns):
    band_struct_larger_k=larger(np.reshape(band_struct[k],(nky,nkz)),3,nky,nkz)
    band_struct_k=interpolate.interp2d(K_list_large[:,0],K_list_large[:,1],band_struct_larger_k)
    for ik in range(len(k_cut_list)):
        
        band_struct_interpolate[k,ik]=band_struct_k(k_cut_list[ik,0],k_cut_list[ik,1])

# %%
fig=plt.figure()
for k in range(nb*ns):
    plt.plot(band_struct_interpolate[k])

#plt.show()

# %%
fig=plt.figure()
pcm= plt.pcolormesh(np.arange(nk_cut),energies,np.log(-np.real(G_retarded_calc).T), 
                    vmin=-10, vmax=4,
                   cmap='hot')

for k in range(nb*ns):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.ylim(0, 14)
plt.savefig(folder_save+'Greaterdedinterp014.pdf')

fig=plt.figure()
pcm= plt.pcolormesh(np.arange(nk_cut),energies,np.real(cur_calc_2.T), 
                    # vmin=-1e-4, vmax=4,
                   cmap=cm.seismic)

for k in range(nb*ns):
    plt.plot(band_struct_interpolate[k],'--',color='black',alpha=0.1)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.ylim(0, 20)
plt.savefig(folder_save+'current.pdf')


fig=plt.figure()
pcm= plt.pcolormesh(np.arange(nk_cut),energies,np.log(-np.real(G_retarded_calc).T), 
                    vmin=-10, vmax=4,
                   cmap='hot')

for k in range(nb*ns):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.ylim(0, 20)
plt.savefig(folder_save+'Gretardeddedinterp.pdf')
# %%
fig=plt.figure()
pcm= plt.pcolormesh(np.arange(nk_cut),energies,np.log(np.real(G_lesser_calc).T), 
                    vmin=-10,vmax=4,
                   cmap='hot')

plt.axhline(7.7,color='green',linestyle='--')
plt.axhline(8.2,color='green',linestyle='--')
for k in range(nb*ns):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.axhline(EF+en_photon/2, color='white', linestyle='--')
plt.ylim(0,14)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')

plt.colorbar(pcm)
plt.savefig(folder_save+'glesserinterp.pdf')
# %%
fig=plt.figure()
pcm= plt.pcolormesh(np.arange(nk_cut),energies,np.log(np.real(sig_l_interp_calc).T), 
                    vmin=-12,vmax=3,
                   cmap='hot')
for k in range(nb*ns):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
plt.savefig(folder_save+'sigl-interp.pdf')
# %%
fig=plt.figure()
pcm= plt.pcolormesh(np.arange(nk_cut),energies,np.log(np.real(sig_l_photon_interp_calc).T), 
                    vmin=-12,vmax=-4,
                   cmap='hot')
for k in range(nb*ns):
    plt.plot(band_struct_interpolate[k],'--',color='white',alpha=0.5)
plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.savefig(folder_save+'sigl_photon-interp.pdf')
plt.colorbar(pcm)



#calculation of the current



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
fig=plt.figure()

x=np.arange(nk_cut)




pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(e_current_dens_list_cut_interp).T,1e-20)), 
                     vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns):
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

fig=plt.figure()
x=np.arange(nk_cut)



pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(e_current_dens_list_cut_interp_photon).T,1e-20)), 
                     vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns):
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


fig=plt.figure()
pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(e_current_dens_list_cut_interp_phonon).T,1e-20)), 
                     vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns):
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



# %%
R_in_cut=interp(R_in,k_cut_list)
R_out_cut=interp(R_out,k_cut_list)
R_in_phonon_cut=interp(R_in_phonon,k_cut_list)
R_out_phonon_cut=interp(R_out_phonon,k_cut_list)
R_in_photon_cut=interp(R_in_photon,k_cut_list)
R_out_photon_cut=interp(R_out_photon,k_cut_list)

# %%

fig=plt.figure()

x=np.arange(nk_cut)
pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_in_cut),1e-20)).T, 
                    vmin=-15,vmax=-3,
                   cmap='hot')
for k in range(nb*ns):
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

fig=plt.figure()

x=np.arange(nk_cut)
pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_out_cut),1e-20)).T, 
                    vmin=-15,vmax=-3,
                   cmap='hot')
for k in range(nb*ns):
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

fig=plt.figure()

x=np.arange(nk_cut)
pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_in_photon_cut),1e-20)).T, 
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
for k in range(nb*ns):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)

plt.savefig(folder_save+'rinphoton.pdf')
#plt.show()

fig=plt.figure()

x=np.arange(nk_cut)
pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_out_photon_cut),1e-20)).T, 
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
for k in range(nb*ns):
    plt.plot(x,band_struct_interpolate[k],'--',color='white',alpha=0.5)

plt.savefig(folder_save+'routphoton.pdf')
#plt.show()


# %%

fig=plt.figure()

x=np.arange(nk_cut)

pcm= plt.pcolormesh(x,energies,np.log(np.maximum(np.real(R_out_cut),1e-20)).T, 
                    vmin=-10,vmax=0,
                   cmap='hot')
for k in range(nb*ns):
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


# %%
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut=(e_current_dens_list)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns):
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


# %%
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut=(e_current_dens_list)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns):
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


# %%
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_phonon=(e_current_dens_list_phonon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_phonon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns):
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
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_phonon=(e_current_dens_list_phonon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_phonon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')
for k in range(nb*ns):
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
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_photon=(e_current_dens_list_photon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_photon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns):
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
fig=plt.figure()

x=np.arange(nk_cut)
e_current_dens_list_cut_photon=(e_current_dens_list_photon)[ind_final]
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,np.log(np.maximum(e_current_dens_list_cut_photon.T,1e-20)), 
                    vmin=-15,vmax=-2,
                   cmap='hot')

for k in range(nb*ns):
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

fig=plt.figure()
x=np.arange(nk_cut)
A=np.log(np.maximum(np.abs(np.real(R_in_phonon)[ind_final]),1e-20))
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-12,vmax=-4,
                   cmap='hot')

for k in range(nb*ns):
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


pcm= plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-12,vmax=-4,
                   cmap='hot')
for k in range(nb*ns):
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

fig=plt.figure()
x=np.arange(nk_cut)
A=np.log(np.maximum(np.abs(np.real(R_in_photon)[ind_final]),1e-20))
#plt.contourf(energies,x,A[:,:,6], cmap='hot')


pcm= plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-12,vmax=-4,
                   cmap='hot')

for k in range(nb*ns):
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


pcm= plt.pcolormesh(x,energies,A[:,:].T, 
                    vmin=-15,vmax=-4,
                   cmap='hot')

plt.xlabel('Path indices')
plt.ylabel('Energies (eV)')
plt.colorbar(pcm)
for k in range(nb*ns):
    plt.plot(x,band_struct[:,ind_final][k],'--',color='white',alpha=0.5)
plt.ylim(0,14)
plt.savefig(folder_save+'routphotonnn.pdf')
#plt.gca().set_aspect("equal")
#plt.show()






fig= plt.figure()

x= np.arange(nk_cut)
A= -np.sum(np.imag(G_retarded_diag_nn), axis=-1)

for i in range(nb*ns):
    plt.plot(x,band_struct[:,ind_final][i],'--',color='white')

pcm= plt.pcolormesh(x, energies, A.T, vmin=0, vmax=10, cmap='hot')
fig.colorbar(pcm)
plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies')
plt.ylim(0, 20)


# plt.ylim(0,20)
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'Band-structure-cut-closest.pdf')
#plt.show()



# %%

fig= plt.figure()

x= np.arange(nk_cut)
A= np.log(np.sum(np.imag(G_lesser_diag_nn), axis=-1))

pcm= plt.pcolormesh(x, energies, A.T, vmin=-10, vmax=4, cmap='hot')
plt.colorbar(pcm)
for i in range(nb*ns):
    plt.plot(x,band_struct[:,ind_final][i],'--',color='white')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies')
plt.ylim(0, 14)

# plt.ylim(0,20)
#plt.gca().set_aspect("equal")
plt.savefig(folder_save+'Delec-cut.pdf')




fig=plt.figure()
x= np.arange(nk_cut)
A= np.sum(np.imag(G_lesser_diag_nn), axis=-1)

pcm= plt.pcolormesh(x, energies, np.log(A.T), vmin=-12, vmax=2, cmap='hot')
fig.colorbar(pcm)
# plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies (eV)')
plt.ylim(0, 14)
plt.savefig(folder_save+'Delec-nn.pdf')
#plt.show()

fig=plt.figure()
x= np.arange(nk_cut)
A= np.log(-np.sum(np.imag(G_greater_diag_nn), axis=-1))

pcm= plt.pcolormesh(x, energies, A.T, vmin=-12, vmax=2, cmap='hot')
fig.colorbar(pcm)
# plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies (eV)')
plt.ylim(0, 20)
plt.savefig(folder_save+'Dos-nn.pdf')
fig, ax= plt.subplots(1, 2, figsize=(15, 7))
fig=plt.figure()
x= np.arange(nk_cut)
A= -np.sum(np.imag(G_greater_diag_nn), axis=-1)

pcm= plt.pcolormesh(x, energies, np.log(A.T), vmin=-12, vmax=2, cmap='hot')
fig.colorbar(pcm)
# plt.title('Band structure obtained after a cut and used closest neighbour')
plt.xlabel('k points indices of the cut')
plt.ylabel('Energies (eV)')
plt.ylim(0, 14)
plt.savefig(folder_save+'Dos014-nn.pdf')

def set_axes_equal(ax):
    x_limits= ax.get_xlim3d()
    y_limits= ax.get_ylim3d()
    z_limits= ax.get_zlim3d()

    x_range= abs(x_limits[1] - x_limits[0])
    y_range= abs(y_limits[1] - y_limits[0])
    z_range= abs(z_limits[1] - z_limits[0])

    plot_radius= 0.5 * max([x_range, y_range, z_range])

    x_middle= np.mean(x_limits)
    y_middle= np.mean(y_limits)
    z_middle= np.mean(z_limits)

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


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


filename= folder_save + 'G_lesser'
fileObject= open(filename, 'wb')

pickle.dump(G_lesser,fileObject)
fileObject.close()


filename= folder_save + 'G_retarded'
fileObject= open(filename, 'wb')

pickle.dump(G_retarded,fileObject)
fileObject.close()


filename= folder_save + 'R_in'
fileObject= open(filename, 'wb')

pickle.dump(R_in,fileObject)
fileObject.close()

filename= folder_save + 'R_out'
fileObject= open(filename, 'wb')

pickle.dump(R_out,fileObject)
fileObject.close()

filename= folder_save + 'R_in_photon'
fileObject= open(filename, 'wb')

pickle.dump(R_in_photon,fileObject)
fileObject.close()

filename= folder_save + 'R_out_photon'
fileObject= open(filename, 'wb')

pickle.dump(R_out_photon,fileObject)
fileObject.close()


filename= folder_save + 'I'
fileObject= open(filename, 'wb')

pickle.dump(I,fileObject)
fileObject.close()

filename= folder_save + 'I_'
fileObject= open(filename, 'wb')

pickle.dump(I_,fileObject)
fileObject.close()

filename= folder_save + 'J'
fileObject= open(filename, 'wb')

pickle.dump(J,fileObject)
fileObject.close()
