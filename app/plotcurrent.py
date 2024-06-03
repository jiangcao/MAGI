# Script realized by Loris CROS during its semester project with Dr. Jiang CAO, during spring 2024. End the 3/06/2024. If questions contact the locros@student.ethz.ch
# This script permits to plot the current in real space calculated before by the script ephphinteraction.py

# %%
import numpy as np
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
import time
from scipy import interpolate
import ipywidgets as widgets
import pickle




# 

# %%

nb=16
nx=10
ny=41
nz=41
folder='/usr/scratch/mont-fort1/sem24f6/MAGI/f2py/QE_local-tests/'


material='graphene'

if material=='graphene':
    hamiltonian='ham_dat'
    alt=[0,1]
    delta_dist=0.4

if material=='hetero':
    hamiltonian='gr_hbn_hr.dat'
    alt=[0,3]
    delta_dist=1.0   


hr,wannier_center,n_range,cell,L = wannierham.load_from_file(fname=folder+hamiltonian,lreorder_axis=True,axis=[3,2,1],nb=nb,nx=nx,ny=ny,nz=nz)
folder_save='/usr/scratch/mont-fort1/sem24f6/MAGI/f2py/QE_local-tests/Save_gr3/'
label='gr0'
# %%



ns = 1
length = 2
nen = 3000 # number of energy points
nky=12
nkz=12
nk=nky*nkz


emin=-15.0
emax=30.0

Lz=L[2]
Ly=L[1]
Lx=L[0]
dkz=2.0*np.pi/Lz / (nkz)
dky=2.0*np.pi/Ly / (nky)
energies = np.linspace(emin,emax,nen)
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


dim_lead = np.ones(2)* nb*ns
temp =  np.ones(2)* 300.0
EF=8
DEF=0.5
mu = np.array( [EF+DEF, EF-DEF] )


# Variables of the phono and photon interaction calculations
Dop2=1e-4
Dac2=np.array([1e-4,1e-4,1e-4])
Ac_fact=np.array([84.9385484086987e10,1120.28324417723e10,649.5300760665194e10])


polarization_direction=np.array([1,0,0])
ph_freq=1e15
en_photon=sc.h*ph_freq/sc.e
dE=energies[1]-energies[0]
nop=int(np.rint(en_photon/(dE)))
J=1e12
V=1e-3
n_photon=J/sc.e*V/sc.c/en_photon

factor=(sc.hbar/sc.m_e)**2/(2*V*sc.epsilon_0*sc.c**2)/(en_photon)*sc.e



epsilon=1e-2




dim_lead = np.ones(2)* nb*ns
temp =  np.ones(2)* 300.0
EF=8
DEF=0.0
mu = np.array( [EF+DEF, EF-DEF] )


# Variables of the phonon and photon interaction calculations
Dop2=1e-4
Dac2=np.array([1e-4,1e-4,1e-4])
Ac_fact=np.array([84.9385484086987e10,1120.28324417723e10,649.5300760665194e10])


polarization_direction=np.array([1,0,0])
ph_freq=1e15
en_photon=sc.h*ph_freq/sc.e
dE=energies[1]-energies[0]
nop=int(np.rint(en_photon/(dE)))
J=1e12
V=1e-3
n_photon=J/sc.e*V/sc.c/en_photon

factor=(sc.hbar/sc.m_e)**2/(2*V*sc.epsilon_0*sc.c**2)/(en_photon)*sc.e



epsilon=1e-2

# %% [markdown]
#######################################################################################
# Some useful functions:
#######################################################################################
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


#######################################################################################
#######################################################################################



#Current matrix between Wannier centers
filename=folder_save+'I'
fileObject2 = open(filename, 'rb')
I = pickle.load(fileObject2)
fileObject2.close()



#Hole current between Wannier centers
filename=folder_save+'I_'
fileObject2 = open(filename, 'rb')
I_ = pickle.load(fileObject2)
fileObject2.close()

#electron current between Wannier centers
filename=folder_save+'J'
fileObject2 = open(filename, 'rb')
J = pickle.load(fileObject2)
fileObject2.close()












def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    x_middle = np.mean(x_limits)
    y_middle = np.mean(y_limits)
    z_middle = np.mean(z_limits)

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])




wannier_center2=np.zeros((3,length*nb))
for i in range(length):
    for j in range(nb):
        wannier_center2[:,j+i*nb]=wannier_center[:,j]+zhat*np.sqrt(np.dot(alpha,alpha))*i



ind=layer(list_coord=wannier_center,delta=delta_dist,L=L, nbr_of_layers=2,ns=length,alt=alt)
bot=wannier_center2[:,np.reshape(ind[:,0],-1)]
top=wannier_center2[:,np.reshape(ind[:,1],-1)]




# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])
mask=np.ones(nb*length,dtype=bool)
mask[np.reshape(ind[:,0],-1)[:nb//2]]= False
mask_opp=np.zeros(nb*length,dtype=bool)
mask_opp[np.reshape(ind[:,0],-1)[:nb//2]]= True
# I1=I/np.max(np.abs(I[:,mask][np.reshape(ind[:,0],-1)[:nb//2],:]))

# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(I1[:,mask][np.reshape(ind[:,0],-1)[:nb//2],:]), vmax=np.max(I1[:,mask][np.reshape(ind[:,0],-1)[:nb//2],:]))

# for i in np.reshape(ind[:,0],-1)[:nb//2]:
#     for j in np.arange(nb*length)[mask]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=I1[i,j]
#         color=cmap(norm(color_value))
#         intens=abs(I1[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)


# plt.title('Current layer1torest')

# set_axes_equal(ax)

# plt.show()




# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])
# mask2=np.ones(nb*length,dtype=bool)
# mask2[np.reshape(ind[:,1],-1)[:nb//2]]= False
# mask2_opp=np.zeros(nb*length,dtype=bool)
# mask2_opp[np.reshape(ind[:,1],-1)[:nb//2]]= True
# I2=I/np.max(np.abs(I[:,mask2][np.reshape(ind[:,1],-1)[:nb//2],:]))
# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(I2[:,mask2][mask2_opp,:]), vmax=np.max(I2[:,mask2][mask2_opp,:]))
# for i in np.reshape(ind[:,1],-1)[:nb//2]:
#     for j in np.arange(nb*length)[mask2]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=I2[i,j]
#         color=cmap(norm(color_value))

#         intens=abs(I2[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax)
# plt.title('Current layer2torest')

# plt.show()



# fig= plt.figure()

# x=np.arange(nb*length)

# #plt.contourf(energies,x,n[:,:,8], cmap='hot')
# pcm = plt.pcolormesh(x,x,np.real(I2).T, 
#                      vmin=-1, vmax=1,
#                    cmap='hot')


# plt.colorbar(pcm)
# plt.show()
# %%
nEF=int((EF-emin)/dE)




wannier_center2=np.zeros((3,length*nb))
for i in range(length):
    for j in range(nb):
        wannier_center2[:,j+i*nb]=wannier_center[:,j]+zhat*np.sqrt(np.dot(alpha,alpha))*i


ind=layer(list_coord=wannier_center,delta=delta_dist,L=L, nbr_of_layers=2,ns=length,alt=alt)
bot=wannier_center2[:,np.reshape(ind[:,0],-1)]
top=wannier_center2[:,np.reshape(ind[:,1],-1)]





mask=np.ones(nb*length,dtype=bool)
mask[np.reshape(ind[:,0],-1)[:nb//2]]= False
mask_opp=np.zeros(nb*length,dtype=bool)
mask_opp[np.reshape(ind[:,0],-1)[:nb//2]]= True
# I_1=I_/np.max(np.abs(I_[:,mask][np.reshape(ind[:,0],-1)[:nb//2],:]))
# J1=J/np.max(np.abs(J[:,mask][np.reshape(ind[:,0],-1)[:nb//2],:]))

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])

# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(I_1[:,mask][mask_opp,:]), vmax=np.max(I_1[:,mask][mask_opp,:]))

# for i in np.reshape(ind[:,0],-1)[:nb//2]:
#     for j in np.arange(nb*length)[mask]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=I_1[i,j]
#         color=cmap(norm(color_value))
#         intens=abs(I_1[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax)        

# plt.title('Current electron layer1torest'+label)
# plt.show()



# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])

# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(J1[:,mask][mask_opp,:]), vmax=np.max(J1[:,mask][mask_opp,:]))

# for i in np.reshape(ind[:,0],-1)[:nb//2]:
#     for j in np.arange(nb*length)[mask]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=J1[i,j]
#         color=cmap(norm(color_value))
#         intens=abs(J1[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax)
# plt.title('Current hole layer1torest'+label)
# plt.show()



mask2=np.ones(nb*length,dtype=bool)
mask2[np.reshape(ind[:,1],-1)[:nb//2]]= False
mask2_opp=np.zeros(nb*length,dtype=bool)
mask2_opp[np.reshape(ind[:,1],-1)[:nb//2]]= True
# I_2=I_/np.max(np.abs(I_[:,mask2][np.reshape(ind[:,1],-1)[:nb//2],:]))
# J2=J/np.max(np.abs(J[:,mask2][np.reshape(ind[:,1],-1)[:nb//2],:]))


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])


# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(I_2[:,mask2][mask2_opp,:]), vmax=np.max(I_2[:,mask2][mask2_opp,:]))
# for i in np.reshape(ind[:,1],-1)[:nb//2]:
#     for j in np.arange(nb*length)[mask2]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=I_2[i,j]
#         color=cmap(norm(color_value))

#         intens=abs(I_2[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax)

# plt.title('Current electron layer2torest'+label)
# plt.show()



# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])


# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(J2[:,mask2][mask2_opp,:]), vmax=np.max(J2[:,mask2][mask2_opp,:]))
# for i in np.reshape(ind[:,1],-1)[:nb//2]:
#     for j in np.arange(nb*length)[mask2]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=J2[i,j]
#         color=cmap(norm(color_value))

#         intens=abs(J2[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax) 
# plt.title('Current hole layer2torest'+label)
# plt.show()



mask3=np.ones(nb*length,dtype=bool)
mask3[np.reshape(ind[:,0],-1)[nb//2:]]= False
mask3_opp=np.zeros(nb*length,dtype=bool)
mask3_opp[np.reshape(ind[:,0],-1)[nb//2:]]= True
# I_3=I_/np.max(np.abs(I_[:,mask3][np.reshape(ind[:,1],-1)[:nb//2],:]))
# J3=J/np.max(np.abs(J[:,mask3][np.reshape(ind[:,1],-1)[:nb//2],:]))


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])


# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(I_3[:,mask3][mask3_opp,:]), vmax=np.max(I_3[:,mask3][mask3_opp,:]))
# for i in np.reshape(ind[:,0],-1)[nb//2:]:
#     for j in np.arange(nb*length)[mask3]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=I_3[i,j]
#         color=cmap(norm(color_value))

#         intens=abs(I_3[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax)
# plt.title('Current electron layer3torest'+label)
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])


# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(J3[:,mask3][mask3_opp,:]), vmax=np.max(J3[:,mask3][mask3_opp,:]))
# for i in np.reshape(ind[:,0],-1)[nb//2:]:
#     for j in np.arange(nb*length)[mask3]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=J3[i,j]
#         color=cmap(norm(color_value))

#         intens=abs(J3[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax) 
# plt.title('Current hole layer3torest'+label)

# plt.show()




















mask4=np.ones(nb*length,dtype=bool)
mask4[np.reshape(ind[:,1],-1)[nb//2:]]= False
mask4_opp=np.zeros(nb*length,dtype=bool)
mask4_opp[np.reshape(ind[:,1],-1)[nb//2:]]= True
I_4=I_/np.max(np.abs(I_[:,mask4][np.reshape(ind[:,1],-1)[:nb//2],:]))
J4=J/np.max(np.abs(J[:,mask4][np.reshape(ind[:,1],-1)[:nb//2],:]))



print('test',np.sum(I[:,mask4][mask,:]+I[:,mask][mask4,:]) )



# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])


# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(I_4[:,mask4][mask4_opp,:]), vmax=np.max(I_4[:,mask4][mask4_opp,:]))
# for i in np.reshape(ind[:,1],-1)[nb//2:]:
#     for j in np.arange(nb*length)[mask4]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=I_4[i,j]
#         color=cmap(norm(color_value))

#         intens=abs(I_4[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax)
# plt.title('Current electron layer4torest'+label)
# plt.show()





# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
# ax.scatter(xs=top[0], ys=top[1], zs=top[2])


# cmap = cm.seismic  # Colormap that goes from blue to red
# norm = colors.Normalize(vmin=np.min(J4[:,mask4][mask4_opp,:]), vmax=np.max(J4[:,mask4][mask4_opp,:]))
# for i in np.reshape(ind[:,1],-1)[nb//2:]:
#     for j in np.arange(nb*length)[mask4]:
#           # Ne pas tracer de vecteur de chaque point à lui-même
#         start_point = wannier_center2[:, i]
#         end_point = wannier_center2[:, j]
#         color_value=J4[i,j]
#         color=cmap(norm(color_value))

#         intens=abs(J4[i,j])
#         ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
#                     end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
#                     color=color,alpha=intens)
# set_axes_equal(ax) 
# plt.title('Current hole layer4torest'+label)
# plt.show()








###############################################################################
#Plot of the current from one layer to the rest
###############################################################################




fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I[:,mask][mask_opp,:]), vmax=np.max(I[:,mask][mask_opp,:]))

for i in np.reshape(ind[:,0],-1)[:nb//2]:
    for j in np.arange(nb*length)[mask]:

        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        
        intens=abs((I/(np.max(np.abs(I[:,mask][mask_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current  layer1torest'+label)
# plt.show()



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I[:,mask2][mask2_opp,:]), vmax=np.max(I[:,mask2][mask2_opp,:]))

for i in np.reshape(ind[:,1],-1)[:nb//2]:
    for j in np.arange(nb*length)[mask2]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        intens=abs((I/(np.max(np.abs(I[:,mask2][mask2_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current  layer2torest'+label)
# plt.show()


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I[:,mask3][mask3_opp,:]), vmax=np.max(I[:,mask3][mask3_opp,:]))

for i in np.reshape(ind[:,0],-1)[nb//2:]:
    for j in np.arange(nb*length)[mask3]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        

        intens=abs((I/(np.max(np.abs(I[:,mask3][mask3_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current  layer3torest'+label)
# plt.show()



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])
vmin=np.min(I[:,mask4][mask4_opp,:])
vmax=np.min(I[:,mask4][mask4_opp,:])
cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I[:,mask4][mask4_opp,:]), vmax=np.max(I[:,mask4][mask4_opp,:]))
print('min',np.min(I[:,mask4][mask4_opp,:]))
print('max',np.max(I[:,mask4][mask4_opp,:]))
for i in np.reshape(ind[:,1],-1)[nb//2:]:
    for j in np.arange(nb*length)[mask4]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        intens=abs((I/(np.max(np.abs(I[:,mask4][mask4_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current layer4torest'+label)
plt.show()







###############################################################################
#Plot of the current hole from one layer to the rest
###############################################################################




fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I_[:,mask][mask_opp,:]), vmax=np.max(I_[:,mask][mask_opp,:]))

for i in np.reshape(ind[:,0],-1)[:nb//2]:
    for j in np.arange(nb*length)[mask]:

        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I_)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        
        intens=abs((I_/(np.max(np.abs(I_[:,mask][mask_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current hole layer1torest'+label)
# plt.show()



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I_[:,mask2][mask2_opp,:]), vmax=np.max(I_[:,mask2][mask2_opp,:]))

for i in np.reshape(ind[:,1],-1)[:nb//2]:
    for j in np.arange(nb*length)[mask2]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I_)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        intens=abs((I_/(np.max(np.abs(I_[:,mask2][mask2_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current hole layer2torest'+label)
# plt.show()


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I_[:,mask3][mask3_opp,:]), vmax=np.max(I_[:,mask3][mask3_opp,:]))

for i in np.reshape(ind[:,0],-1)[nb//2:]:
    for j in np.arange(nb*length)[mask3]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I_)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        

        intens=abs((I_/(np.max(np.abs(I_[:,mask3][mask3_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current hole layer3torest'+label)
# plt.show()



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])
vmin=np.min(I_[:,mask4][mask4_opp,:])
vmax=np.min(I_[:,mask4][mask4_opp,:])
cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(I_[:,mask4][mask4_opp,:]), vmax=np.max(I_[:,mask4][mask4_opp,:]))
print('min',np.min(I_[:,mask4][mask4_opp,:]))
print('max',np.max(I_[:,mask4][mask4_opp,:]))
for i in np.reshape(ind[:,1],-1)[nb//2:]:
    for j in np.arange(nb*length)[mask4]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(I_)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        intens=abs((I_/(np.max(np.abs(I_[:,mask4][mask4_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current hole layer4torest'+label)
plt.show()






###############################################################################
#Plot of the current electron from one layer to the rest
###############################################################################




fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(J[:,mask][mask_opp,:]), vmax=np.max(J[:,mask][mask_opp,:]))

for i in np.reshape(ind[:,0],-1)[:nb//2]:
    for j in np.arange(nb*length)[mask]:

        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(J)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        
        intens=abs((J/(np.max(np.abs(J[:,mask][mask_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current electron layer1torest'+label)
# plt.show()



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(J[:,mask2][mask2_opp,:]), vmax=np.max(J[:,mask2][mask2_opp,:]))

for i in np.reshape(ind[:,1],-1)[:nb//2]:
    for j in np.arange(nb*length)[mask2]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(J)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        intens=abs((J/(np.max(np.abs(J[:,mask2][mask2_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current electron layer2torest'+label)
# plt.show()


fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])

cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(J[:,mask3][mask3_opp,:]), vmax=np.max(J[:,mask3][mask3_opp,:]))

for i in np.reshape(ind[:,0],-1)[nb//2:]:
    for j in np.arange(nb*length)[mask3]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(J)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        

        intens=abs((J/(np.max(np.abs(J[:,mask3][mask3_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current electron layer3torest'+label)
# plt.show()



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xs=bot[0], ys=bot[1], zs=bot[2])
ax.scatter(xs=top[0], ys=top[1], zs=top[2])
vmin=np.min(J[:,mask4][mask4_opp,:])
vmax=np.min(J[:,mask4][mask4_opp,:])
cmap = cm.seismic  # Colormap that goes from blue to red
norm = colors.Normalize(vmin=np.min(J[:,mask4][mask4_opp,:]), vmax=np.max(J[:,mask4][mask4_opp,:]))
print('min',np.min(J[:,mask4][mask4_opp,:]))
print('max',np.max(J[:,mask4][mask4_opp,:]))
for i in np.reshape(ind[:,1],-1)[nb//2:]:
    for j in np.arange(nb*length)[mask4]:
          # Ne pas tracer de vecteur de chaque point à lui-même
        start_point = wannier_center2[:, i]
        end_point = wannier_center2[:, j]
        color_value=(J)[i,j]
        if color_value>0:
            color='red'
        else:
            color='blue'
        intens=abs((J/(np.max(np.abs(J[:,mask4][mask4_opp,:]))))[i,j])
        ax.quiver(start_point[0], start_point[1], start_point[2],  # Début du vecteur
                    end_point[0] - start_point[0], end_point[1] - start_point[1], end_point[2] - start_point[2],  # Direction du vecteur
                    color=color,alpha=intens)
set_axes_equal(ax)        

plt.title('Current electron layer4torest'+label)
plt.show()