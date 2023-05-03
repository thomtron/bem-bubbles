from plyfile import PlyData, PlyElement

import meshcut

import numpy as np
from numpy import genfromtxt

import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True


# NOTE: to use this file, you first have to run a simulation with main.cpp and 
# then compute external potential data with pot-ext.cpp. Some of the parameters
# used here are specifically adopted to the parameters in those two files - if 
# you made changes in those two files, you may have to reproduce them here!

# This script plots the external potential and the pressure field of a specific
# simulation output (here mesh-113.ply generated from main.cpp).



def load_ply(fname):
    plydata = PlyData.read(fname)
    verts_raw = plydata.elements[0].data
    faces_raw = plydata.elements[1].data

    verts = np.empty(shape=(len(verts_raw),3))
    phi = np.empty(shape=(len(verts_raw),))
    psi = np.empty(shape=(len(verts_raw),))
    for i in range(len(verts)):
        verts[i] = np.array([verts_raw[i][0],verts_raw[i][1],verts_raw[i][2]])
        phi[i] = verts_raw[i][3]
        psi[i] = verts_raw[i][4]

    faces = np.empty(shape=(len(faces_raw),3),dtype=np.int32)
    for i in range(len(faces)):
        faces[i] = faces_raw[i][0]

    return verts,faces ,phi,psi


# dimensions of the grid, where the potential values were sampled
# see phi-ext.cpp
N = 801
M = 301

s = (N,M)

# loading the potential
data = genfromtxt("../../build/pot-ext.csv",delimiter=';')
pot = data[:,3].reshape(s)

# loading the time derivative of the potential
data = genfromtxt("../../build/pot_t-ext.csv",delimiter=';')
pot_t = data[:,3].reshape(s)

x = (data[:,0].reshape(s))[:,0]
z = (data[:,2].reshape(s))[0,:]


grad = np.gradient(pot,x,z) # space derivative here!

pressure = 1.0 - (pot_t + 0.5*(grad[0]**2+grad[1]**2))



fig = plt.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)


X,Z = np.meshgrid(x,z)

p_pot = ax1.pcolormesh(x,-z,pot.T,vmin=-0.12,vmax=0.08,cmap='plasma',rasterized=True)

cax = ax1.inset_axes([1.02,0.0, 0.02, 1.0])
cbar_left = fig.colorbar(p_pot,ax=ax1,cax=cax,orientation='vertical')


# pressure field:
p_pressure = ax2.pcolormesh(x,z,pressure.T,cmap='viridis',zorder=1,rasterized=True,vmin=0.7,vmax=1.3)
ax2.contour(x-0.005,z-0.005,pressure.T,[0.9,1,1.1],colors='k')

cax = ax2.inset_axes([1.02,0.0, 0.02, 1.0])
cbar_right = fig.colorbar(p_pressure,ax=ax2,cax=cax,orientation = 'vertical')


# plot the cut of the mesh through the plane on which we sampled potential and pressure:

# origin and normal of cutting plane (see meshcut documentation)
plane_orig = (0, 0, 0)
plane_norm = (0, 1, 0)

plane = meshcut.Plane(plane_orig, plane_norm)

verts,triangles,phi,psi = load_ply('../../build/meshes/mesh-113.ply') # needs to be adapted!
mesh = meshcut.TriangleMesh(verts,triangles)

sections = meshcut.cross_section_mesh(mesh,plane)

for section in sections:
    section = np.append(section,[section[0]],axis=0)
    ax1.plot(section[:,0],section[:,2],'w-',linewidth=2)
    ax1.fill(section[:,0],section[:,2],'k-')

    ax2.plot(section[:,0],section[:,2],'w-',linewidth=2)
    ax2.fill(section[:,0],section[:,2],'k-')


ax1.set_aspect('equal')
ax2.set_aspect('equal')

ax1.set_xlim(-4,4)
ax1.set_ylim(-1.5,1.5)
ax2.set_xlim(-4,4)
ax2.set_ylim(-1.5,1.5)

fontsize = 16

ax1.set_ylabel(r'$z_*$',fontsize=fontsize)
ax2.set_xlabel(r'$x_*$',fontsize=fontsize)
ax2.set_ylabel(r'$z_*$',fontsize=fontsize)

cbar_left.ax.set_ylabel(r'$u_*$',fontsize=fontsize,rotation=0)
cbar_right.ax.set_ylabel(r'$p_*$',fontsize=fontsize,rotation=0)


plt.show()