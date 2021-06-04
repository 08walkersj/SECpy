#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  7 15:25:51 2021

@author: simon
"""

# %% Imports
import numpy as np
from SECpy import RE, SECS
from apexpy import Apex
import datetime as dt
import secsy
CS= secsy.cubedsphere 
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
# import pickle
mpl.rcParams.update({'xtick.labelsize': 5, 'ytick.labelsize': 5, 
                     'xtick.major.pad':0.1, 'xtick.minor.pad':0,
                     'ytick.major.pad':0.1, 'ytick.minor.pad':0,
                     'xtick.major.size': 2.5, 'xtick.minor.size': 1.0,
                     'ytick.major.size': 2.5, 'ytick.minor.size': 1.0,
                     'font.size': 7.5, 'lines.linewidth': 0.5,
                     'lines.markersize': 2.5, 'axes.linewidth': 0.2})
# mpl.rcParams.update(pickle.load(open('/home/simon/Documents/Matplotlib/rc_params.pickle', 'rb')))
# %% Create Subplots
sites=np.array(['LYR', 'BJN', 'TRO', 'HRN', 'ABK', 'SOD', 'NUR', 'BFE', 'NAL',
        'SOR', 'AND', 'KEV', 'KIL', 'MUO', 'PEL', 'OUJ', 'HAN', 'UPS',
        'HLP', 'RVK'], dtype='<U3')
MagLat= np.array([54.61, 55.62, 59.9 , 60.5 , 62.3 , 64.52, 64.94, 66.9 , 67.37,
        68.02, 68.35, 69.02, 69.3 , 69.66, 69.76, 70.54, 74.5 , 77.  ,
        78.2 , 78.92])
MagLon= np.array([18.82, 11.67, 17.35, 24.65, 26.65, 27.23, 10.99, 24.08, 26.63,
        23.53, 18.82, 20.79, 16.03, 18.94, 27.01, 22.22, 19.2 , 15.6 ,
        15.83, 11.95])
"""Create figure"""
fig = plt.figure(constrained_layout=False, num='Spherical Elementary Currents')
gs = fig.add_gridspec(3, 4, width_ratios= [2,2,2,1], height_ratios= [20, 0.1, .9])
lon_centre= 17.7
lat_centre= 68.1
# fig.add_subplot(gs[15:70, 3]) #old dimensions
axes = [fig.add_subplot(gs[:-2,0]), fig.add_subplot(gs[:-2,1]), 
                 fig.add_subplot(gs[:-2,2]),fig.add_subplot(gs[:-2, 3])]
caxes= [fig.add_subplot(gs[-1, i]) for i in range(3)]
# %% Set grid
"""Grid"""
A = Apex(date=dt.datetime(2008, 6, 1, 0, 0, 0))
f1, f2 = A.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
qd_north = f2 / np.linalg.norm(f2)
East, North= qd_north[0], qd_north[1]
Gridproj= CS.CSprojection((lon_centre, lat_centre), [East, North])
node_grid=CS.CSgrid(Gridproj, 3700, 2200, 50., 50.)
Evalproj= CS.CSprojection((lon_centre-0.23, lat_centre+0.23), [East, North])
eval_grid= CS.CSgrid(Evalproj,  3300, 1800, 50., 50.)
node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
lat,lon= eval_grid.lat.flatten(), eval_grid.lon.flatten()
scat=axes[0].scatter(node_grid.xi, node_grid.eta, s=1, color='red')
xlim, ylim = axes[0].get_xlim(), axes[0].get_ylim()
for ax in axes[:-2]:
    for cl in Gridproj.get_projected_coastlines(resolution = '110m'):
        ax.plot(*cl, color='black',zorder=-1)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('xi')
for ax in axes[1:-2]:
    ax.yaxis.set_visible(False)
axes[-2].xaxis.set_visible(False)
axes[0].set_ylabel('eta')
axes[0].set_title('Measured Magnetic Field')
axes[1].set_title('Evaluated Magnetic Field')
axes[2].set_title('Ionospheric Current')
axes[3].set_title(r'Magnetic Meridian ($105^o$)')
axes[4].set_title('Validity', size=6)
fig.suptitle('SECS Solution for One Minute of Data')
# %% Set date
"""Check Date and get Mag locations"""
Date= pd.Timestamp('2000-02-05T22:34')
data_path= '/home/simon/gits/SECpy/SECS_Visualisation/SECS_Data/Scandinavia_HDF.h5' #input data path
# %% Regularisation
"""Regularisiation"""
depth=500
λ1= 1e-23
λ2= λ1*1e2
poles= SECS(node_lons, node_lats, lon, lat, mode='image', image_current_radius=RE-depth*1E3)
Le, Ln=node_grid.get_Le_Ln()
node_f1, node_f2= A.basevectors_qd(node_grid.lat.flatten(), node_grid.lon.flatten(), 110, coords='geo')
e= node_f1/np.linalg.norm(node_f1, axis=0)
L= np.diag(e[0]).dot(Le) + np.diag(e[1]).dot(Ln)
G=poles.G_Matrix(MagLon, MagLat)
GTG= np.dot(G.T, G)
matrix= GTG + λ1*np.identity(GTG.shape[0]) +λ2*np.dot(L.T, L)/np.max(np.abs(np.dot(L.T, L)))
Inverse= poles.Inverse_Matrix(matrix, cond=0)
poles.fitting_matrix['divergence-free']= np.dot(Inverse,G.T)
# %% Meridian
"""Meridian"""
mlat= np.linspace(49, 81, 50)
merid_glat, merid_glon, err= A.apex2geo(alat=mlat, alon= 105, height=0)
f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3= A.basevectors_apex(merid_glat, merid_glon, 0, coords='geo')
MagMlat, MagMlon= A.geo2apex(MagLat, MagLon, 0)
#%% Features
#Colorbars
# mappable0=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-10E3, vmax=10E3,), cmap='seismic')
# cbar0=fig.colorbar(mappable=mappable0, cax=caxes[0], orientation='horizontal', extend='both')
# cbar0.set_label(label='Amplitudes (A)', rotation=0)
# cbar0.ax.locator_params(nbins=7)

mappable1=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=0.012*1e3, vmax= .5*1e3), cmap=mpl.cm.afmhot_r)
cbar1=fig.colorbar(mappable=mappable1, cax=caxes[2], orientation='horizontal', extend='both')
cbar1.set_label(label='Ionospheric Current Magnitude '+r'$Akm^{-1}$', rotation=0)

mappable2=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-4.0E2, vmax=4.0E2,), cmap=mpl.cm.bwr)
cbar2=fig.colorbar(mappable=mappable2, cax=caxes[0], orientation='horizontal', extend='both')
cbar2.set_label(label='Evaluated Br (nT)', rotation=0)

#Add evaluation meridian
axes[1].plot(*Gridproj.geo2cube(merid_glon, merid_glat), color ='blue', alpha=0.6, zorder=200)
axes[2].plot(*Gridproj.geo2cube(merid_glon, merid_glat), color ='blue', alpha=0.6, zorder=200)
axes[-2].plot([0]*len(mlat), mlat)
for ax in axes[:-1]:
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
axes[-1].spines['top'].set_visible(False)
axes[-1].spines['right'].set_visible(False)
#Add y=x line
axes[-1].plot(np.linspace(-1e3,1e3), np.linspace(-1e3, 1e3))
#Resize the figure
DPI=fig.get_dpi()
fig.set_size_inches(2400.0/float(DPI),1220.0/float(DPI))

#%% RUN SECS

DataFrame= pd.read_hdf(data_path, mode='r', where="Date_UTC=Date" )
DataFrame=DataFrame.sort_values(by='glat')
index= DataFrame['Site'].isin(sites)
Btheta, Bphi, Br= -DataFrame['Btheta'][index].values*1e-9,DataFrame['Bphi'][index].values*1e-9,-DataFrame['Br'][index].values*1e-9
if len(Btheta)!= len(sites):
    raise ValueError("not all magnetometers available: ")
scatter=[]
#Currents
CI= poles.Currents(Btheta, Bphi, Br)
CMag=np.sqrt((CI[0]**2) +CI[1]**2)
xi, eta, CurE, CurN= Gridproj.vector_cube_projection(*CI,lon, lat)
CurQ=axes[2].quiver(xi[::3], eta[::3], CurE[::3], CurN[::3], zorder=90,scale=(1e1)/2, color='Black', label='Ionospheric Current', alpha=0.7)
scatter.append(axes[2].scatter(xi, eta, c= CMag*1e3, vmin=0.012*1e3, vmax= .5*1e3, cmap= 'afmhot_r',  zorder=50, alpha=0.5))
#Magnetic field
MI, MT = poles.Magnetic_Field(Btheta, Bphi, Br)
M= (MI[0] +MT[0], MI[1] +MT[1], MI[2] +MT[2])
xi, eta, MagE, MagN= Gridproj.vector_cube_projection(*M[1:],lon, lat)
Mag=axes[1].quiver(xi[::8], eta[::8], MagE[::8], MagN[::8], zorder=80,
                    scale=(1e-5)/8, color='Black', label='Total Field', alpha=0.7)
scatter.append(axes[1].scatter(xi, eta, c= M[0]*1E9, vmin=-4.0E2, vmax=4.0E2,
                                cmap='bwr', zorder=50, alpha=.5))
Measured=[None, None]
xi, eta, Measured[0], Measured[1]= Gridproj.vector_cube_projection(Bphi, -Btheta, MagLon, MagLat)
MagQ=axes[0].quiver(xi, eta, *Measured,zorder=90, 
                     scale = (1e-5)/8, color='red', label='Measurements', alpha=0.7)
axes[0].scatter(xi, eta, c= Br*1E9, vmin=-4.0E2, vmax=4.0E2,
                                cmap='bwr', zorder=50, alpha=.5)
#Model vs Measured
model_I, model_T= poles.Magnetic_Field(Btheta, Bphi, Br, eval_longitude= MagLon, eval_latitude= MagLat)
model= (model_I[0] +model_T[0], model_I[1] +model_T[1], model_I[2] +model_T[2])
model= model[0], *Gridproj.vector_cube_projection(*model[1:], MagLon, MagLat, return_xi_eta=False)
scatter.append(axes[-1].scatter(Br*1e9, model[0]*1e9, label='Br component', color= 'purple'))
scatter.append(axes[-1].scatter(Measured[0]*1e9, model[1]*1e9, label='Be component', color='green'))
scatter.append(axes[-1].scatter(Measured[1]*1e9, model[-1]*1e9, label='Bn component', color='red'))
lim=max(abs(np.append(Measured[0]*1e9, Measured[1]*1e9)))
axes[-1].set_xlim((-lim*1.2, lim*1.2))
axes[-1].set_ylim((-lim*1.2, lim*1.2))
#Meridian Current
meridian= poles.Currents(Btheta, Bphi, Br, eval_longitude= merid_glon, eval_latitude= merid_glat, singularity_limit= 50e3)
poles.Magnetic_Field(Btheta, Bphi, Br, eval_longitude= merid_glon, eval_latitude= merid_glat)
qd_east= g1[0]*meridian[0] +g1[1]*meridian[1]
qd_north= g2[0]*meridian[0] +g2[1]*meridian[1]
MeridCur=axes[-2].quiver([0]*len(mlat), mlat, qd_east, qd_north, scale=1., color='black', label='Meridian Current')

Curqk= axes[2].quiverkey(CurQ, 0.7, 0.17, 2000e-3, r'$2000 \ Akm^{-1}$', labelpos='E',
               coordinates='figure', zorder=500)
Magqk= axes[0].quiverkey(MagQ, 0.45, 0.17, 400e-9, r'$400 \ nT$', labelpos='E',
                   coordinates='figure', zorder=500)
scat.remove()