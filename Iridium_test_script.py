#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 17:41:07 2021

@author: simon
"""

import pandas as pd
from SECpy import SECS, RE
from SECpy import PlottingTools as PT
Geocentric_to_PlateCarree_vector_components, features= PT.Geocentric_to_PlateCarree_vector_components, PT.features
from apexpy import Apex
import datetime as dt
from secsy import cubedsphere as CS
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib as mpl
mpl.rcParams.update({'xtick.labelsize': 5, 'ytick.labelsize': 5, 
                     'xtick.major.pad':0.1, 'xtick.minor.pad':0,
                     'ytick.major.pad':0.1, 'ytick.minor.pad':0,
                     'xtick.major.size': 2.5, 'xtick.minor.size': 1.0,
                     'ytick.major.size': 2.5, 'ytick.minor.size': 1.0,
                     'font.size': 7.5, 'lines.linewidth': 0.5,
                     'lines.markersize': 2.5, 'axes.linewidth': 0.2})

sat_data_path= '/home/simon/BCSS-DAG Dropbox/Simon Walker/Iridium/All_Data.hdf5'
max_date= np.datetime64('2017-09-30T23:59:59.785107456')
min_date= np.datetime64('2009-10-02T00:00:00.682729472')
Sat_Store= pd.HDFStore(sat_data_path, mode='r')

mag_data_path= '/home/simon/gits/SECpy/SECS_Visualisation/SECS_Data/Scandinavia_HDF.h5'
Mag_Store= pd.HDFStore(mag_data_path, mode='r')

MagLat= np.array([54.61, 55.62, 59.9 , 60.5 , 62.3 , 64.52, 64.94, 66.9 , 67.37,
        68.02, 68.35, 69.02, 69.3 , 69.66, 69.76, 70.54, 74.5 , 77.  ,
        78.2 , 78.92])
MagLon= np.array([18.82, 11.67, 17.35, 24.65, 26.65, 27.23, 10.99, 24.08, 26.63,
        23.53, 18.82, 20.79, 16.03, 18.94, 27.01, 22.22, 19.2 , 15.6 ,
        15.83, 11.95])
sites=np.array(['LYR', 'BJN', 'TRO', 'HRN', 'ABK', 'SOD', 'NUR', 'BFE', 'NAL',
        'SOR', 'AND', 'KEV', 'KIL', 'MUO', 'PEL', 'OUJ', 'HAN', 'UPS',
        'HLP', 'RVK'], dtype='<U3')
lon_centre= 17.7
lat_centre= 68.1
"""Create figure"""
fig = plt.figure(constrained_layout=False, num='Spherical Elementary Currents')
gs = fig.add_gridspec(4, 9, width_ratios= [.25, 2,.25, .25, 2,.25, .25, 2, .25], height_ratios=[1, .05, 1, .05], wspace=0.1, hspace=0.2)
lon_centre= 17.7
lat_centre= 68.1
proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
axes_curl = [fig.add_subplot(gs[0,:3], projection= proj), fig.add_subplot(gs[0,3:6], projection= proj), 
                 fig.add_subplot(gs[0,6:], projection= proj)]
axes_div= [fig.add_subplot(gs[2, i*3:(i+1)*3], projection=proj) for i in range(3)]
for ax in axes_curl:
    ax.coastlines()
axes_curl[0].set_title('Magnetic Field')
axes_curl[1].set_title('Pole Amplitude')
axes_curl[2].set_title('Ionospheric Currents')
for ax in axes_div:
    ax.coastlines()
fig.text(0, 0.75, 'Curl-Free')
fig.text(0, 0.25, 'Divergence-Free')
cax_curl= fig.add_subplot(gs[1, 4])
caxes_div= [fig.add_subplot(gs[3, 3*(i)+1]) for i in range(3)]
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

min_lon= min(node_lons)-5
max_lon= max(node_lons)+5
min_lat= min(node_lats)-5
# max_lat= max(node_lats)
max_lat= 90
mlat=np.linspace(49, 81, 50)
merid_glat, merid_glon, err= A.apex2geo(alat=mlat, alon= 105, height=110)
Date= pd.Timestamp('2016-02-08T02:10')
Date1= Date - pd.Timedelta(value=60, unit='m')
Date2= Date + pd.Timedelta(value=60, unit='m')
Sat= Sat_Store.select(key='main', where='Date_UTC>=Date1 & Date_UTC<=Date2 \
                      & glon>= min_lon & glon<=max_lon & glat>=min_lat & glat<=max_lat')
"""Regularisiation"""
depth=500
λ1= 1e-23
λ2= λ1*1e2
poles= SECS(node_lons, node_lats, lon, lat, mode='image', image_current_radius=RE-depth*1E3)
Le, Ln=node_grid.get_Le_Ln()
node_f1, node_f2= A.basevectors_qd(node_grid.lat.flatten(), node_grid.lon.flatten(), 110, coords='geo')
e= node_f1/np.linalg.norm(node_f1, axis=0)
L= np.diag(e[0]).dot(Le.toarray()) + np.diag(e[1]).dot(Ln.toarray())

# Curl free
poles.Fitting_Matrix(Sat.glon.values, Sat.glat.values, eval_radius=RE+800e3, system='curl-free')

# Div free
G=poles.G_Matrix(MagLon, MagLat, system='divergence-free')
GTG= np.dot(G.T, G)
matrix= GTG + λ1*np.identity(GTG.shape[0]) +λ2*np.dot(L.T, L)/np.max(np.abs(np.dot(L.T, L)))
Inverse= poles.Inverse_Matrix(matrix, cond=0)
poles.fitting_matrix['divergence-free']= np.dot(Inverse,G.T)



"""Run SECS Curl"""
# Magnetic Field
Modelled=poles.Magnetic_Field(-Sat.Bn.values, Sat.Be.values, Sat.Bu.values, eval_radius=RE+800e3, system='curl-free')
modQ=axes_curl[0].quiver(lon[::8], lat[::8], 
               *Geocentric_to_PlateCarree_vector_components(Modelled[1][::8], Modelled[2][::8], lat[::8]), 
               transform=ccrs.PlateCarree())
MesQ=axes_curl[0].quiver(Sat.glon.values, Sat.glat.values, 
               *Geocentric_to_PlateCarree_vector_components(Sat.Be.values, Sat.Bn.values, Sat.glat.values),
               transform=ccrs.PlateCarree(), color='red')
# SECS Amplitudes
p=axes_curl[1].pcolormesh(node_grid.lon_mesh, node_grid.lat_mesh, 
                   poles.Amplitude(-Sat.Bn.values, Sat.Be.values, Sat.Bu.values, system='curl-free').reshape(node_grid.shape), cmap='seismic',
                   transform=ccrs.PlateCarree())
# Current
Current= poles.Currents(-Sat.Bn.values, Sat.Be.values, Sat.Bu.values, system='curl-free')
CurQ=axes_curl[2].quiver(lon, lat, *Geocentric_to_PlateCarree_vector_components(Current[0], Current[1], lat),
                         transform=ccrs.PlateCarree(), scale=5e-7)

# Color bars
cbar_curl=fig.colorbar(p, cax=cax_curl, orientation='horizontal')
cbar_curl.set_label('Curl-Free Ionospheric Pole Amplitude')
"""Run SECS Div"""
Mag= Mag_Store.select(key='Main', where='Date_UTC=Date')
index= Mag.Site.isin(sites)
if np.sum(index)!= len(sites):
    raise ValueError('Not all sites are available')
Btheta= -Mag.Btheta[index].values
Bphi= Mag.Bphi[index].values
Br= Mag.Br[index].values
del Mag

# Magnetic Field
Modelled= poles.Magnetic_Field(Btheta, Bphi, Br, eval_radius=RE)
Modelled= np.sum(np.array(Modelled), axis=0)
modQ= axes_div[0].quiver(lon[::8], lat[::8], 
                         *Geocentric_to_PlateCarree_vector_components(Modelled[1][::8], Modelled[2][::8], lat[::8]),
                         transform= ccrs.PlateCarree(), scale=1e3, zorder=100)
MesQ= axes_div[0].quiver(MagLon, MagLat, 
                         *Geocentric_to_PlateCarree_vector_components(Bphi, -Btheta, MagLat),
                         transform=ccrs.PlateCarree(), color='red', scale=1e3, zorder=100)
BrP=axes_div[0].pcolormesh(eval_grid.lon_mesh, eval_grid.lat_mesh, Modelled[0].reshape(eval_grid.shape),
                       transform=ccrs.PlateCarree(), cmap='seismic')

# SECS Amplitudes
p=axes_div[1].pcolormesh(node_grid.lon_mesh, node_grid.lat_mesh, 
                   poles.Amplitude(Btheta, Bphi, Br).reshape(node_grid.shape), 
                   cmap='seismic', transform=ccrs.PlateCarree())
# Currents
Current= poles.Currents(Btheta, Bphi, Br)
CurQ=axes_div[2].quiver(lon, lat, *Geocentric_to_PlateCarree_vector_components(Current[0], Current[1], lat),
                        transform=ccrs.PlateCarree(),scale=5e9)

# Color bars
cbar_mag= fig.colorbar(BrP, cax=caxes_div[0], orientation='horizontal')
cbar_mag.set_label('Modelled Radial Magnetic Field')
cbar_Amp= fig.colorbar(p, cax=caxes_div[1], orientation='horizontal')
cbar_Amp.set_label('Divergence-Free Ionospheric Pole Amplitude')

""" Adjust limits"""
extent= axes_div[1].get_extent()
for ax in axes_curl+ axes_div:
    ax.set_extent(extent, crs=proj)