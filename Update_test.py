#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 16:31:30 2021

@author: simon
"""

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from SECpy import BphiCF, theta, RE, Local2Global, PlottingTools, Deg2Rad, V, SECS
#%% Base functions check
fig= plt.figure()

lon_centre= 19.6
lat_centre= 68.7
proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)

ax= fig.add_subplot(221, projection=proj)
lats, lons= np.meshgrid(np.linspace(70, 90), np.linspace(10, 30))
thetap= theta(80, 20, lats, lons)
Bn= BphiCF(thetap, RE+800e3, RE+110e3)
Bn, Be= Local2Global(Deg2Rad(90-80), Deg2Rad(90-lats), thetap, Deg2Rad(20), Deg2Rad(lons), Bn, 0)
ax.scatter(20, 80, transform= ccrs.PlateCarree())
ax.quiver(lons, lats, 
          *PlottingTools.Geocentric_to_PlateCarree_vector_components(Be, -Bn, lats), 
          transform=ccrs.PlateCarree())
ax.set_extent([10, 30, 79, 81])
ax.set_title('Magnetic Field SECpy')
ax2= fig.add_subplot(222, projection=proj)
ax2.scatter(20, 80, transform= ccrs.PlateCarree())
Cur=V(thetap, RE+110e3)
Curn, Cure= Local2Global(Deg2Rad(90-80), Deg2Rad(90-lats), thetap, Deg2Rad(20), Deg2Rad(lons), 0, Cur)
ax2.quiver(lons, lats, 
          *PlottingTools.Geocentric_to_PlateCarree_vector_components(Cure, -Curn, lats), 
          transform=ccrs.PlateCarree())
ax2.set_extent([10, 30, 79, 81])
ax2.set_title('Current SECpy')
import secsy
ax3= fig.add_subplot(223, projection=proj)
e, n, r=secsy.utils.get_SECS_B_G_matrices(lats, lons, RE+800e3, 80, 20, RI= RE+110e3, current_type='curl_free')
ax3.scatter(20, 80, transform=ccrs.PlateCarree())
ax3.quiver(lons.flatten(), lats.flatten(), 
           *PlottingTools.Geocentric_to_PlateCarree_vector_components(e[:,0], n[:,0], lats.flatten()),
           transform=ccrs.PlateCarree())
ax3.set_extent([10, 30, 79, 81])
ax3.set_title('Magnetic Field secsy')


ax4= fig.add_subplot(224, projection=proj)
e, n=secsy.utils.get_SECS_J_G_matrices(lats, lons, 80, 20, current_type='curl_free')
ax4.scatter(20, 80, transform=ccrs.PlateCarree())
ax4.quiver(lons.flatten(), lats.flatten(), 
           *PlottingTools.Geocentric_to_PlateCarree_vector_components(e[:,0], n[:,0], lats.flatten()),
           transform=ccrs.PlateCarree())
ax4.set_extent([10, 30, 79, 81])
ax4.set_title('Current secsy')
plt.subplots_adjust(wspace=0, hspace=0.2)
#%% Input Check
fig= plt.figure()

lon_centre= 19.6
lat_centre= 68.7
Btheta, Bphi, Br= np.array([1]), np.array([1]), np.array([0])
proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
grid_lons, grid_lats = np.meshgrid(np.linspace(0, 40, 100), np.linspace(70,90, 100))
grid_lons= grid_lons.flatten()
grid_lats= grid_lats.flatten()
pole= SECS(grid_lons, grid_lats, lons.flatten(), lats.flatten())
pole.Fitting_Matrix(np.array([25]), np.array([80.5]), cond=0, eval_radius=RE+800e3, system='curl-free')
Mag= pole.Magnetic_Field(Btheta, Bphi, Br, eval_radius=RE+800e3, system='curl-free', eval_longitude=np.array([25]), eval_latitude=np.array([80.5]))
ax= fig.add_subplot(111, projection=proj)
ax.scatter(grid_lons, grid_lats, transform=ccrs.PlateCarree())
ax.quiver(np.array([25]), np.array([80.5]), 
          *PlottingTools.Geocentric_to_PlateCarree_vector_components(Mag[1], Mag[2], np.array([80.5])),
          transform= ccrs.PlateCarree())
ax.quiver(np.array([25]), np.array([80.5]), *PlottingTools.Geocentric_to_PlateCarree_vector_components(Bphi, -Btheta, np.array([80.5])), 
          transform=ccrs.PlateCarree(), color='red')