#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:04:31 2020
@author: Simon Walker

"""
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
from SECpy import SECS, RE
from SECpy import PlottingTools as PT
Geocentric_to_PlateCarree_vector_components, features= PT.Geocentric_to_PlateCarree_vector_components, PT.features
from apexpy import Apex
import os
import datetime as dt
import cubedsphere as CS
import matplotlib.pyplot as plt
import matplotlib as mpl
from urllib.error import HTTPError
import pylab as plot
params = {'legend.fontsize': 5,
          'legend.handlelength': 2}
plot.rcParams.update(params)
data_path= './SECS_Visualisation/SECS_Data/Scandinavia_HDF.h5' #input data path
sites=np.array(['LYR', 'BJN', 'TRO', 'HRN', 'ABK', 'SOD', 'NUR', 'BFE', 'NAL',
        'SOR', 'AND', 'KEV', 'KIL', 'MUO', 'PEL', 'OUJ', 'HAN', 'UPS',
        'HLP', 'RVK'], dtype='<U3') #List of sites used

"""Create figure"""
fig = plt.figure(constrained_layout=False, num='Spherical Elementary Currents')
gs = fig.add_gridspec(100, 4, width_ratios= [2,2,2,1])
lon_centre= 17.7
lat_centre= 68.1
proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
axes = [fig.add_subplot(gs[:,0], projection= proj), fig.add_subplot(gs[:,1], projection= proj), 
                 fig.add_subplot(gs[:,2], projection= proj),fig.add_subplot(gs[5:70, 3]), fig.add_subplot(gs[75:90, 3], aspect='auto')]
"""Check Date and get Mag locations"""
start_Date= pd.Timestamp('2000-09-21T22:33') #Input timestamp that you want
end_Date= pd.Timestamp('2000-02-05T22:33') #Input timestamp that you want

plots_list=[]
DataFrame= pd.read_hdf(data_path, mode='r', where="Date_UTC=start_Date" )

index= DataFrame['Site'].isin(sites)
MagLat=DataFrame['glat'][index].values
MagLon=DataFrame['glon'][index].values

#Check if magnetometer selection is available
if len(MagLon)!=len(sites):
    raise ValueError("not all magnetometers available!")


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
poles.fitting_matrix= np.dot(Inverse,G.T)


"""Meridian"""
mlat= np.linspace(49, 81, 50)
merid_glat, merid_glon, err= A.apex2geo(alat=mlat, alon= 105, height=0)
f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3= A.basevectors_apex(merid_glat, merid_glon, 0, coords='geo')
proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
MagMlat, MagMlon= A.geo2apex(MagLat, MagLon, 0)

"""Changing features of the plot"""
def PlottingFunction(date, poles):
    """Normal Scales:
        Mag quivers (1e-5)/8
        Cur quivers (1e1)/2
        Merid quivers 2.5
        Amplitudes -10E3 10E3
        Br -6E7*1e-9 6E7*1e-9
        Cur 0.012*1e3 1.5*1e3
        """
    global Btheta, Bphi, Br
    DataFrame= pd.read_hdf(data_path, mode='r', where="Date_UTC=date" )
    index= DataFrame['Site'].isin(sites)
    Btheta, Bphi, Br= -DataFrame['Btheta'][index].values*1e-9,DataFrame['Bphi'][index].values*1e-9,-DataFrame['Br'][index].values*1e-9
    #Check if magnetometer selection is available
    if len(Btheta)!=len(sites):
        raise ValueError("not all magnetometers available!")
    scatter=[]
    #Currents
    CI= poles.Currents(Btheta, Bphi, Br)
    CMag=np.sqrt((CI[0]**2) +CI[1]**2)
    CurE, CurN= Geocentric_to_PlateCarree_vector_components(*CI,lat)
    CurQ=axes[1].quiver(lon[::3], lat[::3], CurE[::3], CurN[::3], zorder=90,scale=(1e0)/2, transform = ccrs.PlateCarree(), color='Black', label='Ionospheric Current', alpha=0.7)
    scatter.append(axes[1].scatter(lon, lat, c= CMag*1e3, vmin=0.001*1e3, vmax= 0.018*1e3, cmap= mpl.cm.afmhot_r, transform= ccrs.PlateCarree(), zorder=50))
    #Magnetic field
    MI, MT = poles.Magnetic_Field(Btheta, Bphi, Br)
    M= (MI[0] +MT[0], MI[1] +MT[1], MI[2] +MT[2])
    MagE, MagN= Geocentric_to_PlateCarree_vector_components(*M[1:],lat)
    Mag=axes[2].quiver(lon[::8], lat[::8], MagE[::8], MagN[::8], zorder=80,
                        scale=(1e-6)/8, transform = ccrs.PlateCarree(), color='Black', label='Total Field', alpha=0.7)
    Measured= Geocentric_to_PlateCarree_vector_components(Bphi, -Btheta, MagLat)
    MagQ=axes[2].quiver(MagLon, MagLat,  *Geocentric_to_PlateCarree_vector_components(Bphi, -Btheta, MagLat),zorder=90, 
                         scale = (1e-6)/8, transform = ccrs.PlateCarree(), color='red', label='Measurements', alpha=0.7)
    scatter.append(axes[2].scatter(lon, lat, c= M[0]*1E9, vmin=-10.0E-9 *1e9, vmax=10.0E-9 *1e9,
                                    cmap=mpl.cm.bwr, transform= ccrs.PlateCarree(), zorder=50))
    #SEC pole amplitudes
    scatter.append(axes[0].scatter(node_lons, node_lats, c= poles.Amplitude(Btheta, Bphi, Br), vmin=-10E1, vmax=10E1, 
                                   cmap=mpl.cm.seismic, transform= ccrs.PlateCarree(), zorder=50, s=1))
    #Model vs Measured
    model_I, model_T= poles.Magnetic_Field(Btheta, Bphi, Br, eval_longitude= MagLon, eval_latitude= MagLat)
    model= (model_I[0] +model_T[0], model_I[1] +model_T[1], model_I[2] +model_T[2])
    model= model[0], *Geocentric_to_PlateCarree_vector_components(*model[1:], MagLat)
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
    MeridCur=axes[-2].quiver([0]*len(mlat), mlat, qd_east, qd_north, scale=1e-1, color='black', label='Meridian Current')
    return CurQ, MagQ, Mag, MeridCur, scatter

"""Add permanent features"""
#Add SEC poles
axes[2].scatter(node_lons, node_lats, transform =ccrs.PlateCarree(), color='red', label= 'SEC Poles', s= 0.5, zorder=50)
for ax in axes[:-2]:
    #Add cartopy features
    ax.spines['geo'].set_visible(False)
    #ax.outline_patch.set_visible(False)
    features(ax, coastlines_only=True)
for ax in axes[-2:]:
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
axes[-2].spines['bottom'].set_visible(False)

#Add magnetometers
for ax in axes[1:-2]:
    ax.scatter(MagLon, MagLat, transform =ccrs.PlateCarree(), marker= '*', color= 'darkorange', label= 'Magnetometer \n Station', zorder=100)
for i, maglat in enumerate(MagMlat):
    if MagMlon[i]-105<0:
        axes[-2].text(x= -0.055*1.1, y= maglat, s='-', color='orange', fontsize=30, zorder=-1)
    else:
        axes[-2].text(x= 0.055*1, y= maglat, s='-', color='orange', fontsize=30, zorder=-1)

try:
    #Add the IMAGE logo
    Image = plt.imread('https://space.fmi.fi/image/www/images/IMAGE_logo_2.png')
    Super = plt.imread('http://supermag.jhuapl.edu/lib/img/supermag-color-white-64.png')
    newax_Im = fig.add_axes([0.01, 0.38, 0.1, 0.1], anchor='NE', zorder=-1)
    newax_Su = fig.add_axes([0.01, 0.45, 0.1, 0.1], anchor='NE', zorder=-1)
    newax_Im.imshow(Image)
    newax_Im.axis('off')
    newax_Su.imshow(Super)
    newax_Su.axis('off')
except HTTPError:
    pass

#Subtitles
axes[1].set_title('Ionospheric Currents')
axes[2].set_title('Magnetic Field on Ground')
axes[0].set_title('SEC POLE Amplitudes')
axes[-2].set_title('Magnetic Meridian 105'+ u'\N{DEGREE SIGN}')
axes[-1].set_title('Validity')

#Axis labels and limits
axes[-2].set_ylabel('Magnetic Latitude')
axes[-1].set_xlabel('Data (nT)')
axes[-1].set_ylabel('Model (nT)')
axes[-2].set_ylim(min(mlat), max(mlat))
axes[-2].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
for ax in axes[:-2]:
    ax.set_xlim((-2067886.9433447798, 1588178.0724498471))
    ax.set_ylim((-2102947.263980966, 2610404.220997629))

#Colorbars
mappable1=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=0.001*1e3, vmax= 0.018*1e3), cmap=mpl.cm.afmhot_r)
cbar1=fig.colorbar(mappable=mappable1, ax=axes[1], orientation='horizontal')
cbar1.set_label(label='Ionospheric Current Magnitude '+r'$Akm^{-1}$', rotation=0)

mappable2=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-10.0E-9 *1e9, vmax=10.0E-9 *1e9), cmap=mpl.cm.bwr)
cbar2=fig.colorbar(mappable=mappable2, ax=axes[2], orientation='horizontal', extend='both')
cbar2.set_label(label='Evaluated Br (nT)', rotation=0)

mappable3=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-10E1, vmax=10E1), cmap=mpl.cm.seismic)
cbar3=fig.colorbar(mappable=mappable3, ax=axes[0], orientation='horizontal')
cbar3.ax.locator_params(nbins=7)
cbar3.set_label(label='SEC Pole Amplitude', rotation=0)

#Add Lat and MagLat circles
for Latitude in np.r_[50:85:5]:
    axes[0].plot(np.linspace(0, 45), [Latitude]*len(np.linspace(0, 45)), transform= ccrs.PlateCarree(), zorder=110, alpha=0.5, linestyle='--', color='blue')
    axes[0].text(np.linspace(0, 45)[0]-6, Latitude, transform= ccrs.PlateCarree(), zorder=110, alpha=1, s= str(Latitude)+ ' '+u'\N{DEGREE SIGN}', horizontalalignment='left', verticalalignment='center', fontsize=10)
for i, maglat in enumerate(np.r_[50:85:5]):
    Latitude, Longitude, err= A.apex2geo(maglat, np.linspace(80,130, 100), height=0)
    axes[1].plot(Longitude, Latitude, transform= ccrs.PlateCarree(), zorder=110, alpha=0.5, linestyle='--', color='blue')
    axes[1].text(Longitude[0]-6, Latitude[0], transform= ccrs.PlateCarree(), zorder=110, alpha=1, s= str(maglat)+ ' '+u'\N{DEGREE SIGN}', horizontalalignment='left', verticalalignment='center', fontsize=10)

#Add evaluation meridian
axes[1].plot(merid_glon, merid_glat, transform= ccrs.PlateCarree(), color ='blue', alpha=0.6, zorder=200)
axes[2].plot(merid_glon, merid_glat, transform= ccrs.PlateCarree(), color ='blue', alpha=0.6, zorder=200)
axes[-2].plot([0]*len(mlat), mlat)

#Add y=x line
axes[-1].plot(np.linspace(-1e3,1e3), np.linspace(-1e3, 1e3))

#Resize the figure
DPI=fig.get_dpi()
fig.set_size_inches(2400.0/float(DPI),1220.0/float(DPI))

plt.pause(1)
sort_index= np.argsort(node_lats)
node_lons2= node_lons[sort_index]
node_lats2= node_lats[sort_index]
"""Runnning the code"""
#Run time series
for j, date in enumerate(pd.date_range(start= start_Date, end= end_Date, freq='min')):
    mlt = A.mlon2mlt(105, date)
    plt.suptitle(str(date)+' UT \n'+str(round(mlt, 2)), fontsize= 15, bbox= {'color':'lightblue', 'pad': 10})
    CurQ, MagQ, Mag, MeridCur, scatter= PlottingFunction(date, poles)
    if j==0:
        #Add legend
        legend_x = -0.6
        legend_y = 1.014
        handle= []
        labels=[]
        for ax in axes[1:]:
            tmp_handles, tmp_labels =ax.get_legend_handles_labels()
            handle+= tmp_handles
            labels+= tmp_labels
        labels, ind= np.unique(labels, return_index=True)
        handles=[]
        for i in ind:
            handles.append(handle[i])
        leg=axes[0].legend(handles, labels, loc='upper left', bbox_to_anchor=(legend_x, legend_y))
        #Add quiver keys
        Curqk= axes[-2].quiverkey(CurQ, 0.45, 0.26, 200e-3, r'$200 \ Akm^{-1}$', labelpos='E',
                       coordinates='figure', zorder=500)
        Magqk= axes[2].quiverkey(MagQ, 0.669, 0.26, 100e-9, r'$100 \ nT$', labelpos='E',
                           coordinates='figure', zorder=500)
    
    plt.pause(.05)
    fig.set_size_inches(2000.0/float(DPI),1000.0/float(DPI))
    break
    plt.savefig('/home/simon/Documents/EGU_2021/Slide_Show/Event_pngs/'+str(j)+'.png')

    amplitudes= poles.Amplitude(Btheta, Bphi, Br)[sort_index]
    columns=['Date_UTC']
    for l in range(len(node_lats2)):
        columns.append(str((node_lons[l], node_lats[l])))
    CurQ.remove()
    MagQ.remove()
    Mag.remove()
    MeridCur.remove()
    for p in scatter:
        p.remove()        