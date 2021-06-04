#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 14:04:31 2020
@author: Simon Walker

"""
# %% Imports
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
from SECpy import RE, SECS
from SECpy import PlottingTools as PT
Geocentric_to_PlateCarree_vector_components, features= PT.Geocentric_to_PlateCarree_vector_components, PT.features
from apexpy import Apex
import os
import datetime as dt
import cubedsphere as CS
import matplotlib.pyplot as plt
import matplotlib as mpl
import shutil
from urllib.error import HTTPError
plt.ioff()

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
gs = fig.add_gridspec(100, 4, width_ratios= [2,2,2,1])
lon_centre= 17.7
lat_centre= 68.1
proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
# fig.add_subplot(gs[15:70, 3]) #old dimensions
axes = [fig.add_subplot(gs[:,0], projection= proj), fig.add_subplot(gs[:,1], projection= proj), 
                 fig.add_subplot(gs[:,2], projection= proj),fig.add_subplot(gs[5:70, 3]), fig.add_subplot(gs[75:90, 3], aspect='auto')]
# %% Set date
"""Check Date and get Mag locations"""
# Date= pd.Timestamp('2000-02-05T22:33')
# Date= pd.Timestamp('2000-10-03T02:40')
# Date= pd.Timestamp('2003-07-16T02:00')
# Date= pd.Timestamp('2003-07-19T02:14')
# Date= pd.Timestamp('2002-02-10T02:31')
# Date=pd.Timestamp('2003-01-28T02:31')
# Date=pd.Timestamp('2001-10-14T02:53')
# Date=pd.Timestamp('2001-10-20T02:37')
# Date=pd.Timestamp('2000-10-10T02:41')
# Date=pd.Timestamp('2002-08-06T02:01')
# Date=pd.Timestamp('2000-06-23T02:18')
# Date=pd.Timestamp('2002-06-06T02:04')
# Date= pd.Timestamp('2002-06-20T02:23')
# Date= pd.Timestamp('2018-04-22T13:00')
# Date= pd.Timestamp('2000-02-05T22:49')

Date= pd.Timestamp('2000-09-21T19:16')
# Date= pd.Timestamp('2000-05-23T20:14')
# subs= pd.read_pickle('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/merged_substormlist.pd')
# subs=subs[(subs['mlon']<120)&(subs['mlon']>95)&(subs.index>'2000')&(subs.AL_peak<-100)&(subs.mlat>0)]
# plots_list=[]
# DataFrame= pd.read_hdf(os.getcwd()+"/Data/Scandinavia/Scandinavia_HDF.h5", mode='r', where="Date_UTC=Date" )
# Extra= pd.read_csv('Data/Extra_converted.csv')
# Extra.Date_UTC= Extra.Date_UTC.astype('datetime64[ns]')
# sites=DataFrame.Site.values
# index= DataFrame['Site'].isin(sites)
# MagLat=DataFrame['glat'][index].values
# MagLon=DataFrame['glon'][index].values
# MagLat= np.array([64.52, 68.35, 54.61, 55.62, 69.3 , 68.02, 64.94, 77.  , 60.5 ,
#         69.76, 78.92, 66.9 , 74.5 , 62.3 , 69.02, 69.66, 78.2 , 67.37,
#         70.54, 59.9 ])
# MagLon= np.array([27.23, 18.82, 18.82, 11.67, 16.03, 23.53, 10.99, 15.6 , 24.65,
#         27.01, 11.95, 24.08, 19.2 , 26.65, 20.79, 18.94, 15.83, 26.63,
#         22.22, 17.35])
# MagLat=np.append(MagLat,[Extra.glat[Extra.Date_UTC==Date].values]*1)
# MagLon=np.append(MagLon,[Extra.glon[Extra.Date_UTC==Date].values]*1)
# MagLon=np.append(MagLon, [-20.2, -22,-8.7,-18.6]*10)

if len(MagLon)!=len(sites):
    raise ValueError("not all magnetometers available!")
# savedir='/Home/siv32/zef014/Documents/New_Examples/'+str(Date.date())+"/"
# savedir='/Home/siv32/zef014/Documents/temperature_spikes_new/'
savedir= '/home/simon/Documents/EGU_2021/Slide_Show/Event_pngs/'
data_path= './SECS_Visualisation/SECS_Data/Scandinavia_HDF.h5' #input data path
#     else:
#         sys.exit()
# %% Set grid
"""Grid"""
A = Apex(date=dt.datetime(2008, 6, 1, 0, 0, 0))
f1, f2 = A.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
qd_north = f2 / np.linalg.norm(f2)
East, North= qd_north[0], qd_north[1]
Gridproj= CS.CSprojection((lon_centre, lat_centre), [East, North])
# node_grid=CS.CSgrid(Gridproj, 4200, 3500, 50., 50.)
node_grid=CS.CSgrid(Gridproj, 3700, 2200, 50., 50.)
Evalproj= CS.CSprojection((lon_centre-0.23, lat_centre+0.23), [East, North])
eval_grid= CS.CSgrid(Evalproj,  3300, 1800, 50., 50.)
node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
lat,lon= eval_grid.lat.flatten(), eval_grid.lon.flatten()
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
#poles.Fitting_Matrix(MagLon, MagLat, cond=0.1**2)
matrix= GTG + λ1*np.identity(GTG.shape[0]) +λ2*np.dot(L.T, L)/np.max(np.abs(np.dot(L.T, L)))
Inverse= poles.Inverse_Matrix(matrix, cond=0)
poles.fitting_matrix= np.dot(Inverse,G.T)
# %% Meridian
"""Meridian"""
mlat= np.linspace(49, 81, 50)
merid_glat, merid_glon, err= A.apex2geo(alat=mlat, alon= 105, height=0)
f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3= A.basevectors_apex(merid_glat, merid_glon, 0, coords='geo')
proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
MagMlat, MagMlon= A.geo2apex(MagLat, MagLon, 0)
# %% Define Plotting function
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
    DataFrame=DataFrame.sort_values(by='glat')
    index= DataFrame['Site'].isin(sites)
    Btheta, Bphi, Br= -DataFrame['Btheta'][index].values*1e-9,DataFrame['Bphi'][index].values*1e-9,-DataFrame['Br'][index].values*1e-9
    # Btheta=np.append(Btheta,[-Extra.Btheta[Extra.Date_UTC==date].values*1e-9]*1)
    # Bphi=np.append(Bphi,[Extra.Bphi[Extra.Date_UTC==date].values*1e-9]*1)
    # Br=np.append(Br,[-Extra.Br[Extra.Date_UTC==date].values*1e-9]*1)
    
    if len(Btheta)!= len(sites):
        raise ValueError("not all magnetometers available: "+str(date))
    scatter=[]
    #Currents
    CI= poles.Currents(Btheta, Bphi, Br)
    CMag=np.sqrt((CI[0]**2) +CI[1]**2)
    CurE, CurN= Geocentric_to_PlateCarree_vector_components(*CI,lat)
    CurQ=axes[1].quiver(lon[::3], lat[::3], CurE[::3], CurN[::3], zorder=90,scale=(1e1)/2, transform = ccrs.PlateCarree(), color='Black', label='Ionospheric Current', alpha=0.7)
    scatter.append(axes[1].scatter(lon, lat, c= CMag*1e3, vmin=0.012*1e3, vmax= .5*1e3, cmap= mpl.cm.afmhot_r, transform= ccrs.PlateCarree(), zorder=50))
    #Magnetic field
    MI, MT = poles.Magnetic_Field(Btheta, Bphi, Br)
    M= (MI[0] +MT[0], MI[1] +MT[1], MI[2] +MT[2])
    MagE, MagN= Geocentric_to_PlateCarree_vector_components(*M[1:],lat)
    Mag=axes[2].quiver(lon[::8], lat[::8], MagE[::8], MagN[::8], zorder=80,
                        scale=(1e-5)/8, transform = ccrs.PlateCarree(), color='Black', label='Total Field', alpha=0.7)
    Measured= Geocentric_to_PlateCarree_vector_components(Bphi, -Btheta, MagLat)
    MagQ=axes[2].quiver(MagLon, MagLat,  *Geocentric_to_PlateCarree_vector_components(Bphi, -Btheta, MagLat),zorder=90, 
                         scale = (1e-5)/8, transform = ccrs.PlateCarree(), color='red', label='Measurements', alpha=0.7)
    scatter.append(axes[2].scatter(lon, lat, c= M[0]*1E9, vmin=-4.0E2, vmax=4.0E2,
                                    cmap=mpl.cm.bwr, transform= ccrs.PlateCarree(), zorder=50))
    #SEC pole amplitudes
    scatter.append(axes[0].scatter(node_lons, node_lats, c= poles.Amplitude(Btheta, Bphi, Br), vmin=-2E3, vmax=2E3, 
                                   cmap=mpl.cm.seismic, transform= ccrs.PlateCarree(), zorder=50, s=10))
    #Model vs Measured
    model_I, model_T= poles.Magnetic_Field(Btheta, Bphi, Br, eval_longitude= MagLon, eval_latitude= MagLat)
    model= (model_I[0] +model_T[0], model_I[1] +model_T[1], model_I[2] +model_T[2])
    model= model[0], *Geocentric_to_PlateCarree_vector_components(*model[1:], MagLat)
    scatter.append(axes[-1].scatter(Br*1e9, model[0]*1e9, label='Br component', color= 'purple'))
    scatter.append(axes[-1].scatter(Measured[0]*1e9, model[1]*1e9, label='Be component', color='green'))
    scatter.append(axes[-1].scatter(Measured[1]*1e9, model[-1]*1e9, label='Bn component', color='red'))
    
    # scatter.append(axes[-1].scatter(Br[33:]*1e9, model[0][33:]*1e9, label='Br component', color= 'purple', marker='x'))
    # scatter.append(axes[-1].scatter(Measured[0][33:]*1e9, model[1][33:]*1e9, label='Be component', color='green', marker='x'))
    # scatter.append(axes[-1].scatter(Measured[1][33:]*1e9, model[-1][33:]*1e9, label='Bn component', color='red', marker='x'))
    lim=max(abs(np.append(Measured[0]*1e9, Measured[1]*1e9)))
    axes[-1].set_xlim((-lim*1.2, lim*1.2))
    axes[-1].set_ylim((-lim*1.2, lim*1.2))
    #Meridian Current
    # start=time.time()
    meridian= poles.Currents(Btheta, Bphi, Br, eval_longitude= merid_glon, eval_latitude= merid_glat, singularity_limit= 50e3)
    poles.Magnetic_Field(Btheta, Bphi, Br, eval_longitude= merid_glon, eval_latitude= merid_glat)
    # end=time.time()
    # print(end-start)
    qd_east= g1[0]*meridian[0] +g1[1]*meridian[1]
    qd_north= g2[0]*meridian[0] +g2[1]*meridian[1]
    MeridCur=axes[-2].quiver([0]*len(mlat), mlat, qd_east, qd_north, scale=1., color='black', label='Meridian Current')
    return CurQ, MagQ, Mag, MeridCur, scatter
# %% Add permanent features
"""Add permanent features"""
#Add SEC poles
axes[2].scatter(node_lons, node_lats, transform =ccrs.PlateCarree(), color='red', label= 'SEC Poles', s= 0.5, zorder=50)
for ax in axes[:-2]:
    #Add carto features
    ax.spines['geo'].set_visible(False)
    # ax.outline_patch.set_visible(False)
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
axes[-2].set_ylabel('Magnetic Latitude'+u' (\N{DEGREE SIGN})')
axes[-1].set_xlabel('Data (nT)')
axes[-1].set_ylabel('Model (nT)')
axes[-2].set_ylim(min(mlat), max(mlat))
axes[-2].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
for ax in axes[:-2]:
    ax.set_xlim((-2067886.9433447798, 1588178.0724498471))
    ax.set_ylim((-2102947.263980966, 2610404.220997629))
#Colorbars
mappable1=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=0.012*1e3, vmax= .5*1e3), cmap=mpl.cm.afmhot_r)
cbar1=fig.colorbar(mappable=mappable1, ax=axes[1], orientation='horizontal')
cbar1.set_label(label='Ionospheric Current Magnitude '+r'$Akm^{-1}$', rotation=0)

mappable2=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-4.0E2, vmax=4.0E2,), cmap=mpl.cm.bwr)
cbar2=fig.colorbar(mappable=mappable2, ax=axes[2], orientation='horizontal', extend='both')
cbar2.set_label(label='Evaluated Br (nT)', rotation=0)

mappable3=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-2E3, vmax=2E3), cmap=mpl.cm.seismic)
cbar3=fig.colorbar(mappable=mappable3, ax=axes[0], orientation='horizontal')
cbar3.ax.locator_params(nbins=7)
cbar3.set_label(label='SEC Pole Amplitude (A)', rotation=0)
#Add Lat and MagLat circles
for Latitude in np.r_[50:85:5]:
    axes[0].plot(np.linspace(0, 45), [Latitude]*len(np.linspace(0, 45)), transform= ccrs.Geodetic(), zorder=110, alpha=0.5, linestyle='--', color='blue')
    axes[0].text(np.linspace(0, 45)[0]-6, Latitude, transform= ccrs.Geodetic(), zorder=110, alpha=1, s= str(Latitude)+ ' '+u'\N{DEGREE SIGN}', horizontalalignment='left', verticalalignment='center', fontsize=10)
for i, maglat in enumerate(np.r_[50:85:5]):
    Latitude, Longitude, err= A.apex2geo(maglat, np.linspace(80,130, 100), height=0)
    axes[1].plot(Longitude, Latitude, transform= ccrs.Geodetic(), zorder=110, alpha=0.5, linestyle='--', color='blue')
    axes[1].text(Longitude[0]-6, Latitude[0], transform= ccrs.Geodetic(), zorder=110, alpha=1, s= str(maglat)+ ' '+u'\N{DEGREE SIGN}', horizontalalignment='left', verticalalignment='center', fontsize=10)
#Add evaluation meridian
axes[1].plot(merid_glon, merid_glat, transform= ccrs.Geodetic(), color ='blue', alpha=0.6, zorder=200)
axes[2].plot(merid_glon, merid_glat, transform= ccrs.Geodetic(), color ='blue', alpha=0.6, zorder=200)
axes[-2].plot([0]*len(mlat), mlat)
#Add y=x line
axes[-1].plot(np.linspace(-1e3,1e3), np.linspace(-1e3, 1e3))
#Resize the figure
DPI=fig.get_dpi()
fig.set_size_inches(2400.0/float(DPI),1220.0/float(DPI))

# manager = plt.get_current_fig_manager()
# manager.window.showMaximized()
# plt.pause(1)
#Run time series
# 1/0
sort_index= np.argsort(node_lats)
node_lons2= node_lons[sort_index]
node_lats2= node_lats[sort_index]
# np.savetxt(savedir+'_node_lons.txt', node_lons2)
# np.savetxt(savedir+'_node_lats.txt', node_lats2)
# breakpoint()
# %% Run loop
# for i, Date, mlon in zip(range(len(subs)),subs.index,subs.mlon):
#     try:
#         DataFrame= pd.read_hdf(os.getcwd()+"/Data/Scandinavia/Scandinavia_HDF.h5", mode='r', where="Date_UTC=Date.round('min')" )
#         # Extra= pd.read_csv('Data/Extra_converted.csv')
#         # Extra.Date_UTC= Extra.Date_UTC.astype('datetime64[ns]')
#         # sites=DataFrame.Site.values
#         index= DataFrame['Site'].isin(sites)
#         MagLat=DataFrame['glat'][index].values
#         MagLon=DataFrame['glon'][index].values
#         # MagLat= np.array([64.52, 68.35, 54.61, 55.62, 69.3 , 68.02, 64.94, 77.  , 60.5 ,
#         #         69.76, 78.92, 66.9 , 74.5 , 62.3 , 69.02, 69.66, 78.2 , 67.37,
#         #         70.54, 59.9 ])
#         # MagLon= np.array([27.23, 18.82, 18.82, 11.67, 16.03, 23.53, 10.99, 15.6 , 24.65,
#         #         27.01, 11.95, 24.08, 19.2 , 26.65, 20.79, 18.94, 15.83, 26.63,
#         #         22.22, 17.35])
#         # MagLat=np.append(MagLat,[Extra.glat[Extra.Date_UTC==Date].values]*1)
#         # MagLon=np.append(MagLon,[Extra.glon[Extra.Date_UTC==Date].values]*1)
#         # MagLon=np.append(MagLon, [-20.2, -22,-8.7,-18.6]*10)
        
#         if len(MagLon)!=len(sites):
#             raise ValueError("not all magnetometers available!")
#     except ValueError:
#         continue
#     savedir= f'/Home/siv32/zef014/Documents/Masters_Presentation/Events/Event{i+1}/'
#     os.makedirs(savedir)

# savedir=f'/home/simon/Documents/EGU_2021/Slide_Show/Event_pngs/{Date.date()}/'
# try:
#     os.makedirs(savedir)
# except FileExistsError:
#     if input('Will overwrite existing directory, continue? (y/n)').lower()!='n':        
#         shutil.rmtree(savedir)
#         os.makedirs(savedir)
# os.makedirs(savedir)
onset= pd.Timestamp('2000-09-21T19:16:32')
# for j, date in enumerate(pd.date_range(start= Date.round('min')- pd.Timedelta(value=30, unit='m'), end= Date.round('min')+ pd.Timedelta(value=30, unit='m'), freq='min')):
for j, date in enumerate(pd.date_range(start= onset.round('min')- pd.Timedelta(value=15, unit='m'), end= Date.round('min')+ pd.Timedelta(value=35, unit='m'), freq='min')):
    diff=np.timedelta64(date-onset, 's').astype('float64')/60
    # date= pd.Timestamp('2000-02-05T22:49')
    # date= pd.Timestamp('2000-02-05T22:34')
    mlt = A.mlon2mlt(105, date)
    if diff<0:
        sign='-'
    else:
        sign='+'
    plt.suptitle(f"Time: {date}\nMLT: {round(mlt, 2)}\nOnset{sign}{abs(round(diff, 2))} (min)",fontsize= 15, bbox= {'color':'lightblue', 'pad': 10})
    # plt.suptitle(str(date)+' UT \n'+str(round(mlt, 2))+' MLT \n'+ 'Onset: '+str(Date) +f',{round(mlon)}'+r'$^\circ$', fontsize= 15, bbox= {'color':'lightblue', 'pad': 10})
    # plt.suptitle(str(date)+' UT \n'+str(round(mlt, 2))+' MLT \n'+ 'Depth: '+str(depth), fontsize= 15, bbox= {'color':'lightblue', 'pad': 10})
    # plt.suptitle('Solved with Ionospheric SEC poles, Image Currents\nand Regularisation', fontsize= 15, y=0.95, bbox= {'color':'lightblue', 'pad': 10})
    try:
        CurQ, MagQ, Mag, MeridCur, scatter= PlottingFunction(date, poles)
    except ValueError:
        continue
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
    Curqk= axes[-2].quiverkey(CurQ, 0.47, 0.26, 2000e-3, r'$2000 \ Akm^{-1}$', labelpos='E',
                   coordinates='figure', zorder=500)
    Magqk= axes[2].quiverkey(MagQ, 0.7, 0.26, 400e-9, r'$400 \ nT$', labelpos='E',
                       coordinates='figure', zorder=500)
    
    # plt.pause(.05)
    fig.set_size_inches(2400.0/float(DPI),1220.0/float(DPI))
    break
    # plt.savefig(savedir+f"{j:02d}"+'.png')

    # amplitudes= poles.Amplitude(Btheta, Bphi, Br)[sort_index]
    columns=['Date_UTC']
    for l in range(len(node_lats2)):
        columns.append(str((node_lons[l], node_lats[l])))
    # dframe= pd.DataFrame(columns=columns)
    # dframe[columns[1:]]= amplitudes
    # dframe['Date_UTC']= date
    # dframe.to_hdf(savedir+'Amplitude.hdf5', mode='a', append=True,format='t', data_columns=True, key='main')
    # np.savetxt(savedir+str(date)+'.txt', amplitudes)
    
    # if j==2:
    # break
    Curqk.remove()
    Magqk.remove()
    CurQ.remove()
    MagQ.remove()
    Mag.remove()
    MeridCur.remove()
    for p in scatter:
        p.remove()        
        
    # import animation
    # animation.saveVideo(savedir, remove_message=True)