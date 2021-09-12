#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:38:32 2020
@author: Simon Walker

"""
#@author: zef014
#%% Import packages
import pandas as pd
import numpy as np
import pyamps
import os
from SECpy import RE, SECS
from apexpy import Apex
import datetime as dt
from secsy.secsy import cubedsphere as CS
import time
from progressbar import progressbar
from Interpolation2d import interpolation2d as interp
import vaex as vx
# import matplotlib.pyplot as plt
import gc
"""CREATING MATRICES"""
def Ultimate_Matrix(Matrices, fitting_matrix):
    #Combines Matrices into one
    return np.concatenate(Matrices)@fitting_matrix
# 1/0
#%% Set Up Matrices
"""Grid"""
#Creating the grid using the cubed sphere method
lon_centre= 17.7
lat_centre= 68.1
A = Apex(date=dt.datetime(2008, 6, 1, 0, 0, 0))
f1, f2 = A.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
qd_north = f2 / np.linalg.norm(f2)
East, North= qd_north[0], qd_north[1]
Gridproj= CS.CSprojection((lon_centre, lat_centre), [East, North])
node_grid=CS.CSgrid(Gridproj, 3700, 2200, 50., 50.)
node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
MagLat= np.array([54.61, 55.62, 59.9 , 60.5 , 62.3 , 64.52, 64.94, 66.9 , 67.37,
        68.02, 68.35, 69.02, 69.3 , 69.66, 69.76, 70.54, 74.5 , 77.  ,
        78.2 , 78.92])
MagLon= np.array([18.82, 11.67, 17.35, 24.65, 26.65, 27.23, 10.99, 24.08, 26.63,
        23.53, 18.82, 20.79, 16.03, 18.94, 27.01, 22.22, 19.2 , 15.6 ,
        15.83, 11.95])

"""Meridian"""
mlat= np.linspace(49, 81, 50)
depth=500
merid_glat, merid_glon, err= A.apex2geo(alat=mlat, alon= 105, height=0)
poles= SECS(node_lons, node_lats, merid_glon, merid_glat, mode='image', image_current_radius=RE-depth*1E3)
f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3= A.basevectors_apex(merid_glat, merid_glon, 0, coords='geo')
"""Regularisiation"""
λ1= 1e-23
λ2= λ1*1e2
Le, Ln=node_grid.get_Le_Ln()
node_f1, node_f2= A.basevectors_qd(node_grid.lat.flatten(), node_grid.lon.flatten(), 110, coords='geo')
e= node_f1/np.linalg.norm(node_f1, axis=0)
L= np.diag(e[0]).dot(Le.toarray()) + np.diag(e[1]).dot(Ln.toarray())
G=poles.G_Matrix(MagLon, MagLat)
GTG= np.dot(G.T, G)
matrix= GTG + λ1*np.identity(GTG.shape[0]) +λ2*np.dot(L.T, L)/np.max(np.abs(np.dot(L.T, L)))
Inverse= poles.Inverse_Matrix(matrix, cond=0)
poles.fitting_matrix['divergence-free']= np.dot(Inverse,G.T)

np.savetxt('/media/simon/Seagate Expansion Drive/SECS/Amplitude/Fitting_Matrix.txt', poles.fitting_matrix['divergence-free'])
#%% Run Evaluation
# 1/0
lon_centre= 17.7
lat_centre= 68.1
A = Apex(date=dt.datetime(2008, 6, 1, 0, 0, 0))
f1, f2 = A.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
qd_north = f2 / np.linalg.norm(f2)
East, North= qd_north[0], qd_north[1]
Gridproj= CS.CSprojection((lon_centre, lat_centre), [East, North])
mlat= np.linspace(49, 81, 50)
# mlat, mlon=A.geo2apex(node_lats, node_lons, 110)
node_grid=CS.CSgrid(Gridproj, 3700, 2200, 50., 50.)
# x, y = node_grid.xi, node_grid.eta
# x, y= np.unique(x), np.unique(y)
node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
sites=np.array(['LYR', 'BJN', 'TRO', 'HRN', 'ABK', 'SOD', 'NUR', 'BFE', 'NAL',
        'SOR', 'AND', 'KEV', 'KIL', 'MUO', 'PEL', 'OUJ', 'HAN', 'UPS',
        'HLP', 'RVK'], dtype='<U3')
fitting_matrix= np.loadtxt('/media/simon/Seagate Expansion Drive/SECS/Amplitude/Fitting_Matrix.txt')
start_time=time.time()
epoch = 2010.
# columns=['Date_UTC', 'MLT']
# columns+=[f'Amplitude{i}' for i in range(len(mlat))]
date_range= pd.date_range(start=  pd.Timestamp('1999-01-01T00:00:00'), end= pd.Timestamp('2020-01-01T00:00:00'), freq='Y')
from multiprocessing.pool import Pool
p=Pool(10)
1/0
def func(dates):
    for k, date in progressbar(enumerate(date_range), max_value= len(dates), prefix='Looping through yearly intervals |'):    
        start= date
        end= date+ pd.DateOffset(years=1)
        DataFrame= pd.read_hdf('/media/simon/Seagate Expansion Drive/SECS_Visualisation/SECS_Data/Scandinavia_HDF.h5', where='Date_UTC>=start & Date_UTC<end')
        DataFrame=DataFrame.sort_values(by='glat')
        index= DataFrame['Site'].isin(sites)
        DataFrame=DataFrame[index]
        dates=[]
        i=0
        for date2 in progressbar(pd.date_range(start=  start, end= end, freq='min'), 
                                 max_value=len(pd.date_range(start=  start, end= end, freq='min')), 
                                 prefix='Collecting data from each minute |', suffix= f'| Current Range: {start.date()}- {end.date()}'):
            index=DataFrame.Date_UTC==date2
            if len(DataFrame[index].Site)!=len(sites):
                continue
            else:
                dates.append(date2)
                i+=1
                Btheta, Bphi, Br= -DataFrame[index].Btheta.values*1e-9,DataFrame[index].Bphi.values*1e-9,-DataFrame[index].Br.values*1e-9
                if i==1:
                    Data= np.concatenate((Btheta, Bphi, Br))
                else:
                    Data= np.c_[Data, np.concatenate((Btheta, Bphi, Br))]
        if len(dates)==0:
            continue
        else:
            #Interpolation
            # tmp_df=pd.DataFrame(columns=np.append(columns,['PC_N_INDEX', 'AL_INDEX', 'BX_GSM', 'BY_GSM', 'BZ_GSM']))
            # mlts = pyamps.mlt_utils.mlon_to_mlt(105., dates, epoch)
            # if Output.T.shape== (len(mlat)*3,):
            #     tmp_df[columns[2:]]= [Output.T]
            # else:
            #     tmp_df[columns[2:]]=pd.DataFrame(Output.T)
            # tmp_df['Date_UTC']= dates
            # tmp_df['MLT']= mlts
            # tmp_df[columns[2:]]= 
            # Omni= pd.read_hdf('/media/simon/Seagate Expansion Drive/Solar_Wind_Data/OMNI_Data_Pandas_shifted.hdf5',
            #                   where= 'Date_UTC>=start & Date_UTC<end', columns=['Date_UTC','PC_N_INDEX', 'AL_INDEX', 'BX_GSE', 'BY_GSM', 'BZ_GSM'])
            # index= Omni.Date_UTC.isin(tmp_df.Date_UTC.values)
            # tmp_df[['PC_N_INDEX', 'AL_INDEX', 'BX_GSM', 'BY_GSM', 'BZ_GSM']]= pd.DataFrame(Omni.loc[index, ['PC_N_INDEX', 'AL_INDEX', 'BX_GSE', 'BY_GSM', 'BZ_GSM']].values)
            # tmp_df.to_hdf('/media/simon/Seagate Expansion Drive/SECS/Magnetic_Only/results_singularity_mod.hdf5','main',mode='a',append=True,format='t', data_columns=True)
            #Full grid
            amplitude= np.dot(fitting_matrix, Data).T.flatten()
            Date_UTC=np.array([np.array(dates).astype('datetime64[m]').astype(float)]*len(node_lons)).flatten()
            vx.from_arrays(Date_UTC= Date_UTC, longitude=np.concatenate([node_lons]*len(dates)), 
                           latitude=np.concatenate([node_lats]*len(dates)), 
                           amplitude= amplitude).export_hdf5(f'/media/simon/Seagate Expansion Drive/SECS/Amplitude/Amplitude_batch{date.year}.hdf5')
            open('/media/simon/Seagate Expansion Drive/SECS/Magnetic_Only/last_save.txt', 'w').write(str(date))
        gc.collect()
end_time=time.time()
Pool.map(np.split(date_range[:-1], 10)+[date_range[-1]], func)

