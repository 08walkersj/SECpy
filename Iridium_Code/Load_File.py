#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 22:57:36 2021

@author: simon
"""

from netCDF4 import Dataset
from astropy import coordinates as coord
from astropy import units
import numpy as np
import pandas as pd
import os
from progressbar import progressbar
from SECpy import RE
hdf_path='/home/simon/BCSS-DAG Dropbox/Simon Walker/Iridium/All_Data.hdf5'
path= '/home/simon/BCSS-DAG Dropbox/Simon Walker/Iridium/'
1/0
max_value= len(os.listdir(path))
count=0
count2=0
for folder in progressbar(os.listdir(path),max_value=max_value):
    if folder.endswith('.hdf5'):
        continue
    print('\n'+folder)
    if folder!='2009' and count<1:
        continue
    else:
        count+=1
    max_value= len(os.listdir(path+folder))
    for file in progressbar(os.listdir(path+folder), max_value=max_value):
        if not file.endswith('.ncdf'): 
            continue
        if file!='20091121Amp_invert.ncdf' and count2<1:
             continue
        # elif count2==0:
        #     count2+=1
        else:
            count2+=1
        data= Dataset(path+folder+'/'+file)
        dimensions= ['datapnt', 'vec_comp', 'plane', 'sv']
        headers= ['time', 'pseudo_sv_num', 'plane_num', 'pos_eci', 'b_eci', 'b_error', 
                  'pseudo_sv_quality', 'data_splice', 'sv_matrix']
        year=file[:4]
        month=file[4:6]
        day= file[6:8]
        time= np.datetime64(f'{year}-{month}-{day}')
        time+=(data['time'][...].data*60*60*1e9).astype('timedelta64[ns]')
        
        cart_pos = coord.CartesianRepresentation(data['pos_eci'][...].data.T[0],
                                                 data['pos_eci'][...].data.T[1],
                                                 data['pos_eci'][...].data.T[2],
                                                 unit = units.m)
        
        cart_B   = coord.CartesianRepresentation(data['b_eci'][...].data.T[0], 
                                                 data['b_eci'][...].data.T[1], 
                                                 data['b_eci'][...].data.T[2],
                                                 unit=units.nT)
        
        gcrs_pos = coord.GCRS(cart_pos, obstime=time)
        itrs_pos = gcrs_pos.transform_to(coord.ITRS(obstime=time))
        
        gcrs_B = coord.GCRS(cart_B, obstime=time)
        itrs_B = gcrs_B.transform_to(coord.ITRS(obstime=time))
        
        
        lat = itrs_pos.spherical.lat.value 
        lon = itrs_pos.spherical.lon.value
        r   = itrs_pos.spherical.distance.value -RE
        
        # convert to enu
        e = np.vstack((-                    np.sin(np.deg2rad(lon)),                      np.cos(np.deg2rad(lon)), np.zeros_like(lon)))
        n = np.vstack((-np.sin(np.deg2rad(lat)) * np.cos(np.deg2rad(lon)), -np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(lon)), np.cos(np.deg2rad(lat) )))
        u = np.vstack(( np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(lon)),  np.cos(np.deg2rad(lat)) * np.sin(np.deg2rad(lon)), np.sin(np.deg2rad(lat) )))
        
        Bx_ecef, By_ecef, Bz_ecef = itrs_B.x.value, itrs_B.y.value, itrs_B.z.value
        Be = e[0] * Bx_ecef + e[1] * By_ecef + e[2] * Bz_ecef
        Bn = n[0] * Bx_ecef + n[1] * By_ecef + n[2] * Bz_ecef
        Bu = u[0] * Bx_ecef + u[1] * By_ecef + u[2] * Bz_ecef
        try:
            b_error= data['b_error'][...].data
        except IndexError:
            b_error= [np.nan]*len(Be)
        pd.DataFrame({'Date_UTC':time, 'glon': lon, 'glat':lat, 'altitude': r, 'Be':Be, 'Bn':Bn, 'Bu':Bn,
                      'B_error':b_error}).to_hdf(hdf_path, key='main', 
                                                                    mode='a',append=True,format='t', 
                                                                    data_columns=True)
                            