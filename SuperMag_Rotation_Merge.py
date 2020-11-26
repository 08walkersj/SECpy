#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 10:55:15 2020
@author: Simon Walker

"""
#@author: zef014
import pandas as pd
import numpy as np
import glob
import sys
def SuperMag2Geo(N, E, D):
    """
    Applies a rotation to the data points using a declination angle
    
    Parameters
    ----------
    N : numpy.ndarray
        Magnetometer measurement in the north direction.
    E : numpy.ndarray
        Magnetometer measurement in the east direction.
    D : numpy.ndarray
        Declination of each data point.

    Returns
    -------
    N2 : numpy.ndarray
        New north vector in the direction of magnetic north.
    E2 : numpy.ndarray
        New east vector in the direction of magnetic east.

    """
    E2= E*np.cos(D) +N*np.sin(D)
    N2= N*np.cos(D)- E*np.sin(D)
    return N2, E2
def declination(glat, glon, radius, year):
    """
    Calculates the Declination angle of the each data point from the CHAOS model
    
    Parameters
    ----------
    glat : numpy.ndarray (degrees)
        geographic lattitude.
    glon : numpy.ndarray (degrees)
        geographic longitude.
    radius : float, optional
        the radius from the centre of the earth to the points.
    year : int
        year.
    Returns
    -------
    numpy.ndarray (radians)
        Declination used to rotate the magnetometer vectors.

    """
    pysymmetry_path = '/Home/siv32/zef014/'
    if pysymmetry_path not in sys.path:
        sys.path.append(pysymmetry_path)
    sys.path.insert(1, '/Home/siv32/zef014/pysymmetry/src/geodesy.py')
    sys.path.append('/Home/siv32/zef014/pysymmetry/src/')
    import geodesy
    colat, radius, _, _ = geodesy.geod2geoc(glat, 0, 0, 0)
    import numpy as np
    from chaosmagpy import load_CHAOS_matfile
    from chaosmagpy.data_utils import mjd2000
    time = mjd2000(year, 1, 1)  # modified Julian date
    # load the CHAOS model
    model = load_CHAOS_matfile(pysymmetry_path + '/pysymmetry/models/' + 'CHAOS-7.mat') 
    B_radius_core, B_theta_core, B_phi_core = model.synth_values_tdep(time, radius, colat, glon)
    B_radius_crust, B_theta_crust, B_phi_crust = model.synth_values_static(radius, colat, glon)
    Br, Bth, Beast = B_radius_core + B_radius_crust, B_theta_core + B_theta_crust, B_phi_core + B_phi_crust
    # convert from gecentric to geodetic:
    glat, h, Bnorth, B_down = geodesy.geoc2geod(colat, radius, Bth, Br)
    return np.arctan(Beast/Bnorth)
def chaos_dec(glat, glon, year, radius= 6371.2e3):
    """

    Parameters
    ----------
    glat : numpy.ndarray (degrees)
        geographic lattitude.
    glon : numpy.ndarray (degrees)
        geographic longitude.
    year : int
        year.
    radius : float, optional
        the radius from the centre of the earth to the points. The default is 6371.2e3 (Radius of Earth).

    Returns
    -------
    numpy.ndarray (radians)
        Declination used to rotate the magnetometer vectors.

    """
    return declination(glat, glon, radius, year)

def SM_Conversion(in_folder, out_folder=False):
    """
    Rotate the magnetometer data from SuperMAG local coordinates to magnetic coordinates using the CHAOS model
    and merge the data sets into a HDF file

    Parameters
    ----------
    in_folder : str
        path to folder containing csv files that contain superMAG data. Ensure that all csv files within the folder contain superMAG data!
    out_folder : str, optional
        Specify a new path for the HDF file. The default is False which saves the HDF file within the same folder as the data.

    Returns
    -------
    None.

    """
    if not out_folder:
        out_folder= in_folder
    if not in_folder.endswith('/'):
        in_folder+='/'
    if not out_folder.endswith('/'):
        out_folder+='/'
    files= glob(in_folder+'.*csv')
    details=pd.read_csv("Super_Mag_stations.csv")
    columns= ['Date_UTC', 'glon', 'glat', 'Bphi', 'Btheta', 'Br', 'Declination', 'Site']
    pd.DataFrame(columns=columns).to_hdf(out_folder+'Data.hdf5', key='Main', format='table', data_columns=True)
    HDFStore= pd.HDFStore(out_folder+'Data.hdf5',)
    for file in files:
        print(f'loading... \n{file}')
        data= pd.read_csv(file)
        year=pd.to_datetime(data.Date_UTC)[0].year
        temp_df= pd.DataFrame(columns=columns)
        temp_df['Date_UTC']= pd.to_datetime(data.Date_UTC)
        temp_df['Site']=  data.IAGA.values
        data_len= len(np.unique(data.IAGA.values))
        print('Rotating Magnetometer Data')
        print('0% Complete')
        for i, site in enumerate(np.unique(data.IAGA.values)):
            glon=details['GEOLON'][details['IAGA']== site].values
            glat=details['GEOLAT'][details['IAGA']== site].values
            index= temp_df['Site']==site
            temp_df.loc[index, 'glon']= glon
            temp_df.loc[index, 'glat']= glat
            dec=chaos_dec(glat, glon, year)
            temp_df.loc[index, 'Declination']= dec* 180/np.pi
            temp_df.loc[index, 'Btheta'], temp_df.loc[index, 'Bphi']= SuperMag2Geo(data.N[index].values, data.E[index].values, dec)
            print(f'{round(i/data_len *100)}% Complete')
        for column in columns[1:-1]:
            exec('temp_df.'+column+'=temp_df.'+column+'.astype(float)')
        temp_df['Br']= data.Z.values
        HDFStore.append('Main', temp_df, format= 'table',data_columns=True)
    print(f"Rotation and Merging of data set complete: \n{out_folder+'Data.hdf5'}")
    HDFStore.close()
