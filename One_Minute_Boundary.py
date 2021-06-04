#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 14:47:17 2021

@author: simon
"""

import pandas as pd
import numpy as np
import os
import traceback
import time
import matplotlib.pyplot as plt
import gc
import matplotlib as mpl
# import pickle
mpl.rcParams.update({'xtick.labelsize': 10, 'ytick.labelsize': 10, 
                     'xtick.major.pad':0.1, 'xtick.minor.pad':0,
                     'ytick.major.pad':0.1, 'ytick.minor.pad':0,
                     'xtick.major.size': 2.5, 'xtick.minor.size': 1.0,
                     'ytick.major.size': 2.5, 'ytick.minor.size': 1.0,
                     'font.size': 15, 'lines.linewidth': 1,
                     'lines.markersize': 4, 'axes.linewidth': 0.2,
                     'legend.fontsize': 7})
def integral(dx, y):
    """
    Finds the area under the curve using the trapezium method

    Parameters
    ----------
    dx : int/float
        Change in x for each data point, assumes consistent spacing of data points.
    y : numpy.ndarray
        y data in order in terms of its x counterpart.

    Returns
    -------
    float
        Area under the curve.

    """
    return dx*(y[0]+np.sum(y[1:-1])+y[-1])
# def main_current(mlat, Currents):
def main_current(Currents):
    """
    Find the boundaries, sheet current density integral, peak sheet current density and the peak location of the sheet current density. 
    Uses the 10% quantile of the data and assumes search for eastward electrojet (reverse the signs of the eastward current before input and the algorithm will find the electrojet properties of the westward electrojet)
    Finds the top 3 strongest electrojet in each profile and will return NAN when the properties cannot be found. 
    i.e. if there is only one electrojet then proprties will be found for that one and in the place of the other 2 there will be NANs.
    
    Parameters
    ----------
    Currents : numpy.ndarray
        Eastward current along the meridian, can take an array of arrays where each array is a new time step.

    Returns
    -------
    equatorward boundary : numpy.ndarray (magnetic lattitude)
        top 3 electrojets equatorward boundary, array of arrays where each array is a new time step.
    poleward boundary : numpy.ndarray (magnetic lattitude)
        top 3 electrojets poleward boundary, array of arrays where each array is a new time step.
    peak locations : numpy.ndarray (magnetic lattitude)
        top 3 electrojets peak sheet current location, array of arrays where each array is a new time step.
    peak sheet current: numpy.ndarray (A/m)
        top 3 electrojets peak sheet current, array of arrays where each array is a new time step.
    total current : numpy.ndarray (Adeg/m )
        top 3 electrojets integrated profile, array of arrays where each array is a new time step.

    """
    mlat= np.linspace(49, 81, 50)    
    quantile= np.vstack(np.quantile(abs(Currents), 0.1, axis=1))
    # Currents= Currents
    # print('1')
    # global top3_lower, index, tmp_index, top3_upper, row, dup, mlats, diff, turning_index, difs, r, dif, peak_index, values, value_index, integrals, pb, eb
    nrows= Currents.shape[0]
    mlats= np.array([mlat]*nrows)
    index= np.concatenate((abs(np.diff(np.sign(Currents-quantile)))==2, np.vstack([True]*nrows)), axis=1)
    index[:,0]= True
    r= np.concatenate([np.vstack(range(Currents.shape[0]))]*50, axis=1)
    # turning_points= []
    top3_lower=[]
    top3_upper=[]
    for row in np.unique(r[index]):
        turning_index=index & (row==r)
        diff= np.zeros(mlats.shape)
        diff[:]=-9.9999e3
        diff[turning_index]=np.append(np.diff(np.cumsum(Currents-quantile, axis=1)[turning_index]),[-9.9999e3])
        top3_lower.append(np.in1d(np.round(abs(diff[row==r]), 15), np.round(np.sort(abs(diff[turning_index & (diff>0)]))[-3:], 15)))
        diff[turning_index]=np.append([-9.9999e3], np.diff(np.cumsum(Currents-quantile, axis=1)[turning_index]))
        top3_upper.append(np.in1d(np.round(abs(diff[row==r]), 15), np.round(np.sort(abs(diff[turning_index & (diff>0)]))[-3:], 15)))
    top3_lower, top3_upper= np.array(top3_lower), np.array(top3_upper)
    dup=top3_lower&top3_upper
    top3_upper[dup]= False
    top3_lower[dup]= False
    dif=np.float128(r.copy())
    dif[:, ::2]=np.append(np.vstack([0]*Currents.shape[0]), np.diff((Currents-quantile)[:, ::2])/(2* np.diff(mlat)[0]), axis=1)
    dif[:,1::2]=np.append(np.diff((Currents-quantile)[:, 1::2])/(2* np.diff(mlat)[0]),np.vstack([0]*Currents.shape[0]), axis=1)
    difs= np.array([arr[arr> np.nanmean(arr)*0.6] for arr in np.array(np.split(abs(dif.flatten()), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    values= np.array([[arr[arr> np.nanmax(arr)*0.4], np.nanmax(arr)] for arr in np.array(np.split(abs(Currents-quantile).flatten(), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    i=0
    integrals=[]
    pb=[]
    eb=[]
    for vals, peaks, d in zip(values[:,0], values[:,1], difs):
        tmp_index= np.in1d(abs(Currents-quantile), vals) |np.in1d(abs(dif), d)
        integrals.append((Currents[tmp_index.reshape(Currents.shape)][1:-1].sum() + (Currents[tmp_index.reshape(Currents.shape)][[0,-1]].sum())/2)*np.diff(mlat)[0])
        pb.append(np.max(mlats[tmp_index.reshape(Currents.shape)]))
        eb.append(np.min(mlats[tmp_index.reshape(Currents.shape)]))
        if i==0:
            peak_index =np.in1d(abs(Currents-quantile), peaks)
            i+=1
        else:
            peak_index |=np.in1d(abs(Currents-quantile), peaks)
    pb, eb = np.array(pb), np.array(eb)
    peak_index=peak_index.reshape(Currents.shape)
    integrals= np.array(integrals)
    total_Current= np.float128(r.copy()[:, :3])
    total_Current[:]=9.999e3
    poleward_boundary= total_Current.copy()
    equatorward_boundary= total_Current.copy()
    peak_values= total_Current.copy()
    peak_locations= total_Current.copy()
    for row in np.unique(r[top3_lower]):
        l= len(integrals[r[top3_lower]==row])
        reorder= np.argsort(-integrals[r[top3_lower]==row])
        total_Current[int(row), :l]= integrals[r[top3_lower]==row][reorder]
        poleward_boundary[int(row), :l]= pb[r[top3_lower]==row][reorder] +  np.diff(np.linspace(49, 81, 50))[0]
        equatorward_boundary[int(row), :l]= eb[r[top3_upper]==row][reorder]
        peak_values[int(row), :l]= values[:,1][r[peak_index]==row][reorder]
        peak_locations[int(row), :l]= mlats[peak_index][r[peak_index]==row][reorder]
    gc.collect()
    peak_values[peak_values==9.999e3]=np.nan
    return equatorward_boundary, poleward_boundary, peak_locations, (peak_values+quantile)*1e3, total_Current
path='/home/simon/gits/SECpy/SECS_Visualisation/SECS_Data/results_singularity_mod_clock2.hdf5'
store= pd.HDFStore(path, mode='r')
Date= pd.Timestamp('2000-02-05T22:34')
data= store.select(key='main', where='Date_UTC=Date')
equatorward_boundary, poleward_boundary, peak_locations, peak_values, total_Current= main_current(data[data.columns[2:52]].values)
East= {}
for i in ['equatorward_boundary', 'poleward_boundary', 'peak_locations', 'peak_values', 'total_Current']:East.update(eval("{i:"+i+'}'))

equatorward_boundary, poleward_boundary, peak_locations, peak_values, total_Current= main_current(-1*data[data.columns[2:52]].values)
West= {}
for i in ['equatorward_boundary', 'poleward_boundary', 'peak_locations', 'peak_values', 'total_Current']:West.update(eval("{i:"+i+'}'))

for i in East:
    East[i][East[i]== 9.999e3]=np.nan
for i in West:
    West[i][West[i]== 9.999e3]=np.nan
fig= plt.figure()
ax= fig.add_subplot(111)
ax.plot(np.linspace(49, 81, 50), data[data.columns[2:52]].values.flatten()*1e3)
xlim= ax.get_xlim()
ylim= ax.get_ylim()
ax.scatter(East['peak_locations'], East['peak_values'], color='green', label='Eastward EJ Peaks')
ax.scatter(West['peak_locations'], -West['peak_values'], color='red', label='Westward EJ Peaks')
ax.plot(xlim, [0]*2, color='black')
for val in East['equatorward_boundary']:
    ax.plot([val]*2, ylim, label='Eastward EJ Boundary', linestyle='--', color='green')
for val in East['poleward_boundary']:
    ax.plot([val]*2, ylim, label='Eastward EJ Boundary', linestyle='--', color='green')
for val in West['equatorward_boundary']:
    ax.plot([val]*2, ylim, label='Westward EJ Boundary', linestyle='--', color='red')
for val in West['poleward_boundary']:
    ax.plot([val]*2, ylim, label='Westward EJ Boundary', linestyle='--', color='red')
ax.plot(xlim, [np.quantile(abs(data[data.columns[2:52]].values), 0.1, axis=1)*1e3]*2, linestyle='--', color='lightblue', label='quantile')
ax.plot(xlim, [-np.quantile(abs(data[data.columns[2:52]].values), 0.1, axis=1)*1e3]*2, linestyle='--', color='lightblue', label='quantile')
handles, labels =ax.get_legend_handles_labels()
labels, index= np.unique(labels, return_index=True)
handles=np.array(handles)[[index]]
leg=ax.legend(handles, labels, loc='upper right')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_frame_on(False)
ax.set_xlabel(r'$^o$mlat')
ax.set_ylabel(r'Eastward Sheet Current Density $Akm^{-1}$')
fig.suptitle('Electojet Detection Algorithm')
#Resize the figure
DPI=fig.get_dpi()
fig.set_size_inches(2400.0/float(DPI),1220.0/float(DPI))