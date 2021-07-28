#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:27:58 2020
@author: Simon Walker

"""
#@author: zef014
import pandas as pd
import numpy as np
import os
import traceback
import time
import matplotlib.pyplot as plt
import gc
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
        poleward_boundary[int(row), :l]= pb[r[top3_lower]==row][reorder]
        equatorward_boundary[int(row), :l]= eb[r[top3_upper]==row][reorder]
        peak_values[int(row), :l]= values[:,1][r[peak_index]==row][reorder]
        peak_locations[int(row), :l]= mlats[peak_index][r[peak_index]==row][reorder]
    gc.collect()
    peak_values[peak_values==9.999e3]=np.nan
    return equatorward_boundary, poleward_boundary, peak_locations, peak_values+quantile, total_Current
def result(results, num_data):
    """
    Used to merge the results from the parallelisation

    Parameters
    ----------
    results : list
        The results from the parallelisation.
    num_data : int
        Number of different data that are returned from the function used in the parallelisation.

    Returns
    -------
    ret : list
        List containing the reordered data that can be called in a way that is the same as if the function used with the parallelisation had been used without.

    """
    ret=[]
    for i in range(num_data):
        ret.append(np.concatenate(np.array(results)[:, i]))
    return ret
if __name__=='__main__':
    #Latitudes used in the original set up
    mlat= np.linspace(49, 81, 50)
    #Setting up the column names for the data frame and HDF file where the electrojet properties will be put
    columns=['Date_UTC', 'MLT']
    for string in ['_Eastward','_Westward']:
        for i in range(1,4):
            for title in ['Equatorward_Boundary', 'Poleward_Boundary','Peak_Lat', 'Peak_Value', 'Current', 'Width']:
                columns.append(title+string+str(i))
    #Loading the store containing the SECS analysis output from evaulating along the meridian at each time step
    store= pd.HDFStore('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/results_singularity_mod.hdf5', mode='r')
    #Loading pool which is used for parallelising this algorithm, allowing python to use more than one core and work faster
    from multiprocessing import Pool
    # 8 cores are selected as defualt adust accordingly
    p= Pool(8)
    #Records the time the algorithm began, this is used to track progress in order to see it will be completed in a sensible amount of time
    start=time.time()
    #Chunks through the data store because it is too large to handle on most computers and less prone to crashes
    for j, data in enumerate(store.select(key='main', chunksize=40_000)):
        #Setting up data frame
        Output= pd.DataFrame(columns=columns)
        #Getting east currents in time order
        Currents= data[data.columns[2:52]].values
        #Getting the times in the data set
        Output.Date_UTC=data.Date_UTC
        #Getting the MLTs of the data set
        Output.MLT= data.MLT
        #Splitting the currents into the first 7 equal sized chunks to be used across 7 cores
        #negative of the currents is used so the westward electrojet properties cane be found
        splits= np.split(-(Currents)[:len(Currents)-(len(Currents)%7)], 7)
        #Any remaining data, if data set length is not a multiple of 7, is placed on the 8th core
        if len(Currents)%7 !=0:
            splits.append(-(Currents)[len(Currents)-(len(Currents)%7):])
        #Peforms the algorithm and places the data in the corresponding collumn
        Output[columns[20:][::6]], Output[columns[20:][1::6]], Output[columns[20:][2::6]], Output[columns[20:][3::6]], Output[columns[20:][4::6]]=result(p.map(main_current, splits), 5)
        #Calculates the width by finding the difference between the boundaries
        Output[columns[20:][5::6]]= Output[columns[20:][1::6]].values- Output[columns[20:][::6]].values
        #Repeats for the eastward electrojet (no negative of the currents)
        splits= np.split((Currents)[:len(Currents)-(len(Currents)%7)], 7)
        if len(Currents)%7 !=0:
            splits.append((Currents)[len(Currents)-(len(Currents)%7):])
        Output[columns[2:20][::6]], Output[columns[2:20][1::6]], Output[columns[2:20][2::6]], Output[columns[2:20][3::6]], Output[columns[2:20][4::6]]=result(p.map(main_current, splits), 5)
        
        Output[columns[2:20][5::6]]= Output[columns[2:20][1::6]].values- Output[columns[2:20][::6]].values
        #Change no data indicator to be a NAN
        Output=Output.replace(9.999e3, np.nan)
        #Change the westward electrojet peaks and integrals to be negative as they should be
        Output[columns[20:][3::6]]=Output[columns[20:][3::6]]*-1
        Output[columns[20:][4::6]]=Output[columns[20:][4::6]]*-1
        Output= pd.concat((Output, data[data.columns[152:]]), axis=1)
        #Convert each column into a float64
        for c in Output.columns[1:]:
            exec('Output.'+c+'=Output.'+c+".astype('float64')")
        #Output into a HDF file
        Output.to_hdf('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/electrojet_boundaries_singularity_mod.hdf5', mode='a', key='main', append=True,format='t', data_columns=True)
        #Record the interation that has been completed
        np.savetxt('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/latest_iterationj.txt', np.array([j]))
        #Record the latest date that has been done
        open('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/latest_iterationd.txt', 'w').write(str(data.Date_UTC.values[-1]))
        #Clear any unnecessary data stored (pandas can some times leave junk on the disk that will reduce the amount of data that can be held and gets worse with each iteration)
        gc.collect()
    end=time.time()
    #Print the average time it takes to find the boudaries for each minute of data
    print('average per minute of data',(end-start)/store.get_storer('main').nrows)
    #Print the total time it took to perform the algorithm on all the data
    print('total time',(end-start))