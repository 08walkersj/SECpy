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
    return dx*(y[0]+np.sum(y[1:-1])+y[-1])
# def main_current(mlat, Currents):
def main_current(Currents):
    mlat= np.linspace(49, 81, 50)
    # Currents= Currents
    # print('1')
    # global top3_lower, index, tmp_index, top3_upper, row, dup, mlats, diff, turning_index, difs, r, dif, peak_index, values, value_index, integrals, pb, eb
    nrows= Currents.shape[0]
    mlats= np.array([mlat]*nrows)
    index= np.concatenate((abs(np.diff(np.sign(Currents)))==2, np.vstack([True]*nrows)), axis=1)
    index[:,0]= True
    r= np.concatenate([np.vstack(range(Currents.shape[0]))]*50, axis=1)
    # turning_points= []
    top3_lower=[]
    top3_upper=[]
    for row in np.unique(r[index]):
        turning_index=index & (row==r)
        # turning_points.append(mlats[turning_index])
        diff= np.zeros(mlats.shape)
        diff[:]=-9.9999e3
        diff[turning_index]=np.append(np.diff(np.cumsum(Currents, axis=1)[turning_index]),[-9.9999e3])
        top3_lower.append(np.in1d(np.round(abs(diff[row==r]), 15), np.round(np.sort(abs(diff[turning_index & (diff>0)]))[-3:], 15)))
        diff[turning_index]=np.append([-9.9999e3], np.diff(np.cumsum(Currents, axis=1)[turning_index]))
        top3_upper.append(np.in1d(np.round(abs(diff[row==r]), 15), np.round(np.sort(abs(diff[turning_index & (diff>0)]))[-3:], 15)))
    top3_lower, top3_upper= np.array(top3_lower), np.array(top3_upper)
    # dup=top3_lower[1::2]&top3_upper[1::2]
    # top3_upper[1::2][dup]= False
    # top3_lower[1::2][dup]= False
    dup=top3_lower&top3_upper
    top3_upper[dup]= False
    top3_lower[dup]= False
    # print('2')
    dif=np.float128(r.copy())
    dif[:, ::2]=np.append(np.vstack([0]*Currents.shape[0]), np.diff(Currents[:, ::2])/(2* np.diff(mlat)[0]), axis=1)
    dif[:,1::2]=np.append(np.diff(Currents[:, 1::2])/(2* np.diff(mlat)[0]),np.vstack([0]*Currents.shape[0]), axis=1)
    difs= np.array([arr[arr> np.nanmean(arr)*0.6] for arr in np.array(np.split(abs(dif.flatten()), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    # Value_index= np.array([[arr> np.nanmax(arr)*0.4, np.nanmax(arr)] for arr in np.array(np.split(abs(Currents).flatten(), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    values= np.array([[arr[arr> np.nanmax(arr)*0.4], np.nanmax(arr)] for arr in np.array(np.split(abs(Currents).flatten(), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    i=0
    integrals=[]
    pb=[]
    eb=[]
    for vals, peaks, d in zip(values[:,0], values[:,1], difs):
        tmp_index= np.in1d(abs(Currents), vals) |np.in1d(abs(dif), d)
        # area= (Currents[tmp_index.reshape(Currents.shape)][1:-1].sum() + (Currents[tmp_index.reshape(Currents.shape)][[0,-1]].sum())/2)*np.diff(mlat)[0]*-1
        integrals.append((Currents[tmp_index.reshape(Currents.shape)][1:-1].sum() + (Currents[tmp_index.reshape(Currents.shape)][[0,-1]].sum())/2)*np.diff(mlat)[0])
        pb.append(np.max(mlats[tmp_index.reshape(Currents.shape)]))
        eb.append(np.min(mlats[tmp_index.reshape(Currents.shape)]))
        if i==0:
            # value_index, peak_index = tmp_index, np.in1d(abs(Currents), peaks)
            peak_index =np.in1d(abs(Currents), peaks)
            i+=1
        else:
            # value_index |= tmp_index
            peak_index |=np.in1d(abs(Currents), peaks)
    # print('3')
    pb, eb = np.array(pb), np.array(eb)
    # value_index, peak_index= value_index.reshape(Currents.shape), peak_index.reshape(Currents.shape)
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
        total_Current[int(row), :l]= integrals[r[top3_lower]==row]
        poleward_boundary[int(row), :l]= pb[r[top3_lower]==row]
        equatorward_boundary[int(row), :l]= eb[r[top3_upper]==row]
        peak_values[int(row), :l]= values[:,1][r[peak_index]==row]
        peak_locations[int(row), :l]= mlats[peak_index][r[peak_index]==row]
    # print('4')
    gc.collect()
    return equatorward_boundary, poleward_boundary, peak_locations, peak_values, total_Current
def result(results, num_data):
    ret=[]
    for i in range(num_data):
        ret.append(np.concatenate(np.array(results)[:, i]))
    return ret
start=time.time()
mlat= np.linspace(49, 81, 50)
columns=['Date_UTC', 'MLT']
for string in ['_Eastward','_Westward']:
    for i in range(1,4):
        for title in ['Equatorward_Boundary', 'Poleward_Boundary','Peak_Lat', 'Peak_Value', 'Current', 'Width']:
            columns.append(title+string+str(i))
# date=pd.Timestamp('2010-03-01 09:28:00')
# end= pd.Timestamp('2010-03-02 09:34:00')
store= pd.HDFStore('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/results_with_clock.hdf5', mode='r')
from multiprocessing import Pool
p= Pool(8)
start=time.time()
for j, data in enumerate(store.select(key='main', chunksize=40_000)):
    print(j)
# for j, data in enumerate(store.select(key='main', chunksize=100)):
    Output= pd.DataFrame(columns=columns)
    # data= pd.read_hdf('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/results_with_clock.hdf5', where= "Date_UTC>=date & Date_UTC<end")
    Currents= data[data.columns[2:52]].values
    Output.Date_UTC=data.Date_UTC
    Output.MLT= data.MLT
    # Output[columns[20:][::6]], Output[columns[20:][1::6]], Output[columns[20:][2::6]], Output[columns[20:][3::6]], Output[columns[20:][4::6]]=equatorward_boundary, poleward_boundary, peak_locations, peak_values, total_Current=main_current(mlat, -(Currents+np.vstack(np.quantile(abs(Currents), 0.1, axis=1))))
    
    splits= np.split(-(Currents+np.vstack(np.quantile(abs(Currents), 0.1, axis=1)))[:len(Currents)-(len(Currents)%7)], 7)
    if len(Currents)%7 !=0:
        splits.append(-(Currents+np.vstack(np.quantile(abs(Currents), 0.1, axis=1)))[len(Currents)-(len(Currents)%7):])
    Output[columns[20:][::6]], Output[columns[20:][1::6]], Output[columns[20:][2::6]], Output[columns[20:][3::6]], Output[columns[20:][4::6]]=result(p.map(main_current, splits), 5)
    
    # Output[columns[20:][::6]], Output[columns[20:][1::6]], Output[columns[20:][2::6]], Output[columns[20:][3::6]], Output[columns[20:][4::6]]=main_current(mlat, -(Currents+np.vstack(np.quantile(abs(Currents), 0.1, axis=1))))
    Output[columns[20:][5::6]]= Output[columns[20:][1::6]].values- Output[columns[20:][::6]].values
    
    splits= np.split((Currents-np.vstack(np.quantile(abs(Currents), 0.1, axis=1)))[:len(Currents)-(len(Currents)%7)], 7)
    if len(Currents)%7 !=0:
        splits.append((Currents-np.vstack(np.quantile(abs(Currents), 0.1, axis=1)))[len(Currents)-(len(Currents)%7):])
    Output[columns[2:20][::6]], Output[columns[2:20][1::6]], Output[columns[2:20][2::6]], Output[columns[2:20][3::6]], Output[columns[2:20][4::6]]=result(p.map(main_current, splits), 5)
    
    # Output[columns[2:20][::6]], Output[columns[2:20][1::6]], Output[columns[2:20][2::6]], Output[columns[2:20][3::6]], Output[columns[2:20][4::6]]=equatorward_boundary, poleward_boundary, peak_locations, peak_values, total_Current=main_current(mlat, Currents+np.vstack(np.quantile(abs(Currents), 0.1, axis=1)))
    # Output[columns[2:20][::6]], Output[columns[2:20][1::6]], Output[columns[2:20][2::6]], Output[columns[2:20][3::6]], Output[columns[2:20][4::6]]=main_current(mlat, Currents+np.vstack(np.quantile(abs(Currents), 0.1, axis=1)))
    Output[columns[2:20][5::6]]= Output[columns[2:20][1::6]].values- Output[columns[2:20][::6]].values
    Output=Output.replace(9.999e3, np.nan)
    Output[columns[20:][3::6]]=Output[columns[20:][3::6]]*-1
    Output[columns[20:][4::6]]=Output[columns[20:][4::6]]*-1
    Output= pd.concat((Output, data[data.columns[152:]]), axis=1)
    for c in Output.columns[1:]:
        exec('Output.'+c+'=Output.'+c+".astype('float128')")
    Output.to_hdf('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/electrojet_boundaries.hdf5', mode='a', key='main', append=True,format='t', data_columns=True)
    np.savetxt('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/latest_iterationj.txt', np.array([j]))
    open('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/latest_iterationd.txt', 'w').write(str(data.Date_UTC.values[-1]))
    gc.collect()
    # print((start-time.time())/data.shape[0])
end=time.time()
print('average per minute of data',(end-start)/store.get_storer('main').nrows)
print('total time',(end-start))

# date=pd.Timestamp('2010-03-01 09:28:00')
# end= pd.Timestamp('2010-03-01 09:31:00')
# data= pd.read_hdf('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/results.hdf5', where= "Date_UTC>=date & Date_UTC<end")
# Currents= data[data.columns[2:52]].values
# equatorward_boundary, poleward_boundary, peak_locations, peak_values, total_Current=main_current((Currents-np.vstack(np.quantile(abs(Currents), 0.01, axis=1))))

# fig, (axes)= plt.subplots(1, 2)
# for ax, cur, eq, pb, m in zip(axes, Currents, equatorward_boundary, poleward_boundary, np.quantile(abs(Currents), 0.01, axis=1)):
#     ax.plot(mlat, cur)
#     ax.plot([eq[eq!=9.999e3]]*2, ax.get_ylim(), color='red', linestyle='--')
#     ax.plot([pb[eq!=9.999e3]]*2, ax.get_ylim(), color='green', linestyle='--')
#     ax.plot(ax.get_xlim(), [m]*2)
