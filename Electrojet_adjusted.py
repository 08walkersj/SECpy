#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 16:04:29 2020
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
def main_current(a):
    Currents, q= a[0]
    mlat= np.linspace(49, 81, 50)    
    quantile= np.vstack(np.quantile(abs(Currents), q, axis=1))
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
        # turning_points.append(mlats[turning_index])
        diff= np.zeros(mlats.shape)
        diff[:]=-9.9999e3
        diff[turning_index]=np.append(np.diff(np.cumsum(Currents-quantile, axis=1)[turning_index]),[-9.9999e3])
        top3_lower.append(np.in1d(np.round(abs(diff[row==r]), 15), np.round(np.sort(abs(diff[turning_index & (diff>0)]))[-3:], 15)))
        diff[turning_index]=np.append([-9.9999e3], np.diff(np.cumsum(Currents-quantile, axis=1)[turning_index]))
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
    dif[:, ::2]=np.append(np.vstack([0]*Currents.shape[0]), np.diff((Currents-quantile)[:, ::2])/(2* np.diff(mlat)[0]), axis=1)
    dif[:,1::2]=np.append(np.diff((Currents-quantile)[:, 1::2])/(2* np.diff(mlat)[0]),np.vstack([0]*Currents.shape[0]), axis=1)
    difs= np.array([arr[arr> np.nanmean(arr)*0.6] for arr in np.array(np.split(abs(dif.flatten()), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    # Value_index= np.array([[arr> np.nanmax(arr)*0.4, np.nanmax(arr)] for arr in np.array(np.split(abs(Currents).flatten(), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    values= np.array([[arr[arr> np.nanmax(arr)*0.4], np.nanmax(arr)] for arr in np.array(np.split(abs(Currents-quantile).flatten(), np.where((top3_lower|top3_upper).flatten())[0]))[1::2]])
    i=0
    integrals=[]
    pb=[]
    eb=[]
    for vals, peaks, d in zip(values[:,0], values[:,1], difs):
        tmp_index= np.in1d(abs(Currents-quantile), vals) |np.in1d(abs(dif), d)
        # area= (Currents[tmp_index.reshape(Currents.shape)][1:-1].sum() + (Currents[tmp_index.reshape(Currents.shape)][[0,-1]].sum())/2)*np.diff(mlat)[0]*-1
        integrals.append((Currents[tmp_index.reshape(Currents.shape)][1:-1].sum() + (Currents[tmp_index.reshape(Currents.shape)][[0,-1]].sum())/2)*np.diff(mlat)[0])
        pb.append(np.max(mlats[tmp_index.reshape(Currents.shape)]))
        eb.append(np.min(mlats[tmp_index.reshape(Currents.shape)]))
        if i==0:
            # value_index, peak_index = tmp_index, np.in1d(abs(Currents), peaks)
            peak_index =np.in1d(abs(Currents-quantile), peaks)
            i+=1
        else:
            # value_index |= tmp_index
            peak_index |=np.in1d(abs(Currents-quantile), peaks)
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
    Jet_num=np.zeros(len(total_Current))
    for row in np.unique(r[top3_lower]):
        l= len(integrals[r[top3_lower]==row])
        Jet_num[row]=l
        reorder= np.argsort(-integrals[r[top3_lower]==row])
        total_Current[int(row), :l]= integrals[r[top3_lower]==row][reorder]
        poleward_boundary[int(row), :l]= pb[r[top3_lower]==row][reorder]
        equatorward_boundary[int(row), :l]= eb[r[top3_upper]==row][reorder]
        peak_values[int(row), :l]= values[:,1][r[peak_index]==row][reorder]
        peak_locations[int(row), :l]= mlats[peak_index][r[peak_index]==row][reorder]
    # print('4')
    gc.collect()
    peak_values[peak_values==9.999e3]=np.nan
    return equatorward_boundary, poleward_boundary, peak_locations, peak_values+quantile, total_Current, np.array(Jet_num)
def result(results, num_data):
    ret=[]
    Jet_num= np.concatenate([res[-1] for res in results])
    results=np.array([res[:-1] for res in results])
    for i in range(num_data):
        ret.append(np.concatenate(np.array(results)[:, i]))
    ret.append(Jet_num)
    return ret
start=time.time()
mlat= np.linspace(49, 81, 50)
columns=['Date_UTC', 'MLT']
for string in ['_Eastward','_Westward']:
    for i in range(1,4):
        for title in ['Equatorward_Boundary', 'Poleward_Boundary','Peak_Lat', 'Peak_Value', 'Current', 'Width']:
            columns.append(title+string+str(i))

profile_store= pd.HDFStore('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/results_with_clock.hdf5', mode='r')
EJ_store= pd.HDFStore('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/electrojet_boundaries.hdf5', mode='r')
from multiprocessing import Pool
p= Pool(8)
for data in EJ_store.select(key='main', chunksize=40_000):
    data2=data[data[data.columns[20:-12]].isnull().any(axis=1)]
    W=data2[data2.Peak_Value_Westward1>-1].Date_UTC.values
    open('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/asjustents_w.txt', 'w').write(str(W))
    start= W[0]
    end= W[-1]
    prof= profile_store.select(key='main', where= "Date_UTC>= start & Date_UTC<=end")
    prof= prof[prof.Date_UTC.isin(W)]
    q=0.1
    Currents= prof[prof.columns[2:52]].values
    splits= [[(split, q)] for split in np.split(-(Currents)[:len(Currents)-(len(Currents)%7)], 7)]
    if len(Currents)%7!=0:
        splits.append([(-(Currents)[len(Currents)-(len(Currents)%7):], q)])
    
    eb, pb, pl, pv, tc, Jet_num= result(p.map(main_current, splits), 5)
    data.loc[data.Date_UTC.isin(W[Jet_num>=2]), data.columns[20:][:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][1:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][2:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][3:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][4:-6:6]]= eb[Jet_num>=2], pb[Jet_num>=2], pl[Jet_num>=2], pv[Jet_num>=2], tc[Jet_num>=2]
    while (Jet_num<2).any():
        q+=0.05
        W= W[Jet_num<2]
        prof= prof[prof.Date_UTC.isin(W)]
        Currents= prof[prof.columns[2:52]].values
        splits= [[(split, q)] for split in np.split(-(Currents)[:len(Currents)-(len(Currents)%7)], 7)]
        if len(Currents)%7 !=0:
            splits.append([(-(Currents)[len(Currents)-(len(Currents)%7):], q)])
        
        eb, pb, pl, pv, tc, Jet_num= result(p.map(main_current, splits), 5)
        if len(Jet_num)!= tc.shape[0]:
            print('fuuhuuuuuuck')
        data.loc[data.Date_UTC.isin(W[Jet_num>=2]), data.columns[20:][:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][1:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][2:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][3:-6:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[20:][4:-6:6]]= eb[Jet_num>=2], pb[Jet_num>=2], pl[Jet_num>=2], pv[Jet_num>=2]*-1, tc[Jet_num>=2]*-1
        if q>= 0.4:
            break
    data2=data[data[data.columns[2:14]].isnull().any(axis=1)]
    W=data2[data2.Peak_Value_Eastward1<1].Date_UTC.values
    start= W[0]
    end= W[-1]
    open('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/asjustents_w.txt', 'w').write(str(W))
    prof= profile_store.select(key='main', where= "Date_UTC>= start & Date_UTC<=end")
    prof= prof[prof.Date_UTC.isin(W)]
    q=0.1
    Currents= prof[prof.columns[2:52]].values
    splits= [[(split, q)] for split in np.split((Currents)[:len(Currents)-(len(Currents)%7)], 7)]
    if len(Currents)%7!=0:
        splits.append([((Currents)[len(Currents)-(len(Currents)%7):], q)])
    
    eb, pb, pl, pv, tc, Jet_num= result(p.map(main_current, splits), 5)
    data.loc[data.Date_UTC.isin(W[Jet_num>=2]), data.columns[2:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[3:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[4:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[5:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[6:20:6]]= eb[Jet_num>=2], pb[Jet_num>=2], pl[Jet_num>=2], pv[Jet_num>=2], tc[Jet_num>=2]
    while (Jet_num<2).any():
        q+=0.05
        W= W[Jet_num<2]
        prof= prof[prof.Date_UTC.isin(W)]
        Currents= prof[prof.columns[2:52]].values
        splits= [[(split, q)] for split in np.split((Currents)[:len(Currents)-(len(Currents)%7)], 7)]
        if len(Currents)%7 !=0:
            splits.append([((Currents)[len(Currents)-(len(Currents)%7):], q)])
        
        eb, pb, pl, pv, tc, Jet_num= result(p.map(main_current, splits), 5)
        if len(Jet_num)!= tc.shape[0]:
            print('fuuhuuuuuuck')
        data.loc[data.Date_UTC.isin(W[Jet_num>=2]), data.columns[2:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[3:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[4:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[5:20:6]], data.loc[data.Date_UTC.isin(W[Jet_num>=2]),data.columns[6:20:6]]= eb[Jet_num>=2], pb[Jet_num>=2], pl[Jet_num>=2], pv[Jet_num>=2], tc[Jet_num>=2]
        if q>= 0.4:
            break
    data=data.replace(9.999e3, np.nan)
    data[columns[20:][5::6]]= data[columns[20:][1::6]].values- data[columns[20:][::6]].values
    data[columns[2:20][5::6]]= data[columns[2:20][1::6]].values- data[columns[2:20][::6]].values
    data.to_hdf('/Home/siv32/zef014/Documents/Masters_Thesis_Code/Masters_Thesis/Data/Meridian_Analysis/electrojet_boundaries_adjusted.hdf5', mode='a', key='main', append=True,format='t', data_columns=True)