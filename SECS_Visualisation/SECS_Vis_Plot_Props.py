#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 18:04:50 2021

@author: simon
"""
import pandas as pd
import numpy as np
def prop_plot(electro_property, season, mlt, path, axis1, axis2):
    store= pd.HDFStore(path, mode='r')
    key= electro_property
    if season!='':
        key+= f'_{season}'
    if mlt!='':
        key+=f'_{mlt}'
    df= store.select(key=key)
    store.close()
    if electro_property== 'Current':
        df.Mean_Eastward_pos=df.Mean_Eastward_pos* (np.pi/180)*((6.371E6+110e3)/1e3)
        df.Mean_Eastward_neg=df.Mean_Eastward_neg* (np.pi/180)*((6.371E6+110e3)/1e3)
        df.Error_Eastward_pos=df.Error_Eastward_pos* (np.pi/180)*((6.371E6+110e3)/1e3)
        df.Error_Eastward_neg=df.Error_Eastward_neg* (np.pi/180)*((6.371E6+110e3)/1e3)
        df.Mean_Westward_pos=df.Mean_Westward_pos* (np.pi/180)*((6.371E6+110e3)/1e3)
        df.Mean_Westward_neg=df.Mean_Westward_neg* (np.pi/180)*((6.371E6+110e3)/1e3)
        df.Error_Westward_pos=df.Error_Westward_pos* (np.pi/180)*((6.371E6+110e3)/1e3)
        df.Error_Westward_neg=df.Error_Westward_neg* (np.pi/180)*((6.371E6+110e3)/1e3)
    #Eastward pos
    East=axis1.plot(df.BY_pos_East, df.Mean_Eastward_pos, marker='o', markersize=1, label='BY +ve', color='green')
    East_errors=axis1.plot(df.BY_pos_East, df.Mean_Eastward_pos+df.Error_Eastward_pos, marker='o', markersize=1, color='red', label='error')
    East_errors.extend(axis1.plot(df.BY_pos_East, df.Mean_Eastward_pos-df.Error_Eastward_pos, marker='o', markersize=1, color='red'))
    #Eastward neg
    East.extend(axis1.plot(abs(df.BY_neg_East), df.Mean_Eastward_neg, marker='o', markersize=1, linestyle='--', label='BY -ve', color='green'))
    East_errors.extend(axis1.plot(abs(df.BY_neg_East), df.Mean_Eastward_neg+df.Error_Eastward_neg, marker='o', markersize=1, linestyle='--', color='red'))
    East_errors.extend(axis1.plot(abs(df.BY_neg_East), df.Mean_Eastward_neg-df.Error_Eastward_neg, marker='o', markersize=1, linestyle='--', color='red'))
    #Westward pos
    West=axis2.plot(df.BY_pos_West, -1*abs(df.Mean_Westward_pos), marker='o',markersize=1, color='green')
    West_errors=axis2.plot(df.BY_pos_West, -1*abs(df.Mean_Westward_pos+df.Error_Westward_pos), marker='o', markersize=1, color='red')
    West_errors.extend(axis2.plot(df.BY_pos_West, -1*abs(df.Mean_Westward_pos-df.Error_Westward_pos), marker='o', markersize=1, color='red'))
    #Westward neg
    West.extend(axis2.plot(abs(df.BY_neg_West), -1*abs(df.Mean_Westward_neg), marker='o', markersize=1, linestyle='--', color='green'))
    West_errors.extend(axis2.plot(abs(df.BY_neg_West), -1*abs(df.Mean_Westward_neg+df.Error_Westward_neg), marker='o', markersize=1, linestyle='--', color='red'))
    West_errors.extend(axis2.plot(abs(df.BY_neg_West), -1*abs(df.Mean_Westward_neg-df.Error_Westward_neg), marker='o', markersize=1, linestyle='--', color='red'))
    return East, East_errors, West, West_errors


if __name__== '__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    fig= plt.figure()
    gs= fig.add_gridspec(2, 1)
    ax= fig.add_subplot(gs[0])
    ax2=fig.add_subplot(gs[1])
    #Setting up broken subplots
    d = .015
    ax.set_xticks([])
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax.tick_params(labeltop=False)  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    prop='Peak_Value'
    season='winter'
    mlt='23'
    plots= np.concatenate(prop_plot(prop, season, mlt, './Plotting_Data/SECS_Visualisation.hdf5', ax, ax2))
    handels, labels= ax.get_legend_handles_labels()
    handels2, labels2= ax2.get_legend_handles_labels()
    ax.legend(handels+handels2, labels+labels2, loc='upper left')
    ax.set_ylabel('Eastward '+prop)
    ax2.set_ylabel('Westward '+prop)