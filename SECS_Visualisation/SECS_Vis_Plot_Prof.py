#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:48:19 2021

@author: simon
"""

import pandas as pd
import numpy as np
def profile_plot(clock_angle, season, mlt, path, axis1, axis2):
    mlat= np.linspace(49, 81, 50)
    colors= ['C0', 'C1', 'C2']
    store= pd.HDFStore(path, mode='r')
    key= 'Profile'
    if season!='':
        key+= f'_{season}'
    if mlt!='':
        key+=f'_{mlt}'
    if clock_angle!='':
        key+= f'_Clock:{clock_angle}'
    df= store.select(key=key)
    store.close()
    if clock_angle=='':
        columns= [str(title)+'_Current_East' for title in mlat]
        east_pos= df.Jet_BY== 'East_pos'
        east_neg= df.Jet_BY== 'East_neg'
        west_pos= df.Jet_BY== 'West_pos'
        west_neg= df.Jet_BY== 'West_neg'
        #Eastward pos
        East= np.concatenate([axis1.plot(mlat, Curr, 
                                         color=color, label=f"IMF BY: +{round(BY, 3)}") for BY, Curr, color in zip(df.BY_GSM[east_pos].values, 
                                                                                    df.loc[east_pos, columns].values,colors)])
        #Eastward neg
        East= np.append(East, np.concatenate([axis1.plot(mlat, Curr, 
                                                         color=color, label=f"IMF BY: {round(BY, 3)}", 
                                                         linestyle='--') for BY, Curr, color in zip(df.BY_GSM[east_neg].values, 
                                                                                    df.loc[east_neg, columns].values,colors)]))
        #Westward pos
        West= np.concatenate([axis2.plot(mlat, Curr, 
                                         color=color, label=f"IMF BY: +{round(BY, 3)}") for BY, Curr, color in zip(df.BY_GSM[west_pos].values, 
                                                                                    df.loc[west_pos, columns].values,colors)])
        #Westward neg
        West= np.append(West, np.concatenate([axis2.plot(mlat, Curr, 
                                                         color=color, label=f"IMF BY: {round(BY, 3)}", 
                                                         linestyle='--') for BY, Curr, color in zip(df.BY_GSM[west_neg].values, 
                                                                                    df.loc[west_neg, columns].values,colors)]))        
        return East, West
    else:
        columns= [str(title)+'_Current_East' for title in mlat]
        east= df.Jet=='East'
        west= df.Jet=='West'
        #Eastward
        East= axis1.plot(mlat, df.loc[east, columns].values.flatten(), color=colors[0], label='East')
        #Westward
        West= axis2.plot(mlat, df.loc[west, columns].values.flatten(), color=colors[0], label='West')

        
        return East, West
    
    
    


if __name__== '__main__':
    import matplotlib.pyplot as plt
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
    season=''
    mlt='12'
    clock_angle='180'
    plots= np.concatenate(profile_plot(clock_angle, season, mlt, './Plotting_Data/SECS_Visualisation.hdf5', ax, ax2))
    handels, labels= ax.get_legend_handles_labels()
    handels2, labels2= ax2.get_legend_handles_labels()
    # ax.legend(handels+handels2, labels+labels2, loc='upper left')
    ax.legend(loc='upper left')

    ax.set_ylabel('Eastward Profile')
    ax2.set_ylabel('Westward Profile')