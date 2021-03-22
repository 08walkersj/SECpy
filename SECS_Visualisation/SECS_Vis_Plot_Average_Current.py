#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 16:26:45 2021

@author: simon
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:48:19 2021

@author: simon
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
def average_current_plot(clock_angle, season, path, axis):
    global df
    mlat= np.linspace(49, 81, 50)
    store= pd.HDFStore(path, mode='r')
    key= 'Average_Current'
    if season!='':
        key+= f'_{season}'
    if clock_angle!='':
        key+= f'_Clock:{clock_angle}'
    df= store.select(key=key)
    store.close()
    mlt= np.arange(0.25, 24, .5)
    scat=axis.scatter(np.array([mlt]*len(mlat)), np.concatenate([np.vstack(mlat)]*len(mlt), axis=1), color='black', s=1/2)
    quiv=axis.quiver(mlt, mlat, df[[f'East_{i}' for i in mlt]].values, df[[f'North_{i}' for i in mlt]].values,  scale=5, headaxislength=0, headlength=0, zorder=2)
    mesh=axis.pcolormesh(np.append(mlt, mlt[0]), mlat, np.append(df[[f'Br_{i}' for i in mlt]].values, df[[f'Br_{i}' for i in mlt]].values[:,0].reshape(50,1), axis=1)*1e9, cmap=mpl.cm.seismic, vmin=-150, vmax=150, zorder=0)
    axis.set_rmax(90)
    return [quiv], [mesh], [scat]
    
    
    


if __name__== '__main__':
    import matplotlib.pyplot as plt
    from ideas import polar
    fig= plt.figure()
    gs= fig.add_gridspec(1, 2, width_ratios=[10, 1])
    ax= polar(fig.add_subplot(gs[0], projection='polar'))
    cax= fig.add_subplot(gs[1])
    season='summer'
    clock_angle='90'
    plots= average_current_plot(clock_angle, season, './Plotting_Data/SECS_Visualisation.hdf5', ax)
    cbar=fig.colorbar(mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-150, vmax=150), cmap=mpl.cm.seismic), cax= cax)
