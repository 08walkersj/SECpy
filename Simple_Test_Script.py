#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 15:12:31 2021

@author: simon
"""

from SECpy import SECS
import numpy as np
node_lons, node_lats= np.meshgrid(np.arange(0, 50, 5), np.arange(60, 80,5))
eval_lons, eval_lats= np.meshgrid(np.arange(5, 45, 2.5), np.arange(65, 75, 2.5))

poles= SECS(node_lons.flatten(), node_lats.flatten(), eval_lons.flatten(), eval_lats.flatten())
Lons= np.array([30.5, 21, 58])
Lats= np.array([64.9, 62, 73])
Btheta= np.array([10, 12, 16])
Bphi= np.array([12, 8, 13])
Br= np.array([-3, -2, -4])
cond=0.05
poles.Fitting_Matrix(Lons, Lats, cond=cond)
Amplitudes= poles.Amplitude(Btheta, Bphi, Br)


import matplotlib.pyplot as plt
fig= plt.figure()
ax= fig.add_subplot(111)
ax.scatter(node_lons, node_lats, s=1, color='red', label='poles')
ax.scatter(Lons, Lats, color='orange', label='measurments')
s=ax.scatter(node_lons, node_lats, c= Amplitudes, cmap= 'seismic')
ax.legend(loc='best')
ax.set_xlabel('longitude')
ax.set_ylabel('lattitude')
cbar=fig.colorbar(s, ax=ax)
cbar.set_label('Pole Amplitudes')
fig.suptitle(f'Cond={cond} SVD Test')