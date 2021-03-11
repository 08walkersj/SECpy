#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 13:48:27 2021

@author: simon
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from SECpy import Vdf, theta
import pickle
lon, lat =np.meshgrid(np.linspace(-20, 20), np.linspace(-20, 20))

fig= plt.figure()
ax= fig.add_subplot(111, projection='3d')

z= Vdf(theta(0, 0, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(10, 0, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(-10, 0, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(0, -10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(0, 10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(10, 0, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(-10, 10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(10, -10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(10, 10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(10, -10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(-10, -10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(-10, -10, lon, lat), 110e3, theta0=np.deg2rad(5))
z+= Vdf(theta(-10, 10, lon, lat), 110e3, theta0=np.deg2rad(5))
ax.plot_surface(lon, lat, z, rstride=1, cstride=1,
                cmap='winter', edgecolor='none')
ax.set_ylabel('latitude')
ax.set_xlabel('longitude')
ax.set_zlabel('Current')
fig.suptitle('Singularity Adjusted')
pickle.dump(fig, open('singularity.pickle', 'wb'))

fig= plt.figure()
ax= fig.add_subplot(111, projection='3d')

z= Vdf(theta(0, 0, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(10, 0, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(-10, 0, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(0, -10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(0, 10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(10, 0, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(-10, 10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(10, -10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(10, 10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(10, -10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(-10, -10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(-10, -10, lon, lat), 110e3, theta0=0)
z+= Vdf(theta(-10, 10, lon, lat), 110e3, theta0=0)
ax.plot_surface(lon, lat, z, rstride=1, cstride=1,
                cmap='winter', edgecolor='none')
ax.set_ylabel('latitude')
ax.set_xlabel('longitude')
ax.set_zlabel('Current')
fig.suptitle('Normal')
pickle.dump(fig, open('normal.pickle', 'wb'))
