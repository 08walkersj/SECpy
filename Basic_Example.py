#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 20 12:25:02 2022

@author: simon
"""

from SECpy import SECS, RE
import numpy as np
import matplotlib.pyplot as plt
lons= np.arange(-20.5, 21, 1)
lats= np.arange(60.5, 81, 1)
node_lon= np.array([0])
node_lat= np.array([70])
lons, lats= np.meshgrid(lons, lats)
poles= SECS(node_lon, node_lat, lons, lats)

fig=plt.figure()
ax=fig.add_subplot(121)
Cur=poles.eval_G_Matrix_J()
ax.quiver(lons, lats, *Cur)
ax.set_title('Divergence Free')

Cur=poles.eval_G_Matrix_J(system='curl-free')
ax2=fig.add_subplot(122)
ax2.quiver(lons, lats, *Cur)
ax2.set_title('Curl Free')