#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  1 22:57:36 2021

@author: simon
"""

from netCDF4 import Dataset
year= 2014
month=8
day= 5
file= f'/home/simon/BCSS-DAG Dropbox/Simon Walker/Iridium/{year:04d}/{year:04d}{month:02d}{day:02d}Amp_invert.ncdf'
data= Dataset(file)
dimensions= ['datapnt', 'vec_comp', 'plane', 'sv']
headers= ['time', 'pseudo_sv_num', 'plane_num', 'pos_eci', 'b_eci', 'b_error', 
          'pseudo_sv_quality', 'data_splice', 'sv_matrix']