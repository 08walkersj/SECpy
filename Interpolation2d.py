#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 17:19:34 2021

@author: simon
"""
import numpy as np
import matplotlib.pyplot as plt
def interpolation2d(x, y, f, xeval, yeval):
    """
    Interpolates linearly a 2d regular grid (x, y) with values f to a 
    new grid or set of values (xeval, yeval) by using the 4 closest coordinates 
    to the evaluation coordinate.

    Parameters
    ----------
    x : numpy.ndarray (1D)
        x- co ordinate. The unique values that describe the x co ordinates that make up the initial grid.
    y : numpy.ndarray (1D)
        y- co ordinate. The unique values that describe the y co ordinates that make up the initial grid.
    f : numpy.ndarray (2D)
        Values to interpolate with dimensions corresponding to the x and y co ordinates and shape (x,y)
    xeval : numpy.ndarray (1D)
        x- co ordinates of the evaluation points. If using a grid then input a flattened xeval.
    yeval : numpy.ndarray (1D)
        y- co ordinates of the evaluation points. If using a grid then input a flattened yeval.

    Returns
    -------
    results : numpy.ndarray(1D)
        Interpolated values that correspond to the cordinates of xeval and yeval.

    """
    x2=np.array([x]*len(xeval))
    y2=np.array([y]*len(yeval))
    diffx= x2-np.vstack(xeval)
    diffxL=np.zeros(diffx.shape)
    diffxL[:]=np.nan
    diffxL[x2<=np.vstack(xeval)]= diffx[x2<=np.vstack(xeval)]
    diffxG=np.zeros(diffx.shape)
    diffxG[:]=np.nan
    diffxG[x2>=np.vstack(xeval)]= diffx[x2>=np.vstack(xeval)]
    diffy= y2-np.vstack(yeval)
    diffyL=np.zeros(diffy.shape)
    diffyL[:]=np.nan
    diffyL[y2<=np.vstack(yeval)]= diffy[y2<=np.vstack(yeval)]
    diffyG=np.zeros(diffy.shape)
    diffyG[:]=np.nan
    diffyG[y2>=np.vstack(yeval)]= diffy[y2>=np.vstack(yeval)]
    x0_ind= np.nanargmin(np.abs(diffxL),axis=1)
    x1_ind= np.nanargmin(np.abs(diffxG),axis=1)
    y0_ind= np.nanargmin(np.abs(diffyL),axis=1)
    y1_ind= np.nanargmin(np.abs(diffyG),axis=1)
    results= np.zeros(xeval.shape)
    results= f[x0_ind, y0_ind]
    dx= xeval-x[x0_ind]
    dy= yeval-y[y0_ind]
    Δx= x[x1_ind]-x[x0_ind]
    Δy= y[y1_ind]-y[y0_ind]
    fracx= (dx/Δx)
    fracy= (dy/Δy)
    ind= ((dx==0)&(dy==0))|((np.isnan(fracx))&(np.isnan(fracy)))
    fracx[np.isnan(fracx)]=0
    fracy[np.isnan(fracy)]=0
    results[~ind]= fracy[~ind]*fracx[~ind]*f[x1_ind, y1_ind][[~ind]] + fracy[~ind]*(1-fracx[~ind])*f[x0_ind,y1_ind][~ind] + fracx[~ind]*(1-fracy[~ind])*f[x1_ind,y0_ind][~ind] + (1-fracx[~ind]-fracy[~ind]+fracx[~ind]*fracy[~ind])*f[x0_ind, y0_ind][~ind]
    return results

#%% Simple Example
if __name__=='__main__':
    def distribution(x, y):
        return np.sin((x+y)*(2*np.pi)/np.max(x+y, axis=1))
    x, y= np.meshgrid(np.arange(0, 50), np.arange(0,50))
    # x= x.flatten()
    # y= y.flatten()
    xeval=np.array([25]*50)
    yeval=np.linspace(0, 48, 50)
    f= distribution(x,y)
    fig= plt.figure()
    ax= fig.add_subplot(121)
    ax2= fig.add_subplot(122)
    s= ax.scatter(x, y, c=f, cmap='seismic')
    s2= ax2.scatter(xeval, yeval, c=interpolation2d(np.arange(0, 50), np.arange(0, 50), f.T, xeval, yeval), vmin=np.min(f), vmax= np.max(f), cmap='seismic')
    ax2.set_ylim(ax.get_ylim())
    ax2.set_xlim(ax.get_xlim())
    cbar= fig.colorbar(s, ax=ax)
    cbar2= fig.colorbar(s2, ax=ax2)
#%% SECS Example
if __name__=='__main__':
    from SECpy import SECS, Example
    node_lons, node_lats, lon, lat= Example.ExampleGrid()
    Btheta, Bphi, Br, MagLon, MagLat= Example.ExampleData()
    ind= (lon>=min(node_lons))&(lon<=max(node_lons))&(lat>=min(node_lats))&(lat<= max(node_lats))
    lon= lon[ind]
    lat= lat[ind]
    poles= SECS(node_lons, node_lats, lon, lat)
    poles.Fitting_Matrix(MagLon, MagLat)
    M= poles.Amplitude(Btheta, Bphi, Br)
    fig=plt.figure()
    ax= fig.add_subplot(121)
    ax2= fig.add_subplot(122)
    s= ax.scatter(node_lons, node_lats, c= M, cmap='seismic')
    s2= ax2.scatter(lon, lat, c=interpolation2d(np.unique(node_lons), np.unique(node_lats), np.vstack([M[node_lons==val] for val in np.unique(node_lons)]), lon, lat), cmap='seismic', vmin=np.min(M), vmax=np.max(M))
    ax2.set_ylim(ax.get_ylim())
    ax2.set_xlim(ax.get_xlim())
    cbar= fig.colorbar(s, ax=ax)
    cbar2= fig.colorbar(s2, ax=ax2)
#%% Cubed Sphere Example
if __name__=='__main__':
    from SECpy import SECS, Example
    import cubedsphere as CS
    from apexpy import Apex
    import datetime as dt
    
    """Grid"""
    lon_centre= 17.7
    lat_centre= 68.1
    A = Apex(date=dt.datetime(2008, 6, 1, 0, 0, 0))
    f1, f2 = A.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
    qd_north = f2 / np.linalg.norm(f2)
    East, North= qd_north[0], qd_north[1]
    Gridproj= CS.CSprojection((lon_centre, lat_centre), [East, North])
    node_grid=CS.CSgrid(Gridproj, 3700, 2200, 50., 50.)
    Evalproj= CS.CSprojection((lon_centre-0.23, lat_centre+0.23), [East, North])
    eval_grid= CS.CSgrid(Evalproj,  3300, 1800, 50., 50.)
    x, y = node_grid.xi, node_grid.eta
    xeval, yeval= eval_grid.xi, eval_grid.eta
    node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
    lat,lon= eval_grid.lat.flatten(), eval_grid.lon.flatten()
    """Interpolate"""    
    Btheta, Bphi, Br, MagLon, MagLat= Example.ExampleData()
    poles= SECS(node_lons, node_lats, lon, lat)
    poles.Fitting_Matrix(MagLon, MagLat)
    M= poles.Amplitude(Btheta, Bphi, Br)
    fig=plt.figure()
    ax= fig.add_subplot(121)
    ax2= fig.add_subplot(122)
    s= ax.scatter(x.flatten(), y.flatten(), c= M, cmap='seismic')
    s2= ax2.scatter(xeval.flatten(), yeval.flatten(), c=interpolation2d(np.unique(x), np.unique(y), np.vstack([M[x.flatten()==val] for val in np.unique(x)]), xeval.flatten(), yeval.flatten()), cmap='seismic', vmin=np.min(M), vmax=np.max(M))
    ax2.set_ylim(ax.get_ylim())
    ax2.set_xlim(ax.get_xlim())
    cbar= fig.colorbar(s, ax=ax)
    cbar2= fig.colorbar(s2, ax=ax2)
#%% Cubed Sphere Cartopy Example
if __name__=='__main__':
    from SECpy import SECS, Example
    from SECpy import PlottingTools as PT
    import cubedsphere as CS
    from apexpy import Apex
    import datetime as dt
    import cartopy.crs as ccrs
    """Grid"""
    lon_centre= 17.7
    lat_centre= 68.1
    A = Apex(date=dt.datetime(2008, 6, 1, 0, 0, 0))
    f1, f2 = A.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
    qd_north = f2 / np.linalg.norm(f2)
    East, North= qd_north[0], qd_north[1]
    Gridproj= CS.CSprojection((lon_centre, lat_centre), [East, North])
    node_grid=CS.CSgrid(Gridproj, 3700, 2200, 50., 50.)
    Evalproj= CS.CSprojection((lon_centre-0.23, lat_centre+0.23), [East, North])
    eval_grid= CS.CSgrid(Evalproj,  3300, 1800, 50., 50.)
    x, y = node_grid.xi, node_grid.eta
    xeval, yeval= eval_grid.xi, eval_grid.eta
    node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
    lat,lon= eval_grid.lat.flatten(), eval_grid.lon.flatten()
    """Figure"""
    proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
    fig=plt.figure()
    ax= fig.add_subplot(121, projection=proj)
    ax2= fig.add_subplot(122, projection=proj)
    PT.features(ax, coastlines_only=True)
    PT.features(ax2, coastlines_only=True)
    """Interpolate"""    
    Btheta, Bphi, Br, MagLon, MagLat= Example.ExampleData()
    poles= SECS(node_lons, node_lats, lon, lat)
    poles.Fitting_Matrix(MagLon, MagLat)
    M= poles.Amplitude(Btheta, Bphi, Br)
    s= ax.scatter(node_lons, node_lats, c= M, cmap='seismic', transform=ccrs.PlateCarree())
    s2= ax2.scatter(lon, lat, 
                    c=interpolation2d(np.unique(x), np.unique(y), np.vstack([M[x.flatten()==val] for val in np.unique(x)]), 
                                      xeval.flatten(), yeval.flatten()), cmap='seismic', 
                    vmin=np.min(M), vmax=np.max(M), transform=ccrs.PlateCarree())
    ax2.set_ylim(ax.get_ylim())
    ax2.set_xlim(ax.get_xlim())
    cbar= fig.colorbar(s, ax=ax)
    cbar2= fig.colorbar(s2, ax=ax2)




