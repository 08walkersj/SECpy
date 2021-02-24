#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 17:05:55 2020
Module for the purpose of using Spherical Elementary Currents (SECs).
The equations referenced in this module are from Vanhamäki, H., & 
Juusola, L. (2020). Introduction to Spherical Elementary Current Systems. 
In Ionospheric Multi-Spacecraft Analysis Tools (pp. 5–33). Springer 
International Publishing. https://doi.org/10.1007/978-3-030-26732-2_2
@author: Simon Walker
"""


"""Module Imports"""
import numpy as np

"""Useful Constants"""
pi = np.pi
RE= 6.371E6
u0= 1.23E-6

"""Errors"""
class ArgumentError(Exception):
     pass
class ModelError(Exception):
    pass

"""Definitions"""
def Deg2Rad(angle):
    """
    Parameters
    ----------
    angle : list/numpy.ndarray/int/float (Degrees)
        Angle to be converted to radians.

    Returns
    -------
    angle : input format (Radians)
        Angle in radians.
    """
    
    return angle*(pi/180)
def Rad2Deg(angle):
    """
    Parameters
    ----------
    angle : list/numpy.ndarray/int/float (Radians)
        Angle to be converted to degrees.

    Returns
    -------
    angle : input format (Degrees)
        Angle in degrees.
    """
    return angle*(180/pi)
def theta(lat_value, long_value, latitude, longitude):
    #Equation 2.16 (Vanhamäki 2020)
    """
    Parameters
    ----------
    lat_value : numpy.ndarray (Degrees)
        Latitude position of the SEC pole(s).
    long_value : numpy.ndarray (Degrees)
        Longitude position of the SEC pole(s).
    latitude : numpy.ndarray (Degrees)
        Latitude position of the evaluation point(s).
    longitude : numpy.ndarray (Degrees)
        Longitude position of the evaluation point(s).

    Returns
    -------
    Colatitude : numpy.ndarray (Radians)
        Theta prime, the colatitude of the evaluation point(s) from the SEC pole(s).
    """
    #theta_P, phi_P: colatitude of SEC pole(s), longitude of SEC pole(s)
    theta_P =Deg2Rad(90- lat_value)
    phi_P =Deg2Rad(long_value)
    #theta0, phi0: colatitude of evaluation point(s), 
    #               longitude of evaluation point(s)
    theta0 =np.vstack(Deg2Rad(90- latitude))
    phi0 =np.vstack(Deg2Rad(longitude))
    cos_theta =np.cos(theta0)*np.cos(theta_P) +np.sin(theta0)*np.sin(theta_P)*np.cos(phi_P- phi0)
    return np.arccos(cos_theta)
def CosC(theta_P, theta0, thetap):
    #Equation 2.19 (Vanhamäki 2020)
    """
    Parameters
    ----------
    theta_P : numpy.ndarray (Radians)
        Colatitude of the SEC pole(s).
    theta0 : numpy.ndarray (Radians)
        Colatitude of the evaluation point(s).
    thetap : numpy.ndarray (Radians)
        Theta prime, the colatitude of the evaluation point(s) from the SEC pole(s).

    Returns
    -------
    CosC : numpy.ndarray
        Variable needed for the conversion of vectors from the local SEC pole system(s) to a global co-ordinate system.
    """
    A= np.cos(theta_P)- np.cos(np.vstack(theta0))*np.cos(thetap)
    B= np.sin(np.vstack(theta0))*np.sin(thetap)
    return A/B
def SinC(theta_P, thetap, phi_P, phi0):
    #Equation 2.20 (Vanhamäki 2020)
    """
    Parameters
    ----------
    theta_P : numpy.ndarray (Radians)
        Colatitude of the SEC pole(s).
    thetap : numpy.ndarray (Radians)
        Theta prime, the colatitude of the evaluation point(s) from the SEC pole(s).
    phi_P : numpy.ndarray (Radians)
        Longitude of the SEC pole(s).
    phi0 : numpy.ndarray (Radians)
        Longitude of the evaluation point(s).

    Returns
    -------
    SinC : numpy.ndarray
        Variable needed for the conversion of vectors from the "local" SEC pole system(s) to a global co-ordinate system.
    """
    A= np.sin(theta_P)*np.sin(phi_P- np.vstack(phi0))
    B= np.sin(thetap)
    return A/B
def Local2Global(theta_P, theta0, thetap, phi_P, phi0, VectorPhi, VectorTheta):
    #Based on a rotation of equation 2.17 and 2.18 (Vanhamäki 2020)
    """
    Parameters
    ----------
    theta_P : numpy.ndarray (Radians)
        Colatitude of the SEC pole(s).
    theta0 : numpy.ndarray (Radians)
        Colatitude of the evaluation point(s).
    thetap : numpy.ndarray (Radians)
        Theta prime, the colatitude of the evaluation point(s) from the SEC pole(s).
    phi_P : numpy.ndarray (Radians)
        Longitude of the SEC pole(s).
    phi0 : numpy.ndarray (Radians)
        Longitude of the evaluation point(s).
    VectorPhi : numpy.ndarray (Vector component)
        Component of the vector in the direction of the phi unit vector in "local" 
        SEC pole system(s) to be converted into a projection in a global system.
    VectorTheta : numpy.ndarray (Vector component)
        Component of the vector in the direction of the theta unit vector in "local" 
        SEC pole system(s) to be converted into a projection in a global system.

    Returns
    -------
    new_VectorTheta : numpy.ndarray (Vector component)
        Theta component of vector converted into a global system.
    new_VectorPhi : numpy.ndarray (Vector component)
        Phi component of vector converted into a global system.
    """
    CosC_= CosC(theta_P, theta0, thetap)
    SinC_= SinC(theta_P, thetap, phi_P, phi0)
    new_VectorTheta= VectorTheta*CosC_ + VectorPhi*SinC_
    new_VectorPhi=-VectorTheta*SinC_ + VectorPhi*CosC_
    return new_VectorTheta, new_VectorPhi
def BrDF(thetap, r, R):
    #Equation 2.13 (Vanhamäki 2020)
    """
    Parameters
    ----------
    thetap : numpy.ndarray (Radians)
        Theta prime, the colatitude of the evaluation point(s) from the SEC pole(s).
    r : float (Meters)
        Radial distance from the centre of the Earth to the evaluation point(s).
    R : float (Meters)
        Radial distance from the centre of the Earth to the SEC pole(s)/ ionospheric layer.

    Raises
    ------
    ValueError
        Error if you try to evaluate the magnetic field at the radius of the SEC pole(s) 
        as there shouldn't be a reason to do this.

    Returns
    -------
    Matrix : numpy.ndarray (Tesla)
        When the matrix is multiplied by the SEC pole amplitude(s)/ 
        scalar(s) will give the radial magnetic field component at the evaluation 
        point(s) due to the SEC pole(s) in the "local" SEC pole system(s).
    """
    s= min([r,R])/max([r,R])
    A= u0/(4*pi*r)
    if r<R:
        B= (1+ s**2 -2*s*np.cos(thetap))**-0.5
        B-=1
    elif r>R:
        B= s/np.sqrt(1+s**2 -2*s*np.cos(thetap))
        B-= s
    else:
        raise ValueError("Invalid radius inputs shouldn't be equal")
    return A*B
def BthetaDF(thetap, r, R):
    #Equation 2.14 (Vanhamäki 2020)
    """
    Parameters
    ----------
    thetap : numpy.ndarray (Radians)
        Theta prime, the colatitude of the evaluation point(s) from the SEC pole(s).
    r : float (Meters)
        Radial distance from the centre of the Earth to the evaluation point(s).
    R : float (Meters)
        Radial distance from the centre of the Earth to the SEC pole(s)/ ionospheric layer.

    Raises
    ------
    ValueError
        Error if you try to evaluate the magnetic field at the radius of the SEC pole(s) 
        as there shouldn't be a reason to do this.

    Returns
    -------
    Matrix : numpy.ndarray (Tesla)
        When the matrix is multiplied by the SEC pole amplitude(s)/ 
        scalar(s) will give the theta component of the magnetic field at the evaluation 
        point(s) due to the SEC pole(s) in the "local" SEC pole system(s).
    """
    s= min([r,R])/max([r,R])
    A= -u0/(4*pi*r*np.sin(thetap))
    if r<R:
        B= (s-np.cos(thetap))*(1+ s**2 -2*s*np.cos(thetap))**-0.5
        B+= np.cos(thetap)
    elif r>R:
        B=(1-s*np.cos(thetap))/np.sqrt(1+s**2 -2*s*np.cos(thetap))
        B-=1
    else:
        raise ValueError('Invalid radius inputs cannot be equal')
    return A*B
def Vdf(thetap, R, theta0= 0):
    #Equation 2.8 (Vanhamäki 2020)
    """
    Parameters
    ----------
    thetap : list/numpy.ndarray (Radians)
        Colatitude of the evaluation point(s) in local SEC pole system(s).
    R : float (Meters)
        Radius from the centre of the Earth to the evaluation point(s).
    theta0 : float (Radians)
        Singularity limit to be used when the evaluation point(s) are too close to the SEC pole(s). The method is based on equations 2.44 from Vanhamäki 2020
        The default is 0.

    Returns
    -------
    Current : numpy.ndarray (Amps)
        Divergence free current around the SEC pole system(s) due 
        to its nature only in the theta unit vector direction in 
        the local co-ordinate system and only once multiplied by the SEC pole amplitude(s).
    """

    A= 1/(4*pi*R)
    if theta0>0:
        Current= np.zeros(thetap.shape)
        Current[thetap<theta0]= np.tan(thetap[thetap<theta0]/2)*(np.tan(theta0/2))**-2
        Current[np.invert(thetap<theta0)]= ((np.tan(thetap[np.invert(thetap<theta0)]/2))**-1)
        return A*Current
    else:
        return A*((np.tan(thetap/2))**-1)
def GridCheck(longitude, latitude, MagLon, MagLat, GridSpacing_km, limit=20, evaluation_altitude=RE, theta_prime=None):
    """
    Parameters
    ----------
    longitude : numpy.ndarray (Degrees)
        Longitude of the SEC poles to check
    latitude : numpy.ndarray (Degrees)
        Latitude of the SEC poles to check
    MagLon : numpy.ndarray (Degrees)
        Longitude of the magnetometers that are being used.
    MagLat : numpy.ndarray (Degrees)
        Latitude of the magnetometers that are being used.
    GridSpacing_km : int/float (kilometers)
        Average spacing between the poles in kilometers.
    limit : int/float, optional
        Percentage of the grid spacing that is the closest the SEC poles can be to magnetometers. The default is 20.
    evaluation_altitude: int/float (kilometers)
        Altitude of the evaluation points/ magnetometers. The default is RE.
    Returns
    -------
    Longitude : numpy.ndarray (Degrees)
        Longitude values of the poles that are problematic.
    Latitude : numpy.ndarray (Degrees)
        Latitude values of the poles that are problematic.
    Index: numpy.ndarray
        Index of the position of the problem poles. Index[0] is the problem measurements, Index[1] is the problem poles

    """
    print('Analysing grid')
    if theta_prime is None:
        theta_prime= theta(latitude, longitude, MagLat, MagLon)
    if len(np.where(GridSpacing_km *1e3 *limit/100>theta_prime*evaluation_altitude)[0])==0:
        print('Grid within requirements')
        return [None]*3
    else:
        print('Grid not optimal returning longitude and latitude values of the poles that are a problem and their index')
        I= np.where(GridSpacing_km *1e3 *limit/100>theta_prime*evaluation_altitude)
        #I[0] is the problem measurements, I[1] is the problem poles
        return longitude[I[1]], latitude[I[1]], I
def ImproveGrid(lon_centre, lat_centre, MagLon, MagLat, GridSpacing_km, length, height, evaluation_altitude=RE, date=None, east=None, north=None, limit=20, movement_limit=5):
    """
    Parameters
    ----------
    lon_centre : int/float (Degrees)
        longitude of the centre of the node grid.
    lat_centre : int/float ((Degrees))
        lattitude of the centre of the node grid.
    MagLon : numpy.ndarray (Degrees)
        Longitude of the magnetometers.
    MagLat : numpy.ndarray (Degrees)
        Latitude of the magnetometers.
    GridSpacing_km : float
        Spacing between nodes in km for the node grid.
    length : int/float
        length of the node grid.
    height : int/float
        height of the node grid.
    date : datetime.datetime
        Date being used.
    limit : int/float, optional
        Percentage of grid spacing that will provide a limit of closeness of magnetometers with the pole nodes. The default is 20.
    evaluation_altitude: int/float (kilometers)
        Altitude of the evaluation points/ magnetometers. The default is RE.
    Returns
    -------
        Change: numpy.ndarray (Degrees)
            Recommended change in the central position of the grid ([[lon_change, lat_change]]).
        Index: numpy.ndarray
            Index for problem poles for recommended grid (empty if no problem poles).
        Num_Bad_Poles: numpy.ndarray
            number of problem poles with each adjustment.
        Changes_Tried: numpy.ndarray
            all adjustments tried.
        Index: numpy.ndarray:
            Index where for each adjusted grid there is a problem pole.
    """
    try:
        import progressbar
        bar=True
    except ModuleNotFoundError:
        bar=False
    if bar:
        prog_bar = progressbar.ProgressBar(max_value=len(np.arange(-5,5,0.1))**2) 
        pbi=0
    if date is not None:
        from apexpy import Apex
    num_minimum=[]
    change=[]
    index=[]
    # from pysymmetry.utils.CSgrid import CSgrid
    from pysymmetry.utils import cubedsphere as CS
    if date is not None:
        A = Apex(date=date)
        f1, f2 = A.basevectors_qd(lat_centre, lon_centre, 0, coords = 'geo')
        qd_north = f2 / np.linalg.norm(f2)
        East, North= qd_north[0], qd_north[1]
    else:
        East, North= east, north
    # node_grid= CSgrid((lon_centre, lat_centre), [East,North], height, length, GridSpacing_km, GridSpacing_km)
    Gridproj= CS.CSprojection((lon_centre, lat_centre), [East, North])
    node_grid=CS.CSgrid(Gridproj, height, length, GridSpacing_km, GridSpacing_km)
    node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten() 
    thetap= theta(node_lats, node_lons, MagLat, MagLon)
    print('\n number of problem poles initially ',len(np.where(GridSpacing_km* 1e3 *limit/100 > thetap*evaluation_altitude)[0]))
    if len(np.where(GridSpacing_km* 1e3 *limit/100 > thetap*evaluation_altitude)[0])==0:
        print('Grid already optimal')
        if bar:
            prog_bar.finish()
        return [None]*5
    print('Looking for grid improvements')
    for dc in np.arange(-movement_limit, movement_limit, 0.1): 
         for d in np.arange(-movement_limit, movement_limit, 0.1):
             if date is not None:
                 A = Apex(date=date)
                 f1, f2 = A.basevectors_qd(lat_centre+dc, lon_centre+d, 0, coords = 'geo')
                 qd_north = f2 / np.linalg.norm(f2)
                 East, North= qd_north[0], qd_north[1]
             else:
                 East, North= east, north
             Gridproj= CS.CSprojection((lon_centre+d, lat_centre+dc), [East, North])
             node_grid=CS.CSgrid(Gridproj, height, length, GridSpacing_km, GridSpacing_km)
             node_lons, node_lats= node_grid.lon.flatten(), node_grid.lat.flatten()
             thetap= theta(node_lats, node_lons, MagLat, MagLon)
             num_minimum.append(len(np.where(GridSpacing_km* 1e3 *limit/100 > thetap*evaluation_altitude)[0]))
             change.append([d,dc])
             index.append([np.where(GridSpacing_km* 1e3 *limit/100 > thetap*evaluation_altitude)])
             if bar:
                 pbi+=1
                 prog_bar.update(pbi)
    if bar:
        prog_bar.finish()
    num_minimum= np.array(num_minimum)
    change= np.array(change)
    index= np.array(index)
    I=np.where(num_minimum== np.min(num_minimum))
    I2= np.where(np.sum(abs(change[I]), axis=1)==np.min(np.sum(abs(change[I]), axis=1)))
    print('Completed! Information returned')
    print('Number of problem poles with this optimisation ',
          num_minimum[I][I2])
    return change[I][I2], index[I][I2], num_minimum, change, index
"""Classes"""
class PlottingTools():
    #Plotting tools that can be used alongisde the SECS analysis
    #Cartopy is used to allow plotting of SECS models onto maps
    def Geocentric_to_PlateCarree_vector_components(east, north, latitude):
        """ convert east north vector components to Plate Carree projection 
    
            parameters
            ----------
            east : array-like
                eastward components
            north : array-like
                westward components
            latitude : array-like
                latitude of each vector
    
            returns
            -------
            east, north in Plate Carree projection
            Preserveres shape and norm
    
        """
        magnitude = np.sqrt(east**2 + north**2)
    
        east_pc = east / np.cos(latitude * np.pi / 180)
    
        magnitude_pc = np.sqrt(east_pc**2 + north**2)
    
        east_pc  = east_pc * magnitude / magnitude_pc
        north_pc = north * magnitude / magnitude_pc
    
        return east_pc, north_pc
    def features(axis, coastlines_only=False, features_list=False, features_kwargs=None):
        """
        Parameters
        ----------
        axis : matplotlib.axes._subplots.AxesSubplot
            Matplotlib axis object that map features are to be added to.

        Returns
        -------
        None.

        """
        import cartopy
        axis.coastlines(resolution='50m', zorder=200)
        if not coastlines_only and not features_list:
            axis.add_feature(cartopy.feature.OCEAN, facecolor='aqua')
            # axis.add_feature(cartopy.feature.COASTLINE, zorder=200)
            
            axis.add_feature(cartopy.feature.BORDERS, linestyle=':')
            axis.add_feature(cartopy.feature.LAKES, alpha=0.5, facecolor= 'mediumseagreen')
            axis.add_feature(cartopy.feature.RIVERS, edgecolor= 'mediumseagreen')
            axis.add_feature(cartopy.feature.LAND, facecolor='green')
        elif coastlines_only and features_list:
            raise ArgumentError('Either coast lines only or a list of features both cannot be inputs')
        elif not coastlines_only and features_list:
            for i, Feature in enumerate(features_list):
                if features_kwargs is None:
                    if Feature.lower()=='borders':
                        axis.add_feature(cartopy.feature.BORDERS, linestyle=':')
                    elif Feature.lower()=='lakes':
                        axis.add_feature(cartopy.feature.LAKES, alpha=0.5, facecolor= 'mediumseagreen')
                    elif Feature.lower()=='rivers':
                        axis.add_feature(cartopy.feature.RIVERS, edgecolor= 'mediumseagreen')
                    elif Feature.lower()=='land':
                        axis.add_feature(cartopy.feature.LAND, facecolor='green')
                    elif Feature.lower()== 'ocean':
                        axis.add_feature(cartopy.feature.OCEAN, facecolor='aqua')
                    elif Feature.lower()=='coastlines':
                        axis.coastlines(resolution='50m', zorder=200)
                else:
                    exec('axis.add_feature(cartopy.feature.'+Feature+',**features_kwargs[i])')
class Example():
    def ExampleGrid():
        """
        Returns
        -------
        node_lons : numpy.ndarray (Degrees)
            SEC pole node longitudes.
        node_lats : numpy.ndarray (Degrees)
            SEC pole node lattitudes.
        lon : numpy.ndarray (Degrees)
            Evaluated node longitudes.
        lat : numpy.ndarray (Degrees)
            Evaluated node lattitudes.
        """
        dlat=0.6
        dlon=1.4
        
        a1= round(58.26/dlat)*dlat
        a2= round(78.92/dlat)*dlat
        b1=round(4.84/dlon)*dlon
        b2= round(36.08/dlon)*dlon
        
        SECSlat,SECSlon = np.meshgrid(np.arange(a1-4.5*dlat,a2+4.5*dlat,dlat), np.arange(b1-4.5*dlon,b2+4.5*dlon,dlon))
        node_lons, node_lats= SECSlon.flatten(),SECSlat.flatten()
        lat,lon = np.meshgrid(np.arange(a1,a2,dlat), np.arange(b1,b2,dlon))
        lat,lon= lat.flatten(),lon.flatten()
        return node_lons, node_lats, lon, lat
    def ExampleData():
        """
        Returns
        -------
        Btheta : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the theta unit vector from the magnetometer(s).
        Bphi : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the phi unit vector from the magnetometer(s).
        Br : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the radial unit vector from the magnetometer(s).
        MagLon : numpy.ndarray (Degrees)
            Longitudinal location of the magnetometers.
        MagLat : numpy.ndarray (Degrees)
            Latitudinal location of the magnetometers.
        """
        #Bnames=['TAR','KAR','UPS','NUR','SOL','DOB','HAN','MEK','OUJ','LYC','RVK','RAN','DON','JCK','PEL','SOD','KIR','MUO','ABK','IVA','KIL','AND','MAS','TRO','KEV','SOR','NOR','BJN','HOP','HOR','LYR','NAL','SUW','BRZ']
        #Date= 23:34 08/11/2015
        BX = -np.array([-14.1,   -6.2,  -11.7,  -32.9,  -21.2,  -43.6,  -61.1,  -85.7, -159.3, -199.9, -336.5, -267.9, -560.3, -465.7, -436.0, -483.1, -673.7, -667.4, -789.5, -688.4, -784.5, -599.0, -714.7, -534.3, -674.2, -321.1, -251.1,   -3.7,   13.0,   37.2,   37.4,   57.3,   -7.1,   -8.2])
        BY = np.array([ -6.2,   52.1,   11.2,   -9.3,   61.0,   37.1,  -25.6,  -31.1,  -40.6,   -3.7,  109.5,  -37.2,  162.0,  105.4,    3.1,  -14.4,   84.1,  109.1,  133.1,   83.7,  158.8,  179.4,  115.1,   78.0,  121.7,   52.2,   46.6,   20.4,   13.6,   17.5,   10.1,   10.2,    6.6,  -0.9])
        BZ = -np.array([-48.2,  -66.9,  -79.6, -101.7, -134.4, -176.1, -144.1, -155.4, -221.3, -255.5, -306.7, -284.0, -284.9, -385.3, -340.3, -270.8, -193.8, -217.3,  -87.3, -149.6,  228.2,  510.4,  376.7,  480.7,  384.6,  492.0,  413.9,  247.0,  180.7,  176.3,  154.6,  159.3,  -23.9,  -32.2])
        MagLat=  np.array([58.26, 59.21, 59.90, 60.50, 61.08, 62.07, 62.30, 62.77, 64.52, 64.61, 64.94, 65.90, 66.11, 66.40, 66.90, 67.37, 67.83, 68.02, 68.35, 68.56, 69.02, 69.30, 69.46, 69.66, 69.76, 70.54, 71.09, 74.50, 76.51, 77.00, 78.20, 78.92, 54.01, 56.17])
        MagLon=  np.array([26.46,  5.24, 17.35, 24.65,  4.84,  9.12, 26.65, 30.97, 27.23, 18.75, 10.98, 26.41, 12.50, 16.98, 24.08, 26.63, 20.42, 23.53, 18.82, 27.29, 20.79, 16.03, 23.70, 18.94, 27.01, 22.22, 25.79, 19.20, 25.01, 15.60, 15.82, 11.95, 23.18, 24.86])
        return BX*1e-9, BY*1e-9, BZ*1e-9, MagLon, MagLat

class SECS:
    def __init__(self, longitude, latitude, eval_longitude, eval_latitude, current_radius=RE+110E3, 
                 mode='default', image_current_radius= None, longitude2=None, latitude2=None):
        """
        Parameters
        ----------
        longitude : list/numpy.ndarray (Degrees)
            Longitude position of the SEC pole(s).
        latitude : list/numpy.ndarray (Degrees)
            Latitude position of the SEC pole(s).
        eval_longitude : list/numpy.ndarray (Degrees)
            Longitude position of the evaluation point(s).
        eval_latitude : list/numpy.ndarray (Degrees)
            Latitude position of the evaluation point(s).
        current_radius : Meters, optional
            The radial distance, from the centre of the Earth, to the current. 
            The default is RE+110E3, the F-layer of the Ionosphere.
        mode : str, optional
            The mode to be used for the SEC model, options are:

            'default'→ which can involve just one set of poles or if a second radial distance is defined
            then a second set of poles will be made (longitude2 and latitude2 can be used to define 
            the location of these poles else the poles will be mirror the first set), the intention 
            of this is to create a second set of poles in the ground at the location of the 
            telluric currents to adjust the model for their influence on the magnetometer measurements.
            
            'image'→ which uses the method described in Juusola, L., Kauristie, K., Vanhamäki, H., 
            Aikio, A., & van de Kamp, M. (2016). Comparison of auroral ionospheric and field-aligned 
            currents derived from Swarm and ground magnetic field measurements. Journal of Geophysical 
            Research A: Space Physics, 121(9), 9256–9283. https://doi.org/10.1002/2016JA022961" 
            This requires just ionospheric poles to be defined and a second set of poles will be infered
            in the ground with a defined relationship with the ionospheric poles this will create a fitting
            matrix that is the same size as if there were just ionospheric poles.
            
            The default is 'default'.
        image_current_radius : Meters, optional
            The radial distance, from the centre of the Earth, to the telluric or image current (for telluric poles). The default is None.
        longitude2 : list/numpy.ndarray (Degrees), optional
            Longitude positon of the telluric SEC pole(s). The default is None.
        latitude2 : list/numpy.ndarray (Degrees), optional
            Latitude position of the telluric SEC pole(s). The default is None.

        Raises
        ------
        ArgumentError
            Will be raised if argument inputs are the wrong type or incorrect.
        
        New Attributes
        --------------
        mode : str
            Model mode.
        R : float (Meters)
            The radial distance, from the centre of the Earth, to the current. 
        RC : float (Meters)
            The radial distance, from the centre of the Earth, to the telluric or image current (for telluric poles).
        (telluric_ /ionospheric_) pole_longitude : numpy.ndarray (Degrees)
            Longitude position of the SEC pole(s).
        (telluric_ /ionospheric_) pole_latitude : numpy.ndarray (Degrees)
            Latitude position of the SEC pole(s).
        eval_longitude : numpy.ndarray (Degrees)
            Longitude position of the evaluation point(s).
        eval_latitude : numpy.ndarray (Degrees)
            Latitude position of the evaluation point(s).
        (telluric_ /ionospheric_) eval_theta_prime : Radians
            The colatitude of the evaluation point(s) from the SEC pole(s).
        Returns
        -------
        None.

        """
        self.mode= mode.lower()
        try:
            self.R= float(current_radius)
            if image_current_radius is not None:
                self.Rc= float(image_current_radius)
            else:
                self.Rc= image_current_radius
        except ValueError:
            raise ArgumentError('Invalid input, data type of current_radius and'+ 
                             'image_current_radius must be of type int or float'+
                             'except are types: ' +str(type(current_radius)) + ' and '
                             +str(type(image_current_radius)))
        if self.Rc is not None and self.mode== 'default' and (longitude2 is None or latitude2 is None):
            print('Creating independent poles to account for telluric currents at image_current_radius. \n'+
                  'if desired method is using image currents that have a specified realtionship with the poles in the ionosphere'+
                  "then set mode='image.'"+ ' Method taken from ' +
                  "Juusola, L., Kauristie, K., Vanhamäki, H., Aikio, A., & van de Kamp, M. (2016). Comparison of auroral ionospheric and field-aligned currents derived from Swarm and ground magnetic field measurements. Journal of Geophysical Research A: Space Physics, 121(9), 9256–9283. https://doi.org/10.1002/2016JA022961")
            print('\n longitude and latitude for ground poles are not defined so copying values from ionospheric poles')
            longitude2, latitude2= longitude, latitude
        if self.mode== 'image' and self.Rc is None:
            raise ArgumentError('image mode selected but image_current_radius is not defined')
        if not ((np.array([str(type(longitude)), str(type(longitude)), str(type(eval_longitude)), str(type(eval_latitude))])== "<class 'list'>") | (np.array([str(type(longitude)), str(type(longitude)), str(type(eval_longitude)), str(type(eval_latitude))])== "<class 'numpy.ndarray'>")).any():
            raise ArgumentError('SEC pole longitude(s), latitude(s), evaluation longitude(s) or evaluation latitudes(s) are of the wrong type must be list or array')
        if self.mode=='default' and self.Rc is None:
            self.pole_longitude= np.array(longitude)
            self.pole_latitude= np.array(latitude)
        elif self.mode=='default':
            self.ionospheric_pole_longitude= np.array(longitude)
            self.ionospheric_pole_latitude= np.array(latitude)
            self.telluric_pole_longitude= np.array(longitude2)
            self.telluric_pole_latitude= np.array(latitude2)
        elif self.mode=='image':
            self.ionospheric_pole_longitude, self.telluric_pole_longitude= [np.array(longitude)]*2
            self.ionospheric_pole_latitude, self.telluric_pole_latitude= [np.array(latitude)]*2
        else:
            raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
        self.eval_longitude= eval_longitude
        self.eval_latitude= eval_latitude
        #Calculate theta prime of the evaluation points
        try:
           self.eval_theta_prime=theta(self.pole_latitude, self.pole_longitude, eval_latitude, eval_longitude)
        except AttributeError:
            if self.mode=='default':
                self.eval_ionospheric_theta_prime=theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
                self.eval_telluric_theta_prime=theta(self.telluric_pole_latitude, self.telluric_pole_longitude, eval_latitude, eval_longitude)
            elif self.mode=='image':
                self.eval_theta_prime=theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
            else:
                raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
    def Fitting_Matrix(self, magnetometer_longitude, magnetometer_latitude, cond=0.05**2, eval_radius= RE):
        """
        Parameters
        ----------
        magnetometer_longitude : numpy.ndarray (Degrees)
            Longitude position of the magnetometer(s).
        magnetometer_latitude : numpy.ndarray (Degrees)
            Latitude position of the magnetometer(s).
        cond : float, optional
            Condition for Singular Value decomposition, singlular values less then cond*max(singular values) are ignored.
            The default is 0.05**2.
        eval_radius : float (Meters), optional
            Radial distance from the centre of the Earth to the magnetometer(s). The default is RE.

        Raises
        ------
        ArgumentError
            Will be raised if argument inputs are the wrong type or incorrect.

        Returns
        -------
        G : numpy.ndarray (Matrix)
            Matrix comprising of properties related to the magnetometer(s) and SEC pole(s), 
            when multiplied by the SEC pole amplitude(s) would match the magnetometer measurement(s)
            if the model is perfectly accurate at the magnetometer locations.
        invG : numpy.ndarray (Matrix)
            Moore-Penrose/sudo inverse of the G matrix using singular value decomposition.

        New Attributes
        --------------
        fitting_matrix : numpy.ndarray (Matrix)
            Matrix that when multiplied by the magnetometer data will give the SEC pole amplitude(s).
        (telluric_ /ionospheric_) theta_prime : numpy.ndarray (Radians)
            The colatitude of the evaluation point(s) from the SEC pole(s).

        """
        #Calculate theta prime of the magnetometers
        try:
            self.theta_prime=theta(self.pole_latitude, self.pole_longitude, magnetometer_latitude, magnetometer_longitude)
        except AttributeError:
            if self.mode=='default':
                self.ionospheric_theta_prime=theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, magnetometer_latitude, magnetometer_longitude)
                self.telluric_theta_prime=theta(self.telluric_pole_latitude, self.telluric_pole_longitude, magnetometer_latitude, magnetometer_longitude)
            elif self.mode=='image':
                self.theta_prime=theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, magnetometer_latitude, magnetometer_longitude)
            else:
                raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
        #Calculate the magnetic field at the magnetometer locations with pole amplitude(s) of value 1
        try:
            G= np.zeros((3*len(magnetometer_latitude), len(self.pole_latitude)))
            G[len(magnetometer_latitude)*2:len(magnetometer_latitude)*3]=BrDF(self.theta_prime, eval_radius, self.R)
            G[0:len(magnetometer_latitude)], G[len(magnetometer_latitude):len(magnetometer_latitude)*2] = Local2Global(Deg2Rad(90- self.pole_latitude), Deg2Rad(90- magnetometer_latitude), self.theta_prime, Deg2Rad(self.pole_longitude), 
                       Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.theta_prime, eval_radius, self.R))
        except AttributeError:
            if self.mode=='default':
                G= np.zeros((3*len(magnetometer_latitude), len(self.ionospheric_pole_latitude)+len(self.telluric_pole_latitude)))
                #Ionospheric and telluric Br
                G[len(magnetometer_latitude)*2:len(magnetometer_latitude)*3]=np.append(BrDF(self.ionospheric_theta_prime, eval_radius, self.R), BrDF(self.telluric_theta_prime, eval_radius, self.Rc), axis=1)
                #Ionospheric Btheta and Bphi
                Gn, Ge= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.ionospheric_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                       Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.ionospheric_theta_prime, eval_radius, self.R))
                #Telluric Btheta and Bphi
                Gn2, Ge2= Local2Global(Deg2Rad(90- self.telluric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.telluric_theta_prime, Deg2Rad(self.telluric_pole_longitude), 
                       Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.telluric_theta_prime, eval_radius, self.Rc))
                G[0:len(magnetometer_latitude)], G[len(magnetometer_latitude):len(magnetometer_latitude)*2] = np.append(Gn, Gn2, axis=1), np.append(Ge, Ge2, axis=1)
            elif self.mode=='image':
                G= np.zeros((3*len(magnetometer_latitude), len(self.ionospheric_pole_latitude)))
                #Ionospheric and telluric Br
                G[len(magnetometer_latitude)*2:len(magnetometer_latitude)*3]=BrDF(self.theta_prime, eval_radius, self.R) -(self.R/self.Rc)* BrDF(self.theta_prime, eval_radius, (self.Rc**2)/self.R)
                #Ionospheric Btheta and Bphi
                Gn, Ge= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                             Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.theta_prime, eval_radius, self.R))
                #Telluric Btheta and Bphi
                Gn2, Ge2= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                             Deg2Rad(magnetometer_longitude), 0, -(self.R/self.Rc)*BthetaDF(self.theta_prime, eval_radius, (self.Rc**2)/self.R))
                G[0:len(magnetometer_latitude)], G[len(magnetometer_latitude):len(magnetometer_latitude)*2] = Gn+Gn2, Ge+Ge2
            else:
                raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
        from scipy.linalg import svd, inv
        #Singular value decomposition of GTG matrix
        U, s, V= svd(np.dot(G.T,G))
        r=len(s)
        #Applying condition
        I=np.where([s<cond*np.max(s)])
        if len(I[0])>0:
            s=s[:I[1][0]]
        U=U[:,:r]
        V=V[:,:r]
        D= np.zeros((len(s), len(s)))
        np.fill_diagonal(D, s)
        invD= np.zeros((r,r))
        invD[:len(s),:len(s)]= inv(D)
        invG= np.dot(np.dot(V.T,invD),U.T)
        #Calculating fitting matrix when multiplied by the magnetometer data will give the SEC pole amplitude(s)
        self.fitting_matrix= np.dot(invG, G.T)
        return G, invG
    def G_Matrix(self, magnetometer_longitude, magnetometer_latitude, eval_radius= RE):
        """
        Parameters
        ----------
        magnetometer_longitude : numpy.ndarray (Degrees)
            Longitude position of the magnetometer(s).
        magnetometer_latitude : numpy.ndarray (Degrees)
            Latitude position of the magnetometer(s).
        eval_radius : float (Meters), optional
            Radial distance from the centre of the Earth to the magnetometer(s). The default is RE.

        Raises
        ------
        ArgumentError
            Will be raised if argument inputs are the wrong type or incorrect.

        Returns
        -------
        G : numpy.ndarray (Matrix)
            Matrix comprising of properties related to the magnetometer(s) and SEC pole(s), 
            when multiplied by the SEC pole amplitude(s) would match the magnetometer measurement(s)
            if the model is perfectly accurate at the magnetometer locations.

        """
        #Calculate theta prime of the magnetometers
        try:
            self.theta_prime=theta(self.pole_latitude, self.pole_longitude, magnetometer_latitude, magnetometer_longitude)
        except AttributeError:
            if self.mode=='default':
                self.ionospheric_theta_prime=theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, magnetometer_latitude, magnetometer_longitude)
                self.telluric_theta_prime=theta(self.telluric_pole_latitude, self.telluric_pole_longitude, magnetometer_latitude, magnetometer_longitude)
            elif self.mode=='image':
                self.theta_prime=theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, magnetometer_latitude, magnetometer_longitude)
            else:
                raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
        #Calculate the magnetic field at the magnetometer locations with pole amplitude(s) of value 1
        try:
            G= np.zeros((3*len(magnetometer_latitude), len(self.pole_latitude)))
            G[len(magnetometer_latitude)*2:len(magnetometer_latitude)*3]=BrDF(self.theta_prime, eval_radius, self.R)
            G[0:len(magnetometer_latitude)], G[len(magnetometer_latitude):len(magnetometer_latitude)*2] = Local2Global(Deg2Rad(90- self.pole_latitude), Deg2Rad(90- magnetometer_latitude), self.theta_prime, Deg2Rad(self.pole_longitude), 
                       Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.theta_prime, eval_radius, self.R))
        except AttributeError:
            if self.mode=='default':
                G= np.zeros((3*len(magnetometer_latitude), len(self.ionospheric_pole_latitude)+len(self.telluric_pole_latitude)))
                #Ionospheric and telluric Br
                G[len(magnetometer_latitude)*2:len(magnetometer_latitude)*3]=np.append(BrDF(self.ionospheric_theta_prime, eval_radius, self.R), BrDF(self.telluric_theta_prime, eval_radius, self.Rc), axis=1)
                #Ionospheric Btheta and Bphi
                Gn, Ge= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.ionospheric_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                       Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.ionospheric_theta_prime, eval_radius, self.R))
                #Telluric Btheta and Bphi
                Gn2, Ge2= Local2Global(Deg2Rad(90- self.telluric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.telluric_theta_prime, Deg2Rad(self.telluric_pole_longitude), 
                       Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.telluric_theta_prime, eval_radius, self.Rc))
                G[0:len(magnetometer_latitude)], G[len(magnetometer_latitude):len(magnetometer_latitude)*2] = np.append(Gn, Gn2, axis=1), np.append(Ge, Ge2, axis=1)
            elif self.mode=='image':
                G= np.zeros((3*len(magnetometer_latitude), len(self.ionospheric_pole_latitude)))
                #Ionospheric and telluric Br
                G[len(magnetometer_latitude)*2:len(magnetometer_latitude)*3]=BrDF(self.theta_prime, eval_radius, self.R) -(self.R/self.Rc)* BrDF(self.theta_prime, eval_radius, (self.Rc**2)/self.R)
                #Ionospheric Btheta and Bphi
                Gn, Ge= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                             Deg2Rad(magnetometer_longitude), 0, BthetaDF(self.theta_prime, eval_radius, self.R))
                #Telluric Btheta and Bphi
                Gn2, Ge2= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- magnetometer_latitude), self.theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                             Deg2Rad(magnetometer_longitude), 0, -(self.R/self.Rc)*BthetaDF(self.theta_prime, eval_radius, (self.Rc**2)/self.R))
                G[0:len(magnetometer_latitude)], G[len(magnetometer_latitude):len(magnetometer_latitude)*2] = Gn+Gn2, Ge+Ge2
            else:
                raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
        return G
    def eval_G_Matrix_B(self, eval_radius= RE, eval_longitude=None, eval_latitude=None):
        """
        When multiplied by the SEC pole amplitudes these matrices will return the magnetic field at the evaluation points.

        Parameters
        ----------
        eval_radius : float (Meters), optional
            Radial distance from the centre of the Earth to the magnetometer(s). The default is RE.
        eval_longitude : numpy.ndarray (Degrees)
            Longitude position of the evaluation point(s). The default is set by user in set up.
        eval_latitude : numpy.ndarray (Degrees)
            Latitude position of the evaluation point(s). The default is set by user in set up.

        Raises
        ------
        ArgumentError
            Raised if not using values in set up and both longitude and latitude are being used.

        Returns (Mode dependent)
        -------
        → mode='default' with no telluric poles
        
        Gr, Ge, Gn : (numpy.ndarray, numpy.ndarray, numpy.ndarray)
            Radial, eastward, northward components of the G matrix at the evaluation point(s).
        
        → mode='image' or 'default' with telluric poles
        
        (Gr, Ge, Gn), (Gr2, Ge2, Gn2) : ((numpy.ndarray, numpy.ndarray, numpy.ndarray), (numpy.ndarray, numpy.ndarray, numpy.ndarray))
            G matrix that will show Ionospheric SEC pole contribution to the magnetic field measurements at the evaluation point(s) when multiplied by SEC pole amplitude,
            G matrix that will show Telluric SEC pole contribution to the magnetic field measurements at the evaluation point(s) when multiplied by SEC pole amplitude.

        """
        if eval_longitude is None and eval_latitude is None:
            eval_longitude= self.eval_longitude
            eval_latitude= self.eval_latitude
            try:
                eval_theta_prime= self.eval_theta_prime
            except AttributeError:
                if self.mode=='image':
                    eval_theta_prime=self.eval_ionospheric_theta_prime
                elif self.mode=='default':
                    ionospheric_eval_theta_prime= self.eval_ionospheric_theta_prime
                    telluric_eval_theta_prime= self.eval_telluric_theta_prime
        elif eval_longitude is None or eval_latitude is None:
            raise ArgumentError('Evaluation latitude and longitude must be defined if not using the values assigned initially')
        else:
            try:
                #Calculate theta prime of the evaluation points
                eval_theta_prime= theta(self.pole_latitude, self.pole_longitude, eval_latitude, eval_longitude)
            except AttributeError:
                if self.mode=='image':
                    eval_theta_prime= theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
                elif self.mode=='default':
                    ionospheric_eval_theta_prime= theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
                    telluric_eval_theta_prime= theta(self.telluric_pole_latitude, self.telluric_pole_longitude, eval_latitude, eval_longitude)
        if self.mode=='default':
            try: 
                self.pole_longitude
                #Radial component of the magnetic field at the evaluation point(s)
                GBr=BrDF(eval_theta_prime, eval_radius, self.R)
                #North and east component of the magnetic field at the evaluation point(s)
                GBn, GBe= Local2Global(Deg2Rad(90- self.pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.pole_longitude), 
                   Deg2Rad(eval_longitude), 0, BthetaDF(eval_theta_prime, eval_radius, self.R))
                
                return  GBr, GBe, GBn
            except AttributeError:
                #Radial component of the ionospheric contribution to the magnetic field at the evaluation point(s)
                GBr= BrDF(ionospheric_eval_theta_prime, eval_radius, self.R)
                #North and east component of the ionospheric contribution to the magnetic field at the evaluation point(s)
                GBn, GBe= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), ionospheric_eval_theta_prime, 
                                     Deg2Rad(self.ionospheric_pole_longitude), Deg2Rad(eval_longitude), 0, BthetaDF(ionospheric_eval_theta_prime, eval_radius, self.R))
                #Radial component of the telluric contribution to the magnetic field at the evaluation point(s)
                GBr2= BrDF(telluric_eval_theta_prime, eval_radius, self.Rc)
                #North and east component of the telluric contribution to the magnetic field at the evaluation point(s)            
                GBn2, GBe2= Local2Global(Deg2Rad(90- self.telluric_pole_latitude), Deg2Rad(90- eval_latitude), telluric_eval_theta_prime, Deg2Rad(self.telluric_pole_longitude), 
                   Deg2Rad(eval_longitude), 0, BthetaDF(telluric_eval_theta_prime, eval_radius, self.Rc))
                return (GBr, GBe, -GBn), (GBr2, GBe2, -GBn2)
        elif self.mode=='image':
            #Radial component of the ionospheric contribution to the magnetic field at the evaluation point(s)
            GBr= BrDF(eval_theta_prime, eval_radius, self.R)
            #North and east component of the ionospheric contribution to the magnetic field at the evaluation point(s)
            GBn, GBe= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
               Deg2Rad(eval_longitude), 0, BthetaDF(eval_theta_prime, eval_radius, self.R))
            #Radial component of the telluric contribution to the magnetic field at the evaluation point(s)
            GBr2= (-self.R/self.Rc)*BrDF(eval_theta_prime, eval_radius, (self.Rc**2)/self.R)
            #North and east component of the telluric contribution to the magnetic field at the evaluation point(s)
            GBn2, GBe2= Local2Global(Deg2Rad(90- self.telluric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.telluric_pole_longitude), 
               Deg2Rad(eval_longitude), 0, (-self.R/self.Rc)*BthetaDF(eval_theta_prime, eval_radius, (self.Rc**2)/self.R))
            return (GBr, GBe, -GBn), (GBr2, GBe2, -GBn2)
        else:
            raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
    def eval_G_Matrix_J(self, eval_radius=None, eval_longitude=None, eval_latitude=None, singularity_limit=0):
        """
        

        When multiplied by the SEC pole amplitudes these matrices will return the sheet current density at the evaluation points.

        Parameters
        ----------
        eval_radius : float (Meters), optional
            Radial distance from the centre of the Earth to the magnetometer(s). The default is set by user in set up.
        eval_longitude : numpy.ndarray (Degrees)
            Longitude position of the evaluation point(s). The default is set by user in set up.
        eval_latitude : numpy.ndarray (Degrees)
            Latitude position of the evaluation point(s). The default is set by user in set up.
        singularity_limit : TYPE, optional
            Limit that signifies a change in formula close to the formula close to the SEC pole. The method is based on equations 2.44 from Vanhamäki 2020.
            The default is 0.

        Raises
        ------
        ArgumentError
            Raised if not using values in set up and both longitude and latitude are being used.

        Returns
        -------
        East_Current_G : numpy.ndarray
            Eastward components of the G matrix at the evaluation point(s) for the sheet current.
        North_Current_G : numpy.ndarray
            Northward components of the G matrix at the evaluation point(s) for the sheet current.

        """
        singularity_limit= singularity_limit/self.R
        if eval_radius is None:
            eval_radius= self.R
        if eval_longitude is None and eval_latitude is None:
            eval_longitude= self.eval_longitude
            eval_latitude= self.eval_latitude
            try:
                eval_theta_prime= self.eval_theta_prime
            except AttributeError:
                eval_theta_prime=self.eval_ionospheric_theta_prime
        elif eval_longitude is None or eval_latitude is None:
            raise ArgumentError('Evaluation latitude and longitude must be defined if not using the values assigned initially')
        else:
            try:
                #Calculate theta prime of the evaluation points
                eval_theta_prime= theta(self.pole_latitude, self.pole_longitude, eval_latitude, eval_longitude)
            except AttributeError:
                eval_theta_prime= theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
        try:
            #North and east component of the current at the evaluation point(s)
            GCn, GCe=Local2Global(Deg2Rad(90- self.pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.pole_longitude), 
                                Deg2Rad(eval_longitude), Vdf(eval_theta_prime, eval_radius, theta0= singularity_limit), 0)
            return GCe, -GCn
        except AttributeError:
            if self.mode=='default':
                #North and east component of the current at the evaluation point(s)
                GCn, GCe=Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                                    Deg2Rad(eval_longitude), Vdf(eval_theta_prime, eval_radius, theta0= singularity_limit), 0)
                return GCe, -GCn
            elif self.mode=='image':
                #North and east component of the current at the evaluation point(s)
                GCn, GCe=Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                                    Deg2Rad(eval_longitude), Vdf(eval_theta_prime, eval_radius, theta0= singularity_limit), 0)
                return GCe, -GCn
        else:
            raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
    def Inverse_Matrix(self, matrix, cond=0.05**2):
        """
        Parameters
        ----------
        matrix : numpy.ndarray (Matrix)
            Matrix to be inverted using the same method in Fitting_Matrix.
        cond : float, optional
            Condition for Singular Value decomposition, singlular values less then cond*max(singular values) are ignored.
            The default is 0.05**2.

        Returns
        -------
        invG : numpy.ndarray (Matrix)
            Moore-Penrose/sudo inverse of the matrix using singular value decomposition.

        """
        from scipy.linalg import svd, inv
        #Singular value decomposition of GTG matrix
        U, s, V= svd(matrix)
        r=len(s)
        #Applying condition
        I=np.where([s<cond*np.max(s)])
        if len(I[0])>0:
            s=s[:I[1][0]]
        U=U[:,:r]
        V=V[:,:r]
        D= np.zeros((len(s), len(s)))
        np.fill_diagonal(D, s)
        invD= np.zeros((r,r))
        invD[:len(s),:len(s)]= inv(D)
        invG= np.dot(np.dot(V.T,invD),U.T)
        return invG
    def Magnetic_Field(self, Btheta, Bphi, Br, eval_radius= RE, eval_longitude=None, eval_latitude=None):
        """
        Parameters
        ----------
        Btheta : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the theta unit vector from the magnetometer(s).
        Bphi : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the phi unit vector from the magnetometer(s).
        Br : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the radial unit vector from the magnetometer(s).
        eval_radius : float (Meters), optional
            Radius from the centre of the Earth to the evaluation point(s). The default is RE.
        eval_longitude : numpy.ndarray (Degrees), optional
            Longitude position of the evaluation point(s). The default is None.
        eval_latitude : numpy.ndarray (Degrees), optional
            Latitude position of the evaluation point(s). The default is None.

        Raises
        ------
        ArgumentError
            Will be raised if argument inputs are the wrong type or incorrect.
        ModelError
            Will be rasied if there is no "fitting_matrix" that is an attribute of the class

        Returns (Mode dependent)
        -------
        → mode='default' with no telluric poles
        
        Br, Be, Bn : (numpy.ndarray, numpy.ndarray, numpy.ndarray) (Tesla)
            Radial, eastward, northward components of the magnetic field at the evaluation point(s).
        
        → mode='image' or 'default' with telluric poles
        
        (Br, Be, Bn), (Br2, Be2, Bn2) : ((numpy.ndarray, numpy.ndarray, numpy.ndarray), (numpy.ndarray, numpy.ndarray, numpy.ndarray)) (Tesla)
            Ionospheric SEC pole contribution to the magnetic field measurements at the evaluation point(s),
            Telluric SEC pole contribution to the magnetic field measurements at the evaluation point(s).

        """
        try:
            #Calculate the amplitudes of the SEC pole(s)
            m= np.dot(self.fitting_matrix, np.concatenate((Btheta, Bphi, Br)).T)
        except AttributeError:
            raise ModelError('Fitting matrix could not be found ensure you have successfully run the "Fitting_Matrix" function or assigned a "fitting_matrix" as an attribute to this class')
        if eval_longitude is None and eval_latitude is None:
            eval_longitude= self.eval_longitude
            eval_latitude= self.eval_latitude
            try:
                eval_theta_prime= self.eval_theta_prime
            except AttributeError:
                if self.mode=='image':
                    eval_theta_prime=self.eval_ionospheric_theta_prime
                elif self.mode=='default':
                    ionospheric_eval_theta_prime= self.eval_ionospheric_theta_prime
                    telluric_eval_theta_prime= self.eval_telluric_theta_prime
        elif eval_longitude is None or eval_latitude is None:
            raise ArgumentError('Evaluation latitude and longitude must be defined if not using the values assigned initially')
        else:
            try:
                #Calculate theta prime of the evaluation points
                eval_theta_prime= theta(self.pole_latitude, self.pole_longitude, eval_latitude, eval_longitude)
            except AttributeError:
                if self.mode=='image':
                    eval_theta_prime= theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
                elif self.mode=='default':
                    ionospheric_eval_theta_prime= theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
                    telluric_eval_theta_prime= theta(self.telluric_pole_latitude, self.telluric_pole_longitude, eval_latitude, eval_longitude)
        if self.mode=='default':
            try: 
                self.pole_longitude
                #Radial component of the magnetic field at the evaluation point(s)
                Br=m*BrDF(eval_theta_prime, eval_radius, self.R)
                #North and east component of the magnetic field at the evaluation point(s)
                Bn, Be= Local2Global(Deg2Rad(90- self.pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.pole_longitude), 
                   Deg2Rad(eval_longitude), 0, m*BthetaDF(eval_theta_prime, eval_radius, self.R))
                
                return  np.sum(Br, 1), np.sum(Be, 1), np.sum(-Bn, 1)
            except AttributeError:
                #Radial component of the ionospheric contribution to the magnetic field at the evaluation point(s)
                Br= m[:len(self.ionospheric_pole_latitude)]*BrDF(ionospheric_eval_theta_prime, eval_radius, self.R)
                #North and east component of the ionospheric contribution to the magnetic field at the evaluation point(s)
                Bn, Be= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), ionospheric_eval_theta_prime, 
                                     Deg2Rad(self.ionospheric_pole_longitude), Deg2Rad(eval_longitude), 0, m[:len(self.ionospheric_pole_latitude)]*BthetaDF(ionospheric_eval_theta_prime, eval_radius, self.R))
                #Radial component of the telluric contribution to the magnetic field at the evaluation point(s)
                Br2= m[len(self.ionospheric_pole_latitude):]*BrDF(telluric_eval_theta_prime, eval_radius, self.Rc)
                #North and east component of the telluric contribution to the magnetic field at the evaluation point(s)            
                Bn2, Be2= Local2Global(Deg2Rad(90- self.telluric_pole_latitude), Deg2Rad(90- eval_latitude), telluric_eval_theta_prime, Deg2Rad(self.telluric_pole_longitude), 
                   Deg2Rad(eval_longitude), 0, m[len(self.ionospheric_pole_latitude):]*BthetaDF(telluric_eval_theta_prime, eval_radius, self.Rc))
                return (np.sum(Br, 1), np.sum(Be, 1), np.sum(-Bn, 1)), (np.sum(Br2, 1), np.sum(Be2, 1), np.sum(-Bn2, 1))
        elif self.mode=='image':
            #Radial component of the ionospheric contribution to the magnetic field at the evaluation point(s)
            Br= m*BrDF(eval_theta_prime, eval_radius, self.R)
            #North and east component of the ionospheric contribution to the magnetic field at the evaluation point(s)
            Bn, Be= Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
               Deg2Rad(eval_longitude), 0, m*BthetaDF(eval_theta_prime, eval_radius, self.R))
            #Radial component of the telluric contribution to the magnetic field at the evaluation point(s)
            Br2= (-self.R/self.Rc)*m*BrDF(eval_theta_prime, eval_radius, (self.Rc**2)/self.R)
            #North and east component of the telluric contribution to the magnetic field at the evaluation point(s)
            Bn2, Be2= Local2Global(Deg2Rad(90- self.telluric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.telluric_pole_longitude), 
               Deg2Rad(eval_longitude), 0, (-self.R/self.Rc)*m*BthetaDF(eval_theta_prime, eval_radius, (self.Rc**2)/self.R))
            return (np.sum(Br, 1), np.sum(Be, 1), np.sum(-Bn, 1)), (np.sum(Br2, 1), np.sum(Be2, 1), np.sum(-Bn2, 1))
        else:
            raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
    def Currents(self, Btheta, Bphi, Br, eval_radius=None, eval_longitude=None, eval_latitude=None, singularity_limit=0):
        """
        Parameters
        ----------
        Btheta : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the theta unit vector from the magnetometer(s).
        Bphi : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the phi unit vector from the magnetometer(s).
        Br : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the radial unit vector from the magnetometer(s).
        eval_radius : float (Meters), optional
            Radius from the centre to the Earth of the evaluation point(s). The default is RE.
        eval_longitude : numpy.ndarray (Degrees), optional
            Longitude position of the evaluation point(s). The default is None.
        eval_latitude : numpy.ndarray (Degrees), optional
            Latitude position of the evaluation point(s). The default is None.
        singularity_limit : float (Meters), optional
            Singularity limit to be used when the evaluation point(s) are too close to the SEC pole(s). The method is based on equation 2.44 from Vanhamäki 2020
            The default is 0 which means ordinary SECS will be used if the singularity limit is not user defined.
            
        Raises
        ------
        ArgumentError
            Will be raised if argument inputs are the wrong type or incorrect.
        ModelError
            Will be rasied if there is no "fitting_matrix" that is an attribute of the class

        Returns
        -------
        Ce, Cn : (numpy.ndarray, numpy.ndarray) (Amps)
            Eastward and northward current components at the evaluation point(s).

        """
        singularity_limit= singularity_limit/self.R
        try:
            #Calculate the amplitudes of the SEC pole(s)
            m= np.dot(self.fitting_matrix, np.concatenate((Btheta, Bphi, Br)).T)
        except AttributeError:
            raise ModelError('Fitting matrix could not be found ensure you have successfully run the "Fitting_Matrix" function or assigned a "fitting_matrix" as an attribute to this class')
        if eval_radius is None:
            eval_radius= self.R
        if eval_longitude is None and eval_latitude is None:
            eval_longitude= self.eval_longitude
            eval_latitude= self.eval_latitude
            try:
                eval_theta_prime= self.eval_theta_prime
            except AttributeError:
                eval_theta_prime=self.eval_ionospheric_theta_prime
        elif eval_longitude is None or eval_latitude is None:
            raise ArgumentError('Evaluation latitude and longitude must be defined if not using the values assigned initially')
        else:
            try:
                #Calculate theta prime of the evaluation points
                eval_theta_prime= theta(self.pole_latitude, self.pole_longitude, eval_latitude, eval_longitude)
            except AttributeError:
                eval_theta_prime= theta(self.ionospheric_pole_latitude, self.ionospheric_pole_longitude, eval_latitude, eval_longitude)
        try:
            #North and east component of the current at the evaluation point(s)
            Cn, Ce=Local2Global(Deg2Rad(90- self.pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.pole_longitude), 
                                Deg2Rad(eval_longitude), m*Vdf(eval_theta_prime, eval_radius, theta0= singularity_limit), 0)
            return np.sum(Ce, 1), np.sum(-Cn, 1)
        except AttributeError:
            if self.mode=='default':
                #North and east component of the current at the evaluation point(s)
                Cn, Ce=Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                                    Deg2Rad(eval_longitude), m[:len(self.ionospheric_pole_latitude)]*Vdf(eval_theta_prime, eval_radius, theta0= singularity_limit), 0)
                return np.sum(Ce, 1), np.sum(-Cn, 1)
            elif self.mode=='image':
                #North and east component of the current at the evaluation point(s)
                Cn, Ce=Local2Global(Deg2Rad(90- self.ionospheric_pole_latitude), Deg2Rad(90- eval_latitude), eval_theta_prime, Deg2Rad(self.ionospheric_pole_longitude), 
                                    Deg2Rad(eval_longitude), m*Vdf(eval_theta_prime, eval_radius, theta0= singularity_limit), 0)
                return np.sum(Ce, 1), np.sum(-Cn, 1)
        else:
            raise ArgumentError("Invalid mode selected use either 'default' or 'image'")
    def Amplitude(self, Btheta, Bphi, Br):
        """
        Parameters
        ----------
        Btheta : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the theta unit vector from the magnetometer(s).
        Bphi : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the phi unit vector from the magnetometer(s).
        Br : numpy.ndarray (Tesla)
            Component of the magnetic field measurement(s) in the radial unit vector from the magnetometer(s).

        Raises
        ------
        ArgumentError
            Will be raised if argument inputs are the wrong type or incorrect.

        Returns
        -------
        Amplitude : numpy.ndarray
            Amplitude(s) of the SEC pole(s).
        """
        if self.mode=='image':
            return np.dot(self.fitting_matrix, np.concatenate((Btheta, Bphi, Br)).T)
        elif self.mode=='default':
            if self.Rc is None:
                return np.dot(self.fitting_matrix, np.concatenate((Btheta, Bphi, Br)).T)
            else:
                return np.dot(self.fitting_matrix, np.concatenate((Btheta, Bphi, Br)).T)[:len(self.ionospheric_pole_latitude)], np.dot(self.fitting_matrix, np.concatenate((Btheta, Bphi, Br)).T)[len(self.ionospheric_pole_latitude):]
        else:
            raise ArgumentError("Invalid mode selected use either 'default' or 'image'")

"""Working Example"""
if __name__== '__main__':
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import matplotlib as mpl
    #Retrieve the SEC pole grid and the evaluation grid
    node_lons, node_lats, eval_lon, eval_lat= Example.ExampleGrid()
    #Retrieve the magnetometer locations and the magnetometer data
    Btheta, Bphi, Br, MagLon, MagLat= Example.ExampleData()
    #Create a SEC pole object
    poles= SECS(node_lons, node_lats, eval_lon, eval_lat)
    #Calculate a fitting matrix
    poles.Fitting_Matrix(MagLon, MagLat)
    #Find the magnetic field at the points on the evaluation grid at RE
    Mag= poles.Magnetic_Field(Btheta, Bphi, Br)
    #Find the current at the points on the evaluation grid at the SEC pole altitude
    Cur= poles.Currents(Btheta, Bphi, Br)
    #Create figure and subplots
    fig = plt.figure(num= 'SEC Example -Grid And Data Taken From Vanhamäki (2020)', constrained_layout=False)
    gs = fig.add_gridspec(1, 4, width_ratios= [2,2,2,1])
    lon_centre, lat_centre= (min(node_lons) + max(node_lons))/2, (min(node_lats) + max(node_lats))/2
    proj=ccrs.LambertConformal(central_longitude=lon_centre, central_latitude=lat_centre)
    axes = [fig.add_subplot(gs[0], projection= proj), fig.add_subplot(gs[1], projection= proj), 
                     fig.add_subplot(gs[2], projection= proj), fig.add_subplot(gs[3])]
    #Add cartopy features to the subplots
    for ax in axes[:-1]:
        ax.scatter(node_lons, node_lats, transform =ccrs.PlateCarree(), color='red', label= 'SEC Poles', s= 0.5, zorder=50)
        PlottingTools.features(ax)
    #Plot magnetometer locations
    for ax in axes[1:-1]:
        ax.scatter(MagLon, MagLat, transform =ccrs.PlateCarree(), marker= '*', color= 'darkorange', label= 'Magnetometer \n Station', zorder=100)
    #Plot magnetometer data
    axes[-2].set_title('Magnetic Field on Ground') 
    Je_m, Jn_m = PlottingTools.Geocentric_to_PlateCarree_vector_components(Bphi, -Btheta, MagLat) #Conversion of vectors to platecarree for plotting
    MagQ=  axes[-2].quiver(MagLon, MagLat,  Je_m, Jn_m,zorder=100, transform = ccrs.PlateCarree(), scale = (1e-5)/8, color='red', label='Measurements')
    Magqk= axes[-2].quiverkey(MagQ, 0.69, 0.71, 200e-9, r'$200 \ nT$', labelpos='E',
                       coordinates='figure')
    #Plot the magnetic field at the points on the evaluation grid
    axes[-2].quiver(eval_lon, eval_lat, *PlottingTools.Geocentric_to_PlateCarree_vector_components(*Mag[1:],eval_lat), zorder=100,scale=(1e-5)/8, transform = ccrs.PlateCarree(), color='Black', label='Total Field', alpha=0.7)
    #Plot the radial component of the magnetic field at the points on the evaluation grid using a matplotlib colour map
    axes[-2].scatter(eval_lon, eval_lat, c= Mag[0], vmin=min(Mag[0]), vmax=max(Mag[0]), cmap=mpl.cm.bwr, transform= ccrs.PlateCarree(), zorder=50)
    #Plot the current at the points on the evaluation grid
    axes[1].set_title('Ionospheric Currents')
    CurQ=  axes[1].quiver(eval_lon, eval_lat, *PlottingTools.Geocentric_to_PlateCarree_vector_components(*Cur,eval_lat), zorder=100,scale=(1e1)/2, transform = ccrs.PlateCarree(), color='Black', label='Ionospheric Current', alpha=0.7)
    Curqk= axes[-1].quiverkey(CurQ, 0.45, 0.71, 500e-3, r'$500 \ Akm^{-1}$', labelpos='E',
                       coordinates='figure')
    #Plot the magnitudes of currents using a matplotlib colour map
    CurMag=np.sqrt((Cur[0]**2) +Cur[1]**2)
    axes[1].scatter(eval_lon, eval_lat, c= CurMag, vmin=min(CurMag), vmax= 0.6*max(CurMag), cmap= mpl.cm.afmhot_r, transform= ccrs.PlateCarree(), zorder=50)
    #Plot the amplitudes of the SEC poles using a matplotlib colour map
    axes[0].set_title('SEC POLE Amplitudes')
    scalars= poles.Amplitude(Btheta, Bphi, Br)
    axes[0].scatter(node_lons, node_lats, c= scalars, vmin=-45000, vmax=45000, cmap=mpl.cm.jet, transform= ccrs.PlateCarree(), zorder=50)
    mappable=mpl.cm.ScalarMappable(norm = mpl.colors.Normalize(vmin=-45000, vmax=45000), cmap=mpl.cm.jet)
    cbar=plt.colorbar(mappable=mappable)
    cbar.set_label('Ionospheric SEC Pole Amplitudes', rotation=270, labelpad=30)
    #Create a plot comparing the measurements to the model for validity
    model= poles.Magnetic_Field(Btheta, Bphi, Br, eval_longitude= MagLon, eval_latitude= MagLat)
    model= model[0], *PlottingTools.Geocentric_to_PlateCarree_vector_components(*model[1:], MagLat) #Conversion of vectors to platecarree
    axes[-1].plot(np.linspace(-1.1*np.max(np.abs(model))*1e9, 1.1*np.max(np.abs(model))*1e9, 100), np.linspace(-1.1*np.max(np.abs(model))*1e9, 1.1*np.max(np.abs(model))*1e9, 100))
    axes[-1].set_aspect('auto')
    axes[-1].scatter(Br*1e9, model[0]*1e9, label='Br component')
    axes[-1].scatter(Je_m*1e9, model[1]*1e9, label='Be component')
    axes[-1].scatter(Jn_m*1e9, model[-1]*1e9, label='Bn component')
    axes[-1].set_title('Validity')
    axes[-1].set_xlabel('Data (nT)')
    axes[-1].set_ylabel('Model (nT)', y=0.15)
    legend_x = -0.6
    legend_y = 1.014
    handle= []
    labels=[]
    for ax in axes[1:]:
        tmp_handles, tmp_labels =ax.get_legend_handles_labels()
        handle+= tmp_handles
        labels+= tmp_labels
    labels, index= np.unique(labels, return_index=True)
    handles=[]
    for i in index:
        handles.append(handle[i])
    leg=axes[0].legend(handles, labels, loc='upper left', bbox_to_anchor=(legend_x, legend_y))
    try:
        #Add the IMAGE logo
        im = plt.imread('https://space.fmi.fi/image/www/images/IMAGE_logo_2.png')
        newax = fig.add_axes([0, 0.3, 0.1, 0.1], anchor='NE', zorder=-1)
        newax.imshow(im)
        newax.axis('off')
    except:
        pass