#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 13:04:35 2021

@author: simon
"""
import numpy as np
from SECpy import BrDF, BthetaDF, RE, Rad2Deg
import matplotlib.pyplot as plt
u0= 1.23E-6
def new_BrDF(thetap, r, R, singularity_limit):
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
    A= u0/(4*np.pi*r)
    B= np.zeros(thetap.shape)
    if r<R:
        B[thetap>=singularity_limit]= (1+ s**2 -2*s*np.cos(thetap[thetap>=singularity_limit]))**-0.5
        B[thetap<singularity_limit]= (1+ s**2 -2*s*(np.cos(singularity_limit)**2)*(1/np.cos(thetap[thetap<singularity_limit])))**-0.5

        B-=1
    elif r>R:
        B[thetap>=singularity_limit]= s/np.sqrt(1+s**2 -2*s*np.cos(thetap[thetap>=singularity_limit]))
        B[thetap<singularity_limit]=s/np.sqrt(1+s**2 -2*s*(np.cos(singularity_limit)**2)*(1/np.cos(thetap[thetap<singularity_limit])))
        B-= s
    else:
        raise ValueError("Invalid radius inputs shouldn't be equal")

    return A*B
def new_BthetaDF(thetap, r, R, singularity_limit):
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
    A= -u0/(4*np.pi*r)
    B= np.zeros(thetap.shape)
    if r<R:
        B[thetap>=singularity_limit]= ((s-np.cos(thetap[thetap>=singularity_limit]))*(1+ s**2 -2*s*np.cos(thetap[thetap>=singularity_limit]))**-0.5)
        B[thetap>=singularity_limit]+= np.cos(thetap[thetap>=singularity_limit])
        B[thetap>=singularity_limit]/=np.sin(thetap[thetap>=singularity_limit])
        B[thetap<singularity_limit]= ((s- (np.cos(singularity_limit)**2)/np.cos(thetap[thetap<singularity_limit]))*(1+ s**2 -2*s*(np.cos(singularity_limit)**2)/np.cos(thetap[thetap<singularity_limit]))**-0.5)
        B[thetap<singularity_limit]+= (np.cos(singularity_limit)**2)/np.cos(thetap[thetap<singularity_limit])
        B[thetap<singularity_limit]/=(np.sin(singularity_limit)**2)/np.sin(thetap[thetap<singularity_limit])
    elif r>R:
        B[thetap>=singularity_limit]=((1-s*np.cos(thetap[thetap>=singularity_limit]))/np.sqrt(1+s**2 -2*s*np.cos(thetap[thetap>=singularity_limit])))
        B[thetap>=singularity_limit]/=np.sin(thetap[thetap>=singularity_limit])
        B[thetap<singularity_limit]=((1-s*(np.cos(singularity_limit)**2)/np.cos(thetap[thetap<singularity_limit]))/np.sqrt(1+s**2 -2*s*(np.cos(singularity_limit)**2)/np.cos(thetap[thetap<singularity_limit])))
        B[thetap<singularity_limit]/= (np.sin(singularity_limit)**2)/np.sin(thetap[thetap<singularity_limit])
#(np.sin(singularity_limit)**2)/np.sin(thetap[thetap<singularity_limit])
        B-=1
    else:
        raise ValueError('Invalid radius inputs cannot be equal')
    return A*B
rm= RE
RI= RE+110e3
fig= plt.figure()
ax= fig.add_subplot(121)
ax2= fig.add_subplot(122)
thetap= np.linspace(0, np.pi/4)
ax.plot(Rad2Deg(thetap), BrDF(thetap, rm, RI), label='old')
ax.plot(Rad2Deg(thetap), new_BrDF(thetap, rm, RI, np.deg2rad(10)), label='new')
ax.set_title('Br')
ax.set_xlabel("theta'")
ax.set_ylabel('Br (T)/ I0')
ax.legend(loc='best')
ax2.plot(Rad2Deg(thetap), BthetaDF(thetap, rm, RI), label='old')
ax2.plot(Rad2Deg(thetap), new_BthetaDF(thetap, rm, RI, np.deg2rad(10)), label='new')
ax2.set_title('Btheta')
ax2.set_xlabel("theta'")
ax2.set_ylabel('Btheta (T)/ I0')
ax2.legend(loc='best')