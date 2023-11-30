# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 10:55:08 2023

@author: dscho
"""

import numpy as np
import matplotlib.pyplot as plt


def temp_warm(lats_in):
    '''
    Create a temperature profile for modern day "warm" earth.

    Parameters
    ----------
    lats_in : Numpy array
        Array of latitudes in degrees where temperature is required.
        0 corresponds to the south pole, 180 to the north.

    Returns
    -------
    temp : Numpy array
        Temperature in Celcius.
    '''

    # Set initial temperature curve
    T_warm = np.array([-47, -19, -11, 1, 9, 14, 19, 23, 25, 25,
                       23, 19, 14, 9, 1, -11, -19, -47])

    # Get base grid:
    npoints = T_warm.size
    dlat = 180 / npoints  # Latitude spacing.
    lats = np.linspace(dlat/2., 180-dlat/2., npoints)  # Lat cell centers.

    # Fit a parabola to the above values
    coeffs = np.polyfit(lats, T_warm, 2)

    # Now, return fitting sampled at "lats".
    temp = coeffs[2] + coeffs[1]*lats_in + coeffs[0] * lats_in**2

    return temp


def insolation(S0, lats):
    '''
    Given a solar constant (`S0`), calculate average annual, longitude-averaged
    insolation values as a function of latitude.
    Insolation is returned at position `lats` in units of W/m^2.

    Parameters
    ----------
    S0 : float
        Solar constant (1370 for typical Earth conditions.)
    lats : Numpy array
        Latitudes to output insolation. Following the grid standards set in
        the diffusion program, polar angle is defined from the south pole.
        In other words, 0 is the south pole, 180 the north.

    Returns
    -------
    insolation : numpy array
        Insolation returned over the input latitudes.
    '''

    # Constants:
    max_tilt = 23.5   # tilt of earth in degrees

    # Create an array to hold insolation:
    insolation = np.zeros(lats.size)

    #  Daily rotation of earth reduces solar constant by distributing the sun
    #  energy all along a zonal band
    dlong = 0.01  # Use 1/100 of a degree in summing over latitudes
    angle = np.cos(np.pi/180. * np.arange(0, 360, dlong))
    angle[angle < 0] = 0
    total_solar = S0 * angle.sum()
    S0_avg = total_solar / (360/dlong)

    # Accumulate normalized insolation through a year.
    # Start with the spin axis tilt for every day in 1 year:
    tilt = [max_tilt * np.cos(2.0*np.pi*day/365) for day in range(365)]

    # Apply to each latitude zone:
    for i, lat in enumerate(lats):
        # Get solar zenith; do not let it go past 180. Convert to latitude.
        zen = lat - 90. + tilt
        zen[zen > 90] = 90
        # Use zenith angle to calculate insolation as function of latitude.
        insolation[i] = S0_avg * np.sum(np.cos(np.pi/180. * zen)) / 365.

    # Average over entire year; multiply by S0 amplitude:
    insolation = S0_avg * insolation / 365

    return insolation


def gen_grid(npoints=18):
    '''
    
    Creating the grid for pole-to-pole latitudes
    
    Parameters
    ----------
    npoints : TYPE, optional
        DESCRIPTION. number of points on the grid

    Returns
    -------
    dlat : float
        latitude spacing in degrees
    lats : numpy array
        latitudes in degrees
    edge : numpy array
        latitude bin edges in degrees

    '''
    dlat = 180/ npoints
    lats =  np.linspace(dlat/2, 180-dlat/2, npoints)
    edge = np.linspace(0, 180, npoints+1)
    
    return dlat, lats, edge

def snowearth(npoints=18, dt = 1, tstop = 10000, radearth=6357000, lambda_heat=100 , dlat = 10, S0=1370):
    '''
    aye docstring
    
    Parameters
    --------
    npoints: integer, 18
        number of latitude points
    
    '''
    dlat, lats, edges = gen_grid(npoints)
    
    delta_t = dt
    delta_y = 2*np.pi*radearth * (dlat/360)
    
    #change units lambda
    lambda_heat_year = lambda_heat * 60 * 60 * 24 * 365
    
    identity_mat = np.identity(npoints)
    
    #create values for original temperature line
    T_warm = temp_warm(lats)
    T_warm_init = T_warm
    
    #initialize A grid
    A_mat = -2*np.identity(len(lats))
    
    # creating tridiagonal matrix
    for i in range(1, len(A_mat)-1):
        A_mat[i, i+1]= 1
        A_mat[i, i-1]= 1
        
    A_mat[0 ,1] = 2
    A_mat[-1, -2] = 2
    
    A_mat = (1/(delta_y**2))*A_mat
   
    #L = I - deltat*lambda*A
    L_mat = identity_mat - (lambda_heat_year*delta_t*A_mat)
    
    #invert L and matrix multiply to create values for diffusion
    for i in range(0,tstop):
        L_inv = np.linalg.inv(L_mat)
        T_warm = np.matmul(L_inv,T_warm)
        
    
    #L = I - deltat*lambda*A
    #f = deltat*lambda*B*T_warm_ * 
    
    B = np.zeros((npoints,npoints))
    B[np.arange(npoints-1), np.arange(npoints-1)+1] = 1
    B[np.arange(npoints-1)+1, np.arange(npoints-1)] = -1
    B[0, :] = B[-1, :] = 0
    
    #set the surface area of the side of each latitude ring at bin center
    Axz = np.pi * ((radearth+50.0)**2 - radearth**2) * np.sin(np.pi/180.*lats)
    
    #dAxz/dlat, doesn't change 
    dAxz = np.matmul(B, Axz) / (Axz * 4 * delta_y**2)
    
    #get total number of steps
    nsteps = int(tstop/delta_t)
    
    insol = insolation(S0, lats)
    
    for i in range(nsteps):
        
        spherecorr = lambda_heat_year*delta_t*np.matmul(B, T_warm)*dAxz
        
        
        T_warm += spherecorr
        
        
        radiative = (1-0.3)
        
        T_warm += delta_t*insol
        
        T_warm = np.matmul(L_inv, T_warm)
        
        
    
    
    plt.plot(lats,T_warm)
    plt.plot(lats,T_warm_init)

    return lats, T_warm








    
snowearth()