# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:29:57 2023

@author: dscho
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

plt.style.use("fivethirtyeight")

#spacestep in mm
#timestep in seconds




def run_heat(dt = 0.02, dx = 0.2, csquare=1, xmax = 1.0, tmax=0.2):

    '''
    
    
    Returns
    ----
    
    x: numpy vector
        - array of position locations
        
    t: numpy vector
        - array of time points
        
    temp: numpy 2D array
        - temperature as a function of time & space
        
    xmax, tmax : float, default to 1 and 0.2
    
    '''
    
    #Set constant r
    r = csquare * dt/dx**2
    
    #Create space and time grids
    x = nparange(0, xmax+dx, dx)
    t = nparange(0, tmax+dt, dt)
    #Save number of points
    M, N = x.size, t.size
    
    #Temp solution array
    temp = np.zeros([M,N])

    temp[0, :] = 0
    temp[-1, :] = 0
    temp[:, 0] = 4*x - 4*(x**2)
    
    
    #Solution to equation
    for j in range(0, N-1):
        for i in range(1, M-1):
            temp[i, j+1] = (1-(2*r))*temp[i, j] + \
                r*(temp[i+1, j] + temp[i-1, j])
                
    return x, t, temp
    
    
    

#map = axes.pcolor(time, x, heat, cmap='seismic', vmin=-25, vmax=25)
#plt.colorbar(map, ax=axes, label='Temperature ($C$)')
