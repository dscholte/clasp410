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


def run_heat(dt, dx, csquare, xmax, tmax):

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
    
    if (dt > ((dx**2)/(2*csquare))):
        raise ValueError('Stability Criterion not met' + 
                         f'dt={dt:6.2f}; dx={dx:6.2f}; csquare={csquare}')
    
    #Set constant r
    r = csquare * dt/dx**2
    
    #Create space and time grids
    x = np.arange(0, xmax+dx, dx)
    t = np.arange(0, tmax+dt, dt)
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


def run_neumann(dt, dx, csquare, xmax, tmax):

    
    if (dt > ((dx**2)/(2*csquare))):
        raise ValueError('Stability Criterion not met' + 
                         f'dt={dt:6.2f}; dx={dx:6.2f}; csquare={csquare}')
    
    #Set constant r
    r = csquare * dt/dx**2
    
    #Create space and time grids
    x = np.arange(0, xmax+dx, dx)
    t = np.arange(0, tmax+dt, dt)
    #Save number of points
    M, N = x.size, t.size
    
    #Temp solution array
    temp = np.zeros([M,N])

    temp[0, :] = 0
    temp[-1, :] = 0
    temp[:, 0] = 4*x - 4*(x**2)
    
    
    #Solution to equation
    for j in range(0, N-1):
        temp[0,j]=temp[1,j]
        temp[-1,j]=temp[-2,j]
        for i in range(1, M-1):
            temp[i, j+1] = (1-(2*r))*temp[i, j] + \
                r*(temp[i+1, j] + temp[i-1, j])
                
    return x, t, temp

#---------------------------------------------------------------
#Validation Code Question 1
    
#Dirichlet
x, t, temp = run_heat(0.02, 0.2, 1, 1, .2)
fig, axes = plt.subplots(1, 1)

map = axes.pcolor(t, x, temp, cmap='inferno', vmin=0, vmax=1)

plt.colorbar(map, ax=axes, label='Temperature ($C$)')
plt.title('Dirichlet Boundary Condition')
plt.show()

#Neumann
x1, t1, temp1 = run_neumann(0.0002, 0.02, .025, 1.0, 2)
fig, axes = plt.subplots(1, 1)

map = axes.pcolor(t1, x1, temp1, cmap='inferno', vmin=0, vmax=1)


plt.colorbar(map, ax=axes, label='Temperature ($C$)')
plt.title('Neumann Boundary Condition')
plt.show()
    
    
#-----------------------------------------------------------------    
#Question 2
#csquare mm^2/s to #m^2/day


t_kanger = np.array([-19.7, -21.0, -17., -8.4, 2.3, 8.4,
10.7, 8.5, 3.1, -6.0, -12.0, -16.9])
    
def temp_kanger(t):
    '''
    For an array of times in days, return timeseries of temperature for
    Kangerlussuaq, Greenland.
    '''
    t_amp = (t_kanger - t_kanger.mean()).max()
    return t_amp*np.sin(np.pi/180 * t - np.pi/2) + t_kanger.mean()
    

def greenland(dt, dx, csquare, xmax, tmax):
    
    landc2 = csquare*(1/1000000)*24*60*60
    
    if (dt > ((dx**2)/(2*landc2))):
        raise ValueError('Stability Criterion not met' + 
                         f'dt={dt:6.2f}; dx={dx:6.2f}; csquare={csquare}')
    
    
    #Set constant r
    r = landc2 * dt/dx**2
    
    #Create space and time grids
    x = np.arange(0, xmax+dx, dx)
    t = np.arange(0, tmax+dt, dt)
    #Save number of points
    M, N = x.size, t.size
    
    #Temp solution array
    temp = np.zeros([M,N])

    temp[0, :] = 5
    temp[-1, :] = temp_kanger(t)
    temp[1:-1,0] = 0
    
    
    #Solution to equation
    for j in range(0, N-1):
        for i in range(1, M-1):
            temp[i, j+1] = (1-(2*r))*temp[i, j] + \
                r*(temp[i+1, j] + temp[i-1, j])
                
    return x, t, temp


#Question 2 Results
xland, tland, templand = greenland(.5, 1.0, 0.25, 100, 1000)
fig, axes = plt.subplots(1, 1)

map = axes.pcolor(tland, xland, templand, cmap='seismic', vmin=-25, vmax=25)


plt.colorbar(map, ax=axes, label='Temperature ($C$)')
plt.title('Greenland')
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    