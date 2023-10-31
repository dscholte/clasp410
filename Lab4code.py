# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:29:57 2023

@author: dscho
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

plt.style.use("fivethirtyeight")


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
#x1, t1, temp1 = run_neumann(0.0002, 0.02, .025, 1.0, 2)
#fig, axes = plt.subplots(1, 1)

#map = axes.pcolor(t1, x1, temp1, cmap='inferno', vmin=0, vmax=1)


#plt.colorbar(map, ax=axes, label='Temperature ($C$)')
#plt.title('Neumann Boundary Condition')
#plt.show()
    
    
#-----------------------------------------------------------------    
#Question 2

t_kanger = np.array([-19.7, -21.0, -17., -8.4, 2.3, 8.4,
10.7, 8.5, 3.1, -6.0, -12.0, -16.9])
 
# initializing K
ques3addahalf = 0.5
ques3addone = 1
ques3addthree = 3
 
# using list comprehension
# adding K to each element
addhalf = [x + ques3addahalf for x in t_kanger]
addone = [x + ques3addone for x in t_kanger]
addthree = [x + ques3addthree for x in t_kanger]


def temp_kanger(t,temps1):
    '''
    For an array of times in days, return timeseries of temperature for
    Kangerlussuaq, Greenland.
    '''
    t_amp = (temps1 - temps1.mean()).max()
    return t_amp*np.sin(np.pi/180 * t - np.pi/2) + temps1.mean()
    

def greenland(dt, dx, csquare, xmax, tmax, climate=temp_kanger(t,t_kanger)):
    
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
    temp[-1, :] = temp_kanger(t,t_kanger)
    temp[1:-1,0] = 0
    
    
    #Solution to equation
    for j in range(0, N-1):
        for i in range(1, M-1):
            temp[i, j+1] = (1-(2*r))*temp[i, j] + \
                r*(temp[i+1, j] + temp[i-1, j])
                
    return x, t, temp


def plot_temp(x, time, temp, axes, xlabel='Time ($s$)', title='',
              ylabel='Distance ($m$)', clabel=r'Temperature ($^{\circ} C$)',
              cmap='inferno', inverty=False, **kwargs):
    '''
    Add a pcolor plot of the heat equation to `axes`. Add a color bar.

    Parameters
    ----------
    x
        Array of position values
    time
        Array of time values 
    temp
        array of temperature solution
    xlabel, ylabel, title, clabel
        Axes labels and titles
    cmap : inferno
        Matplotlib colormap name.
    '''

    map = axes.pcolor(time, x, temp, cmap=cmap, **kwargs)
    plt.colorbar(map, ax=axes, label=clabel)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title)
    

def profile():
    
    dt = 10
    dx = 1.0
    years = 50
    
    x, time, temp = greenland(dt, dx, 0.25, 100, years*365, temp_kanger(t,t_kanger))
    
    maxtemp = np.abs(temp).max()
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    plot_temp(x, time/365., temp, ax, xlabel='Time (Years)', ylabel='Depth ($m$)',
              cmap='seismic', vmin=-maxtemp, vmax=maxtemp,
              clabel='Temperature', 
              title='Question 2')
    fig.tight_layout()
    
    # Create fianl temp profile
    fig, ax2 = plt.subplots(1, 1, figsize=(10, 8))
    ax2.plot(temp[:, int(-365/dt):].min(axis=1), x, color='blue', label='Winter')
    ax2.plot(temp[:, int(-365/dt):].max(axis=1), x, '--', color='red', label='Summer')

    
    ax2.legend(loc='best')
    ax2.set_ylabel('Depth (m)')
    ax2.set_xlabel('Temperature (degC)')
    ax2.set_title('Ground Temperature')
    ax2.set_xlim([-7, 8])
    ax2.set_ylim([-5, 105])
    fig.tight_layout()
    
    
profile()
    
    #Question 3
    
def three():
    
    dt = 10
    dx = 1.0
    nyear = 50
    xhalf, timehalf, temphalf = greenland(dt, dx, 0.25, 100, nyear*365, temp_kanger(t,addhalf))    
    xone, timeone, tempone = greenland(dt, dx, 0.25, 100, nyear*365, temp_kanger(t,addone))    
    xthree, timethree, tempthree = greenland(dt, dx, 0.25, 100, nyear*365, temp_kanger(t,addthree))    

    fig, ax3 = plt.subplots(1, 1, figsize=(10, 8))
    ax3.plot(temphalf[:, int(-365/dt):].min(axis=1), xhalf, color='blue', label='Winter')
    ax3.plot(temphalf[:, int(-365/dt):].max(axis=1), xhalf, '--', color='red', label='Summer')
    fig
    
three()
    
