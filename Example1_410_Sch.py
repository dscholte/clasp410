# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 11:12:39 2023
@author: dscho

Coffee Example
"""

import numpy as np
import matplotlib.pyplot as plt

# Create a time array:
tfinal, tstep = 600, 1
time = np.arange(0, tfinal, tstep)

def tempsolve(time, k=1./300, T_env=25, T_init=90):
    '''
    This function takes an array of times and returns an array of 
    temperatures corresponding to each time.
    -------
    Parameters
    ==========
    time = numpy array of times
        Array of time inputs for which you want corresponding temps
    

    '''
    
    temp  = T_env + ((T_init-T_env)* np.exp(-k * time))
    
    return temp

def time_to_temp(T_targ, k = 1/300, T_env=20, T_init=90):
    
    '''
    '''
    return (-1/k) * np.log((T_targ - T_env)/(T_init - T_env))

T_cream = tempsolve(time, T_init=85)
T_nocream = tempsolve(time, T_init=90)

# Time until drinkable
t_cream = time_to_temp(60, T_init = 85) #Add cream right away
t_nocream = time_to_temp(60, T_init = 90) #Add cream at 60
t_smart = time_to_temp(65, T_init = 90) #Add cream once at 65
    
#Create figure and axis
fig, ax = plt.subplots(1,1)

ax.plot(time, T_nocream, label="No cream until 60")
ax.plot(time, T_cream, label="Cream immediately")

    
ax.axvline(t_nocream, c='red', ls='--', label='No Cream: T=60')
ax.axvline(t_cream, c='blue', ls='--', label='Added cream at 60')
ax.axvline(t_smart, c='green', ls='--', label='No Cream: T=60')

ax.legend(loc='best')


    
    
    
    
    
    
    