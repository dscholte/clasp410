# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 10:20:04 2023

@author: dscho

aye coffee problem eulers
"""


import numpy as np
import matplotlib.pyplot as plt



def hot(T, k, Tenv):
    '''
    Given current temp (T), cooling rate (k) and environmental temp (Tenv), 
    return the derivative of the temperature of the cooling body
    '''
    dt = -k*(T-Tenv)
    return dt



def euler(T_init, deltat=60, tstart = 0, tstop = 600, k = 1/300, 
          Tenv = 21, Tstop = 60):
    '''
    solve the cooling equaiton given the time array, 'time', with times 
    evenly spaced by some deltat in units of seconds.
    
    Also require T_init an initial temperature of the coffee in degrees celsius
    '''
    
    #Create time array
    time  = np.arange(tstart, tstop, deltat)
    
    # cReate solution array
    T_sol = np.zeros(time.size)
    T_sol[0]= T_init
    
    #iterate through solver
    for i in range(time.size-1):
        T_sol[i+1] = T_sol[i]+ deltat * hot(T_sol[i], k=k, Tenv = Tenv)
        
        if T_sol[i+1] <= Tstop:
            T_sol = T_sol[:i+2]
            time = time[:i+2]
            break
        
    return time, T_sol
        






