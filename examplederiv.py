# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 09:53:24 2023

@author: dscho

This file explores different numerical approximations for first derivatives.

"""
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('fivethirtyeight')

def derivative(fx, dx=1):
    '''
    take the forward difference first derivative of 'fx'
    Final element uses a backwards difference that is the same as second to 
    last element.
    '''
    
    
    #create result array of zeros
    dfdx = np.zeros(fx.size)
    
    dfdx[:-1] = (fx[1:] - fx[:-1]) / dx
    
    dfdx[-1] = (fx[-1] - fx[-2]) / dx
    
    return dfdx

dx = 0.5
x = np.arange(0,4*np.pi,0.5)
fx = np.sin(x)
dfdx_sol = np.cos(x)
dfdx = derivative(fx, dx = dx)

#Plot

plt.plot(x, dfdx_sol, label = 'Analytical Solution')
plt.plot(x, dfdx, label = 'Numerical approxi')
plt.legend()


dx_all = np.linspace(0.001,1,100)
err = []

for dx in dx_all:
    x = np.arange(0,4*np.pi,0.5)
    fx = np.sin(x)
    dfdx_sol = np.cos(x)
    dfdx = derivative(fx, dx = dx)
    
    
    # calculate error:
    err.append(np.max(np.abs(dfdx_sol - dfdx)))

fig, ax = plt.subplots(1,1)
ax.plot(dx_all, err)
    
    
    
    
    
    
    
    





