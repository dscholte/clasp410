# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:12:48 2023

@author: dscho





"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
from scipy.integrate import odeint
from scipy.integrate import solve_ivp, ode

y0 = [10,1]

t = np.linspace(0,50,num=1000)

alpha = 1.1
beta = 0.4
delta = 0.1
gamma = 0.4

params = [alpha,beta,delta,gamma]


def simulation(t, variables, alpha,beta,delta,gamma):
    
    # prey population
    x = variables[0]
    
    # predator population
    y = variables[1]
    
    #derivatives
    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - gamma * y
    
    return ([dxdt,dydt])

dT=0.1

#Solve ODE
y = solve_ivp(simulation, [0,100], [10,1], args = params, method='DOP853',
              max_step = dT)


f, (ax1,ax2) = plt.subplots(2)

time, N1, N2 = y.t, y.y[0, :], y.y[1, :]

line1, = ax1.plot(time,N1, color = "b")

line2, = ax2.plot(time,N2, color = "r")

ax1.set_ylabel("Fish")
ax2.set_ylabel("Bears")
ax2.set_xlabel("Time")
    
    
plt.show()
    
    
    
    
    
    
    
    
    
    
    
    