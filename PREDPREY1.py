# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:12:48 2023

@author: dscho





"""

import numpy as np
import matplotlib.pyplot as plt

plt.style.use("fivethirtyeight")

from scipy.integrate import solve_ivp
from dataclasses import dataclass


@dataclass
class modelrundata:
    times: list[float]
    P_sol: list[float]
    label: str

    def __post_init__(self):
        if len(self.P_sol) != len(self.times):
            raise ValueError("Lists are not the same length! Doesn't make sense!")

t = np.linspace(0, 50, num=1000)

alpha = 1
beta = 2
delta = 1
gamma = 3

params = [alpha, beta, delta, gamma]

N1initial = 0.3
N2initial = 0.6

#Competition equations
def dNdt_comp(t, variables, alpha, beta, delta, gamma):
    # prey population
    x = variables[0]

    # predator population
    y = variables[1]
    
    dn1dt = (alpha * x * (1 - x)) - (beta * x * y)
    dn2dt = (delta * y * (1 - y)) - (gamma * x * y)
    
    return (dn1dt, dn2dt)

#Predator-Prey equations
def simulation(t, variables, alpha, beta, delta, gamma):
    # prey population
    x = variables[0]

    # predator population
    y = variables[1]

    # derivatives
    dxdt = (alpha * x) - (beta * x * y)
    dydt = (-delta * y) + (gamma * x * y)

    return (dxdt, dydt)


def taylor_step(prey, pred, dt, a, b, c, d) -> float:
    """
    Take First Order Taylor Step  Forward given a dt and initial vals

    """

    
    results = (simulation(dt, [prey, pred], a, b, c, d))
    prey_deriv = dt * results[0]
    pred_deriv = dt *  results[1]

    return prey + prey_deriv, pred + pred_deriv



def euler_method(
    N_prey_init: float,
    N_pred_init: float,
    a,
    b,
    c,
    d,
    dt: float = 1,
    t_stop: float = 600,
) -> modelrundata:
    """
    Define Euler Method
    """

    ## Create Initial Arrays ##
    times = np.arange(0, t_stop, dt)
    solution_prey = np.zeros(len(times))
    solution_pred = np.zeros(len(times))

    ## Set Initial Conditions ##
    solution_prey[0] = N_prey_init
    solution_pred[0] = N_pred_init

    ## Run Euler Method with taylor_step function##
    for index, value in enumerate(times[0:-1]):
        solution_prey[index + 1], solution_pred[index + 1] = taylor_step(
            solution_prey[index], solution_pred[index], dt, a, b, c, d
        )

    ## Return Results ##
    return times,solution_prey,solution_pred


dT = 1
t_stop = 100

#Running Euler Method pred-prey
time, N1, N2 = euler_method(N2initial, N1initial, 1, 2, 1, 3,0.05,t_stop)
f, (ax1, ax2) = plt.subplots(2)

ax1.plot(time, N1, color="b")

ax2.plot(time, N2, color="r")

ax1.set_ylabel("Fish")
ax2.set_ylabel("Bears")
ax2.set_xlabel("Time")







'''
Solve the Lotka-Volterra competition and predator/prey equations using
Scipy's ODE class and the adaptive step 8th order solver.
Parameters
----------
simulation : function
A python function that takes `time`, [`N1`, `N2`] as inputs and
returns the time derivative of N1 and N2.

N1_init, N2_init : float
Initial conditions for `N1` and `N2`, ranging from (0,1]
                                                    
dT : float, default=10
Largest timestep allowed in years.

t_final : float, default=100
Integrate until this value is reached, in years.

args=params : float, alpha,beta,gamma,delta
Lotka-Volterra coefficient values

Returns
-------
time : Numpy array
Time elapsed in years.

N1, N2 : Numpy arrays
Normalized population density solutions.
'''

# Solve Competition ODE with ivp
zeta = solve_ivp(dNdt_comp, [0, t_stop], [N1initial,N2initial], args=params, 
              method="DOP853", max_step=dT)

# Solve Predator-Prey ODE with ivp
y = solve_ivp(simulation, [0, t_stop], [N1initial,N2initial], args=params, 
              method="DOP853", max_step=dT)



f, (ax1) = plt.subplots(1,1, figsize = (10,6))

timepp, N1pp, N2pp = y.t, y.y[0, :], y.y[1, :]

(line1,) = ax1.plot(timepp, N1pp, color="b")

(line2,) = ax1.plot(timepp, N2pp, color="r")

ax1.set_title('Predator-Prey RK8')
ax1.set_ylabel("Moose")



plt.show()

v, (ax3) = plt.subplots(1, 1, figsize = (10,6))

timecomp, N1comp, N2comp = zeta.t, zeta.y[0, :], zeta.y[1, :]

(line1,) = ax3.plot(timecomp, N1comp, color="b")

(line2,) = ax3.plot(timecomp, N2comp, color="r")

ax3.set_title('Competition RK8')
ax3.set_ylabel("hey")



plt.show()