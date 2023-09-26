# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:12:48 2023

@author: dscho





"""

import numpy as np
import matplotlib.pyplot as plt

plt.style.use("fivethirtyeight")
from scipy.integrate import odeint
from scipy.integrate import solve_ivp, ode
from dataclasses import dataclass


@dataclass
class modelrundata:
    times: list[float]
    P_sol: list[float]
    label: str

    def __post_init__(self):
        if len(self.P_sol) != len(self.times):
            raise ValueError("Lists are not the same length! Doesn't make sense!")


y0 = [10, 1]

t = np.linspace(0, 50, num=1000)

alpha = 1.1
beta = 0.4
delta = 0.1
gamma = 0.4

params = [alpha, beta, delta, gamma]


def simulation(t, variables, alpha, beta, delta, gamma):
    # prey population
    x = variables[0]

    # predator population
    y = variables[1]

    # derivatives
    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - gamma * y

    return (dxdt, dydt)


def taylor_step(prey, pred, dt, a, b, c, d) -> float:
    """
    Take First Order Taylor Step  Forward given a dt and initial vals

    """

    prey_deriv, pred_deriv = dt * (simulation(dt, [prey, pred], a, b, c, d))
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

    ## Run Euler Method ##
    for index, value in enumerate(times[0:-1]):
        solution_prey[index + 1], solution_pred[index + 1] = taylor_step(
            solution_prey[index], solution_pred[index], dt, a, b, c, d
        )

    ## Return Results ##
    return times,solution_prey,solution_pred


dT = 0.1
time, N1, N2 = euler_method(100,5,1.1,0.4,0.1,0.4)
f, (ax1, ax2) = plt.subplots(2)

ax1.plot(time, N1, color="b")

ax2.plot(time, N2, color="r")

ax1.set_ylabel("Fish")
ax2.set_ylabel("Bears")
ax2.set_xlabel("Time")
# Solve ODE
y = solve_ivp(simulation, [0, 100], [10, 1], args=params, method="DOP853", max_step=dT)


f, (ax1, ax2) = plt.subplots(2)

time, N1, N2 = y.t, y.y[0, :], y.y[1, :]

(line1,) = ax1.plot(time, N1, color="b")

(line2,) = ax2.plot(time, N2, color="r")

ax1.set_ylabel("Fish")
ax2.set_ylabel("Bears")
ax2.set_xlabel("Time")


plt.show()
