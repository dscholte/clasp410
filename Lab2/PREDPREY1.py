# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:12:48 2023

@author: dscho, manishrv





"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from dataclasses import dataclass
import os

plt.style.use("fivethirtyeight")


@dataclass
class modelrundata:
    times: list[float]
    prey: list[float]
    pred: list[float]
    label: str

    def __post_init__(self):
        if len(self.prey) != len(self.times) or len(self.pred) != len(self.times):
            raise ValueError("Lists are not the same length! Doesn't make sense!")


# Competition equations
def dNdt_comp(t, variables, alpha, beta, delta, gamma):
    # prey population
    x = variables[0]

    # predator population
    y = variables[1]

    dn1dt = (alpha * x * (1 - x)) - (beta * x * y)
    dn2dt = (delta * y * (1 - y)) - (gamma * x * y)

    return (dn1dt, dn2dt)


# Predator-Prey equations
def simulation(t, variables, alpha, beta, delta, gamma):
    # prey population
    x = variables[0]

    # predator population
    y = variables[1]

    # derivatives
    dxdt = (alpha * x) - (beta * x * y)
    dydt = (-delta * y) + (gamma * x * y)

    return (dxdt, dydt)


def taylor_step(deriv_func, prey, pred, dt, a, b, c, d) -> float:
    """
    Take First Order Taylor Step  Forward given a dt and initial vals

    """

    results = deriv_func(dt, [prey, pred], a, b, c, d)
    prey_deriv = dt * results[0]
    pred_deriv = dt * results[1]

    return prey + prey_deriv, pred + pred_deriv


def euler_method(
    deriv_func,
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
            deriv_func, solution_prey[index], solution_pred[index], dt, a, b, c, d
        )

    ## Return Results ##
    return times, solution_prey, solution_pred


def make_plotting_folder():
    if not os.path.exists("Lab2/Plots"):
        os.mkdir("Lab2/Plots")
    return "Lab2/Plots"


def plot(model_data: modelrundata) -> None:
    fig, axes = plt.subplots(1)
    axes.plot(model_data.times, model_data.prey, color="b")
    axes.plot(model_data.times, model_data.pred, color="r")
    axes.set_ylabel("Population")
    axes.set_xlabel("Time (Years)")
    plot_dir = make_plotting_folder()
    plot_name = model_data.label + "_Plot.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    # plt.show(block=False)


def main():
    ## Set Model Parameters ##
    alpha = 1
    beta = 2
    delta = 1
    gamma = 3
    N1initial = 2
    N2initial = 6
    dT = 0.005
    t_stop = 1000

    ##################### Running Predator - Prey with two methods ############################################

    # Run Euler Method
    times, prey, pred = euler_method(
        simulation, N1initial, N2initial, alpha, beta, gamma, delta, dT, t_stop
    )
    predpreyeulerdata = modelrundata(times, prey, pred, "Predator_Prey_Euler")

    # Run RK8 Method
    predprey_ode = solve_ivp(
        simulation,
        [0, t_stop],
        [N1initial, N2initial],
        args=[alpha, beta, delta, gamma],
        method="DOP853",
        max_step=dT,
    )
    timepp, N1pp, N2pp = predprey_ode.t, predprey_ode.y[0, :], predprey_ode.y[1, :]
    predpreyRK8modeldata = modelrundata(timepp, N1pp, N2pp, "Predator_Prey_RK8")

    ##################### Running Competition with two methods############################################

    # Run Euler Method
    times, prey, pred = euler_method(
        dNdt_comp, N1initial, N2initial, alpha, beta, gamma, delta, dT, t_stop
    )
    compeulerdata = modelrundata(times, prey, pred, "Competition_Euler")

    # Run RK8 Method
    comp_ode = solve_ivp(
        dNdt_comp,
        [0, t_stop],
        [N1initial, N2initial],
        args=[alpha, beta, delta, gamma],
        method="DOP853",
        max_step=dT,
    )

    timecomp, N1comp, N2comp = comp_ode.t, comp_ode.y[0, :], comp_ode.y[1, :]
    compRK8modeldata = modelrundata(timecomp, N1comp, N2comp, "Competition_RK8")

    ## Plot Single Line Preliminary Analysis Model Data ##
    plot(predpreyeulerdata)
    plot(predpreyRK8modeldata)
    plot(compRK8modeldata)
    plot(compeulerdata)


if __name__ == "__main__":
    main()
