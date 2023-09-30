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
import json

plt.style.use("fivethirtyeight")


@dataclass
class modelrundata:
    times: list[float]
    N1: list[float]
    N2: list[float]
    label: str

    def __post_init__(self):
        if len(self.N1) != len(self.times) or len(self.N2) != len(self.times):
            raise ValueError("Lists are not the same length! Doesn't make sense!")


# Competition equations
def dNdt_comp(t, variables, alpha, beta, delta, gamma):
    # N1 population
    N1 = variables[0]

    # N2 population
    N2 = variables[1]

    dn1dt = (alpha * N1 * (1 - N1)) - (beta * N1 * N2)
    dn2dt = (delta * N2 * (1 - N2)) - (gamma * N1 * N2)

    return (dn1dt, dn2dt)


# Predator - Prey equations
def simulation(t, variables, alpha, beta, delta, gamma):
    # Prey population
    N1 = variables[0]

    # Pred population
    N2 = variables[1]

    # derivatives
    dn1dt = (alpha * N1) - (beta * N1 * N2)
    dn2dt = (-delta * N2) + (gamma * N1 * N2)

    return (dn1dt, dn2dt)


def taylor_step(deriv_func, N1, N2, dt, a, b, c, d) -> float:
    """
    Take First Order Taylor Step  Forward given a dt and initial vals

    """

    results = deriv_func(dt, [N1, N2], a, b, c, d)
    N1_step = N1 + (dt * results[0])
    N2_step = N2 + (dt * results[1])

    return N1_step, N2_step


def euler_method(
    deriv_func,
    N_N1_init: float,
    N_N2_init: float,
    a: float,
    b: float,
    c: float,
    d: float,
    dt: float = 1,
    t_stop: float = 600,
) -> modelrundata:
    """
    Run Euler Method

    """

    ## Create Initial Arrays ##
    times = np.arange(0, t_stop, dt)
    solution_N1 = np.zeros(len(times))
    solution_N2 = np.zeros(len(times))

    ## Set Initial Conditions ##
    solution_N1[0] = N_N1_init
    solution_N2[0] = N_N2_init

    ## Run Euler Method with taylor_step function ##
    for index, value in enumerate(times[0:-1]):
        solution_N1[index + 1], solution_N2[index + 1] = taylor_step(
            deriv_func, solution_N1[index], solution_N2[index], dt, a, b, c, d
        )

    ## Return Results ##
    return times, solution_N1, solution_N2


def make_plotting_folder(save_id: str, iter_label: str):
    if not os.path.exists("Lab2"):
        os.mkdir("Lab2")
    if not os.path.exists("Lab2/Plots"):
        os.mkdir("Lab2/Plots")
    if not os.path.exists("Lab2/Plots/" + save_id):
        os.mkdir("Lab2/Plots/" + save_id)
    if not os.path.exists("Lab2/Plots/" + save_id + "/" + iter_label):
        os.mkdir("Lab2/Plots/" + save_id + "/" + iter_label)
    return "Lab2/Plots/" + save_id + "/" + iter_label


def single_plot(model_data: modelrundata, save_id: str, iter_label: str) -> None:
    fig, axes = plt.subplots(1)
    axes.plot(model_data.times, model_data.N1, color="b", label="N1 Population")
    axes.plot(model_data.times, model_data.N2, color="r", label="N2 Population")
    axes.legend(loc="upper left")
    axes.set_ylabel("Population")
    axes.set_xlabel("Time (Years)")
    plot_dir = make_plotting_folder(save_id, iter_label)
    plot_name = "Single_" + model_data.label + "_Plot.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    plt.close()


def comparison_plots(
    euler_model_data: modelrundata,
    RK8_model_data: modelrundata,
    save_id: str,
    iter_label: str,
    label: str,
) -> None:
    fig, axes = plt.subplots(1)
    axes.plot(
        euler_model_data.times,
        euler_model_data.N1,
        color="b",
        linestyle="dotted",
        label="Euler N1 Population",
    )
    axes.plot(
        euler_model_data.times,
        euler_model_data.N2,
        color="r",
        linestyle="dotted",
        label="Euler N2 Population",
    )
    axes.plot(
        RK8_model_data.times,
        RK8_model_data.N1,
        color="b",
        linestyle="solid",
        label="RK8 N1 Population",
    )
    axes.plot(
        RK8_model_data.times,
        RK8_model_data.N2,
        color="r",
        linestyle="solid",
        label="Rk8 N2 Population",
    )
    axes.legend(loc="upper left")
    axes.set_ylabel("Population")
    axes.set_xlabel("Time (Years)")
    plot_dir = make_plotting_folder(save_id, iter_label)
    plot_name = "Comparison_RK8_Euler" + label + "_Plot.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    plt.close()


def phase_plot(model_data: modelrundata, save_id: str, iter_label: str) -> None:
    fig, axes = plt.subplots(1)
    axes.plot(model_data.N1, model_data.N2, color="b", label="Phase - Prey & Pred")
    axes.set_ylabel("Predator Species")
    axes.set_xlabel("Prey Species")
    plot_dir = make_plotting_folder(save_id, iter_label)
    plot_name = "Phase_" + model_data.label + "_Plot.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    plt.close()


def run_models(model_params):
    save_id = model_params["save_id"]
    iter_label = model_params["iter_label"]
    alpha = model_params["alpha"]
    beta = model_params["beta"]
    delta = model_params["delta"]
    gamma = model_params["gamma"]
    N1initial = model_params["N1initial"]
    N2initial = model_params["N2initial"]
    t_step_comp = model_params["t_step_comp"]
    t_step_preypred = model_params["t_step_preypred"]
    t_stop = model_params["t_stop"]

    ##################### Run Predator - Prey ODE Solvers ############################################

    # Run Euler Method
    times, N1, N2 = euler_method(
        simulation,
        N1initial,
        N2initial,
        alpha,
        beta,
        delta,
        gamma,
        t_step_preypred,
        t_stop,
    )
    ppeulerdata = modelrundata(times, N1, N2, "Prey_Pred_Euler")

    # Run RK8 Method
    N2de = solve_ivp(
        simulation,
        [0, t_stop],
        [N1initial, N2initial],
        args=[alpha, beta, delta, gamma],
        method="DOP853",
        max_step=t_step_preypred,
    )
    timepp, N1pp, N2pp = N2de.t, N2de.y[0, :], N2de.y[1, :]
    ppRK8data = modelrundata(timepp, N1pp, N2pp, "Prey_Pred_RK8")

    ##################### Run Competition ODE Solvers ############################################

    # Run Euler Method
    times, N1, N2 = euler_method(
        dNdt_comp,
        N1initial,
        N2initial,
        alpha,
        beta,
        delta,
        gamma,
        t_step_comp,
        t_stop,
    )
    compeulerdata = modelrundata(times, N1, N2, "Competition_Euler")

    # Run RK8 Method
    comp_ode = solve_ivp(
        dNdt_comp,
        [0, t_stop],
        [N1initial, N2initial],
        args=[alpha, beta, delta, gamma],
        method="DOP853",
        max_step=t_step_comp,
    )

    timecomp, N1comp, N2comp = comp_ode.t, comp_ode.y[0, :], comp_ode.y[1, :]
    compRK8data = modelrundata(timecomp, N1comp, N2comp, "Competition_RK8")

    ##################### Create Outputs ############################################

    ## Export Paramaters File ##
    export_dir = make_plotting_folder(save_id, iter_label)
    with open(os.path.join(export_dir, "parameters.json"), "w") as f:
        json.dump(model_params, f, indent=4)

    ## Plot Single Line Preliminary Analysis Model Data ##
    single_plot(ppeulerdata, save_id, iter_label)
    single_plot(ppRK8data, save_id, iter_label)
    single_plot(compeulerdata, save_id, iter_label)
    single_plot(compRK8data, save_id, iter_label)

    ## Plot Comparison Plots ##
    comparison_plots(ppeulerdata, ppRK8data, save_id, iter_label, "PredPrey")
    comparison_plots(compeulerdata, compRK8data, save_id, iter_label, "Comparison")

    ## Plot Phase Diagrams ##
    phase_plot(ppeulerdata, save_id, iter_label)
    phase_plot(ppRK8data, save_id, iter_label)
    return None


def main():
    """
    Let's start going through the lab! We'll go through each question in the lab.
    """

    default_model_params = {
        "save_id": "Basic",
        "iter_label": "General",
        "alpha": 1,
        "beta": 2,
        "delta": 1,
        "gamma": 3,
        "N1initial": 0.3,
        "N2initial": 0.6,
        "t_step_comp": 1,
        "t_step_preypred": 0.05,
        "t_stop": 100,
    }
    ## Question 1: The first thing was to show we could reproduce the figures given! ##
    model_params = default_model_params.copy()
    run_models(model_params)

    ## Question 1 Part 2: We need to vary the timestep and see what happens! ##
    model_params["save_id"] = "varytimestep"
    for time_step in np.arange(1, 4, 0.5):
        model_params["iter_label"] = "ts_" + str(time_step)
        model_params["t_step_comp"] = time_step
        model_params["t_step_preypred"] = time_step
        run_models(model_params)

    ## Question 1 Part 2: We need to vary the timestep and see what happens! ##
    model_params = default_model_params.copy()
    model_params["save_id"] = "varytimestep"
    for time_step in np.arange(1, 4, 0.5):
        model_params["iter_label"] = "ts_" + str(time_step)
        model_params["t_step_comp"] = time_step
        model_params["t_step_preypred"] = time_step
        run_models(model_params)

    ## Question 2 & 3: Vary Initial Conditions!! ##
    model_params = default_model_params.copy()
    model_params["save_id"] = "VaryInitConds!"
    for N1_step in np.arange(0.1, 1, 0.1):
        model_params["iter_label"] = "N1_" + str(round(N1_step, 2))
        model_params["N1initial"] = N1_step
        run_models(model_params)
    for N2_step in np.arange(0.1, 1, 0.1):
        model_params["iter_label"] = "N2_" + str(round(N2_step, 2))
        model_params["N2initial"] = N2_step
        run_models(model_params)


if __name__ == "__main__":
    main()
