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


# Object that holds model results
@dataclass
class modelrundata:
    """
    Class to hold all the information required from the ODE Solver.
    Variables
    ---------
    times
        It's a list of the times at each point
    N1
        The population of the N1 Species (Prey for the PreyPred Sim) at each time (t) in the times list
    N2
        The population of the N2 Species (Pred for the PreyPred Sim) at each time (t) in the times list
    label
        Informative label for the plots

    """

    times: list[float]
    N1: list[float]
    N2: list[float]
    label: str

    def __post_init__(self):
        """
        Sanity Check when putting data in this class"""
        if len(self.N1) != len(self.times) or len(self.N2) != len(self.times):
            raise ValueError("Lists are not the same length! Doesn't make sense!")


# Competition ODEs
def dNdt_comp(
    t: float,
    variables: list[float],
    alpha: float,
    beta: float,
    delta: float,
    gamma: float,
) -> tuple[float, float]:
    """ "
    The Competition ODEs
    Parameters
    ----------
    variables
        List of the N1 and N2 Populations
    alpha, beta, delta, gamma
        Input variables to the Competition ODEs
    Returns
    -------
    dn1dt, dn2dt
        Derivative at the next time step for N1, N2
    """
    # N1 population
    N1 = variables[0]

    # N2 population
    N2 = variables[1]

    dn1dt = (alpha * N1 * (1 - N1)) - (beta * N1 * N2)
    dn2dt = (delta * N2 * (1 - N2)) - (gamma * N1 * N2)

    return (dn1dt, dn2dt)


# Predator - Prey equations
def simulation(
    t: float,
    variables: list[float],
    alpha: float,
    beta: float,
    delta: float,
    gamma: float,
) -> tuple[float, float]:
    """ "
    The Predator Prey ODEs
    Parameters
    ----------
    t
        time (Unused)
    variables
        List of the N1 (Prey) and N2 (Pred) Populations
    alpha, beta, delta, gamma
        Input variables to the Competition ODEs
    Returns
    -------
    dn1dt, dn2dt
        Derivative at the next time step for N1, N2
    """
    # Prey population
    N1 = variables[0]

    # Pred population
    N2 = variables[1]

    # derivatives
    dn1dt = (alpha * N1) - (beta * N1 * N2)
    dn2dt = (-delta * N2) + (gamma * N1 * N2)

    return (dn1dt, dn2dt)


def taylor_step(
    deriv_func, N1: float, N2: float, dt: float, a: float, b: float, c: float, d: float
) -> float:
    """
    Take First Order Taylor Step Forward given a dt and initial vals
    Parameters
    ----------
    deriv_func
        The Derivative Function to use (Either Compeition ODE or Predator-Prey ODE)
    N1
        Initial N1 Population
    N2
        Initial N2 Population
    dt
        Timestep (Max time step for RK8)
    a,b,c,d
        Known Parameters for ODEs
    Returns
    -------
    float
        The taylor step forward in time
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
    Calls taylor step function
    Taylor function calls derivative of either equation

    Parameters
    ----------
    deriv_func
        The Derivative Function to use (Either Compeition ODE or Predator-Prey ODE)
    N_N1_init
        Initial N1 Population
    N_N2_init
        Initial N2 Population
    dt
        Timestep (Max time step for RK8)
    t_stop
        End of time (Default is a 100 years)
    a,b,c,d
        Known Parameters for ODEs
    Returns
    -------
    times, N1, N2
        The output from the euler method, time, populations for both species
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

    ## Return Results , times and populations ##
    return times, solution_N1, solution_N2


# Creates folders for plots
def make_plotting_folder(save_id: str, iter_label: str) -> str:
    """
    Creates/Checks the required plotting folders
    Parameters
    ----------
    save_id
        First level Folder for what you are varying
    iter_label
        Second Level Folder for iterating through initial conditions
    Returns
    -------
    The plotting directory you are dumping into.
    """
    if not os.path.exists("Lab2"):
        os.mkdir("Lab2")
    if not os.path.exists("Lab2/Plots"):
        os.mkdir("Lab2/Plots")
    if not os.path.exists("Lab2/Plots/" + save_id):
        os.mkdir("Lab2/Plots/" + save_id)
    if not os.path.exists("Lab2/Plots/" + save_id + "/" + iter_label):
        os.mkdir("Lab2/Plots/" + save_id + "/" + iter_label)
    return "Lab2/Plots/" + save_id + "/" + iter_label


# Plots model data and save it to directory
def single_plot(model_data: modelrundata, save_id: str, iter_label: str) -> None:
    """
    Plots one model data in a specific folder
    Parameters
    ----------
    model_data
        The model data
    save_id, iter_label
        Where to save plot
    Returns
    -------
    None
    """
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


# compares runga-kutta vs euler
def comparison_plots(
    euler_model_data: modelrundata,
    RK8_model_data: modelrundata,
    save_id: str,
    iter_label: str,
    label: str,
    title: str,
) -> None:
    """
        Plots  model data against each other in a specific folder
    Parameters
    ----------
    euler model_data, RK8_model_data
        The model data for the models to compare
    save_id, iter_label
        Where to save plot
    label
        What type, should be "Competition" or "Predator Prey"
    title
        Title of plot
    Returns
    -------
    None
    """
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
    plot_dir = make_plotting_folder(
        save_id, iter_label
    )  # make sure you have plotting folder
    plot_name = "Comparison_RK8_Euler" + label + "_Plot.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    plt.close()


def phase_plot(model_data: modelrundata, save_id: str, iter_label: str) -> None:
    """
        Plots one model data of Prey Population against Predator Population
    Parameters
    ----------
    model_data
        The model data
    save_id, iter_label
        Where to save plot
    Returns
    -------
    None
    """
    fig, axes = plt.subplots(1)
    axes.plot(model_data.N1, model_data.N2, color="b", label="Phase - Prey & Pred")
    axes.set_ylabel("Predator Species")
    axes.set_xlabel("Prey Species")
    plot_dir = make_plotting_folder(save_id, iter_label)
    plot_name = "Phase_" + model_data.label + "_Plot.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    plt.close()


# Takes in parameters a,b,c,d and others and runs euler and RK8
def run_models(
    model_params,
) -> tuple[modelrundata, modelrundata, modelrundata, modelrundata]:
    """
    Runs the Euler and RK8 method for both competition and predator prey
        Parameters
    ----------
    model_params
        Holds all model_params, listed here
    save_id
        What folder to save in
    iter_label
        What sub folder to save in
    alpha, beta, delta, gamma
        ODE Parameters
    N1initial, N2initial
        Intial Population for both species
    t_step_comp
        Time Step For Competition
    t_step_preypred
        Time Step for Predator Prey
    t_stop
        What time to stop at (Usually 100 years)
    Returns
    -------
    Tuple (4 Items)
        All of the model run data of all 4 options
    """
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
    comparison_plots(
        ppeulerdata,
        ppRK8data,
        save_id,
        iter_label,
        "PredPrey",
        "Lotka-Volterra Predator-Prey Model",
    )
    comparison_plots(
        compeulerdata,
        compRK8data,
        save_id,
        iter_label,
        "Comparison",
        "Lotka-Volterra Competition Model",
    )

    ## Plot Phase Diagrams ##
    phase_plot(ppeulerdata, save_id, iter_label)
    phase_plot(ppRK8data, save_id, iter_label)
    return ppeulerdata, ppRK8data, compeulerdata, compRK8data


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
    # for time_step in np.arange(1, 4, 0.5):
    #     model_params["iter_label"] = "ts_" + str(time_step)
    #     model_params["t_step_comp"] = time_step
    #     model_params["t_step_preypred"] = time_step
    #     run_models(model_params)

    ## Question 2 & 3: Vary Initial Conditions!! ##
    model_params = default_model_params.copy()
    model_params["save_id"] = "VaryInitConds!"
    modelsruns = []
    for N1_step in np.arange(0.1, 1, 0.1):
        model_params["iter_label"] = "N1_" + str(round(N1_step, 2))
        model_params["N1initial"] = N1_step
        results1 = run_models(model_params)
        modelsruns.append(results1)
    for N2_step in np.arange(0.1, 1, 0.1):
        model_params["iter_label"] = "N2_" + str(round(N2_step, 2))
        model_params["N2initial"] = N2_step
        run_models(model_params)
    population_prey_max = []
    population_pred_max = []
    for i in modelsruns:
        population_prey_max.append(max(i[1].N1))
        population_pred_max.append(max(i[1].N2))
    N1_Values = np.arange(0.1, 1, 0.1)
    fig, axes = plt.subplots(1)
    axes.plot(N1_Values, population_prey_max, color="b", label="N1 Max Population")
    axes.plot(N1_Values, population_pred_max, color="r", label="N2 Max Population")
    axes.legend(loc="upper left")
    axes.set_ylabel("Population")
    axes.set_xlabel("N1")
    plot_dir = make_plotting_folder(
        model_params["save_id"], "Overall"
    )  # make sure you have plotting folder
    plot_name = "Max_Population_Plot.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    plt.close()
    # plot_dir = make_plotting_folder(save_id, iter_label)

    # Pedator Prey ModelComparoisons
    colors = ["b", "g", "r", "c", "m", "y", "k", "#A52A2A", "#FF7F00"]
    fig, axes = plt.subplots(1)
    for counter, i in enumerate(modelsruns):
        axes.plot(
            i[1].N1, i[1].N2, color=colors[counter], label=round(N1_Values[counter], 2)
        )
    axes.set_ylabel("Predator Species")
    axes.set_xlabel("Prey Species")
    axes.legend()
    plot_dir = make_plotting_folder(
        model_params["save_id"], "Overall"
    )  # make sure you have plotting folder
    plot_name = "AllPreyPredPhaseDiagrams.pdf"
    plt.savefig(os.path.join(plot_dir, plot_name), format="pdf", bbox_inches="tight")
    plt.close()

    plt.show()


if __name__ == "__main__":
    main()
