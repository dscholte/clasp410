# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 09:29:57 2023

@author: dscho
"""

import numpy as np
import matplotlib.pyplot as plt
from enum import Enum  # Enum ships with python standard library

plt.style.use("fivethirtyeight")


class ModelMode(Enum):
    hot_rods = 1
    greenland = 2


def run_heat(
    dt,
    dx,
    csquare,
    xmax,
    tmax,
    model_mode=ModelMode.hot_rods,
    is_hotrods_neumann_or_greenland_add_temps=0,
):
    """
    Parameters
    ----------
    xmax, tmax : float
        default to 1 and 0.2
    dt
        Time Step
    dx
        Spatial Step
    csquare
        Thermal Diffusivity Squared
    is_neumann
        Is the boundary conditions dirichlet (False) or neumann (True)
    Returns
    ----
    x: numpy vector
        - array of position locations
    t: numpy vector
        - array of time points
    temp: numpy 2D array
        - temperature as a function of time & space

    """

    if dt > ((dx**2) / (2 * csquare)):
        raise ValueError(
            "Stability Criterion not met"
            + f"dt={dt:6.2f}; dx={dx:6.2f}; csquare={csquare}"
        )

    # Set constant r
    r = csquare * dt / dx**2

    # Create space and time grids
    x = np.arange(0, xmax + dx, dx)
    t = np.arange(0, tmax + dt, dt)
    # Save number of points
    M, N = x.size, t.size

    # Temp solution array
    temp = np.zeros([M, N])
    if model_mode == ModelMode.hot_rods:
        temp[0, :] = 0
        temp[-1, :] = 0
        temp[:, 0] = 4 * x - 4 * (x**2)
    elif model_mode == ModelMode.greenland:
        temp[0, :] = 5
        temp[-1, :] = temp_kanger(t, is_hotrods_neumann_or_greenland_add_temps)
        temp[1:-1, 0] = 0

    # Solution to equation
    for j in range(0, N - 1):
        for i in range(1, M - 1):
            if (
                model_mode == ModelMode.hot_rods
                and is_hotrods_neumann_or_greenland_add_temps == 1
            ):
                temp[0, j] = temp[1, j]
                temp[-1, j] = temp[-2, j]
            temp[i, j + 1] = (1 - (2 * r)) * temp[i, j] + r * (
                temp[i + 1, j] + temp[i - 1, j]
            )

    return x, t, temp


def apply_dirichlet_and_validate_model():
    ## Validate Solver - Hot Rods ##
    size_of_rod = 1  # Meters
    spatial_step = 0.2  # Meters
    amount_of_time = 0.2  # Seconds
    time_step = 0.02  # Seconds
    thermal_diffusivity_squared = 1  # Constant
    x, t, temp = run_heat(
        time_step,
        spatial_step,
        thermal_diffusivity_squared,
        size_of_rod,
        amount_of_time,
        ModelMode.hot_rods,
        False,
    )

    sol10p3 = [
        [0.000000, 0.640000, 0.960000, 0.960000, 0.640000, 0.000000],
        [0.000000, 0.480000, 0.800000, 0.800000, 0.480000, 0.000000],
        [0.000000, 0.400000, 0.640000, 0.640000, 0.400000, 0.000000],
        [0.000000, 0.320000, 0.520000, 0.520000, 0.320000, 0.000000],
        [0.000000, 0.260000, 0.420000, 0.420000, 0.260000, 0.000000],
        [0.000000, 0.210000, 0.340000, 0.340000, 0.210000, 0.000000],
        [0.000000, 0.170000, 0.275000, 0.275000, 0.170000, 0.000000],
        [0.000000, 0.137500, 0.222500, 0.222500, 0.137500, 0.000000],
        [0.000000, 0.111250, 0.180000, 0.180000, 0.111250, 0.000000],
        [0.000000, 0.090000, 0.145625, 0.145625, 0.090000, 0.000000],
        [0.000000, 0.072812, 0.117813, 0.117813, 0.072812, 0.000000],
    ]
    # Convert to an array and transpose it to get correct ordering:
    sol10p3 = np.array(sol10p3).transpose()
    assert (sol10p3 - temp < 0.00001).all()
    plot_temp(
        x,
        t,
        temp,
        "Temperature (degC)",
        "Dirichlet Boundary Condition",
        "Depth (m)",
        "Temperature ($C$)",
        "inferno",
        False,
    )


# Neumann
def apply_neumann():
    x1, t1, temp1 = run_heat(0.0002, 0.02, 0.025, 1.0, 2, ModelMode.hot_rods, True)
    plot_temp(
        x1,
        t1,
        temp1,
        "Temperature (degC)",
        "Neumann Boundary Condition",
        "Depth (m)",
        "Temperature ($C$)",
        "inferno",
        False,
    )


def plot_greenland_ground_profile(title, x, t, temp):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.plot(
        temp[:, int(-365 / (t[1] - t[0])) :].min(axis=1),
        x,
        color="blue",
        label="Winter",
    )
    ax.plot(
        temp[:, int(-365 / (t[1] - t[0])) :].max(axis=1),
        x,
        "--",
        color="red",
        label="Summer",
    )
    ax.legend(loc="best")
    ax.set_ylabel("Depth (m)")
    ax.set_xlabel("Temperature (degC)")
    ax.set_title(title)
    ax.set_xlim([-7, 8])
    ax.set_ylim([-5, 105])
    fig.tight_layout()


def temp_kanger(t, warming):
    """
    For an array of times in days, return timeseries of temperature for
    Kangerlussuaq, Greenland.
    """
    t_kanger = np.array(
        [-19.7, -21.0, -17.0, -8.4, 2.3, 8.4, 10.7, 8.5, 3.1, -6.0, -12.0, -16.9]
    )
    t_amp = (t_kanger - t_kanger.mean()).max()

    return t_amp * np.sin(np.pi / 180 * t - np.pi / 2) + t_kanger.mean() + warming


def plot_temp(
    x,
    time,
    temp,
    xlabel="Time ($s$)",
    title="",
    ylabel="Distance ($m$)",
    clabel=r"Temperature ($^{\circ} C$)",
    cmap="inferno",
    inverty=False,
    **kwargs,
):
    """
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
    """
    fig, axes = plt.subplots(1, 1, figsize=(10, 8))
    map = axes.pcolor(time, x, temp, cmap=cmap, **kwargs)
    plt.colorbar(map, ax=axes, label=clabel)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title)


def run_general_permafrost_model():
    dt = 10
    dx = 1.0
    xmax = 100
    years = 50
    c_square = 0.25
    landc2 = c_square * (1 / 1000000) * 24 * 60 * 60
    add_temps = 0
    x, time, temp = run_heat(
        dt,
        dx,
        landc2,
        xmax,
        years * 365,
        ModelMode.greenland,
        add_temps,
    )
    maxtemp = np.abs(temp).max()

    plot_temp(
        x,
        time / 365.0,
        temp,
        xlabel="Time (Years)",
        ylabel="Depth ($m$)",
        cmap="seismic",
        vmin=-maxtemp,
        vmax=maxtemp,
        clabel="Temperature",
        title="Question 2",
    )
    plot_greenland_ground_profile("Ground Temperature", x, time, temp)
    return


# ------------------------------------------------
# Question 3


def run_general_permafrost_model_with_ghg_effect():
    dt = 10
    dx = 1.0
    nyear = 50
    xmax = 100
    c_square = 0.25
    landc2 = c_square * (1 / 1000000) * 24 * 60 * 60
    add_temps = 0.5
    xhalf, timehalf, temphalf = run_heat(
        dt,
        dx,
        landc2,
        xmax,
        nyear * 365,
        ModelMode.greenland,
        add_temps,
    )
    add_temps = 1
    xone, timeone, tempone = run_heat(
        dt,
        dx,
        landc2,
        xmax,
        nyear * 365,
        ModelMode.greenland,
        add_temps,
    )
    add_temps = 3
    xthree, timethree, tempthree = run_heat(
        dt,
        dx,
        landc2,
        xmax,
        nyear * 365,
        ModelMode.greenland,
        add_temps,
    )
    plot_greenland_ground_profile(
        "Ground Temperature warming 0.5degC", xhalf, timehalf, temphalf
    )
    plot_greenland_ground_profile(
        "Ground Temperature warming 1degC", xone, timeone, tempone
    )
    plot_greenland_ground_profile(
        "Ground Temperature warming 3degC", xthree, timethree, tempthree
    )
    return


def main():
    # ---------------------------------------------------------------
    # Question 1
    apply_dirichlet_and_validate_model()

    # Neumann for funzies
    apply_neumann()

    # -----------------------------------------------------------------
    # Question 2
    run_general_permafrost_model()

    # -----------------------------------------------------------------
    # Question 3
    run_general_permafrost_model_with_ghg_effect()
    plt.show()


if __name__ == "__main__":
    main()
