# -*- coding: utf-8 -*-

"""This file contains the functions that generate the current densities for the simulation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

# Importing functions
from modules.settings_modules import EIS_parameters

# General edition
mpl.rcParams['font.family'] = 'cmr10'  # 'cmr10' for English characters and 'DejaVu Serif' for French ones
plt.rcParams['axes.formatter.use_mathtext'] = True
plt.rcParams['axes.linewidth'] = 2.0
colors = mpl.colormaps['tab10']


# __________________________________________________Current densities___________________________________________________

def step_current(t, parameters):
    """Represents a step change in current density.
    The current starts at 0 and smoothly stabilizes at i_ini_step A.m-2 in delta_t_load seconds.
    Around t_switch seconds, the current increases smoothly and stabilizes at i_final_step A.m-2 in delta_t_load seconds.
    This is a C∞ function, which is advantageous for enhancing the overall stability of the results.

    Parameters:
    ----------
    t : float
        Time in seconds.
    parameters : dict
        A dictionary containing the parameters for the current density function.

    Returns:
    -------
    i_fc : float
        The step current density at time t.
    """

    # Initialisation
    t0_step, tf_step, delta_t_load_step, delta_t_dyn_step = parameters['t_step']  # (s, s, s, s). It is the initial, final,
    #                                                                     loading and dynamic time for display.
    i_ini_step, i_final_step = parameters['i_step']  # (A.m-2, A.m-2). It is the initial and final current density values.
    t_switch = tf_step // 2  # The current density value changes around this time.

    # Step current density
    i_fc = i_ini_step * (1.0 + math.tanh(4 * (t - 2 * (delta_t_load_step / 2)) / delta_t_load_step)) / 2 + \
           + (i_final_step - i_ini_step) * (1.0 + math.tanh(4 * (t - t_switch - (delta_t_load_step / 2)) / delta_t_load_step)) / 2

    return i_fc


def polarization_current(t, parameters):
    """Represents a current density used for creating a polarization curve.
    Starting from 0, the current density increases by the value of delta_i_pola every delta_t, following C∞ step current
    increments, until it reaches i_max_pola. Each increment lasts for delta_t_load_pola seconds. After each increment,
    there is a pause of delta_t_break_pola seconds to allow the stack to reach equilibrium.

    Parameters:
    ----------
    t : float
        Time in seconds.
    parameters : dict
        A dictionary containing the parameters for the current density function.

    Returns:
    -------
    i_fc : float
        The polarization current density at time t.
    """

    # Initialisation
    delta_t_load_pola, delta_t_break_pola, delta_i_pola, delta_t_ini_pola = parameters['delta_pola']
    #                                 (s, s, A.m-2, s). It is the loading time, the breaking time, the current density
    #                                 step, and the initial breaking time.
    i_max_pola = parameters['i_max_pola']  # A.m-2. It is the maximum current density for the polarization curve.
    delta_t = delta_t_load_pola + delta_t_break_pola  # s. It is the time of one load.
    tf = delta_t_ini_pola + int(i_max_pola / delta_i_pola + 1) * delta_t  # s. It is the duration of this polarization current.
    n = int(tf / delta_t)  # . It is the number of loads made for this polarization current.

    # Current density for the polarization curve
    i_fc = 0  # A.m-2. Initialisation of the current density.
    if t < delta_t_ini_pola:  # It is the initial break for having homogeneity inside the cell
        i_fc = 0  # before starting the measurements.
    else:
        for i in range(n):
            t_switch = delta_t * i  # The current density value changes around this time.
            i_fc += delta_i_pola * (1.0 + math.tanh(4 * (t - delta_t_ini_pola - delta_t - t_switch - (delta_t_load_pola / 2)) /
                                             delta_t_load_pola)) / 2

    return i_fc


def EIS_current(t, parameters):
    """
    Represents a current density used for creating an EIS curve and Bode diagrams.
    The current density is first equilibrated at i_EIS A.m-2 from 0 to t0_EIS seconds using a step increase.
    Then, a sinusoidal perturbation is added to the current density. This perturbation has an amplitude of
    (ratio_EIS * i_EIS) A.m-2 and a frequency of f[n_inf] Hz.

    Parameters:
    ----------
    t : float
        Time in seconds.
    parameters : dict
        A dictionary containing the parameters for the current density function.

    Returns:
    -------
    i_fc : float
        The polarization current density at time t.
    """

    # Initialisation
    i_EIS, ratio_EIS = parameters['i_EIS'], parameters['ratio_EIS']  # (A/m², ). i_EIS is the current for which a
    #                                                                  ratio_EIS perturbation is added.
    t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = parameters['t_EIS']  # It is the initial
    #         EIS time after stack equilibrium, a list of time parameters which gives the beginning of each frequency
    #         change, the final time, a list of time parameters which gives the estimated time for reaching equilibrium
    #         at each frequency, and a list of time parameters which gives the estimated time for measuring the voltage
    #         response at each frequency.
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = parameters['f_EIS']  # It is the power of the initial
    #         frequency: f_min_EIS = 10**f_power_min_EIS, the power of the final frequency, the number of frequencies
    #         tested and the number of points calculated per specific period.
    f = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)  # It is a list of all the frequency tested,
    #                                                      ranged logarithmically.

    # Current density for the EIS curve
    if t < t0_EIS:
        delta_t_ini = t0_EIS / 4  # s. It is the required time for elevating i_fc from 0 to i_EIS without starving the
        #                              cell.
        i_fc = i_EIS * (1.0 + math.tanh(4 * (t - 2 * (delta_t_ini / 2)) / delta_t_ini)) / 2
    else:
        n_inf = np.where(t_new_start_EIS <= t)[0][-1]  # It is the number of frequency changes which has been made so far.
        i_disruption = (ratio_EIS * i_EIS) * math.cos(2 * math.pi * f[n_inf] * t)
        i_fc = i_EIS + i_disruption

    return i_fc


# _________________________________________________Test of the program__________________________________________________

if __name__ == "__main__":
    fig, ax = plt.subplots(1, 3, figsize=(18, 6))

    # Tests for step_current curve:
    #   step_current parameters
    t_step = 0, 500, 30, 10  # (s, s, s, s). It is the initial, final, loading and dynamic time for display.
    i_step = 0.4e4, 0.8e4  # (A.m-2, A.m-2). It is the initial and final current density values.
    parameters = {'t_step': t_step, 'i_step': i_step}
    #   Display
    n = 10000
    t = np.linspace(t_step[0], t_step[1], n)
    i_fc = np.zeros(n)
    for i in range(n):
        i_fc[i] = step_current(t[i], parameters) / 10000  # Conversion in A/cm²
    ax[0].plot(t, i_fc, color=colors(0), label=r'$\mathregular{i_{fc}}$ (step)')
    ax[0].set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax[0].set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax[0].set_title('The step current density behaviour over time')
    ax[0].legend(loc='best')
    plt.show()

    # Tests for polarization curves:
    #   polarization_current parameters
    delta_pola = 30, 30, 0.1e4, 1 * 60  # (s, s, A.m-2, s). It is the loading time, the breaking time, the current
    #                                     density step, and the initial breaking time.
    i_max_pola = 1e4  # A.m-2. The maximum current density for the polarization curve.
    parameters = {'delta_pola': delta_pola, 'i_max_pola': i_max_pola}

    delta_t_load, delta_t_break_pola, delta_i_pola, delta_t_ini_pola = delta_pola
    t0 = 0
    tf = delta_t_ini_pola + int(i_max_pola / delta_i_pola + 1) * (delta_t_load + delta_t_break_pola)
    #   Display
    n = 10000
    t = np.linspace(t0, tf, n)
    i_fc = np.zeros(n)
    for i in range(n):
        i_fc[i] = polarization_current(t[i], parameters) / 10000  # Conversion in A/cm²
    ax[1].plot(t, i_fc, color=colors(1), label=r'$\mathregular{i_{fc}}$ (pola)')
    ax[1].set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax[1].set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax[1].set_title('The current density behaviour over time\nfor a polarization curve')
    ax[1].legend(loc='best')
    plt.show()

    # Tests for EIS curve:
    #   EIS_current parameters
    i_EIS, ratio_EIS = 1.0e4, 5 / 100  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
    f_EIS = -1.0, 1.0, 3, 50  # Frequency parameters for the EIS_current density function. It is a tuple containing the
    #                           power of the initial frequency 'f_power_min_EIS': f_min_EIS = 10**f_power_min_EIS, the
    #                           power of the final frequency 'f_power_max_EIS', the number of frequencies tested
    #                           'nb_f_EIS' and the number of points calculated per specific period 'nb_points_EIS'.
    t_EIS = EIS_parameters(f_EIS)  # It is the EIS parameters.
    parameters = {'i_EIS': i_EIS, 'ratio_EIS': ratio_EIS, 'f_EIS': f_EIS, 't_EIS': t_EIS}
    #   Display
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS
    t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
    f = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)  # It is a list of all the frequency tested.
    #        Step current for reaching i_EIS
    n = 1000
    t = np.linspace(0, t0_EIS, n)
    i_fc = np.zeros(n)
    for i in range(n):
        i_fc[i] = EIS_current(t[i], parameters) / 10000  # Conversion in A/cm²
    ax[2].plot(t, i_fc, color=colors(2))
    #        EIS current density
    for i in range(len(t_new_start_EIS)):  # The EIS curve is displayed for each frequency change in order to reduce the
        #                               number of points to calculate.
        t0_EIS_temp = t_new_start_EIS[i]
        tf_EIS_temp = t_new_start_EIS[i] + delta_t_break_EIS[i] + delta_t_measurement_EIS[i]
        n_inf = np.where(t_new_start_EIS <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes
        #                                                      which has been made so far.
        n = int(f[n_inf] * (tf_EIS_temp - t0_EIS_temp) * nb_points_EIS)  # It is the number of points to calculate.
        t = np.linspace(t0_EIS_temp, tf_EIS_temp, n)  # It is the time interval for the portion of the current density
        #                                              at a given frequency.
        i_fc = np.zeros(n)
        for j in range(n):
            i_fc[j] = EIS_current(t[j], parameters) / 10000  # Conversion in A/cm²
        ax[2].plot(t, i_fc, color=colors(2))
    ax[2].set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax[2].set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax[2].set_title('The current density behaviour over time\nfor an EIS curve')
    ax[2].legend([r'$\mathregular{i_{fc}}$ (EIS)'], loc='best')
    plt.show()
