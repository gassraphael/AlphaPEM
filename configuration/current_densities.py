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
from calibration.experimental_values import pola_exp_values_calibration

# General edition
mpl.rcParams['font.family'] = 'cmr10'  # 'cmr10' for English characters and 'DejaVu Serif' for French ones
plt.rcParams['axes.formatter.use_mathtext'] = True
plt.rcParams['axes.linewidth'] = 2.0
colors = mpl.colormaps['tab10']


# __________________________________________________Current densities___________________________________________________

def step_current(t, parameters):
    """ This function represents a step current density experiment. For the first delta_t_ini_step seconds, the current
    density is set to i_ini A.m-2 to allow the internal states of the fuel cell to stabilise. Then, the current density
    increases from i_ini to i_step A.m-2 in a step change over delta_t_load seconds. Finally, the current density
    remains at i_step A.m-2.
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

    # Extraction of the parameters
    #   Initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_step = parameters['step_current_parameters']['delta_t_ini_step'] # (s).
    #   Initial current density used for the stabilisation of the internal states.
    i_ini = 1.0e4  # (A.m-2). This is the standard value for the initialisation.
    #   Loading time for the step current density function, from 0 to i_step.
    delta_t_load_step = parameters['step_current_parameters']['delta_t_load_step'] # (s).
    #   Current density for the step current density function.
    i_step = parameters['step_current_parameters']['i_step'] # (A.m-2).

    # Step current density
    return i_ini + (i_step - i_ini) * (1.0 + math.tanh(4 * (t - delta_t_ini_step - (delta_t_load_step / 2)) /
                                                       (delta_t_load_step / 2))) / 2


def polarization_current(t, parameters):
    """ This function represents a current density used for creating a polarization curve. For the first
    delta_t_ini_step seconds, the current density is set to i_ini A.m-2 to allow the internal states of the fuel cell to
    stabilise. Then, the current density increases by the value of delta_i_pola every delta_t, following C∞ step current
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

    # Extraction of the parameters
    #   Initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola = parameters['pola_current_parameters']['delta_t_ini_pola'] # (s).
    #   Loading time for one step current of the polarisation current density function.
    delta_t_load_pola =parameters['pola_current_parameters']['delta_t_load_pola'] # (s).
    #   Breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola = parameters['pola_current_parameters']['delta_t_break_pola'] # (s).
    #   Current density step for the polarisation current density function.
    delta_i_pola = parameters['pola_current_parameters']['delta_i_pola'] # (A.m-2).
    #   Maximum current density for the polarization curve.
    i_max_pola = parameters['pola_current_parameters']['i_max_pola'] # (A.m-2).

    # Calculation of the time parameters
    #   Time of one step.
    delta_t = delta_t_load_pola + delta_t_break_pola  # (s).
    #   Duration of this polarization curve.
    tf = delta_t_ini_pola + int(i_max_pola / delta_i_pola + 1) * delta_t  # (s).
    #   Number of loads made for this polarization curve.
    n = int(tf / delta_t)

    # Current density for the polarization curve
    i_fc = 0  # A.m-2. Initialisation of the current density.
    for i in range(n):
        i_fc += delta_i_pola * (1.0 + math.tanh(4 * (t - delta_t_ini_pola - i * delta_t - (delta_t_load_pola / 2)) /
                                                (delta_t_load_pola / 2))) / 2
    return i_fc


def polarization_current_for_calibration(t, parameters):
    """This function represents a current density used for creating a polarisation curve dedicated to the calibration of
    a specific fuel cell. The principle is similar to the polarization_current function, but it uses experimental values
    for the current density load.

    Parameters:
    ----------
    t : float
        Time in seconds.
    parameters : dict
        A dictionary containing the parameters for the current density function.

    Returns:
    -------
    i_fc : float
        The polarisation current density at time t.
    """

    # Extraction of the parameters
    #   Initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola_cali = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']  # (s).
    #   Loading time for one step current of the polarisation current density function.
    delta_t_load_pola_cali = parameters['pola_current_for_cali_parameters']['delta_t_load_pola_cali']  # (s).
    #   Breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola_cali = parameters['pola_current_for_cali_parameters']['delta_t_break_pola_cali']  # (s).
    type_fuel_cell = parameters['type_fuel_cell']  # The fuel cell for which the calibration is performed.
    i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(type_fuel_cell)  # (A.m-2, V). It is the experimental
    #                                                           current density and voltage values for the calibration.

    # Calculation of the time parameters
    #   Time of one step.
    delta_t = delta_t_load_pola_cali + delta_t_break_pola_cali  # (s).

    # Current density for the polarization curve used for calibration
    i_fc = 0  # A.m-2. Initialisation of the current density.
    for e in range(len(i_exp_cali_t)):
        if e == 0:
            i_fc += i_exp_cali_t[0] * (1.0 + math.tanh(4 * (t - delta_t_ini_pola_cali - (delta_t_load_pola_cali / 2)) /
                                                       (delta_t_load_pola_cali / 2))) / 2
        else:
            delta_i_exp_cali = (i_exp_cali_t[e] - i_exp_cali_t[e-1]) # (A.m-2). It is the difference between the
            #                                                  current density at the current step and the previous one.
            i_fc += delta_i_exp_cali * (1.0 + math.tanh(4 * (t - delta_t_ini_pola_cali - e * delta_t - (delta_t_load_pola_cali / 2)) /
                                                        (delta_t_load_pola_cali / 2))) / 2
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
        delta_t_ini = 3*60 # s. It is the required time for elevating i_fc from 0 to i_EIS without starving the cell.
        i_fc = i_EIS * (1.0 + math.tanh(4 * (t - 2 * (delta_t_ini / 2)) / delta_t_ini)) / 2
    else:
        n_inf = np.where(t_new_start_EIS <= t)[0][-1]  # It is the number of frequency changes which has been made so far.
        i_disruption = (ratio_EIS * i_EIS) * math.cos(2 * math.pi * f[n_inf] * t)
        i_fc = i_EIS + i_disruption

    return i_fc


# _________________________________________________Test of the program__________________________________________________

if __name__ == "__main__":
    fig, ax = plt.subplots(1, 4, figsize=(24, 6))

    # Tests for step_current curve:
    #   Step current parameters
    delta_t_ini_step = 120*60 # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_step = 3 * 60 # (s). Loading time for the step current density function, from 0 to i_step.
    delta_t_break_step = 30 * 60 # (s). Time at i_step current density for the stabilization of the internal states.
    i_step = 1.5e4 # (A.m-2). Current density for the step current density function.
    step_current_parameters = {'delta_t_ini_step': delta_t_ini_step, 'delta_t_load_step': delta_t_load_step,
                               'delta_t_break_step': delta_t_break_step, 'i_step': i_step}
    parameters = {'step_current_parameters': step_current_parameters}
    #   Display
    n = 10000
    t = np.linspace(0, delta_t_ini_step + delta_t_load_step + delta_t_break_step, n)
    i_fc_step = np.zeros(n)
    for i in range(n):
        i_fc_step[i] = step_current(t[i], parameters) / 10000  # Conversion in A/cm²
    ax[0].plot(t, i_fc_step, color=colors(0), label=r'$\mathregular{i_{fc}}$ (step)')
    ax[0].set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax[0].set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax[0].set_title('The step current density behaviour over time')
    ax[0].legend(loc='best')
    plt.show()


    # Tests for polarization curves:
    #   Polarization current parameters
    delta_t_ini_pola = 120 * 60 # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_pola = 30 # (s). Loading time for one step current of the polarisation current density function.
    delta_t_break_pola = 15 * 60 # (s). Breaking time for one step current, for the stabilisation of the internal states.
    delta_i_pola = 0.1e4 # (A.m-2). Current density step for the polarisation current density function.
    i_max_pola = 2.5e4 # (A.m-2). Maximum current density for the polarization curve.
    pola_current_parameters = {'delta_t_ini_pola': delta_t_ini_pola, 'delta_t_load_pola': delta_t_load_pola,
                               'delta_t_break_pola': delta_t_break_pola, 'delta_i_pola': delta_i_pola,
                               'i_max_pola': i_max_pola}
    parameters = {'pola_current_parameters': pola_current_parameters}
    #   Display
    n = 10000
    tf = delta_t_ini_pola + int(i_max_pola / delta_i_pola) * (delta_t_load_pola + delta_t_break_pola)
    t = np.linspace(0, tf, n)
    i_fc_pola = np.zeros(n)
    for i in range(n):
        i_fc_pola[i] = polarization_current(t[i], parameters) / 1e4  # Conversion in A/cm²
    ax[1].plot(t, i_fc_pola, color=colors(1), label=r'$\mathregular{i_{fc}}$ (pola)')
    ax[1].set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax[1].set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax[1].set_title('The current density behaviour over time\nfor a polarization curve')
    ax[1].legend(loc='best')
    plt.show()


    # Tests for calibration curves:
    #   Polarization current for calibration parameters
    delta_t_ini_pola_cali = 120 * 60 # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_pola_cali = 30 # (s). Loading time for one step current of the polarisation current density function.
    delta_t_break_pola_cali = 15 * 60 # (s). Breaking time for one step current, for the stabilisation of the internal states.
    pola_current_for_cali_parameters = {'delta_t_ini_pola_cali': delta_t_ini_pola_cali,
                                        'delta_t_load_pola_cali': delta_t_load_pola_cali,
                                        'delta_t_break_pola_cali': delta_t_break_pola_cali}
    parameters = {'pola_current_for_cali_parameters': pola_current_for_cali_parameters, 'type_fuel_cell': "EH-31_2.0"}
    i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(parameters['type_fuel_cell'])  # (A.m-2, V). Experimental
    #                                                            current density and voltage values for the calibration.
    #   Display
    n = 10000
    tf = delta_t_ini_pola_cali + len(i_exp_cali_t) * (delta_t_load_pola_cali + delta_t_break_pola_cali)
    t = np.linspace(0, tf, n)
    i_fc_cali = np.zeros(n)
    for i in range(n):
        i_fc_cali[i] = polarization_current_for_calibration(t[i], parameters) / 1e4  # Conversion in A/cm²
    ax[2].plot(t, i_fc_cali, color=colors(2), label=r'$\mathregular{i_{fc}}$ (pola)')
    ax[2].set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax[2].set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax[2].set_title('The current density behaviour over time\nfor a polarization curve (calibration)')
    ax[2].legend(loc='best')
    plt.show()


    # Tests for EIS curve:
    #   EIS_current parameters
    i_EIS, ratio_EIS = 1.0e4, 5 / 100  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
    f_EIS = -3, 5, 90, 50  # Frequency parameters for the EIS_current density function. It is a tuple containing the
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
    i_fc_EIS = np.zeros(n)
    for i in range(n):
        i_fc_EIS[i] = EIS_current(t[i], parameters) / 10000  # Conversion in A/cm²
    ax[3].plot(t, i_fc_EIS, color=colors(3))
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
        i_fc_EIS = np.zeros(n)
        for j in range(n):
            i_fc_EIS[j] = EIS_current(t[j], parameters) / 10000  # Conversion in A/cm²
        ax[3].plot(t, i_fc_EIS, color=colors(3))
    ax[3].set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax[3].set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax[3].set_title('The current density behaviour over time\nfor an EIS curve')
    ax[3].legend([r'$\mathregular{i_{fc}}$ (EIS)'], loc='best')
    plt.show()
