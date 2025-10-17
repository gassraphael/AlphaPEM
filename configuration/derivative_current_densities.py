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

def dstep_currentdt(t, parameters):
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
    return 2 * i_ini * (1.0 - math.tanh(4 * (t - (delta_t_load_step / 2)) / (delta_t_load_step / 2)) ** 2) + \
           2 * (i_step - i_ini) * (1.0 - math.tanh(4 * (t - delta_t_ini_step - (delta_t_load_step / 2)) /
                                                   (delta_t_load_step / 2)) ** 2)


def dpolarization_currentdt(t, parameters):
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

    raise ValueError("The function 'dpolarization_currentdt' has been disabled.")

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


def dpolarization_current_for_calibrationdt(t, parameters):
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

    raise ValueError("The function 'dpolarization_current_for_calibrationdt' has been disabled.")

    # Extraction of the parameters
    #   Initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola_cali = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']  # (s).
    #   Loading time for one step current of the polarisation current density function.
    delta_t_load_pola_cali = parameters['pola_current_for_cali_parameters']['delta_t_load_pola_cali']  # (s).
    #   Breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola_cali = parameters['pola_current_for_cali_parameters']['delta_t_break_pola_cali']  # (s).
    type_fuel_cell, voltage_zone = parameters['type_fuel_cell'], parameters['voltage_zone']  # The fuel cell for which the calibration is performed.
    i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(type_fuel_cell, voltage_zone)  # (A.m-2, V). It is the experimental
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


def dEIS_currentdt(t, parameters):
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

    raise ValueError("The function 'dEIS_currentdt' has been disabled.")

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