# -*- coding: utf-8 -*-

"""This module contains some of the required functions for the main.py file.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import time
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Importing constants' value and functions
from configuration.settings import current_density_parameters, operating_inputs_function, physical_parameters, \
                                   computing_parameters
from model.AlphaPEM import AlphaPEM
from calibration.experimental_values import pola_exp_values, pola_exp_values_calibration

# _____________________________________________________Main modules_____________________________________________________

def figures_preparation(type_current, type_display):
    """ This function create the required figures and axes according to the type_current and type_display.

    Parameters
    ----------
    type_current : str
        Type of current density function.
    type_display : str
        Type of display.

    Returns
    -------
    fig1 : matplotlib.figure.Figure
        Figure for the first plot.
    ax1 : matplotlib.axes._subplots.AxesSubplot
        Axes for the first plot.
    fig2 : matplotlib.figure.Figure
        Figure for the second plot.
    ax2 : matplotlib.axes._subplots.AxesSubplot
        Axes for the second plot.
    """

    mpl.rcParams['font.family'] = 'cmr10'  # 'cmr10' for English characters and 'DejaVu Serif' for French ones
    mpl.rcParams['axes.formatter.use_mathtext'] = True  # For the scientific notation
    mpl.rcParams['lines.linewidth'] = 2.0
    mpl.rcParams['lines.markersize'] = 5.0

    if type_display == "no_display":
        fig1, ax1 = None, None
        fig2, ax2 = None, None
        fig3, ax3 = None, None

    # For the step current
    if type_current == "step":
        if type_display == "multiple":  # saving instruction is directly implemented within AlphaPEM.Display here.
            mpl.rcParams['font.size'] = 18  # Font size for all text
            fig1, ax1 = None, None  # Here, additional plots are unnecessary
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
        elif type_display == "synthetic":
            mpl.rcParams['font.size'] = 13  # Font size for all text
            fig1, ax1 = plt.subplots(3, 3, figsize=(14, 14))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
            plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)

    # For the polarization curve
    elif type_current == "polarization":
        if type_display == "multiple":
            mpl.rcParams['font.size'] = 11  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3, figsize=(14, 4.7))
            fig2, ax2 = plt.subplots(1, 3, figsize=(14, 4.7))
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
            plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        elif type_display == "synthetic":
            mpl.rcParams['font.size'] = 18  # Font size for all text
            fig1, ax1 = plt.subplots(figsize=(8, 8))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary

    # For the polarization curve used for the calibration
    elif type_current == "polarization_for_cali":
        if type_display == "multiple":
            mpl.rcParams['font.size'] = 11  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3, figsize=(14, 4.7))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
            plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        elif type_display == "synthetic":
            mpl.rcParams['font.size'] = 18  # Font size for all text
            fig1, ax1 = plt.subplots(figsize=(8, 8))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary

    # For the EIS curve
    elif type_current == "EIS":
        if type_display == "multiple":
            mpl.rcParams['font.size'] = 18  # Font size for all text
            fig1, ax1 = plt.subplots(figsize=(8, 8))
            fig2, ax2 = plt.subplots(figsize=(8, 8))
            fig3, ax3 = plt.subplots(figsize=(8, 8))
        elif type_display == "synthetic":
            mpl.rcParams['font.size'] = 13  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3, figsize=(14, 4.7))
            fig2, ax2 = None, None  # Here, additional plots are unnecessary
            fig3, ax3 = None, None  # Here, additional plots are unnecessary
            plt.subplots_adjust(left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)

    return fig1, ax1, fig2, ax2, fig3, ax3


def launch_AlphaPEM_for_step_current(type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3, type_fuel_cell_4,
                                     type_current, type_auxiliary, type_control_1, type_control_2, type_control_3,
                                     type_control_4, type_purge, type_display, type_plot):
    """Launch the AlphaPEM simulator for a step current density and display the results.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration (computing parameter).
    type_current : str
        Type of current density function (computing parameter).
    type_auxiliary : str
        Type of auxiliary system (computing parameter).
    type_control : str
        Type of control system (computing parameter).
    type_purge : str
        Type of purge system (computing parameter).
    type_display : str
        Type of display (computing parameter).
    type_plot : str
        Type of plot (computing parameter).
    """

    # Starting time
    start_time = time.time()

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(type_current, type_display)

    # Retrieving parameters from the settings.py file
    #   Imposed inputs
    (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
     current_density) = current_density_parameters(type_current)
    #   Operating conditions
    T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, pola_current_parameters_1 = \
        operating_inputs_function(pola_current_parameters, type_fuel_cell_1)
    if type_fuel_cell_2 is not None:
        T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, pola_current_parameters_2 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_2)
    if type_fuel_cell_3 is not None:
        T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3, pola_current_parameters_3 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_3)
    if type_fuel_cell_4 is not None:
        T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4, pola_current_parameters_4 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_4)
    #   Physical parameters
    (Hcl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, i0_c_ref, kappa_co, kappa_c,
     a_slim, b_slim, a_switch, C_scl) = physical_parameters(type_fuel_cell_1)
    #   Computing parameters
    n_gdl, t_purge, step_current_parameters = computing_parameters(step_current_parameters, type_current, Hgdl, Hcl)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if type_plot == "dynamic":

        # Certain conditions must be met.
        if type_fuel_cell_2 is not None or type_fuel_cell_3 is not None or type_fuel_cell_4 is not None:
            raise ValueError('dynamic plot is not currently intended for use with different inputs.')
        if type_current == "step" and type_display == "multiple":
            raise ValueError('dynamic plot is not thought to be used with step current and multiple display.' +
                             'There would be too much plots to handle.')

        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = None
        #           Extraction of the parameters
        tf_step = step_current_parameters['delta_t_ini_step'] + step_current_parameters['delta_t_load_step'] + \
                  step_current_parameters['delta_t_break_step']  # (s).
        delta_t_dyn_step = step_current_parameters['delta_t_dyn_step']  # (s).
        #           Calculation
        n = int(tf_step / delta_t_dyn_step)  # It is the plot update number.
        time_interval = [0, delta_t_dyn_step]  # (s). It is the initial time interval.

        # Dynamic simulation
        for i in range(n):
            Simulator_1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                                  step_current_parameters, pola_current_parameters_1, pola_current_for_cali_parameters,
                                  i_EIS, ratio_EIS,
                                  t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, epsilon_mc,
                                  epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl,
                                  n_gdl, t_purge, type_fuel_cell_1, type_current, type_auxiliary, type_control_1,
                                  type_purge, type_display, type_plot, initial_variable_values, time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                t0_interval = Simulator_1.variables['t'][-1]
                tf_interval = (i + 2) * delta_t_dyn_step
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator_1.solver_variable_names:
                initial_variable_values.append(Simulator_1.variables[x][-1])

            # Display
            if type_display != "no_display":
                Simulator_1.Display(ax1, ax2, ax3)

    else:  # elif type_plot == "fixed":

        # Certain conditions must be met.
        if type_current == "step" and \
                (type_fuel_cell_2 is not None or type_fuel_cell_3 is not None or type_fuel_cell_4 is not None):
            raise ValueError('step current contains too much information for a common plot with different inputs.')

        # Simulation
        Simulator_1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                               step_current_parameters, pola_current_parameters_1, pola_current_for_cali_parameters,
                               i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                               epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl,
                               n_gdl, t_purge, type_fuel_cell_1, type_current, type_auxiliary, type_control_1,
                               type_purge, type_display, type_plot)
        if type_fuel_cell_2 is not None:
            Simulator_2 = AlphaPEM(current_density, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2,
                                   step_current_parameters, pola_current_parameters_2, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch,
                                   C_scl, n_gdl, t_purge, type_fuel_cell_2, type_current, type_auxiliary,
                                   type_control_2, type_purge, type_display, type_plot)
        if type_fuel_cell_3 is not None:
            Simulator_3 = AlphaPEM(current_density, T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3,
                                   step_current_parameters, pola_current_parameters_3, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch,
                                   C_scl, n_gdl, t_purge, type_fuel_cell_3, type_current, type_auxiliary,
                                   type_control_3, type_purge, type_display, type_plot)
        if type_fuel_cell_4 is not None:
            Simulator_4 = AlphaPEM(current_density, T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4,
                                   step_current_parameters, pola_current_parameters_4, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch,
                                   C_scl, n_gdl, t_purge, type_fuel_cell_4, type_current, type_auxiliary,
                                   type_control_4, type_purge, type_display, type_plot)
        # Display
        if type_display != "no_display":
            Simulator_1.Display(ax1, ax2, ax3)
            if type_fuel_cell_2 is not None:
                Simulator_2.Display(ax1, ax2, ax3)
            if type_fuel_cell_3 is not None:
                Simulator_3.Display(ax1, ax2, ax3)
            if type_fuel_cell_4 is not None:
                Simulator_4.Display(ax1, ax2, ax3)

    # Plot saving
    Simulator_1.Save_plot(fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)


def launch_AlphaPEM_for_polarization_current(type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3, type_fuel_cell_4,
                                             type_current, type_auxiliary, type_control_1, type_control_2,
                                             type_control_3, type_control_4, type_purge, type_display, type_plot):
    """Launch the AlphaPEM simulator for a polarization current density and display the results.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration (computing parameter).
    type_current : str
        Type of current density function (computing parameter).
    type_auxiliary : str
        Type of auxiliary system (computing parameter).
    type_control : str
        Type of control system (computing parameter).
    type_purge : str
        Type of purge system (computing parameter).
    type_display : str
        Type of display (computing parameter).
    type_plot : str
        Type of plot (computing parameter).
    """

    # Starting time
    start_time = time.time()

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(type_current, type_display)

    # Retrieving parameters from the settings.py file
    #   Imposed inputs
    (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
     current_density) = current_density_parameters(type_current)
    #   Operating conditions
    T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, pola_current_parameters_1 = \
        operating_inputs_function(pola_current_parameters, type_fuel_cell_1)
    if type_fuel_cell_2 is not None:
        T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, pola_current_parameters_2 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_2)
    if type_fuel_cell_3 is not None:
        T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3, pola_current_parameters_3 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_3)
    if type_fuel_cell_4 is not None:
        T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4, pola_current_parameters_4 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_4)
    #   Physical parameters
    (Hcl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, i0_c_ref, kappa_co, kappa_c,
     a_slim, b_slim, a_switch, C_scl) = physical_parameters(type_fuel_cell_1)
    #   Computing parameters
    n_gdl, t_purge, step_current_parameters = computing_parameters(step_current_parameters, type_current, Hgdl, Hcl)

    # Condition to fill for the comparison with experimental values
    if type_fuel_cell_1 is not None and type_fuel_cell_1 != "manual_setup" and \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":  # Experimental points are accessible
        i_exp_t_1, U_exp_t_1 = pola_exp_values(type_fuel_cell_1)
        if pola_current_parameters_1['i_max_pola'] < i_exp_t_1[-1]:
            raise ValueError('The given maximum current density of the polarization curve i_max_pola_1 is lower than the '
                             'maximum current density of the experimental values. Please increase it.')
    if type_fuel_cell_2 is not None and type_fuel_cell_2 != "manual_setup" and \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":  # Experimental points are accessible
        i_exp_t_2, U_exp_t_2 = pola_exp_values(type_fuel_cell_2)
        if pola_current_parameters_2['i_max_pola'] < i_exp_t_2[-1]:
            raise ValueError('The given maximum current density of the polarization curve i_max_pola_2 is lower than the '
                             'maximum current density of the experimental values. Please increase it.')
    if type_fuel_cell_3 is not None and type_fuel_cell_3 != "manual_setup" and \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":  # Experimental points are accessible
        i_exp_t_3, U_exp_t_3 = pola_exp_values(type_fuel_cell_3)
        if pola_current_parameters_3['i_max_pola'] < i_exp_t_3[-1]:
            raise ValueError('The given maximum current density of the polarization curve i_max_pola_3 is lower than the '
                             'maximum current density of the experimental values. Please increase it.')
    if type_fuel_cell_4 is not None and type_fuel_cell_4 != "manual_setup" and \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":  # Experimental points are accessible
        i_exp_t_4, U_exp_t_4 = pola_exp_values(type_fuel_cell_4)
        if pola_current_parameters_4['i_max_pola'] < i_exp_t_4[-1]:
            raise ValueError('The given maximum current density of the polarization curve i_max_pola_4 is lower than the '
                             'maximum current density of the experimental values. Please increase it.')

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if type_plot == "dynamic":

        # Certain conditions must be met.
        if type_fuel_cell_2 is not None or type_fuel_cell_3 is not None or type_fuel_cell_4 is not None:
            raise ValueError('dynamic plot is not currently intended for use with different inputs.')

        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = None
        #           Extraction of the parameters
        delta_t_ini_pola = pola_current_parameters_1['delta_t_ini_pola']  # (s).
        delta_t_load_pola = pola_current_parameters_1['delta_t_load_pola']  # (s).
        delta_t_break_pola = pola_current_parameters_1['delta_t_break_pola']  # (s).
        delta_i_pola = pola_current_parameters_1['delta_i_pola']  # (A.m-2).
        i_max_pola_1 = pola_current_parameters_1['i_max_pola_1']  # (A.m-2).
        #           Calculation
        delta_t_pola = delta_t_load_pola + delta_t_break_pola  # s. It is the time of one load.
        tf = delta_t_ini_pola + int(
            i_max_pola_1 / delta_i_pola) * delta_t_pola  # s. It is the polarization current duration.
        n = int(tf / delta_t_pola)  # It is the plot update number.
        time_interval = [0, delta_t_ini_pola + delta_t_pola]  # It is the initial time interval.

        # Dynamic simulation
        for i in range(n):
            Simulator_1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                                   step_current_parameters, pola_current_parameters_1, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_1, type_current, type_auxiliary,
                                   type_control_1, type_purge, type_display, type_plot, initial_variable_values,
                                   time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                t0_interval = Simulator_1.variables['t'][-1]
                tf_interval = delta_t_ini_pola + (i + 2) * delta_t_pola
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator_1.solver_variable_names:
                initial_variable_values.append(Simulator_1.variables[x][-1])

            # Display
            if type_display != "no_display":
                Simulator_1.Display(ax1, ax2, ax3)

    else:  # elif type_plot == "fixed":
        # Simulation
        Simulator_1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                               step_current_parameters, pola_current_parameters_1, pola_current_for_cali_parameters,
                               i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                               epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch,
                               C_scl, n_gdl, t_purge, type_fuel_cell_1, type_current, type_auxiliary, type_control_1,
                               type_purge, type_display, type_plot)
        if type_fuel_cell_2 is not None:
            Simulator_2 = AlphaPEM(current_density, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2,
                                   step_current_parameters, pola_current_parameters_2, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_2, type_current, type_auxiliary,
                                   type_control_2, type_purge, type_display, type_plot)
        if type_fuel_cell_3 is not None:
            Simulator_3 = AlphaPEM(current_density, T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3,
                                   step_current_parameters, pola_current_parameters_3, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_3, type_current, type_auxiliary,
                                   type_control_3, type_purge, type_display, type_plot)
        if type_fuel_cell_4 is not None:
            Simulator_4 = AlphaPEM(current_density, T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4,
                                   step_current_parameters, pola_current_parameters_4, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_4, type_current, type_auxiliary,
                                   type_control_4, type_purge, type_display, type_plot)

        # Display
        if type_display != "no_display":
            Simulator_1.Display(ax1, ax2, ax3)
            if type_fuel_cell_2 is not None:
                Simulator_2.Display(ax1, ax2, ax3)
            if type_fuel_cell_3 is not None:
                Simulator_3.Display(ax1, ax2, ax3)
            if type_fuel_cell_4 is not None:
                Simulator_4.Display(ax1, ax2, ax3)

    # Plot saving
    Simulator_1.Save_plot(fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)


def launch_AlphaPEM_for_polarization_current_for_calibration(type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3,
                                                             type_fuel_cell_4, type_current, type_auxiliary,
                                                             type_control_1, type_control_2, type_control_3,
                                                             type_control_4, type_purge, type_display, type_plot):
    """Launch the AlphaPEM simulator for a polarization current density made for the calibration of the undetermined
    parameters, and display the results.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration (computing parameter).
    type_current : str
        Type of current density function (computing parameter).
    type_auxiliary : str
        Type of auxiliary system (computing parameter).
    type_control : str
        Type of control system (computing parameter).
    type_purge : str
        Type of purge system (computing parameter).
    type_display : str
        Type of display (computing parameter).
    type_plot : str
        Type of plot (computing parameter).
    """

    # Starting time
    start_time = time.time()

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(type_current, type_display)

    # Retrieving parameters from the settings.py file
    #   Imposed inputs
    (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
     current_density) = current_density_parameters(type_current)
    #   Operating conditions
    T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, pola_current_parameters_1 = \
        operating_inputs_function(pola_current_parameters, type_fuel_cell_1)
    if type_fuel_cell_2 is not None:
        T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, pola_current_parameters_2 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_2)
    if type_fuel_cell_3 is not None:
        T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3, pola_current_parameters_3 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_3)
    if type_fuel_cell_4 is not None:
        T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4, pola_current_parameters_4 = \
            operating_inputs_function(pola_current_parameters, type_fuel_cell_4)
    #   Physical parameters
    (Hcl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, i0_c_ref, kappa_co, kappa_c,
     a_slim, b_slim, a_switch, C_scl) = physical_parameters(type_fuel_cell_1)
    #   Computing parameters
    n_gdl, t_purge, step_current_parameters = computing_parameters(step_current_parameters, type_current,
                                                                             Hgdl, Hcl)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if type_plot == "dynamic":

        # Certain conditions must be met.
        if type_fuel_cell_2 is not None or type_fuel_cell_3 is not None or type_fuel_cell_4 is not None:
            raise ValueError('dynamic plot is not currently intended for use with different inputs.')
        if type_current == "polarization_for_cali":
            raise ValueError('calibration should not use dynamic plot, as it is not intended for real-time display.')

        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = None
        #           Extraction of the parameters
        delta_t_ini_pola_cali = pola_current_parameters_1['delta_t_ini_pola_cali']  # (s).
        delta_t_load_pola_cali = pola_current_parameters_1['delta_t_load_pola_cali']  # (s).
        delta_t_break_pola_cali = pola_current_parameters_1['delta_t_break_pola_cali']  # (s).
        i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(type_fuel_cell_1)  # (A.m-2, V).
        #           Calculation
        delta_t_pola_cali = delta_t_load_pola_cali + delta_t_break_pola_cali  # s. It is the time of one load.
        tf = delta_t_ini_pola_cali + len(
            i_exp_cali_t) * delta_t_pola_cali  # s. It is the polarization current duration.
        n = int(tf / delta_t_pola_cali)  # It is the plot update number.
        time_interval = [0, delta_t_ini_pola_cali + delta_t_pola_cali]  # It is the initial time interval.

        # Dynamic simulation
        for i in range(n):
            Simulator_1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                                   step_current_parameters, pola_current_parameters_1, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_1, type_current, type_auxiliary,
                                   type_control_1, type_purge, type_display, type_plot, initial_variable_values,
                                   time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                t0_interval = Simulator_1.variables['t'][-1]
                tf_interval = delta_t_ini_pola_cali + (i + 2) * delta_t_pola_cali
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator_1.solver_variable_names:
                initial_variable_values.append(Simulator_1.variables[x][-1])

            # Display
            if type_display != "no_display":
                Simulator_1.Display(ax1, ax2, ax3)

    else:  # elif type_plot == "fixed":

        # Certain conditions must be met.
        if type_current == "polarization_for_cali" and (type_fuel_cell_1 == "manual_setup" or \
                                                        type_auxiliary != "forced-convective_cathode_with_flow-through_anode"):
            raise ValueError('polarization current for calibration should be done with experimental data.')

        # Simulation
        Simulator_1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                               step_current_parameters, pola_current_parameters_1, pola_current_for_cali_parameters,
                               i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                               epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch,
                               C_scl, n_gdl, t_purge, type_fuel_cell_1, type_current, type_auxiliary, type_control_1,
                               type_purge, type_display, type_plot)
        if type_fuel_cell_2 is not None:
            Simulator_2 = AlphaPEM(current_density, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2,
                                   step_current_parameters, pola_current_parameters_2, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_2, type_current, type_auxiliary,
                                   type_control_2, type_purge, type_display, type_plot)
        if type_fuel_cell_3 is not None:
            Simulator_3 = AlphaPEM(current_density, T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3,
                                   step_current_parameters, pola_current_parameters_3, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_3, type_current, type_auxiliary,
                                   type_control_3, type_purge, type_display, type_plot)
        if type_fuel_cell_4 is not None:
            Simulator_4 = AlphaPEM(current_density, T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4,
                                   step_current_parameters, pola_current_parameters_4, pola_current_for_cali_parameters,
                                   i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                   epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim,
                                   a_switch, C_scl, n_gdl, t_purge, type_fuel_cell_4, type_current, type_auxiliary,
                                   type_control_4, type_purge, type_display, type_plot)

        # Display
        if type_display != "no_display":
            Simulator_1.Display(ax1, ax2, ax3)
            if type_fuel_cell_2 is not None:
                Simulator_2.Display(ax1, ax2, ax3)
            if type_fuel_cell_3 is not None:
                Simulator_3.Display(ax1, ax2, ax3)
            if type_fuel_cell_4 is not None:
                Simulator_4.Display(ax1, ax2, ax3)

    # Plot saving
    Simulator_1.Save_plot(fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)


def launch_AlphaPEM_for_EIS_current(type_fuel_cell, type_current, type_auxiliary, type_control, type_purge,
                                    type_display, type_plot):
    """Launch the AlphaPEM simulator for an EIS current density and display the results.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration (computing parameter).
    type_current : str
        Type of current density function (computing parameter).
    type_auxiliary : str
        Type of auxiliary system (computing parameter).
    type_control : str
        Type of control system (computing parameter).
    type_purge : str
        Type of purge system (computing parameter).
    type_display : str
        Type of display (computing parameter).
    type_plot : str
        Type of plot (computing parameter).
    """

    # Starting time
    start_time = time.time()

    # Check if the type_current is valid
    if type_plot != "dynamic":
        raise ValueError('EIS has to be plot with a dynamic type_plot setting, '
                         'because max_step has to be adjusted at each frequency.')

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(type_current, type_display)

    # Retrieving parameters from the settings.py file
    #   Imposed inputs
    (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
     current_density) = current_density_parameters(type_current)
    #   Operating conditions
    T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, pola_current_parameters = \
        operating_inputs_function(pola_current_parameters, type_fuel_cell)
    #   Physical parameters
    (Hcl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, i0_c_ref, kappa_co, kappa_c,
     a_slim, b_slim, a_switch, C_scl) = physical_parameters(type_fuel_cell)
    #   Computing parameters
    n_gdl, t_purge, step_current_parameters = computing_parameters(step_current_parameters, type_current, Hgdl, Hcl)

    # Initialization
    #       Calculation of the plot update number (n) and the initial time interval (time_interval).
    initial_variable_values = None
    t0_EIS, t_new_start, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS  # These are used for EIS max_step
    #                                                                    actualization.
    f = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)  # It is a list of all the frequency tested.
    n = len(t_new_start)  # It is the plot update number.
    time_interval = [0, t0_EIS]  # It is the initial time interval.

    #       A preliminary simulation run is necessary to equilibrate the internal variables of the cell at i_EIS
    #       prior to initiating the EIS.
    Simulator = AlphaPEM(current_density, T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des,
                         step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS,
                         ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, epsilon_mc,
                         epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl, n_gdl, t_purge,
                         type_fuel_cell, type_current, type_auxiliary, type_control, type_purge,"no_display",
                         type_plot, initial_variable_values, time_interval)

    # time_interval actualization
    t0_EIS_temp = t0_EIS  # It is the initial time for 1 EIS point.
    tf_EIS_temp = t_new_start[0] + delta_t_break_EIS[0] + delta_t_measurement_EIS[0]  # It is the final time for
    #                                                                                  1 EIS point.
    n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which has been
    #                                                      made.
    time_interval = [t0_EIS_temp, tf_EIS_temp]

    # Recovery of the internal states from the end of the preceding simulation.
    initial_variable_values = []
    for x in Simulator.solver_variable_names:
        initial_variable_values.append(Simulator.variables[x][-1])

    if type_display == "multiple":
        print("A display bug prevents the dynamic updating of the graphs, as it appears that too much data is "
              "involved. However, the data is correctly calculated, and the appropriate plots are saved in the "
              "'results' folder. This display bug does not occur when using a 'synthetic' type_display.")

    # Dynamic simulation
    for i in range(n):
        Simulator = AlphaPEM(current_density, T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des,
                             step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters,
                             i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                             epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl,
                             n_gdl, t_purge, type_fuel_cell, type_current, type_auxiliary, type_control, type_purge,
                             type_display, type_plot, initial_variable_values, time_interval)

        # time_interval actualization
        if i < (n - 1):  # The final simulation does not require actualization.
            t0_EIS_temp = Simulator.variables['t'][-1]  # It is the initial time for 1 EIS point.
            tf_EIS_temp = t_new_start[i + 1] + delta_t_break_EIS[i + 1] + delta_t_measurement_EIS[i + 1]  # It
            #                                                               is the final time for 1 EIS point.
            n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which
            #                                                      has been made.
            time_interval = [t0_EIS_temp, tf_EIS_temp]  # It is the time interval for 1 EIS point.

        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values = []
        for x in Simulator.solver_variable_names:
            initial_variable_values.append(Simulator.variables[x][-1])

        # Display
        if type_display != "no_display":
            Simulator.Display(ax1, ax2, ax3)

    # Plot saving
    Simulator.Save_plot(fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)
