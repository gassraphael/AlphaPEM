# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import time
import numpy as np

# Importing constants' value and functions
from configuration.settings import current_density_parameters, operating_inputs, physical_parameters, \
                                   computing_parameters
from model.AlphaPEM import AlphaPEM
from modules.main_modules import figures_preparation, plot_saving

# __________________________________________________AlphaPEM settings___________________________________________________
"""
- Users can select various preconfigured configurations for execution.
- Adjustments to these configurations can be made within setting.py or current_densities.py and their associated files.
- Selecting different type_fuel_cell during a single run results in simultaneous plots for various configurations.
"""

if __name__ == '__main__':
    # Starting time
    start_time = time.time()

    # Fuel cell possibilities: "EH-31_1.5"(2021), "EH-31_2.0"(2021), "EH-31_2.25"(2021), "EH-31_2.5"(2021),
    #                          "BX_1.0"(2015), "BX_1.35"(2015), "LF"(2010), "manual_setup".
    # This parameter includes the fuel cell used in the model and the corresponding operating conditions.
    # - EH-31 is a fuel cell developed by EH GROUP. 1.5, 2.0, 2.25 and 2.5 corresponds to the different pressure options.
    # - BX corresponds to the fuel cell used in Biao Xie's work: https://doi.org/10.1016/j.ijheatmasstransfer.2022.122705.
    # 1.0 and 1.35 corresponds to the different pressure options.
    # - LF corresponds to the fuel cell used in Linhao Fan work: http://dx.doi.org/10.1016/j.enconman.2017.08.034.
    type_fuel_cell_1 = "EH-31_2.0"
    type_fuel_cell_2 = None
    type_fuel_cell_3 = None
    type_fuel_cell_4 = None
    # Current density possibilities: "step", "polarization", "EIS".
    type_current = "polarization"
    # Auxiliary system possibilities: "closed_anode_with_recirculation", "opened_anode", "no_auxiliary".
    type_auxiliary = "opened_anode"
    # Control strategy for the operating inputs: "Phi_des", "no_control".
    type_control_1 = "no_control"
    type_control_2 = "no_control"
    type_control_3 = "no_control"
    type_control_4 = "no_control"
    # Purges possibilities: "constant_purge", "periodic_purge", "no_purge".
    type_purge = "no_purge"
    # Display possibilities: "multiple", "synthetic", "no_display".
    type_display = "synthetic"
    # Plot possibilities: "dynamic", "fixed". Using dynamic plot option enables real-time figure updates during program
    # execution, albeit at the cost of decreased program speed.
    type_plot = "fixed"

# ___________________________________Retrieving parameters from the settings.py file____________________________________
    """This should remain unaltered for regular program usage."""

    # Imposed inputs
    t_step, i_step, delta_pola, i_EIS, ratio_EIS, f_EIS, t_EIS, current_density = \
        current_density_parameters(type_current)
    # Operating conditions
    Tfc_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, i_max_pola_1 = \
        operating_inputs(type_fuel_cell_1)
    if type_fuel_cell_2 is not None:
        Tfc_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, i_max_pola_2 = \
            operating_inputs(type_fuel_cell_2)
    if type_fuel_cell_3 is not None:
        Tfc_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3, i_max_pola_3 = \
            operating_inputs(type_fuel_cell_3)
    if type_fuel_cell_4 is not None:
        Tfc_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4, i_max_pola_4 = \
            operating_inputs(type_fuel_cell_4)
    # Physical parameters
    Hcl, epsilon_mc, tau, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hgc, Wgc, Lgc, Aact, e, Re, i0_c_ref, kappa_co, \
        kappa_c, a_slim, b_slim, a_switch, C_dl = physical_parameters(type_fuel_cell_1)
    # Computing parameters
    max_step, n_gdl, t_purge = computing_parameters(type_current, Hgdl, Hcl)

# __________________________________________________________Main________________________________________________________
    """This section is dedicated to ensuring the proper execution of the simulator, considering all the various 
    possibilities including real-time figure updates and simultaneous plotting for different configurations. 
    This should remain unaltered for regular program usage.
    """

    # Check if the type_current is valid
    if type_current != "step" and type_current != "polarization" and type_current != "EIS":
        raise ValueError('You have to specify a type_current which is accepted.')

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(type_current, type_display)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if type_plot == "dynamic":

        # Check if the type_fuel_cell and type_current are valid
        if type_fuel_cell_2 is not None or type_fuel_cell_3 is not None or type_fuel_cell_4 is not None:
            raise ValueError('dynamic plot is not currently intended for use with different inputs.')
        if type_current == "step" and type_display == "multiple":
            raise ValueError('dynamic plot is not thought to be used with step current and multiple display.' +
                             'There would be too much plots to handle.')

        # Initialization
        #       ... of the plot update number (n) and the initial time interval (time_interval)
        initial_variable_values = None
        if type_current == "step":
            t0_step, tf_step, delta_t_load_step, delta_t_dyn_step = t_step
            n = int(tf_step / delta_t_dyn_step)  # It is the plot update number.
            time_interval = [0, delta_t_dyn_step]  # It is the initial time interval.
        elif type_current == "polarization":
            delta_t_load_pola, delta_t_break_pola, delta_i_pola, delta_t_ini_pola = delta_pola
            delta_t = delta_t_load_pola + delta_t_break_pola  # s. It is the time of one load.
            tf = delta_t_ini_pola + int(i_max_pola_1 / delta_i_pola + 1) * delta_t  # s. It is the polarization current
            #                                                                         duration.
            n = int(tf / delta_t)  # It is the plot update number.
            time_interval = [0, delta_t_ini_pola + delta_t]  # It is the initial time interval.
        elif type_current == "EIS":
            t0_EIS, t_new_start, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
            f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS  # These are used for EIS max_step actualization.
            f = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)  # It is a list of all the frequency tested.
            n = len(t_new_start)  # It is the plot update number.
            time_interval = [0, t0_EIS]  # It is the initial time interval.

        #       A preliminary simulation run is necessary to equilibrate the internal variables of the cell at i_EIS
        #       prior to initiating the EIS.
        if type_current == "EIS":
            Simulator1 = AlphaPEM(current_density, Tfc_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                                  t_step, i_step, i_max_pola_1, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                  Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref,
                                  kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge,
                                  type_fuel_cell_1, type_current, type_auxiliary, type_control_1, type_purge,
                                  "no_display", type_plot, initial_variable_values, time_interval)

            # time_interval actualization
            t0_EIS_temp = t0_EIS  # It is the initial time for 1 EIS point.
            tf_EIS_temp = t_new_start[0] + delta_t_break_EIS[0] + delta_t_measurement_EIS[0]  # It is the final time for
            #                                                                                  1 EIS point.
            n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which has been
            #                                                      made.
            max_step = 1 / (f[n_inf] * nb_points_EIS)  # max_step is actualized according to the current frequency
            #                                        for increased calculation
            time_interval = [t0_EIS_temp, tf_EIS_temp]

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator1.solver_variable_names:
                initial_variable_values.append(Simulator1.variables[x][-1])

            if type_display == "multiple":
                print("A display bug prevents the dynamic updating of the graphs, as it appears that too much data is "
                      "involved. However, the data is correctly calculated, and the appropriate plots are saved in the "
                      "'results' folder. This display bug does not occur when using a 'synthetic' type_display.")

        # Dynamic simulation
        for i in range(n):
            Simulator1 = AlphaPEM(current_density, Tfc_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                                  t_step, i_step, i_max_pola_1, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                  Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref,
                                  kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge,
                                  type_fuel_cell_1, type_current, type_auxiliary, type_control_1, type_purge,
                                  type_display, type_plot, initial_variable_values, time_interval)

            # time_interval actualization
            if i < (n - 1):  # The final simulation does not require actualization.
                if type_current == "step":
                    t0_interval = Simulator1.variables['t'][-1]
                    tf_interval = (i + 2) * delta_t_dyn_step
                    time_interval = [t0_interval, tf_interval]  # Reset of the time interval
                elif type_current == "polarization":
                    t0_interval = Simulator1.variables['t'][-1]
                    tf_interval = delta_t_ini_pola + (i + 2) * delta_t
                    time_interval = [t0_interval, tf_interval]  # Reset of the time interval
                elif type_current == "EIS":
                    t0_EIS_temp = Simulator1.variables['t'][-1]  # It is the initial time for 1 EIS point.
                    tf_EIS_temp = t_new_start[i + 1] + delta_t_break_EIS[i + 1] + delta_t_measurement_EIS[i + 1]  # It
                    #                                                               is the final time for 1 EIS point.
                    n_inf = np.where(t_new_start <= t0_EIS_temp)[0][-1]  # It is the number of frequency changes which
                    #                                                      has been made.
                    max_step = 1 / (f[n_inf] * nb_points_EIS)  # max_step is actualized according to the current frequency
                    #                                        for increased calculation
                    time_interval = [t0_EIS_temp, tf_EIS_temp]  # It is the time interval for 1 EIS point.

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = []
            for x in Simulator1.solver_variable_names:
                initial_variable_values.append(Simulator1.variables[x][-1])

            # Display
            if type_display != "no_display":
                Simulator1.Display(ax1, ax2, ax3)

    else:  # elif type_plot == "fixed":

        # Certain conditions must be met.
        if type_current == "step" and \
                (type_fuel_cell_2 is not None or type_fuel_cell_3 is not None or type_fuel_cell_4 is not None):
            raise ValueError('step current contains too much information for a common plot with different inputs.')
        elif type_current == "EIS":
            raise ValueError('EIS has to be plot with a dynamic type_plot setting, '
                             'because max_step has to be adjusted at each frequency.')

        # Simulation
        Simulator1 = AlphaPEM(current_density, Tfc_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, t_step,
                              i_step, i_max_pola_1, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl,
                              Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co,
                              kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell_1,
                              type_current, type_auxiliary, type_control_1, type_purge, type_display, type_plot)
        if type_fuel_cell_2 is not None:
            Simulator2 = AlphaPEM(current_density, Tfc_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2,
                                  t_step, i_step, i_max_pola_2, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                  Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref,
                                  kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge,
                                  type_fuel_cell_2, type_current, type_auxiliary, type_control_2, type_purge,
                                  type_display, type_plot)
        if type_fuel_cell_3 is not None:
            Simulator3 = AlphaPEM(current_density, Tfc_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3,
                                  t_step, i_step, i_max_pola_3, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                  Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref,
                                  kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge,
                                  type_fuel_cell_3, type_current, type_auxiliary, type_control_3, type_purge,
                                  type_display, type_plot)
        if type_fuel_cell_4 is not None:
            Simulator4 = AlphaPEM(current_density, Tfc_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4,
                                  t_step, i_step, i_max_pola_4, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl,
                                  Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref,
                                  kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl, max_step, n_gdl, t_purge,
                                  type_fuel_cell_4, type_current, type_auxiliary, type_control_4, type_purge,
                                  type_display, type_plot)

        # Display
        if type_display != "no_display":
            Simulator1.Display(ax1, ax2, ax3)
            if type_fuel_cell_2 is not None:
                Simulator2.Display(ax1, ax2, ax3)
            if type_fuel_cell_3 is not None:
                Simulator3.Display(ax1, ax2, ax3)
            if type_fuel_cell_4 is not None:
                Simulator4.Display(ax1, ax2, ax3)

    # Plot saving
    plot_saving(type_fuel_cell_1, type_current, type_display, fig1, fig2, fig3)

    # Ending time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)
