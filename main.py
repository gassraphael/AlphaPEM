# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
import copy
from configuration.settings import current_density_parameters, operating_inputs_function, physical_parameters, \
                                   computing_parameters
from modules.main_modules import (select_nth_elements, launch_AlphaPEM_for_step_current,
                                  launch_AlphaPEM_for_polarization_current,
                                  launch_AlphaPEM_for_polarization_current_for_calibration,
                                  launch_AlphaPEM_for_EIS_current)

# __________________________________________________AlphaPEM settings___________________________________________________
"""
- Users can select various preconfigured configurations for execution.
- Adjustments to these configurations can be made within setting.py or current_densities.py and their associated files.
- Selecting different type_fuel_cell during a single run results in simultaneous plots for various configurations.
"""

if __name__ == '__main__':
    # Fuel cell possibilities: "ZSW-GenStack"(2022), "ZSW-GenStack_Pa_1.61_Pc_1.41", "ZSW-GenStack_Pa_2.01_Pc_1.81",
    #                          "ZSW-GenStack_Pa_2.4_Pc_2.2", "ZSW-GenStack_Pa_2.8_Pc_2.6", "ZSW-GenStack_T_62",
    #                          "ZSW-GenStack_T_76", "ZSW-GenStack_T_84", "EH-31_1.5"(2021), "EH-31_2.0"(2021), "EH-31_2.25"(2021),
    #                          "EH-31_2.5"(2021), "manual_setup".
    # This parameter includes the fuel cell used in the model and the corresponding operating conditions.
    # - GenStack is a fuel cell developed in open source by ZSW (https://zenodo.org/records/14223364).
    # - EH-31 is a fuel cell developed by EH GROUP. 1.5, 2.0, 2.25 and 2.5 corresponds to the different pressure options.
    type_fuel_cell_1 = "EH-31_2.25"
    type_fuel_cell_2 = None
    type_fuel_cell_3 = None
    type_fuel_cell_4 = None
    # Current density possibilities: "step", "polarization", "polarization_for_cali", "EIS".
    type_current = "polarization"
    # Calibration zone: "before_voltage_drop", "full".
    if type_current == "polarization" or type_current == "polarization_for_cali":
        voltage_zone = "before_voltage_drop"
    else:
        voltage_zone = None
    # Auxiliary system possibilities: "forced-convective_cathode_with_anodic_recirculation",
    #                                 "forced-convective_cathode_with_flow-through_anode", "no_auxiliary".
    type_auxiliary = "forced-convective_cathode_with_flow-through_anode"
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

# __________________________________________________________Main________________________________________________________
    """This section is dedicated to ensuring the proper execution of the simulator, considering all the various 
    possibilities including real-time figure updates and simultaneous plotting for different configurations. 
    This should remain unaltered for regular program usage.
    """

    # Retrieving parameters from the settings.py file
    #   Imposed inputs
    (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
     current_density) = current_density_parameters(type_current)
    #   Operating conditions
    T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, y_H2_in_1, pola_current_parameters_1 = \
        operating_inputs_function(copy.deepcopy(pola_current_parameters), type_fuel_cell_1, voltage_zone)
    T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, y_H2_in_2, pola_current_parameters_2 = \
            operating_inputs_function(copy.deepcopy(pola_current_parameters), type_fuel_cell_2, voltage_zone)
    T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3, y_H2_in_3, pola_current_parameters_3 = \
            operating_inputs_function(copy.deepcopy(pola_current_parameters), type_fuel_cell_3, voltage_zone)
    T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4, y_H2_in_4, pola_current_parameters_4 = \
            operating_inputs_function(copy.deepcopy(pola_current_parameters), type_fuel_cell_4, voltage_zone)
    #   Physical parameters
    (Hacl, Hccl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_cl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc,
     Lgc, Vsm_a, Vsm_c, Vem_a, Vem_c, A_T_a, A_T_c, Aact, n_cell, e, Re, i0_d_c_ref, i0_h_c_ref, kappa_co, kappa_c,
     a_slim, b_slim, a_switch, C_scl) = physical_parameters(type_fuel_cell_1)
    #   Computing parameters
    n_gdl, n_mpl, t_purge, rtol, atol, step_current_parameters = computing_parameters(copy.deepcopy(step_current_parameters), Hgdl, Hmpl, Hacl, type_fuel_cell_1)

    # Initialize the operating inputs and parameters dictionaries.
    operating_inputs = {'current_density': current_density, 'T_des': [None, T_des_1, T_des_2, T_des_3, T_des_4],
                        'Pa_des': [None, Pa_des_1, Pa_des_2, Pa_des_3, Pa_des_4],
                        'Pc_des': [None, Pc_des_1, Pc_des_2, Pc_des_3, Pc_des_4],
                        'Sa': [None, Sa_1, Sa_2, Sa_3, Sa_4],
                        'Sc': [None, Sc_1, Sc_2, Sc_3, Sc_4],
                        'Phi_a_des': [None, Phi_a_des_1, Phi_a_des_2, Phi_a_des_3, Phi_a_des_4],
                        'Phi_c_des': [None, Phi_c_des_1, Phi_c_des_2, Phi_c_des_3, Phi_c_des_4],
                        'y_H2_in': [None, y_H2_in_1, y_H2_in_2, y_H2_in_3, y_H2_in_4]}
    current_parameters = {'step_current_parameters': step_current_parameters,
                          'pola_current_parameters': [None, pola_current_parameters_1, pola_current_parameters_2,
                                                      pola_current_parameters_3, pola_current_parameters_4],
                          'pola_current_for_cali_parameters': pola_current_for_cali_parameters,
                          'i_EIS': i_EIS, 'ratio_EIS': ratio_EIS, 't_EIS': t_EIS, 'f_EIS': f_EIS}
    accessible_physical_parameters = {'Aact': Aact, 'n_cell': n_cell, 'Hagc': Hagc, 'Hcgc': Hcgc, 'Wagc': Wagc,
                                      'Wcgc': Wcgc, 'Lgc': Lgc, 'Vsm_a': Vsm_a, 'Vsm_c': Vsm_c, 'Vem_a': Vem_a,
                                      'Vem_c': Vem_c, 'A_T_a': A_T_a, 'A_T_c': A_T_c}
    undetermined_physical_parameters = {'Hgdl': Hgdl, 'Hmpl': Hmpl, 'Hmem': Hmem, 'Hacl': Hacl, 'Hccl': Hccl,
                                        'epsilon_gdl': epsilon_gdl, 'epsilon_cl': epsilon_cl,
                                        'epsilon_mpl': epsilon_mpl, 'epsilon_mc': epsilon_mc, 'epsilon_c': epsilon_c,
                                        'e': e, 'Re': Re, 'i0_d_c_ref': i0_d_c_ref, 'i0_h_c_ref': i0_h_c_ref,
                                        'kappa_co': kappa_co, 'kappa_c': kappa_c, 'a_slim': a_slim, 'b_slim': b_slim,
                                        'a_switch': a_switch, 'C_scl': C_scl}
    computing_parameters = {'n_gdl': n_gdl, 'n_mpl': n_mpl, 't_purge': t_purge, 'rtol': rtol, 'atol': atol,
                            'type_fuel_cell': [None, type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3, type_fuel_cell_4],
                            'type_current': type_current, 'voltage_zone': voltage_zone, 'type_auxiliary': type_auxiliary,
                            'type_control': [None, type_control_1, type_control_2, type_control_3, type_control_4],
                            'type_purge': type_purge, 'type_display': type_display, 'type_plot': type_plot}

    # Check if the type_current is valid and launch the simulation
    if type_current == "step":
        Simulator = launch_AlphaPEM_for_step_current(select_nth_elements(operating_inputs,1),
                                         select_nth_elements(current_parameters, 1),
                                         accessible_physical_parameters, undetermined_physical_parameters,
                                         select_nth_elements(computing_parameters, 1))
    elif type_current == "polarization":
        Simulator = launch_AlphaPEM_for_polarization_current(operating_inputs, current_parameters, accessible_physical_parameters,
                                                 undetermined_physical_parameters, computing_parameters)
    elif type_current == "polarization_for_cali":
        Simulator = launch_AlphaPEM_for_polarization_current_for_calibration(operating_inputs, current_parameters,
                                                                 accessible_physical_parameters,
                                                                 undetermined_physical_parameters, computing_parameters)
    elif type_current == "EIS":
        Simulator = launch_AlphaPEM_for_EIS_current(select_nth_elements(operating_inputs,1),
                                        select_nth_elements(current_parameters, 1),
                                        accessible_physical_parameters, undetermined_physical_parameters,
                                        select_nth_elements(computing_parameters, 1))
    else:
        raise ValueError('You have to specify a type_current which is accepted.')

