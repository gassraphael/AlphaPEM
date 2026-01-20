# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
import cProfile
import copy
from configuration.settings import calculate_current_density_parameters, calculate_operating_inputs, calculate_physical_parameters, \
                                   calculate_computing_parameters
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

def main():
    # Fuel cell possibilities: "ZSW-GenStack"(2022), "ZSW-GenStack_Pa_1.61_Pc_1.41", "ZSW-GenStack_Pa_2.01_Pc_1.81",
    #                          "ZSW-GenStack_Pa_2.4_Pc_2.2", "ZSW-GenStack_Pa_2.8_Pc_2.6", "ZSW-GenStack_T_62",
    #                          "ZSW-GenStack_T_76", "ZSW-GenStack_T_84", "EH-31_1.5"(2021), "EH-31_2.0"(2021),
    #                          "EH-31_2.25"(2021), "EH-31_2.5"(2021), "manual_setup".
    # This parameter includes the fuel cell used in the model and the corresponding operating conditions.
    # - GenStack is a fuel cell developed in open source by ZSW (https://zenodo.org/records/14223364).
    # - EH-31 is a fuel cell developed by EH GROUP. 1.5, 2.0, 2.25 and 2.5 corresponds to the different pressure options.
    type_fuel_cell_1 = "ZSW-GenStack"
    type_fuel_cell_2 = None
    type_fuel_cell_3 = None
    type_fuel_cell_4 = None
    type_fuel_cell_5 = None
    # Current density possibilities: "step", "polarization", "polarization_for_cali", "EIS".
    type_current = "step"
    # Calibration zone : "before_voltage_drop", "full".
    # (only for "polarization" and "polarization_for_cali" currend densities
    voltage_zone = "full"
    # Auxiliary system possibilities: "forced-convective_cathode_with_anodic_recirculation",
    #                                 "forced-convective_cathode_with_flow-through_anode", "no_auxiliary".
    type_auxiliary = "no_auxiliary"
    # Control strategy for the operating inputs: "Phi_des", "no_control".
    type_control_1 = "no_control"
    type_control_2 = "no_control"
    type_control_3 = "no_control"
    type_control_4 = "no_control"
    type_control_5 = "no_control"
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
     current_density) = calculate_current_density_parameters(type_current)
    #   Operating conditions
    T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, y_H2_in_1, pola_current_parameters_1 = \
        calculate_operating_inputs(copy.deepcopy(pola_current_parameters), type_fuel_cell_1, voltage_zone)
    T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, y_H2_in_2, pola_current_parameters_2 = \
        calculate_operating_inputs(copy.deepcopy(pola_current_parameters), type_fuel_cell_2, voltage_zone)
    T_des_3, Pa_des_3, Pc_des_3, Sa_3, Sc_3, Phi_a_des_3, Phi_c_des_3, y_H2_in_3, pola_current_parameters_3 = \
        calculate_operating_inputs(copy.deepcopy(pola_current_parameters), type_fuel_cell_3, voltage_zone)
    T_des_4, Pa_des_4, Pc_des_4, Sa_4, Sc_4, Phi_a_des_4, Phi_c_des_4, y_H2_in_4, pola_current_parameters_4 = \
        calculate_operating_inputs(copy.deepcopy(pola_current_parameters), type_fuel_cell_4, voltage_zone)
    T_des_5, Pa_des_5, Pc_des_5, Sa_5, Sc_5, Phi_a_des_5, Phi_c_des_5, y_H2_in_5, pola_current_parameters_5 = \
        calculate_operating_inputs(copy.deepcopy(pola_current_parameters), type_fuel_cell_5, voltage_zone)
    #   Physical parameters
    (Hacl, Hccl, IC, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc,
     Lgc, nb_channel_in_gc, Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem, Aact, nb_cell, e, Re, i0_d_c_ref,
     i0_h_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl) \
        = calculate_physical_parameters(type_fuel_cell_1)
    #   Computing parameters
    nb_gc, nb_gdl, nb_mpl, t_purge, rtol, atol = calculate_computing_parameters(step_current_parameters)

    # Initialize the operating inputs and parameters dictionaries.
    operating_inputs = {'current_density': current_density,
                        'T_des': [None, T_des_1, T_des_2, T_des_3, T_des_4, T_des_5],
                        'Pa_des': [None, Pa_des_1, Pa_des_2, Pa_des_3, Pa_des_4, Pa_des_5],
                        'Pc_des': [None, Pc_des_1, Pc_des_2, Pc_des_3, Pc_des_4, Pc_des_5],
                        'Sa': [None, Sa_1, Sa_2, Sa_3, Sa_4, Sa_5],
                        'Sc': [None, Sc_1, Sc_2, Sc_3, Sc_4, Sc_5],
                        'Phi_a_des': [None, Phi_a_des_1, Phi_a_des_2, Phi_a_des_3, Phi_a_des_4, Phi_a_des_5],
                        'Phi_c_des': [None, Phi_c_des_1, Phi_c_des_2, Phi_c_des_3, Phi_c_des_4, Phi_c_des_5],
                        'y_H2_in': [None, y_H2_in_1, y_H2_in_2, y_H2_in_3, y_H2_in_4, y_H2_in_5]}
    current_parameters = {'step_current_parameters': step_current_parameters,
                          'pola_current_parameters': [None, pola_current_parameters_1, pola_current_parameters_2,
                                                      pola_current_parameters_3, pola_current_parameters_4,
                                                      pola_current_parameters_5],
                          'pola_current_for_cali_parameters': pola_current_for_cali_parameters,
                          'i_EIS': i_EIS, 'ratio_EIS': ratio_EIS, 't_EIS': t_EIS, 'f_EIS': f_EIS}
    accessible_physical_parameters = {'Aact': Aact, 'nb_cell': nb_cell, 'Hagc': Hagc, 'Hcgc': Hcgc, 'Wagc': Wagc,
                                      'Wcgc': Wcgc, 'Lgc': Lgc, 'nb_channel_in_gc': nb_channel_in_gc, 'Ldist': Ldist,
                                      'Lm': Lm, 'A_T_a': A_T_a, 'A_T_c': A_T_c, 'Vasm': Vasm, 'Vcsm': Vcsm,
                                      'Vaem': Vaem, 'Vcem': Vcem}
    undetermined_physical_parameters = {'Hgdl': Hgdl, 'Hmpl': Hmpl, 'Hmem': Hmem, 'Hacl': Hacl,
                                        'Hccl': Hccl, 'epsilon_gdl': epsilon_gdl, 'epsilon_mpl': epsilon_mpl, 'IC': IC,
                                        'epsilon_c': epsilon_c, 'e': e, 'Re': Re, 'i0_d_c_ref': i0_d_c_ref,
                                        'i0_h_c_ref': i0_h_c_ref, 'kappa_co': kappa_co, 'kappa_c': kappa_c,
                                        'a_slim': a_slim, 'b_slim': b_slim, 'a_switch': a_switch, 'C_scl': C_scl}
    computing_parameters = {'nb_gc': nb_gc, 'nb_gdl': nb_gdl, 'nb_mpl': nb_mpl, 't_purge': t_purge,
                            'rtol': rtol, 'atol': atol,
                            'type_fuel_cell': [None, type_fuel_cell_1, type_fuel_cell_2, type_fuel_cell_3,
                                               type_fuel_cell_4, type_fuel_cell_5],
                            'type_current': type_current, 'voltage_zone': voltage_zone,
                            'type_auxiliary': type_auxiliary,
                            'type_control': [None, type_control_1, type_control_2, type_control_3, type_control_4,
                                             type_control_5],
                            'type_purge': type_purge, 'type_display': type_display, 'type_plot': type_plot}

    # Check if the type_current is valid and launch the simulation
    if type_current == "step":
        Simulator = launch_AlphaPEM_for_step_current(select_nth_elements(operating_inputs, 1),
                                                     select_nth_elements(current_parameters, 1),
                                                     accessible_physical_parameters, undetermined_physical_parameters,
                                                     select_nth_elements(computing_parameters, 1))
    elif type_current == "polarization":
        Simulator = launch_AlphaPEM_for_polarization_current(operating_inputs, current_parameters,
                                                             accessible_physical_parameters,
                                                             undetermined_physical_parameters, computing_parameters)
    elif type_current == "polarization_for_cali":
        Simulator = launch_AlphaPEM_for_polarization_current_for_calibration(operating_inputs, current_parameters,
                                                                             accessible_physical_parameters,
                                                                             undetermined_physical_parameters,
                                                                             computing_parameters)
    elif type_current == "EIS":
        Simulator = launch_AlphaPEM_for_EIS_current(select_nth_elements(operating_inputs, 1),
                                                    select_nth_elements(current_parameters, 1),
                                                    accessible_physical_parameters, undetermined_physical_parameters,
                                                    select_nth_elements(computing_parameters, 1))
    else:
        raise ValueError('You have to specify a type_current which is accepted.')

    return Simulator

if __name__ == '__main__':
    main()
    # cProfile.run('main()', sort='ncalls')  # Profiling sorted by number of calls

