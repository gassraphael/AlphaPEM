# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from model.flows import calculate_flows
from model.cell_voltage import calculate_eta_c_intermediate_values
from model.heat_transfer import calculate_heat_transfers
from model.control import control_operating_conditions
from modules.dif_eq_modules import calculate_dif_eq_int_values
from model.dif_eq_MEA import (calculate_dyn_dissoved_water_evolution_inside_MEA, calculate_dyn_liquid_water_evolution_inside_MEA,
                              calculate_dyn_vapor_evolution_inside_MEA, calculate_dyn_H2_O2_N2_evolution_inside_MEA,
                              calculate_dyn_voltage_evolution, calculate_dyn_temperature_evolution_inside_MEA)
from model.dif_eq_GC_manifold import (calculate_dyn_gas_evolution_inside_gas_channel,
                                      calculate_dyn_temperature_evolution_inside_gas_channel,
                                      calculate_dyn_manifold_pressure_and_humidity_evolution)
from model.dif_eq_auxiliaries import (calculate_dyn_air_compressor_evolution, calculate_dyn_humidifier_evolution,
                                      calculate_dyn_throttle_area_controler)

# ______________________Objective function to solve. It gives the system of differential equations______________________

def dydt(t, y, operating_inputs, parameters, solver_variable_names, control_variables):
    """This function gives the system of differential equations to solve.

    Parameters
    ----------
    t : float
        Time (s).
    y : numpy.ndarray
        Numpy list of the solver variables.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    solver_variable_names : list
        Names of the solver variables.
    control_variables : dict
        Variables controlled by the user.

    Returns
    -------
    dydt : list
        List containing the derivative of the solver variables.
    """

    # Creation of the dif_eq dictionary. It is an intermediate calculation to simplify the writing of the code.
    dif_eq = {('d' + key + ' / dt'): None for key in solver_variable_names}
    # Creation of the solver_variables dict. It is an intermediate calculation to simplify the writing of the code.
    solver_variables = {}
    for index, key in enumerate(solver_variable_names):
        solver_variables[key] = y[index]

    # Modifications of the operating conditions in real time, if required.
    if parameters["type_control"] != "no_control":
        control_operating_conditions(t, solver_variables, operating_inputs, parameters, control_variables)

    if parameters['type_auxiliary'] != "no_auxiliary":
        raise ValueError("Currently, the simulation with auxiliaries are in reconstruction. Please set 'type_auxiliary' to 'no_auxiliary'.")

    # Intermediate values
    i_fc = operating_inputs['current_density'](t, parameters)
    dif_eq_int_values = calculate_dif_eq_int_values(t, solver_variables, control_variables, i_fc, operating_inputs,
                                                    parameters)
    eta_c_intermediate_values = calculate_eta_c_intermediate_values(solver_variables, operating_inputs, parameters)

    # Calculation of the flows
    matter_flows_dico = calculate_flows(t, solver_variables, control_variables, i_fc, operating_inputs, parameters)
    heat_flows_dico = calculate_heat_transfers(solver_variables, i_fc, operating_inputs, parameters, **matter_flows_dico)

    # Calculation of the dynamic evolutions
    #       Inside the MEA
    calculate_dyn_dissoved_water_evolution_inside_MEA(dif_eq, **parameters, **matter_flows_dico)
    calculate_dyn_liquid_water_evolution_inside_MEA(dif_eq, solver_variables, **parameters, **matter_flows_dico)
    calculate_dyn_vapor_evolution_inside_MEA(dif_eq, solver_variables, **parameters, **matter_flows_dico)
    calculate_dyn_H2_O2_N2_evolution_inside_MEA(dif_eq, solver_variables, **parameters, **matter_flows_dico)
    calculate_dyn_voltage_evolution(dif_eq, i_fc, **solver_variables, **operating_inputs, **parameters,
                                    **eta_c_intermediate_values)
    calculate_dyn_temperature_evolution_inside_MEA(dif_eq, **parameters, **dif_eq_int_values, **heat_flows_dico)
    #       Inside the gaz channels and the manifolds
    calculate_dyn_gas_evolution_inside_gas_channel(dif_eq, **parameters, **matter_flows_dico)
    calculate_dyn_temperature_evolution_inside_gas_channel(dif_eq, **parameters)
    if parameters['type_auxiliary'] != "no_auxiliary":
        calculate_dyn_manifold_pressure_and_humidity_evolution(dif_eq, **operating_inputs, **parameters, **matter_flows_dico)
    #       Inside the auxiliaries
    if parameters['type_auxiliary'] != "no_auxiliary":
        calculate_dyn_air_compressor_evolution(dif_eq, **solver_variables, **parameters)
        calculate_dyn_humidifier_evolution(dif_eq, **solver_variables, **parameters, **dif_eq_int_values)
        calculate_dyn_throttle_area_controler(dif_eq, solver_variables, **operating_inputs, **parameters,
                                              **matter_flows_dico)

    # dif_eq is converted to dydt because the solver requires an ordered list to work
    dydt = np.array([dif_eq['d' + key + ' / dt'] for key in solver_variable_names])
    return dydt