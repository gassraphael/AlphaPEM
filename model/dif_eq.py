# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model. It is a master file which converts
the 1D+1D+1D model to several 1D models in order to ease the coding.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

from configuration.settings import E0
from model.flows_1D_GC_manifold import calculate_flows_1D_GC_manifold
# Importing constants' value and functions
from model.velocity import calculate_velocity_evolution
from model.flows_1D_MEA import calculate_flows_1D_MEA
from model.cell_voltage import calculate_1D_GC_current_density, calculate_C_O2_Pt
from model.heat_transfer import calculate_heat_transfers
from modules.dif_eq_modules import calculate_dif_eq_int_values
from model.dif_eq_1D_MEA import (calculate_dyn_dissoved_water_evolution_inside_MEA, calculate_dyn_liquid_water_evolution_inside_MEA,
                                 calculate_dyn_vapor_evolution_inside_MEA, calculate_dyn_H2_O2_N2_evolution_inside_MEA,
                                 calculate_dyn_voltage_evolution, calculate_dyn_temperature_evolution_inside_MEA)
from model.dif_eq_1D_GC_manifold import (calculate_dyn_gas_evolution_inside_gas_channel,
                                         calculate_dyn_temperature_evolution_inside_gas_channel,
                                         calculate_dyn_manifold_pressure_and_humidity_evolution,
                                         calculate_dyn_liq_evolution_inside_gas_channel)
from model.dif_eq_auxiliaries import (calculate_dyn_air_compressor_evolution, calculate_dyn_humidifier_evolution,
                                      calculate_dyn_throttle_area_controler)


# ______________________Objective function to solve. It gives the system of differential equations______________________

def dydt(t, y, operating_inputs, parameters, solver_variable_names):
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

    Returns
    -------
    dydt : list
        List containing the derivative of the solver variables.
    """

    # Creation of the solver_variables and dif_eq dictionaries. They are intermediates to simplify the code.
    nb_gc = parameters['nb_gc']
    solver_variables_1D_cell = [None] * (nb_gc + 1)  # First index is None to be consistent with GC numbering starting at 1
    solver_variable_1D_manifold = {} # Dictionary corresponding to the variable inside the gas channel and manifold 1D line
    solver_variable_auxiliary = {} # Dictionary corresponding to the auxiliary variables
    dif_eq_1D_cell = [None] * (nb_gc + 1) # First index is None to be consistent with GC numbering starting at 1
    dif_eq_1D_manifold = {} # Dictionary corresponding to the differential equations inside the gas channel 1D line
    dif_eq_auxiliary = {} # Dictionary corresponding to the differential equations of the auxiliary variables
    for i in range(nb_gc): # Each dictionary in this loop corresponds to one gas channel node
        solver_variables_1D_cell[i + 1] = {} # Dictionary corresponding to the variable inside the GC and the MEA at the GC node i
        dif_eq_1D_cell[i + 1] = {} # Dictionary corresponding to the differential equations inside the GC and the MEA at the GC node i
        for index, variable in enumerate(solver_variable_names[0]): # Concern only the MEA and GC variables
            solver_variables_1D_cell[i + 1][variable] = y[index + i * len(solver_variable_names[0])]
            dif_eq_1D_cell[i + 1]['d' + variable + ' / dt'] = None
    if parameters['type_auxiliary'] == "forced-convective_cathode_with_flow-through_anode" or \
            parameters['type_auxiliary'] == "forced-convective_cathode_with_anodic_recirculation":
        for index, variable in enumerate(solver_variable_names[1]): # Concern only the manifold variables
            solver_variable_1D_manifold[f'{variable}'] = y[index + nb_gc * len(solver_variable_names[0])]
            dif_eq_1D_manifold['d' + f'{variable}' + ' / dt'] = None
        for index, variable in enumerate(solver_variable_names[2]):  # Concern only the auxiliary variables
            solver_variable_auxiliary[f'{variable}'] = y[index + nb_gc * len(solver_variable_names[0]) + len(solver_variable_names[1])]
            dif_eq_auxiliary['d' + f'{variable}' + ' / dt'] = None

    # Conditions to pursue the calculations
    if solver_variables_1D_cell[1]['eta_c'] > E0:                                                                        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        raise ValueError("The cathode overpotential is higher than the open circuit voltage at time t = " + str(t) +
                         " s. It means that the voltage is negative, which is not possible.")

    # Intermediate values

    dif_eq_int_values = [None] + [calculate_dif_eq_int_values(t, solver_variables_1D_cell[i], operating_inputs,
                                                              parameters) for i in range(1, nb_gc + 1)]

    # Calculate the local current density at each node of the GC.
    i_fc_cell = operating_inputs['current_density'](t, parameters)
    i_fc = calculate_1D_GC_current_density(i_fc_cell, solver_variables_1D_cell, parameters)

    # Calculation of the oxygen concentration at the platinum surface in the cathode catalyst layer
    C_O2_Pt = [None] + [calculate_C_O2_Pt(i_fc[i], solver_variables_1D_cell[i], parameters)
                        for i in range(1, nb_gc + 1)]

    # Calculation of the velocities inside the GC and the manifolds
    v_a, v_c, Pa_in, Pc_in = calculate_velocity_evolution(solver_variables_1D_cell, i_fc_cell, operating_inputs, parameters)

    # Calculation of the flows. These elements are lists of dictionaries. Each index corresponds to one GC node.
    # Each element of the dictionary corresponds to one flow inside this GC node.
    flows_1D_MEA = [None] + [calculate_flows_1D_MEA(solver_variables_1D_cell[i], i_fc[i], v_a[i], v_c[i],
                                                    operating_inputs, parameters)
                             for i in range(1, nb_gc + 1)]
    flows_1D_GC_manifold = calculate_flows_1D_GC_manifold(solver_variables_1D_cell, solver_variable_1D_manifold,
                                                          solver_variable_auxiliary, i_fc_cell, v_a, v_c, Pa_in, Pc_in,
                                                          operating_inputs, parameters)
    for dif_eq, flow_names in {'agc_agdl': ('Jl', 'Jv', 'J_H2'), 'cgdl_cgc': ('Jl', 'Jv', 'J_O2')}.items():
        for flow_name in flow_names:
            flows_1D_GC_manifold[flow_name][dif_eq] = [None]
            for i in range(1, nb_gc + 1):
                flows_1D_GC_manifold[flow_name][dif_eq].append(flows_1D_MEA[i][flow_name][dif_eq])
    heat_flows_global = [None] + [calculate_heat_transfers(solver_variables_1D_cell[i], i_fc[i], operating_inputs,
                                                           parameters, **flows_1D_MEA[i])
                                  for i in range(1, nb_gc + 1)]

    # Calculation of the dynamic evolutions
    for i in range(1, nb_gc + 1):
    #       Inside the MEA
        calculate_dyn_dissoved_water_evolution_inside_MEA(dif_eq_1D_cell[i], solver_variables_1D_cell[i], **parameters,
                                                          **flows_1D_MEA[i])
        calculate_dyn_liquid_water_evolution_inside_MEA(dif_eq_1D_cell[i], solver_variables_1D_cell[i], **parameters,
                                                        **flows_1D_MEA[i])
        calculate_dyn_vapor_evolution_inside_MEA(dif_eq_1D_cell[i], solver_variables_1D_cell[i], **parameters,
                                                 **flows_1D_MEA[i])
        calculate_dyn_H2_O2_N2_evolution_inside_MEA(dif_eq_1D_cell[i], solver_variables_1D_cell[i], **parameters,
                                                    **flows_1D_MEA[i])
        calculate_dyn_voltage_evolution(dif_eq_1D_cell[i], i_fc[i], C_O2_Pt[i], **solver_variables_1D_cell[i],
                                        **operating_inputs, **parameters, **dif_eq_int_values[i])                       # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        calculate_dyn_temperature_evolution_inside_MEA(dif_eq_1D_cell[i], **parameters, **dif_eq_int_values[i],
                                                       **heat_flows_global[i])
    #       Inside the gaz channels and the manifolds
    calculate_dyn_gas_evolution_inside_gas_channel(dif_eq_1D_cell, solver_variables_1D_cell, **parameters,
                                                   **flows_1D_GC_manifold)
    calculate_dyn_liq_evolution_inside_gas_channel(dif_eq_1D_cell, **operating_inputs, **parameters,
                                                   **flows_1D_GC_manifold)
    calculate_dyn_temperature_evolution_inside_gas_channel(dif_eq_1D_cell, **parameters)
    if parameters['type_auxiliary'] != "no_auxiliary":
        calculate_dyn_manifold_pressure_and_humidity_evolution(dif_eq_1D_manifold, **operating_inputs, **parameters,
                                                               **flows_1D_GC_manifold)
    #       Inside the auxiliaries
    if parameters['type_auxiliary'] != "no_auxiliary":
        calculate_dyn_air_compressor_evolution(dif_eq_auxiliary, **solver_variable_auxiliary, **parameters)
        calculate_dyn_humidifier_evolution(dif_eq_auxiliary, **solver_variable_auxiliary, **parameters, **dif_eq_int_values)
        calculate_dyn_throttle_area_controler(dif_eq_auxiliary, solver_variable_auxiliary, **operating_inputs, **parameters,
                                              **flows_1D_GC_manifold)

    # all the dif_eq dictionaries are converted because the solver requires an ordered list to work
    dif_eq_global = [dif_eq_1D_cell[i]['d' + key + ' / dt'] for i in range(1, nb_gc + 1)
                     for key in solver_variable_names[0]]
    if parameters['type_auxiliary'] == "forced-convective_cathode_with_flow-through_anode" or \
            parameters['type_auxiliary'] == "forced-convective_cathode_with_anodic_recirculation":
        dif_eq_global += [dif_eq_1D_manifold['d' + key + ' / dt'] for key in solver_variable_names[1]]
        dif_eq_global += [dif_eq_auxiliary['d' + key + ' / dt'] for key in solver_variable_names[2]]
    return dif_eq_global