# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the cell voltage. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing constants' value and functions
from configuration.settings import F, R, E0, Pref
from modules.transitory_functions import average, k_H2, k_O2, sigma_p_eff


# _____________________________________________________Cell voltage_____________________________________________________

def calculate_eta_c_intermediate_values(solver_variables, operating_inputs, parameters):
    """This function calculates the intermediate values needed for the calculation of the cathode overpotential dynamic
    evolution.

    Parameters
    ----------
    solver_variables : dict
        The dictionary containing the variables calculated by the solver.
    operating_inputs : dict
        The dictionary containing the operating inputs.
    parameters : dict
        The dictionary containing the parameters.

    Returns
    -------
    dict
        The dictionary containing the crossover current density i_n at time t, and the liquid water induced voltage drop
        function f_drop at time t.
    """

    # Extraction of the variables
    s_ccl, lambda_mem = solver_variables['s_ccl'], solver_variables['lambda_mem']
    C_H2_acl, C_O2_ccl = solver_variables['C_H2_acl'], solver_variables['C_O2_ccl']
    T_acl, T_mem, T_ccl = solver_variables['T_acl'], solver_variables['T_mem'], solver_variables['T_ccl']
    # Extraction of the operating inputs and the parameters
    Pc_des = operating_inputs['Pc_des']
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    kappa_co, kappa_c = parameters['kappa_co'], parameters['kappa_c']
    a_slim, b_slim, a_switch = parameters['a_slim'], parameters['b_slim'], parameters['a_switch']

    # The crossover current density i_n
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
    i_n = i_H2 + i_O2

    # The liquid water induced voltage drop function f_drop
    slim = a_slim * (Pc_des / 1e5) + b_slim
    s_switch = a_switch * slim
    f_drop = 0.5 * (1.0 - math.tanh((4 * s_ccl - 2 * slim - 2 * s_switch) / (slim - s_switch)))

    return {'i_n': i_n, 'f_drop': f_drop}


def calculate_cell_voltage(variables, operating_inputs, parameters):
    """This function calculates the cell voltage at each time step.

    Parameters
    ----------
    variables : dict
        The dictionary containing the variables calculated by the solver.
    operating_inputs : dict
        The dictionary containing the operating inputs.
    parameters : dict
        The dictionary containing the parameters.

    Returns
    -------
    Ucell_t : list
        The cell voltage at each time step.
    """

    # Extraction of the variables
    t, lambda_mem_t, lambda_ccl_t = variables['t'], variables['lambda_mem'], variables['lambda_ccl']
    C_H2_acl_t, C_O2_ccl_t, eta_c_t = variables['C_H2_acl'], variables['C_O2_ccl'], variables['eta_c']
    T_acl_t, T_mem_t, T_ccl_t = variables['T_acl'], variables['T_mem'], variables['T_ccl']
    # Extraction of the operating inputs and the parameters
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    epsilon_mc, Re, kappa_co = parameters['epsilon_mc'], parameters['Re'], parameters['kappa_co']

    # Initialisation
    n = len(t)
    Ucell_t = [0] * n

    # Loop for having Ucell_t at each time step
    for i in range(n):

        # Recovery of the already calculated variable values at each time step
        lambda_mem, lambda_ccl = lambda_mem_t[i], lambda_ccl_t[i]
        C_H2_acl, C_O2_ccl = C_H2_acl_t[i], C_O2_ccl_t[i]
        T_acl, T_mem, T_ccl = T_acl_t[i], T_mem_t[i], T_ccl_t[i]
        eta_c = eta_c_t[i]

        # Current density value at this time step
        i_fc = operating_inputs['current_density'](t[i], parameters)

        # The equilibrium potential
        Ueq = E0 - 8.5e-4 * (T_ccl - 298.15) + R * T_ccl / (2 * F) * (math.log(R * T_acl * C_H2_acl / Pref) +
                                                                  0.5 * math.log(R * T_ccl * C_O2_ccl / Pref))

        # The crossover current density
        T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                              weights=[Hacl/(Hacl + Hmem + Hccl), Hmem/(Hacl + Hmem + Hccl), Hccl/(Hacl + Hmem + Hccl)])
        i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
        i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
        i_n = i_H2 + i_O2

        # The proton resistance
        #       The proton resistance at the membrane : Rmem
        Rmem = Hmem / sigma_p_eff('mem', lambda_mem, T_mem)
        #       The proton resistance at the cathode catalyst layer : Rccl
        Rccl = Hccl / sigma_p_eff('ccl', lambda_ccl, T_ccl, epsilon_mc)
        #       The total proton resistance
        Rp = Rmem + Rccl  # its value is around [4-7]e-6 ohm.m².

        # The cell voltage
        Ucell_t[i] = Ueq - eta_c - (i_fc + i_n) * (Rp + Re)
    return Ucell_t
