# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the cell voltage. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import F, R, E0, Pref
from modules.transitory_functions import k_H2, k_O2


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
    # Extraction of the operating inputs and the parameters
    Tfc, Pc_des = operating_inputs['Tfc'], operating_inputs['Pc_des']
    Hmem = parameters['Hmem']
    i0_c_ref, kappa_co, kappa_c = parameters['i0_c_ref'], parameters['kappa_co'], parameters['kappa_c']
    a_slim, b_slim, a_switch = parameters['a_slim'], parameters['b_slim'], parameters['a_switch']

    # The crossover current density i_n
    i_H2 = 2 * F * R * Tfc / Hmem * C_H2_acl * k_H2(lambda_mem, Tfc, kappa_co)
    i_O2 = 4 * F * R * Tfc / Hmem * C_O2_ccl * k_O2(lambda_mem, Tfc, kappa_co)
    i_n = i_H2 + i_O2

    # The liquid water induced voltage drop function f_drop
    slim = a_slim * (Pc_des / 1e5) + b_slim
    s_switch = a_switch * slim
    f_drop = 0.5 * (1.0 - np.tanh((4 * s_ccl - 2 * slim - 2 * s_switch) / (slim - s_switch)))

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
    # Extraction of the operating inputs and the parameters
    Tfc = operating_inputs['Tfc']
    Hmem, Hcl, epsilon_mc, tau = parameters['Hmem'], parameters['Hcl'], parameters['epsilon_mc'], parameters['tau']
    Re, kappa_co = parameters['Re'], parameters['kappa_co']

    # Initialisation
    n = len(t)
    Ucell_t = [0] * n

    # Loop for having Ucell_t at each time step
    for i in range(n):

        # Recovery of the already calculated variable values at each time step
        lambda_mem, lambda_ccl = lambda_mem_t[i], lambda_ccl_t[i]
        C_H2_acl, C_O2_ccl = C_H2_acl_t[i], C_O2_ccl_t[i]
        eta_c = eta_c_t[i]

        # Current density value at this time step
        i_fc = operating_inputs['current_density'](t[i], parameters)

        # The equilibrium potential
        Ueq = E0 - 8.5e-4 * (Tfc - 298.15) + R * Tfc / (2 * F) * (np.log(R * Tfc * C_H2_acl / Pref) +
                                                                  0.5 * np.log(R * Tfc * C_O2_ccl / Pref))

        # The crossover current density
        i_H2 = 2 * F * R * Tfc / Hmem * C_H2_acl * k_H2(lambda_mem, Tfc, kappa_co)
        i_O2 = 4 * F * R * Tfc / Hmem * C_O2_ccl * k_O2(lambda_mem, Tfc, kappa_co)
        i_n = i_H2 + i_O2

        # The proton resistance
        #       The proton resistance at the membrane : Rmem
        if lambda_mem >= 1:
            Rmem = Hmem / ((0.5139 * lambda_mem - 0.326) * np.exp(1268 * (1 / 303.15 - 1 / Tfc)))
        else:
            Rmem = Hmem / (0.1879 * np.exp(1268 * (1 / 303.15 - 1 / Tfc)))
        #       The proton resistance at the cathode catalyst layer : Rccl
        if lambda_ccl >= 1:
            Rccl = 1 / 3 * 1 / (epsilon_mc / tau) * \
                   Hcl / ((0.5139 * lambda_ccl - 0.326) * np.exp(1268 * (1 / 303.15 - 1 / Tfc)))
        else:
            Rccl = 1 / 3 * 1 / (epsilon_mc / tau) * \
                   Hcl / (0.1879 * np.exp(1268 * (1 / 303.15 - 1 / Tfc)))
        #       The total proton resistance
        Rp = Rmem + Rccl  # its value is around [4-7]e-6 ohm.m².

        # The cell voltage
        Ucell_t[i] = Ueq - eta_c - (i_fc + i_n) * (Rp + Re)
    return Ucell_t
