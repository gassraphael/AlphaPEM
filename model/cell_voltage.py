# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the cell voltage. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing constants' value and functions
from configuration.settings import F, R, E0, Pref_eq
from modules.transitory_functions import average, k_H2, k_O2, sigma_p_eff, R_T_O2_Pt, a_c


# _____________________________________________________Cell voltage_____________________________________________________

def calculate_C_O2_Pt(i_fc, s_ccl, lambda_ccl, C_O2_ccl, T_ccl, Hccl, IC, **kwargs):
    """This function calculates the oxygen concentration at the platinum surface in the cathode catalyst layer.
    Parameters
    ----------
    i_fc : float
        The current density (A/m²).
    s_ccl : float
        The saturation in the cathode catalyst layer (-).
    lambda_ccl : float
        The membrane water content in the cathode catalyst layer (-).
    C_O2_ccl : float
        The oxygen concentration in the cathode catalyst layer (mol/m³).
    T_ccl : float
        The temperature in the cathode catalyst layer (K).
    Hccl : float
        The thickness of the cathode catalyst layer (m).
    IC : float
        The ionomer content in the cathode catalyst layer (-).

    Returns
    -------
    C_O2_Pt : float
        The oxygen concentration at the platinum surface in the cathode catalyst layer (mol/m³).

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return C_O2_ccl - i_fc / (4 * F * Hccl) * R_T_O2_Pt(s_ccl, lambda_ccl, T_ccl, Hccl, IC) / a_c(lambda_ccl, T_ccl, Hccl, IC)


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
    t = variables['t']
    s_ccl_t, lambda_mem_t, lambda_ccl_t = variables['s_ccl'], variables['lambda_mem'], variables['lambda_ccl']
    C_H2_acl_t, C_O2_ccl_t, eta_c_t = variables['C_H2_acl'], variables['C_O2_ccl'], variables['eta_c']
    T_acl_t, T_mem_t, T_ccl_t = variables['T_acl'], variables['T_mem'], variables['T_ccl']
    # Extraction of the operating inputs and the parameters
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    IC, Re, kappa_co = parameters['IC'], parameters['Re'], parameters['kappa_co']

    # Initialisation
    n = len(t)
    Ucell_t = [0] * n

    # Loop for having Ucell_t at each time step
    for i in range(n):

        # Recovery of the already calculated variable values at each time step
        s_ccl, lambda_mem, lambda_ccl = s_ccl_t[i], lambda_mem_t[i], lambda_ccl_t[i]
        C_H2_acl, C_O2_ccl = C_H2_acl_t[i], C_O2_ccl_t[i]
        T_acl, T_mem, T_ccl = T_acl_t[i], T_mem_t[i], T_ccl_t[i]
        eta_c = eta_c_t[i]

        # Current density value at this time step
        i_fc = operating_inputs['current_density'](t[i], parameters)

        C_O2_Pt = calculate_C_O2_Pt(i_fc, s_ccl, lambda_ccl, C_O2_ccl, T_ccl, Hccl, IC)

        # The equilibrium potential
        Ueq = E0 - 8.5e-4 * (T_ccl - 298.15) + R * T_ccl / (2 * F) * (math.log(R * T_acl * C_H2_acl / Pref_eq) +
                                                                      0.5 * math.log(R * T_ccl * C_O2_Pt / Pref_eq))

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
        Rccl = Hccl / sigma_p_eff('ccl', lambda_ccl, T_ccl, Hcl=Hccl, IC=IC)
        #       The total proton resistance
        Rp = Rmem + Rccl  # its value is around [4-7]e-6 ohm.m².

        # The cell voltage
        Ucell_t[i] = Ueq - eta_c - (i_fc + i_n) * (Rp + Re)

        # print('i_fc = ', i_fc, 'Ucell_t = ', Ucell_t[i], 'eta_c = ', eta_c, 'C_O2_ccl = ', C_O2_ccl, 'C_O2_Pt = ', C_O2_Pt)
    return Ucell_t
