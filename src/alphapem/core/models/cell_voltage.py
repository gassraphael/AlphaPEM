# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the cell voltage. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing constants' value and functions
from alphapem.utils.physics_constants import F, R, E0, Pref_eq
from alphapem.utils.maths_functions import average
from alphapem.core.modules.heat_modules import sigma_p_eff
from alphapem.core.modules.flows_1D_MEA_modules import k_H2, k_O2


# _____________________________________________________Cell voltage_____________________________________________________

def calculate_cell_voltage(i_fc, C_O2_Pt, sv, parameters):
    """This function calculates the cell voltage in volt.

    Parameters
    ----------
    i_fc : float
        The current density (A/m²).
    C_O2_Pt : float
        The oxygen concentration at the platinum surface in the cathode catalyst layer (mol/m³).
    sv : dict
        The dictionary containing the variables calculated by the solver.
    parameters : dict
        The dictionary containing the parameters.

    Returns
    -------
    Ucell_t : list
        The cell voltage in volt.
    """

    # Extraction of the variables
    s_ccl, lambda_mem, lambda_ccl = sv['s_ccl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_acl, C_O2_ccl = sv['C_H2_acl'], sv['C_O2_ccl']
    eta_c = sv['eta_c']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']
    # Extraction of the operating inputs and the parameters
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    Re, kappa_co = parameters['Re'], parameters['kappa_co']

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
    Rccl = Hccl / sigma_p_eff('ccl', lambda_ccl, T_ccl, Hcl=Hccl)
    #       The total proton resistance
    Rp = Rmem + Rccl  # its value is around [4-7]e-6 ohm.m².

    # The cell voltage
    Ucell = Ueq - eta_c - (i_fc + i_n) * (Rp + Re)

    return Ucell
