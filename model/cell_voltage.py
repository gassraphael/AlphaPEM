# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the cell voltage. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math
import numpy as np
from scipy.optimize import least_squares

# Importing constants' value and functions
from configuration.settings import F, R, E0, Pref_eq
from modules.transitory_functions import average, k_H2, k_O2, sigma_p_eff, R_T_O2_Pt, a_c


# _____________________________________________________Cell voltage_____________________________________________________

def calculate_1D_GC_current_density(i_fc_cell, sv, parameters):
    """This function calculates the local current density distribution in the 1D direction of the GC.

    Parameters
    ----------
    i_fc_cell : float
        Fuel cell current density at time t (A.m-2).
    sv : dict
        Variables calculated by the solver. They correspond to the cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    i_fc : list
        Local current density distribution in the 1D direction of the GC (A.m-2).

    """

    # Extraction of the operating inputs and the parameters
    nb_gc = parameters['nb_gc']
    # Extraction of the variables
    C_O2_ccl = [None] + [sv[i]['C_O2_ccl'] for i in range(1, nb_gc + 1)]

    # Residual function for least_squares solver applied on the local current density
    def residuals(variable_guessed):

        # Recovery of the guessed variable values
        U_cell_guessed = variable_guessed[0]
        i_fc_guessed = [None] + [variable_guessed[i] for i in range(1, nb_gc + 1)]
        C_O2_Pt_guessed = [None] + [variable_guessed[i] for i in range(nb_gc + 1, 2 * nb_gc + 1)]

        # Residuals: difference between
        res_U_cell = [calculate_cell_voltage(i_fc_guessed[i], C_O2_Pt_guessed[i], sv[i], parameters) - U_cell_guessed
                      for i in range(1, nb_gc + 1)]
        res_i_fc = i_fc_cell - average(i_fc_guessed[1:])
        res_C_O2_Pt = [calculate_C_O2_Pt(i_fc_guessed[i], sv[i], parameters)  - C_O2_Pt_guessed[i]
                       for i in range(1, nb_gc + 1)]

        return res_U_cell + [res_i_fc] + res_C_O2_Pt

    # Calculation of the 1D GC current density by solving the system of equations defined by the residuals function using least squares solver
    #       Initial guesses, bounds
    x0 = [calculate_cell_voltage(i_fc_cell, C_O2_ccl[1], sv[1], parameters)] + \
         [i_fc_cell] * nb_gc + \
         [C_O2_ccl[i] for i in range(1, nb_gc + 1)]  # Initial guesses for the least square solver.
    #       Solver call
    sol = least_squares(residuals, x0, method='lm')

    #       Check for convergence
    if not sol.success:
        raise RuntimeError(f"Convergence failed in calculate_1D_GC_current_density: {sol.message}")
    #      Extract the results
    i_fc = [None] + [float(sol.x[i]) for i in range(1, nb_gc + 1)]

    return i_fc


def calculate_C_O2_Pt(i_fc, sv_1D, parameters):
    """This function calculates the oxygen concentration at the platinum surface in the cathode catalyst layer.
    Parameters
    ----------
    i_fc : float
        The current density (A/m²).
    sv_1D : dict
        The dictionary containing the variables calculated by the solver.
    parameters : dict
        The dictionary containing the parameters.

    Returns
    -------
    C_O2_Pt : float
        The oxygen concentration at the platinum surface in the cathode catalyst layer (mol/m³).

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    # Extraction of the variables
    s_ccl, lambda_ccl = sv_1D['s_ccl'], sv_1D['lambda_ccl']
    C_O2_ccl = sv_1D['C_O2_ccl']
    T_ccl = sv_1D['T_ccl']
    # Extraction of the operating inputs and the parameters
    Hccl, K_O2_ad_Pt = parameters['Hccl'], parameters['K_O2_ad_Pt']

    C_O2_Pt = C_O2_ccl - i_fc / (4 * F * Hccl) * R_T_O2_Pt(s_ccl, lambda_ccl, T_ccl, Hccl, K_O2_ad_Pt) / a_c(lambda_ccl, T_ccl, Hccl)

    if C_O2_Pt <= 0:
        raise ValueError("Calculated C_O2_Pt is non-physical (negative or zero). Check input parameters and variables.")

    return C_O2_Pt


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
