# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the current distribution through the GC.
It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
from scipy.optimize import least_squares
from alphapem.core.models.cell_voltage import calculate_cell_voltage
from alphapem.core.modules.cell_voltage_modules import calculate_C_O2_Pt

# Importing constants' value and functions
from alphapem.utils.maths_functions import average


# _________________________________________________Current distribution_________________________________________________
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
