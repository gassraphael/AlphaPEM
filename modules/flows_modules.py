# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import epsilon_cl, theta_c_gdl, theta_c_cl, R
from modules.transitory_functions import nu_l


# _____________________________________________________Flow modules_____________________________________________________

def flows_int_values(sv, operating_inputs, parameters):
    """This functions calculates intermediate values for the flows calculation.

    Parameters
    ----------
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    Pagc : float
        Global pressure in the anode gas channel (Pa).
    Pcgc : float
        Global pressure in the cathode gas channel (Pa).
    s_agdl_agdl : list
        Mean value of the saturated liquid water variable in the anode GDL between two adjacent GDL nodes.
    s_agdl_acl : float
        Mean value of the saturated liquid water variable between the last GDL node and the anode CL.
    s_ccl_cgdl : float
        Mean value of the saturated liquid water variable between the cathode CL and the first GDL node.
    s_cgdl_cgdl : list
        Mean value of the saturated liquid water variable in the cathode GDL between two adjacent GDL nodes.
    epsilon_mean : float
        Mean value of the porosity in the GDL and the CL.
    theta_c_mean : float
        Mean value of the contact angle in the GDL and the CL.
    lambda_acl_mem : float
        Mean value of the dissolved water variable between the anode CL and the membrane.
    lambda_mem_ccl : float
        Mean value of the dissolved water variable between the membrane and the cathode CL.
    Pagc_agdl : float
        Mean value of the pressure between the anode gas channel and the first GDL node.
    Pagdl_agdl : list
        Mean value of the pressure in the anode GDL between two adjacent GDL nodes.
    Pagdl_acl : float
        Mean value of the pressure between the last GDL node and the anode CL.
    Pccl_cgdl : float
        Mean value of the pressure between the cathode CL and the first GDL node.
    Pcgdl_cgdl : list
        Mean value of the pressure in the cathode GDL between two adjacent GDL nodes.
    Pcgdl_cgc : float
        Mean value of the pressure between the last GDL node and the cathode gas channel.
    nu_l : float
        Liquid water kinematic viscosity (m².s-1).
    """

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv['C_v_agc'], sv['C_v_acl'], sv['C_v_ccl'], sv['C_v_cgc']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv['C_H2_agc'], sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_O2_cgc']
    C_N2 = sv['C_N2']
    # Extraction of the operating inputs and the parameters
    Tfc = operating_inputs['Tfc']
    epsilon_gdl, n_gdl = parameters['epsilon_gdl'], parameters['n_gdl']

    # Pressures in the stack
    Pagc = (C_v_agc + C_H2_agc) * R * Tfc
    Pagdl = [(sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}']) * R * Tfc for i in range(1, n_gdl + 1)]
    Pacl = (C_v_acl + C_H2_acl) * R * Tfc
    Pccl = (C_v_ccl + C_O2_ccl + C_N2) * R * Tfc
    Pcgdl = [(sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2) * R * Tfc for i in range(1, n_gdl + 1)]
    Pcgc = (C_v_cgc + C_O2_cgc + C_N2) * R * Tfc

    # Mean values ...
    #       ... of the saturated liquid water variable
    s_agdl_agdl = [None] + [sv[f's_agdl_{i}'] / 2 + sv[f's_agdl_{i + 1}'] / 2 for i in range(1, n_gdl)]
    s_agdl_acl = sv[f's_agdl_{n_gdl}'] / 2 + s_acl / 2
    s_ccl_cgdl = s_ccl / 2 + sv['s_cgdl_1'] / 2
    s_cgdl_cgdl = [None] + [sv[f's_cgdl_{i}'] / 2 + sv[f's_cgdl_{i + 1}'] / 2 for i in range(1, n_gdl)]
    #       ... of the porosity and the contact angle
    epsilon_mean = epsilon_gdl / 2 + epsilon_cl / 2
    theta_c_mean = theta_c_gdl / 2 + theta_c_cl / 2
    #       ... of the dissolved water variable
    lambda_acl_mem = lambda_acl / 2 + lambda_mem / 2
    lambda_mem_ccl = lambda_mem / 2 + lambda_ccl / 2
    #       ... of the pressure
    Pagc_agdl = Pagc / 2 + Pagdl[0] / 2
    Pagdl_agdl = [None] + [Pagdl[i] / 2 + Pagdl[i + 1] / 2 for i in range(0, n_gdl - 1)]
    Pagdl_acl = Pagdl[-1] / 2 + Pacl / 2
    Pccl_cgdl = Pccl / 2 + Pcgdl[0] / 2
    Pcgdl_cgdl = [None] + [Pcgdl[i] / 2 + Pcgdl[i + 1] / 2 for i in range(0, n_gdl - 1)]
    Pcgdl_cgc = Pcgdl[n_gdl - 1] / 2 + Pcgc / 2

    return (Pagc, Pcgc, s_agdl_agdl, s_agdl_acl, s_ccl_cgdl, s_cgdl_cgdl, epsilon_mean, theta_c_mean, lambda_acl_mem,
            lambda_mem_ccl, Pagc_agdl, Pagdl_agdl, Pagdl_acl, Pccl_cgdl, Pcgdl_cgdl, Pcgdl_cgc, nu_l(Tfc))
