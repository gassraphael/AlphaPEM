# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the heat transfer calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
from scipy.stats import hmean

# Importing constants' value and functions
from configuration.settings import epsilon_cl
from modules.transitory_functions import k_th_eff


# _____________________________________________________Flow modules_____________________________________________________

def heat_transfer_int_values(sv, parameters):
    """This functions calculates intermediate values for the heat calculation.

    Parameters
    ----------
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    tuple
        A tuple containing the average effective thermal diffusivity values between different layers:
        - k_th_eff_agc_agdl : Effective thermal diffusivity between the AGC and the first GDL layer.
        - k_th_eff_agdl_agdl : List of effective thermal diffusivities between adjacent GDL layers on the anode side.
        - k_th_eff_agdl_acl : Effective thermal diffusivity between the last GDL layer and the anode catalyst layer.
        - k_th_eff_acl_mem : Effective thermal diffusivity between the anode catalyst layer and the membrane.
        - k_th_eff_mem_ccl : Effective thermal diffusivity between the membrane and the cathode catalyst layer.
        - k_th_eff_ccl_cgdl : Effective thermal diffusivity between the cathode catalyst layer and the first GDL layer.
        - k_th_eff_cgdl_cgdl : List of effective thermal diffusivities between adjacent GDL layers on the cathode side.
        - k_th_eff_cgdl_cgc : Effective thermal diffusivity between the last GDL layer and the CGC.
    """

    # Extraction of the variables
    C_v_acl, C_v_ccl, s_acl, s_ccl = sv['C_v_acl'], sv['C_v_ccl'], sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_acl, C_O2_ccl, C_N2 = sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_N2']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']

    # Extraction of the operating inputs and the parameters
    Hgdl, Hcl, Hmem = parameters['Hgdl'], parameters['Hcl'], parameters['Hmem']
    epsilon_mc, tau, epsilon_gdl = parameters['epsilon_mc'], parameters['tau'], parameters['epsilon_gdl']
    n_gdl = parameters['n_gdl']

    # Weighted harmonic means of the effective thermal diffusivity
    k_th_eff_agc_agdl = k_th_eff('agdl', sv[f'T_agdl_{1}'], C_v=sv[f'C_v_agdl_{1}'], s=sv[f's_agdl_{1}'],
                                 C_H2=sv[f'C_H2_agdl_{1}'], epsilon=epsilon_gdl)

    k_th_eff_agdl_agdl = [None] + [hmean([k_th_eff('agdl', sv[f'T_agdl_{i}'], C_v=sv[f'C_v_agdl_{i}'],
                                                   s=sv[f's_agdl_{i}'], C_H2=sv[f'C_H2_agdl_{i}'], epsilon=epsilon_gdl),
                                          k_th_eff('agdl', sv[f'T_agdl_{i + 1}'], C_v=sv[f'C_v_agdl_{i + 1}'],
                                                   s=sv[f's_agdl_{i+1}'], C_H2=sv[f'C_H2_agdl_{i+1}'], epsilon=epsilon_gdl)])
                                   for i in range(1, n_gdl)]

    k_th_eff_agdl_acl = hmean([k_th_eff('agdl', sv[f'T_agdl_{n_gdl}'], C_v=sv[f'C_v_agdl_{n_gdl}'],
                                        s=sv[f's_agdl_{n_gdl}'], C_H2=sv[f'C_H2_agdl_{n_gdl}'], epsilon=epsilon_gdl),
                               k_th_eff('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl,
                                        C_H2=C_H2_acl, epsilon=epsilon_cl, epsilon_mc=epsilon_mc)],
                              weights=[Hgdl / 2, Hcl / 2])

    k_th_eff_acl_mem = hmean([k_th_eff('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl,
                                       C_H2=C_H2_acl, epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
                              k_th_eff('mem', T_mem, lambdaa=lambda_mem)],
                             weights=[Hcl / 2, Hmem / 2])

    k_th_eff_mem_ccl = hmean([k_th_eff('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl, C_O2=C_O2_ccl,
                                       C_N2=C_N2, epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
                              k_th_eff('mem', T_mem, lambdaa=lambda_mem)],
                             weights=[Hcl / 2, Hmem / 2])

    k_th_eff_ccl_cgdl = hmean([k_th_eff('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl, C_O2=C_O2_ccl,
                                        C_N2=C_N2, epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
                               k_th_eff('cgdl', sv['T_cgdl_1'], C_v=sv['C_v_cgdl_1'], s=sv['s_cgdl_1'],
                                        C_O2=sv[f'C_O2_cgdl_1'], C_N2=C_N2, epsilon=epsilon_gdl)],
                              weights=[Hcl / 2, Hgdl / 2])

    k_th_eff_cgdl_cgdl = [None] + [hmean([k_th_eff('cgdl', sv[f'T_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                                                   s=sv[f's_cgdl_{i}'], C_O2=sv[f'C_O2_cgdl_{i}'], C_N2=C_N2,
                                                   epsilon=epsilon_gdl),
                                          k_th_eff('cgdl', sv[f'T_cgdl_{i + 1}'], C_v=sv[f'C_v_cgdl_{i + 1}'],
                                                   s=sv[f's_cgdl_{i+1}'], C_O2=sv[f'C_O2_cgdl_{i+1}'], C_N2=C_N2,
                                                   epsilon=epsilon_gdl)])
                                   for i in range(1, n_gdl)]

    k_th_eff_cgdl_cgc = k_th_eff('cgdl', sv[f'T_cgdl_{n_gdl}'], C_v=sv[f'C_v_cgdl_{n_gdl}'],
                                 s=sv[f's_cgdl_{n_gdl}'], C_O2=sv[f'C_O2_cgdl_{n_gdl}'], C_N2=C_N2, epsilon=epsilon_gdl)

    return (k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_acl, k_th_eff_acl_mem, k_th_eff_mem_ccl,
            k_th_eff_ccl_cgdl, k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc)
