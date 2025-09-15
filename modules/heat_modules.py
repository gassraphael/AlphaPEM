# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the heat transfer calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from modules.transitory_functions import hmean, k_th_eff


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
    k_th_eff_agc_agdl : float
        Effective thermal diffusivity between the AGC and the first GDL layer (J.m-1.s-1.K-1).
    k_th_eff_agdl_agdl : list of floats
        List of effective thermal diffusivities between adjacent GDL layers on the anode side (J.m-1.s-1.K-1).
    k_th_eff_agdl_ampl : float
        Effective thermal diffusivity between the last GDL layer and the anode microporous layer (J.m-1.s-1.K-1).
    k_th_eff_ampl_ampl : list of floats
        List of effective thermal diffusivities between adjacent microporous layers on the anode side (J.m-1.s-1.K-1).
    k_th_eff_ampl_acl : float
        Effective thermal diffusivity between the anode microporous layer and the anode catalyst layer (J.m-1.s-1.K-1).
    k_th_eff_acl_mem : float
        Effective thermal diffusivity between the anode catalyst layer and the membrane (J.m-1.s-1.K-1).
    k_th_eff_mem_ccl : float
        Effective thermal diffusivity between the membrane and the cathode catalyst layer (J.m-1.s-1.K-1).
    k_th_eff_ccl_cmpl : float
        Effective thermal diffusivity between the cathode catalyst layer and the cathode microporous layer (J.m-1.s-1.K-1).
    k_th_eff_cmpl_cmpl : list of floats
        List of effective thermal diffusivities between adjacent microporous layers on the cathode side (J.m-1.s-1.K-1).
    k_th_eff_cmpl_cgdl : float
        Effective thermal diffusivity between the cathode microporous layer and the first GDL layer (J.m-1.s-1.K-1).
    k_th_eff_cgdl_cgdl : list of floats
        List of effective thermal diffusivities between adjacent GDL layers on the cathode side (J.m-1.s-1.K-1).
    k_th_eff_cgdl_cgc : float
        Effective thermal diffusivity between the last GDL layer and the CGC (J.m-1.s-1.K-1).
    """

    # Extraction of the variables
    C_v_acl, C_v_ccl = sv['C_v_acl'], sv['C_v_ccl']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_acl, C_O2_ccl, C_N2_a, C_N2_c = sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_N2_a'], sv['C_N2_c']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']

    # Extraction of the operating inputs and the parameters
    Hgdl, Hmpl, Hacl, Hccl = parameters['Hgdl'], parameters['Hmpl'], parameters['Hacl'], parameters['Hccl']
    Hmem, epsilon_mc, epsilon_gdl = parameters['Hmem'], parameters['epsilon_mc'], parameters['epsilon_gdl']
    epsilon_cl, epsilon_mpl = parameters['epsilon_cl'], parameters['epsilon_mpl']
    epsilon_c, n_gdl, n_mpl = parameters['epsilon_c'], parameters['n_gdl'], parameters['n_mpl']

    # Weighted harmonic means of the effective thermal diffusivity
    k_th_eff_agc_agdl = k_th_eff('agdl', sv[f'T_agdl_{1}'], C_v=sv[f'C_v_agdl_{1}'], s=sv[f's_agdl_{1}'],
                                 C_H2=sv[f'C_H2_agdl_{1}'], C_N2=C_N2_a, epsilon=epsilon_gdl, epsilon_c=epsilon_c)

    k_th_eff_agdl_agdl = [None] + [hmean([k_th_eff('agdl', sv[f'T_agdl_{i}'], C_v=sv[f'C_v_agdl_{i}'],
                                                   s=sv[f's_agdl_{i}'], C_H2=sv[f'C_H2_agdl_{i}'], C_N2=C_N2_a, epsilon=epsilon_gdl,
                                                     epsilon_c=epsilon_c),
                                          k_th_eff('agdl', sv[f'T_agdl_{i + 1}'], C_v=sv[f'C_v_agdl_{i + 1}'],
                                                   s=sv[f's_agdl_{i + 1}'], C_H2=sv[f'C_H2_agdl_{i + 1}'], C_N2=C_N2_a,
                                                   epsilon=epsilon_gdl, epsilon_c=epsilon_c)])
                                   for i in range(1, n_gdl)]

    k_th_eff_agdl_ampl = hmean([k_th_eff('agdl', sv[f'T_agdl_{n_gdl}'], C_v=sv[f'C_v_agdl_{n_gdl}'],
                                                  s=sv[f's_agdl_{n_gdl}'], C_H2=sv[f'C_H2_agdl_{n_gdl}'], C_N2=C_N2_a,
                                           epsilon=epsilon_gdl, epsilon_c=epsilon_c),
                                        k_th_eff('ampl', sv[f'T_ampl_{1}'], C_v=sv[f'C_v_ampl_{1}'],
                                                 s=sv[f's_ampl_{1}'], C_H2=sv[f'C_H2_ampl_{1}'], C_N2=C_N2_a,
                                                 epsilon=epsilon_mpl)],
                                   weights=[(Hgdl / n_gdl) / 2, (Hmpl / n_mpl) / 2])

    k_th_eff_ampl_ampl = [None] + [hmean([k_th_eff('ampl', sv[f'T_ampl_{i}'], C_v=sv[f'C_v_ampl_{i}'],
                                                     s=sv[f's_ampl_{i}'], C_H2=sv[f'C_H2_ampl_{i}'], C_N2=C_N2_a,
                                                     epsilon=epsilon_mpl),
                                            k_th_eff('ampl', sv[f'T_ampl_{i + 1}'], C_v=sv[f'C_v_ampl_{i + 1}'],
                                                     s=sv[f's_ampl_{i + 1}'], C_H2=sv[f'C_H2_ampl_{i + 1}'], C_N2=C_N2_a,
                                                     epsilon=epsilon_mpl)])
                                   for i in range(1, n_mpl)]

    k_th_eff_ampl_acl = hmean([k_th_eff('ampl', sv[f'T_ampl_{n_mpl}'], C_v=sv[f'C_v_ampl_{n_mpl}'],
                                                  s=sv[f's_ampl_{n_mpl}'], C_H2=sv[f'C_H2_ampl_{n_mpl}'],
                                                 C_N2=C_N2_a, epsilon=epsilon_mpl),
                                 k_th_eff('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl,
                                          C_H2=C_H2_acl, C_N2=C_N2_a, epsilon=epsilon_cl, epsilon_mc=epsilon_mc)],
                                weights=[(Hmpl / n_mpl) / 2, Hacl / 2])

    k_th_eff_acl_mem = hmean([k_th_eff('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl,
                                       C_H2=C_H2_acl, C_N2=C_N2_a, epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
                                      k_th_eff('mem', T_mem, lambdaa=lambda_mem)],
                               weights=[Hacl / 2, Hmem / 2])

    k_th_eff_mem_ccl = hmean([k_th_eff('mem', T_mem, lambdaa=lambda_mem),
                                      k_th_eff('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl,
                                                C_O2=C_O2_ccl, C_N2=C_N2_c, epsilon=epsilon_cl, epsilon_mc=epsilon_mc)],
                               weights=[Hmem / 2, Hccl / 2])

    k_th_eff_ccl_cmpl = hmean([k_th_eff('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl, C_O2=C_O2_ccl,
                                                C_N2=C_N2_c, epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
                                       k_th_eff('cmpl', sv[f'T_cmpl_{1}'], C_v=sv[f'C_v_cmpl_{1}'],
                                                 s=sv[f's_cmpl_{1}'], C_O2=sv[f'C_O2_cmpl_{1}'], C_N2=C_N2_c,
                                                epsilon=epsilon_mpl)],
                                weights=[Hccl / 2, (Hmpl / n_mpl) / 2])

    k_th_eff_cmpl_cmpl = [None] + [hmean([k_th_eff('cmpl', sv[f'T_cmpl_{i}'], C_v=sv[f'C_v_cmpl_{i}'],
                                                   s=sv[f's_cmpl_{i}'], C_O2=sv[f'C_O2_cmpl_{i}'], C_N2=C_N2_c,
                                                   epsilon=epsilon_mpl),
                                            k_th_eff('cmpl', sv[f'T_cmpl_{i + 1}'], C_v=sv[f'C_v_cmpl_{i + 1}'],
                                                    s=sv[f's_cmpl_{i + 1}'], C_O2=sv[f'C_O2_cmpl_{i + 1}'], C_N2=C_N2_c,
                                                    epsilon=epsilon_mpl)])
                                      for i in range(1, n_mpl)]

    k_th_eff_cmpl_cgdl = hmean([k_th_eff('cmpl', sv[f'T_cmpl_{n_mpl}'], C_v=sv[f'C_v_cmpl_{n_mpl}'],
                                                  s=sv[f's_cmpl_{n_mpl}'], C_O2=sv[f'C_O2_cmpl_{n_mpl}'], C_N2=C_N2_c,
                                                  epsilon=epsilon_mpl),
                                        k_th_eff('cgdl', sv['T_cgdl_1'], C_v=sv['C_v_cgdl_1'], s=sv['s_cgdl_1'],
                                          C_O2=sv[f'C_O2_cgdl_1'], C_N2=C_N2_c, epsilon=epsilon_gdl, epsilon_c=epsilon_c)],
                                weights=[(Hmpl / n_mpl) / 2, (Hgdl / n_gdl) / 2])

    k_th_eff_cgdl_cgdl = [None] + [hmean([k_th_eff('cgdl', sv[f'T_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                                                   s=sv[f's_cgdl_{i}'], C_O2=sv[f'C_O2_cgdl_{i}'], C_N2=C_N2_c,
                                                   epsilon=epsilon_gdl, epsilon_c=epsilon_c),
                                            k_th_eff('cgdl', sv[f'T_cgdl_{i + 1}'], C_v=sv[f'C_v_cgdl_{i + 1}'],
                                                   s=sv[f's_cgdl_{i + 1}'], C_O2=sv[f'C_O2_cgdl_{i + 1}'], C_N2=C_N2_c,
                                                   epsilon=epsilon_gdl, epsilon_c=epsilon_c)])
                                   for i in range(1, n_gdl)]

    k_th_eff_cgdl_cgc = k_th_eff('cgdl', sv[f'T_cgdl_{n_gdl}'], C_v=sv[f'C_v_cgdl_{n_gdl}'],
                                 s=sv[f's_cgdl_{n_gdl}'], C_O2=sv[f'C_O2_cgdl_{n_gdl}'], C_N2=C_N2_c, epsilon=epsilon_gdl,
                                 epsilon_c=epsilon_c)

    return (k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_ampl, k_th_eff_ampl_ampl, k_th_eff_ampl_acl,
            k_th_eff_acl_mem, k_th_eff_mem_ccl, k_th_eff_ccl_cmpl, k_th_eff_cmpl_cmpl, k_th_eff_cmpl_cgdl,
            k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc)
