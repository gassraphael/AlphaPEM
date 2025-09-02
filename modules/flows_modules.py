# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import R
from modules.transitory_functions import average, Dcap, Da_eff, Dc_eff, h_a, h_c, D


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
    lambda_acl_mem : float
        Water content in the ACL and the membrane (kg/kg).
    """

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv['C_v_agc'], sv['C_v_acl'], sv['C_v_ccl'], sv['C_v_cgc']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv['C_H2_agc'], sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_O2_cgc']
    C_N2_a, C_N2_c = sv['C_N2_a'], sv['C_N2_c']
    T_agc, T_acl, T_mem, T_ccl, T_cgc = sv['T_agc'], sv['T_acl'], sv['T_mem'], sv['T_ccl'], sv['T_cgc']
    # Extraction of the operating inputs and the parameters
    epsilon_gdl, epsilon_cl = parameters['epsilon_gdl'], parameters['epsilon_cl']
    epsilon_mpl, epsilon_c = parameters['epsilon_mpl'], parameters['epsilon_c']
    e, Hacl, Hccl, Hmem = parameters['e'], parameters['Hacl'], parameters['Hccl'], parameters['Hmem']
    Hgdl, Hmpl, Wagc, Wcgc = parameters['Hgdl'], parameters['Hmpl'], parameters['Wagc'], parameters['Wcgc']
    Hagc, Hcgc, n_gdl, n_mpl = parameters['Hagc'], parameters['Hcgc'], parameters['n_gdl'], parameters['n_mpl']

    # Transitory parameter
    H_gdl_node = Hgdl / n_gdl
    H_mpl_node = Hmpl / n_mpl

    # Pressures in the stack
    Pagc = (C_v_agc + C_H2_agc + C_N2_a) * R * T_agc
    Pagdl = [None] + [(sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}'] + C_N2_a) * R * sv[f'T_agdl_{i}'] for i in range(1, n_gdl + 1)]
    Pampl = [None] + [(sv[f'C_v_ampl_{i}'] + sv[f'C_H2_ampl_{i}'] + C_N2_a) * R * sv[f'T_ampl_{i}'] for i in range(1, n_mpl + 1)]
    Pacl = (C_v_acl + C_H2_acl + C_N2_a) * R * T_acl
    Pccl = (C_v_ccl + C_O2_ccl + C_N2_c) * R * T_ccl
    Pcmpl = [None] + [(sv[f'C_v_cmpl_{i}'] + sv[f'C_O2_cmpl_{i}'] + C_N2_c) * R * sv[f'T_cmpl_{i}'] for i in range(1, n_mpl + 1)]
    Pcgdl = [None] + [(sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2_c) * R * sv[f'T_cgdl_{i}'] for i in range(1, n_gdl + 1)]
    Pcgc = (C_v_cgc + C_O2_cgc + C_N2_c) * R * T_cgc

    # Weighted mean values ...
    #       ... of the water content
    lambda_acl_mem = average([lambda_acl, lambda_mem], weights = [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    lambda_mem_ccl = average([lambda_mem, lambda_ccl], weights = [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])
    #       ... of the diffusion coefficient of water in the membrane
    D_acl_mem = average([D(lambda_acl), D(lambda_mem)], weights = [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_mem_ccl = average([D(lambda_mem), D(lambda_ccl)], weights = [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])
    #       ... of the capillary coefficient
    D_cap_agdl_agdl = [None] + [average([Dcap('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                         Dcap('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, n_gdl)]
    D_cap_agdl_ampl = average([Dcap('gdl', sv[f's_agdl_{n_gdl}'], sv[f'T_agdl_{n_gdl}'], epsilon_gdl, e,
                                   epsilon_c=epsilon_c),
                                     Dcap('mpl', sv['s_ampl_1'], sv['T_ampl_1'], epsilon_mpl, e)],
                             weights=[H_gdl_node / (H_gdl_node + H_mpl_node), H_mpl_node / (H_gdl_node + H_mpl_node)])
    D_cap_ampl_ampl = [None] + [average([Dcap('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], epsilon_mpl, e),
                                         Dcap('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], epsilon_mpl, e)])
                                for i in range(1, n_mpl)]
    D_cap_ampl_acl = average([Dcap('mpl', sv[f's_ampl_{n_mpl}'], sv[f'T_ampl_{n_mpl}'], epsilon_mpl, e),
                                    Dcap('cl', s_acl, T_acl, epsilon_cl, e)],
                             weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    D_cap_cgdl_cgdl = [None] + [average([Dcap('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                        Dcap('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, n_gdl)]
    D_cap_cmpl_cgdl = average([Dcap('mpl', sv[f's_cmpl_{n_mpl}'], sv[f'T_cmpl_{n_mpl}'], epsilon_mpl, e),
                                     Dcap('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], epsilon_gdl, e,
                                           epsilon_c=epsilon_c)],
                            weights = [H_mpl_node / (H_gdl_node + H_mpl_node), H_gdl_node / (H_gdl_node + H_mpl_node)])
    D_cap_cmpl_cmpl = [None] + [average([Dcap('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], epsilon_mpl, e),
                                        Dcap('mpl', sv[f's_cmpl_{i+1}'], sv[f'T_cmpl_{i+1}'], epsilon_mpl, e)])
                                for i in range(1, n_mpl)]
    D_cap_ccl_cmpl = average([Dcap('cl', s_ccl, T_ccl, epsilon_cl, e),
                                    Dcap('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], epsilon_mpl, e)],
                             weights=[Hccl / (H_mpl_node + Hccl), H_mpl_node / (H_mpl_node + Hccl)])
    #       ... of the effective diffusion coefficient between the gas channel and the gas diffusion layer
    ha_Da_eff_agc_agdl = average([h_a(Pagc, T_agc, Wagc, Hagc) * Hagc,
                                        Da_eff('gdl', sv['s_agdl_1'], sv['T_agdl_1'], Pagdl[1], epsilon_gdl,
                                             epsilon_c = epsilon_c)],
                               weights = [Hagc / (Hagc + H_gdl_node), H_gdl_node / (Hagc + H_gdl_node)])
    hc_Dc_eff_cgdl_cgc = average([h_c(Pcgc, T_cgc, Wcgc, Hcgc) * Hcgc,
                                        Dc_eff('gdl', sv[f's_cgdl_{n_gdl}'], sv[f'T_cgdl_{n_gdl}'], Pcgdl[n_gdl],
                                             epsilon_gdl, epsilon_c = epsilon_c)],
                                 weights = [Hcgc / (Hcgc + H_gdl_node), H_gdl_node / (Hcgc + H_gdl_node)])
    #       ... of the effective diffusion coefficient
    Da_eff_agdl_agdl = [None] + [average([Da_eff('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], Pagdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                          Da_eff('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], Pagdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, n_gdl)]
    Da_eff_agdl_ampl = average([Da_eff('gdl', sv[f's_agdl_{n_gdl}'], sv[f'T_agdl_{n_gdl}'], Pagdl[n_gdl],
                                          epsilon_gdl, epsilon_c = epsilon_c),
                                      Da_eff('mpl', sv['s_ampl_1'], sv['T_ampl_1'], Pampl[1], epsilon_mpl)],
                               weights = [H_gdl_node / (H_gdl_node + H_mpl_node), H_mpl_node / (H_gdl_node + H_mpl_node)])
    Da_eff_ampl_ampl = [None] + [average([Da_eff('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], Pampl[i],
                                              epsilon_mpl),
                                        Da_eff('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], Pampl[i+1],
                                              epsilon_mpl)]) for i in range(1, n_mpl)]
    Da_eff_ampl_acl = average([Da_eff('mpl', sv[f's_ampl_{n_mpl}'], sv[f'T_ampl_{n_mpl}'], Pampl[n_mpl],
                                            epsilon_mpl),
                                     Da_eff('cl', s_acl, T_acl, Pacl, epsilon_cl)],
                              weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    Dc_eff_cgdl_cgdl = [None] + [average([Dc_eff('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], Pcgdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                        Dc_eff('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], Pcgdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, n_gdl)]
    Dc_eff_cmpl_cgdl = average([Dc_eff('mpl', sv[f's_cmpl_{n_mpl}'], sv[f'T_cmpl_{n_mpl}'], Pcmpl[n_mpl],
                                              epsilon_mpl),
                                      Dc_eff('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], Pcgdl[1],
                                          epsilon_gdl, epsilon_c = epsilon_c)],
                               weights = [H_mpl_node / (H_gdl_node + H_mpl_node), H_gdl_node / (H_gdl_node + H_mpl_node)])
    Dc_eff_cmpl_cmpl = [None] + [average([Dc_eff('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], Pcmpl[i],
                                              epsilon_mpl),
                                        Dc_eff('mpl', sv[f's_cmpl_{i+1}'], sv[f'T_cmpl_{i+1}'], Pcmpl[i+1],
                                              epsilon_mpl)]) for i in range(1, n_mpl)]
    Dc_eff_ccl_cmpl = average([Dc_eff('cl', s_ccl, T_ccl, Pccl, epsilon_cl),
                                     Dc_eff('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], Pcmpl[1], epsilon_mpl)],
                              weights=[Hccl / (H_mpl_node + Hccl), H_mpl_node / (H_mpl_node + Hccl)])
    #       ... of the temperature
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])

    return (H_gdl_node, H_mpl_node, Pagc, Pcgc, lambda_acl_mem, lambda_mem_ccl, D_acl_mem, D_mem_ccl, D_cap_agdl_agdl,
            D_cap_agdl_ampl, D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_cgdl_cgdl, D_cap_cmpl_cgdl, D_cap_cmpl_cmpl,
            D_cap_ccl_cmpl, ha_Da_eff_agc_agdl, hc_Dc_eff_cgdl_cgc, Da_eff_agdl_agdl, Da_eff_agdl_ampl,
            Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_cgdl_cgdl, Dc_eff_cmpl_cgdl, Dc_eff_cmpl_cmpl, Dc_eff_ccl_cmpl,
            T_acl_mem_ccl)
