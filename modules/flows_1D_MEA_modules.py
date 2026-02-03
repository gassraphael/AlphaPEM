# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import R, M_H2, M_O2, M_N2, M_H2O
from modules.transitory_functions import (hmean, average, Dcap, Da_eff, Dc_eff, D_lambda, D_lambda_eff, D_EOD,
                                          D_EOD_eff, epsilon_cl, Pcap)


# _____________________________________________________Flow modules_____________________________________________________

def flows_1D_MEA_int_values(sv, i_fc, parameters):
    """This functions calculates intermediate values for the flows calculation.

    Parameters
    ----------
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    i_fc : float
        Current density of the fuel cell (A/m²).
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
    C_N2_agc, C_N2_cgc = sv['C_N2_agc'], sv['C_N2_cgc']
    T_agc, T_acl, T_mem, T_ccl, T_cgc = sv['T_agc'], sv['T_acl'], sv['T_mem'], sv['T_ccl'], sv['T_cgc']
    # Extraction of the operating inputs and the parameters
    epsilon_gdl, epsilon_mpl, epsilon_c = parameters['epsilon_gdl'], parameters['epsilon_mpl'], parameters['epsilon_c']
    e, Hacl, Hccl, Hmem = parameters['e'], parameters['Hacl'], parameters['Hccl'], parameters['Hmem']
    Hgdl, Hmpl, Wagc, Wcgc = parameters['Hgdl'], parameters['Hmpl'], parameters['Wagc'], parameters['Wcgc']
    nb_gdl, nb_mpl = parameters['nb_gdl'], parameters['nb_mpl']

    # Transitory parameter
    H_gdl_node = Hgdl / nb_gdl
    H_mpl_node = Hmpl / nb_mpl

    # Pressures in the stack
    Pagc = (C_v_agc + C_H2_agc + C_N2_agc) * R * T_agc
    Pagdl = [None] + [(sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}'] + C_N2_agc) * R * sv[f'T_agdl_{i}'] for i in range(1, nb_gdl + 1)]
    Pampl = [None] + [(sv[f'C_v_ampl_{i}'] + sv[f'C_H2_ampl_{i}'] + C_N2_agc) * R * sv[f'T_ampl_{i}'] for i in range(1, nb_mpl + 1)]
    Pacl = (C_v_acl + C_H2_acl + C_N2_agc) * R * T_acl
    Pccl = (C_v_ccl + C_O2_ccl + C_N2_cgc) * R * T_ccl
    Pcmpl = [None] + [(sv[f'C_v_cmpl_{i}'] + sv[f'C_O2_cmpl_{i}'] + C_N2_cgc) * R * sv[f'T_cmpl_{i}'] for i in range(1, nb_mpl + 1)]
    Pcgdl = [None] + [(sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2_cgc) * R * sv[f'T_cgdl_{i}'] for i in range(1, nb_gdl + 1)]
    Pcgc = (C_v_cgc + C_O2_cgc + C_N2_cgc) * R * T_cgc

    # Capillary pressures in the stack
    Pcap_agdl = Pcap('gdl', sv['s_agdl_1'], sv['T_agdl_1'], epsilon_gdl, epsilon_c=epsilon_c)
    Pcap_cgdl = Pcap('gdl', sv[f's_cgdl_{nb_gdl}'], sv[f'T_cgdl_{nb_gdl}'], epsilon_gdl, epsilon_c=epsilon_c)

    # Densities in the GC
    rho_agc = C_H2_agc * M_H2 + C_v_agc * M_H2O + C_N2_agc * M_N2
    rho_cgc = C_O2_cgc * M_O2 + C_v_cgc * M_H2O + C_N2_cgc * M_N2

    # Weighted mean values ...
    #       ... of the EOD flow of water in the membrane
    D_eff_EOD_acl_mem = hmean([D_EOD_eff(i_fc, lambda_acl, T_acl, Hacl), D_EOD(lambda_mem)],
                          weights = [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_eff_EOD_mem_ccl = hmean([D_EOD(lambda_mem), D_EOD_eff(i_fc, lambda_ccl, T_ccl, Hccl)],
                          weights = [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])
    #       ... of the diffusion coefficient of water in the membrane
    D_lambda_eff_acl_mem = hmean([D_lambda_eff(lambda_acl, T_acl, Hacl), D_lambda(lambda_mem)],
                          weights = [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_lambda_eff_mem_ccl = hmean([D_lambda(lambda_mem), D_lambda_eff(lambda_ccl, T_ccl, Hccl)],
                          weights = [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])
    #       ... of the capillary coefficient
    D_cap_agdl_agdl = [None] + [hmean([Dcap('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                         Dcap('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, nb_gdl)]
    D_cap_agdl_ampl = hmean([Dcap('gdl', sv[f's_agdl_{nb_gdl}'], sv[f'T_agdl_{nb_gdl}'], epsilon_gdl, e,
                                   epsilon_c=epsilon_c),
                                  Dcap('mpl', sv['s_ampl_1'], sv['T_ampl_1'], epsilon_mpl, e)],
                             weights=[H_gdl_node / (H_gdl_node + H_mpl_node), H_mpl_node / (H_gdl_node + H_mpl_node)])
    D_cap_ampl_ampl = [None] + [hmean([Dcap('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], epsilon_mpl, e),
                                         Dcap('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], epsilon_mpl, e)])
                                for i in range(1, nb_mpl)]
    D_cap_ampl_acl = hmean([Dcap('mpl', sv[f's_ampl_{nb_mpl}'], sv[f'T_ampl_{nb_mpl}'], epsilon_mpl, e),
                                    Dcap('cl', s_acl, T_acl, epsilon_cl(lambda_acl, T_acl, Hacl), e)],
                             weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    D_cap_ccl_cmpl = hmean([Dcap('cl', s_ccl, T_ccl, epsilon_cl(lambda_ccl, T_ccl, Hccl), e),
                            Dcap('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], epsilon_mpl, e)],
                           weights=[Hccl / (Hccl + H_mpl_node), H_mpl_node / (Hccl + H_mpl_node)])
    D_cap_cmpl_cmpl = [None] + [hmean([Dcap('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], epsilon_mpl, e),
                                       Dcap('mpl', sv[f's_cmpl_{i + 1}'], sv[f'T_cmpl_{i + 1}'], epsilon_mpl, e)])
                                for i in range(1, nb_mpl)]
    D_cap_cmpl_cgdl = hmean([Dcap('mpl', sv[f's_cmpl_{nb_mpl}'], sv[f'T_cmpl_{nb_mpl}'], epsilon_mpl, e),
                             Dcap('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], epsilon_gdl, e,
                                  epsilon_c=epsilon_c)],
                            weights=[H_mpl_node / (H_mpl_node + H_gdl_node), H_gdl_node / (H_mpl_node + H_gdl_node)])
    D_cap_cgdl_cgdl = [None] + [hmean([Dcap('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                        Dcap('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, nb_gdl)]
    #       ... of the effective diffusion coefficient
    Da_eff_agdl_agdl = [None] + [hmean([Da_eff('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], Pagdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                          Da_eff('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], Pagdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, nb_gdl)]
    Da_eff_agdl_ampl = hmean([Da_eff('gdl', sv[f's_agdl_{nb_gdl}'], sv[f'T_agdl_{nb_gdl}'], Pagdl[nb_gdl],
                                          epsilon_gdl, epsilon_c = epsilon_c),
                                      Da_eff('mpl', sv['s_ampl_1'], sv['T_ampl_1'], Pampl[1], epsilon_mpl)],
                               weights = [H_gdl_node / (H_gdl_node + H_mpl_node), H_mpl_node / (H_gdl_node + H_mpl_node)])
    Da_eff_ampl_ampl = [None] + [hmean([Da_eff('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], Pampl[i],
                                              epsilon_mpl),
                                        Da_eff('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], Pampl[i+1],
                                              epsilon_mpl)]) for i in range(1, nb_mpl)]
    Da_eff_ampl_acl = hmean([Da_eff('mpl', sv[f's_ampl_{nb_mpl}'], sv[f'T_ampl_{nb_mpl}'], Pampl[nb_mpl],
                                            epsilon_mpl),
                                     Da_eff('cl', s_acl, T_acl, Pacl, epsilon_cl(lambda_acl, T_acl, Hacl))],
                              weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    Dc_eff_ccl_cmpl = hmean([Dc_eff('cl', s_ccl, T_ccl, Pccl, epsilon_cl(lambda_ccl, T_ccl, Hccl)),
                               Dc_eff('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], Pcmpl[1], epsilon_mpl)],
                              weights=[Hccl / (H_mpl_node + Hccl), H_mpl_node / (H_mpl_node + Hccl)])
    Dc_eff_cmpl_cmpl = [None] + [hmean([Dc_eff('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], Pcmpl[i],
                                                 epsilon_mpl),
                                          Dc_eff('mpl', sv[f's_cmpl_{i + 1}'], sv[f'T_cmpl_{i + 1}'], Pcmpl[i + 1],
                                                 epsilon_mpl)]) for i in range(1, nb_mpl)]
    Dc_eff_cmpl_cgdl = hmean([Dc_eff('mpl', sv[f's_cmpl_{nb_mpl}'], sv[f'T_cmpl_{nb_mpl}'], Pcmpl[nb_mpl],
                                       epsilon_mpl),
                                    Dc_eff('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], Pcgdl[1], epsilon_gdl,
                                           epsilon_c=epsilon_c)],
                               weights=[H_mpl_node / (H_mpl_node + H_gdl_node), H_gdl_node / (H_mpl_node + H_gdl_node)])
    Dc_eff_cgdl_cgdl = [None] + [hmean([Dc_eff('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], Pcgdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                        Dc_eff('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], Pcgdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, nb_gdl)]
    #       ... of the temperature
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])

    return (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc, D_eff_EOD_acl_mem,
            D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem, D_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl,
            D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl,
            Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl,
            Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl)
