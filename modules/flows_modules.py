# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import R, F
from modules.transitory_functions import hmean, average, interpolate, Dcap, Da_eff, Dc_eff, h_a, h_c, D


# _____________________________________________________Flow modules_____________________________________________________

def flows_int_values(sv, i_fc, operating_inputs, parameters):
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
    C_v_acl, C_v_ccl = sv['C_v_acl'], sv['C_v_ccl']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_acl, C_O2_ccl, C_N2_a, C_N2_c = sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_N2_a'], sv['C_N2_c']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']
    # Extraction of the operating inputs and the parameters
    T_des = operating_inputs['T_des']
    epsilon_gdl, epsilon_cl = parameters['epsilon_gdl'], parameters['epsilon_cl']
    epsilon_mpl, epsilon_c = parameters['epsilon_mpl'], parameters['epsilon_c']
    e, Hacl, Hccl, Hmem = parameters['e'], parameters['Hacl'], parameters['Hccl'], parameters['Hmem']
    Hgdl, Hmpl, Wagc, Wcgc = parameters['Hgdl'], parameters['Hmpl'], parameters['Wagc'], parameters['Wcgc']
    Hagc, Hcgc, Htl = parameters['Hagc'], parameters['Hcgc'], parameters['Htl']
    epsilon_atl, epsilon_ctl = parameters['epsilon_atl'], parameters['epsilon_ctl']
    nb_gc, nb_gdl, nb_tl, nb_mpl = parameters['nb_gc'], parameters['nb_gdl'], parameters['nb_tl'], parameters['nb_mpl']

    # Transitory parameter
    H_gdl_node = Hgdl / nb_gdl
    H_mpl_node = Hmpl / nb_mpl
    H_tl_node = Htl / nb_tl

    # Pressures in the stack
    Pagc = [None] + [(sv[f'C_v_agc_{i}'] + sv[f'C_H2_agc_{i}'] + C_N2_a) * R * sv[f'T_agc_{i}'] for i in range(1, nb_gc + 1)]
    Pagdl = [None] + [(sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}'] + C_N2_a) * R * sv[f'T_agdl_{i}'] for i in range(1, nb_gdl + 1)]
    Patl = [None] + [(sv[f'C_v_atl_{i}'] + sv[f'C_H2_atl_{i}'] + C_N2_a) * R * sv[f'T_atl_{i}'] for i in range(1, nb_tl + 1)]
    Pampl = [None] + [(sv[f'C_v_ampl_{i}'] + sv[f'C_H2_ampl_{i}'] + C_N2_a) * R * sv[f'T_ampl_{i}'] for i in range(1, nb_mpl + 1)]
    Pacl = (C_v_acl + C_H2_acl + C_N2_a) * R * T_acl
    Pccl = (C_v_ccl + C_O2_ccl + C_N2_c) * R * T_ccl
    Pcmpl = [None] + [(sv[f'C_v_cmpl_{i}'] + sv[f'C_O2_cmpl_{i}'] + C_N2_c) * R * sv[f'T_cmpl_{i}'] for i in range(1, nb_mpl + 1)]
    Pctl = [None] + [(sv[f'C_v_ctl_{i}'] + sv[f'C_O2_ctl_{i}'] + C_N2_c) * R * sv[f'T_ctl_{i}'] for i in range(1, nb_tl + 1)]
    Pcgdl = [None] + [(sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2_c) * R * sv[f'T_cgdl_{i}'] for i in range(1, nb_gdl + 1)]
    Pcgc = [None] + [(sv[f'C_v_cgc_{i}'] + sv[f'C_O2_cgc_{i}'] + C_N2_c) * R * sv[f'T_cgc_{i}'] for i in range(1, nb_gc + 1)]

    # Weighted mean values ...
    #       ... of the EOD flow of water in the membrane
    J_EOD_acl_mem = 2.5 / 22 * i_fc / F * interpolate([lambda_acl, lambda_mem], [Hacl, Hmem])
    J_EOD_mem_ccl = 2.5 / 22 * i_fc / F * interpolate([lambda_mem, lambda_ccl], [Hmem, Hccl])
    #       ... of the diffusion coefficient of water in the membrane
    D_acl_mem = hmean([D(lambda_acl), D(lambda_mem)], weights = [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_mem_ccl = hmean([D(lambda_mem), D(lambda_ccl)], weights = [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])
    #       ... of the capillary coefficient
    D_cap_agdl_agdl = [None] + [hmean([Dcap('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                         Dcap('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, nb_gdl)]
    D_cap_agdl_atl = hmean([Dcap('gdl', sv[f's_agdl_{nb_gdl}'], sv[f'T_agdl_{nb_gdl}'], epsilon_gdl, e,
                                   epsilon_c=epsilon_c),
                                  Dcap('atl', sv['s_atl_1'], sv['T_atl_1'], epsilon_atl[1], e,
                                       epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=1)],
                             weights=[H_gdl_node / (H_gdl_node + H_tl_node), H_tl_node / (H_gdl_node + H_tl_node)])
    D_cap_atl_atl = [None] + [hmean([Dcap('atl', sv[f's_atl_{i}'], sv[f'T_atl_{i}'], epsilon_atl[i], e,
                                          epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=i),
                                       Dcap('atl', sv[f's_atl_{i + 1}'], sv[f'T_atl_{i + 1}'], epsilon_atl[i+1], e,
                                          epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=i + 1)])
                                for i in range(1, nb_tl)]
    D_cap_atl_ampl = hmean([Dcap('atl', sv[f's_atl_{nb_tl}'], sv[f'T_atl_{nb_tl}'], epsilon_atl[nb_tl], e,
                                  epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=nb_tl),
                             Dcap('mpl', sv['s_ampl_1'], sv['T_ampl_1'], epsilon_mpl, e)],
                            weights=[H_tl_node / (H_tl_node + H_mpl_node), H_mpl_node / (H_tl_node + H_mpl_node)])
    D_cap_ampl_ampl = [None] + [hmean([Dcap('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], epsilon_mpl, e),
                                         Dcap('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], epsilon_mpl, e)])
                                for i in range(1, nb_mpl)]
    D_cap_ampl_acl = hmean([Dcap('mpl', sv[f's_ampl_{nb_mpl}'], sv[f'T_ampl_{nb_mpl}'], epsilon_mpl, e),
                                    Dcap('cl', s_acl, T_acl, epsilon_cl, e)],
                             weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    D_cap_cgdl_cgdl = [None] + [hmean([Dcap('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                        Dcap('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, nb_gdl)]
    D_cap_cmpl_ctl = hmean([Dcap('mpl', sv[f's_cmpl_{nb_mpl}'], sv[f'T_cmpl_{nb_mpl}'], epsilon_mpl, e),
                                  Dcap('ctl', sv['s_ctl_1'], sv['T_ctl_1'], epsilon_ctl[1], e, epsilon_c=epsilon_c,
                                       n_tl=nb_tl, Htl=Htl, node=1)],
                            weights = [H_mpl_node / (H_mpl_node + H_tl_node), H_tl_node / (H_mpl_node + H_tl_node)])
    D_cap_ctl_ctl = [None] + [hmean([Dcap('ctl', sv[f's_ctl_{i}'], sv[f'T_ctl_{i}'], epsilon_ctl[i], e,
                                          epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=i),
                                       Dcap('ctl', sv[f's_ctl_{i + 1}'], sv[f'T_ctl_{i + 1}'], epsilon_ctl[i+1],
                                            e, epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=i + 1)])
                                for i in range(1, nb_tl)]
    D_cap_ctl_cgdl = hmean([Dcap('ctl', sv[f's_ctl_{nb_tl}'], sv[f'T_ctl_{nb_tl}'], epsilon_ctl[nb_tl], e,
                                        epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=nb_tl),
                             Dcap('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], epsilon_gdl, e,
                                  epsilon_c=epsilon_c)],
                            weights=[H_tl_node / (H_tl_node + H_gdl_node), H_gdl_node / (H_tl_node + H_gdl_node)])
    D_cap_cmpl_cmpl = [None] + [hmean([Dcap('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], epsilon_mpl, e),
                                        Dcap('mpl', sv[f's_cmpl_{i+1}'], sv[f'T_cmpl_{i+1}'], epsilon_mpl, e)])
                                for i in range(1, nb_mpl)]
    D_cap_ccl_cmpl = hmean([Dcap('cl', s_ccl, T_ccl, epsilon_cl, e),
                                    Dcap('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], epsilon_mpl, e)],
                             weights=[Hccl / (Hccl + H_mpl_node), H_mpl_node / (Hccl + H_mpl_node)])
    #       ... of the effective diffusion coefficient
    Da_eff_agdl_agdl = [None] + [hmean([Da_eff('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], Pagdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                          Da_eff('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], Pagdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, nb_gdl)]
    Da_eff_agdl_atl = hmean([Da_eff('gdl', sv[f's_agdl_{nb_gdl}'], sv[f'T_agdl_{nb_gdl}'], Pagdl[nb_gdl],
                                          epsilon_gdl, epsilon_c = epsilon_c),
                                      Da_eff('tl', sv['s_atl_1'], sv['T_atl_1'], Patl[1], epsilon_atl[1],
                                             epsilon_c = epsilon_c, n_tl=nb_tl, Htl=Htl, node=1)],
                               weights = [H_gdl_node / (H_gdl_node + H_tl_node), H_tl_node / (H_gdl_node + H_tl_node)])
    Da_eff_atl_atl = [None] + [hmean([Da_eff('tl', sv[f's_atl_{i}'], sv[f'T_atl_{i}'], Patl[i], epsilon_atl[i],
                                             epsilon_c = epsilon_c, n_tl=nb_tl, Htl=Htl, node=i),
                                        Da_eff('tl', sv[f's_atl_{i + 1}'], sv[f'T_atl_{i + 1}'], Patl[i + 1],
                                               epsilon_atl[i + 1], epsilon_c = epsilon_c, n_tl=nb_tl, Htl=Htl, node=i + 1)])
                               for i in range(1, nb_tl)]
    Da_eff_atl_ampl = hmean([Da_eff('tl', sv[f's_atl_{nb_tl}'], sv[f'T_atl_{nb_tl}'], Patl[nb_tl],
                                     epsilon_atl[nb_tl], epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=nb_tl),
                              Da_eff('mpl', sv['s_ampl_1'], sv['T_ampl_1'], Pampl[1], epsilon_mpl)],
                             weights=[H_tl_node / (H_tl_node + H_mpl_node), H_mpl_node / (H_tl_node + H_mpl_node)])
    Da_eff_ampl_ampl = [None] + [hmean([Da_eff('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], Pampl[i],
                                              epsilon_mpl),
                                        Da_eff('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], Pampl[i+1],
                                              epsilon_mpl)]) for i in range(1, nb_mpl)]
    Da_eff_ampl_acl = hmean([Da_eff('mpl', sv[f's_ampl_{nb_mpl}'], sv[f'T_ampl_{nb_mpl}'], Pampl[nb_mpl],
                                            epsilon_mpl),
                                     Da_eff('cl', s_acl, T_acl, Pacl, epsilon_cl)],
                              weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    Dc_eff_ccl_cmpl = hmean([Dc_eff('cl', s_ccl, T_ccl, Pccl, epsilon_cl),
                               Dc_eff('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], Pcmpl[1], epsilon_mpl)],
                              weights=[Hccl / (H_mpl_node + Hccl), H_mpl_node / (H_mpl_node + Hccl)])
    Dc_eff_cmpl_cmpl = [None] + [hmean([Dc_eff('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], Pcmpl[i],
                                                 epsilon_mpl),
                                          Dc_eff('mpl', sv[f's_cmpl_{i + 1}'], sv[f'T_cmpl_{i + 1}'], Pcmpl[i + 1],
                                                 epsilon_mpl)]) for i in range(1, nb_mpl)]
    Dc_eff_cmpl_ctl = hmean([Dc_eff('mpl', sv[f's_cmpl_{nb_mpl}'], sv[f'T_cmpl_{nb_mpl}'], Pcmpl[nb_mpl],
                                       epsilon_mpl),
                                   Dc_eff('tl', sv['s_ctl_1'], sv['T_ctl_1'], Pctl[1], epsilon_ctl[1],
                                          epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=1)],
                               weights=[H_mpl_node / (H_mpl_node + H_tl_node), H_tl_node / (H_mpl_node + H_tl_node)])
    Dc_eff_ctl_ctl = [None] + [hmean([Dc_eff('tl', sv[f's_ctl_{i}'], sv[f'T_ctl_{i}'], Pctl[i], epsilon_ctl[i],
                                             epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=i),
                                      Dc_eff('tl', sv[f's_ctl_{i + 1}'], sv[f'T_ctl_{i + 1}'], Pctl[i + 1],
                                             epsilon_ctl[i+1], epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=i + 1)])
                               for i in range(1, nb_tl)]
    Dc_eff_ctl_cgdl = hmean([Dc_eff('tl', sv[f's_ctl_{nb_tl}'], sv[f'T_ctl_{nb_tl}'], Pctl[nb_tl],
                                     epsilon_ctl[nb_tl], epsilon_c=epsilon_c, n_tl=nb_tl, Htl=Htl, node=nb_tl),
                              Dc_eff('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], Pcgdl[1],
                                     epsilon_gdl, epsilon_c=epsilon_c)],
                             weights=[H_tl_node / (H_tl_node + H_gdl_node), H_gdl_node / (H_tl_node + H_gdl_node)])
    Dc_eff_cgdl_cgdl = [None] + [hmean([Dc_eff('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], Pcgdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                        Dc_eff('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], Pcgdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, nb_gdl)]
    #       ... of the temperature
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])

    return (H_gdl_node, H_tl_node, H_mpl_node, Pagc, Pcgc, J_EOD_acl_mem, J_EOD_mem_ccl, D_acl_mem, D_mem_ccl,
            D_cap_agdl_agdl, D_cap_agdl_atl, D_cap_atl_atl, D_cap_atl_ampl, D_cap_ampl_ampl, D_cap_ampl_acl,
            D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_ctl, D_cap_ctl_ctl, D_cap_ctl_cgdl, D_cap_cgdl_cgdl,
            Da_eff_agdl_agdl, Da_eff_agdl_atl, Da_eff_atl_atl, Da_eff_atl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl,
            Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl, Dc_eff_cmpl_ctl, Dc_eff_ctl_ctl, Dc_eff_ctl_cgdl, Dc_eff_cgdl_cgdl,
            T_acl_mem_ccl)
