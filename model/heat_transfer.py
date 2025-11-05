# -*- coding: utf-8 -*-

"""This file represents all the heat transfers occuring inside the fuel cell system.
It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import F, delta_s_HOR, delta_s_ORR
from modules.transitory_functions import d_dx, sigma_p_eff, sigma_e_eff, delta_h_liq, delta_h_abs
from modules.heat_modules import heat_transfer_int_values


# ____________________________________________________Heat transfers____________________________________________________

def calculate_heat_transfers(sv, i_fc, operating_inputs, parameters, S_abs, Sl, **kwargs):
    """This function calculates the heat transfers occurring inside the fuel cell system.

    Parameters
    ----------
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    i_fc : float
        Fuel cell current density at time t (A.m-2).
    parameters : dict
        Parameters of the fuel cell model.
    S_abs : dict
        Water absorption rates from the CL to the membrane (mol.m-3.s-1).
    Sl : dict
        Liquid water absorption rates (mol.m-3.s-1).

    Returns
    -------
    dict
        Heat transfers occuring inside the fuel cell system.
    """

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of the variables
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    s_acl, s_ccl, eta_c = sv['s_acl'], sv['s_ccl'], sv['eta_c']

    # Extraction of the operating inputs and parameters
    T_des = operating_inputs['T_des']
    epsilon_mc, epsilon_gdl, epsilon_cl = parameters['epsilon_mc'], parameters['epsilon_gdl'], parameters['epsilon_cl']
    epsilon_mpl, epsilon_c = parameters['epsilon_mpl'], parameters['epsilon_c']
    nb_gdl, nb_mpl = parameters['nb_gdl'], parameters['nb_mpl']
    Hmem, Hgdl, Hmpl = parameters['Hmem'], parameters['Hgdl'], parameters['Hmpl']
    Hacl, Hccl = parameters['Hacl'], parameters['Hccl']

    # Intermediate values
    (Hgdl_node, Hmpl_node, k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_ampl, k_th_eff_ampl_ampl,
     k_th_eff_ampl_acl, k_th_eff_acl_mem, k_th_eff_mem_ccl, k_th_eff_ccl_cmpl, k_th_eff_cmpl_cmpl, k_th_eff_cmpl_cgdl,
     k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc) = heat_transfer_int_values(sv, parameters)

    # ______________________________________________Heat flows (J.m-2.s-1)______________________________________________

    # Anode side
    T_agc_mean = T_des
    T_cgc_mean = T_des
    Jt_agc_agdl = - k_th_eff_agc_agdl * d_dx(y_minus = T_agc_mean, y_plus = sv['T_agdl_1'],
                                             dx = Hgdl_node / 2)
    Jt_agdl_agdl = {f'agdl_agdl_{i}': -k_th_eff_agdl_agdl[i] * d_dx(y_minus = sv[f'T_agdl_{i}'], y_plus = sv[f'T_agdl_{i+1}'],
                                                                    dx = Hgdl_node / 2)
                    for i in range(1, nb_gdl)}

    Jt_agdl_ampl = - k_th_eff_agdl_ampl * d_dx(y_minus = sv[f'T_agdl_{nb_gdl}'], y_plus = sv['T_ampl_1'],
                                             dx_minus = Hgdl_node / 2, dx_plus = Hmpl_node / 2)
    Jt_ampl_ampl = {f'ampl_ampl_{i}': -k_th_eff_ampl_ampl[i] * d_dx(y_minus = sv[f'T_ampl_{i}'], y_plus = sv[f'T_ampl_{i+1}'],
                                                                    dx = Hmpl_node / 2)
                    for i in range(1, nb_mpl)}
    Jt_ampl_acl = - k_th_eff_ampl_acl * d_dx(y_minus = sv[f'T_ampl_{nb_mpl}'], y_plus = T_acl,
                                             dx_minus = Hmpl_node / 2, dx_plus = Hacl / 2)

    # Membrane side
    Jt_acl_mem = - k_th_eff_acl_mem * d_dx(y_minus = T_acl, y_plus = T_mem,
                                           dx_minus = Hacl / 2, dx_plus = Hmem / 2)
    Jt_mem_ccl = - k_th_eff_mem_ccl * d_dx(y_minus = T_mem, y_plus = T_ccl,
                                           dx_minus = Hmem / 2, dx_plus = Hccl / 2)

    # Cathode side
    Jt_ccl_cmpl = - k_th_eff_ccl_cmpl * d_dx(y_minus = T_ccl, y_plus = sv['T_cmpl_1'],
                                             dx_minus = Hccl / 2, dx_plus = Hmpl_node / 2)
    Jt_cmpl_cmpl = {f'cmpl_cmpl_{i}': -k_th_eff_cmpl_cmpl[i] * d_dx(y_minus = sv[f'T_cmpl_{i}'], y_plus = sv[f'T_cmpl_{i+1}'],
                                                                    dx = Hmpl_node / 2)
                    for i in range(1, nb_mpl)}
    Jt_cmpl_cgdl = - k_th_eff_cmpl_cgdl * d_dx(y_minus = sv[f'T_cmpl_{nb_mpl}'], y_plus = sv['T_cgdl_1'],
                                             dx_minus = Hmpl_node, dx_plus = Hgdl_node)
    Jt_cgdl_cgdl = {f'cgdl_cgdl_{i}': -k_th_eff_cgdl_cgdl[i] * d_dx(y_minus = sv[f'T_cgdl_{i}'], y_plus = sv[f'T_cgdl_{i+1}'],
                                                                    dx = Hgdl_node / 2)
                    for i in range(1, nb_gdl)}

    Jt_cgdl_cgc = - k_th_eff_cgdl_cgc * d_dx(y_minus = sv[f'T_cgdl_{nb_gdl}'], y_plus = T_cgc_mean,
                                             dx = Hgdl_node / 2)

    Jt = {'agc_agdl': Jt_agc_agdl, **Jt_agdl_agdl, 'agdl_ampl': Jt_agdl_ampl, **Jt_ampl_ampl, 'ampl_acl': Jt_ampl_acl,
          'acl_mem': Jt_acl_mem, 'mem_ccl': Jt_mem_ccl, 'ccl_cmpl': Jt_ccl_cmpl, **Jt_cmpl_cmpl,
          'cmpl_cgdl': Jt_cmpl_cgdl, **Jt_cgdl_cgdl, 'cgdl_cgc': Jt_cgdl_cgc}

    # ____________________________________________Heat generated (J.m-3.s-1)____________________________________________

    # The heat dissipated by the electrochemical reaction 2*H2 + O2 -> 2*H2O, in J.m-3.s-1.
    #    It is given by the sum of Peltier and activation heats [vetterFreeOpenReference2019].
    S_r_acl = i_fc / (2 * F * Hacl) # mol.m-3.s-1. It is the amount of hydrogen consumed at the ACL.
    S_r_ccl = i_fc / (4 * F * Hccl) # mol.m-3.s-1. It is the amount of oxygen consumed at the CCL.
    Q_r = {'acl': S_r_acl * T_acl * (-delta_s_HOR),
           'ccl': S_r_ccl * T_ccl * (-delta_s_ORR) + i_fc * eta_c / Hccl}

    # The heat dissipated by the absorption of water from the CL to the membrane, in J.m-3.s-1.
    Q_sorp = {'acl': S_abs['acl'] * (-delta_h_abs(T_acl)),
              'ccl': S_abs['ccl'] * (-delta_h_abs(T_ccl))}

    # The heat dissipated by the liquefaction of vapor water, in J.m-3.s-1.
    Q_liq = {**{f'agdl_{i}': Sl['agdl'][i] * (- delta_h_liq(sv[f'T_agdl_{i}'])) for i in range(1, nb_gdl + 1)},
             **{f'cgdl_{i}': Sl['cgdl'][i] * (- delta_h_liq(sv[f'T_cgdl_{i}'])) for i in range(1, nb_gdl + 1)},
             **{f'ampl_{i}': Sl['ampl'][i] * (- delta_h_liq(sv[f'T_ampl_{i}'])) for i in range(1, nb_mpl + 1)},
             **{f'cmpl_{i}': Sl['cmpl'][i] * (- delta_h_liq(sv[f'T_cmpl_{i}'])) for i in range(1, nb_mpl + 1)},
             'acl': Sl['acl'] * (-delta_h_liq(T_acl)),
             'ccl': Sl['ccl'] * (-delta_h_liq(T_ccl))             }

    # The heat dissipated by the ionic currents (Joule heating + Ohm's law), in J.m-3.s-1.
    Q_p = {'mem': i_fc ** 2 / sigma_p_eff('mem', lambda_mem, T_mem),
           'ccl': i_fc ** 2 / (3 * sigma_p_eff('ccl', lambda_ccl, T_ccl, epsilon_mc))}

    # The heat dissipated by the electric currents (Joule heating + Ohm's law), in J.m-3.s-1.
    Q_e = {**{f'agdl_{i}': i_fc ** 2 / sigma_e_eff('gdl', epsilon_gdl, epsilon_c=epsilon_c) for i in range(1, nb_gdl + 1)},
           **{f'ampl_{i}': i_fc ** 2 / sigma_e_eff('mpl', epsilon_mpl) for i in range(1, nb_mpl + 1)},
           'acl': i_fc ** 2 / sigma_e_eff('cl', epsilon_cl, epsilon_mc=epsilon_mc),
           'ccl': i_fc ** 2 / (3 * sigma_e_eff('cl', epsilon_cl, epsilon_mc=epsilon_mc)),
           **{f'cmpl_{i}': i_fc ** 2 / sigma_e_eff('mpl', epsilon_mpl) for i in range(1, nb_mpl + 1)},
           **{f'cgdl_{i}': i_fc ** 2 / sigma_e_eff('gdl', epsilon_gdl, epsilon_c=epsilon_c) for i in range(1, nb_gdl + 1)}}

    return {'Jt': Jt, 'Q_r': Q_r, 'Q_sorp': Q_sorp, 'Q_liq': Q_liq, 'Q_p': Q_p, 'Q_e': Q_e}