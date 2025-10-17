# -*- coding: utf-8 -*-

"""This file represents all the heat transfers occuring inside the fuel cell system.
It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import F, delta_s_HOR, delta_s_ORR
from modules.transitory_functions import sigma_p_eff, sigma_e_eff, delta_h_liq, delta_h_abs
from modules.heat_modules import heat_transfer_int_values


# ____________________________________________________Heat transfers____________________________________________________

def calculate_heat_transfers(sv, i_fc, operating_inputs, parameters, S_abs_acl, S_abs_ccl, Sl_agdl, Sl_atl, Sl_ampl, Sl_acl, Sl_ccl,
                             Sl_cmpl, Sl_ctl, Sl_cgdl, **kwargs):
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
    S_abs_acl : float
        Water absorption in the ACL (mol.m-3.s-1).
    S_abs_ccl : float
        Water absorption in the CCL (mol.m-3.s-1).
    Sl_agdl : list
        Liquid water generated through vapor water liquefaction in the AGDL (mol.m-3.s-1).
    Sl_atl : list
        Liquid water generated through vapor water liquefaction in the ATL (mol.m-3.s-1).
    Sl_ampl : list
        Liquid water generated through vapor water liquefaction in the AMPL (mol.m-3.s-1).
    Sl_acl : float
        Liquid water generated through vapor water liquefaction in the ACL (mol.m-3.s-1).
    Sl_ccl : float
        Liquid water generated through vapor water liquefaction in the CCL (mol.m-3.s-1).
    Sl_cmpl : list
        Liquid water generated through vapor water liquefaction in the CMPL (mol.m-3.s-1).
    Sl_ctl : list
        Liquid water generated through vapor water liquefaction in the CTL (mol.m-3.s-1).
    Sl_cgdl : list
        Liquid water generated through vapor water liquefaction in the CGDL (mol.m-3.s-1).

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
    n_gdl, n_mpl = parameters['n_gdl'], parameters['n_mpl']
    Hmem, Hgdl, Hmpl = parameters['Hmem'], parameters['Hgdl'], parameters['Hmpl']
    Hacl, Hccl, Htl, n_tl = parameters['Hacl'], parameters['Hccl'], parameters['Htl'], parameters['n_tl']
    epsilon_atl, epsilon_ctl = parameters['epsilon_atl'], parameters['epsilon_ctl']

    # Intermediate values
    (k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_atl, k_th_eff_atl_atl, k_th_eff_atl_ampl, k_th_eff_ampl_ampl,
     k_th_eff_ampl_acl, k_th_eff_acl_mem, k_th_eff_mem_ccl, k_th_eff_ccl_cmpl, k_th_eff_cmpl_cmpl, k_th_eff_cmpl_ctl,
     k_th_eff_ctl_ctl, k_th_eff_ctl_cgdl, k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc) \
            =  heat_transfer_int_values(sv, parameters)

    # ______________________________________________Heat flows (J.m-2.s-1)______________________________________________

    # Anode side
    T_agc_mean = T_des
    T_cgc_mean = T_des
    Jt_agc_agdl = - 2 * k_th_eff_agc_agdl * (sv['T_agdl_1'] - T_agc_mean) / (Hgdl / n_gdl)
    Jt_agdl_agdl = {f'agdl_agdl_{i}': -k_th_eff_agdl_agdl[i] * (sv[f'T_agdl_{i+1}'] - sv[f'T_agdl_{i}']) / (Hgdl/n_gdl)
                    for i in range(1, n_gdl)}

    Jt_agdl_atl = - 2 * k_th_eff_agdl_atl * (sv['T_atl_1'] - sv[f'T_agdl_{n_gdl}']) / (Hgdl / n_gdl + Htl / n_tl)
    Jt_atl_atl = {f'atl_atl_{i}': -k_th_eff_atl_atl[i] * (sv[f'T_atl_{i + 1}'] - sv[f'T_atl_{i}']) / (Htl / n_tl)
                  for i in range(1, n_tl)}
    Jt_atl_ampl = - 2 * k_th_eff_atl_ampl * (sv['T_ampl_1'] - sv[f'T_atl_{n_tl}']) / (Htl / n_tl + Hmpl / n_mpl)
    Jt_ampl_ampl = {f'ampl_ampl_{i}': -k_th_eff_ampl_ampl[i] * (sv[f'T_ampl_{i+1}'] - sv[f'T_ampl_{i}']) / (Hmpl/n_mpl)
                    for i in range(1, n_mpl)}
    Jt_ampl_acl = - 2 * k_th_eff_ampl_acl * (T_acl - sv[f'T_ampl_{n_mpl}']) / (Hmpl / n_mpl + Hacl)

    # Membrane side
    Jt_acl_mem = - 2 * k_th_eff_acl_mem * (T_mem - T_acl) / (Hacl + Hmem)
    Jt_mem_ccl = - 2 * k_th_eff_mem_ccl * (T_ccl - T_mem) / (Hmem + Hccl)

    # Cathode side
    Jt_ccl_cmpl = - 2 * k_th_eff_ccl_cmpl * (sv['T_cmpl_1'] - T_ccl) / (Hccl + Hmpl / n_mpl)
    Jt_cmpl_cmpl = {f'cmpl_cmpl_{i}': -k_th_eff_cmpl_cmpl[i] * (sv[f'T_cmpl_{i+1}'] - sv[f'T_cmpl_{i}']) / (Hmpl/n_mpl)
                    for i in range(1, n_mpl)}
    Jt_cmpl_ctl = - 2 * k_th_eff_cmpl_ctl * (sv['T_ctl_1'] - sv[f'T_cmpl_{n_mpl}']) / (Htl / n_tl + Hgdl / n_gdl)
    Jt_ctl_ctl = {f'ctl_ctl_{i}': -k_th_eff_ctl_ctl[i] * (sv[f'T_ctl_{i + 1}'] - sv[f'T_ctl_{i}']) / (Htl / n_tl)
                  for i in range(1, n_tl)}
    Jt_ctl_cgdl = - 2 * k_th_eff_ctl_cgdl * (sv['T_cgdl_1'] - sv[f'T_ctl_{n_tl}']) / (Htl / n_tl + Hgdl / n_gdl)
    Jt_cgdl_cgdl = {f'cgdl_cgdl_{i}': -k_th_eff_cgdl_cgdl[i] * (sv[f'T_cgdl_{i+1}'] - sv[f'T_cgdl_{i}']) / (Hgdl/n_gdl)
                    for i in range(1, n_gdl)}

    Jt_cgdl_cgc = - 2 * k_th_eff_cgdl_cgc * (T_cgc_mean - sv[f'T_cgdl_{n_gdl}']) / (Hgdl / n_gdl)

    Jt = {'agc_agdl': Jt_agc_agdl, **Jt_agdl_agdl, 'agdl_atl': Jt_agdl_atl, **Jt_atl_atl, 'atl_ampl': Jt_atl_ampl,
          **Jt_ampl_ampl, 'ampl_acl': Jt_ampl_acl, 'acl_mem': Jt_acl_mem, 'mem_ccl': Jt_mem_ccl, 'ccl_cmpl': Jt_ccl_cmpl,
          **Jt_cmpl_cmpl, 'cmpl_ctl': Jt_cmpl_ctl, **Jt_ctl_ctl, 'ctl_cgdl': Jt_ctl_cgdl, **Jt_cgdl_cgdl,
          'cgdl_cgc': Jt_cgdl_cgc}

    # ____________________________________________Heat generated (J.m-3.s-1)____________________________________________

    # The heat dissipated by the electrochemical reaction 2*H2 + O2 -> 2*H2O, in J.m-3.s-1.
    #    It is given by the sum of Peltier and activation heats [vetterFreeOpenReference2019].
    S_r_acl = i_fc / (2 * F * Hacl) # mol.m-3.s-1. It is the amount of hydrogen consumed at the ACL.
    S_r_ccl = i_fc / (4 * F * Hccl) # mol.m-3.s-1. It is the amount of oxygen consumed at the CCL.
    Q_r = {'acl': S_r_acl * T_acl * (-delta_s_HOR),
           'ccl': S_r_ccl * T_ccl * (-delta_s_ORR) + i_fc * eta_c / Hccl}

    # The heat dissipated by the absorption of water from the CL to the membrane, in J.m-3.s-1.
    Q_sorp = {'acl': S_abs_acl * (-delta_h_abs(T_acl)),
              'ccl': S_abs_ccl * (-delta_h_abs(T_ccl))}

    # The heat dissipated by the liquefaction of vapor water, in J.m-3.s-1.
    Q_liq = {**{f'agdl_{i}': Sl_agdl[i] * (- delta_h_liq(sv[f'T_agdl_{i}'])) for i in range(1, n_gdl + 1)},
             **{f'cgdl_{i}': Sl_cgdl[i] * (- delta_h_liq(sv[f'T_cgdl_{i}'])) for i in range(1, n_gdl + 1)},
             **{f'atl_{i}': Sl_atl[i] * (- delta_h_liq(sv[f'T_atl_{i}'])) for i in range(1, n_tl + 1)},
             **{f'ctl_{i}': Sl_ctl[i] * (- delta_h_liq(sv[f'T_ctl_{i}'])) for i in range(1, n_tl + 1)},
             **{f'ampl_{i}': Sl_ampl[i] * (- delta_h_liq(sv[f'T_ampl_{i}'])) for i in range(1, n_mpl + 1)},
             **{f'cmpl_{i}': Sl_cmpl[i] * (- delta_h_liq(sv[f'T_cmpl_{i}'])) for i in range(1, n_mpl + 1)},
             'acl': Sl_acl * (-delta_h_liq(T_acl)),
             'ccl': Sl_ccl * (-delta_h_liq(T_ccl))             }

    # The heat dissipated by the ionic currents (Joule heating + Ohm's law), in J.m-3.s-1.
    Q_p = {'mem': i_fc ** 2 / sigma_p_eff('mem', lambda_mem, T_mem),
           'ccl': i_fc ** 2 / (3 * sigma_p_eff('ccl', lambda_ccl, T_ccl, epsilon_mc))}

    # The heat dissipated by the electric currents (Joule heating + Ohm's law), in J.m-3.s-1.
    Q_e = {**{f'agdl_{i}': i_fc ** 2 / sigma_e_eff('gdl', epsilon_gdl, epsilon_c=epsilon_c) for i in range(1, n_gdl + 1)},
           **{f'atl_{i}': i_fc ** 2 / sigma_e_eff('atl', epsilon_atl[i], epsilon_c=epsilon_c, n_tl=n_tl, Htl=Htl, node=i) for i in range(1, n_tl + 1)},
           **{f'ampl_{i}': i_fc ** 2 / sigma_e_eff('mpl', epsilon_mpl) for i in range(1, n_mpl + 1)},
           'acl': i_fc ** 2 / sigma_e_eff('cl', epsilon_cl, epsilon_mc=epsilon_mc),
           'ccl': i_fc ** 2 / (3 * sigma_e_eff('cl', epsilon_cl, epsilon_mc=epsilon_mc)),
           **{f'cmpl_{i}': i_fc ** 2 / sigma_e_eff('mpl', epsilon_mpl) for i in range(1, n_mpl + 1)},
           **{f'ctl_{i}': i_fc ** 2 / sigma_e_eff('ctl', epsilon_ctl[i], epsilon_c=epsilon_c, n_tl=n_tl, Htl=Htl, node=i) for i in range(1, n_tl + 1)},
           **{f'cgdl_{i}': i_fc ** 2 / sigma_e_eff('gdl', epsilon_gdl, epsilon_c=epsilon_c) for i in range(1, n_gdl + 1)}}

    return {'Jt': Jt, 'Q_r': Q_r, 'Q_sorp': Q_sorp, 'Q_liq': Q_liq, 'Q_p': Q_p, 'Q_e': Q_e}