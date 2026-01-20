# -*- coding: utf-8 -*-

"""This file represents all the matter flows inside the fuel cell system. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import rho_mem, M_eq, F, R, ECSA_0
from model.auxiliaries import auxiliaries
from modules.transitory_functions import (interpolate, d_dx, Dcap, h_a, h_c, lambda_eq, gamma_sorp, Svl, k_H2, k_O2,
                                          epsilon_cl)
from modules.flows_modules import flows_int_values


# ________________________________________________________Flows_________________________________________________________

def calculate_flows(t, sv, control_variables, i_fc, operating_inputs, parameters):
    """This function calculates the flows inside the fuel cell system.

    Parameters
    ----------
    t : float
        Time (s).
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    control_variables : dict
        Variables controlled by the user.
    i_fc : float
        Fuel cell current density at time t (A.m-2).
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    dict
        Flows inside the fuel cell system.
    """

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of the variables
    C_v_acl, C_v_ccl = sv['C_v_acl'], sv['C_v_ccl']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_acl, C_O2_ccl = sv['C_H2_acl'], sv['C_O2_ccl']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']

    # Extraction of the operating inputs and parameters
    T_des = operating_inputs['T_des']
    Aact, Hmem, Hacl, Hccl = parameters['Aact'], parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    Wagc, Wcgc, Hagc, Hcgc = parameters['Wagc'], parameters['Wcgc'], parameters['Hagc'], parameters['Hcgc']
    epsilon_gdl, epsilon_mpl, epsilon_c = parameters['epsilon_gdl'], parameters['epsilon_mpl'], parameters['epsilon_c']
    e, kappa_co = parameters['e'], parameters['kappa_co']
    nb_gc, nb_gdl, nb_mpl = parameters['nb_gc'], parameters['nb_gdl'], parameters['nb_mpl']

    # Intermediate values
    (H_gdl_node, H_mpl_node, Pagc, Pcgc, D_eff_EOD_acl_mem, D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem,
     D_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl, D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl,
     D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl, Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl,
     Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl, Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl) = \
        flows_int_values(sv, i_fc, operating_inputs, parameters)
    C_N2_a_mean = (sum(sv[f'C_N2_agc_{i}'] for i in range(1, nb_gc + 1)) / nb_gc)
    C_N2_c_mean = (sum(sv[f'C_N2_cgc_{i}'] for i in range(1, nb_gc + 1)) / nb_gc)

    # ________________________________________Dissolved water flows (mol.m-2.s-1)_______________________________________

    # Anode side
    J_lambda_acl_mem = D_eff_EOD_acl_mem * interpolate([lambda_acl, lambda_mem], [Hacl, Hmem]) - \
                       rho_mem / M_eq * D_lambda_eff_acl_mem * d_dx(y_minus = lambda_acl, y_plus = lambda_mem,
                                                                    dx_minus = Hacl / 2, dx_plus = Hmem / 2)
    # Cathode side
    J_lambda_mem_ccl = D_eff_EOD_mem_ccl * interpolate([lambda_mem, lambda_ccl], [Hmem, Hccl]) - \
                       rho_mem / M_eq * D_lambda_eff_mem_ccl * d_dx(y_minus = lambda_mem, y_plus = lambda_ccl,
                                                                    dx_minus = Hmem / 2, dx_plus = Hccl / 2)

    # _________________________________________Liquid water flows (kg.m-2.s-1)__________________________________________

    # Anode side
    s_agc = 0  # Dirichlet boundary condition (taken at the agc/agdl border).
    Jl_agc_agdl = - Dcap('gdl', sv['s_agdl_1'], sv['T_agdl_1'], epsilon_gdl, e, epsilon_c=epsilon_c) * \
                    d_dx(y_minus = s_agc, y_plus = sv['s_agdl_1'], dx = H_gdl_node / 2)
    Jl_agdl_agdl = [None] + [- D_cap_agdl_agdl[i] * d_dx(y_minus = sv[f's_agdl_{i}'], y_plus = sv[f's_agdl_{i + 1}'],
                                                         dx = H_gdl_node / 2)
                             for i in range(1, nb_gdl)]
    Jl_agdl_ampl = - D_cap_agdl_ampl * d_dx(y_minus = sv[f's_agdl_{nb_gdl}'], y_plus = sv['s_ampl_1'],
                                          dx_minus = H_gdl_node / 2, dx_plus = H_mpl_node / 2)
    Jl_ampl_ampl = [None] + [- D_cap_ampl_ampl[i] * d_dx(y_minus = sv[f's_ampl_{i}'], y_plus = sv[f's_ampl_{i + 1}'],
                                                         dx = H_mpl_node / 2)
                             for i in range(1, nb_mpl)]
    Jl_ampl_acl = - D_cap_ampl_acl * d_dx(y_minus = sv[f's_ampl_{nb_mpl}'], y_plus = s_acl,
                                          dx_minus = H_mpl_node / 2, dx_plus = Hacl / 2)

    # Cathode side
    s_cgc = 0  # Dirichlet boundary condition (taken at the cgc/cgdl border).
    Jl_ccl_cmpl = - D_cap_ccl_cmpl * d_dx(y_minus = s_ccl, y_plus = sv['s_cmpl_1'],
                                          dx_minus = Hccl / 2, dx_plus = H_mpl_node / 2)
    Jl_cmpl_cmpl = [None] + [- D_cap_cmpl_cmpl[i] * d_dx(y_minus = sv[f's_cmpl_{i}'], y_plus = sv[f's_cmpl_{i + 1}'],
                                                         dx = H_mpl_node / 2)
                            for i in range(1, nb_mpl)]
    Jl_cmpl_cgdl = - D_cap_cmpl_cgdl * d_dx(y_minus = sv[f's_cmpl_{nb_mpl}'], y_plus = sv['s_cgdl_1'],
                                          dx_minus = H_mpl_node / 2, dx_plus = H_gdl_node / 2)
    Jl_cgdl_cgdl = [None] + [- D_cap_cgdl_cgdl[i] * d_dx(y_minus = sv[f's_cgdl_{i}'], y_plus = sv[f's_cgdl_{i + 1}'],
                                                         dx = H_gdl_node / 2)
                             for i in range(1, nb_gdl)]
    Jl_cgdl_cgc = - Dcap('gdl', sv[f's_cgdl_{nb_gdl}'], sv[f'T_cgdl_{nb_gdl}'], epsilon_gdl, e, epsilon_c=epsilon_c) * \
                  d_dx(y_minus = sv[f's_cgdl_{nb_gdl}'], y_plus = s_cgc, dx = H_gdl_node / 2)

    # _____________________________________________Vapor flows (mol.m-2.s-1)____________________________________________

    # Convective vapor flows
    #   Anode side
    Jv_agc_agdl = [None] + [h_a(Pagc[i], T_des, Wagc, Hagc) * (sv[f'C_v_agc_{i}'] - sv['C_v_agdl_1'])
                            for i in range(1, nb_gc + 1)]
    #   Cathode side
    Jv_cgdl_cgc = [None] + [h_c(Pcgc[i], T_des, Wcgc, Hcgc) * (sv[f'C_v_cgdl_{nb_gdl}'] - sv[f'C_v_cgc_{i}'])
                            for i in range(1, nb_gc + 1)]

    # Conductive vapor flows
    #   Anode side
    Jv_agdl_agdl = [None] + [- Da_eff_agdl_agdl[i] * d_dx(y_minus = sv[f'C_v_agdl_{i}'], y_plus = sv[f'C_v_agdl_{i + 1}'],
                                                          dx = H_gdl_node / 2)
                             for i in range(1, nb_gdl)]
    Jv_agdl_ampl = - Da_eff_agdl_ampl * d_dx(y_minus = sv[f'C_v_agdl_{nb_gdl}'], y_plus = sv['C_v_ampl_1'],
                                           dx_minus = H_gdl_node / 2, dx_plus = H_mpl_node / 2)
    Jv_ampl_ampl = [None] + [- Da_eff_ampl_ampl[i] * d_dx(y_minus = sv[f'C_v_ampl_{i}'], y_plus = sv[f'C_v_ampl_{i + 1}'],
                                                          dx = H_mpl_node / 2)
                             for i in range(1, nb_mpl)]
    Jv_ampl_acl = - Da_eff_ampl_acl * d_dx(y_minus = sv[f'C_v_ampl_{nb_mpl}'], y_plus = C_v_acl,
                                           dx_minus = H_mpl_node / 2, dx_plus = Hacl / 2)

    #   Cathode side
    Jv_ccl_cmpl = - Dc_eff_ccl_cmpl * d_dx(y_minus = C_v_ccl, y_plus = sv['C_v_cmpl_1'],
                                           dx_minus = Hccl / 2, dx_plus = H_mpl_node / 2)
    Jv_cmpl_cmpl = [None] + [- Dc_eff_cmpl_cmpl[i] * d_dx(y_minus = sv[f'C_v_cmpl_{i}'], y_plus = sv[f'C_v_cmpl_{i + 1}'],
                                                          dx = H_mpl_node / 2)
                             for i in range(1, nb_mpl)]
    Jv_cmpl_cgdl = - Dc_eff_cmpl_cgdl * d_dx(y_minus = sv[f'C_v_cmpl_{nb_mpl}'], y_plus = sv['C_v_cgdl_1'],
                                           dx_minus = H_mpl_node / 2, dx_plus = H_gdl_node / 2)
    Jv_cgdl_cgdl = [None] + [- Dc_eff_cgdl_cgdl[i] * d_dx(y_minus = sv[f'C_v_cgdl_{i}'], y_plus = sv[f'C_v_cgdl_{i + 1}'],
                                                          dx = H_gdl_node / 2)
                             for i in range(1, nb_gdl)]

    # ______________________________H2 and O2 flows (mol.m-2.s-1 for J, mol.m-3.s-1 for S)______________________________

    # Hydrogen and oxygen consumption
    #   Anode side
    S_H2_reac = i_fc / (2 * F * Hacl)
    S_H2_cros = R * T_acl_mem_ccl / (Hmem * Hacl) * (k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl +
                                                    2 * k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl)
    #   Cathode side
    S_O2_reac = i_fc / (4 * F * Hccl)
    S_O2_cros = R * T_acl_mem_ccl / (Hmem * Hccl) * (k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl +
                                                    1 / 2 * k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl)

    # Conductive-convective H2 and O2 flows
    #   Anode side
    J_H2_agc_agdl = [None] + [h_a(Pagc[i], T_des, Wagc, Hagc) * (sv[f'C_H2_agc_{i}'] - sv['C_H2_agdl_1'])
                              for i in range(1, nb_gc + 1)]
    #   Cathode side
    J_O2_cgdl_cgc = [None] + [h_c(Pcgc[i], T_des, Wcgc, Hcgc) * (sv[f'C_O2_cgdl_{nb_gdl}'] - sv[f'C_O2_cgc_{i}'])
                              for i in range(1, nb_gc + 1)]

    # Conductive H2 and O2 flows
    #   Anode side
    J_H2_agdl_agdl = [None] + [- Da_eff_agdl_agdl[i] * d_dx(y_minus = sv[f'C_H2_agdl_{i}'], y_plus = sv[f'C_H2_agdl_{i+1}'],
                                                            dx = H_gdl_node / 2)
                               for i in range(1, nb_gdl)]
    J_H2_agdl_ampl = - Da_eff_agdl_ampl * d_dx(y_minus = sv[f'C_H2_agdl_{nb_gdl}'], y_plus = sv['C_H2_ampl_1'],
                                             dx_minus = H_gdl_node / 2, dx_plus = H_mpl_node / 2)
    J_H2_ampl_ampl = [None] + [- Da_eff_ampl_ampl[i] * d_dx(y_minus = sv[f'C_H2_ampl_{i}'], y_plus = sv[f'C_H2_ampl_{i + 1}'],
                                                            dx = H_mpl_node / 2)
                               for i in range(1, nb_mpl)]
    J_H2_ampl_acl = - Da_eff_ampl_acl * d_dx(y_minus = sv[f'C_H2_ampl_{nb_mpl}'], y_plus = C_H2_acl,
                                             dx_minus = H_mpl_node / 2, dx_plus = Hacl / 2)

    #   Cathode side
    J_O2_ccl_cmpl = - Dc_eff_ccl_cmpl * d_dx(y_minus = C_O2_ccl, y_plus = sv['C_O2_cmpl_1'],
                                             dx_minus = Hccl / 2, dx_plus = H_mpl_node / 2)
    J_O2_cmpl_cmpl = [None] + [- Dc_eff_cmpl_cmpl[i] * d_dx(y_minus = sv[f'C_O2_cmpl_{i}'], y_plus = sv[f'C_O2_cmpl_{i + 1}'],
                                                            dx = H_mpl_node / 2)
                               for i in range(1, nb_mpl)]
    J_O2_cmpl_cgdl = - Dc_eff_cmpl_cgdl * d_dx(y_minus = sv[f'C_O2_cmpl_{nb_mpl}'], y_plus = sv['C_O2_cgdl_1'],
                                             dx_minus = H_mpl_node / 2, dx_plus = H_gdl_node / 2)
    J_O2_cgdl_cgdl = [None] + [- Dc_eff_cgdl_cgdl[i] * d_dx(y_minus = sv[f'C_O2_cgdl_{i}'], y_plus = sv[f'C_O2_cgdl_{i + 1}'],
                                                            dx = H_gdl_node / 2)
                               for i in range(1, nb_gdl)]

    # __________________________________________Water generated (mol.m-3.s-1)___________________________________________

    # Water produced in the membrane at the CL through the chemical reaction and crossover
    #   Anode side
    Sp_acl = 2 * k_O2(lambda_mem, T_mem, kappa_co) * R * T_acl_mem_ccl / (Hmem * Hacl) * C_O2_ccl
    #   Cathode side
    Sp_ccl = i_fc / (2 * F * Hccl) + k_H2(lambda_mem, T_mem, kappa_co) * R * T_acl_mem_ccl / (Hmem * Hccl) * C_H2_acl

    # Water absorption in the CL:
    #   Anode side
    S_abs_acl = gamma_sorp(C_v_acl, s_acl, lambda_acl, T_acl, Hacl) * rho_mem / M_eq * \
                 (lambda_eq(C_v_acl, s_acl, T_acl) - lambda_acl)
    #   Cathode side
    S_abs_ccl = gamma_sorp(C_v_ccl, s_ccl, lambda_ccl, T_ccl, Hccl) * rho_mem / M_eq * \
                 (lambda_eq(C_v_ccl, s_ccl, T_ccl) - lambda_ccl)

    # Liquid water generated through vapor condensation or degenerated through evaporation
    #   Anode side
    Sl_agdl = [None] + [Svl(element='anode', s=sv[f's_agdl_{i}'], C_v=sv[f'C_v_agdl_{i}'],
                            Ctot=sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}'] + C_N2_a_mean,
                            T=sv[f'T_agdl_{i}'], epsilon=epsilon_gdl) for i in range(1, nb_gdl + 1)]
    Sl_ampl = [None] + [Svl(element='anode', s=sv[f's_ampl_{i}'], C_v=sv[f'C_v_ampl_{i}'],
                            Ctot=sv[f'C_v_ampl_{i}'] + sv[f'C_H2_ampl_{i}'] + C_N2_a_mean,
                            T=sv[f'T_ampl_{i}'], epsilon=epsilon_mpl) for i in range(1, nb_mpl + 1)]
    Sl_acl = Svl(element='anode', s=s_acl, C_v=C_v_acl, Ctot=C_v_acl + C_H2_acl + C_N2_a_mean, T=T_acl,
                 epsilon=epsilon_cl(sv['lambda_acl'], sv['T_acl'], Hacl))
    #   Cathode side
    Sl_ccl = Svl(element='cathode', s=s_ccl, C_v=C_v_ccl, Ctot=C_v_ccl + C_O2_ccl + C_N2_c_mean, T=T_ccl,
                 epsilon=epsilon_cl(sv['lambda_ccl'], sv['T_ccl'], Hccl))
    Sl_cmpl = [None] + [Svl(element='cathode', s=sv[f's_cmpl_{i}'], C_v=sv[f'C_v_cmpl_{i}'],
                            Ctot=sv[f'C_v_cmpl_{i}'] + sv[f'C_O2_cmpl_{i}'] + C_N2_c_mean,
                            T=sv[f'T_cmpl_{i}'], epsilon=epsilon_mpl) for i in range(1, nb_mpl + 1)]
    Sl_cgdl = [None] + [Svl(element='cathode', s=sv[f's_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                            Ctot=sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2_c_mean,
                            T=sv[f'T_cgdl_{i}'], epsilon=epsilon_gdl) for i in range(1, nb_gdl + 1)]

    # Vapor generated through liquid water evaporation or degenerated through condensation
    #   Anode side
    Sv_agdl = [None] + [-x for x in Sl_agdl[1:]]
    Sv_ampl = [None] + [-x for x in Sl_ampl[1:]]
    Sv_acl = - Sl_acl
    #   Cathode side
    Sv_ccl = - Sl_ccl
    Sv_cmpl = [None] + [-x for x in Sl_cmpl[1:]]
    Sv_cgdl = [None] + [-x for x in Sl_cgdl[1:]]

    # ____________________________________________Auxiliary flows (mol.s-1)_____________________________________________

    auxiliary_flows_dico = auxiliaries(t, sv, control_variables, i_fc, Jv_agc_agdl, Jv_cgdl_cgc, J_H2_agc_agdl,
                                       J_O2_cgdl_cgc, operating_inputs, parameters)

    # _____________________________________Assemble and return the flow dictionary______________________________________

    return {**auxiliary_flows_dico,
            'Jv': {**auxiliary_flows_dico.get('Jv', {}),
                   'agc_agdl': Jv_agc_agdl, 'agdl_agdl': Jv_agdl_agdl, 'agdl_ampl': Jv_agdl_ampl,
                   'ampl_ampl': Jv_ampl_ampl, 'ampl_acl': Jv_ampl_acl,
                   'ccl_cmpl': Jv_ccl_cmpl, 'cmpl_cmpl': Jv_cmpl_cmpl, 'cmpl_cgdl': Jv_cmpl_cgdl,
                   'cgdl_cgdl': Jv_cgdl_cgdl, 'cgdl_cgc': Jv_cgdl_cgc},
            'Jl': {'agc_agdl': Jl_agc_agdl, 'agdl_agdl': Jl_agdl_agdl, 'agdl_ampl': Jl_agdl_ampl,
                   'ampl_ampl': Jl_ampl_ampl, 'ampl_acl': Jl_ampl_acl,
                   'ccl_cmpl': Jl_ccl_cmpl, 'cmpl_cmpl': Jl_cmpl_cmpl, 'cmpl_cgdl': Jl_cmpl_cgdl,
                   'cgdl_cgdl': Jl_cgdl_cgdl, 'cgdl_cgc': Jl_cgdl_cgc},
            'J_lambda': {'acl_mem': J_lambda_acl_mem, 'mem_ccl': J_lambda_mem_ccl},
            'J_H2': {**auxiliary_flows_dico.get('J_H2', {}),
                     'agc_agdl': J_H2_agc_agdl, 'agdl_agdl': J_H2_agdl_agdl, 'agdl_ampl': J_H2_agdl_ampl,
                     'ampl_ampl': J_H2_ampl_ampl, 'ampl_acl': J_H2_ampl_acl},
            'J_O2': {**auxiliary_flows_dico.get('J_O2', {}),
                     'ccl_cmpl': J_O2_ccl_cmpl, 'cmpl_cmpl': J_O2_cmpl_cmpl, 'cmpl_cgdl': J_O2_cmpl_cgdl,
                     'cgdl_cgdl': J_O2_cgdl_cgdl, 'cgdl_cgc': J_O2_cgdl_cgc},
            'S_abs': {'acl': S_abs_acl, 'ccl': S_abs_ccl},
            'Sp' :{'acl': Sp_acl, 'ccl': Sp_ccl},
            'S_H2' :{'reac': S_H2_reac, 'cros': S_H2_cros}, 'S_O2': {'reac': S_O2_reac, 'cros': S_O2_cros},
            'Sv' : {'agdl': Sv_agdl, 'ampl': Sv_ampl, 'acl': Sv_acl, 'ccl': Sv_ccl, 'cmpl': Sv_cmpl, 'cgdl': Sv_cgdl},
            'Sl': {'agdl': Sl_agdl, 'ampl': Sl_ampl, 'acl': Sl_acl, 'ccl': Sl_ccl, 'cmpl': Sl_cmpl, 'cgdl': Sl_cgdl}}