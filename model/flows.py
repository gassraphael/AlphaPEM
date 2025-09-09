# -*- coding: utf-8 -*-

"""This file represents all the matter flows inside the fuel cell system. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import rho_mem, M_eq, F, R
from model.auxiliaries import auxiliaries
from modules.transitory_functions import Dcap, lambda_eq, gamma_sorp, Svl, k_H2, k_O2
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
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv['C_v_agc'], sv['C_v_acl'], sv['C_v_ccl'], sv['C_v_cgc']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv['C_H2_agc'], sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_O2_cgc']
    C_N2_a, C_N2_c = sv['C_N2_a'], sv['C_N2_c']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']

    # Extraction of the operating inputs and parameters
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    epsilon_gdl, epsilon_cl = parameters['epsilon_gdl'], parameters['epsilon_cl']
    epsilon_mpl, epsilon_c = parameters['epsilon_mpl'], parameters['epsilon_c']
    e, kappa_co, n_gdl, n_mpl = parameters['e'], parameters['kappa_co'], parameters['n_gdl'], parameters['n_mpl']

    # Intermediate values
    (H_gdl_node, H_mpl_node, Pagc, Pcgc, lambda_acl_mem, lambda_mem_ccl, D_acl_mem, D_mem_ccl, D_cap_agdl_agdl,
     D_cap_agdl_ampl, D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_cgdl_cgdl, D_cap_cmpl_cgdl, D_cap_cmpl_cmpl,
     D_cap_ccl_cmpl, ha_Da_eff_agc_agdl, hc_Dc_eff_cgdl_cgc, Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl,
     Da_eff_ampl_acl, Dc_eff_cgdl_cgdl, Dc_eff_cmpl_cgdl, Dc_eff_cmpl_cmpl, Dc_eff_ccl_cmpl, T_acl_mem_ccl) = \
        flows_int_values(sv, operating_inputs, parameters)

    # Inlet and outlet flows
    (Jv_a_in, Jv_a_out, Jv_c_in, Jv_c_out, J_H2_in, J_H2_out, J_O2_in, J_O2_out, J_N2_a_in, J_N2_a_out, J_N2_c_in,
     J_N2_c_out, Wasm_in, Wasm_out, Waem_in, Waem_out, Wcsm_in, Wcsm_out, Wcem_in, Wcem_out, Ware, Wv_asm_in,
     Wv_aem_out, Wv_csm_in, Wv_cem_out) = auxiliaries(t, sv, control_variables, i_fc, operating_inputs, parameters)

    # ________________________________________Dissolved water flows (mol.m-2.s-1)_______________________________________

    # Anode side
    J_lambda_acl_mem = 2.5 / 22 * i_fc / F * lambda_acl_mem - \
                       2 * rho_mem / M_eq * D_acl_mem * (lambda_mem - lambda_acl) / (Hmem + Hacl)
    # Cathode side
    J_lambda_mem_ccl = 2.5 / 22 * i_fc / F * lambda_mem_ccl - \
                       2 * rho_mem / M_eq * D_mem_ccl * (lambda_ccl - lambda_mem) / (Hmem + Hccl)

    # _________________________________________Liquid water flows (kg.m-2.s-1)__________________________________________

    # Anode side
    s_agc = 0  # Dirichlet boundary condition (taken at the agc/agdl border).
    Jl_agc_agdl = - 2 * Dcap('gdl', sv['s_agdl_1'], sv['T_agdl_1'], epsilon_gdl, e, epsilon_c=epsilon_c) * \
                    (sv['s_agdl_1'] - s_agc) / H_gdl_node
    Jl_agdl_agdl = [None] + [- D_cap_agdl_agdl[i] * (sv[f's_agdl_{i + 1}'] - sv[f's_agdl_{i}']) / H_gdl_node
                             for i in range(1, n_gdl)]
    Jl_agdl_ampl = - 2 * D_cap_agdl_ampl * (sv['s_ampl_1'] - sv[f's_agdl_{n_gdl}']) / (H_gdl_node + H_mpl_node)
    Jl_ampl_ampl = [None] + [- D_cap_ampl_ampl[i] * (sv[f's_ampl_{i + 1}'] - sv[f's_ampl_{i}']) / H_mpl_node
                             for i in range(1, n_mpl)]
    Jl_ampl_acl = - 2 * D_cap_ampl_acl * (s_acl - sv[f's_ampl_{n_mpl}']) / (H_mpl_node + Hacl)

    # Cathode side
    s_cgc = 0  # Dirichlet boundary condition (taken at the cgc/cgdl border).
    Jl_ccl_cmpl = - 2 * D_cap_ccl_cmpl * (sv['s_cmpl_1'] - s_ccl) / (Hccl + H_mpl_node)
    Jl_cmpl_cmpl = [None] + [- D_cap_cmpl_cmpl[i] * (sv[f's_cmpl_{i + 1}'] - sv[f's_cmpl_{i}']) / H_mpl_node
                            for i in range(1, n_mpl)]
    Jl_cmpl_cgdl = - 2 * D_cap_cmpl_cgdl * (sv['s_cgdl_1'] - sv[f's_cmpl_{n_mpl}']) / (H_mpl_node + H_gdl_node)
    Jl_cgdl_cgdl = [None] + [- D_cap_cgdl_cgdl[i] * (sv[f's_cgdl_{i + 1}'] - sv[f's_cgdl_{i}']) / H_gdl_node
                             for i in range(1, n_gdl)]
    Jl_cgdl_cgc = - 2 * Dcap('gdl', sv[f's_cgdl_{n_gdl}'], sv[f'T_cgdl_{n_gdl}'], epsilon_gdl, e,
                             epsilon_c=epsilon_c) * (s_cgc - sv[f's_cgdl_{n_gdl}']) / H_gdl_node

    # _____________________________________________Vapor flows (mol.m-2.s-1)____________________________________________

    # Convective vapor flows
    #   Anode side
    Jv_agc_agdl = - 2 * ha_Da_eff_agc_agdl * (sv['C_v_agdl_1'] - C_v_agc) / H_gdl_node
    #   Cathode side
    Jv_cgdl_cgc = - 2 * hc_Dc_eff_cgdl_cgc * (C_v_cgc - sv[f'C_v_cgdl_{n_gdl}']) / H_gdl_node

    # Conductive vapor flows
    #   Anode side
    Jv_agdl_agdl = [None] + [- Da_eff_agdl_agdl[i] * (sv[f'C_v_agdl_{i + 1}'] - sv[f'C_v_agdl_{i}']) / H_gdl_node
                             for i in range(1, n_gdl)]
    Jv_agdl_ampl = - 2 * Da_eff_agdl_ampl * (sv['C_v_ampl_1'] - sv[f'C_v_agdl_{n_gdl}']) / (H_gdl_node + H_mpl_node)
    Jv_ampl_ampl = [None] + [- Da_eff_ampl_ampl[i] * (sv[f'C_v_ampl_{i + 1}'] - sv[f'C_v_ampl_{i}']) / H_mpl_node
                             for i in range(1, n_mpl)]
    Jv_ampl_acl = - 2 * Da_eff_ampl_acl * (C_v_acl - sv[f'C_v_ampl_{n_mpl}']) / (H_mpl_node + Hacl)

    #   Cathode side
    Jv_ccl_cmpl = - 2 * Dc_eff_ccl_cmpl * (sv['C_v_cmpl_1'] - C_v_ccl) / (Hccl + H_mpl_node)
    Jv_cmpl_cmpl = [None] + [- Dc_eff_cmpl_cmpl[i] * (sv[f'C_v_cmpl_{i + 1}'] - sv[f'C_v_cmpl_{i}']) / H_mpl_node
                             for i in range(1, n_mpl)]
    Jv_cmpl_cgdl = - 2 * Dc_eff_cmpl_cgdl * (sv['C_v_cgdl_1'] - sv[f'C_v_cmpl_{n_mpl}']) / (H_mpl_node + H_gdl_node)
    Jv_cgdl_cgdl = [None] + [- Dc_eff_cgdl_cgdl[i] * (sv[f'C_v_cgdl_{i + 1}'] - sv[f'C_v_cgdl_{i}']) / H_gdl_node
                             for i in range(1, n_gdl)]

    # __________________________________________H2 and O2 flows (mol.m-2.s-1)___________________________________________

    # Hydrogen and oxygen consumption
    #   Anode side
    S_H2_acl = - i_fc / (2 * F * Hacl) - \
               R * T_acl_mem_ccl / (Hmem * Hacl) * (k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl +
                                         2 * k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl)
    #   Cathode side
    S_O2_ccl = - i_fc / (4 * F * Hccl) - \
               R * T_acl_mem_ccl / (Hmem * Hccl) * (k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl +
                                         1 / 2 * k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl)

    # Conductive-convective H2 and O2 flows
    #   Anode side
    J_H2_agc_agdl = - 2 * ha_Da_eff_agc_agdl * (sv['C_H2_agdl_1'] - C_H2_agc) / H_gdl_node
    #   Cathode side
    J_O2_cgdl_cgc = - 2 * hc_Dc_eff_cgdl_cgc * (C_O2_cgc - sv[f'C_O2_cgdl_{n_gdl}']) / H_gdl_node

    # Conductive H2 and O2 flows
    #   Anode side
    J_H2_agdl_agdl = [None] + [- Da_eff_agdl_agdl[i] * (sv[f'C_H2_agdl_{i+1}'] - sv[f'C_H2_agdl_{i}']) / H_gdl_node
                               for i in range(1, n_gdl)]
    J_H2_agdl_ampl = - 2 * Da_eff_agdl_ampl * (sv['C_H2_ampl_1'] - sv[f'C_H2_agdl_{n_gdl}']) / (H_gdl_node + H_mpl_node)
    J_H2_ampl_ampl = [None] + [- Da_eff_ampl_ampl[i] * (sv[f'C_H2_ampl_{i + 1}'] - sv[f'C_H2_ampl_{i}']) / H_mpl_node
                               for i in range(1, n_mpl)]
    J_H2_ampl_acl = - 2 * Da_eff_ampl_acl * (C_H2_acl - sv[f'C_H2_ampl_{n_mpl}']) / (H_mpl_node + Hacl)

    #   Cathode side
    J_O2_ccl_cmpl = - 2 * Dc_eff_ccl_cmpl * (sv['C_O2_cmpl_1'] - C_O2_ccl) / (Hccl + H_mpl_node)
    J_O2_cmpl_cmpl = [None] + [- Dc_eff_cmpl_cmpl[i] * (sv[f'C_O2_cmpl_{i + 1}'] - sv[f'C_O2_cmpl_{i}']) / H_mpl_node
                               for i in range(1, n_mpl)]
    J_O2_cmpl_cgdl = - 2 * Dc_eff_cmpl_cgdl * (sv['C_O2_cgdl_1'] - sv[f'C_O2_cmpl_{n_mpl}']) / (H_mpl_node + H_gdl_node)
    J_O2_cgdl_cgdl = [None] + [- Dc_eff_cgdl_cgdl[i] * (sv[f'C_O2_cgdl_{i + 1}'] - sv[f'C_O2_cgdl_{i}']) / H_gdl_node
                               for i in range(1, n_gdl)]

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
                            Ctot=sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}'] + C_N2_a,
                            T=sv[f'T_agdl_{i}'], epsilon=epsilon_gdl) for i in range(1, n_gdl + 1)]
    Sl_ampl = [None] + [Svl(element='anode', s=sv[f's_ampl_{i}'], C_v=sv[f'C_v_ampl_{i}'],
                            Ctot=sv[f'C_v_ampl_{i}'] + sv[f'C_H2_ampl_{i}'] + C_N2_a,
                            T=sv[f'T_ampl_{i}'], epsilon=epsilon_mpl) for i in range(1, n_mpl + 1)]
    Sl_acl = Svl(element='anode', s=s_acl, C_v=C_v_acl, Ctot=C_v_acl + C_H2_acl + C_N2_a, T=T_acl, epsilon=epsilon_cl)
    #   Cathode side
    Sl_ccl = Svl(element='cathode', s=s_ccl, C_v=C_v_ccl, Ctot=C_v_ccl + C_O2_ccl + C_N2_c, T=T_ccl, epsilon=epsilon_cl)
    Sl_cmpl = [None] + [Svl(element='cathode', s=sv[f's_cmpl_{i}'], C_v=sv[f'C_v_cmpl_{i}'],
                            Ctot=sv[f'C_v_cmpl_{i}'] + sv[f'C_O2_cmpl_{i}'] + C_N2_c,
                            T=sv[f'T_cmpl_{i}'], epsilon=epsilon_mpl) for i in range(1, n_mpl + 1)]
    Sl_cgdl = [None] + [Svl(element='cathode', s=sv[f's_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                            Ctot=sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2_c,
                            T=sv[f'T_cgdl_{i}'], epsilon=epsilon_gdl) for i in range(1, n_gdl + 1)]

    # Vapor generated through liquid water evaporation or degenerated through condensation
    #   Anode side
    Sv_agdl = [None] + [-x for x in Sl_agdl[1:]]
    Sv_ampl = [None] + [-x for x in Sl_ampl[1:]]
    Sv_acl = - Sl_acl
    #   Cathode side
    Sv_ccl = - Sl_ccl
    Sv_cmpl = [None] + [-x for x in Sl_cmpl[1:]]
    Sv_cgdl = [None] + [-x for x in Sl_cgdl[1:]]

    return {'Jv_a_in': Jv_a_in, 'Jv_a_out': Jv_a_out, 'Jv_c_in': Jv_c_in, 'Jv_c_out': Jv_c_out, 'J_H2_in': J_H2_in,
            'J_H2_out': J_H2_out, 'J_O2_in': J_O2_in, 'J_O2_out': J_O2_out, 'J_N2_a_in': J_N2_a_in,
            'J_N2_a_out': J_N2_a_out, 'J_N2_c_in': J_N2_c_in, 'J_N2_c_out': J_N2_c_out, 'Jv_agc_agdl': Jv_agc_agdl,
            'Jv_agdl_agdl': Jv_agdl_agdl, 'Jv_agdl_ampl': Jv_agdl_ampl, 'Jv_ampl_ampl': Jv_ampl_ampl,
            'Jv_ampl_acl': Jv_ampl_acl, 'S_abs_acl': S_abs_acl, 'S_abs_ccl': S_abs_ccl, 'Jv_ccl_cmpl': Jv_ccl_cmpl,
            'Jv_cmpl_cmpl': Jv_cmpl_cmpl, 'Jv_cmpl_cgdl': Jv_cmpl_cgdl, 'Jv_cgdl_cgdl': Jv_cgdl_cgdl,
            'Jv_cgdl_cgc': Jv_cgdl_cgc, 'Jl_agc_agdl': Jl_agc_agdl, 'Jl_agdl_agdl': Jl_agdl_agdl,
            'Jl_agdl_ampl': Jl_agdl_ampl, 'Jl_ampl_ampl': Jl_ampl_ampl, 'Jl_ampl_acl': Jl_ampl_acl,
            'J_lambda_acl_mem': J_lambda_acl_mem, 'J_lambda_mem_ccl': J_lambda_mem_ccl, 'Jl_ccl_cmpl': Jl_ccl_cmpl,
            'Jl_cmpl_cmpl': Jl_cmpl_cmpl, 'Jl_cmpl_cgdl': Jl_cmpl_cgdl, 'Jl_cgdl_cgdl': Jl_cgdl_cgdl,
            'Jl_cgdl_cgc': Jl_cgdl_cgc, 'Sp_acl': Sp_acl, 'Sp_ccl': Sp_ccl, 'J_H2_agc_agdl': J_H2_agc_agdl,
            'J_H2_agdl_agdl': J_H2_agdl_agdl, 'J_H2_agdl_ampl': J_H2_agdl_ampl, 'J_H2_ampl_ampl': J_H2_ampl_ampl,
            'J_H2_ampl_acl': J_H2_ampl_acl, 'J_O2_ccl_cmpl': J_O2_ccl_cmpl, 'J_O2_cmpl_cmpl': J_O2_cmpl_cmpl,
            'J_O2_cmpl_cgdl': J_O2_cmpl_cgdl, 'J_O2_cgdl_cgdl': J_O2_cgdl_cgdl, 'J_O2_cgdl_cgc': J_O2_cgdl_cgc,
            'S_H2_acl': S_H2_acl, 'S_O2_ccl': S_O2_ccl, 'Sv_agdl': Sv_agdl, 'Sv_ampl': Sv_ampl, 'Sv_acl': Sv_acl,
            'Sv_ccl': Sv_ccl, 'Sv_cmpl': Sv_cmpl, 'Sv_cgdl': Sv_cgdl, 'Sl_agdl': Sl_agdl, 'Sl_ampl': Sl_ampl,
            'Sl_acl': Sl_acl, 'Sl_ccl': Sl_ccl, 'Sl_cmpl': Sl_cmpl, 'Sl_cgdl': Sl_cgdl, 'Pagc': Pagc, 'Pcgc': Pcgc,
            'Wasm_in': Wasm_in, 'Wasm_out': Wasm_out, 'Waem_in': Waem_in, 'Waem_out': Waem_out, 'Wcsm_in': Wcsm_in,
            'Wcsm_out': Wcsm_out, 'Wcem_in': Wcem_in, 'Wcem_out': Wcem_out, 'Ware': Ware, 'Wv_asm_in': Wv_asm_in,
            'Wv_aem_out': Wv_aem_out, 'Wv_csm_in': Wv_csm_in, 'Wv_cem_out': Wv_cem_out}