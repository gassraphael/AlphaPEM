# -*- coding: utf-8 -*-

"""This file represents all the flows inside the fuel cell system. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import Kshape, rho_mem, M_eq, epsilon_cl, theta_c_gdl, gamma_cond, gamma_evap, F, R
from model.auxiliaries import auxiliaries
from modules.transitory_functions import lambda_eq, gamma_sorp, D, Svl, sigma, K0, Da_eff, Dc_eff, h_a, h_c, k_H2, k_O2
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
    C_N2 = sv['C_N2']

    # Extraction of the operating inputs and parameters
    Tfc, Pc_des = operating_inputs['Tfc'], operating_inputs['Pc_des']
    Hgdl, Hmem, Hcl = parameters['Hgdl'], parameters['Hmem'], parameters['Hcl']
    Hgc, Wgc = parameters['Hgc'], parameters['Wgc']
    epsilon_gdl, epsilon_c = parameters['epsilon_gdl'], parameters['epsilon_c']
    e, kappa_co, n_gdl = parameters['e'], parameters['kappa_co'], parameters['n_gdl']
    a_slim, b_slim = parameters['a_slim'], parameters['b_slim']

    # Intermediate values
    (Pagc, Pcgc, s_agdl_agdl, s_agdl_acl, s_ccl_cgdl, s_cgdl_cgdl, epsilon_mean, theta_c_mean, lambda_acl_mem,
     lambda_mem_ccl, Pagc_agdl, Pagdl_agdl, Pagdl_acl, Pccl_cgdl, Pcgdl_cgdl, Pcgdl_cgc, nu_l) \
        = flows_int_values(sv, operating_inputs, parameters)

    # Inlet and outlet flows
    Jv_a_in, Jv_a_out, Jv_c_in, Jv_c_out, J_H2_in, J_H2_out, J_O2_in, J_O2_out, J_N2_in, J_N2_out, \
        Wasm_in, Wasm_out, Waem_in, Waem_out, Wcsm_in, Wcsm_out, Wcem_in, Wcem_out, Ware, \
        Wv_asm_in, Wv_aem_out, Wv_csm_in, Wv_cem_out \
        = auxiliaries(t, sv, control_variables, i_fc, operating_inputs, parameters)

    # ________________________________________Dissolved water flows (mol.m-2.s-1)_______________________________________

    # Anode side
    J_lambda_mem_acl = 2.5 / 22 * i_fc / F * lambda_acl_mem - 2 * rho_mem / M_eq * \
               D(lambda_acl_mem) * (lambda_mem - lambda_acl) / (Hmem + Hcl)
    # Cathode side
    J_lambda_mem_ccl = 2.5 / 22 * i_fc / F * lambda_mem_ccl - 2 * rho_mem / M_eq * \
               D(lambda_mem_ccl) * (lambda_ccl - lambda_mem) / (Hmem + Hcl)

    # _________________________________________Liquid water flows (kg.m-2.s-1)__________________________________________

    # Anode side
    Jl_agdl_agdl = [None] + [0] * (n_gdl - 1)
    for i in range(1, n_gdl):
        Jl_agdl_agdl[i] = - sigma(Tfc) * K0(epsilon_gdl, epsilon_c, epsilon_gdl) / nu_l * abs(np.cos(theta_c_gdl)) * \
                          (epsilon_gdl / K0(epsilon_gdl, epsilon_c, epsilon_gdl)) ** 0.5 * \
                          (s_agdl_agdl[i] ** e + 1e-7) * (1.417 - 4.24 * s_agdl_agdl[i] + 3.789 * s_agdl_agdl[i] ** 2) * \
                          (sv[f's_agdl_{i + 1}'] - sv[f's_agdl_{i}']) / (Hgdl / n_gdl)
    Jl_agdl_acl = - 2 * sigma(Tfc) * K0(epsilon_mean, epsilon_c, epsilon_gdl) / nu_l * abs(np.cos(theta_c_mean)) * \
                  (epsilon_mean / K0(epsilon_mean, epsilon_c, epsilon_gdl)) ** 0.5 * \
                  (s_agdl_acl ** e + 1e-7) * (1.417 - 4.24 * s_agdl_acl + 3.789 * s_agdl_acl ** 2) * \
                  (s_acl - sv[f's_agdl_{n_gdl}']) / (Hgdl / n_gdl + Hcl)

    # Cathode side
    Jl_cgdl_cgdl = [None] + [0] * (n_gdl - 1)
    for i in range(1, n_gdl):
        Jl_cgdl_cgdl[i] = - sigma(Tfc) * K0(epsilon_gdl, epsilon_c, epsilon_gdl) / nu_l * abs(np.cos(theta_c_gdl)) * \
                          (epsilon_gdl / K0(epsilon_gdl, epsilon_c, epsilon_gdl)) ** 0.5 * \
                          (s_cgdl_cgdl[i] ** e + 1e-7) * (1.417 - 4.24 * s_cgdl_cgdl[i] + 3.789 * s_cgdl_cgdl[i] ** 2) * \
                          (sv[f's_cgdl_{i + 1}'] - sv[f's_cgdl_{i}']) / (Hgdl / n_gdl)
    Jl_ccl_cgdl = - 2 * sigma(Tfc) * K0(epsilon_mean, epsilon_c, epsilon_gdl) / nu_l * abs(np.cos(theta_c_mean)) * \
                  (epsilon_mean / K0(epsilon_mean, epsilon_c, epsilon_gdl)) ** 0.5 * \
                  (s_ccl_cgdl ** e + 1e-7) * (1.417 - 4.24 * s_ccl_cgdl + 3.789 * s_ccl_cgdl ** 2) * \
                  (sv['s_cgdl_1'] - s_ccl) / (Hgdl / n_gdl + Hcl)

    # _____________________________________________Vapor flows (mol.m-2.s-1)____________________________________________

    # Convective vapor flows
    #   Anode side
    Jv_agc_agdl = h_a(Pagc_agdl, Tfc, Wgc, Hgc) * (C_v_agc - sv['C_v_agdl_1'])
    #   Cathode side
    Jv_cgdl_cgc = h_c(Pcgdl_cgc, Tfc, Wgc, Hgc) * (sv[f'C_v_cgdl_{n_gdl}'] - C_v_cgc)

    # Conductive vapor flows
    #   Anode side
    Jv_agdl_agdl = [None] + [0] * (n_gdl - 1)
    for i in range(1, n_gdl):
        Jv_agdl_agdl[i] = - Da_eff(s_agdl_agdl[i], epsilon_gdl, Pagdl_agdl[i], Tfc, epsilon_c, epsilon_gdl) * \
                          (sv[f'C_v_agdl_{i + 1}'] - sv[f'C_v_agdl_{i}']) / (Hgdl / n_gdl)
    Jv_agdl_acl = - 2 * Da_eff(s_agdl_acl, epsilon_mean, Pagdl_acl, Tfc, epsilon_c, epsilon_gdl) * \
                  (C_v_acl - sv[f'C_v_agdl_{n_gdl}']) / (Hgdl / n_gdl + Hcl)

    #   Cathode side
    Jv_cgdl_cgdl = [None] + [0] * (n_gdl - 1)
    for i in range(1, n_gdl):
        Jv_cgdl_cgdl[i] = - Dc_eff(s_cgdl_cgdl[i], epsilon_gdl, Pcgdl_cgdl[i], Tfc, epsilon_c, epsilon_gdl) * \
                          (sv[f'C_v_cgdl_{i + 1}'] - sv[f'C_v_cgdl_{i}']) / (Hgdl / n_gdl)
    Jv_ccl_cgdl = - 2 * Dc_eff(s_ccl_cgdl, epsilon_mean, Pccl_cgdl, Tfc, epsilon_c, epsilon_gdl) * \
                  (sv['C_v_cgdl_1'] - C_v_ccl) / (Hgdl / n_gdl + Hcl)

    # __________________________________________H2 and O2 flows (mol.m-2.s-1)___________________________________________

    # Hydrogen and oxygen consumption
    #   Anode side
    S_H2_acl = - i_fc / (2 * F * Hcl) - \
               R * Tfc / (Hmem * Hcl) * (k_H2(lambda_mem, Tfc, kappa_co) * C_H2_acl +
                                         2 * k_O2(lambda_mem, Tfc, kappa_co) * C_O2_ccl)
    #   Cathode side
    S_O2_ccl = - i_fc / (4 * F * Hcl) - \
               R * Tfc / (Hmem * Hcl) * (k_O2(lambda_mem, Tfc, kappa_co) * C_O2_ccl +
                                         1 / 2 * k_H2(lambda_mem, Tfc, kappa_co) * C_H2_acl)

    # Conductive-convective H2 and O2 flows
    #   Anode side
    J_H2_agc_agdl = h_a(Pagc_agdl, Tfc, Wgc, Hgc) * (C_H2_agc - sv['C_H2_agdl_1'])
    #   Cathode side
    J_O2_cgdl_cgc = h_c(Pcgdl_cgc, Tfc, Wgc, Hgc) * (sv[f'C_O2_cgdl_{n_gdl}'] - C_O2_cgc)

    # Conductive H2 and O2 flows
    #   Anode side
    J_H2_agdl_agdl = [None] + [0] * (n_gdl - 1)
    for i in range(1, n_gdl):
        J_H2_agdl_agdl[i] = - Da_eff(s_agdl_agdl[i], epsilon_gdl, Pagdl_agdl[i], Tfc, epsilon_c, epsilon_gdl) * \
                            (sv[f'C_H2_agdl_{i + 1}'] - sv[f'C_H2_agdl_{i}']) / (Hgdl / n_gdl)
    J_H2_agdl_acl = - 2 * Da_eff(s_agdl_acl, epsilon_mean, Pagdl_acl, Tfc, epsilon_c, epsilon_gdl) * \
                    (C_H2_acl - sv[f'C_H2_agdl_{n_gdl}']) / (Hgdl / n_gdl + Hcl)

    #   Cathode side
    J_O2_cgdl_cgdl = [None] + [0] * (n_gdl - 1)
    for i in range(1, n_gdl):
        J_O2_cgdl_cgdl[i] = - Dc_eff(s_cgdl_cgdl[i], epsilon_gdl, Pcgdl_cgdl[i], Tfc, epsilon_c, epsilon_gdl) * \
                            (sv[f'C_O2_cgdl_{i + 1}'] - sv[f'C_O2_cgdl_{i}']) / (Hgdl / n_gdl)
    J_O2_ccl_cgdl = - 2 * Dc_eff(s_ccl_cgdl, epsilon_mean, Pccl_cgdl, Tfc, epsilon_c, epsilon_gdl) * \
                    (sv['C_O2_cgdl_1'] - C_O2_ccl) / (Hgdl / n_gdl + Hcl)

    # __________________________________________Water generated (mol.m-3.s-1)___________________________________________

    # Water produced in the membrane at the CL through the chemical reaction and crossover
    #   Anode side
    Sp_acl = 2 * k_O2(lambda_mem, Tfc, kappa_co) * R * Tfc / (Hmem * Hcl) * C_O2_ccl
    #   Cathode side
    Sp_ccl = i_fc / (2 * F * Hcl) + k_H2(lambda_mem, Tfc, kappa_co) * R * Tfc / (Hmem * Hcl) * C_H2_acl

    # Water sorption in the CL:
    #   Anode side
    S_sorp_acl = gamma_sorp(C_v_acl, s_acl, lambda_acl, Tfc, Hcl, Kshape) * rho_mem / M_eq * \
                 (lambda_eq(C_v_acl, s_acl, Tfc, Kshape) - lambda_acl)
    #   Cathode side
    S_sorp_ccl = gamma_sorp(C_v_ccl, s_ccl, lambda_ccl, Tfc, Hcl, Kshape) * rho_mem / M_eq * \
                 (lambda_eq(C_v_ccl, s_ccl, Tfc, Kshape) - lambda_ccl)

    # Liquid water generated through vapor condensation or degenerated through evaporation
    #   Anode side
    Sl_agdl = [None] + [Svl(sv[f's_agdl_{i}'], sv[f'C_v_agdl_{i}'], sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}'],
                            epsilon_gdl, Tfc, gamma_cond, gamma_evap) for i in range(1, n_gdl + 1)]
    Sl_acl = Svl(s_acl, C_v_acl, C_v_acl + C_H2_acl, epsilon_cl, Tfc, gamma_cond, gamma_evap)
    #   Cathode side
    Sl_cgdl = [None] + [Svl(sv[f's_cgdl_{i}'], sv[f'C_v_cgdl_{i}'], sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2,
                            epsilon_gdl, Tfc, gamma_cond, gamma_evap) for i in range(1, n_gdl + 1)]
    Sl_ccl = Svl(s_ccl, C_v_ccl, C_v_ccl + C_O2_ccl + C_N2, epsilon_cl, Tfc, gamma_cond, gamma_evap)

    # Vapor generated through liquid water evaporation or degenerated through condensation
    #   Anode side
    Sv_agdl = [None] + [-x for x in Sl_agdl[1:]]
    Sv_acl = - Sl_acl
    #   Cathode side
    Sv_cgdl = [None] + [-x for x in Sl_cgdl[1:]]
    Sv_ccl = - Sl_ccl

    return {'Jv_a_in': Jv_a_in, 'Jv_a_out': Jv_a_out, 'Jv_c_in': Jv_c_in, 'Jv_c_out': Jv_c_out, 'J_H2_in': J_H2_in,
            'J_H2_out': J_H2_out, 'J_O2_in': J_O2_in, 'J_O2_out': J_O2_out, 'J_N2_in': J_N2_in, 'J_N2_out': J_N2_out,
            'Jv_agc_agdl': Jv_agc_agdl, 'Jv_agdl_agdl': Jv_agdl_agdl, 'Jv_agdl_acl': Jv_agdl_acl,
            'S_sorp_acl': S_sorp_acl, 'S_sorp_ccl': S_sorp_ccl, 'Jv_ccl_cgdl': Jv_ccl_cgdl,
            'Jv_cgdl_cgdl': Jv_cgdl_cgdl, 'Jv_cgdl_cgc': Jv_cgdl_cgc, 'Jl_agdl_agdl': Jl_agdl_agdl,
            'Jl_agdl_acl': Jl_agdl_acl, 'J_lambda_mem_acl': J_lambda_mem_acl, 'J_lambda_mem_ccl': J_lambda_mem_ccl,
            'Jl_ccl_cgdl': Jl_ccl_cgdl, 'Jl_cgdl_cgdl': Jl_cgdl_cgdl, 'Sp_acl': Sp_acl, 'Sp_ccl': Sp_ccl,
            'J_H2_agc_agdl': J_H2_agc_agdl, 'J_H2_agdl_agdl': J_H2_agdl_agdl, 'J_H2_agdl_acl': J_H2_agdl_acl,
            'J_O2_ccl_cgdl': J_O2_ccl_cgdl, 'J_O2_cgdl_cgdl': J_O2_cgdl_cgdl, 'J_O2_cgdl_cgc': J_O2_cgdl_cgc,
            'S_H2_acl': S_H2_acl, 'S_O2_ccl': S_O2_ccl, 'Sv_agdl': Sv_agdl, 'Sv_acl': Sv_acl, 'Sv_ccl': Sv_ccl,
            'Sv_cgdl': Sv_cgdl, 'Sl_agdl': Sl_agdl, 'Sl_acl': Sl_acl, 'Sl_ccl': Sl_ccl, 'Sl_cgdl': Sl_cgdl,
            'Pagc': Pagc, 'Pcgc': Pcgc, 'Wasm_in': Wasm_in, 'Wasm_out': Wasm_out, 'Waem_in': Waem_in,
            'Waem_out': Waem_out, 'Wcsm_in': Wcsm_in, 'Wcsm_out': Wcsm_out, 'Wcem_in': Wcem_in, 'Wcem_out': Wcem_out,
            'Ware': Ware, 'Wv_asm_in': Wv_asm_in, 'Wv_aem_out': Wv_aem_out, 'Wv_csm_in': Wv_csm_in,
            'Wv_cem_out': Wv_cem_out}