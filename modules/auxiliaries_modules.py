# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the auxiliaries flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import Text, Pext, Phi_ext, M_H2, M_O2, M_N2, M_H2O, y_O2_ext, R, F
from modules.transitory_functions import average, interpolate, Psat, C_v_sat, k_H2, k_O2


# _________________________________________________Auxiliaries modules__________________________________________________

def auxiliaries_int_values_which_are_commun_with_dif_eq(t, sv, operating_inputs, parameters):
    """This functions calculates intermediate values for the auxiliaries flows calculation.

    Parameters
    ----------
    t : float
        Time (s).
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    operating_inputs : dict
        Operating inputs of the fuel cell model.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    k_purge : float
        Purge coefficient. It is equal to 1 if the purge is active and 0 otherwise.
    Abp_a : float
        Area of the back pressure valve in the anode external manifold (m²).
    Abp_c : float
        Area of the back pressure valve in the cathode external manifold (m²).
    """

    # Extraction of the variables
    lambda_mem = sv['lambda_mem']
    C_H2_acl, C_O2_ccl, C_N2_a, C_N2_c =  sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_N2_a'], sv['C_N2_c']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']
    Pasm_in, Pasm_in_re, Pasm, Pasm_out = sv.get('Pasm_in', None), sv.get('Pasm_in_re', None), sv.get('Pasm', None), sv['Pasm_out']
    Paem, Paem_in, Paem_out_re, Paem_out = sv.get('Paem', None), sv['Paem_in'], sv.get('Paem_out_re', None), sv.get('Paem_out', None)
    Pcsm_in, Pcsm, Pcsm_out = sv.get('Pcsm_in', None), sv.get('Pcsm', None), sv['Pcsm_out']
    Pcem_in, Pcem, Pcem_out = sv['Pcem_in'], sv.get('Pcem', None), sv.get('Pcem_out', None)
    Phi_asm_in, Phi_asm_in_re, Phi_asm, Phi_asm_out = sv.get('Phi_asm_in', None), sv.get('Phi_asm_in_re', None), sv.get('Phi_asm', None), sv['Phi_asm_out']
    Phi_aem_in, Phi_aem, Phi_aem_out_re, Phi_aem_out = sv['Phi_aem_in'], sv.get('Phi_aem', None), sv.get('Phi_aem_out_re', None), sv.get('Phi_aem_out', None)
    Phi_csm_in, Phi_csm, Phi_csm_out = sv.get('Phi_csm_in', None), sv.get('Phi_csm', None), sv['Phi_csm_out']
    Phi_cem_in, Phi_cem, Phi_cem_out = sv['Phi_cem_in'], sv.get('Phi_cem', None), sv.get('Phi_cem_out', None)
    Abp_a, Abp_c = sv.get('Abp_a', None), sv.get('Abp_c', None)
    # Extraction of the operating inputs and the parameters
    T_des, y_H2_in, Pa_des = operating_inputs['T_des'], operating_inputs['y_H2_in'], operating_inputs['Pa_des']
    Phi_a_des, Phi_c_des = operating_inputs['Phi_a_des'], operating_inputs['Phi_c_des']
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    A_T_a, A_T_c = parameters['A_T_a'], parameters['A_T_c']
    kappa_co, n_gc, t_purge, type_purge = parameters['kappa_co'], parameters['n_gc'], parameters['t_purge'], parameters['type_purge']

    # Physical quantities outside the stack
    #       Molar masses
    M = {}
    M['ext'] = Phi_ext * Psat(Text) / Pext * M_H2O + \
               y_O2_ext * (1 - Phi_ext * Psat(Text) / Pext) * M_O2 + \
               (1 - y_O2_ext) * (1 - Phi_ext * Psat(Text) / Pext) * M_N2
    M['H2_N2_in'] = y_H2_in * M_H2 + (1 - y_H2_in) * M_N2

    # Physical quantities inside the stack
    #       Pressures
    P = {}
    for i in range(1, n_gc + 1):
        P[f'agc_{i}'] = (sv[f'C_v_agc_{i}'] + sv[f'C_H2_agc_{i}'] + C_N2_a) * R * sv[f'T_agc_{i}']
        P[f'cgc_{i}'] = (sv[f'C_v_cgc_{i}'] + sv[f'C_O2_cgc_{i}'] + C_N2_c) * R * sv[f'T_cgc_{i}']
    #       Humidities
    Phi = {}
    for i in range(1, n_gc + 1):
        Phi[f'agc_{i}'] = sv[f'C_v_agc_{i}'] / C_v_sat(sv[f'T_agc_{i}'])
        Phi[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] / C_v_sat(sv[f'T_cgc_{i}'])
    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_O2 = {}
    y_H2 = {}
    for i in range(1, n_gc + 1):
        y_H2[f'agc_{i}'] = sv[f'C_H2_agc_{i}'] / (sv[f'C_H2_agc_{i}'] + C_N2_a)
        y_O2[f'cgc_{i}'] = sv[f'C_O2_cgc_{i}'] / (sv[f'C_O2_cgc_{i}'] + C_N2_c)
    #       Molar masses
    for i in range(1, n_gc + 1):
        M[f'agc_{i}'] = sv[f'C_v_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2O + \
                        sv[f'C_H2_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2 + \
                        C_N2_a * R * T_des / P[f'agc_{i}'] * M_N2
        M[f'cgc_{i}'] = Phi[f'cgc_{i}'] * Psat(T_des) / P[f'cgc_{i}'] * M_H2O + \
                        y_O2[f'cgc_{i}'] * (1 - Phi[f'cgc_{i}'] * Psat(T_des) / P[f'cgc_{i}']) * M_O2 + \
                        (1 - y_O2[f'cgc_{i}']) * (1 - Phi[f'cgc_{i}'] * Psat(T_des) / P[f'cgc_{i}']) * M_N2
    #       Density of the gas mixture.
    rho = {}
    for i in range(1, n_gc + 1):
        rho[f'agc_{i}'] = P[f'agc_{i}'] / (R * sv[f'T_agc_{i}']) * M[f'agc_{i}']
        rho[f'cgc_{i}'] = P[f'cgc_{i}'] / (R * sv[f'T_cgc_{i}']) * M[f'cgc_{i}']
    #       Internal current density
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                            weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl),
                                     Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
    i_n = i_H2 + i_O2

    # Physical quantities in the auxiliary system
    if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation" or \
       parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
        # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the EM
        y_H2_asm_out = (Pasm_out - Phi_asm_out * Psat(T_des) - C_N2_a * R * T_des) / (Pasm_out - Phi_asm_out * Psat(T_des))
        y_H2_aem_in = (Paem_in - Phi_aem_in * Psat(T_des) - C_N2_a * R * T_des) / (Paem_in - Phi_aem_in * Psat(T_des))
        y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_a * R * T_des) / (Paem - Phi_aem * Psat(T_des))
        y_H2_aem_out = (Paem_out - Phi_aem_out * Psat(T_des) - C_N2_a * R * T_des) / (Paem_out - Phi_aem_out * Psat(T_des))
        y_O2_csm_out = (Pcsm_out - Phi_csm_out * Psat(T_des) - C_N2_c * R * T_des) / (Pcsm_out - Phi_csm_out * Psat(T_des))
        y_O2_cem_in = (Pcem_in - Phi_cem_in * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem_in - Phi_cem_in * Psat(T_cgc))
        y_O2_cem = (Pcem - Phi_cem * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem - Phi_cem * Psat(T_cgc))
        y_O2_cem_out = (Pcem_out - Phi_cem_out * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem_out - Phi_cem_out * Psat(T_cgc))

        # Molar masses
        if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation":
            Masm_in = Phi_asm_in * Psat(T_des) / Pasm_in * M_H2O + \
                      (1 - Phi_asm_in * Psat(T_des) / Pasm_in) * M_H2
            Masm_in_re = Phi_asm_in_re * Psat(T_des) / Pasm_in_re * M_H2O + \
            (1 - Phi_asm_in_re * Psat(T_des) / Pasm_in_re) * M_H2
            Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
            Masm_out = Phi_asm_out * Psat(T_des) / Pasm_out * M_H2O + \
                       (1 - Phi_asm_out * Psat(T_des) / Pasm_out) * M_H2
            Maem_in = Phi_aem_in * Psat(T_des) / Paem_in * M_H2O + \
                      (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_H2
            Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
            Maem_out = Phi_aem_out * Psat(T_des) / Paem_out * M_H2O + \
                       (1 - Phi_aem_out * Psat(T_des) / Paem_out) * M_H2
            Maem_out_re = Phi_aem_out_re * Psat(T_des) / Paem_out_re * M_H2O + \
                          (1 - Phi_aem_out_re * Psat(T_des) / Paem_out_re) * M_H2
        else:  # parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
            Masm_in = Phi_asm_in * Psat(T_des) / Pasm_in * M_H2O + \
                      y_H2_in * (1 - Phi_asm_in * Psat(T_des) / Pasm_in) * M_H2 + \
                      (1 - y_H2_in) * (1 - Phi_asm_in * Psat(T_des) / Pasm_in) * M_N2
            Masm_in_re = None
            Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 + \
                   (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
            Masm_out = Phi_asm_out * Psat(T_des) / Pasm_out * M_H2O + \
                       y_H2_in * (1 - Phi_asm_out * Psat(T_des) / Pasm_out) * M_H2 + \
                       (1 - y_H2_in) * (1 - Phi_asm_out * Psat(T_des) / Pasm_out) * M_N2
            Maem_in = Phi_aem_in * Psat(T_des) / Paem_in * M_H2O + \
                      y_H2_aem_in * (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_H2 + \
                      (1 - y_H2_aem_in) * (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_N2
            Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 + \
                   (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
            Maem_out = Phi_aem_out * Psat(T_des) / Paem_out * M_H2O + \
                       y_H2_aem_out * (1 - Phi_aem_out * Psat(T_des) / Paem_out) * M_H2 + \
                       (1 - y_H2_aem_out) * (1 - Phi_aem_out * Psat(T_des) / Paem_out) * M_N2
            Maem_out_re = None
        # Molar masses at the cathode side
        Mcsm_in = Phi_csm_in * Psat(T_des) / Pcsm_in * M_H2O + \
                  y_O2_ext * (1 - Phi_csm_in * Psat(T_des) / Pcsm_in) * M_O2 + \
                  (1 - y_O2_ext) * (1 - Phi_csm_in * Psat(T_des) / Pcsm_in) * M_N2
        Mcsm = Phi_csm * Psat(T_des) / Pcsm * M_H2O + \
               y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 + \
               (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
        Mcsm_out = Phi_csm_out * Psat(T_des) / Pcsm_out * M_H2O + \
                   y_O2_ext * (1 - Phi_csm_out * Psat(T_des) / Pcsm_out) * M_O2 + \
                   (1 - y_O2_ext) * (1 - Phi_csm_out * Psat(T_des) / Pcsm_out) * M_N2
        Mcem_in = Phi_cem_in * Psat(T_cgc) / Pcem_in * M_H2O + \
                  y_O2_cem_in * (1 - Phi_cem_in * Psat(T_cgc) / Pcem_in) * M_O2 + \
                  (1 - y_O2_cem_in) * (1 - Phi_cem_in * Psat(T_cgc) / Pcem_in) * M_N2
        Mcem = Phi_cem * Psat(T_des) / Pcem * M_H2O + \
               y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 + \
               (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2
        Mcem_out = Phi_cem_out * Psat(T_cgc) / Pcem_out * M_H2O + \
                   y_O2_cem_out * (1 - Phi_cem_out * Psat(T_cgc) / Pcem_out) * M_O2 + \
                   (1 - y_O2_cem_out) * (1 - Phi_cem_out * Psat(T_cgc) / Pcem_out) * M_N2

        # Density of the gas mixture.
        if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation":
            rho_asm_in_re = Pasm_in_re / (R * T_des) * Masm_in_re
            rho_aem_out_re = Paem_out_re / (R * T_des) * Maem_out_re
        else: # parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
            rho_asm_in_re, rho_aem_out_re = [None] * 2
        rho_a_in = Pa_des / (R * T_des) * M_H2
        rho_asm_in = Pasm_in / (R * T_des) * Masm_in
        rho_asm = Pasm / (R * T_des) * Masm
        rho_asm_out = Pasm_out / (R * T_des) * Masm_out
        rho_aem_in = Paem_in / (R * T_des) * Maem_in
        rho_aem = Paem / (R * T_des) * Maem
        rho_aem_out = Paem_out / (R * T_des) * Maem_out
        rho_csm_in = Pcsm_in / (R * T_des) * Mcsm_in
        rho_csm = Pcsm / (R * T_des) * Mcsm
        rho_csm_out = Pcsm_out / (R * T_des) * Mcsm_out
        rho_cgc = Pcgc / (R * T_cgc) * Mcgc
        rho_cem_in = Pcem_in / (R * T_cgc) * Mcem_in
        rho_cem = Pcem / (R * T_cgc) * Mcem
        rho_cem_out = Pcem_out / (R * T_cgc) * Mcem_out

        # Purge
        if type_purge == "no_purge":
            k_purge = 0
        elif type_purge == "constant_purge":
            k_purge = 1
        elif type_purge == "periodic_purge":
            purge_time, delta_purge = t_purge
            if (t - int(t / (purge_time + delta_purge)) * (purge_time + delta_purge)) <= purge_time:
                k_purge = 1
            else:
                k_purge = 0
        else:
            raise ValueError("The type_purge variable should be correctly referenced.")
        # Back pressure valve area
        if Abp_a > A_T_a:
            Abp_a = A_T_a
        elif Abp_a < 0:
            Abp_a = 0
        if Abp_c > A_T_c:
            Abp_c = A_T_c
        elif Abp_c < 0:
            Abp_c = 0

    else:  # parameters["type_auxiliary"] == "no_auxiliary"
        y_O2['csm_out'] = y_O2_ext
        y_O2['cem_in'] = y_O2[f'cgc_{n_gc}']
        M['aem_in'] = Phi_aem_in * Psat(T_des) / Paem_in * M_H2O + \
                  (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_H2
        M['cem_in'] = Phi_cem_in * Psat(sv[f'T_cgc_{i}']) / Pcem_in * M_H2O + \
                      y_O2['cem_in'] * (1 - Phi_cem_in * Psat(sv[f'T_cgc_{i}']) / Pcem_in) * M_O2 + \
                      (1 - y_O2['cem_in']) * (1 - Phi_cem_in * Psat(sv[f'T_cgc_{i}']) / Pcem_in) * M_N2
        M['asm_out'] = Phi_a_des * Psat(T_des) / Pasm_out * M_H2O + \
                   (1 - Phi_a_des * Psat(T_des) / Pasm_out) * M_H2
        M['csm_out'] = Phi_c_des * Psat(T_des) / Pcsm_out * M_H2O + \
                   y_O2['csm_out'] * (1 - Phi_c_des * Psat(T_des) / Pcsm_out) * M_O2 + \
                   (1 - y_O2['csm_out']) * (1 - Phi_c_des * Psat(T_des) / Pcsm_out) * M_N2
        rho['asm_out'] = Pasm_out / (R * T_des) * M['asm_out']
        rho['aem_in'] = Paem_in / (R * T_des) * M['aem_in']
        rho['csm_out'] = Pcsm_out / (R * T_des) * M['csm_out']
        rho['cem_in'] = Pcem_in / (R * T_des) * M['cem_in']
        k_purge, Abp_a, Abp_c = [None] * 3

    return (i_n, P, Phi, y_H2, y_O2, M, rho, k_purge, Abp_a, Abp_c)


def auxiliaries_int_values(sv, parameters, P, Phi, y_O2):

    # Extraction of the variables
    v_asm_in, v_asm_in_re, v_asm, v_asm_out = sv.get('v_asm_in', None), sv.get('v_asm_in_re', None), sv.get('v_asm', None), sv['v_asm_out']
    v_aem_in, v_aem = sv['v_aem_in'], sv.get('v_aem', None)
    v_aem_out, v_aem_out_re, v_a_ext = sv.get('v_aem_out', None), sv.get('v_aem_out_re', None), sv['v_a_ext']
    v_csm_in, v_csm, v_csm_out = sv.get('v_csm_in', None), sv.get('v_csm', None), sv['v_csm_out']
    v_cem_in, v_cem, v_cem_out, v_c_ext = sv['v_cem_in'], sv.get('v_cem', None), sv.get('v_cem_out', None), sv['v_c_ext']
    Pasm_in, Pasm_in_re, Pasm, Pasm_out = sv.get('Pasm_in', None), sv.get('Pasm_in_re', None), sv.get('Pasm', None), sv['Pasm_out']
    Paem_in, Paem_out, Paem, Paem_out_re, Pa_ext = sv['Paem_in'], sv.get('Paem_out', None), sv.get('Paem', None), sv.get('Paem_out_re', None), sv['Pa_ext']
    Pcsm_in, Pcsm, Pcsm_out = sv.get('Pcsm_in', None), sv.get('Pcsm', None), sv['Pcsm_out']
    Pcem_in, Pcem, Pcem_out, Pc_ext = sv['Pcem_in'], sv.get('Pcem', None), sv.get('Pcem_out', None), sv['Pc_ext']
    Phi_asm_in, Phi_asm_in_re, Phi_asm = sv.get('Phi_asm_in', None), sv.get('Phi_asm_in_re', None), sv.get('Phi_asm', None)
    Phi_asm_out, Phi_aem_in, Phi_aem = sv['Phi_asm_out'], sv['Phi_aem_in'], sv.get('Phi_aem', None)
    Phi_aem_out, Phi_aem_out_re, Phi_a_ext = sv.get('Phi_aem_out', None), sv.get('Phi_aem_out_re', None), sv['Phi_a_ext']
    Phi_csm_in, Phi_csm, Phi_csm_out = sv.get('Phi_csm_in', None), sv.get('Phi_csm', None), sv['Phi_csm_out']
    Phi_cem_in, Phi_cem, Phi_cem_out, Phi_c_ext = sv['Phi_cem_in'], sv.get('Phi_cem', None), sv.get('Phi_cem_out', None), sv['Phi_c_ext']

    # Extraction of the operating inputs and the parameters
    Wagc, Wcgc, Lgc = parameters['Wagc'], parameters['Wcgc'], parameters['Lgc']
    Lm, L_endplate, L_man_gc = parameters['Lm'], parameters['L_endplate'], parameters['L_man_gc']
    n_gc, type_auxiliary = parameters['n_gc'], parameters['type_auxiliary']

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
         type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # Mean velocities in the different sections of the auxiliary in m.s-1
        v_asm_in_to_asm = interpolate([v_asm_in, v_asm], [L_endplate, Lm])
        v_asm_to_asm_out = interpolate([v_asm, v_asm_out], [Lm, L_man_gc])
        v_asm_out_to_agc = interpolate([v_asm_out, v_agc], [L_man_gc, Lgc / n_gc])
        v_agc_to_aem_in = interpolate([v_agc, v_aem_in], [Lgc / n_gc, L_man_gc])
        v_aem_in_to_aem = interpolate([v_aem_in, v_aem], [L_man_gc, Lm])
        v_asm_in_re_to_asm = interpolate([v_asm_in_re, v_asm], [L_endplate, Lm])
        v_aem_to_aem_out_re = interpolate([v_aem, v_aem_out_re], [Lm, L_endplate])
        v_aem_to_aem_out = interpolate([v_aem, v_aem_out], [Lm, L_endplate])
        v_aem_out_to_ext = None
        v_csm_in_to_csm = interpolate([v_csm_in, v_csm], [L_endplate, Lm])
        v_csm_to_csm_out = interpolate([v_csm, v_csm_out], [Lm, L_man_gc])
        v_csm_out_to_cgc = interpolate([v_csm_out, v_cgc], [L_man_gc, Lgc / n_gc])
        v_cgc_to_cem_in = interpolate([v_cgc, v_cem_in], [Lgc / n_gc, L_man_gc])
        v_cem_in_to_cem = interpolate([v_cem_in, v_cem], [L_man_gc, Lm])
        v_cem_to_cem_out = interpolate([v_cem, v_cem_out], [Lm, L_endplate])
        v_cem_out_to_ext = None
        # Mean pressures in the different sections of the auxiliary in Pa
        Pasm_in_to_asm = interpolate([Pasm_in, Pasm], [L_endplate, Lm])
        Pasm_to_asm_out = interpolate([Pasm, Pasm_out], [Lm, L_man_gc])
        Pasm_out_to_agc = interpolate([Pasm_out, Pagc], [L_man_gc, Lgc / n_gc])
        Pagc_to_aem_in = interpolate([Pagc, Paem_in], [Lgc / n_gc, L_man_gc])
        Paem_out_to_ext = None
        Paem_in_to_aem = interpolate([Paem_in, Paem], [L_man_gc, Lm])
        Pasm_in_re_to_asm = interpolate([Pasm_in_re, Pasm], [L_endplate, Lm])
        Paem_to_aem_out_re = interpolate([Paem, Paem_out_re], [Lm, L_endplate])
        Paem_to_aem_out = interpolate([Paem, Paem_out], [Lm, L_endplate])
        Pcsm_in_to_csm = interpolate([Pcsm_in, Pcsm], [L_endplate, Lm])
        Pcsm_to_csm_out = interpolate([Pcsm, Pcsm_out], [Lm, L_man_gc])
        Pcsm_out_to_cgc = interpolate([Pcsm_out, Pcgc], [L_man_gc, Lgc / n_gc])
        Pcgc_to_cem_in = interpolate([Pcgc, Pcem_in], [Lgc / n_gc, L_man_gc])
        Pcem_in_to_cem = interpolate([Pcem_in, Pcem], [L_man_gc, Lm])
        Pcem_to_cem_out = interpolate([Pcem, Pcem_out], [Lm, L_endplate])
        Pcem_out_to_ext = None
        # Mean humidities in the different sections of the auxiliary
        Phi_asm_in_to_asm = interpolate([Phi_asm_in, Phi_asm], [L_endplate, Lm])
        Phi_asm_to_asm_out = interpolate([Phi_asm, Phi_asm_out], [Lm, L_endplate])
        Phi_asm_out_to_agc = interpolate([Phi_asm_out, Phi_agc], [L_endplate, Lgc / n_gc])
        Phi_agc_to_aem_in = interpolate([Phi_agc, Phi_aem_in], [Lgc / n_gc, L_endplate])
        Phi_aem_in_to_aem = interpolate([Phi_aem_in, Phi_aem], [L_endplate, Lm])
        Phi_asm_in_re_to_asm = interpolate([Phi_asm_in_re, Phi_asm], [L_endplate, Lm])
        Phi_aem_to_aem_out_re = interpolate([Phi_aem, Phi_aem_out_re], [Lm, L_endplate])
        Phi_aem_to_aem_out = interpolate([Phi_aem, Phi_aem_out], [Lm, L_endplate])
        Phi_aem_out_to_ext = None
        Phi_csm_in_to_csm = interpolate([Phi_csm_in, Phi_csm], [L_endplate, Lm])
        Phi_csm_to_csm_out = interpolate([Phi_csm, Phi_csm_out], [Lm, L_endplate])
        Phi_csm_out_to_cgc = interpolate([Phi_csm_out, Phi_cgc], [L_endplate, Lgc / n_gc])
        Phi_cgc_to_cem_in = interpolate([Phi_cgc, Phi_cem_in], [Lgc / n_gc, L_endplate])
        Phi_cem_in_to_cem = interpolate([Phi_cem_in, Phi_cem], [L_endplate, Lm])
        Phi_cem_to_cem_out = interpolate([Phi_cem, Phi_cem_out], [Lm, L_endplate])
        Phi_cem_out_to_ext = None

    else:  # type_auxiliary == "no_auxiliary"
        C_v_agc_to_agc = [None] + [(sv[f'C_v_agc_{i}'] + sv[f'C_v_agc_{i + 1}']) / 2 for i in range(1, n_gc)]
        C_v_cgc_to_cgc = [None] + [(sv[f'C_v_cgc_{i}'] + sv[f'C_v_cgc_{i + 1}']) / 2 for i in range(1, n_gc)]
        C_H2_agc_to_agc = [None] + [(sv[f'C_H2_agc_{i}'] + sv[f'C_H2_agc_{i + 1}']) / 2 for i in range(1, n_gc)]
        C_O2_cgc_to_cgc = [None] + [(sv[f'C_O2_cgc_{i}'] + sv[f'C_O2_cgc_{i + 1}']) / 2 for i in range(1, n_gc)]
        v_asm_out_to_agc = interpolate([v_asm_out, sv['v_agc_1']], [L_man_gc, Lgc / n_gc])
        v_agc_to_agc = [None] + [(sv[f'v_agc_{i}'] + sv[f'v_agc_{i + 1}']) / 2 for i in range(1, n_gc)]
        v_agc_to_aem_in = interpolate([sv[f'v_agc_{n_gc}'], v_aem_in], [Lgc / n_gc, L_man_gc])
        v_aem_out_to_ext = interpolate([v_aem_in, v_a_ext], [L_man_gc, Lgc / n_gc])
        v_csm_out_to_cgc = interpolate([v_csm_out, sv['v_cgc_1']], [L_man_gc, Lgc / n_gc])
        v_cgc_to_cgc = [None] + [(sv[f'v_cgc_{i}'] + sv[f'v_cgc_{i + 1}']) / 2 for i in range(1, n_gc)]
        v_cgc_to_cem_in = interpolate([sv[f'v_cgc_{n_gc}'], v_cem_in], [Lgc / n_gc, L_man_gc])
        v_cem_out_to_ext = interpolate([v_cem_in, v_c_ext], [L_man_gc, Lgc / n_gc])
        Pasm_out_to_agc = interpolate([Pasm_out, P['agc_1']], [L_man_gc, Lgc / n_gc])
        Pagc_to_aem_in = interpolate([P[f'agc_{n_gc}'], Paem_in], [Lgc / n_gc, L_man_gc])
        Paem_out_to_ext = interpolate([Paem_in, Pa_ext], [L_man_gc, Lgc / n_gc])
        Pcsm_out_to_cgc = interpolate([Pcsm_out, P['cgc_1']], [L_man_gc, Lgc / n_gc])
        Pcgc_to_cem_in = interpolate([P[f'cgc_{n_gc}'], Pcem_in], [Lgc / n_gc, L_man_gc])
        Pcem_out_to_ext = interpolate([Pcem_in, Pc_ext], [L_man_gc, Lgc / n_gc])
        Phi_asm_out_to_agc = interpolate([Phi_asm_out, Phi['agc_1']], [L_man_gc, Lgc / n_gc])
        Phi_agc_to_aem_in = interpolate([Phi[f'agc_{n_gc}'], Phi_aem_in], [Lgc / n_gc, L_man_gc])
        Phi_aem_out_to_ext = interpolate([Phi_aem_in, Phi_a_ext], [L_man_gc, Lgc / n_gc] )                     # Boundary condition: at the exit, flow is isothermal and pressure is constant through space, so humidity is constant through space.
        Phi_csm_out_to_cgc = interpolate([Phi_csm_out, Phi['cgc_1']], [L_man_gc, Lgc / n_gc])
        Phi_cgc_to_cem_in = interpolate([Phi[f'cgc_{n_gc}'], Phi_cem_in], [Lgc / n_gc, L_man_gc])
        Phi_cem_out_to_ext = interpolate([Phi_cem_in, Phi_c_ext], [L_man_gc, Lgc / n_gc])                      # Boundary condition: at the exit, flow is isothermal and pressure is constant through space, so humidity is constant through space.
        y_O2_csm_out_to_cgc = interpolate([y_O2['csm_out'], y_O2['cgc_1']], [L_man_gc, Lgc / n_gc])
        y_O2_cgc_to_cem_in = interpolate([y_O2[f'cgc_{n_gc}'], y_O2['cem_in']], [Lgc / n_gc, L_man_gc])
        (v_asm_in_to_asm, v_asm_to_asm_out, v_aem_in_to_aem, v_asm_in_re_to_asm, v_aem_to_aem_out_re, v_aem_to_aem_out,
         v_csm_in_to_csm, v_csm_to_csm_out, v_cem_in_to_cem, v_cem_to_cem_out) = [None] * 10
        (Pasm_in_to_asm, Pasm_to_asm_out, Paem_in_to_aem, Pasm_in_re_to_asm, Paem_to_aem_out_re, Paem_to_aem_out,
         Pcsm_in_to_csm, Pcsm_to_csm_out, Pcem_in_to_cem, Pcem_to_cem_out) = [None] * 10
        (Phi_asm_in_to_asm, Phi_asm_to_asm_out, Phi_aem_in_to_aem, Phi_asm_in_re_to_asm, Phi_aem_to_aem_out_re,
         Phi_aem_to_aem_out, Phi_csm_in_to_csm, Phi_csm_to_csm_out, Phi_cem_in_to_cem, Phi_cem_to_cem_out) = [None] * 10

    return (C_v_agc_to_agc, C_v_cgc_to_cgc, C_H2_agc_to_agc, C_O2_cgc_to_cgc, v_asm_in_to_asm, v_asm_to_asm_out,
            v_asm_out_to_agc, v_agc_to_agc, v_agc_to_aem_in, v_aem_in_to_aem,  v_asm_in_re_to_asm, v_aem_to_aem_out_re,
            v_aem_to_aem_out, v_aem_out_to_ext, v_csm_in_to_csm, v_csm_to_csm_out, v_csm_out_to_cgc, v_cgc_to_cgc,
            v_cgc_to_cem_in, v_cem_in_to_cem, v_cem_to_cem_out, v_cem_out_to_ext, Pasm_in_to_asm, Pasm_to_asm_out,
            Pasm_out_to_agc, Pagc_to_aem_in, Paem_in_to_aem, Pasm_in_re_to_asm, Paem_to_aem_out_re, Paem_to_aem_out,
            Paem_out_to_ext, Pcsm_in_to_csm, Pcsm_to_csm_out, Pcsm_out_to_cgc, Pcgc_to_cem_in, Pcem_in_to_cem,
            Pcem_to_cem_out, Pcem_out_to_ext, Phi_asm_in_to_asm, Phi_asm_to_asm_out, Phi_asm_out_to_agc,
            Phi_agc_to_aem_in, Phi_aem_in_to_aem, Phi_asm_in_re_to_asm, Phi_aem_to_aem_out_re, Phi_aem_to_aem_out,
            Phi_aem_out_to_ext, Phi_csm_in_to_csm, Phi_csm_to_csm_out, Phi_csm_out_to_cgc, Phi_cgc_to_cem_in,
            Phi_cem_in_to_cem, Phi_cem_to_cem_out, Phi_cem_out_to_ext, y_O2_csm_out_to_cgc, y_O2_cgc_to_cem_in)