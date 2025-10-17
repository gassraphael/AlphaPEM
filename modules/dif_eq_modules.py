# -*- coding: utf-8 -*-

"""This module is used to determine intermediate values for the calculation of the differential equations
and to implement integration events.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import Text, Pext, Phi_ext, M_H2, M_O2, M_N2, M_H2O, y_O2_ext, R, F
from modules.transitory_functions import average, Psat, C_v_sat, mu_mixture_gases, k_H2, k_O2, calculate_rho_Cp0


# ____________________________________________Differential equations modules____________________________________________

def calculate_dif_eq_int_values(t, sv, operating_inputs, parameters, Ware, **kwarks):
    """This functions calculates intermediate values for the calculation of the differential equations

        Parameters
        ----------
        t : float
            Time (s).
        sv : dict
            Variables calculated by the solver. They correspond to the fuel cell internal states.
            sv is a contraction of solver_variables for enhanced readability.
        operating_inputs : dict
            Operating inputs of the fuel cell.
        parameters : dict
            Parameters of the fuel cell model.
        Wasm_ext_to_in : float
            Mass flow entering inside the anode supply manifold (kg.s-1).
        Wcsm_ext_to_in : float
            Mass flow entering inside the cathode supply manifold (kg.s-1).
        Ware : float
            Mass flow inside the anode recirculation loop (kg.s-1).
        Ja_in : float
            Molar flow entering inside the anode gas channel (mol.m-2.s-1).
        Jc_in : float
            Molar flow entering inside the cathode gas channel (mol.m-2.s-1).

        Returns
        -------
        Mext : float
            Molar mass of the ambient air outside the stack (kg/mol).
        M_H2_N2_in : float
            Molar mass of the inlet gas at the anode side (H2/N2 mixture) (kg/mol).
        i_n : float
            Internal current density (A/m²).
        Masm : float
            Molar mass of all the gas species in the anode supply manifold (kg/mol).
        Maem : float
            Molar mass of all the gas species in the anode external manifold (kg/mol).
        Mcsm : float
            Molar mass of all the gas species in the cathode supply manifold (kg/mol).
        Mcem : float
            Molar mass of all the gas species in the cathode external manifold (kg/mol).
        rho_Cp0 : dict
            Volumetric heat capacity of each component in the stack (J.m-3.K-1).
        """

    # Extraction of the variables
    C_v_acl, C_v_ccl = sv['C_v_acl'], sv['C_v_ccl']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_acl, C_O2_ccl = sv['C_H2_acl'], sv['C_O2_ccl']
    C_N2_a, C_N2_c = sv['C_N2_a'], sv['C_N2_c']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']
    Pasm_in_re, Pasm_in, Pasm, Pasm_out = sv.get('Pasm_in_re', None), sv.get('Pasm_in', None), sv.get('Pasm', None), sv['Pasm_out']
    Paem_in, Paem, Paem_out, Paem_out_re, Pa_ext = sv['Paem_in'], sv.get('Paem', None), sv.get('Paem_out', None), sv.get('Paem_out_re', None), sv['Pa_ext']
    Pcsm_in, Pcsm, Pcsm_out = sv.get('Pcsm_in', None), sv.get('Pcsm', None), sv['Pcsm_out']
    Pcem_in, Pcem, Pcem_out, Pc_ext = sv['Pcem_in'], sv.get('Pcem', None), sv.get('Pcem_out', None), sv['Pc_ext']
    Phi_asm_in_re, Phi_asm_in, Phi_asm = sv.get('Phi_asm_in_re', None), sv.get('Phi_asm_in', None), sv.get('Phi_asm', None)
    Phi_asm_out, Phi_aem_in, Phi_aem = sv['Phi_asm_out'], sv['Phi_aem_in'], sv.get('Phi_aem', None)
    Phi_aem_out, Phi_aem_out_re, Phi_a_ext = sv.get('Phi_aem_out', None), sv.get('Phi_aem_out_re', None), sv['Phi_a_ext']
    Phi_csm_in, Phi_csm, Phi_csm_out = sv.get('Phi_csm_in', None), sv.get('Phi_csm', None), sv['Phi_csm_out']
    Phi_cem_in, Phi_cem, Phi_cem_out, Phi_c_ext = sv['Phi_cem_in'], sv.get('Phi_cem', None), sv.get('Phi_cem_out', None), sv['Phi_c_ext']

    # Extraction of the operating inputs and the parameters
    T_des, y_H2_in = operating_inputs['T_des'], operating_inputs['y_H2_in']
    Lgc, Lm, L_endplate, L_man_gc = parameters['Lgc'], parameters['Lm'], parameters['L_endplate'], parameters['L_man_gc']
    A_T_a, A_T_c, Hagc, Wagc = parameters['A_T_a'], parameters['A_T_c'], parameters['Hagc'], parameters['Wagc']
    Hmem, Hacl = parameters['Hmem'], parameters['Hacl']
    Hccl, epsilon_gdl, epsilon_cl = parameters['Hccl'], parameters['epsilon_gdl'], parameters['epsilon_cl']
    epsilon_mpl, kappa_co, epsilon_mc = parameters['epsilon_mpl'], parameters['kappa_co'], parameters['epsilon_mc']
    epsilon_atl, epsilon_ctl = parameters['epsilon_atl'], parameters['epsilon_ctl']
    Htl, n_gc, n_gdl, n_mpl, n_tl = parameters['Htl'], parameters['n_gc'], parameters['n_gdl'], parameters['n_mpl'], parameters['n_tl']
    t_purge, type_auxiliary, type_purge = parameters['t_purge'], parameters['type_auxiliary'], parameters['type_purge']

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
    #       Total concentration
    C_tot = {}
    for i in range(1, n_gc + 1):
        C_tot[f'agc_{i}'] = sv[f'C_v_agc_{i}'] + sv[f'C_H2_agc_{i}'] + C_N2_a
        C_tot[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] + sv[f'C_O2_cgc_{i}'] + C_N2_c

    #       Humidities
    Phi = {}
    for i in range(1, n_gc + 1):
        Phi[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] / C_v_sat(sv[f'T_cgc_{i}'])

    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_O2 = {}
    for i in range(1, n_gc + 1):
        y_O2[f'cgc_{i}'] = sv[f'C_O2_cgc_{i}'] / (sv[f'C_O2_cgc_{i}'] + C_N2_c)

    #       Molar masses
    for i in range(1, n_gc + 1):
        M[f'agc_{i}'] = sv[f'C_v_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2O + \
               sv[f'C_H2_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2 + \
               C_N2_a * R * T_des / P[f'agc_{i}'] * M_N2
        M[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] * R * T_des / P[f'cgc_{i}'] * M_H2O + \
                        sv[f'C_O2_cgc_{i}'] * R * T_des / P[f'cgc_{i}'] * M_O2 + \
                        C_N2_c * R * T_des / P[f'cgc_{i}'] * M_N2

    #       Internal current density
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
    i_n = i_H2 + i_O2

    #       Volumetric heat capacity (J.m-3.K-1)
    rho_Cp0 = {
        **{f'agdl_{i}': calculate_rho_Cp0('agdl', sv[f'T_agdl_{i}'], C_v=sv[f'C_v_agdl_{i}'],
                                          s=sv[f's_agdl_{i}'], C_H2=sv[f'C_H2_agdl_{i}'], C_N2=C_N2_a, epsilon=epsilon_gdl)
           for i in range(1, n_gdl + 1)},
        **{f'atl_{i}': calculate_rho_Cp0('atl', sv[f'T_atl_{i}'], C_v=sv[f'C_v_atl_{i}'], s=sv[f's_atl_{i}'],
                                         C_H2=sv[f'C_H2_atl_{i}'], C_N2=C_N2_a, epsilon=epsilon_atl[i], n_tl=n_tl,
                                         Htl=Htl, node=i)
           for i in range(1, n_tl + 1)},
        **{f'ampl_{i}': calculate_rho_Cp0('ampl', sv[f'T_ampl_{i}'], C_v=sv[f'C_v_ampl_{i}'],
                                          s=sv[f's_ampl_{i}'], C_H2=sv[f'C_H2_ampl_{i}'], C_N2=C_N2_a, epsilon=epsilon_mpl)
           for i in range(1, n_mpl + 1)},
        'acl': calculate_rho_Cp0('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl, C_N2=C_N2_a, C_H2=C_H2_acl,
                                 epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
        'mem': calculate_rho_Cp0('mem', T_mem, lambdaa=lambda_mem),
        'ccl': calculate_rho_Cp0('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl, C_O2=C_O2_ccl, C_N2=C_N2_c,
                                 epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
        **{f'cmpl_{i}': calculate_rho_Cp0('cmpl', sv[f'T_cmpl_{i}'], C_v=sv[f'C_v_cmpl_{i}'],
                                          s=sv[f's_cmpl_{i}'], C_O2=sv[f'C_O2_cmpl_{i}'], C_N2=C_N2_c, epsilon=epsilon_mpl)
           for i in range(1, n_mpl + 1)},
        **{f'ctl_{i}': calculate_rho_Cp0('ctl', sv[f'T_ctl_{i}'], C_v=sv[f'C_v_ctl_{i}'], s=sv[f's_ctl_{i}'],
                                         C_O2=sv[f'C_O2_ctl_{i}'], C_N2=C_N2_c, epsilon=epsilon_ctl[i], n_tl=n_tl,
                                         Htl=Htl, node=i)
           for i in range(1, n_tl + 1)},
        **{f'cgdl_{i}': calculate_rho_Cp0('cgdl', sv[f'T_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                                          s=sv[f's_cgdl_{i}'], C_O2=sv[f'C_O2_cgdl_{i}'], C_N2=C_N2_c, epsilon=epsilon_gdl)
           for i in range(1, n_gdl + 1)}
        }

    # Physical quantities inside the auxiliary system
    if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation" or \
       parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
        # Lengths
        Lman_to_endplate = (Lm + L_endplate) / 2 # Length from the manifold node to the endplate node.
        Lman_to_man_gc = (Lm + L_man_gc) / 2 # Length from the manifold node to the manifold_to_gas_channel node.
        Lman_gc_to_gc = (L_man_gc + Lgc / n_gc) / 2 # Length from the manifold_to_gas_channel node to the first gas channel node.
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

        # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the EM
        y_H2_aem_in = (Paem_in - Phi_aem_in * Psat(T_des) - C_N2_a * R * T_des) / (Paem_in - Phi_aem_in * Psat(T_des))
        y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_a * R * T_des) / (Paem - Phi_aem * Psat(T_des))
        y_H2_aem_out = (Paem_out - Phi_aem_out * Psat(T_des) - C_N2_a * R * T_des) / (Paem_out - Phi_aem_out * Psat(T_des))
        y_O2_cem_in = (Pcem_in - Phi_cem_in * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem_in - Phi_cem_in * Psat(T_cgc))
        y_O2_cem = (Pcem - Phi_cem * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem - Phi_cem * Psat(T_cgc))
        y_O2_cem_out = (Pcem_out - Phi_cem_out * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem_out - Phi_cem_out * Psat(T_cgc))

        # Molar masses at the anode side
        if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation":
            M['asm_in_re'] = Phi_asm_in_re * Psat(T_des) / Pasm_in_re * M_H2O + \
                         (1 - Phi_asm_in_re * Psat(T_des) / Pasm_in_re) * M_H2
            M['asm_in'] = Phi_asm_in * Psat(T_des) / Pasm_in * M_H2O + \
                      (1 - Phi_asm_in * Psat(T_des) / Pasm_in) * M_H2
            M['asm'] = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
            M['asm_out'] = Phi_asm_out * Psat(T_des) / Pasm_out * M_H2O + \
                       (1 - Phi_asm_out * Psat(T_des) / Pasm_out) * M_H2
            M['aem_in'] = Phi_aem_in * Psat(T_des) / Paem_in * M_H2O + \
                      (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_H2
            M['aem'] = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
            M['aem_out'] = Phi_aem_out * Psat(T_des) / Paem_out * M_H2O + \
                       (1 - Phi_aem_out * Psat(T_des) / Paem_out) * M_H2
            M['aem_out_re'] = Phi_aem_out_re * Psat(T_des) / Paem_out_re * M_H2O + \
                          (1 - Phi_aem_out_re * Psat(T_des) / Paem_out_re) * M_H2
        else: #parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
            M['asm_in'] = Phi_asm_in * Psat(T_des) / Pasm_in * M_H2O + \
                      y_H2_in * (1 - Phi_asm_in * Psat(T_des) / Pasm_in) * M_H2 + \
                      (1 - y_H2_in) * (1 - Phi_asm_in * Psat(T_des) / Pasm_in) * M_N2
            M['asm'] = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 + \
                   (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
            M['asm_out'] = Phi_asm_out * Psat(T_des) / Pasm_out * M_H2O + \
                       y_H2_in * (1 - Phi_asm_out * Psat(T_des) / Pasm_out) * M_H2 + \
                       (1 - y_H2_in) * (1 - Phi_asm_out * Psat(T_des) / Pasm_out) * M_N2
            M['aem_in'] = Phi_aem_in * Psat(T_des) / Paem_in * M_H2O + \
                      y_H2_aem_in * (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_H2 + \
                      (1 - y_H2_aem_in) * (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_N2
            M['aem'] = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 + \
                   (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
            M['aem_out'] = Phi_aem_out * Psat(T_des) / Paem_out * M_H2O + \
                       y_H2_aem_out * (1 - Phi_aem_out * Psat(T_des) / Paem_out) * M_H2 + \
                       (1 - y_H2_aem_out) * (1 - Phi_aem_out * Psat(T_des) / Paem_out) * M_N2
        # Molar masses at the cathode side
        M['csm_in'] = Phi_csm_in * Psat(T_des) / Pcsm_in * M_H2O + \
                  y_O2_ext * (1 - Phi_csm_in * Psat(T_des) / Pcsm_in) * M_O2 + \
                  (1 - y_O2_ext) * (1 - Phi_csm_in * Psat(T_des) / Pcsm_in) * M_N2
        M['csm'] = Phi_csm * Psat(T_des) / Pcsm * M_H2O + \
               y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 + \
               (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
        M['csm_out'] = Phi_csm_out * Psat(T_des) / Pcsm_out * M_H2O + \
                   y_O2_ext * (1 - Phi_csm_out * Psat(T_des) / Pcsm_out) * M_O2 + \
                   (1 - y_O2_ext) * (1 - Phi_csm_out * Psat(T_des) / Pcsm_out) * M_N2
        M['cem_in'] = Phi_cem_in * Psat(T_cgc) / Pcem_in * M_H2O + \
                  y_O2_cem_in * (1 - Phi_cem_in * Psat(T_cgc) / Pcem_in) * M_O2 + \
                  (1 - y_O2_cem_in) * (1 - Phi_cem_in * Psat(T_cgc) / Pcem_in) * M_N2
        M['cem'] = Phi_cem * Psat(T_des) / Pcem * M_H2O + \
               y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 + \
               (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2
        M['cem_out'] = Phi_cem_out * Psat(T_cgc) / Pcem_out * M_H2O + \
                   y_O2_cem_out * (1 - Phi_cem_out * Psat(T_cgc) / Pcem_out) * M_O2 + \
                   (1 - y_O2_cem_out) * (1 - Phi_cem_out * Psat(T_cgc) / Pcem_out) * M_N2

        # Density/concentration of the gas mixture.
        C_tot_a_in = Pasm_in / (R * T_des)
        rho_asm_in = Pasm_in / (R * T_des) * Masm_in
        rho_asm = Pasm / (R * T_des) * Masm
        rho_asm_out = Pasm_out / (R * T_des) * Masm_out
        rho_agc = P[f'agc_{i}'] / (R * sv[f'T_agc_{i}']) * Magc
        rho_aem_in = Paem_in / (R * T_des) * Maem_in
        rho_aem = Paem / (R * T_des) * Maem
        rho_aem_out = Paem_out / (R * T_des) * Maem_out
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            rho_asm_in_re = Pasm_in_re / (R * T_des) * Masm_in_re
            rho_aem_out_re = Paem_out_re / (R * T_des) * Maem_out_re
        else:
            rho_asm_in_re, rho_aem_out_re = None, None
        rho_a_ext = Pext / (R * T_des) * Maem_out
        C_tot_a_ext = Pext / (R * T_des)                                                                                # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed.
        C_tot_c_in = Pcsm_in / (R * T_des)
        rho_csm_in = Pcsm_in / (R * T_des) * Mcsm_in
        rho_csm = Pcsm / (R * T_des) * Mcsm
        rho_csm_out = Pcsm_out / (R * T_des) * Mcsm_out
        rho_cgc = P[f'cgc_{i}'] / (R * sv[f'T_cgc_{i}']) * Mcgc
        rho_cem_in = Pcem_in / (R * T_cgc) * Mcem_in
        rho_cem = Pcem / (R * T_cgc) * Mcem
        rho_cem_out = Pcem_out / (R * T_cgc) * Mcem_out
        rho_c_ext = Pext / (R * T_des) * Mcem_out
        C_tot_c_ext = Pext * Mcem_out / (R * T_des)                                                                     # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed.

        # Vapor ratio over the gas mixture.
        x_H2O_v_asm_in_re = Phi_asm_in_re * Psat(T_des) / Pasm_in_re
        x_H2O_v_asm_in = Phi_asm_in * Psat(T_des) / Pasm_in
        x_H2O_v_asm = Phi_asm * Psat(T_des) / Pasm
        x_H2O_v_asm_out = Phi_asm_out * Psat(T_des) / Pasm_out
        x_H2O_v_agc = C_v_agc / (C_v_agc + C_H2_agc + C_N2_a)
        x_H2O_v_aem_in = Phi_aem_in * Psat(T_des) / Paem_in
        x_H2O_v_aem = Phi_aem * Psat(T_des) / Paem
        x_H2O_v_aem_out = Phi_aem_out * Psat(T_des) / Paem_out
        x_H2O_v_aem_out_re = Phi_aem_out_re * Psat(T_des) / Paem_out_re
        x_H2O_v_a_ext = Phi_a_ext * Psat(T_des) / Pext
        x_H2O_v_csm_in = Phi_csm_in * Psat(T_des) / Pcsm_in
        x_H2O_v_csm = Phi_csm * Psat(T_des) / Pcsm
        x_H2O_v_csm_out = Phi_csm_out * Psat(T_des) / Pcsm_out
        x_H2O_v_cgc = C_v_cgc / (C_v_cgc + C_O2_cgc + C_N2_c)
        x_H2O_v_cem_in = Phi_cem_in * Psat(T_des) / Pcem_in
        x_H2O_v_cem = Phi_cem * Psat(T_des) / Pcem
        x_H2O_v_cem_out = Phi_cem_out * Psat(T_des) / Pcem_out
        x_H2O_v_c_ext = Phi_c_ext * Psat(T_des) / Pext

        # Molar fraction of H2 in the dry gas mixture (H2/N2)
        y_H2_agc = C_H2_agc / (C_H2_agc + C_N2_a)
        y_O2_cgc = C_O2_cgc / (C_O2_cgc + C_N2_c)

        # Dynamic viscosity of the gas mixture at the anode side.
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            mu_gaz_asm_in_re = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_asm_in_re, 1 - x_H2O_v_asm_in_re], T_des)
            mu_gaz_asm_in = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_asm_in, 1 - x_H2O_v_asm_in], T_des)
            mu_gaz_asm = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_asm, 1 - x_H2O_v_asm], T_des)
            mu_gaz_asm_out = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_asm_out, 1 - x_H2O_v_asm_out], T_des)
            mu_gaz_agc = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_agc, 1 - x_H2O_v_agc], T_agc)
            mu_gaz_aem_in = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_aem_in, 1 - x_H2O_v_aem_in], T_des)
            mu_gaz_aem = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_aem, 1 - x_H2O_v_aem], T_des)
            mu_gaz_aem_out = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_aem_out, 1 - x_H2O_v_aem_out], T_des)
            mu_gaz_aem_out_re = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_aem_out_re, 1 - x_H2O_v_aem_out_re], T_des)
            mu_gaz_a_ext = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_a_ext, 1 - x_H2O_v_a_ext], T_des)
        else:  # type_auxiliary == "forced-convective_cathode_with_flow-through_anode"
            mu_gaz_asm_in = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                            [x_H2O_v_asm_in, y_H2_in * (1 - x_H2O_v_asm_in), (1 - y_H2_in) * (1 - x_H2O_v_asm_in)],
                                             T_des)
            mu_gaz_asm = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                          [x_H2O_v_asm, y_H2_in * (1 - x_H2O_v_asm), (1 - y_H2_in) * (1 - x_H2O_v_asm)],
                                          T_des)
            mu_gaz_asm_out = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                              [x_H2O_v_asm_out, y_H2_in * (1 - x_H2O_v_asm_out), (1 - y_H2_in) * (1 - x_H2O_v_asm_out)],
                                              T_des)
            mu_gaz_agc = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                          [x_H2O_v_agc, y_H2_agc * (1 - x_H2O_v_agc),
                                           (1 - y_H2_agc) * (1 - x_H2O_v_agc)], T_agc)
            mu_gaz_aem_in = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                             [x_H2O_v_aem_in, y_H2_aem_in * (1 - x_H2O_v_aem_in), (1 - y_H2_aem_in) * (1 - x_H2O_v_aem_in)],
                                             T_des)
            mu_gaz_aem = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                          [x_H2O_v_aem, y_H2_aem * (1 - x_H2O_v_aem),
                                           (1 - y_H2_aem) * (1 - x_H2O_v_aem)], T_des)
            mu_gaz_aem_out = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                              [x_H2O_v_aem_out, y_H2_aem_out * (1 - x_H2O_v_aem_out), (1 - y_H2_aem_out) * (1 - x_H2O_v_aem_out)],
                                              T_des)
            mu_gaz_a_ext = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                            [x_H2O_v_a_ext, y_H2_aem_out * (1 - x_H2O_v_a_ext), (1 - y_H2_aem_out) * (1 - x_H2O_v_a_ext)],
                                            T_des)
        # Dynamic viscosity of the gas mixture at the cathode side.
        mu_gaz_csm_in = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                        [x_H2O_v_csm_in, y_O2_ext * (1 - x_H2O_v_csm_in), (1 - y_O2_ext) * (1 - x_H2O_v_csm_in)],
                                         T_des)
        mu_gaz_csm = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                      [x_H2O_v_csm, y_O2_ext * (1 - x_H2O_v_csm), (1 - y_O2_ext) * (1 - x_H2O_v_csm)],
                                      T_des)
        mu_gaz_csm_out = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                          [x_H2O_v_csm_out, y_O2_ext * (1 - x_H2O_v_csm_out), (1 - y_O2_ext) * (1 - x_H2O_v_csm_out)],
                                          T_des)
        mu_gaz_cgc = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                      [x_H2O_v_cgc, y_O2_cgc * (1 - x_H2O_v_cgc), (1 - y_O2_cgc) * (1 - x_H2O_v_cgc)],
                                      T_cgc)
        mu_gaz_cem_in = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                         [x_H2O_v_cem_in, y_O2_cem_in * (1 - x_H2O_v_cem_in), (1 - y_O2_cem_in) * (1 - x_H2O_v_cem_in)],
                                         T_des)
        mu_gaz_cem = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                      [x_H2O_v_cem, y_O2_cem * (1 - x_H2O_v_cem), (1 - y_O2_cem) * (1 - x_H2O_v_cem)],
                                      T_des)
        mu_gaz_cem_out = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                          [x_H2O_v_cem_out, y_O2_cem_out * (1 - x_H2O_v_cem_out), (1 - y_O2_cem_out) * (1 - x_H2O_v_cem_out)],
                                          T_des)
        mu_gas_c_ext = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                        [x_H2O_v_c_ext, y_O2_cem_out * (1 - x_H2O_v_c_ext),
                                         (1 - y_O2_cem_out) * (1 - x_H2O_v_c_ext)],
                                        T_des)

        # Boundary velocities
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            v_re = Ware / rho_aem_out_re / A_T_a
        else:  # type_auxiliary == "forced-convective_cathode_with_flow-through_anode"
            v_re = None

    else:  # parameters["type_auxiliary"] == "no_auxiliary"
        # Lengths
        Lman_gc_to_gc = (L_man_gc + Lgc / n_gc) / 2 # Length from the manifold_to_gas_channel node to the first gas channel node.
        # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2).
        y_O2_csm_out = y_O2_ext
        y_O2_cem_in = sv[f'C_O2_cgc_{n_gc}'] / (sv[f'C_O2_cgc_{n_gc}'] + C_N2_c)
        # Molar masses
        M['asm_out'] = Phi_asm_out * Psat(T_des) / Pasm_out * M_H2O + (1 - Phi_asm_out * Psat(T_des) / Pasm_out) * M_H2
        M['aem_in'] = Phi_aem_in * Psat(T_des) / Paem_in * M_H2O + (1 - Phi_aem_in * Psat(T_des) / Paem_in) * M_H2
        M['csm_out'] = Phi_csm_out * Psat(T_des) / Pcsm_out * M_H2O + \
                   y_O2_csm_out * (1 - Phi_csm_out * Psat(T_des) / Pcsm_out) * M_O2 + \
                   (1 - y_O2_csm_out) * (1 - Phi_csm_out * Psat(T_des) / Pcsm_out) * M_N2
        M['cem_in'] = Phi_cem_in * Psat(sv[f'T_cgc_{i}']) / Pcem_in * M_H2O + \
                  y_O2_cem_in * (1 - Phi_cem_in * Psat(sv[f'T_cgc_{i}']) / Pcem_in) * M_O2 + \
                  (1 - y_O2_cem_in) * (1 - Phi_cem_in * Psat(sv[f'T_cgc_{i}']) / Pcem_in) * M_N2
        # Density of the gas mixture.
        rho = {}
        rho['asm_out'] = Pasm_out / (R * T_des) * M['asm_out']
        for i in range(1, n_gc + 1):
            rho[f'agc_{i}'] = P[f'agc_{i}'] / (R * sv[f'T_agc_{i}']) * M[f'agc_{i}']
        rho['aem_in'] = Paem_in / (R * T_des) * M['aem_in']
        rho['a_ext'] = Pa_ext / (R * T_des) * M['aem_in']
        rho['csm_out'] = Pcsm_out / (R * T_des) * M['csm_out']
        for i in range(1, n_gc + 1):
            rho[f'cgc_{i}'] = P[f'cgc_{i}'] / (R * sv[f'T_cgc_{i}']) * M[f'cgc_{i}']
        rho['cem_in'] = Pcem_in / (R * T_des) * M['cem_in']
        rho['c_ext'] = Pc_ext / (R * T_des) * M['cem_in']
        # Concentration of the gas mixture.
        C_tot['a_in'] = Pasm_out / (R * T_des)                                                                             #Boundary condition: at the inlet, .............
        C_tot['asm_out'] = Pasm_out / (R * T_des)
        C_tot['aem_in'] = Paem_in / (R * T_des)
        C_tot['a_ext'] = Pa_ext / (R * T_des)                                                                              # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed too.
        C_tot['c_in'] = Pcsm_out / (R * T_des)
        C_tot['csm_out'] = Pcsm_out / (R * T_des)                                                                          #Boundary condition: at the inlet, .............
        C_tot['cem_in'] = Pcem_in / (R * T_des)
        C_tot['c_ext'] = Pc_ext / (R * T_des)                                                                              # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed too.
        # Vapor ratio over the gas mixture.
        x_H2O_v = {}
        x_H2O_v['asm_out'] = Phi_asm_out * Psat(T_des) / Pasm_out
        for i in range(1, n_gc + 1):
            x_H2O_v[f'agc_{i}'] = sv[f'C_v_agc_{i}'] / (sv[f'C_v_agc_{i}'] + sv[f'C_H2_agc_{i}'] + C_N2_a)
        x_H2O_v['aem_in'] = Phi_aem_in * Psat(T_des) / Paem_in
        x_H2O_v['a_ext'] = Phi_a_ext * Psat(T_des) / Pa_ext
        x_H2O_v['csm_out'] = Phi_csm_out * Psat(T_des) / Pcsm_out
        for i in range(1, n_gc + 1):
            x_H2O_v[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] / (sv[f'C_v_cgc_{i}'] + sv[f'C_O2_cgc_{i}'] + C_N2_c)
        x_H2O_v['cem_in'] = Phi_cem_in * Psat(T_des) / Pcem_in
        x_H2O_v['c_ext'] = Phi_c_ext * Psat(T_des) / Pc_ext
        # Dynamic viscosity of the gas mixture.
        mu_gaz = {}
        mu_gaz['asm_out'] = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v['asm_out'], 1 - x_H2O_v['asm_out']], T_des)
        for i in range(1, n_gc + 1):
            mu_gaz[f'agc_{i}'] = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v[f'agc_{i}'], 1 - x_H2O_v[f'agc_{i}']], sv[f'T_agc_{i}'])
        mu_gaz['aem_in'] = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v['aem_in'], 1 - x_H2O_v['aem_in']], T_des)
        mu_gaz['a_ext'] = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v['a_ext'], 1 - x_H2O_v['a_ext']], T_des)
        mu_gaz['csm_out'] = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                          [x_H2O_v['csm_out'], y_O2_csm_out * (1 - x_H2O_v['csm_out']), (1 - y_O2_csm_out) * (1 - x_H2O_v['csm_out'])],
                                          T_des)
        for i in range(1, n_gc + 1):
            mu_gaz[f'cgc_{i}'] = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                          [x_H2O_v[f'cgc_{i}'], y_O2[f'cgc_{i}'] * (1 - x_H2O_v[f'cgc_{i}']), (1 - y_O2[f'cgc_{i}']) * (1 - x_H2O_v[f'cgc_{i}'])],
                                             sv[f'T_cgc_{i}'])
        mu_gaz['cem_in'] = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                         [x_H2O_v['cem_in'], y_O2_cem_in * (1 - x_H2O_v['cem_in']), (1 - y_O2_cem_in) * (1 - x_H2O_v['cem_in'])],
                                         T_des)
        mu_gaz['c_ext'] = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                         [x_H2O_v['c_ext'], y_O2_cem_in * (1 - x_H2O_v['c_ext']), (1 - y_O2_cem_in) * (1 - x_H2O_v['c_ext'])],
                                         T_des)
        # Set to None the variables not used when "no_auxiliary" system is considered
        v_re, Lman_to_endplate, Lman_to_man_gc, k_purge = [None] * 4

    return {'i_n': i_n, 'rho_Cp0': rho_Cp0, 'v_re': v_re, 'Lman_to_endplate': Lman_to_endplate,
            'Lman_to_man_gc': Lman_to_man_gc, 'Lman_gc_to_gc': Lman_gc_to_gc, 'k_purge': k_purge, 'M': M, 'rho': rho,
            'C_tot': C_tot, 'mu_gaz': mu_gaz, 'P': P}


def desired_flows(dif_eq, solver_variables, control_variables, i_fc, di_fcdt, operating_inputs, parameters, i_n,
                  **kwargs):
    """
    This function calculates the desired flow for the air compressor and the humidifiers. These desired flow are
    different from the real ones as the corresponding machines takes time to reach the desired values.

    Parameters
    ----------
    solver_variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    control_variables : dict
        Variables controlled by the user.
    i_fc : float
        Fuel cell current density (A/m²).
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    i_n : float
        Internal current density (A/m²).

    Returns
    -------
    Wcp_des : float
        Desired air compressor flow rate (kg/s).
    Wa_inj_des : float
        Desired humidifier flow rate at the anode side (kg/s).
    Wc_inj_des : float
        Desired humidifier flow rate at the cathode side (kg/s).
    """

    # Extraction of the variables
    Pacp_des, Pasm, Pcsm, Wcp = solver_variables['Pacp_des'], solver_variables.get('Pasm', None), solver_variables.get('Pcsm', None), solver_variables.get('Wcp', None)
    Pasm_out, Pccp_des, Pcsm_out = solver_variables['Pasm_out'], solver_variables['Pccp_des'], solver_variables['Pcsm_out']
    # Extraction of the operating inputs and the parameters
    T_des, Sa, Sc = operating_inputs['T_des'], operating_inputs['Sa'], operating_inputs['Sc']
    y_H2_in = operating_inputs['y_H2_in']
    Phi_a_des, Phi_c_des = control_variables['Phi_a_des'], control_variables['Phi_c_des']
    Aact, n_cell, type_auxiliary = parameters['Aact'], parameters['n_cell'], parameters['type_auxiliary']

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # The desired air compressor volume flow rate (mol.s-1)
        Wacp_des = 1 / y_H2_in * Sa * (i_fc + i_n) / (2 * F) * (n_cell * Aact)
        Wccp_des = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * (i_fc + i_n) / (4 * F) * (n_cell * Aact)

        # The desired humidifier volume flow rate at the anode side Wa_v_inj_des (mol.s-1)
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            Prd = Pasm
            Wacp_des = 1 / y_H2_in * Sa * (i_fc + i_n) / (2 * F) * (n_cell * Aact)
            Wa_inj_des = (Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des)) /
                          (1 - Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des))) * Wacp_des)
        else:  # type_auxiliary == "forced-convective_cathode_with_anodic_recirculation"
            Wa_inj_des = 0

        # The desired humidifier volume flow rate at the cathode side Wc_inj_des (mol.s-1)
        Pcp = Pcsm
        Wv_hum_in = Phi_ext * Psat(Text) / Pext * Wcp  # Vapor flow rate from the outside
        Wc_v_des = Phi_c_des * Psat(T_des) / Pcp * Wcp  # Desired vapor flow rate
        Wc_inj_des = Wc_v_des - Wv_hum_in  # Desired humidifier flow rate

    else:  # elif type_auxiliary == "no_auxiliary":
        # At the anode side
        # Wacp_des = Sa * (i_fc + i_n) / (2 * F) * (n_cell * Aact)
        Wacp_des = Sa * 1.0e4 / (2 * F) * (n_cell * Aact)                                                               # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dWacp_desdt = Sa * di_fcdt / (2 * F) * (n_cell * Aact)

        Wa_inj_des = (Phi_a_des * Psat(T_des) / Pacp_des) / (1 - Phi_a_des * Psat(T_des) / Pacp_des) * Wacp_des         # Desired vapor flow rate
        dWa_inj_desdt = 0                                                                                               # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # At the cathode side
        # Wccp_des = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * (i_fc + i_n) / (4 * F) * (n_cell * Aact)
        Wccp_des = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * 1.0e4 / (4 * F) * (n_cell * Aact)         # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dWccp_desdt = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * di_fcdt / (4 * F) * (n_cell * Aact)

        Wc_v_hum_in = Phi_ext * Psat(Text) / Pext * Wccp_des  # Vapor flow rate from the outside
        dWc_v_hum_indt = Phi_ext * Psat(Text) / Pext * dWccp_desdt  # Vapor flow rate from the outside
        Wc_v_des = (Phi_c_des * Psat(T_des) / Pccp_des) / (1 - Phi_c_des * Psat(T_des) / Pccp_des) * Wccp_des           # Desired vapor flow rate
        dWc_v_desdt = 0                                                                                                 # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        Wc_inj_des = Wc_v_des - Wc_v_hum_in  # Desired humidifier flow rate
        dWc_inj_desdt = dWc_v_desdt - dWc_v_hum_indt  # Desired humidifier flow rate

    return Wacp_des, dWacp_desdt, Wa_inj_des, dWa_inj_desdt, Wccp_des, dWccp_desdt, Wc_inj_des, dWc_inj_desdt


# ______________________________________Function which gives the integration event______________________________________

def event_negative(t, y, operating_inputs, parameters, solver_variable_names, control_variables):
    """This function creates an event that will be checked at each step of solve_ivp integration. The integration stops
    if one of the crucial variables (C_v, lambda, C_O2, C_H2) becomes negative (or smaller than 1e-5).

    Parameters
    ----------
    t : float
        Time (s).
    y : numpy.ndarray
        Numpy list of the solver variables.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    solver_variable_names : list
        Names of the solver variables.
    control_variables : dict
        Variables controlled by the user.

    Returns
    -------
    The difference between the minimum value of the crucial variables and 1e-5.
    """

    negative_solver_variables = {} # Dictionary to store the crucial variables
    for index, key in enumerate(solver_variable_names):
        if (key.startswith("C_v_")) or (key.startswith("lambda_")) or \
                (key.startswith("C_O2_")) or (key.startswith("C_H2_")):
            negative_solver_variables[key] = y[index]
    return min(negative_solver_variables.values()) - 1e-5  # 1e-5 is a control parameter to stop the program before
    #                                                        having negative values.
