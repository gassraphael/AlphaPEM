# -*- coding: utf-8 -*-

"""This file represents all the flows passing through the auxiliaries. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import F, R, Text, Pext, Phi_ext, y_O2_ext, M_H2O
from modules.transitory_functions import Psat
from modules.auxiliaries_modules import auxiliaries_int_values_which_are_commun_with_dif_eq, auxiliaries_int_values


# ______________________________________________________Auxiliaries_____________________________________________________

def auxiliaries(t, sv, control_variables, i_fc, operating_inputs, parameters):
    """This function calculates the flows passing through the auxiliaries.

    Parameters
    ----------
    t : float
        Time (s).
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
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
    Jv_a_in : float
        Vapor flow at the anode gas channel inlet (mol.m-2.s-1).
    Jv_a_out : float
        Vapor flow at the anode gas channel outlet (mol.m-2.s-1).
    Jv_c_in : float
        Vapor flow at the cathode gas channel inlet (mol.m-2.s-1).
    Jv_c_out : float
        Vapor flow at the cathode gas channel outlet (mol.m-2.s-1).
    J_H2_in : float
        H2 flow at the anode gas channel inlet (mol.m-2.s-1).
    J_H2_out : float
        H2 flow at the anode gas channel outlet (mol.m-2.s-1).
    J_O2_in : float
        O2 flow at the cathode gas channel inlet (mol.m-2.s-1).
    J_O2_out : float
        O2 flow at the cathode gas channel outlet (mol.m-2.s-1).
    J_N2_a_in : float
        N2 flow at the anode gas channel inlet (mol.m-2.s-1).
    J_N2_a_out : float
        N2 flow at the anode gas channel outlet (mol.m-2.s-1).
    J_N2_c_in : float
        N2 flow at the cathode gas channel inlet (mol.m-2.s-1).
    J_N2_c_out : float
        N2 flow at the cathode gas channel outlet (mol.m-2.s-1).
    Wasm_ext_to_in : float
        Anode side supply manifold inlet flow (kg.s-1).
    Wasm_out : float
        Anode side supply manifold outlet flow (kg.s-1).
    Waem_in : float
        Anode side external manifold inlet flow (kg.s-1).
    Waem_out_to_ext : float
        Anode side external manifold outlet flow (kg.s-1).
    Wcsm_ext_to_in : float
        Cathode side supply manifold inlet flow (kg.s-1).
    Wcsm_out : float
        Cathode side supply manifold outlet flow (kg.s-1).
    Wcem_in : float
        Cathode side external manifold inlet flow (kg.s-1).
    Wcem_out_to_ext : float
        Cathode side external manifold outlet flow (kg.s-1).
    Ware : float
        Anode side recirculation flow (kg.s-1).
    Wv_asm_ext_to_in : float
        Vapor flow at the anode supply manifold inlet (mol.s-1).
    Wv_aem_out_to_ext : float
        Vapor flow at the anode external manifold outlet (mol.s-1).
    Wv_csm_ext_to_in : float
        Vapor flow at the cathode supply manifold inlet (mol.s-1).
    Wv_cem_out_to_ext : float
        Vapor flow at the cathode external manifold outlet (mol.s-1).
    """

    # __________________________________________________Preliminaries___________________________________________________

    # Extraction of the variables
    v_aem_in, v_aem, v_aem_out, v_a_ext = sv['v_aem_in'], sv.get('v_aem', None), sv.get('v_aem_out', None), sv['v_a_ext']
    v_cem_in, v_cem, v_cem_out, v_c_ext = sv['v_cem_in'], sv.get('v_cem', None), sv.get('v_cem_out', None), sv['v_c_ext']
    Pasm, Pasm_out, Paem_in = sv.get('Pasm', None), sv['Pasm_out'], sv['Paem_in']
    Paem_out, Paem, Paem_out_re, Pa_ext = sv.get('Paem_out', None), sv.get('Paem', None), sv.get('Paem_out_re', None), sv['Pa_ext']
    Pcsm, Pcsm_out, Pcem_in, Pcem, Pcem_out, Pc_ext = sv.get('Pcsm', None), sv['Pcsm_out'], sv['Pcem_in'], sv.get('Pcem', None), sv.get('Pcem_out', None), sv['Pc_ext']
    Phi_asm_out, Phi_aem_in, Phi_aem = sv['Phi_asm_out'], sv['Phi_aem_in'], sv.get('Phi_aem', None)
    Phi_aem_out, Phi_aem_out_re, Phi_a_ext = sv.get('Phi_aem_out', None), sv.get('Phi_aem_out_re', None), sv['Phi_a_ext']
    Phi_csm_out, Phi_cem_in, Phi_cem_out, Phi_c_ext = sv['Phi_csm_out'], sv['Phi_cem_in'], sv.get('Phi_cem_out', None), sv['Phi_c_ext']
    Wcp, Wa_inj, Wc_inj = sv.get('Wcp', None), sv.get('Wa_inj', None), sv.get('Wc_inj', None)

    # Extraction of the operating inputs and the parameters
    T_des = operating_inputs['T_des']
    Sa, Sc, y_H2_in = operating_inputs['Sa'], operating_inputs['Sc'], operating_inputs['y_H2_in']
    Aact, n_cell, Hagc, Hcgc = parameters['Aact'], parameters['n_cell'], parameters['Hagc'], parameters['Hcgc']
    Wagc, Wcgc = parameters['Wagc'], parameters['Wcgc']
    A_T_a, A_T_c = parameters['A_T_a'], parameters['A_T_c']
    n_gc, type_auxiliary = parameters['n_gc'], parameters['type_auxiliary']

    # Intermediate values
        # Commun intermediate values with dif_eq_modules.py which allows to avoid redundant calculations
    (i_n, P, Phi, y_H2, y_O2, M, rho, k_purge, Abp_a, Abp_c) = auxiliaries_int_values_which_are_commun_with_dif_eq(t, sv, operating_inputs, parameters)
        # Other intermediate values
    (C_v_agc_to_agc, C_v_cgc_to_cgc, C_H2_agc_to_agc, C_O2_cgc_to_cgc, v_asm_in_to_asm, v_asm_to_asm_out,
            v_asm_out_to_agc, v_agc_to_agc, v_agc_to_aem_in, v_aem_in_to_aem,  v_asm_in_re_to_asm, v_aem_to_aem_out_re,
            v_aem_to_aem_out, v_aem_out_to_ext, v_csm_in_to_csm, v_csm_to_csm_out, v_csm_out_to_cgc, v_cgc_to_cgc,
            v_cgc_to_cem_in, v_cem_in_to_cem, v_cem_to_cem_out, v_cem_out_to_ext, Pasm_in_to_asm, Pasm_to_asm_out,
            Pasm_out_to_agc, Pagc_to_aem_in, Paem_in_to_aem, Pasm_in_re_to_asm, Paem_to_aem_out_re, Paem_to_aem_out,
            Paem_out_to_ext, Pcsm_in_to_csm, Pcsm_to_csm_out, Pcsm_out_to_cgc, Pcgc_to_cem_in, Pcem_in_to_cem,
            Pcem_to_cem_out, Pcem_out_to_ext, Phi_asm_in_to_asm, Phi_asm_to_asm_out, Phi_asm_out_to_agc,
            Phi_agc_to_aem_in, Phi_aem_in_to_aem, Phi_asm_in_re_to_asm, Phi_aem_to_aem_out_re, Phi_aem_to_aem_out,
            Phi_aem_out_to_ext, Phi_csm_in_to_csm, Phi_csm_to_csm_out, Phi_csm_out_to_cgc, Phi_cgc_to_cem_in,
            Phi_cem_in_to_cem, Phi_cem_to_cem_out, Phi_cem_out_to_ext, y_O2_csm_out_to_cgc, y_O2_cgc_to_cem_in) = auxiliaries_int_values(sv, parameters, P, Phi, y_O2)

    # _________________________________________Inlet and outlet global flows____________________________________________
    """Global flows here refer to flows that integrate all the chemical species circulating together.
    Slight differences are to be noted in the expression of these flows depending on the type of auxiliary selected.
    """

    # Anode flow through the auxiliaries in mol.s-1
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
         type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        Wasm_in_to_asm = rho_asm_in_to_asm * v_asm_in_to_asm * A_T_a
        Wasm_to_asm_out = rho_asm_to_asm_out * v_asm_to_asm_out * Hagc * Wagc
        Wasm_out_to_agc = rho_asm_out_to_agc * v_asm_out_to_agc * Hagc * Wagc
        Wagc_to_aem_in = rho_agc_to_aem_in * v_agc_to_aem_in * Hagc * Wagc
        Waem_in_to_aem = rho_aem_in_to_aem * v_aem_in_to_aem * Hagc * Wagc
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            Ware = Maem_out_re * (Paem_out_re / (Paem_out_re - Phi_aem_out_re * Psat(T_des))) * \
                   (Sa - 1) * (i_fc + i_n) / (2 * F) * (n_cell * Aact) # The pump exactly compensates the pressure drop.
            Wasm_in_re_to_asm = rho_asm_in_re_to_asm * v_asm_in_re_to_asm * A_T_a
            Waem_to_aem_out_re = rho_aem_to_aem_out_re * v_aem_to_aem_out_re * A_T_a
            Waem_to_aem_out = k_purge * rho_aem_to_aem_out * v_aem_to_aem_out * A_T_a
            Waem_out_to_ext = k_purge * rho_aem_out_to_ext * v_aem_out * A_T_a
        else: # type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            Ware = None
            Wasm_in_re_to_asm = None
            Waem_to_aem_out_re = None
            Waem_to_aem_out = rho_aem_to_aem_out * v_aem_to_aem_out * Abp_a
            Waem_out_to_ext = rho_aem_out_to_ext * v_aem_out * Abp_a
    else:  # elif type_auxiliary == "no_auxiliary" (only 1 cell):
        Wasm_out_to_agc = Pasm_out_to_agc * v_asm_out_to_agc * Hagc * Wagc / (R * T_des)
        Wagc_to_aem_in = Pagc_to_aem_in * v_agc_to_aem_in * Hagc * Wagc / (R * T_des)
        Waem_out_to_ext = Paem_out_to_ext * v_aem_out_to_ext * Hagc * Wagc / (R * T_des)
        Wa_ext = Pa_ext * v_a_ext * Hagc * Wagc / (R * T_des)                                                           # Boundary condition: at the exit, flow is isothermal, chemical composition remains unchanged and pressure is constant through space, so velocity is constant through space.
        Wasm_in_to_asm, Ware, Wasm_in_re_to_asm, Wasm_to_asm_out = [None] *4
        Waem_in_to_aem, Waem_to_aem_out_re, Waem_to_aem_out = [None] *3
    # Anode flow entering/leaving the stack in mol.m-2.s-1
    Ja_in = Wasm_out_to_agc / (Hagc * Wagc)
    Ja_out = Wagc_to_aem_in / (Hagc * Wagc)

    # Cathode flow through the auxiliaries in mol.s-1
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        Wcsm_in_to_csm = rho_csm_in_to_csm * v_csm_in_to_csm * A_T_c
        Wcsm_to_csm_out = rho_csm_to_csm_out * v_csm_to_csm_out * Hcgc * Wcgc
        Wcsm_out_to_cgc = rho_csm_out_to_cgc * v_csm_out_to_cgc * Hcgc * Wcgc
        Wcgc_to_cem_in = rho_cgc_to_cem_in * v_cgc_to_cem_in * Hcgc * Wcgc
        Wcem_in_to_cem = rho_cem_in_to_cem * v_cem_in_to_cem * Hcgc * Wcgc
        Wcem_to_cem_out = rho_cem_to_cem_out * v_cem_to_cem_out * Abp_c
        Wcem_out_to_ext = rho_cem_out_to_ext * v_cem_out * Abp_c
    else:  # elif type_auxiliary == "no_auxiliary" (only 1 cell):
        Wcsm_out_to_cgc = Pcsm_out_to_cgc * v_csm_out_to_cgc * Hcgc * Wcgc / (R * T_des)
        Wcgc_to_cem_in = Pcgc_to_cem_in * v_cgc_to_cem_in * Hcgc * Wcgc / (R * T_des)
        Wcem_out_to_ext = Pcem_out_to_ext * v_cem_out_to_ext * Hcgc * Wcgc / (R * T_des)
        Wc_ext = Pc_ext * v_c_ext * Hcgc * Wcgc / (R * T_des)                                                           # Boundary condition: at the exit, flow is isothermal, chemical composition remains unchanged and pressure is constant through space, so velocity is constant through space.
        Wcsm_in_to_csm, Wcsm_to_csm_out, Wcem_in_to_cem, Wcem_to_cem_out, Wcsm_out = [None] * 5

    # Cathode flow entering/leaving the stack in mol.m-2.s-1
    Jc_in = Wcsm_out_to_cgc / (Hcgc * Wcgc)
    Jc_out = Wcgc_to_cem_in / (Hcgc * Wcgc)

    # ________________________________________Inlet and outlet specific flows___________________________________________
    """Specific flows here refer to flows that integrate only a single chemical species within the ensemble of species
     circulating together. For example, only the water vapor flow within the ensemble of hydrogen and water vapor.
     """

    # Vapor flows at the GC (mol.m-2.s-1)
    Jv_a_in = Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc * Ja_in
    Jv_agc_agc = [None] + [C_v_agc_to_agc[i] * v_agc_to_agc[i] for i in range(1, n_gc)]
    Jv_a_out = Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in * Ja_out
    Jv_c_in = Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc * Jc_in
    Jv_cgc_cgc = [None] + [C_v_cgc_to_cgc[i] * v_cgc_to_cgc[i] for i in range(1, n_gc)]
    Jv_c_out = Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in * Jc_out

    # H2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        J_H2_in = y_H2['asm_out'] * (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        J_H2_agc_agc = None
        J_H2_out = y_H2_agc * (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out
    else:  # elif type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or type_auxiliary == "no_auxiliary":
        J_H2_in = (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        J_H2_agc_agc = [None] + [C_H2_agc_to_agc[i] * v_agc_to_agc[i] for i in range(1, n_gc)]
        J_H2_out = (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out

    # O2 flows at the GC (mol.m-2.s-1)
    J_O2_in = y_O2_csm_out_to_cgc * (1 - Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc) * Jc_in
    J_O2_cgc_cgc = [None] + [C_O2_cgc_to_cgc[i] * v_cgc_to_cgc[i] for i in range(1, n_gc)]
    J_O2_out = y_O2_cgc_to_cem_in * (1 - Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in) * Jc_out

    # N2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        J_N2_a_in = (1 - y_H2['asm_out']) * (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        J_N2_a_out = (1 - y_H2_agc) * (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out
    else:  # elif type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or type_auxiliary == "no_auxiliary":
        J_N2_a_in = 0
        J_N2_a_out = 0
    J_N2_c_in = (1 - y_O2_csm_out_to_cgc) * (1 - Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc) * Jc_in
    J_N2_c_out = (1 - y_O2_cgc_to_cem_in) * (1 - Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in) * Jc_out

    # Vapor flows at the manifold (mol.s-1)
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        Wv_asm_in_to_asm = Phi_asm_in_to_asm * Psat(T_des) / Pasm_in_to_asm * Wasm_in_to_asm
        Wv_asm_to_asm_out = Phi_asm_to_asm_out * Psat(T_des) / Pasm_to_asm_out * Wasm_to_asm_out
        Wv_asm_out_to_agc = Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc * Wasm_out_to_agc
        Wv_agc_to_aem_in = Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in * Wagc_to_aem_in
        Wv_aem_in_to_aem = Phi_aem_in_to_aem * Psat(T_des) / Paem_in_to_aem * Waem_in_to_aem
        Wv_aem_to_aem_out = Phi_aem_to_aem_out * Psat(T_des) / Paem_to_aem_out * Waem_to_aem_out
        Wv_aem_out_to_ext = Phi_aem_out * Psat(T_des) / Paem_out * Waem_out_to_ext
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            # At the anode side
            Wv_asm_ext_to_in = 0
            Wv_asm_in_re_to_asm = Phi_asm_in_re_to_asm * Psat(T_des) / Pasm_in_re_to_asm * Wasm_in_re_to_asm
            Wv_aem_to_aem_out_re = Phi_aem_to_aem_out_re * Psat(T_des) / Paem_to_aem_out_re * Waem_to_aem_out_re
            Wv_are = Phi_aem_out_re * Psat(T_des) / Paem_out_re * (Ware / M['aem_out_re']) # The pump exactly compensates the pressure drop.
        else: # type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            # At the anode side
            Wv_asm_ext_to_in = Wa_inj / M_H2O
            Wv_asm_in_re_to_asm = None
            Wv_aem_to_aem_out_re = None
            Wv_are = None
        # At the cathode side
        Wv_csm_ext_to_in = Phi_ext * Psat(Text) / Pext * (Wcp / M['ext']) + Wc_inj / M_H2O
        Wv_csm_in_to_csm = Phi_csm_in_to_csm * Psat(T_des) / Pcsm_in_to_csm * Wcsm_in_to_csm
        Wv_csm_to_csm_out = Phi_csm_to_csm_out * Psat(T_des) / Pcsm_to_csm_out * Wcsm_to_csm_out
        Wv_csm_out_to_cgc = Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc * Wcsm_out_to_cgc
        Wv_cgc_to_cem_in = Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in * Wcgc_to_cem_in
        Wv_cem_in_to_cem = Phi_cem_in_to_cem * Psat(T_des) / Pcem_in_to_cem * Wcem_in_to_cem
        Wv_cem_to_cem_out = Phi_cem_to_cem_out * Psat(T_des) / Pcem_to_cem_out * Wcem_to_cem_out
        Wv_cem_out_to_ext = Phi_cem_out * Psat(T_des) / Pcem_out * Wcem_out_to_ext
    else:  # elif type_auxiliary == "no_auxiliary":
        Wv_asm_out_to_agc = Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc * Wasm_out_to_agc
        Wv_agc_to_aem_in = Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in * Wagc_to_aem_in
        Wv_aem_out_to_ext = Phi_aem_out_to_ext * Psat(T_des) / Paem_out_to_ext * Waem_out_to_ext
        Wv_a_ext = Phi_a_ext * Psat(T_des) / Pa_ext * Wa_ext                                                            # Boundary condition: at the exit, flow is isothermal, chemical composition remains unchanged and pressure is constant through space, so humidity is constant through space.
        Wv_csm_out_to_cgc = Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc * Wcsm_out_to_cgc
        Wv_cgc_to_cem_in = Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in * Wcgc_to_cem_in
        Wv_cem_out_to_ext = Phi_cem_out_to_ext * Psat(T_des) / Pcem_out_to_ext * Wcem_out_to_ext
        Wv_c_ext = Phi_c_ext * Psat(T_des) / Pc_ext * Wc_ext                                                            # Boundary condition: at the exit, flow is isothermal, chemical composition remains unchanged and pressure is constant through space, so humidity is constant through space.
        Wv_asm_ext_to_in, Wv_asm_in_to_asm, Wv_asm_in_re_to_asm, Wv_asm_to_asm_out, Wv_aem_in_to_aem = [None] * 5
        Wv_aem_to_aem_out, Wv_aem_to_aem_out_re, Wv_are, Wv_csm_ext_to_in, Wv_csm_in_to_csm = [None] * 5
        Wv_csm_to_csm_out, Wv_cem_in_to_cem, Wv_cem_to_cem_out = [None] * 3

    return {'Ja_in': Jv_a_in, 'Jc_in': Jc_in, 'Jv_a_in': Jv_a_in, 'Jv_agc_agc': Jv_agc_agc, 'Jv_a_out': Jv_a_out,
            'Jv_c_in': Jv_c_in, 'Jv_cgc_cgc': Jv_cgc_cgc, 'Jv_c_out': Jv_c_out, 'J_H2_in': J_H2_in,
            'J_H2_agc_agc': J_H2_agc_agc, 'J_H2_out': J_H2_out, 'J_O2_in': J_O2_in, 'J_O2_cgc_cgc': J_O2_cgc_cgc,
            'J_O2_out': J_O2_out, 'J_N2_a_in': J_N2_a_in, 'J_N2_a_out': J_N2_a_out, 'J_N2_c_in': J_N2_c_in,
            'J_N2_c_out': J_N2_c_out, 'Ware': Ware, 'Wasm_in_re_to_asm': Wasm_in_re_to_asm,
            'Wasm_in_to_asm': Wasm_in_to_asm, 'Wasm_to_asm_out': Wasm_to_asm_out, 'Wasm_out_to_agc': Wasm_out_to_agc,
            'Wagc_to_aem_in': Wagc_to_aem_in, 'Waem_in_to_aem': Waem_in_to_aem, 'Waem_to_aem_out': Waem_to_aem_out,
            'Waem_to_aem_out_re': Waem_to_aem_out_re, 'Waem_out_to_ext': Waem_out_to_ext,
            'Wcsm_in_to_csm': Wcsm_in_to_csm, 'Wcsm_to_csm_out': Wcsm_to_csm_out, 'Wcsm_out_to_cgc': Wcsm_out_to_cgc,
            'Wcgc_to_cem_in': Wcgc_to_cem_in, 'Wcem_in_to_cem': Wcem_in_to_cem, 'Wcem_to_cem_out': Wcem_to_cem_out,
            'Wcem_out_to_ext': Wcem_out_to_ext, 'Wv_asm_ext_to_in': Wv_asm_ext_to_in, 'Wv_asm_in_to_asm': Wv_asm_in_to_asm,
            'Wv_asm_in_re_to_asm': Wv_asm_in_re_to_asm, 'Wv_asm_to_asm_out': Wv_asm_to_asm_out,
            'Wv_asm_out_to_agc': Wv_asm_out_to_agc, 'Wv_agc_to_aem_in': Wv_agc_to_aem_in,
            'Wv_aem_in_to_aem': Wv_aem_in_to_aem, 'Wv_aem_to_aem_out_re': Wv_aem_to_aem_out_re,
            'Wv_aem_to_aem_out': Wv_aem_to_aem_out, 'Wv_aem_out_to_ext': Wv_aem_out_to_ext, 'Wv_are': Wv_are,
            'Wv_a_ext': Wv_a_ext, 'Wv_csm_ext_to_in': Wv_csm_ext_to_in,
            'Wv_csm_in_to_csm': Wv_csm_in_to_csm, 'Wv_csm_to_csm_out': Wv_csm_to_csm_out,
            'Wv_csm_out_to_cgc': Wv_csm_out_to_cgc, 'Wv_cgc_to_cem_in': Wv_cgc_to_cem_in,
            'Wv_cem_in_to_cem': Wv_cem_in_to_cem, 'Wv_cem_to_cem_out': Wv_cem_to_cem_out,
            'Wv_cem_out_to_ext': Wv_cem_out_to_ext, 'Wv_c_ext': Wv_c_ext}
