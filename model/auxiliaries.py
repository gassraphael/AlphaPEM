# -*- coding: utf-8 -*-

"""This file represents all the flows passing through the auxiliaries. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import F, R, Text, Pext, Phi_ext, y_O2_ext, M_H2O
from modules.transitory_functions import Psat
from modules.auxiliaries_modules import auxiliaries_int_values_which_are_commun_with_dif_eq
from model.velocity import calculate_velocity_evolution, desired_flows


# ______________________________________________________Auxiliaries_____________________________________________________

def auxiliaries(t, sv, control_variables, i_fc, Jv_agc_agdl, Jv_cgdl_cgc, J_H2_agc_agdl, J_O2_cgdl_cgc,
                operating_inputs, parameters):
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
    Jv_agc_agdl : list
        Vapor flow between the AGC and the AGDL at each GC node (mol.m-2.s-1).
    Jv_cgdl_cgc : list
        Vapor flow between the CGDL and the CGC at each GC node (mol.m-2.s-1).
    J_H2_agc_agdl : list
        H2 flow between the AGC and the AGDL at each GC node (mol.m-2.s-1).
    J_O2_cgdl_cgc : list
        O2 flow between the CGDL and the CGC at each GC node (mol.m-2.s-1).
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    Jv : dict
        Vapor flow between the different layers (mol.m-2.s-1).
    J_H2 : dict
        Hydrogen flow between the different layers (mol.m-2.s-1).
    J_O2 : dict
        Oxygen flow between the different layers (mol.m-2.s-1).
    J_N2 : dict
        Nitrogen flow between the different layers (mol.m-2.s-1).
    W : dict
        Global flows through the auxiliaries (mol.s-1).
    W_v : dict
        Vapor flows through the auxiliaries (mol.s-1).
    v_a_in : float
        Velocity evolution at the inlet of the anode (m.s-1).
    v_c_in : float
        Velocity evolution at the inlet of the cathode (m.s-1).
    Pa_in : float
        Inlet pressure at the anode side (Pa).
    Pc_in : float
        Inlet pressure at the cathode side (Pa).
    """

    # __________________________________________________Preliminaries___________________________________________________

    # Extraction of the variables
    Pasm, Paem, Pcsm, Pcem = sv.get('Pasm', None), sv.get('Paem', None), sv.get('Pcsm', None), sv.get('Pcem', None)
    Phi_asm, Phi_aem = sv.get('Phi_asm', None), sv.get('Phi_aem', None)
    Phi_csm, Phi_cem = sv.get('Phi_csm', None), sv.get('Phi_cem', None)
    Wcp, Wa_inj, Wc_inj = sv.get('Wcp', None), sv.get('Wa_inj', None), sv.get('Wc_inj', None)

    # Extraction of the operating inputs and the parameters
    T_des, Phi_a_des, Phi_c_des = operating_inputs['T_des'], operating_inputs['Phi_a_des'], operating_inputs['Phi_c_des']
    Sa, Sc, y_H2_in = operating_inputs['Sa'], operating_inputs['Sc'], operating_inputs['y_H2_in']
    Aact, nb_cell, Hagc, Hcgc = parameters['Aact'], parameters['nb_cell'], parameters['Hagc'], parameters['Hcgc']
    Wagc, Wcgc, nb_channel_in_gc = parameters['Wagc'], parameters['Wcgc'], parameters['nb_channel_in_gc']
    A_T_a, A_T_c = parameters['A_T_a'], parameters['A_T_c']
    nb_gc, type_auxiliary = parameters['nb_gc'], parameters['type_auxiliary']

    # Intermediate values
        # Commun intermediate values with dif_eq_modules.py which allows to avoid redundant calculations
    P, Phi, y_H2, y_O2, M, rho, k_purge, Abp_a, Abp_c, mu_gaz, i_n = \
        auxiliaries_int_values_which_are_commun_with_dif_eq(t, sv, operating_inputs, parameters)
    v_a, v_c, Pa_in, Pc_in = calculate_velocity_evolution(sv, control_variables, i_fc, i_n, Jv_agc_agdl, Jv_cgdl_cgc, J_H2_agc_agdl,
                                            J_O2_cgdl_cgc, operating_inputs, parameters, mu_gaz)
    W_des = desired_flows(sv, control_variables, i_fc, i_n, Pa_in, Pc_in, operating_inputs, parameters)

    # _________________________________________Inlet and outlet global flows____________________________________________
    """Global flows here refer to flows that integrate all the chemical species circulating together.
    Slight differences are to be noted in the expression of these flows depending on the type of auxiliary selected.
    """

    # Anode flow through the auxiliaries in mol.s-1
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
         type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # Wa_in = rho_asm_in_to_asm * v_a * A_T_a
        # Wasm_to_asm_out = rho_asm_to_asm_out * v_a * Hagc * Wagc
        # Wasm_out_to_agc = rho_asm_out_to_agc * v_a * Hagc * Wagc
        # Wagc_to_aem_in = rho_agc_to_aem_in * v_a * Hagc * Wagc
        # Waem_in_to_aem = rho_aem_in_to_aem * v_a * Hagc * Wagc
        # if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":                                     # Attention: prévoir un débit minimal pour la pompe, comme les débits entrants.
        #     Ware = Maem_out_re * (Paem_out_re / (Paem_out_re - Phi_aem_out_re * Psat(T_des))) * \
        #            (Sa - 1) * i_fc / (2 * F) * (nb_cell * Aact)                                                  # The pump exactly compensates the pressure drop.
        #     Wasm_in_re_to_asm = rho_asm_in_re_to_asm * v_a * A_T_a
        #     Waem_to_aem_out_re = rho_aem_to_aem_out_re * v_a * A_T_a
        #     Waem_to_aem_out = k_purge * rho_aem_to_aem_out * v_a * A_T_a
        #     Wa_out = k_purge * rho_aem_out_to_ext * v_a * A_T_a
        # else: # type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        #     Ware = None
        #     Wasm_in_re_to_asm = None
        #     Waem_to_aem_out_re = None
        #     Waem_to_aem_out = rho_aem_to_aem_out * v_a * Abp_a
        #     Wa_out = rho_aem_out_to_ext * v_a * Abp_a
    else:  # elif type_auxiliary == "no_auxiliary" (only 1 cell):
        Wa_in = W_des['H2'] + W_des['H2O_inj_a']  # This expression is also present in calculate_velocity_evolution.
        Wa_out = P[f'agc_{nb_gc}'] / (R * T_des) * v_a[nb_gc] * Hagc * Wagc * nb_cell * nb_channel_in_gc
        Ware, Wasm_to_agc, Wagc_to_aem = [None] * 3
    # Anode flow entering/leaving the stack in mol.m-2.s-1
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # Ja_in = 0
        # Ja_out = 0
    else:  # elif type_auxiliary == "no_auxiliary" (only 1 cell):
        Ja_in = Wa_in / (Hagc * Wagc) / nb_cell / nb_channel_in_gc # This expression is also present in calculate_velocity_evolution.
        Ja_out = Wa_out / (Hagc * Wagc) / nb_cell / nb_channel_in_gc

    # Cathode flow through the auxiliaries in mol.s-1
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # Wc_in = rho_csm_in_to_csm * v_c * A_T_c
        # Wcsm_to_csm_out = rho_csm_to_csm_out * v_c * Hcgc * Wcgc
        # Wcsm_out_to_cgc = rho_csm_out_to_cgc * v_c * Hcgc * Wcgc
        # Wcgc_to_cem_in = rho_cgc_to_cem_in * v_c * Hcgc * Wcgc
        # Wcem_in_to_cem = rho_cem_in_to_cem * v_c * Hcgc * Wcgc
        # Wcem_to_cem_out = rho_cem_to_cem_out * v_c * Abp_c
        # Wc_out = rho_cem_out_to_ext * v_c * Abp_c
    else:  # elif type_auxiliary == "no_auxiliary" (only 1 cell):
        Wc_in = W_des['dry_air'] + W_des['H2O_inj_c']  # This expression is also present in calculate_velocity_evolution.
        Wc_out = P[f'cgc_{nb_gc}'] / (R * T_des) * v_c[nb_gc] * Hcgc * Wcgc * nb_cell * nb_channel_in_gc
        Wcsm_to_cgc, Wcgc_to_cem = [None] * 2

    # Cathode flow entering/leaving the stack in mol.m-2.s-1
    Jc_in = Wc_in / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc  # This expression is also present in calculate_velocity_evolution.
    Jc_out = Wc_out / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc

    # ________________________________________Inlet and outlet specific flows___________________________________________
    """Specific flows here refer to flows that integrate only a single chemical species within the ensemble of species
     circulating together. For example, only the water vapor flow within the ensemble of hydrogen and water vapor.
     """

    # Vapor flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # Jv_agc_in = Phi_asm * Psat(T_des) / Pasm * Ja_in
    else:  # elif type_auxiliary == "no_auxiliary":
        Jv_agc_in = Phi_a_des * Psat(T_des) / Pa_in * Ja_in
    Jv_agc_agc = [None] + [sv[f'C_v_agc_{i}'] * v_a[i] for i in range(1, nb_gc)]
    Jv_agc_out = sv[f'C_v_agc_{nb_gc}'] * R * T_des / P[f'agc_{nb_gc}'] * Ja_out
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # Jv_cgc_in = Phi_csm * Psat(T_des) / Pcsm * Jc_in
    else:  # elif type_auxiliary == "no_auxiliary":
        Jv_cgc_in = Phi_c_des * Psat(T_des) / Pc_in * Jc_in
    Jv_cgc_cgc = [None] + [sv[f'C_v_cgc_{i}'] * v_c[i] for i in range(1, nb_gc)]
    Jv_cgc_out = sv[f'C_v_cgc_{nb_gc}'] * R * T_des / P[f'cgc_{nb_gc}'] * Jc_out

    # H2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # J_H2_agc_in = y_H2['asm_out'] * (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        # J_H2_agc_agc = None
        # J_H2_agc_out = y_H2_agc * (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out
    else:  # elif type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or type_auxiliary == "no_auxiliary":
        J_H2_agc_in = (1 - Phi_a_des * Psat(T_des) / Pa_in) * Ja_in
        J_H2_agc_agc = [None] + [sv[f'C_H2_agc_{i}'] * v_a[i] for i in range(1, nb_gc)]
        J_H2_agc_out = sv[f'C_H2_agc_{nb_gc}'] * R * T_des / P[f'agc_{nb_gc}'] * Ja_out

    # O2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # J_O2_cgc_in = y_O2_csm * (1 - Phi_csm * Psat(T_des) / Pcsm) * Jc_in
    else:  # elif type_auxiliary == "no_auxiliary":
        J_O2_cgc_in = y_O2_ext * (1 - Phi_c_des * Psat(T_des) / Pc_in) * Jc_in
    J_O2_cgc_cgc = [None] + [sv[f'C_O2_cgc_{i}'] * v_c[i] for i in range(1, nb_gc)]
    J_O2_cgc_out = sv[f'C_O2_cgc_{nb_gc}'] * R * T_des / P[f'cgc_{nb_gc}'] * Jc_out

    # N2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # J_N2_agc_in = (1 - y_H2['asm_out']) * (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        # J_N2_agc_out = (1 - y_H2_agc) * (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out
        # J_N2_cgc_in = (1 - y_O2_csm_out_to_cgc) * (1 - Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc) * Jc_in
        # J_N2_cgc_out = (1 - y_O2_cgc_to_cem_in) * (1 - Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in) * Jc_out
    else:  # elif type_auxiliary == "no_auxiliary":
        J_N2_agc_in = 0
        J_N2_agc_agc = [None] + [0] * (nb_gc - 1)
        J_N2_agc_out = 0
        J_N2_cgc_in = (1 - y_O2_ext) * (1 - Phi_c_des * Psat(T_des) / Pc_in) * Jc_in
        J_N2_cgc_cgc = [None] + [sv[f'C_N2_cgc_{i}'] * v_c[i] for i in range(1, nb_gc)]
        J_N2_cgc_out = sv[f'C_N2_cgc_{nb_gc}'] * R * T_des / P[f'cgc_{nb_gc}'] * Jc_out

    # Vapor flows at the manifold (mol.s-1)
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # Wv_asm_in_to_asm = Phi_asm_in_to_asm * Psat(T_des) / Pasm_in_to_asm * Wa_in
        # Wv_asm_to_asm_out = Phi_asm_to_asm_out * Psat(T_des) / Pasm_to_asm_out * Wasm_to_asm_out
        # Wv_asm_out_to_agc = Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc * Wasm_out_to_agc
        # Wv_agc_to_aem_in = Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in * Wagc_to_aem_in
        # Wv_aem_in_to_aem = Phi_aem_in_to_aem * Psat(T_des) / Paem_in_to_aem * Waem_in_to_aem
        # Wv_aem_to_aem_out = Phi_aem_to_aem_out * Psat(T_des) / Paem_to_aem_out * Waem_to_aem_out
        # Wv_a_out = Phi_aem_out * Psat(T_des) / Paem_out * Wa_out
        # if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
        #     # At the anode side
        #     Wv_asm_ext_to_in = 0
        #     Wv_asm_in_re_to_asm = Phi_asm_in_re_to_asm * Psat(T_des) / Pasm_in_re_to_asm * Wasm_in_re_to_asm
        #     Wv_aem_to_aem_out_re = Phi_aem_to_aem_out_re * Psat(T_des) / Paem_to_aem_out_re * Waem_to_aem_out_re
        #     Wv_are = Phi_aem_out_re * Psat(T_des) / Paem_out_re * (Ware / M['aem_out_re']) # The pump exactly compensates the pressure drop.
        # else: # type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        #     # At the anode side
        #     Wv_asm_ext_to_in = Wa_inj / M_H2O
        #     Wv_asm_in_re_to_asm = None
        #     Wv_aem_to_aem_out_re = None
        #     Wv_are = None
        # # At the cathode side
        # Wv_csm_ext_to_in = Phi_ext * Psat(Text) / Pext * (Wcp / M['ext']) + Wc_inj / M_H2O
        # Wv_csm_in_to_csm = Phi_csm_in_to_csm * Psat(T_des) / Pcsm_in_to_csm * Wc_in
        # Wv_csm_to_csm_out = Phi_csm_to_csm_out * Psat(T_des) / Pcsm_to_csm_out * Wcsm_to_csm_out
        # Wv_csm_out_to_cgc = Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc * Wcsm_out_to_cgc
        # Wv_cgc_to_cem_in = Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in * Wcgc_to_cem_in
        # Wv_cem_in_to_cem = Phi_cem_in_to_cem * Psat(T_des) / Pcem_in_to_cem * Wcem_in_to_cem
        # Wv_cem_to_cem_out = Phi_cem_to_cem_out * Psat(T_des) / Pcem_to_cem_out * Wcem_to_cem_out
        # Wv_c_out = Phi_cem_out * Psat(T_des) / Pcem_out * Wc_out
    else:  # elif type_auxiliary == "no_auxiliary":
        Wv_are, Wv_a_in, Wv_asm_to_agc, Wv_agc_to_aem, Wv_a_out = [None] * 5
        Wv_c_in, Wv_csm_to_cgc, Wv_cgc_to_cem, Wv_c_out = [None] * 4

    return {'Jv': {'agc_in': Jv_agc_in, 'agc_agc': Jv_agc_agc, 'agc_out': Jv_agc_out,
                   'cgc_in': Jv_cgc_in, 'cgc_cgc': Jv_cgc_cgc, 'cgc_out': Jv_cgc_out},
            'J_H2': {'agc_in': J_H2_agc_in, 'agc_agc': J_H2_agc_agc, 'agc_out': J_H2_agc_out},
            'J_O2': {'cgc_in': J_O2_cgc_in, 'cgc_cgc': J_O2_cgc_cgc, 'cgc_out': J_O2_cgc_out},
            'J_N2': {'agc_in': J_N2_agc_in, 'agc_agc': J_N2_agc_agc, 'agc_out': J_N2_agc_out,
                     'cgc_in': J_N2_cgc_in, 'cgc_cgc': J_N2_cgc_cgc, 'cgc_out': J_N2_cgc_out},
            'W': {'are': Ware, 'a_in': Wa_in, 'asm_to_agc': Wasm_to_agc, 'agc_to_aem': Wagc_to_aem, 'a_out': Wa_out,
                  'c_in': Wc_in, 'csm_to_cgc': Wcsm_to_cgc, 'cgc_to_cem': Wcgc_to_cem, 'c_out': Wc_out},
            'W_v': {'a_in': Wv_a_in, 'are': Wv_are, 'asm_to_agc': Wv_asm_to_agc, 'agc_to_aem': Wv_agc_to_aem,
                    'a_out': Wv_a_out, 'c_in': Wv_c_in, 'csm_to_cgc': Wv_csm_to_cgc, 'cgc_to_cem': Wv_cgc_to_cem,
                    'c_out': Wv_c_out},
            'v_a_in': v_a[0], 'v_c_in': v_c[0], 'Pa_in': Pa_in, 'Pc_in': Pc_in}