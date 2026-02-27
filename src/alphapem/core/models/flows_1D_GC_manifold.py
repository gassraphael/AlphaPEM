# -*- coding: utf-8 -*-

"""This file represents all the flows passing through the auxiliaries. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from alphapem.utils.physics_constants import R, y_O2_ext, K_v_liq_gas, D_liq_dif
from alphapem.utils.maths_functions import d_dx
from alphapem.utils.physics_functions import rho_H2O_l, Psat
from alphapem.core.modules.flows_1D_GC_manifold_modules import flow_1D_GC_manifold_int_values
from alphapem.core.models.velocity import desired_flows


# ______________________________________________________Auxiliaries_____________________________________________________

def calculate_flows_1D_GC_manifold(sv_1D_cell, sv_1D_manifold, sv_auxiliary, i_fc_cell, v_a, v_c, Pa_in,
                                   Pc_in, operating_inputs, parameters):
    """This function calculates the flows passing through the auxiliaries.

    Parameters
    ----------
    sv_1D_cell : dict
        Variables calculated by the solver. They correspond to the cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    sv_1D_manifold : dict
        Variables calculated by the solver. They correspond to the manifold internal states.
        sv is a contraction of solver_variables for enhanced readability.
    sv_auxiliary : dict
        Variables calculated by the solver. They correspond to the auxiliary internal states.
        sv is a contraction of solver_variables for enhanced readability.
    i_fc_cell : float
        Fuel cell current density at time t (A.m-2).
    v_a : list
        Velocity evolution at the anode side (m.s-1).
    v_c : list
        Velocity evolution at the cathode side (m.s-1).
    Pa_in : float
        Inlet pressure at the anode side (Pa).
    Pc_in : float
        Inlet pressure at the cathode side (Pa).
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

    # Extraction of the operating inputs and the parameters
    T_des, Phi_a_des, Phi_c_des = operating_inputs['T_des'], operating_inputs['Phi_a_des'], operating_inputs['Phi_c_des']
    Sa, Sc, y_H2_in = operating_inputs['Sa'], operating_inputs['Sc'], operating_inputs['y_H2_in']
    Aact, nb_cell, Hagc, Hcgc = parameters['Aact'], parameters['nb_cell'], parameters['Hagc'], parameters['Hcgc']
    Wagc, Wcgc, Lgc, nb_channel_in_gc = parameters['Wagc'], parameters['Wcgc'], parameters['Lgc'], parameters['nb_channel_in_gc']
    A_T_a, A_T_c = parameters['A_T_a'], parameters['A_T_c']
    nb_gc, type_auxiliary = parameters['nb_gc'], parameters['type_auxiliary']

    # Extraction of the variables
    Pasm, Paem = sv_1D_manifold.get('Pasm', None), sv_1D_manifold.get('Paem', None)
    Pcsm, Pcem = sv_1D_manifold.get('Pcsm', None), sv_1D_manifold.get('Pcem', None)
    Phi_asm, Phi_aem = sv_1D_manifold.get('Phi_asm', None), sv_1D_manifold.get('Phi_aem', None)
    Phi_csm, Phi_cem = sv_1D_manifold.get('Phi_csm', None), sv_1D_manifold.get('Phi_cem', None)
    Wcp, Wa_inj = sv_auxiliary.get('Wcp', None), sv_auxiliary.get('Wa_inj', None)
    Wc_inj = sv_auxiliary.get('Wc_inj', None)
    C_v_agc = [None] + [sv_1D_cell[i]['C_v_agc'] for i in range(1, nb_gc + 1)]
    C_v_cgc = [None] + [sv_1D_cell[i]['C_v_cgc'] for i in range(1, nb_gc + 1)]
    s_agc = [None] + [sv_1D_cell[i]['s_agc'] for i in range(1, nb_gc + 1)] + [None]                                             # Adding a placeholder for the boundary condition at the outlet
    s_cgc = [None] + [sv_1D_cell[i]['s_cgc'] for i in range(1, nb_gc + 1)] + [None]                                             # Adding a placeholder for the boundary condition at the outlet
    C_H2_agc = [None] + [sv_1D_cell[i]['C_H2_agc'] for i in range(1, nb_gc + 1)]
    C_O2_cgc = [None] + [sv_1D_cell[i]['C_O2_cgc'] for i in range(1, nb_gc + 1)]
    C_N2_agc = [None] + [sv_1D_cell[i]['C_N2_agc'] for i in range(1, nb_gc + 1)]
    C_N2_cgc = [None] + [sv_1D_cell[i]['C_N2_cgc'] for i in range(1, nb_gc + 1)]

    # Intermediate values
    (P_agc, P_cgc, Phi_agc, Phi_cgc, y_H2_agc, y_O2_cgc, M_agc, M_cgc, M_ext, M_H2_N2_in, rho_agc, rho_cgc, k_purge,
     Abp_a, Abp_c, mu_gaz_agc, mu_gaz_cgc) = flow_1D_GC_manifold_int_values(sv_1D_cell, sv_auxiliary, operating_inputs,
                                                                            parameters)                                 # Some of them will remain useless ?!
    W_des = desired_flows(sv_1D_cell, i_fc_cell, Pa_in, Pc_in, operating_inputs, parameters)

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
        #            (Sa - 1) * i_fc_cell / (2 * F) * (nb_cell * Aact)                                                  # The pump exactly compensates the pressure drop.
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
        Wa_out = P_agc[nb_gc] / (R * T_des) * v_a[nb_gc] * Hagc * Wagc * nb_cell * nb_channel_in_gc
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
        Wc_out = P_cgc[nb_gc] / (R * T_des) * v_c[nb_gc] * Hcgc * Wcgc * nb_cell * nb_channel_in_gc
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
    Jv_agc_agc = [None] + [C_v_agc[i] * v_a[i] for i in range(1, nb_gc)]
    Jv_agc_out = C_v_agc[nb_gc] * R * T_des / P_agc[nb_gc] * Ja_out
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # Jv_cgc_in = Phi_csm * Psat(T_des) / Pcsm * Jc_in
    else:  # elif type_auxiliary == "no_auxiliary":
        Jv_cgc_in = Phi_c_des * Psat(T_des) / Pc_in * Jc_in
    Jv_cgc_cgc = [None] + [C_v_cgc[i] * v_c[i] for i in range(1, nb_gc)]
    Jv_cgc_out = C_v_cgc[nb_gc] * R * T_des / P_cgc[nb_gc] * Jc_out

    # Liquid water flows at the GC (kg.m-2.s-1)
    #   At the anode side
    s_agc[nb_gc + 1] = 0 # Boundary condition at the outlet of the anode GC: no liquid water at the outlet.
    Jl_agc_agc_conv = [None] + [rho_H2O_l(T_des) * K_v_liq_gas * v_a[i] * s_agc[i] for i in range(1, nb_gc + 1)]
    Jl_agc_agc_dif = [None] + [- D_liq_dif * d_dx(y_minus = s_agc[i], y_plus = s_agc[i + 1], dx = (Lgc / nb_gc) / 2)
                               for i in range(1, nb_gc + 1)]
    Jl_agc_agc = [None] + [Jl_agc_agc_conv[i] + Jl_agc_agc_dif[i] for i in range(1, nb_gc + 1)]
    Jl_agc_out = Jl_agc_agc_conv[-1] + Jl_agc_agc_dif[-1]
    #   At the cathode side
    s_cgc[nb_gc + 1] = 0 # Boundary condition at the outlet of the cathode GC: no liquid water at the outlet.
    Jl_cgc_cgc_conv = [None] + [rho_H2O_l(T_des) * K_v_liq_gas * v_c[i] * s_cgc[i] for i in range(1, nb_gc + 1)]
    Jl_cgc_cgc_dif = [None] + [- D_liq_dif * d_dx(y_minus = s_cgc[i], y_plus = s_cgc[i + 1], dx = (Lgc / nb_gc) / 2)
                               for i in range(1, nb_gc + 1)]
    Jl_cgc_cgc = [None] + [Jl_cgc_cgc_conv[i] + Jl_cgc_cgc_dif[i] for i in range(1, nb_gc + 1)]
    Jl_cgc_out = Jl_cgc_cgc_conv[-1] + Jl_cgc_cgc_dif[-1]

    # H2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # J_H2_agc_in = y_H2['asm_out'] * (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        # J_H2_agc_agc = None
        # J_H2_agc_out = y_H2_agc * (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out
    else:  # elif type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or type_auxiliary == "no_auxiliary":
        J_H2_agc_in = (1 - Phi_a_des * Psat(T_des) / Pa_in) * Ja_in
        J_H2_agc_agc = [None] + [C_H2_agc[i] * v_a[i] for i in range(1, nb_gc)]
        J_H2_agc_out = C_H2_agc[nb_gc] * R * T_des / P_agc[nb_gc] * Ja_out

    # O2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # J_O2_cgc_in = y_O2_csm * (1 - Phi_csm * Psat(T_des) / Pcsm) * Jc_in
    else:  # elif type_auxiliary == "no_auxiliary":
        J_O2_cgc_in = y_O2_ext * (1 - Phi_c_des * Psat(T_des) / Pc_in) * Jc_in
    J_O2_cgc_cgc = [None] + [C_O2_cgc[i] * v_c[i] for i in range(1, nb_gc)]
    J_O2_cgc_out = C_O2_cgc[nb_gc] * R * T_des / P_cgc[nb_gc] * Jc_out

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
        J_N2_cgc_cgc = [None] + [C_N2_cgc[i] * v_c[i] for i in range(1, nb_gc)]
        J_N2_cgc_out = C_N2_cgc[nb_gc] * R * T_des / P_cgc[nb_gc] * Jc_out

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
            'Jl': {'agc_agc': Jl_agc_agc, 'agc_out': Jl_agc_out, 'cgc_cgc': Jl_cgc_cgc, 'cgc_out': Jl_cgc_out},
            'J_H2': {'agc_in': J_H2_agc_in, 'agc_agc': J_H2_agc_agc, 'agc_out': J_H2_agc_out},
            'J_O2': {'cgc_in': J_O2_cgc_in, 'cgc_cgc': J_O2_cgc_cgc, 'cgc_out': J_O2_cgc_out},
            'J_N2': {'agc_in': J_N2_agc_in, 'agc_agc': J_N2_agc_agc, 'agc_out': J_N2_agc_out,
                     'cgc_in': J_N2_cgc_in, 'cgc_cgc': J_N2_cgc_cgc, 'cgc_out': J_N2_cgc_out},
            'W': {'are': Ware, 'a_in': Wa_in, 'asm_to_agc': Wasm_to_agc, 'agc_to_aem': Wagc_to_aem, 'a_out': Wa_out,
                  'c_in': Wc_in, 'csm_to_cgc': Wcsm_to_cgc, 'cgc_to_cem': Wcgc_to_cem, 'c_out': Wc_out},
            'W_v': {'a_in': Wv_a_in, 'are': Wv_are, 'asm_to_agc': Wv_asm_to_agc, 'agc_to_aem': Wv_agc_to_aem,
                    'a_out': Wv_a_out, 'c_in': Wv_c_in, 'csm_to_cgc': Wv_csm_to_cgc, 'cgc_to_cem': Wv_cgc_to_cem,
                    'c_out': Wv_c_out}}