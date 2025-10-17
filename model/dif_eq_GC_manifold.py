# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import Pext, R
from modules.transitory_functions import d_dx, d2_dx2, Psat


# ____________________________________________________Main functions____________________________________________________


def calculate_dyn_gas_evolution_inside_gas_channel(dif_eq, Hagc, Hcgc, Lgc, n_gc, type_auxiliary, Jv_a_in, Jv_agc_agc,
                                                   Jv_a_out, Jv_c_in, Jv_cgc_cgc, Jv_c_out, Jv_agc_agdl, Jv_cgdl_cgc,
                                                   J_H2_in, J_H2_agc_agc, J_H2_out, J_O2_in, J_O2_cgc_cgc,
                                                   J_O2_out, J_N2_a_in, J_N2_a_out, J_N2_c_in, J_N2_c_out,
                                                   J_H2_agc_agdl, J_O2_cgdl_cgc, **kwargs):
    """This function calculates the dynamic evolution of the vapor, hydrogen, oxygen and nitrogen gases in the gas
    channels.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    Hagc : float
        Thickness of the anode gas channel (m).
    Hcgc : float
        Thickness of the cathode gas channel (m).
    Lgc : float
        Length of the gas channel (m).
    type_auxiliary : str
        Type of auxiliary components used in the fuel cell system.
    Jv_a_in : float
        Water vapor flow at the anode inlet (mol.m-2.s-1).
    Jv_a_out : float
        Water vapor flow at the anode outlet (mol.m-2.s-1).
    Jv_c_in : float
        Water vapor flow at the cathode inlet (mol.m-2.s-1).
    Jv_c_out : float
        Water vapor flow at the cathode outlet (mol.m-2.s-1).
    Jv_agc_agdl : float
        Water vapor flow between the anode gas channel and the anode GDL (mol.m-2.s-1).
    Jv_cgdl_cgc : float
        Water vapor flow between the cathode GDL and the cathode gas channel (mol.m-2.s-1).
    J_H2_in : float
        Hydrogen flow at the anode inlet (mol.m-2.s-1).
    J_H2_out : float
        Hydrogen flow at the anode outlet (mol.m-2.s-1).
    J_O2_in : float
        Oxygen flow at the cathode inlet (mol.m-2.s-1).
    J_O2_out : float
        Oxygen flow at the cathode outlet (mol.m-2.s-1).
    J_N2_a_in : float
        Nitrogen flow at the anode inlet (mol.m-2.s-1).
    J_N2_a_out : float
        Nitrogen flow at the anode outlet (mol.m-2.s-1).
    J_N2_c_in : float
        Nitrogen flow at the cathode inlet (mol.m-2.s-1).
    J_N2_c_out : float
        Nitrogen flow at the cathode outlet (mol.m-2.s-1).
    J_H2_agc_agdl : float
        Hydrogen flow between the anode gas channel and the anode GDL (mol.m-2.s-1).
    J_O2_cgdl_cgc : float
        Oxygen flow between the cathode GDL and the cathode gas channel (mol.m-2.s-1).
    """

    # At the anode side, inside the AGC
    Jv_agc_agdl = [0] * (n_gc + 1)                                                                                      # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    J_H2_agc_agdl = [0] * (n_gc + 1)                                                                                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Jv_cgdl_cgc = [0] * (n_gc + 1)                                                                                      # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    J_O2_cgdl_cgc = [0] * (n_gc + 1)                                                                                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if n_gc == 1:
        dif_eq['dC_v_agc_1 / dt'] = (Jv_a_in - Jv_a_out) / Lgc - Jv_agc_agdl[1] / Hagc
    elif n_gc == 2:
        dif_eq['dC_v_agc_1 / dt'] = (Jv_a_in - Jv_agc_agc[1]) / (Lgc / n_gc) - Jv_agc_agdl[1] / Hagc
        dif_eq['dC_v_agc_2 / dt'] = (Jv_agc_agc[1] - Jv_a_out) / (Lgc / n_gc) - Jv_agc_agdl[2] / Hagc
    else: # n_gc > 2:
        dif_eq['dC_v_agc_1 / dt'] = (Jv_a_in - Jv_agc_agc[1]) / (Lgc / n_gc) - Jv_agc_agdl[1] / Hagc
        for i in range(2, n_gc):
            dif_eq[f'dC_v_agc_{i} / dt'] = (Jv_agc_agc[i - 1] - Jv_agc_agc[i]) / (Lgc / n_gc) - Jv_agc_agdl[i] / Hagc
        dif_eq[f'dC_v_agc_{n_gc} / dt'] = (Jv_agc_agc[n_gc - 1] - Jv_a_out) / (Lgc / n_gc) - Jv_agc_agdl[n_gc] / Hagc

    if n_gc == 1:
        dif_eq['dC_H2_agc_1 / dt'] = (J_H2_in - J_H2_out) / Lgc - J_H2_agc_agdl[1] / Hagc
    elif n_gc == 2:
        dif_eq['dC_H2_agc_1 / dt'] = (J_H2_in - J_H2_agc_agc[1]) / (Lgc / n_gc) - J_H2_agc_agdl[1] / Hagc
        dif_eq['dC_H2_agc_2 / dt'] = (J_H2_agc_agc[1] - J_H2_out) / (Lgc / n_gc) - J_H2_agc_agdl[2] / Hagc
    else: # n_gc > 2:
        dif_eq['dC_H2_agc_1 / dt'] = (J_H2_in - J_H2_agc_agc[1]) / (Lgc / n_gc) - J_H2_agc_agdl[1] / Hagc
        for i in range(2, n_gc):
            dif_eq[f'dC_H2_agc_{i} / dt'] = (J_H2_agc_agc[i - 1] - J_H2_agc_agc[i]) / (Lgc / n_gc) - J_H2_agc_agdl[i] / Hagc
        dif_eq[f'dC_H2_agc_{n_gc} / dt'] = (J_H2_agc_agc[n_gc - 1] - J_H2_out) / (Lgc / n_gc) - J_H2_agc_agdl[n_gc] / Hagc

    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dC_N2_a / dt'] = (J_N2_a_in - J_N2_a_out) / Lgc  # Test bench: simulated H2 recirculation which leads to N2 in the anode.
    else:
        dif_eq['dC_N2_a / dt'] = 0

    # At the cathode side, inside the CGC
    if n_gc == 1:
        dif_eq['dC_v_cgc_1 / dt'] = (Jv_c_in - Jv_c_out) / Lgc + Jv_cgdl_cgc[1] / Hcgc
    elif n_gc == 2:
        dif_eq['dC_v_cgc_1 / dt'] = (Jv_c_in - Jv_cgc_cgc[1]) / (Lgc / n_gc) + Jv_cgdl_cgc[1] / Hcgc
        dif_eq['dC_v_cgc_2 / dt'] = (Jv_cgc_cgc[1] - Jv_c_out) / (Lgc / n_gc) + Jv_cgdl_cgc[2] / Hcgc
    else: # n_gc > 2:
        dif_eq['dC_v_cgc_1 / dt'] = (Jv_c_in - Jv_cgc_cgc[1]) / (Lgc / n_gc) + Jv_cgdl_cgc[1] / Hcgc
        for i in range(2, n_gc):
            dif_eq[f'dC_v_cgc_{i} / dt'] = (Jv_cgc_cgc[i - 1] - Jv_cgc_cgc[i]) / (Lgc / n_gc) + Jv_cgdl_cgc[i] / Hcgc
        dif_eq[f'dC_v_cgc_{n_gc} / dt'] = (Jv_cgc_cgc[n_gc - 1] - Jv_c_out) / (Lgc / n_gc) + Jv_cgdl_cgc[n_gc] / Hcgc

    if n_gc == 1:
        dif_eq['dC_O2_cgc_1 / dt'] = (J_O2_in - J_O2_out) / Lgc + J_O2_cgdl_cgc[1] / Hcgc
    elif n_gc == 2:
        dif_eq['dC_O2_cgc_1 / dt'] = (J_O2_in - J_O2_cgc_cgc[1]) / (Lgc / n_gc) + J_O2_cgdl_cgc[1] / Hcgc
        dif_eq['dC_O2_cgc_2 / dt'] = (J_O2_cgc_cgc[1] - J_O2_out) / (Lgc / n_gc) + J_O2_cgdl_cgc[2] / Hcgc
    else: # n_gc > 2:
        dif_eq['dC_O2_cgc_1 / dt'] = (J_O2_in - J_O2_cgc_cgc[1]) / (Lgc / n_gc) + J_O2_cgdl_cgc[1] / Hcgc
        for i in range(2, n_gc):
            dif_eq[f'dC_O2_cgc_{i} / dt'] = (J_O2_cgc_cgc[i - 1] - J_O2_cgc_cgc[i]) / (Lgc / n_gc) + J_O2_cgdl_cgc[i] / Hcgc
        dif_eq[f'dC_O2_cgc_{n_gc} / dt'] = (J_O2_cgc_cgc[n_gc - 1] - J_O2_out) / (Lgc / n_gc) + J_O2_cgdl_cgc[n_gc] / Hcgc

    dif_eq['dC_N2_c / dt'] = (J_N2_c_in - J_N2_c_out) / Lgc


def calculate_dyn_temperature_evolution_inside_gas_channel(dif_eq, n_gc, **kwarks):
    """
    This function calculates the dynamic evolution of the temperature in the fuel cell.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    """

    # At the anode side, inside the AGC
    for i in range(1, n_gc + 1):
        dif_eq[f'dT_agc_{i} / dt'] = 0                                                                                  # Dirichlet boundary condition. T_agc is initialized to T_fc and remains constant.
    # At the cathode side, inside the CGC
    for i in range(1, n_gc + 1):
        dif_eq[f'dT_cgc_{i} / dt'] = 0                                                                                  # Dirichlet boundary condition. T_cgc is initialized to T_fc and remains constant.


def calculate_dyn_manifold_pressure_and_humidity_evolution(dif_eq, T_des, n_cell, Hagc, Hcgc, Wagc, Wcgc, Lgc, Vasm,
                                                           Vcsm, Vaem, Vcem, V_endplate_a, V_endplate_c, V_man_agc,
                                                           V_man_cgc, n_gc, type_auxiliary, Ware, Wasm_in_re_to_asm,
                                                           Wasm_in_to_asm, Wasm_to_asm_out, Wasm_out_to_agc,
                                                           Wagc_to_aem_in, Waem_in_to_aem,  Waem_to_aem_out,
                                                           Waem_to_aem_out_re, Waem_out_to_ext, Wcsm_in_to_csm,
                                                           Wcsm_to_csm_out, Wcsm_out_to_cgc, Wcgc_to_cem_in, Wcem_in_to_cem,
                                                           Wcem_to_cem_out, Wcem_out_to_ext, Wv_asm_ext_to_in,
                                                           Wv_asm_in_to_asm, Wv_asm_in_re_to_asm, Wv_asm_to_asm_out,
                                                           Wv_asm_out_to_agc, Wv_agc_to_aem_in, Wv_aem_in_to_aem,
                                                           Wv_aem_to_aem_out_re, Wv_aem_to_aem_out, Wv_aem_out_to_ext,
                                                           Wv_are, Wv_a_ext, Wv_csm_ext_to_in, Wv_csm_in_to_csm,
                                                           Wv_csm_to_csm_out, Wv_csm_out_to_cgc, Wv_cgc_to_cem_in,
                                                           Wv_cem_in_to_cem, Wv_cem_to_cem_out, Wv_cem_out_to_ext,
                                                           Wv_c_ext, **kwargs):
    """This function calculates the dynamic evolution of the pressure and humidity inside the manifolds.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    T_des : float
        Fuel cell temperature (K).
    Hagc : float
        Thickness of the anode gas channel (m).
    Hcgc : float
        Thickness of the cathode gas channel (m).
    Wagc : float
        Width of the anode gas channel (m).
    Wcgc : float
        Width of the cathode gas channel (m).
    type_auxiliary : str
        Type of auxiliary components used in the fuel cell model.
    Jv_a_in : float
        Water vapor flow at the anode inlet (mol.m-2.s-1).
    Jv_a_out : float
        Water vapor flow at the anode outlet (mol.m-2.s-1).
    Jv_c_in : float
        Water vapor flow at the cathode inlet (mol.m-2.s-1).
    Jv_c_out : float
        Water vapor flow at the cathode outlet (mol.m-2.s-1).
    Wasm_out : float
        Flow at the anode supply manifold outlet (kg.s-1).
    Waem_in : float
        Flow at the anode exhaust manifold inlet (kg.s-1).
    Waem_out_to_ext : float
        Flow at the anode exhaust manifold outlet (kg.s-1).
    Wcsm_out : float
        Flow at the cathode supply manifold outlet (kg.s-1).
    Wcem_in : float
        Flow at the cathode exhaust manifold inlet (kg.s-1).
    Wcem_out_to_ext : float
        Flow at the cathode exhaust manifold outlet (kg.s-1).
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

    # Pressure evolution inside the manifolds
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # At the anode side
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            dif_eq['dPasm_in_re / dt'] = (Ware - Wasm_in_re_to_asm) / V_endplate_a * R * T_des
            dif_eq['dPasm / dt'] = (Wasm_in_to_asm + Wasm_in_re_to_asm - n_cell * Wasm_to_asm_out) / Vasm * R * T_des
            dif_eq['dPaem / dt'] = (n_cell * Waem_in_to_aem - Waem_to_aem_out - Waem_to_aem_out_re) / Vaem * R * T_des
            dif_eq['dPaem_out_re / dt'] = (Waem_to_aem_out_re - Ware) / V_endplate_a * R * T_des
        else: # type_auxiliary == "forced-convective_cathode_with_flow-through_anode"
            dif_eq['dPasm / dt'] = (Wasm_in_to_asm - n_cell * Wasm_to_asm_out) / Vasm * R * T_des
            dif_eq['dPaem / dt'] = (n_cell * Waem_in_to_aem - Waem_to_aem_out) / Vaem * R * T_des
        dif_eq['dPasm_out / dt'] = (Wasm_to_asm_out - Wasm_out_to_agc) / V_man_agc * R * T_des
        dif_eq['dPaem_in / dt'] = (Wagc_to_aem_in - Waem_in_to_aem) / V_man_agc * R * T_des
        dif_eq['dPaem_out / dt'] = (Waem_to_aem_out - Waem_out_to_ext) / V_endplate_a * R * T_des
        # At the cathode side
        dif_eq['dPcsm / dt'] = (Wcsm_in_to_csm - n_cell * Wcsm_to_csm_out) / Vcsm * R * T_des
        dif_eq['dPcsm_out / dt'] = (Wcsm_to_csm_out - Wcsm_out_to_cgc) / V_man_cgc * R * T_des
        dif_eq['dPcem_in / dt'] = (Wcgc_to_cem_in - Wcem_in_to_cem) / V_man_cgc * R * T_des
        dif_eq['dPcem / dt'] = (n_cell * Wcem_in_to_cem - Wcem_to_cem_out) / Vcem * R * T_des
        dif_eq['dPcem_out / dt'] = (Wcem_to_cem_out - Wcem_out_to_ext) / V_endplate_c * R * T_des
    else:  # elif type_auxiliary == "no_auxiliary":
        dif_eq['dPaem_in / dt'] = (Wagc_to_aem_in - Waem_out_to_ext) / V_man_agc * R * T_des
        dif_eq['dPa_ext / dt'] = 0                                                                                      # Boundary condition - Dirichlet condition.
        dif_eq['dPcem_in / dt'] = (Wcgc_to_cem_in - Wcem_out_to_ext) / V_man_cgc * R * T_des
        dif_eq['dPc_ext / dt'] = 0                                                                                      # Boundary condition - Dirichlet condition.

    # Humidity evolution inside the manifolds
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # At the anode side
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            dif_eq['dPhi_asm_in_re / dt'] = (Wv_are - Wv_asm_in_re_to_asm) / V_endplate_a * R * T_des / Psat(T_des)
            dif_eq['dPhi_asm / dt'] = (Wv_asm_in_to_asm + Wv_asm_in_re_to_asm - n_cell * Wv_asm_to_asm_out) / Vasm * R * T_des / Psat(T_des)
            dif_eq['dPhi_aem / dt'] = (n_cell * Wv_aem_in_to_aem - Wv_aem_to_aem_out_re - Wv_aem_to_aem_out) / Vaem * R * T_des / Psat(T_des)
            dif_eq['dPhi_aem_out_re / dt'] = (Wv_aem_to_aem_out_re - Wv_are) / V_endplate_a * R * T_des / Psat(T_des)
        else: # type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            dif_eq['dPhi_asm / dt'] = (Wv_asm_in_to_asm - n_cell * Wv_asm_to_asm_out) / Vasm * R * T_des / Psat(T_des)
            dif_eq['dPhi_aem / dt'] = (n_cell * Wv_aem_in_to_aem - Wv_aem_to_aem_out) / Vaem * R * T_des / Psat(T_des)
        dif_eq['dPhi_asm_in / dt'] = (Wv_asm_ext_to_in - Wv_asm_in_to_asm) / V_endplate_a * R * T_des / Psat(T_des)
        dif_eq['dPhi_asm_out / dt'] = (Wv_asm_to_asm_out - Wv_asm_out_to_agc) / V_man_agc * R * T_des / Psat(T_des)
        dif_eq['dPhi_aem_in / dt'] = (Wv_agc_to_aem_in - Wv_aem_in_to_aem) / V_man_agc * R * T_des / Psat(T_des)
        dif_eq['dPhi_aem_out / dt'] = (Wv_aem_to_aem_out - Wv_aem_out_to_ext) / V_endplate_a * R * T_des / Psat(T_des)
        # At the cathode side
        dif_eq['dPhi_csm_in / dt'] = (Wv_csm_ext_to_in - Wv_csm_in_to_csm) / V_endplate_c * R * T_des / Psat(T_des)
        dif_eq['dPhi_csm / dt'] = (Wv_csm_in_to_csm - n_cell * Wv_csm_to_csm_out) / Vcsm * R * T_des / Psat(T_des)
        dif_eq['dPhi_csm_out / dt'] = (Wv_csm_to_csm_out - Wv_csm_out_to_cgc) / V_man_cgc * R * T_des / Psat(T_des)
        dif_eq['dPhi_cem_in / dt'] = (Wv_cgc_to_cem_in - Wv_cem_in_to_cem) / V_man_cgc * R * T_des / Psat(T_des)
        dif_eq['dPhi_cem / dt'] = None
        dif_eq['dPhi_cem_out / dt'] = (Wv_cem_to_cem_out - Wv_cem_out_to_ext) / V_endplate_c * R * T_des / Psat(T_des)
    else:  # elif type_auxiliary == "no_auxiliary":
        dif_eq['dPhi_asm_out / dt'] = 0                                                                                 # Boundary condition - Dirichlet condition.
        dif_eq['dPhi_aem_in / dt'] = (Wv_agc_to_aem_in - Wv_aem_out_to_ext) / V_man_agc * R * T_des / Psat(T_des)
        dif_eq['dPhi_a_ext / dt'] = (Wv_aem_out_to_ext - Wv_a_ext) / (Hagc * Wagc * Lgc / n_gc) * R * T_des / Psat(T_des)      # Boundary condition: at the exit, flow is isothermal, chemical composition remains unchanged and pressure is constant through space, so humidity is constant through space.
        dif_eq['dPhi_csm_out / dt'] = 0                                                                                 # Boundary condition - Dirichlet condition.
        dif_eq['dPhi_cem_in / dt'] = (Wv_cgc_to_cem_in - Wv_cem_out_to_ext) / V_man_cgc * R * T_des / Psat(T_des)
        dif_eq['dPhi_c_ext / dt'] = (Wv_cem_out_to_ext - Wv_c_ext) / (Hcgc * Wcgc * Lgc / n_gc) * R * T_des / Psat(T_des)      # Boundary condition: at the exit, flow is isothermal, chemical composition remains unchanged and pressure is constant through space, so humidity is constant through space.


def calculate_dC_totdx(C_tot, type_auxiliary, Lman_gc_to_gc):
    dC_totdx = {}
    if type_auxiliary == "no_auxiliary":
        # At the anode side
        dC_totdx['asm_out'] = d_dx(y_minus=C_tot['asm_out'], y_0=None, y_plus=C_tot['agc_1'], dx_minus=Lman_gc_to_gc) # Boundary condition: pressure is constant through space, and the inlet flow is isotherm, so the total concentration is constant through space.
        dC_totdx['a_ext'] = d_dx(y_minus=C_tot['aem_in'], y_0=None, y_plus=C_tot['a_ext'], dx_minus=Lman_gc_to_gc)  # Boundary condition: at the exit, pressure and temperature are fixed, so the total concentration is fixed.
        # At the cathode side
        dC_totdx['csm_out'] = d_dx(y_minus=C_tot['csm_out'], y_0=None, y_plus=C_tot['cgc_1'], dx_minus=Lman_gc_to_gc) # Boundary condition: pressure is constant through space, and the inlet flow is isotherm, so the total concentration is constant through space.
        dC_totdx['c_ext'] = d_dx(y_minus=C_tot['cem_in'], y_0=None, y_plus=C_tot['c_ext'], dx_minus=Lman_gc_to_gc)  # Boundary condition: at the exit, pressure and temperature are fixed, so the total concentration is fixed.
    return dC_totdx


def calculate_d2C_totdx2(C_tot, type_auxiliary, Lman_gc_to_gc):
    d2C_totdx2 = {}
    if type_auxiliary == "no_auxiliary":
        # At the anode side
        d2C_totdx2['asm_out'] = d2_dx2(y_minus=C_tot['asm_out'], y_0=C_tot['asm_out'], y_plus=C_tot['agc_1'],
                                       dx_minus=Lman_gc_to_gc)                                                        # Boundary condition: pressure is constant through space, and the inlet flow is isotherm, so the total concentration is constant through space.
        d2C_totdx2['a_ext'] = d2_dx2(y_minus=C_tot['aem_in'], y_0=C_tot['a_ext'], y_plus=C_tot['a_ext'],
                                     dx_minus=Lman_gc_to_gc)                                                        # Boundary condition: at the exit, pressure and temperature are fixed, so the total concentration is fixed.
        # At the cathode side
        d2C_totdx2['csm_out'] = d2_dx2(y_minus=C_tot['csm_out'], y_0=C_tot['csm_out'], y_plus=C_tot['cgc_1'],
                                       dx_minus=Lman_gc_to_gc)                                                        # Boundary condition: pressure is constant through space, and the inlet flow is isotherm, so the total concentration is constant through space.
        d2C_totdx2['c_ext'] = d2_dx2(y_minus=C_tot['cem_in'], y_0=C_tot['c_ext'], y_plus=C_tot['c_ext'],
                                     dx_minus=Lman_gc_to_gc)                                                        # Boundary condition: at the exit, pressure and temperature are fixed, so the total concentration is fixed.
    return d2C_totdx2


def calculate_dPdx(P, Lgc, n_gc, type_auxiliary, Lman_gc_to_gc):
    dPdx = {}
    if type_auxiliary == "no_auxiliary":
        # At the anode side
        dPdx['asm_out'] = d_dx(y_minus=P['asm_out'], y_0=None, y_plus=P['agc_1'], dx_minus=Lman_gc_to_gc)               # Boundary condition: pressure is constant through space.
        if n_gc == 1:
            dPdx['agc_1'] = d_dx(y_minus=P['asm_out'], y_0=None, y_plus=P['aem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            dPdx['agc_1'] = d_dx(y_minus=P['asm_out'], y_0=P['agc_1'], y_plus=P['agc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            dPdx['agc_2'] = d_dx(y_minus=P['agc_1'], y_0=P['agc_2'], y_plus=P['aem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else: # n_gc > 2
            dPdx['agc_1'] = d_dx(y_minus=P['asm_out'], y_0=P['agc_1'], y_plus=P['agc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                dPdx[f'agc_{i}'] = d_dx(y_minus=P[f'agc_{i - 1}'], y_0=None, y_plus=P[f'agc_{i + 1}'],
                                        dx_minus=Lgc / n_gc)
            dPdx[f'agc_{n_gc}'] = d_dx(y_minus=P[f'agc_{n_gc - 1}'], y_0=P[f'agc_{n_gc}'], y_plus=P['aem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        dPdx['aem_in'] = d_dx(y_minus=P[f'agc_{n_gc}'], y_0=None, y_plus=P['a_ext'], dx_minus=Lman_gc_to_gc)
        dPdx['a_ext'] = d_dx(y_minus=P['aem_in'], y_0=None, y_plus=P['a_ext'], dx_minus=Lman_gc_to_gc)                  # Boundary condition: at the exit, pressure is fixed.
        # At the cathode side
        dPdx['csm_out'] = d_dx(y_minus=P['csm_out'], y_0=None, y_plus=P['cgc_1'], dx_minus=Lman_gc_to_gc)               # Boundary condition: pressure is constant through space.
        if n_gc == 1:
            dPdx['cgc_1'] = d_dx(y_minus=P['csm_out'], y_0=None, y_plus=P['cem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            dPdx['cgc_1'] = d_dx(y_minus=P['csm_out'], y_0=P['cgc_1'], y_plus=P['cgc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            dPdx['cgc_2'] = d_dx(y_minus=P['cgc_1'], y_0=P['cgc_2'], y_plus=P['cem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else: # n_gc > 2
            dPdx['cgc_1'] = d_dx(y_minus=P['csm_out'], y_0=P['cgc_1'], y_plus=P['cgc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                dPdx[f'cgc_{i}'] = d_dx(y_minus=P[f'cgc_{i - 1}'], y_0=None, y_plus=P[f'cgc_{i + 1}'],
                                           dx_minus=Lgc / n_gc)
            dPdx[f'cgc_{n_gc}'] = d_dx(y_minus=P[f'cgc_{n_gc - 1}'], y_0=P[f'cgc_{n_gc}'], y_plus=P['cem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        dPdx['cem_in'] = d_dx(y_minus=P[f'cgc_{n_gc}'], y_0=None, y_plus=P['c_ext'], dx_minus=Lman_gc_to_gc)
        dPdx['c_ext'] = d_dx(y_minus=P['cem_in'], y_0=None, y_plus=P['c_ext'], dx_minus=Lman_gc_to_gc)              # Boundary condition: at the exit, pressure is fixed.
    return dPdx


def calculate_dmudx(Lgc, n_gc, type_auxiliary, Lman_gc_to_gc, mu):
    dmudx = {}
    if type_auxiliary == "no_auxiliary":
        # At the anode side
        dmudx['asm_out'] = d_dx(y_minus=mu['asm_out'], y_0=None, y_plus=mu['agc_1'], dx_minus=Lman_gc_to_gc)          # Boundary condition: at the entry, temperature and gas composition are fixed, so the dynamic viscosity is fixed.
        if n_gc == 1:
            dmudx['agc_1'] = d_dx(y_minus=mu['asm_out'], y_0=None, y_plus=mu['aem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            dmudx['agc_1'] = d_dx(y_minus=mu['asm_out'], y_0=mu['agc_1'], y_plus=mu['agc_2'],
                                  dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            dmudx['agc_2'] = d_dx(y_minus=mu['agc_1'], y_0=mu['agc_2'], y_plus=mu['aem_in'],
                                  dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else: # n_gc > 2
            dmudx['agc_1'] = d_dx(y_minus=mu['asm_out'], y_0=mu['agc_1'], y_plus=mu['agc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                dmudx[f'agc_{i}'] = d_dx(y_minus=mu[f'agc_{i - 1}'], y_0=None, y_plus=mu[f'agc_{i + 1}'],
                                        dx_minus=Lgc / n_gc)
            dmudx[f'agc_{n_gc}'] = d_dx(y_minus=mu[f'agc_{n_gc - 1}'], y_0=mu[f'agc_{n_gc}'], y_plus=mu['aem_in'],
                                       dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        dmudx['aem_in'] = d_dx(y_minus=mu[f'agc_{n_gc}'], y_0=None, y_plus=mu['a_ext'], dx_minus=Lman_gc_to_gc)
        dmudx['a_ext'] = d_dx(y_minus=mu['aem_in'], y_0=None, y_plus=mu['a_ext'], dx_minus=Lman_gc_to_gc)           # Boundary condition: at the exit, temperature and gas composition are fixed, so the dynamic viscosity is fixed.
        # At the cathode side
        dmudx['csm_out'] = d_dx(y_minus=mu['csm_out'], y_0=None, y_plus=mu['cgc_1'], dx_minus=Lman_gc_to_gc)          # Boundary condition: at the entry, temperature and gas composition are fixed, so the dynamic viscosity is fixed.
        if n_gc == 1:
            dmudx['cgc_1'] = d_dx(y_minus=mu['csm_out'], y_0=None, y_plus=mu['cem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            dmudx['cgc_1'] = d_dx(y_minus=mu['csm_out'], y_0=mu['cgc_1'], y_plus=mu['cgc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            dmudx['cgc_2'] = d_dx(y_minus=mu['cgc_1'], y_0=mu['cgc_2'], y_plus=mu['cem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else: # n_gc > 2
            dmudx['cgc_1'] = d_dx(y_minus=mu['csm_out'], y_0=mu['cgc_1'], y_plus=mu['cgc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                dmudx[f'cgc_{i}'] = d_dx(y_minus=mu[f'cgc_{i - 1}'], y_0=None, y_plus=mu[f'cgc_{i + 1}'],
                                           dx_minus=Lgc / n_gc)
            dmudx[f'cgc_{n_gc}'] = d_dx(y_minus=mu[f'cgc_{n_gc - 1}'], y_0=mu[f'cgc_{n_gc}'], y_plus=mu['cem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        dmudx['cem_in'] = d_dx(y_minus=mu[f'cgc_{n_gc}'], y_0=None, y_plus=mu['c_ext'], dx_minus=Lman_gc_to_gc)
        dmudx['c_ext'] = d_dx(y_minus=mu['cem_in'], y_0=None, y_plus=mu['c_ext'], dx_minus=Lman_gc_to_gc)           # Boundary condition: at the exit, temperature and gas composition are fixed, so the dynamic viscosity is fixed.
    return dmudx


def calculate_dvdx(dif_eq, v, T_des, Lgc, n_gc, type_auxiliary, Lman_gc_to_gc, C_tot, dC_totdx):
    dvdx = {}
    if type_auxiliary == "no_auxiliary":
        # At the anode side
        dC_totdt_asm_out = dif_eq['dPasm_out / dt'] / (R * T_des)                                                       # True as long as the flow is isothermal.
        dvdx['asm_out'] = - v['asm_out'] / C_tot['asm_out'] * dC_totdx['asm_out'] - dC_totdt_asm_out / C_tot['asm_out'] # Obtained from the molar balance differential equation, knowing dif_eq['dPcsm_out / dt'].
        if n_gc == 1:
            dvdx['agc_1'] = d_dx(y_minus=v['asm_out'], y_0=None, y_plus=v['aem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            dvdx['agc_1'] = d_dx(y_minus=v['asm_out'], y_0=v['agc_1'], y_plus=v['agc_2'],
                                  dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            dvdx['agc_2'] = d_dx(y_minus=v['agc_1'], y_0=v['agc_2'], y_plus=v['aem_in'],
                                  dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else:  # n_gc > 2
            dvdx['agc_1'] = d_dx(y_minus=v['asm_out'], y_0=v['agc_1'], y_plus=v['agc_2'],
                                  dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                dvdx[f'agc_{i}'] = d_dx(y_minus=v[f'agc_{i - 1}'], y_0=None, y_plus=v[f'agc_{i + 1}'],
                                         dx_minus=Lgc / n_gc)
            dvdx[f'agc_{n_gc}'] = d_dx(y_minus=v[f'agc_{n_gc - 1}'], y_0=v[f'agc_{n_gc}'], y_plus=v['aem_in'],
                                        dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        dvdx['aem_in'] = d_dx(y_minus=v[f'agc_{n_gc}'], y_0=None, y_plus=v['a_ext'], dx_minus=Lman_gc_to_gc)
        dvdx['a_ext'] = - v['a_ext'] / C_tot['a_ext'] * dC_totdx['a_ext']                                               # Obtained from the boundary condition C_tot_a_ext = cte (fixed pressure and temperature at the outlet), injected in the molar balance differential equation.
        # At the cathode side
        dC_totdt_csm_out = dif_eq['dPcsm_out / dt'] / (R * T_des)                                                       # True as long as the flow is isothermal.
        dvdx['csm_out'] = - v['csm_out'] / C_tot['csm_out'] * dC_totdx['csm_out'] - dC_totdt_csm_out / C_tot['csm_out'] # Obtained from the molar balance differential equation, knowing dif_eq['dPcsm_out / dt'].
        if n_gc == 1:
            dvdx['cgc_1'] = d_dx(y_minus=v['csm_out'], y_0=None, y_plus=v['cem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            dvdx['cgc_1'] = d_dx(y_minus=v['csm_out'], y_0=v['cgc_1'], y_plus=v['cgc_2'],
                                  dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            dvdx['cgc_2'] = d_dx(y_minus=v['cgc_1'], y_0=v['cgc_2'], y_plus=v['cem_in'],
                                  dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else:  # n_gc > 2
            dvdx['cgc_1'] = d_dx(y_minus=v['csm_out'], y_0=v['cgc_1'], y_plus=v['cgc_2'],
                                  dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                dvdx[f'cgc_{i}'] = d_dx(y_minus=v[f'cgc_{i - 1}'], y_0=None, y_plus=v[f'cgc_{i + 1}'],
                                         dx_minus=Lgc / n_gc)
            dvdx[f'cgc_{n_gc}'] = d_dx(y_minus=v[f'cgc_{n_gc - 1}'], y_0=v[f'cgc_{n_gc}'], y_plus=v['cem_in'],
                                        dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        dvdx['cem_in'] = d_dx(y_minus=v[f'cgc_{n_gc}'], y_0=None, y_plus=v['c_ext'], dx_minus=Lman_gc_to_gc)
        dvdx['c_ext'] = - v['c_ext'] / C_tot['c_ext'] * dC_totdx['c_ext']                                               # Obtained from the boundary condition C_tot_c_ext = cte (fixed pressure and temperature at the outlet), injected in the molar balance differential equation.
    return dvdx


def calculate_d2vdx2(dif_eq, v, T_des, Lgc, n_gc, type_auxiliary, Lman_gc_to_gc, C_tot, dC_totdx, dvdx):
    d2vdx2 = {}
    d2C_totdx2 = calculate_d2C_totdx2(C_tot, type_auxiliary, Lman_gc_to_gc)
    if type_auxiliary == "no_auxiliary":
        # At the anode side
        dC_totdt_asm_out = dif_eq['dPasm_out / dt'] / (R * T_des)                                                       # True as long as the flow is isothermal.
        d2vdx2['asm_out'] = - (2 * dC_totdx['asm_out'] * dvdx['asm_out'] + v['asm_out'] * d2C_totdx2['asm_out']
                               + dC_totdt_asm_out) / C_tot['asm_out']                                                   # Obtained from the molar balance differential equation, knowing dif_eq['dPcsm_out / dt'].
        if n_gc == 1:
            d2vdx2['agc_1'] = d2_dx2(y_minus=v['asm_out'], y_0=v['agc_1'], y_plus=v['aem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            d2vdx2['agc_1'] = d2_dx2(y_minus=v['asm_out'], y_0=v['agc_1'], y_plus=v['agc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            d2vdx2['agc_2'] = d2_dx2(y_minus=v['agc_1'], y_0=v['agc_2'], y_plus=v['aem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else:  # n_gc > 2
            d2vdx2['agc_1'] = d2_dx2(y_minus=v['asm_out'], y_0=v['agc_1'], y_plus=v['agc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                d2vdx2[f'agc_{i}'] = d2_dx2(y_minus=v[f'agc_{i - 1}'], y_0=v[f'agc_{i}'], y_plus=v[f'agc_{i + 1}'],
                                        dx_minus=Lgc / n_gc)
            d2vdx2[f'agc_{n_gc}'] = d2_dx2(y_minus=v[f'agc_{n_gc - 1}'], y_0=v[f'agc_{n_gc}'], y_plus=v['aem_in'],
                                       dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        d2vdx2['aem_in'] = d2_dx2(y_minus = v[f'agc_{n_gc}'], y_0 = v['aem_in'], y_plus = v['a_ext'], dx_minus = Lman_gc_to_gc)
        d2vdx2['a_ext'] = - (2 * dC_totdx['a_ext'] * dvdx['a_ext'] + v['a_ext'] * d2C_totdx2['a_ext']) / C_tot['a_ext'] # Obtained from the boundary condition C_tot_a_ext = cte (fixed pressure and temperature at the outlet), injected in the molar balance differential equation.
        # At the cathode side
        dC_totdt_csm_out = dif_eq['dPcsm_out / dt'] / (R * T_des)                                                       # True as long as the flow is isothermal.
        d2vdx2['csm_out'] = - (2 * dC_totdx['csm_out'] * dvdx['csm_out'] + v['csm_out'] * d2C_totdx2['csm_out']
                               + dC_totdt_csm_out) / C_tot['csm_out']                                                   # Obtained from the molar balance differential equation, knowing dif_eq['dPcsm_out / dt'].
        if n_gc == 1:
            d2vdx2['cgc_1'] = d2_dx2(y_minus=v['csm_out'], y_0=v['cgc_1'], y_plus=v['cem_in'], dx_minus=Lman_gc_to_gc)
        elif n_gc == 2:
            d2vdx2['cgc_1'] = d2_dx2(y_minus=v['csm_out'], y_0=v['cgc_1'], y_plus=v['cgc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            d2vdx2['cgc_2'] = d2_dx2(y_minus=v['cgc_1'], y_0=v['cgc_2'], y_plus=v['cem_in'],
                                 dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        else:  # n_gc > 2
            d2vdx2['cgc_1'] = d2_dx2(y_minus=v['csm_out'], y_0=v['cgc_1'], y_plus=v['cgc_2'],
                                 dx_minus=Lman_gc_to_gc, dx_plus=Lgc / n_gc)
            for i in range(2, n_gc):
                d2vdx2[f'cgc_{i}'] = d2_dx2(y_minus=v[f'cgc_{i - 1}'], y_0=v[f'cgc_{i}'], y_plus=v[f'cgc_{i + 1}'],
                                        dx_minus=Lgc / n_gc)
            d2vdx2[f'cgc_{n_gc}'] = d2_dx2(y_minus=v[f'cgc_{n_gc - 1}'], y_0=v[f'cgc_{n_gc}'], y_plus=v['cem_in'],
                                       dx_minus=Lgc / n_gc, dx_plus=Lman_gc_to_gc)
        d2vdx2['cem_in'] = d2_dx2(y_minus = v[f'cgc_{n_gc}'], y_0 = v['cem_in'], y_plus = v['c_ext'], dx_minus = Lman_gc_to_gc)
        d2vdx2['c_ext'] = - (2 * dC_totdx['c_ext'] * dvdx['c_ext'] + v['c_ext'] * d2C_totdx2['c_ext']) / C_tot['c_ext'] # Obtained from the boundary condition C_tot_c_ext = cte (fixed pressure and temperature at the outlet), injected in the molar balance differential equation.
    return d2vdx2


def calculate_dyn_velocity_evolution(dif_eq, sv, T_des, Pa_des, Lgc, Hagc, Hcgc, n_gc, type_auxiliary,
                                     Lman_to_endplate, Lman_to_man_gc, Lman_gc_to_gc, v_re, k_purge, rho, C_tot, mu_gaz,
                                     P, Jv_agc_agdl, Jv_cgdl_cgc, J_H2_agc_agdl, J_O2_cgdl_cgc, **kwargs):

    # Boundary conditions:
    #   - At the inlet:
    #       - for the anode if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
    #           - for asm_in:
    #               - Dirichlet boundary condition: dv/dx = 0.
    #           - for asm_in_re:
    #               - Neumann boundary condition: dP/dx = 0. The pump exactly compensates the pressure drop.
    #               - Dirichlet boundary condition: v = v_re.
    #       - for the cathode, or for the anode in other scenarios:
    #           - Dirichlet boundary condition: v = v_in.
    #   - At the outlet:
    #       - if type_auxiliary != "no_auxiliary":
    #           - for aem_out_re:
    #               - Neumann boundary condition: dP/dx = 0. The pump exactly compensates the pressure drop.
    #               - Dirichlet boundary condition: v = v_re.
    #           - for aem_out et cem_out:
    #               - if A_bp > 0 or k_purge == 1:
    #                 - Dirichlet boundary condition: P = Pext.
    #               - if A_bp = 0 or k_purge == 0:
    #                 - Dirichlet boundary condition: v = 0.
    #       - if type_auxiliary == "no_auxiliary":
    #           - Dirichlet boundary condition: P = P_des.

    # Preliminary calculations
    v_asm_in_re, v_asm_in, v_asm = sv.get('v_asm_in_re', None), sv.get('v_asm_in', None), sv.get('v_asm', None)
    v_asm_out, v_aem_in, v_aem = sv['v_asm_out'], sv['v_aem_in'], sv.get('v_aem', None)
    v_aem_out, v_aem_out_re, v_a_ext = sv.get('v_aem_out', None), sv.get('v_aem_out_re', None), sv['v_a_ext']
    v_csm_in, v_csm, v_csm_out = sv.get('v_csm_in', None), sv.get('v_csm', None), sv['v_csm_out']
    v_cem_in, v_cem, v_cem_out, v_c_ext = sv['v_cem_in'], sv.get('v_cem', None), sv.get('v_cem_out', None), sv['v_c_ext']
    Pasm_in_re, Pasm_in, Pasm = sv.get('Pasm_in_re', None), sv.get('Pasm_in', None), sv.get('Pasm', None)
    Pasm_out, Paem_in, Paem, Paem_out = sv['Pasm_out'], sv['Paem_in'], sv.get('Paem', None), sv.get('Paem_out', None)
    Paem_out_re, Pa_ext, Pcsm_in, Pcsm= sv.get('Paem_out_re', None), sv['Pa_ext'], sv.get('Pcsm_in', None), sv.get('Pcsm', None)
    Pcsm_out, Pcem_in , Pcem, Pcem_out, Pc_ext = sv['Pcsm_out'], sv['Pcem_in'], sv.get('Pcem', None), sv.get('Pcem_out', None), sv['Pc_ext']
    Abp_a, Abp_c = sv.get('Abp_a', None), sv.get('Abp_c', None)

    v = {'asm_in_re': v_asm_in_re, 'asm_in': v_asm_in, 'asm': v_asm, 'asm_out': v_asm_out, 'aem_in': v_aem_in,
         'aem': v_aem, 'aem_out': v_aem_out, 'aem_out_re': v_aem_out_re, 'a_ext': v_a_ext, 'csm_in': v_csm_in,
         'csm': v_csm, 'csm_out': v_csm_out, 'cem_in': v_cem_in, 'cem': v_cem, 'cem_out': v_cem_out, 'c_ext': v_c_ext}
    for i in range(1, n_gc + 1):
        v[f'agc_{i}'], v[f'cgc_{i}'] = sv[f'v_agc_{i}'], sv[f'v_cgc_{i}']
    P.update({'asm_in_re': Pasm_in_re, 'asm_in': Pasm_in, 'asm': Pasm, 'asm_out': Pasm_out, 'aem_in': Paem_in,
              'aem': Paem, 'aem_out': Paem_out, 'aem_out_re': Paem_out_re, 'a_ext': Pa_ext, 'csm_in': Pcsm_in,
              'csm': Pcsm, 'csm_out': Pcsm_out, 'cem_in': Pcem_in, 'cem': Pcem, 'cem_out': Pcem_out, 'c_ext': Pc_ext})
    S_dif = {} # Diffusion term
    for i in range(1, n_gc + 1):
        # S_dif[f'agc_{i}'] = - (Jv_agc_agdl[i] + J_H2_agc_agdl[i]) / Hagc
        # S_dif[f'cgc_{i}'] =   (Jv_cgdl_cgc[i] + J_O2_cgdl_cgc[i]) / Hcgc
        S_dif[f'agc_{i}'] = 0                                                                                           # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        S_dif[f'cgc_{i}'] = 0                                                                                           # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dC_totdx = calculate_dC_totdx(C_tot, type_auxiliary, Lman_gc_to_gc)
    dPdx = calculate_dPdx(P, Lgc, n_gc, type_auxiliary, Lman_gc_to_gc)
    dmudx = calculate_dmudx(Lgc, n_gc, type_auxiliary, Lman_gc_to_gc, mu_gaz)
    dvdx = calculate_dvdx(dif_eq, v, T_des, Lgc, n_gc, type_auxiliary, Lman_gc_to_gc, C_tot, dC_totdx)
    d2vdx2 = calculate_d2vdx2(dif_eq, v, T_des, Lgc, n_gc, type_auxiliary, Lman_gc_to_gc, C_tot, dC_totdx, dvdx)

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
        dif_eq['dv_asm_in_re / dt'] = - d_dx(Pasm_in_re, Pasm_in_re, Pasm, Lman_to_endplate/2) / rho['asm_in_re'] - \
                                        v_asm_in_re * d_dx(v_re, v_asm_in_re, v_asm, Lman_to_endplate/2) + \
                                        mu_gaz['asm_in_re'] / rho['asm_in_re'] * d2_dx2(v_re, v_asm_in_re, v_asm, Lman_to_endplate/2)
        dif_eq['dv_aem_out_re / dt'] = - d_dx(Paem, Paem_out_re, Paem_out_re, Lman_to_endplate/2) / rho['aem_out_re'] - \
                                        v_aem_out_re * d_dx(v_aem, v_aem_out_re, v_re, Lman_to_endplate/2) + \
                                        mu_gaz['aem_out_re'] / rho['aem_out_re'] * d2_dx2(v_aem, v_aem_out_re, v_re, Lman_to_endplate/2)

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # At the anode side
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            dif_eq['dv_asm_in / dt'] = - d_dx(Pa_des, Pasm_in, Pasm, Lman_to_endplate / 2) / rho['asm_in'] - \
                                         v_asm_in * d_dx(v_asm_in, v_asm_in, v_asm_out, Lman_to_endplate / 2) + \
                                         mu_gaz['asm_in'] / rho['asm_in'] * d2_dx2(v_asm_in, v_asm_in, v_asm_out, Lman_to_endplate / 2)
        else: # type_auxiliary == "forced-convective_cathode_with_flow-through_anode"
            dif_eq['dv_asm_in / dt'] = - d_dx(Pasm_in, Pasm_in, Pasm, Lman_to_endplate/2) / rho['asm_in'] - \
                                         v_asm_in * d_dx(v_a_in, v_asm_in, v_asm_out, Lman_to_endplate/2) + \
                                         mu_gaz['asm_in'] / rho['asm_in'] * d2_dx2(v_a_in, v_asm_in, v_asm_out, Lman_to_endplate/2)
        dif_eq['dv_asm / dt'] = - d_dx(Pasm_in, Pasm, Pasm_out, Lman_to_endplate/2, Lman_to_man_gc/2) / rho['asm'] - \
                                  v_asm * d_dx(v_asm_in, v_asm, v_asm_out, Lman_to_endplate/2, Lman_to_man_gc/2) + \
                                  mu_gaz['asm'] / rho['asm'] * d2_dx2(v_asm_in, v_asm, v_asm_out, Lman_to_endplate/2, Lman_to_man_gc/2)
        dif_eq['dv_asm_out / dt'] = - d_dx(Pasm, Pasm_out, Pagc, Lman_to_man_gc/2, Lman_gc_to_gc) / rho['asm_out'] - \
                                      v_asm_out * d_dx(v_asm, v_asm_out, v_agc, Lman_to_man_gc/2, Lman_gc_to_gc) + \
                                      mu_gaz['asm_out'] / rho['asm_out'] * d2_dx2(v_asm, v_asm_out, v_agc, Lman_to_man_gc/2, Lman_gc_to_gc)
        dif_eq['dv_agc / dt'] = - d_dx(Pasm_out, Pagc, Paem_in, Lman_gc_to_gc) / rho['agc'] - \
                                  v_agc * d_dx(v_asm_out, v_agc, v_aem_in, Lman_gc_to_gc) + \
                                  mu_gaz['agc'] / rho['agc'] * d2_dx2(v_asm_out, v_agc, v_aem_in, Lman_gc_to_gc)
        dif_eq['dv_aem_in / dt'] = - d_dx(Pagc, Paem_in, Paem, Lman_gc_to_gc, Lman_to_man_gc/2) / rho['aem_in'] - \
                                     v_aem_in * d_dx(v_agc, v_aem_in, v_aem, Lman_gc_to_gc, Lman_to_man_gc/2) + \
                                     mu_gaz['aem_in'] / rho['aem_in'] * d2_dx2(v_agc, v_aem_in, v_aem, Lman_gc_to_gc, Lman_to_man_gc/2)
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            if Abp_a > 0: # The valve is open which impacts the boundary conditions at the outlet.
                dif_eq['dv_aem / dt'] = - d_dx(Paem_in, Paem, Paem_out, Lman_to_man_gc / 2, Lman_to_endplate / 2) / rho['aem'] - \
                                        v_aem * d_dx(v_aem_in, v_aem, v_aem_out, Lman_to_man_gc / 2, Lman_to_endplate / 2) + \
                                        mu_gaz['aem'] / rho['aem'] * d2_dx2(v_aem_in, v_aem, v_aem_out, Lman_to_man_gc / 2, Lman_to_endplate / 2)
                dif_eq['dv_aem_out / dt'] = - d_dx(Paem, Paem_out, Pext, Lman_to_endplate/2) / rho['aem_out'] - \
                                              v_aem_out * d_dx(v_aem, v_aem_out, v_aem_out, Lman_to_endplate/2) + \
                                              mu_gaz['aem_out'] / rho['aem_out'] * d2_dx2(v_aem, v_aem_out, v_aem_out, Lman_to_endplate/2)
            else: # Abp_a == 0. The valve is closed which impacts the boundary conditions at the outlet.
                dif_eq['dv_aem / dt'] = - d_dx(Paem_in, Paem, Paem_out_re, Lman_to_man_gc / 2, Lman_to_endplate / 2) / rho['aem'] - \
                                        v_aem * d_dx(v_aem_in, v_aem, v_aem_out_re, Lman_to_man_gc / 2, Lman_to_endplate / 2) + \
                                        mu_gaz['aem'] / rho['aem'] * d2_dx2(v_aem_in, v_aem, v_aem_out_re, Lman_to_man_gc / 2, Lman_to_endplate / 2)
                dif_eq['dv_aem_out / dt'] = - d_dx(Paem, Paem_out, Paem_out, Lman_to_endplate/2) / rho['aem_out'] - \
                                              v_aem_out * d_dx(v_aem, v_aem_out, 0, Lman_to_endplate/2) + \
                                              mu_gaz['aem_out'] / rho['aem_out'] * d2_dx2(v_aem, v_aem_out, 0, Lman_to_endplate/2)
        else: # type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            dif_eq['dv_aem / dt'] = - d_dx(Paem_in, Paem, Paem_out, Lman_to_man_gc / 2, Lman_to_endplate / 2) / rho['aem'] - \
                                    v_aem * d_dx(v_aem_in, v_aem, v_aem_out, Lman_to_man_gc / 2, Lman_to_endplate / 2) + \
                                    mu_gaz['aem'] / rho['aem'] * d2_dx2(v_aem_in, v_aem, v_aem_out, Lman_to_man_gc / 2, Lman_to_endplate / 2)
            if k_purge == 1: # The valve is open which impacts the boundary conditions at the outlet.
                dif_eq['dv_aem_out / dt'] = - d_dx(Paem, Paem_out, Pext, Lman_to_endplate/2) / rho['aem_out'] - \
                                              v_aem_out * d_dx(v_aem, v_aem_out, v_aem_out, Lman_to_endplate/2) + \
                                              mu_gaz['aem_out'] / rho['aem_out'] * d2_dx2(v_aem, v_aem_out, v_aem_out, Lman_to_endplate/2)
            else: # k_purge == 0. The valve is closed which impacts the boundary conditions at the outlet.
                dif_eq['dv_aem_out / dt'] = - d_dx(Paem, Paem_out, Paem_out, Lman_to_endplate/2) / rho['aem_out'] - \
                                              v_aem_out * d_dx(v_aem, v_aem_out, 0, Lman_to_endplate/2) + \
                                              mu_gaz['aem_out'] / rho['aem_out'] * d2_dx2(v_aem, v_aem_out, 0, Lman_to_endplate/2)
        # At the cathode side
        dif_eq['dv_csm_in / dt'] = - d_dx(Pcsm_in, Pcsm_in, Pcsm, Lman_to_endplate/2) / rho['csm_in'] - \
                                     v_csm_in * d_dx(v_c_in, v_csm_in, v_csm_out, Lman_to_endplate/2) + \
                                     mu_gaz['csm_in'] / rho['csm_in'] * d2_dx2(v_c_in, v_csm_in, v_csm_out, Lman_to_endplate/2)
        dif_eq['dv_csm / dt'] = - d_dx(Pcsm_in, Pcsm, Pcsm_out, Lman_to_endplate/2, Lman_to_man_gc/2) / rho['csm'] - \
                                  v_csm * d_dx(v_csm_in, v_csm, v_csm_out, Lman_to_endplate/2, Lman_to_man_gc/2) + \
                                  mu_gaz['csm'] / rho['csm'] * d2_dx2(v_csm_in, v_csm, v_csm_out, Lman_to_endplate/2, Lman_to_man_gc/2)
        dif_eq['dv_csm_out / dt'] = - d_dx(Pcsm, Pcsm_out, Pcgc, Lman_to_man_gc/2, Lman_gc_to_gc) / rho['csm_out'] - \
                                      v_csm_out * d_dx(v_csm, v_csm_out, Lman_to_man_gc/2, Lman_gc_to_gc) + \
                                      mu_gaz['csm_out'] / rho['csm_out'] * d2_dx2(v_csm, v_csm_out, v_cgc, Lman_to_man_gc/2, Lman_gc_to_gc)
        dif_eq['dv_cgc / dt'] = - d_dx(Pcsm_out, Pcgc, Pcem_in, Lman_gc_to_gc) / rho['cgc'] - \
                                  v_cgc * d_dx(v_csm_out, v_cgc, v_cem_in, Lman_gc_to_gc) + \
                                  mu_gaz['cgc'] / rho['cgc'] * d2_dx2(v_csm_out, v_cgc, v_cem_in, Lman_gc_to_gc)
        dif_eq['dv_cem_in / dt'] = - d_dx(Pcgc, Pcem_in, Pcem, Lman_gc_to_gc, Lman_to_man_gc/2) / rho['cem_in'] - \
                                    v_cem_in * d_dx(v_cgc, v_cem_in, v_cem, Lman_gc_to_gc, Lman_to_man_gc/2) + \
                                    mu_gaz['cem_in'] / rho['cem_in'] * d2_dx2(v_cgc, v_cem_in, v_cem, Lman_gc_to_gc, Lman_to_man_gc/2)
        dif_eq['dv_cem / dt'] = - d_dx(Pcem_in, Pcem, Pcem_out, Lman_to_man_gc/2, Lman_to_endplate/2) / rho['cem'] - \
                                  v_cem * d_dx(v_cem_in, v_cem, v_cem_out, Lman_to_man_gc/2, Lman_to_endplate/2) + \
                                  mu_gaz['cem'] / rho['cem'] * d2_dx2(v_cem_in, v_cem, v_cem_out, Lman_to_man_gc/2, Lman_to_endplate/2)
        if Abp_c > 0:  # The valve is open which impacts the boundary conditions at the outlet.
            dif_eq['dv_cem_out / dt'] = - d_dx(Pcem, Pcem_out, Pext, Lman_to_endplate/2) / rho['cem_out'] - \
                                          v_cem_out * d_dx(v_cem, v_cem_out, v_cem_out, Lman_to_endplate/2) + \
                                          mu_gaz['cem_out'] / rho['cem_out'] * d2_dx2(v_cem, v_cem_out, v_cem_out, Lman_to_endplate/2)
        else:  # Abp_c == 0 The valve is closed which impacts the boundary conditions at the outlet.
            dif_eq['dv_cem_out / dt'] = - d_dx(Pcem, Pcem_out, Pcem_out, Lman_to_endplate/2) / rho['cem_out'] - \
                                        v_cem_out * d_dx(v_cem, v_cem_out, 0, Lman_to_endplate/2) + \
                                        mu_gaz['cem_out'] / rho['cem_out'] * d2_dx2(v_cem, v_cem_out, 0, Lman_to_endplate/2)
    else: # type_auxiliary == "no_auxiliary"
        for key in ['asm_out', 'agc', 'aem_in', 'a_ext', 'csm_out', 'cgc', 'cem_in', 'c_ext']:
            if key in ['agc', 'cgc']: # Diffusion only occurs in the gas channels.
                for i in range(1, n_gc + 1):
                    key_i = f"{key}_{i}"
                    dif_eq['dv_' + key_i + ' / dt'] = - dPdx[key_i] / rho[key_i] + \
                                                      (-v[key_i] + 4/3 * 1/rho[key_i] * dmudx[key_i]) * dvdx[key_i] + \
                                                      4/3 * mu_gaz[key_i] / rho[key_i] * d2vdx2[key_i] + \
                                                      - S_dif[key_i] * v[key_i] / rho[key_i]
            else:
                dif_eq['dv_' + key + ' / dt'] = - dPdx[key] / rho[key] + \
                                                (-v[key] + 4/3 * 1/rho[key] * dmudx[key]) * dvdx[key] + \
                                                4/3 * mu_gaz[key] / rho[key] * d2vdx2[key]