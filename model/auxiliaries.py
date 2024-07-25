# -*- coding: utf-8 -*-

"""This file represents all the flows passing through the auxiliaries. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import F, R, Text, Pext, Phi_ext, yO2_ext, gamma, gamma_H2, M_H2, M_H2O, n_cell, A_T, \
    Ksm_in, Ksm_out, Kem_in, C_D
from modules.transitory_functions import Psat
from modules.auxiliaries_modules import auxiliaries_int_values


# ______________________________________________________Auxiliaries_____________________________________________________

def auxiliaries(t, solver_variables, control_variables, i_fc, operating_inputs, parameters):
    """This function calculates the flows passing through the auxiliaries.

    Parameters
    ----------
    t : float
        Time (s).
    solver_variables : dict
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
    J_N2_in : float
        N2 flow at the cathode gas channel inlet (mol.m-2.s-1).
    J_N2_out : float
        N2 flow at the cathode gas channel outlet (mol.m-2.s-1).
    Wasm_in : float
        Anode side supply manifold inlet flow (kg.s-1).
    Wasm_out : float
        Anode side supply manifold outlet flow (kg.s-1).
    Waem_in : float
        Anode side external manifold inlet flow (kg.s-1).
    Waem_out : float
        Anode side external manifold outlet flow (kg.s-1).
    Wcsm_in : float
        Cathode side supply manifold inlet flow (kg.s-1).
    Wcsm_out : float
        Cathode side supply manifold outlet flow (kg.s-1).
    Wcem_in : float
        Cathode side external manifold inlet flow (kg.s-1).
    Wcem_out : float
        Cathode side external manifold outlet flow (kg.s-1).
    Ware : float
        Anode side recirculation flow (kg.s-1).
    Wv_asm_in : float
        Vapor flow at the anode supply manifold inlet (mol.s-1).
    Wv_aem_out : float
        Vapor flow at the anode external manifold outlet (mol.s-1).
    Wv_csm_in : float
        Vapor flow at the cathode supply manifold inlet (mol.s-1).
    Wv_cem_out : float
        Vapor flow at the cathode external manifold outlet (mol.s-1).
    """

    # __________________________________________________Preliminaries___________________________________________________

    # Extraction of the variables
    Pasm, Paem, Pcsm = solver_variables['Pasm'], solver_variables['Paem'], solver_variables['Pcsm']
    Pcem, Phi_asm, Phi_aem = solver_variables['Pcem'], solver_variables['Phi_asm'], solver_variables['Phi_aem']
    Phi_csm, Phi_cem = solver_variables['Phi_csm'], solver_variables['Phi_cem']
    Wcp, Wa_inj, Wc_inj = solver_variables['Wcp'], solver_variables['Wa_inj'], solver_variables['Wc_inj']

    # Extraction of the operating inputs and the parameters
    Tfc, Pa_des, Pc_des = operating_inputs['Tfc'], operating_inputs['Pa_des'], operating_inputs['Pc_des']
    Sa, Sc = operating_inputs['Sa'], operating_inputs['Sc']
    Phi_a_des, Phi_c_des = control_variables['Phi_a_des'], control_variables['Phi_c_des']
    Aact, Hgc, Wgc = parameters['Aact'], parameters['Hgc'], parameters['Wgc']
    type_auxiliary = parameters['type_auxiliary']

    # Intermediate values
    Mext, Pagc, Pcgc, Phi_agc, Phi_cgc, y_cgc, Magc, Mcgc, Pr_aem, Pr_cem, \
        Maem, Masm, Mcem, Mcsm, k_purge, Abp_a, Abp_c, i_n = \
        auxiliaries_int_values(t, solver_variables, operating_inputs, parameters)

    # _________________________________________Inlet and outlet global flows____________________________________________
    """Global flows here refer to flows that integrate all the chemical species circulating together.
    Slight differences are to be noted in the expression of these flows depending on the type of auxiliary selected.
    """

    # At the anode side
    if type_auxiliary == "closed_anode_with_recirculation":
        # Anode inlet
        Wasm_in = Ksm_in * (Pa_des - Pasm)  # kg.s-1
        Wasm_out = Ksm_out * (Pasm - Pagc)  # kg.s-1
        Ja_in = Wasm_out / (Hgc * Wgc * Masm)  # mol.m-2.s-1
        # Anode outlet
        Waem_in = Kem_in * (Pagc - Paem)  # kg.s-1
        Ware = n_cell * Maem * (Paem / (Paem - Phi_aem * Psat(Tfc))) * (Sa - 1) * (i_fc + i_n) / (
                    2 * F) * Aact  # kg.s-1
        Waem_out = k_purge * C_D * A_T * Paem / np.sqrt(R * Tfc) * Pr_aem ** (1 / gamma_H2) * \
                   np.sqrt(Magc * 2 * gamma_H2 / (gamma_H2 - 1) * (1 - Pr_aem ** ((gamma_H2 - 1) / gamma_H2)))  # kg.s-1
        Ja_out = Waem_in / (Hgc * Wgc * Magc)  # mol.m-2.s-1

    elif type_auxiliary == "opened_anode":
        # Anode inlet
        Wrd = n_cell * M_H2 * Sa * (i_fc + i_n) / (2 * F) * Aact  # kg.s-1
        Wasm_in = Wrd + Wa_inj  # kg.s-1
        Wasm_out = Ksm_out * (Pasm - Pagc)  # kg.s-1
        Ja_in = Wasm_out / (Hgc * Wgc * Masm)  # mol.m-2.s-1
        # Anode outlet
        Waem_in = Kem_in * (Pagc - Paem)  # kg.s-1
        Ware = 0  # kg.s-1
        Waem_out = C_D * Abp_a * Paem / np.sqrt(R * Tfc) * Pr_aem ** (1 / gamma_H2) * \
                   np.sqrt(Magc * 2 * gamma_H2 / (gamma_H2 - 1) * (1 - Pr_aem ** ((gamma_H2 - 1) / gamma_H2)))
        # kg.s-1
        Ja_out = Waem_in / (Hgc * Wgc * Magc)  # mol.m-2.s-1

    else:  # elif type_auxiliary == "no_auxiliary" (only 1 cell):
        # Anode inlet
        Wasm_in, Wasm_out = 0, 0  # kg.s-1
        Ja_in = (1 + Phi_a_des * Psat(Tfc) / (Pagc - Phi_a_des * Psat(Tfc))) * \
                Sa * (i_fc + i_n) / (2 * F) * Aact / (Hgc * Wgc)  # mol.m-2.s-1
        # Anode outlet
        Waem_in, Ware, Waem_out = 0, 0, 0  # kg.s-1
        Ja_out = Kem_in * (Pagc - Pa_des) / (Hgc * Wgc * Magc)  # mol.m-2.s-1

    # At the cathode side
    if type_auxiliary == "closed_anode_with_recirculation" or type_auxiliary == "opened_anode":
        # Cathode inlet         
        Wcsm_in = Wcp + Wc_inj  # kg.s-1
        Wcsm_out = Ksm_out * (Pcsm - Pcgc)  # kg.s-1
        Jc_in = Wcsm_out / (Hgc * Wgc * Mcsm)  # mol.m-2.s-1
        # Cathode outlet
        Wcem_in = Kem_in * (Pcgc - Pcem)  # kg.s-1
        Wcem_out = C_D * Abp_c * Pcem / np.sqrt(R * Tfc) * Pr_cem ** (1 / gamma) * \
                   np.sqrt(Mcgc * 2 * gamma / (gamma - 1) * (1 - Pr_cem ** ((gamma - 1) / gamma)))  # kg.s-1
        Jc_out = Wcem_in / (Hgc * Wgc * Mcgc)  # mol.m-2.s-1

    else:  # elif type_auxiliary == "no_auxiliary" (only 1 cell):
        # Cathode inlet   
        Wcsm_in, Wcsm_out = 0, 0  # kg.s-1
        Jc_in = (1 + Phi_c_des * Psat(Tfc) / (Pcgc - Phi_c_des * Psat(Tfc))) * \
                1 / yO2_ext * Sc * (i_fc + i_n) / (4 * F) * Aact / (Hgc * Wgc)  # mol.m-2.s-1
        # Cathode outlet
        Wcem_in, Wcem_out = 0, 0  # kg.s-1
        Jc_out = Kem_in * (Pcgc - Pc_des) / (Hgc * Wgc * Mcgc)  # mol.m-2.s-1

    # ________________________________________Inlet and outlet specific flows___________________________________________
    """Specific flows here refer to flows that integrate only a single chemical species within the ensemble of species
     circulating together. For example, only the water vapor flow within the ensemble of hydrogen and water vapor.
     """

    # Vapor flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "closed_anode_with_recirculation" or type_auxiliary == "opened_anode":
        Jv_a_in = Phi_asm * Psat(Tfc) / Pasm * Ja_in
        Jv_c_in = Phi_csm * Psat(Tfc) / Pcsm * Jc_in
    else:  # elif type_auxiliary == "no_auxiliary":
        Jv_a_in = Phi_a_des * Psat(Tfc) / Pagc * Ja_in
        Jv_c_in = Phi_c_des * Psat(Tfc) / Pcgc * Jc_in
    Jv_a_out = Phi_agc * Psat(Tfc) / Pagc * Ja_out
    Jv_c_out = Phi_cgc * Psat(Tfc) / Pcgc * Jc_out

    # H2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "closed_anode_with_recirculation" or type_auxiliary == "opened_anode":
        J_H2_in = (1 - Phi_asm * Psat(Tfc) / Pasm) * Ja_in
    else:  # elif type_auxiliary == "no_auxiliary":
        J_H2_in = (1 - Phi_a_des * Psat(Tfc) / Pagc) * Ja_in
    J_H2_out = (1 - Phi_agc * Psat(Tfc) / Pagc) * Ja_out

    # O2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "closed_anode_with_recirculation" or type_auxiliary == "opened_anode":
        J_O2_in = yO2_ext * (1 - Phi_csm * Psat(Tfc) / Pcsm) * Jc_in
    else:  # elif type_auxiliary == "no_auxiliary":
        J_O2_in = yO2_ext * (1 - Phi_c_des * Psat(Tfc) / Pcgc) * Jc_in
    J_O2_out = y_cgc * (1 - Phi_cgc * Psat(Tfc) / Pcgc) * Jc_out

    # N2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == "closed_anode_with_recirculation" or type_auxiliary == "opened_anode":
        J_N2_in = (1 - yO2_ext) * (1 - Phi_csm * Psat(Tfc) / Pcsm) * Jc_in
    else:  # elif type_auxiliary == "no_auxiliary":
        J_N2_in = (1 - yO2_ext) * (1 - Phi_c_des * Psat(Tfc) / Pcgc) * Jc_in
    J_N2_out = (1 - y_cgc) * (1 - Phi_cgc * Psat(Tfc) / Pcgc) * Jc_out

    # Vapor flows at the manifold (mol.s-1)
    if type_auxiliary == "closed_anode_with_recirculation":
        Wv_asm_in = Phi_aem * Psat(Tfc) / Paem * (Ware / Maem)
        Wv_aem_out = Phi_aem * Psat(Tfc) / Paem * (Waem_out / Maem)
        Wv_csm_in = Phi_ext * Psat(Text) / Pext * (Wcp / Mext) + Wc_inj / M_H2O
        Wv_cem_out = Phi_cem * Psat(Tfc) / Pcem * (Wcem_out / Mcem)
    elif type_auxiliary == "opened_anode":
        Wv_asm_in = Wa_inj / M_H2O
        Wv_aem_out = Phi_aem * Psat(Tfc) / Paem * (Waem_out / Maem)
        Wv_csm_in = Phi_ext * Psat(Text) / Pext * (Wcp / Mext) + Wc_inj / M_H2O
        Wv_cem_out = Phi_cem * Psat(Tfc) / Pcem * (Wcem_out / Mcem)
    else:  # elif type_auxiliary == "no_auxiliary":
        Wv_asm_in, Wv_aem_out, Wv_csm_in, Wv_cem_out = [0] * 4

    return Jv_a_in, Jv_a_out, Jv_c_in, Jv_c_out, \
        J_H2_in, J_H2_out, J_O2_in, J_O2_out, J_N2_in, J_N2_out, \
        Wasm_in, Wasm_out, Waem_in, Waem_out, Wcsm_in, Wcsm_out, Wcem_in, Wcem_out, Ware, \
        Wv_asm_in, Wv_aem_out, Wv_csm_in, Wv_cem_out
