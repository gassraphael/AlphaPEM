# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the auxiliaries flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import Text, Pext, Phi_ext, M_H2, M_O2, M_N2, M_H2O, y_O2_ext, R, F, A_T
from modules.transitory_functions import average, Psat, C_v_sat, k_H2, k_O2


# _________________________________________________Auxiliaries modules__________________________________________________

def auxiliaries_int_values(t, solver_variables, operating_inputs, parameters):
    """This functions calculates intermediate values for the auxiliaries flows calculation.

    Parameters
    ----------
    t : float
        Time (s).
    solver_variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell model.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    Mext : float
        Molar mass of the ambient air outside the stack (kg/mol).
    Pagc : float
        Global pressure in the anode gas channel (Pa).
    Pcgc : float
        Global pressure in the cathode gas channel (Pa).
    Phi_agc : float
        Relative humidity in the anode gas channel.
    Phi_cgc : float
        Relative humidity in the cathode gas channel.
    y_H2_agc : float
        Hydrogen ratio in the dry anode gas mixture (H2/N2) injected at the inlet.
    y_O2_cgc : float
        Oxygen ratio in dry air in the cathode gas channel.
    Magc : float
        Molar mass of all the gas species in the anode gas channel (kg/mol).
    Mcgc : float
        Molar mass of all the gas species in the cathode gas channel (kg/mol).
    Pr_aem : float
        Pressure ratio in the anode external manifold.
    Pr_cem : float
        Pressure ratio in the cathode external manifold.
    Maem : float
        Molar mass of all the gas species in the anode external manifold (kg/mol).
    Masm : float
        Molar mass of all the gas species in the anode supply manifold (kg/mol).
    Mcem : float
        Molar mass of all the gas species in the cathode external manifold (kg/mol).
    Mcsm : float
        Molar mass of all the gas species in the cathode supply manifold (kg/mol).
    k_purge : float
        Purge coefficient. It is equal to 1 if the purge is active and 0 otherwise.
    Abp_a : float
        Area of the back pressure valve in the anode external manifold (m²).
    Abp_c : float
        Area of the back pressure valve in the cathode external manifold (m²).
    i_n : float
        Internal current density (A/m²).
    """

    # Extraction of the variables
    C_v_agc, C_v_cgc = solver_variables['C_v_agc'], solver_variables['C_v_cgc']
    lambda_mem = solver_variables['lambda_mem']
    C_H2_agc, C_H2_acl = solver_variables['C_H2_agc'], solver_variables['C_H2_acl']
    C_O2_ccl, C_O2_cgc = solver_variables['C_O2_ccl'], solver_variables['C_O2_cgc']
    C_N2_a, C_N2_c = solver_variables['C_N2_a'], solver_variables['C_N2_c']
    T_acl, T_mem, T_ccl = solver_variables['T_acl'], solver_variables['T_mem'], solver_variables['T_ccl']
    Pasm, Paem, Pcsm = solver_variables['Pasm'], solver_variables['Paem'], solver_variables['Pcsm']
    Pcem, Phi_asm, Phi_aem = solver_variables['Pcem'], solver_variables['Phi_asm'], solver_variables['Phi_aem']
    Phi_csm, Phi_cem = solver_variables['Phi_csm'], solver_variables['Phi_cem']
    Abp_a, Abp_c = solver_variables['Abp_a'], solver_variables['Abp_c']
    # Extraction of the operating inputs and the parameters
    T_des, y_H2_in = operating_inputs['T_des'], operating_inputs['y_H2_in']
    Hmem, Hacl, Hccl, kappa_co = parameters['Hmem'], parameters['Hacl'], parameters['Hccl'], parameters['kappa_co']
    t_purge, type_purge = parameters['t_purge'], parameters['type_purge']

    # Molar mass of the ambient air
    Mext = Phi_ext * Psat(Text) / Pext * M_H2O + \
           y_O2_ext * (1 - Phi_ext * Psat(Text) / Pext) * M_O2 + \
           (1 - y_O2_ext) * (1 - Phi_ext * Psat(Text) / Pext) * M_N2
    M_H2_N2_in = y_H2_in * M_H2 + (1 - y_H2_in) * M_N2

    # Physical quantities inside the stack
    #       Pressures
    Pagc = (C_v_agc + C_H2_agc + C_N2_a) * R * T_des
    Pcgc = (C_v_cgc + C_O2_cgc + C_N2_c) * R * T_des
    #       Humidities
    Phi_agc = C_v_agc / C_v_sat(T_des)
    Phi_cgc = C_v_cgc / C_v_sat(T_des)
    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_H2_agc = C_H2_agc / (C_H2_agc + C_N2_a)
    y_O2_cgc = C_O2_cgc / (C_O2_cgc + C_N2_c)
    #       Molar masses
    Magc = C_v_agc * R * T_des / Pagc * M_H2O + \
           C_H2_agc * R * T_des / Pagc * M_H2 + \
           C_N2_a * R * T_des / Pagc * M_N2
    Mcgc = Phi_cgc * Psat(T_des) / Pcgc * M_H2O + \
           y_O2_cgc * (1 - Phi_cgc * Psat(T_des) / Pcgc) * M_O2 + \
           (1 - y_O2_cgc) * (1 - Phi_cgc * Psat(T_des) / Pcgc) * M_N2

    # Internal current density
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
    i_n = i_H2 + i_O2

    # Physical quantities in the auxiliary system
    if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation" or \
       parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
        # Pressure ratios
        Pr_aem = (Pext / Paem)
        Pr_cem = (Pext / Pcem)
        # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the EM
        y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_a * R * T_des) / (Paem - Phi_aem * Psat(T_des))
        y_O2_cem = (Pcem - Phi_cem * Psat(T_des) - C_N2_c * R * T_des) / (Pcem - Phi_cem * Psat(T_des))
        # Molar masses
        if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation":
            Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
            Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
        else: # parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode"
            Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 + \
                   (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
            Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 + \
                   (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
        Mcem = Phi_cem * Psat(T_des) / Pcem * M_H2O + \
               y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 + \
               (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2
        Mcsm = Phi_csm * Psat(T_des) / Pcsm * M_H2O + \
               y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 + \
               (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
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
        if Abp_a > A_T:
            Abp_a = A_T
        elif Abp_a < 0:
            Abp_a = 0
        if Abp_c > A_T:
            Abp_c = A_T
        elif Abp_c < 0:
            Abp_c = 0
    else:  # parameters["type_auxiliary"] == "no_auxiliary"
        Pr_aem, Pr_cem, Maem, Masm, Mcem, Mcsm, k_purge, Abp_a, Abp_c = [0] * 9

    return (Mext, M_H2_N2_in, Pagc, Pcgc, Phi_agc, Phi_cgc, y_H2_agc, y_O2_cgc, Magc, Mcgc, Pr_aem, Pr_cem, Maem, Masm, Mcem, Mcsm,
            k_purge, Abp_a, Abp_c, i_n)
