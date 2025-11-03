# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the auxiliaries flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import Text, Pext, Phi_ext, M_H2, M_O2, M_N2, M_H2O, y_O2_ext, R, F
from modules.transitory_functions import Psat, C_v_sat, mu_mixture_gases


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
    C_H2_acl, C_O2_ccl =  sv['C_H2_acl'], sv['C_O2_ccl']
    Pasm, Paem, Pcsm, Pcem = sv.get('Pasm', None), sv.get('Paem', None), sv.get('Pcsm', None), sv.get('Pcem', None)
    Phi_asm, Phi_aem = sv.get('Phi_asm', None), sv.get('Phi_aem', None)
    Phi_csm, Phi_cem = sv.get('Phi_csm', None), sv.get('Phi_cem', None)
    Abp_a, Abp_c = sv.get('Abp_a', None), sv.get('Abp_c', None)
    # Extraction of the operating inputs and the parameters
    T_des, y_H2_in, Pa_des = operating_inputs['T_des'], operating_inputs['y_H2_in'], operating_inputs['Pa_des']
    Lgc, A_T_a, A_T_c = parameters['Lgc'], parameters['A_T_a'], parameters['A_T_c']
    kappa_co, nb_gc, t_purge, type_purge = parameters['kappa_co'], parameters['nb_gc'], parameters['t_purge'], parameters['type_purge']

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
    for i in range(1, nb_gc + 1):
        P[f'agc_{i}'] = (sv[f'C_v_agc_{i}'] + sv[f'C_H2_agc_{i}'] + sv[f'C_N2_agc_{i}']) * R * sv[f'T_agc_{i}']
        P[f'cgc_{i}'] = (sv[f'C_v_cgc_{i}'] + sv[f'C_O2_cgc_{i}'] + sv[f'C_N2_cgc_{i}']) * R * sv[f'T_cgc_{i}']
    #       Humidities
    Phi = {}
    for i in range(1, nb_gc + 1):
        Phi[f'agc_{i}'] = sv[f'C_v_agc_{i}'] / C_v_sat(sv[f'T_agc_{i}'])
        Phi[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] / C_v_sat(sv[f'T_cgc_{i}'])
    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_O2 = {}
    y_H2 = {}
    for i in range(1, nb_gc + 1):
        y_H2[f'agc_{i}'] = sv[f'C_H2_agc_{i}'] / (sv[f'C_H2_agc_{i}'] + sv[f'C_N2_agc_{i}'])
        y_O2[f'cgc_{i}'] = sv[f'C_O2_cgc_{i}'] / (sv[f'C_O2_cgc_{i}'] + sv[f'C_N2_cgc_{i}'])
    #       Molar masses
    for i in range(1, nb_gc + 1):
        M[f'agc_{i}'] = sv[f'C_v_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2O + \
                        sv[f'C_H2_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2 + \
                        sv[f'C_N2_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_N2
        M[f'cgc_{i}'] = Phi[f'cgc_{i}'] * Psat(T_des) / P[f'cgc_{i}'] * M_H2O + \
                        y_O2[f'cgc_{i}'] * (1 - Phi[f'cgc_{i}'] * Psat(T_des) / P[f'cgc_{i}']) * M_O2 + \
                        (1 - y_O2[f'cgc_{i}']) * (1 - Phi[f'cgc_{i}'] * Psat(T_des) / P[f'cgc_{i}']) * M_N2
    #       Density of the gas mixture.
    rho = {}
    for i in range(1, nb_gc + 1):
        rho[f'agc_{i}'] = P[f'agc_{i}'] / (R * sv[f'T_agc_{i}']) * M[f'agc_{i}']
        rho[f'cgc_{i}'] = P[f'cgc_{i}'] / (R * sv[f'T_cgc_{i}']) * M[f'cgc_{i}']

    #       Vapor ratio over the gas mixture.
    x_H2O_v = {}
    for i in range(1, nb_gc + 1):
        x_H2O_v[f'agc_{i}'] = sv[f'C_v_agc_{i}'] / (sv[f'C_v_agc_{i}'] + sv[f'C_H2_agc_{i}'] + sv[f'C_N2_agc_{i}'])
    for i in range(1, nb_gc + 1):
        x_H2O_v[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] / (sv[f'C_v_cgc_{i}'] + sv[f'C_O2_cgc_{i}'] + sv[f'C_N2_cgc_{i}'])

    #       Dynamic viscosity of the gas mixture.
    mu_gaz = {}
    for i in range(1, nb_gc + 1):
        mu_gaz[f'agc_{i}'] = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v[f'agc_{i}'], 1 - x_H2O_v[f'agc_{i}']],
                                              sv[f'T_agc_{i}'])
    for i in range(1, nb_gc + 1):
        mu_gaz[f'cgc_{i}'] = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                              [x_H2O_v[f'cgc_{i}'], y_O2[f'cgc_{i}'] * (1 - x_H2O_v[f'cgc_{i}']),
                                               (1 - y_O2[f'cgc_{i}']) * (1 - x_H2O_v[f'cgc_{i}'])],
                                              sv[f'T_cgc_{i}'])

    # Physical quantities in the auxiliary system
    if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation" or \
       parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
        pass
        # # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the EM
        # y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_a * R * T_des) / (Paem - Phi_aem * Psat(T_des))
        # y_O2_cem = (Pcem - Phi_cem * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem - Phi_cem * Psat(T_cgc))
        #
        # # Molar masses
        # if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation":
        #     Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
        #            (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
        #     Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
        #            (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
        # else:  # parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
        #     Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
        #            y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 + \
        #            (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
        #     Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
        #            y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 + \
        #            (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
        # # Molar masses at the cathode side
        # Mcsm = Phi_csm * Psat(T_des) / Pcsm * M_H2O + \
        #        y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 + \
        #        (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
        # Mcem = Phi_cem * Psat(T_des) / Pcem * M_H2O + \
        #        y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 + \
        #        (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2
        #
        # # Density of the gas mixture.
        # rho_asm = Pasm / (R * T_des) * Masm
        # rho_aem = Paem / (R * T_des) * Maem
        # rho_csm = Pcsm / (R * T_des) * Mcsm
        # rho_cgc = Pcgc / (R * T_cgc) * Mcgc
        # rho_cem = Pcem / (R * T_cgc) * Mcem
        #
        # # Purge
        # if type_purge == "no_purge":
        #     k_purge = 0
        # elif type_purge == "constant_purge":
        #     k_purge = 1
        # elif type_purge == "periodic_purge":
        #     purge_time, delta_purge = t_purge
        #     if (t - int(t / (purge_time + delta_purge)) * (purge_time + delta_purge)) <= purge_time:
        #         k_purge = 1
        #     else:
        #         k_purge = 0
        # else:
        #     raise ValueError("The type_purge variable should be correctly referenced.")
        # # Back pressure valve area
        # if Abp_a > A_T_a:
        #     Abp_a = A_T_a
        # elif Abp_a < 0:
        #     Abp_a = 0
        # if Abp_c > A_T_c:
        #     Abp_c = A_T_c
        # elif Abp_c < 0:
        #     Abp_c = 0

    else:  # parameters["type_auxiliary"] == "no_auxiliary"
        k_purge, Abp_a, Abp_c = [None] * 3

    return P, Phi, y_H2, y_O2, M, rho, k_purge, Abp_a, Abp_c, mu_gaz