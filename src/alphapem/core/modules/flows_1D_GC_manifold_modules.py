# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the auxiliaries flows calculation.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries

# Importing constants' value and functions
from alphapem.utils.physics_constants import Text, Pext, Phi_ext, M_H2, M_O2, M_N2, M_H2O, y_O2_ext, R
from alphapem.utils.physics_functions import Psat, C_v_sat, mu_mixture_gases


# _________________________________________________Auxiliaries modules__________________________________________________

def flow_1D_GC_manifold_int_values(sv_1D_cell, sv_auxiliary, operating_inputs, parameters):
    """This functions calculates intermediate values for the auxiliaries flows calculation.

    Parameters
    ----------
    sv_1D_cell : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    sv_auxiliary : dict
        Variables calculated by the auxiliary system. They correspond to the auxiliary system internal states.
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

    # Extraction of the operating inputs and the parameters
    T_des, y_H2_in, Pa_des = operating_inputs['T_des'], operating_inputs['y_H2_in'], operating_inputs['Pa_des']
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    kappa_co, nb_gc, t_purge, type_purge = parameters['kappa_co'], parameters['nb_gc'], parameters['t_purge'], parameters['type_purge']
    # Extraction of the variables
    Abp_a, Abp_c = sv_auxiliary.get('Abp_a', None), sv_auxiliary.get('Abp_c', None)
    C_v_agc = [None] + [sv_1D_cell[i]['C_v_agc'] for i in range(1, nb_gc + 1)]
    C_v_cgc = [None] + [sv_1D_cell[i]['C_v_cgc'] for i in range(1, nb_gc + 1)]
    C_H2_agc = [None] + [sv_1D_cell[i]['C_H2_agc'] for i in range(1, nb_gc + 1)]
    C_O2_cgc = [None] + [sv_1D_cell[i]['C_O2_cgc'] for i in range(1, nb_gc + 1)]
    C_N2_agc = [None] + [sv_1D_cell[i]['C_N2_agc'] for i in range(1, nb_gc + 1)]
    C_N2_cgc = [None] + [sv_1D_cell[i]['C_N2_cgc'] for i in range(1, nb_gc + 1)]
    T_agc = [None] + [sv_1D_cell[i]['T_agc'] for i in range(1, nb_gc + 1)]
    T_cgc = [None] + [sv_1D_cell[i]['T_cgc'] for i in range(1, nb_gc + 1)]

    # Physical quantities outside the stack
    #       Molar masses
    M_ext = Phi_ext * Psat(Text) / Pext * M_H2O + \
               y_O2_ext * (1 - Phi_ext * Psat(Text) / Pext) * M_O2 + \
               (1 - y_O2_ext) * (1 - Phi_ext * Psat(Text) / Pext) * M_N2
    M_H2_N2_in = y_H2_in * M_H2 + (1 - y_H2_in) * M_N2

    # Physical quantities inside the stack
    #       Pressures
    P_agc = [None]
    P_cgc = [None]
    for i in range(1, nb_gc + 1):
        P_agc += [(C_v_agc[i] + C_H2_agc[i] + C_N2_agc[i]) * R * T_agc[i]]
        P_cgc += [(C_v_cgc[i] + C_O2_cgc[i] + C_N2_cgc[i]) * R * T_cgc[i]]
    #       Humidities
    Phi_agc = [None]
    Phi_cgc = [None]
    for i in range(1, nb_gc + 1):
        Phi_agc += [C_v_agc[i] / C_v_sat(T_agc[i])]
        Phi_cgc += [C_v_cgc[i] / C_v_sat(T_cgc[i])]
    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_O2_cgc = [None]
    y_H2_agc = [None]
    for i in range(1, nb_gc + 1):
        y_H2_agc += [C_H2_agc[i] / (C_H2_agc[i] + C_N2_agc[i])]
        y_O2_cgc += [C_O2_cgc[i] / (C_O2_cgc[i] + C_N2_cgc[i])]
    #       Molar masses
    M_agc = [None]
    M_cgc = [None]
    for i in range(1, nb_gc + 1):
        M_agc += [C_v_agc[i] * R * T_des / P_agc[i] * M_H2O + \
                        C_H2_agc[i] * R * T_des / P_agc[i] * M_H2 + \
                        C_N2_agc[i] * R * T_des / P_agc[i] * M_N2]
        M_cgc += [Phi_cgc[i] * Psat(T_des) / P_cgc[i] * M_H2O + \
                        y_O2_cgc[i] * (1 - Phi_cgc[i] * Psat(T_des) / P_cgc[i]) * M_O2 + \
                        (1 - y_O2_cgc[i]) * (1 - Phi_cgc[i] * Psat(T_des) / P_cgc[i]) * M_N2]
    #       Density of the gas mixture.
    rho_agc = [None]
    rho_cgc = [None]
    for i in range(1, nb_gc + 1):
        rho_agc += [P_agc[i] / (R * T_agc[i]) * M_agc[i]]
        rho_cgc += [P_cgc[i] / (R * T_cgc[i]) * M_cgc[i]]

    #       Vapor ratio over the gas mixture.
    x_H2O_v_agc = [None]
    x_H2O_v_cgc = [None]
    for i in range(1, nb_gc + 1):
        x_H2O_v_agc += [C_v_agc[i] / (C_v_agc[i] + C_H2_agc[i] + C_N2_agc[i])]
    for i in range(1, nb_gc + 1):
        x_H2O_v_cgc += [C_v_cgc[i] / (C_v_cgc[i] + C_O2_cgc[i] + C_N2_cgc[i])]

    #       Dynamic viscosity of the gas mixture.
    mu_gaz_agc = [None]
    mu_gaz_cgc = [None]
    for i in range(1, nb_gc + 1):
        mu_gaz_agc += [mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_agc[i], 1 - x_H2O_v_agc[i]], T_agc[i])]
    for i in range(1, nb_gc + 1):
        mu_gaz_cgc += [mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                              [x_H2O_v_cgc[i], y_O2_cgc[i] * (1 - x_H2O_v_cgc[i]),
                                               (1 - y_O2_cgc[i]) * (1 - x_H2O_v_cgc[i])], T_cgc[i])]

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

    return (P_agc, P_cgc, Phi_agc, Phi_cgc, y_H2_agc, y_O2_cgc, M_agc, M_cgc, M_ext, M_H2_N2_in, rho_agc, rho_cgc,
            k_purge, Abp_a, Abp_c, mu_gaz_agc, mu_gaz_cgc)