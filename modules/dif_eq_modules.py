# -*- coding: utf-8 -*-

"""This module is used to determine intermediate values for the calculation of the differential equations
and to implement integration events.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import Text, Pext, Phi_ext, M_H2, M_O2, M_N2, M_H2O, y_O2_ext, R
from modules.transitory_functions import Psat, C_v_sat, mu_mixture_gases, calculate_rho_Cp0


# ____________________________________________Differential equations modules____________________________________________

def calculate_dif_eq_int_values(t, sv, control_variables, i_fc, operating_inputs, parameters):
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
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']
    Pasm, Paem, Pcsm, Pcem = sv.get('Pasm', None), sv.get('Paem', None), sv.get('Pcsm', None), sv.get('Pcem', None)
    Phi_asm, Phi_aem = sv.get('Phi_asm', None), sv.get('Phi_aem', None)
    Phi_csm, Phi_cem = sv.get('Phi_csm', None), sv.get('Phi_cem', None)

    # Extraction of the operating inputs and the parameters
    T_des, y_H2_in = operating_inputs['T_des'], operating_inputs['y_H2_in']
    Lgc, nb_channel_in_gc, Lm = parameters['Lgc'], parameters['nb_channel_in_gc'], parameters['Lm']
    Hccl, epsilon_gdl, epsilon_cl = parameters['Hccl'], parameters['epsilon_gdl'], parameters['epsilon_cl']
    epsilon_mpl, kappa_co, epsilon_mc = parameters['epsilon_mpl'], parameters['kappa_co'], parameters['epsilon_mc']
    epsilon_atl, epsilon_ctl = parameters['epsilon_atl'], parameters['epsilon_ctl']
    Htl, nb_gc, nb_gdl, nb_mpl, nb_tl = parameters['Htl'], parameters['nb_gc'], parameters['nb_gdl'], parameters['nb_mpl'], parameters['nb_tl']
    t_purge, type_auxiliary, type_purge = parameters['t_purge'], parameters['type_auxiliary'], parameters['type_purge']

    # Calculation of intermediate values
    C_N2_a_mean = (sum(sv[f'C_N2_agc_{i}'] for i in range(1, nb_gc + 1)) / nb_gc)
    C_N2_c_mean = (sum(sv[f'C_N2_cgc_{i}'] for i in range(1, nb_gc + 1)) / nb_gc)

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
    #       Total concentration
    C_tot = {}
    for i in range(1, nb_gc + 1):
        C_tot[f'agc_{i}'] = sv[f'C_v_agc_{i}'] + sv[f'C_H2_agc_{i}'] + sv[f'C_N2_agc_{i}']
        C_tot[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] + sv[f'C_O2_cgc_{i}'] + sv[f'C_N2_cgc_{i}']

    #       Humidities
    Phi = {}
    for i in range(1, nb_gc + 1):
        Phi[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] / C_v_sat(sv[f'T_cgc_{i}'])

    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_O2 = {}
    for i in range(1, nb_gc + 1):
        y_O2[f'cgc_{i}'] = sv[f'C_O2_cgc_{i}'] / (sv[f'C_O2_cgc_{i}'] + sv[f'C_N2_cgc_{i}'])

    #       Molar masses
    for i in range(1, nb_gc + 1):
        M[f'agc_{i}'] = sv[f'C_v_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2O + \
               sv[f'C_H2_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_H2 + \
               sv[f'C_N2_agc_{i}'] * R * T_des / P[f'agc_{i}'] * M_N2
        M[f'cgc_{i}'] = sv[f'C_v_cgc_{i}'] * R * T_des / P[f'cgc_{i}'] * M_H2O + \
                        sv[f'C_O2_cgc_{i}'] * R * T_des / P[f'cgc_{i}'] * M_O2 + \
                        sv[f'C_N2_cgc_{i}'] * R * T_des / P[f'cgc_{i}'] * M_N2

    #       Density of the gas mixture.
    rho = {}
    for i in range(1, nb_gc + 1):
        rho[f'agc_{i}'] = P[f'agc_{i}'] / (R * sv[f'T_agc_{i}']) * M[f'agc_{i}']
    for i in range(1, nb_gc + 1):
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


    #       Volumetric heat capacity (J.m-3.K-1)
    rho_Cp0 = {
        **{f'agdl_{i}': calculate_rho_Cp0('agdl', sv[f'T_agdl_{i}'], C_v=sv[f'C_v_agdl_{i}'],
                                          s=sv[f's_agdl_{i}'], C_H2=sv[f'C_H2_agdl_{i}'], C_N2=C_N2_a_mean, epsilon=epsilon_gdl)
           for i in range(1, nb_gdl + 1)},
        **{f'atl_{i}': calculate_rho_Cp0('atl', sv[f'T_atl_{i}'], C_v=sv[f'C_v_atl_{i}'], s=sv[f's_atl_{i}'],
                                         C_H2=sv[f'C_H2_atl_{i}'], C_N2=C_N2_a_mean, epsilon=epsilon_atl[i], n_tl=nb_tl,
                                         Htl=Htl, node=i)
           for i in range(1, nb_tl + 1)},
        **{f'ampl_{i}': calculate_rho_Cp0('ampl', sv[f'T_ampl_{i}'], C_v=sv[f'C_v_ampl_{i}'],
                                          s=sv[f's_ampl_{i}'], C_H2=sv[f'C_H2_ampl_{i}'], C_N2=C_N2_a_mean, epsilon=epsilon_mpl)
           for i in range(1, nb_mpl + 1)},
        'acl': calculate_rho_Cp0('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl, C_N2=C_N2_a_mean, C_H2=C_H2_acl,
                                 epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
        'mem': calculate_rho_Cp0('mem', T_mem, lambdaa=lambda_mem),
        'ccl': calculate_rho_Cp0('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl, C_O2=C_O2_ccl, C_N2=C_N2_c_mean,
                                 epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
        **{f'cmpl_{i}': calculate_rho_Cp0('cmpl', sv[f'T_cmpl_{i}'], C_v=sv[f'C_v_cmpl_{i}'],
                                          s=sv[f's_cmpl_{i}'], C_O2=sv[f'C_O2_cmpl_{i}'], C_N2=C_N2_c_mean, epsilon=epsilon_mpl)
           for i in range(1, nb_mpl + 1)},
        **{f'ctl_{i}': calculate_rho_Cp0('ctl', sv[f'T_ctl_{i}'], C_v=sv[f'C_v_ctl_{i}'], s=sv[f's_ctl_{i}'],
                                         C_O2=sv[f'C_O2_ctl_{i}'], C_N2=C_N2_c_mean, epsilon=epsilon_ctl[i], n_tl=nb_tl,
                                         Htl=Htl, node=i)
           for i in range(1, nb_tl + 1)},
        **{f'cgdl_{i}': calculate_rho_Cp0('cgdl', sv[f'T_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                                          s=sv[f's_cgdl_{i}'], C_O2=sv[f'C_O2_cgdl_{i}'], C_N2=C_N2_c_mean, epsilon=epsilon_gdl)
           for i in range(1, nb_gdl + 1)}
        }



    # Physical quantities inside the auxiliary system
    if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation" or \
       parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
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
        y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_a * R * T_des) / (Paem - Phi_aem * Psat(T_des))
        y_O2_cem = (Pcem - Phi_cem * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem - Phi_cem * Psat(T_cgc))

        # Molar masses at the anode side
        if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation":
            M['asm'] = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
            M['aem'] = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
        else: #parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
            M['asm'] = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 + \
                   (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
            M['aem'] = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 + \
                   (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
        # Molar masses at the cathode side
        M['csm'] = Phi_csm * Psat(T_des) / Pcsm * M_H2O + \
               y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 + \
               (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
        M['cem'] = Phi_cem * Psat(T_des) / Pcem * M_H2O + \
               y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 + \
               (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2

        # Density/concentration of the gas mixture.
        C_tot_a_in = Pasm_in / (R * T_des)
        rho_asm = Pasm / (R * T_des) * Masm
        rho_agc = P[f'agc_{i}'] / (R * sv[f'T_agc_{i}']) * Magc
        rho_aem = Paem / (R * T_des) * Maem
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            rho_asm_in_re = Pasm_in_re / (R * T_des) * Masm_in_re
            rho_aem_out_re = Paem_out_re / (R * T_des) * Maem_out_re
        else:
            rho_asm_in_re, rho_aem_out_re = None, None
        rho_a_ext = Pext / (R * T_des) * Maem_out
        C_tot_a_ext = Pext / (R * T_des)                                                                                # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed.
        C_tot_c_in = Pcsm_in / (R * T_des)
        rho_csm = Pcsm / (R * T_des) * Mcsm
        rho_cgc = P[f'cgc_{i}'] / (R * sv[f'T_cgc_{i}']) * Mcgc
        rho_cem = Pcem / (R * T_cgc) * Mcem
        rho_c_ext = Pext / (R * T_des) * Mcem_out
        C_tot_c_ext = Pext * Mcem_out / (R * T_des)                                                                     # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed.

        # Vapor ratio over the gas mixture.
        x_H2O_v_asm = Phi_asm * Psat(T_des) / Pasm
        x_H2O_v_agc = C_v_agc / (C_v_agc + C_H2_agc + C_N2_a)
        x_H2O_v_aem = Phi_aem * Psat(T_des) / Paem
        x_H2O_v_a_ext = Phi_a_ext * Psat(T_des) / Pext
        x_H2O_v_csm = Phi_csm * Psat(T_des) / Pcsm
        x_H2O_v_cgc = C_v_cgc / (C_v_cgc + C_O2_cgc + C_N2_c)
        x_H2O_v_cem = Phi_cem * Psat(T_des) / Pcem
        x_H2O_v_c_ext = Phi_c_ext * Psat(T_des) / Pext

        # Molar fraction of H2 in the dry gas mixture (H2/N2)
        y_H2_agc = C_H2_agc / (C_H2_agc + C_N2_a)
        y_O2_cgc = C_O2_cgc / (C_O2_cgc + C_N2_c)

        # Dynamic viscosity of the gas mixture at the anode side.
        if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
            mu_gaz_asm = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_asm, 1 - x_H2O_v_asm], T_des)
            mu_gaz_agc = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_agc, 1 - x_H2O_v_agc], T_agc)
            mu_gaz_aem = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_aem, 1 - x_H2O_v_aem], T_des)
            mu_gaz_a_ext = mu_mixture_gases(['H2O_v', 'H2'], [x_H2O_v_a_ext, 1 - x_H2O_v_a_ext], T_des)
        else:  # type_auxiliary == "forced-convective_cathode_with_flow-through_anode"
            mu_gaz_asm = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                          [x_H2O_v_asm, y_H2_in * (1 - x_H2O_v_asm), (1 - y_H2_in) * (1 - x_H2O_v_asm)],
                                          T_des)
            mu_gaz_agc = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                          [x_H2O_v_agc, y_H2_agc * (1 - x_H2O_v_agc),
                                           (1 - y_H2_agc) * (1 - x_H2O_v_agc)], T_agc)
            mu_gaz_aem = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                          [x_H2O_v_aem, y_H2_aem * (1 - x_H2O_v_aem),
                                           (1 - y_H2_aem) * (1 - x_H2O_v_aem)], T_des)
            mu_gaz_a_ext = mu_mixture_gases(['H2O_v', 'H2', 'N2'],
                                            [x_H2O_v_a_ext, y_H2_aem_out * (1 - x_H2O_v_a_ext), (1 - y_H2_aem_out) * (1 - x_H2O_v_a_ext)],
                                            T_des)
        # Dynamic viscosity of the gas mixture at the cathode side.
        mu_gaz_csm = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                      [x_H2O_v_csm, y_O2_ext * (1 - x_H2O_v_csm), (1 - y_O2_ext) * (1 - x_H2O_v_csm)],
                                      T_des)
        mu_gaz_cgc = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                      [x_H2O_v_cgc, y_O2_cgc * (1 - x_H2O_v_cgc), (1 - y_O2_cgc) * (1 - x_H2O_v_cgc)],
                                      T_cgc)
        mu_gaz_cem = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                      [x_H2O_v_cem, y_O2_cem * (1 - x_H2O_v_cem), (1 - y_O2_cem) * (1 - x_H2O_v_cem)],
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
        # Set to None the variables not used when "no_auxiliary" system is considered
        v_re, Lman_to_endplate, Lman_to_man_gc, k_purge = [None] * 4

    return {'rho_Cp0': rho_Cp0, 'v_re': v_re, 'k_purge': k_purge, 'rho': rho, 'C_tot': C_tot, 'mu_gaz': mu_gaz, 'P': P}


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
