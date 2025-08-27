# -*- coding: utf-8 -*-

"""This module is used to determine intermediate values for the calculation of the differential equations
and to implement integration events.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import Text, Pext, Phi_ext, n_cell, M_H2, M_O2, M_N2, M_H2O, y_O2_ext, R, F
from modules.transitory_functions import average, Psat, C_v_sat, k_H2, k_O2, calculate_rho_Cp0


# ____________________________________________Differential equations modules____________________________________________

def dif_eq_int_values(sv, operating_inputs, control_variables, parameters):
    """This functions calculates intermediate values for the calculation of the differential equations

        Parameters
        ----------
        sv : dict
            Variables calculated by the solver. They correspond to the fuel cell internal states.
            sv is a contraction of solver_variables for enhanced readability.
        operating_inputs : dict
            Operating inputs of the fuel cell.
        control_variables : dict
            Variables controlled by the user.
        parameters : dict
            Parameters of the fuel cell model.

        Returns
        -------
        Mext : float
            Molar mass of the ambient air outside the stack (kg/mol).
        M_H2_N2_in : float
            Molar mass of the inlet gas at the anode side (H2/N2 mixture) (kg/mol).
        Pagc : float
            Global pressure in the anode gas channel (Pa).
        Pcgc : float
            Global pressure in the cathode gas channel (Pa).
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
    C_v_agc, C_v_ampl, C_v_acl = sv['C_v_agc'], sv['C_v_ampl'], sv['C_v_acl']
    C_v_ccl, C_v_cmpl, C_v_cgc = sv['C_v_ccl'], sv['C_v_cmpl'], sv['C_v_cgc']
    s_ampl, s_acl, s_ccl, s_cmpl = sv['s_ampl'], sv['s_acl'], sv['s_ccl'], sv['s_cmpl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_agc, C_H2_ampl, C_H2_acl = sv['C_H2_agc'], sv['C_H2_ampl'], sv['C_H2_acl']
    C_O2_ccl, C_O2_cmpl, C_O2_cgc = sv['C_O2_ccl'], sv['C_O2_cmpl'], sv['C_O2_cgc']
    C_N2_a, C_N2_c = sv['C_N2_a'], sv['C_N2_c']
    T_agc, T_ampl, T_acl, T_mem = sv['T_agc'], sv['T_ampl'], sv['T_acl'], sv['T_mem']
    T_ccl, T_cmpl, T_cgc = sv['T_ccl'], sv['T_cmpl'], sv['T_cgc']
    Pasm, Paem, Pcsm, Pcem = sv['Pasm'], sv['Paem'], sv['Pcsm'], sv['Pcem']
    # Extraction of the operating inputs and the parameters
    T_des, y_H2_in = operating_inputs['T_des'], operating_inputs['y_H2_in']
    Phi_a_des, Phi_c_des = control_variables['Phi_a_des'], control_variables['Phi_c_des']
    Hmem, Hacl, Hccl = parameters['Hmem'], parameters['Hacl'], parameters['Hccl']
    epsilon_gdl, epsilon_cl, epsilon_mpl = parameters['epsilon_gdl'], parameters['epsilon_cl'], parameters['epsilon_mpl']
    kappa_co, epsilon_mc = parameters['kappa_co'], parameters['epsilon_mc']
    n_gdl = parameters['n_gdl']

    # Physical quantities outside the stack
    # Molar masses
    Mext = Phi_ext * Psat(Text) / Pext * M_H2O + \
           y_O2_ext * (1 - Phi_ext * Psat(Text) / Pext) * M_O2 + \
           (1 - y_O2_ext) * (1 - Phi_ext * Psat(Text) / Pext) * M_N2
    M_H2_N2_in = y_H2_in * M_H2 + (1 - y_H2_in) * M_N2

    # Physical quantities inside the stack
    #       Pressures
    Pagc = (C_v_agc + C_H2_agc + C_N2_a) * R * T_agc
    Pcgc = (C_v_cgc + C_O2_cgc + C_N2_c) * R * T_cgc
    #       Humidities
    Phi_agc = C_v_agc / C_v_sat(T_agc)
    Phi_cgc = C_v_cgc / C_v_sat(T_cgc)
    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_H2_agc = C_H2_agc / (C_H2_agc + C_N2_a)
    y_O2_cgc = C_O2_cgc / (C_O2_cgc + C_N2_c)
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
        'ampl': calculate_rho_Cp0('ampl', T_ampl, C_v=C_v_ampl, s=s_ampl, C_H2=C_H2_ampl, C_N2=C_N2_a, epsilon=epsilon_mpl),
        'acl': calculate_rho_Cp0('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl, C_N2=C_N2_a, C_H2=C_H2_acl,
                                 epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
        'mem': calculate_rho_Cp0('mem', T_mem, lambdaa=lambda_mem),
        'ccl': calculate_rho_Cp0('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl, C_O2=C_O2_ccl, C_N2=C_N2_c,
                                 epsilon=epsilon_cl, epsilon_mc=epsilon_mc),
        'cmpl': calculate_rho_Cp0('cmpl', T_cmpl, C_v=C_v_cmpl, s=s_cmpl, C_O2=C_O2_cmpl, C_N2=C_N2_c,
                                 epsilon=epsilon_mpl),
        **{f'cgdl_{i}': calculate_rho_Cp0('cgdl', sv[f'T_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                                          s=sv[f's_cgdl_{i}'], C_O2=sv[f'C_O2_cgdl_{i}'], C_N2=C_N2_c, epsilon=epsilon_gdl)
           for i in range(1, n_gdl + 1)}
        }

    # Physical quantities inside the auxiliary system
    if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation" or \
       parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode":
        # Pressure
        Pp = Pasm
        # Humidities
        Phi_aem = Phi_agc * Paem / Pagc
        Phi_asm = Phi_aem * Pp / Paem
        Phi_cem = Phi_cgc * Pcem / Pcgc
        # Molar masses
        if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation":
            Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
            Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
        else:  # parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode"
            Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O + \
                   y_H2_in * (1 - Phi_a_des * Psat(T_des) / Pasm) * M_H2 + \
                   (1 - y_H2_in) * (1 - Phi_a_des * Psat(T_des) / Pasm) * M_N2
            Maem = Phi_aem * Psat(T_des) / Paem * M_H2O + \
                   y_H2_agc * (1 - Phi_aem * Psat(T_des) / Paem) * M_O2 + \
                   (1 - y_H2_agc) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
        Mcsm = Phi_c_des * Psat(T_des) / Pcsm * M_H2O + \
               y_O2_ext * (1 - Phi_c_des * Psat(T_des) / Pcsm) * M_O2 + \
               (1 - y_O2_ext) * (1 - Phi_c_des * Psat(T_des) / Pcsm) * M_N2
        Mcem = Phi_cem * Psat(T_des) / Pcem * M_H2O + \
               y_O2_cgc * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 + \
               (1 - y_O2_cgc) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2
    else:  # parameters["type_auxiliary"] == "no_auxiliary"
        Masm, Maem, Mcsm, Mcem = [0] * 4

    return Mext, M_H2_N2_in, Pagc, Pcgc, i_n, Masm, Maem, Mcsm, Mcem, rho_Cp0


def desired_flows(solver_variables, control_variables, i_n, i_fc, operating_inputs, parameters, Mext, M_H2_N2_in):
    """
    This function calculates the desired flow for the air compressor and the humidifiers. These desired flow are
    different from the real ones as the corresponding machines takes time to reach the desired values.

    Parameters
    ----------
    solver_variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    control_variables : dict
        Variables controlled by the user.
    i_n : float
        Internal current density (A/m²).
    i_fc : float
        Fuel cell current density (A/m²).
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    Mext : float
        Molar mass of the ambient air outside the stack (kg/mol).
    M_H2_N2_in : float
        Molar mass of the inlet gas at the anode side (H2/N2 mixture) (kg/mol).

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
    Pasm, Pcsm, Wcp = solver_variables['Pasm'], solver_variables['Pcsm'], solver_variables['Wcp']
    # Extraction of the operating inputs and the parameters
    T_des, Sa, Sc = operating_inputs['T_des'], operating_inputs['Sa'], operating_inputs['Sc']
    y_H2_in = operating_inputs['y_H2_in']
    Phi_a_des, Phi_c_des = control_variables['Phi_a_des'], control_variables['Phi_c_des']
    Aact, type_auxiliary = parameters['Aact'], parameters['type_auxiliary']

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # The desired air compressor flow rate Wcp_des (kg.s-1)
        Wcp_des = n_cell * Mext * Pext / (Pext - Phi_ext * Psat(Text)) * \
                  1 / y_O2_ext * Sc * (i_fc + i_n) / (4 * F) * Aact

        # The desired humidifier flow rate at the anode side Wa_v_inj_des (kg.s-1)
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            Prd = Pasm
            Wrd = n_cell * M_H2_N2_in / y_H2_in * Sa * (i_fc + i_n) / (2 * F) * Aact
            Wa_inj_des = (M_H2O * Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des)) /
                          (1 - Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des))) * (Wrd / M_H2_N2_in))
        else:  # type_auxiliary == "forced-convective_cathode_with_anodic_recirculation"
            Wa_inj_des = 0

        # The desired humidifier flow rate at the cathode side Wc_inj_des (kg.s-1)
        Pcp = Pcsm
        Wv_hum_in = M_H2O * Phi_ext * Psat(Text) / Pext * (Wcp / Mext)  # Vapor flow rate from the outside
        Wc_v_des = M_H2O * Phi_c_des * Psat(T_des) / Pcp * (Wcp / Mext)  # Desired vapor flow rate
        Wc_inj_des = Wc_v_des - Wv_hum_in  # Desired humidifier flow rate

    else:  # elif type_auxiliary == "no_auxiliary":
        Wcp_des, Wa_inj_des, Wc_inj_des = 0, 0, 0

    return Wcp_des, Wa_inj_des, Wc_inj_des


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
