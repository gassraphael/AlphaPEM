# -*- coding: utf-8 -*-

"""This module is used to determine intermediate values for the calculation of the differential equations
and to implement integration events.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value and functions
from configuration.settings import Text, Pext, Phi_ext, n_cell, M_H2, M_O2, M_N2, M_H2O, yO2_ext, R, F
from modules.transitory_functions import Psat, C_v_sat, k_H2, k_O2, rho_H2O


# ____________________________________________Differential equations modules____________________________________________

def dif_eq_int_values(solver_variables, control_variables, operating_inputs, parameters):
    """This functions calculates intermediate values for the calculation of the differential equations

        Parameters
        ----------
        solver_variables : dict
            Variables calculated by the solver. They correspond to the fuel cell internal states.
        control_variables : dict
            Variables controlled by the user.
        operating_inputs : dict
            Operating inputs of the fuel cell.
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
        rho_H2O(Tfc) : float
            Density of water vapor at the fuel cell temperature (kg/m³).
        """

    # Extraction of the variables
    C_v_agc, C_v_cgc = solver_variables['C_v_agc'], solver_variables['C_v_cgc']
    C_H2_agc, C_H2_acl = solver_variables['C_H2_agc'], solver_variables['C_H2_acl']
    C_O2_ccl, C_O2_cgc, C_N2 = solver_variables['C_O2_ccl'], solver_variables['C_O2_cgc'], solver_variables['C_N2']
    lambda_mem, Pasm = solver_variables['lambda_mem'], solver_variables['Pasm']
    Paem, Pcsm, Pcem = solver_variables['Paem'], solver_variables['Pcsm'], solver_variables['Pcem']
    # Extraction of the operating inputs and the parameters
    Tfc, Phi_c_des = operating_inputs['Tfc'], control_variables['Phi_c_des']
    Hmem, kappa_co = parameters['Hmem'], parameters['kappa_co']

    # Physical quantities outside the stack
    # Molar masses
    Mext = Phi_ext * Psat(Text) / Pext * M_H2O + \
           yO2_ext * (1 - Phi_ext * Psat(Text) / Pext) * M_O2 + \
           (1 - yO2_ext) * (1 - Phi_ext * Psat(Text) / Pext) * M_N2

    # Physical quantities inside the stack
    #       Pressures
    Pagc = (C_v_agc + C_H2_agc) * R * Tfc
    Pcgc = (C_v_cgc + C_O2_cgc + C_N2) * R * Tfc
    #       Humidities
    Phi_agc = C_v_agc / C_v_sat(Tfc)
    Phi_cgc = C_v_cgc / C_v_sat(Tfc)
    #       Oxygen ratio in dry air
    y_c = C_O2_cgc / (C_O2_cgc + C_N2)
    #       Internal current density
    i_H2 = 2 * F * R * Tfc / Hmem * C_H2_acl * k_H2(lambda_mem, Tfc, kappa_co)
    i_O2 = 4 * F * R * Tfc / Hmem * C_O2_ccl * k_O2(lambda_mem, Tfc, kappa_co)
    i_n = i_H2 + i_O2

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
        Masm = Phi_asm * Psat(Tfc) / Pasm * M_H2O + \
               (1 - Phi_asm * Psat(Tfc) / Pasm) * M_H2
        Maem = Phi_aem * Psat(Tfc) / Paem * M_H2O + \
               (1 - Phi_aem * Psat(Tfc) / Paem) * M_H2
        Mcsm = Phi_c_des * Psat(Tfc) / Pcsm * M_H2O + \
               yO2_ext * (1 - Phi_c_des * Psat(Tfc) / Pcsm) * M_O2 + \
               (1 - yO2_ext) * (1 - Phi_c_des * Psat(Tfc) / Pcsm) * M_N2
        Mcem = Phi_cem * Psat(Tfc) / Pcem * M_H2O + \
               y_c * (1 - Phi_cem * Psat(Tfc) / Pcem) * M_O2 + \
               (1 - y_c) * (1 - Phi_cem * Psat(Tfc) / Pcem) * M_N2
    else:  # parameters["type_auxiliary"] == "no_auxiliary"
        Masm, Maem, Mcsm, Mcem = [0] * 4

    return Mext, Pagc, Pcgc, i_n, Masm, Maem, Mcsm, Mcem, rho_H2O(Tfc)


def desired_flows(solver_variables, control_variables, i_n, i_fc, operating_inputs, parameters, Mext):
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
    Tfc, Sa, Sc = operating_inputs['Tfc'], operating_inputs['Sa'], operating_inputs['Sc']
    Phi_a_des, Phi_c_des = control_variables['Phi_a_des'], control_variables['Phi_c_des']
    Aact, type_auxiliary = parameters['Aact'], parameters['type_auxiliary']

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # Intermediate values
        Prd = Pasm
        Pcp = Pcsm

        # The desired air compressor flow rate Wcp_des (kg.s-1)
        Wcp_des = n_cell * Mext * Pext / (Pext - Phi_ext * Psat(Text)) * \
                  1 / yO2_ext * Sc * (i_fc + i_n) / (4 * F) * Aact

        # The desired humidifier flow rate at the anode side Wa_v_inj_des (kg.s-1)
        Wrd = n_cell * M_H2 * Sa * (i_fc + i_n) / (2 * F) * Aact
        Wa_inj_des = (M_H2O * Phi_a_des * Psat(Tfc) / (Prd + Phi_a_des * Psat(Tfc)) /
                      (1 - Phi_a_des * Psat(Tfc) / (Prd + Phi_a_des * Psat(Tfc))) * (Wrd / M_H2))

        # The desired humidifier flow rate at the cathode side Wc_inj_des (kg.s-1)
        Wv_hum_in = M_H2O * Phi_ext * Psat(Text) / Pext * (Wcp / Mext)  # Vapor flow rate from the outside
        Wc_v_des = M_H2O * Phi_c_des * Psat(Tfc) / Pcp * (Wcp / Mext)  # Desired vapor flow rate
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
