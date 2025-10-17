# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from configuration.settings import F, Text, Pext, Phi_ext, y_O2_ext, tau_cp, tau_hum, R, Kp_acp, Kd_acp, Kp_ccp, Kd_ccp, Kp_T, Kd_T
from modules.transitory_functions import Psat


# ____________________________________________________Main functions____________________________________________________


def calculate_dyn_air_compressor_evolution(dif_eq, Pacp_des, Pasm_out, Pccp_des, Pcsm_out, type_auxiliary, **kwargs):
    """This function calculates the dynamic evolution of the air compressor.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    """
    if type_auxiliary == "no_auxiliary":
        # At the anode side
        dif_eq['dPasm_out / dt'] = (Pacp_des - Pasm_out) / tau_cp
        # At the cathode side
        dif_eq['dPcsm_out / dt'] = (Pccp_des - Pcsm_out) / tau_cp


def calculate_dyn_air_compressor_controler_evolution(dif_eq, Wacp_des, dWacp_desdt, Wa_inj_des, dWa_inj_desdt,
                                                     Wccp_des, dWccp_desdt, Wc_inj_des, dWc_inj_desdt, Pasm_out,
                                                     Pcsm_out, v_asm_out, v_csm_out, T_des, Hagc, Hcgc, Wagc, Wcgc,
                                                     type_auxiliary, **kwargs):

    if type_auxiliary == "no_auxiliary":
        # At the anode side
        Wacp_to_asm_out = Pasm_out * v_asm_out * Hagc * Wagc / (R * T_des)                                                     # Volume flow rate at the anode compressor outlet (mol.s-1)
        dWacp_to_asm_outdt = (v_asm_out * dif_eq['dPasm_out / dt'] + Pasm_out * dif_eq['dv_asm_out / dt']) * Hagc * Wagc / (R * T_des)  # Derivative of the anode volume flow rate (mol.s-2). This is calculated assuming an isothermal flow.

        Wa_des = Wacp_des + Wa_inj_des
        # dWa_desdt = dWacp_desdt + dWa_inj_desdt
        dWa_desdt = 0                                                                                                   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        dif_eq['dPacp_des / dt'] = Kp_acp * (Wa_des - Wacp_to_asm_out) + Kd_acp * (dWa_desdt - dWacp_to_asm_outdt)                      # PD controller for the anode compressor outlet pressure
        # dif_eq['dPacp_des / dt'] = 0                                                                                    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # At the cathode side
        Wccp_to_csm_out = Pcsm_out * v_csm_out * Hcgc * Wcgc / (R * T_des)                                                     # Volume flow rate at the cathode compressor outlet (mol.s-1)
        dWccp_to_csm_outdt = (v_csm_out * dif_eq['dPcsm_out / dt'] + Pcsm_out * dif_eq['dv_csm_out / dt']) * Hcgc * Wcgc / (R * T_des)  # Derivative of the cathode volume flow rate (mol.s-2). This is calculated assuming an isothermal flow.

        Wc_des = Wccp_des + Wc_inj_des  # Desired volume flow rate at the entry of the cathode side (mol.s-1)
        # dWc_desdt = dWccp_desdt + dWc_inj_desdt
        dWc_desdt = 0                                                                                                   # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        dif_eq['dPccp_des / dt'] = Kp_ccp * (Wc_des - Wccp_to_csm_out) + Kd_ccp * (dWc_desdt - dWccp_to_csm_outdt)                      # PD controller for the cathode compressor outlet pressure
        # dif_eq['dPccp_des / dt'] = 0                                                                                  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


def calculate_dyn_humidifier_evolution(dif_eq, Wa_inj_des, Wc_inj_des, Wa_inj, Wc_inj, type_auxiliary, **kwargs):
    # Anode and cathode humidifiers evolution
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
        dif_eq['dWa_inj / dt'] = 0
        dif_eq['dWc_inj / dt'] = (Wc_inj_des - Wc_inj) / tau_hum  # Estimation at the first order.
    elif type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dWa_inj / dt'] = (Wa_inj_des - Wa_inj) / tau_hum  # Estimation at the first order.
        dif_eq['dWc_inj / dt'] = (Wc_inj_des - Wc_inj) / tau_hum  # Estimation at the first order.


def calculate_dyn_throttle_area_controler(dif_eq, sv, Pa_des, Pc_des, A_T_a, A_T_c, type_auxiliary, Pagc, Pcgc, **kwargs):
    """This function calculates the dynamic evolution of the throttle area inside the anode and cathode auxiliaries.
    This function has to be executed after 'calculate_dyn_vapor_evolution' and 'calculate_dyn_H2_O2_N2_evolution'.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    sv : dict
        Dictionary containing the solver variables.
    Pa_des : float
        Desired pressure inside the anode gas channel (Pa).
    Pc_des : float
        Desired pressure inside the cathode gas channel (Pa).
    type_auxiliary : str
        Type of auxiliary components used in the fuel cell model.
    Pagc : float
        Pressure inside the anode gas channel (Pa).
    Pcgc : float
        Pressure inside the cathode gas channel (Pa).
    """

    # Extraction of the variables
    T_agc, T_cgc, Abp_a, Abp_c = sv['T_agc'], sv['T_cgc'], sv.get('Abp_a', None), sv.get('Abp_c', None)

    # Calculation of the pressure derivative inside the gas channels
    dPagcdt = (dif_eq['dC_v_agc / dt'] + dif_eq['dC_H2_agc / dt'] + dif_eq['dC_N2_a / dt']) * R * T_agc
    dPcgcdt = (dif_eq['dC_v_cgc / dt'] + dif_eq['dC_O2_cgc / dt'] + dif_eq['dC_N2_c / dt']) * R * T_cgc

    # Throttle area evolution inside the anode auxiliaries
    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dAbp_a / dt'] = - Kp_T * (Pa_des - Pagc) + Kd_T * dPagcdt  # PD controller
        if Abp_a > A_T_a and dif_eq['dAbp_a / dt'] > 0:  # The throttle area cannot be higher than the maximum value
            dif_eq['dAbp_a / dt'] = 0
        elif Abp_a < 0 and dif_eq['dAbp_a / dt'] < 0:  # The throttle area cannot be lower than 0
            dif_eq['dAbp_a / dt'] = 0

    # Throttle area evolution inside the cathode auxiliaries
    dif_eq['dAbp_c / dt'] = - Kp_T * (Pc_des - Pcgc) + Kd_T * dPcgcdt  # PD controller
    if Abp_c > A_T_c and dif_eq['dAbp_c / dt'] > 0:  # The throttle area cannot be higher than the maximum value
        dif_eq['dAbp_c / dt'] = 0
    elif Abp_c < 0 and dif_eq['dAbp_c / dt'] < 0:  # The throttle area cannot be lower than 0
        dif_eq['dAbp_c / dt'] = 0
