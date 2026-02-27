# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from alphapem.utils.physics_constants import tau_hum, R, Kp_T, Kd_T


# ____________________________________________________Main functions____________________________________________________


def calculate_dyn_air_compressor_evolution(dif_eq, Pacp_des, Pasm_out, Pccp_des, Pcsm_out, type_auxiliary, **kwargs):
    """This function calculates the dynamic evolution of the air compressor.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    """
    pass
    # dif_eq['dPasm_out / dt'] = (Pacp_des - Pasm_out) / tau_cp # Estimation at the first order.


def calculate_dyn_humidifier_evolution(dif_eq, Wa_inj, Wc_inj, type_auxiliary, Wa_inj_des, Wc_inj_des, **kwargs):
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
