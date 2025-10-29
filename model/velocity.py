# -*- coding: utf-8 -*-

"""This file represents the calculation of the velocity over time. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________
# Importing the necessary libraries
import math
from scipy.optimize import least_squares

# Importing constants' value and functions
from configuration.settings import R, F, Text, Pext, Phi_ext, y_O2_ext
from modules.transitory_functions import Psat


# ________________________________________________________Velocity______________________________________________________


def calculate_inlet_pressure_evolution(sv, v_a, v_c, Jv_agc_agdl, Jv_cgdl_cgc, J_H2_agc_agdl, J_O2_cgdl_cgc, rho,
                                       mu_gaz, Pa_des, Pc_des, Hagc, Hcgc, Wagc, Wcgc, Lgc, nb_gc,  nb_gdl,
                                       type_auxiliary, **kwarks):
    """
        Calculate inlet pressures at anode and cathode using the Hagen-Poiseuille relation.

        Parameters
        ----------
        sv : dict
            Solver state variables.
        v_a : float
            Anode gas velocity (m/s).
        v_c : float
            Cathode gas velocity (m/s).
        Jv_agc_agdl : list
            Vapor flow between the AGC and the AGDL at each GC node (mol.m-2.s-1).
        Jv_cgdl_cgc : list
            Vapor flow between the CGDL and the CGC at each GC node (mol.m-2.s-1).
        J_H2_agc_agdl : list
            H2 flow between the AGC and the AGDL at each GC node (mol.m-2.s-1).
        J_O2_cgdl_cgc : list
            O2 flow between the CGDL and the CGC at each GC node (mol.m-2.s-1).
        rho: dict
            Gas densities at anode and cathode (kg/m³).
        mu_gaz: dict
            Gas dynamic viscosities at anode and cathode (Pa·s).
        Pa_des : float
            Desired anode pressure (Pa).
        Pc_des : float
            Desired cathode pressure (Pa).
        Hagc : float
            Anode gas channel thickness (m).
        Hcgc : float
            Cathode gas channel thickness (m).
        Wagc : float
            Anode gas channel width (m).
        Wcgc : float
            Cathode gas channel width (m).
        Lgc : float
            Length of the gas channel (m).
        nb_gc : int
            Number of gas channel nodes.
        nb_gdl : int
            Number of gas diffusion layer nodes.
        type_auxiliary : str
            Auxiliary configuration identifier.

        Returns
        -------
        Pa_in : float
            Anode inlet pressure (Pa).
        Pc_in : float
            Cathode inlet pressure (Pa).
    """

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        pass
        # L_node_gc = Lgc / nb_gc                                                                                         # Length of one gas channel node (m).
        # Dh_agc = math.sqrt(4 * Hagc * Wagc / math.pi)                                                                   # Hydraulic diameter of the anode gas channel (m).
        # Dh_cgc = math.sqrt(4 * Hcgc * Wcgc / math.pi)                                                                   # Hydraulic diameter of the cathode gas channel (m).
        # Delta_Pagc = [None] + [0.0] * nb_gc
        # Delta_Pcgc = [None] + [0.0] * nb_gc
        #
        # # At the anode side
        # for i in range(1, nb_gc + 1):
        #     Re_a_i = rho[f'agc_{i}'] * v_a * Dh_agc / mu_gaz[f'agc_{i}']                                                # Reynolds number at anode node i.
        #     f_D_a_i = 0.3164 * Re_a_i ** (-1 / 4)                                                                       # Blasius correlation for the Darcy friction factor at anode node i.
        #     Delta_Pagc[i] = f_D_a_i * L_node_gc / Dh_agc * rho[f'agc_{i}'] * v_a ** 2 / 2                               # Darcy-Weisbach equation for the pressure drop at anode node i.
        # Pa_in = Pa_des + sum(Delta_Pagc[1:])
        #
        # # At the cathode side
        # for i in range(1, nb_gc + 1):
        #     Re_c_i = rho[f'cgc_{i}'] * v_c * Dh_cgc / mu_gaz[f'cgc_{i}']                                                # Reynolds number at cathode node i.
        #     f_D_c_i = 0.3164 * Re_c_i ** (-1 / 4)                                                                       # Blasius correlation for the Darcy friction factor at cathode node i.
        #     Delta_Pcgc[i] = f_D_c_i * L_node_gc / Dh_cgc * rho[f'cgc_{i}'] * v_c ** 2 / 2                               # Darcy-Weisbach equation for the pressure drop at cathode node i.
        # Pc_in = Pc_des + sum(Delta_Pcgc[1:])

    if type_auxiliary == "no_auxiliary":
        # Pa_in = Pa_des + 8 * math.pi * mu_a_times_Lgc / (Hagc * Wagc) * v_a                                             # Hagen-Poiseuille equation for the pressure drop.
        # Pc_in = Pc_des + 8 * math.pi * mu_c_times_Lgc / (Hcgc * Wcgc) * v_c                                             # Hagen-Poiseuille equation for the pressure drop.

        L_node_gc = Lgc / nb_gc                                                                                         # Length of one gas channel node (m).
        Delta_Pagc = [None] + [0.0] * nb_gc
        Delta_Pcgc = [None] + [0.0] * nb_gc
        C_tot_agdl = sv['C_v_agdl_1'] + sv['C_H2_agdl_1']
        C_tot_cgdl = sv[f'C_v_cgdl_{nb_gdl}'] + sv[f'C_O2_cgdl_{nb_gdl}'] + sv['C_N2_c']
        J_tot_agc_agdl = [None] + [0.0] * nb_gc
        J_tot_cgdl_cgc = [None] + [0.0] * nb_gc
        for i in range(1, nb_gc + 1):
            J_tot_agc_agdl[i] = Jv_agc_agdl[i] + J_H2_agc_agdl[i]                                                       # Total molar flow from the AGC to the AGDL (mol.m-2.s-1).
            J_tot_cgdl_cgc[i] = Jv_cgdl_cgc[i] + J_O2_cgdl_cgc[i]                                                       # Total molar flow from the CGDL to the CGC (mol.m-2.s-1).


        # At the anode side
        for i in range(1, nb_gc + 1):
            Delta_Pagc[i] = 8 * math.pi * mu_gaz[f'agc_{i}'] * L_node_gc / (Hagc * Wagc) * \
                            (v_a + J_tot_agc_agdl[i] / C_tot_agdl)                                                      # Modified Hagen-Poiseuille equation, considering the convective vapor flow from the GC to the GDL.
        Pa_in = Pa_des + sum(Delta_Pagc[1:])

        # At the cathode side
        for i in range(1, nb_gc + 1):
            Delta_Pcgc[i] = 8 * math.pi * mu_gaz[f'cgc_{i}'] * L_node_gc / (Hcgc * Wcgc) * \
                            (v_c - J_tot_cgdl_cgc[i] / C_tot_cgdl)                                                      # Modified Hagen-Poiseuille equation, considering the convective vapor flow from the GC to the GDL.
        Pc_in = Pc_des + sum(Delta_Pcgc[1:])

    return Pa_in, Pc_in


def calculate_velocity_evolution(sv, control_variables, i_fc, Jv_agc_agdl, Jv_cgdl_cgc, J_H2_agc_agdl, J_O2_cgdl_cgc,
                                 operating_inputs, parameters, rho, mu_gaz):
    """
        Calculate the gas velocities at the anode and cathode.
        This function finds the velocities v_a and v_c that make the pressures computed from the Hagen-Poiseuille
        relation equal to the pressures inferred from the ideal gas law applied to the desired molar flows.
        A bounded nonlinear least squares solver is used.

        Parameters
        ----------
        sv : dict
            Solver state variables (passed to `desired_flows`).
        control_variables : dict
            Control inputs (passed to `desired_flows`).
        i_fc : int
            Fuel cell current density at time t.
        Jv_agc_agdl : list
            Vapor flow between the AGC and the AGDL at each GC node (mol.m-2.s-1).
        Jv_cgdl_cgc : list
            Vapor flow between the CGDL and the CGC at each GC node (mol.m-2.s-1).
        J_H2_agc_agdl : list
            H2 flow between the AGC and the AGDL at each GC node (mol.m-2.s-1).
        J_O2_cgdl_cgc : list
            O2 flow between the CGDL and the CGC at each GC node (mol.m-2.s-1).
        operating_inputs : dict
            Operating conditions.
        parameters : dict
            Model parameters.
        rho : dict
            Gas densities at anode and cathode.
        mu_gaz : dict
            Gas dynamic viscosities at anode and cathode.

        Returns
        -------
        v_a : float
            Anode gas velocity (m/s).
        v_c : float
            Cathode gas velocity (m/s).

        Raises
        ------
        RuntimeError
            If the least squares solver does not converge.
        """

    T_des= operating_inputs['T_des']
    Hagc, Hcgc, Wagc, Wcgc = parameters['Hagc'], parameters['Hcgc'], parameters['Wagc'], parameters['Wcgc']
    nb_channel_in_gc, nb_cell, type_auxiliary = parameters['nb_channel_in_gc'], parameters['nb_cell'], parameters['type_auxiliary']

    # Residual function for least_squares
    def residuals(v):
        v_a, v_c = float(v[0]), float(v[1])
        # Pa_in and Pc_in from Hagen-Poiseuille pressure-drop relation
        Pa_in, Pc_in = calculate_inlet_pressure_evolution(sv, v_a, v_c, Jv_agc_agdl, Jv_cgdl_cgc, J_H2_agc_agdl,
                                                          J_O2_cgdl_cgc, rho, mu_gaz, **operating_inputs, **parameters)
        # Desired flows
        W_des = desired_flows(sv, control_variables, i_fc, Pa_in, Pc_in, operating_inputs, parameters)

        if type_auxiliary == "no_auxiliary":
            # At the anode side
            Wa_in = W_des['H2'] + W_des['H2O_inj_a'] # This expression is also present in the auxiliaries model.
            Ja_in = Wa_in / (Hagc * Wagc) / nb_cell / nb_channel_in_gc # This expression is also present in the auxiliaries model.
            Pa_from_flow = Ja_in / v_a * R * T_des  # Ideal gas law to get the velocity from the inlet molar flow.

            # At the cathode side
            Wc_in = W_des['dry_air'] + W_des['H2O_inj_c'] # This expression is also present in the auxiliaries model.
            Jc_in = Wc_in / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc # This expression is also present in the auxiliaries model.
            Pc_from_flow = Jc_in / v_c * R * T_des  # Ideal gas law to get the velocity from the inlet molar flow.

            # Residuals: difference between pressures from ideal gas law and Hagen-Poiseuille relation
            res_a = Pa_from_flow - Pa_in
            res_c = Pc_from_flow - Pc_in
            return [res_a, res_c]

    # Bounded solve in [1e-8, 100]
    x0 = [1.0, 1.0]
    lb = [1e-8, 1e-8]
    ub = [100.0, 100.0]
    sol = least_squares(residuals, x0, bounds=(lb, ub), xtol=1e-10, ftol=1e-10, max_nfev=2000, method='trf')

    if not sol.success:
        raise RuntimeError(f"Convergence failed in calculate_velocity_evolution: {sol.message}")

    v_a, v_c = float(sol.x[0]), float(sol.x[1])
    return v_a, v_c


def desired_flows(solver_variables, control_variables, i_fc, Pa_in, Pc_in, operating_inputs, parameters):
    """
    This function calculates the desired flow for the air compressor and the humidifiers. These desired flow are
    different from the real ones as the corresponding machines takes time to reach the desired values.

    Parameters
    ----------
    solver_variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    control_variables : dict
        Variables controlled by the user.
    i_fc : float
        Fuel cell current density (A/m²).
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    Wcp_des : float
        Desired air compressor flow rate (kg/s).
    W_H2O_inj_a_des : float
        Desired humidifier flow rate at the anode side (kg/s).
    W_H2O_inj_c_des : float
        Desired humidifier flow rate at the cathode side (kg/s).
    """

    # Extraction of the variables
    Pasm, Pcsm, Wcp = solver_variables.get('Pasm', None), solver_variables.get('Pcsm', None), solver_variables.get('Wcp', None)
    # Extraction of the operating inputs and the parameters
    T_des, Sa, Sc = operating_inputs['T_des'], operating_inputs['Sa'], operating_inputs['Sc']
    y_H2_in = operating_inputs['y_H2_in']
    Phi_a_des, Phi_c_des = control_variables['Phi_a_des'], control_variables['Phi_c_des']
    Aact, nb_cell, type_auxiliary = parameters['Aact'], parameters['nb_cell'], parameters['type_auxiliary']

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # The desired air compressor volume flow rate (mol.s-1)                                                         # Warning: considérer le débit minimal du compresseur !
        W_H2_des = 1 / y_H2_in * Sa * i_fc / (2 * F) * (nb_cell * Aact)
        Wacp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc, W_H2_des)                                         # Adjust the desired compressor flow to ensure a minimum flow is maintained.
        W_air_ext_des = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * i_fc / (4 * F) * (nb_cell * Aact)
        Wccp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc, W_air_ext_des)                                         # Adjust the desired compressor flow to ensure a minimum flow is maintained.

        # The desired humidifier volume flow rate at the anode side Wa_v_inj_des (mol.s-1)                              # Warning: considérer le débit minimal du compresseur !
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            Prd = Pasm
            W_H2_des = 1 / y_H2_in * Sa * i_fc / (2 * F) * (nb_cell * Aact)
            W_H2O_inj_a_des = (Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des)) /
                          (1 - Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des))) * W_H2_des)
        else:  # type_auxiliary == "forced-convective_cathode_with_anodic_recirculation"
            W_H2O_inj_a_des = 0

        # The desired humidifier volume flow rate at the cathode side W_H2O_inj_c_des (mol.s-1)
        Pcp = Pcsm
        Wv_hum_in = Phi_ext * Psat(Text) / Pext * Wcp  # Vapor flow rate from the outside
        W_H20_c_des = Phi_c_des * Psat(T_des) / Pcp * Wcp  # Desired vapor flow rate
        W_H2O_inj_c_des = W_H20_c_des - Wv_hum_in  # Desired humidifier flow rate

    else:  # elif type_auxiliary == "no_auxiliary":
        # At the anode side
        W_H2_des = Sa * i_fc / (2 * F) * (nb_cell * Aact)
        W_H2O_inj_a_des = (Phi_a_des * Psat(T_des) / (Pa_in - Phi_a_des * Psat(T_des))) * W_H2_des

        # At the cathode side
        W_dry_air_des = 1 / y_O2_ext * Sc * i_fc / (4 * F) * (nb_cell * Aact)
        W_H2O_inj_c_des = (Phi_c_des * Psat(T_des) / (Pc_in - Phi_c_des * Psat(T_des))) * W_dry_air_des

    return {'H2': W_H2_des, 'dry_air': W_dry_air_des, 'H2O_inj_a': W_H2O_inj_a_des, 'H2O_inj_c': W_H2O_inj_c_des}


def adjust_compressor_flow_with_minimum(i_fc, Wcp_des):
    """
    Adjusts the desired compressor flow rate to ensure a minimum flow is maintained, based on the current density and
     model parameters.

    Parameters
    ----------
    i_fc : float
        Actual fuel cell current density (A/m²).
    parameters : dict
        Model parameters, must include 'pola_current_parameters' with 'delta_i_pola'.
    Wcp_des : float
        Desired compressor flow rate (mol/s).

    Returns
    -------
    Wcp_des_adjusted : float
        Adjusted compressor flow rate (mol/s) ensuring the minimum flow.
    """

    # Parameters for minimum current density adjustment
    i_cp_min = 0.3e4  # (A/m²) Minimum current density for compressor flow.
    delta_i_load_step = 0.01e4  # (A/m²) Minimum current density step for reaching the minimum compressor flow.

    if i_fc <= i_cp_min + 3 * delta_i_load_step:
        Wcp_des_adjusted = (
            Wcp_des * i_cp_min / i_fc * (1.0 + math.tanh(4 * (i_fc - (delta_i_load_step / 2)) / (delta_i_load_step / 2))) / 2
            + Wcp_des * (1 - i_cp_min / i_fc) * (1.0 + math.tanh(4 * (i_fc - i_cp_min - (delta_i_load_step / 2)) / (delta_i_load_step / 2))) / 2
        )
        return Wcp_des_adjusted
    else: # For higher current densities, the compressor flow is not adjusted, and so it is faster to return the original value.
        return Wcp_des