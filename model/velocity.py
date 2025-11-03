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

    T_des, Pa_des, Pc_des = operating_inputs['T_des'], operating_inputs['Pa_des'], operating_inputs['Pc_des']
    Hagc, Hcgc, Wagc, Wcgc = parameters['Hagc'], parameters['Hcgc'], parameters['Wagc'], parameters['Wcgc']
    Lgc, Ldist = parameters['Lgc'], parameters['Ldist']
    nb_channel_in_gc, nb_cell, type_auxiliary = parameters['nb_channel_in_gc'], parameters['nb_cell'], parameters['type_auxiliary']
    nb_gc, nb_gdl = parameters['nb_gc'], parameters['nb_gdl']

    # Intermediate calculation
    L_node_gc = Lgc / nb_gc                                                                                             # Length of one gas channel node (m).
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        Pa_ext = Pext
        Pc_ext = Pext
    else:  # elif type_auxiliary == "no_auxiliary":
        Pa_ext = Pa_des
        Pc_ext = Pc_des
    C_N2_a_mean = (sum(sv[f'C_N2_agc_{i}'] for i in range(1, nb_gc + 1)) / nb_gc)
    C_N2_c_mean = (sum(sv[f'C_N2_cgc_{i}'] for i in range(1, nb_gc + 1)) / nb_gc)

    # Calculation of the boundary conditions at the GC/GDL interface
    C_tot_agdl = sv['C_v_agdl_1'] + sv['C_H2_agdl_1'] + C_N2_a_mean
    C_tot_cgdl = sv[f'C_v_cgdl_{nb_gdl}'] + sv[f'C_O2_cgdl_{nb_gdl}'] + C_N2_c_mean
    J_tot_agc_agdl = [None] + [0.0] * nb_gc
    J_tot_cgdl_cgc = [None] + [0.0] * nb_gc
    for i in range(1, nb_gc + 1):
        J_tot_agc_agdl[i] = Jv_agc_agdl[i] + J_H2_agc_agdl[i]                                                           # Total molar flow from the AGC to the AGDL (mol.m-2.s-1).
        J_tot_cgdl_cgc[i] = Jv_cgdl_cgc[i] + J_O2_cgdl_cgc[i]                                                           # Total molar flow from the CGDL to the CGC (mol.m-2.s-1).

    # Residual function for least_squares solver applied on the inlet molar flows
    def residuals(J_in_guessed):
        # Intermediate values
        J_a_in_guesses, J_c_in_guessed = float(J_in_guessed[0]), float(J_in_guessed[1])                                 # Inlet molar flows at anode and cathode (mol/s).
        J_a, J_c = [J_a_in_guesses] + [0] * (nb_gc + 1), [J_c_in_guessed] + [0] * (nb_gc + 1)                           # J[0] is the cell inlet molar flow, J[i] is the outlet molar flow at node i (upwind numerical scheme), and J[-1] is the cell outlet molar flow.
        P_a, P_c = [0] * (nb_gc + 1) + [Pa_ext], [0] * (nb_gc + 1) + [Pc_ext]                                           # P[0] is the inlet pressure, P[i] is the pressure at the outlet of node i, and P[-1] is the outlet pressure.
        v_a, v_c = [0] * (nb_gc + 2), [0] * (nb_gc + 2)                                                                 # v[0] is the inlet velocity, v[i] is the velocity at the outlet of node i, and v[-1] is the outlet velocity.

        # Continuity equation are used for calculating the molar flows at all points along the GC at stationary state.
        for i in range(1, nb_gc + 1):
            J_a[i] = J_a[i-1] - J_tot_agc_agdl[i] * L_node_gc / Hagc
            J_c[i] = J_c[i-1] + J_tot_cgdl_cgc[i] * L_node_gc / Hcgc
        J_a[-1] = J_a[nb_gc]                                                                                            # Inside the distributor, there are no mass transfer with the GDL.
        J_c[-1] = J_c[nb_gc]                                                                                            # Inside the distributor, there are no mass transfer with the GDL.

        # Velocities at the outlets of the GCs (m/s).
        v_a[-1] = J_a[-1] / P_a[-1] * R * T_des
        v_c[-1] = J_c[-1] / P_c[-1] * R * T_des

        # Backward calculation of pressures and velocities along the GC using modified Hagen-Poiseuille equation
        P_a[nb_gc] = P_a[-1] + 8 * math.pi * mu_gaz[f'agc_{nb_gc}'] * Ldist / (Hagc * Wagc) * v_a[-1]                   # Hagen-Poiseuille equation inside the distributor.
        v_a[nb_gc] = J_a[nb_gc] / P_a[nb_gc] * R * T_des                                                                # It is considered that mu_gaz[-1] at static state is closed to mu_gaz[f'cgc_{nb_gc}'] in order to simplify the calculations.
        P_c[nb_gc] = P_c[-1] + 8 * math.pi * mu_gaz[f'cgc_{nb_gc}'] * Ldist / (Hcgc * Wcgc) * v_c[-1]
        v_c[nb_gc] = J_c[nb_gc] / P_c[nb_gc] * R * T_des
        for i in range(nb_gc, 0, -1):
            # At the node side
            P_a[i-1] = P_a[i] + 8 * math.pi * mu_gaz[f'agc_{i}'] * L_node_gc / (Hagc * Wagc) * \
                        (v_a[i] + J_tot_agc_agdl[i] / C_tot_agdl)                                                       # Modified Hagen-Poiseuille equation, considering the convective vapor flow from the GC to the GDL.
            v_a[i-1] = J_a[i-1] / P_a[i-1] * R * T_des                                                                  # Velocity calculated using the known molar flow and pressure at node i-1.
            # At the cathode side
            P_c[i-1] = P_c[i] + 8 * math.pi * mu_gaz[f'cgc_{i}'] * L_node_gc / (Hcgc * Wcgc) * \
                        (v_c[i] - J_tot_cgdl_cgc[i] / C_tot_cgdl)                                                       # It is considered that mu_gaz at static state is closed to mu_gaz[f'cgc_{i}'] in order to simplify the calculations.
            v_c[i-1] = J_c[i-1] / P_c[i-1] * R * T_des

        # Desired molar flows at anode and cathode
        if type_auxiliary == "no_auxiliary":
            W_des_calculated = desired_flows(sv, control_variables, i_fc, P_a[0], P_c[0], operating_inputs, parameters)
            Wa_in_calculated = W_des_calculated['H2'] + W_des_calculated['H2O_inj_a']  # This expression is also present in auxiliaries.py.
            J_a_in_calculated = Wa_in_calculated / (Hagc * Wagc) / nb_cell / nb_channel_in_gc
            Wc_in_calculated = W_des_calculated['dry_air'] + W_des_calculated['H2O_inj_c']  # This expression is also present in calculate_velocity_evolution.
            J_c_in_calculated = Wc_in_calculated / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc

        # Residuals: difference between
        res_a = J_a_in_calculated - J_a_in_guesses
        res_c = J_c_in_calculated - J_c_in_guessed
        return [res_a, res_c]

    # Calculation of the initial molar flow using least squares and pressure drop relation
    #       Initial guesses, bounds
    v_medium = 10 # Initial guess for the velocity (m/s).
    v_max = 100  # Maximum velocity (m/s)
    P_max = 4e5  # Maximum pressure (Pa)
    x0 = [v_medium * Pa_des / (R * T_des), v_medium * Pc_des / (R * T_des)]
    lb = [1e-8, 1e-8]
    ub = [v_max * P_max / (R * T_des), v_max * P_max / (R * T_des)]
    #       Solver call
    sol = least_squares(residuals, x0, bounds=(lb, ub), method='trf')
    #       Check for convergence
    if not sol.success:
        raise RuntimeError(f"Convergence failed in calculate_velocity_evolution: {sol.message}")
    #      Extract initial flow rates
    J_a_in, J_c_in = float(sol.x[0]), float(sol.x[1])

    # Calculation of the velocity profiles along the gas channels using the found inlet molar flows
    #       Intermediate values
    J_a, J_c = [J_a_in] + [0] * (nb_gc + 1), [J_c_in] + [0] * (nb_gc + 1)                                               # J[0] is the cell inlet molar flow, J[i] is the outlet molar flow at node i (upwind numerical scheme), and J[-1] is the cell outlet molar flow.
    P_a, P_c = [0] * (nb_gc + 1) + [Pa_ext], [0] * (nb_gc + 1) + [Pc_ext]                                               # P[0] is the inlet pressure, P[i] is the pressure at the outlet of node i, and P[-1] is the outlet pressure.
    v_a, v_c = [0] * (nb_gc + 2), [0] * (nb_gc + 2)                                                                     # v[0] is the inlet velocity, v[i] is the velocity at the outlet of node i, and v[-1] is the outlet velocity.

    #       Continuity equation are used for calculating the molar flows at all points along the GC at stationary state.
    for i in range(1, nb_gc + 1):
        J_a[i] = J_a[i - 1] - J_tot_agc_agdl[i] * L_node_gc / Hagc
        J_c[i] = J_c[i - 1] + J_tot_cgdl_cgc[i] * L_node_gc / Hcgc
    J_a[-1] = J_a[nb_gc]                                                                                                # Inside the distributor, there are no mass transfer with the GDL.
    J_c[-1] = J_c[nb_gc]                                                                                                # Inside the distributor, there are no mass transfer with the GDL.

    #       Velocities at the outlets of the GCs (m/s).
    v_a[-1] = J_a[-1] / P_a[-1] * R * T_des
    v_c[-1] = J_c[-1] / P_c[-1] * R * T_des

    #       Backward calculation of pressures and velocities along the GC using modified Hagen-Poiseuille equation
    #       It is considered that mu_gaz at static state is closed to mu_gaz[f'cgc_{i}'] in order to simplify the calculations.
    P_a[nb_gc] = P_a[-1] + 8 * math.pi * mu_gaz[f'agc_{nb_gc}'] * Ldist / (Hagc * Wagc) * v_a[-1]                       # Hagen-Poiseuille equation inside the distributor.
    v_a[nb_gc] = J_a[nb_gc] / P_a[nb_gc] * R * T_des                                                                    # It is considered that mu_gaz[-1] at static state is closed to mu_gaz[f'cgc_{nb_gc}'] in order to simplify the calculations.
    P_c[nb_gc] = P_c[-1] + 8 * math.pi * mu_gaz[f'cgc_{nb_gc}'] * Ldist / (Hcgc * Wcgc) * v_c[-1]
    v_c[nb_gc] = J_c[nb_gc] / P_c[nb_gc] * R * T_des
    for i in range(nb_gc, 0, -1):
        #       At the node side
        P_a[i - 1] = P_a[i] + 8 * math.pi * mu_gaz[f'agc_{i}'] * L_node_gc / (Hagc * Wagc) * \
                     (v_a[i] + J_tot_agc_agdl[i] / C_tot_agdl)
        v_a[i - 1] = J_a[i - 1] / P_a[i - 1] * R * T_des                                                                # Velocity calculated using the known molar flow and pressure at node i-1.
        #       At the cathode side
        P_c[i - 1] = P_c[i] + 8 * math.pi * mu_gaz[f'cgc_{i}'] * L_node_gc / (Hcgc * Wcgc) * \
                     (v_c[i] - J_tot_cgdl_cgc[i] / C_tot_cgdl)
        v_c[i - 1] = J_c[i - 1] / P_c[i - 1] * R * T_des

    return v_a, v_c, P_a[0], P_c[0]


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