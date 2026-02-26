# -*- coding: utf-8 -*-

"""This file represents the calculation of the velocity over time. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________
# Importing the necessary libraries
import math
from scipy.optimize import least_squares

# Importing constants' value and functions
from alphapem.utils.physics_constants import R, F, Text, Pext, Phi_ext, y_O2_ext, M_H2O
from alphapem.utils.maths_functions import average
from alphapem.domain.modules.flows_1D_MEA_modules import h_a, h_c, k_H2, k_O2
from alphapem.utils.physics_functions import Psat, mu_mixture_gases


# ________________________________________________________Velocity______________________________________________________

def calculate_velocity_evolution(sv, i_fc_cell, operating_inputs, parameters):
    """
        Calculate the gas velocities at the anode and cathode.
        This function finds the velocities v_a and v_c that make the pressures computed from the Hagen-Poiseuille
        relation equal to the pressures inferred from the ideal gas law applied to the desired molar flows.
        A bounded nonlinear least squares solver is used.

        Parameters
        ----------
        sv : dict
            Solver state variables (passed to `desired_flows`).
        i_fc_cell : int
            Fuel cell current density at time t.
        operating_inputs : dict
            Operating conditions.
        parameters : dict
            Model parameters.

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

    # Extraction of the operating inputs and the parameters
    T_des, Pa_des, Pc_des = operating_inputs['T_des'], operating_inputs['Pa_des'], operating_inputs['Pc_des']
    Hagc, Hcgc, Wagc, Wcgc = parameters['Hagc'], parameters['Hcgc'], parameters['Wagc'], parameters['Wcgc']
    Lgc, Ldist = parameters['Lgc'], parameters['Ldist']
    Hmem, Hacl, Hccl, kappa_co = parameters['Hmem'], parameters['Hacl'], parameters['Hccl'], parameters['kappa_co']
    nb_channel_in_gc, nb_cell, type_auxiliary = parameters['nb_channel_in_gc'], parameters['nb_cell'], parameters['type_auxiliary']
    nb_gc, nb_gdl = parameters['nb_gc'], parameters['nb_gdl']
    # Extraction of the variables
    C_v_agc = [None] + [sv[i]['C_v_agc'] for i in range(1, nb_gc + 1)]
    C_v_agdl_1 = [None] + [sv[i]['C_v_agdl_1'] for i in range(1, nb_gc + 1)]
    C_v_cgdl_nb_gdl = [None] + [sv[i][f'C_v_cgdl_{nb_gdl}'] for i in range(1, nb_gc + 1)]
    C_v_cgc = [None] + [sv[i]['C_v_cgc'] for i in range(1, nb_gc + 1)]
    C_H2_agc = [None] + [sv[i]['C_H2_agc'] for i in range(1, nb_gc + 1)]
    C_H2_agdl_1 = [None] + [sv[i]['C_H2_agdl_1'] for i in range(1, nb_gc + 1)]
    C_O2_cgdl_nb_gdl = [None] + [sv[i][f'C_O2_cgdl_{nb_gdl}'] for i in range(1, nb_gc + 1)]
    C_O2_cgc = [None] + [sv[i]['C_O2_cgc'] for i in range(1, nb_gc + 1)]
    C_N2_agc = [None] + [sv[i]['C_N2_agc'] for i in range(1, nb_gc + 1)]
    C_N2_cgc = [None] + [sv[i]['C_N2_cgc'] for i in range(1, nb_gc + 1)]
    T_agc = [None] + [sv[i]['T_agc'] for i in range(1, nb_gc + 1)]
    T_cgc = [None] + [sv[i]['T_cgc'] for i in range(1, nb_gc + 1)]

    # Intermediate calculation
    #       Length of one gas channel node (m).
    L_node_gc = Lgc / nb_gc
    #       Pressures (Pa).
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
            type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        Pa_ext = Pext
        Pc_ext = Pext
    else:  # elif type_auxiliary == "no_auxiliary":
        Pa_ext = Pa_des
        Pc_ext = Pc_des
    Pagc = [None] + [(C_v_agc[i] + C_H2_agc[i] + C_N2_agc[i]) * R * T_agc[i]
                     for i in range(1, nb_gc + 1)]
    Pcgc = [None] + [(C_v_cgc[i] + C_O2_cgc[i] + C_N2_cgc[i]) * R * T_cgc[i]
                     for i in range(1, nb_gc + 1)]
    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
    y_O2 = {}
    y_H2 = {}
    for i in range(1, nb_gc + 1):
        y_H2[f'agc_{i}'] = C_H2_agc[i] / (C_H2_agc[i] + C_N2_agc[i])
        y_O2[f'cgc_{i}'] = C_O2_cgc[i] / (C_O2_cgc[i] + C_N2_cgc[i])
    #       Vapor ratio over the gas mixture.
    x_H2O_v = {}
    for i in range(1, nb_gc + 1):
        x_H2O_v[f'agc_{i}'] = C_v_agc[i] / (C_v_agc[i] + C_H2_agc[i] + C_N2_agc[i])
    for i in range(1, nb_gc + 1):
        x_H2O_v[f'cgc_{i}'] = C_v_cgc[i] / (C_v_cgc[i] + C_O2_cgc[i] + C_N2_cgc[i])
    #       Dynamic viscosity of the gas mixture.
    mu_gaz = {}
    for i in range(1, nb_gc + 1):
        mu_gaz[f'agc_{i}'] = mu_mixture_gases(['H2O_v', 'H2'],
                                              [x_H2O_v[f'agc_{i}'], 1 - x_H2O_v[f'agc_{i}']], T_agc[i])
    for i in range(1, nb_gc + 1):
        mu_gaz[f'cgc_{i}'] = mu_mixture_gases(['H2O_v', 'O2', 'N2'],
                                              [x_H2O_v[f'cgc_{i}'], y_O2[f'cgc_{i}'] * (1 - x_H2O_v[f'cgc_{i}']),
                                               (1 - y_O2[f'cgc_{i}']) * (1 - x_H2O_v[f'cgc_{i}'])], T_cgc[i])

    # Calculation of the boundary conditions at the GC/GDL interface
    C_tot_agdl = C_v_agdl_1[i] + C_H2_agdl_1[i] + C_N2_agc[i]
    C_tot_cgdl = C_v_cgdl_nb_gdl[i] + C_O2_cgdl_nb_gdl[i] + C_N2_cgc[i]
    J_tot_agc_agdl = [None] + [0.0] * nb_gc
    J_tot_cgdl_cgc = [None] + [0.0] * nb_gc
    for i in range(1, nb_gc + 1):
        Jv_agc_agdl = h_a(Pagc[i], T_des, Wagc, Hagc) * (C_v_agc[i] - C_v_agdl_1[i])                                    # This equation is also calcultaed in flows.py.
        J_H2_agc_agdl = h_a(Pagc[i], T_des, Wagc, Hagc) * (C_H2_agc[i] - C_H2_agdl_1[i])                                # This equation is also calcultaed in flows.py.
        Jv_cgdl_cgc = h_c(Pcgc[i], T_des, Wcgc, Hcgc) * (C_v_cgdl_nb_gdl[i] - C_v_cgc[i])                                # This equation is also calcultaed in flows.py.
        J_O2_cgdl_cgc = h_c(Pcgc[i], T_des, Wcgc, Hcgc) * (C_O2_cgdl_nb_gdl[i] - C_O2_cgc[i])                            # This equation is also calcultaed in flows.py.
        Jl_agc_agdl = 0.0 / M_H2O                                                                                       # Should be added later, knowing that it requires the knowledge of v_agc...
        Jl_cgdl_cgc = 0.0 / M_H2O                                                                                       # Should be added later, knowing that it requires the knowledge of v_cgc...
        J_tot_agc_agdl[i] = Jv_agc_agdl + J_H2_agc_agdl                                                                 # Total molar flow from the AGC to the AGDL (mol.m-2.s-1).
        J_tot_cgdl_cgc[i] = Jv_cgdl_cgc + J_O2_cgdl_cgc                                                                 # Total molar flow from the CGDL to the CGC (mol.m-2.s-1).

    # Residual function for least_squares solver applied on the inlet molar flows
    def residuals(J_in_guessed):
        # Intermediate values
        J_a_in_guessed, J_c_in_guessed = float(J_in_guessed[0]), float(J_in_guessed[1])                                 # Inlet molar flows at anode and cathode (mol/s).
        J_a, J_c = [J_a_in_guessed] + [0] * (nb_gc + 1), [J_c_in_guessed] + [0] * (nb_gc + 1)                           # J[0] is the cell inlet molar flow, J[i] is the outlet molar flow at node i (upwind numerical scheme), and J[-1] is the cell outlet molar flow.
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
            W_des_calculated = desired_flows(sv, i_fc_cell, P_a[0], P_c[0], operating_inputs, parameters)
            Wa_in_calculated = W_des_calculated['H2'] + W_des_calculated['H2O_inj_a']                                   # This expression is also present in auxiliaries.py.
            J_a_in_calculated = Wa_in_calculated / (Hagc * Wagc) / nb_cell / nb_channel_in_gc
            Wc_in_calculated = W_des_calculated['dry_air'] + W_des_calculated['H2O_inj_c']                              # This expression is also present in calculate_velocity_evolution.
            J_c_in_calculated = Wc_in_calculated / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc

        # Residuals: difference between
        res_a = J_a_in_calculated - J_a_in_guessed
        res_c = J_c_in_calculated - J_c_in_guessed
        return [res_a, res_c]

    # Calculation of the initial molar flow using least squares and pressure drop relation
    #       Initial guesses, bounds
    v_medium = 10 # Initial guess for the velocity (m/s).
    x0 = [v_medium * Pa_des / (R * T_des), v_medium * Pc_des / (R * T_des)]
    #       Solver call
    sol = least_squares(residuals, x0, method='lm')
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


def desired_flows(solver_variables, i_fc_cell, Pa_in, Pc_in, operating_inputs, parameters):
    """
    This function calculates the desired flow for the air compressor and the humidifiers. These desired flow are
    different from the real ones as the corresponding machines takes time to reach the desired values.

    Parameters
    ----------
    solver_variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    i_fc_cell : float
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

    # Extraction of the operating inputs and the parameters
    T_des, Sa, Sc = operating_inputs['T_des'], operating_inputs['Sa'], operating_inputs['Sc']
    Phi_a_des, Phi_c_des = operating_inputs['Phi_a_des'], operating_inputs['Phi_c_des']
    y_H2_in = operating_inputs['y_H2_in']
    kappa_co = parameters['kappa_co']
    Hacl, Hmem, Hccl, Aact = parameters['Hacl'], parameters['Hmem'], parameters['Hccl'], parameters['Aact']
    nb_gc, nb_cell, type_auxiliary = parameters['nb_gc'], parameters['nb_cell'], parameters['type_auxiliary']
    # Extraction of the variables
    # Pasm, Pcsm, Wcp = solver_variables.get('Pasm', None), solver_variables.get('Pcsm', None), solver_variables.get('Wcp', None)
    C_H2_acl = [None] + [solver_variables[i]['C_H2_acl'] for i in range(1, nb_gc + 1)]
    C_O2_ccl = [None] + [solver_variables[i]['C_O2_ccl'] for i in range(1, nb_gc + 1)]
    lambda_mem = [None] + [solver_variables[i]['lambda_mem'] for i in range(1, nb_gc + 1)]
    T_acl = [None] + [solver_variables[i]['T_acl'] for i in range(1, nb_gc + 1)]
    T_mem = [None] + [solver_variables[i]['T_mem'] for i in range(1, nb_gc + 1)]
    T_ccl = [None] + [solver_variables[i]['T_ccl'] for i in range(1, nb_gc + 1)]

    # Physical quantities inside the stack
    #       The crossover current density i_n
    i_n = []
    for i in range(1, nb_gc + 1):
        T_acl_mem_ccl = average([T_acl[i], T_mem[i], T_ccl[i]],
                                weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl),
                                         Hccl / (Hacl + Hmem + Hccl)])
        i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl[i] * k_H2(lambda_mem[i], T_mem[i], kappa_co)
        i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl[i] * k_O2(lambda_mem[i], T_mem[i], kappa_co)
        i_n += [i_H2 + i_O2]

    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        # The desired air compressor volume flow rate (mol.s-1)                                                         # Warning: considérer le débit minimal du compresseur !
        W_H2_des = 1 / y_H2_in * Sa * i_fc_cell / (2 * F) * (nb_cell * Aact)
        Wacp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc_cell, W_H2_des)                                         # Adjust the desired compressor flow to ensure a minimum flow is maintained.
        W_air_ext_des = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * i_fc_cell / (4 * F) * (nb_cell * Aact)
        Wccp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc_cell, W_air_ext_des)                                         # Adjust the desired compressor flow to ensure a minimum flow is maintained.

        # The desired humidifier volume flow rate at the anode side Wa_v_inj_des (mol.s-1)                              # Warning: considérer le débit minimal du compresseur !
        if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            Prd = Pasm
            W_H2_des = 1 / y_H2_in * Sa * i_fc_cell / (2 * F) * (nb_cell * Aact)
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
        W_H2_des = Sa * (i_fc_cell + max(i_n)) / (2 * F) * (nb_cell * Aact)
        W_H2O_inj_a_des = (Phi_a_des * Psat(T_des) / (Pa_in - Phi_a_des * Psat(T_des))) * W_H2_des

        # At the cathode side
        W_dry_air_des = 1 / y_O2_ext * Sc * (i_fc_cell + max(i_n)) / (4 * F) * (nb_cell * Aact)
        W_H2O_inj_c_des = (Phi_c_des * Psat(T_des) / (Pc_in - Phi_c_des * Psat(T_des))) * W_dry_air_des

    return {'H2': W_H2_des, 'dry_air': W_dry_air_des, 'H2O_inj_a': W_H2O_inj_a_des, 'H2O_inj_c': W_H2O_inj_c_des}


def adjust_compressor_flow_with_minimum(i_fc_cell, Wcp_des):
    """
    Adjusts the desired compressor flow rate to ensure a minimum flow is maintained, based on the current density and
     model parameters.

    Parameters
    ----------
    i_fc_cell : float
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

    if i_fc_cell <= i_cp_min + 3 * delta_i_load_step:
        Wcp_des_adjusted = (
            Wcp_des * i_cp_min / i_fc_cell * (1.0 + math.tanh(4 * (i_fc_cell - (delta_i_load_step / 2)) / (delta_i_load_step / 2))) / 2
            + Wcp_des * (1 - i_cp_min / i_fc_cell) * (1.0 + math.tanh(4 * (i_fc_cell - i_cp_min - (delta_i_load_step / 2)) / (delta_i_load_step / 2))) / 2
        )
        return Wcp_des_adjusted
    else: # For higher current densities, the compressor flow is not adjusted, and so it is faster to return the original value.
        return Wcp_des