# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing constants' value and functions
from configuration.settings import C_O2ref, alpha_c, tau_cp, tau_hum, rho_mem, M_eq, epsilon_cl, F, R, M_H2O, \
    n_cell, Vsm, Vem, A_T, Kp, Kd
from model.flows import calculate_flows
from model.cell_voltage import calculate_eta_c_intermediate_values
from model.heat_transfer import calculate_heat_transfers
from model.control import control_operating_conditions
from modules.transitory_functions import rho_H2O_l, Psat
from modules.dif_eq_modules import dif_eq_int_values, desired_flows


# ______________________Objective function to solve. It gives the system of differential equations______________________

def dydt(t, y, operating_inputs, parameters, solver_variable_names, control_variables):
    """This function gives the system of differential equations to solve.

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
    dydt : list
        List containing the derivative of the solver variables.
    """

    # Creation of the dif_eq dictionary. It is an intermediate calculation to simplify the writing of the code.
    dif_eq = {('d' + key + ' / dt'): 0 for key in solver_variable_names}
    # Creation of the solver_variables dict. It is an intermediate calculation to simplify the writing of the code.
    solver_variables = {}
    for index, key in enumerate(solver_variable_names):
        solver_variables[key] = y[index]

    # Modifications of the operating conditions in real time, if required.
    if parameters["type_control"] != "no_control":
        control_operating_conditions(t, solver_variables, operating_inputs, parameters, control_variables)

    # Intermediate values
    i_fc = operating_inputs['current_density'](t, parameters)
    Mext, Pagc, Pcgc, i_n, Masm, Maem, Mcsm, Mcem, rho_Cp0 = dif_eq_int_values(solver_variables,
                                                                        operating_inputs, control_variables, parameters)
    Wcp_des, Wa_inj_des, Wc_inj_des = desired_flows(solver_variables, control_variables, i_n, i_fc, operating_inputs,
                                                    parameters, Mext)
    eta_c_intermediate_values = calculate_eta_c_intermediate_values(solver_variables, operating_inputs, parameters)

    # Calculation of the flows
    matter_flows_dico = calculate_flows(t, solver_variables, control_variables, i_fc, operating_inputs, parameters)
    heat_flows_dico = calculate_heat_transfers(solver_variables, i_fc, parameters, **matter_flows_dico)

    # Calculation of the dynamic evolutions
    #       Inside the cell
    calculate_dyn_dissoved_water_evolution(dif_eq, **parameters, **matter_flows_dico)
    calculate_dyn_liquid_water_evolution(dif_eq, solver_variables, **parameters, **matter_flows_dico)
    calculate_dyn_vapor_evolution(dif_eq, solver_variables, **parameters, **matter_flows_dico)
    calculate_dyn_H2_O2_N2_evolution(dif_eq, solver_variables, **parameters, **matter_flows_dico)
    calculate_dyn_voltage_evolution(dif_eq, i_fc, **solver_variables, **operating_inputs, **parameters,
                                    **eta_c_intermediate_values)
    calculate_dyn_temperature_evolution(dif_eq, rho_Cp0, **parameters, **heat_flows_dico)
    #       Inside the auxiliary components
    calculate_dyn_manifold_pressure_and_humidity_evolution(dif_eq, Masm, Maem, Mcsm, Mcem, **solver_variables,
                                                           **operating_inputs, **parameters, **matter_flows_dico)
    calculate_dyn_air_compressor_and_humidifier_evolution(dif_eq, Wcp_des, Wa_inj_des, Wc_inj_des,
                                                          **solver_variables, **parameters)
    calculate_dyn_throttle_area_evolution(dif_eq, Pagc, Pcgc, **solver_variables, **operating_inputs, **parameters)

    # dif_eq is converted to dydt because the solver requires an ordered list to work
    dydt = []
    for key in solver_variable_names:
        dydt.append(dif_eq['d' + key + ' / dt'])
    return dydt


# ___________________________Elementary functions which gives specific differential equations___________________________

def calculate_dyn_dissoved_water_evolution(dif_eq, Hmem, Hcl, epsilon_mc, S_abs_acl, S_abs_ccl, J_lambda_acl_mem,
                                           J_lambda_mem_ccl, Sp_acl, Sp_ccl, **kwargs):
    """
    This function calculates the dynamic evolution of the dissolved water in the membrane and the catalyst layers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    Hmem : float
        Thickness of the membrane (m).
    Hcl : float
        Thickness of the catalyst layer (m).
    epsilon_mc : float
        Volume fraction of ionomer in the catalyst layer.
    S_abs_acl : float
        Water absorption in the anode catalyst layer (mol.m-3.s-1)
    S_abs_ccl : float
        Water absorption in the cathode catalyst layer (mol.m-3.s-1)
    J_lambda_acl_mem : float
        Dissolved water flow between the anode catalyst layer and the membrane (mol.m-2.s-1)
    J_lambda_mem_ccl : float
        Dissolved water flow between the membrane and the cathode catalyst layer (mol.m-2.s-1)
    Sp_acl : float
        Water produced in the membrane at the ACL through the chemical reaction and crossover (mol.m-3.s-1)
    Sp_ccl : float
        Water produced in the membrane at the CCL through the chemical reaction and crossover (mol.m-3.s-1)
    """

    dif_eq['dlambda_acl / dt'] = M_eq / (rho_mem * epsilon_mc) * (-J_lambda_acl_mem / Hcl + S_abs_acl + Sp_acl)
    dif_eq['dlambda_mem / dt'] = M_eq / rho_mem * (J_lambda_acl_mem - J_lambda_mem_ccl) / Hmem
    dif_eq['dlambda_ccl / dt'] = M_eq / (rho_mem * epsilon_mc) * (J_lambda_mem_ccl / Hcl + S_abs_ccl + Sp_ccl)


def calculate_dyn_liquid_water_evolution(dif_eq, sv, Hgdl, Hcl, epsilon_gdl, n_gdl, Jl_agc_agdl, Jl_agdl_agdl, Jl_agdl_acl,
                                         Jl_ccl_cgdl, Jl_cgdl_cgdl, Jl_cgdl_cgc, Sl_agdl, Sl_acl, Sl_ccl, Sl_cgdl, **kwargs):
    """
    This function calculates the dynamic evolution of the liquid water in the gas diffusion and catalyst layers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    Hgdl : float
        Thickness of the gas diffusion layer (m).
    Hcl : float
        Thickness of the catalyst layer (m).
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    Jl_agc_agdl : float
        Liquid water flow between the anode gas channel border and the anode GDL (kg.m-2.s-1).
    Jl_agdl_agdl : list
        Liquid water flow between two nodes of the anode GDL (kg.m-2.s-1).
    Jl_agdl_acl : float
        Liquid water flow between the anode GDL and the anode CL (kg.m-2.s-1).
    Jl_ccl_cgdl : float
        Liquid water flow between the cathode CL and the cathode GDL (kg.m-2.s-1).
    Jl_cgdl_cgdl : list
        Liquid water flow between two nodes of the cathode GDL (kg.m-2.s-1).
    Jl_cgdl_cgc : float
        Liquid water flow between the cathode GDL and the cathode gas channel border (kg.m-2.s-1).
    Sl_agdl : list
        Liquid water produced in the anode GDL (kg.m-3.s-1).
    Sl_acl : float
        Liquid water produced in the anode CL (kg.m-3.s-1).
    Sl_ccl : float
        Liquid water produced in the cathode CL (kg.m-3.s-1).
    Sl_cgdl : list
        Liquid water produced in the cathode GDL (kg.m-3.s-1).
    """

    # At the anode side
    #       Inside the AGDL
    dif_eq['ds_agdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_agdl_1']) * epsilon_gdl) * \
                                 ((Jl_agc_agdl - Jl_agdl_agdl[1]) / (Hgdl / n_gdl) + M_H2O * Sl_agdl[1])
    for i in range(2, n_gdl):
        dif_eq[f'ds_agdl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_agdl_{i}']) * epsilon_gdl) * \
                                      ((Jl_agdl_agdl[i - 1] - Jl_agdl_agdl[i]) / (Hgdl / n_gdl) +
                                       M_H2O * Sl_agdl[i])
    dif_eq[f'ds_agdl_{n_gdl} / dt'] = 1 / (rho_H2O_l(sv[f'T_agdl_{n_gdl}']) * epsilon_gdl) * \
                                      ((Jl_agdl_agdl[n_gdl - 1] - Jl_agdl_acl) / (Hgdl / n_gdl) +
                                       M_H2O * Sl_agdl[n_gdl])
    #      Inside the ACL
    dif_eq['ds_acl / dt'] = 1 / (rho_H2O_l(sv['T_acl']) * epsilon_cl) * (Jl_agdl_acl / Hcl + M_H2O * Sl_acl)

    # At the cathode side
    #       Inside the CCL
    dif_eq['ds_ccl / dt'] = 1 / (rho_H2O_l(sv['T_ccl']) * epsilon_cl) * (- Jl_ccl_cgdl / Hcl + M_H2O * Sl_ccl)
    #       Inside the CGDL
    dif_eq['ds_cgdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_cgdl_1']) * epsilon_gdl) * \
                               ((Jl_ccl_cgdl - Jl_cgdl_cgdl[1]) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[1])
    for i in range(2, n_gdl):
        dif_eq[f'ds_cgdl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_cgdl_{i}']) * epsilon_gdl) * \
                                      ((Jl_cgdl_cgdl[i - 1] - Jl_cgdl_cgdl[i]) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[i])
    dif_eq[f'ds_cgdl_{n_gdl} / dt'] = 1 / (rho_H2O_l(sv[f'T_cgdl_{n_gdl}']) * epsilon_gdl) * \
                                     ((Jl_cgdl_cgdl[n_gdl - 1] - Jl_cgdl_cgc) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[n_gdl])


def calculate_dyn_vapor_evolution(dif_eq, sv, Hgdl, Hcl, Hgc, Lgc, epsilon_gdl, n_gdl, Jv_a_in, Jv_a_out, Jv_c_in,
                                  Jv_c_out, Jv_agc_agdl, Jv_agdl_agdl, Jv_agdl_acl, S_abs_acl, S_abs_ccl, Jv_ccl_cgdl,
                                  Jv_cgdl_cgdl, Jv_cgdl_cgc, Sv_agdl, Sv_acl, Sv_ccl, Sv_cgdl, **kwargs):
    """This function calculates the dynamic evolution of the vapor in the gas channels, the gas diffusion layers and the
    catalyst layers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    Hgdl : float
        Thickness of the gas diffusion layer (m).
    Hcl : float
        Thickness of the catalyst layer (m).
    Hgc : float
        Thickness of the gas channel (m).
    Lgc : float
        Length of the gas channel (m).
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    Jv_a_in : float
        Water vapor flow at the anode inlet (mol.m-2.s-1).
    Jv_a_out : float
        Water vapor flow at the anode outlet (mol.m-2.s-1).
    Jv_c_in : float
        Water vapor flow at the cathode inlet (mol.m-2.s-1).
    Jv_c_out : float
        Water vapor flow at the cathode outlet (mol.m-2.s-1).
    Jv_agc_agdl : float
        Water vapor flow between the anode gas channel and the anode GDL (mol.m-2.s-1).
    Jv_agdl_agdl : list
        Water vapor flow between two nodes of the anode GDL (mol.m-2.s-1).
    Jv_agdl_acl : float
        Water vapor flow between the anode GDL and the anode CL (mol.m-2.s-1).
    S_abs_acl: float
        Water vapor absorption in the anode CL (mol.m-3.s-1).
    S_abs_ccl: float
        Water vapor absorption in the cathode CL (mol.m-3.s-1).
    Jv_ccl_cgdl : float
        Water vapor flow between the cathode CL and the cathode GDL (mol.m-2.s-1).
    Jv_cgdl_cgdl : list
        Water vapor flow between two nodes of the cathode GDL (mol.m-2.s-1).
    Jv_cgdl_cgc : float
        Water vapor flow between the cathode GDL and the cathode gas channel (mol.m-2.s-1).
    Sv_agdl : list
        Water vapor produced in the anode GDL (mol.m-3.s-1).
    Sv_acl : float
        Water vapor produced in the anode CL (mol.m-3.s-1).
    Sv_ccl : float
        Water vapor produced in the cathode CL (mol.m-3.s-1).
    Sv_cgdl : list
        Water vapor produced in the cathode GDL (mol.m-3.s-1).
    """

    # At the anode side
    #       Inside the AGC
    dif_eq['dC_v_agc / dt'] = (Jv_a_in - Jv_a_out) / Lgc - Jv_agc_agdl / Hgc
    #       Inside the AGDL
    dif_eq['dC_v_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                 ((Jv_agc_agdl - Jv_agdl_agdl[1]) / (Hgdl / n_gdl) + Sv_agdl[1])
    for i in range(2, n_gdl):
        dif_eq[f'dC_v_agdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{i}'])) * \
                                        ((Jv_agdl_agdl[i - 1] - Jv_agdl_agdl[i]) / (Hgdl / n_gdl) + Sv_agdl[i])
    dif_eq[f'dC_v_agdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{n_gdl}'])) * \
                                        ((Jv_agdl_agdl[n_gdl - 1] - Jv_agdl_acl) / (Hgdl / n_gdl) + Sv_agdl[n_gdl])
    #       Inside the ACL
    dif_eq['dC_v_acl / dt'] = 1 / (epsilon_cl * (1 - sv['s_acl'])) * (Jv_agdl_acl / Hcl - S_abs_acl + Sv_acl)

    # At the cathode side
    #       Inside the CCL
    dif_eq['dC_v_ccl / dt'] = 1 / (epsilon_cl * (1 - sv['s_ccl'])) * (- Jv_ccl_cgdl / Hcl - S_abs_ccl + Sv_ccl)
    #       Inside the CGDL
    dif_eq['dC_v_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * \
                                 ((Jv_ccl_cgdl - Jv_cgdl_cgdl[1]) / (Hgdl / n_gdl) + Sv_cgdl[1])
    for i in range(2, n_gdl):
        dif_eq[f'dC_v_cgdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{i}'])) * \
                                        ((Jv_cgdl_cgdl[i - 1] - Jv_cgdl_cgdl[i]) / (Hgdl / n_gdl) + Sv_cgdl[i])
    dif_eq[f'dC_v_cgdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{n_gdl}'])) * \
                                        ((Jv_cgdl_cgdl[n_gdl - 1] - Jv_cgdl_cgc) / (Hgdl / n_gdl) + Sv_cgdl[n_gdl])
    #       Inside the CGC
    dif_eq['dC_v_cgc / dt'] = (Jv_c_in - Jv_c_out) / Lgc + Jv_cgdl_cgc / Hgc


def calculate_dyn_H2_O2_N2_evolution(dif_eq, sv, Hgdl, Hcl, Hgc, Lgc, epsilon_gdl, n_gdl, J_H2_in, J_H2_out, J_O2_in,
                                     J_O2_out, J_N2_in, J_N2_out, J_H2_agc_agdl, J_H2_agdl_agdl, J_H2_agdl_acl,
                                     J_O2_ccl_cgdl, J_O2_cgdl_cgdl, J_O2_cgdl_cgc, S_H2_acl, S_O2_ccl, **kwargs):
    """This function calculates the dynamic evolution of the hydrogen, oxygen and nitrogen in the gas channels, the gas
    diffusion layers and the catalyst layers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    Hgdl : float
        Thickness of the gas diffusion layer (m).
    Hcl : float
        Thickness of the catalyst layer (m).
    Hgc : float
        Thickness of the gas channel (m).
    Lgc : float
        Length of the gas channel (m).
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    J_H2_in : float
        Hydrogen flow at the anode inlet (mol.m-2.s-1).
    J_H2_out : float
        Hydrogen flow at the anode outlet (mol.m-2.s-1).
    J_O2_in : float
        Oxygen flow at the cathode inlet (mol.m-2.s-1).
    J_O2_out : float
        Oxygen flow at the cathode outlet (mol.m-2.s-1).
    J_N2_in : float
        Nitrogen flow at the cathode inlet (mol.m-2.s-1).
    J_N2_out : float
        Nitrogen flow at the cathode outlet (mol.m-2.s-1).
    J_H2_agc_agdl : float
        Hydrogen flow between the anode gas channel and the anode GDL (mol.m-2.s-1).
    J_H2_agdl_agdl : list
        Hydrogen flow between two nodes of the anode GDL (mol.m-2.s-1).
    J_H2_agdl_acl : float
        Hydrogen flow between the anode GDL and the anode CL (mol.m-2.s-1).
    J_O2_ccl_cgdl : float
        Oxygen flow between the cathode CL and the cathode GDL (mol.m-2.s-1).
    J_O2_cgdl_cgdl : list
        Oxygen flow between two nodes of the cathode GDL (mol.m-2.s-1).
    J_O2_cgdl_cgc : float
        Oxygen flow between the cathode GDL and the cathode gas channel (mol.m-2.s-1).
    S_H2_acl : float
        Hydrogen consumed in the anode CL (mol.m-3.s-1).
    S_O2_ccl : float
        Oxygen consumed in the cathode CL (mol.m-3.s-1).
    """

    # At the anode side
    #      Inside the AGC
    dif_eq['dC_H2_agc / dt'] = (J_H2_in - J_H2_out) / Lgc - J_H2_agc_agdl / Hgc
    #      Inside the AGDL
    dif_eq['dC_H2_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                  (J_H2_agc_agdl - J_H2_agdl_agdl[1]) / (Hgdl / n_gdl)
    for i in range(2, n_gdl):
        dif_eq[f'dC_H2_agdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{i}'])) * \
                                         (J_H2_agdl_agdl[i - 1] - J_H2_agdl_agdl[i]) / (Hgdl / n_gdl)
    dif_eq[f'dC_H2_agdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{n_gdl}'])) * \
                                         (J_H2_agdl_agdl[n_gdl - 1] - J_H2_agdl_acl) / (Hgdl / n_gdl)
    #      Inside the ACL
    dif_eq['dC_H2_acl / dt'] = 1 / (epsilon_cl * (1 - sv['s_acl'])) * (J_H2_agdl_acl / Hcl + S_H2_acl)

    # At the cathode side
    #      Inside the CCL
    dif_eq['dC_O2_ccl / dt'] = 1 / (epsilon_cl * (1 - sv['s_ccl'])) * (-J_O2_ccl_cgdl / Hcl + S_O2_ccl)
    #      Inside the CGDL
    dif_eq['dC_O2_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * \
                                  (J_O2_ccl_cgdl - J_O2_cgdl_cgdl[1]) / (Hgdl / n_gdl)
    for i in range(2, n_gdl):
        dif_eq[f'dC_O2_cgdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{i}'])) * \
                                         (J_O2_cgdl_cgdl[i - 1] - J_O2_cgdl_cgdl[i]) / (Hgdl / n_gdl)
    dif_eq[f'dC_O2_cgdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{n_gdl}'])) * \
                                         (J_O2_cgdl_cgdl[n_gdl - 1] - J_O2_cgdl_cgc) / (Hgdl / n_gdl)
    #      Inside the CGC
    dif_eq['dC_O2_cgc / dt'] = (J_O2_in - J_O2_out) / Lgc + J_O2_cgdl_cgc / Hgc

    # Inside the whole cell
    dif_eq['dC_N2 / dt'] = (J_N2_in - J_N2_out) / Lgc


def calculate_dyn_temperature_evolution(dif_eq, rho_Cp0, Hgdl, Hcl, Hmem, n_gdl, Jt, Q_r, Q_sorp, Q_liq, Q_p, Q_e,
                                        **kwargs):
    """
    This function calculates the dynamic evolution of the temperature in the fuel cell.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    rho_Cp0 : dict
        Volumetric heat capacity of the different components of the fuel cell system, in J.m-3.K-1.
    Hgdl : float
        Thickness of the gas diffusion layer, in m.
    Hcl : float
        Thickness of the catalyst layer, in m.
    Hmem : float
        Thickness of the membrane, in m.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    Jt : dict
        Heat flows occuring inside the fuel cell system, J.m-2.s-1.
    Q_r : dict
        Heat dissipated by the electrochemical reaction 2*H2 + O2 -> 2*H2O, in J.m-3.s-1.
    Q_sorp : dict
        Heat dissipated by the absorption of water from the CL to the membrane, in J.m-3.s-1.
    Q_liq : dict
        Heat dissipated by the evaporation of liquid water, in J.m-3.s-1.
    Q_p : dict
        Heat dissipated by the ionic currents (Joule heating + Ohm's law), in J.m-3.s-1.
    Q_e : dict
        Heat dissipated by the electric currents (Joule heating + Ohm's law), in J.m-3.s-1.
    """

    # At the anode side
    #       Inside the AGC
    dif_eq['dT_agc / dt'] = 0  # Dirichlet boundary condition. T_agc is initialized to T_fc and remains constant.
    #       Inside the AGDL
    dif_eq['dT_agdl_1 / dt'] = (1 / rho_Cp0['agdl_1']) * \
                               ( (Jt['agc_agdl'] - Jt['agdl_agdl_1']) / (Hgdl / n_gdl)  +
                                 Q_liq['agdl_1'] + Q_e['agdl_1'] )
    for i in range(2, n_gdl):
        dif_eq[f'dT_agdl_{i} / dt'] = (1 / rho_Cp0[f'agdl_{i}']) * \
                                      ( (Jt[f'agdl_agdl_{i - 1}'] - Jt[f'agdl_agdl_{i}']) / (Hgdl / n_gdl) +
                                        Q_liq[f'agdl_{i}'] + Q_e[f'agdl_{i}'] )
    dif_eq[f'dT_agdl_{n_gdl} / dt'] = (1 / rho_Cp0[f'agdl_{n_gdl}']) * \
                                      ( (Jt[f'agdl_agdl_{n_gdl - 1}'] - Jt['agdl_acl']) / (Hgdl / n_gdl) +
                                        Q_liq[f'agdl_{n_gdl}'] + Q_e[f'agdl_{n_gdl}'] )
    #      Inside the ACL
    dif_eq['dT_acl / dt'] = (1 / rho_Cp0['acl']) * \
                            ( (Jt['agdl_acl'] - Jt['acl_mem']) / Hcl +
                              Q_r['acl'] + Q_sorp['acl'] + Q_liq['acl'] + Q_e['acl'] )

    # Inside the membrane
    dif_eq['dT_mem / dt'] = (1 / rho_Cp0['mem']) * \
                            ( (Jt['acl_mem'] - Jt['mem_ccl']) / Hmem +
                              Q_p['mem'] )

    # At the cathode side
    #       Inside the CCL
    dif_eq['dT_ccl / dt'] = (1 / rho_Cp0['ccl']) * \
                            ( (Jt['mem_ccl'] - Jt['ccl_cgdl']) / Hcl +
                              Q_r['ccl'] + Q_sorp['acl'] + Q_liq['ccl'] + Q_p['ccl'] + Q_e['ccl'] )
    #       Inside the CGDL
    dif_eq['dT_cgdl_1 / dt'] = (1 / rho_Cp0['cgdl_1']) * \
                               ( (Jt['ccl_cgdl'] - Jt['cgdl_cgdl_1']) / (Hgdl / n_gdl) +
                                 Q_liq['cgdl_1'] + Q_e['cgdl_1'] )
    for i in range(2, n_gdl):
        dif_eq[f'dT_cgdl_{i} / dt'] = (1 / rho_Cp0[f'cgdl_{i}']) * \
                                      ( (Jt[f'cgdl_cgdl_{i - 1}'] - Jt[f'cgdl_cgdl_{i}']) / (Hgdl / n_gdl) +
                                        Q_liq[f'cgdl_{i}'] + Q_e[f'cgdl_{i}'] )
    dif_eq[f'dT_cgdl_{n_gdl} / dt'] = (1 / rho_Cp0[f'cgdl_{n_gdl}']) * \
                                      ( (Jt[f'cgdl_cgdl_{n_gdl - 1}'] - Jt['cgdl_cgc']) / (Hgdl / n_gdl) +
                                        Q_liq[f'cgdl_{n_gdl}'] + Q_e[f'cgdl_{n_gdl}'] )
    #       Inside the CCG
    dif_eq['dT_cgc / dt'] = 0  # Dirichlet boundary condition. T_cgc is initialized to T_fc and remains constant.


def calculate_dyn_voltage_evolution(dif_eq, i_fc, C_O2_ccl, T_ccl, eta_c, Hcl, i0_c_ref, kappa_c, C_scl, i_n, f_drop,
                                    **kwargs):
    """This function calculates the dynamic evolution of the cell overpotential eta_c.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    i_fc : float
        Fuel cell current density (A.m-2).
    C_O2_ccl : float
        Oxygen concentration in the cathode catalyst layer (mol.m-3).
    T_ccl : float
        Fuel cell temperature in the cathode catalyst layer (K).
    eta_c : float
        Cell overpotential (V).
    Hcl : float
        Thickness of the catalyst layer (m).
    i0_c_ref : float
        Reference exchange current density at the cathode (A.m-2).
    kappa_c : float
        Overpotential correction exponent.
    C_scl : float
        Volumetric space-charge layer capacitance (F.m-3).
    i_n : float
        Crossover current density (A.m-2).
    f_drop : float
        Liquid water induced voltage drop function.
    """

    dif_eq['deta_c / dt'] = 1 / (C_scl * Hcl) * ((i_fc + i_n) - i0_c_ref * (C_O2_ccl / C_O2ref) ** kappa_c *
                                                 math.exp(f_drop * alpha_c * F / (R * T_ccl) * eta_c))


def calculate_dyn_manifold_pressure_and_humidity_evolution(dif_eq, Masm, Maem, Mcsm, Mcem, T_des, Hgc, Wgc,
                                                           type_auxiliary, Jv_a_in, Jv_a_out, Jv_c_in, Jv_c_out,
                                                           Wasm_in, Wasm_out, Waem_in, Waem_out, Wcsm_in, Wcsm_out,
                                                           Wcem_in, Wcem_out, Ware, Wv_asm_in, Wv_aem_out, Wv_csm_in,
                                                           Wv_cem_out, **kwargs):
    """This function calculates the dynamic evolution of the pressure and humidity inside the manifolds.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    Masm : float
        Molar mass of all the gaseous species inside the anode supply manifold (kg.mol-1).
    Maem : float
        Molar mass of all the gaseous species inside the anode exhaust manifold (kg.mol-1).
    Mcsm : float
        Molar mass of all the gaseous species inside the cathode supply manifold (kg.mol-1).
    Mcem : float
        Molar mass of all the gaseous species inside the cathode exhaust manifold (kg.mol-1).
    T_des : float
        Fuel cell temperature (K).
    Hgc : float
        Thickness of the gas channel (m).
    Wgc : float
        Width of the gas channel (m).
    type_auxiliary : str
        Type of auxiliary components used in the fuel cell model.
    Jv_a_in : float
        Water vapor flow at the anode inlet (mol.m-2.s-1).
    Jv_a_out : float
        Water vapor flow at the anode outlet (mol.m-2.s-1).
    Jv_c_in : float
        Water vapor flow at the cathode inlet (mol.m-2.s-1).
    Jv_c_out : float
        Water vapor flow at the cathode outlet (mol.m-2.s-1).
    Wasm_in : float
        Flow at the anode supply manifold inlet (kg.s-1).
    Wasm_out : float
        Flow at the anode supply manifold outlet (kg.s-1).
    Waem_in : float
        Flow at the anode exhaust manifold inlet (kg.s-1).
    Waem_out : float
        Flow at the anode exhaust manifold outlet (kg.s-1).
    Wcsm_in : float
        Flow at the cathode supply manifold inlet (kg.s-1).
    Wcsm_out : float
        Flow at the cathode supply manifold outlet (kg.s-1).
    Wcem_in : float
        Flow at the cathode exhaust manifold inlet (kg.s-1).
    Wcem_out : float
        Flow at the cathode exhaust manifold outlet (kg.s-1).
    Ware : float
        Anode side recirculation flow (kg.s-1).
    Wv_asm_in : float
        Vapor flow at the anode supply manifold inlet (mol.s-1).
    Wv_aem_out : float
        Vapor flow at the anode external manifold outlet (mol.s-1).
    Wv_csm_in : float
        Vapor flow at the cathode supply manifold inlet (mol.s-1).
    Wv_cem_out : float
        Vapor flow at the cathode external manifold outlet (mol.s-1).
    """

    # Pressure evolution inside the manifolds
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
        dif_eq['dPasm / dt'] = (Wasm_in + Ware - n_cell * Wasm_out) / (Vsm * Masm) * R * T_des
        dif_eq['dPaem / dt'] = (n_cell * Waem_in - Ware - Waem_out) / (Vem * Maem) * R * T_des
        dif_eq['dPcsm / dt'] = (Wcsm_in - n_cell * Wcsm_out) / (Vsm * Mcsm) * R * T_des
        dif_eq['dPcem / dt'] = (n_cell * Wcem_in - Wcem_out) / (Vem * Mcem) * R * T_des
    elif type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dPasm / dt'] = (Wasm_in - n_cell * Wasm_out) / (Vsm * Masm) * R * T_des
        dif_eq['dPaem / dt'] = (n_cell * Waem_in - Waem_out) / (Vem * Maem) * R * T_des
        dif_eq['dPcsm / dt'] = (Wcsm_in - n_cell * Wcsm_out) / (Vsm * Mcsm) * R * T_des
        dif_eq['dPcem / dt'] = (n_cell * Wcem_in - Wcem_out) / (Vem * Mcem) * R * T_des
    else:  # elif type_auxiliary == "no_auxiliary":
        dif_eq['dPasm / dt'], dif_eq['dPaem / dt'], dif_eq['dPcsm / dt'], dif_eq['dPcem / dt'] = 0, 0, 0, 0

    # Humidity evolution inside the manifolds
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
        dif_eq['dPhi_asm / dt'] = (Wv_asm_in - Jv_a_in * Hgc * Wgc * n_cell) / Vsm * R * T_des / Psat(T_des)
        dif_eq['dPhi_aem / dt'] = (Jv_a_out * Hgc * Wgc * n_cell - Wv_asm_in - Wv_aem_out) / Vem * R * T_des / Psat(T_des)
        dif_eq['dPhi_csm / dt'] = (Wv_csm_in - Jv_c_in * Hgc * Wgc * n_cell) / Vsm * R * T_des / Psat(T_des)
        dif_eq['dPhi_cem / dt'] = (Jv_c_out * Hgc * Wgc * n_cell - Wv_cem_out) / Vem * R * T_des / Psat(T_des)
    elif type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dPhi_asm / dt'] = (Wv_asm_in - Jv_a_in * Hgc * Wgc * n_cell) / Vsm * R * T_des / Psat(T_des)
        dif_eq['dPhi_aem / dt'] = (Jv_a_out * Hgc * Wgc * n_cell - Wv_aem_out) / Vem * R * T_des / Psat(T_des)
        dif_eq['dPhi_csm / dt'] = (Wv_csm_in - Jv_c_in * Hgc * Wgc * n_cell) / Vsm * R * T_des / Psat(T_des)
        dif_eq['dPhi_cem / dt'] = (Jv_c_out * Hgc * Wgc * n_cell - Wv_cem_out) / Vem * R * T_des / Psat(T_des)
    else:  # elif type_auxiliary == "no_auxiliary":
        dif_eq['dPhi_asm / dt'], dif_eq['dPhi_aem / dt'], dif_eq['dPhi_csm / dt'], dif_eq['dPhi_cem / dt'] = 0, 0, 0, 0


def calculate_dyn_air_compressor_and_humidifier_evolution(dif_eq, Wcp_des, Wa_inj_des, Wc_inj_des, type_auxiliary, Wcp,
                                                          Wa_inj, Wc_inj, **kwargs):
    """This function calculates the dynamic evolution of the air compressor and the humidifiers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    Wcp_des : float
        Desired air compressor flow rate (kg.s-1).
    Wa_inj_des : float
        Desired anode humidifier flow rate (kg.s-1).
    Wc_inj_des : float
        Desired cathode humidifier flow rate (kg.s-1).
    type_auxiliary : str
        Type of auxiliary components used in the fuel cell model.
    Wcp : float
        Air compressor flow rate (kg.s-1).
    Wa_inj : float
        Anode humidifier flow rate (kg.s-1).
    Wc_inj : float
        Cathode humidifier flow rate (kg.s-1).
    """

    # Air compressor evolution
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
        dif_eq['dWcp / dt'] = (Wcp_des - Wcp) / tau_cp  # Estimation at the first order.
    elif type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dWcp / dt'] = (Wcp_des - Wcp) / tau_cp  # Estimation at the first order.
    else:  # elif type_auxiliary == "no_auxiliary":
        dif_eq['dWcp / dt'] = 0

    # Anode and cathode humidifiers evolution
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation":
        dif_eq['dWa_inj / dt'] = 0
        dif_eq['dWc_inj / dt'] = (Wc_inj_des - Wc_inj) / tau_hum  # Estimation at the first order.
    elif type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dWa_inj / dt'] = (Wa_inj_des - Wa_inj) / tau_hum  # Estimation at the first order.
        dif_eq['dWc_inj / dt'] = (Wc_inj_des - Wc_inj) / tau_hum  # Estimation at the first order.
    else:  # elif type_auxiliary == "no_auxiliary":
        dif_eq['dWa_inj / dt'], dif_eq['dWc_inj / dt'] = 0, 0


def calculate_dyn_throttle_area_evolution(dif_eq, Pagc, Pcgc, T_agc, T_cgc, Abp_a, Abp_c, Pa_des, Pc_des,
                                          type_auxiliary, **kwargs):
    """This function calculates the dynamic evolution of the throttle area inside the anode and cathode auxiliaries.
    This function has to be executed after 'calculate_dyn_vapor_evolution' and 'calculate_dyn_H2_O2_N2_evolution'.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    Pagc : float
        Pressure inside the anode gas channel (Pa).
    Pcgc : float
        Pressure inside the cathode gas channel (Pa).
    T_agc : float
        Fuel cell temperature in the anode gas channel (K).
    T_cgc : float
        Fuel cell temperature in the cathode gas channel (K).
    Abp_a : float
        Throttle area inside the anode auxiliaries (m2).
    Abp_c : float
        Throttle area inside the cathode auxiliaries (m2).
    Pa_des : float
        Desired pressure inside the anode gas channel (Pa).
    Pc_des : float
        Desired pressure inside the cathode gas channel (Pa).
    type_auxiliary : str
        Type of auxiliary components used in the fuel cell model.
    """

    # Calculation of the pressure derivative inside the gas channels
    dPagcdt = (dif_eq['dC_v_agc / dt'] + dif_eq['dC_H2_agc / dt']) * R * T_agc
    dPcgcdt = (dif_eq['dC_v_cgc / dt'] + dif_eq['dC_O2_cgc / dt'] + dif_eq['dC_N2 / dt']) * R * T_cgc

    # Throttle area evolution inside the anode auxiliaries
    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dAbp_a / dt'] = - Kp * (Pa_des - Pagc) + Kd * dPagcdt  # PD controller
        if Abp_a > A_T and dif_eq['dAbp_a / dt'] > 0:  # The throttle area cannot be higher than the maximum value
            dif_eq['dAbp_a / dt'] = 0
        elif Abp_a < 0 and dif_eq['dAbp_a / dt'] < 0:  # The throttle area cannot be lower than 0
            dif_eq['dAbp_a / dt'] = 0
    else:  # elif type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
           #      type_auxiliary == "no_auxiliary":
        dif_eq['dAbp_a / dt'] = 0  # The throttle area is not considered

    # Throttle area evolution inside the cathode auxiliaries
    if type_auxiliary == "forced-convective_cathode_with_anodic_recirculation" or \
       type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
        dif_eq['dAbp_c / dt'] = - Kp * (Pc_des - Pcgc) + Kd * dPcgcdt  # PD controller
        if Abp_c > A_T and dif_eq['dAbp_c / dt'] > 0:  # The throttle area cannot be higher than the maximum value
            dif_eq['dAbp_c / dt'] = 0
        elif Abp_c < 0 and dif_eq['dAbp_c / dt'] < 0:  # The throttle area cannot be lower than 0
            dif_eq['dAbp_c / dt'] = 0
    else:  # elif type_auxiliary == "no_auxiliary":
        dif_eq['dAbp_a / dt'] = 0  # The throttle area is not considered
