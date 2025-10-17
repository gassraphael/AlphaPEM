# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing constants' value and functions
from configuration.settings import C_O2ref, alpha_c, rho_mem, M_eq, F, R, M_H2O
from modules.transitory_functions import rho_H2O_l


# ____________________________________________________Main functions____________________________________________________


def calculate_dyn_dissoved_water_evolution_inside_MEA(dif_eq, Hmem, Hacl, Hccl, epsilon_mc, S_abs_acl, S_abs_ccl,
                                                      J_lambda_acl_mem, J_lambda_mem_ccl, Sp_acl, Sp_ccl, **kwargs):
    """
    This function calculates the dynamic evolution of the dissolved water in the membrane and the catalyst layers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    Hmem : float
        Thickness of the membrane (m).
    Hacl : float
        Thickness of the anode catalyst layer (m).
    Hccl : float
        Thickness of the cathode catalyst layer (m).
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

    dif_eq['dlambda_acl / dt'] = M_eq / (rho_mem * epsilon_mc) * (-J_lambda_acl_mem / Hacl + S_abs_acl + Sp_acl)
    dif_eq['dlambda_mem / dt'] = M_eq / rho_mem * (J_lambda_acl_mem - J_lambda_mem_ccl) / Hmem
    dif_eq['dlambda_ccl / dt'] = M_eq / (rho_mem * epsilon_mc) * (J_lambda_mem_ccl / Hccl + S_abs_ccl + Sp_ccl)


def calculate_dyn_liquid_water_evolution_inside_MEA(dif_eq, sv, Hgdl, Hmpl, Hacl, Hccl, epsilon_gdl, epsilon_mpl,
                                                    epsilon_cl, epsilon_atl, epsilon_ctl, Htl, n_gdl, n_tl, n_mpl,
                                                    Jl_agc_agdl, Jl_agdl_agdl, Jl_agdl_atl, Jl_atl_atl, Jl_atl_ampl,
                                                    Jl_ampl_ampl, Jl_ampl_acl, Jl_ccl_cmpl, Jl_cmpl_cmpl, Jl_cmpl_ctl,
                                                    Jl_ctl_ctl, Jl_ctl_cgdl, Jl_cgdl_cgdl, Jl_cgdl_cgc, Sl_agdl, Sl_atl,
                                                    Sl_ampl, Sl_acl, Sl_ccl, Sl_cmpl, Sl_ctl, Sl_cgdl, **kwargs):
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
    Hmpl : float
        Thickness of the microporous layer (m).
    Hacl : float
        Thickness of the anode catalyst layer (m).
    Hccl : float
        Thickness of the cathode catalyst layer (m).
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_mpl : float
        Anode/cathode MPL porosity.
    epsilon_cl : float
        Anode/cathode CL porosity.
    epsilon_atl : list
        Anode transition layer porosities.
    epsilon_ctl : list
        Cathode transition layer porosities.
    Htl : float
        Thickness of the anode transition layer (m).
    n_gdl : int
        Number of model nodes placed inside each GDL.
    n_tl : int
        Number of model nodes placed inside the anode transition layer.
    n_mpl : int
        Number of model nodes placed inside each MPL.
    Jl_agc_agdl : float
        Liquid water flow between the anode gas channel border and the anode GDL (kg.m-2.s-1).
    Jl_agdl_agdl : list
        Liquid water flow between two nodes of the anode GDL (kg.m-2.s-1).
    Jl_agdl_atl : float
        Liquid water flow between the last node of the anode GDL and the first node of the transition layer (kg.m-2.s-1).
    Jl_atl_atl : list
        Liquid water flow between two nodes of the anode transition layer (kg.m-2.s-1).
    Jl_atl_ampl : float
        Liquid water flow between the anode transition layer and the anode MPL (kg.m-2.s-1).
    Jl_ampl_ampl : list
        Liquid water flow between two nodes of the anode MPL (kg.m-2.s-1).
    Jl_ampl_acl : float
        Liquid water flow between the anode microporous layer and the anode catalyst layer (kg.m-2.s-1).
    Jl_ccl_cmpl : list
        Liquid water flow between the cathode catalyst layer and the cathode microporous layer (kg.m-2.s-1).
    Jl_cmpl_cmpl: list
        Liquid water flow between two nodes of the cathode MPL (kg.m-2.s-1).
    Jl_cmpl_ctl : float
        Liquid water flow between the last node of the cathode microporous layer and the first node of the cathode transition layer (kg.m-2.s-1).
    Jl_ctl_cgdl : float
        Liquid water flow between the last node of the cathode transition layer and the first node of the cathode GDL (kg.m-2.s-1).
    Jl_cgdl_cgdl : list
        Liquid water flow between two nodes of the cathode GDL (kg.m-2.s-1).
    Jl_cgdl_cgc : float
        Liquid water flow between the cathode GDL and the cathode gas channel border (kg.m-2.s-1).
    Sl_agdl : list
        Liquid water produced in the anode GDL (kg.m-3.s-1).
    Sl_atl : list
        Liquid water produced in the anode transition layer (kg.m-3.s-1).
    Sl_ampl : list
        Liquid water produced in the anode MPL (kg.m-3.s-1).
    Sl_acl : float
        Liquid water produced in the anode CL (kg.m-3.s-1).
    Sl_ccl : float
        Liquid water produced in the cathode CL (kg.m-3.s-1).
    Sl_cmpl : list
        Liquid water produced in the cathode MPL (kg.m-3.s-1).
    Sl_ctl : list
        Liquid water produced in the cathode transition layer (kg.m-3.s-1).
    Sl_cgdl : list
        Liquid water produced in the cathode GDL (kg.m-3.s-1).
    """

    # At the anode side
    #       Inside the AGDL
    if n_gdl == 1:
        dif_eq['ds_agdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_agdl_1']) * epsilon_gdl) * \
                                   ((Jl_agc_agdl - Jl_agdl_atl) / Hgdl + M_H2O * Sl_agdl[1])
    elif n_gdl == 2:
        dif_eq['ds_agdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_agdl_1']) * epsilon_gdl) * \
                                   ((Jl_agc_agdl - Jl_agdl_agdl[1]) / (Hgdl / n_gdl) + M_H2O * Sl_agdl[1])
        dif_eq['ds_agdl_2 / dt'] = 1 / (rho_H2O_l(sv['T_agdl_2']) * epsilon_gdl) * \
                                   ((Jl_agdl_agdl[1] - Jl_agdl_atl) / (Hgdl / n_gdl) + M_H2O * Sl_agdl[2])
    else: # n_gdl > 2
        dif_eq['ds_agdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_agdl_1']) * epsilon_gdl) * \
                                     ((Jl_agc_agdl - Jl_agdl_agdl[1]) / (Hgdl / n_gdl) + M_H2O * Sl_agdl[1])
        for i in range(2, n_gdl):
            dif_eq[f'ds_agdl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_agdl_{i}']) * epsilon_gdl) * \
                                          ((Jl_agdl_agdl[i - 1] - Jl_agdl_agdl[i]) / (Hgdl / n_gdl) +
                                           M_H2O * Sl_agdl[i])
        dif_eq[f'ds_agdl_{n_gdl} / dt'] = 1 / (rho_H2O_l(sv[f'T_agdl_{n_gdl}']) * epsilon_gdl) * \
                                          ((Jl_agdl_agdl[n_gdl - 1] - Jl_agdl_atl) / (Hgdl / n_gdl) +
                                           M_H2O * Sl_agdl[n_gdl])
    #      Inside the ATL
    if n_tl == 2:
        dif_eq['ds_atl_1 / dt'] = 1 / (rho_H2O_l(sv['T_atl_1']) * epsilon_atl[1]) * \
                                   ((Jl_agdl_atl - Jl_atl_atl[1]) / (Htl / n_tl) + M_H2O * Sl_atl[1])
        dif_eq['ds_atl_2 / dt'] = 1 / (rho_H2O_l(sv['T_atl_2']) * epsilon_atl[2]) * \
                                   ((Jl_atl_atl[1] - Jl_atl_ampl) / (Htl / n_tl) + M_H2O * Sl_atl[2])
    else: # n_tl > 2
        dif_eq['ds_atl_1 / dt'] = 1 / (rho_H2O_l(sv['T_atl_1']) * epsilon_atl[1]) * \
                                    ((Jl_agdl_atl - Jl_atl_atl[1]) / (Htl / n_tl) + M_H2O * Sl_atl[1])
        for i in range(2, n_tl):
            dif_eq[f'ds_atl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_atl_{i}']) * epsilon_atl[i]) * \
                                          ((Jl_atl_atl[i - 1] - Jl_atl_atl[i]) / (Htl / n_tl) + M_H2O * Sl_atl[i])
        dif_eq[f'ds_atl_{n_tl} / dt'] = 1 / (rho_H2O_l(sv[f'T_atl_{n_tl}']) * epsilon_atl[n_tl]) * \
                                            ((Jl_atl_atl[n_tl - 1] - Jl_atl_ampl) / (Htl / n_tl) + M_H2O * Sl_atl[n_tl])
    #      Inside the AMPL
    if n_mpl == 1:
        dif_eq['ds_ampl_1 / dt'] = 1 / (rho_H2O_l(sv['T_ampl_1']) * epsilon_mpl) * \
                                   ((Jl_atl_ampl - Jl_ampl_acl) / Hmpl + M_H2O * Sl_ampl[1])
    elif n_mpl == 2:
        dif_eq['ds_ampl_1 / dt'] = 1 / (rho_H2O_l(sv['T_ampl_1']) * epsilon_mpl) * \
                                   ((Jl_atl_ampl - Jl_ampl_ampl[1]) / (Hmpl / n_mpl) + M_H2O * Sl_ampl[1])
        dif_eq['ds_ampl_2 / dt'] = 1 / (rho_H2O_l(sv['T_ampl_2']) * epsilon_mpl) * \
                                   ((Jl_ampl_ampl[1] - Jl_ampl_acl) / (Hmpl / n_mpl) + M_H2O * Sl_ampl[2])
    else: # n_mpl > 2
        dif_eq['ds_ampl_1 / dt'] = 1 / (rho_H2O_l(sv['T_ampl_1']) * epsilon_mpl) * \
                                    ((Jl_atl_ampl - Jl_ampl_ampl[1]) / (Hmpl / n_mpl) + M_H2O * Sl_ampl[1])
        for i in range(2, n_mpl):
            dif_eq[f'ds_ampl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_ampl_{i}']) * epsilon_mpl) * \
                                          ((Jl_ampl_ampl[i - 1] - Jl_ampl_ampl[i]) / (Hmpl / n_mpl) +
                                           M_H2O * Sl_ampl[i])
        dif_eq[f'ds_ampl_{n_mpl} / dt'] = 1 / (rho_H2O_l(sv[f'T_ampl_{n_mpl}']) * epsilon_mpl) * \
                                            ((Jl_ampl_ampl[n_mpl - 1] - Jl_ampl_acl) / (Hmpl / n_mpl) +
                                                M_H2O * Sl_ampl[n_mpl])
    #      Inside the ACL
    dif_eq['ds_acl / dt'] = 1 / (rho_H2O_l(sv['T_acl']) * epsilon_cl) * (Jl_ampl_acl / Hacl + M_H2O * Sl_acl)

    # At the cathode side
    #       Inside the CCL
    dif_eq['ds_ccl / dt'] = 1 / (rho_H2O_l(sv['T_ccl']) * epsilon_cl) * (- Jl_ccl_cmpl / Hccl + M_H2O * Sl_ccl)
    #       Inside the CMPL
    if n_mpl == 1:
        dif_eq['ds_cmpl_1 / dt'] = 1 / (rho_H2O_l(sv['T_cmpl_1']) * epsilon_mpl) * \
                                   ((Jl_ccl_cmpl - Jl_cmpl_ctl) / Hmpl + M_H2O * Sl_cmpl[1])
    elif n_mpl == 2:
        dif_eq['ds_cmpl_1 / dt'] = 1 / (rho_H2O_l(sv['T_cmpl_1']) * epsilon_mpl) * \
                                   ((Jl_ccl_cmpl - Jl_cmpl_cmpl[1]) / (Hmpl / n_mpl) + M_H2O * Sl_cmpl[1])
        dif_eq['ds_cmpl_2 / dt'] = 1 / (rho_H2O_l(sv['T_cmpl_2']) * epsilon_mpl) * \
                                   ((Jl_cmpl_cmpl[1] - Jl_cmpl_ctl) / (Hmpl / n_mpl) + M_H2O * Sl_cmpl[2])
    else: # n_mpl > 2
        dif_eq['ds_cmpl_1 / dt'] = 1 / (rho_H2O_l(sv['T_cmpl_1']) * epsilon_mpl) * \
                                    ((Jl_ccl_cmpl - Jl_cmpl_cmpl[1]) / (Hmpl / n_mpl) + M_H2O * Sl_cmpl[1])
        for i in range(2, n_mpl):
            dif_eq[f'ds_cmpl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_cmpl_{i}']) * epsilon_mpl) * \
                                          ((Jl_cmpl_cmpl[i - 1] - Jl_cmpl_cmpl[i]) / (Hmpl / n_mpl) +
                                           M_H2O * Sl_cmpl[i])
        dif_eq[f'ds_cmpl_{n_mpl} / dt'] = 1 / (rho_H2O_l(sv[f'T_cmpl_{n_mpl}']) * epsilon_mpl) * \
                                            ((Jl_cmpl_cmpl[n_mpl - 1] - Jl_cmpl_ctl) / (Hmpl / n_mpl) +
                                                M_H2O * Sl_cmpl[n_mpl])
    #       Inside the CTL
    if n_tl == 2:
        dif_eq['ds_ctl_1 / dt'] = 1 / (rho_H2O_l(sv['T_ctl_1']) * epsilon_ctl[1]) * \
                                   ((Jl_cmpl_ctl - Jl_ctl_ctl[1]) / (Htl / n_tl) + M_H2O * Sl_ctl[1])
        dif_eq['ds_ctl_2 / dt'] = 1 / (rho_H2O_l(sv['T_ctl_2']) * epsilon_ctl[2]) * \
                                   ((Jl_ctl_ctl[1] - Jl_ctl_cgdl) / (Htl / n_tl) + M_H2O * Sl_ctl[2])
    else: # n_tl > 2
        dif_eq['ds_ctl_1 / dt'] = 1 / (rho_H2O_l(sv['T_ctl_1']) * epsilon_ctl[1]) * \
                                    ((Jl_cmpl_ctl - Jl_ctl_ctl[1]) / (Htl / n_tl) + M_H2O * Sl_ctl[1])
        for i in range(2, n_tl):
            dif_eq[f'ds_ctl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_ctl_{i}']) * epsilon_ctl[i]) * \
                                          ((Jl_ctl_ctl[i - 1] - Jl_ctl_ctl[i]) / (Htl / n_tl) + M_H2O * Sl_ctl[i])
        dif_eq[f'ds_ctl_{n_tl} / dt'] = 1 / (rho_H2O_l(sv[f'T_ctl_{n_tl}']) * epsilon_ctl[n_tl]) * \
                                            ((Jl_ctl_ctl[n_tl - 1] - Jl_ctl_cgdl) / (Htl / n_tl) + M_H2O * Sl_ctl[n_tl])
    #       Inside the CGDL
    if n_gdl == 1:
        dif_eq['ds_cgdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_cgdl_1']) * epsilon_gdl) * \
                                 ((Jl_ctl_cgdl - Jl_cgdl_cgc) / Hgdl + M_H2O * Sl_cgdl[1])
    elif n_gdl == 2:
        dif_eq['ds_cgdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_cgdl_1']) * epsilon_gdl) * \
                                   ((Jl_ctl_cgdl - Jl_cgdl_cgdl[1]) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[1])
        dif_eq['ds_cgdl_2 / dt'] = 1 / (rho_H2O_l(sv['T_cgdl_2']) * epsilon_gdl) * \
                                   ((Jl_cgdl_cgdl[1] - Jl_cgdl_cgc) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[2])
    else:
        dif_eq['ds_cgdl_1 / dt'] = 1 / (rho_H2O_l(sv['T_cgdl_1']) * epsilon_gdl) * \
                                   ((Jl_ctl_cgdl - Jl_cgdl_cgdl[1]) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[1])
        for i in range(2, n_gdl):
            dif_eq[f'ds_cgdl_{i} / dt'] = 1 / (rho_H2O_l(sv[f'T_cgdl_{i}']) * epsilon_gdl) * \
                                          ((Jl_cgdl_cgdl[i - 1] - Jl_cgdl_cgdl[i]) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[i])
        dif_eq[f'ds_cgdl_{n_gdl} / dt'] = 1 / (rho_H2O_l(sv[f'T_cgdl_{n_gdl}']) * epsilon_gdl) * \
                                     ((Jl_cgdl_cgdl[n_gdl - 1] - Jl_cgdl_cgc) / (Hgdl / n_gdl) + M_H2O * Sl_cgdl[n_gdl])


def calculate_dyn_vapor_evolution_inside_MEA(dif_eq, sv, Hgdl, Hmpl, Hacl, Hccl, epsilon_gdl, epsilon_cl, epsilon_mpl,
                                             epsilon_atl, epsilon_ctl, Htl, n_gdl, n_tl, n_mpl, Jv_agc_agdl,
                                             Jv_agdl_agdl, Jv_agdl_atl, Jv_atl_atl, Jv_atl_ampl, Jv_ampl_ampl,
                                             Jv_ampl_acl, S_abs_acl, S_abs_ccl, Jv_ccl_cmpl, Jv_cmpl_cmpl, Jv_cmpl_ctl,
                                             Jv_ctl_ctl, Jv_ctl_cgdl, Jv_cgdl_cgdl, Jv_cgdl_cgc, Sv_agdl, Sv_atl,
                                             Sv_ampl, Sv_acl, Sv_ccl, Sv_cmpl, Sv_ctl, Sv_cgdl, **kwargs):
    """This function calculates the dynamic evolution of the vapor in the gas diffusion layers, the microporous layers,
    and the catalyst layers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    Hgdl : float
        Thickness of the gas diffusion layer (m).
    Hmpl : float
        Thickness of the microporous layer (m).
    Hacl : float
        Thickness of the anode catalyst layer (m).
    Hccl : float
        Thickness of the cathode catalyst layer (m).
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_cl : float
        Anode/cathode CL porosity.
    epsilon_mpl : float
        Anode/cathode MPL porosity.
    epsilon_atl : list
        Anode transition layer porosities.
    epsilon_ctl : list
        Cathode transition layer porosities.
    Htl : float
        Thickness of the anode transition layer (m).
    n_gdl : int
        Number of model nodes placed inside each GDL.
    n_tl : int
        Number of model nodes placed inside the anode transition layer.
    n_mpl : int
        Number of model nodes placed inside each MPL.
    Jv_agc_agdl : float
        Water vapor flow between the anode gas channel and the anode GDL (mol.m-2.s-1).
    Jv_agdl_agdl : float
        Water vapor flow between two nodes of the anode GDL (mol.m-2.s-1).
    Jv_agdl_atl : float
        Water vapor flow between the last node of the anode GDL and the first node of the anode transition layer (mol.m-2.s-1).
    Jv_atl_atl : list
        Water vapor flow between two nodes of the anode transition layer (mol.m-2.s-1).
    Jv_atl_ampl : float
        Water vapor flow between the anode transition layer and the anode MPL (mol.m-2.s-1).
    Jv_ampl_ampl : list
        Water vapor flow between two nodes of the anode MPL (mol.m-2.s-1).
    Jv_ampl_acl : list
        Water vapor flow between the anode microporous layer and the anode catalyst layer (mol.m-2.s-1).
    S_abs_acl: float
        Water vapor absorption in the anode CL (mol.m-3.s-1).
    S_abs_ccl: float
        Water vapor absorption in the cathode CL (mol.m-3.s-1).
    Jv_ccl_cmpl : list
        Water vapor flow between the cathode catalyst layer and the cathode microporous layer (mol.m-2.s-1).
    Jv_cmpl_cmpl: list
        Water vapor flow between two nodes of the cathode MPL (mol.m-2.s-
    Jv_cmpl_ctl : float
        Water vapor flow between the last node of the cathode microporous layer and the first node of the cathode transition layer (mol.m-2.s-1).
    Jv_ctl_ctl : list
        Water vapor flow between two nodes of the cathode transition layer (mol.m-2.s-1).
    Jv_ctl_cgdl : float
        Water vapor flow between the last node of the cathode transition layer and the first node of the cathode GDL (mol.m-2.s-1).
    Jv_cgdl_cgdl : list
        Water vapor flow between two nodes of the cathode GDL (mol.m-2.s-1).
    Jv_cgdl_cgc : float
        Water vapor flow between the cathode GDL and the cathode gas channel (mol.m-2.s-1).
    Sv_agdl : list
        Water vapor produced in the anode GDL (mol.m-3.s-1).
    Sv_atl : list
        Water vapor produced in the anode transition layer (mol.m-3.s-1).
    Sv_ampl : float
        Water vapor produced in the anode microporous layer (mol.m-3.s-1).
    Sv_acl : float
        Water vapor produced in the anode CL (mol.m-3.s-1).
    Sv_ccl : float
        Water vapor produced in the cathode CL (mol.m-3.s-1).
    Sv_cmpl : float
        Water vapor produced in the cathode microporous layer (mol.m-3.s-1).
    Sv_ctl : list
        Water vapor produced in the cathode transition layer (mol.m-3.s-1).
    Sv_cgdl : list
        Water vapor produced in the cathode GDL (mol.m-3.s-1).
    """

    # At the anode side
    #       Inside the AGDL
    n_gc = len(Jv_agdl_agdl)  # Number of gas channel nodes
    Jv_agc_agdl_avg = sum(Jv_agc_agdl[1:n_gc+1]) / n_gc # Average vapor flow from the gas channel to the GDL
    if n_gdl == 1:
        dif_eq['dC_v_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                   ((Jv_agc_agdl_avg - Jv_agdl_atl) / Hgdl + Sv_agdl[1])
    elif n_gdl == 2:
        dif_eq['dC_v_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                   ((Jv_agc_agdl_avg - Jv_agdl_agdl[1]) / (Hgdl / n_gdl) + Sv_agdl[1])
        dif_eq['dC_v_agdl_2 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_2'])) * \
                                   ((Jv_agdl_agdl[1] - Jv_agdl_atl) / (Hgdl / n_gdl) + Sv_agdl[2])
    else: # n_gdl > 2
        dif_eq['dC_v_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                     ((Jv_agc_agdl_avg - Jv_agdl_agdl[1]) / (Hgdl / n_gdl) + Sv_agdl[1])
        for i in range(2, n_gdl):
            dif_eq[f'dC_v_agdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{i}'])) * \
                                            ((Jv_agdl_agdl[i - 1] - Jv_agdl_agdl[i]) / (Hgdl / n_gdl) + Sv_agdl[i])
        dif_eq[f'dC_v_agdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{n_gdl}'])) * \
                                            ((Jv_agdl_agdl[n_gdl - 1] - Jv_agdl_atl) / (Hgdl / n_gdl) + Sv_agdl[n_gdl])
    #       Inside the ATL
    if n_tl == 2:
        dif_eq['dC_v_atl_1 / dt'] = 1 / (epsilon_atl[1] * (1 - sv['s_atl_1'])) * \
                                     ((Jv_agdl_atl - Jv_atl_atl[1]) / (Htl / n_tl) + Sv_atl[1])
        dif_eq['dC_v_atl_2 / dt'] = 1 / (epsilon_atl[2] * (1 - sv['s_atl_2'])) * \
                                     ((Jv_atl_atl[1] - Jv_atl_ampl) / (Htl / n_tl) + Sv_atl[2])
    else: # n_tl > 2
        dif_eq['dC_v_atl_1 / dt'] = 1 / (epsilon_atl[1] * (1 - sv['s_atl_1'])) * \
                                        ((Jv_agdl_atl - Jv_atl_atl[1]) / (Htl / n_tl) + Sv_atl[1])
        for i in range(2, n_tl):
            dif_eq[f'dC_v_atl_{i} / dt'] = 1 / (epsilon_atl[i] * (1 - sv[f's_atl_{i}'])) * \
                                            ((Jv_atl_atl[i - 1] - Jv_atl_atl[i]) / (Htl / n_tl) + Sv_atl[i])
        dif_eq[f'dC_v_atl_{n_tl} / dt'] = 1 / (epsilon_atl[n_tl] * (1 - sv[f's_atl_{n_tl}'])) * \
                                            ((Jv_atl_atl[n_tl - 1] - Jv_atl_ampl) / (Htl / n_tl) + Sv_atl[n_tl])
    #       Inside the AMPL
    if n_mpl == 1:
        dif_eq['dC_v_ampl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_1'])) * ((Jv_atl_ampl - Jv_ampl_acl) / Hmpl +
                                                                                 Sv_ampl[1])
    elif n_mpl == 2:
        dif_eq['dC_v_ampl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_1'])) * \
                                   ((Jv_atl_ampl - Jv_ampl_ampl[1]) / (Hmpl / n_mpl) + Sv_ampl[1])
        dif_eq['dC_v_ampl_2 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_2'])) * \
                                   ((Jv_ampl_ampl[1] - Jv_ampl_acl) / (Hmpl / n_mpl) + Sv_ampl[2])
    else: # n_mpl > 2
        dif_eq['dC_v_ampl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_1'])) * \
                                    ((Jv_atl_ampl - Jv_ampl_ampl[1]) / (Hmpl / n_mpl) + Sv_ampl[1])
        for i in range(2, n_mpl):
            dif_eq[f'dC_v_ampl_{i} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_ampl_{i}'])) * \
                                            ((Jv_ampl_ampl[i - 1] - Jv_ampl_ampl[i]) / (Hmpl / n_mpl) + Sv_ampl[i])
        dif_eq[f'dC_v_ampl_{n_mpl} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_ampl_{n_mpl}'])) * \
                                            ((Jv_ampl_ampl[n_mpl - 1] - Jv_ampl_acl) / (Hmpl / n_mpl) + Sv_ampl[n_mpl])
    #       Inside the ACL
    dif_eq['dC_v_acl / dt'] = 1 / (epsilon_cl * (1 - sv['s_acl'])) * (Jv_ampl_acl / Hacl - S_abs_acl + Sv_acl)

    # At the cathode side
    #       Inside the CCL
    dif_eq['dC_v_ccl / dt'] = 1 / (epsilon_cl * (1 - sv['s_ccl'])) * (- Jv_ccl_cmpl / Hccl - S_abs_ccl + Sv_ccl)
    #       Inside the CMPL
    if n_mpl == 1:
        dif_eq['dC_v_cmpl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_1'])) * ((Jv_ccl_cmpl - Jv_cmpl_ctl) / Hmpl +
                                                                                 Sv_cmpl[1])
    elif n_mpl == 2:
        dif_eq['dC_v_cmpl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_1'])) * \
                                   ((Jv_ccl_cmpl - Jv_cmpl_cmpl[1]) / (Hmpl / n_mpl) + Sv_cmpl[1])
        dif_eq['dC_v_cmpl_2 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_2'])) * \
                                   ((Jv_cmpl_cmpl[1] - Jv_cmpl_ctl) / (Hmpl / n_mpl) + Sv_cmpl[2])
    else: # n_mpl > 2
        dif_eq['dC_v_cmpl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_1'])) * \
                                    ((Jv_ccl_cmpl - Jv_cmpl_cmpl[1]) / (Hmpl / n_mpl) + Sv_cmpl[1])
        for i in range(2, n_mpl):
            dif_eq[f'dC_v_cmpl_{i} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_cmpl_{i}'])) * \
                                            ((Jv_cmpl_cmpl[i - 1] - Jv_cmpl_cmpl[i]) / (Hmpl / n_mpl) + Sv_cmpl[i])
        dif_eq[f'dC_v_cmpl_{n_mpl} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_cmpl_{n_mpl}'])) * \
                                            ((Jv_cmpl_cmpl[n_mpl - 1] - Jv_cmpl_ctl) / (Hmpl / n_mpl) + Sv_cmpl[n_mpl])
    #       Inside the CTL
    if n_tl == 2:
        dif_eq['dC_v_ctl_1 / dt'] = 1 / (epsilon_ctl[1] * (1 - sv['s_ctl_1'])) * \
                                     ((Jv_cmpl_ctl - Jv_ctl_ctl[1]) / (Htl / n_tl) + Sv_ctl[1])
        dif_eq['dC_v_ctl_2 / dt'] = 1 / (epsilon_ctl[2] * (1 - sv['s_ctl_2'])) * \
                                     ((Jv_ctl_ctl[1] - Jv_ctl_cgdl) / (Htl / n_tl) + Sv_ctl[2])
    else: # n_tl > 2
        dif_eq['dC_v_ctl_1 / dt'] = 1 / (epsilon_ctl[1] * (1 - sv['s_ctl_1'])) * \
                                        ((Jv_cmpl_ctl - Jv_ctl_ctl[1]) / (Htl / n_tl) + Sv_ctl[1])
        for i in range(2, n_tl):
            dif_eq[f'dC_v_ctl_{i} / dt'] = 1 / (epsilon_ctl[i] * (1 - sv[f's_ctl_{i}'])) * \
                                            ((Jv_ctl_ctl[i - 1] - Jv_ctl_ctl[i]) / (Htl / n_tl) + Sv_ctl[i])
        dif_eq[f'dC_v_ctl_{n_tl} / dt'] = 1 / (epsilon_ctl[n_tl] * (1 - sv[f's_ctl_{n_tl}'])) * \
                                            ((Jv_ctl_ctl[n_tl - 1] - Jv_ctl_cgdl) / (Htl / n_tl) + Sv_ctl[n_tl])
    #       Inside the CGDL
    Jv_cgdl_cgc_avg = sum(Jv_cgdl_cgc[1:n_gc + 1]) / n_gc  # Average vapor flow from the GDL to the gas channel
    if n_gdl == 1:
        dif_eq['dC_v_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * ((Jv_ctl_cgdl - Jv_cgdl_cgc_avg) / Hgdl +
                                                                                 Sv_cgdl[1])
    elif n_gdl == 2:
        dif_eq['dC_v_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * \
                                   ((Jv_ctl_cgdl - Jv_cgdl_cgdl[1]) / (Hgdl / n_gdl) + Sv_cgdl[1])
        dif_eq['dC_v_cgdl_2 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_2'])) * \
                                   ((Jv_cgdl_cgdl[1] - Jv_cgdl_cgc_avg) / (Hgdl / n_gdl) + Sv_cgdl[2])
    else: # n_gdl > 2
        dif_eq['dC_v_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * \
                                     ((Jv_ctl_cgdl - Jv_cgdl_cgdl[1]) / (Hgdl / n_gdl) + Sv_cgdl[1])
        for i in range(2, n_gdl):
            dif_eq[f'dC_v_cgdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{i}'])) * \
                                            ((Jv_cgdl_cgdl[i - 1] - Jv_cgdl_cgdl[i]) / (Hgdl / n_gdl) + Sv_cgdl[i])
        dif_eq[f'dC_v_cgdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{n_gdl}'])) * \
                                            ((Jv_cgdl_cgdl[n_gdl - 1] - Jv_cgdl_cgc_avg) / (Hgdl / n_gdl) + Sv_cgdl[n_gdl])


def calculate_dyn_H2_O2_N2_evolution_inside_MEA(dif_eq, sv, Hgdl, Hmpl, Hacl, Hccl, epsilon_gdl, epsilon_cl,
                                                epsilon_mpl, epsilon_atl, epsilon_ctl, Htl, n_gdl, n_tl, n_mpl,
                                                J_H2_agc_agdl, J_H2_agdl_agdl, J_H2_agdl_atl, J_H2_atl_atl,
                                                J_H2_atl_ampl, J_H2_ampl_ampl, J_H2_ampl_acl, J_O2_ccl_cmpl,
                                                J_O2_cmpl_cmpl, J_O2_cmpl_ctl, J_O2_ctl_ctl, J_O2_ctl_cgdl,
                                                J_O2_cgdl_cgdl, J_O2_cgdl_cgc, S_H2_acl, S_O2_ccl, **kwargs):
    """This function calculates the dynamic evolution of the hydrogen, oxygen and nitrogen in the gas diffusion layers,
    the microporous layers, and the catalyst layers.

    Parameters
    ----------
    dif_eq : dict
        Dictionary used for saving the differential equations.
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    Hgdl : float
        Thickness of the gas diffusion layer (m).
    Hmpl : float
        Thickness of the microporous layer (m).
    Hacl : float
        Thickness of the anode catalyst layer (m).
    Hccl : float
        Thickness of the cathode catalyst layer (m).
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_cl : float
        Anode/cathode CL porosity.
    epsilon_mpl : float
        Anode/cathode MPL porosity.
    epsilon_atl : list
        Anode transition layer porosities.
    epsilon_ctl : list
        Cathode transition layer porosities.
    Htl : float
        Thickness of the anode transition layer (m).
    n_gdl : int
        Number of model nodes placed inside each GDL.
    n_tl : int
        Number of model nodes placed inside the anode transition layer.
    n_mpl : int
        Number of model nodes placed inside each MPL.
    J_H2_agc_agdl : float
        Hydrogen flow between the anode gas channel and the anode GDL (mol.m-2.s-1).
    J_H2_agdl_agdl : list
        Hydrogen flow between two nodes of the anode GDL (mol.m-2.s-1).
    J_H2_agdl_atl : float
        Hydrogen flow between the last node of the anode GDL and the first node of the anode transition layer (mol.m-2.s-1).
    J_H2_atl_atl : list
        Hydrogen flow between two nodes of the anode transition layer (mol.m-2.s-1).
    J_H2_atl_ampl : float
        Hydrogen flow between the anode transition layer and the anode MPL (mol.m-2.s-1).
    J_H2_ampl_ampl : list
        Hydrogen flow between two nodes of the anode MPL (mol.m-2.s-1
    J_H2_ampl_acl : float
        Hydrogen flow between the anode microporous layer and the anode catalyst layer (mol.m-2.s-1).
    J_O2_ccl_cmpl : float
        Oxygen flow between the cathode catalyst layer and the cathode microporous layer (mol.m-2.s-1).
    J_O2_cmpl_cmpl: list
        Oxygen flow between two nodes of the cathode MPL (mol.m-2.s-1).
    J_O2_cmpl_ctl : float
        Oxygen flow between the last node of the cathode microporous layer and the first node of the cathode transition layer (mol.m-2.s-1).
    J_O2_ctl_ctl : list
        Oxygen flow between two nodes of the cathode transition layer (mol.m-2.s-
    J_O2_ctl_cgdl : float
        Oxygen flow between the last node of the cathode transition layer and the first node of the cathode GDL (mol.m-2.s-1).
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
    #      Inside the AGDL
    n_gc = len(J_H2_agc_agdl)  # Number of gas channel nodes
    J_H2_agc_agdl_avg = sum(J_H2_agc_agdl[1:n_gc + 1]) / n_gc  # Average H2 flow from the gas channel to the GDL
    if n_gdl == 1:
        dif_eq['dC_H2_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                      (J_H2_agc_agdl_avg - J_H2_agdl_atl) / Hgdl
    elif n_gdl == 2:
        dif_eq['dC_H2_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                      (J_H2_agc_agdl_avg - J_H2_agdl_agdl[1]) / (Hgdl / n_gdl)
        dif_eq['dC_H2_agdl_2 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_2'])) * \
                                      (J_H2_agdl_agdl[1] - J_H2_agdl_atl) / (Hgdl / n_gdl)
    else: # n_gdl > 2
        dif_eq['dC_H2_agdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_agdl_1'])) * \
                                      (J_H2_agc_agdl_avg - J_H2_agdl_agdl[1]) / (Hgdl / n_gdl)
        for i in range(2, n_gdl):
            dif_eq[f'dC_H2_agdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{i}'])) * \
                                             (J_H2_agdl_agdl[i - 1] - J_H2_agdl_agdl[i]) / (Hgdl / n_gdl)
        dif_eq[f'dC_H2_agdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_agdl_{n_gdl}'])) * \
                                             (J_H2_agdl_agdl[n_gdl - 1] - J_H2_agdl_atl) / (Hgdl / n_gdl)
    #      Inside the ATL
    if n_tl == 2:
        dif_eq['dC_H2_atl_1 / dt'] = 1 / (epsilon_atl[1] * (1 - sv['s_atl_1'])) * \
                                      (J_H2_agdl_atl - J_H2_atl_atl[1]) / (Htl / n_tl)
        dif_eq['dC_H2_atl_2 / dt'] = 1 / (epsilon_atl[2] * (1 - sv['s_atl_2'])) * \
                                      (J_H2_atl_atl[1] - J_H2_atl_ampl) / (Htl / n_tl)
    else: # n_tl > 2
        dif_eq['dC_H2_atl_1 / dt'] = 1 / (epsilon_atl[1] * (1 - sv['s_atl_1'])) * \
                                        (J_H2_agdl_atl - J_H2_atl_atl[1]) / (Htl / n_tl)
        for i in range(2, n_tl):
            dif_eq[f'dC_H2_atl_{i} / dt'] = 1 / (epsilon_atl[i] * (1 - sv[f's_atl_{i}'])) * \
                                             (J_H2_atl_atl[i - 1] - J_H2_atl_atl[i]) / (Htl / n_tl)
        dif_eq[f'dC_H2_atl_{n_tl} / dt'] = 1 / (epsilon_atl[n_tl] * (1 - sv[f's_atl_{n_tl}'])) * \
                                             (J_H2_atl_atl[n_tl - 1] - J_H2_atl_ampl) / (Htl / n_tl)
    #      Inside the AMPL
    if n_mpl == 1:
        dif_eq['dC_H2_ampl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_1'])) * (J_H2_atl_ampl - J_H2_ampl_acl) / Hmpl
    elif n_mpl == 2:
        dif_eq['dC_H2_ampl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_1'])) * \
                                      (J_H2_atl_ampl - J_H2_ampl_ampl[1]) / (Hmpl / n_mpl)
        dif_eq['dC_H2_ampl_2 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_2'])) * \
                                      (J_H2_ampl_ampl[1] - J_H2_ampl_acl) / (Hmpl / n_mpl)
    else: # n_mpl > 2
        dif_eq['dC_H2_ampl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_ampl_1'])) * \
                                      (J_H2_atl_ampl - J_H2_ampl_ampl[1]) / (Hmpl / n_mpl)
        for i in range(2, n_mpl):
            dif_eq[f'dC_H2_ampl_{i} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_ampl_{i}'])) * \
                                             (J_H2_ampl_ampl[i - 1] - J_H2_ampl_ampl[i]) / (Hmpl / n_mpl)
        dif_eq[f'dC_H2_ampl_{n_mpl} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_ampl_{n_mpl}'])) * \
                                             (J_H2_ampl_ampl[n_mpl - 1] - J_H2_ampl_acl) / (Hmpl / n_mpl)
    #      Inside the ACL
    dif_eq['dC_H2_acl / dt'] = 1 / (epsilon_cl * (1 - sv['s_acl'])) * (J_H2_ampl_acl / Hacl + S_H2_acl)

    # At the cathode side
    #      Inside the CCL
    dif_eq['dC_O2_ccl / dt'] = 1 / (epsilon_cl * (1 - sv['s_ccl'])) * (-J_O2_ccl_cmpl / Hccl + S_O2_ccl)
    #      Inside the CMPL
    if n_mpl == 1:
        dif_eq['dC_O2_cmpl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_1'])) * (J_O2_ccl_cmpl - J_O2_cmpl_ctl) / Hmpl
    elif n_mpl == 2:
        dif_eq['dC_O2_cmpl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_1'])) * \
                                      (J_O2_ccl_cmpl - J_O2_cmpl_cmpl[1]) / (Hmpl / n_mpl)
        dif_eq['dC_O2_cmpl_2 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_2'])) * \
                                      (J_O2_cmpl_cmpl[1] - J_O2_cmpl_ctl) / (Hmpl / n_mpl)
    else: # n_mpl > 2
        dif_eq['dC_O2_cmpl_1 / dt'] = 1 / (epsilon_mpl * (1 - sv['s_cmpl_1'])) * \
                                        (J_O2_ccl_cmpl - J_O2_cmpl_cmpl[1]) / (Hmpl / n_mpl)
        for i in range(2, n_mpl):
            dif_eq[f'dC_O2_cmpl_{i} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_cmpl_{i}'])) * \
                                             (J_O2_cmpl_cmpl[i - 1] - J_O2_cmpl_cmpl[i]) / (Hmpl / n_mpl)
        dif_eq[f'dC_O2_cmpl_{n_mpl} / dt'] = 1 / (epsilon_mpl * (1 - sv[f's_cmpl_{n_mpl}'])) * \
                                             (J_O2_cmpl_cmpl[n_mpl - 1] - J_O2_cmpl_ctl) / (Hmpl / n_mpl)
    #      Inside the CTL
    if n_tl == 2:
        dif_eq['dC_O2_ctl_1 / dt'] = 1 / (epsilon_ctl[1] * (1 - sv['s_ctl_1'])) * \
                                      (J_O2_cmpl_ctl - J_O2_ctl_ctl[1]) / (Htl / n_tl)
        dif_eq['dC_O2_ctl_2 / dt'] = 1 / (epsilon_ctl[2] * (1 - sv['s_ctl_2'])) * \
                                      (J_O2_ctl_ctl[1] - J_O2_ctl_cgdl) / (Htl / n_tl)
    else: # n_tl > 2
        dif_eq['dC_O2_ctl_1 / dt'] = 1 / (epsilon_ctl[1] * (1 - sv['s_ctl_1'])) * \
                                        (J_O2_cmpl_ctl - J_O2_ctl_ctl[1]) / (Htl / n_tl)
        for i in range(2, n_tl):
            dif_eq[f'dC_O2_ctl_{i} / dt'] = 1 / (epsilon_ctl[i] * (1 - sv[f's_ctl_{i}'])) * \
                                             (J_O2_ctl_ctl[i - 1] - J_O2_ctl_ctl[i]) / (Htl / n_tl)
        dif_eq[f'dC_O2_ctl_{n_tl} / dt'] = 1 / (epsilon_ctl[n_tl] * (1 - sv[f's_ctl_{n_tl}'])) * \
                                             (J_O2_ctl_ctl[n_tl - 1] - J_O2_ctl_cgdl) / (Htl / n_tl)
    #      Inside the CGDL
    J_O2_cgdl_cgc_avg = sum(J_O2_cgdl_cgc[1:n_gc + 1]) / n_gc  # Average O2 flow from the GDL to the gas channel
    if n_gdl == 1:
        dif_eq['dC_O2_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * (J_O2_ctl_cgdl - J_O2_cgdl_cgc_avg) / Hgdl
    elif n_gdl == 2:
        dif_eq['dC_O2_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * \
                                      (J_O2_ctl_cgdl - J_O2_cgdl_cgdl[1]) / (Hgdl / n_gdl)
        dif_eq['dC_O2_cgdl_2 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_2'])) * \
                                      (J_O2_cgdl_cgdl[1] - J_O2_cgdl_cgc_avg) / (Hgdl / n_gdl)
    else:
        dif_eq['dC_O2_cgdl_1 / dt'] = 1 / (epsilon_gdl * (1 - sv['s_cgdl_1'])) * \
                                      (J_O2_ctl_cgdl - J_O2_cgdl_cgdl[1]) / (Hgdl / n_gdl)
        for i in range(2, n_gdl):
            dif_eq[f'dC_O2_cgdl_{i} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{i}'])) * \
                                             (J_O2_cgdl_cgdl[i - 1] - J_O2_cgdl_cgdl[i]) / (Hgdl / n_gdl)
        dif_eq[f'dC_O2_cgdl_{n_gdl} / dt'] = 1 / (epsilon_gdl * (1 - sv[f's_cgdl_{n_gdl}'])) * \
                                             (J_O2_cgdl_cgdl[n_gdl - 1] - J_O2_cgdl_cgc_avg) / (Hgdl / n_gdl)


def calculate_dyn_temperature_evolution_inside_MEA(dif_eq, Hgdl, Hmpl, Hacl, Hccl, Hmem, Htl, n_gdl, n_tl,
                                                   n_mpl, rho_Cp0, Jt, Q_r, Q_sorp, Q_liq, Q_p, Q_e, **kwargs):
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
    Hmpl : float
        Thickness of the microporous layer, in m.
    Hacl : float
        Thickness of the anode catalyst layer, in m.
    Hccl : float
        Thickness of the cathode catalyst layer, in m.
    Hmem : float
        Thickness of the membrane, in m.
    Htl : float
        Thickness of the transition layer, in m.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    n_tl : int
        Number of model nodes placed inside the transition layer.
    n_mpl : int
        Number of model nodes placed inside each MPL.
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
    #       Inside the AGDL
    if n_gdl == 1:
        dif_eq['dT_agdl_1 / dt'] = (1 / rho_Cp0['agdl_1']) * ( (Jt['agc_agdl'] - Jt['agdl_atl']) / Hgdl +
                                                               Q_liq['agdl_1'] + Q_e['agdl_1'] )
    elif n_gdl == 2:
        dif_eq['dT_agdl_1 / dt'] = (1 / rho_Cp0['agdl_1']) * ( (Jt['agc_agdl'] - Jt['agdl_agdl_1']) / (Hgdl / n_gdl) +
                                                               Q_liq['agdl_1'] + Q_e['agdl_1'] )
        dif_eq['dT_agdl_2 / dt'] = (1 / rho_Cp0['agdl_2']) * ( (Jt['agdl_agdl_1'] - Jt['agdl_atl']) / (Hgdl / n_gdl) +
                                                               Q_liq['agdl_2'] + Q_e['agdl_2'] )
    else: # n_gdl > 2
        dif_eq['dT_agdl_1 / dt'] = (1 / rho_Cp0['agdl_1']) * ( (Jt['agc_agdl'] - Jt['agdl_agdl_1']) / (Hgdl / n_gdl)  +
                                                               Q_liq['agdl_1'] + Q_e['agdl_1'] )
        for i in range(2, n_gdl):
            dif_eq[f'dT_agdl_{i} / dt'] = (1 / rho_Cp0[f'agdl_{i}']) * \
                                          ( (Jt[f'agdl_agdl_{i - 1}'] - Jt[f'agdl_agdl_{i}']) / (Hgdl / n_gdl) +
                                            Q_liq[f'agdl_{i}'] + Q_e[f'agdl_{i}'] )
        dif_eq[f'dT_agdl_{n_gdl} / dt'] = (1 / rho_Cp0[f'agdl_{n_gdl}']) * \
                                          ( (Jt[f'agdl_agdl_{n_gdl - 1}'] - Jt['agdl_atl']) / (Hgdl / n_gdl) +
                                            Q_liq[f'agdl_{n_gdl}'] + Q_e[f'agdl_{n_gdl}'] )
    #       Inside the ATL
    if n_tl == 2:
        dif_eq['dT_atl_1 / dt'] = (1 / rho_Cp0['atl_1']) * ((Jt['agdl_atl'] - Jt['atl_atl_1']) / (Htl / n_tl) +
                                                              Q_liq['atl_1'] + Q_e['atl_1'])
        dif_eq['dT_atl_2 / dt'] = (1 / rho_Cp0['atl_2']) * ((Jt['atl_atl_1'] - Jt['atl_ampl']) / (Htl / n_tl) +
                                                              Q_liq['atl_2'] + Q_e['atl_2'])
    else: # n_tl > 2
        dif_eq['dT_atl_1 / dt'] = (1 / rho_Cp0['atl_1']) * ((Jt['agdl_atl'] - Jt['atl_atl_1']) / (Htl / n_tl) +
                                                                Q_liq['atl_1'] + Q_e['atl_1'])
        for i in range(2, n_tl):
            dif_eq[f'dT_atl_{i} / dt'] = (1 / rho_Cp0[f'atl_{i}']) * \
                                          ( (Jt[f'atl_atl_{i - 1}'] - Jt[f'atl_atl_{i}']) / (Htl / n_tl) +
                                            Q_liq[f'atl_{i}'] + Q_e[f'atl_{i}'] )
        dif_eq[f'dT_atl_{n_tl} / dt'] = (1 / rho_Cp0[f'atl_{n_tl}']) * \
                                          ( (Jt[f'atl_atl_{n_tl - 1}'] - Jt['atl_ampl']) / (Htl / n_tl) +
                                            Q_liq[f'atl_{n_tl}'] + Q_e[f'atl_{n_tl}'] )
    #      Inside the AMPL
    if n_mpl == 1:
        dif_eq['dT_ampl_1 / dt'] = (1 / rho_Cp0['ampl_1']) * ( (Jt['atl_ampl'] - Jt['ampl_acl']) / Hmpl +
                                                             Q_liq['ampl_1'] + Q_e['ampl_1'] )
    elif n_mpl == 2:
        dif_eq['dT_ampl_1 / dt'] = (1 / rho_Cp0['ampl_1']) * ( (Jt['atl_ampl'] - Jt['ampl_ampl_1']) / (Hmpl / n_mpl) +
                                                               Q_liq['ampl_1'] + Q_e['ampl_1'] )
        dif_eq['dT_ampl_2 / dt'] = (1 / rho_Cp0['ampl_2']) * ( (Jt['ampl_ampl_1'] - Jt['ampl_acl']) / (Hmpl / n_mpl) +
                                                               Q_liq['ampl_2'] + Q_e['ampl_2'] )
    else: # n_mpl > 2
        dif_eq['dT_ampl_1 / dt'] = (1 / rho_Cp0['ampl_1']) * ( (Jt['atl_ampl'] - Jt['ampl_ampl_1']) / (Hmpl / n_mpl) +
                                                                Q_liq['ampl_1'] + Q_e['ampl_1'] )
        for i in range(2, n_mpl):
            dif_eq[f'dT_ampl_{i} / dt'] = (1 / rho_Cp0[f'ampl_{i}']) * \
                                          ( (Jt[f'ampl_ampl_{i - 1}'] - Jt[f'ampl_ampl_{i}']) / (Hmpl / n_mpl) +
                                            Q_liq[f'ampl_{i}'] + Q_e[f'ampl_{i}'] )
        dif_eq[f'dT_ampl_{n_mpl} / dt'] = (1 / rho_Cp0[f'ampl_{n_mpl}']) * \
                                            ( (Jt[f'ampl_ampl_{n_mpl - 1}'] - Jt['ampl_acl']) / (Hmpl / n_mpl) +
                                              Q_liq[f'ampl_{n_mpl}'] + Q_e[f'ampl_{n_mpl}'] )
    #      Inside the ACL
    dif_eq['dT_acl / dt'] = (1 / rho_Cp0['acl']) * \
                            ( (Jt['ampl_acl'] - Jt['acl_mem']) / Hacl +
                             Q_r['acl'] + Q_sorp['acl'] + Q_liq['acl'] + Q_e['acl'] )

    # Inside the membrane
    dif_eq['dT_mem / dt'] = (1 / rho_Cp0['mem']) * \
                            ( (Jt['acl_mem'] - Jt['mem_ccl']) / Hmem + Q_p['mem'] )

    # At the cathode side
    #       Inside the CCL
    dif_eq['dT_ccl / dt'] = (1 / rho_Cp0['ccl']) * \
                            ( (Jt['mem_ccl'] - Jt['ccl_cmpl']) / Hccl +
                              Q_r['ccl'] + Q_sorp['ccl'] + Q_liq['ccl'] + Q_p['ccl'] + Q_e['ccl'] )
    #      Inside the CMPL
    if n_mpl == 1:
        dif_eq['dT_cmpl_1 / dt'] = (1 / rho_Cp0['cmpl_1']) * ((Jt['ccl_cmpl'] - Jt['cmpl_ctl']) / Hmpl +
                                                              Q_liq['cmpl_1'] + Q_e['cmpl_1'])
    elif n_mpl == 2:
        dif_eq['dT_cmpl_1 / dt'] = (1 / rho_Cp0['cmpl_1']) * \
                                   ( (Jt['ccl_cmpl'] - Jt['cmpl_cmpl_1']) / (Hmpl / n_mpl) +
                                     Q_liq['cmpl_1'] + Q_e['cmpl_1'] )
        dif_eq['dT_cmpl_2 / dt'] = (1 / rho_Cp0['cmpl_2']) * \
                                   ( (Jt['cmpl_cmpl_1'] - Jt['cmpl_ctl']) / (Hmpl / n_mpl) +
                                     Q_liq['cmpl_2'] + Q_e['cmpl_2'] )
    else: # n_mpl > 2
        dif_eq['dT_cmpl_1 / dt'] = (1 / rho_Cp0['cmpl_1']) * \
                                   ( (Jt['ccl_cmpl'] - Jt['cmpl_cmpl_1']) / (Hmpl / n_mpl) +
                                     Q_liq['cmpl_1'] + Q_e['cmpl_1'] )
        for i in range(2, n_mpl):
            dif_eq[f'dT_cmpl_{i} / dt'] = (1 / rho_Cp0[f'cmpl_{i}']) * \
                                          ( (Jt[f'cmpl_cmpl_{i - 1}'] - Jt[f'cmpl_cmpl_{i}']) / (Hmpl / n_mpl) +
                                            Q_liq[f'cmpl_{i}'] + Q_e[f'cmpl_{i}'] )
        dif_eq[f'dT_cmpl_{n_mpl} / dt'] = (1 / rho_Cp0[f'cmpl_{n_mpl}']) * \
                                          ( (Jt[f'cmpl_cmpl_{n_mpl - 1}'] - Jt['cmpl_ctl']) / (Hmpl / n_mpl) +
                                            Q_liq[f'cmpl_{n_mpl}'] + Q_e[f'cmpl_{n_mpl}'] )
    #       Inside the CTL
    if n_tl == 2:
        dif_eq['dT_ctl_1 / dt'] = (1 / rho_Cp0['ctl_1']) * ((Jt['cmpl_ctl'] - Jt['ctl_ctl_1']) / (Htl / n_tl) +
                                    Q_liq['ctl_1'] + Q_e['ctl_1'])
        dif_eq['dT_ctl_2 / dt'] = (1 / rho_Cp0['ctl_2']) * ((Jt['ctl_ctl_1'] - Jt['ctl_cgdl']) / (Htl / n_tl) +
                                    Q_liq['ctl_2'] + Q_e['ctl_2'])
    else: # n_tl > 2
        dif_eq['dT_ctl_1 / dt'] = (1 / rho_Cp0['ctl_1']) * ((Jt['cmpl_ctl'] - Jt['ctl_ctl_1']) / (Htl / n_tl) +
                                    Q_liq['ctl_1'] + Q_e['ctl_1'])
        for i in range(2, n_tl):
            dif_eq[f'dT_ctl_{i} / dt'] = (1 / rho_Cp0[f'ctl_{i}']) * \
                                          ( (Jt[f'ctl_ctl_{i - 1}'] - Jt[f'ctl_ctl_{i}']) / (Htl / n_tl) +
                                            Q_liq[f'ctl_{i}'] + Q_e[f'ctl_{i}'] )
        dif_eq[f'dT_ctl_{n_tl} / dt'] = (1 / rho_Cp0[f'ctl_{n_tl}']) * \
                                          ( (Jt[f'ctl_ctl_{n_tl - 1}'] - Jt['ctl_cgdl']) / (Htl / n_tl) +
                                            Q_liq[f'ctl_{n_tl}'] + Q_e[f'ctl_{n_tl}'] )
    #       Inside the CGDL
    if n_gdl == 1:
        dif_eq['dT_cgdl_1 / dt'] = (1 / rho_Cp0['cgdl_1']) * ( (Jt['ctl_cgdl'] - Jt['cgdl_cgc']) / Hgdl +
                                                               Q_liq['cgdl_1'] + Q_e['cgdl_1'] )
    elif n_gdl == 2:
        dif_eq['dT_cgdl_1 / dt'] = (1 / rho_Cp0['cgdl_1']) * \
                                   ( (Jt['ctl_cgdl'] - Jt['cgdl_cgdl_1']) / (Hgdl / n_gdl) +
                                     Q_liq['cgdl_1'] + Q_e['cgdl_1'] )
        dif_eq['dT_cgdl_2 / dt'] = (1 / rho_Cp0['cgdl_2']) * \
                                   ( (Jt['cgdl_cgdl_1'] - Jt['cgdl_cgc']) / (Hgdl / n_gdl) +
                                     Q_liq['cgdl_2'] + Q_e['cgdl_2'] )
    else: # n_gdl > 2
        dif_eq['dT_cgdl_1 / dt'] = (1 / rho_Cp0['cgdl_1']) * \
                                   ( (Jt['ctl_cgdl'] - Jt['cgdl_cgdl_1']) / (Hgdl / n_gdl) +
                                     Q_liq['cgdl_1'] + Q_e['cgdl_1'] )
        for i in range(2, n_gdl):
            dif_eq[f'dT_cgdl_{i} / dt'] = (1 / rho_Cp0[f'cgdl_{i}']) * \
                                          ( (Jt[f'cgdl_cgdl_{i - 1}'] - Jt[f'cgdl_cgdl_{i}']) / (Hgdl / n_gdl) +
                                            Q_liq[f'cgdl_{i}'] + Q_e[f'cgdl_{i}'] )
        dif_eq[f'dT_cgdl_{n_gdl} / dt'] = (1 / rho_Cp0[f'cgdl_{n_gdl}']) * \
                                          ( (Jt[f'cgdl_cgdl_{n_gdl - 1}'] - Jt['cgdl_cgc']) / (Hgdl / n_gdl) +
                                            Q_liq[f'cgdl_{n_gdl}'] + Q_e[f'cgdl_{n_gdl}'] )


def calculate_dyn_voltage_evolution(dif_eq, i_fc, C_O2_ccl, T_ccl, eta_c, Hccl, i0_d_c_ref, i0_h_c_ref, kappa_c, C_scl,
                                    i_n, f_drop, **kwargs):
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
    Hccl : float
        Thickness of the cathode catalyst layer (m).
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

    # dif_eq['deta_c / dt'] = 1 / (C_scl * Hccl) * ((i_fc + i_n) - i0_c_ref * (C_O2_ccl / C_O2ref) ** kappa_c *
    #                                              math.exp(f_drop * alpha_c * F / (R * T_ccl) * eta_c))

    dif_eq['deta_c / dt'] = 1 / (C_scl * Hccl) * ((i_fc + i_n) - (i0_d_c_ref ** f_drop * i0_h_c_ref ** (1 - f_drop)) *
                                                  (C_O2_ccl / C_O2ref) ** kappa_c * math.exp(alpha_c * F / (R * T_ccl) * eta_c))


