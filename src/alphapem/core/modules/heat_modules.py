# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the heat transfer calculation.
"""
import math
from functools import lru_cache

from alphapem.core.modules.cell_voltage_modules import epsilon_mc, epsilon_cl
from alphapem.core.modules.flows_1D_MEA_modules import fv
from alphapem.utils.physics_constants import (sigma_e_gdl, sigma_e_mpl, sigma_e_cl, M_H2O, M_H2, M_N2, M_O2, k_th_gdl,
                                              k_th_mpl, k_th_mem, k_th_cl, rho_gdl, Cp_gdl, rho_mpl, Cp_mpl, rho_cl,
                                              Cp_cl, rho_mem, Cp_mem)
# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from alphapem.utils.maths_functions import hmean, average
from alphapem.utils.physics_functions import mu_gaz, rho_H2O_l


# _________________________________________________Heat transfer modules________________________________________________

def heat_transfer_int_values(sv, parameters):
    """This functions calculates intermediate values for the heat calculation.

    Parameters
    ----------
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    k_th_eff_agc_agdl : float
        Effective thermal diffusivity between the AGC and the first GDL layer (J.m-1.s-1.K-1).
    k_th_eff_agdl_agdl : list of floats
        List of effective thermal diffusivities between adjacent GDL layers on the anode side (J.m-1.s-1.K-1).
    k_th_eff_agdl_atl : float
        Effective thermal diffusivity between the last GDL layer and the anode transport layer (J.m-1.s-1.K-1).
    k_th_eff_atl_atl : list
        List of effective thermal diffusivities between adjacent transport layers on the anode side (J.m-1.s-1.K-1).
    k_th_eff_atl_ampl : float
        Effective thermal diffusivity between the anode transport layer and the first microporous layer (J.m-1.s-1.K-1).
    k_th_eff_ampl_ampl : list of floats
        List of effective thermal diffusivities between adjacent microporous layers on the anode side (J.m-1.s-1.K-1).
    k_th_eff_ampl_acl : float
        Effective thermal diffusivity between the anode microporous layer and the anode catalyst layer (J.m-1.s-1.K-1).
    k_th_eff_acl_mem : float
        Effective thermal diffusivity between the anode catalyst layer and the membrane (J.m-1.s-1.K-1).
    k_th_eff_mem_ccl : float
        Effective thermal diffusivity between the membrane and the cathode catalyst layer (J.m-1.s-1.K-1).
    k_th_eff_ccl_cmpl : float
        Effective thermal diffusivity between the cathode catalyst layer and the cathode microporous layer (J.m-1.s-1.K-1).
    k_th_eff_cmpl_cmpl : list of floats
        List of effective thermal diffusivities between adjacent microporous layers on the cathode side (J.m-1.s-1.K-1).
    k_th_eff_cmpl_ctl : float
        Effective thermal diffusivity between the last microporous layer and the cathode transport layer (J.m-1.s-1.K-1).
    k_th_eff_ctl_ctl : list of floats
        List of effective thermal diffusivities between adjacent transport layers on the cathode side (J.m-1.s-1.K-1).
    k_th_eff_ctl_cgdl : float
        Effective thermal diffusivity between the cathode transport layer and the first GDL layer (J.m-1.s-1.K-1).
    k_th_eff_cgdl_cgdl : list of floats
        List of effective thermal diffusivities between adjacent GDL layers on the cathode side (J.m-1.s-1.K-1).
    k_th_eff_cgdl_cgc : float
        Effective thermal diffusivity between the last GDL layer and the CGC (J.m-1.s-1.K-1).
    """

    # Extraction of the variables
    C_v_acl, C_v_ccl = sv['C_v_acl'], sv['C_v_ccl']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_acl, C_O2_ccl = sv['C_H2_acl'], sv['C_O2_ccl']
    C_N2_agc, C_N2_cgc = sv['C_N2_agc'], sv['C_N2_cgc']
    T_acl, T_mem, T_ccl = sv['T_acl'], sv['T_mem'], sv['T_ccl']

    # Extraction of the operating inputs and the parameters
    Hgdl, Hmpl, Hacl, Hccl = parameters['Hgdl'], parameters['Hmpl'], parameters['Hacl'], parameters['Hccl']
    Hmem, epsilon_gdl, epsilon_mpl = parameters['Hmem'], parameters['epsilon_gdl'], parameters['epsilon_mpl']
    epsilon_c, nb_gc, nb_gdl, nb_mpl = parameters['epsilon_c'], parameters['nb_gc'], parameters['nb_gdl'], parameters['nb_mpl']

    # Calculation of intermediate values
    Hgdl_node = Hgdl / nb_gdl
    Hmpl_node = Hmpl / nb_mpl

    # Weighted harmonic means of the effective thermal diffusivity
    k_th_eff_agc_agdl = k_th_eff('agdl', sv[f'T_agdl_{1}'], C_v=sv[f'C_v_agdl_{1}'], s=sv[f's_agdl_{1}'],
                                 C_H2=sv[f'C_H2_agdl_{1}'], C_N2=C_N2_agc, epsilon=epsilon_gdl, epsilon_c=epsilon_c)

    k_th_eff_agdl_agdl = [None] + [hmean([k_th_eff('agdl', sv[f'T_agdl_{i}'], C_v=sv[f'C_v_agdl_{i}'],
                                                   s=sv[f's_agdl_{i}'], C_H2=sv[f'C_H2_agdl_{i}'], C_N2=C_N2_agc,
                                                   epsilon=epsilon_gdl, epsilon_c=epsilon_c),
                                          k_th_eff('agdl', sv[f'T_agdl_{i + 1}'], C_v=sv[f'C_v_agdl_{i + 1}'],
                                                   s=sv[f's_agdl_{i + 1}'], C_H2=sv[f'C_H2_agdl_{i + 1}'], C_N2=C_N2_agc,
                                                   epsilon=epsilon_gdl, epsilon_c=epsilon_c)])
                                   for i in range(1, nb_gdl)]

    k_th_eff_agdl_ampl = hmean([k_th_eff('agdl', sv[f'T_agdl_{nb_gdl}'], C_v=sv[f'C_v_agdl_{nb_gdl}'],
                                               s=sv[f's_agdl_{nb_gdl}'], C_H2=sv[f'C_H2_agdl_{nb_gdl}'], C_N2=C_N2_agc,
                                               epsilon=epsilon_gdl, epsilon_c=epsilon_c),
                                     k_th_eff('ampl', sv[f'T_ampl_{1}'], C_v=sv[f'C_v_ampl_{1}'],
                                              s=sv[f's_ampl_{1}'], C_H2=sv[f'C_H2_ampl_{1}'], C_N2=C_N2_agc,
                                              epsilon=epsilon_mpl)],
                                   weights=[Hgdl_node / 2, Hmpl_node / 2])

    k_th_eff_ampl_ampl = [None] + [hmean([k_th_eff('ampl', sv[f'T_ampl_{i}'], C_v=sv[f'C_v_ampl_{i}'],
                                                   s=sv[f's_ampl_{i}'], C_H2=sv[f'C_H2_ampl_{i}'], C_N2=C_N2_agc,
                                                   epsilon=epsilon_mpl),
                                          k_th_eff('ampl', sv[f'T_ampl_{i + 1}'], C_v=sv[f'C_v_ampl_{i + 1}'],
                                                   s=sv[f's_ampl_{i + 1}'], C_H2=sv[f'C_H2_ampl_{i + 1}'], C_N2=C_N2_agc,
                                                   epsilon=epsilon_mpl)])
                                   for i in range(1, nb_mpl)]

    k_th_eff_ampl_acl = hmean([k_th_eff('ampl', sv[f'T_ampl_{nb_mpl}'], C_v=sv[f'C_v_ampl_{nb_mpl}'],
                                                  s=sv[f's_ampl_{nb_mpl}'], C_H2=sv[f'C_H2_ampl_{nb_mpl}'],
                                                 C_N2=C_N2_agc, epsilon=epsilon_mpl),
                                 k_th_eff('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl,
                                          C_H2=C_H2_acl, C_N2=C_N2_agc, Hcl = Hacl)],
                                weights=[Hmpl_node / 2, Hacl / 2])

    k_th_eff_acl_mem = hmean([k_th_eff('acl', T_acl, C_v=C_v_acl, s=s_acl, lambdaa=lambda_acl,
                                       C_H2=C_H2_acl, C_N2=C_N2_agc, Hcl = Hacl),
                                      k_th_eff('mem', T_mem, lambdaa=lambda_mem)],
                               weights=[Hacl / 2, Hmem / 2])

    k_th_eff_mem_ccl = hmean([k_th_eff('mem', T_mem, lambdaa=lambda_mem),
                                      k_th_eff('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl,
                                                C_O2=C_O2_ccl, C_N2=C_N2_cgc, Hcl = Hccl)],
                               weights=[Hmem / 2, Hccl / 2])

    k_th_eff_ccl_cmpl = hmean([k_th_eff('ccl', T_ccl, C_v=C_v_ccl, s=s_ccl, lambdaa=lambda_ccl, C_O2=C_O2_ccl,
                                                C_N2=C_N2_cgc, Hcl = Hccl),
                                       k_th_eff('cmpl', sv[f'T_cmpl_{1}'], C_v=sv[f'C_v_cmpl_{1}'],
                                                 s=sv[f's_cmpl_{1}'], C_O2=sv[f'C_O2_cmpl_{1}'], C_N2=C_N2_cgc,
                                                epsilon=epsilon_mpl)],
                                weights=[Hccl / 2, Hmpl_node / 2])

    k_th_eff_cmpl_cmpl = [None] + [hmean([k_th_eff('cmpl', sv[f'T_cmpl_{i}'], C_v=sv[f'C_v_cmpl_{i}'],
                                                   s=sv[f's_cmpl_{i}'], C_O2=sv[f'C_O2_cmpl_{i}'], C_N2=C_N2_cgc,
                                                   epsilon=epsilon_mpl),
                                            k_th_eff('cmpl', sv[f'T_cmpl_{i + 1}'], C_v=sv[f'C_v_cmpl_{i + 1}'],
                                                    s=sv[f's_cmpl_{i + 1}'], C_O2=sv[f'C_O2_cmpl_{i + 1}'], C_N2=C_N2_cgc,
                                                    epsilon=epsilon_mpl)])
                                      for i in range(1, nb_mpl)]

    k_th_eff_cmpl_cgdl = hmean([k_th_eff('cmpl', sv[f'T_cmpl_{nb_mpl}'], C_v=sv[f'C_v_cmpl_{nb_mpl}'],
                                               s=sv[f's_cmpl_{nb_mpl}'], C_O2=sv[f'C_O2_cmpl_{nb_mpl}'], C_N2=C_N2_cgc,
                                               epsilon=epsilon_mpl),
                                     k_th_eff('cgdl', sv['T_cgdl_1'], C_v=sv['C_v_cgdl_1'], s=sv['s_cgdl_1'],
                                              C_O2=sv[f'C_O2_cgdl_1'], C_N2=C_N2_cgc, epsilon=epsilon_gdl,
                                              epsilon_c=epsilon_c)],
                               weights=[Hmpl_node / 2, Hgdl_node / 2])

    k_th_eff_cgdl_cgdl = [None] + [hmean([k_th_eff('cgdl', sv[f'T_cgdl_{i}'], C_v=sv[f'C_v_cgdl_{i}'],
                                                   s=sv[f's_cgdl_{i}'], C_O2=sv[f'C_O2_cgdl_{i}'], C_N2=C_N2_cgc,
                                                   epsilon=epsilon_gdl, epsilon_c=epsilon_c),
                                          k_th_eff('cgdl', sv[f'T_cgdl_{i + 1}'], C_v=sv[f'C_v_cgdl_{i + 1}'],
                                                   s=sv[f's_cgdl_{i + 1}'], C_O2=sv[f'C_O2_cgdl_{i + 1}'], C_N2=C_N2_cgc,
                                                   epsilon=epsilon_gdl, epsilon_c=epsilon_c)])
                                   for i in range(1, nb_gdl)]

    k_th_eff_cgdl_cgc = k_th_eff('cgdl', sv[f'T_cgdl_{nb_gdl}'], C_v=sv[f'C_v_cgdl_{nb_gdl}'],
                                 s=sv[f's_cgdl_{nb_gdl}'], C_O2=sv[f'C_O2_cgdl_{nb_gdl}'], C_N2=C_N2_cgc,
                                 epsilon=epsilon_gdl, epsilon_c=epsilon_c)

    return (Hgdl_node, Hmpl_node, k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_ampl, k_th_eff_ampl_ampl,
            k_th_eff_ampl_acl, k_th_eff_acl_mem, k_th_eff_mem_ccl, k_th_eff_ccl_cmpl, k_th_eff_cmpl_cmpl,
            k_th_eff_cmpl_cgdl, k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def sigma_p_eff(element, lambdaa, T, Hcl=None):
    """This function calculates the effective proton conductivity, in Ω-1.m-1, in either the membrane or the CCL.

    Parameters
    ----------
    element : str
        Specifies the element for which the proton conductivity is calculated.
        Must be either 'mem' (membrane) or 'ccl' (cathode catalyst layer).
    lambdaa : float
        Water content in the membrane.
    T : float
        Temperature in K.
    Hcl : float, optional
        Thickness of the CL layer.

    Returns
    -------
    float
        Proton conductivity in Ω-1.m-1.
    """

    lambda_transition = 1

    if element == 'mem': # The proton conductivity at the membrane
        sigma_p_eff_low = 0.1879 * math.exp(1268 * (1 / 303.15 - 1 / T))
        sigma_p_eff_high = (0.5139 * lambdaa - 0.326) * math.exp(1268 * (1 / 303.15 - 1 / T))
    elif element == 'ccl': # The effective proton conductivity at the cathode catalyst layer
        sigma_p_eff_low = epsilon_mc(lambdaa, T, Hcl) * 0.1879 * math.exp(1268 * (1 / 303.15 - 1 / T))
        sigma_p_eff_high = epsilon_mc(lambdaa, T, Hcl) * (0.5139 * lambdaa - 0.326) * math.exp(1268 * (1 / 303.15 - 1 / T))
    else:
        raise ValueError("The element should be either 'mem' or 'ccl'.")

    # Transition function between low and high lambda
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_transition - lambdaa)))  # transition function

    return w * sigma_p_eff_low + (1 - w) * sigma_p_eff_high  # interpolation between sigma_p_eff value at low and high lambda.


@lru_cache(maxsize=None) # Cache the results to optimize performance
def sigma_e_eff(element, epsilon=None, epsilon_c=None, lambda_cl=None, T_cl=None, Hcl=None):
    """This function calculates the effective electrical conductivity, in Ω-1.m-1, in either the GDL, the MPL or the CL,
    considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the proton conductivity is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.
    lambda_cl : float, optional
        Water content in the CL.
    T_cl : float, optional
        Temperature inside the CL in K.
    Hcl : float, optional
        Thickness of the CL layer.

    Returns
    -------
    float
        Effective electrical conductivity in Ω-1.m-1.
    """
    if element == 'gdl': # The effective electrical conductivity at the GDL
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta3 = 4.04
        else:
            beta3 = 4.40
        return (1 - epsilon) * sigma_e_gdl * math.exp(beta3 * epsilon_c) # Using the volume fraction of conductive material.

    elif element == 'mpl': # The effective electrical conductivity at the MPL
        return (1 - epsilon) * sigma_e_mpl # Using the volume fraction of conductive material.

    elif element == 'cl': # The effective electrical conductivity at the CL
        return (1 - epsilon_cl(lambda_cl, T_cl, Hcl) - epsilon_mc(lambda_cl, T_cl, Hcl)) * sigma_e_cl # Using the volume fraction of conductive material.

    else:
        raise ValueError("The element should be either 'gdl', 'mpl' or 'cl'.")


@lru_cache(maxsize=None) # Cache the results to optimize performance
def k_th(component, T):
    """This function calculates the thermal conductivity of fluids, in J.m-1.s-1.K-1, as a function of the
    temperature.

    Parameters
    ----------
    component : str
        Specifies the gas for which the thermal conductivity is calculated.
        Must be either 'H2O_l' (liquid water), 'H2O_v' (vapor), 'H2' (hydrogen), 'O2' (hydrogen), or 'N2' (nitrogen).
    T : float
        Temperature in K.

    Returns
    -------
    float
        Thermal conductivity of the selected fluid in J.m-1.s-1.K-1.

    Notes
    -----
    Source : Carl L. Yaws - Manuel 2014 - Transport properties of chemicals and hydrocarbons
    (https://www.sciencedirect.com/book/9780323286589/transport-properties-of-chemicals-and-hydrocarbons)"""

    if component == 'H2O_l':  # For T >= 273.16 and T <= 633.15 K.
        return -0.2987 + 4.7054e-3 * T - 5.6209e-6 * T ** 2
    elif component == 'H2O_v': # For T >= 150 K and T <= 1500 k.
        return 5.6199e-3 + 1.5699e-5 * T + 1.0106e-7 * T ** 2 - 2.4282e-11 * T ** 3
    elif component == 'H2': # For T >= 14 K and T <= 1500 K.
        return 1.0979e-2 + 6.6411e-4 * T - 3.4378e-7 * T ** 2 + 9.7283e-11 * T ** 3
    elif component == 'O2': # For T >= 80 K and T <= 2000 K.
        return 1.5475e-4 + 9.4153e-5 * T - 2.7529e-8 * T ** 2 + 5.2069e-12 * T ** 3
    elif component == 'N2': # For T >= 63 K and T <= 1500 K.
        return -2.2678e-4 + 1.0275e-4 * T - 6.0151e-8 * T ** 2 + 2.2332e-11 * T ** 3
    else:
        raise ValueError("The element should be either 'H2O_l', 'H2O_v', 'H2', 'O2' or 'N2'.")


def k_th_gaz_mixture(k_th_g, mu_g, x, M):
    """This function calculates the thermal conductivity of a gas mixture, in J.m-1.s-1.K-1. The Lindsay–Bromley
    (Wassiljewa) method is used.

    Parameters
    ----------
    k_th_g : list
        Thermal conductivities of each pure gas component, in J.m-1.s-1.K-1, at the same temperature.
    mu_g : list
        Viscosity of each pure gas component, in Pa.s, at the same temperature.
    x : list
        Mole fractions of each gas component in the mixture (must sum to 1).
    M : list
        Molar masses of each gas component (in kg.mol-1).

    Returns
    -------
    lambda_mix : float
        Thermal conductivity of the gas mixture, in J.m-1.s-1.K-1.

    Notes
    -----
    Source : [wuMathematicalModelingTransient2009] and [polingPropertiesGasesLiquids2001]"""

    total_x = 0.0
    for xi in x:
        total_x += xi
    if abs(total_x - 1.0) > 1e-6:
        raise ValueError("The sum of the molar fractions should be 1.")

    n = len(k_th_g)
    epsilon_TS = 0.85 # Value suggested by Tandon and Saxena in 1965.

    # Calculation of A_W using Maxon and Saxena suggestion.
    A_W = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(1.0)
            else:
                val = (epsilon_TS * (1 + (mu_g[i] / mu_g[j]) ** 0.5 * (M[j] / M[i]) ** 0.25) ** 2) / \
                      (8 * (1 + M[i] / M[j])) ** 0.5
                row.append(val)
        A_W.append(row)

    # Calculation of the thermal conductivity of the gas mixture.
    k_th_gaz_mixture = 0.0
    for i in range(n):
        prod_x_A_w = 0.0
        for j in range(n):
            prod_x_A_w += x[j] * A_W[i][j]
        k_th_gaz_mixture += x[i] * k_th_g[i] / prod_x_A_w

    return k_th_gaz_mixture


@lru_cache(maxsize=None) # Cache the results to optimize performance
def k_th_eff(element, T, C_v=None, s=None, lambdaa=None, C_H2=None, C_O2=None, C_N2=None, epsilon=None, Hcl = None,
             epsilon_c=None):
    """This function calculates the effective thermal conductivity, in J.m-1.s-1.K-1, in either the GDL, the MPL, the CL
    or the membrane. A weighted harmonic average is used for characterizing the conductivity of each material in a layer,
    instead of a weighted arithmetic average. The physical meaning is that all the heat energy is forced to pass through
    all the material, as a series resistance network, instead of a parallel one
    [pharoahEffectiveTransportCoefficients2006].

    Parameters
    ----------
    element : str
        Specifies the element for which the proton conductivity is calculated.
        Must be either 'agdl' (anode gas diffusion layer), 'cgdl' (cathode gas diffusion layer),
        'acl' (anode catalyst layer), 'ccl' (cathode catalyst layer) or 'mem' (membrane).
    T : float
        Temperature in K.
    C_v : float
        Water concentration variable in mol.m-3.
    s : float
        Liquid water saturation variable.
    lambdaa : float
        Water content in the membrane.
    C_H2 : float
        Concentration of hydrogen in the AGDL or ACL.
    C_O2 : float
        Concentration of oxygen in the CGDL or CCL.
    C_N2 : float
        Concentration of nitrogen in the CGDL or CCL.
    epsilon : float
        Porosity.
    Hcl : float
        Thickness of the CL layer.
    epsilon_c : float
        Compression ratio of the GDL.

    Returns
    -------
    float
        Effective thermal conductivity in J.m-1.s-1.K-1."""

    if element in ('agdl', 'cgdl'): # The effective thermal conductivity at the GDL
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta3 = 4.04
        else:
            beta3 = 4.40
        if element == 'agdl': # The thermal conductivity of the gas mixture in the AGDL
            sum_C_v_C_H2_C_N2 = C_v + C_H2 + C_N2
            x_v, x_h2, x_n2 = C_v / sum_C_v_C_H2_C_N2, C_H2 / sum_C_v_C_H2_C_N2, C_N2 / sum_C_v_C_H2_C_N2
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('H2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('H2', T), mu_gaz('N2', T)],
                                        [x_v, x_h2, x_n2],
                                        [M_H2O, M_H2, M_N2])
        else:                 # The thermal conductivity of the gas mixture in the CGDL
            sum_C_v_C_O2_C_N2 = C_v + C_O2 + C_N2
            x_v, x_o2, x_n2 = C_v / sum_C_v_C_O2_C_N2, C_O2 / sum_C_v_C_O2_C_N2, C_N2 / sum_C_v_C_O2_C_N2
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('O2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('O2', T), mu_gaz('N2', T)],
                                        [x_v, x_o2, x_n2],
                                        [M_H2O, M_O2, M_N2])
        return hmean([k_th_gdl * math.exp(beta3 * epsilon_c), k_th('H2O_l', T), k_th_gaz],
                     weights=[1 - epsilon, epsilon * s, epsilon * (1 - s)])

    elif element in ('ampl', 'cmpl'): # The effective thermal conductivity at the GDL
        if element == 'ampl': # The thermal conductivity of the gas mixture in the AGDL
            sum_C_v_C_H2_C_N2 = C_v + C_H2 + C_N2
            x_v, x_h2, x_n2 = C_v / sum_C_v_C_H2_C_N2, C_H2 / sum_C_v_C_H2_C_N2, C_N2 / sum_C_v_C_H2_C_N2
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('H2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('H2', T), mu_gaz('N2', T)],
                                        [x_v, x_h2, x_n2],
                                        [M_H2O, M_H2, M_N2])
        else:                 # The thermal conductivity of the gas mixture in the CGDL
            sum_C_v_C_O2_C_N2 = C_v + C_O2 + C_N2
            x_v, x_o2, x_n2 = C_v / sum_C_v_C_O2_C_N2, C_O2 / sum_C_v_C_O2_C_N2, C_N2 / sum_C_v_C_O2_C_N2
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('O2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('O2', T), mu_gaz('N2', T)],
                                        [x_v, x_o2, x_n2],
                                        [M_H2O, M_O2, M_N2])
        return hmean([k_th_mpl, k_th('H2O_l', T), k_th_gaz],
                     weights=[1 - epsilon, epsilon * s, epsilon * (1 - s)])

    elif element in ('acl', 'ccl'): # The effective thermal conductivity at the CL
        fv_val = fv(lambdaa, T)
        epsilon_mc_val = epsilon_mc(lambdaa, T, Hcl)
        epsilon_cl_val = epsilon_cl(lambdaa, T, Hcl)
        k_th_eff_mem = hmean([k_th_mem, k_th('H2O_l', T)],
                             weights=[1 - fv_val, fv_val]) # The effective thermal conductivity at the membrane
        if element == 'acl':  # The thermal conductivity of the gas mixture in the ACL
            sum_C_v_C_H2_C_N2 = C_v + C_H2 + C_N2
            x_v, x_h2, x_n2 = C_v / sum_C_v_C_H2_C_N2, C_H2 / sum_C_v_C_H2_C_N2, C_N2 / sum_C_v_C_H2_C_N2
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('H2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('H2', T), mu_gaz('N2', T)],
                                        [x_v, x_h2, x_n2],
                                        [M_H2O, M_H2, M_N2])
        else:  # The thermal conductivity of the gas mixture in the CCL
            sum_C_v_C_O2_C_N2 = C_v + C_O2 + C_N2
            x_v, x_o2, x_n2 = C_v / sum_C_v_C_O2_C_N2, C_O2 / sum_C_v_C_O2_C_N2, C_N2 / sum_C_v_C_O2_C_N2
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('O2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('O2', T), mu_gaz('N2', T)],
                                        [x_v, x_o2, x_n2],
                                        [M_H2O, M_O2, M_N2])
        return hmean([k_th_cl, k_th_eff_mem, k_th('H2O_l', T), k_th_gaz],
                     weights=[1 - epsilon_cl_val - epsilon_mc_val, epsilon_mc_val, epsilon_cl_val * s, epsilon_cl_val * (1-s)])

    elif element == 'mem': # The effective thermal conductivity at the membrane
        fv_val = fv(lambdaa, T)
        return hmean([k_th_mem, k_th('H2O_l', T)],
                     weights=[1 - fv_val, fv_val])

    else:
        raise ValueError("The element should be either 'agdl', 'cgdl', 'ampl', 'cmpl', 'acl', 'ccl' or 'mem'.")


@lru_cache(maxsize=None) # Cache the results to optimize performance
def Cp0(component, T):
    """This function calculates the specific heat capacity of fluids, in J.kg-1.K-1, as a function of the
    temperature.

    Parameters
    ----------
    component : str
        Specifies the gas for which the specific heat capacity is calculated.
        Must be either 'H2O_l' (liquid water), 'H2O_v' (vapor), 'H2' (hydrogen), 'O2' (oxygen), or 'N2' (nitrogen).
    T : float
        Temperature in K.

    Returns
    -------
    float
        Specific heat capacity of the selected fluid in J.kg-1.K-1.

    Notes
    -----
    Source : Chase, M. W. (1998). NIST-JANAF Thermochemical Tables, 4th edition"""

    if component == 'H2O_l':  # For T >= 298 and T <= 500 K.
        return 1/M_H2O * (- 203.6060 + 1523.290 * (T/1000) - 3196.413 * (T/1000)**2 + 2474.455 * (T/1000)**3 +
                         3.855326 / (T/1000)**2)
    elif component == 'H2O_v':  # For T = 350 K. I failed to find a proper equation at the good range of temperature.
        return 1880
    elif component == 'H2':  # For T >= 298 K and T <= 1000 K.
        return 1/M_H2 * (33.066178 - 11.363417 * (T/1000) + 11.432816 * (T/1000)**2 - 2.772874 * (T/1000)**3 -
                         0.158558 / (T/1000)**2)
    elif component == 'O2':  # For T >= 100 K and T <= 700 K.
        return 1/M_O2 * (31.32234 - 20.23531 * (T/1000) + 57.86644 * (T/1000)**2 - 36.50624 * (T/1000)**3 -
                         0.007374 / (T/1000)**2)
    elif component == 'N2':  # For T >= 100 K and T <= 500 K.
        return 1/M_N2 * (28.98641 + 1.853978 * (T/1000) - 9.647459 * (T/1000)**2 + 16.63537 * (T/1000)**3 +
                         0.000117 / (T/1000)**2)
    else:
        raise ValueError("The element should be either 'H2O_l', 'H2O_v', 'H2', 'O2' or 'N2'.")


@lru_cache(maxsize=None) # Cache the results to optimize performance
def h0(component, T):
    """This function calculates the standard enthalpy of fluids, in J.mol-1, as a function of the temperature.
    The variation of the enthalpy of reaction with temperature is given by Kirchhoff's Law of Thermochemistry.

    Parameters
    ----------
    component : str
        Specifies the gas for which the specific heat capacity is calculated.
        Must be either 'H2O_l' (liquid water) or 'H2O_v' (vapor).
    T : float
        Temperature in K.

    Returns
    -------
    float
        Standard enthalpy of the selected fluid in J.mol-1.

    Notes
    -----
    Source : Chase, M. W. (1998). NIST-JANAF Thermochemical Tables, 4th edition"""

    if component == 'H2O_l':  # For T >= 298 and T <= 500 K.
        return (-285.83 - 203.6060 * (T/1000) + 1523.290 * (T/1000)**2 / 2 - 3196.413 * (T/1000)**3 / 3 +
                2474.455 * (T/1000)**4 / 4 - 3.855326 / (T/1000) - 256.5478 + 285.8304) * 1e3
    elif component == 'H2O_v':  # For T = 298.15 K. I failed to find a proper equation at the good range of temperature.
        return -241.83 * 1e3 + Cp0('H2O_v', T) * M_H2O * (T - 298.15)
    else:
        raise ValueError("The element should be either 'H2O_l' or 'H2O_v'")


@lru_cache(maxsize=None) # Cache the results to optimize performance
def calculate_rho_Cp0(element, T, C_v=None, s=None, lambdaa=None, C_H2=None, C_O2=None, C_N2=None, epsilon=None,
                      Hcl = None):
    """This function calculates the volumetric heat capacity, in J.m-3.K-1, in either the GDL, the MPL, the CL or
    the membrane.

    Parameters
    ----------
    element : str
        Specifies the element for which the volumetric heat capacity is calculated.
        Must be either 'agdl' (anode gas diffusion layer), 'cgdl' (cathode gas diffusion layer),
        'acl' (anode catalyst layer), 'ccl' (cathode catalyst layer) or 'mem' (membrane).
    T : float
        Temperature in K.
    C_v : float
        Water concentration variable in mol.m-3.
    s : float
        Liquid water saturation variable.
    lambdaa : float
        Water content in the membrane.
    C_H2 : float
        Concentration of hydrogen in the AGDL or ACL.
    C_O2 : float
        Concentration of oxygen in the CGDL or CCL.
    C_N2 : float
        Concentration of nitrogen in the CGDL or CCL.
    epsilon : float
        Porosity.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Volumetric heat capacity in J.m-3.K-1."""

    if element in ('agdl', 'cgdl', 'ampl', 'cmpl'):  # The volumetric heat capacity at the GDL
        if element in ('agdl', 'ampl'):  # In the anode
            sum_C_v_C_H2_C_N2 = C_v + C_H2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_H2 * C_H2 * Cp0('H2', T),
                                     M_N2 * C_N2 * Cp0('N2', T)],
                                    weights=[C_v / sum_C_v_C_H2_C_N2, C_H2 / sum_C_v_C_H2_C_N2, C_N2 / sum_C_v_C_H2_C_N2])
        else:  # In the cathode
            sum_C_v_C_O2_C_N2 = C_v + C_O2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_O2 * C_O2 * Cp0('O2', T),
                                     M_N2 * C_N2 * Cp0('N2', T)],
                                    weights=[C_v / sum_C_v_C_O2_C_N2, C_O2 / sum_C_v_C_O2_C_N2, C_N2 / sum_C_v_C_O2_C_N2])
        if element in ('agdl', 'cgdl'): # In the GDLs
            return average([rho_gdl * Cp_gdl, rho_H2O_l(T) * Cp0('H2O_l', T), rho_Cp0_gaz],
                           weights=[1 - epsilon, epsilon * s, epsilon * (1 - s)])
        else: # In the MPLs
            return average([rho_mpl * Cp_mpl, rho_H2O_l(T) * Cp0('H2O_l', T), rho_Cp0_gaz],
                           weights=[1 - epsilon, epsilon * s, epsilon * (1 - s)])

    elif element == 'acl' or element == 'ccl':  # The volumetric heat capacity at the CL
        epsilon_mc_val = epsilon_mc(lambdaa, T, Hcl)
        epsilon_cl_val = epsilon_cl(lambdaa, T, Hcl)
        if element == 'acl':  # The heat capacity of the gas mixture in the ACL
            sum_C_v_C_H2_C_N2 = C_v + C_H2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_H2 * C_H2 * Cp0('H2', T),
                                     M_N2 * C_N2 * Cp0('N2', T)],
                                    weights=[C_v / sum_C_v_C_H2_C_N2, C_H2 / sum_C_v_C_H2_C_N2, C_N2 / sum_C_v_C_H2_C_N2])
        else:  # The heat capacity of the gas mixture in the CCL
            sum_C_v_C_O2_C_N2 = C_v + C_O2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_O2 * C_O2 * Cp0('O2', T),
                                     M_N2 * C_N2 * Cp0('N2', T)],
                                    weights=[C_v / sum_C_v_C_O2_C_N2, C_O2 / sum_C_v_C_O2_C_N2, C_N2 / sum_C_v_C_O2_C_N2])
        return average([rho_cl * Cp_cl, rho_mem * Cp_mem, rho_H2O_l(T) * Cp0('H2O_l', T), rho_Cp0_gaz],
                       weights=[1 - epsilon_cl_val - epsilon_mc_val, epsilon_mc_val, epsilon_cl_val * s, epsilon_cl_val * (1 - s)])

    elif element == 'mem':  # The volumetric heat capacity at the membrane
        fv_val = fv(lambdaa, T)
        return average([rho_mem * Cp_mem, rho_H2O_l(T) * Cp0('H2O_l', T)],
                       weights=[1 - fv_val, fv_val])

    else:
        raise ValueError("The element should be either 'agdl', 'cgdl', 'ampl', 'cmpl', 'acl', 'ccl' or 'mem'.")


def delta_h_liq(T):
    """This function computes the molar enthalpy of liquefaction of water at a given temperature, in J.mol-1.
       It is calculated as the difference in molar enthalpy between liquid water (H2O_l) and water vapor (H2O_v).

        Parameters
        ----------
        T : float
            Temperature in K.

        Returns
        -------
        delta_h_liq : float
            Molar enthalpy of liquefaction in J.mol-1.

        Notes
        -----
        This value should be close to -42 000 J.mol-1 [vetterFreeOpenReference2019].
    """

    return h0('H2O_l', T) - h0('H2O_v', T)


def delta_h_abs(T):
    """This function computes the molar enthalpy of absorption of water at a given temperature, in J.mol-1.
    This reaction is exothermic.

        Parameters
        ----------
        T : float
            Temperature in K.

        Returns
        -------
        delta_h_sorp : float
            Molar enthalpy of absorption in the CL in J.mol-1.

        Notes
        -----
        For Nafion, the enthalpy of absorption is almost equal to that of liquefaction [vetterFreeOpenReference2019].
    """

    return delta_h_liq(T)
