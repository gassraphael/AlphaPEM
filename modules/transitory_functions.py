# -*- coding: utf-8 -*-

"""This module contains transitory functions which all have a specific physical meaning for modeling the PEM fuel cell.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np
import math

# Importing constants' value
from configuration.settings import (M_eq, rho_mem, theta_c_gdl, theta_c_cl, M_H2, M_O2, M_N2, M_H2O, R, Kshape,
                                    epsilon_p, alpha_p, k_th_gdl, k_th_cl, k_th_mem, Cp_gdl, Cp_cl, Cp_mem, rho_gdl,
                                    rho_cl, sigma_e_gdl, sigma_e_cl)


# _________________________________________________Transitory functions_________________________________________________

def hmean(terms, weights=None):
    """
    Calculate the weighted harmonic mean of a list of terms with corresponding weights.
    It is more efficient to express this function in the code than calling hmean from scipy.stats.

    Parameters
    ----------
    terms (list of float):
        The terms to calculate the harmonic mean for.
    weights (list of float):
        The weights corresponding to each term. If None, uniform weights are assumed.

    Returns
    -------
    float:
        The weighted harmonic mean.
    """
    if weights is None:
        weights = [1] * len(terms)  # Assign equal weights if not provided

    if len(terms) != len(weights):
        raise ValueError("The length of terms and weights must be the same.")

    # Calculate the weighted harmonic mean
    weighted_sum = sum((w / t) for w, t in zip(weights, terms) if t != 0)
    total_weight = sum(weights)

    if weighted_sum == 0:
        return float('inf')  # Avoid division by zero

    return total_weight / weighted_sum


def average(terms, weights=None):
    """
    Calculate the weighted arithmetic mean of a list of terms with corresponding weights.
    It is more efficient to express this function in the code than calling average from numpy.

    Parameters
    ----------
    terms (list of float):
        The terms to calculate the average for.
    weights (list of float, optional):
        The weights corresponding to each term. If None, uniform weights are assumed.

    Returns
    -------
    float:
        The weighted arithmetic mean.
    """
    if weights is None:
        # If no weights are provided, use uniform weights
        weights = [1] * len(terms)

    if len(terms) != len(weights):
        raise ValueError("The length of terms and weights must be the same.")

    # Calculate the weighted arithmetic mean
    weighted_sum = sum(w * t for w, t in zip(weights, terms))
    total_weight = sum(weights)

    if total_weight == 0:
        return float('nan')  # Avoid division by zero

    return weighted_sum / total_weight


def rho_H2O_l(T):
    """This function calculates the water density, in kg.m-3, as a function of the temperature.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Water density in kg.m-3.
    """
    T_Celsius = T - 273.15
    return ((999.83952 + 16.945176 * T_Celsius - 7.9870401e-3 * T_Celsius ** 2 - 46.170461e-6 * T_Celsius ** 3
             + 105.56302e-9 * T_Celsius ** 4 - 280.54253e-12 * T_Celsius ** 5) / (1 + 16.879850e-3 * T_Celsius))


def nu_l(T):
    """This function calculates the liquid water kinematic viscosity, in m².s-1, as a function of the temperature.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Liquid water kinematic viscosity in m².s-1.
    """
    mu_l = 2.414 * 10 ** (-5 + 247.8 / (T - 140.0))  # Pa.s. It is the liquid water dynamic viscosity.
    return mu_l / rho_H2O_l(T)

def mu_gaz(component, T):
    """This function calculates the dynamic viscosity of different gases, in Pa.s, as a function of the temperature.

    Parameters
    ----------
    component : str
        Specifies the gas for which the dynamic viscosity is calculated.
        Must be either 'H2O_v' (vapor), 'H2' (hydrogen), 'O2' (hydrogen), or 'N2' (nitrogen).
    T : float
        Temperature in K.

    Returns
    -------
    float
        Dynamic viscosity of the selected gas in Pa.s.

    Notes
    -----
    Source : Carl L. Yaws - Manuel 2014 - Transport properties of chemicals and hydrocarbons
    (https://www.sciencedirect.com/book/9780323286589/transport-properties-of-chemicals-and-hydrocarbons)"""

    if component == 'H2O_v': # For T >= 150 K and T <= 1500 k.
        return (22.8211 + 1.7387e-1 * T + 3.2465e-4 * T ** 2 - 1.4334e-7 * T ** 3) * 10**-7
    elif component == 'H2': # For T >= 15 K and T <= 1500 K.
        return (1.7611 + 3.4165e-1 * T - 1.8368e-4 * T ** 2 + 5.1147e-8 * T ** 3) * 10**-7
    elif component == 'O2': # For T >= 54 K and T <= 1500 K.
        return (-4.9433 + 8.0673e-1 * T - 4.0416e-4 * T ** 2 + 1.0111e-7 * T ** 3) * 10**-7
    elif component == 'N2': # For T >= 63 K and T <= 1970 K.
        return (4.4656 + 6.3814e-1 * T - 2.6596e-4 * T ** 2 + 5.4113e-8 * T ** 3) * 10**-7
    else:
        raise ValueError("The element should be either 'H2O_v', 'H2', 'O2' or 'N2'.")


def Psat(T):
    """This function calculates the saturated partial pressure of vapor, in Pa, as a function of the temperature.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Saturated partial pressure of vapor in Pa.
    """
    Tcelsius = T - 273.15
    return 101325 * 10 ** (-2.1794 + 0.02953 * Tcelsius - 9.1837e-5 * Tcelsius ** 2 + 1.4454e-7 * Tcelsius ** 3)


def C_v_sat(T):
    """This function calculates the saturated vapor concentration for a perfect gas, in mol.m-3, as a function of the
    temperature.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Saturated vapor concentration for a perfect gas in mol.m-3.
    """
    return Psat(T) / (R * T)


def Dcap(element, s, T, epsilon, e, epsilon_c=None):
    """ This function calculates the capillary coefficient at the GDL or the CL and at the anode, in kg.m.s-1,
    considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the capillary coefficient is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.
    epsilon : float
        Porosity.
    e : float
        Capillary exponent.
    epsilon_c : float, optional
        Compression ratio of the GDL.
    """

    if element == 'gdl':
        if epsilon_c==None:
            raise ValueError("For the GDL, epsilon_c must be provided.")
        return sigma(T) * K0(element, epsilon, epsilon_c) / nu_l(T) * abs(math.cos(theta_c_gdl)) * \
               (epsilon / K0(element, epsilon, epsilon_c)) ** 0.5 * (s ** e + 1e-7) * (1.417 - 4.24 * s + 3.789 * s**2)

    elif element == 'cl':
        return sigma(T) * K0(element, epsilon) / nu_l(T) * abs(math.cos(theta_c_cl)) * \
               (epsilon / K0(element, epsilon)) ** 0.5 * (s ** e + 1e-7) * (1.417 - 4.24 * s + 3.789 * s**2)

    else:
        raise ValueError("The element should be either 'gdl' or 'cl'.")


def Da(P, T):
    """This function calculates the diffusion coefficient at the anode, in m².s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Diffusion coefficient at the anode in m².s-1.
    """
    return 1.644e-4 * (T / 333) ** 2.334 * (101325 / P)


def Dc(P, T):
    """This function calculates the diffusion coefficient at the cathode, in m².s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Diffusion coefficient at the cathode in m².s-1.
    """
    return 3.242e-5 * (T / 333) ** 2.334 * (101325 / P)


def Da_eff(element, s, T, P, epsilon, epsilon_c=None, tau=None):
    """This function calculates the effective diffusion coefficient at the GDL or the CL and at the anode, in m².s-1,
    considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the effective diffusion coefficient is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.
    tau : float, optional
        Pore structure coefficient in the CL. Required if element is 'cl'.

    Returns
    -------
    float
        Effective diffusion coefficient at the anode in m².s-1.
    """

    if element == 'gdl': # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        if epsilon_c==None:
            raise ValueError("For the GDL, epsilon_c must be provided.")
        # According to the GDL porosity, the GDL compression effect is different.
        if 0.55 <= epsilon < 0.67:
            beta2 = -1.59
        elif 0.67 <= epsilon < 0.8:
            beta2 = -0.90
        else:
            raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                             "epsilon_gdl should be between 0.55 and 0.8.")
        return epsilon * ((epsilon - epsilon_p) / (1 - epsilon_p)) ** alpha_p * math.exp(beta2 * epsilon_c) * (1 - s) ** 2 * Da(P, T)

    elif element == 'cl': # The effective diffusion coefficient at the CL using Bruggeman model.
        if tau==None:
            raise ValueError("For the CL, tau must be provided.")
        return epsilon ** tau * (1 - s) ** tau * Da(P, T)

    else:
        raise ValueError("The element should be either 'gdl' or 'cl'.")



def Dc_eff(element, s, T, P, epsilon, epsilon_c=None, tau=None):
    """This function calculates the effective diffusion coefficient at the GDL or the CL and at the cathode, in m².s-1,
    considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the effective diffusion coefficient is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.
    tau : float, optional
        Pore structure coefficient in the CL. Required if element is 'cl'.

    Returns
    -------
    float
        Effective diffusion coefficient at the cathode in m².s-1.
    """

    if element == 'gdl': # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        if epsilon_c==None:
            raise ValueError("For the GDL, epsilon_c must be provided.")
        # According to the GDL porosity, the GDL compression effect is different.
        if 0.55 <= epsilon < 0.67:
            beta2 = -1.59
        elif 0.67 <= epsilon < 0.8:
            beta2 = -0.90
        else:
            raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                             "epsilon_gdl should be between 0.55 and 0.8.")
        return epsilon * ((epsilon - epsilon_p) / (1 - epsilon_p)) ** alpha_p * math.exp(beta2 * epsilon_c) * (1 - s) ** 2 * Dc(P, T)

    elif element == 'cl': # The effective diffusion coefficient at the CL using Bruggeman model.
        if tau==None:
            raise ValueError("For the CL, tau must be provided.")
        return epsilon ** tau * (1 - s) ** tau * Dc(P, T)

    else:
        raise ValueError("The element should be either 'gdl' or 'cl'.")


def h_a(P, T, Wgc, Hgc):
    """This function calculates the effective convective-conductive mass transfer coefficient at the anode, in m.s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    Wgc : float
        Width of the gas channel in m.
    Hgc : float
        Thickness of the gas channel in m.

    Returns
    -------
    float
        Effective convective-conductive mass transfer coefficient at the anode in m.s-1.
    """
    Sh = 0.9247 * math.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.
    return Sh * Da(P, T) / Hgc


def h_c(P, T, Wgc, Hgc):
    """This function calculates the effective convective-conductive mass transfer coefficient at the cathode, in m.s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    Wgc : float
        Width of the gas channel in m.
    Hgc : float
        Thickness of the gas channel in m.

    Returns
    -------
    float
        Effective convective-conductive mass transfer coefficient at the cathode in m.s-1.
    """
    Sh = 0.9247 * math.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.
    return Sh * Dc(P, T) / Hgc


def lambda_eq(C_v, s, T):
    """This function calculates the equilibrium water content in the membrane. Hinatsu's expression modified with
    Bao's formulation has been selected.

    Parameters
    ----------
    C_v : float
        Water concentration variable in mol.m-3.
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Equilibrium water content in the membrane.
    """
    a_w = C_v / C_v_sat(T) + 2 * s  # water activity
    return 0.5 * (0.300 + 10.8 * a_w - 16.0 * a_w ** 2 + 14.1 * a_w ** 3) * (1 - math.tanh(100 * (a_w - 1))) \
        + 0.5 * (9.2 + 8.6 * (1 - math.exp(-Kshape * (a_w - 1)))) * (1 + math.tanh(100 * (a_w - 1)))


def D(lambdaa):
    """This function calculates the diffusion coefficient of water in the membrane, in m².s-1.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.

    Returns
    -------
    float
        Diffusion coefficient of water in the membrane in m².s-1.
    """
    return 4.1e-10 * (lambdaa / 25.0) ** 0.15 * (1.0 + math.tanh((lambdaa - 2.5) / 1.4))


def fv(lambdaa, T):
    """This function calculates the water volume fraction of the membrane.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.

    Returns
    -------
    float
        Water volume fraction of the membrane.
    """

    return (lambdaa * M_H2O / rho_H2O_l(T)) / (M_eq / rho_mem + lambdaa * M_H2O / rho_H2O_l(T))


def gamma_sorp(C_v, s, lambdaa, T, Hcl):
    """This function calculates the sorption rate of water in the membrane, in s-1.

    Parameters
    ----------
    C_v : float
        Water concentration variable in mol.m-3.
    s : float
        Liquid water saturation variable.
    lambdaa : float
        Water content in the membrane.
    T : float
        Temperature in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Sorption rate of water in the membrane in s-1.
    """

    if lambda_eq(C_v, s, T) >= lambdaa:  # absorption
        return (1.14e-5 * fv(lambdaa, T)) / Hcl * math.exp(2416 * (1 / 303 - 1 / T))
    else:  #                               desorption
        return (4.59e-5 * fv(lambdaa, T)) / Hcl * math.exp(2416 * (1 / 303 - 1 / T))


def Svl(s, C_v, Ctot, T, epsilon, gamma_cond, gamma_evap):
    """This function calculates the phase transfer rate of water condensation or evaporation, in mol.m-3.s-1.

    Parameters
    ----------
    s : float
        Liquid water saturation variable.
    C_v : float
        Water concentration variable in mol.m-3.
    Ctot : float
        Total gas concentration in mol.m-3.
    T : float
        Temperature in K.
    epsilon : float
        Porosity.
    gamma_cond : float
        Overall condensation rate constant for water in s-1.
    gamma_evap : float
        Overall evaporation rate constant for water in Pa-1.s-1.

    Returns
    -------
    float
        Phase transfer rate of water condensation or evaporation in mol.m-3.s-1.
    """

    if C_v > C_v_sat(T):  # condensation
        return gamma_cond * epsilon * (1 - s) * (C_v / Ctot) * (C_v - C_v_sat(T))
    else:  # evaporation
        return -gamma_evap * epsilon * s * rho_H2O_l(T) / M_H2O * R * T * (C_v_sat(T) - C_v)


def sigma(T):
    """This function calculates the water surface tension, in N.m-1, as a function of the temperature.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Water surface tension in N.m-1.
    """
    return 235.8e-3 * ((647.15 - T) / 647.15) ** 1.256 * (1 - 0.625 * (647.15 - T) / 647.15)


def K0(element, epsilon, epsilon_c=None):
    """This function calculates the intrinsic permeability, in m², considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the intrinsic permeability is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.

    Returns
    -------
    float
        Intrinsic permeability in m².
    """

    if element == 'gdl':
        if epsilon_c==None:
            raise ValueError("For the GDL, epsilon_c must be provided.")
        # According to the GDL porosity, the GDL compression effect is different.
        if 0.55 <= epsilon < 0.67:
            beta1 = -3.60
        elif 0.67 <= epsilon < 0.8:
            beta1 = -2.60
        else:
            raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                             "epsilon_gdl should be between 0.55 and 0.8.")
        return epsilon / (8 * math.log(epsilon) ** 2) * (epsilon - epsilon_p) ** (alpha_p + 2) * \
            4.6e-6 ** 2 / ((1 - epsilon_p) ** alpha_p * ((alpha_p + 1) * epsilon - epsilon_p) ** 2) * math.exp(beta1 * epsilon_c)

    elif element == 'cl':
        return epsilon / (8 * math.log(epsilon) ** 2) * (epsilon - epsilon_p) ** (alpha_p + 2) * \
            4.6e-6 ** 2 / ((1 - epsilon_p) ** alpha_p * ((alpha_p + 1) * epsilon - epsilon_p) ** 2)

    else:
        raise ValueError("The element should be either 'gdl' or 'cl'.")


def k_H2(lambdaa, T, kappa_co):
    """This function calculates the permeability coefficient of the membrane for hydrogen, in mol.m−1.s−1.Pa−1.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.
    T : float
        Temperature in K.
    kappa_co : float
        Crossover correction coefficient in mol.m-1.s-1.Pa-1.

    Returns
    -------
    float
        Permeability coefficient of the membrane for hydrogen in mol.m−1.s−1.Pa−1.
    """

    # Initialisation of the constants
    E_H2_v = 2.1e4  # J.mol-1. It is the activation energy of H2 for crossover in the under saturated membrane.
    E_H2_l = 1.8e4  # J.mol-1. It is the activation energy of H2 for crossover in the liquide-equilibrated membrane.
    Tref = 303.15  # K.

    # Calculation of the permeability coefficient of the membrane for hydrogen
    if lambdaa < 17.6:
        return kappa_co * (0.29 + 2.2 * fv(lambdaa, T)) * 1e-14 * math.exp(E_H2_v / R * (1 / Tref - 1 / T))
    else:
        return kappa_co * 1.8 * 1e-14 * math.exp(E_H2_l / R * (1 / Tref - 1 / T))


def k_O2(lambdaa, T, kappa_co):
    """This function calculates the permeability coefficient of the membrane for oxygen, in mol.m−1.s−1.Pa−1.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.
    T : float
        Temperature in K.
    kappa_co : float
        Crossover correction coefficient in mol.m-1.s-1.Pa-1.

    Returns
    -------
    float
        Permeability coefficient of the membrane for oxygen in mol.m−1.s−1.Pa−1.
    """

    # Initialisation of the constants
    E_O2_v = 2.2e4  # J.mol-1. It is the activation energy of oxygen for crossover in the under saturated membrane.
    E_O2_l = 2.0e4  # J.mol-1. It is the activation energy of oxygen for crossover in the liquide-equilibrated membrane.
    Tref = 303.15  # K.

    # Calculation of the permeability coefficient of the membrane for oxygen
    if lambdaa < 17.6:
        return kappa_co * (0.11 + 1.9 * fv(lambdaa, T)) * 1e-14 * math.exp(E_O2_v / R * (1 / Tref - 1 / T))
    else:
        return kappa_co * 1.2 * 1e-14 * math.exp(E_O2_l / R * (1 / Tref - 1 / T))


def sigma_p_eff(element, lambdaa, T, epsilon_mc=None):
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
    epsilon_mc : float
        Volume fraction of ionomer in the CCL.

    Returns
    -------
    float
        Proton conductivity in Ω-1.m-1.
    """
    if element == 'mem': # The proton conductivity at the membrane
        if lambdaa >= 1:
            return (0.5139 * lambdaa - 0.326) * math.exp(1268 * (1 / 303.15 - 1 / T))
        else:
            return 0.1879 * math.exp(1268 * (1 / 303.15 - 1 / T))
    elif element == 'ccl': # The effective proton conductivity at the cathode catalyst layer
        if epsilon_mc==None:
            raise ValueError("For the CCL, epsilon_mc must be provided.")
        if lambdaa >= 1:
            return epsilon_mc * (0.5139 * lambdaa - 0.326) * math.exp(1268 * (1 / 303.15 - 1 / T))
        else:
            return epsilon_mc * 0.1879 * math.exp(1268 * (1 / 303.15 - 1 / T))
    else:
        raise ValueError("The element should be either 'mem' or 'ccl'.")


def sigma_e_eff(element, epsilon, epsilon_c=None, epsilon_mc=None, tau=None):
    """This function calculates the effective electrical conductivity, in Ω-1.m-1, in either the GDL or the CL,
    considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the proton conductivity is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    epsilon : float
        Porosity.
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
    tau : float
        Pore structure coefficient in the CL.

    Returns
    -------
    float
        Effective electrical conductivity in Ω-1.m-1.
    """
    if element == 'gdl': # The effective electrical conductivity at the GDL
        if epsilon_c==None:
            raise ValueError("For the GDL, epsilon_c must be provided.")
        # According to the GDL porosity, the GDL compression effect is different.
        if 0.55 <= epsilon < 0.67:
            beta3 = 4.04
        elif 0.67 <= epsilon < 0.8:
            beta3 = 4.40
        else:
            raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                             "epsilon_gdl should be between 0.55 and 0.8.")
        return (1 - epsilon) * sigma_e_gdl * math.exp(beta3 * epsilon_c) # Using the volume fraction of conductive material.

    elif element == 'cl': # The effective electrical conductivity at the CL
        if epsilon_mc==None or tau==None:
            raise ValueError("For the CL, epsilon_mc and tau must be provided.")
        return (1 - epsilon - epsilon_mc  ) * sigma_e_cl # Using the volume fraction of conductive material.

    else:
        raise ValueError("The element should be either 'gdl' or 'cl'.")


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

    if abs(sum(x) - 1.0) > 1e-6:
        raise ValueError("The sum of the molar fractions should be 1.")

    n = len(k_th_g)
    A_W = np.zeros((n, n)) # Interaction coefficient from Wassiljewa equation.
    epsilon_TS = 0.85 # Value suggested by Tandon and Saxena in 1965.

    # Calculation of A_W using Maxon and Saxena suggestion.
    for i in range(n):
        for j in range(n):
            if i == j:
                A_W[i, j] = 1.0
            else:
                A_W[i, j] = (epsilon_TS * (1 + (mu_g[i] / mu_g[j])**0.5 * (M[j] / M[i]) ** 0.25) ** 2) /  \
                            (8 * (1 + M[i] / M[j]))**0.5

    # Calculation of the thermal conductivity of the gas mixture.
    k_th_gaz_mixture = 0.0
    for i in range(n):
        k_th_gaz_mixture += x[i] * k_th_g[i] / sum([x[j] * A_W[i, j] for j in range(n)])

    return k_th_gaz_mixture


def k_th_eff(element, T, C_v=None, s=None, lambdaa=None, C_H2=None, C_O2=None, C_N2=None, epsilon=None, epsilon_mc=None,
             epsilon_c=None):
    """This function calculates the effective thermal conductivity, in J.m-1.s-1.K-1, in either the GDL, the CL or the
    membrane. A weighted harmonic average is used for characterizing the conductivity of each material in a layer,
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
    epsilon_mc : float
        Volume fraction of ionomer in the CL.

    Returns
    -------
    float
        Effective thermal conductivity in J.m-1.s-1.K-1."""

    if element == 'agdl' or element == 'cgdl': # The effective thermal conductivity at the GDL
        if C_v==None or s==None or epsilon==None or epsilon_c==None:
            raise ValueError("For the GDL, C_v, s, epsilon and epsilon_c must be provided.")
        # According to the GDL porosity, the GDL compression effect is different.
        if 0.55 <= epsilon < 0.67:
            beta3 = 4.04
        elif 0.67 <= epsilon < 0.8:
            beta3 = 4.40
        else:
            raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                             "epsilon_gdl should be between 0.55 and 0.8.")
        if element == 'agdl': # The thermal conductivity of the gas mixture in the AGDL
            if C_H2 == None:
                raise ValueError("For the AGDL, C_H2 must be provided.")
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('H2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('H2', T)],
                                        [C_v / (C_v + C_H2), C_H2 / (C_v + C_H2)],
                                        [M_H2O, M_H2])
        else:                 # The thermal conductivity of the gas mixture in the CGDL
            if C_O2 == None or C_N2 == None:
                raise ValueError("For the CGDL, C_O2 and C_N2 must be provided.")
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('O2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('O2', T), mu_gaz('N2', T)],
                                        [C_v / (C_v + C_O2 + C_N2), C_O2 / (C_v + C_O2 + C_N2), C_N2 / (C_v + C_O2 + C_N2)],
                                        [M_H2O, M_O2, M_N2])
        return hmean([k_th_gdl * math.exp(beta3 * epsilon_c), k_th('H2O_l', T), k_th_gaz],
                     weights=[1 - epsilon, epsilon * s, epsilon * (1 - s)])

    elif element == 'acl' or element == 'ccl': # The effective thermal conductivity at the CL
        if C_v==None or lambdaa==None or s==None or epsilon==None or epsilon_mc==None:
            raise ValueError("For the CL, C_v, lambdaa, s, epsilon, epsilon_mc and tau must be provided.")
        k_th_eff_mem = hmean([k_th_mem, k_th('H2O_l', T)],
                             weights=[1 - fv(lambdaa, T), fv(lambdaa, T)]) # The effective thermal conductivity at the
        #                                                                    membrane
        if element == 'acl':  # The thermal conductivity of the gas mixture in the ACL
            if C_H2 == None:
                raise ValueError("For the ACL, C_H2 must be provided.")
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('H2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('H2', T)],
                                        [C_v / (C_v + C_H2), C_H2 / (C_v + C_H2)],
                                        [M_H2O, M_H2])
        else:  # The thermal conductivity of the gas mixture in the CCL
            if C_O2 == None or C_N2 == None:
                raise ValueError("For the CCL, C_O2 and C_N2 must be provided.")
            k_th_gaz = k_th_gaz_mixture([k_th('H2O_v', T), k_th('O2', T), k_th('N2', T)],
                                        [mu_gaz('H2O_v', T), mu_gaz('O2', T), mu_gaz('N2', T)],
                                        [C_v / (C_v + C_O2 + C_N2), C_O2 / (C_v + C_O2 + C_N2), C_N2 / (C_v + C_O2 + C_N2)],
                                        [M_H2O, M_O2, M_N2])
        return hmean([k_th_cl, k_th_eff_mem, k_th('H2O_l', T), k_th_gaz],
                     weights=[1 - epsilon - epsilon_mc, epsilon_mc, epsilon * s, epsilon * (1-s)])

    elif element == 'mem': # The effective thermal conductivity at the membrane
        if lambdaa==None:
            raise ValueError("For the membrane, lambdaa must be provided.")
        return hmean([k_th_mem, k_th('H2O_l', T)],
                     weights=[1 - fv(lambdaa, T), fv(lambdaa, T)])

    else:
        raise ValueError("The element should be either 'agdl', 'cgdl', 'acl', 'ccl' or 'mem'.")


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


def calculate_rho_Cp0(element, T, C_v=None, s=None, lambdaa=None, C_H2=None, C_O2=None, C_N2=None, epsilon=None, epsilon_mc=None):
    """This function calculates the volumetric heat capacity, in J.m-3.K-1, in either the GDL, the CL or the membrane.

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
    epsilon_mc : float
        Volume fraction of ionomer in the CL.

    Returns
    -------
    float
        Volumetric heat capacity in J.m-3.K-1."""

    if element == 'agdl' or element == 'cgdl':  # The volumetric heat capacity at the GDL
        if C_v is None or s is None or epsilon is None:
            raise ValueError("For the GDL, C_v, s and epsilon must be provided.")
        if element == 'agdl':  # The heat capacity of the gas mixture in the AGDL
            if C_H2 is None:
                raise ValueError("For the AGDL, C_H2 must be provided.")
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_H2 * C_H2 * Cp0('H2', T)],
                                    weights=[C_v / (C_v + C_H2), C_H2 / (C_v + C_H2)])
        else:  # The heat capacity of the gas mixture in the CGDL
            if C_O2 is None or C_N2 is None:
                raise ValueError("For the CGDL, C_O2 and C_N2 must be provided.")
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_O2 * C_O2 * Cp0('O2', T),
                                     M_N2 * C_N2 * Cp0('N2', T)],
                                    weights=[C_v / (C_v + C_O2 + C_N2), C_O2 / (C_v + C_O2 + C_N2), C_N2 / (C_v + C_O2 + C_N2)])
        return average([rho_gdl * Cp_gdl, rho_H2O_l(T) * Cp0('H2O_l', T), rho_Cp0_gaz],
                          weights=[1 - epsilon, epsilon * s, epsilon * (1 - s)])

    elif element == 'acl' or element == 'ccl':  # The volumetric heat capacity at the CL
        if C_v is None or lambdaa is None or s is None or epsilon is None or epsilon_mc is None:
            raise ValueError("For the CL, C_v, lambdaa, s, epsilon, and epsilon_mc must be provided.")
        if element == 'acl':  # The heat capacity of the gas mixture in the ACL
            if C_H2 is None:
                raise ValueError("For the ACL, C_H2 must be provided.")
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_H2 * C_H2 * Cp0('H2', T)],
                                    weights=[C_v / (C_v + C_H2), C_H2 / (C_v + C_H2)])
        else:  # The heat capacity of the gas mixture in the CCL
            if C_O2 is None or C_N2 is None:
                raise ValueError("For the CCL, C_O2 and C_N2 must be provided.")
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0('H2O_v', T), M_O2 * C_O2 * Cp0('O2', T),
                                     M_N2 * C_N2 * Cp0('N2', T)],
                                    weights=[C_v / (C_v + C_O2 + C_N2), C_O2 / (C_v + C_O2 + C_N2), C_N2 / (C_v + C_O2 + C_N2)])
        return average([rho_cl * Cp_cl, rho_mem * Cp_mem, rho_H2O_l(T) * Cp0('H2O_l', T), rho_Cp0_gaz],
                          weights=[1 - epsilon - epsilon_mc, epsilon_mc, epsilon * s, epsilon * (1 - s)])

    elif element == 'mem':  # The volumetric heat capacity at the membrane
        if lambdaa is None:
            raise ValueError("For the membrane, lambdaa must be provided.")
        return average([rho_mem * Cp_mem, rho_H2O_l(T) * Cp0('H2O_l', T)],
                          weights=[1 - fv(lambdaa, T), fv(lambdaa, T)])

    else:
        raise ValueError("The element should be either 'agdl', 'cgdl', 'acl', 'ccl' or 'mem'.")


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