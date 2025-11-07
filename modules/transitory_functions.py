# -*- coding: utf-8 -*-

"""This module contains transitory functions which all have a specific physical meaning for modeling the PEM fuel cell.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import math

# Importing constants' value
from functools import lru_cache
from configuration.settings import (M_eq, rho_mem, Dp_mpl, Dp_cl, theta_c_gdl, theta_c_mpl, theta_c_cl, gamma_cond,
                                    gamma_evap, M_H2, M_O2, M_N2, M_H2O, R, F, Kshape, epsilon_p, alpha_p,
                                    tau_mpl, tau_cl, r_s_gdl, r_s_mpl, r_s_cl, k_th_gdl, k_th_mpl, k_th_cl, k_th_mem,
                                    Cp_gdl, Cp_mpl, Cp_cl, Cp_mem, rho_gdl, rho_mpl, rho_cl, sigma_e_gdl, sigma_e_mpl,
                                    sigma_e_cl)


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

    n = len(terms)
    if weights is None:
        weights = [1] * n  # Assign equal weights if not provided

    if len(weights) != n:
        raise ValueError("The length of terms and weights must be the same.")

    # Calculate the weighted harmonic mean
    weighted_sum = 0
    total_weight = 0
    for w, t in zip(weights, terms):
        if t != 0:
            weighted_sum += w / t
        total_weight += w

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
    n = len(terms)
    if weights is None:
        total_weight = n
        weighted_sum = 0.0
        for t in terms:
            weighted_sum += t
    else:
        if n != len(weights):
            raise ValueError("The length of terms and weights must be the same.")
        total_weight = 0.0
        weighted_sum = 0.0
        for i in range(n):
            w = weights[i]
            total_weight += w
            weighted_sum += w * terms[i]

    if total_weight == 0:
        return float('nan')
    return weighted_sum / total_weight


def interpolate(terms, distances):
    """
    Fast inverse distance interpolation for exactly 2 points.

    Parameters
    ----------
    terms : list of float
        The values at each node ([y1, y2]).
    distances : list of float
        The distances from each node to the interpolation point ([d1, d2]).

    Returns
    -------
    float
        The interpolated value.
    """
    if len(terms) != 2 or len(distances) != 2:
        raise ValueError("This function only supports interpolation with 2 points.")
    y1, y2 = terms
    d1, d2 = distances
    if d1 == 0: return y1
    if d2 == 0: return y2
    return (d2 * y1 + d1 * y2) / (d1 + d2)



def d_dx(y_minus, y_plus, dx = None, dx_minus = None, dx_plus = None):
    """
    Computes the centered first derivative (second order) with different steps to the left and right.

    Parameters
    ----------
    y_minus : float
        Value at the left point (i-1).
    y_plus : float
        Value at the right point (i+1).
    dx : float
        Step between (i-1) and (i+1) when dx_minus = dx_plus.
    dx_minus : float
        Step between (i-1) and i.
    dx_plus : float
        Step between i and (i+1).

    Returns
    -------
    float
        Approximation of the first derivative at i.
    """

    # Case of uniform grid spacing
    if dx is None:
        if dx_minus is None or dx_plus is None:
            raise ValueError("Either dx or both dx_minus and dx_plus must be provided.")
    else:
        if dx == 0:
            raise ValueError("dx must be non-zero.")
        return (y_plus - y_minus) / (2.0 * dx)

    # Case of non-uniform grid spacing (dx is None and dx_minus and dx_plus are provided)
    if dx_minus <= 0 or dx_plus <= 0:
        raise ValueError("dx_minus and dx_plus must be positive non-zero values.")
    y_0 = interpolate([y_minus, y_plus], [dx_minus, dx_plus])
    return (y_plus * dx_minus**2 + y_0 * (dx_plus**2 - dx_minus**2) - y_minus * dx_plus**2)  / \
           (dx_minus * dx_plus * (dx_minus + dx_plus))


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


def mu_mixture_gases(components, x, T):
    """This function calculates the dynamic viscosity of a gas mixture, in Pa.s, as a function of the temperature.

    Parameters
    ----------
    components : list of str
        List of gas components in the mixture. Each component must be either 'H2O_v' (vapor), 'H2' (hydrogen),
        'O2' (oxygen), or 'N2' (nitrogen).
    x : list of float
        List of mole fractions corresponding to each gas component in the mixture.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Dynamic viscosity of the gas mixture in Pa.s.

    Notes
    -----
    A simple mixture law is used here to calculate the dynamic viscosity of the gas mixture.
    """

    # Calculate the dynamic viscosities of each gas component in Pa.s.
    mu_values = [mu_gaz(comp, T) for comp in components]

    # Calculate the molar mass of the gas mixture in kg/mol.
    M_mix = 0.0
    for j in range(len(components)):
        M_j = M_H2O if components[j] == 'H2O_v' else M_H2 if components[j] == 'H2' else M_O2 if components[j] == 'O2' \
              else M_N2 if components[j] == 'N2' else None
        M_mix += M_j * x[j]

    inv_mu_mix = 0.0
    for j, mu_j in enumerate(mu_values):
        M_j = M_H2O if components[j] == 'H2O_v' else M_H2 if components[j] == 'H2' else M_O2 if components[j] == 'O2' else M_N2
        c_j = M_j * x[j] / M_mix
        inv_mu_mix += c_j / mu_j

    return 1 / inv_mu_mix



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

    K0_value = K0(element, epsilon, epsilon_c)
    if element == 'gdl':
        theta_c_value = theta_c_gdl
    elif element == 'mpl':
        theta_c_value = theta_c_mpl
    elif element == 'cl':
        theta_c_value = theta_c_cl
    else:
        raise ValueError("The element should be either 'gdl', 'mpl' or 'cl'.")

    return sigma(T) * K0_value / nu_l(T) * abs(math.cos(theta_c_value)) * \
           (epsilon / K0_value) ** 0.5 * (s ** e + 1e-7) * (1.417 - 4.24 * s + 3.789 * s ** 2)


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


def Da_eff(element, s, T, P, epsilon, epsilon_c=None):
    """This function calculates the effective diffusion coefficient at the GDL, TL, MPL or the CL and at the anode,
     in m².s-1, considering GDL compression.

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

    Returns
    -------
    float
        Effective diffusion coefficient at the anode in m².s-1.
    """

    if element == 'gdl': # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta2 = -1.59
        else:
            beta2 = -0.90
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p)) ** alpha_p)
        return epsilon / tau_gdl * math.exp(beta2 * epsilon_c) * (1 - s) ** r_s_gdl * Da(P, T)

    elif element == 'mpl': # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s) ** r_s_mpl * Da(P, T)

    elif element == 'cl': # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_cl * (1 - s) ** r_s_cl * Da(P, T)

    else:
        raise ValueError("The element should be either 'gdl', 'tl', 'mpl' or 'cl'.")



def Dc_eff(element, s, T, P, epsilon, epsilon_c=None):
    """This function calculates the effective diffusion coefficient at the GDL, MPL, TL or the CL and at the cathode,
     in m².s-1, considering GDL compression.

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

    Returns
    -------
    float
        Effective diffusion coefficient at the cathode in m².s-1.
    """

    if element == 'gdl': # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta2 = -1.59
        else:
            beta2 = -0.90
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p)) ** alpha_p)
        return epsilon / tau_gdl * math.exp(beta2 * epsilon_c) * (1 - s) ** r_s_gdl * Dc(P, T)

    elif element == 'mpl': # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s) ** r_s_mpl * Dc(P, T)

    elif element == 'cl': # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_cl * (1 - s) ** r_s_cl * Dc(P, T)

    else:
        raise ValueError("The element should be either 'gdl', 'tl', 'mpl' or 'cl'.")


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


def lambda_v_eq(a_w):
    """This function calculates the equilibrium water content in the membrane from the vapor phase. Hinatsu's expression
     has been selected.

    Parameters
    ----------
    a_w: float
        Water activity.

    Returns
    -------
    float
        Equilibrium water content in the membrane from the vapor phase.
    """
    return 0.300 + 10.8 * a_w - 16.0 * a_w ** 2 + 14.1 * a_w ** 3


def lambda_l_eq(T):
    """This function calculates the equilibrium water content in the membrane from the liquid phase.
    Hinatsu's expression has been selected. It is valid for N-form membranes for 25 to 100 °C.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Equilibrium water content in the membrane from the liquid phase.
    """
    return 10.0 * 1.84e-2 * (T - 273.15) + 9.90e-4 * (T - 273.15)**2


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
    return 0.5 * lambda_v_eq(a_w)                                          * (1 - math.tanh(100 * (a_w - 1))) + \
           0.5 * (lambda_v_eq(1) + ((lambda_l_eq(T) - lambda_v_eq(1)) / 2) * (1 - math.exp(-Kshape * (a_w - 1)))) * \
                                                                             (1 + math.tanh(100 * (a_w - 1)))


def D_lambda(lambdaa):
    """This function calculates the diffusion coefficient of water in the bulk membrane, in m².s-1.

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


def D_lambda_eff(lambdaa, epsilon_mc):
    """This function calculates the effective diffusion coefficient of water in the catalyst layers, in m².s-1.

    Parameters
    ----------
    lambdaa : float
        Water content in the catalyst layer.

    Returns
    -------
    float
        Effective diffusion coefficient of water in the catalyst layer in m².s-1.
    """
    return epsilon_mc / tau_cl * D_lambda(lambdaa)


def D_EOD(i_fc):
    """This function calculates the electro-osmotic drag diffusion coefficient of water in the membrane, in mol.m-2.s-1.

    Parameters
    ----------
    i_fc : float
        Fuel cell current density in A.m-2.

    Returns
    -------
    float
        Electro-osmotic drag diffusion coefficient of water in the membrane in mol.m-2.s-1.
    """
    return 2.5 / 22 * i_fc / F


def D_EOD_eff(i_fc, epsilon_mc):
    """This function calculates the effective electro-osmotic drag diffusion coefficient of water in the catalyst layers,
    in mol.m-2.s-1.

    Parameters
    ----------
    i_fc : float
        Fuel cell current density in A.m-2.

    Returns
    -------
    float
        Effective electro-osmotic drag diffusion coefficient of water in the catalyst layer in mol.m-2.s-1.
    """
    return epsilon_mc / tau_cl * D_EOD(i_fc)


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

    fv_value = fv(lambdaa, T)
    gamma_abs = (1.14e-5 * fv_value) / Hcl * math.exp(2416 * (1 / 303 - 1 / T))
    gamma_des = (4.59e-5 * fv_value) / Hcl * math.exp(2416 * (1 / 303 - 1 / T))

    # Transition function between absorption and desorption
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_eq(C_v, s, T) - lambdaa))) # transition function

    return w * gamma_abs + (1 - w) * gamma_des # interpolation between absorption and desorption


def Svl(element, s, C_v, Ctot, T, epsilon):
    """This function calculates the phase transfer rate of water condensation or evaporation, in mol.m-3.s-1.
    It is positive for condensation and negative for evaporation.

    Parameters
    ----------
    element : str
        Specifies the element for which the phase transfer rate is calculated.
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

    Returns
    -------
    float
        Phase transfer rate of water condensation or evaporation in mol.m-3.s-1.
    """

    # Calculation of the total and partial pressures
    Ptot = Ctot * R * T # Total pressure
    P_v = C_v * R * T # Partial pressure of vapor

    # Determination of the diffusion coefficient at the anode or the cathode
    if element == 'anode':
        D_value = Da(Ptot, T)  # Diffusion coefficient at the anode
    else:  # element == 'cathode'
        D_value = Dc(Ptot, T)  # Diffusion coefficient at the cathode

    Svl_cond = gamma_cond * M_H2O / (R * T) * epsilon * (1 - s) * D_value * Ptot * math.log((Ptot - Psat(T)) / (Ptot - P_v))
    Svl_evap = gamma_evap * M_H2O / (R * T) * epsilon * s * D_value * Ptot * math.log((Ptot - Psat(T)) / (Ptot - P_v))

    # Transition function between condensation and evaporation
    K_transition = 3e-3 # This is a constant that defines the sharpness of the transition between two states.  The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (Psat(T) - P_v))) # transition function

    return w * Svl_evap + (1 - w) * Svl_cond # interpolation between condensation and evaporation


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


@lru_cache(maxsize=None) # Cache the results to optimize performance
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

    Sources
    -------
    1. Qin Chen 2020 - Two-dimensional multi-physics modeling of porous transport layer in polymer electrolyte membrane
    electrolyzer for water splitting - for the Blake-Kozeny equation.
    2. M.L. Stewart 2005 - A study of pore geometry effects on anisotropy in hydraulic permeability using the
    lattice-Boltzmann method - for the Blake-Kozeny equation.
    """

    if element == 'gdl':
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta1 = -3.60
        else:
            beta1 = -2.60
        return epsilon / (8 * math.log(epsilon) ** 2) * (epsilon - epsilon_p) ** (alpha_p + 2) * \
            4.6e-6 ** 2 / ((1 - epsilon_p) ** alpha_p * ((alpha_p + 1) * epsilon - epsilon_p) ** 2) * math.exp(beta1 * epsilon_c)

    elif element == 'mpl':
        return (Dp_mpl**2 / 150) * (epsilon**3 / ((1-epsilon)**2)) # Using the Blake-Kozeny equation

    elif element == 'cl':
        return (Dp_cl**2 / 150) * (epsilon**3 / ((1-epsilon)**2)) # Using the Blake-Kozeny equation

    else:
        raise ValueError("The element should be either 'gdl', 'mpl' or 'cl'.")


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
    E_H2_l = 1.8e4  # J.mol-1. It is the activation energy of H2 for crossover in the liquid-equilibrated membrane.
    Tref = 303.15  # K.

    # Calculation of the permeability coefficient of the membrane for hydrogen
    k_H2_d = kappa_co * (0.29 + 2.2 * fv(lambdaa, T)) * 1e-14 * math.exp(E_H2_v / R * (1 / Tref - 1 / T))
    k_H2_l = kappa_co * 1.8 * 1e-14 * math.exp(E_H2_l / R * (1 / Tref - 1 / T))

    # Transition function between under-saturated and liquid-saturated states
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_l_eq(T) - lambdaa)))  # transition function

    return w * k_H2_d + (1 - w) * k_H2_l  # interpolation between under-saturated and liquid-equilibrated H2 crossover


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
    E_O2_l = 2.0e4  # J.mol-1. It is the activation energy of oxygen for crossover in the liquid-equilibrated membrane.
    Tref = 303.15  # K.

    # Calculation of the permeability coefficient of the membrane for oxygen
    k_O2_v = kappa_co * (0.11 + 1.9 * fv(lambdaa, T)) * 1e-14 * math.exp(E_O2_v / R * (1 / Tref - 1 / T))
    k_O2_l = kappa_co * 1.2 * 1e-14 * math.exp(E_O2_l / R * (1 / Tref - 1 / T))

    # Transition function between under-saturated and liquid-saturated states
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_l_eq(T) - lambdaa)))  # transition function

    return w * k_O2_v + (1 - w) * k_O2_l  # interpolation between under-saturated and liquid-equilibrated O2 crossover


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

    lambda_transition = 1

    if element == 'mem': # The proton conductivity at the membrane
        sigma_p_eff_low = 0.1879 * math.exp(1268 * (1 / 303.15 - 1 / T))
        sigma_p_eff_high = (0.5139 * lambdaa - 0.326) * math.exp(1268 * (1 / 303.15 - 1 / T))
    elif element == 'ccl': # The effective proton conductivity at the cathode catalyst layer
        sigma_p_eff_low = epsilon_mc * 0.1879 * math.exp(1268 * (1 / 303.15 - 1 / T))
        sigma_p_eff_high = epsilon_mc * (0.5139 * lambdaa - 0.326) * math.exp(1268 * (1 / 303.15 - 1 / T))
    else:
        raise ValueError("The element should be either 'mem' or 'ccl'.")

    # Transition function between low and high lambda
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_transition - lambdaa)))  # transition function

    return w * sigma_p_eff_low + (1 - w) * sigma_p_eff_high  # interpolation between sigma_p_eff value at low and high lambda.


@lru_cache(maxsize=None) # Cache the results to optimize performance
def sigma_e_eff(element, epsilon, epsilon_c=None, epsilon_mc=None):
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
    epsilon_mc : float, optional
        Volume fraction of ionomer in the CL.

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
        return (1 - epsilon - epsilon_mc) * sigma_e_cl # Using the volume fraction of conductive material.

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


def k_th_eff(element, T, C_v=None, s=None, lambdaa=None, C_H2=None, C_O2=None, C_N2=None, epsilon=None, epsilon_mc=None,
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
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
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
                     weights=[1 - epsilon - epsilon_mc, epsilon_mc, epsilon * s, epsilon * (1-s)])

    elif element == 'mem': # The effective thermal conductivity at the membrane
        fv_val = fv(lambdaa, T)
        return hmean([k_th_mem, k_th('H2O_l', T)],
                     weights=[1 - fv_val, fv_val])

    else:
        raise ValueError("The element should be either 'agdl', 'cgdl', 'ampl', 'cmpl', 'acl', 'ccl' or 'mem'.")


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


def calculate_rho_Cp0(element, T, C_v=None, s=None, lambdaa=None, C_H2=None, C_O2=None, C_N2=None, epsilon=None,
                      epsilon_mc=None):
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
    epsilon_mc : float
        Volume fraction of ionomer in the CL.

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
                          weights=[1 - epsilon - epsilon_mc, epsilon_mc, epsilon * s, epsilon * (1 - s)])

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