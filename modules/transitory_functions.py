# -*- coding: utf-8 -*-

"""This module contains transitory functions which all have a specific physical meaning for modeling the PEM fuel cell.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np

# Importing constants' value
from configuration.settings import M_eq, rho_mem, M_H2O, R


# _________________________________________________Transitory functions_________________________________________________

def rho_H2O(T):
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
    return ((999.83952 + 16.945176 * (T - 273.15) - 7.9870401e-3 * (T - 273.15) ** 2 - 46.170461e-6 * (T - 273.15) ** 3
             + 105.56302e-9 * (T - 273.15) ** 4 - 280.54253e-12 * (T - 273.15) ** 5) /
            (1 + 16.879850e-3 * (T - 273.15)))


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
    return mu_l / rho_H2O(T)


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
    return 101325 * 10 ** (-2.1794 + 0.02953 * (T - 273.15) - 9.1837e-5 * (T - 273.15) ** 2 +
                           1.4454e-7 * (T - 273.15) ** 3)


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


def Da_eff(s, epsilon, P, T, epsilon_c, epsilon_gdl):
    """This function calculates the effective diffusion coefficient at the anode, in m².s-1, considering GDL
    compression.
    Remark: it is considered here that the compression of the stack has a similar effect on the GDL and the CL, which
    may be wrong. This is why two porosities are considered in the parameters of this function: epsilon and epsilon_gdl.

    Parameters
    ----------
    s : float
        Liquid water saturation variable.
    epsilon : float
        Porosity.
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    epsilon_c : float
        Compression ratio of the GDL.
    epsilon_gdl : float
        Porosity of the GDL.

    Returns
    -------
    float
        Effective diffusion coefficient at the anode in m².s-1.
    """

    # According to the GDL porosity, the GDL compression effect is different.
    if 0.55 <= epsilon_gdl < 0.67:
        beta2 = -1.59
    elif 0.67 <= epsilon_gdl < 0.8:
        beta2 = -0.90
    else:
        raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                         "epsilon_gdl should be between 0.55 and 0.8.")

    return epsilon * ((epsilon - 0.11) / (1 - 0.11)) ** 0.785 * np.exp(beta2 * epsilon_c) * (1 - s) ** 2 * Da(P, T)


def Dc_eff(s, epsilon, P, T, epsilon_c, epsilon_gdl):
    """This function calculates the effective diffusion coefficient at the cathode, in m².s-1, considering GDL
    compression.
    Remark: it is considered here that the compression of the stack has a similar effect on the GDL and the CL, which
    may be wrong. This is why two porosities are considered in the parameters of this function: epsilon and epsilon_gdl.

    Parameters
    ----------
    s : float
        Liquid water saturation variable.
    epsilon : float
        Porosity.
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    epsilon_c : float
        Compression ratio of the GDL.
    epsilon_gdl : float
        Porosity of the GDL.

    Returns
    -------
    float
        Effective diffusion coefficient at the cathode in m².s-1.
    """

    # According to the GDL porosity, the GDL compression effect is different.
    if 0.55 <= epsilon_gdl < 0.67:
        beta2 = -1.59
    elif 0.67 <= epsilon_gdl < 0.8:
        beta2 = -0.90
    else:
        raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                         "epsilon_gdl should be between 0.55 and 0.8.")

    return epsilon * ((epsilon - 0.11) / (1 - 0.11)) ** 0.785 * np.exp(beta2 * epsilon_c) * (1 - s) ** 2 * Dc(P, T)


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
    Sh = 0.9247 * np.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.
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
    Sh = 0.9247 * np.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.
    return Sh * Dc(P, T) / Hgc


def lambda_eq(C_v, s, T, Kshape):
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
    Kshape : float
        Mathematical factor governing lambda_eq smoothing

    Returns
    -------
    float
        Equilibrium water content in the membrane.
    """
    a_w = C_v / C_v_sat(T) + 2 * s  # water activity
    return 0.5 * (0.300 + 10.8 * a_w - 16.0 * a_w ** 2 + 14.1 * a_w ** 3) * (1 - np.tanh(100 * (a_w - 1))) \
        + 0.5 * (9.2 + 8.6 * (1 - np.exp(-Kshape * (a_w - 1)))) * (1 + np.tanh(100 * (a_w - 1)))


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
    return 4.1e-10 * (lambdaa / 25.0) ** 0.15 * (1.0 + np.tanh((lambdaa - 2.5) / 1.4))


def gamma_sorp(C_v, s, lambdaa, T, Hcl, Kshape):
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
    Kshape : float
        Mathematical factor governing lambda_eq smoothing

    Returns
    -------
    float
        Sorption rate of water in the membrane in s-1.
    """

    fv = (lambdaa * M_H2O / rho_H2O(T)) / (M_eq / rho_mem + lambdaa * M_H2O / rho_H2O(T))  # water volume fraction of
    #                                                                                       the membrane
    if lambda_eq(C_v, s, T, Kshape) >= lambdaa:  # type_flow = absorption
        return (1.14e-5 * fv) / Hcl * np.exp(2416 * (1 / 303 - 1 / T))
    else:  # type_flow = desorption
        return (4.59e-5 * fv) / Hcl * np.exp(2416 * (1 / 303 - 1 / T))


def Svl(s, C_v, Ctot, epsilon, T, gamma_cond, gamma_evap):
    """This function calculates the phase transfer rate of water condensation or evaporation, in mol.m-3.s-1.

    Parameters
    ----------
    s : float
        Liquid water saturation variable.
    C_v : float
        Water concentration variable in mol.m-3.
    Ctot : float
        Total gas concentration in mol.m-3.
    epsilon : float
        Porosity.
    T : float
        Temperature in K.
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
        return -gamma_evap * epsilon * s * rho_H2O(T) / M_H2O * R * T * (C_v_sat(T) - C_v)


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


def K0(epsilon, epsilon_c, epsilon_gdl):
    """This function calculates the intrinsic permeability, in m², considering GDL compression.
    Remark: it is considered here that the compression of the stack has a similar effect on the GDL and the CL, which
    may be wrong. This is why two porosities are considered in the parameters of this function: epsilon and epsilon_gdl.

    Parameters
    ----------
    epsilon : float
        Porosity.
    epsilon_c : float
        Compression ratio of the GDL.
    epsilon_gdl : float
        Porosity of the GDL.

    Returns
    -------
    float
        Intrinsic permeability in m².
    """

    # According to the GDL porosity, the GDL compression effect is different.
    if 0.55 <= epsilon_gdl < 0.67:
        beta1 = -3.60
    elif 0.67 <= epsilon_gdl < 0.8:
        beta1 = -2.60
    else:
        raise ValueError("In order to calculate the effects of the GDL compression on its structure, "
                         "epsilon_gdl should be between 0.55 and 0.8.")

    return epsilon / (8 * np.log(epsilon) ** 2) * (epsilon - 0.11) ** (0.785 + 2) * \
        4.6e-6 ** 2 / ((1 - 0.11) ** 0.785 * ((0.785 + 1) * epsilon - 0.11) ** 2) * np.exp(beta1 * epsilon_c)


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
        fv = (lambdaa * M_H2O / rho_H2O(T)) / (M_eq / rho_mem + lambdaa * M_H2O / rho_H2O(T))
        return kappa_co * (0.29 + 2.2 * fv) * 1e-14 * np.exp(E_H2_v / R * (1 / Tref - 1 / T))
    else:
        return kappa_co * 1.8 * 1e-14 * np.exp(E_H2_l / R * (1 / Tref - 1 / T))


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
        fv = (lambdaa * M_H2O / rho_H2O(T)) / (M_eq / rho_mem + lambdaa * M_H2O / rho_H2O(T))
        return kappa_co * (0.11 + 1.9 * fv) * 1e-14 * np.exp(E_O2_v / R * (1 / Tref - 1 / T))
    else:
        return kappa_co * 1.2 * 1e-14 * np.exp(E_O2_l / R * (1 / Tref - 1 / T))
