# -*- coding: utf-8 -*-

"""This module contains physical functions which are used for modeling the PEM fuel cell."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value
from functools import lru_cache
from alphapem.utils.physics_constants import M_H2, M_O2, M_N2, M_H2O, R


# __________________________________________________Physical functions__________________________________________________


@lru_cache(maxsize=None) # Cache the results to optimize performance
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


@lru_cache(maxsize=None) # Cache the results to optimize performance
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


@lru_cache(maxsize=None) # Cache the results to optimize performance
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


@lru_cache(maxsize=None) # Cache the results to optimize performance
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


@lru_cache(maxsize=None) # Cache the results to optimize performance
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
