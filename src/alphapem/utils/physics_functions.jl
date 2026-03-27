# -*- coding: utf-8 -*-

"""This module contains physical functions which are used for modeling the PEM fuel cell."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value
include(joinpath(@__DIR__, "physics_constants.jl"))


# __________________________________________________Physical functions__________________________________________________

"""
    rho_H2O_l(T)

This function calculates the water density, in kg.m-3, as a function of the temperature.

# Arguments
- `T`: Temperature in K.

# Returns
- `rho_H2O_l`: Water density in kg.m-3.
"""
function rho_H2O_l(T)
    T_Celsius = T - 273.15
    return ((999.83952 + 16.945176 * T_Celsius - 7.9870401e-3 * T_Celsius^2 - 46.170461e-6 * T_Celsius^3 +
             105.56302e-9 * T_Celsius^4 - 280.54253e-12 * T_Celsius^5) / (1 + 16.879850e-3 * T_Celsius))
end


"""
    nu_l(T)

This function calculates the liquid water kinematic viscosity, in m².s-1, as a function of the temperature.

# Arguments
- `T`: Temperature in K.

# Returns
- `nu_l`: Liquid water kinematic viscosity in m².s-1.
"""
function nu_l(T)
    mu_l = 2.414 * 10^(-5 + 247.8 / (T - 140.0))  # Pa.s. It is the liquid water dynamic viscosity.
    return mu_l / rho_H2O_l(T)
end


"""
    Psat(T)

This function calculates the saturated partial pressure of vapor, in Pa, as a function of the temperature.

# Arguments
- `T`: Temperature in K.

# Returns
- `Psat`: Saturated partial pressure of vapor in Pa.
"""
function Psat(T)
    Tcelsius = T - 273.15
    return 101325 * 10^(-2.1794 + 0.02953 * Tcelsius - 9.1837e-5 * Tcelsius^2 + 1.4454e-7 * Tcelsius^3)
end


"""
    C_v_sat(T)

This function calculates the saturated vapor concentration for a perfect gas, in mol.m-3, as a function of the
temperature.

# Arguments
- `T`: Temperature in K.

# Returns
- `C_v_sat`: Saturated vapor concentration for a perfect gas in mol.m-3.
"""
function C_v_sat(T)
    return Psat(T) / (R * T)
end


"""
    molar_mass(component)

This function returns the molar mass of a given gas component, in kg.mol-1.

# Arguments
- `component::String`: Specifies the gas for which the molar mass is returned.
  Must be either `"H2O_v"` (vapor), `"H2"` (hydrogen), `"O2"` (oxygen), or `"N2"` (nitrogen).

# Returns
- `M_gas`: Molar mass of the selected gas component in kg.mol-1.
"""
function molar_mass(component::String)
    if component == "H2O_v"
        return M_H2O
    elseif component == "H2"
        return M_H2
    elseif component == "O2"
        return M_O2
    elseif component == "N2"
        return M_N2
    else
        throw(ArgumentError("The element should be either 'H2O_v', 'H2', 'O2' or 'N2'."))
    end
end


"""
    mu_gaz(component, T)

This function calculates the dynamic viscosity of different gases, in Pa.s, as a function of the temperature.

# Arguments
- `component::String`: Specifies the gas for which the dynamic viscosity is calculated.
  Must be either `"H2O_v"` (vapor), `"H2"` (hydrogen), `"O2"` (oxygen), or `"N2"` (nitrogen).
- `T`: Temperature in K.

# Returns
- `mu_gaz`: Dynamic viscosity of the selected gas in Pa.s.

# Notes
Source : Carl L. Yaws - Manuel 2014 - Transport properties of chemicals and hydrocarbons
(https://www.sciencedirect.com/book/9780323286589/transport-properties-of-chemicals-and-hydrocarbons)
"""
function mu_gaz(component::String, T)

    if component == "H2O_v"  # For T >= 150 K and T <= 1500 k.
        return (22.8211 + 1.7387e-1 * T + 3.2465e-4 * T^2 - 1.4334e-7 * T^3) * 1e-7
    elseif component == "H2"  # For T >= 15 K and T <= 1500 K.
        return (1.7611 + 3.4165e-1 * T - 1.8368e-4 * T^2 + 5.1147e-8 * T^3) * 1e-7
    elseif component == "O2"  # For T >= 54 K and T <= 1500 K.
        return (-4.9433 + 8.0673e-1 * T - 4.0416e-4 * T^2 + 1.0111e-7 * T^3) * 1e-7
    elseif component == "N2"  # For T >= 63 K and T <= 1970 K.
        return (4.4656 + 6.3814e-1 * T - 2.6596e-4 * T^2 + 5.4113e-8 * T^3) * 1e-7
    else
        throw(ArgumentError("The element should be either 'H2O_v', 'H2', 'O2' or 'N2'."))
    end
end


"""
    mu_mixture_gases(components, x, T)

This function calculates the dynamic viscosity of a gas mixture, in Pa.s, as a function of the temperature.

# Arguments
- `components::Vector{String}`: List of gas components in the mixture. Each component must be either `"H2O_v"`
  (vapor), `"H2"` (hydrogen), `"O2"` (oxygen), or `"N2"` (nitrogen).
- `x::Vector`: List of mole fractions corresponding to each gas component in the mixture.
- `T`: Temperature in K.

# Returns
- `mu_mixture_gases`: Dynamic viscosity of the gas mixture in Pa.s.

# Notes
A simple mixture law is used here to calculate the dynamic viscosity of the gas mixture.
"""
function mu_mixture_gases(components::Vector, x::Vector, T)

    # Calculate the dynamic viscosities of each gas component in Pa.s.
    mu_values = [mu_gaz(comp, T) for comp in components]

    # Calculate the molar mass of the gas mixture in kg/mol.
    M_mix = 0.0
    @inbounds for j in eachindex(components)
        M_j = molar_mass(components[j])
        M_mix += M_j * x[j]
    end

    inv_mu_mix = 0.0
    @inbounds for j in eachindex(mu_values)
        mu_j = mu_values[j]
        M_j = molar_mass(components[j])
        c_j = M_j * x[j] / M_mix
        inv_mu_mix += c_j / mu_j
    end

    return 1 / inv_mu_mix
end
