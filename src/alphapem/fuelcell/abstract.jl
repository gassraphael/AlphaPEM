"""
Abstract type for all fuel cell models.
"""

abstract type AbstractFuelCell end

# Abstract interface: all concrete fuel cell models must implement these methods

"""
    physical_params(fc::AbstractFuelCell)::PhysicalParams
Return the physical parameters of the fuel cell.
"""
function physical_params(fc::AbstractFuelCell)::PhysicalParams
    throw(MethodError(physical_params, (fc,)))
end

"""
    operating_conditions(fc::AbstractFuelCell, type_fuel_cell::String)::OperatingConditions
Return the operating parameters of the fuel cell.
"""
function operating_conditions(fc::AbstractFuelCell, type_fuel_cell::String)::OperatingConditions
    throw(MethodError(operating_conditions, (fc, type_fuel_cell)))
end

"""
    pola_exp_data(fc::AbstractFuelCell, type_fuel_cell::Symbol, voltage_zone::Symbol)::PolaExperimentalData
Return the polarization experimental data of the fuel cell.
"""
function pola_exp_data(fc::AbstractFuelCell, type_fuel_cell::Symbol, voltage_zone::Symbol)::PolaExperimentalData
    throw(MethodError(experimental_values, (fc, type_fuel_cell, voltage_zone)))
end

"""
    numerical_params(fc::AbstractFuelCell)::NumericalParams
Return the numerical parameters of the fuel cell.
"""
function numerical_params(fc::AbstractFuelCell)::NumericalParams
    throw(MethodError(numerical_params, (fc,)))
end