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
    experimental_values(fc::AbstractFuelCell, type_fuel_cell::String, voltage_zone::String)::ExperimentalValues
Return the experimental values of the fuel cell.
"""
function experimental_values(fc::AbstractFuelCell, type_fuel_cell::String, voltage_zone::String)::ExperimentalValues
    throw(MethodError(experimental_values, (fc, type_fuel_cell, voltage_zone)))
end

"""
    numerical_params(fc::AbstractFuelCell)::NumericalParams
Return the numerical parameters of the fuel cell.
"""
function numerical_params(fc::AbstractFuelCell)::NumericalParams
    throw(MethodError(numerical_params, (fc,)))
end