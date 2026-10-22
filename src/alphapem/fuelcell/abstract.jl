"""
Abstract type for all fuel cell models.
"""

abstract type AbstractFuelCell end

using ...Config: PhysicalParams, OperatingConditions, PolaExperimentalData

# Abstract interface: all concrete fuel cell models must implement these methods

"""
    physical_parameters(fc::AbstractFuelCell)::PhysicalParams
Return the physical parameters of the fuel cell.
"""
function physical_parameters(fc::AbstractFuelCell)::PhysicalParams
    return PhysicalParams()
end

"""
    operating_conditions(fc::AbstractFuelCell, type_fuel_cell::Symbol)::OperatingConditions
Return the operating parameters of the fuel cell.
"""
function operating_conditions(fc::AbstractFuelCell, type_fuel_cell::Symbol)::OperatingConditions
    return OperatingConditions()
end

"""
    pola_exp_data(fc::AbstractFuelCell, type_fuel_cell::Symbol, voltage_zone::Symbol)::PolaExperimentalData
Return the polarization experimental data of the fuel cell.
"""
function pola_exp_data(fc::AbstractFuelCell, type_fuel_cell::Symbol, voltage_zone::Symbol)::PolaExperimentalData
    return PolaExperimentalData()
end

"""
    pola_exp_data_calibration(fc::AbstractFuelCell, type_fuel_cell::Symbol, voltage_zone::Symbol)::PolaExperimentalData
Return the polarization experimental data used for calibration.
"""
function pola_exp_data_calibration(fc::AbstractFuelCell, type_fuel_cell::Symbol, voltage_zone::Symbol)::PolaExperimentalData
    return PolaExperimentalData()
end

"""
    undetermined_parameter_bounds(fc::AbstractFuelCell, voltage_zone::Symbol = :full) -> Vector{Tuple{Symbol, Float64, Float64}}
Return the list of undetermined parameter bounds for the fuel cell.
"""
function undetermined_parameter_bounds(fc::AbstractFuelCell, voltage_zone::Symbol = :full)::Vector{Tuple{Symbol, Float64, Float64}}
    # PARAMETER_BOUNDS lists every parameter that may fluctuate, not the ones left undetermined:
    # the selection is made per fuel cell in the `data` directory.
    error("undetermined_parameter_bounds is not implemented for $(typeof(fc)).")
end


