"""
Abstract type for all fuel cell models.
"""

abstract type AbstractFuelCell end

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
    undetermined_parameters(fc::AbstractFuelCell, voltage_zone::Symbol) -> Vector{Tuple{Symbol, Float64, Float64}}
Return the list of undetermined parameters for the fuel cell.
"""
function undetermined_parameters(fc::AbstractFuelCell, voltage_zone::Symbol = :full)::Vector{Tuple{Symbol, Float64, Float64}}
    # Generic implementation using default bounds from Config
    # This can be overridden by specific fuel cell types
    
    # We use a standard set of parameters that are typically undetermined in PEMFC models
    params = [
        (:Hacl,          5e-6, 20e-6),
        (:Hccl,          5e-6, 20e-6),
        (:Hmem,          5e-6, 50e-6),
        (:Hgdl,        100e-6, 150e-6),
        (:Hmpl,         40e-6, 100e-6),
        (:epsilon_gdl,    0.5, 0.9),
        (:e,              3.0, 5.0),
        (:Re,            5e-8, 5e-6),
        (:i0_c_ref,       0.1, 80.0),
        (:kappa_co,      0.01, 40.0),
        (:kappa_c,       0.25, 4.0),
    ]
    
    if voltage_zone == :full
        push!(params, (:K_O2_ad_Pt, 0.1, 10.0))
    end
    
    return params
end

