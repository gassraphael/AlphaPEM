"""
Abstract type for all fuel cell models.
"""

abstract type AbstractFuelCell end

using ...Config: PhysicalParams, OperatingConditions, PolaExperimentalData, UNDETERMINED_PARAMETER_BOUNDS

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
    # Generic implementation: delegate to the canonical bounds defined in
    # `fuel_cell_parameters.jl` (UNDETERMINED_PARAMETER_BOUNDS).
    params = Tuple{Symbol, Float64, Float64}[]
    for (name, (minv, maxv, _kind)) in UNDETERMINED_PARAMETER_BOUNDS
        # Preserve historical behaviour: `:K_O2_ad_Pt` was only present in the
        # undetermined list for the full voltage zone.
        if voltage_zone != :full && name == :K_O2_ad_Pt
            continue
        end
        push!(params, (name, minv, maxv))
    end
    return params
end


