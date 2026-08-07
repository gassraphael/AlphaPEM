# -*- coding: utf-8 -*-

"""
    DefaultFuelCell

A generic fallback FuelCell using default parameters, for unknown or manual types.
"""

mutable struct DefaultFuelCell <: AbstractFuelCell
    physical_parameters::PhysicalParams
    operating_conditions::OperatingConditions
    pola_exp_data::PolaExperimentalData
    pola_exp_data_cali::PolaExperimentalData
end

function DefaultFuelCell()
    return DefaultFuelCell(
        PhysicalParams(),
        OperatingConditions(),
        PolaExperimentalData(),
        PolaExperimentalData()
    )
end


"""
    GenericFuelCell

Generic concrete fuel cell type populated entirely from data files.
Used by the factory for any stack defined under `data/`.
"""
mutable struct GenericFuelCell <: AbstractFuelCell
    physical_parameters::PhysicalParams
    operating_conditions::OperatingConditions
    pola_exp_data::PolaExperimentalData
    pola_exp_data_cali::PolaExperimentalData
    undetermined_parameters::Dict{Symbol, Vector{Tuple{Symbol, Float64, Float64}}}
end

function GenericFuelCell(; physical_parameters::PhysicalParams,
                           operating_conditions::OperatingConditions,
                           pola_exp_data::PolaExperimentalData,
                           pola_exp_data_cali::PolaExperimentalData,
                           undetermined_parameters::Dict{Symbol, Vector{Tuple{Symbol, Float64, Float64}}})
    return GenericFuelCell(
        physical_parameters,
        operating_conditions,
        pola_exp_data,
        pola_exp_data_cali,
        undetermined_parameters
    )
end


"""
    undetermined_parameters(fc::GenericFuelCell, voltage_zone::Symbol = :full) -> Vector{Tuple{Symbol, Float64, Float64}}

Return the undetermined parameters stored in the fuel cell for the requested
voltage zone. The parameters are loaded once from the stack file when the fuel
cell is constructed.
"""
function undetermined_parameters(fc::GenericFuelCell, voltage_zone::Symbol = :full)::Vector{Tuple{Symbol, Float64, Float64}}
    voltage_zone in (:full, :before_voltage_drop) ||
        throw(ArgumentError("voltage_zone must be :full or :before_voltage_drop (got $voltage_zone)"))

    return fc.undetermined_parameters[voltage_zone]
end

