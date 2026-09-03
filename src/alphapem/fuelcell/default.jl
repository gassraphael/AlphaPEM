# -*- coding: utf-8 -*-

"""
    GenericFuelCell

Generic concrete fuel cell type populated entirely from data files.
Used by the factory for any stack defined under `data/`.
"""
const DiscretizationParams = Dict{Symbol, Vector{Tuple{Symbol, Float64, Float64}}}
const UndeterminedParameterSets = Dict{Symbol, DiscretizationParams}

mutable struct GenericFuelCell <: AbstractFuelCell
    physical_parameters::PhysicalParams
    operating_conditions::OperatingConditions
    pola_exp_data::PolaExperimentalData
    pola_exp_data_cali::PolaExperimentalData
    undetermined_parameter_bounds::DiscretizationParams
end

function GenericFuelCell(; physical_parameters::PhysicalParams,
                           operating_conditions::OperatingConditions,
                           pola_exp_data::PolaExperimentalData,
                           pola_exp_data_cali::PolaExperimentalData,
                           undetermined_parameter_bounds::DiscretizationParams)
    return GenericFuelCell(
        physical_parameters,
        operating_conditions,
        pola_exp_data,
        pola_exp_data_cali,
        undetermined_parameter_bounds
    )
end


"""
    physical_parameters(fc::GenericFuelCell) -> PhysicalParams

Return the physical parameters stored in the fuel cell. The parameters are
selected according to the GC spatial resolution at construction time.
"""
physical_parameters(fc::GenericFuelCell)::PhysicalParams = fc.physical_parameters


"""
    undetermined_parameter_bounds(fc::GenericFuelCell, voltage_zone::Symbol = :full) -> Vector{Tuple{Symbol, Float64, Float64}}

Return the undetermined parameter bounds stored in the fuel cell for the requested
voltage zone. The discretisation-aware set is selected at construction time.
"""
function undetermined_parameter_bounds(fc::GenericFuelCell, voltage_zone::Symbol = :full)::Vector{Tuple{Symbol, Float64, Float64}}
    voltage_zone in (:full, :before_voltage_drop) ||
        throw(ArgumentError("voltage_zone must be :full or :before_voltage_drop (got $voltage_zone)"))

    return fc.undetermined_parameter_bounds[voltage_zone]
end

