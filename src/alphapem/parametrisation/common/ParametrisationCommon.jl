# -*- coding: utf-8 -*-

"""
    ParametrisationCommon

Common utilities for AlphaPEM parameter identification, shared by calibration and
validity analysis modules.

This module handles parameter bounds definitions, physical parameter mapping,
sampling logic, and results export.
"""
module ParametrisationCommon

using Random
using Dates
using Printf
using YAML

using AlphaPEM.Config: PhysicalParams, PARAMETER_METADATA
using AlphaPEM.Fuelcell: create_fuelcell, DefaultFuelCell, undetermined_parameters

export ParameterBound,
       ParameterBounds,
       bounds_for_fuel_cell,
       new_PhysicalParams_from_sample,
       get_reference_config,
       export_parameter_bounds,
       export_calibrated_params

include("types.jl")
include("bounds.jl")
include("mapping.jl")
include("export.jl")

end # module ParametrisationCommon
