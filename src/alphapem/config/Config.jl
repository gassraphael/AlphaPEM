# -*- coding: utf-8 -*-

"""
    AlphaPEM.Config

This module is the configuration entry point for AlphaPEM,
and provides access to the main configuration parameters.

Modules:
    - current_parameters: Structures for current density parameters
    - fuel_cell_parameters: Structures for physical, operating, numerical, and experimental parameters
    - simulation_config: Structure and validation for simulation configuration parameters

Exports:
    - AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
    - AbstractFuelCellParams, PhysicalParams, OperatingConditions, PolaExperimentalData, NumericalParams
    - SimulationConfig, validate_config
"""
module Config

include("current_parameters.jl")
include("fuel_cell_parameters.jl")
include("simulation_config.jl")

export AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
export AbstractFuelCellParams, PhysicalParams, OperatingConditions, PolaExperimentalData, NumericalParams
export SimulationConfig, validate_config

end  # module Config

