# -*- coding: utf-8 -*-

"""
    AlphaPEM.Fuelcell

This module is the entry point for fuel cell models and their parameters in AlphaPEM.
It provides access to the main types and functions related to fuel cell modeling.

Fuel-cell data (physical parameters, nominal operating conditions, undetermined
parameter bounds and experimental polarization curves) are loaded from files
under `data/`. Adding a new stack only requires new data files; no source code
in `src/` needs to be modified.

Modules:
    - abstract: Abstract types and interface for fuel cell models
    - default: Default and generic concrete fuel cell types
    - factory: Factory function to instantiate fuel cells from data files
    - data_loader: Utilities for loading stack and polarization data
"""
module Fuelcell

using ..Config: PhysicalParams, OperatingConditions, PolaExperimentalData, NumericalParams

include("data_loader.jl")
include("abstract.jl")
include("default.jl")
include("factory.jl")

export create_fuelcell, DefaultFuelCell, GenericFuelCell, AbstractFuelCell,
       undetermined_parameters, physical_parameters, operating_conditions,
       pola_exp_data, pola_exp_data_calibration

end  # module Fuelcell

