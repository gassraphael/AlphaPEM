# -*- coding: utf-8 -*-

"""
    AlphaPEM.Fuelcell

This module is the entry point for fuel cell models and their parameters in AlphaPEM.
It provides access to the main types and functions related to fuel cell modeling.

Modules:
    - abstract: Abstract types for fuel cell models
    - eh31: EH-31 fuel cell model implementation
    - zsw: ZSW fuel cell model implementation
"""
module Fuelcell

include("../config/fuel_cell_parameters.jl")   # PhysicalParams, OperatingConditions, etc.

include("abstract.jl")
include("zsw.jl")
include("eh31.jl")
include("factory.jl")
include("default.jl")

export create_fuelcell, DefaultFuelCell, AbstractFuelCell, EH31FuelCell, ZSWFuelCell

end  # module Fuelcell

