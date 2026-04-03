# -*- coding: utf-8 -*-

"""
    AlphaPEM.Currents

This module is the entry point for current density profiles and related calculations in AlphaPEM.
It provides access to the main types and functions for current density modeling and analysis.

Modules:
    - abstract: Abstract types for current profiles
    - step: Step current profile implementation
    - polarization: Polarization current profile implementation
    - polarization_for_cali: Polarization profile for calibration
    - eis: Electrochemical Impedance Spectroscopy (EIS) profile implementation
    - factory: Factory functions for current profiles
"""
module Currents

using ..Config: AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams

include("abstract.jl")
include("step.jl")
include("polarization.jl")
include("polarization_for_cali.jl")
include("eis.jl")
include("factory.jl")

export AbstractCurrent, StepCurrent, PolarizationCurrent, PolarizationCalibrationCurrent, EISCurrent, create_current
export AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams

end  # module Currents

