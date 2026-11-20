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
    - polarization_ascending: Single ascending sweep used by the extended local sensitivity analysis
    - eis: Electrochemical Impedance Spectroscopy (EIS) profile implementation
    - factory: Factory functions for current profiles
"""
module Currents

using ..Config: AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams,
                AscendingPolarizationParams, EISParams

include("abstract.jl")
include("step.jl")
include("polarization.jl")
include("polarization_for_cali.jl")
include("polarization_ascending.jl")
include("eis.jl")
include("factory.jl")

export AbstractCurrent, StepCurrent, PolarizationCurrent, PolarizationCalibrationCurrent,
       AscendingPolarizationCurrent, EISCurrent, create_current
export AbstractCurrentParams, StepParams, PolarizationParams, PolarizationCalibrationParams,
       AscendingPolarizationParams, EISParams
export solver_tstops, solver_dtmax

end  # module Currents

