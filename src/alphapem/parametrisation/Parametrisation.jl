# -*- coding: utf-8 -*-

"""
    AlphaPEM.Parametrisation

This module provides parameter calibration tools used to optimize undetermined
fuel-cell parameters against experimental polarization curves.

Modules:
    - calibration_modules: Bounds, parameter updates, error calculation, and reporting helpers
"""
module Parametrisation

include("calibration_modules.jl")

export parameter_bounds_for_calibration,
       parameters_for_calibration,
       update_undetermined_parameters,
       calculate_simulation_error,
       print_calibration_results,
       save_calibration_results

end  # module Parametrisation

