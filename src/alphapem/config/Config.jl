# -*- coding: utf-8 -*-

"""
    AlphaPEM.Config

This module contains configuration files for fuel cell parameters, current density profiles,
and experimental data used throughout the AlphaPEM package.

Modules:
    - CurrentDensities: Current density profile functions (step, polarization, EIS)
    - Parameters: Fuel cell physical and operating parameters
    - ParametersSpecific: Fuel cell stored physical and operating parameters
    - PolaExpValues: Experimental polarization curve data
"""
module Config

module CurrentDensities
include("current_densities.jl")

export step_current,
       polarization_current,
       polarization_current_for_calibration,
       EIS_current

end  # module CurrentDensities

module Parameters
include("parameters.jl")

export calculate_current_density_parameters,
       calculate_operating_inputs,
       calculate_physical_parameters,
       calculate_computing_parameters

end  # module Parameters

module ParametersSpecific
include("parameters_specific.jl")

export stored_operating_inputs,
       stored_physical_parameters

end  # module ParametersSpecific

module PolaExpValues
include("pola_exp_values.jl")

export pola_exp_values,
       pola_exp_values_calibration,
       plot_experimental_polarisation_curve

end  # module PolaExpValues

using .CurrentDensities: step_current,
                         polarization_current,
                         polarization_current_for_calibration,
                         EIS_current

using .Parameters: calculate_current_density_parameters,
                   calculate_operating_inputs,
                   calculate_physical_parameters,
                   calculate_computing_parameters

using .ParametersSpecific: stored_operating_inputs,
                           stored_physical_parameters

using .PolaExpValues: pola_exp_values,
                      pola_exp_values_calibration,
                      plot_experimental_polarisation_curve

export CurrentDensities,
       Parameters,
       ParametersSpecific,
       PolaExpValues,
       step_current,
       polarization_current,
       polarization_current_for_calibration,
       EIS_current,
       calculate_current_density_parameters,
       calculate_operating_inputs,
       calculate_physical_parameters,
       calculate_computing_parameters,
       stored_operating_inputs,
       stored_physical_parameters,
       pola_exp_values,
       pola_exp_values_calibration,
       plot_experimental_polarisation_curve

end  # module Config

