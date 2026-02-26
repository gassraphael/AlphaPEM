# -*- coding: utf-8 -*-

"""AlphaPEM Configuration Module

This module contains configuration files for fuel cell parameters, current density profiles,
and experimental data used throughout the AlphaPEM package.

Modules:
    - current_densities: Current density profile functions (step, polarization, EIS)
    - specific_parameters: Fuel cell physical and operating parameters
    - pola_exp_values: Experimental polarization curve data
"""

from alphapem.config.current_densities import (
    step_current,
    polarization_current,
    polarization_current_for_calibration,
    EIS_current
)

from alphapem.config.parameters_specific import (
    stored_operating_inputs,
    stored_physical_parameters
)

from alphapem.config.pola_exp_values import (
    pola_exp_values,
    pola_exp_values_calibration,
    plot_experimental_polarisation_curve
)

__all__ = [
    "step_current",
    "polarization_current",
    "polarization_current_for_calibration",
    "EIS_current",
    "stored_operating_inputs",
    "stored_physical_parameters",
    "pola_exp_values",
    "pola_exp_values_calibration",
    "plot_experimental_polarisation_curve",
]

