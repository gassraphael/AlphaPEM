# -*- coding: utf-8 -*-

"""AlphaPEM Parameter Calibration Module

This module provides automated calibration tools for determining undetermined physical parameters
of fuel cell systems using genetic algorithms (PyGAD). The calibration process fits simulation
results to experimental polarization curves by optimizing parameters like layer thicknesses,
porosities, and electrochemical constants.

Key Features:
    - Genetic algorithm-based parameter optimization
    - Multi-objective calibration across different operating conditions
    - Support for multiple fuel cell types (ZSW-GenStack, EH-31)
    - Parallel processing capabilities for cluster computing

Modules:
    - calibration: Main calibration execution script
    - calibration_modules: Helper functions for bounds, parameter updates, and error calculation

Important Note:
    This module uses PyGAD (Python Genetic Algorithm), which is licensed under BSD-3-Clause.
    See LICENSE-BSD-3-CLAUSE in this directory for PyGAD's license terms.
    The rest of AlphaPEM is licensed under GPLv3 (see LICENSE at project root).

Usage:
    For local execution:
        python -m alphapem.parametrisation.calibration

    For HPC cluster execution:
        Use scripts/run_calibration_cluster.sh
        qsub scripts/run_calibration_cluster.sh calibration.py
"""

from alphapem.parametrisation.calibration_modules import (
    parameter_bounds_for_calibration,
    parameters_for_calibration,
    update_undetermined_parameters,
    calculate_simulation_error,
    print_calibration_results,
    save_calibration_results
)

__all__ = [
    "parameter_bounds_for_calibration",
    "parameters_for_calibration",
    "update_undetermined_parameters",
    "calculate_simulation_error",
    "print_calibration_results",
    "save_calibration_results",
]
