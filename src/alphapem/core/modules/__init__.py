# -*- coding: utf-8 -*-

"""AlphaPEM Computational Modules

This module provides access to all the physical simulation modules used for modeling
different aspects of the PEMFC system including:
- Fluid flow dynamics (1D manifold and MEA regions)
- Heat transfer and thermal management
- Electrochemical reactions and voltage calculations
- Differential equations solver
- Display management

These modules implement the finite-volume method and handle the physics-based calculations
for the PEMFC simulator.
"""

__all__ = [
    # Main simulation modules
    "run_simulation_modules.py",
    "dif_eq_modules",
    "flows_1D_MEA_modules",
    "flows_1D_GC_manifold_modules",
    "heat_modules",
    "cell_voltage_modules",
    # Display and GUI modules
    "display_modules",
]

