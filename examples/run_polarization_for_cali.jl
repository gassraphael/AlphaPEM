# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :polarization_for_cali
You can modify any SimulationConfig parameter below.
"""

include("../src/alphapem/config/simulation_config.jl")
include("../src/alphapem/application/run_simulation.jl")

using .SimulationConfigModule

cfg = SimulationConfig(
    type_fuel_cells = [:ZSW_GenStack],
    type_current = :polarization_for_cali,
    voltage_zone = :full,
    type_auxiliary = :no_auxiliary,
    type_purge = :no_purge,
    type_display = :synthetic,
    type_plot = :fixed
)

run_simulation(cfg)
