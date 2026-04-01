# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :EIS
You can modify any SimulationConfig parameter below.
"""

include("../src/alphapem/config/simulation_config.jl")
include("../src/alphapem/application/run_simulation.jl")
include("../src/alphapem/config/current_parameters.jl")

using .SimulationConfigModule
using .Config: EISParams

current_params = EISParams(
    i_EIS = 1.0e4,        # (A/m²). Parameters for the EIS curve.
    ratio = 5.0 / 100.0,  # (.). Parameters for the EIS curve.
    f_power_min = -3.0,   # (.). Power of the minimum frequency for the EIS current density function.
    f_power_max = 5.0,    # (.). Power of the maximum frequency for the EIS current density function.
    nb_f = 90,            # (.). Number of frequencies tested for the EIS current density function.
    nb_points = 50,       # (.). Number of points calculated per specific period for the EIS current density function.
)

cfg = SimulationConfig(
    type_fuel_cells = [:ZSW_GenStack],
    type_current = current_params,
    voltage_zone = :full,
    type_auxiliary = :no_auxiliary,
    type_purge = :no_purge,
    type_display = :synthetic,
    type_plot = :dynamic
)

run_simulation(cfg)
