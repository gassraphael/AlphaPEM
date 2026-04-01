# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :step
You can modify any SimulationConfig parameter below.
"""

include("../src/alphapem/config/simulation_config.jl")
include("../src/alphapem/application/run_simulation.jl")
include("../src/alphapem/config/current_parameters.jl")

using .SimulationConfigModule
using .Config: StepParams

current_params = StepParams(
    delta_t_ini = 30.0 * 60.0,   # (s). Initial time at zero current density for the stabilisation of the internal states (standard value).
    delta_t_load = 30.0,         # (s). Loading time for the step current density function, from 0 to i_step.
    delta_t_break = 2.0 * 60.0,  # (s). Time at i_step current density for the stabilisation of the internal states.
    i_ini = 1.0e4,               # (A.m-2). Initial current density for the step current density function.
    i_step = 2.0e4               # (A.m-2). Current density for the step current density function.
)

cfg = SimulationConfig(
    type_fuel_cells = [:ZSW_GenStack],
    type_current = current_params,
    voltage_zone = :full,
    type_auxiliary = :no_auxiliary,
    type_purge = :no_purge,
    type_display = :synthetic,
    type_plot = :fixed
)

run_simulation(cfg)
