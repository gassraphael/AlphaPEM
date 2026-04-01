# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :polarization
You can modify any SimulationConfig parameter below.
"""

include("../src/alphapem/config/simulation_config.jl")
include("../src/alphapem/application/run_simulation.jl")
include("../src/alphapem/config/current_parameters.jl")

using .SimulationConfigModule
using .Config: PolarizationParams

current_params = PolarizationParams(
    delta_t_ini = 120.0 * 60.0,  # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_i = 0.05e4,            # (A.m-2). Current density step for the polarisation current density function.
    v_load = 0.01e4,             # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break = 15.0 * 60.0, # (s). Breaking time for one step current, for the stabilisation of the internal states.
    i_max = 3.0e4                # Maximum current (default value, can be overridden by experimental current values if provided).
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
