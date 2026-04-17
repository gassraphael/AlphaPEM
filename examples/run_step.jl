# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :step
You can modify any SimulationConfig parameter below.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# --- Load AlphaPEM from the project environment ---
using AlphaPEM.Config: SimulationConfig, StepParams
using AlphaPEM.Application: run_simulation

current_params = StepParams(
    delta_t_ini = 30.0 * 60.0,   # (s). Initial time at zero current density for the stabilisation of the internal states (standard value).
    delta_t_load = 30.0,         # (s). Loading time for the step current density function, from 0 to i_step.
    delta_t_break = 2.0 * 60.0,  # (s). Time at i_step current density for the stabilisation of the internal states.
    i_ini = 1.0e4,               # (A.m-2). Initial current density for the step current density function.
    i_step = 1.5e4               # (A.m-2). Current density for the step current density function.
)

cfg = SimulationConfig(
    type_fuel_cell = :ZSW_GenStack,
    type_current = current_params,
    voltage_zone = :full, # :before_voltage_drop, :full.
    type_auxiliary = :no_auxiliary, # :forced_convective_cathode_with_anodic_recirculation, :forced_convective_cathode_with_flow_through_anode, :no_auxiliary.
    type_flow = :counter_flow, # :co_flow, :counter_flow.
    type_purge = :no_purge, # :constant_purge, :periodic_purge, :no_purge.
    type_display = :synthetic, # :multiple, :synthetic, :no_display.
    display_timing = :postrun, # :live, :postrun.
)

start_time = time() # Starting time
run_simulation(cfg)
algo_time = time() - start_time # Ending time
println("Time of the algorithm in second : ", algo_time)
