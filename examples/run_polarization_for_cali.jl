# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :polarization_for_cali
You can modify any SimulationConfig parameter below.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Config: SimulationConfig, PolarizationCalibrationParams
using AlphaPEM.Application: run_simulation

current_params = PolarizationCalibrationParams(
    delta_t_ini= 120.0 * 60.0,    # (s). Initial time at zero current density for the stabilisation of the internal states.
    v_load = 0.01e4,              # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break = 10.0 * 60.0,  # (s). Breaking time for one step current, for the stabilisation of the internal states.
    i_exp = Float64[]             # Experimental current density values (A/m²) for calibration
)

cfg = SimulationConfig(
    type_fuel_cell = :ZSW_GenStack,
    type_current = current_params,
    voltage_zone = :full, # :before_voltage_drop, :full.
    type_auxiliary = :no_auxiliary, # :forced_convective_cathode_with_anodic_recirculation, :forced_convective_cathode_with_flow_through_anode, :no_auxiliary.
    type_purge = :no_purge, # :constant_purge, :periodic_purge, :no_purge.
    type_display = :synthetic, # :multiple, :synthetic, :no_display.
    display_timing = :postrun # :live, :postrun.
)

start_time = time() # Starting time
run_simulation(cfg)
algo_time = time() - start_time # Ending time
println("Time of the algorithm in second : ", algo_time)
