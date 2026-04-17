# -*- coding: utf-8 -*-

"""
Example: faster EIS smoke run with live plotting.

This profile keeps the same EIS workflow as `run_EIS.jl` but uses:
- fewer frequencies,
- a high-frequency band,
- fewer points per period.

It is intended for quick validation of the EIS plotting pipeline.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Config: SimulationConfig, EISParams
using AlphaPEM.Application: run_simulation

current_params = EISParams(
    i_EIS = 1.0e4,        # (A/m²)
    ratio = 5.0 / 100.0,  # (-)
    f_power_min = 2.0,    # 10^2 Hz
    f_power_max = 3.0,    # 10^3 Hz
    nb_f = 3,             # minimal number of frequencies for a smoke run
    nb_points = 20,       # lighter Fourier grid
)

cfg = SimulationConfig(
    type_fuel_cell = :ZSW_GenStack,
    type_current = current_params,
    voltage_zone = :full,
    type_auxiliary = :no_auxiliary,
    type_purge = :no_purge,
    type_display = :synthetic,
    display_timing = :live,
)

start_time = time()
run_simulation(cfg)
algo_time = time() - start_time
println("Time of the algorithm in second : ", algo_time)
