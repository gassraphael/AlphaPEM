# -*- coding: utf-8 -*-

# run_simple_sensitivity_analysis.jl
#
# Example script: simple sensitivity analysis of AlphaPEM's PhysicalParams.
#
# This script activates the project, spawns the Distributed workers, broadcasts the
# sensitivity-analysis code to them, and launches the analysis.
#
# Usage:
#   julia --project=. examples/run_simple_sensitivity_analysis.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using Distributed
using AlphaPEM.Config: SimulationConfig

# Start workers if not already initialized
if nprocs() == 1
    addprocs(max(1, Sys.CPU_THREADS - 1))
end

# Broadcast the project activation and the sensitivity-analysis code (module imports,
# helpers and worker task logic) to all processes, including the master.
@everywhere begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)
    include(joinpath(@__DIR__, "..", "src", "alphapem", "parametrisation", "simple_sensitivity_analysis.jl"))
end

# ─────────────────────────────────────────────────────────────────────────────
# EXECUTION
# ─────────────────────────────────────────────────────────────────────────────

# Default configuration (can be edited below).
base_config = SimulationConfig(
    type_fuel_cell = :ZSW_nominal,
    voltage_zone   = :full,
    numerical_parameters = NumericalParams(nb_gc = 5)
)

run_simple_sensitivity_analysis(
    base_config;
    variation_pct      = 5.0,
    region_thresholds  = (0.5, 1.5),   # A/cm²: activation / ohmic / mass-transport
)
