# -*- coding: utf-8 -*-

# run_calibration.jl
#
# Example script: parameter calibration for AlphaPEM using Genetic Algorithms.
#
# This script set up and launch a calibration process
# to identify the best undetermined physical parameters for a given fuel cell
# stack by matching experimental polarization curves.
#
# Usage:
#   julia --project=. examples/run_calibration.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Config: SimulationConfig, PolarizationCalibrationParams, GAConfig
using AlphaPEM.Parametrisation
using AlphaPEM.Parametrisation.Calibration
using Logging
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# THREADING
# ─────────────────────────────────────────────────────────────────────────────

const PARALLEL  = true   # true  → multi-threaded population evaluation
const N_THREADS = 0      # 0 → use all available cores

if PARALLEL
    n_desired = N_THREADS == 0 ? Sys.CPU_THREADS : min(N_THREADS, Sys.CPU_THREADS)
    if Threads.nthreads() < n_desired
        julia_bin = joinpath(Sys.BINDIR, "julia")
        project   = Base.active_project()
        script    = @__FILE__
        @info "Re-launching with $n_desired thread(s) (currently $(Threads.nthreads()))…"
        exit(run(`$julia_bin --threads=$n_desired --project=$project $script $ARGS`).exitcode)
    end
end

@info "Running with $(Threads.nthreads()) thread(s) (parallel simulation enabled)."

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

# ── Experimental conditions ───────────────────────────────────────────────────
# Each SimulationConfig represents one experimental dataset to calibrate on.
# All conditions must share the same fuel cell type (same parameters to identify).
# Add or remove blocks below to calibrate on more or fewer datasets.

calibration_conditions = [
    SimulationConfig(
        type_fuel_cell = :ZSW_GenStack,
        voltage_zone   = :before_voltage_drop,                # :full or :before_voltage_drop
    ),
    # SimulationConfig(
    #     type_fuel_cell = :ZSW_GenStack_Pa_1_61_Pc_1_41,
    #     voltage_zone   = :before_voltage_drop,
    # ),
]

# ── Genetic Algorithm settings ────────────────────────────────────────────────

ga_cfg = GAConfig(
    num_generations = 10_000,                       # 1000+ recommended for high precision
    pop_size        = 128,                      # 128+ recommended
    target_error    = 1/100,                    # Stop if RMSE < 1%
)

# ── Calibration config ────────────────────────────────────────────────────────

calib_cfg = CalibrationConfig(
    simulation_configs = calibration_conditions,
    ga_config          = ga_cfg,
    parallel           = PARALLEL,
    output_dir         = "results/calibration/ZSW_GenStack",
)

# ─────────────────────────────────────────────────────────────────────────────
# EXECUTION
# ─────────────────────────────────────────────────────────────────────────────

println("="^72)
println("  AlphaPEM — Parameter Calibration (Genetic Algorithm)")
println("="^72)
println()
ref = first(calibration_conditions)
@info "Starting calibration for $(ref.type_fuel_cell) — $(length(calibration_conditions)) condition(s)..."

# Run the calibration
result = calibrate(calib_cfg)

# ─────────────────────────────────────────────────────────────────────────────
# RESULTS
# ─────────────────────────────────────────────────────────────────────────────

println()
println("─"^72)
println("  Calibration Results")
println("─"^72)
@printf("  Final RMSE       : %.4f %%\n", result.min_rmse)
@printf("  Execution Time   : %.2f hours\n", result.execution_time / 3600)
println("  Results saved to : ", calib_cfg.output_dir)

println()
println("="^72)
println("  Calibration process complete.")
println("="^72)
