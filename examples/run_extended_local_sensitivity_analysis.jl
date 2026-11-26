# -*- coding: utf-8 -*-

# run_extended_local_sensitivity_analysis.jl
#
# Example script: extended local sensitivity and parameter-subset selection of AlphaPEM.
#
# This script only sets the analysis parameters and displays the outcome. The analysis itself is
# coordinated step by step by `run_extended_local_sensitivity_analysis`
# (src/alphapem/parametrisation/extended_local_sensitivity_analysis.jl).
#
# Long runs:
#   - A checkpoint is written regularly. If the run is interrupted, set RESUME = true and launch the
#     script again to continue where it stopped.
#   - When the simulations are finished, the complete result is archived. Set LOAD_SAVED_RESULT = true
#     to repeat the analysis (thresholds, figures, selection) without simulating again.
#
# Usage:
#   julia --project=. examples/run_extended_local_sensitivity_analysis.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# ─────────────────────────────────────────────────────────────────────────────
# THREADING
# ─────────────────────────────────────────────────────────────────────────────

const PARALLEL  = true
const N_THREADS = 0   # 0 → use all available cores

if PARALLEL
    n_desired = N_THREADS == 0 ? Sys.CPU_THREADS : min(N_THREADS, Sys.CPU_THREADS)
    if Threads.nthreads() < n_desired
        julia_bin = joinpath(Sys.BINDIR, "julia")
        project   = Base.active_project()
        script    = @__FILE__
        @info "Re-launching with $n_desired thread(s) (currently $(Threads.nthreads()))..."
        exit(run(`$julia_bin --threads=$n_desired --project=$project $script $ARGS`).exitcode)
    end
end

@info "Running with $(Threads.nthreads()) thread(s)."

# ─────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ─────────────────────────────────────────────────────────────────────────────

using AlphaPEM.Application: wait_for_interactive_figures!
using AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

# Sobol balance is guaranteed for a power-of-two sample size. A non-power-of-two edit is rounded up automatically.
const N_CANDIDATES = 2^9
const CURRENT_DENSITY_GRID = Float64[
    0.0, 0.03, 0.05, 0.1, 0.2, 0.6, 1.0, 1.2, 1.5, 1.7, 1.9, 2.1, 2.4, 2.7, 3.0
]
const OPERATING_BOUNDS_OVERRIDES = Dict{Symbol, Tuple{Float64, Float64, Symbol}}()

# Resume an interrupted run from its checkpoint (the screening is not repeated).
const RESUME = false
# Skip the simulations: reload the archived result of a previous run and only repeat the analysis.
const LOAD_SAVED_RESULT = false

cfg = ExtendedLocalSensitivityConfig(
    n_candidates = N_CANDIDATES,                                 # Sobol candidates generated before experimental screening.
    current_density_grid = CURRENT_DENSITY_GRID,                 # Shared current-density points for every operating case (A/cm²).
    current_density_zones = copy(DEFAULT_CURRENT_DENSITY_ZONES), # Magnitude summaries by current range; copy before customizing zone boundaries.
    pressure_gap = 0.2e5,                                        # Exact pressure difference: Pa_des - Pc_des (Pa).
    operating_bounds_overrides = OPERATING_BOUNDS_OVERRIDES,     # Optional partial bound overrides; empty keeps framework defaults.
    base_operating_conditions = nothing,                         # Use selected fuel-cell nominal conditions as the base; swept fields vary, others stay nominal.
    sensitivity_threshold = 0.1,                                 # FIGURES AND SUBSET SELECTION: minimum median voltage-sensitivity magnitude for a parameter to count as non-negligible; a collinear pair needs both parameters above it, and a retained parameter must reach it to be selected.
    collinearity_threshold = 0.707,                              # FIGURES: collinearity index above which a pair is highlighted (angle below 45°).
    psi_th = 0.9,                                                # SUBSET SELECTION: median collinearity above which the less sensitive parameter of a pair is removed from the candidates.
    parallel = PARALLEL,                                         # Run independent simulations concurrently when several Julia threads exist.
    resume = RESUME,                                             # Continue from the checkpoint of an interrupted run.
    checkpoint_frequency = 10,                                   # Save the checkpoint every this many finished work units (0 disables it).
    save_result = true,                                          # Archive the complete result at the end so the analysis can be repeated without simulating.
    screening_datasets = [
        ScreeningDataset(:ZSW_nominal, 2024, :full),             # Candidate must match this experimental polarization curve.
        ScreeningDataset(:ZSW_Pa_1_61_Pc_1_41, 2024, :full),     # Candidate must also match this curve.
    ],                                                           # A candidate is retained only if it passes screening on every listed dataset.
)
validate_config(cfg)

# ─────────────────────────────────────────────────────────────────────────────
# RUN
# ─────────────────────────────────────────────────────────────────────────────

println("="^72)
println("  AlphaPEM — Extended Local Sensitivity Analysis")
println("="^72)

# Either simulate (or resume), or reload the archive and repeat only the analysis.
outcome = LOAD_SAVED_RESULT ?
    reanalyze_extended_local_sensitivity(cfg) :
    run_extended_local_sensitivity_analysis(cfg)

# ─────────────────────────────────────────────────────────────────────────────
# DISPLAY SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

result, analysis = outcome.result, outcome.analysis
println()
println("─"^72)
println("  Results summary")
println("─"^72)
println("  Accepted candidates : $(length(result.candidates)) / $(length(outcome.records))")
println("  Operating cases     : $(length(result.operating_cases))")
println("  Execution time      : $(round(result.execution_time / 3600; digits = 2)) h")
println("  Output directory    : $(cfg.output_dir)")
isempty(outcome.fixed_parameters) ||
    println("  Fixed parameters    : $(join(string.(outcome.fixed_parameters), ", "))")
println()
println("  Collinear pairs above thresholds : $(size(analysis.pairs, 1))")
println("  Removed by the collinearity filter (psi_th = $(cfg.psi_th)):")
for row in eachrow(analysis.collinearity_filter.removed)
    println("    $(rpad(row.removed, 16)) collinear with $(rpad(row.kept, 16)) (psi = $(round(row.psi_median; digits = 2)))")
end
println()
println("  Selected parameters (median sensitivity >= $(cfg.sensitivity_threshold)): $(length(analysis.selection.selected))")
println("    ", join(string.(analysis.selection.selected), ", "))
println("="^72)

# In scripts, keep the browser page alive after all figures are complete.
wait_for_interactive_figures!(
    "Sensitivity figures are complete. Press Ctrl+C to close the browser plot server.")
