# -*- coding: utf-8 -*-

# run_sobol_sensitivity_analysis.jl
#
# Example script: Sobol global sensitivity analysis of AlphaPEM.
#
# This script runs a variance-based sensitivity analysis directly on AlphaPEM
# simulations (no surrogate model). It follows the workflow of the student notebook
# 03_sobol_analysis_AlphaPEM.ipynb:
#   1. Generate Sobol samples for physical + operating parameters.
#   2. Run AlphaPEM for each sample and extract the polarization curve.
#   3. Impute missing curves with KNN when simulations fail.
#   4. Compute AUC in three regions: activation, ohmic, mass transport.
#   5. Compute Sobol S1/ST (and optional S2) indices.
#   6. Export CSV/YAML results and generate plots.
#
# Usage:
#   julia --project=. examples/run_sobol_sensitivity_analysis.jl

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
        @info "Re-launching with $n_desired thread(s)..."
        exit(run(`$julia_bin --threads=$n_desired --project=$project $script $ARGS`).exitcode)
    end
end

@info "Running with $(Threads.nthreads()) thread(s)."

# ─────────────────────────────────────────────────────────────────────────────
# IMPORTS
# ─────────────────────────────────────────────────────────────────────────────

using AlphaPEM.Parametrisation.SobolSensitivityAnalysis
using AlphaPEM.Config: OperatingConditionConstraint

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

# NOTE: The total number of AlphaPEM simulations is N * (2 + n_params) for S1/ST
# and much larger for S2. With the default N = 10_000 and ~20 input parameters,
# this can be ~220 000 simulations. Adjust N according to available compute time.

# Example: fix some operating conditions to their nominal values so they are not
# sampled by Sobol. This reduces the number of input parameters and simulations.
excluded_oc = [:Sa, :Sc, :y_H2_in, :i_min_stoich]

# Example: define custom operating-condition constraints. The default constraint
# enforces Pc_des = Pa_des - 0.5e5. You can add or replace constraints here.
my_constraints = [
    OperatingConditionConstraint(
        target = :Pc_des,
        sources = [:Pa_des],
        fn = d -> d[:Pa_des] - 0.5e5
    ),
]

cfg = SobolAnalysisConfig(
    fuel_cell_type              = :ZSW_nominal,
    voltage_zone                = :full,
    N                           = 1024,          # Start small; increase once validated
    second_order                = false,         # Set true only if compute budget allows
    seed                        = 42,
    region_thresholds           = (0.4, 1.6),    # A/cm²: activation / ohmic / mass-transport
    include_operating_conditions = true,
    operating_condition_constraints = my_constraints,
    excluded_operating_conditions = excluded_oc,
    parallel                    = PARALLEL,
    max_run_time_s              = 60.0,
    knn_k                       = 10,
    output_dir                  = "results/sobol_sensitivity",
    save_curves                 = true,
)

# ─────────────────────────────────────────────────────────────────────────────
# RUN
# ─────────────────────────────────────────────────────────────────────────────

println("=" ^ 72)
println("  AlphaPEM — Sobol Global Sensitivity Analysis")
println("=" ^ 72)

result = run_sobol_analysis(cfg)

# ─────────────────────────────────────────────────────────────────────────────
# DISPLAY SUMMARY
# ─────────────────────────────────────────────────────────────────────────────

println()
println("─" ^ 72)
println("  Results summary")
println("─" ^ 72)
println("  Total time: $(round(result.execution_time; digits=1)) s")
println("  Output directory: $(result.config.output_dir)")
println()
println("  Top parameters per region (ST index):")
for r in result.regions
    df = result.sobol_indices[r.name]
    df_sorted = sort(df, :ST, rev = true)
    top = df_sorted.parameter[1:min(5, nrow(df_sorted))]
    println("    $(rpad(string(r.name), 15)) : $(join(top, ", "))")
end
println("=" ^ 72)
