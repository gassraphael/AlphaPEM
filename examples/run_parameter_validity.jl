# -*- coding: utf-8 -*-

# run_parameter_validity.jl
#
# Example script: parameter validity region analysis for AlphaPEM.
#
# This script runs the complete pipeline to identify the sub-region of
# AlphaPEM's undetermined-parameter space that reliably produces physical
# polarization curves.
#
# Workflow:
#   STEP 1 — Draw configurations by Latin Hypercube Sampling
#   STEP 2 — Simulate each configuration with AlphaPEM and classify the curve
#   STEP 3 — Restrict the valid region via IRD methods (PRIM, MaxBox, etc. — requires R + IRD package)
#   STEP 4 — Export all results (CSV, YAML bounds, text reports)
#
# All output files are written to:
#   <project_root>/results/model_validity/   (fixed, see VALIDITY_OUTPUT_DIR)
#
# Usage:
#   julia --project=. examples/run_parameter_validity.jl
#
# For IRD methods (STEP 3), R must be installed and the IRD package cloned:
#   git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/IRD_method_2023
#   See README.md § Installation — step 8 for full instructions.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# ─────────────────────────────────────────────────────────────────────────────
# THREADING
#
# These two constants must stay at the very top of the script, BEFORE any
# heavy `using` statement.  Julia's thread count is fixed at startup; if more
# threads are needed the script re-executes itself immediately with the correct
# --threads flag, before paying the cost of loading any package.
# ─────────────────────────────────────────────────────────────────────────────

const PARALLEL  = true   # true  → multi-threaded batch simulation
                         # false → sequential (useful for debugging)

const N_THREADS = 0      # Number of Julia threads to request.
                         #   0  → use all available cores (Sys.CPU_THREADS)
                         #   N  → use exactly min(N, Sys.CPU_THREADS) threads

# ── Auto-restart with the right thread count ─────────────────────────────────
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

@info "Running with $(Threads.nthreads()) thread(s)$(PARALLEL ? " — parallel batch enabled" : " — sequential mode")."

# ─────────────────────────────────────────────────────────────────────────────
# IMPORTS  (after the restart check — only reached with the correct thread count)
# ─────────────────────────────────────────────────────────────────────────────

using AlphaPEM.Config: NumericalParams
using AlphaPEM.Parametrisation.ValidParameterRegion
using Logging
using Printf

# ─────────────────────────────────────────────────────────────────────────────
# LOGGING FILTER
#
# During batch simulation the IDA solver may trigger "Safety stop" warnings
# whenever a configuration exceeds the fuel cell's physical operating limit and
# the simulation is terminated gracefully.  These warnings are expected, already
# captured in the :invalid classification, and only pollute the progress bar
# display.  The filter below silences them while keeping all other
# Info / Warn / Error messages intact.
# ─────────────────────────────────────────────────────────────────────────────

struct _SuppressSafetyWarnings <: AbstractLogger
    base::AbstractLogger
end
Logging.min_enabled_level(l::_SuppressSafetyWarnings) = Logging.min_enabled_level(l.base)
Logging.shouldlog(l::_SuppressSafetyWarnings, args...)  = Logging.shouldlog(l.base, args...)
function Logging.handle_message(l::_SuppressSafetyWarnings, level, msg, args...; kwargs...)
    level == Logging.Warn && startswith(string(msg), "Safety stop") && return
    Logging.handle_message(l.base, level, msg, args...; kwargs...)
end
Logging.catch_exceptions(l::_SuppressSafetyWarnings) = Logging.catch_exceptions(l.base)

# ─────────────────────────────────────────────────────────────────────────────
# BUILD CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

criteria_cfg = ValidityCriteriaConfig(
    voltage_range           = (0.0, 1.23),   # V — acceptable starting-voltage range
    monotonic_threshold     = 0.005,         # V — max upward bump tolerated (5 mV)
    apply_start_range       = true,
    apply_monotonicity      = true,
    apply_positive_voltages = true,
)

analysis_cfg = ValidityAnalysisConfig(
    fuel_cell_type         = :ZSW_GenStack,
    voltage_zone           = :before_voltage_drop,  # :before_voltage_drop, :full.
    n_samples              = 10_000,                # Total number of configurations to simulate (LHS samples).
    validation_criteria    = criteria_cfg,
    parallel               = PARALLEL,              # ← driven by the constant above
    save_curves            = true,                  # Set to true to save polarization curves
    reuse_from             = "results/model_validity/2026.06.02 - 10000 samples - before voltage drop - V1",               # Set to "path/to/previous/run" to reuse curves. ex: "results/model_validity/2026.06.02 - 10000 samples - before voltage drop - V1"
    hyperbox_finder_method = [:PRIM, :MaxBox],                 # Vector of IRD methods: :PRIM, :MaxBox
    max_run_time_s         = 30.0,       # Maximum simulation runtime (seconds)
)

# IRD configuration (required — STEP 3 is no longer optional)
# Validate the methods and build the IRD config
allowed = Set([:PRIM, :MaxBox])
methods_vec = if isa(analysis_cfg.hyperbox_finder_method, Symbol)
    # Single method: convert to vector
    if !(analysis_cfg.hyperbox_finder_method in allowed)
        error("Invalid hyperbox_finder_method: $(analysis_cfg.hyperbox_finder_method). Use :PRIM, :MaxBox, or a vector of these.")
    end
    [analysis_cfg.hyperbox_finder_method]
else
    # Vector of methods: validate each one
    isempty(analysis_cfg.hyperbox_finder_method) &&
        error("hyperbox_finder_method cannot be empty. Provide at least one IRD method (:PRIM, :MaxBox, or others).")
    for m in analysis_cfg.hyperbox_finder_method
        if !(m in allowed)
            error("Invalid hyperbox_finder_method: $m. Use :PRIM or :MaxBox only.")
        end
    end
    analysis_cfg.hyperbox_finder_method
end

ird_cfg = IRDConfig(
    ird_package_dir   = "external/IRD_method_2023/irdpackage",
    methods           = methods_vec,
    probability_range = (0.8, 1.0),
    seed              = 42,
)

# ─────────────────────────────────────────────────────────────────────────────
# RUN  (full pipeline in one call)
# ─────────────────────────────────────────────────────────────────────────────

println("="^72)
println("  AlphaPEM — Parameter Validity Region Analysis")
println("="^72)
println()

result = with_logger(_SuppressSafetyWarnings(current_logger())) do
    run_validity_analysis(analysis_cfg, ird_cfg)
end

# ─────────────────────────────────────────────────────────────────────────────
# DISPLAY RESULTS
# ─────────────────────────────────────────────────────────────────────────────

s = result.validation_summary
valid_pct = round(s.valid_count / max(1, s.total_simulations) * 100; digits = 1)

println()
println("─"^72)
println("  Classification results")
println("─"^72)
println("  Simulations : ", s.total_simulations)
println("  Valid       : ", s.valid_count, "  (", valid_pct, " %)")
println("  Invalid     : ", s.invalid_count)
println("  Failed      : ", s.failed_count)

println()
println("─"^72)
println("  Restricted region analysis (IRD methods)")
println("─"^72)

# Display results for each method
for r in result.prim_results
    println("  $(rpad(string(r.method), 10)) │ Precision: $(round(r.precision; digits=3))  Coverage: $(round(r.coverage; digits=3))  Valid in box: $(r.n_valid_inside)/$(r.n_inside_box)")
end

println()
println("─"^72)
println("  Key output files")
println("─"^72)

# Display essential files
if haskey(result.output_files, :parameter_classification_csv)
    println("  Classifications : ", basename(result.output_files[:parameter_classification_csv]))
end
if haskey(result.output_files, :generated_curves_csv)
    println("  Curves          : ", basename(result.output_files[:generated_curves_csv]))
end
if haskey(result.output_files, :bounds_restricted_yaml)
    println("  Bounds          : ", basename(result.output_files[:bounds_restricted_yaml]))
end
if haskey(result.output_files, :final_report_txt)
    println("  Report          : ", basename(result.output_files[:final_report_txt]))
end

println()
println("="^72)
println("  Pipeline steps")
println("="^72)
println("    ✅  STEP 1 — LHS sampling                         (complete)")
println("    ✅  STEP 2 — Batch simulation & classification    (complete)")
println("    ✅  STEP 3 — IRD methods analysis (PRIM/MaxBox)   (complete)")
println("    ✅  STEP 4 — Export results                       (complete)")
println("="^72)
