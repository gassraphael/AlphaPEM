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
#   STEP 3 — (optional) Restrict the valid region via PRIM (requires R + IRD package)
#   STEP 4 — Export all results (CSV, YAML bounds, text reports)
#
# All output files are written to:
#   <project_root>/results/model_validity/   (fixed, see VALIDITY_OUTPUT_DIR)
#
# Usage:
#   julia --project=. --threads=auto examples/run_parameter_validity.jl
#
# For PRIM (STEP 3), R must be installed and the IRD package cloned:
#   git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/
#   See README.md § Installation — step 8 for full instructions.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Parametrisation.ValidParameterRegion

# ─────────────────────────────────────────────────────────────────────────────
# USER SETTINGS  (edit here)
# ─────────────────────────────────────────────────────────────────────────────

# Set to true to run PRIM (STEP 3).
# Requires R + the IRD package cloned at the default IRD_PACKAGE_DIR
# (see README § Installation).
const RUN_PRIM = lowercase(get(ENV, "RUN_PRIM", "false")) in ("1", "true", "yes")


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
    fuel_cell_type       = :ZSW_GenStack,
    voltage_zone         = :full,
    n_samples            = 100,
    sampling_seed        = 42,
    validation_criteria  = criteria_cfg,
    parallel             = true,
    checkpoint_interval  = 100,
)

# PRIM configuration (only used when RUN_PRIM = true)
prim_cfg = RUN_PRIM ? PRIMConfig(
    ird_package_dir   = "external/supplementary_2023_ird/irdpackage",
    methods           = [:PRIM, :MaxBox],
    probability_range = (0.8, 1.0),
    seed              = 42,
) : nothing

# ─────────────────────────────────────────────────────────────────────────────
# RUN  (full pipeline in one call)
# ─────────────────────────────────────────────────────────────────────────────

println("="^72)
println("  AlphaPEM — Parameter Validity Region Analysis")
println("="^72)
println()

result = run_validity_analysis(analysis_cfg, prim_cfg)

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
s.failed_count > 0 && println("  Failed      : ", s.failed_count)

println()
println("─"^72)
println("  Generated files")
println("─"^72)
for (k, v) in sort(collect(result.output_files); by = x -> string(x[1]))
    println("  $(rpad(string(k), 38)) → $v")
end

if RUN_PRIM && !isempty(result.prim_results)
    println()
    println("─"^72)
    println("  PRIM bounds summary")
    println("─"^72)
    for r in result.prim_results
        println("  Method    : ", r.method)
        println("  Precision : ", round(r.precision; digits = 3))
        println("  Coverage  : ", round(r.coverage;  digits = 3))
        println()
    end
end

println()
println("="^72)
println("  Pipeline steps")
println("="^72)
println("    ✅  STEP 1 — LHS sampling                         (complete)")
println("    ✅  STEP 2 — Batch simulation & classification    (complete)")
if RUN_PRIM
    println("    ✅  STEP 3 — PRIM restriction via R irdpackage   (complete)")
else
    println("    -- STEP 3 — PRIM restriction via R irdpackage    (skipped: RUN_PRIM=false)")
end
println("    ✅  STEP 4 — Export results                       (complete)")
println("="^72)

