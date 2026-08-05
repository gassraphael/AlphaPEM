# -*- coding: utf-8 -*-

# run_cv_extraction.jl
#
# Example script: extract physical parameters from a cyclic voltammetry (CV)
# curve for a ZSW fuel cell stack.
#
# Usage:
#   julia --project=. examples/run_cv_extraction.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Parametrisation.CVExtraction
using Printf

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

# Geometric electrode area reported in the experimental file (cm^2)
area_cm2 = 283.87

# Path to the experimental CV file (absolute path from the current working directory)
cv_file = abspath(joinpath("data", "experimental", "ZSW", "2026", "cv", "cv_cell_10.txt"))

# Extraction configuration
# The default limits are taken from the initial file this program was inspired by.
# All potential limits are in V vs the reference electrode used in the experiment.
cfg = CVExtractionConfig(
    area_cm2 = area_cm2,                    # geometric electrode area (cm²)
    cycle_for_extraction = 3,               # use the 3rd CV cycle (0 = mean cycle)
    double_layer_limit_min = 0.30,          # V vs reference electrode
    double_layer_limit_max = 0.50,          # V vs reference electrode
    ohmic_drop_limit_min = 0.35,            # V vs reference electrode
    ohmic_drop_limit_max = 0.45,            # V vs reference electrode
    ecsa_limit_min = 0.09,                  # V vs reference electrode
    ecsa_limit_max = 0.50,                  # V vs reference electrode
    conversion_factor_uc_cm2 = 210.0,       # μC·cm⁻², H monolayer charge on Pt
    covering_degree = 0.77,                 # empirical H coverage on Pt (0–1)
)

# ---------------------------------------------------------------------------
# EXECUTION
# ---------------------------------------------------------------------------

println("="^72)
println("  AlphaPEM - CV Parameter Extraction")
println("="^72)
println()

@info "Extracting parameters from: $cv_file"
result = extract_cv_parameters(cv_file, cfg)

# ---------------------------------------------------------------------------
# RESULTS
# ---------------------------------------------------------------------------

println()
println("-"^72)
println("  Extracted parameters for $(result.cv_file_name)")
println("-"^72)
@printf("  ECSA (H2 adsorption) : %.2f cm2 Pt / cm2 electrode\n", result.ecsa_adsorption_cm2)
@printf("  ECSA (H2 desorption) : %.2f cm2 Pt / cm2 electrode\n", result.ecsa_desorption_cm2)
@printf("  H2 crossover         : %.4e A cm-2\n", result.crossover_a_cm2)
@printf("  Double-layer cap.    : %.4e F cm-2\n", result.dlc_f_cm2)
@printf("  Ohmic-drop slope     : %.4e A V-1 cm-2\n", result.ohmic_drop_slope)
@printf("  Scan rate            : %.4f V s-1\n", result.scan_rate_vs)
println()
println("="^72)
println("  CV extraction complete.")
println("="^72)

# Plot original and corrected CV side by side
using CairoMakie
fig = plot_cv_analysis(result; title = "CV Analysis - $(result.cv_file_name)")

plot_dir = joinpath(@__DIR__, "..", "results", "cv_extraction")
mkpath(plot_dir)
plot_path = joinpath(plot_dir, "$(result.cv_file_name)_cv_analysis.png")
save(plot_path, fig)
@info "CV analysis plot saved to: $plot_path"
