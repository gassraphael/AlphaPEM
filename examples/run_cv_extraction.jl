# -*- coding: utf-8 -*-

# run_cv_extraction.jl
#
# Example script: extract physical parameters from cyclic voltammetry (CV)
# curves for a ZSW fuel cell stack. All CV files in a directory are processed,
# a report is written, and the middle cell (by alphabetical order) is plotted.
#
# Usage:
#   julia --project=. examples/run_cv_extraction.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Parametrisation.CVExtraction
using Printf
using Statistics: mean

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

# Geometric electrode area reported in the experimental files (cm^2)
area_cm2 = 283.87

# Directory containing the experimental CV files
cv_dir = abspath(joinpath("data", "ZSW", "2026", "cv"))

# Extraction configuration
# These are the settings of the ZSW Matlab tool this module is ported from.
# All potential limits are in V vs the reference electrode used in the experiment.
cfg = CVExtractionConfig(
    area_cm2 = area_cm2,                    # geometric electrode area (cm²)
    cycle_for_extraction = 0,               # 0 = mean cycle used (what the ZSW tool reports)
    ignore_cycles_for_mean = [1],           # skip the first, not-yet-conditioned cycle
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

function _cell_number(filename::String)::Int
    m = match(r"cell_(\d+)", filename)
    return m === nothing ? 0 : parse(Int, m.captures[1])
end

cv_files = sort(
    filter(f -> endswith(lowercase(f), ".txt"), readdir(cv_dir));
    by = _cell_number,
)
if isempty(cv_files)
    error("No .txt CV files found in $cv_dir")
end

results = Dict{String, CVExtractionResult}()
for file in cv_files
    path = joinpath(cv_dir, file)
    @info "Extracting parameters from: $path"
    results[file] = extract_cv_parameters(path, cfg)
end
result_values = collect(values(results))

# ---------------------------------------------------------------------------
# RESULTS
# ---------------------------------------------------------------------------

println()
println("-"^72)
println("  Extracted parameters summary")
println("-"^72)

param_names = [
    "ECSA adsorption (cm2 Pt / cm2 electrode)",
    "ECSA desorption (cm2 Pt / cm2 electrode)",
    "H2 crossover (A cm-2)",
    "Double-layer cap. (F cm-2)",
    "Ohmic-drop slope (A V-1 cm-2)",
    "Scan rate (V s-1)",
]

values = [
    [r.ecsa_adsorption_cm2 for r in result_values],
    [r.ecsa_desorption_cm2 for r in result_values],
    [r.crossover_a_cm2 for r in result_values],
    [r.dlc_f_cm2 for r in result_values],
    [r.ohmic_drop_slope for r in result_values],
    [r.scan_rate_vs for r in result_values],
]

for (name, vals) in zip(param_names, values)
    @printf("  %-45s : mean = %.4e\n", name, mean(vals))
end

println()
println("="^72)
println("  CV extraction complete.")
println("="^72)

# Plot the middle cell (numerical order)
using CairoMakie

middle_file = cv_files[(length(cv_files) + 1) ÷ 2]
middle_result = results[middle_file]
fig = plot_cv_analysis(middle_result; title = "CV Analysis - $(middle_result.cv_file_name)")

plot_dir = joinpath(@__DIR__, "..", "results", "cv_extraction")
mkpath(plot_dir)
plot_path = joinpath(plot_dir, "$(middle_result.cv_file_name)_cv_analysis.png")
save(plot_path, fig)
@info "CV analysis plot saved to: $plot_path"
display(fig)

# Write report
report_path = joinpath(plot_dir, "$(middle_result.cv_file_name)_cv_analysis_report.txt")
open(report_path, "w") do io
    println(io, "AlphaPEM - CV Parameter Extraction Report")
    println(io, "="^72)
    println(io)
    println(io, "Directory: $cv_dir")
    println(io, "Number of cells: $(length(cv_files))")
    println(io)
    println(io, "Mean values")
    println(io, "-"^72)
    for (name, vals) in zip(param_names, values)
        @printf(io, "  %-45s : %.4e\n", name, mean(vals))
    end
    println(io)
    println(io, "Per-cell values")

    # Build a table whose column widths adapt to the longest header/value so
    # that long parameter names do not push the numbers out of alignment.
    short_names = [
        "Cell",
        "ECSA ads",
        "ECSA des",
        "H2 crossover",
        "DLC",
        "Ohmic slope",
        "Scan rate",
    ]
    header = short_names
    rows = Vector{String}[]
    for file in cv_files
        r = results[file]
        push!(rows, [
            r.cv_file_name,
            @sprintf("%.4e", r.ecsa_adsorption_cm2),
            @sprintf("%.4e", r.ecsa_desorption_cm2),
            @sprintf("%.4e", r.crossover_a_cm2),
            @sprintf("%.4e", r.dlc_f_cm2),
            @sprintf("%.4e", r.ohmic_drop_slope),
            @sprintf("%.4e", r.scan_rate_vs),
        ])
    end

    widths = [max(length(header[j]), maximum((length(r[j]) for r in rows); init = 0))
              for j in 1:length(header)]

    sep(l, m, r) = l * join((repeat("─", w + 2) for w in widths), m) * r
    fmt_row(cells) = "│" * join((" " * rpad(cells[j], widths[j]) * " │") for j in 1:length(header))

    println(io, sep("┌", "┬", "┐"))
    println(io, fmt_row(header))
    println(io, sep("├", "┼", "┤"))
    for r in rows
        println(io, fmt_row(r))
    end
    println(io, sep("└", "┴", "┘"))
end
@info "CV extraction report saved to: $report_path"
