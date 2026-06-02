# -*- coding: utf-8 -*-

"""
    ResultsExport

Module for saving validity-analysis and PRIM outputs to disk.

Handled formats:
- **CSV** — classified configurations (parameter columns + `classification`).
- **YAML** — parameter bounds, original and PRIM-restricted, in a format compatible
  with downstream sampling and re-use across runs.
- **Plain text** — validation summary and comparative bounds report.

# Exports
- `ExportConfig`: Output directory and file-naming options
- `ValidationSummary`: Aggregated statistics from a batch classification run
- `export_classified_configurations`: Write classified configs to CSV
- `export_parameter_bounds`: Write parameter bounds to YAML
- `export_prim_report`: Write a human-readable comparison report
- `export_validation_summary`: Write the batch-run statistics to a text file
- `generate_comparison_table`: Format a side-by-side table of original vs. restricted bounds
"""
module ResultsExport

using CSV
using DataFrames
using Dates
using Printf
using YAML

export ExportConfig,
       ValidationSummary,
       export_classified_configurations,
       export_parameter_bounds,
       export_prim_report,
       export_validation_summary,
       generate_comparison_table,
       load_restricted_bounds

# ─────────────────────────────────────────────────────────────────────────────
# DATA STRUCTURES
# ─────────────────────────────────────────────────────────────────────────────

"""
    ExportConfig

File-naming and directory options for all export functions.

# Fields
- `output_dir::String`: Root directory for all output files.  Created if absent.
  Default `"results/model_validity"`.
- `prefix::String`: Optional prefix prepended to every file name. Default `""`.
- `timestamp::Bool`: Append a `YYYYMMDD_HHMMSS` timestamp to file names.  Default `true`.
- `overwrite::Bool`: Silently overwrite existing files without prompting.  Default `false`.
"""
Base.@kwdef struct ExportConfig
    output_dir::String  = "results/model_validity"
    prefix::String      = ""
    timestamp::Bool     = true
    overwrite::Bool     = false
end


"""
    ValidationSummary

Aggregated statistics produced at the end of a batch simulation + classification run.

# Fields
- `total_simulations::Int`: Number of configurations attempted.
- `valid_count::Int`: Configurations whose polarization curve passed all enabled criteria.
- `invalid_count::Int`: Configurations that failed at least one criterion.
- `failed_count::Int`: Configurations where the ODE solver diverged or errored.
- `criteria_breakdown::Dict{Symbol,Int}`: Per-criterion failure count
  (keys: `:start_voltage`, `:monotonicity`, `:positive_voltage`).
- `execution_time::Float64`: Total wall-clock time of the batch run (seconds).
- `timestamp::String`: ISO-8601 datetime string of when the run completed.
"""
struct ValidationSummary
    total_simulations::Int
    valid_count::Int
    invalid_count::Int
    failed_count::Int
    criteria_breakdown::Dict{Symbol, Int}
    execution_time::Float64
    timestamp::String
end


# ─────────────────────────────────────────────────────────────────────────────
# CSV EXPORT
# ─────────────────────────────────────────────────────────────────────────────

"""
    export_classified_configurations(data, cfg = ExportConfig())::String

Write a `DataFrame` of classified configurations to a CSV file.

Expected columns:
- `sample_id` — unique integer identifier
- One column per sampled parameter (e.g. `Hacl`, `Re`, …)
- `classification` — `"valid"`, `"invalid"`, or `"failed"`
- Optional: `validation_details`, `start_in_range`, `is_monotonic`,
  `has_positive_voltages`, `error_message`

# Returns
The absolute path of the created CSV file.
"""
function export_classified_configurations(data::DataFrame,
                                           cfg::ExportConfig = ExportConfig())::String
    mkpath(cfg.output_dir)

    # Build file name.
    ts      = cfg.timestamp ? "_" * Dates.format(Dates.now(), "yyyymmdd_HHMMSS") : ""
    prefix  = isempty(cfg.prefix) ? "" : cfg.prefix * "_"
    fname   = prefix * "classified_configurations" * ts * ".csv"
    fpath   = abspath(joinpath(cfg.output_dir, fname))

    # Guard against accidental overwrite.
    if isfile(fpath) && !cfg.overwrite
        @warn "File already exists and overwrite=false; skipping: $fpath"
        return fpath
    end

    CSV.write(fpath, data)
    return fpath
end


# ─────────────────────────────────────────────────────────────────────────────
# YAML EXPORT
# ─────────────────────────────────────────────────────────────────────────────

"""
    export_parameter_bounds(bounds, filepath; method = :PRIM, metadata = Dict())::Nothing

Write parameter bounds to a YAML file.

Generated structure:

```yaml
metadata:
  method: PRIM
  timestamp: 2026-05-22T14:30:00
  fuel_cell_type: ZSW_GenStack
parameters:
  Hacl: {min: 5.0e-6, max: 1.2e-5, unit: m}
  Re:   {min: 1.0e-7, max: 3.0e-6, unit: "Ω·m²"}
  # ...
```

This format is compatible with downstream YAML readers and can be passed directly to
future sampling / calibration workflows.
"""
function export_parameter_bounds(bounds::Dict{Symbol, Tuple{Float64, Float64}},
                                 filepath::String;
                                 method::Symbol = :PRIM,
                                 metadata::Dict = Dict())::Nothing
    mkpath(dirname(filepath))

    # Build the YAML structure metadata (used only for comments below)
    ts = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    meta_block = merge(
        Dict("method" => string(method), "timestamp" => ts),
        Dict(string(k) => v for (k, v) in metadata)
    )

     # Write YAML manually for deterministic key ordering and consistent float format
     open(filepath, "w") do io
         println(io, "metadata:")
         for k in sort(collect(keys(meta_block)))
             v = meta_block[k]
             println(io, "  $(k): $(v)")
         end
        println(io, "parameters:")
        for name in sort(collect(keys(bounds)))
            lo, hi = bounds[Symbol(name)]
            println(io, "  $(name):")
            println(io, "    min: $(lo)")
            println(io, "    max: $(hi)")
        end
    end
    return nothing
end


# ─────────────────────────────────────────────────────────────────────────────
# TEXT REPORTS
# ─────────────────────────────────────────────────────────────────────────────

"""
    export_prim_report(original_bounds, restricted_bounds, filepath;
                       prim_metrics = Dict())::Nothing

Write a human-readable report comparing original and PRIM-restricted bounds.

The report includes:
- Side-by-side table: original vs. restricted `[min, max]` for each parameter
- Relative interval shrinkage (%) per parameter
- PRIM precision and coverage
- Random Forest AUC, precision, recall
"""
function export_prim_report(original_bounds::Dict{Symbol, Tuple{Float64, Float64}},
                            restricted_bounds::Dict{Symbol, Tuple{Float64, Float64}},
                            filepath::String;
                            prim_metrics::Dict = Dict())::Nothing
    mkpath(dirname(filepath))
    table = generate_comparison_table(original_bounds, restricted_bounds)
    open(filepath, "w") do io
        println(io, "="^72)
        println(io, "  AlphaPEM — PRIM Valid-Region Restriction Report")
        println(io, "="^72)
        println(io, "  Generated : ", Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))
         if !isempty(prim_metrics)
             println(io, "-"^72)
             println(io, "  PRIM metrics")
             println(io, "-"^72)
             for k in sort(collect(keys(prim_metrics)))
                 v = prim_metrics[k]
                 @printf(io, "    %-22s : %g\n", string(k), v)
             end
         end
        println(io, "-"^72)
        println(io, table)
        println(io, "="^72)
    end
    return nothing
end


"""
    export_validation_summary(summary::ValidationSummary, filepath::String)::Nothing

Write a plain-text summary of a batch classification run.

Reported items: total / valid / invalid / failed counts, percentages,
per-criterion failure counts, total execution time, and run timestamp.
"""
function export_validation_summary(summary::ValidationSummary,
                                   filepath::String)::Nothing
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        println(io, "="^60)
        println(io, "  AlphaPEM — Batch Validity Classification Summary")
        println(io, "="^60)
        println(io, "  Timestamp  : ", summary.timestamp)
        @printf(io, "  Elapsed    : %.2f s\n", summary.execution_time)
        println(io, "-"^60)
        @printf(io, "  Total     : %d\n", summary.total_simulations)
        @printf(io, "  Valid     : %d  (%.1f %%)\n",
                summary.valid_count,
                100.0 * summary.valid_count / max(1, summary.total_simulations))
        @printf(io, "  Invalid   : %d  (%.1f %%)\n",
                summary.invalid_count,
                100.0 * summary.invalid_count / max(1, summary.total_simulations))
        @printf(io, "  Failed    : %d  (%.1f %%)\n",
                summary.failed_count,
                100.0 * summary.failed_count / max(1, summary.total_simulations))
         println(io, "-"^60)
         println(io, "  Per-criterion failures (among invalid/failed):")
         for k in sort(collect(keys(summary.criteria_breakdown)))
             v = summary.criteria_breakdown[k]
             @printf(io, "    %-22s : %d\n", string(k), v)
         end
        println(io, "="^60)
    end
    return nothing
end


# ─────────────────────────────────────────────────────────────────────────────
# UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

"""
    generate_comparison_table(original_bounds, restricted_bounds)::String

Return a formatted ASCII table comparing original and PRIM-restricted bounds.

Example output:
```
┌──────────────┬────────────┬────────────┬────────────┬────────────┬──────────┐
│  Parameter   │  Orig.Min  │  Orig.Max  │  Rest.Min  │  Rest.Max  │  Shrink  │
├──────────────┼────────────┼────────────┼────────────┼────────────┼──────────┤
│     Hacl     │  5.00e-06  │  1.50e-05  │  7.00e-06  │  1.20e-05  │  44.4 %  │
│     Hmem     │  5.00e-06  │  3.00e-05  │  1.00e-05  │  2.50e-05  │  25.0 %  │
└──────────────┴────────────┴────────────┴────────────┴────────────┴──────────┘
```
"""
function generate_comparison_table(original_bounds::Dict{Symbol, Tuple{Float64, Float64}},
                                   restricted_bounds::Dict{Symbol, Tuple{Float64, Float64}})::String
    params = sort(collect(keys(original_bounds)))

    # Column headers
    H = ("Parameter", "Orig Min", "Orig Max", "Rest Min", "Rest Max", "Shrink %")

    rows = Tuple[]
    for name in params
        orig_lo, orig_hi = get(original_bounds,    name, (NaN, NaN))
        rest_lo, rest_hi = get(restricted_bounds,  name, (NaN, NaN))

        orig_span = orig_hi - orig_lo
        rest_span = rest_hi - rest_lo
        shrink = if orig_span > 0.0
            round((1.0 - rest_span / orig_span) * 100.0; digits = 1)
        else
            0.0
        end

        fmt(x) = isnan(x) ? "N/A" : @sprintf("%.3e", x)
        push!(rows, (string(name),
                     fmt(orig_lo), fmt(orig_hi),
                     fmt(rest_lo), fmt(rest_hi),
                     @sprintf("%.1f %%", shrink)))
    end

    # Compute column widths
    widths = [max(length(H[j]), maximum(length(r[j]) for r in rows; init = 0))
              for j in 1:6]

    make_sep(l, m, r) =
        l * join([repeat("─", w + 2) for w in widths], m) * r

    fmt_row(cells) =
        "│" * join((" " * rpad(cells[j], widths[j]) * " │") for j in 1:6)

    buf = IOBuffer()
    println(buf, make_sep("┌", "┬", "┐"))
    println(buf, fmt_row(H))
    println(buf, make_sep("├", "┼", "┤"))
    for r in rows
        println(buf, fmt_row(r))
    end
    println(buf, make_sep("└", "┴", "┘"))

    return String(take!(buf))
end


# ─────────────────────────────────────────────────────────────────────────────
# YAML LOADER (downstream helper)
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_restricted_bounds(yaml_path::String)
        -> Dict{Symbol, Tuple{Float64, Float64}}

Read a YAML file produced by `export_parameter_bounds` and return a
`Dict` mapping each parameter name to its `(min, max)` interval.

This helper is the inverse of `export_parameter_bounds` and is intended to
inject PRIM-restricted bounds into a subsequent calibration workflow:

```julia
bounds = load_restricted_bounds("results/model_validity/restricted_bounds_PRIM.yaml")
# → Dict(:Hacl => (7.0e-6, 1.2e-5), :Re => (2.0e-7, 1.8e-6), ...)
```

The file must follow the format written by `export_parameter_bounds`:
```yaml
metadata:
  method: PRIM
  ...
parameters:
  Hacl:
    min: 7.0e-6
    max: 1.2e-5
  Re:
    min: 2.0e-7
    max: 1.8e-6
```
"""
function load_restricted_bounds(yaml_path::String)::Dict{Symbol, Tuple{Float64, Float64}}
    isfile(yaml_path) ||
        error("Bounds YAML file not found: $yaml_path")

    data   = YAML.load_file(yaml_path)
    params = get(data, "parameters", nothing)
    params === nothing &&
        error("Unexpected YAML format — missing 'parameters' key in: $yaml_path")

    bounds = Dict{Symbol, Tuple{Float64, Float64}}()
    for (name, val) in params
        lo = Float64(val["min"])
        hi = Float64(val["max"])
        bounds[Symbol(name)] = (lo, hi)
    end
    return bounds
end

end # module ResultsExport


