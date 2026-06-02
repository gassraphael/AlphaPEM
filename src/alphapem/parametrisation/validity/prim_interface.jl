# -*- coding: utf-8 -*-

"""
    PRIMInterface

Julia interface to the R `irdpackage` for PRIM-based valid-region analysis.

This module wraps calls to the IRD (Interpretable Regional Descriptors) R package via
`RCall.jl`.  The package provides the PRIM (Patient Rule Induction Method) and MaxBox
algorithms, which find axis-aligned hyperboxes in the parameter space where the
AlphaPEM model is most likely to produce valid polarization curves.

## Pipeline overview

1. A Random Forest is trained on the classified configurations to predict P(valid).
2. PRIM (and optionally MaxBox) peels and pastes the box to maximise precision while
   keeping the reference configuration inside the box.
3. Bounds are post-processed to tighten them without reducing the validity rate.
4. Results are saved as YAML files in `output_dir`.

## R package

The IRD package must be cloned locally before use:

```bash
git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/
```

The relevant sub-directory is `external/supplementary_2023_ird/irdpackage`.
It is loaded via `devtools::load_all()` through `RCall.jl`.

# Exports
- `PRIMConfig`: All options for one PRIM analysis run
- `PRIMResult`: Structured output of one PRIM/MaxBox run
- `run_prim_analysis`: Execute the full PRIM pipeline and return results
"""
module PRIMInterface

using RCall
using CSV
using DataFrames
using Dates
using Printf
using YAML

export PRIMConfig,
       PRIMResult,
       run_prim_analysis

# ─────────────────────────────────────────────────────────────────────────────
# DATA STRUCTURES
# ─────────────────────────────────────────────────────────────────────────────

"""
    PRIMConfig

All options needed to run one PRIM analysis via the R `irdpackage`.

# Fields
- `data_path::String`: Path to the CSV of classified configurations.
  **Leave empty** when calling via `find_valid_region` — it is set automatically.
- `ird_package_dir::String`: Path to the local `irdpackage` folder.
  Default: `"external/supplementary_2023_ird/irdpackage"`.
- `reference_config_path::String`: Path to a YAML file with the reference (known-valid)
  configuration used as *x_interest* in PRIM.
  **Leave empty** when calling via `find_valid_region` — it is set automatically.
- `probability_range::Tuple{Float64,Float64}`: Target P(valid) interval for the box.
  Default `(0.8, 1.0)` — keep only regions estimated ≥ 80 % valid.
- `methods::Vector{Symbol}`: IRD methods to run. Supported: `:PRIM`, `:MaxBox`.
- `categorical_params::Vector{Symbol}`: Parameters to treat as categorical (e.g., `[:e]`).
- `seed::Int`: Random seed for the Random Forest and PRIM. Default `42`.
- `output_dir::String`: Output directory. **Leave empty** when calling via
  `find_valid_region` — it is fixed to `VALIDITY_OUTPUT_DIR` automatically.
- `target_column::String`: Name of the classification column. Default `"classification"`.
- `positive_class::String`: Label of the valid class. Default `"valid"`.

# Example
```julia
cfg = PRIMConfig(
    data_path             = "results/model_validity/classified_configs.csv",
    reference_config_path = "results/model_validity/reference_ZSW.yaml",
    output_dir            = "results/model_validity",
)
```
"""
Base.@kwdef struct PRIMConfig
    data_path::String                           = ""   # set internally by find_valid_region
    ird_package_dir::String                     = "external/supplementary_2023_ird/irdpackage"
    reference_config_path::String               = ""   # set internally by find_valid_region
    probability_range::Tuple{Float64, Float64}  = (0.8, 1.0)
    methods::Vector{Symbol}                     = [:PRIM, :MaxBox]
    categorical_params::Vector{Symbol}          = [:e]
    seed::Int                                   = 42
    output_dir::String                          = ""   # set internally by find_valid_region
    target_column::String                       = "classification"
    positive_class::String                      = "valid"
end


"""
    PRIMResult

Structured output of one completed PRIM or MaxBox run.

# Fields
- `method::Symbol`: Which algorithm was used (`:PRIM` or `:MaxBox`).
- `restricted_bounds::Dict{Symbol,Tuple{Float64,Float64}}`: PRIM-derived bounds per parameter.
- `original_bounds::Dict{Symbol,Tuple{Float64,Float64}}`: Original (pre-PRIM) bounds.
- `precision::Float64`: Share of valid configurations inside the found box (0–1).
- `coverage::Float64`: Share of all configurations retained by the box (0–1).
- `n_inside_box::Int`: Total number of configurations inside the box.
- `n_valid_inside::Int`: Valid configurations inside the box.
- `rf_metrics::Dict{String,Float64}`: Random Forest performance (keys: `"AUC"`,
  `"precision"`, `"recall"`).
- `output_files::Dict{Symbol,String}`: Paths to generated files
  (keys: `:bounds_yaml`, `:report_txt`).
"""
struct PRIMResult
    method::Symbol
    restricted_bounds::Dict{Symbol, Tuple{Float64, Float64}}
    original_bounds::Dict{Symbol, Tuple{Float64, Float64}}
    precision::Float64
    coverage::Float64
    n_inside_box::Int
    n_valid_inside::Int
    rf_metrics::Dict{String, Float64}
    output_files::Dict{Symbol, String}
end


# ─────────────────────────────────────────────────────────────────────────────
# PRIVATE HELPERS
# ─────────────────────────────────────────────────────────────────────────────

"""
    _filter_r_stderr(r_stderr)::String

Filter R stderr to remove informational messages and keep only actual errors/warnings.
Removes lines related to:
- Package loading and patching
- PRIM algorithm trace messages (peeling/pasting iterations)
- Informational messages from IRD helpers
- Final completion messages
"""
function _filter_r_stderr(r_stderr::String)::String
    # List of patterns for informational messages that should be filtered
    informational_patterns = [
        r"^\[irdpackage\]",                    # Package loading messages
        r"^peeling iteration",                 # PRIM peeling iterations
        r"^pasting iteration",                 # PRIM pasting iterations
        r"^Ignoring x_interest keys",          # Unused keys info
        r"^The `desired_class`",               # Class assignment confirmation
        r"^Saved RF metrics:",                 # Saved files info
        r"^run_prim\.R — done\.",              # Completion message
    ]

    lines = split(r_stderr, "\n")
    filtered_lines = String[]

    for line in lines
        # Keep empty lines and lines that are NOT in the informational patterns
        is_informational = any(
            occursin(pattern, line)
            for pattern in informational_patterns
        )
        if !is_informational && !isempty(strip(line))
            push!(filtered_lines, line)
        end
    end

    return join(filtered_lines, "\n")
end


"""
    _find_rscript()::String

Return the absolute path to the `Rscript` executable used by RCall.jl.
Falls back to `Sys.which("Rscript")` if the RCall query fails.
"""
function _find_rscript()::String
    # Primary: ask the embedded R session where its executables live
    rscript = try
        r_bin = rcopy(String, R"R.home('bin')")
        joinpath(r_bin, "Rscript")
    catch
        ""
    end

    if !isempty(rscript) && isfile(rscript)
        return rscript
    end

    # Fallback: search PATH
    rscript_which = Sys.which("Rscript")
    isempty(rscript_which) &&
        error("Rscript not found. Ensure R is installed and in PATH, " *
              "or that Julia's RCall.jl is properly configured.")
    return rscript_which
end


"""
    _parse_bounds_yaml(yaml_path)
        -> (bounds::Dict{Symbol,Tuple{Float64,Float64}},
            cat_values::Dict{Symbol,Vector{Float64}})

Parse a YAML produced by `write_ird_yaml()` in `ird_helpers.R`.

Expected format (list of parameter entries):
```yaml
parameters:
  - name: Hacl
    type: continuous
    low: 5.0e-6
    high: 1.5e-5
    fixed: False
  - name: e
    type: categorical
    values: [3, 4, 5]
    fixed: False
```

Returns:
- `bounds`    — interval `(min, max)` for every parameter.
  Categorical params use `(minimum(values), maximum(values))`.
- `cat_values` — for categorical params only: sorted Float64 vector of allowed values.
"""
function _parse_bounds_yaml(yaml_path::String)
    data = YAML.load_file(yaml_path)

    haskey(data, "parameters") ||
        error("Unexpected bounds YAML format — missing 'parameters' key in: $yaml_path")

    bounds     = Dict{Symbol, Tuple{Float64, Float64}}()
    cat_values = Dict{Symbol, Vector{Float64}}()

    for p in data["parameters"]
        name = Symbol(p["name"])
        ptype = get(p, "type", "continuous")

        if ptype == "categorical"
            vals = sort(Float64[Float64(v) for v in p["values"]])
            cat_values[name] = vals
            bounds[name]     = (vals[1], vals[end])
        else
            lo = Float64(p["low"])
            hi = Float64(p["high"])
            bounds[name] = (lo, hi)
        end
    end

    return bounds, cat_values
end


"""
    _parse_rf_metrics(rf_txt_path)::Dict{String,Float64}

Parse the `RF_metrics_<run_name>.txt` file saved by `run_prim.R`.

Extracts `classif.auc`, `classif.precision`, `classif.recall` from the
"[Test set]" block (first table encountered). Returns zeros for any metric
that cannot be parsed.
"""
function _parse_rf_metrics(rf_txt_path::String)::Dict{String, Float64}
    result = Dict{String, Float64}(
        "AUC"       => 0.0,
        "precision" => 0.0,
        "recall"    => 0.0,
    )
    isfile(rf_txt_path) || return result

    lines = readlines(rf_txt_path)
    # Find the header line after "[Test set]" and parse the following value line
    in_test_block = false
    header_line   = ""
    for (i, line) in enumerate(lines)
        stripped = strip(line)
        if stripped == "[Test set]"
            in_test_block = true
            continue
        end
        if in_test_block && startswith(stripped, "classif.")
            header_line = stripped
            # The next non-empty line contains the values
            for j in (i + 1):min(i + 5, length(lines))
                val_line = strip(lines[j])
                isempty(val_line) && continue
                headers = split(header_line)
                values  = split(val_line)
                length(headers) == length(values) || break
                for (h, v) in zip(headers, values)
                    parsed = tryparse(Float64, v)
                    parsed === nothing && continue
                    if h == "classif.auc"       result["AUC"]       = parsed end
                    if h == "classif.precision" result["precision"] = parsed end
                    if h == "classif.recall"    result["recall"]    = parsed end
                end
                break
            end
            break   # Only use the [Test set] block
        end
        # Stop if we leave the test block
        if in_test_block && stripped == "[Full dataset]"
            break
        end
    end
    return result
end


"""
    _compute_box_statistics(data_path, bounds, cat_values, positive_class)
        -> (precision, coverage, n_inside, n_valid_inside)

Read the classified configurations CSV and evaluate how many rows fall inside
the box defined by `bounds`.

For continuous parameters, a row is inside if `lo ≤ val ≤ hi`.
For categorical parameters (listed in `cat_values`), a row is inside if its
rounded-integer value is within the allowed set.
"""
function _compute_box_statistics(data_path::String,
                                  bounds::Dict{Symbol, Tuple{Float64, Float64}},
                                  cat_values::Dict{Symbol, Vector{Float64}},
                                  positive_class::String)
    data = CSV.read(data_path, DataFrame)
    n_total = nrow(data)

    bound_names = collect(keys(bounds))

    function _in_box(row)::Bool
        for name in bound_names
            col = String(name)
            (col in names(data)) || continue   # column not in data → skip
            lo, hi = bounds[name]
            val = Float64(row[col])
            if haskey(cat_values, name)
                # Categorical: accept if rounded value is in the allowed set
                rv = round(val)
                in_set = any(v -> isapprox(rv, v; atol = 0.5), cat_values[name])
                in_set || return false
            else
                (lo <= val <= hi) || return false
            end
        end
        return true
    end

    n_inside       = 0
    n_valid_inside = 0
    for row in eachrow(data)
        if _in_box(row)
            n_inside += 1
            if string(row.classification) == positive_class
                n_valid_inside += 1
            end
        end
    end

    precision = n_inside > 0 ? n_valid_inside / n_inside : 0.0
    coverage  = n_total  > 0 ? n_inside       / n_total  : 0.0
    return precision, coverage, n_inside, n_valid_inside
end


"""
    _original_bounds_from_data(data_path, param_names)
        -> Dict{Symbol,Tuple{Float64,Float64}}

Compute the observed (min, max) range of each parameter column in the classified CSV.
This represents the original parameter bounds before PRIM restriction.
"""
function _original_bounds_from_data(data_path::String,
                                     param_names::Vector{Symbol})::Dict{Symbol, Tuple{Float64, Float64}}
    data   = CSV.read(data_path, DataFrame)
    bounds = Dict{Symbol, Tuple{Float64, Float64}}()
    for name in param_names
        col = String(name)
        if col in names(data)
            vals = skipmissing(data[!, col])
            bounds[name] = (minimum(vals), maximum(vals))
        end
    end
    return bounds
end


# ─────────────────────────────────────────────────────────────────────────────
# MAIN ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_prim_analysis(cfg::PRIMConfig)::Vector{PRIMResult}

Execute the full IRD/PRIM pipeline and return one `PRIMResult` per requested method.

Steps performed:
1. Check that the IRD package directory and required R scripts are present.
2. Locate the bundled `run_prim.R` and `ird_helpers.R` scripts.
3. Call `Rscript run_prim.R` as a subprocess with appropriate arguments.
4. Parse the generated `IRD_bounds_<method>_<run_name>.yaml` files.
5. Parse the `RF_metrics_<run_name>.txt` file.
6. Compute precision/coverage from the classified configurations CSV.
7. Return one `PRIMResult` per method.

# Example
```julia
cfg = PRIMConfig(
    data_path             = "results/model_validity/classified_configs.csv",
    reference_config_path = "results/model_validity/reference_ZSW.yaml",
)
results = run_prim_analysis(cfg)
for r in results
    @info "\$(r.method): precision=\$(round(r.precision; digits=3)), " *
          "coverage=\$(round(r.coverage; digits=3))"
end
```
"""
function run_prim_analysis(cfg::PRIMConfig)::Vector{PRIMResult}
    # ── 1. Validate inputs ────────────────────────────────────────────────────
    isfile(cfg.data_path) ||
        error("Classified data CSV not found: $(cfg.data_path)")
    isfile(cfg.reference_config_path) ||
        error("Reference config YAML not found: $(cfg.reference_config_path)")

    # ── 2. Check IRD package directory ────────────────────────────────────────
    isdir(cfg.ird_package_dir) ||
        error("""
        IRD package directory not found: $(cfg.ird_package_dir)

        Please run the following commands before using run_prim_analysis:

            git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/

        Also ensure that R is installed and available in your PATH, and that the
        required R packages are installed (devtools, mlr3, mlr3learners, mlr3pipelines,
        iml, ranger, yaml, jsonlite, data.table, optparse).
        """)

    # ── 3. Locate bundled R scripts ───────────────────────────────────────────
    r_dir       = joinpath(@__DIR__, "R")
    run_prim_r  = joinpath(r_dir, "run_prim.R")
    helpers_r   = joinpath(r_dir, "ird_helpers.R")

    isfile(run_prim_r) ||
        error("Bundled R script not found: $run_prim_r")
    isfile(helpers_r)  ||
        error("Bundled R helpers not found: $helpers_r")

    # ── 4. Build Rscript command ──────────────────────────────────────────────
    mkpath(cfg.output_dir)

    rscript    = _find_rscript()
    run_name   = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    methods_s  = join([string(m) for m in cfg.methods], ",")
    cat_s      = join([string(p) for p in cfg.categorical_params], ",")
    range_s    = "$(cfg.probability_range[1]),$(cfg.probability_range[2])"

    cmd = Cmd([
        rscript, run_prim_r,
        "--data",                   abspath(cfg.data_path),
        "--target",                 cfg.target_column,
        "--positive",               cfg.positive_class,
        "--xinterest",              abspath(cfg.reference_config_path),
        "--range",                  range_s,
        "--methods",                methods_s,
        "--outdir",                 abspath(cfg.output_dir),
        "--categorical_overrides",  cat_s,
        "--run_name",               run_name,
        "--seed",                   string(cfg.seed),
        "--ird_pkg_dir",            abspath(cfg.ird_package_dir),
        "--helpers_path",           abspath(helpers_r),
    ])

    # ── 5. Run the R script ───────────────────────────────────────────────────
    @info "Running PRIM analysis via Rscript…  (methods: $methods_s)"
    stdout_buf = IOBuffer()
    stderr_buf = IOBuffer()
    try
        run(pipeline(cmd; stdout = stdout_buf, stderr = stderr_buf); wait = true)
    catch e
        r_stdout = String(take!(stdout_buf))
        r_stderr = String(take!(stderr_buf))
        err_msg = sprint(showerror, e)
        error("""
        PRIM Rscript execution failed.

        Command:
          $(join(cmd.exec, " "))

        Julia error:
          $err_msg

        R stdout:
        $r_stdout

        R stderr:
        $r_stderr
        """)
    end

    r_stdout = String(take!(stdout_buf))
    r_stderr = String(take!(stderr_buf))

    if !isempty(r_stdout)
        @info "R stdout:\n$r_stdout"
    end

    # Filter and display only significant stderr messages (actual errors/warnings)
    r_stderr_filtered = _filter_r_stderr(r_stderr)
    if !isempty(r_stderr_filtered)
        @warn "R stderr:\n$r_stderr_filtered"
    end

    # ── 6. Parse RF metrics (shared across all methods in this run) ───────────
    rf_txt      = joinpath(cfg.output_dir, "RF_metrics_$(run_name).txt")
    rf_metrics  = _parse_rf_metrics(rf_txt)

    # ── 7. Build PRIMResult for each method ───────────────────────────────────
    results = PRIMResult[]

    for method in cfg.methods
        method_str = string(method)
        bounds_yaml = joinpath(cfg.output_dir,
                               "IRD_bounds_$(method_str)_$(run_name).yaml")
        report_txt  = joinpath(cfg.output_dir,
                               "IRD_report_$(method_str)_$(run_name).txt")

        if !isfile(bounds_yaml)
            @warn "PRIM output YAML not found for method=$method_str; skipping."
            continue
        end

        # Parse the YAML bounds produced by write_ird_yaml()
        restricted_bounds, cat_values = _parse_bounds_yaml(bounds_yaml)

        # Compute original (pre-PRIM) bounds from the observed data ranges
        original_bounds = _original_bounds_from_data(
            cfg.data_path, collect(keys(restricted_bounds))
        )

        # Compute precision / coverage directly from the classified dataset
        precision, coverage, n_inside, n_valid_inside = _compute_box_statistics(
            cfg.data_path, restricted_bounds, cat_values, cfg.positive_class
        )

        output_files = Dict{Symbol, String}(
            :bounds_yaml => bounds_yaml,
            :report_txt  => report_txt,
            :rf_metrics  => rf_txt,
        )

        push!(results, PRIMResult(
            method,
            restricted_bounds,
            original_bounds,
            precision,
            coverage,
            n_inside,
            n_valid_inside,
            rf_metrics,
            output_files,
        ))

        @info @sprintf(
            "  %s — precision=%.3f, coverage=%.3f  (%d/%d valid inside box)",
            method_str, precision, coverage, n_valid_inside, n_inside
        )
    end

    isempty(results) &&
        @warn "run_prim_analysis: no results were produced. Check R stderr above."

    return results
end

end # module PRIMInterface

