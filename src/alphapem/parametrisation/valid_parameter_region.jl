# -*- coding: utf-8 -*-

"""
    ValidParameterRegion

Main module for the validity-region analysis of AlphaPEM undetermined parameters.

## Purpose

The undetermined-parameter space of AlphaPEM is large and contains configurations
that produce non-physical polarization curves (negative voltages, etc.).  This
module provides a complete pipeline to identify a compact, interpretable sub-region
where simulations are reliable, so that subsequent calibration, sampling, and
surrogate-modelling efforts can be focused there.

## Workflow

When the user calls `run_validity_analysis`, the following steps are executed:

| Step | Action                                                                               |
|------|--------------------------------------------------------------------------------------|
| 1    | Define parameter bounds and draw configurations by Latin Hypercube Sampling          |
| 2    | Simulate each configuration with AlphaPEM and classify the polarization curve        |
| 3    | *(optional)* Restrict the valid region via PRIM/MaxBox (requires R + IRD package)   |
| 4    | Export results: classified CSV, bounds YAML, validation summary, PRIM report         |

`run_validity_analysis` orchestrates all four steps in a single call.
Steps 1–2–4 always run; Step 3 is enabled by passing a `PRIMConfig`.

## Sub-modules

| Sub-module               | Responsibility                                         |
|--------------------------|--------------------------------------------------------|
| `ValidityCriteria`       | Per-curve validity checks (voltage range, monotonicity)|
| `ConfigurationSampling`  | LHS sample generation and parameter-struct mapping     |
| `PRIMInterface`          | Julia ↔ R bridge for `irdpackage` (PRIM, MaxBox)      |
| `ResultsExport`          | CSV / YAML / text-report export                        |

## Credits

This work is inspired by the master's thesis:
  *Sensitivity Analysis and Surrogate Modeling of PEM Fuel Cells*,
  by Nathaly Vergel Serrano, Dejvis Toptani, and Camila Bermudez Valderrama,
  supervised by Dr. Giuseppe Casalicchio and Fiona Ewald,
  in collaboration with Luis Winkler (ZSW Ulm) and Prof. Herbert Palm (UAS Munich).
  Statistical Consulting, Master in Statistics and Data Science,
  Institute of Statistics, Ludwig-Maximilians-Universität München.
  Repository: https://github.com/nathaly-vergel/Official-Sensitivity-Analysis-and-Surrogate-Modeling-of-PEM-Fuel-Cells

The IRD/PRIM framework used downstream relies on:
  https://github.com/slds-lmu/supplementary_2023_ird

## Example

```julia
using AlphaPEM.Parametrisation.ValidParameterRegion

cfg = ValidityAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    n_samples      = 2000,
    voltage_zone   = :full,
    output_dir     = "results/model_validity",
    hyperbox_finder_method = :PRIM,
)

results = run_validity_analysis(cfg)   # ← full pipeline (Steps 1-2-4)
println("Restricted bounds saved to: ", results.output_files[:bounds_yaml])
```
"""
module ValidParameterRegion

# ─── Standard library ────────────────────────────────────────────────────────
using Printf
using Dates
using LinearAlgebra: BLAS

# ─── External packages ───────────────────────────────────────────────────────
using DataFrames
using CSV
using ProgressMeter: Progress, next!, finish!

# ─── AlphaPEM simulation components ──────────────────────────────────────────
using AlphaPEM.Config:  SimulationConfig, PolarizationParams, NumericalParams
using AlphaPEM.Fuelcell: create_fuelcell
using AlphaPEM.Currents: create_current
# Import the AlphaPEM simulator struct under an alias to avoid name collision
# with the parent module (AlphaPEM the package).
import AlphaPEM.Core.Models: AlphaPEM as AlphaPEMSimulator, simulate_model!

# ─── Fixed output directory (not user-configurable) ──────────────────────────
"""
    VALIDITY_OUTPUT_DIR

Absolute path to the fixed output directory for all validity-analysis results.
Resolves to `<project_root>/results/model_validity` regardless of the working directory.
"""
const VALIDITY_OUTPUT_DIR = abspath(joinpath(@__DIR__, "..", "..", "..", "results", "model_validity"))

# ─── Sub-modules ─────────────────────────────────────────────────────────────
include(joinpath(@__DIR__, "validity/validity_criteria.jl"))
include(joinpath(@__DIR__, "validity/configuration_sampling.jl"))
include(joinpath(@__DIR__, "validity/prim_interface.jl"))
include(joinpath(@__DIR__, "validity/results_export.jl"))

using .ValidityCriteria
using .ConfigurationSampling
using .PRIMInterface
using .ResultsExport

# Re-export sub-module types into this namespace for convenient access
using .ValidityCriteria:      ValidityCriteriaConfig, ValidationResult,
                               classify_polarization_curve
using .ConfigurationSampling: ParameterBound, ParameterBounds, SamplingConfig,
                               bounds_for_fuel_cell, generate_lhs_samples,
                               apply_bounds_to_params, get_reference_config
using .PRIMInterface:         PRIMConfig, PRIMResult,
                               run_prim_analysis
using .ResultsExport:         ExportConfig, ValidationSummary,
                               export_classified_configurations,
                               export_parameter_bounds,
                               export_prim_report,
                               export_validation_summary,
                               generate_comparison_table,
                               load_restricted_bounds

# Public API
export ValidityCriteria,
       ConfigurationSampling,
       PRIMInterface,
       ResultsExport,
       # Fixed output path
       VALIDITY_OUTPUT_DIR,
       # Top-level workflow
       ValidityAnalysisConfig,
       ValidityAnalysisResult,
       run_validity_analysis,
       generate_test_samples,
       classify_batch_simulations,
       find_valid_region,
       # ValidityCriteria types & functions
       ValidityCriteriaConfig,
       ValidationResult,
       classify_polarization_curve,
       check_start_voltage_range,
       check_monotonicity,
       check_positive_voltages,
       # ConfigurationSampling types & functions
       ParameterBound,
       ParameterBounds,
       SamplingConfig,
       bounds_for_fuel_cell,
       generate_lhs_samples,
       apply_bounds_to_params,
       get_reference_config,
       # PRIMInterface types & functions
       PRIMConfig,
       PRIMResult,
       run_prim_analysis,
       # ResultsExport types & functions
       ExportConfig,
       ValidationSummary,
       export_classified_configurations,
       export_parameter_bounds,
       export_prim_report,
       export_validation_summary,
       generate_comparison_table,
       load_restricted_bounds

# ─────────────────────────────────────────────────────────────────────────────
# TOP-LEVEL CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

"""
    ValidityAnalysisConfig

Master configuration for the complete validity-analysis pipeline.

Consolidates options from all sub-modules so that a full run can be
launched with a single struct.

# Fields
- `fuel_cell_type::Symbol`: Fuel-cell type to analyse. Default `:ZSW_GenStack`.
- `voltage_zone::Symbol`: `:before_voltage_drop` or `:full`. Default `:full`.
- `n_samples::Int`: Number of LHS configurations to simulate. Default `10000`.
- `sampling_seed::Int`: Random seed for reproducible LHS. Default `42`.
- `validation_criteria::ValidityCriteriaConfig`: Thresholds for the curve classifiers.
- `polarization_params::PolarizationParams`: Current profile used for batch simulations.
  Defaults to fast-batch settings (short stabilisation times) optimised for validity
  testing rather than high-fidelity calibration.
- `nb_gc::Int`: Number of gas-channel nodes used in batch runs. Default `1` (fastest).
- `parallel::Bool`: Use multi-threading (`Threads.@threads`) for batch simulation.
  Default `true`.
- `checkpoint_interval::Int`: Save intermediate results every N simulations.
  Set to `0` to disable. Default `100`.
- `save_curves::Bool`: Save polarization curves to CSV. Default `true`.
- `reuse_from::Union{String, Nothing}": Path to previous run directory to reuse curves.
- `hyperbox_finder_method::Union{Symbol, Nothing}`: Hyperbox-finder method to run.
  Supported values: `:PRIM` (default) and `:MaxBox`.
  Set to `nothing` to skip the hyperbox-finder step entirely.
  Default `:PRIM`.

# Example
```julia
cfg = ValidityAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    n_samples      = 500,       # quick test
    parallel       = false,
    hyperbox_finder_method   = :PRIM,    # or `nothing` to skip, or `:MaxBox` to select MaxBox
)
```
"""
Base.@kwdef struct ValidityAnalysisConfig
    fuel_cell_type::Symbol              = :ZSW_GenStack
    voltage_zone::Symbol                = :full
    n_samples::Int                      = 2000
    sampling_seed::Int                  = 42
    validation_criteria::ValidityCriteriaConfig = ValidityCriteriaConfig()
    polarization_params::PolarizationParams = PolarizationParams()
    nb_gc::Int                          = 1   # Minimum spatial resolution for speed
    parallel::Bool                      = true
    save_curves::Bool                   = true   # Save polarization curves to curves.csv
    reuse_from::Union{String, Nothing}  = nothing  # Path to previous run directory to reuse curves
    hyperbox_finder_method::Union{Symbol, Nothing} = :PRIM # e.g. :PRIM or :MaxBox; set to `nothing` to skip
end


"""
    ValidityAnalysisResult

Complete output of a validity-analysis run.

# Fields
- `config::ValidityAnalysisConfig`: The configuration used for this run.
- `original_bounds::Dict{Symbol,Tuple{Float64,Float64}}`: Prior parameter bounds.
- `restricted_bounds::Dict{Symbol,Tuple{Float64,Float64}}`: PRIM-restricted bounds.
- `validation_summary::ValidationSummary`: Classification statistics.
- `prim_results::Vector{PRIMResult}`: One entry per IRD method that was run.
- `output_files::Dict{Symbol,String}`: Paths to all generated files.
- `execution_time::Float64`: Total wall-clock time in seconds.
"""
struct ValidityAnalysisResult
    config::ValidityAnalysisConfig
    original_bounds::Dict{Symbol, Tuple{Float64, Float64}}
    restricted_bounds::Dict{Symbol, Tuple{Float64, Float64}}
    validation_summary::ValidationSummary
    prim_results::Vector{PRIMResult}
    output_files::Dict{Symbol, String}
    execution_time::Float64
end


# ─────────────────────────────────────────────────────────────────────────────
# ORCHESTRATION FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

"""
    _generate_run_directory(cfg::ValidityAnalysisConfig) -> String

Generate a timestamped and auto-indexed output directory for the current run.

Directory structure: `<VALIDITY_OUTPUT_DIR>/<YYYYMMDD>_<NNN>_s<NSAMPLES>/`
where:
  - YYYYMMDD = current date
  - NNN = auto-incremented index if multiple runs on the same day
  - NSAMPLES = number of samples in `cfg`

All output files for this run will be written to the returned directory.

# Example
```julia
dir = _generate_run_directory(cfg)  # → "results/model_validity/20260601_001_s0020/"
```
"""
function _generate_run_directory(cfg::ValidityAnalysisConfig)::String
    date_str = Dates.format(Dates.today(), "yyyymmdd")

    # Find the next available index for this date
    base_dir = VALIDITY_OUTPUT_DIR
    mkpath(base_dir)

    index = 1
    while isdir(joinpath(base_dir, @sprintf("%s_%03d_s%05d", date_str, index, cfg.n_samples)))
        index += 1
    end

    run_dir = joinpath(base_dir, @sprintf("%s_%03d_s%05d", date_str, index, cfg.n_samples))
    mkpath(run_dir)
    return run_dir
end


"""
    run_validity_analysis(cfg::ValidityAnalysisConfig,
                          prim_cfg::Union{PRIMConfig,Nothing} = nothing)
        ::ValidityAnalysisResult

Execute the complete validity-analysis pipeline (Steps 1–4) and return all results.

## Pipeline

| Step | Action                                                                                  |
|------|-----------------------------------------------------------------------------------------|
| 1    | Draw `cfg.n_samples` configurations by Latin Hypercube Sampling                         |
| 2    | Simulate each configuration with AlphaPEM and classify the polarization curve           |
| 3    | *(optional)* Restrict bounds via PRIM/MaxBox if `prim_cfg` is provided                 |
| 4    | Export results: classified CSV, curves CSV, bounds YAML, validation summary, PRIM report|

## Output Organization

All results are automatically organised in a timestamped directory:
```
results/model_validity/<YYYYMMDD>_<NNN>_s<NSAMPLES>/
├── summary.txt               (validation summary)
├── configurations.csv        (classified parameter combinations)
├── curves.csv               (polarization curves, if save_curves=true)
├── bounds_prior.yaml        (original parameter bounds)
├── bounds_restricted_*.yaml (PRIM-restricted bounds, if PRIM is run)
└── ...
```

## Arguments
- `cfg::ValidityAnalysisConfig` — master configuration (sampling, simulation, criteria, …).
  - Set `save_curves=true` to save U(I) curves for each configuration.
  - Set `reuse_from="/path/to/previous/run"` to skip simulation and reuse curves.
- `prim_cfg::Union{PRIMConfig,Nothing}` — PRIM options.  Pass `nothing` (default) to skip
  PRIM analysis.

## Returns
A `ValidityAnalysisResult` with all intermediate outputs and paths to generated files.

## Reuse curves from a previous run with different classification criteria
cfg = ValidityAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    n_samples      = 500,
    reuse_from     = "results/model_validity/20260601_001_s0500/",
    validation_criteria = ValidityCriteriaConfig(monotonic_threshold=0.01),
)
result = run_validity_analysis(cfg)
```
"""
function run_validity_analysis(cfg::ValidityAnalysisConfig,
                                prim_cfg::Union{PRIMConfig, Nothing} = nothing
                                )::ValidityAnalysisResult
    overall_start = time()

    # ── Create or reuse run directory ─────────────────────────────────────────
    run_dir = if cfg.reuse_from !== nothing
        # Reuse existing run directory
        cfg.reuse_from
    else
        # Generate new timestamped directory
        _generate_run_directory(cfg)
    end

    mkpath(run_dir)
    output_files = Dict{Symbol, String}()

    @info "Run directory: $run_dir"

    # ── STEP 1: LHS sampling ──────────────────────────────────────────────────
    @info "STEP 1 — Generating $(cfg.n_samples) LHS samples…"
    X, pb = generate_test_samples(cfg)
    @info @sprintf("  → Sample matrix: %d × %d", size(X, 1), size(X, 2))

    # Extract original bounds from ParameterBounds
    orig_bounds = Dict{Symbol, Tuple{Float64, Float64}}(
        b.name => (b.min, b.max) for b in pb.bounds
    )

    # Export original bounds to YAML
    bounds_path = abspath(joinpath(run_dir, "bounds_prior.yaml"))
    export_parameter_bounds(orig_bounds, bounds_path;
                            method   = :prior,
                            metadata = Dict("fuel_cell_type" => string(cfg.fuel_cell_type),
                                            "voltage_zone"   => string(cfg.voltage_zone)))
    output_files[:bounds_prior_yaml] = bounds_path
    @info "  → Prior bounds: $bounds_path"

    # ── STEP 2: Batch simulation + classification ─────────────────────────────
    if cfg.reuse_from !== nothing
        @info "STEP 2 — Loading cached curves from $(cfg.reuse_from)…"
        try
            curves_path = joinpath(cfg.reuse_from, "curves.csv")
            data = classify_batch_simulations(X, pb, cfg, run_dir; curves_df = _load_curves_from_csv(curves_path))
            output_files[:curves_csv] = curves_path
            @info "  → Using curves from: $curves_path"
        catch e
            @warn "Failed to load cached curves: $e. Running fresh simulation…"
            result_tuple = classify_batch_simulations(X, pb, cfg, run_dir)
            if result_tuple isa NamedTuple
                data, summary = result_tuple.data, result_tuple.summary
            else
                data = result_tuple
                summary = nothing
            end
        end
    else
        @info "STEP 2 — Batch simulation and classification…"
        result_tuple = classify_batch_simulations(X, pb, cfg, run_dir)
        if result_tuple isa NamedTuple
            data, summary = result_tuple.data, result_tuple.summary
        else
            data = result_tuple
            summary = nothing
        end
    end

    # Export classified CSV
    csv_path = abspath(joinpath(run_dir, "configurations.csv"))
    CSV.write(csv_path, data)
    output_files[:configurations_csv] = csv_path
    @info "  → Configurations: $csv_path"

    # Extract summary stats from data
    classifications = data.classification
    n_valid   = count(==("valid"),   classifications)
    n_invalid = count(==("invalid"), classifications)
    n_failed  = count(==("failed"),  classifications)

    criteria_breakdown = Dict{Symbol, Int}(
        :start_voltage    => count(x -> x === false, data.start_in_range),
        :monotonicity     => count(x -> x === false, data.is_monotonic),
        :positive_voltage => count(x -> x === false, data.has_positive_voltages),
    )

    summary = ValidationSummary(
        cfg.n_samples,
        n_valid,
        n_invalid,
        n_failed,
        criteria_breakdown,
        0.0,  # elapsed time, already computed in classify_batch_simulations
        Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    )

    # Export validation summary
    summary_path = abspath(joinpath(run_dir, "summary.txt"))
    export_validation_summary(summary, summary_path)
    output_files[:summary] = summary_path
    @info "  → Summary: $summary_path"

    # ── STEP 3: PRIM (optional) ───────────────────────────────────────────────
    prim_results       = PRIMResult[]
    restricted_bounds  = orig_bounds          # default: unchanged when PRIM is skipped

    if prim_cfg !== nothing
        @info "STEP 3 — Running PRIM analysis…"
        ref_config   = get_reference_config(cfg.fuel_cell_type)
        prim_results = find_valid_region(data, ref_config, prim_cfg, run_dir)

        for r in prim_results
            mstr = string(r.method)

            # Export restricted bounds YAML
            bounds_path = abspath(joinpath(run_dir, "bounds_restricted_$(mstr).yaml"))
            export_parameter_bounds(r.restricted_bounds, bounds_path;
                                    method   = r.method,
                                    metadata = Dict(
                                        "fuel_cell_type" => string(cfg.fuel_cell_type),
                                        "precision"      => r.precision,
                                        "coverage"       => r.coverage,
                                    ))
            output_files[Symbol("bounds_restricted_$(mstr)_yaml")] = bounds_path
            @info "  → Restricted bounds ($(mstr)): $bounds_path"

            # Export PRIM comparison report
            report_path = abspath(joinpath(run_dir, "report_prim_$(mstr).txt"))
            prim_metrics = merge(
                r.rf_metrics,
                Dict("box_precision" => r.precision,
                     "box_coverage"  => r.coverage,
                     "n_inside_box"  => Float64(r.n_inside_box),
                     "n_valid_inside"=> Float64(r.n_valid_inside))
            )
            export_prim_report(orig_bounds, r.restricted_bounds, report_path;
                               prim_metrics = prim_metrics)
            output_files[Symbol("report_prim_$(mstr)")] = report_path
            @info "  → PRIM report ($(mstr)): $report_path"
        end

        # Use the first successful method's bounds as the canonical restricted result
        if !isempty(prim_results)
            restricted_bounds = first(prim_results).restricted_bounds
        end
    else
        @info "STEP 3 — Skipped (pass a PRIMConfig to enable PRIM analysis)."
    end

    overall_elapsed = time() - overall_start
    @info @sprintf("Pipeline complete in %.1f s.", overall_elapsed)

    return ValidityAnalysisResult(
        cfg,
        orig_bounds,
        restricted_bounds,
        summary,
        prim_results,
        output_files,
        overall_elapsed,
    )
end


"""
    generate_test_samples(cfg::ValidityAnalysisConfig)

Step 1 of the pipeline: define parameter bounds and draw LHS samples.

Returns `(samples::Matrix{Float64}, bounds::ParameterBounds)`.
"""
function generate_test_samples(cfg::ValidityAnalysisConfig)
    pb = bounds_for_fuel_cell(cfg.fuel_cell_type, cfg.voltage_zone)
    s_cfg = SamplingConfig(
        n_samples = cfg.n_samples,
        method = :lhs,
        seed = cfg.sampling_seed,
        include_reference = true,
    )
    X = generate_lhs_samples(pb, s_cfg)
    return X, pb
end


"""
    classify_batch_simulations(samples, bounds, cfg, run_dir; curves_df=nothing)
        -> (DataFrame, ValidationSummary) | DataFrame

Step 2 of the pipeline: simulate all sampled configurations with AlphaPEM and
classify each polarization curve.

Each row in `samples` is mapped to a `PhysicalParams` struct via
`apply_bounds_to_params`, injected into the reference fuel cell, and simulated
with a fast polarization profile.  The resulting curve is classified by
`ValidityCriteria.classify_polarization_curve`.

If `curves_df` is provided (reusing cached curves), the function skips simulation
and re-classifies the cached curves with the current validation criteria.

# Arguments
- `samples::Matrix{Float64}`: Sample matrix of size `(n, p)` (one configuration per row).
- `bounds::ParameterBounds`: Column-to-parameter mapping and fuel-cell metadata.
- `cfg::ValidityAnalysisConfig`: Master configuration (parallel, criteria, …).
- `run_dir::String`: Output directory for saving curves (if `cfg.save_curves=true`).
- `curves_df::Union{DataFrame,Nothing}`: Pre-computed curves (for reuse). If provided,
  skip simulation. Default: `nothing` (simulate).

# Returns
A `DataFrame` with one row per configuration:
  - Parameter columns (from `samples`)
  - Classification columns: `classification`, `validation_details`, `start_in_range`,
    `is_monotonic`, `has_positive_voltages`, `error_message`
  - (optional) Curve columns: `current_density`, `voltage` (flattened by sample_id)

Also returns `ValidationSummary` if newly simulated. If reusing curves, returns only the DataFrame.
"""
function classify_batch_simulations(samples::Matrix{Float64},
                                     bounds::ParameterBounds,
                                     cfg::ValidityAnalysisConfig,
                                     run_dir::String;
                                     curves_df::Union{DataFrame, Nothing} = nothing)
    n_samples = size(samples, 1)
    param_names = [b.name for b in bounds.bounds]

    # Reference physical parameters: supply all fixed (non-sampled) fields.
    base_params = get_reference_config(bounds.fuel_cell_type)

    # Numerical setup for batch runs: minimum spatial resolution for speed.
    num_params = NumericalParams(nb_gc = cfg.nb_gc)

    # If reusing curves, skip simulation and just re-classify
    if curves_df !== nothing
        @info "Re-classifying $(n_samples) cached curves with current criteria…"

        # Re-classify each configuration using the cached curves
        classifications  = Vector{Symbol}(undef, n_samples)
        val_details      = Vector{Union{String, Missing}}(undef, n_samples)
        start_in_range   = Vector{Union{Bool, Missing}}(undef, n_samples)
        is_monotonic     = Vector{Union{Bool, Missing}}(undef, n_samples)
        has_positive_v   = Vector{Union{Bool, Missing}}(undef, n_samples)
        error_messages   = Vector{Union{String, Missing}}(undef, n_samples)

        for i in 1:n_samples
            # Extract the I and U arrays for this sample
            sample_curves = filter(row -> row.sample_id == i, curves_df)
            if nrow(sample_curves) == 0
                classifications[i] = :failed
                val_details[i] = "no cached curve"
                start_in_range[i] = missing
                is_monotonic[i] = missing
                has_positive_v[i] = missing
                error_messages[i] = "sample_id not found in cache"
                continue
            end

            I_arr = Vector{Float64}(sample_curves.current_density)
            U_arr = Vector{Float64}(sample_curves.voltage)

            # Re-classify using the vector-based overload: Ucell, ifc
            vr = classify_polarization_curve(U_arr, I_arr, cfg.validation_criteria)

            _n2m(x) = x === nothing ? missing : x
            classifications[i] = vr.classification
            val_details[i]     = vr.details
            start_in_range[i]  = _n2m(vr.start_in_range)
            is_monotonic[i]    = _n2m(vr.is_monotonic)
            has_positive_v[i]  = _n2m(vr.has_positive_voltages)
            error_messages[i]  = missing
        end

        # Build results dataframe (no new curves to save)
        data = _build_results_dataframe(
            samples, param_names,
            classifications, val_details,
            start_in_range, is_monotonic, has_positive_v, error_messages
        )

        return data
    end

    # === FRESH SIMULATION MODE ===

    # Pre-allocate result vectors.
    classifications  = Vector{Symbol}(undef, n_samples)
    val_details      = Vector{Union{String, Missing}}(undef, n_samples)
    start_in_range   = Vector{Union{Bool, Missing}}(undef, n_samples)
    is_monotonic     = Vector{Union{Bool, Missing}}(undef, n_samples)
    has_positive_v   = Vector{Union{Bool, Missing}}(undef, n_samples)
    error_messages   = Vector{Union{String, Missing}}(undef, n_samples)

    # If saving curves, prepare a vector to collect (sample_id, I, U) tuples
    curves_collection = cfg.save_curves ? Tuple{Int, Vector{Float64}, Vector{Float64}}[] : nothing

    start_time = time()

    if cfg.parallel && Threads.nthreads() > 1
        # ── Parallel branch ─────────────────────────────────────────────────
        # Restrict BLAS to 1 thread per Julia thread to avoid over-subscription:
        # with N Julia threads each spawning M BLAS threads, you would get N×M
        # OS threads competing for the same CPU cores.
        BLAS.set_num_threads(1)

        prog = Progress(n_samples;
                        desc   = "Batch simulation ($(Threads.nthreads()) threads): ",
                        barlen = 40,
                        color  = :cyan)
        prog_lock = ReentrantLock()

        Threads.@threads for i in 1:n_samples
            cl, det, sir, im, hpv, err, I_arr, U_arr = _simulate_one_configuration(
                samples[i, :], bounds, base_params,
                cfg.polarization_params, num_params, cfg.validation_criteria
            )
            classifications[i] = cl
            val_details[i]     = det
            start_in_range[i]  = sir
            is_monotonic[i]    = im
            has_positive_v[i]  = hpv
            error_messages[i]  = err

            if cfg.save_curves && U_arr !== nothing
                lock(prog_lock) do
                    push!(curves_collection, (i, I_arr, U_arr))
                end
            end

            lock(prog_lock) do
                next!(prog)
            end
        end
        finish!(prog)

    else
        # ── Sequential branch ────────────────────────────────────────────────
        prog = Progress(n_samples;
                        desc   = "Batch simulation (sequential): ",
                        barlen = 40,
                        color  = :cyan)

        for i in 1:n_samples
            cl, det, sir, im, hpv, err, I_arr, U_arr = _simulate_one_configuration(
                samples[i, :], bounds, base_params,
                cfg.polarization_params, num_params, cfg.validation_criteria
            )
            classifications[i] = cl
            val_details[i]     = det
            start_in_range[i]  = sir
            is_monotonic[i]    = im
            has_positive_v[i]  = hpv
            error_messages[i]  = err

            if cfg.save_curves && U_arr !== nothing
                push!(curves_collection, (i, I_arr, U_arr))
            end

            next!(prog)
        end
        finish!(prog)
    end

    elapsed = time() - start_time

    # Save curves to CSV if requested
    if cfg.save_curves && curves_collection !== nothing && !isempty(curves_collection)
        curves_path = abspath(joinpath(run_dir, "curves.csv"))
        _save_curves_to_csv(curves_collection, curves_path)
        @info "  → Curves: $curves_path"
    end

    # Build the final DataFrame.
    data = _build_results_dataframe(
        samples, param_names,
        classifications, val_details,
        start_in_range, is_monotonic, has_positive_v, error_messages
    )

    # Compute aggregate statistics.
    n_valid   = count(==(:valid),   classifications)
    n_invalid = count(==(:invalid), classifications)
    n_failed  = count(==(:failed),  classifications)

    criteria_breakdown = Dict{Symbol, Int}(
        :start_voltage    => count(x -> x === false, start_in_range),
        :monotonicity     => count(x -> x === false, is_monotonic),
        :positive_voltage => count(x -> x === false, has_positive_v),
    )

    summary = ValidationSummary(
        n_samples,
        n_valid,
        n_invalid,
        n_failed,
        criteria_breakdown,
        elapsed,
        Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    )

    @info @sprintf(
        "Batch complete — %d valid (%.1f %%), %d invalid, %d failed  [%.1f s]",
        n_valid, 100.0 * n_valid / n_samples,
        n_invalid, n_failed, elapsed
    )

    return (data = data, summary = summary)
end


"""
    _simulate_one_configuration(sample, bounds, base_params, polar_params, num_params, val_cfg)
      -> (classification, details, start_in_range, is_monotonic, has_positive_v, error_msg,
          I_array, U_array)

Run a single AlphaPEM polarization simulation with the given sampled parameters
and return its classification tuple along with the polarization curve arrays.

Any exception raised by the solver is caught and recorded as `:failed`.

# Returns
A tuple:
  - classification::Symbol (:valid, :invalid, or :failed)
  - details::String
  - start_in_range::Union{Bool, Missing}
  - is_monotonic::Union{Bool, Missing}
  - has_positive_voltages::Union{Bool, Missing}
  - error_message::Union{String, Missing}
  - I_array::Union{Vector{Float64}, Nothing} — current density (A/m²), or nothing if failed
  - U_array::Union{Vector{Float64}, Nothing} — cell voltage (V), or nothing if failed
"""
function _simulate_one_configuration(sample::Vector{Float64},
                                      bounds::ParameterBounds,
                                      base_params,
                                      polar_params::PolarizationParams,
                                      num_params::NumericalParams,
                                      val_cfg::ValidityCriteriaConfig)
    try
        # Map the sample to PhysicalParams (applies bounds clamping + EH-31 constraint).
        modified_params = apply_bounds_to_params(sample, bounds, base_params)

        # Build a fuel-cell object and inject the modified physical parameters.
        fc = create_fuelcell(bounds.fuel_cell_type, bounds.voltage_zone)
        fc.physical_parameters = modified_params

        # Assemble a minimal SimulationConfig (no display, no plots).
        sim_cfg = SimulationConfig(
            type_fuel_cell       = bounds.fuel_cell_type,
            type_current         = polar_params,
            numerical_parameters = num_params,
            voltage_zone         = bounds.voltage_zone,
            type_display         = :no_display,
            display_timing       = :postrun,
        )

        # Create current profile and simulator, then run.
        cd   = create_current(polar_params, fc)
        simu = AlphaPEMSimulator(fc, cd, sim_cfg)
        simulate_model!(simu)

        # Classify the resulting polarization curve.
        vr = classify_polarization_curve(simu, val_cfg)

        # Extract current and voltage arrays from simulation results
        I_arr = Vector{Float64}(simu.results.I)
        U_arr = Vector{Float64}(simu.results.U)

        # `ValidationResult` uses `nothing` for disabled criteria; normalise to `missing`.
        _n2m(x) = x === nothing ? missing : x

        return (
            vr.classification,
            vr.details,
            _n2m(vr.start_in_range),
            _n2m(vr.is_monotonic),
            _n2m(vr.has_positive_voltages),
            missing,
            I_arr,
            U_arr
        )

    catch e
        # Solver divergence, numerical overflow, or any other exception.
        return (:failed, "simulation exception", missing, missing, missing,
                sprint(showerror, e), nothing, nothing)
    end
end


"""
    _build_results_dataframe(samples, param_names, classifications, ...) -> DataFrame

Assemble the classified configurations into a `DataFrame`.

Columns: `sample_id`, one column per sampled parameter, then
`classification`, `validation_details`, `start_in_range`, `is_monotonic`,
`has_positive_voltages`, `error_message`.
"""
function _build_results_dataframe(samples::Matrix{Float64},
                                   param_names::Vector{Symbol},
                                   classifications::Vector{Symbol},
                                   val_details::Vector{Union{String, Missing}},
                                   start_in_range::Vector{Union{Bool, Missing}},
                                   is_monotonic::Vector{Union{Bool, Missing}},
                                   has_positive_v::Vector{Union{Bool, Missing}},
                                   error_messages::Vector{Union{String, Missing}})::DataFrame
    n = size(samples, 1)
    df = DataFrame()
    df.sample_id = 1:n
    for (j, name) in enumerate(param_names)
        df[!, name] = samples[:, j]
    end
    df.classification       = [string(c) for c in classifications]
    df.validation_details   = val_details
    df.start_in_range       = start_in_range
    df.is_monotonic         = is_monotonic
    df.has_positive_voltages = has_positive_v
    df.error_message        = error_messages
    return df
end


# ─────────────────────────────────────────────────────────────────────────────

"""
    find_valid_region(classified_data, reference_config, prim_cfg::PRIMConfig, run_dir::String)
        -> Vector{PRIMResult}

Step 3 of the pipeline: run PRIM via the R `irdpackage` and return the restricted bounds.

## Arguments
- `classified_data::DataFrame` — output of `classify_batch_simulations` (columns:
  parameter names + `classification`).
- `reference_config` — either:
  - a `String` path to an already existing YAML reference file, **or**
  - a `Dict{Symbol,<:Any}` / `NamedTuple` of parameter → value mappings.
    In this case the reference YAML is written to `run_dir/reference_config.yaml`.
- `prim_cfg::PRIMConfig` — PRIM options (methods, bounds on probability, seed, …).
- `run_dir::String` — output directory for PRIM results.

## Returns
`Vector{PRIMResult}` — one entry per IRD method that was successfully run.

## Notes
- The classified CSV used as input to R is written to `run_dir/configurations_for_prim.csv`.
- All output files go to `run_dir`.
"""
function find_valid_region(classified_data::DataFrame,
                            reference_config,
                            prim_cfg::PRIMConfig,
                            run_dir::String)::Vector{PRIMResult}
    mkpath(run_dir)

    # IRD expects a binary target (valid / invalid). Failed simulations are
    # conservatively re-labeled as invalid for the PRIM learning step.
    data_for_prim = copy(classified_data)
    if :classification in propertynames(data_for_prim)
        n_failed = count(==("failed"), data_for_prim.classification)
        if n_failed > 0
            data_for_prim.classification =
                [c == "failed" ? "invalid" : c for c in data_for_prim.classification]
            @info "PRIM preprocessing: re-labeled $n_failed failed simulations as invalid."
        end
    end

    # ── 1. Write classified CSV (only param + classification cols) ────────────
    # Determine which columns are parameter columns (drop diagnostics)
    diag_cols = Set(["sample_id", "validation_details", "start_in_range",
                     "is_monotonic", "has_positive_voltages", "error_message"])
    keep_cols = [c for c in names(data_for_prim) if !(c in diag_cols)]

    csv_path = abspath(joinpath(run_dir, "configurations_for_prim.csv"))
    CSV.write(csv_path, data_for_prim[:, keep_cols])
    @info "Classified data for PRIM written to: $csv_path"

    # ── 2. Resolve reference config YAML path ─────────────────────────────────
    ref_yaml_path = if reference_config isa String
        abspath(reference_config)
    else
        # Write the reference mapping to a temporary YAML file
        ref_path = abspath(joinpath(run_dir, "reference_config.yaml"))
        _write_reference_yaml(reference_config, ref_path)
        @info "Reference config written to: $ref_path"
        ref_path
    end

    # ── 3. Build updated PRIMConfig with the correct data / reference paths ───
    updated_cfg = PRIMConfig(
        data_path             = csv_path,
        ird_package_dir       = prim_cfg.ird_package_dir,
        reference_config_path = ref_yaml_path,
        probability_range     = prim_cfg.probability_range,
        methods               = prim_cfg.methods,
        categorical_params    = prim_cfg.categorical_params,
        seed                  = prim_cfg.seed,
        output_dir            = run_dir,
        target_column         = prim_cfg.target_column,
        positive_class        = prim_cfg.positive_class,
    )

    # ── 4. Delegate to PRIMInterface ──────────────────────────────────────────
    return PRIMInterface.run_prim_analysis(updated_cfg)
end


"""
    _save_curves_to_csv(curves_collection, filepath)

Save polarization curves to a CSV file.

# Arguments
- `curves_collection::Vector{Tuple{Int, Vector, Vector}}` — vector of (sample_id, I_array, U_array) tuples.
- `filepath::String` — output CSV file path.

# CSV Structure
```
sample_id,current_density,voltage
1,0.0,1.23
1,100.0,1.20
1,200.0,1.15
...
```

Each configuration can have multiple I-U points. Rows are grouped by sample_id.
"""
function _save_curves_to_csv(curves_collection::Vector{Tuple{Int, Vector{Float64}, Vector{Float64}}},
                              filepath::String)::Nothing
    mkpath(dirname(filepath))

    # Build a flat structure: for each (sample_id, I_array, U_array), create one row per point
    all_rows = []
    for (sample_id, I_arr, U_arr) in curves_collection
        if length(I_arr) == length(U_arr)
            for (i_val, u_val) in zip(I_arr, U_arr)
                push!(all_rows, Dict(:sample_id => sample_id, :current_density => i_val, :voltage => u_val))
            end
        end
    end

    # Convert to DataFrame and write
    curves_df = DataFrame(all_rows)
    CSV.write(filepath, curves_df)
    return nothing
end


"""
    _load_curves_from_csv(filepath) -> DataFrame

Load polarization curves from a CSV file.

Returns a `DataFrame` with columns: `sample_id`, `current_density`, `voltage`.
"""
function _load_curves_from_csv(filepath::String)::DataFrame
    return CSV.read(filepath, DataFrame)
end


"""
    _write_reference_yaml(reference_config, filepath)

Serialise a reference configuration (Dict, NamedTuple, or PhysicalParams-like struct)
to a YAML file suitable for use as `x_interest` by the R `irdpackage`.

The YAML format is a flat mapping of parameter name to scalar value, e.g.:
```yaml
Hacl: 9.0e-6
Hccl: 1.2e-5
e: 4
Re: 1.0e-6
```
"""
function _write_reference_yaml(reference_config, filepath::String)::Nothing
    mkpath(dirname(filepath))
    open(filepath, "w") do io
        if reference_config isa Dict
            for (k, v) in sort(collect(reference_config); by = x -> string(x[1]))
                println(io, "$(k): $(v)")
            end
        elseif reference_config isa NamedTuple
            for k in keys(reference_config)
                println(io, "$(k): $(reference_config[k])")
            end
        else
            # Assume struct with fieldnames (e.g. PhysicalParams)
            for f in fieldnames(typeof(reference_config))
                println(io, "$(f): $(getfield(reference_config, f))")
            end
        end
    end
    return nothing
end

end # module ValidParameterRegion
