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

# Example
```julia
cfg = ValidityAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    n_samples      = 500,       # quick test
    parallel       = false,
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
    checkpoint_interval::Int            = 100
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
    run_validity_analysis(cfg::ValidityAnalysisConfig,
                          prim_cfg::Union{PRIMConfig,Nothing} = nothing)
        ::ValidityAnalysisResult

Execute the complete validity-analysis pipeline (Steps 3–6) and return all results.

## Pipeline

| Step | Action                                                                                  |
|------|-----------------------------------------------------------------------------------------|
| 1    | Draw `cfg.n_samples` configurations by Latin Hypercube Sampling                         |
| 2    | Simulate each configuration with AlphaPEM and classify the polarization curve           |
| 3    | *(optional)* Restrict bounds via PRIM/MaxBox if `prim_cfg` is provided                 |
| 4    | Export results: classified CSV, bounds YAML, validation summary, PRIM report            |

## Arguments
- `cfg::ValidityAnalysisConfig` — master configuration (sampling, simulation, criteria, …).
- `prim_cfg::Union{PRIMConfig,Nothing}` — PRIM options.  Pass `nothing` (default) to skip
  Step 5 and return only the classified dataset + original bounds.

## Returns
A `ValidityAnalysisResult` with all intermediate outputs and paths to generated files.

## Example
```julia
using AlphaPEM.Parametrisation.ValidParameterRegion

cfg = ValidityAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    n_samples      = 2000,
    output_dir     = "results/model_validity",
)

# Steps 1+2+4 only (without PRIM):
result = run_validity_analysis(cfg)

# Full pipeline with PRIM:
prim_cfg = PRIMConfig(
    ird_package_dir       = "external/supplementary_2023_ird/irdpackage",
    reference_config_path = "results/model_validity/reference_config.yaml",
    output_dir            = "results/model_validity",
)
result = run_validity_analysis(cfg, prim_cfg)

println("Restricted bounds: ", result.output_files[:restricted_bounds_PRIM_yaml])
```
"""
function run_validity_analysis(cfg::ValidityAnalysisConfig,
                                prim_cfg::Union{PRIMConfig, Nothing} = nothing
                                )::ValidityAnalysisResult
    overall_start = time()
    mkpath(VALIDITY_OUTPUT_DIR)
    output_files = Dict{Symbol, String}()

    # ── STEP 1: LHS sampling ──────────────────────────────────────────────────
    @info "STEP 1 — Generating $(cfg.n_samples) LHS samples…"
    X, pb = generate_test_samples(cfg)
    @info @sprintf("  → Sample matrix: %d × %d", size(X, 1), size(X, 2))

    # Extract original bounds from ParameterBounds
    orig_bounds = Dict{Symbol, Tuple{Float64, Float64}}(
        b.name => (b.min, b.max) for b in pb.bounds
    )

    # Export original bounds to YAML
    orig_bounds_path = abspath(joinpath(VALIDITY_OUTPUT_DIR, "original_bounds.yaml"))
    export_parameter_bounds(orig_bounds, orig_bounds_path;
                            method   = :prior,
                            metadata = Dict("fuel_cell_type" => string(cfg.fuel_cell_type),
                                            "voltage_zone"   => string(cfg.voltage_zone)))
    output_files[:original_bounds_yaml] = orig_bounds_path
    @info "  → Original bounds: $orig_bounds_path"

    # ── STEP 2: Batch simulation + classification ─────────────────────────────
    @info "STEP 2 — Batch simulation and classification…"
    data, summary = classify_batch_simulations(X, pb, cfg)

    # Export classified CSV
    export_cfg = ExportConfig(output_dir = VALIDITY_OUTPUT_DIR, timestamp = true, overwrite = true)
    csv_path   = export_classified_configurations(data, export_cfg)
    output_files[:classified_csv] = csv_path
    @info "  → Classified CSV: $csv_path"

    # Export validation summary
    summary_path = abspath(joinpath(VALIDITY_OUTPUT_DIR, "validation_summary.txt"))
    export_validation_summary(summary, summary_path)
    output_files[:validation_summary] = summary_path
    @info "  → Validation summary: $summary_path"

    # ── STEP 3: PRIM (optional) ───────────────────────────────────────────────
    prim_results       = PRIMResult[]
    restricted_bounds  = orig_bounds          # default: unchanged when PRIM is skipped

    if prim_cfg !== nothing
        @info "STEP 3 — Running PRIM analysis…"
        ref_config   = get_reference_config(cfg.fuel_cell_type)
        prim_results = find_valid_region(data, ref_config, prim_cfg)

        for r in prim_results
            mstr = string(r.method)

            # Export restricted bounds YAML
            bounds_path = abspath(joinpath(VALIDITY_OUTPUT_DIR,
                                           "restricted_bounds_$(mstr).yaml"))
            export_parameter_bounds(r.restricted_bounds, bounds_path;
                                    method   = r.method,
                                    metadata = Dict(
                                        "fuel_cell_type" => string(cfg.fuel_cell_type),
                                        "precision"      => r.precision,
                                        "coverage"       => r.coverage,
                                    ))
            output_files[Symbol("restricted_bounds_$(mstr)_yaml")] = bounds_path
            @info "  → Restricted bounds ($(mstr)): $bounds_path"

            # Export PRIM comparison report
            report_path = abspath(joinpath(VALIDITY_OUTPUT_DIR,
                                           "prim_report_$(mstr).txt"))
            prim_metrics = merge(
                r.rf_metrics,
                Dict("box_precision" => r.precision,
                     "box_coverage"  => r.coverage,
                     "n_inside_box"  => Float64(r.n_inside_box),
                     "n_valid_inside"=> Float64(r.n_valid_inside))
            )
            export_prim_report(orig_bounds, r.restricted_bounds, report_path;
                               prim_metrics = prim_metrics)
            output_files[Symbol("prim_report_$(mstr)")] = report_path
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
    classify_batch_simulations(samples, bounds, cfg) -> (DataFrame, ValidationSummary)

Step 2 of the pipeline: simulate all sampled configurations with AlphaPEM and
classify each polarization curve.

Each row in `samples` is mapped to a `PhysicalParams` struct via
`apply_bounds_to_params`, injected into the reference fuel cell, and simulated
with a fast polarization profile.  The resulting curve is classified by
`ValidityCriteria.classify_polarization_curve`.

# Arguments
- `samples::Matrix{Float64}`: Sample matrix of size `(n, p)` (one configuration per row).
- `bounds::ParameterBounds`: Column-to-parameter mapping and fuel-cell metadata.
- `cfg::ValidityAnalysisConfig`: Master configuration (parallel, checkpoint, criteria, …).

# Returns
A named tuple `(data::DataFrame, summary::ValidationSummary)` where:
- `data` has one row per configuration with parameter columns + classification columns.
- `summary` holds aggregate statistics (counts, criterion breakdown, elapsed time).
"""
function classify_batch_simulations(samples::Matrix{Float64},
                                     bounds::ParameterBounds,
                                     cfg::ValidityAnalysisConfig)
    n_samples = size(samples, 1)
    param_names = [b.name for b in bounds.bounds]

    # Reference physical parameters: supply all fixed (non-sampled) fields.
    base_params = get_reference_config(bounds.fuel_cell_type)

    # Numerical setup for batch runs: minimum spatial resolution for speed.
    num_params = NumericalParams(nb_gc = cfg.nb_gc)

    # Pre-allocate result vectors.
    classifications  = Vector{Symbol}(undef, n_samples)
    val_details      = Vector{Union{String, Missing}}(undef, n_samples)
    start_in_range   = Vector{Union{Bool, Missing}}(undef, n_samples)
    is_monotonic     = Vector{Union{Bool, Missing}}(undef, n_samples)
    has_positive_v   = Vector{Union{Bool, Missing}}(undef, n_samples)
    error_messages   = Vector{Union{String, Missing}}(undef, n_samples)

    # Checkpointing setup.
    do_checkpoints = cfg.checkpoint_interval > 0
    ckpt_dir = joinpath(VALIDITY_OUTPUT_DIR, "checkpoints")
    do_checkpoints && mkpath(ckpt_dir)

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
            cl, det, sir, im, hpv, err = _simulate_one_configuration(
                samples[i, :], bounds, base_params,
                cfg.polarization_params, num_params, cfg.validation_criteria
            )
            classifications[i] = cl
            val_details[i]     = det
            start_in_range[i]  = sir
            is_monotonic[i]    = im
            has_positive_v[i]  = hpv
            error_messages[i]  = err

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
            cl, det, sir, im, hpv, err = _simulate_one_configuration(
                samples[i, :], bounds, base_params,
                cfg.polarization_params, num_params, cfg.validation_criteria
            )
            classifications[i] = cl
            val_details[i]     = det
            start_in_range[i]  = sir
            is_monotonic[i]    = im
            has_positive_v[i]  = hpv
            error_messages[i]  = err

            next!(prog)

            # Periodic checkpoint: flush partial results to CSV.
            if do_checkpoints && i % cfg.checkpoint_interval == 0
                partial = _build_results_dataframe(
                    samples[1:i, :], param_names,
                    classifications[1:i], val_details[1:i],
                    start_in_range[1:i], is_monotonic[1:i],
                    has_positive_v[1:i], error_messages[1:i]
                )
                ckpt_path = joinpath(ckpt_dir, @sprintf("checkpoint_%06d.csv", i))
                CSV.write(ckpt_path, partial)
                @info @sprintf("  Checkpoint saved: %s  (%d valid / %d sim so far)",
                               basename(ckpt_path),
                               count(==(:valid), classifications[1:i]), i)
            end
        end
        finish!(prog)
    end

    elapsed = time() - start_time

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
      -> (classification, details, start_in_range, is_monotonic, has_positive_v, error_msg)

Run a single AlphaPEM polarization simulation with the given sampled parameters
and return its classification tuple.

Any exception raised by the solver is caught and recorded as `:failed`.
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

        # `ValidationResult` uses `nothing` for disabled criteria; normalise to `missing`.
        _n2m(x) = x === nothing ? missing : x

        return (
            vr.classification,
            vr.details,
            _n2m(vr.start_in_range),
            _n2m(vr.is_monotonic),
            _n2m(vr.has_positive_voltages),
            missing
        )

    catch e
        # Solver divergence, numerical overflow, or any other exception.
        return (:failed, "simulation exception", missing, missing, missing, sprint(showerror, e))
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
    find_valid_region(classified_data, reference_config, prim_cfg::PRIMConfig)
        -> Vector{PRIMResult}

Step 3 of the pipeline: run PRIM via the R `irdpackage` and return the restricted bounds.

## Arguments
- `classified_data::DataFrame` — output of `classify_batch_simulations` (columns:
  parameter names + `classification`).
- `reference_config` — either:
  - a `String` path to an already existing YAML reference file, **or**
  - a `Dict{Symbol,<:Any}` / `NamedTuple` of parameter → value mappings.
    In this case the reference YAML is written to
    `VALIDITY_OUTPUT_DIR/reference_config.yaml` automatically.
- `prim_cfg::PRIMConfig` — PRIM options (methods, bounds on probability, seed, …).

## Returns
`Vector{PRIMResult}` — one entry per IRD method that was successfully run.

## Notes
- The classified CSV used as input to R is written to
  `VALIDITY_OUTPUT_DIR/classified_for_prim.csv` (only parameter + classification
  columns; diagnostic columns are dropped).
- All output files go to `VALIDITY_OUTPUT_DIR` regardless of any path set in `prim_cfg`.
"""
function find_valid_region(classified_data::DataFrame,
                            reference_config,
                            prim_cfg::PRIMConfig)::Vector{PRIMResult}
    mkpath(VALIDITY_OUTPUT_DIR)

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

    csv_path = abspath(joinpath(VALIDITY_OUTPUT_DIR, "classified_for_prim.csv"))
    CSV.write(csv_path, data_for_prim[:, keep_cols])
    @info "Classified data written to: $csv_path"

    # ── 2. Resolve reference config YAML path ─────────────────────────────────
    ref_yaml_path = if reference_config isa String
        abspath(reference_config)
    else
        # Write the reference mapping to a temporary YAML file
        ref_path = abspath(joinpath(VALIDITY_OUTPUT_DIR, "reference_config.yaml"))
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
        output_dir            = VALIDITY_OUTPUT_DIR,
        target_column         = prim_cfg.target_column,
        positive_class        = prim_cfg.positive_class,
    )

    # ── 4. Delegate to PRIMInterface ──────────────────────────────────────────
    return PRIMInterface.run_prim_analysis(updated_cfg)
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
