# -*- coding: utf-8 -*-

# ============================================================
# RESULT ARCHIVE (SAVE AND RELOAD)
# ============================================================

# Bump when the layout of the archive changes incompatibly.
const EXTENDED_LOCAL_RESULT_FORMAT_VERSION = 1

"""
    extended_local_result_path(cfg) -> String

Default location of the final result archive: `extended_local_sensitivity_result.jld2` in
`cfg.output_dir`.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Configuration providing `output_dir`.

Returns
-------
String
    Path of the archive.
"""
extended_local_result_path(cfg::ExtendedLocalSensitivityConfig)::String =
    joinpath(cfg.output_dir, "extended_local_sensitivity_result.jld2")


"""
    save_extended_local_sensitivity_result(result; path=extended_local_result_path(result.config))
        -> Bool

Archive every computed quantity of a finished sensitivity run.

In plain terms: the simulations are by far the most expensive part of the analysis. This file keeps
their complete output (all sensitivity blocks, validity masks, screening records, candidates,
operating cases, configuration and timings), so that the analysis stages (collinearity, selection,
figures) can be repeated or modified later with
[`load_extended_local_sensitivity_result`](@ref), without simulating again.

The file is written to a temporary name and then renamed, so an interrupted write never leaves a
corrupt archive. A failure only emits a warning, since it must not discard a finished run.

Parameters
----------
result : ExtendedLocalSimulationResult
    Result returned by [`run_local_sensitivity_simulations`](@ref).
path : AbstractString, optional
    Destination file. Default is [`extended_local_result_path`](@ref).

Returns
-------
Bool
    `true` if the archive was written.
"""
function save_extended_local_sensitivity_result(
    result::ExtendedLocalSimulationResult;
    path::AbstractString = extended_local_result_path(result.config),
)::Bool
    mkpath(dirname(path))
    tmp_path = path * ".tmp"
    try
        JLD2.jldsave(tmp_path; format_version = EXTENDED_LOCAL_RESULT_FORMAT_VERSION, result)
        mv(tmp_path, path; force = true)
        return true
    catch e
        @warn "Failed to save the extended local sensitivity result archive: $e"
        rm(tmp_path; force = true)
        return false
    end
end


"""
    load_extended_local_sensitivity_result(path) -> ExtendedLocalSimulationResult

Read an archive written by [`save_extended_local_sensitivity_result`](@ref).

Parameters
----------
path : AbstractString
    Archive file.

Returns
-------
ExtendedLocalSimulationResult
    The stored result, ready for [`analyze_extended_local_sensitivity`](@ref).

Raises
------
ArgumentError
    If the file does not exist, is not an archive of this analysis, or has an unsupported version.
"""
function load_extended_local_sensitivity_result(path::AbstractString)::ExtendedLocalSimulationResult
    isfile(path) || throw(ArgumentError("No result archive at $path"))
    data = try
        JLD2.load(path)
    catch e
        e isa InterruptException && rethrow()
        throw(ArgumentError("$path is not a readable JLD2 archive"))
    end
    get(data, "format_version", nothing) == EXTENDED_LOCAL_RESULT_FORMAT_VERSION &&
        get(data, "result", nothing) isa ExtendedLocalSimulationResult ||
        throw(ArgumentError("$path is not a supported extended local sensitivity archive"))
    return data["result"]
end


# ============================================================
# SUMMARY TABLE
# ============================================================

"""
    sensitivity_summary_table(collinearity, params, collinearity_filter, selection) -> DataFrame

One-row-per-parameter summary of the analysis outcome.

Parameters
----------
collinearity : CollinearitySummary
    Result of [`collinearity_summary`](@ref), providing the median sensitivity magnitudes.
params : Vector{ParameterBound}
    Analysed parameters in summary order.
collinearity_filter : CollinearityFilterResult
    Result of [`filter_collinear_parameters`](@ref).
selection : NamedTuple
    Result of [`select_parameters_by_sensitivity`](@ref).

Returns
-------
DataFrame
    Columns `parameter`, `median_sensitivity`, `removed_as_collinear` (with `collinear_with`, the
    parameter kept instead, or `missing`) and `selected`.
"""
function sensitivity_summary_table(
    collinearity::CollinearitySummary,
    params::Vector{ParameterBound},
    collinearity_filter::CollinearityFilterResult,
    selection,
)::DataFrame
    names = Symbol[p.name for p in params]
    partner = Dict(zip(collinearity_filter.removed.removed, collinearity_filter.removed.kept))
    return DataFrame(
        parameter = names,
        median_sensitivity = collinearity.sensitivity_median,
        removed_as_collinear = [name in keys(partner) for name in names],
        collinear_with = Union{Symbol,Missing}[get(partner, name, missing) for name in names],
        selected = [name in selection.selected for name in names])
end


# ============================================================
# TABLES
# ============================================================

"""
    export_sensitivity_tables(output_dir, pairs, summary) -> Nothing

Write `collinear_pairs.csv` and `sensitivity_summary.csv` in `output_dir`.

Parameters
----------
output_dir : AbstractString
    Destination directory, created if needed.
pairs : DataFrame
    Result of [`collinear_pairs`](@ref).
summary : DataFrame
    Result of [`sensitivity_summary_table`](@ref).

Returns
-------
Nothing
"""
function export_sensitivity_tables(output_dir::AbstractString, pairs::DataFrame, summary::DataFrame)::Nothing
    mkpath(output_dir)
    CSV.write(joinpath(output_dir, "collinear_pairs.csv"), pairs)
    CSV.write(joinpath(output_dir, "sensitivity_summary.csv"), summary)
    return nothing
end



# ============================================================
# YAML SUMMARY
# ============================================================

"""
    write_extended_local_summary_yaml(output_dir, result, analysis) -> String

Write `extended_local_sensitivity_summary.yaml`, a human-readable record of the run.

The file gathers what is needed to reproduce and interpret the analysis: software versions, seed and
configuration hash, sizes of the design (candidates, parameters, operating cases), simulation
outcomes, thresholds, and the parameter selection (fixed, removed as collinear, selected).

Parameters
----------
output_dir : AbstractString
    Destination directory, created if needed.
result : ExtendedLocalSimulationResult
    Finished sensitivity simulations.
analysis : NamedTuple
    Outcome of the analysis stage, with the fields `collinearity_filter` and `selection`.

Returns
-------
String
    Path of the written file.
"""
function write_extended_local_summary_yaml(
    output_dir::AbstractString,
    result::ExtendedLocalSimulationResult,
    analysis,
)::String
    cfg = result.config

    # 1. Reproducibility: versions, seed and configuration hash
    alphapem_version = try
        string(pkgversion(@__MODULE__))
    catch
        "unknown"
    end
    metadata = Dict{String, Any}(
        "julia_version" => string(VERSION),
        "alphapem_version" => alphapem_version,
        "seed" => cfg.seed,
        "config_sha256" => bytes2hex(SHA.sha256(codeunits(_struct_to_string(cfg)))),
        "timestamp" => Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"),
        "execution_time_s" => round(result.execution_time; digits = 2),
    )

    # 2. Design sizes and simulation outcomes
    screening = Dict{String, Any}(
        "n_candidates" => length(result.screening_records),
        "n_accepted" => length(result.candidates),
        "norm" => string(cfg.screening_norm),
        "datasets" => ["$(d.fuel_cell_type)_$(d.year)_$(d.voltage_zone)" for d in cfg.screening_datasets],
    )
    design = Dict{String, Any}(
        "fuel_cell_type" => string(cfg.fuel_cell_type),
        "n_parameters" => length(result.params),
        "n_levels" => cfg.n_levels,
        "n_operating_cases" => length(result.operating_cases),
        "current_density_grid_A_cm2" => cfg.current_density_grid,
    )
    simulations = Dict{String, Any}(
        "block_status_counts" => Dict(string(k) => v for (k, v) in result.run_status_counts),
        "run_errors" => Dict(string(k) => v for (k, v) in result.run_errors),
    )

    # 3. Thresholds and parameter selection
    thresholds = Dict{String, Any}(
        "sensitivity_threshold" => cfg.sensitivity_threshold,
        "collinearity_threshold" => cfg.collinearity_threshold,
        "psi_th" => cfg.psi_th,
    )
    removed = analysis.collinearity_filter.removed
    selection = Dict{String, Any}(
        "analysed" => [string(p.name) for p in result.params],
        "removed_as_collinear" => [Dict("parameter" => string(r.removed), "kept" => string(r.kept),
                                        "median_psi" => r.psi_median) for r in eachrow(removed)],
        "below_sensitivity_threshold" => string.(analysis.selection.below_threshold),
        "selected" => string.(analysis.selection.selected),
    )

    # 4. Write the file
    mkpath(output_dir)
    path = joinpath(output_dir, "extended_local_sensitivity_summary.yaml")
    YAML.write_file(path, Dict{String, Any}(
        "metadata" => metadata, "design" => design, "screening" => screening,
        "simulations" => simulations, "thresholds" => thresholds, "selection" => selection))
    return path
end
