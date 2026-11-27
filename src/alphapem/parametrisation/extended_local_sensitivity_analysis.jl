# -*- coding: utf-8 -*-

"""
    ExtendedLocalSensitivityAnalysis

Extended local sensitivity and parameter identifiability analysis of AlphaPEM.

Local sensitivities (normalised finite differences) are evaluated at many points
sampled globally in the parameter space with a Sobol sequence, then analysed through
their magnitude and collinearity.

The analysis is a **practical** (data-driven), sensitivity-based identifiability analysis, as
opposed to a *structural* one: a parameter is deemed identifiable if it (1) has an observable
effect on the predicted outputs and (2) that effect cannot be reproduced by other parameters.
The conclusions depend on the model structure, the parameter bounds and the chosen measurements,
and do not by themselves establish the physical importance of a parameter.

`run_extended_local_sensitivity_analysis` coordinates the complete analysis, step by step, from the
Sobol candidates to the parameter selection; `analyze_extended_local_sensitivity` repeats the
analysis stage alone on a stored result. `run_local_sensitivity_simulations` executes the accepted
parameter candidates against the operating-case library and returns the signed local-sensitivity
blocks and their validity masks.

# References
- A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of
  pem fuel cell models—part i: Sensitivity analysis and parameter identifiability,"
  Journal of The Electrochemical Society, vol. 167, p. 044504, Mar. 2020.
  https://doi.org/10.1149/1945-7111/ab7091.
"""
module ExtendedLocalSensitivityAnalysis

using CSV
using DataFrames
using Dates
using JLD2
using LinearAlgebra: BLAS
using Logging
using ProgressMeter: Progress, finish!, next!
using Random
using QuasiMonteCarlo
using SHA
using YAML
using Statistics: mean, median, std

using AlphaPEM.Config: NumericalParams, PARAMETER_BOUNDS, PARAMETER_METADATA, PhysicalParams,
                       SimulationConfig, PolarizationCalibrationParams, AscendingPolarizationParams,
                       OperatingConditions,
                       OPERATING_CONDITIONS_BOUNDS, OperatingConditionConstraint,
                       exact_pressure_gap_constraint
using AlphaPEM.Fuelcell: create_fuelcell
using AlphaPEM.Application: run_simulation
using AlphaPEM.Application: serve_interactive_figures!, wait_for_interactive_figures!
using AlphaPEM.Core.Models: _polarization_points_cali, time_history, derived_outputs
using AlphaPEM.Core.Models.PlotHelpers: export_figure_bytes
using AlphaPEM.Currents: AscendingPolarizationCurrent
import WGLMakie
using ..Calibration
using ..Calibration: ParameterBound, ParameterBounds

include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "definitions", "types.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "definitions", "current_density_summaries.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "definitions", "operating_cases.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "definitions", "outputs.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "definitions", "parameters.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "results", "hashing.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "results", "checkpoint.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "simulation", "warning_logger.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "simulation", "sensitivity.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "analysis", "statistics.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "simulation", "sampling.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "analysis", "convergence.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "plots", "sensitivity_plots.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "simulation", "sensitivity_engine.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "analysis", "collinearity.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "analysis", "multi_output.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "analysis", "parameter_selection.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "results", "export.jl"))

export ExtendedLocalSensitivityConfig,
       ExtendedLocalSensitivityResult,
       ExtendedLocalSimulationResult,
       CurrentDensityZone,
       DEFAULT_CURRENT_DENSITY_ZONES,
       OperatingCase,
       generate_operating_case_library,
       ascending_polarization_params,
       VoltageCurve,
       voltage_curve_from_history,
       extract_voltage_curve,
       summarize_current_density_sensitivities,
       ScreeningDataset,
       ScreeningRecord,
       OutputSpec,
       SUPPORTED_OUTPUTS,
       validate_config,
       all_parameter_bounds,
       select_parameters,
       to_unit,
       from_unit,
       PARAMETER_LATEX,
       equidistant_levels,
       extended_local_sensitivity_matrix,
       output_sensitivity_magnitude,
       sensitivity_statistics,
       collinearity_index,
       collinearity_indices,
       median_collinearity,
       influential_mask,
       generate_unit_design,
       nominal_unit_point,
       unit_design_to_physical,
       screening_error,
       screen_candidates,
       alphapem_screening_simulator,
       convergence_sample_counts,
       SensitivityConvergenceSnapshot,
       candidate_rms_magnitudes,
       convergence_snapshot,
       SensitivityConvergencePlot,
       create_sensitivity_convergence_plot,
       update_sensitivity_convergence_plot!,
       SENSITIVITY_INTERPRETATION_WARNING,
       VoltageSensitivityPlot,
       create_voltage_sensitivity_plot,
       add_voltage_sensitivity_candidate!,
       finish_voltage_sensitivity_plot!,
       COLLINEARITY_INTERPRETATION_WARNING,
       CollinearityPlot,
       create_collinearity_plot,
       finish_collinearity_plot!,
       serve_sensitivity_plots!,
       FigureAutosave,
       autosave_figure!,
       candidate_collinearity,
       CollinearitySummary,
       collinearity_summary,
       collinear_pairs,
       CollinearityFilterResult,
       filter_collinear_parameters,
       select_parameters_by_sensitivity,
       extended_local_result_path,
       save_extended_local_sensitivity_result,
       load_extended_local_sensitivity_result,
       sensitivity_summary_table,
       analyze_extended_local_sensitivity,
       run_extended_local_sensitivity_analysis,
       reanalyze_extended_local_sensitivity,
       run_candidate_screening,
       restored_screening_records,
       create_sensitivity_plots,
       fill_sensitivity_plots!,
       finish_sensitivity_plots!,
       export_sensitivity_tables,
       write_extended_local_summary_yaml,
       FUTURE_OUTPUTS,
       MeasurementSet,
       DEFAULT_MEASUREMENT_SETS,
       combine_outputs,
       compare_measurement_sets,
       run_local_sensitivity_simulations,
       alphapem_local_sensitivity_simulator


# ============================================================
# ANALYSIS COORDINATOR
# ============================================================

"""
    run_extended_local_sensitivity_analysis(cfg; interactive_figures=true) -> NamedTuple

Run the complete extended local sensitivity analysis.

The stages are: choice of the analysed parameters, screening of Sobol candidates against
experimental polarization curves, generation of the operating cases, local-sensitivity simulations,
archiving of the result, and analysis ([`analyze_extended_local_sensitivity`](@ref)).

Two files make a long run safe and reusable:

- **Checkpoint** (`extended_local_sensitivity_checkpoint.jld2`, in `output_dir`): written
  regularly during the simulations (every `cfg.checkpoint_frequency` finished work units). If the run
  is interrupted, call this function again with `cfg.resume = true` to continue where it stopped; the
  screening is then not repeated. The checkpoint is deleted once the final archive is written.
- **Result archive** (`extended_local_sensitivity_result.jld2`, written when `cfg.save_result` is
  `true`): the complete output of the simulations, reloadable with
  [`load_extended_local_sensitivity_result`](@ref) to repeat the analysis without recomputing.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Analysis configuration.
interactive_figures : Bool, optional
    Serve the figures in the browser while the simulations run. Default is `true`.

Returns
-------
NamedTuple
    `result` (simulation result), `analysis` (see [`analyze_extended_local_sensitivity`](@ref)),
    `records` (screening records) and `fixed_parameters` (parameters excluded from the analysis).

Raises
------
ArgumentError
    If the configuration is invalid or no candidate passes the screening.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function run_extended_local_sensitivity_analysis(
    cfg::ExtendedLocalSensitivityConfig;
    interactive_figures::Bool = true,
)
    # 0. Preliminary checks
    validate_config(cfg)

    # 1. Analysed parameters (integer parameters are fixed and reported separately)
    params, fixed_parameters = select_parameters(cfg)

    # 2. Candidate screening, skipped when a resumed checkpoint already holds it
    records = restored_screening_records(cfg)
    if records === nothing
        records = run_candidate_screening(cfg, params)
    else
        @info "Screening records restored from the checkpoint" count = length(records)
    end
    accepted_count = count(r -> r.status == :accepted, records)
    accepted_count > 0 || throw(ArgumentError("No candidate passed the screening"))

    # 3. Operating-case library and budget warning
    operating_cases = generate_operating_case_library(cfg)
    @info "Projected local-sensitivity run count; all operating cases are used." accepted_parameter_points = accepted_count parameter_count = length(params) finite_difference_levels = cfg.n_levels operating_case_count = length(operating_cases) projected_simulations = accepted_count * length(params) * cfg.n_levels * length(operating_cases)

    # 4. Figures, filled live while the simulations run
    plots = create_sensitivity_plots(params, accepted_count, cfg; interactive = interactive_figures)

    # 5. Local-sensitivity simulations (checkpointed; archived when finished)
    result = run_local_sensitivity_simulations(
        records, params, operating_cases, cfg;
        on_convergence_update = snapshot ->
            update_sensitivity_convergence_plot!(plots.convergence, snapshot),
        on_candidate_update = (_, magnitudes) ->
            add_voltage_sensitivity_candidate!(plots.voltage, magnitudes))

    # 6. Analysis of the finished result
    analysis = analyze_extended_local_sensitivity(result; plots)

    return (result = result, analysis = analysis, records = records,
            fixed_parameters = fixed_parameters)
end


"""
    analyze_extended_local_sensitivity(result; output_dir=result.config.output_dir, plots=nothing,
                                       figures=true, interactive_figures=true) -> NamedTuple

Analyse the sensitivities of a finished run: collinearity, parameter selection, figures and tables.

This stage needs no simulation. It can therefore be repeated on a result reloaded with
load_extended_local_sensitivity_result, for instance after changing the thresholds
`sensitivity_threshold`, `collinearity_threshold` or `psi_th` of the configuration.

Its outputs, written in `output_dir`, are the figures (`convergence.pdf`, `voltage_sensitivity.pdf`,
`collinearity.pdf`), `collinear_pairs.csv` (strongly collinear pairs) and `sensitivity_summary.csv`
(one row per parameter: median sensitivity, collinearity removal and final selection) and
`extended_local_sensitivity_summary.yaml` (versions, seed, hash, design sizes, outcomes, selection).

Parameters
----------
result : ExtendedLocalSimulationResult
    Finished sensitivity simulations.
output_dir : AbstractString, optional
    Destination of the figures and tables. Default is the configured `output_dir`.
plots : NamedTuple or nothing, optional
    Figures already filled while the simulations ran. If `nothing`, figures are created and filled
    from `result`.
figures : Bool, optional
    Create the figures; `false` skips every figure. Default is `true`.
interactive_figures : Bool, optional
    Serve newly created figures in the browser. Default is `true`.

Returns
-------
NamedTuple
    `collinearity`, `pairs`, `collinearity_filter`, `selection`, `summary`, `summary_yaml` (path of
    the YAML summary) and `plots`.
"""
function analyze_extended_local_sensitivity(
    result::ExtendedLocalSimulationResult;
    output_dir::AbstractString = result.config.output_dir,
    plots = nothing,
    figures::Bool = true,
    interactive_figures::Bool = true,
)
    cfg = result.config
    params = result.params

    # 1. Figures: created and filled from the result unless supplied by a live run
    if !figures
        plots = nothing
    elseif plots === nothing
        plots = create_sensitivity_plots(params, length(result.candidates), cfg;
                                         output_dir, interactive = interactive_figures)
        fill_sensitivity_plots!(plots, result)
    end

    # 2. Collinearity of the sensitivity vectors, summarised over all candidates
    collinearity = collinearity_summary(result)
    pairs = collinear_pairs(
        collinearity, params; collinearity_threshold = cfg.collinearity_threshold,
        sensitivity_threshold = cfg.sensitivity_threshold)
    plots === nothing || finish_sensitivity_plots!(plots, result, collinearity)

    # 3. Parameter-subset selection: collinearity filter, then sensitivity threshold
    collinearity_filter = filter_collinear_parameters(collinearity, params; psi_th = cfg.psi_th)
    selection = select_parameters_by_sensitivity(
        collinearity, params, collinearity_filter.retained;
        sensitivity_threshold = cfg.sensitivity_threshold)

    # 4. Tables
    summary = sensitivity_summary_table(collinearity, params, collinearity_filter, selection)
    export_sensitivity_tables(output_dir, pairs, summary)
    summary_yaml = write_extended_local_summary_yaml(
        output_dir, result, (collinearity_filter = collinearity_filter, selection = selection))

    return (collinearity = collinearity, pairs = pairs, collinearity_filter = collinearity_filter,
            selection = selection, summary = summary, summary_yaml = summary_yaml, plots = plots)
end

"""
    reanalyze_extended_local_sensitivity(cfg; path=extended_local_result_path(cfg),
                                         interactive_figures=true) -> NamedTuple

Repeat the analysis stage on an archived result, without simulating again.

The simulations come from the archive written by a previous
run_extended_local_sensitivity_analysis; the settings of the analysis stage (thresholds
`sensitivity_threshold`, `collinearity_threshold`, `psi_th`, `output_dir`) come from `cfg`, so they can
be changed freely.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Configuration whose analysis settings are applied.
path : AbstractString, optional
    Archive to read. Default is [`extended_local_result_path`](@ref).
interactive_figures : Bool, optional
    Serve the figures in the browser. Default is `true`.

Returns
-------
NamedTuple
    Same fields as [`run_extended_local_sensitivity_analysis`](@ref).

Raises
------
ArgumentError
    If the archive is missing or unsupported.
"""
function reanalyze_extended_local_sensitivity(
    cfg::ExtendedLocalSensitivityConfig;
    path::AbstractString = extended_local_result_path(cfg),
    interactive_figures::Bool = true,
)
    # 1. Stored simulations, with the analysis settings of the current configuration
    validate_config(cfg)
    stored = load_extended_local_sensitivity_result(path)
    result = ExtendedLocalSimulationResult(
        cfg, stored.params, stored.screening_records, stored.candidates, stored.theta,
        stored.operating_cases, stored.sensitivities, stored.validity_masks, stored.statuses,
        stored.unit_times, stored.run_status_counts, stored.run_errors, stored.execution_time)

    # 2. Analysis stage
    analysis = analyze_extended_local_sensitivity(result; interactive_figures)

    _, fixed_parameters = select_parameters(cfg)
    return (result = result, analysis = analysis, records = result.screening_records,
            fixed_parameters = fixed_parameters)
end

end  # module ExtendedLocalSensitivityAnalysis
