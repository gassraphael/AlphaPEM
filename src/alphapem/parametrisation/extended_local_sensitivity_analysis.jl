# -*- coding: utf-8 -*-

"""
    ExtendedLocalSensitivityAnalysis

Extended local sensitivity and parameter identifiability analysis of AlphaPEM.

Local sensitivities (normalised finite differences) are evaluated at many points
sampled globally in the parameter space with a Sobol sequence, then analysed through
their magnitude, collinearity and Hessian-based identifiability.

The analysis is a **practical** (data-driven), sensitivity-based identifiability analysis, as
opposed to a *structural* one: a parameter is deemed identifiable if it (1) has an observable
effect on the predicted outputs and (2) that effect cannot be reproduced by other parameters.
The conclusions depend on the model structure, the parameter bounds and the chosen measurements,
and do not by themselves establish the physical importance of a parameter.

`run_local_sensitivity_simulations` executes the accepted parameter candidates against the
operating-case library and returns the signed local-sensitivity blocks and their validity masks.
The result is the simulation stage consumed by the later summary and identifiability analyses.

# References
- A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of
  pem fuel cell models—part i: Sensitivity analysis and parameter identifiability,"
  Journal of The Electrochemical Society, vol. 167, p. 044504, Mar. 2020.
  https://doi.org/10.1149/1945-7111/ab7091.
- P. Affonso Nobrega, *marapendi* (PERSEE – Mines Paris – PSL University), MIT licence,
  <https://github.com/affonsop/marapendi>.

The module also hosts the simulation-infrastructure helpers shared with
`SobolSensitivityAnalysis` (warning logger, hashing, checkpoint I/O).
"""
module ExtendedLocalSensitivityAnalysis

using DataFrames
using JLD2
using LinearAlgebra: BLAS
using Logging
using ProgressMeter: Progress, finish!, next!
using Random
using QuasiMonteCarlo
using SHA
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
using AlphaPEM.Currents: AscendingPolarizationCurrent
import WGLMakie
using ..Calibration
using ..Calibration: ParameterBound, ParameterBounds

include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "warning_logger.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "hashing.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "checkpoint_io.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "types.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "current_density_summaries.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "operating_cases.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "outputs.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "parameters.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "sensitivity.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "statistics.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "sampling.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "convergence.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "sensitivity_plots.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "sensitivity_engine.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "collinearity.jl"))

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
       candidate_collinearity,
       CollinearitySummary,
       collinearity_summary,
       collinear_pairs,
       run_local_sensitivity_simulations,
       alphapem_local_sensitivity_simulator

end  # module ExtendedLocalSensitivityAnalysis
