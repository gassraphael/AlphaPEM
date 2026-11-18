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
using Logging
using Random
using QuasiMonteCarlo
using SHA
using Statistics: median, std

using AlphaPEM.Config: NumericalParams, PARAMETER_BOUNDS, PARAMETER_METADATA, PhysicalParams,
                       SimulationConfig, PolarizationCalibrationParams
using AlphaPEM.Fuelcell: create_fuelcell
using AlphaPEM.Application: run_simulation
using AlphaPEM.Core.Models: _polarization_points_cali
using ..Calibration
using ..Calibration: ParameterBound, ParameterBounds

include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "warning_logger.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "hashing.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "checkpoint_io.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "types.jl"))       # step 1: config / result types
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "parameters.jl"))  # step 2: parameters, lin/log scaling
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "sensitivity.jl"))  # step 3.1: normalised sensitivities (eq. 7)
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "statistics.jl"))   # step 3.2: magnitude median / std (eq. 8-9), collinearity (eq. 10)
include(joinpath(@__DIR__, "sensitivity_analysis", "extended_local", "sampling.jl"))     # step 4: scrambled Sobol design + screening against the user's datasets

export ExtendedLocalSensitivityConfig,
       ExtendedLocalSensitivityResult,
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
       alphapem_screening_simulator

end  # module ExtendedLocalSensitivityAnalysis
