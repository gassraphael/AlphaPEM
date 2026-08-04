# -*- coding: utf-8 -*-

"""
    SobolSensitivityAnalysis

Global sensitivity analysis of AlphaPEM using the Sobol method.

This module performs variance-based sensitivity analysis directly on AlphaPEM
simulations (no surrogate model). It follows the workflow of the student notebook
`03_sobol_analysis_AlphaPEM.ipynb`:

1. Generate Sobol samples for the input parameters.
2. Run AlphaPEM for each sample and extract the polarization curve.
3. Impute missing curves with KNN when simulations fail.
4. Compute the area under the curve (AUC) in three characteristic regions:
   activation, ohmic, and mass transport.
5. Compute Sobol first-order (S1), total-order (ST), and optional second-order (S2)
   sensitivity indices.
6. Export results (CSV, YAML) and generate visualisations.

# Example
```julia
using AlphaPEM.Parametrisation.SobolSensitivityAnalysis

cfg = SobolAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    voltage_zone   = :full,
    N              = 1024,
    second_order   = false,
)

result = run_sobol_analysis(cfg)
```
"""
module SobolSensitivityAnalysis

using Printf
using Dates
using Statistics
using LinearAlgebra: BLAS
using DataFrames
using CSV
using ProgressMeter
using GlobalSensitivity
using QuasiMonteCarlo
using CairoMakie

using AlphaPEM.Config: SimulationConfig, PolarizationParams, NumericalParams, PhysicalParams, PARAMETER_METADATA,
                       OPERATING_CONDITIONS_BOUNDS, OperatingConditionConstraint, default_operating_condition_constraints
using AlphaPEM.Fuelcell: create_fuelcell
using AlphaPEM.Currents: create_current
using AlphaPEM.Application: run_simulation
import AlphaPEM.Core.Models: AlphaPEM as AlphaPEMSimulator, simulate_model!, _polarization_points

using ..ParametrisationCommon: ParameterBound, ParameterBounds, bounds_for_fuel_cell, new_PhysicalParams_from_sample

export SobolAnalysisConfig,
       PolarizationRegion,
       SobolAnalysisResult,
       run_sobol_analysis,
       build_input_parameters,
       build_regions,
       region_auc,
       compute_regional_aucs,
       generate_sobol_design_matrices,
       is_valid_operating_conditions,
       plot_sobol_indices,
       plot_sobol_ranking

include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_types.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_regions.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_sampling.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_simulation.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_knn.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_analysis.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_export.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_plots.jl"))


"""
    run_sobol_analysis(cfg::SobolAnalysisConfig)

Run the complete Sobol sensitivity analysis pipeline.

Returns a `SobolAnalysisResult`.
"""
function run_sobol_analysis(cfg::SobolAnalysisConfig)::SobolAnalysisResult
    overall_start = time()
    mkpath(cfg.output_dir)

    @info "=" ^ 72
    @info "AlphaPEM — Sobol Global Sensitivity Analysis"
    @info "=" ^ 72
    @info "Fuel cell : $(cfg.fuel_cell_type)"
    @info "Voltage zone: $(cfg.voltage_zone)"
    @info "Sobol base size N: $(cfg.N)"
    @info "Second order: $(cfg.second_order)"
    @info "Region thresholds: $(cfg.region_thresholds) A/cm²"

    # Step 1: define input parameters and regions
    params = build_input_parameters(cfg)
    regions = build_regions(cfg.region_thresholds)
    @info "Number of input parameters: $(length(params))"

    # Step 2: generate Sobol design matrices
    @info "Generating Sobol design matrices..."
    A, B = generate_sobol_design_matrices(cfg, params)
    all_points = GlobalSensitivity.fuse_designs(A, B; second_order = cfg.second_order)
    n_evals = size(all_points, 2)
    @info "Total model evaluations required: $(n_evals)"

    # Step 3: run simulations
    @info "Running AlphaPEM simulations..."
    df_curves = run_sobol_simulations(cfg, params, all_points)

    # Step 4: KNN imputation for failed simulations
    n_imputed = impute_missing_curves!(df_curves, params; k = cfg.knn_k)
    @info "Imputed $(n_imputed) missing curves."

    # Step 5: compute Sobol indices
    @info "Computing Sobol indices..."
    sobol_indices = compute_sobol_indices(cfg, params, A, B, df_curves, regions)

    # Step 6: export results
    @info "Exporting results to $(cfg.output_dir)..."
    result = SobolAnalysisResult(
        cfg,
        regions,
        [p.name for p in params],
        sobol_indices,
        Dict{Symbol, String}(),
        time() - overall_start
    )
    output_files = export_sobol_results(result, cfg.output_dir)
    result = SobolAnalysisResult(
        cfg,
        regions,
        [p.name for p in params],
        sobol_indices,
        output_files,
        time() - overall_start
    )

    # Step 7: plots
    try
        plot_sobol_indices(result; index_type = :ST, top_k = 10)
        plot_sobol_ranking(result; index_type = :ST, top_k = 10)
    catch e
        @warn "Failed to generate plots: $e"
    end

    @info "Analysis complete in $(round(result.execution_time; digits=1)) s."
    return result
end

end # module SobolSensitivityAnalysis
