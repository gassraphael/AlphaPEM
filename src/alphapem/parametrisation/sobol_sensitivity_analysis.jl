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
    fuel_cell_type = :ZSW_nominal,
    year           = 2024,
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
using Logging
using DataFrames
using CSV
using ProgressMeter
using GlobalSensitivity
using QuasiMonteCarlo
using Random
using SHA
using JLD2
using CairoMakie

using AlphaPEM.Config: SimulationConfig, PolarizationParams, NumericalParams, PhysicalParams, PARAMETER_METADATA,
                       OPERATING_CONDITIONS_BOUNDS, OperatingConditionConstraint, default_operating_condition_constraints
using AlphaPEM.Fuelcell: create_fuelcell
using AlphaPEM.Currents: create_current
using AlphaPEM.Application: run_simulation
import AlphaPEM.Core.Models: AlphaPEM as AlphaPEMSimulator, simulate_model!, _polarization_points

using ..ParametrisationCommon: ParameterBound, ParameterBounds, bounds_for_fuel_cell, new_PhysicalParams_from_sample
using ..ValidParameterRegion: ValidityCriteria

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
       plot_sobol_ranking,
       add_confidence_intervals,
       build_sobol_region_summary,
       build_sobol_summary_table,
       select_top_features,
       run_sobol_convergence_analysis,
       plot_sobol_index_convergence,
       plot_top_k_rankings_across_regions,
       plot_sobol_heatmap

include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_types.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_regions.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_sampling.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_analysis.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_checkpoint.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_simulation.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_knn.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_summary.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "sobol_convergence.jl"))
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

    # Step 2: try to resume from an existing checkpoint
    checkpoint = nothing
    elapsed_before = 0.0
    if cfg.resume
        checkpoint = _load_checkpoint(cfg.output_dir, cfg)
    end

    if checkpoint !== nothing
        A = checkpoint.A
        B = checkpoint.B
        params = checkpoint.params
        regions = checkpoint.regions
        df_existing = checkpoint.df_curves
        elapsed_before = checkpoint.elapsed_time
        @info "Resumed from checkpoint; $(nrow(df_existing)) sample(s) already simulated."
    else
        # Step 2b: generate fresh Sobol design matrices
        @info "Generating Sobol design matrices..."
        A, B = generate_sobol_design_matrices(cfg, params)
        df_existing = nothing
    end

    all_points = _fuse_designs_bootstrap(A, B; second_order = cfg.second_order, nboot = _sobol_nboot(cfg.N))
    n_evals = size(all_points, 2)
    @info "Total model evaluations required: $(n_evals)"

    # Save the initial checkpoint so that a fresh run can be resumed even if it
    # is interrupted before the first batch of simulations completes.
    if cfg.checkpoint_frequency > 0 && df_existing === nothing
        empty_df = DataFrame(
            sample_id = Int[],
            config_hash = String[],
            status = Symbol[],
            ifc = Vector{Union{Vector{Float64}, Nothing}}(),
            Ucell = Vector{Union{Vector{Float64}, Nothing}}(),
        )
        for p in params
            empty_df[!, p.name] = Float64[]
        end
        _save_checkpoint(cfg.output_dir, cfg, A, B, params, regions, empty_df, 0.0)
    end

    # Step 3: run simulations
    @info "Running AlphaPEM simulations..."
    function _on_checkpoint(df_partial)
        _save_checkpoint(
            cfg.output_dir, cfg, A, B, params, regions, df_partial,
            time() - overall_start + elapsed_before
        )
    end
    df_curves = run_sobol_simulations(
        cfg, params, all_points;
        df_existing = df_existing,
        on_checkpoint = _on_checkpoint
    )

    # Step 4: KNN imputation for failed simulations
    n_total = nrow(df_curves)
    failed_before = findall(df_curves.status .!= :ok)
    n_failed = length(failed_before)
    n_imputed = impute_missing_curves!(df_curves, params; k = cfg.knn_k)
    imputed_idx = findall(df_curves.status .== :imputed)
    imputation_report = Dict{Symbol, Any}(
        :n_total => n_total,
        :n_failed_before_imputation => n_failed,
        :n_imputed => n_imputed,
        :imputation_rate => n_total > 0 ? n_imputed / n_total : 0.0,
        :sample_id_imputed => df_curves.sample_id[imputed_idx],
        :sample_id_failed => df_curves.sample_id[failed_before],
    )
    @info "Imputed $(n_imputed)/$(n_failed) missing curves."

    # Step 5: build regional AUC output matrix and impute missing entries
    @info "Building regional AUC output matrix..."
    Y = build_output_matrix(df_curves, regions, all_points)
    auc_report = impute_missing_aucs!(Y, df_curves, params, regions; k = cfg.knn_k)
    merge!(imputation_report, auc_report)

    if any(isnan, Y)
        missing_samples = Int[col for col in 1:size(Y, 2) if any(isnan, Y[:, col])]
        error("Regional AUCs remain missing after KNN imputation for sample(s) $(missing_samples); cannot compute Sobol indices.")
    end

    @info "$(auc_report[:n_auc_imputed]) missing regional AUC entry(ies) imputed across $(auc_report[:n_auc_missing]) sample(s)."

    # Step 6: compute Sobol indices
    @info "Computing Sobol indices..."
    sobol_indices = compute_sobol_indices(cfg, params, A, B, Y, regions)

    # Clean up the checkpoint now that the simulations have succeeded
    _remove_checkpoint(cfg.output_dir)

    # Step 6: export results
    @info "Exporting results to $(cfg.output_dir)..."
    total_execution_time = time() - overall_start + elapsed_before
    result = SobolAnalysisResult(
        cfg,
        regions,
        [p.name for p in params],
        sobol_indices,
        Dict{Symbol, String}(),
        total_execution_time,
        imputation_report
    )
    output_files = export_sobol_results(result, cfg.output_dir; curves_df = df_curves)
    result = SobolAnalysisResult(
        cfg,
        regions,
        [p.name for p in params],
        sobol_indices,
        output_files,
        total_execution_time,
        imputation_report
    )

    # Step 7: plots
    try
        plot_sobol_indices(result; index_type = :ST, top_k = cfg.top_k)
        plot_sobol_ranking(result; index_type = :ST, top_k = cfg.top_k)
        plot_top_k_rankings_across_regions(result; index_type = :ST, top_k = cfg.top_k)
        if cfg.second_order
            plot_sobol_heatmap(result)
        end
    catch e
        @warn "Failed to generate plots: $e"
    end

    @info "Analysis complete in $(round(result.execution_time; digits=1)) s."
    return result
end

end # module SobolSensitivityAnalysis
