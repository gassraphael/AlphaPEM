# -*- coding: utf-8 -*-

"""
    SobolAnalysisConfig

Configuration for the Sobol global sensitivity analysis of AlphaPEM.

# Fields
- `fuel_cell_type::Symbol`: Fuel-cell type to analyse. Default `:ZSW_nominal`.
- `voltage_zone::Symbol`: `:before_voltage_drop` or `:full`. Default `:full`.
- `N::Int`: Base size of the Sobol sequence. Default `10_000`.
- `second_order::Bool`: Compute second-order indices S2. Default `true`.
- `seed::Int`: Random seed for reproducibility. Default `42`.
- `region_thresholds::Tuple{Float64, Float64}`: Current-density thresholds (A/cm²)
  separating activation / ohmic / mass-transport regions. Default `(0.4, 1.6)`.
- `include_operating_conditions::Bool`: Include operating conditions as inputs.
  Default `true`.
- `parameter_bounds::Union{ParameterBounds, Nothing}`: Optional custom bounds.
  If `nothing`, bounds are inferred from the fuel cell.
- `polarization_params::PolarizationParams`: Current profile used for simulations.
- `numerical_params::NumericalParams`: Numerical settings for batch runs.
- `parallel::Bool`: Use multi-threading. Default `true`.
- `max_run_time_s::Float64`: Maximum simulation runtime per curve (s). Default `60.0`.
- `knn_k::Int`: Number of neighbours for KNN imputation. Default `10`.
- `top_k::Int`: Number of top parameters shown in plots. Default `13`.
- `output_dir::String`: Directory for results. Default `"results/sobol_sensitivity"`.
- `operating_condition_constraints::Vector{OperatingConditionConstraint}`: Constraints
  applied after sampling operating conditions. Supports `:(=)`, `:(>=)` and `:(<=)`
  constraints. Default: `Pc_des >= Pa_des - 0.5e5`.
- `excluded_operating_conditions::Vector{Symbol}`: Operating conditions that are kept at
  their nominal value and excluded from the Sobol sampling. Default: `Symbol[]`.
- `save_curves::Bool`: Save raw polarization curves. Default `true`.
"""
Base.@kwdef struct SobolAnalysisConfig
    fuel_cell_type::Symbol              = :ZSW_nominal
    year::Union{Int,Nothing}            = 2024
    voltage_zone::Symbol                = :full
    N::Int                              = 10_000
    second_order::Bool                  = true
    seed::Int                           = 42
    region_thresholds::Tuple{Float64, Float64} = (0.4, 1.6)
    include_operating_conditions::Bool  = true
    parameter_bounds::Union{ParameterBounds, Nothing} = nothing
    operating_condition_constraints::Vector{OperatingConditionConstraint} = default_operating_condition_constraints()
    excluded_operating_conditions::Vector{Symbol} = Symbol[]
    polarization_params::PolarizationParams = PolarizationParams(di_step = 0.05e4)
    numerical_params::NumericalParams   = NumericalParams(nb_gc = 1)
    parallel::Bool                      = true
    max_run_time_s::Float64             = 60.0
    knn_k::Int                          = 10
    top_k::Int                          = 13
    output_dir::String                  = abspath(joinpath(@__DIR__, "..", "..", "..", "..", "results", "sobol_sensitivity"))
    save_curves::Bool                   = true
end


"""
    PolarizationRegion

Named tuple-like container for a characteristic region of the polarization curve.

# Fields
- `name::Symbol`: Region name (`:activation`, `:ohmic`, `:mass_transport`).
- `i_min::Float64`: Lower current-density bound (A/cm²).
- `i_max::Float64`: Upper current-density bound (A/cm²).
"""
struct PolarizationRegion
    name::Symbol
    i_min::Float64
    i_max::Float64
end


"""
    SobolAnalysisResult

Complete output of a Sobol sensitivity analysis run.

# Fields
- `config::SobolAnalysisConfig`: Configuration used.
- `regions::Vector{PolarizationRegion}`: Analysed regions.
- `param_names::Vector{Symbol}`: Input parameter names in order.
- `sobol_indices::Dict{Symbol, DataFrame}`: Maps region name to a DataFrame of
  S1, ST (and S2 if requested) with confidence intervals.
- `output_files::Dict{Symbol, String}`: Paths to generated files.
- `execution_time::Float64`: Total wall-clock time (s).
- `imputation_report::Dict{Symbol, Any}`: Report on KNN imputation (total samples,
  failed samples, imputed samples, imputation rate, list of imputed sample ids).
"""
struct SobolAnalysisResult
    config::SobolAnalysisConfig
    regions::Vector{PolarizationRegion}
    param_names::Vector{Symbol}
    sobol_indices::Dict{Symbol, DataFrame}
    output_files::Dict{Symbol, String}
    execution_time::Float64
    imputation_report::Dict{Symbol, Any}
end
