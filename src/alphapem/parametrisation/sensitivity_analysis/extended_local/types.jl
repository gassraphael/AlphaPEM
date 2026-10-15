# -*- coding: utf-8 -*-

"""
    SUPPORTED_OUTPUTS

Predicted outputs `ỹ` that can be analysed (Goshtasbi et al. 2020, eq. 3). Only the cell
voltage is available for now; HFR and membrane water crossover are not exposed by AlphaPEM yet.
"""
const SUPPORTED_OUTPUTS = (:voltage,)


"""
    ExtendedLocalSensitivityConfig

Configuration of the extended local sensitivity and identifiability analysis.

# Fuel cell and parameters
- `fuel_cell_type::Symbol`: Fuel-cell type to analyse. Default `:ZSW_nominal`.
- `year::Union{Int,Nothing}`: Data year. Default `2024`.
- `voltage_zone::Symbol`: `:before_voltage_drop` or `:full`. Default `:full`.
- `nb_gc::Int`: Number of gas-channel nodes used for the runs. Default `5`.
- `numerical_params::NumericalParams`: Other numerical settings for the runs.
- `parameter_subset::Union{Vector{Symbol},Nothing}`: Restrict the analysis to these parameters. If `nothing`,
   all parameters are analysed from the fuel cell default configuration.
- `parameter_bounds::Union{ParameterBounds,Nothing}`: Custom bounds. If `nothing`, bounds are
  inferred from the fuel cell default configuration.
- `integer_policy::Symbol`: `:fix` keeps integer parameters at their nominal value (excluded
  from the analysis), `:round` rounds the perturbed levels. Default `:fix`.

# Outputs
- `outputs::Vector{Symbol}`: Predicted outputs, among `SUPPORTED_OUTPUTS`. Default `[:voltage]`.

# Local sensitivity (Goshtasbi et al. 2020, eq. 7-10)
- `fd_scheme::Symbol`: `:equidistant` (`n_levels` equidistant levels, article) or `:forward`
  (single forward step, marapendi). Default `:equidistant`.
- `n_levels::Int`: Number of equidistant levels per parameter. Default `7`.
- `eps_theta::Float64`: Forward step in the unit hypercube. Default `1e-2` (must stay above
  the solver tolerance noise).
- `collinearity_floor::Float64`: Lower bound of the sensitivity-norm product in the collinearity
  index. Default `1e-12`.
- `collinearity_threshold::Float64`: Collinearity above which a pair is highlighted. Default `0.707`.
- `sensitivity_threshold::Float64`: Sensitivity above which a parameter is deemed identifiable. Default `0.1`.

# Sampling and screening
- `n_candidates::Int`: Number of candidate Sobol points. Default `256`.
- `n_target::Int`: Number of accepted points after which screening stops. Default `55`.
- `scramble::Bool`: Randomise the Sobol sequence with `seed`. Default `true`.
- `include_nominal::Bool`: Add the nominal parameter set as a candidate. Default `false`.
- `screening_norm::Symbol`: `:linf` (article) or `:rmse` (marapendi). Default `:linf`.
- `screening_threshold_V::Float64`: Maximum voltage error (V) for a candidate to be accepted. Default `0.3`.
- `seed::Int`: Random seed for reproducibility. Default `42`.

# Identifiability (Goshtasbi et al. 2020, eq. 14-15)
- `identifiability_matrix::Symbol`: `:single_point` (article), `:median` (marapendi) or
  `:all_points`. Default `:single_point`.
- `lambda_threshold::Float64`: Minimum Hessian eigenvalue for an identifiable subset. Default `1.0`.

# Execution
- `parallel::Bool`: Use multi-threading. Default `true`.
- `max_run_time_s::Float64`: Maximum runtime per simulation (s). Default `60.0`.
- `max_simulations::Int`: Budget guard: the run is refused above this number of simulations
  unless `force` is `true`. Default `500_000`.
- `force::Bool`: Bypass the budget guard. Default `false`.
- `resume::Bool`: Resume from an existing checkpoint in `output_dir`. Default `false`.
- `checkpoint_frequency::Int`: Save a checkpoint every `n` completed work units (`0` disables). Default `10`.
- `output_dir::String`: Directory for results. Default `results/extended_local_sensitivity`.
- `save_raw::Bool`: Export the raw sensitivity blocks. Default `false`.
- `top_k::Int`: Number of top parameters shown in plots. Default `13`.
"""
Base.@kwdef struct ExtendedLocalSensitivityConfig
    # Fuel cell and parameters
    fuel_cell_type::Symbol                 = :ZSW_nominal
    year::Union{Int,Nothing}               = 2024
    voltage_zone::Symbol                   = :full
    nb_gc::Int                             = 5
    numerical_params::NumericalParams      = NumericalParams(nb_gc = nb_gc)
    parameter_subset::Union{Vector{Symbol},Nothing}  = nothing
    parameter_bounds::Union{ParameterBounds,Nothing} = nothing
    integer_policy::Symbol                 = :fix
    # Outputs
    outputs::Vector{Symbol}                = [:voltage]
    # Local sensitivity
    fd_scheme::Symbol                      = :equidistant
    n_levels::Int                          = 7
    eps_theta::Float64                     = 1e-2
    collinearity_floor::Float64            = 1e-12
    collinearity_threshold::Float64        = 0.707
    sensitivity_threshold::Float64         = 0.1
    # Sampling and screening
    n_candidates::Int                      = 256
    n_target::Int                          = 55
    scramble::Bool                         = true
    include_nominal::Bool                  = false
    screening_norm::Symbol                 = :linf
    screening_threshold_V::Float64         = 0.3
    seed::Int                              = 42
    # Identifiability
    identifiability_matrix::Symbol         = :single_point
    lambda_threshold::Float64              = 1.0
    # Execution
    parallel::Bool                         = true
    max_run_time_s::Float64                = 60.0
    max_simulations::Int                   = 500_000
    force::Bool                            = false
    resume::Bool                           = false
    checkpoint_frequency::Int              = 10
    output_dir::String                     = abspath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "results", "extended_local_sensitivity"))
    save_raw::Bool                         = false
    top_k::Int                             = 13
end


"""
    validate_config(cfg::ExtendedLocalSensitivityConfig)

Check the consistency of a configuration and return `cfg`.

Throws an `ArgumentError` describing the first invalid field.
"""
function validate_config(cfg::ExtendedLocalSensitivityConfig)::ExtendedLocalSensitivityConfig
    _check(cond, msg) = cond || throw(ArgumentError(msg))
    _check(cfg.voltage_zone in (:full, :before_voltage_drop),
           "voltage_zone must be :full or :before_voltage_drop (got $(cfg.voltage_zone))")
    _check(cfg.integer_policy in (:fix, :round), "integer_policy must be :fix or :round")
    _check(!isempty(cfg.outputs) && all(o -> o in SUPPORTED_OUTPUTS, cfg.outputs),
           "outputs must be a non-empty subset of $(SUPPORTED_OUTPUTS)")
    _check(cfg.fd_scheme in (:equidistant, :forward), "fd_scheme must be :equidistant or :forward")
    _check(cfg.n_levels >= 3, "n_levels must be at least 3")
    _check(0.0 < cfg.eps_theta < 1.0, "eps_theta must lie in (0, 1)")
    _check(cfg.collinearity_floor > 0.0, "collinearity_floor must be positive")
    _check(0.0 < cfg.collinearity_threshold < 1.0, "collinearity_threshold must lie in (0, 1)")
    _check(cfg.sensitivity_threshold >= 0.0, "sensitivity_threshold must be non-negative")
    _check(cfg.n_candidates >= 1, "n_candidates must be positive")
    _check(1 <= cfg.n_target <= cfg.n_candidates, "n_target must lie in [1, n_candidates]")
    _check(cfg.screening_norm in (:linf, :rmse), "screening_norm must be :linf or :rmse")
    _check(cfg.screening_threshold_V > 0.0, "screening_threshold_V must be positive")
    _check(cfg.identifiability_matrix in (:single_point, :median, :all_points),
           "identifiability_matrix must be :single_point, :median or :all_points")
    _check(cfg.lambda_threshold > 0.0, "lambda_threshold must be positive")
    _check(cfg.max_run_time_s > 0.0, "max_run_time_s must be positive")
    _check(cfg.max_simulations >= 1, "max_simulations must be positive")
    _check(cfg.checkpoint_frequency >= 0, "checkpoint_frequency must be non-negative")
    _check(cfg.top_k >= 1, "top_k must be positive")
    _check(cfg.parameter_subset === nothing || !isempty(cfg.parameter_subset),
           "parameter_subset must not be empty")
    return cfg
end


"""
    OutputSpec

Description of one predicted output `ỹ` (Goshtasbi et al. 2020, eq. 3).

# Fields
- `name::Symbol`: Output name (see `SUPPORTED_OUTPUTS`).
- `extract::Function`: Maps a finished simulation to a fixed-size `Vector{Float64}`.
- `mask::Function`: Selects the current-density points kept for this output.
- `weight::Float64`: Weight applied when outputs are concatenated (`1.0` = article).
"""
struct OutputSpec
    name::Symbol
    extract::Function
    mask::Function
    weight::Float64
end


"""
    ExtendedLocalSensitivityResult

Complete output of an extended local sensitivity analysis run.

# Fields
- `config::ExtendedLocalSensitivityConfig`: Configuration used.
- `param_names::Vector{Symbol}`: Analysed parameters, in order.
- `case_labels::Vector{String}`: Operating cases, in order.
- `theta::Matrix{Float64}`: Accepted sample points in the unit hypercube (`n_params × n_points`).
- `S_n::Array{Float64,3}`: Normalised sensitivities (`n_points × n_params × n_measurements`).
- `output_columns::Dict{Symbol,UnitRange{Int}}`: Measurement columns of `S_n` for each output.
- `sensitivity_median::Vector{Float64}`: Median sensitivity magnitude `s̄_i` (eq. 8).
- `sensitivity_std::Vector{Float64}`: Standard deviation `s*_i` of the magnitudes (eq. 9).
- `collinearity_median::Matrix{Float64}`: Median collinearity index `ψ_ij` (eq. 10).
- `ranking::Dict{Symbol,Vector{Int}}`: Parameter rankings (`:magnitude`, `:orthogonal`).
- `min_eigenvalues::Dict{Symbol,Vector{Float64}}`: Smallest Hessian eigenvalue per subset size.
- `selected_parameters::Vector{Symbol}`: Identifiable subset retained.
- `screening_report::Dict{Symbol,Any}`: Candidates, acceptance rate and failures.
- `output_files::Dict{Symbol,String}`: Paths to generated files.
- `execution_time::Float64`: Total wall-clock time (s).
"""
struct ExtendedLocalSensitivityResult
    config::ExtendedLocalSensitivityConfig
    param_names::Vector{Symbol}
    case_labels::Vector{String}
    theta::Matrix{Float64}
    S_n::Array{Float64,3}
    output_columns::Dict{Symbol,UnitRange{Int}}
    sensitivity_median::Vector{Float64}
    sensitivity_std::Vector{Float64}
    collinearity_median::Matrix{Float64}
    ranking::Dict{Symbol,Vector{Int}}
    min_eigenvalues::Dict{Symbol,Vector{Float64}}
    selected_parameters::Vector{Symbol}
    screening_report::Dict{Symbol,Any}
    output_files::Dict{Symbol,String}
    execution_time::Float64
end
