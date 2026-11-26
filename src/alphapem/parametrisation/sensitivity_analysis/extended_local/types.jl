# -*- coding: utf-8 -*-

"""
    SUPPORTED_OUTPUTS

Predicted outputs `ỹ` that can be analysed, as defined in equation 3 of the reference below. Only
the cell voltage is available for now; HFR and membrane water crossover are not exposed by AlphaPEM.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
const SUPPORTED_OUTPUTS = (:voltage,)

"""
    OPERATING_CONDITIONS_NAMES_FOR_PARAMETRISATION

Names of the independent operating conditions varied in the factorial library. Each condition takes its lower
bound, midpoint, and upper bound, forming the complete `3^6 = 729` combinations. `Pc_des` is
not included because it is derived from `Pa_des` and the configured pressure gap.
"""
const OPERATING_CONDITIONS_NAMES_FOR_PARAMETRISATION = (
    :T_des, :Pa_des, :Sa, :Sc, :Phi_a_des, :Phi_c_des
)
const DEFAULT_CURRENT_DENSITY_GRID = Float64[
    0.0, 0.03, 0.05, 0.1, 0.2, 0.6, 1.0, 1.2, 1.5, 1.7, 1.9, 2.1, 2.4, 2.7, 3.0
]

"""
    CurrentDensityZone(name, lower, upper; include_lower=true, include_upper=false)

Named current-density interval used to summarize sensitivities.

Bounds are in A/cm². Inclusion flags make adjacent zones unambiguous (for example, one zone can
include a shared boundary while the next excludes it).

# Fields
- `name::Symbol`: Stable name used in summary tables.
- `lower::Float64`: Lower current-density bound (A/cm²).
- `upper::Float64`: Upper current-density bound (A/cm²); `Inf` means unbounded.
- `include_lower::Bool`: Whether the lower bound is part of the zone.
- `include_upper::Bool`: Whether the upper bound is part of the zone.
"""
struct CurrentDensityZone
    name::Symbol
    lower::Float64
    upper::Float64
    include_lower::Bool
    include_upper::Bool
end

"""
    CurrentDensityZone(name, lower, upper; include_lower=true, include_upper=false)

Create a current-density interval, converting its bounds to `Float64`.

Parameters
----------
name : Symbol
    Zone name used in sensitivity summary rows.
lower, upper : Real
    Zone bounds in A/cm²; `upper` may be `Inf`.
include_lower, include_upper : Bool
    Whether each bound is included. Defaults are `true` and `false`, respectively.

Returns
-------
CurrentDensityZone
    Named interval; validity and overlaps are checked when the run configuration is validated.
"""
function CurrentDensityZone(
    name::Symbol,
    lower::Real,
    upper::Real;
    include_lower::Bool = true,
    include_upper::Bool = false,
)::CurrentDensityZone
    return CurrentDensityZone(name, Float64(lower), Float64(upper), include_lower, include_upper)
end

"""
    DEFAULT_CURRENT_DENSITY_ZONES

Default non-overlapping current-density zones in A/cm². The boundary 0.5 belongs to the low zone,
and 1.5 belongs to the medium zone.
"""
const DEFAULT_CURRENT_DENSITY_ZONES = CurrentDensityZone[
    CurrentDensityZone(:low_current_density_or_kinetic, 0.0, 0.5; include_upper = true),
    CurrentDensityZone(:medium_current_density_or_ohmic, 0.5, 1.5; include_lower = false, include_upper = true),
    CurrentDensityZone(:high_current_density_or_mass_transport, 1.5, Inf; include_lower = false),
]


"""
    ScreeningDataset(fuel_cell_type, year, voltage_zone)

One experimental polarisation curve used for the initial screening of the Sobol candidates,
e.g. `ScreeningDataset(:ZSW_Pa_1_61_Pc_1_41, 2024, :full)`.
"""
struct ScreeningDataset
    fuel_cell_type::Symbol
    year::Union{Int,Nothing}
    voltage_zone::Symbol
end


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

# Operating cases
- `operating_bounds_overrides::Dict{Symbol,Tuple{Float64,Float64,Symbol}}`: Partial overrides
  for `OPERATING_CONDITIONS_BOUNDS`. Unspecified bounds retain their defaults.
- `base_operating_conditions::Union{Nothing,OperatingConditions}`: Starting operating conditions
  for every generated case. If `nothing`, use the selected fuel cell's nominal conditions;
  otherwise use this value as the base. The six parameterized conditions are replaced by their
  design levels, while fields not parameterized (e.g. `y_H2_in`, `i_min_stoich`) retain the base
  values.
- `pressure_gap::Float64`: Exact `Pa_des - Pc_des` difference in pascals. Default `0.2e5`.
- `current_density_grid::Vector{Float64}`: Shared current-density grid in A/cm². Default is
  Goshtasbi's 15-point grid.
- `current_density_zones::Vector{CurrentDensityZone}`: Named, non-overlapping current-density
  intervals used for sensitivity summaries. The default zones are low `0–0.5`, medium `(0.5, 1.5]`,
  and high `>1.5 A/cm²`.
- `delta_t_ini::Float64`: Stabilisation time at the initial 1 A/cm² current (s).
- `v_load::Float64`: Loading rate for initialization and sweep ramps (A/m²/s).
- `delta_t_break::Float64`: Stabilisation time above 0.3 A/cm² (s).
- `delta_t_break_lower_0_3::Float64`: Stabilisation time up to 0.3 A/cm², excluding OCV (s).
- `delta_t_break_OCV::Float64`: Stabilisation time at OCV after returning from 1 A/cm² (s).
- `delta_t_measurement::Float64`: Trailing voltage averaging-window duration (s).

# Outputs
- `outputs::Vector{Symbol}`: Predicted outputs, among `SUPPORTED_OUTPUTS`. Default `[:voltage]`.

# Local sensitivity (equations 7-10 of the reference below)
- `n_levels::Int`: Number of equidistant levels `n_s` per parameter in `[0, 1]` (eq. 7). Default `7`.
- `collinearity_floor::Float64`: Lower bound of the sensitivity-norm product in the collinearity
  index. Default `1e-12`.
- `collinearity_threshold::Float64`: Collinearity above which a pair is highlighted. Default `0.707`.
- `sensitivity_threshold::Float64`: Minimum median sensitivity magnitude for a parameter to be considered influential.
   Default `0.1`.

# Parameter-subset selection
- `psi_th::Float64`: Collinearity-filter threshold: when the median collinearity of a pair exceeds
  it, the parameter with the smaller median sensitivity is removed. Default `0.9` (Goshtasbi et al.).

# Sampling and screening
- `n_candidates::Int`: Requested number of Sobol candidates. If it is not a power of 2, the design
  generator rounds it up to the next power of 2 and reports the effective count. All candidates are
  screened; the number accepted is whatever screening yields (no target). Default `2^9` (`512`).
- `screening_datasets::Vector{ScreeningDataset}`: Experimental polarisation curves used to screen the
  candidates. A candidate is accepted if it is within `screening_threshold_V` on every dataset.
  Default: `ZSW_nominal` (2024, `:full`).
- `include_nominal::Bool`: Add the nominal parameter set as a candidate. Default `false`.
- `screening_norm::Symbol`: `:linf` (article) or `:rmse` (marapendi). Default `:linf`.
- `screening_threshold_V::Float64`: Maximum voltage error (V) for a candidate to be accepted. Default `0.3`.
- `seed::Int`: Random seed for reproducibility. Default `42`.

# Identifiability (equations 14-15 of the reference below)
- `identifiability_matrix::Symbol`: `:single_point` (article), `:median` (marapendi) or
  `:all_points`. Default `:single_point`.
- `lambda_threshold::Float64`: Minimum Hessian eigenvalue for an identifiable subset. Default `1.0`.

# Execution
- `parallel::Bool`: Use multi-threading. Default `true`.
- `max_run_time_s::Float64`: Maximum runtime per simulation (s). Default `60.0`.
- `resume::Bool`: Resume from an existing checkpoint in `output_dir`. Default `false`.
- `checkpoint_frequency::Int`: Save a checkpoint every `n` completed work units (`0` disables). Default `10`.
- `output_dir::String`: Directory for results. Default `results/extended_local_sensitivity`.
- `save_raw::Bool`: Export the raw sensitivity blocks. Default `false`.
- `top_k::Int`: Number of top parameters shown in plots. Default `13`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
Base.@kwdef struct ExtendedLocalSensitivityConfig
    # Fuel cell and parameters
    fuel_cell_type::Symbol                 = :ZSW_nominal      # Fuel-cell model used for nominal properties and screening.
    year::Union{Int,Nothing}               = 2024              # Experimental-data year; `nothing` selects the default year.
    voltage_zone::Symbol                   = :full             # Polarization-curve region used by the analysis.
    nb_gc::Int                             = 5                 # Number of gas-channel nodes in each simulation.
    numerical_params::NumericalParams      = NumericalParams(nb_gc = nb_gc) # Solver tolerances, mesh sizes and runtime settings.
    parameter_subset::Union{Vector{Symbol},Nothing}  = nothing # Parameter names to analyse; `nothing` selects all eligible parameters.
    parameter_bounds::Union{ParameterBounds,Nothing} = nothing # Optional analysis bounds; `nothing` uses fuel-cell defaults.
    integer_policy::Symbol                 = :fix              # Fix integer parameters at nominal (`:fix`) or round perturbations (`:round`).
    # Factorial operating-case library
    operating_bounds_overrides::Dict{Symbol, Tuple{Float64, Float64, Symbol}} =
        Dict{Symbol, Tuple{Float64, Float64, Symbol}}()        # Per-condition bound overrides; omitted conditions keep default bounds.
    base_operating_conditions::Union{Nothing, OperatingConditions} = nothing # Base for cases: nominal if `nothing`; swept fields are replaced, others retained.
    pressure_gap::Float64                   = 0.2e5            # Enforce Pa_des - Pc_des = this pressure difference (Pa).
    current_density_grid::Vector{Float64}   = copy(DEFAULT_CURRENT_DENSITY_GRID) # Requested current points (A/cm²), shared by every case.
    current_density_zones::Vector{CurrentDensityZone} = copy(DEFAULT_CURRENT_DENSITY_ZONES) # Named, non-overlapping current intervals for summaries.
    # Ascending current profile (defaults: AlphaPEM polarization durations)
    delta_t_ini::Float64                    = AscendingPolarizationParams().delta_t_ini             # Stabilization duration at the initial 1 A/cm² current (s).
    v_load::Float64                         = AscendingPolarizationParams().v_load                  # Current ramp rate (A/m²/s).
    delta_t_break::Float64                  = AscendingPolarizationParams().delta_t_break           # Stabilization duration above 0.3 A/cm² (s).
    delta_t_break_lower_0_3::Float64        = AscendingPolarizationParams().delta_t_break_lower_0_3 # Stabilization duration at or below 0.3 A/cm² (s).
    delta_t_break_OCV::Float64              = AscendingPolarizationParams().delta_t_break_OCV       # Stabilization duration at OCV after returning from 1 A/cm² (s).
    delta_t_measurement::Float64            = AscendingPolarizationParams().delta_t_measurement     # Trailing voltage averaging-window duration (s).
    # Outputs
    outputs::Vector{Symbol}                = [:voltage]   # Model outputs to extract; currently only voltage is supported.
    # Local sensitivity
    n_levels::Int                          = 7             # Equidistant parameter levels used to estimate each local sensitivity.
    collinearity_floor::Float64            = 1e-12         # Minimum norm product used in collinearity calculations.
    collinearity_threshold::Float64        = 0.707         # Collinearity threshold for highlighting parameter pairs.
    sensitivity_threshold::Float64         = 0.1           # Minimum median sensitivity magnitude for an influential parameter.
    # Parameter-subset selection
    psi_th::Float64                        = 0.9           # Pairs with median collinearity above it lose their less sensitive parameter.
    # Sampling and screening
    n_candidates::Int                      = 2^9           # Sobol candidates generated before experimental screening.
    screening_datasets::Vector{ScreeningDataset} = [ScreeningDataset(:ZSW_nominal, 2024, :full)] # Experimental curves every candidate must match.
    include_nominal::Bool                  = false         # Add the nominal parameter set as a screening candidate.
    screening_norm::Symbol                 = :linf         # Voltage-error metric used during screening (`:linf` or `:rmse`).
    screening_threshold_V::Float64         = 0.3           # Maximum allowed screening error (V).
    seed::Int                              = 42            # Seed for repeatable Sobol candidate generation.
    # Identifiability
    identifiability_matrix::Symbol         = :single_point # Sensitivity matrix assembly mode for identifiability.
    lambda_threshold::Float64              = 1.0           # Minimum Hessian eigenvalue for an identifiable subset.
    # Execution
    parallel::Bool                         = true          # Run independent work units concurrently when multiple Julia threads exist.
    max_run_time_s::Float64                = 120.0         # Wall-clock limit for one simulation run (s).
    resume::Bool                           = false         # Resume analysis from a compatible checkpoint in `output_dir`.
    checkpoint_frequency::Int              = 10            # Save progress every this many work units; 0 disables periodic saves.
    output_dir::String                     = abspath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "results", "extended_local_sensitivity")) # Directory for checkpoints and analysis outputs.
    save_raw::Bool                         = false         # Export raw per-point sensitivity blocks in addition to summaries.
    top_k::Int                             = 13            # Number of leading parameters shown in ranking plots.
end


"""
    validate_config(cfg::ExtendedLocalSensitivityConfig)

Check the consistency of an extended local sensitivity configuration.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Configuration to validate.

Returns
-------
ExtendedLocalSensitivityConfig
    The unchanged configuration if valid.

Raises
------
ArgumentError
    If a configuration field is invalid.
"""
function validate_config(cfg::ExtendedLocalSensitivityConfig)::ExtendedLocalSensitivityConfig
    _check(cond, msg) = cond || throw(ArgumentError(msg))

    # Validate the selected fuel-cell settings and outputs
    _check(cfg.voltage_zone in (:full, :before_voltage_drop),
           "voltage_zone must be :full or :before_voltage_drop (got $(cfg.voltage_zone))")
    _check(cfg.integer_policy in (:fix, :round), "integer_policy must be :fix or :round")
    _check(!isempty(cfg.outputs) && all(o -> o in SUPPORTED_OUTPUTS, cfg.outputs),
           "outputs must be a non-empty subset of $(SUPPORTED_OUTPUTS)")

    # Validate factorial operating-case settings
    _check(isfinite(cfg.pressure_gap) && cfg.pressure_gap > 0.0,
           "pressure_gap must be finite and positive")
    _check(length(cfg.current_density_grid) >= 2,
           "current_density_grid must contain at least two points")
    _check(all(isfinite, cfg.current_density_grid) &&
           all(i -> i >= 0.0, cfg.current_density_grid),
           "current_density_grid must contain finite, non-negative values")
    _check(all(diff(cfg.current_density_grid) .> 0.0),
           "current_density_grid must be strictly increasing")
    _check(!isempty(cfg.current_density_zones), "current_density_zones must not be empty")
    zone_names = [zone.name for zone in cfg.current_density_zones]
    _check(length(unique(zone_names)) == length(zone_names),
           "current_density_zones must have unique names")
    for zone in cfg.current_density_zones
        _check(zone.name != :all_current_densities,
               ":all_current_densities is reserved for the full-grid summary")
        _check(isfinite(zone.lower) && zone.lower >= 0.0 && !isnan(zone.upper) &&
               zone.upper > zone.lower && zone.upper > 0.0,
               "current-density zone $(zone.name) must have lower >= 0 and upper > lower")
    end
    ordered_zones = sort(cfg.current_density_zones; by = zone -> (zone.lower, zone.upper))
    for (previous, next) in zip(ordered_zones, ordered_zones[2:end])
        overlaps_interior = next.lower < previous.upper
        overlaps_at_boundary = next.lower == previous.upper &&
                               previous.include_upper && next.include_lower
        _check(!overlaps_interior && !overlaps_at_boundary,
               "current-density zones $(previous.name) and $(next.name) must not overlap")
    end

    # Validate the ascending current profile
    _check(all(isfinite, (cfg.delta_t_ini, cfg.v_load, cfg.delta_t_break,
                          cfg.delta_t_break_lower_0_3, cfg.delta_t_break_OCV,
                          cfg.delta_t_measurement)),
           "ascending profile durations and loading rate must be finite")
    _check(cfg.v_load > 0.0, "v_load must be positive")
    _check(cfg.delta_t_measurement >= 0.0, "delta_t_measurement must be non-negative")
    _check(min(cfg.delta_t_break, cfg.delta_t_break_lower_0_3) >= cfg.delta_t_measurement,
           "stabilisation times must be at least delta_t_measurement")
    _check(cfg.delta_t_break_OCV >= 0.0,
           "delta_t_break_OCV must be non-negative")
    _check(!iszero(first(cfg.current_density_grid)) ||
           cfg.delta_t_break_OCV >= cfg.delta_t_measurement,
           "delta_t_break_OCV must be at least delta_t_measurement when OCV is included in current_density_grid")
    for (name, bound) in cfg.operating_bounds_overrides
        _check(haskey(OPERATING_CONDITIONS_BOUNDS, name),
               "operating_bounds_overrides contains unknown operating condition $name")
        lo, hi, kind = bound
        _check(isfinite(lo) && isfinite(hi) && lo < hi,
               "operating bounds for $name must be finite and have min < max")
        _check(kind == OPERATING_CONDITIONS_BOUNDS[name][3],
               "operating bound type for $name must remain $(OPERATING_CONDITIONS_BOUNDS[name][3])")
    end

    # Validate local-sensitivity settings
    _check(cfg.n_levels >= 3, "n_levels must be at least 3")
    _check(cfg.collinearity_floor > 0.0, "collinearity_floor must be positive")
    _check(0.0 < cfg.collinearity_threshold < 1.0, "collinearity_threshold must lie in (0, 1)")
    _check(isfinite(cfg.sensitivity_threshold) && cfg.sensitivity_threshold >= 0.0,
           "sensitivity_threshold must be finite and non-negative")
    _check(0.0 < cfg.psi_th < 1.0, "psi_th must lie in (0, 1)")

    # Validate sampling and screening settings
    _check(cfg.n_candidates >= 1, "n_candidates must be positive")
    _check(cfg.screening_norm in (:linf, :rmse), "screening_norm must be :linf or :rmse")
    _check(cfg.screening_threshold_V > 0.0, "screening_threshold_V must be positive")
    _check(!isempty(cfg.screening_datasets), "screening_datasets must not be empty")
    for d in cfg.screening_datasets
        _check(d.voltage_zone in (:full, :before_voltage_drop),
               "screening dataset $(d.fuel_cell_type): voltage_zone must be :full or :before_voltage_drop")
    end

    # Validate identifiability and execution limits
    _check(cfg.identifiability_matrix in (:single_point, :median, :all_points),
           "identifiability_matrix must be :single_point, :median or :all_points")
    _check(cfg.lambda_threshold > 0.0, "lambda_threshold must be positive")
    _check(cfg.max_run_time_s > 0.0, "max_run_time_s must be positive")
    _check(cfg.checkpoint_frequency >= 0, "checkpoint_frequency must be non-negative")
    _check(cfg.top_k >= 1, "top_k must be positive")

    # Validate the optional parameter subset
    _check(cfg.parameter_subset === nothing || !isempty(cfg.parameter_subset),
           "parameter_subset must not be empty")
    return cfg
end


"""
    OutputSpec

Description of one predicted output `ỹ` (equation 3 of the reference below).

# Fields
- `name::Symbol`: Output name (see `SUPPORTED_OUTPUTS`).
- `extract::Function`: Maps a finished simulation to a fixed-size `Vector{Float64}`.
- `mask::Function`: Selects the current-density points kept for this output.
- `weight::Float64`: Weight applied when outputs are concatenated (`1.0` = article).

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
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
