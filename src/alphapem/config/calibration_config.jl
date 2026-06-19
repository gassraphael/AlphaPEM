# -*- coding: utf-8 -*-

"""
    CalibrationConfig

Configuration structures for Genetic Algorithm-based parameter calibration.
"""

"""
    GAConfig

Genetic Algorithm hyperparameters.

# Fields
- `num_generations::Int`: Number of generations. Default `1000`.
- `pop_size::Int`: Population size. Default `128`.
- `mating_rate::Float64`: Fraction of population selected as parents. Default `0.2`.
- `mutation_genes::Int`: Number of genes mutated per individual (exactly). Default `1`.
- `elitism::Int`: Number of best individuals carried to next generation. Default `1`.
- `target_error::Float64`: Stop if RMSE < this value (%). Default `0.005` (0.5%).
- `seed::Union{Int, Nothing}`: Random seed for reproducibility. Set to an `Int` to fix results, or `nothing` for stochastic runs. Default `nothing`.
"""
Base.@kwdef struct GAConfig
    num_generations::Int              = 1000
    pop_size::Int                     = 128
    mating_rate::Float64              = 0.2
    mutation_genes::Int               = 1
    elitism::Int                      = 1
    target_error::Float64             = 0.5/100
    seed::Union{Int, Nothing}         = nothing  # nothing = stochastic
end


"""
    CalibrationConfig

Complete configuration for calibration workflow. Supports one or several experimental
conditions; all conditions must share the same fuel cell type and voltage zone, as they
calibrate the same physical parameter set.

# Fields
- `simulation_configs::Vector{SimulationConfig}`: One entry per experimental condition
  (e.g. different operating temperatures or pressures). Each carries its own
  `operating_conditions`; the first entry's `type_fuel_cell` and `voltage_zone` govern
  the parameter bounds for all conditions.
- `ga_config::GAConfig`: Genetic algorithm hyperparameters. Default `GAConfig()`.
- `parallel::Bool`: Use multi-threading for population evaluation. Default `true`.
- `initial_population_file::Union{String, Nothing}`: Path to YAML file for warm-start
  (checkpoint, calibrated_bounds, or final_population). Default `nothing`.
- `output_dir::String`: Directory for results and logs. Default `"results/calibration"`.
- `save_frequency::Int`: Frequency (generations) to save checkpoints. Default `1`.
"""
Base.@kwdef struct CalibrationConfig
    simulation_configs::Vector{SimulationConfig}   = [SimulationConfig(
                                                            type_current   = PolarizationCalibrationParams(),
                                                            type_display   = :no_display)]
    ga_config::GAConfig                            = GAConfig()
    parallel::Bool                                 = true
    initial_population_file::Union{String,Nothing} = nothing
    output_dir::String                             = "results/calibration"
    save_frequency::Int                            = 1
end


"""
    CalibrationResult

Results of a calibration run.

# Fields
- `config::CalibrationConfig`: Configuration used.
- `best_params::PhysicalParams`: Calibrated physical parameters.
- `best_fitness::Float64`: Final fitness (RMSE %).
- `min_rmse::Float64`: Minimum RMSE achieved (%).
- `history::Vector{Float64}`: Best RMSE per generation.
- `execution_time::Float64`: Total wall-clock time (seconds).
"""
struct CalibrationResult
    config::CalibrationConfig
    best_params::PhysicalParams
    best_fitness::Float64
    min_rmse::Float64
    history::Vector{Float64}
    execution_time::Float64
end
