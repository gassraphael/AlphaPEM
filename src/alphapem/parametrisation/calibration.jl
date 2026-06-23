# -*- coding: utf-8 -*-

"""
    Calibration

Core calibration engine for AlphaPEM undetermined parameters using Genetic Algorithms.
"""
module Calibration

using PyCall
using Random
using Dates
using Printf

using AlphaPEM.Config: PolarizationCalibrationParams, PhysicalParams
using AlphaPEM.Config: CalibrationConfig, CalibrationResult, GAConfig
using AlphaPEM.Fuelcell: create_fuelcell
using AlphaPEM.Currents: create_current
import AlphaPEM.Core.Models: AlphaPEM, simulate_model!

# Import common parametrisation helpers
using ..ParametrisationCommon: bounds_for_fuel_cell, new_PhysicalParams_from_sample,
                               get_reference_config, export_parameter_bounds, export_calibrated_params

export CalibrationConfig,
       CalibrationResult,
       calibrate

# Include helper module for internal functions
include("calibration/helpers.jl")
using .CalibrationHelpers

# ________________________________________________________________─────────────
# CORE CALIBRATION
# ________________________________________________________________─────────────

"""
    calibrate(cfg::CalibrationConfig) -> CalibrationResult

Launch the GA-based parameter calibration for AlphaPEM.
"""
function calibrate(cfg::CalibrationConfig)::CalibrationResult

    # 0. Initialisation
    start_time = time() # Record the initial time for performance measurement
    mkpath(cfg.output_dir) # Initialize the output directory for logs and results
    ga_config = cfg.ga_config # Extract GA hyperparameter settings

    # 1. Setup bounds and reference configuration
    # All conditions share the same fuel cell type → same parameter bounds
    ref_cfg = first(cfg.simulation_configs)
    parameter_bounds = bounds_for_fuel_cell(ref_cfg.type_fuel_cell, ref_cfg.voltage_zone) # Retrieve parameter constraints
    base_params = get_reference_config(ref_cfg.type_fuel_cell) # Load baseline physical parameters

    lower_bounds = [bound.min for bound in parameter_bounds.bounds] # Define the lower boundary vector for GA
    upper_bounds = [bound.max for bound in parameter_bounds.bounds] # Define the upper boundary vector for GA

    # 2. Prepare fuel cells and current profiles — one entry per experimental condition
    fuel_cells       = [create_fuelcell(sc.type_fuel_cell, sc.voltage_zone) for sc in cfg.simulation_configs]
    current_profiles = [create_current(PolarizationCalibrationParams(), fc) for fc in fuel_cells]

    # 3. Define fitness function
    # When parallel=true, we use a batch fitness function evaluated with Threads.@threads in Julia.
    # PyGAD's parallel_processing uses Python threads which cannot parallelize Julia code (GIL-equivalent);
    # instead, fitness_batch_size=pop_size sends the whole population at once to a Julia batch function.
    if cfg.parallel
        jl_fitness_func_batch = (ga_instance, population, population_idx) -> CalibrationHelpers._fitness_function_batch(
            ga_instance, population, population_idx, parameter_bounds, base_params, fuel_cells, current_profiles, cfg.simulation_configs
        )
        fitness_func = py"lambda f: lambda ga_instance, population, population_idx: list(f(ga_instance, population, population_idx))"(jl_fitness_func_batch)
    else
        jl_fitness_func = (ga_instance, solution, solution_idx) -> CalibrationHelpers._fitness_function(
            solution, parameter_bounds, base_params, fuel_cells, current_profiles, cfg.simulation_configs
        )
        # Wrap Julia function in a Python lambda to provide __code__ attribute for PyGAD inspection.
        fitness_func = py"lambda f: lambda ga_instance, solution, solution_idx: f(ga_instance, solution, solution_idx)"(jl_fitness_func)
    end

    # 4. Setup callback (on_generation)
    history = Float64[] # Initialize history buffer for RMSE tracking
    last_save_time = Ref(time()) # Track the last checkpoint save time

    # Wrap Julia callback in a Python lambda to provide __code__ attribute for PyGAD inspection.
    on_generation = py"lambda f: lambda ga_instance: f(ga_instance)"(
        ga -> CalibrationHelpers._on_generation(ga, history, ga_config, cfg, last_save_time, parameter_bounds, base_params)
    )

    # 5. Build initial population
    num_params = length(lower_bounds)
    rng = ga_config.seed === nothing ? Random.default_rng() : MersenneTwister(ga_config.seed) # Initialize random number generator

    initial_population = if cfg.initial_population_file !== nothing # Handle warm-start population loading
        loaded = CalibrationHelpers._load_warm_start_population(
            cfg.initial_population_file, parameter_bounds, ga_config.pop_size, lower_bounds, upper_bounds, rng) # Attempt to load file
        if loaded !== nothing
            @info "Warm-start: loaded $(length(loaded)) individuals from $(cfg.initial_population_file)" # Confirm loading
            loaded # Use loaded population
        else
            @warn "Failed to load warm-start, using random population" # Issue warning on failure
            [lower_bounds .+ rand(rng, num_params) .* (upper_bounds .- lower_bounds) for _ in 1:ga_config.pop_size] # Fallback to random
        end
    else
        [lower_bounds .+ rand(rng, num_params) .* (upper_bounds .- lower_bounds) for _ in 1:ga_config.pop_size] # Initialize random population
    end

    # Convert initial population to Python-compatible matrix (list of lists)
    initial_population_py = [collect(ind) for ind in initial_population]

    if cfg.parallel
        @info "Calibration: parallel population evaluation enabled ($(Threads.nthreads()) threads)." # Log threading status
    else
        @info "Calibration: sequential population evaluation." # Log serial status
    end

    # 6. Run optimization with PyGAD
    pygad = pyimport("pygad") # Import PyGAD

    # Gene space: list of dicts with 'low' and 'high' for each gene
    gene_space = [Dict("low" => lower_bounds[i], "high" => upper_bounds[i]) for i in 1:num_params]

    ga_instance = pygad.GA(
        num_generations        = ga_config.num_generations,    # Maximum number of generations
        sol_per_pop            = ga_config.pop_size,           # Population size
        num_parents_mating     = ga_config.num_parents_mating, # Parents selected per generation
        mutation_num_genes     = ga_config.mutation_num_genes, # Number of genes to mutate.
        keep_parents           = ga_config.elitism,            # Keep the best solution of the previous generation. 1 elite should be enough.
        parent_selection_type  = "rws",                        # Parent selection type. "rws" means roulette wheel selection. Best found.
        crossover_type         = "single_point",               # Crossover type. "single_point" means single point crossover. Best found.
        mutation_type          = "random",                     # Uniform random mutation
        fitness_func           = fitness_func,                 # Fitness function (maximized)
        num_genes              = num_params,                   # Number of parameters to optimize
        gene_space             = gene_space,                   # Per-gene bounds
        initial_population     = initial_population_py,        # Warm-start or random initial population
        on_generation          = on_generation,                # Callback called after each generation
        stop_criteria          = "reach_$(1.0 / ga_config.target_error)", # Stop if fitness = 1/RMSE exceeds 1/target_error
        random_seed            = ga_config.seed === nothing ? nothing : ga_config.seed, # Reproducibility seed
        fitness_batch_size     = cfg.parallel ? ga_config.pop_size : nothing # Batch mode: whole population sent at once for Julia-side threading
    )

    ga_instance.run() # Execute the genetic algorithm
    println() # Finalize progress line

    best_sol, best_fitness_py, _ = ga_instance.best_solution()
    best_fitness = 1.0 / best_fitness_py # Convert 1/RMSE back to RMSE
    optimized_genes = Float64.(best_sol) # Extract best individual found (global best)
    best_params = new_PhysicalParams_from_sample(optimized_genes, parameter_bounds, base_params) # Convert genes to physical parameters

    # 7. Save results
    final_params = Dict{Symbol, Float64}( # Prepare final calibrated parameters dictionary
        bound.name => optimized_genes[i] for (i, bound) in enumerate(parameter_bounds.bounds)
    )
    export_calibrated_params(final_params, joinpath(cfg.output_dir, "calibrated_bounds.yaml"); # Export calibrated parameters
                             method = :calibration,
                             metadata = Dict("fuel_cell_type" => string(ref_cfg.type_fuel_cell),
                                             "rmse_percent" => best_fitness))

    elapsed_time = time() - start_time # Compute total execution time
    result = CalibrationResult(cfg, best_params, best_fitness, best_fitness, history, elapsed_time) # Construct result object

    # Extract full final population and fitness values from PyGAD
    final_population_matrix = Float64.(ga_instance.population) # Shape: (pop_size, num_genes)
    final_fitness_values    = Float64.(ga_instance.last_generation_fitness) # Shape: (pop_size,)
    final_population_list   = [collect(final_population_matrix[i, :]) for i in 1:size(final_population_matrix, 1)]
    final_fitness_list      = [1.0 / v for v in final_fitness_values] # Convert 1/RMSE back to RMSE

    CalibrationHelpers._save_final_results(result, cfg.output_dir, final_population_list, final_fitness_list) # Persist final data to disk
    CalibrationHelpers._plot_calibration_results(result, cfg.output_dir) # Generate results figure

    for checkpoint_file in ("calibration_checkpoint.yaml", "calibration_checkpoint_population.yaml") # Clean up temporary checkpoint files
        checkpoint_path = joinpath(cfg.output_dir, checkpoint_file)
        if isfile(checkpoint_path) # Check if checkpoint exists
            rm(checkpoint_path) # Remove temporary checkpoint file
        end
    end

    return result # Return the final calibration results
end

end # module Calibration
