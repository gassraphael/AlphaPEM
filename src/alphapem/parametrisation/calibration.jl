# -*- coding: utf-8 -*-

"""
    Calibration

Core calibration engine for AlphaPEM undetermined parameters using Genetic Algorithms.
"""
module Calibration

using Evolutionary
using Random
using Dates
using Printf

using AlphaPEM.Config: PolarizationCalibrationParams, PhysicalParams
using AlphaPEM.Config: CalibrationConfig, CalibrationResult, GAConfig
using AlphaPEM.Fuelcell: create_fuelcell
using AlphaPEM.Currents: create_current
import AlphaPEM.Core.Models: AlphaPEM, simulate_model!

# Import sampling helpers
using ..ValidParameterRegion.ConfigurationSampling: bounds_for_fuel_cell, new_PhysicalParams_from_sample, get_reference_config
using ..ValidParameterRegion.ResultsExport: export_parameter_bounds

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
    best_parameters_ever = copy(lower_bounds) # Initialize tracker for the best discovered parameters
    best_fitness_ever = Inf # Initialize best fitness tracker with infinity
    lock_capture = ReentrantLock() # Define lock for thread-safe updates to shared best results

    function objective(gene_values::Vector{Float64}) # Objective function to minimize RMSE
        fit = CalibrationHelpers._fitness_function(gene_values, parameter_bounds, base_params, fuel_cells, current_profiles, cfg.simulation_configs) # Compute simulation error

        lock(lock_capture) do # Safely update the global best individual
            if fit < best_fitness_ever
                best_fitness_ever = fit
                best_parameters_ever .= gene_values
            end
        end

        return fit # Return the computed fitness value
    end

    # 4. Configure GA
    # Custom mutation: mutates exactly `mutation_genes` genes per individual with uniform random resampling
    num_genes = length(lower_bounds) # Determine the number of genes in the individual
    num_mutated_genes = clamp(ga_config.mutation_genes, 1, num_genes) # Restrict mutations within valid bounds

    function uniform_fixed_gene_mutation(recombinant::Vector{Float64}; rng=Random.default_rng()) # Custom mutation operator
        mutation_indices = randperm(rng, num_genes)[1:num_mutated_genes] # Select random genes to mutate
        for i in mutation_indices
            recombinant[i] = lower_bounds[i] + rand(rng) * (upper_bounds[i] - lower_bounds[i]) # Apply uniform random resampling
        end
        return recombinant # Return the mutated recombinant
    end

    genetic_algorithm_options = GA( # Instantiate the Genetic Algorithm operator
        populationSize = ga_config.pop_size, # Set population size
        selection = roulette, # Use roulette wheel selection
        crossover = singlepoint, # Use single-point crossover
        mutation = uniform_fixed_gene_mutation, # Use the custom fixed-gene mutation
        epsilon = ga_config.elitism # Apply elitism to preserve best individuals
    )

    # 5. Setup callback
    history = Float64[] # Initialize history buffer for RMSE tracking
    last_save_time = time() # Track the last checkpoint save time

    callback = (record) -> begin # Callback triggered at each GA generation
        current_best_rmse = record.value # Extract current best RMSE from record
        push!(history, current_best_rmse) # Log current RMSE to history

        generation = record.iteration # Get current generation index
        @info @sprintf("Generation %d/%d: Best RMSE = %.4f %%", generation, ga_config.num_generations, best_fitness_ever) # Log progress

        if (cfg.save_frequency > 0 && generation % cfg.save_frequency == 0) || (time() - last_save_time > 300) # Check for save triggers
            CalibrationHelpers._save_intermediate(best_parameters_ever, best_fitness_ever, generation, parameter_bounds, base_params, cfg.output_dir) # Save checkpoint
            last_save_time = time() # Reset last save timer
        end

        return current_best_rmse <= ga_config.target_error # Return true if target error is reached
    end

    # 6. Run optimization
    rng = ga_config.seed === nothing ? Random.default_rng() : MersenneTwister(ga_config.seed) # Initialize random number generator

    constraints = Evolutionary.BoxConstraints(lower_bounds, upper_bounds) # Define search space constraints

    initial_population = if cfg.initial_population_file !== nothing # Handle warm-start population loading
        loaded = CalibrationHelpers._load_warm_start_population(
            cfg.initial_population_file, parameter_bounds, ga_config.pop_size, lower_bounds, upper_bounds, rng) # Attempt to load file
        if loaded !== nothing
            @info "Warm-start: loaded $(length(loaded)) individuals from $(cfg.initial_population_file)" # Confirm loading
            loaded # Use loaded population
        else
            @warn "Failed to load warm-start, using random population" # Issue warning on failure
            [lower_bounds .+ rand(rng, length(lower_bounds)) .* (upper_bounds .- lower_bounds) for _ in 1:ga_config.pop_size] # Fallback to random
        end
    else
        [lower_bounds .+ rand(rng, length(lower_bounds)) .* (upper_bounds .- lower_bounds) for _ in 1:ga_config.pop_size] # Initialize random population
    end

    parallelization = cfg.parallel ? :thread : :serial # Determine parallelization mode
    if cfg.parallel
        @info "Calibration: parallel population evaluation enabled ($(Threads.nthreads()) threads)." # Log threading status
    else
        @info "Calibration: sequential population evaluation." # Log serial status
    end

    optimization_result = Evolutionary.optimize(objective, constraints, genetic_algorithm_options, initial_population, # Execute optimization
                 Evolutionary.Options(
                     iterations = ga_config.num_generations, # Set maximum iterations
                     callback = callback, # Set progress callback
                     rng = rng, # Set random number generator
                     parallelization = parallelization # Set evaluation mode
                 ))

    optimized_genes = best_parameters_ever # Extract best individual found (global best, not just last generation)
    best_fitness = best_fitness_ever # Extract best fitness achieved (global best)
    best_params = new_PhysicalParams_from_sample(optimized_genes, parameter_bounds, base_params) # Convert genes to physical parameters

    # 7. Save results
    final_bounds = Dict{Symbol, Tuple{Float64, Float64}}( # Prepare final results dictionary
        bound.name => (optimized_genes[i], optimized_genes[i]) for (i, bound) in enumerate(parameter_bounds.bounds)
    )
    export_parameter_bounds(final_bounds, joinpath(cfg.output_dir, "calibrated_bounds.yaml"); # Export calibrated parameters
                            method = :calibration,
                            metadata = Dict("fuel_cell_type" => string(ref_cfg.type_fuel_cell),
                                            "rmse_percent" => best_fitness))

    elapsed_time = time() - start_time # Compute total execution time
    @info "Calibration complete in $(round(elapsed_time/3600, digits=2)) hours." # Log total time
    @info "Final RMSE: $(round(best_fitness, digits=4)) %" # Log final accuracy

    result = CalibrationResult(cfg, best_params, best_fitness, best_fitness, history, elapsed_time) # Construct result object

    CalibrationHelpers._save_final_results(result, cfg.output_dir, [optimized_genes], [best_fitness]) # Persist final data to disk

    checkpoint_path = joinpath(cfg.output_dir, "calibration_checkpoint.yaml") # Define path to temporary checkpoint
    if isfile(checkpoint_path) # Check if checkpoint exists
        rm(checkpoint_path) # Remove temporary checkpoint file
        @info "Temporary checkpoint file removed." # Log removal
    end

    return result # Return the final calibration results
end

end # module Calibration
