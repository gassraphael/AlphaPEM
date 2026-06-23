# -*- coding: utf-8 -*-

"""
    CalibrationHelpers

Internal helper functions for the calibration module.
"""
module CalibrationHelpers

using Statistics: mean
using Dates
using Printf
using YAML
using CairoMakie

import AlphaPEM.Core.Models: AlphaPEM, simulate_model!, _polarization_points_cali, _calculate_rmse,
                             plot_polarization_curve_for_cali, plot_polarization_curve
import AlphaPEM.Config: PolaExperimentalData, PolarizationCalibrationParams, PolarizationParams, SimulationConfig
import AlphaPEM.Fuelcell: create_fuelcell
import AlphaPEM.Currents: create_current
using ...Calibration: CalibrationConfig
using ...ParametrisationCommon: export_parameter_bounds, export_calibrated_params, new_PhysicalParams_from_sample

export _fitness_function,
       _fitness_function_batch,
       _on_generation,
       calculate_simulation_error,
       _save_intermediate,
       _load_warm_start_population,
       _save_final_results,
       _plot_calibration_results

"""
    _fitness_function(solution, parameter_bounds, base_params, fuel_cells, current_profiles, simulation_configs) -> Float64

Evaluate the negative Root Mean Square Error (RMSE) for a given parameter gene vector across all operating conditions.
PyGAD maximizes fitness, so we return -RMSE.
"""
function _fitness_function(solution,
                           parameter_bounds,
                           base_params,
                           fuel_cells,
                           current_profiles,
                           simulation_configs)::Float64
    gene_values = Float64.(solution) # Convert PyObject to Julia Float64 vector
    try
        physical_params = new_PhysicalParams_from_sample(gene_values, parameter_bounds, base_params) # Map normalized genes to physical parameters
        rmse_values = zeros(Float64, length(fuel_cells)) # Initialize array for individual condition RMSEs

        for i in eachindex(fuel_cells) # Iterate through each fuel cell operating condition
            fuel_cell = deepcopy(fuel_cells[i]) # Create a local copy of the fuel cell model
            fuel_cell.physical_parameters = physical_params # Apply the current candidate parameters

            simulation = AlphaPEM(fuel_cell, current_profiles[i], simulation_configs[i]) # Initialize the simulation engine
            simulate_model!(simulation) # Execute the numerical simulation

            experimental_data = fuel_cell.pola_exp_data_cali # Retrieve the calibration experimental data
            if isempty(experimental_data.i_exp) || isempty(experimental_data.U_exp) # Check if experimental data is valid
                rmse_values[i] = 100.0 # Assign high error penalty for missing data
            else
                rmse_values[i] = calculate_simulation_error(simulation, experimental_data) # Compute RMSE for this specific condition
            end
        end

        return 1.0 / (mean(rmse_values) + 1e-6) # Return 1/RMSE (PyGAD maximizes, RWS requires positive fitness)

    catch e # Catch any simulation or numerical errors
        println("\nAn error occurred during the evaluation of the solution.")
        params = [string(b.name, ": ", gene_values[i]) for (i, b) in enumerate(parameter_bounds.bounds)]
        println("Attempted parameters: " * join(params, " | "))
        println("Exception : ", e)
        println("Refusing this solution and continuing the optimization.\n")
        return 1e-6 # Return a near-zero fitness value on failure (RWS-compatible)
    end
end

"""
    _fitness_function_batch(ga_instance, population, population_idx, parameter_bounds, base_params, fuel_cells, current_profiles, simulation_configs) -> Vector{Float64}

Batch version of the fitness function for parallel evaluation of the entire population.
Uses Threads.@threads to distribute individual evaluations.
"""
function _fitness_function_batch(ga_instance,
                                 population,
                                 population_idx,
                                 parameter_bounds,
                                 base_params,
                                 fuel_cells,
                                 current_profiles,
                                 simulation_configs)::Vector{Float64}
    n = size(population, 1)
    fitness_values = Vector{Float64}(undef, n)
    @sync for i in 1:n
        Threads.@spawn begin
            fitness_values[i] = _fitness_function(
                population[i, :], parameter_bounds, base_params, fuel_cells, current_profiles, simulation_configs
            )
        end
    end
    return fitness_values
end

"""
    _on_generation(ga_instance, history, ga_config, cfg, last_save_time, parameter_bounds, base_params)

Callback function executed after each generation of the Genetic Algorithm.
Handles progress logging and periodic checkpoint saving.
"""
function _on_generation(ga_instance, history, ga_config, cfg, last_save_time, parameter_bounds, base_params)
    # PyGAD maximizes fitness, so fitness = 1/RMSE (positive, compatible with RWS)
    best_sol, best_fitness_py, _ = ga_instance.best_solution()
    current_best_rmse =  1.0 / best_fitness_py # Convert 1/RMSE back to RMSE
    push!(history, current_best_rmse) # Log current RMSE to history

    generation = ga_instance.generations_completed # Get current generation index
    msg = @sprintf("Generation %d/%d: Best RMSE = %.2f %%", generation, ga_config.num_generations, current_best_rmse)
    print("\r[ Info: ", msg) # Log progress in-place
    flush(stdout)

    if (cfg.save_frequency > 0 && generation % cfg.save_frequency == 0) || (time() - last_save_time[] > 300) # Check for save triggers
        best_parameters = Float64.(best_sol)
        _save_intermediate(ga_instance, best_parameters, current_best_rmse, generation, parameter_bounds, base_params, cfg.output_dir) # Save checkpoint
        last_save_time[] = time() # Reset last save timer
    end
end

"""
    calculate_simulation_error(simulation, experimental_data) -> Float64

Calculate the RMSE (%) between the simulation results and the experimental polarization data.

For calibration protocols, extracts the simulated polarization points using the canonical sampling
infrastructure. No interpolation needed: the simulation exactly follows the experimental current profile,
so averaged simulation points correspond directly to experimental measurements.
"""
function calculate_simulation_error(simulation::AlphaPEM, experimental_data::PolaExperimentalData)::Float64
    _, Ucell_sim = _polarization_points_cali(simulation.outputs, simulation.current_density, experimental_data.i_exp)
    return _calculate_rmse(Ucell_sim, experimental_data.U_exp)
end

"""
    _save_intermediate(ga_instance, best_gene_values, best_fitness, generation, parameter_bounds, base_params, output_dir)

Save the current best individual's calibrated parameter values and the full population as checkpoints for potential warm-start recovery.
"""
function _save_intermediate(ga_instance, best_gene_values, best_fitness, generation, parameter_bounds, base_params, output_dir)
    best_params_dict = Dict{Symbol, Float64}( # Prepare calibrated parameters dictionary for export
        bound.name => best_gene_values[i] for (i, bound) in enumerate(parameter_bounds.bounds)
    )

    try
        export_calibrated_params(best_params_dict, joinpath(output_dir, "calibration_checkpoint.yaml"); # Export checkpoint to YAML
                                 method = :checkpoint, # Mark as an intermediate checkpoint
                                 metadata = Dict("generation" => generation, # Log the current generation index
                                                 "rmse_percent" => best_fitness)) # Log the current best accuracy
    catch e # Handle filesystem or export errors
        @warn "Failed to save checkpoint: $e" # Log warning message
    end

    # Save the full current population for warm-start recovery
    try
        population = Float64.(ga_instance.population) # Extract current population matrix from PyGAD
        fitness_values = Float64.(ga_instance.last_generation_fitness) # Extract last generation fitness values
        population_data = [
            Dict("individual" => i, "params" => collect(population[i, :]), "rmse" => -fitness_values[i])
            for i in 1:size(population, 1)
        ]
        YAML.write_file(joinpath(output_dir, "calibration_checkpoint_population.yaml"), population_data) # Write full population
    catch e
        @warn "Failed to save checkpoint population: $e"
    end
end

"""
    _load_warm_start_population(file, parameter_bounds, pop_size, lower_bounds, upper_bounds, rng) -> Union{Vector{Vector{Float64}}, Nothing}

Load a warm-start population from a YAML checkpoint or a final population file.
"""
function _load_warm_start_population(file::String, parameter_bounds, pop_size::Int, lower_bounds, upper_bounds, rng)
    try
        data = YAML.load_file(file) # Parse the YAML source file
        individuals = Float64[][] # Initialize container for loaded individuals

        if isa(data, Dict) && haskey(data, "parameters") # Handle checkpoint/calibrated_bounds format
            params_dict = data["parameters"] # Extract the parameter data dictionary
            param_index_map = Dict(string(bound.name) => i for (i, bound) in enumerate(parameter_bounds.bounds)) # Map names to vector indices

            individual = zeros(Float64, length(parameter_bounds.bounds)) # Initialize a single candidate vector
            for (param_name, param_data) in params_dict # Iterate through each parameter entry
                if haskey(param_index_map, param_name) # Check if parameter exists in current model
                    idx = param_index_map[param_name] # Retrieve the target index
                    # Support both new single-value format and legacy min/max format
                    if isa(param_data, Dict)
                        individual[idx] = Float64(param_data["min"]) # Legacy format: use min field
                    else
                        individual[idx] = Float64(param_data) # New format: direct scalar value
                    end
                else
                    @warn "Parameter $param_name from checkpoint not found in current bounds" # Warn if model mismatch
                end
            end
            push!(individuals, individual) # Add the reconstructed individual to the list

        elseif isa(data, Vector) # Handle final_population list format
            for entry in data # Iterate through list entries
                if isa(entry, Dict) && haskey(entry, "params") # Check for expected parameter key
                    push!(individuals, copy(entry["params"])) # Load the gene vector copy
                end
            end
        else # Handle unsupported formats
            @warn "Unrecognized file format in $file" # Log error for unknown structure
            return nothing # Return nothing to signal failure
        end

        if isempty(individuals) # Check if any individuals were successfully retrieved
            @warn "No individuals loaded from $file" # Log warning for empty results
            return nothing # Signal empty population
        end

        # Fill to pop_size with random individuals if the loaded set is smaller than required
        num_loaded_individuals = length(individuals) # Count currently loaded candidates
        if num_loaded_individuals < pop_size # Check if expansion is needed
            for _ in (num_loaded_individuals+1):pop_size # Fill remaining population slots
                random_individual = lower_bounds .+ rand(rng, length(lower_bounds)) .* (upper_bounds .- lower_bounds) # Generate random candidate
                push!(individuals, random_individual) # Append to population
            end
        elseif num_loaded_individuals > pop_size # Truncate if loaded set exceeds requested size
            individuals = individuals[1:pop_size] # Keep only the required number of candidates
        end

        return individuals # Return the finalized population vector

    catch e # Handle general loading errors
        @warn "Failed to load warm-start population from $file: $e" # Log failure details
        return nothing # Signal failure
    end
end

"""
    _save_final_results(result, output_dir, final_population, final_fitness)

Save the comprehensive calibration report and the final population gene vectors to YAML files.
"""
function _save_final_results(result, output_dir::String, final_population, final_fitness)
    report = Dict( # Construct the final report dictionary
        "metadata" => Dict(
            "fuel_cell_types" => unique([string(sc.type_fuel_cell) for sc in result.config.simulation_configs]), # Log all unique fuel cell types
            "voltage_zones" => unique([string(sc.voltage_zone) for sc in result.config.simulation_configs]), # Log all unique voltage zones
            "date" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), # Timestamp the calibration completion
            "execution_time_seconds" => result.execution_time, # Record total execution duration in seconds
            "execution_time_formatted" => @sprintf("%.2f hours", result.execution_time / 3600) # Formatted execution time for readability
        ),
        "results" => Dict(
            "final_rmse_percent" => result.best_fitness, # Log the final best RMSE achieved
            "min_rmse_achieved" => result.min_rmse, # Log the minimum RMSE throughout history
            "generations_completed" => length(result.history) # Record total number of GA generations
        ),
        "config" => Dict(
            "pop_size" => result.config.ga_config.pop_size, # Log the population size hyperparameter
            "num_generations" => result.config.ga_config.num_generations, # Log the generation limit hyperparameter
            "num_parents_mating" => result.config.ga_config.num_parents_mating, # Log parents selected for mating
            "mutation_num_genes" => result.config.ga_config.mutation_num_genes, # Log number of genes to mutate
            "elitism" => result.config.ga_config.elitism, # Log the elitism count hyperparameter
            "target_error" => result.config.ga_config.target_error, # Log the target error hyperparameter
            "seed" => result.config.ga_config.seed # Record the random seed used for reproducibility
        ),
        "history" => result.history # Include the full RMSE history vector
    )

    try
        YAML.write_file(joinpath(output_dir, "calibration_report.yaml"), report) # Write the report to disk
    catch e # Handle reporting write errors
        @warn "Failed to save calibration report: $e" # Log failure warning
    end

    if !isempty(final_population) # Check if there is population data to save
        population_data = [ # Construct population summary list
            Dict("individual" => i, "params" => final_population[i], "rmse" => final_fitness[i])
            for i in 1:length(final_population) # Map each individual to its genes and fitness
        ]
        try
            YAML.write_file(joinpath(output_dir, "final_population.yaml"), population_data) # Write population data to disk
        catch e # Handle population write errors
            @warn "Failed to save final population: $e" # Log failure warning
        end
    end
end

"""
    _plot_calibration_results(result, output_dir)

Generate and save a figure showing the experimental vs. simulated polarization curves
for all calibration conditions using the best identified parameters.
"""
function _plot_calibration_results(result, output_dir::String)
    configs = result.config.simulation_configs
    n_configs = length(configs)
    n_configs == 0 && return nothing

    @info "Generating calibration results figure for $n_configs condition(s)..."

    # 1. Setup Figure
    # Adjust figure size based on the number of subplots
    fig = Figure(size = (800, 400 * n_configs))

    try
        for (i, sc) in enumerate(configs)
            # 2. Re-simulate each condition with the best identified parameters
            # Use the specified voltage zone for the final plot
            full_sc = SimulationConfig(
                type_fuel_cell = sc.type_fuel_cell,
                voltage_zone   = sc.voltage_zone,
                display_timing = :postrun # Ensure we get discretized points
            )

            fc = create_fuelcell(full_sc.type_fuel_cell, full_sc.voltage_zone)
            fc.physical_parameters = result.best_params # Use calibrated parameters

            # Use the complete polarization current profile
            cd = create_current(PolarizationParams(), fc)

            # Execute simulation
            sim = AlphaPEM(fc, cd, full_sc)
            simulate_model!(sim)

            # 3. Plot on a new axis
            ax = Axis(fig[i, 1])
            plot_polarization_curve(sim.outputs, fc, cd, full_sc, ax)

            # Add a subtitle or prefix to the axis title if multiple conditions exist
            if n_configs > 1
                ax.title = "Condition $i: " * ax.title[]
            end
        end

        # 4. Save the figure
        out_path = joinpath(output_dir, "calibration_polarization_curves.png")
        save(out_path, fig)
    catch e
        @warn "Failed to generate calibration plots: $e"
        # Don't rethrow, as this is a post-processing step that shouldn't crash the calibration
    end

    return nothing
end

end # module CalibrationHelpers
