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

import AlphaPEM.Core.Models: AlphaPEM, simulate_model!, _polarization_points_cali, _calculate_rmse
import AlphaPEM.Config: PolaExperimentalData
using ...Calibration: CalibrationConfig
using ...ValidParameterRegion.ResultsExport: export_parameter_bounds
using ...ValidParameterRegion.ConfigurationSampling: new_PhysicalParams_from_sample

export _fitness_function,
       calculate_simulation_error,
       _save_intermediate,
       _load_warm_start_population,
       _save_final_results

"""
    _fitness_function(gene_values, parameter_bounds, base_params, fuel_cells, current_profiles, simulation_configs) -> Float64

Evaluate the Root Mean Square Error (RMSE) for a given parameter gene vector across all operating conditions.
"""
function _fitness_function(gene_values::Vector{Float64},
                           parameter_bounds,
                           base_params,
                           fuel_cells,
                           current_profiles,
                           simulation_configs)::Float64
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

        return mean(rmse_values) # Return the average RMSE across all conditions

    catch e # Catch any simulation or numerical errors
        println("\nAn error occurred during the evaluation of the solution.")
        params = [string(b.name, ": ", gene_values[i]) for (i, b) in enumerate(parameter_bounds.bounds)]
        println("Attempted parameters: " * join(params, " | "))
        println("Exception : ", e)
        println("Refusing this solution and continuing the optimization.\n")
        return 1e6 # Return a large penalty value on failure
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
    _save_intermediate(best_gene_values, best_fitness, generation, parameter_bounds, base_params, output_dir)

Save the current best individual's gene values as a checkpoint for potential warm-start recovery.
"""
function _save_intermediate(best_gene_values, best_fitness, generation, parameter_bounds, base_params, output_dir)
    intermediate_bounds = Dict{Symbol, Tuple{Float64, Float64}}( # Prepare temporary bounds dictionary for export
        bound.name => (best_gene_values[i], best_gene_values[i]) for (i, bound) in enumerate(parameter_bounds.bounds)
    )

    try
        export_parameter_bounds(intermediate_bounds, joinpath(output_dir, "calibration_checkpoint.yaml"); # Export checkpoint to YAML
                                method = :checkpoint, # Mark as an intermediate checkpoint
                                metadata = Dict("generation" => generation, # Log the current generation index
                                                "rmse_percent" => best_fitness)) # Log the current best accuracy
    catch e # Handle filesystem or export errors
        @warn "Failed to save checkpoint: $e" # Log warning message
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
                    individual[idx] = param_data["min"] # Set gene value from minimum field
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
            "mating_rate" => result.config.ga_config.mating_rate, # Log the mating rate hyperparameter
            "mutation_genes" => result.config.ga_config.mutation_genes, # Log the number of mutated genes hyperparameter
            "elitism" => result.config.ga_config.elitism, # Log the elitism count hyperparameter
            "seed" => result.config.ga_config.seed # Record the random seed used for reproducibility
        ),
        "history" => result.history # Include the full RMSE history vector
    )

    try
        YAML.write_file(joinpath(output_dir, "calibration_report.yaml"), report) # Write the report to disk
        @info "Calibration report saved to $(joinpath(output_dir, "calibration_report.yaml"))" # Log successful save
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
            @info "Final population saved to $(joinpath(output_dir, "final_population.yaml"))" # Log successful save
        catch e # Handle population write errors
            @warn "Failed to save final population: $e" # Log failure warning
        end
    end
end

end # module CalibrationHelpers
