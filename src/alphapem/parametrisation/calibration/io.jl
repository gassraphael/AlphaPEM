"""
    _params_dict_for_export(best_params, parameter_bounds) -> Dict{Symbol, Real}

Build a dictionary of calibrated parameters for export. Angular parameters such as `theta_c_gdl` are converted from radians to degrees for readability.
"""
function _params_dict_for_export(best_params, parameter_bounds)
    return Dict{Symbol, Real}(
        bound.name => if bound.name == :theta_c_gdl
            rad2deg(getfield(best_params, bound.name))
        else
            getfield(best_params, bound.name)
        end
        for bound in parameter_bounds.bounds
    )
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
    _generate_calibration_output_dir(cfg::CalibrationConfig) -> String

Generate a timestamped output directory for a calibration run, mirroring the convention
used in `run_parameter_validity.jl`.

The directory name includes the run date, fuel-cell type, voltage zone, number of
experimental conditions, and an auto-incremented version suffix to avoid overwriting
previous runs.
"""
function _generate_calibration_output_dir(cfg::CalibrationConfig)::String
    ref_cfg = first(cfg.simulation_configs)
    date_str = Dates.format(Dates.today(), "yyyy.mm.dd")
    fc_str   = string(ref_cfg.type_fuel_cell)
    zone_str = string(ref_cfg.voltage_zone)
    n_cond   = length(cfg.simulation_configs)

    base_dir = cfg.output_dir
    mkpath(base_dir)

    base_dirname = "$date_str - $fc_str - $zone_str - $n_cond condition$(n_cond > 1 ? "s" : "")"

    version = 1
    while isdir(joinpath(base_dir, "$base_dirname - V$version"))
        version += 1
    end

    run_dir = joinpath(base_dir, "$base_dirname - V$version")
    mkpath(run_dir)
    return run_dir
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
