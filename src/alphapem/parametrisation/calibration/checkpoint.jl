"""
    _on_generation(ga_instance, history, ga_config, cfg, last_save_time, parameter_bounds, base_params)

Callback function executed after each generation of the Genetic Algorithm.
Handles progress logging and periodic checkpoint saving.
"""
function _on_generation(ga_instance, history, ga_config, cfg, last_save_time, parameter_bounds, base_params)
    # PyGAD maximizes fitness, so fitness = 1/RMSE (positive, compatible with RWS)
    best_sol, best_fitness_py, _ = ga_instance.best_solution()
    current_best_rmse = 1.0 / pyconvert(Float64, best_fitness_py) # Convert 1/RMSE back to RMSE
    push!(history, current_best_rmse) # Log current RMSE to history

    generation = pyconvert(Int, ga_instance.generations_completed) # Get current generation index
    msg = @sprintf("Generation %d/%d: Best RMSE = %.2f %%", generation, ga_config.num_generations, current_best_rmse)
    print("\r[ Info: ", msg) # Log progress in-place
    flush(stdout)

    if (cfg.save_frequency > 0 && generation % cfg.save_frequency == 0) || (time() - last_save_time[] > 300) # Check for save triggers
        best_parameters = pyconvert(Vector{Float64}, best_sol)
        _save_intermediate(ga_instance, best_parameters, current_best_rmse, generation, parameter_bounds, base_params, cfg.output_dir) # Save checkpoint
        last_save_time[] = time() # Reset last save timer
    end
end

"""
    _save_intermediate(ga_instance, best_gene_values, best_fitness, generation, parameter_bounds, base_params, output_dir)

Save the current best individual's calibrated parameter values and the full population as checkpoints for potential warm-start recovery.
"""
function _save_intermediate(ga_instance, best_gene_values, best_fitness, generation, parameter_bounds, base_params, output_dir)
    best_params = new_PhysicalParams_from_sample(Float64.(best_gene_values), parameter_bounds, base_params)
    best_params_dict = _params_dict_for_export(best_params, parameter_bounds)

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
        population = pyconvert(Matrix{Float64}, ga_instance.population) # Extract current population matrix from PyGAD
        fitness_values = pyconvert(Vector{Float64}, ga_instance.last_generation_fitness) # Extract last generation fitness values
        population_data = [
            Dict("individual" => i, "params" => collect(population[i, :]), "rmse" => 1.0 / fitness_values[i])
            for i in 1:size(population, 1)
        ]
        YAML.write_file(joinpath(output_dir, "calibration_checkpoint_population.yaml"), population_data) # Write full population
    catch e
        @warn "Failed to save checkpoint population: $e"
    end
end
