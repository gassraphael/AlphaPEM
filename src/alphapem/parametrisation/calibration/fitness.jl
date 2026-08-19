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
    gene_values = solution isa Py ? pyconvert(Vector{Float64}, solution) : Float64.(solution) # Convert PyObject to Julia Float64 vector
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
    jl_population = pyconvert(Matrix{Float64}, population)  # copy Python array to Julia before threading
    n = size(jl_population, 1)
    fitness_values = Vector{Float64}(undef, n)
    @sync for i in 1:n
        Threads.@spawn begin
            fitness_values[i] = _fitness_function(
                jl_population[i, :], parameter_bounds, base_params, fuel_cells, current_profiles, simulation_configs
            )
        end
    end
    return fitness_values
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
    if isempty(Ucell_sim) # If simulation failed to produce valid points, return a high error penalty
        return 100.0
    end
    return _calculate_rmse(Ucell_sim, experimental_data.U_exp)
end
