# -*- coding: utf-8 -*-

"""
This file is designated for executing the undetermined parameters' calibration.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using PyCall
using Dates
using Statistics

# Keep matplotlib and pygad via Python interop
const plt = pyimport("matplotlib.pyplot")
const pygad = pyimport("pygad")
const np = pyimport("numpy")

# Importing constants' value and functions
include(joinpath(@__DIR__, "../core/models/AlphaPEM.jl"))
include(joinpath(@__DIR__, "../config/current_densities.jl"))
include(joinpath(@__DIR__, "calibration_modules.jl"))

# ________________________________________________Determined parameters_________________________________________________
"""
Users can select here the fuel cell to calibrate. The determined parameters of the experimental stack, along with the
different operating conditions, are imported from the calibration_modules.jl.
"""

# Fuel cell possibilities: "ZSW-GenStack"(2022), "ZSW-GenStack_Pa_1.61_Pc_1.41", "ZSW-GenStack_Pa_2.01_Pc_1.81",
#                          "ZSW-GenStack_Pa_2.4_Pc_2.2", "ZSW-GenStack_Pa_2.8_Pc_2.6", "ZSW-GenStack_T_62",
#                          "ZSW-GenStack_T_76", "ZSW-GenStack_T_84", "EH-31_1.5"(2021), "EH-31_2.0"(2021),
#                          "EH-31_2.25"(2021), "EH-31_2.5"(2021)
type_fuel_cell_1::String = "ZSW-GenStack_Pa_2.01_Pc_1.81"
type_fuel_cell_2::String = "ZSW-GenStack_Pa_2.8_Pc_2.6"

# Calibration zone: "before_voltage_drop", "full".
voltage_zone::String = "full"

(operating_inputs_1, current_parameters, accessible_physical_parameters, undetermined_physical_parameters,
 model_parameters, computing_parameters_1, i_exp_1, U_exp_1) =
    parameters_for_calibration(type_fuel_cell_1, voltage_zone)

(operating_inputs_2, current_parameters, accessible_physical_parameters, undetermined_physical_parameters,
 model_parameters, computing_parameters_2, i_exp_2, U_exp_2) =
    parameters_for_calibration(type_fuel_cell_2, voltage_zone)

# _________________________________________________Calibration settings_________________________________________________
"""
Users can select here the parameters for the Genetic Algorithm. The experimental fuel cell data are stored in the files
calibration_modules.jl and experimental_values.jl.
The parameters employed for the Genetic Algorithm here have proven to be effective, though not necessarily optimal.
"""

# Parameter bounds
varbound, gene_space = parameter_bounds_for_calibration(type_fuel_cell_1, voltage_zone, operating_inputs_1, operating_inputs_2)

# PyGAD parameters for the genetic algorithm:
# Number of generations:
num_generations::Int = 500  # 1000 generations are good, 2000 generations are better.
#                             10 generation of 128 elements takes approximatively ??min on my computer (16 CPU cores).
#                             100 generations of 128 elements take approximatively 10h on UR cluster (1 node of 32 CPU cores).
# Initial population (one solution means a member of the population):
# 1) random initial population.
initial_population = nothing  # It is the initial population, which can be loaded from a file.
# 2) custom initial population.
# initial_population = pygad.load(filename="results/ZSW-GenStack/parameter_calibration_1").population # Remark: it seems that the
#                                 file should be loaded from the same computer and in the same folder it is created ...?
sol_per_pop::Int = 128  # It is the population size. It should be between 100 and 200 for a good compromise between speed and
#                         precision. Select a multiple of the available number of CPU cores for optimal multiprocessing.
num_genes::Int = length(varbound)  # Number of genes in the solution. It is the number of undetermined parameters.
# Number of solutions to be selected as parents in the mating pool.
num_parents_mating::Int = Int(floor(0.2 * sol_per_pop))  # Here it is 20% of the population, the best found (has to be an integer).
# Number of genes to mutate.
mutation_num_genes::Int = 1  # It is the number of undetermined parameters to mutate, here taken to 1. It is the best found.

# _______________________________________________________Functions______________________________________________________

"""Calculates the fitness value for the genetic algorithm.

This function calculates first the maximum error between the simulated polarization curves and the experimental
ones. Then, it returns the inverse of this error as the fitness value to maximize by the genetic algorithm.

# Arguments
- `ga_instance`: The instance of the genetic algorithm.
- `solution::Vector{Float64}`: The list of undetermined parameters to calibrate.
- `solution_idx::Int`: The index of the solution in the population.

# Returns
- `fitness::Float64`: The fitness value to maximize, which is the inverse of the maximum error between the simulated
  and experimental polarization curves.
"""
function pola_points(ga_instance, solution::Vector{Float64}, solution_idx::Int)::Float64

    # Extraction of the undetermined parameters from the solution
    solution_of_undetermined_physical_parameters = update_undetermined_parameters(
        type_fuel_cell_1, solution, varbound, deepcopy(undetermined_physical_parameters)
    )

    try
        # Calculation of the model polarization curve
        # Create the simulators
        simulator_1 = AlphaPEM(accessible_physical_parameters, solution_of_undetermined_physical_parameters, model_parameters)
        simulator_2 = AlphaPEM(accessible_physical_parameters, solution_of_undetermined_physical_parameters, model_parameters)

        # Simulate the model with the current parameters and the operating conditions of the experimental points
        simulate_model!(simulator_1, operating_inputs_1, current_parameters, computing_parameters_1)
        simulate_model!(simulator_2, operating_inputs_2, current_parameters, computing_parameters_2)

        # Calculation of the simulation error between the simulated and experimental polarization curves
        sim_error = calculate_simulation_error(simulator_1, U_exp_1, i_exp_1, simulator_2, U_exp_2, i_exp_2)

        if sim_error == 0  # If the error is zero, it means the simulation is perfect.
            sim_error = 1e-6
        end

        # Calculation of the fitness value to maximize
        fitness = 1.0 / sim_error  # pygad maximise the fitness value, not the error value, so the inverse is taken.
        return fitness

    catch e  # If an error occurs during the simulation, return a very low fitness value to refuse this solution
        println("\nAn error occurred during the evaluation of the solution.")
        params = [string(k, ": ", v) for (k, v) in pairs(solution_of_undetermined_physical_parameters)]
        println("Attempted parameters: " * join(params, " | "))
        println("Exception : ", e)
        println("Refusing this solution and continuing the optimization.\n")
        return 1e-12  # Very low fitness value to refuse this solution
    end
end

# Global variable to track fitness
last_fitness::Float64 = 0.0

"""Function to display the generation number and the fitness value, and to save the GA instance.

# Arguments
- `ga_instance`: An instance of the PyGAD library, which is used to perform the optimization.

# Returns
- `Nothing`: This function does not return any value. It is used as a callback function to display information during
  the optimization process and to save the GA instance.
"""
function callback_generation(ga_instance)::Nothing
    global last_fitness
    idx = py"int"(np.argmax(ga_instance.last_generation_fitness))  # Get the index of the best solution in the last generation.
    fitness = py"float"(ga_instance.last_generation_fitness[idx])  # Get the fitness value of the best solution.
    println("Generation = $(ga_instance.generations_completed) / $num_generations = $(ga_instance.generations_completed / num_generations * 100)%")
    println("Fitness    = $fitness")
    println("Change     = $(fitness - last_fitness)")
    global last_fitness = fitness
    ga_instance.save(filename="parameter_calibration_ongoing")  # Save the GA instance.
    return nothing
end

# ____________________________________________________Main program______________________________________________________
"""
This section is dedicated to ensuring the proper execution of the calibration.
In this program, multi processing is permitted.
This section should remain unaltered for regular program usage.
"""

function run_calibration()::Nothing

    # Starting time
    start_time = time()

    # pYGAD execution
    ga_instance = pygad.GA(
        num_generations=num_generations,  # Number of generations.
        num_parents_mating=num_parents_mating,  # Number of solutions to be selected as parents in the mating pool.
        fitness_func=pola_points,  # Function to maximize.
        initial_population=initial_population,  # Initial population. It can be None for a random population.
        sol_per_pop=sol_per_pop,  # Population size. It is the number of solutions in the population.
        num_genes=num_genes,  # Number of genes in the solution.
        gene_space=gene_space,  # Bounds of the undetermined parameters.
        parent_selection_type="rws",  # Parent selection type. "rws" means roulette wheel selection. Best found.
        keep_parents=1,  # Keep the best solution of the previous generation. 1 elite should be enough.
        crossover_type="single_point",  # Crossover type. "single_point" means single point crossover. Best found.
        mutation_type="random",  # Mutation type. "random" means random mutation.
        mutation_num_genes=mutation_num_genes,  # Number of genes to mutate.
        stop_criteria=["reach_2"],  # Stop when the error is less than 0.5%
        parallel_processing=["process", Sys.cpu_info() |> length],  # Use multiprocessing with available CPU cores.
        on_generation=callback_generation,  # Callback function to display generation number and fitness value.
    )
    ga_instance.run()  # Run the genetic algorithm.

    # Recovery of the results
    # ga_instance.plot_fitness() # Plot the fitness value of the best solution at each generation.
    idx = py"int"(np.argmax(ga_instance.last_generation_fitness))  # Get the index of the best solution in the last generation.
    solution = ga_instance.population[idx]  # Get the best solution from the last generation.
    solution_fitness = py"float"(ga_instance.last_generation_fitness[idx])  # Get the fitness value of the best solution.
    sim_error = 1.0 / solution_fitness  # The error is the inverse of the fitness value.
    undetermined_physical_parameters_calibrated = update_undetermined_parameters(
        type_fuel_cell_1, solution, varbound, deepcopy(undetermined_physical_parameters)
    )

    # Print of the parameter calibration results
    convergence = [float(1.0 / f) for f in ga_instance.best_solutions_fitness]
    print_calibration_results(convergence, ga_instance, solution, varbound, sim_error)

    # Save the data in files
    save_calibration_results(convergence, ga_instance, solution, varbound, sim_error, type_fuel_cell_1 * " and " * type_fuel_cell_2)

    # Calculate, display and save the calibrated and experimental polarization curve (from the complete experimental points)
    #       Calculate the calibrated polarization curve
    operating_inputs_final_1 = deepcopy(operating_inputs_1)
    computing_parameters_final_1 = deepcopy(computing_parameters_1)
    operating_inputs_final_1["current_density"] = polarization_current
    computing_parameters_final_1["type_current"] = "polarization"

    operating_inputs_final_2 = deepcopy(operating_inputs_2)
    computing_parameters_final_2 = deepcopy(computing_parameters_2)
    operating_inputs_final_2["current_density"] = polarization_current
    computing_parameters_final_2["type_current"] = "polarization"

    simulator_final_1 = AlphaPEM(accessible_physical_parameters, undetermined_physical_parameters_calibrated, model_parameters)
    simulator_final_2 = AlphaPEM(accessible_physical_parameters, undetermined_physical_parameters_calibrated, model_parameters)
    simulate_model!(simulator_final_1, operating_inputs_final_1, current_parameters, computing_parameters_final_1)
    simulate_model!(simulator_final_2, operating_inputs_final_2, current_parameters, computing_parameters_final_2)

    #       Display the calibrated and experimental polarization curve
    fig, ax = plt.subplots(figsize=(8, 8))
    Display(simulator_final_1, ax)
    Display(simulator_final_2, ax)

    #       Save the calibrated and experimental polarization curve
    last_underscore_idx = findlast('_', type_fuel_cell_1)
    subfolder_name = (last_underscore_idx !== nothing) ? type_fuel_cell_1[1:(last_underscore_idx - 1)] : type_fuel_cell_1
    Saving_instructions(simulator_final_1, "results", subfolder_name, "final_calibration_1.pdf", fig)

    # Algorithm time
    algo_time = time() - start_time
    println("Time of the algorithm in hours : ", algo_time / 3600)

    return nothing
end

# Run the calibration if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_calibration()
end

