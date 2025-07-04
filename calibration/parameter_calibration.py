# -*- coding: utf-8 -*-

"""
This file is designated for executing the undetermined parameters' calibration.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import sys
import os
import time
import pygad
import matplotlib.pyplot as plt
import numpy as np

# Importing constants' value and functions
sys.path.append(os.path.abspath('..'))  # Add parent (root) folder path to list of module search paths
from model.AlphaPEM import AlphaPEM
from modules.calibration_modules import determined_parameters, calculate_simulation_error, print_calibration_results, \
                                        save_calibration_results
from modules.main_modules import saving_instructions

# _________________________________________________Calibration settings_________________________________________________
"""
Users can select here the fuel cell to calibrate, along with the undetermined parameters to modify, and the parameters 
for the Genetic Algorithm. The experimental fuel cell data are stored in the files calibration_modules.py and 
experimental_values.py.
The parameters employed for the Genetic Algorithm here have proven to be effective, though not necessarily optimal.
"""

# Fuel cell possibilities: "EH-31_1.5"(2021), "EH-31_2.0"(2021), "EH-31_2.25"(2021), "EH-31_2.5"(2021),
#                          "BX_1.0"(2015), "BX_1.35"(2015), "LF"(2010), "manual_setup".
type_fuel_cell_1 = "EH-31_2.0"
type_fuel_cell_2 = "EH-31_2.25"

# Parameter bounds
#       Fuel cell physical parameters
epsilon_gdl_min, epsilon_gdl_max = 0.55, 0.80  # It is the anode/cathode GDL porosity, without units.
epsilon_mc_min, epsilon_mc_max = 0.15, 0.40  # It is the volume fraction of ionomer in the CL.
tau_min, tau_max = 1.0, 4.0  # It is the pore structure coefficient, without units.
epsilon_c_min, epsilon_c_max = 0.15, 0.30  # It is the compression ratio of the GDL.
#       Constants based on the interaction between water and the structure
e_min, e_max = 3, 5  # It is the capillary exponent, and should be an int number.
#       Voltage polarization
Re_min, Re_max = 5e-7, 5e-6  # ohm.m². It is the electron conduction resistance of the circuit.
i0_c_ref_min, i0_c_ref_max = 1e-3, 80  # A.m-2.It is the reference exchange current density at the cathode.
kappa_co_min, kappa_co_max = 0.01, 40  # A.m-2. It is the crossover correction coefficient.
kappa_c_min, kappa_c_max = 0.25, 4  # It is the overpotential correction exponent.
#       The bounds on liquid saturation coefficients are constrained to facilitate calibration.
a_slim_min, a_slim_max = 0.0, 0.2  # It is one of the limit liquid saturation coefficients.
b_slim_min, b_slim_max = 0.0, 0.4  # It is one of the limit liquid saturation coefficients.
a_switch_min, a_switch_max = 0.5, 1.0  # It is one of the limit liquid saturation coefficients.
#       Undetermined parameter which is not considered yet (require the use of EIS curves to be calibrated)
C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
#       Bounds gathering and type
varbound = [[epsilon_gdl_min, epsilon_gdl_max, 'real'],
            [epsilon_mc_min, epsilon_mc_max, 'real'],
            [tau_min, tau_max, 'real'],
            [epsilon_c_min, epsilon_c_max, 'real'],
            [e_min, e_max, 'int'],
            [Re_min, Re_max, 'real'],
            [i0_c_ref_min, i0_c_ref_max, 'real'],
            [kappa_co_min, kappa_co_max, 'real'],
            [kappa_c_min, kappa_c_max, 'real'],
            [a_slim_min, a_slim_max, 'real'],
            [b_slim_min, b_slim_max, 'real'],
            [a_switch_min, a_switch_max, 'real']]
gene_space = [] # List used to define the bounds of the undetermined parameters for pygad.
for i in range(len(varbound)):
    min_val, max_val, type_val = varbound[i]
    if type_val == 'int':
        gene_space.append({'low': min_val, 'high': max_val, 'step': 1})
    else:
        gene_space.append({'low': min_val, 'high': max_val})

# PyGAD parameters for the genetic algorithm:
    # Number of generations:
num_generations = 1 # It should be between 1000 and 1500, depending on the population_size,
    #                 for a good compromise between speed and precision.
    # Initial population (one solution means a member of the population):
        # 1) custom initial population.
initial_population = None # It is the initial population, which can be loaded from a file.
# initial_population = pygad.load(filename="results/EH-31/parameter_calibration_1").population
        # 2) random initial population.
sol_per_pop = 16 # It is the population size. It should be between 100 and 200 for a good compromise between speed and
    #              precision. Select a multiple of the available number of CPU cores for optimal multiprocessing.
num_genes = len(varbound) # Number of genes in the solution. It is the number of undetermined parameters.
    # Number of solutions to be selected as parents in the mating pool.
num_parents_mating = int(0.2 * sol_per_pop) # Here it is 20% of the population, the best found (has to be an integer).
    # Number of genes to mutate.
mutation_num_genes = 1 # It is the number of undetermined parameters to mutate, here taken to 1. It is the best found.

# ________________________________________________Determined parameters_________________________________________________
"""
The determined parameters of the experimental stack, along with the different operating conditions, are imported from 
the calibration_module.py.
This section should remain unaltered for regular program usage.
"""

T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, i_max_pola_1, Aact, Hmem, Hcl, Hgdl, Hgc, Wgc, Lgc, \
    type_auxiliary, type_control, type_purge, type_display, type_plot, type_current, current_density, t_step, i_step, \
    delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, t_purge, max_step, n_gdl, i_exp1, U_exp1 \
    = determined_parameters(type_fuel_cell_1)

T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, i_max_pola_2, Aact, Hmem, Hcl, Hgdl, Hgc, Wgc, Lgc, \
    type_auxiliary, type_control, type_purge, type_display, type_plot, type_current, current_density, t_step, i_step, \
    delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, t_purge, max_step, n_gdl, i_exp2, U_exp2 \
    = determined_parameters(type_fuel_cell_2)

# _______________________________________________________Functions______________________________________________________

def pola_points(ga_instance, solution, solution_idx): # Function to maximize.
    """
    This function calculates first the maximum error between the simulated polarization curves and the experimental
    ones. Then, it returns the inverse of this error as the fitness value to maximize by the genetic algorithm.

    Parameters
    ----------
    ga_instance : pygad.GA
        The instance of the genetic algorithm.
    solution : list
        The list of undetermined parameters to calibrate.
    solution_idx : int
        The index of the solution in the population.

    Returns
    -------
    fitness : float
        The fitness value to maximize, which is the inverse of the maximum error between the simulated and experimental
        polarization curves.
    """

    # Extraction of the undetermined parameters from the solution
    epsilon_gdl, epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch = solution

    # Calculation of the model polarization curve
    Simulator1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, t_step,
                          i_step, i_max_pola_1, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_scl, max_step, n_gdl, t_purge, type_fuel_cell_1, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)
    Simulator2 = AlphaPEM(current_density, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, t_step,
                          i_step, i_max_pola_2, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_scl, max_step, n_gdl, t_purge, type_fuel_cell_2, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)

    # Calculation of the simulation error between the simulated and experimental polarization curves
    sim_error = calculate_simulation_error(Simulator1, U_exp1, i_exp1, Simulator2, U_exp2, i_exp2)
    if sim_error == 0:  # If the error is zero, it means the simulation is perfect.
        sim_error = 1e-6

    # Calculation of the fitness value to maximize
    fitness = 1.0 / sim_error # pygad maximise the fitness value, not the error value, so the inverse is taken.
    return fitness

last_fitness = 0
def callback_generation(ga_instance):
    """
    Function to display the generation number and the fitness value, and to save the GA instance.

    Parameters
    ----------
    ga_instance : PyGAD object
        An instance of the PyGAD library, which is used to perform the optimization.

    Returns
    -------
    None
        This function does not return any value. It is used as a callback function to display information during the
    optimization process and to save the GA instance..
    """
    global last_fitness
    idx = np.argmax(ga_instance.last_generation_fitness)  # Get the index of the best solution in the last generation.
    fitness = ga_instance.last_generation_fitness[idx]  # Get the fitness value of the best solution.
    print(f"Generation = {ga_instance.generations_completed} / {num_generations+1} "
          f"= {ga_instance.generations_completed/(num_generations+1)*100:.2f}%")
    print(f"Fitness    = {fitness}")
    print(f"Change     = {fitness - last_fitness}")
    last_fitness = fitness
    ga_instance.save(filename="parameter_calibration_ongoing")  # Save the GA instance.

# ____________________________________________________Main program______________________________________________________
"""
This section is dedicated to ensuring the proper execution of the calibration. 
In this program, multi processing is permitted.
In certain instances, uncommon parameter combinations may result in errors if the chosen max_step value is excessively 
high, which sometimes crashes the program. If this occurs, it is restarted to avoid wasting cluster access time.
Warning: It may result in an infinite loop. A maximum execution time must be imposed at the cluster level.
This section should remain unaltered for regular program usage.
"""
if __name__ == '__main__':
    # Starting time
    start_time = time.time()

    # pYGAD execution
    ga_instance = pygad.GA( # 1 iteration takes 720s = 12min = 0.2h inside LIS cluster with 80 nodes.
    #                         800 iterations takes 160h = 6.67j.
        num_generations=num_generations, # Number of generations.
        num_parents_mating=num_parents_mating, # Number of solutions to be selected as parents in the mating pool.
        fitness_func=pola_points, # Function to maximize.
        initial_population=initial_population, # Initial population. It can be None for a random population.
        sol_per_pop=sol_per_pop, # Population size. It is the number of solutions in the population.
        num_genes=num_genes, # Number of genes in the solution.
        gene_space=gene_space, # Bounds of the undetermined parameters.
        parent_selection_type="rws", # Parent selection type. "rws" means roulette wheel selection. Best found.
        keep_parents=1, # Keep the best solution of the previous generation. 1 elite should be enough.
        crossover_type="single_point", # Crossover type. "single_point" means single point crossover. Best found.
        mutation_type="random", # Mutation type. "random" means random mutation.
        mutation_num_genes=mutation_num_genes, # Number of genes to mutate.
        stop_criteria=["reach_2"], # Stop when the error is less than 0.5%
        parallel_processing = ["process", os.cpu_count()],  # Use multiprocessing with the number of CPU cores available.
        on_generation = callback_generation # Callback function to display the generation number and the fitness value.
    )
    ga_instance.run() # Run the genetic algorithm.

    # Recovery of the results
    ga_instance.plot_fitness() # Plot the fitness value of the best solution at each generation.
    idx = np.argmax(ga_instance.last_generation_fitness) # Get the index of the best solution in the last generation.
    solution = ga_instance.population[idx] # Get the best solution from the last generation.
    solution_fitness = ga_instance.last_generation_fitness[idx] # Get the fitness value of the best solution.
    epsilon_gdl, epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch = solution
    sim_error = 1.0 / solution_fitness  # The error is the inverse of the fitness value.

    # Print of the parameter calibration results
    convergence = [float(1.0 / f) for f in ga_instance.best_solutions_fitness]
    print_calibration_results(convergence, ga_instance, epsilon_gdl, epsilon_mc, tau, epsilon_c, e,
                              Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error)

    # Save the data in files
    save_calibration_results(convergence, ga_instance, epsilon_gdl, epsilon_mc, tau, epsilon_c, e,
                             Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error, type_fuel_cell_1)

    # Calculate, display and save the calibrated and experimental polarization curve
    #       Calculate the calibrated polarization curve
    Simulator1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, t_step,
                          i_step, i_max_pola_1, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_scl, max_step, n_gdl, t_purge, type_fuel_cell_1, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)
    Simulator2 = AlphaPEM(current_density, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, t_step,
                          i_step, i_max_pola_2, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_scl, max_step, n_gdl, t_purge, type_fuel_cell_2, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)
    #       Display the calibrated and experimental polarization curve
    fig, ax = plt.subplots(figsize=(8, 8))
    Simulator1.Display(ax)
    Simulator2.Display(ax)
    #       Save the calibrated and experimental polarization curve
    subfolder_name = type_fuel_cell_1[:type_fuel_cell_1.rfind('_')] if type_fuel_cell_1.rfind('_') != -1 \
        else type_fuel_cell_1
    saving_instructions("results", subfolder_name, "final_calibration_1.pdf", fig)

    # Algorithm time
    algo_time = time.time() - start_time
    print('Time of the algorithm in hours :', algo_time/3600)