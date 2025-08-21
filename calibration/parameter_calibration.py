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
from configuration.current_densities import step_current, polarization_current
from modules.calibration_modules import parameters_for_calibration, calculate_simulation_error, print_calibration_results, \
                                        save_calibration_results

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
Hcl_min, Hcl_max = 8e-6, 2e-5  # m. It is the thickness of the CL.
Hmem_min, Hmem_max = 1.5e-5, 5e-5  # m. It is the thickness of the membrane.
epsilon_gdl_min, epsilon_gdl_max = 0.50, 0.90  # It is the anode/cathode GDL porosity, without units.
epsilon_mc_min, epsilon_mc_max = 0.15, 0.40  # It is the volume fraction of ionomer in the CL.
# epsilon_c_min, epsilon_c_max = 0.15, 0.30  # It is the compression ratio of the GDL.
epsilon_c = 0.20  # It is the compression ratio of the GDL.
#       Constants based on the interaction between water and the structure
e_min, e_max = 3, 5  # It is the capillary exponent, and should be an int number.
#       Voltage polarization
i0_c_ref_min, i0_c_ref_max = 1e-3, 80  # A.m-2.It is the reference exchange current density at the cathode.
kappa_co_min, kappa_co_max = 0.01, 40  # A.m-2. It is the crossover correction coefficient.
# kappa_co = 1  # A.m-2. It is the crossover correction coefficient.
kappa_c_min, kappa_c_max = 0.25, 4  # It is the overpotential correction exponent.
#       The bounds on liquid saturation coefficients are constrained to facilitate calibration.
# a_slim_min, a_slim_max = 0.0, 0.2  # It is one of the limit liquid saturation coefficients.
# b_slim_min, b_slim_max = 0.05, 0.4  # It is one of the limit liquid saturation coefficients.
# a_switch_min, a_switch_max = 0.5, 0.95  # It is one of the limit liquid saturation coefficients.
a_slim = 0.05553  # It is one of the limit liquid saturation coefficients.
b_slim = 0.10514  # It is one of the limit liquid saturation coefficients.
a_switch = 0.63654  # It is one of the limit liquid saturation coefficients.
#       Undetermined parameter which is not considered yet (require the use of EIS curves to be calibrated)
C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
#       Bounds gathering and type
varbound = [[Hcl_min, Hcl_max, 'real'],
            [Hmem_min, Hmem_max, 'real'],
            [epsilon_gdl_min, epsilon_gdl_max, 'real'],
            [epsilon_mc_min, epsilon_mc_max, 'real'],
            [e_min, e_max, 'int'],
            [i0_c_ref_min, i0_c_ref_max, 'real'],
            [kappa_co_min, kappa_co_max, 'real'],
            [kappa_c_min, kappa_c_max, 'real']]
gene_space = [] # List used to define the bounds of the undetermined parameters for pygad.
for i in range(len(varbound)):
    min_val, max_val, type_val = varbound[i]
    if type_val == 'int':
        gene_space.append({'low': min_val, 'high': max_val, 'step': 1})
    else:
        gene_space.append({'low': min_val, 'high': max_val})

# PyGAD parameters for the genetic algorithm:
    # Number of generations:
num_generations = 1000 # 1000 generations should be enough.
    #                   10 generation of 128 elements takes approximatively 100min on my computer (16 CPU cores).
    #                   1000 generations of 128 elements take less than 80h on UR cluster (1 node of 32 CPU cores).
    # Initial population (one solution means a member of the population):
        # 1) random initial population.
initial_population = None # It is the initial population, which can be loaded from a file.
        # 2) custom initial population.
# initial_population = pygad.load(filename="results/EH-31/parameter_calibration_1").population # Remark: it seems that the
#                                 file should be loaded from the same computer and in the same folder it is created ...?
sol_per_pop = 128 # It is the population size. It should be between 100 and 200 for a good compromise between speed and
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

(T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, step_current_parameters, pola_current_parameters,
 pola_current_for_cali_parameters_1, Aact, Hgdl, Hgc, Wgc, Lgc, estimated_undetermined_parameters_for_initialisation,
 type_auxiliary, type_control, type_purge, type_display, type_plot, type_current, current_density, i_EIS, ratio_EIS,
 t_EIS, f_EIS, t_purge, n_gdl, i_exp_1, U_exp_1) \
    = parameters_for_calibration(type_fuel_cell_1)

(T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, step_current_parameters, pola_current_parameters,
 pola_current_for_cali_parameters_2, Aact, Hgdl, Hgc, Wgc, Lgc, estimated_undetermined_parameters_for_initialisation,
 type_auxiliary, type_control, type_purge, type_display, type_plot, type_current, current_density, i_EIS, ratio_EIS,
 t_EIS, f_EIS, t_purge, n_gdl,  i_exp_2, U_exp_2) \
    = parameters_for_calibration(type_fuel_cell_2)

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
    # Hcl, Hmem, epsilon_gdl, epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch = solution
    Hcl, Hmem, epsilon_gdl, epsilon_mc, e, i0_c_ref, kappa_co, kappa_c = solution

    # Calculation of the model polarization curve
    Simulator_1 = AlphaPEM(current_density, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                           step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters_1, i_EIS,
                           ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, epsilon_mc,
                           epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl, n_gdl, t_purge,
                           type_fuel_cell_1, type_current, type_auxiliary, type_control, type_purge, type_display,
                           type_plot, initial_variable_values_1)
    Simulator_2 = AlphaPEM(current_density, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2,
                           step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters_2, i_EIS,
                           ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl, epsilon_mc,
                           epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl, n_gdl, t_purge,
                           type_fuel_cell_2, type_current, type_auxiliary, type_control, type_purge, type_display,
                           type_plot, initial_variable_values_2)

    # Calculation of the simulation error between the simulated and experimental polarization curves
    sim_error = calculate_simulation_error(Simulator_1, U_exp_1, i_exp_1, Simulator_2, U_exp_2, i_exp_2)
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
    print(f"Generation = {ga_instance.generations_completed} / {num_generations} "
          f"= {ga_instance.generations_completed/(num_generations)*100:.2f}%")
    print(f"Fitness    = {fitness}")
    print(f"Change     = {fitness - last_fitness}")
    last_fitness = fitness
    ga_instance.save(filename="parameter_calibration_ongoing")  # Save the GA instance.

# ____________________________________________________Main program______________________________________________________
"""
This section is dedicated to ensuring the proper execution of the calibration. 
In this program, multi processing is permitted.
This section should remain unaltered for regular program usage.
"""
if __name__ == '__main__':
    # Starting time
    start_time = time.time()

    # Initialization of the simulator in order to equilibrate the internal states of the fuel cell
    Hcl_ini = estimated_undetermined_parameters_for_initialisation['Hcl']
    Hmem_ini = estimated_undetermined_parameters_for_initialisation['Hmem']
    epsilon_gdl_ini = estimated_undetermined_parameters_for_initialisation['epsilon_gdl']
    epsilon_mc_ini = estimated_undetermined_parameters_for_initialisation['epsilon_mc']
    # epsilon_c_ini = estimated_undetermined_parameters_for_initialisation['epsilon_c']
    e_ini = estimated_undetermined_parameters_for_initialisation['e']
    i0_c_ref_ini = estimated_undetermined_parameters_for_initialisation['i0_c_ref']
    kappa_co_ini = estimated_undetermined_parameters_for_initialisation['kappa_co']
    kappa_c_ini = estimated_undetermined_parameters_for_initialisation['kappa_c']
    a_slim_ini = estimated_undetermined_parameters_for_initialisation['a_slim']
    b_slim_ini = estimated_undetermined_parameters_for_initialisation['b_slim']
    a_switch_ini = estimated_undetermined_parameters_for_initialisation['a_switch']
    C_scl_ini = estimated_undetermined_parameters_for_initialisation['C_scl']
    Simulator_ini_1 = AlphaPEM(step_current, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                               step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters_1,
                               i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem_ini, Hcl_ini, Hgc, Wgc, Lgc,
                               epsilon_gdl_ini, epsilon_mc_ini, epsilon_c, e_ini, i0_c_ref_ini, kappa_co_ini,
                               kappa_c_ini, a_slim_ini, b_slim_ini, a_switch_ini, C_scl_ini, n_gdl, t_purge,
                               type_fuel_cell_1, 'step', type_auxiliary, type_control, type_purge,
                               type_display, type_plot)
    Simulator_ini_2 = AlphaPEM(step_current, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2,
                               step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters_2,
                               i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem_ini, Hcl_ini, Hgc, Wgc, Lgc,
                               epsilon_gdl_ini, epsilon_mc_ini, epsilon_c, e_ini, i0_c_ref_ini, kappa_co_ini,
                               kappa_c_ini, a_slim_ini, b_slim_ini, a_switch_ini, C_scl_ini, n_gdl, t_purge,
                               type_fuel_cell_2, 'step', type_auxiliary, type_control, type_purge,
                               type_display, type_plot)
    #   Recovery of the internal states from the end of the preceding simulation.
    initial_variable_values_1, initial_variable_values_2 = [], []
    for x in Simulator_ini_1.solver_variable_names:
        initial_variable_values_1.append(Simulator_ini_1.variables[x][-1])
        initial_variable_values_2.append(Simulator_ini_2.variables[x][-1])

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
    # ga_instance.plot_fitness() # Plot the fitness value of the best solution at each generation.
    idx = np.argmax(ga_instance.last_generation_fitness) # Get the index of the best solution in the last generation.
    solution = ga_instance.population[idx] # Get the best solution from the last generation.
    solution_fitness = ga_instance.last_generation_fitness[idx] # Get the fitness value of the best solution.
    # epsilon_gdl, epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch = solution
    Hcl, Hmem, epsilon_gdl, epsilon_mc, e, i0_c_ref, kappa_co, kappa_c = solution
    sim_error = 1.0 / solution_fitness  # The error is the inverse of the fitness value.

    # Print of the parameter calibration results
    convergence = [float(1.0 / f) for f in ga_instance.best_solutions_fitness]
    print_calibration_results(convergence, ga_instance, Hcl, Hmem, epsilon_gdl, epsilon_mc, epsilon_c, e,
                              i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error)

    # Save the data in files
    save_calibration_results(convergence, ga_instance, Hcl, Hmem, epsilon_gdl, epsilon_mc, epsilon_c, e,
                             i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error, type_fuel_cell_1)

    # Calculate, display and save the calibrated and experimental polarization curve (from the complete experimental points)
    #       Calculate the calibrated polarization curve
    Simulator_final_1 = AlphaPEM(polarization_current, T_des_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1,
                                 step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters_1,
                                 i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                 epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl,
                                 n_gdl, t_purge, type_fuel_cell_1, 'polarization', type_auxiliary,
                                 type_control, type_purge, type_display, type_plot, initial_variable_values_1)
    Simulator_final_2 = AlphaPEM(polarization_current, T_des_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2,
                                 step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters_2,
                                 i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc, Lgc, epsilon_gdl,
                                 epsilon_mc, epsilon_c, e, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl,
                                 n_gdl, t_purge, type_fuel_cell_2, 'polarization', type_auxiliary,
                                 type_control, type_purge, type_display, type_plot, initial_variable_values_2)
    #       Display the calibrated and experimental polarization curve
    fig, ax = plt.subplots(figsize=(8, 8))
    Simulator_final_1.Display(ax)
    Simulator_final_2.Display(ax)
    #       Save the calibrated and experimental polarization curve
    subfolder_name = type_fuel_cell_1[:type_fuel_cell_1.rfind('_')] if type_fuel_cell_1.rfind('_') != -1 \
        else type_fuel_cell_1
    Simulator_final_1.Saving_instructions("results", subfolder_name, "final_calibration_1.pdf", fig)

    # Algorithm time
    algo_time = time.time() - start_time
    print('Time of the algorithm in hours :', algo_time/3600)