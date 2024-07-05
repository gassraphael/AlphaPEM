# -*- coding: utf-8 -*-

"""
This file is designated for executing the undetermined parameters' calibration.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import sys
import os
import time
import matplotlib.pyplot as plt
import numpy as np
from geneticalgorithm2 import GeneticAlgorithm2 as ga
from geneticalgorithm2 import AlgorithmParams

# Importing constants' value and functions
sys.path.append(os.path.abspath('..'))  # Add parent (root) folder path to list of module search paths
from model.AlphaPEM import AlphaPEM
from modules.calibration_modules import determined_parameters, calculate_simulation_error, \
    print_calibration_results, save_calibration_results
from modules.main_modules import saving_instructions

# _________________________________________________Calibration settings_________________________________________________
"""
Users can select the experimental fuel cell data for calibration from calibration_modules.py and experimental_values.py,
along with the undetermined parameters to modify, and the parameters for the Genetic Algorithm. 
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
i0_c_ref_min, i0_c_ref_max = 1e-3, 5e2  # A.m-2.It is the reference exchange current density at the cathode.
kappa_co_min, kappa_co_max = 0.01, 40  # A.m-2. It is the crossover correction coefficient.
kappa_c_min, kappa_c_max = 0, 100  # It is the overpotential correction exponent.
#       The bounds on liquid saturation coefficients are constrained to facilitate calibration.
a_slim_min, a_slim_max = 0.0, 0.2  # It is one of the limit liquid saturation coefficients.
b_slim_min, b_slim_max = 0.0, 0.4  # It is one of the limit liquid saturation coefficients.
a_switch_min, a_switch_max = 0.5, 1.0  # It is one of the limit liquid saturation coefficients.
#       Bounds gathering
varbound = np.array([[epsilon_gdl_min, epsilon_gdl_max],
                     [epsilon_mc_min, epsilon_mc_max], [tau_min, tau_max],
                     [epsilon_c_min, epsilon_c_max], [e_min, e_max], [Re_min, Re_max],
                     [i0_c_ref_min, i0_c_ref_max], [kappa_co_min, kappa_co_max],
                     [kappa_c_min, kappa_c_max], [a_slim_min, a_slim_max],
                     [b_slim_min, b_slim_max], [a_switch_min, a_switch_max]])
#       Bounds type
vartype = np.array(['real', 'real', 'real', 'real', 'int', 'real', 'real', 'real', 'real', 'real', 'real', 'real'])

# Undetermined parameter which is not considered yet (require the use of EIS curves to be calibrated)
C_dl = 2e7  # F.m-3. It is the volumetric double layer capacitance.

# GeneticAlgorithm parameters
algorithm_param = AlgorithmParams(  # 1 iteration takes 720s = 12min = 0.2h inside LIS cluster with 80 nodes.
    #                                 800 iterations takes 160h = 6.67j.
    max_num_iteration=800,  # It should be between 1000 and 1500, depending on the population_size, for a good compromise
    #                         between speed and precision.
    population_size=160,  # It should be between 100 and 200 for a good compromise between speed and precision.
    #                       Select a multiple of the available number of CPU cores for optimal multiprocessing.
    mutation_probability=0.33 / 12,  # Best found.
    elit_ratio=1 / 160,  # 1 elite should be enough. The value given in population_size should then be reused here.
    parents_portion=0.2,  # Best found. Default is 0.3.
    crossover_type='one_point',  # Best found. Default is 'uniform'.
    mutation_type='uniform_by_x',  # Best found. Default is 'uniform_by_center'.
    selection_type='roulette',  # Best found and default setting.
    max_iteration_without_improv=None)

# ________________________________________________Determined parameters_________________________________________________
"""
The determined parameters of the experimental stack, along with the different operating conditions, are imported from 
the calibration_module.py.
This section should remain unaltered for regular program usage.
"""

Tfc_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, i_pola_1, Aact, Hmem, Hcl, Hgdl, Hgc, Wgc, Lgc, \
    type_auxiliary, type_control, type_purge, type_display, type_plot, type_current, current_density, t_step, i_step, \
    delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, t_purge, max_step, n_gdl, i_exp1, U_exp1 \
    = determined_parameters(type_fuel_cell_1)

Tfc_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, i_pola_2, Aact, Hmem, Hcl, Hgdl, Hgc, Wgc, Lgc, \
    type_auxiliary, type_control, type_purge, type_display, type_plot, type_current, current_density, t_step, i_step, \
    delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, t_purge, max_step, n_gdl, i_exp2, U_exp2 \
    = determined_parameters(type_fuel_cell_2)


# _________________________________________________Function to minimize_________________________________________________

def pola_points(x):
    """
    This function calculates the maximum error between the simulated polarization curves and the experimental ones.
    It is the function to minimize using the genetic algorithm.

    Parameters
    ----------
    x: numpy.ndarray
        Undetermined parameters proposed by the genetic algorithm.

    Returns
    -------
    sim_error: float
        Maximum error between the simulated polarization curves and the experimental ones.
    """

    # Undetermined parameters proposed by the genetic algorithm.
    epsilon_gdl, epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch = x

    # Calculation of the model polarization curve
    Simulator1 = AlphaPEM(current_density, Tfc_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, t_step,
                          i_step, i_pola_1, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell_1, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)
    Simulator2 = AlphaPEM(current_density, Tfc_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, t_step,
                          i_step, i_pola_2, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell_2, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)

    # Calculation of the error between the simulated polarization curves and the experimental ones
    sim_error = calculate_simulation_error(Simulator1, U_exp1, i_exp1, Simulator2, U_exp2, i_exp2)
    return sim_error


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

    # GeneticAlgorithm execution
    unfinished = True  # As long as this value remains true, the execution of the GeneticAlgorithm restart.
    while unfinished:
        try:
            model = ga(dimension=12, variable_type=vartype, variable_boundaries=varbound,
                       algorithm_parameters=algorithm_param)
            # Multi processing is permitted.
            model.run(function=pola_points, set_function=ga.set_function_multiprocess(pola_points, n_jobs=-1))
        except ValueError:
            print("\nA ValueError has been detected due to a parameter combination that is challenging to satisfy "
                  "with the current max_step value, which is too high in the model. The program is relaunched.")
        except TypeError:
            print("\nTypeError has been detected due to a parameter combination that is challenging to satisfy "
                  "with the current max_step value, which is too high in the model. The program is relaunched.")
        except AssertionError:
            print("\nAssertionError has been detected due to a parameter combination that is challenging to satisfy "
                  "with the current max_step value, which is too high in the model. The program is relaunched.")
        else:  # Loop exit
            unfinished = False

    # Recovery of the results
    convergence = model.report
    solution = model.result
    epsilon_gdl, epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch = \
        solution['variable'][0], solution['variable'][1], solution['variable'][2], solution['variable'][3], \
            solution['variable'][4], solution['variable'][5], solution['variable'][6], solution['variable'][7], \
            solution['variable'][8], solution['variable'][9], solution['variable'][10], solution['variable'][11]
    sim_error = solution['score']

    # Print of the parameter calibration results
    print_calibration_results(convergence, epsilon_gdl, epsilon_mc, tau, epsilon_c, e,
                              Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error)

    # Save data in a text file
    save_calibration_results(convergence, epsilon_gdl, epsilon_mc, tau, epsilon_c, e,
                             Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error, type_fuel_cell_1)

    # Calculate, display and save the calibrated and experimental polarization curve
    #       Calculate the calibrated polarization curve
    Simulator1 = AlphaPEM(current_density, Tfc_1, Pa_des_1, Pc_des_1, Sa_1, Sc_1, Phi_a_des_1, Phi_c_des_1, t_step,
                          i_step, i_pola_1, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell_1, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)
    Simulator2 = AlphaPEM(current_density, Tfc_2, Pa_des_2, Pc_des_2, Sa_2, Sc_2, Phi_a_des_2, Phi_c_des_2, t_step,
                          i_step, i_pola_2, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, Aact, Hgdl, Hmem, Hcl, Hgc, Wgc,
                          Lgc, epsilon_gdl, tau, epsilon_mc, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim,
                          b_slim, a_switch, C_dl, max_step, n_gdl, t_purge, type_fuel_cell_2, type_current,
                          type_auxiliary, type_control, type_purge, type_display, type_plot)
    #       Display the calibrated and experimental polarization curve
    fig, ax = plt.subplots(figsize=(6, 6))
    Simulator1.Display(ax)
    Simulator2.Display(ax)
    #       Save the calibrated and experimental polarization curve
    subfolder_name = type_fuel_cell_1[:type_fuel_cell_1.rfind('_')] if type_fuel_cell_1.rfind('_') != -1 \
        else type_fuel_cell_1
    saving_instructions("results", subfolder_name, "step_current_syn_1.pdf", fig)

    # Algorithm time
    algo_time = time.time() - start_time
    print('Time of the algorithm in second :', algo_time)
