# -*- coding: utf-8 -*-

"""This module contains some of the functions required for the parameter calibration.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import os
from colorama import Fore, Style
import numpy as np

# Importing functions
from configuration.current_densities import polarization_current_for_calibration
from calibration.experimental_values import pola_exp_values_calibration


# _________________________________________________Calibration modules__________________________________________________

def parameters_for_calibration(type_fuel_cell):
    """This function is used to determine the parameters of the fuel cell model for the calibration when a registered
    type_fuel_cell is considered.

    Parameters
    ----------
    type_fuel_cell : str
            Type of fuel cell configuration.

    Returns
    -------
    T_des : float
            Desired fuel cell temperature in Kelvin.
    Pa_des : float
        Desired anode pressure in Pascal.
    Pc_des : float
        Desired cathode pressure in Pascal.
    Sa : float
        Stoichiometric ratio of hydrogen.
    Sc : float
        Stoichiometric ratio of oxygen.
    Phi_a_des : float
        Desired anode relative humidity.
    Phi_c_des : float
        Desired cathode relative humidity.
    i_max_pola : float
        Maximum current density for the polarization curve.
    Aact : float
        Active area of the cell in m².
    Hmem : float
        Thickness of the membrane in m.
    Hcl : float
        Thickness of the catalyst layer in m.
    Hgc : float
        Thickness of the gas channel in m.
    Wgc : float
        Width of the gas channel in m.
    Lgc : float
        Length of the gas channel in m.
    type_auxiliary : str
        Type of auxiliary system.
    type_control : str
        Type of control system.
    type_purge : str
        Type of purge system.
    type_display : str
        Type of display.
    type_plot : str
        Type of plot.
    type_current : str
        Type of current density function.
    current_density : function
        Current density evolution over time. It is a function of time and parameters dictionary.
    step_current_parameters : dict
        Parameters for the step current density. It is a dictionary containing:
        - 'delta_t_ini_step': the initial time (in seconds) at zero current density for the stabilisation of the
        internal states,
        - 'delta_t_load_step': the loading time (in seconds) for the step current density function, from 0 to
        i_step,
        - 'delta_t_break_step': the time (in seconds) at i_step current density for the stabilisation of the
        internal states,
        - 'i_step': the current density (in A.m-2) for the step current density function,
        - 'delta_t_dyn_step': the time (in seconds) for dynamic display of the step current density function.
    pola_current_parameters : dict
        Parameters for the polarization current density. It is a dictionary containing:
        - 'delta_t_ini_pola': the initial time (in seconds) at zero current density for the stabilisation of the
        internal states,
        - 'delta_t_load_pola': the loading time (in seconds) for one step current of the polarisation current
        density function,
        - 'delta_t_break_pola': the breaking time (in seconds) for one step current, for the stabilisation of the
        internal states,
        - 'delta_i_pola': the current density step (in A.m-2) for the polarisation current density function.
        - 'i_max_pola': the maximum current density (in A.m-2) for the polarization curve.
    pola_current_for_cali_parameters : dict
        Parameters for the polarization current density for calibration. It is a dictionary containing:
        - 'delta_t_ini_pola_cali': the initial time (in seconds) at zero current density for the stabilisation of
        the internal states,
        - 'delta_t_load_pola_cali': the loading time (in seconds) for one step current of the polarisation current
        density function,
        - 'delta_t_break_pola_cali': the breaking time (in seconds) for one step current, for the stabilisation of
        the internal states.
    i_EIS : float
        Current for which a ratio_EIS perturbation is added.
    ratio_EIS : float
        Value of the perturbation on the current density for building the EIS curve.
    t_EIS : tuple
        EIS parameters. It is a tuple containing the initial EIS time after stack equilibrium 't0_EIS', a list of time
        parameters which gives the beginning of each frequency change 't_new_start', the final time 'tf_EIS', a list of
        time parameters which gives the estimated time for reaching equilibrium at each frequency 'delta_t_break_EIS',
        and a list of time parameters which gives the estimated time for measuring the voltage response at each
        frequency 'delta_t_measurement_EIS'.
    f_EIS : tuple
        EIS parameters. It is a tuple containing the power of the initial frequency
        'f_power_min': f_min = 10**f_power_min, the power of the final frequency 'f_power_max', the number of
        frequencies tested 'nb_f' and the number of points calculated per specific period 'nb_points'.
    t_purge : tuple
        Time parameters for purging the system.
        It is the purge time interval 'purge_time' and the time between two purges 'delta_purge'.
    n_gdl : int
        Number of points considered in the GDL.
    i_exp : numpy.ndarray
        Experimental values of the current density.
    U_exp : numpy.ndarray
        Experimental values of the voltage.
    """

    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
        # Given values by the author
        #       Operating inputs
        T_des = 74 + 273.15  # K. It is the temperature of the fuel cell.
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        if type_fuel_cell == "EH-31_1.5":
            Pa_des, Pc_des = 1.5e5, 1.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        elif type_fuel_cell == "EH-31_2.0":
            Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        elif type_fuel_cell == "EH-31_2.25":
            Pa_des, Pc_des = 2.25e5, 2.25e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        else:  # type_fuel_cell == "EH-31_2.5":
            Pa_des, Pc_des = 2.5e5, 2.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        #       Fuel cell physical parameters
        Aact = 8.5e-3  # m². It is the active area of the catalyst layer.
        Wgc = 4.5e-4  # m. It is the width of the gas channel.
        Lgc = 9.67  # m. It is the length of the gas channel.

        # Extrapolated physical parameters
        Hgdl = 2e-4  # m. It is the thickness of the gas diffusion layer.
        Hmpl = 3e-5  # m. It is the thickness of the microporous layer.
        epsilon_mpl = 0.4  # It is the porosity of the microporous layer.
        Hgc = 5e-4  # m. It is the thickness of the gas channel.

        # Estimated undetermined parameters for the initialisation
        #   Gas diffusion layer
        epsilon_gdl = 0.7943  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.2  # It is the compression ratio of the GDL.
        #   Catalyst layer
        Hcl = 8e-6  # m. It is the thickness of the anode or cathode catalyst layer.
        epsilon_mc = 0.2111  # It is the volume fraction of ionomer in the CL.
        #   Membrane
        Hmem = 1.5e-5  # m. It is the thickness of the membrane.
        #   Interaction parameters between water and PEMFC structure
        e = 3.0  # It is the capillary exponent
        #   Voltage polarization
        i0_c_ref = 14.86  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co = 1  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 0.6386  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.05553, 0.10514, 0.63654  # It is the limit liquid saturation coefficients.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
        estimated_undetermined_parameters_for_initialisation = {'epsilon_gdl': epsilon_gdl, 'epsilon_c': epsilon_c,
                                                                'Hcl': Hcl, 'epsilon_mc': epsilon_mc, 'Hmem': Hmem,
                                                                'e': e, 'i0_c_ref': i0_c_ref, 'kappa_co': kappa_co,
                                                                'kappa_c': kappa_c, 'a_slim': a_slim, 'b_slim': b_slim,
                                                                'a_switch': a_switch, 'C_scl': C_scl}

        # Algorithm parameters for polarization curve generation
        type_auxiliary = "forced-convective_cathode_with_flow-through_anode"
        type_control = "no_control"
        type_purge = "no_purge"
        type_display = "no_display"
        type_plot = "fixed"
        type_current = "polarization_for_cali"
        current_density = polarization_current_for_calibration
        delta_t_ini_step = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_step = 1e-15  # (s). Loading time for the step current density function, from 0 to i_step.
        delta_t_break_step = 0  # (s). Time at i_step current density for the stabilisation of the internal states.
        i_step = 0  # (A.m-2). Current density for the step current density function.
        step_current_parameters = {'delta_t_ini_step': delta_t_ini_step, 'delta_t_load_step': delta_t_load_step,
                                   'delta_t_break_step': delta_t_break_step, 'i_step': i_step}
        delta_t_ini_pola = 30 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        delta_i_pola = 0.05e4  # (A.m-2). Current density step for the polarisation current density function.
        i_max_pola = 1.7e4  # (A.m-2). It is the maximum current density for the polarization curve.
        pola_current_parameters = {'delta_t_ini_pola': delta_t_ini_pola, 'delta_t_load_pola': delta_t_load_pola,
                                   'delta_t_break_pola': delta_t_break_pola, 'delta_i_pola': delta_i_pola,
                                   'i_max_pola': i_max_pola}
        delta_t_ini_pola_cali = 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola_cali = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola_cali = 10 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        pola_current_for_cali_parameters = {'delta_t_ini_pola_cali': delta_t_ini_pola_cali,
                                            'delta_t_load_pola_cali': delta_t_load_pola_cali,
                                            'delta_t_break_pola_cali': delta_t_break_pola_cali}
        i_EIS, ratio_EIS = np.nan, np.nan  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
        f_EIS, t_EIS = np.nan, np.nan  # It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        n_gdl = int(Hgdl / Hcl / 4)  # It is the number of model points placed inside each GDL.

    elif type_fuel_cell == "LF":
        # Given values by the author
        #       Operating inputs
        T_des = 80 + 273.15  # K. It is the temperature of the fuel cell.
        Pa_des, Pc_des = 101325, 101325  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 2.0, 1.5  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.84, 0.59  # It is the desired relative humidity.
        #       Fuel cell physical parameters
        Hmem = 5.08e-5  # m. It is the thickness of the membrane.
        Hcl = 1e-5  # m. It is the thickness of the anode or cathode catalyst layer.
        Hgdl = 4.2e-4  # m. It is the thickness of the gas diffusion layer.
        Hgc = 1e-3  # m. It is the thickness of the gas channel.
        Wgc = 8e-4  # m. It is the width of the gas channel.

        # Extrapolated physical parameters
        Aact = 0.0025  # m². It is the active area of the catalyst layer.
        Lgc = 1.6  # m. It is the length of the gas channel.
        Hmpl = 3e-5  # m. It is the thickness of the microporous layer.
        epsilon_mpl = 0.4  # It is the porosity of the microporous layer.

        # Estimated undetermined parameters for the initialisation
        # Catalyst layer
        epsilon_mc = 0.27  # It is the volume fraction of ionomer in the CL.
        # Gas diffusion layer
        epsilon_gdl = 0.6  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.21  # It is the compression ratio of the GDL.
        # Interaction parameters between water and PEMFC structure
        e = 3.0  # It is the capillary exponent
        # Voltage polarization
        i0_c_ref = 10  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co = 25  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 1.5  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0, 1, 1  # It is the limit liquid saturation coefficients.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.
        estimated_undetermined_parameters_for_initialisation = {'epsilon_gdl': epsilon_gdl, 'epsilon_c': epsilon_c,
                                                                'Hcl': Hcl, 'epsilon_mc': epsilon_mc, 'Hmem': Hmem,
                                                                'e': e, 'i0_c_ref': i0_c_ref, 'kappa_co': kappa_co,
                                                                'kappa_c': kappa_c, 'a_slim': a_slim, 'b_slim': b_slim,
                                                                'a_switch': a_switch, 'C_scl': C_scl}

        # Algorithm parameters for polarization curve generation
        type_auxiliary = "forced-convective_cathode_with_anodic_recirculation"
        type_control = "no_control"
        type_purge = "no_purge"
        type_display = "no_display"
        type_plot = "fixed"
        type_current = "polarization_for_cali"
        current_density = polarization_current_for_calibration
        delta_t_ini_step = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_step = 0  # (s). Loading time for the step current density function, from 0 to i_step.
        delta_t_break_step = 0  # (s). Time at i_step current density for the stabilisation of the internal states.
        i_step = 0  # (A.m-2). Current density for the step current density function.
        step_current_parameters = {'delta_t_ini_step': delta_t_ini_step, 'delta_t_load_step': delta_t_load_step,
                                   'delta_t_break_step': delta_t_break_step, 'i_step': i_step}
        delta_t_ini_pola = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        delta_i_pola = 0.1e4  # (A.m-2). Current density step for the polarisation current density function.
        pola_current_parameters = {'delta_t_ini_pola': delta_t_ini_pola, 'delta_t_load_pola': delta_t_load_pola,
                                   'delta_t_break_pola': delta_t_break_pola, 'delta_i_pola': delta_i_pola}
        delta_t_ini_pola_cali = 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
        delta_t_load_pola_cali = 30  # (s). Loading time for one step current of the polarisation current density function.
        delta_t_break_pola_cali = 10 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
        pola_current_for_cali_parameters = {'delta_t_ini_pola_cali': delta_t_ini_pola_cali,
                                            'delta_t_load_pola_cali': delta_t_load_pola_cali,
                                            'delta_t_break_pola_cali': delta_t_break_pola_cali}
        i_EIS, ratio_EIS = np.nan, np.nan  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
        f_EIS, t_EIS = np.nan, np.nan  # It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        n_gdl = int(Hgdl / Hcl / 4)  # It is the number of model points placed inside each GDL.

    else:
        ValueError("A correct type_fuel_cell should be given.")

    # Characteristic points of the experimental polarization curve
    i_exp, U_exp = pola_exp_values_calibration(type_fuel_cell)

    return (T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, step_current_parameters, pola_current_parameters,
            pola_current_for_cali_parameters, Aact, Hgdl, Hmpl, Hgc, Wgc, Lgc, epsilon_mpl,
            estimated_undetermined_parameters_for_initialisation, type_auxiliary, type_control, type_purge,
            type_display, type_plot, type_current, current_density, i_EIS, ratio_EIS, t_EIS, f_EIS, t_purge, n_gdl,
            i_exp, U_exp)


def calculate_simulation_error(Simulator_1, U_exp_1, i_exp_1, Simulator_2, U_exp_2, i_exp_2):
    """This function is used to calculate the simulation maximal error between the experimental and the simulated
    polarization curves. Two simulations on different operating conditions and on the same stack, and so two set of
    experimental data, are considered as it is the minimum amount of data which is required for the calibration.

    Parameters
    ----------
    Simulator_1 : AlphaPEM object
        PEM simulator which contains the simulation results for the first simulation.
    U_exp_1 : numpy.ndarray
        Experimental values of the voltage for the first simulation.
    i_exp_1 : numpy.ndarray
        Experimental values of the current density for the first simulation.
    Simulator_2 : AlphaPEM object
        PEM simulator which contains the simulation results for the second simulation.
    U_exp_2 : numpy.ndarray
        Experimental values of the voltage for the second simulation.
    i_exp_2 : numpy.ndarray
        Experimental values of the current density for the second simulation.

    Returns
    -------
    sim_error : float
        Maximum error between the experimental and the simulated polarization curves in percentage.
    """

    # Recovery of ifc_1
    t1 = np.array(Simulator_1.variables['t'])
    n1 = len(t1)
    ifc_t_1 = np.zeros(n1)
    for i in range(n1):  # Creation of ifc_t
        ifc_t_1[i] = Simulator_1.operating_inputs['current_density'](t1[i], Simulator_1.parameters)
    # Recovery of ifc_2
    t2 = np.array(Simulator_2.variables['t'])
    n2 = len(t2)
    ifc_t_2 = np.zeros(n2)
    for i in range(n2):  # Creation of ifc_t
        ifc_t_2[i] = Simulator_2.operating_inputs['current_density'](t2[i], Simulator_2.parameters)

    # Polarisation curve point recovery after stack stabilisation for Simulator1
    #   Extraction of the parameters
    #       The initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola_cali_1 = Simulator_1.parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']  # (s).
    #       The loading time for one step current of the polarisation current density function.
    delta_t_load_pola_cali_1 = Simulator_1.parameters['pola_current_for_cali_parameters']['delta_t_load_pola_cali']  # (s).
    #       The breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola_cali_1 = Simulator_1.parameters['pola_current_for_cali_parameters']['delta_t_break_pola_cali']  # (s).
    #   Calculation
    nb_loads1 = len(i_exp_1)  # Number of load which are made
    delta_t_cali_1 = delta_t_load_pola_cali_1 + delta_t_break_pola_cali_1  # s. It is the time of one load.
    ifc_discretized1 = np.zeros(nb_loads1)
    Ucell_discretized1 = np.zeros(nb_loads1)
    for i in range(nb_loads1):
        t_load_1 = delta_t_ini_pola_cali_1 + (i + 1) * delta_t_cali_1 # time for measurement
        idx1 = (np.abs(t1 - t_load_1)).argmin()  # the corresponding index
        ifc_discretized1[i] = ifc_t_1[idx1]  # the last value at the end of each load
        Ucell_discretized1[i] = Simulator_1.variables['Ucell'][idx1]  # the last value at the end of each load
    # Polarisation curve point recovery after stack stabilisation for Simulator2
    #   Extraction of the parameters
    #       The initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola_cali_2 = Simulator_2.parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']  # (s).
    #       The loading time for one step current of the polarisation current density function.
    delta_t_load_pola_cali_2 = Simulator_2.parameters['pola_current_for_cali_parameters']['delta_t_load_pola_cali']  # (s).
    #       The breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola_cali_2 = Simulator_2.parameters['pola_current_for_cali_parameters']['delta_t_break_pola_cali']  # (s).
    #   Calculation
    nb_loads2 = len(i_exp_2)  # Number of load which are made
    delta_t_cali_2 = delta_t_load_pola_cali_2 + delta_t_break_pola_cali_2  # s. It is the time of one load.
    ifc_discretized2 = np.zeros(nb_loads2)
    Ucell_discretized2 = np.zeros(nb_loads2)
    for i in range(nb_loads2):
        t_load_2 = delta_t_ini_pola_cali_2 + (i + 1) * delta_t_cali_2 # time for measurement
        idx2 = (np.abs(t2 - t_load_2)).argmin()  # the corresponding index
        ifc_discretized2[i] = ifc_t_2[idx2]  # the last value at the end of each load
        Ucell_discretized2[i] = Simulator_2.variables['Ucell'][idx2]  # the last value at the end of each load

    # Distance between the simulated and the experimental polarization curves.
    sim_error = (np.max(np.abs(Ucell_discretized1 - U_exp_1) / U_exp_1 * 100)
                     + np.max(np.abs(Ucell_discretized2 - U_exp_2) / U_exp_2 * 100)) / 2  # in %.

    return sim_error


def print_calibration_results(convergence, ga_instance, Hcl, Hmem, epsilon_gdl, epsilon_mc, epsilon_c, e, i0_c_ref,
                              kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error):
    """This function is used to print the calibration results.

    Parameters
    ----------
    convergence : dict
        A dictionary generated by the GeneticAlgorithm2 model's report method. It contains information about the
        convergence of the genetic algorithm used for optimizing the undetermined parameters in the calibration files.
    ga_instance : PyGAD object
        An instance of the PyGAD library, which is used to perform the optimization.
    Hcl : float
        Thickness of the catalyst layer in m.
    Hmem : float
        Thickness of the membrane in m.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
    epsilon_c : float
        Compression ratio of the GDL.
    e : float
        Capillary exponent.
    i0_c_ref : float
        Reference exchange current density at the cathode in A.m-2.
    kappa_co : float
        Crossover correction coefficient in mol.m-1.s-1.Pa-1.
    kappa_c : float
        Overpotential correction exponent.
    a_slim : float
        One of the limit liquid saturation coefficients: the slop of slim function.
    b_slim : float
        One of the limit liquid saturation coefficients: the intercept of slim function.
    a_switch : float
        One of the limit liquid saturation coefficients: the slop of s_switch function.
    sim_error : float
        Maximum error between the experimental and the simulated polarization curves in percentage.
    """

    print("The convergence is:\n", convergence)
    print("\nThe optimized Hcl: ", Hcl)
    print("The optimized Hmem: ", Hmem)
    print("The optimized epsilon_gdl: ", epsilon_gdl)
    print("The optimized epsilon_mc: ", epsilon_mc)
    print("The optimized epsilon_c: ", epsilon_c)
    print("The optimized e: ", e)
    print("The optimized i0_c_ref: ", i0_c_ref)
    print("The optimized kappa_co: ", kappa_co)
    print("The optimized kappa_c: ", kappa_c)
    print("The optimized a_slim: ", a_slim)
    print("The optimized b_slim: ", b_slim)
    print("The optimized a_switch: ", a_switch)
    print(Fore.RED + "\nThe max simulation error is: ", sim_error, "%")
    print(Style.RESET_ALL)
    if ga_instance.best_solution_generation != -1:
        print(f"Best fitness value reached after {ga_instance.best_solution_generation} generations.")


def save_calibration_results(convergence, ga_instance, Hcl, Hmem, epsilon_gdl, epsilon_mc, epsilon_c, e, i0_c_ref,
                             kappa_co, kappa_c, a_slim, b_slim, a_switch, sim_error, type_fuel_cell):
    """This function is used to save in a text file and a PyGAD file the calibration results.

    Parameters
    ----------
    convergence : dict
        A dictionary generated by the GeneticAlgorithm2 model's report method. It contains information about the
        convergence of the genetic algorithm used for optimizing the undetermined parameters in the calibration files.
    ga_instance : PyGAD object
        An instance of the PyGAD library, which is used to perform the optimization.
    Hcl : float
        Thickness of the catalyst layer in m.
    Hmem : float
        Thickness of the membrane in m.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
    epsilon_c : float
        Compression ratio of the GDL.
    e : float
        Capillary exponent.
    i0_c_ref : float
        Reference exchange current density at the cathode in A.m-2.
    kappa_co : float
        Crossover correction coefficient in mol.m-1.s-1.Pa-1.
    kappa_c : float
        Overpotential correction exponent.
    a_slim : float
        One of the limit liquid saturation coefficients: the slop of slim function.
    b_slim : float
        One of the limit liquid saturation coefficients: the intercept of slim function.
    a_switch : float
        One of the limit liquid saturation coefficients: the slop of s_switch function.
    sim_error : float
        Maximum error between the experimental and the simulated polarization curves in percentage.
    type_fuel_cell : str
        Type of fuel cell configuration.

    Returns
    -------
    None
        The function saves the calibration results in a text file and a PyGAD file.
    """

    root_folder, filename = "results", "parameter_calibration_1.txt"
    subfolder_name = type_fuel_cell[:type_fuel_cell.rfind('_')] if type_fuel_cell.rfind('_') != -1 else type_fuel_cell
    counter = 1
    # Create the folder if necessary
    folder_name = os.path.join(root_folder, subfolder_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    # Create the file without erasing the previous ones
    while os.path.isfile(os.path.join(folder_name, filename)):
        counter += 1
        filename = "parameter_calibration_" + str(counter) + ".txt"
    # Write information
    file_path = os.path.join(folder_name, filename)
    with open(file_path, "w") as file:
        file.write("The convergence is: " + str(convergence) +
                   "\nThe optimized Hcl: " + str(Hcl) +
                   "\nThe optimized Hmem: " + str(Hmem) +
                   "\nThe optimized epsilon_gdl: " + str(epsilon_gdl) +
                   "\nThe optimized epsilon_mc: " + str(epsilon_mc) +
                   "\nThe optimized epsilon_c: " + str(epsilon_c) +
                   "\nThe optimized e: " + str(e) +
                   "\nThe optimized i0_c_ref: " + str(i0_c_ref) +
                   "\nThe optimized kappa_co: " + str(kappa_co) +
                   "\nThe optimized kappa_c: " + str(kappa_c) +
                   "\nThe optimized a_slim: " + str(a_slim) +
                   "\nThe optimized b_slim: " + str(b_slim) +
                   "\nThe optimized a_switch: " + str(a_switch) +
                   "\nThe max simulation error is: " + str(sim_error) + "%" +
                   "\nHere the algorithm works with RG2&3 (global)")
        if ga_instance.best_solution_generation != -1:
            file.write(f"\nBest fitness value reached after {ga_instance.best_solution_generation} generations.")
    # Save the GA instance (the name is without extension).
    ga_instance.save(filename=os.path.join(folder_name, "parameter_calibration_" + str(counter)))
    # Delete the ongoing GA file because the calibration ended.
    if os.path.isfile('parameter_calibration_ongoing.pkl'):
        os.remove('parameter_calibration_ongoing.pkl')