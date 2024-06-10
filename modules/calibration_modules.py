# -*- coding: utf-8 -*-

"""This module contains some of the functions required for the parameter calibration.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import os
from colorama import Fore, Style
import numpy as np
from scipy.interpolate import interp1d

# Importing functions
from configuration.current_densities import polarization_current
from calibration.experimental_values import pola_exp_values


# _________________________________________________Calibration modules__________________________________________________

def determined_parameters(type_fuel_cell):
    """This function is used to determine the parameters of the fuel cell model for the calibration when a registered
    type_fuel_cell is considered.

    Parameters
    ----------
    type_fuel_cell : str
            Type of fuel cell configuration.

    Returns
    -------
    Tfc : float
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
    i_pola : float
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
    t_step : tuple
        Time parameters for the step_current density function.
        It is a tuple containing the initial time 't0_step', final time 'tf_step', loading time 'delta_t_load' and
        dynamic time for display 'delta_t_dyn'.
    i_step : tuple
        Current parameters for the step_current density function.
        It is a tuple containing the initial and final current density value 'i_ini' and 'i_final'.
    delta_pola : tuple
        Parameters for the polarization curve. It is a tuple containing the loading time
        'delta_t_load', the breaking time 'delta_t_break', the current density step 'delta_i', and the initial
        breaking time 'delta_t_ini'.
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
    max_step : float
            Maximum time step for the solver.
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
        Tfc = 74 + 273.15  # K. It is the temperature of the fuel cell.
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        if type_fuel_cell == "EH-31_1.5":
            Pa_des, Pc_des = 1.5e5, 1.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
            i_pola = 2.3e4
        elif type_fuel_cell == "EH-31_2.0":
            Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
            i_pola = 2.5e4
        elif type_fuel_cell == "EH-31_2.25":
            Pa_des, Pc_des = 2.25e5, 2.25e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
            i_pola = 2.85e4
        else:  # type_fuel_cell == "EH-31_2.5":
            Pa_des, Pc_des = 2.5e5, 2.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
            i_pola = 3.05e4
        #       Fuel cell physical parameters
        Aact = 8.5e-3  # m². It is the active area of the catalyst layer.
        Wgc = 4.5e-4  # m. It is the width of the gas channel.
        Lgc = 9.67  # m. It is the length of the gas channel.

        # Extrapolated physical parameters
        Hmem = 2e-5  # m. It is the thickness of the membrane.
        Hcl = 1e-5  # m. It is the thickness of the anode or cathode catalyst layer.
        Hgdl = 2e-4  # m. It is the thickness of the gas diffusion layer.
        Hgc = 5e-4  # m. It is the thickness of the gas channel.

        # Algorithm parameters for polarization curve generation
        type_auxiliary = "opened_anode"
        type_control = "no_control"
        type_purge = "no_purge"
        type_display = "multiple"
        type_plot = "final"
        type_current = "polarization"
        current_density = polarization_current
        t_step = np.nan, np.nan, np.nan, np.nan  # It is the time parameters for the step_current density function.
        i_step = np.nan, np.nan  # It is the current parameters for the step_current density function.
        delta_pola = 30, 30, 0.1e4, 1 * 60  # It is the parameters for the polarization curve.
        i_EIS, ratio_EIS = np.nan, np.nan  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
        f_EIS, t_EIS = np.nan, np.nan  # It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        max_step = 0.05  # It is good enough for having graphs without instabilities.
        n_gdl = int(Hgdl / Hcl / 2)  # It is the number of model points placed inside each GDL.

    elif type_fuel_cell == "BX_1.0" or type_fuel_cell == "BX_1.35":
        # Given values by the author
        #       Operating inputs
        Tfc = 80 + 273.15  # K. It is the temperature of the fuel cell.
        Sa, Sc = 1.5, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        if type_fuel_cell == "BX_1.0":
            Pa_des, Pc_des = 1.0 * 101325, 1.0 * 101325  # Pa. It is the desired pressures of the fuel gasS.
            Phi_a_des, Phi_c_des = 0.25, 0.25  # It is the desired relative humidity.
            i_pola = 2.1e4
        else:  # type_fuel_cell == "BX_1.35":
            Pa_des, Pc_des = 1.35 * 101325, 1.35 * 101325  # Pa. It is the desired pressures of the fuel gas.
            Phi_a_des, Phi_c_des = 1.0, 0.5  # It is the desired relative humidity.
            i_pola = 1.6e4
        #       Fuel cell physical parameters
        Aact = 0.005  # m². It is the active area of the catalyst layer.
        Hmem = 1e-5  # m. It is the thickness of the membrane.
        Hcl = 1e-5  # m. It is the thickness of the anode or cathode catalyst layer.
        Hgdl = 2e-4  # m. It is the thickness of the gas diffusion layer.
        Hgc = 5e-4  # m. It is the thickness of the gas channel.
        Wgc = 8e-4  # m. It is the width of the gas channel.

        # Extrapolated physical parameters
        Lgc = 3.0  # m. It is the length of the gas channel.

        # Algorithm parameters for polarization curve generation
        type_auxiliary = "opened_anode"
        type_control = "no_control"
        type_purge = "no_purge"
        type_display = "no_display"
        type_plot = "final"
        type_current = "polarization"
        current_density = polarization_current
        t_step = np.nan  # It is the time parameters for the step_current density function.
        i_step = np.nan, np.nan  # It is the current parameters for the step_current density function.
        delta_pola = 30, 30, 0.1e4, 60 * 60  # It is the parameters for the polarization curve.
        i_EIS, ratio_EIS = np.nan, np.nan  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
        f_EIS, t_EIS = np.nan, np.nan  # It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        max_step = 0.05  # It is good enough for having graphs without instabilities.
        n_gdl = int(Hgdl / Hcl / 2)  # It is the number of model points placed inside each GDL.

    elif type_fuel_cell == "LF":
        # Given values by the author
        #       Operating inputs
        Tfc = 80 + 273.15  # K. It is the temperature of the fuel cell.
        Pa_des, Pc_des = 101325, 101325  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 2.0, 1.5  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.84, 0.59  # It is the desired relative humidity.
        i_pola = 1.45e4
        #       Fuel cell physical parameters
        Hmem = 5.08e-5  # m. It is the thickness of the membrane.
        Hcl = 1e-5  # m. It is the thickness of the anode or cathode catalyst layer.
        Hgdl = 4.2e-4  # m. It is the thickness of the gas diffusion layer.
        Hgc = 1e-3  # m. It is the thickness of the gas channel.
        Wgc = 8e-4  # m. It is the width of the gas channel.

        # Extrapolated physical parameters
        Aact = 0.0025  # m². It is the active area of the catalyst layer.
        Lgc = 1.6  # m. It is the length of the gas channel.

        # Algorithm parameters for polarization curve generation
        type_auxiliary = "closed_anode"
        type_control = "no_control"
        type_purge = "no_purge"
        type_display = "no_display"
        type_plot = "final"
        type_current = "polarization"
        current_density = polarization_current
        t_step = np.nan  # It is the time parameters for the step_current density function.
        i_step = np.nan, np.nan  # It is the current parameters for the step_current density function.
        delta_pola = 30, 30, 0.1e4, 60 * 60  # It is the parameters for the polarization curve.
        i_EIS, ratio_EIS = np.nan, np.nan  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
        f_EIS, t_EIS = np.nan, np.nan  # It is the EIS parameters.
        t_purge = 0.6, 15  # s It is the purge time and the distance between two purges.
        max_step = 0.05  # It is good enough for having graphs without instabilities.
        n_gdl = int(Hgdl / Hcl / 2)  # It is the number of model points placed inside each GDL.

    else:
        ValueError("A correct type_fuel_cell should be given.")

    # Characteristic points of the experimental polarization curve
    i_exp, U_exp = pola_exp_values(type_fuel_cell)

    return (Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, i_pola, Aact, Hmem, Hcl, Hgdl, Hgc, Wgc, Lgc,
            type_auxiliary, type_control, type_purge, type_display, type_plot, type_current, current_density, t_step,
            i_step, delta_pola, i_EIS, ratio_EIS, t_EIS, f_EIS, t_purge, max_step, n_gdl, i_exp, U_exp)


def calculate_simulation_error(Simulator1, U_exp1, i_exp1, Simulator2, U_exp2, i_exp2):
    """This function is used to calculate the simulation maximal error between the experimental and the simulated
    polarization curves. Two simulations on different operating conditions and on the same stack, and so two set of
    experimental data, are considered as it is the minimum amount of data which is required for the calibration.

    Parameters
    ----------
    Simulator1 : AlphaPEM object
        PEM simulator which contains the simulation results for the first simulation.
    U_exp1 : numpy.ndarray
        Experimental values of the voltage for the first simulation.
    i_exp1 : numpy.ndarray
        Experimental values of the current density for the first simulation.
    Simulator2 : AlphaPEM object
        PEM simulator which contains the simulation results for the second simulation.
    U_exp2 : numpy.ndarray
        Experimental values of the voltage for the second simulation.
    i_exp2 : numpy.ndarray
        Experimental values of the current density for the second simulation.

    Returns
    -------
    sim_error : float
        Maximum error between the experimental and the simulated polarization curves in percentage.
    """

    # Recovery of ifc_1
    t1 = np.array(Simulator1.variables['t'])
    n1 = len(t1)
    ifc_t_1 = np.zeros(n1)
    for i in range(n1):  # Creation of ifc_t and conversion in A/cm²
        ifc_t_1[i] = Simulator1.operating_inputs['current_density'](t1[i], Simulator1.parameters) / 1e4
    # Recovery of ifc_2
    t2 = np.array(Simulator2.variables['t'])
    n2 = len(t2)
    ifc_t_2 = np.zeros(n2)
    for i in range(n2):  # Creation of ifc_t and conversion in A/cm²
        ifc_t_2[i] = Simulator2.operating_inputs['current_density'](t2[i], Simulator2.parameters) / 1e4

    # Polarisation curve point recovery after stack stabilisation for Simulator1
    delta_t_load, delta_t_break, delta_i, delta_t_ini = Simulator1.parameters['delta_pola']
    nb_loads1 = int(Simulator1.parameters['i_pola'] / delta_i + 1)  # Number of load which are made
    ifc_discretized1 = np.zeros(nb_loads1)
    Ucell_discretized1 = np.zeros(nb_loads1)
    for i in range(nb_loads1):
        t_load = delta_t_ini + (i + 1) * (delta_t_load + delta_t_break) - delta_t_break / 10  # time for the measurement
        idx1 = (np.abs(t1 - t_load)).argmin()  # the corresponding index
        ifc_discretized1[i] = ifc_t_1[idx1]  # the last value at the end of each load
        Ucell_discretized1[i] = Simulator1.variables['Ucell'][idx1]  # the last value at the end of each load
    # Polarisation curve point recovery after stack stabilisation for Simulator2
    delta_t_load, delta_t_break, delta_i, delta_t_ini = Simulator2.parameters['delta_pola']
    nb_loads2 = int(Simulator2.parameters['i_pola'] / delta_i + 1)  # Number of load which are made
    ifc_discretized2 = np.zeros(nb_loads2)
    Ucell_discretized2 = np.zeros(nb_loads2)
    for i in range(nb_loads2):
        t_load = delta_t_ini + (i + 1) * (delta_t_load + delta_t_break) - delta_t_break / 10  # time for the measurement
        idx2 = (np.abs(t2 - t_load)).argmin()  # the corresponding index
        ifc_discretized2[i] = ifc_t_2[idx2]  # the last value at the end of each load
        Ucell_discretized2[i] = Simulator2.variables['Ucell'][idx2]  # the last value at the end of each load

    # Interpolation of experimental points to match model points for Simulator 1
    i_fc_reduced1 = ifc_discretized1[(ifc_discretized1 >= i_exp1[0]) & (ifc_discretized1 <= i_exp1[-1])]
    Ucell_reduced1 = Ucell_discretized1[(ifc_discretized1 >= i_exp1[0]) & (ifc_discretized1 <= i_exp1[-1])]
    U_exp_interpolated1 = interp1d(i_exp1, U_exp1, kind='linear')(i_fc_reduced1)
    # Interpolation of experimental points to match model points for Simulator 2
    i_fc_reduced2 = ifc_discretized2[(ifc_discretized2 >= i_exp2[0]) & (ifc_discretized2 <= i_exp2[-1])]
    Ucell_reduced2 = Ucell_discretized2[(ifc_discretized2 >= i_exp2[0]) & (ifc_discretized2 <= i_exp2[-1])]
    U_exp_interpolated2 = interp1d(i_exp2, U_exp2, kind='linear')(i_fc_reduced2)

    # Distance between the simulated and the experimental polarization curves.
    if np.isnan(Ucell_reduced1).any() or np.isnan(Ucell_reduced2).any():
        sim_error = 1e5  # If the mass loss arrives before i_exp[-1], some value would be equal to NaN.
    else:
        sim_error = (np.max(np.abs(Ucell_reduced1 - U_exp_interpolated1) / U_exp_interpolated1 * 100)
                     + np.max(np.abs(Ucell_reduced2 - U_exp_interpolated2) / U_exp_interpolated2 * 100)) / 2  # in %.

    return sim_error


def print_calibration_results(convergence, epsilon_gdl, epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c,
                              a_slim, b_slim, a_switch, sim_error):
    """This function is used to print the calibration results.

    Parameters
    ----------
    convergence : dict
        A dictionary generated by the GeneticAlgorithm2 model's report method. It contains information about the
        convergence of the genetic algorithm used for optimizing the undetermined parameters in the calibration files.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
    tau : float
        Pore structure coefficient.
    epsilon_c : float
        Compression ratio of the GDL.
    e : float
        Capillary exponent.
    Re : float
        Electron conduction resistance of the circuit in ohm.m².
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
    print("\nThe optimized epsilon_gdl: ", epsilon_gdl)
    print("The optimized epsilon_mc: ", epsilon_mc)
    print("The optimized tau: ", tau)
    print("The optimized epsilon_c: ", epsilon_c)
    print("The optimized e: ", e)
    print("The optimized Re: ", Re)
    print("The optimized i0_c_ref: ", i0_c_ref)
    print("The optimized kappa_co: ", kappa_co)
    print("The optimized kappa_c: ", kappa_c)
    print("The optimized a_slim: ", a_slim)
    print("The optimized b_slim: ", b_slim)
    print("The optimized a_switch: ", a_switch)
    print(Fore.RED + "\nThe max simulation error is: ", sim_error, "%")
    print(Style.RESET_ALL)


def save_calibration_results(convergence, epsilon_gdl, epsilon_mc, tau, epsilon_c, e, Re, i0_c_ref, kappa_co, kappa_c,
                             a_slim, b_slim, a_switch, sim_error, type_fuel_cell):
    """This function is used to save in a text file the calibration results.

    Parameters
    ----------
    convergence : dict
        A dictionary generated by the GeneticAlgorithm2 model's report method. It contains information about the
        convergence of the genetic algorithm used for optimizing the undetermined parameters in the calibration files.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
    tau : float
        Pore structure coefficient.
    epsilon_c : float
        Compression ratio of the GDL.
    e : float
        Capillary exponent.
    Re : float
        Electron conduction resistance of the circuit in ohm.m².
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
                   "\nThe optimized epsilon_gdl: " + str(epsilon_gdl) +
                   "\nThe optimized epsilon_mc: " + str(epsilon_mc) +
                   "\nThe optimized tau: " + str(tau) +
                   "\nThe optimized epsilon_c: " + str(epsilon_c) +
                   "\nThe optimized e: " + str(e) +
                   "\nThe optimized Re: " + str(Re) +
                   "\nThe optimized i0_c_ref: " + str(i0_c_ref) +
                   "\nThe optimized kappa_co: " + str(kappa_co) +
                   "\nThe optimized kappa_c: " + str(kappa_c) +
                   "\nThe optimized a_slim: " + str(a_slim) +
                   "\nThe optimized b_slim: " + str(b_slim) +
                   "\nThe optimized a_switch: " + str(a_switch) +
                   "\nThe max simulation error is: " + str(sim_error) + "%" +
                   "\nHere the algorithm works with RG2&3 (global)")
