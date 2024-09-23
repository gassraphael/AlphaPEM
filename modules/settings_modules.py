# -*- coding: utf-8 -*-

"""This modul contains some of the required functions for the settings.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np


# ___________________________________________________Settings modules___________________________________________________

def stored_operating_inputs(type_fuel_cell):
    """This function gives the operating inputs which correspond to the given type_fuel_cell.

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
    i_max_pola : float
        Maximum current density for the polarization curve.
    """

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5":
        Tfc = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 1.5e5, 1.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        i_max_pola = 3.0e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "EH-31_2.0":
        Tfc = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        i_max_pola = 3.0e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "EH-31_2.25":
        Tfc = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.25e5, 2.25e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        i_max_pola = 3.0e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "EH-31_2.5":
        Tfc = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.5e5, 2.5e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        i_max_pola = 3.0e4  # A.m-2. It is the maximum current density for the polarization curve.

    # For LF fuel cell
    elif type_fuel_cell == "LF":
        Tfc = 80 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 101325, 101325  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 2.0, 1.5  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.84, 0.59  # It is the desired relative humidity.
        i_max_pola = 1.6e4  # A.m-2. It is the maximum current density for the polarization curve.

    # For other fuel cells
    else:
        raise ValueError('the type_fuel_cell given is not valid.')

    return Tfc, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, i_max_pola


def stored_physical_parameters(type_fuel_cell):
    """This function gives the physical parameters which correspond to the given type_fuel_cell.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.

    Returns
    -------
    Hcl : float
        Thickness of the catalyst layer in m.
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
    tau : float
        Pore structure coefficient.
    Hmem : float
            Thickness of the membrane in m.
    Hgdl : float
            Thickness of the gas diffusion layer in m.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_c : float
        Compression ratio of the GDL.
    Hgc : float
        Thickness of the gas channel in m.
    Wgc : float
        Width of the gas channel in m.
    Lgc : float
        Length of the gas channel in m.
    Aact : float
            Active area of the cell in m².
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
    C_dl : float
        Volumetric double layer capacitance in F.m-3.
    """

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
        # Catalyst layer
        Aact = 8.5e-3  # m². It is the active area of the catalyst layer.
        Hcl = 1e-5  # m. It is the thickness of the anode or cathode catalyst layer.
        epsilon_mc = 0.399  # It is the volume fraction of ionomer in the CL.
        tau = 1.016  # It is the pore structure coefficient, without units.
        # Membrane
        Hmem = 2e-5  # m. It is the thickness of the membrane.
        # Gas diffusion layer
        Hgdl = 2e-4  # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl = 0.701  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.271  # It is the compression ratio of the GDL.
        # Gas channel
        Hgc = 5e-4  # m. It is the thickness of the gas channel.
        Wgc = 4.5e-4  # m. It is the width of the gas channel.
        Lgc = 9.67  # m. It is the length of the gas channel.
        # Interaction parameters between water and PEMFC structure
        e = 5.0  # It is the capillary exponent
        # Voltage polarization
        Re = 5.70e-07  # ohm.m². It is the electron conduction resistance of the circuit.
        i0_c_ref = 2.79  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co = 27.2  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 1.61  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.05553, 0.10514, 0.63654  # It is the limit liquid saturation coefficients.
        C_dl = 2e7  # F.m-3. It is the volumetric double layer capacitance.

    # For LF fuel cell
    elif type_fuel_cell == "LF":
        # Catalyst layer
        Hcl = 1e-5  # m. It is the thickness of the anode or cathode catalyst layer.
        epsilon_mc = 0.27  # It is the volume fraction of ionomer in the CL.
        tau = 1.2  # It is the pore structure coefficient, without units.
        # Membrane
        Hmem = 5.08e-5  # m. It is the thickness of the membrane.
        # Gas diffusion layer
        Hgdl = 4.2e-4  # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl = 0.6  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.21  # It is the compression ratio of the GDL.
        # Gas channel.
        Hgc = 1e-3  # m. It is the thickness of the gas channel.
        Wgc = 8e-4  # m. It is the width of the gas channel.
        Wrib = 7e-4  # m. It is the rib of the gas channel.
        L1gc = 0.1  # m. It is the length of 1 channel.
        Ngc = 16  # . It is the number of channels in the gas channel. It is taken to have Aact = 25cm².
        Lgc = (L1gc + Wrib) * Ngc  # m. It is the length of the gas channel.
        Aact = (L1gc + 2 * Wrib) * (Wgc + Wrib) * Ngc + (L1gc + 2 * Wrib) * 7e-4  # m². It is the CL active area.
        # Interaction parameters between water and PEMFC structure
        e = 3.0  # It is the capillary exponent
        # Voltage polarization
        Re = 1e-6  # ohm.m². It is the electron conduction resistance of the circuit.
        i0_c_ref = 10  # A.m-2.It is the reference exchange current density at the cathode.
        kappa_co = 25  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 1.5  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0, 1, 1  # It is the limit liquid saturation coefficients.
        C_dl = 2e7  # F.m-3. It is the volumetric double layer capacitance.

    # For other fuel cells
    else:
        raise ValueError('the type_input given is not valid.')

    return Hcl, epsilon_mc, tau, Hmem, Hgdl, epsilon_gdl, epsilon_c, \
        Hgc, Wgc, Lgc, Aact, e, Re, i0_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_dl


def EIS_parameters(f_EIS):
    """This function gives the time parameters for the EIS_current density function.

    Parameters
    ----------
    f_EIS : tuple
        EIS parameters. It is a tuple containing the power of the initial frequency 'f_power_min_EIS':
        f_min_EIS = 10**f_power_min_EIS, the power of the final frequency 'f_power_max_EIS', the number of frequencies
        tested 'nb_f_EIS' and the number of points calculated per specific period 'nb_points_EIS'.

    Returns
    -------
    t_EIS : tuple
        EIS parameters. It is a tuple containing the initial EIS time after stack equilibrium 't0_EIS', a list of time
        parameters which gives the beginning of each frequency change 't_new_start_EIS', the final time 'tf_EIS', a list
        of time parameters which gives the estimated time for reaching equilibrium at each frequency
        'delta_t_break_EIS', and a list of time parameters which gives the estimated time for measuring the voltage
        response at each frequency 'delta_t_measurement_EIS'.
    """

    # Initialisation
    #       Frequencies
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS  # They are the frequency parameters for the EIS
    #                                                                    simulation.
    f = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)  # It is the tested frequencies
    nb_period_break_EIS, nb_period_measurement_EIS = 50, 50  # They are the number of temporal periods which are used
    #                                                          for break and for measurement. It is more accurate to use
    #                                                          periods than time as the frequency range is big.
    #       Time parameters
    delta_t_break_EIS = np.array([])  # It is the estimated time for reaching equilibrium at each frequency.
    delta_t_measurement_EIS = np.array([])  # It is the estimated time for measuring the voltage response.

    # Time parameters calculation
    t0_EIS = 1 / f[0]  # s. It is the simulation starting time.
    #                   [0, t0_EIS] is used to let the stack equilibrate to i_EIS.
    t_new_start_EIS = np.array([t0_EIS])  # It is a list of time parameters which gives the beginning of each frequency
    #                                   change.
    for i in range(nb_f_EIS):  # The goal is to measure nb_f_EIS periods of the signal in order to have precise enough values.
        T_i = 1 / (f[i]) # s. It is the period of the signal.
        if i < (nb_f_EIS - 1):
            delta_t_break_EIS = np.concatenate((delta_t_break_EIS, [nb_period_break_EIS * T_i]))
            delta_t_measurement_EIS = np.concatenate((delta_t_measurement_EIS, [nb_period_measurement_EIS * T_i]))
            next_start_EIS = t_new_start_EIS[i] + delta_t_break_EIS[i] + delta_t_measurement_EIS[i]
            t_new_start_EIS = np.concatenate((t_new_start_EIS, [next_start_EIS]))
        else:
            delta_t_break_EIS = np.concatenate((delta_t_break_EIS, [nb_period_break_EIS * T_i]))
            delta_t_measurement_EIS = np.concatenate((delta_t_measurement_EIS, [nb_period_measurement_EIS * T_i]))
            tf_EIS = t_new_start_EIS[-1] + delta_t_break_EIS[-1] + delta_t_measurement_EIS[-1]  # s. It is the
            #                                                                                     simulation ending time

    t_EIS = t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS
    return t_EIS
