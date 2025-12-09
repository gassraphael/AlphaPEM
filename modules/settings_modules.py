# -*- coding: utf-8 -*-

"""This modul contains some of the required functions for the settings.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np


# ___________________________________________________Settings modules___________________________________________________

def stored_operating_inputs(type_fuel_cell, voltage_zone):
    """This function gives the operating inputs which correspond to the given type_fuel_cell.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.
    voltage_zone : str
        Zone of the polarization curve which is considered.

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
    """

    # For the ZSW Generic Stack fuel cell
    if type_fuel_cell == "ZSW-GenStack":
        T_des = 68 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.2e5, 2.0e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.500e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.700e4 # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41":
        T_des = 68 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 1.61e5, 1.41e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7 # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.200e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 0.700e4 # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81":
        T_des = 68 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.01e5, 1.81e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.500e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.300e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2":
        T_des = 68 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.4e5, 2.2e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.500e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.900e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6":
        T_des = 68 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.8e5, 2.6e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.500e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.900e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "ZSW-GenStack_T_62":
        T_des = 62 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.2e5, 2.0e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.500e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.500e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "ZSW-GenStack_T_76":
        T_des = 76 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.2e5, 2.0e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.500e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.100e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "ZSW-GenStack_T_84":
        T_des = 84 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.2e5, 2.0e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.000e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 0.700e4  # A.m-2. It is the maximum current density for the polarization curve.

    # For EH-31 fuel cell
    elif type_fuel_cell == "EH-31_1.5":
        T_des = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 1.5e5, 1.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        y_H2_in = 1 # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.300e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.700e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "EH-31_2.0":
        T_des = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        y_H2_in = 1  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.500e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.300e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "EH-31_2.25":
        T_des = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.25e5, 2.25e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        y_H2_in = 1  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 2.800e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.700e4  # A.m-2. It is the maximum current density for the polarization curve.
    elif type_fuel_cell == "EH-31_2.5":
        T_des = 74 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des, Pc_des = 2.5e5, 2.5e5  # Pa. It is the desired pressures of the fuel gas.
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        y_H2_in = 1  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full":
            i_max_pola = 3.000e4  # A.m-2. It is the maximum current density for the polarization curve.
        elif voltage_zone == "before_voltage_drop":
            i_max_pola = 1.600e4  # A.m-2. It is the maximum current density for the polarization curve.

    # For other fuel cells
    else:
        raise ValueError('the type_fuel_cell given is not valid.')

    return T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, i_max_pola


def stored_physical_parameters(type_fuel_cell):
    """This function gives the physical parameters which correspond to the given type_fuel_cell.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.

    Returns
    -------
    Hacl : float
        Thickness of the anode catalyst layer in m.
    Hccl : float
        Thickness of the cathode catalyst layer in m.
    epsilon_mc : float
        Volume fraction of ionomer in the CL.
    Hmem : float
            Thickness of the membrane in m.
    Hgdl : float
            Thickness of the gas diffusion layer in m.
    epsilon_gdl : float
        Anode/cathode GDL porosity.
    epsilon_c : float
        Compression ratio of the GDL.
    Hagc : float
        Thickness of the anode gas channel in m.
    Hcgc : float
        Thickness of the cathode gas channel in m.
    Wagc : float
        Width of the anode gas channel in m.
    Wcgc : float
        Width of the cathode gas channel in m.
    Lgc : float
        Length of the gas channel in m.
    Aact : float
            Active area of the cell in m².
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
    C_dl : float
        Volumetric double layer capacitance in F.m-3.
    """

    # For the ZSW Generic Stack fuel cell
    if type_fuel_cell == "ZSW-GenStack" or type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" or \
            type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" or type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" or \
            type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" or type_fuel_cell == "ZSW-GenStack_T_62" or \
            type_fuel_cell == "ZSW-GenStack_T_76" or type_fuel_cell == "ZSW-GenStack_T_84":
        # Global
        Aact = 283.87e-4  # m². It is the MEA active area.
        nb_cell = 26  # . It is the number of cell in the stack.
        # Catalyst layer
        Hacl = 9.76346510970484e-6  # m. It is the thickness of the anode catalyst layer.
        Hccl = 5.307345509355617e-6  # m. It is the thickness of the cathode catalyst layer.
        epsilon_mc = 0.48352845256040067  # It is the volume fraction of ionomer in the CL.
        # Membrane
        Hmem = 10.684275696739961e-06  # m. It is the thickness of the membrane.
        # Gas diffusion layer
        Hgdl = 103.92479846707957e-6  # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl = 0.8184931467542469  # It is the anode/cathode GDL porosity.
        epsilon_cl = 0.5  # It is the porosity of the catalyst layer, without units.
        epsilon_c = 0.2  # It is the compression ratio of the GDL.
        #   Microporous layer
        Hmpl = 43.98306893354156e-6  # m. It is the thickness of the microporous layer.
        epsilon_mpl = 0.425  # It is the porosity of the microporous layer.
        # Gas channel
        Hagc = 230e-6  # m. It is the thickness of the anode gas channel.
        Hcgc = 300e-6  # m. It is the thickness of the cathode gas channel.
        Wagc = 430e-6  # m. It is the width of the anode gas channel.
        Wcgc = 532e-6  # m. It is the width of the cathode gas channel.
        Lgc = 246.2e-3  # m. It is the length of one channel in the bipolar plate.
        nb_channel_in_gc = 105 # . It is the number of channels in the bipolar plate.
        Ldist = 71.1e-3 # m. It is the length of the distributor, which is the volume between the gas channel and the manifold.
        #   Auxiliaries
        Lm = 25.8e-3  # m. It is the length of the manifold.
        A_T_a = 9.01e-4  # m². It is the inlet/exhaust anode manifold throttle area
        A_T_c = 22.61e-4  # m². It is the inlet/exhaust cathode manifold throttle area
        Vasm, Vcsm = Lm * A_T_a, Lm * A_T_c  # m3. It is the supply manifold volume.
        Vaem, Vcem = Vasm, Vcsm  # m-3. It is the exhaust manifold volume.
        # Interaction parameters between water and PEMFC structure
        e = 4.0  # It is the capillary exponent
        # Voltage polarization
        Re = 6.635540741264315e-8  # ohm.m². It is the electron conduction resistance of the circuit.
        i0_d_c_ref = 3.1021792383083717  # A.m-2. It is the dry reference exchange current density at the cathode.
        i0_h_c_ref = 1.0  # A.m-2. It is the fully humidified reference exchange current density at the cathode.
        kappa_co = 3.6088595047482706 # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 0.6404199981185146  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.05553, 0.10514, 0.63654  # It is the limit liquid saturation coefficients.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.

    # For EH-31 fuel cell
    elif type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
        # Global
        Aact = 85e-4  # m². It is the active area of the catalyst layer.
        nb_cell = 1  # . It is the number of cell in the stack.
        # Catalyst layer
        Hacl = 8.593e-6  # m. It is the thickness of the anode catalyst layer.
        Hccl = Hacl  # m. It is the thickness of the cathode catalyst layer.
        epsilon_mc = 0.3986  # It is the volume fraction of ionomer in the CL.
        # Membrane
        Hmem = 16.06e-6  # m. It is the thickness of the membrane.
        # Gas diffusion layer
        Hgdl = 200e-6  # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl = 0.5002  # It is the anode/cathode GDL porosity.
        epsilon_cl = 0.25  # It is the porosity of the catalyst layer, without units.
        epsilon_c = 0.2  # It is the compression ratio of the GDL.
        #   Microporous layer
        Hmpl = 30e-6  # m. It is the thickness of the microporous layer.
        epsilon_mpl = 0.4  # It is the porosity of the microporous layer.
        # Gas channel
        Hagc = 500e-6  # m. It is the thickness of the anode gas channel.
        Hcgc = Hagc  # m. It is the thickness of the cathode gas channel.
        Wagc = 450e-6  # m. It is the width of the anode gas channel.
        Wcgc = Wagc  # m. It is the width of the cathode gas channel.
        Lgc = 144e-3  # m. It is the length of one channel in the bipolar plate.
        nb_channel_in_gc = 67  # . It is the number of channels in the bipolar plate.
        Ldist = 5e-2  # m. It is the estimated length of the distributor, which is the volume between the gas channel and the manifold.
        #   Auxiliaries
        Lm = 2.03  # m. It is the length of the manifold.
        A_T_a = 11.8e-4  # m². It is the inlet/exhaust anode manifold throttle area
        A_T_c = 34.4e-4  # m². It is the inlet/exhaust cathode manifold throttle area
        Vasm, Vcsm = Lm * A_T_a, Lm * A_T_c  # m3. It is the supply manifold volume.
        Vaem, Vcem = Vasm, Vcsm  # m-3. It is the exhaust manifold volume.
        # Interaction parameters between water and PEMFC structure
        e = 4.0  # It is the capillary exponent
        # Voltage polarization
        Re = 1e-06  # ohm.m². It is the electron conduction resistance of the circuit.
        i0_d_c_ref = 14.43  # A.m-2. It is the dry reference exchange current density at the cathode.
        i0_h_c_ref = 1.0  # A.m-2. It is the fully humidified reference exchange current density at the cathode.
        kappa_co = 30.42  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 0.4152  # It is the overpotential correction exponent.
        a_slim, b_slim, a_switch = 0.05553, 0.10514, 0.82  # It is the limit liquid saturation coefficients.
        C_scl = 20e6  # F.m-3. It is the volumetric space-charge layer capacitance.

    # For other fuel cells
    else:
        raise ValueError('the type_input given is not valid.')

    return (Hacl, Hccl, epsilon_mc, Hmem, Hgdl, epsilon_gdl, epsilon_cl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc,
            Wcgc, Lgc, nb_channel_in_gc, Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem, Aact, nb_cell, e, Re, i0_d_c_ref,
            i0_h_c_ref, kappa_co, kappa_c, a_slim, b_slim, a_switch, C_scl)


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
    t0_EIS = 120*60  # s. It is the simulation starting time. [0, t0_EIS] is used to let the stack equilibrate to i_EIS.
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
