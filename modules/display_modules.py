# -*- coding: utf-8 -*-

"""This module is used to accurately plot the figures.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatter, FormatStrFormatter
from numpy.fft import fft, fftfreq
from scipy.interpolate import interp1d

# Importing constants' value and functions
from configuration.settings import F, R, E0, Pref
from modules.transitory_functions import average, Psat, C_v_sat, k_H2, k_O2
from calibration.experimental_values import (pola_exp_values, plot_experimental_polarisation_curve,
                                             pola_exp_values_calibration)

# General edition
colors = mpl.colormaps['tab10']


# __________________________________________________Polarisation curve__________________________________________________

def plot_polarisation_curve(variables, operating_inputs, parameters, ax, show=True):
    """
    This function plots the model polarisation curve, and compare it to the experimental one (if it exists). The
    polarisation curve is a classical representation of the cell performances, showing the cell voltage as a function
    of the current density.
    To generate it, the current density is increased step by step, and the cell voltage is recorded at each step.
    The time for which this point is captured is determined using the following approach: at the beginning of each load,
    a delta_t_load_pola time is needed to raise the current density to its next value. Subsequently, a delta_t_break_pola
    time is observed to ensure the dynamic stability of the stack's variables before initiating a new load. Finally,
    each polarisation point is recorded at the end of each delta_t_break_pola time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the polarisation curve will be plotted.
    show : bool, optional
        If True, the polarisation curve will be displayed. If False, it will not be displayed.
    """

    # Extraction of the variables
    t, Ucell_t = np.array(variables['t']), np.array(variables['Ucell'])
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    pola_current_parameters = parameters['pola_current_parameters']
    delta_t_ini_pola = pola_current_parameters['delta_t_ini_pola']
    delta_t_load_pola = pola_current_parameters['delta_t_load_pola']
    delta_t_break_pola = pola_current_parameters['delta_t_break_pola']
    delta_i_pola, i_max_pola = pola_current_parameters['delta_i_pola'], pola_current_parameters['i_max_pola']
    type_fuel_cell, type_current = parameters['type_fuel_cell'], parameters['type_current']
    type_auxiliary, type_control = parameters['type_auxiliary'], parameters['type_control']
    type_plot = parameters['type_plot']
    # Extraction of the experimental current density and voltage values.
    i_exp_t, U_exp_t = pola_exp_values(type_fuel_cell)  # (A.m-2, V).

    if type_plot == "fixed":
        # Creation of ifc_t
        n = len(t)
        ifc_t = np.zeros(n)
        for i in range(n):
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²

        # Recovery of ifc and Ucell from the model after each stack stabilisation
        nb_loads = int(i_max_pola / delta_i_pola)  # Number of loads which are made
        ifc_discretized = np.zeros(nb_loads + 1) # One point is taken at ifc = 0, before the first load.
        Ucell_discretized = np.zeros(nb_loads + 1) # One point is taken at ifc = 0, before the first load.
        for i in range(nb_loads+1):
            t_load = delta_t_ini_pola + i * (delta_t_load_pola + delta_t_break_pola) # time for measurement
            idx = (np.abs(t - t_load)).argmin()  # the corresponding index
            ifc_discretized[i] = ifc_t[idx]  # the last value at the end of each load
            Ucell_discretized[i] = Ucell_t[idx]  # the last value at the end of each load

        # Plot the experimental polarization curve and calculate the simulation error compared with experimental data
        if type_fuel_cell != "manual_setup" and \
           type_auxiliary == "forced-convective_cathode_with_flow-through_anode":  # Experimental points are accessible
            # Plot of the experimental polarization curve
            i_exp_t = i_exp_t / 1e4  # Conversion in A/cm²
            plot_experimental_polarisation_curve(type_fuel_cell, i_exp_t, U_exp_t, ax)
            # Calculate the simulation error compared with experimental data
            #       Experimental points are interpolated to correspond to the model points
            Ucell_interpolated = interp1d(ifc_discretized, Ucell_discretized, kind='linear')(i_exp_t)
            sim_error = calculate_simulation_error(Ucell_interpolated, U_exp_t)
        else:
            sim_error = None

        # Plot the model polarisation curve
        plot_specific_line(ifc_discretized, Ucell_discretized, type_fuel_cell, type_current, type_auxiliary,
                           type_control, sim_error, ax)
        plot_pola_instructions(type_fuel_cell, ax, show)

    else:  # type_plot == "dynamic"
        # Plot of the polarisation curve produced by the model
        idx = (np.abs(t - t[-1])).argmin()  # index for polarisation measurement
        ifc = np.array(current_density(t[idx], parameters) / 1e4)  # time for polarisation measurement
        Ucell = np.array(Ucell_t[idx])  # voltage measurement
        ax.plot(ifc, Ucell, 'og', markersize=2)

    # Add the common instructions for the plot
    ax.set_xlabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                  labelpad=3)
    ax.set_ylabel(r'$\mathbf{Cell}$ $\mathbf{voltage}$ $\mathbf{U_{cell}}$ $\mathbf{\left( V \right)}$', labelpad=3)
    if type_plot == "fixed":
        ax.legend(loc='best')

def plot_polarisation_curve_for_cali(variables, operating_inputs, parameters, ax):
    """
    This function plots the model polarisation curve, and compare it to the experimental one. The
    polarisation curve is a classical representation of the cell performances, showing the cell voltage as a function
    of the current density.
    To generate it, the current density is increased step by step, and the cell voltage is recorded at each step.
    The time for which this point is captured is determined using the following approach: at the beginning of each load,
    a delta_t_load_pola time is needed to raise the current density to its next value. Subsequently, a delta_t_break_pola
    time is observed to ensure the dynamic stability of the stack's variables before initiating a new load. Finally,
    each polarisation point is recorded at the end of each delta_t_break_pola time.


    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the polarisation curve will be plotted.
    """

    # Extraction of the variables
    t, Ucell_t = np.array(variables['t']), np.array(variables['Ucell'])
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    pola_current_for_cali_parameters = parameters['pola_current_for_cali_parameters']
    delta_t_ini_pola_cali = pola_current_for_cali_parameters['delta_t_ini_pola_cali']
    delta_t_load_pola_cali = pola_current_for_cali_parameters['delta_t_load_pola_cali']
    delta_t_break_pola_cali = pola_current_for_cali_parameters['delta_t_break_pola_cali']
    type_fuel_cell, type_current = parameters['type_fuel_cell'], parameters['type_current']
    type_auxiliary, type_control = parameters['type_auxiliary'], parameters['type_control']
    type_plot = parameters['type_plot']
    # Extraction of the experimental current density and voltage values for the calibration.
    i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(parameters['type_fuel_cell'])  # (A.m-2, V).

    # Creation of ifc_t
    n = len(t)
    ifc_t = np.zeros(n)
    for i in range(n):
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²

    # Recovery of ifc and Ucell from the model after each stack stabilisation
    nb_loads = len(i_exp_cali_t)  # Number of loads which are made
    delta_t_cali = delta_t_load_pola_cali + delta_t_break_pola_cali  # s. It is the time of one load.
    ifc_discretized = np.zeros(nb_loads)
    Ucell_discretized = np.zeros(nb_loads)
    for i in range(nb_loads):
        t_load = delta_t_ini_pola_cali + (i + 1) * delta_t_cali # time for measurement
        idx = (np.abs(t - t_load)).argmin()  # the corresponding index
        ifc_discretized[i] = ifc_t[idx]  # the last value at the end of each load
        Ucell_discretized[i] = Ucell_t[idx]  # the last value at the end of each load

    # Plot the experimental polarization curve
    i_exp_cali_t = i_exp_cali_t / 1e4  # Conversion in A/cm²
    plot_experimental_polarisation_curve(type_fuel_cell, i_exp_cali_t, U_exp_cali_t, ax)

    # Plot the model polarisation curve
    sim_error = calculate_simulation_error(Ucell_discretized, U_exp_cali_t) # Calculate the simulation error
    plot_specific_line(ifc_discretized, Ucell_discretized, type_fuel_cell, type_current, type_auxiliary,
                       type_control, sim_error, ax)
    plot_pola_instructions(type_fuel_cell, ax)

    # Add the common instructions for the plot
    ax.set_xlabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                  labelpad=3)
    ax.set_ylabel(r'$\mathbf{Cell}$ $\mathbf{voltage}$ $\mathbf{U_{cell}}$ $\mathbf{\left( V \right)}$', labelpad=3)
    if type_plot == "fixed":
        ax.legend(loc='best')


# _______________________________________________________EIS curves_____________________________________________________

def make_Fourier_transformation(variables, operating_inputs, parameters):
    """
    This function calculates the Fourier transformation of both cell voltage and current density. It will be used to
    display the Nyquist and Bode diagrams.
    To generate it at each frequency change, the cell voltage and the current density are recorded. The time for which
    these points are captured is determined using the following approach: at the beginning of each frequency change, a
    delta_t_break_EIS time is observed to ensure the dynamic stability of the stack's variables. Subsequently, a
    delta_t_measurement_EIS time is needed to record the cell voltage and the current density.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    dict
        Dictionary containing the Fourier transformation (FT) of the cell voltage and the current density, all amplitude
        values of the cell voltage calculated by the FT, the amplitude of the cell voltage at the frequency of the
        perturbation, all frequency values used vy the FT, the frequency of the perturbation, and the number of points
        used in the FT.
    """

    # Extraction of the variables
    t, Ucell_t = np.array(variables['t']), np.array(variables['Ucell'])
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    t_EIS = parameters['t_EIS']

    # Creation of ifc
    ifc_t = np.zeros(len(t))
    for i in range(len(t)):
        ifc_t[i] = current_density(t[i], parameters)

    # Identify the areas where Ucell and ifc can be measured for the EIS: after equilibrium and at each frequency change
    t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
    n_inf = np.where(t_new_start_EIS <= t[0])[0][-1]  # The number of frequency changes which has been mad so far.
    Ucell_EIS_measured = Ucell_t[np.where((t > (t[0] + delta_t_break_EIS[n_inf])) &
                                          (t < (t[0] + delta_t_break_EIS[n_inf] + delta_t_measurement_EIS[n_inf])))]
    ifc_EIS_measured = ifc_t[np.where((t > (t[0] + delta_t_break_EIS[n_inf])) &
                                      (t < (t[0] + delta_t_break_EIS[n_inf] + delta_t_measurement_EIS[n_inf])))]

    # Determination of the Fourier transformation
    N = Ucell_EIS_measured.size  # Number of points used for the Fourier transformation
    Ucell_Fourier = fft(Ucell_EIS_measured)  # Ucell Fourier transformation
    ifc_Fourier = fft(ifc_EIS_measured)  # ifc Fourier transformation
    A_period_t = np.concatenate(
        ([np.abs(Ucell_Fourier)[0] / N], np.abs(Ucell_Fourier[1:N // 2]) * 2 / N))  # Recovery of
    #                                                                             all amplitude values calculated by fft
    A = max(A_period_t[1:])  # Amplitude at the frequency of the perturbation
    freq_t = fftfreq(N)[:N // 2]  # Recovery of all frequency values used by fft
    f = freq_t[np.argmax(A_period_t == A)]  # Recovery of the studied frequency

    return {'Ucell_Fourier': Ucell_Fourier, 'ifc_Fourier': ifc_Fourier, 'A_period_t': A_period_t, 'A': A,
            'freq_t': freq_t, 'f': f, 'N': N}


def plot_EIS_curve_Nyquist(parameters, Fourier_results, ax):
    """
    This function is used to plot the Nyquist diagram of the EIS curves.

    Parameters
    ----------
    parameters : dict
        Parameters of the fuel cell model.
    Fourier_results : dict
        Dictionary containing the Fourier transformation (FT) of the cell voltage and the current density, all amplitude
        values of the cell voltage calculated by the FT, the amplitude of the cell voltage at the frequency of the
        perturbation, all frequency values used vy the FT, the frequency of the perturbation, and the number of points
        used in the FT.
    ax : matplotlib.axes.Axes
        Axes on which the Nyquist diagram will be plotted.
    """

    # Extraction of the parameters
    i_EIS, ratio_EIS, type_fuel_cell = parameters['i_EIS'], parameters['ratio_EIS'], parameters['type_fuel_cell']
    # Extraction of the Fourier results
    Ucell_Fourier, ifc_Fourier = Fourier_results['Ucell_Fourier'], Fourier_results['ifc_Fourier']
    f_Fourier = Fourier_results['f']
    A_period_t, A, N = Fourier_results['A_period_t'], Fourier_results['A'], Fourier_results['N']

    # Calculation of the real and imaginary component of the impedance for each period
    Z0 = A / (ratio_EIS * (-i_EIS)) * 1e7  # Impedance of the perturbation in mΩ.cm². The sign of i is inverted to
    #                  comply with the standards of EIS, which measure a device under load rather than a current source.
    theta_U_t = np.angle(Ucell_Fourier[0:N // 2])  # Recovery of all dephasing values calculated by fft
    theta_i_t = np.angle(ifc_Fourier[0:N // 2])  # Recovery of all dephasing values calculated by fft
    theta_U = theta_U_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation
    theta_i = theta_i_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation
    Z_real = Z0 * np.cos(theta_U - theta_i)  # Real component of the impedance for each period
    Z_imag = Z0 * np.sin(theta_U - theta_i)  # Imaginary component of the impedance for each period

    # Plot the Nyquist diagram
    ax.plot(Z_real, -Z_imag, 'o', color=colors(0), label='Nyquist diagram')
    ax.set_xlabel(r'$\mathbf{Z_{real}}$ $\mathbf{(m\Omega.cm^{2})}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{-Z_{imag}}$ $\mathbf{(m\Omega.cm^{2})}$', labelpad=3)
    #       Plot instructions
    plot_EIS_Nyquist_instructions(type_fuel_cell, f_Fourier, Z_real, -Z_imag, ax)


def plot_EIS_curve_Bode_amplitude(parameters, Fourier_results, ax):
    """This function is used to plot the amplitude Bode diagram of the EIS curves.

    Parameters
    ----------
    parameters : dict
        Parameters of the fuel cell model.
    Fourier_results : dict
        Dictionary containing the Fourier transformation (FT) of the cell voltage and the current density, all amplitude
        values of the cell voltage calculated by the FT, the amplitude of the cell voltage at the frequency of the
        perturbation, all frequency values used vy the FT, the frequency of the perturbation, and the number of points
        used in the FT.
    ax : matplotlib.axes.Axes
        Axes on which the amplitude Bode diagram will be plotted.

    """

    # Extraction of the parameters
    i_EIS, ratio_EIS, f_EIS = parameters['i_EIS'], parameters['ratio_EIS'], parameters['f_EIS']
    type_fuel_cell = parameters['type_fuel_cell']
    # Extraction of the Fourier results
    A, f = Fourier_results['A'], Fourier_results['f']

    # Calculation of the impedance of the perturbation
    Z0 = A / (ratio_EIS * (-i_EIS)) * 1e7  # in mΩ.cm². The sign of i is inverted to comply with the standards of EIS,
    #                                        which measure a device under load rather than a current source.

    # Plot the amplitude Bode diagram
    ax.plot(f, np.abs(Z0), 'o', color=colors(1), label='Amplitude Bode diagram')
    ax.set_xlabel(r'$\mathbf{Frequency}$ $\mathbf{(Hz,}$ $\mathbf{logarithmic}$ $\mathbf{scale)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Impedance}$ $\mathbf{amplitude}$ $\mathbf{(m\Omega.cm^{2})}$', labelpad=3)
    #   Plot instructions
    plot_Bode_amplitude_instructions(f_EIS, type_fuel_cell, ax)


def plot_EIS_curve_Bode_angle(parameters, Fourier_results, ax):
    """This function is used to plot the angle Bode diagram. It only works with an entry signal made with a cosinus
    (not a sinus).

    Parameters
    ----------
    Fourier_results : dict
        Dictionary containing the Fourier transformation (FT) of the cell voltage and the current density, all amplitude
        values of the cell voltage calculated by the FT, the amplitude of the cell voltage at the frequency of the
        perturbation, all frequency values used vy the FT, the frequency of the perturbation, and the number of points
        used in the FT.
    ax : matplotlib.axes.Axes
        Axes on which the angle Bode diagram will be plotted.
    """

    # Extraction of the parameters
    f_EIS, type_fuel_cell = parameters['f_EIS'], parameters['type_fuel_cell']
    # Extraction of the Fourier results
    Ucell_Fourier, ifc_Fourier = Fourier_results['Ucell_Fourier'], Fourier_results['ifc_Fourier']
    A_period_t, A = Fourier_results['A_period_t'], Fourier_results['A']
    f, N = Fourier_results['f'], Fourier_results['N']

    # Calculation of the dephasing values at the frequency of the perturbation
    theta_U_t = np.angle(Ucell_Fourier[0:N // 2])  # Recovery of all dephasing values calculated by fft
    theta_i_t = np.angle(ifc_Fourier[0:N // 2]) + np.pi  # Recovery of all dephasing values calculated by fft.
    #                                                    An angle of pi is added to comply with the standards of EIS,
    #                                                    which measure a device under load rather than a current source.
    theta_U = theta_U_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation
    theta_i = theta_i_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation
    phi_U_i = ((theta_U - theta_i) * 180 / np.pi) % 360 # Dephasing between Ucell and ifc with a value between 0 and 360
    if phi_U_i > 180:
        phi_U_i -= 360 # To have a value between -180 and 180

    # Plot the angle Bode diagram
    ax.plot(f, phi_U_i, 'o', color=colors(2), label='Angle Bode diagram')
    ax.set_xlabel(r'$\mathbf{Frequency}$ $\mathbf{(Hz,}$ $\mathbf{logarithmic}$ $\mathbf{scale)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Phase}$ $\mathbf{(^\circ)}$', labelpad=3)
    #   Plot instructions
    plot_Bode_phase_instructions(f_EIS, type_fuel_cell, ax)


def plot_EIS_curve_tests(variables, operating_inputs, parameters, Fourier_results):
    """This function is used to test the accuracy of the EIS results. It compares the reconstructed Ucell_Fourier(t)
    from the Fourier transformation with the current density ifc(t), and displays Ucell(t) given by the model with the
    reconstructed Ucell_Fourier(t).

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    Fourier_results : dict
        Dictionary containing the Fourier transformation (FT) of the cell voltage and the current density, all amplitude
        values of the cell voltage calculated by the FT, the amplitude of the cell voltage at the frequency of the
        perturbation, all frequency values used vy the FT, the frequency of the perturbation, and the number of points
        used in the FT.
    """

    # Extraction of the variables
    t, Ucell_t = np.array(variables['t']), variables['Ucell']
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    i_EIS, ratio_EIS = parameters['i_EIS'], parameters['ratio_EIS']
    t_EIS, f_EIS = parameters['t_EIS'], parameters['f_EIS']
    # Extraction of the Fourier results
    Ucell_Fourier, ifc_Fourier = Fourier_results['Ucell_Fourier'], Fourier_results['ifc_Fourier']
    A_period_t, A = Fourier_results['A_period_t'], Fourier_results['A']
    f, N = Fourier_results['f'], Fourier_results['N']

    # Reconstructed Ucell with a cosinus form, and comparison of its form with the current density one.
    t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS
    n_inf = np.where(t_new_start_EIS <= t[0])[0][-1]  # The number of frequency changes which has been made.
    f_current = np.logspace(f_power_min_EIS, f_power_max_EIS, num=nb_f_EIS)
    theta_U_t = np.angle(Ucell_Fourier[0:N // 2])  # Recovery of all dephasing values calculated by fft
    theta_i_t = np.angle(ifc_Fourier[0:N // 2])  # Recovery of all dephasing values calculated by fft
    theta_U = theta_U_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation
    theta_i = theta_i_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation
    print("Ucell:", round(A_period_t[0], 4), ' + ', round(A, 6), " * np.cos(2*np.pi*", round(f, 4), "*t + ",
          round(theta_U, 4), "). ")
    print("Current:", i_EIS, ' + ', ratio_EIS * i_EIS, " * np.cos(2*np.pi*", round(f_current[n_inf], 4), "*t + ",
          round(theta_i, 4), "). \n")

    # Display ifc(t)
    plt.figure(3)
    plt.subplot(2, 1, 1)
    #   Creation of ifc_t
    n = len(t)
    ifc_t = np.zeros(n)
    for i in range(n):  # Conversion in A/cm²
        ifc_t[i] = current_density(t[i], parameters) / 1e4
    #   Plot of ifc_t
    plt.plot(t, ifc_t, color='blue', label='ifc')
    plt.xlabel('Time (s)')
    plt.ylabel('Current density (A/cm²)')
    plt.title('The current density\nbehaviour over time')

    # Display Ucell(t) and compare it with the reconstructed Ucell_Fourier(t) from the Fourier transformation
    plt.subplot(2, 1, 2)
    Ucell_Fourier = A_period_t[0] + A * np.cos(2 * np.pi * f * t + theta_U)
    plt.plot(t, Ucell_t, color='blue', label='Ucell')
    plt.plot(t, Ucell_Fourier, color='black', label='Ucell_Fourier')
    plt.xlabel('Time (s)')
    plt.ylabel('Cell voltage (V)')
    plt.title('The cell voltage\nbehaviour over time')

# ___________________________________________________Internal variables_________________________________________________

def plot_ifc(variables, operating_inputs, parameters, ax):
    """This function plots the current density as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    n : int
        Number of points used to plot the current density.
    ax : matplotlib.axes.Axes
        Axes on which the current density will be plotted.
    """

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]

    # Plot the current density: ifc
    n = len(t)
    ifc_t = np.zeros(n)
    for i in range(n):  # Creation of ifc_t
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
    ax.plot(t, ifc_t, color=colors(0), label=r'$\mathregular{i_{fc}}$')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                  labelpad=3)
    ax.legend([r'$\mathregular{i_{fc}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)


def plot_J(variables, parameters, ax):
    """This function plots the sorption and dissolved water flows as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the flows will be plotted.
    """
    # Extraction of the operating inputs and the parameters
    Hacl, Hccl = parameters['Hacl'], parameters['Hccl']
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    S_abs_acl_t = np.array(variables['S_abs_acl'])[mask]
    S_abs_ccl_t = np.array(variables['S_abs_ccl'])[mask]
    J_lambda_acl_mem_t = np.array(variables['J_lambda_acl_mem'])[mask]
    J_lambda_mem_ccl_t = np.array(variables['J_lambda_mem_ccl'])[mask]

    # Plot the sorption and dissolved water flows: J
    J_abs_acl, J_abs_ccl = S_abs_acl_t * Hacl, S_abs_ccl_t * Hccl  # Conversion in mol.m⁻².s⁻¹ for comparison
    ax.plot(t, J_abs_acl, color=colors(2))
    ax.plot(t, J_lambda_acl_mem_t, color=colors(3))
    ax.plot(t, J_abs_ccl, color=colors(4))
    ax.plot(t, J_lambda_mem_ccl_t, color=colors(7))
    ax.legend([r'$\mathregular{J_{abs,acl}}$', r'$\mathregular{J_{\lambda,mem,acl}}$', r'$\mathregular{J_{abs,ccl}}$',
               r'$\mathregular{J_{\lambda,mem,ccl}}$'], loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Flows}$ $\mathbf{J}$ $\mathbf{\left( mol.m^{-2}.s^{-1} \right)}$', labelpad=3)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # Plot instructions
    plot_general_instructions(ax)


def plot_C_v(variables, parameters, ax):
    """This function plots the vapor concentrations at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the vapor concentration will be plotted.
    """

    # Extraction of the parameter
    n_gdl, n_mpl, type_current, type_plot = parameters['n_gdl'], parameters['n_mpl'], parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    C_v_agc_t = np.array(variables['C_v_agc'])[mask]
    C_v_agdl_t = np.array(variables[f'C_v_agdl_{int(np.ceil(n_gdl / 2))}'])[mask]
    C_v_ampl_t = np.array(variables[f'C_v_ampl_{int(np.ceil(n_mpl / 2))}'])[mask]
    C_v_acl_t = np.array(variables['C_v_acl'])[mask]
    C_v_ccl_t = np.array(variables['C_v_ccl'])[mask]
    C_v_cmpl_t = np.array(variables[f'C_v_cmpl_{int(np.ceil(n_mpl / 2))}'])[mask]
    C_v_cgdl_t = np.array(variables[f'C_v_cgdl_{int(np.ceil(n_gdl / 2))}'])[mask]
    C_v_cgc_t = np.array(variables['C_v_cgc'])[mask]
    T_ccl = np.array(variables['T_ccl'])[mask]

    # Plot the vapor concentrations at different spatial localisations Cv
    C_v_sat_ccl_t = np.array([C_v_sat(T) for T in T_ccl])
    ax.plot(t, C_v_agc_t, color=colors(0))
    ax.plot(t, C_v_agdl_t, color=colors(1))
    ax.plot(t, C_v_ampl_t, color=colors(2))
    ax.plot(t, C_v_acl_t, color=colors(3))
    ax.plot(t, C_v_ccl_t, color=colors(5))
    ax.plot(t, C_v_cmpl_t, color=colors(6))
    ax.plot(t, C_v_cgdl_t, color=colors(7))
    ax.plot(t, C_v_cgc_t, color=colors(8))
    ax.plot(t, C_v_sat_ccl_t, color='k')
    ax.legend([r'$\mathregular{C_{v,agc}}$', r'$\mathregular{C_{v,agdl}}$', r'$\mathregular{C_{v,ampl}}$',
               r'$\mathregular{C_{v,acl}}$', r'$\mathregular{C_{v,ccl}}$', r'$\mathregular{C_{v,cmpl}}$',
               r'$\mathregular{C_{v,cgdl}}$', r'$\mathregular{C_{v,cgc}}$', r'$\mathregular{C_{v,sat,ccl}}$'], loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r"$\mathbf{Vapor}$ $\mathbf{concentration}$ $\mathbf{C_{v}}$ $\mathbf{\left( mol.m^{-3} \right)}$",
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_lambda(variables, operating_inputs, parameters, ax):
    """This function plots the water content at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the water content will be plotted.
    """

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    pola_current_parameters, type_current = parameters['pola_current_parameters'], parameters['type_current']
    type_plot = parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    lambda_acl_t = np.array(variables['lambda_acl'])[mask]
    lambda_mem_t = np.array(variables['lambda_mem'])[mask]
    lambda_ccl_t = np.array(variables['lambda_ccl'])[mask]

    # Plot the water content at different spatial localisations: lambda
    if type_current == "polarization":
        n = len(t)
        ifc_t = np.zeros(n)
        for i in range(n):  # Creation of i_fc
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        # Recovery of the internal states from the model after each stack stabilisation
        delta_t_ini_pola = pola_current_parameters['delta_t_ini_pola']
        delta_t_load_pola = pola_current_parameters['delta_t_load_pola']
        delta_t_break_pola = pola_current_parameters['delta_t_break_pola']
        nb_loads = int(pola_current_parameters['i_max_pola'] / pola_current_parameters['delta_i_pola'])  # Number of loads which are made
        ifc_discretized_t, lambda_acl_discretized_t = np.zeros(nb_loads), np.zeros(nb_loads)
        lambda_mem_discretized_t, lambda_ccl_discretized_t = np.zeros(nb_loads), np.zeros(nb_loads)
        for i in range(nb_loads):
            t_load = delta_t_ini_pola + (i + 1) * (delta_t_load_pola + delta_t_break_pola)  # time for measurement
            idx = (np.abs(t - t_load)).argmin()  # the corresponding index
            ifc_discretized_t[i] = ifc_t[idx]  # the last value at the end of each load
            lambda_acl_discretized_t[i] = lambda_acl_t[idx]  # the last value at the end of each load
            lambda_mem_discretized_t[i] = lambda_mem_t[idx]  # the last value at the end of each load
            lambda_ccl_discretized_t[i] = lambda_ccl_t[idx]  # the last value at the end of each load
        ax.scatter(ifc_discretized_t, lambda_acl_discretized_t, marker='o', color=colors(2))
        ax.scatter(ifc_discretized_t, lambda_mem_discretized_t, marker='o', color=colors(3))
        ax.scatter(ifc_discretized_t, lambda_ccl_discretized_t, marker='o', color=colors(4))
        ax.set_xlabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                      labelpad=3)
    else:
        ax.plot(t, lambda_acl_t, color=colors(3))
        ax.plot(t, lambda_mem_t, color=colors(4))
        ax.plot(t, lambda_ccl_t, color=colors(5))
        ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Water}$ $\mathbf{content}$ $\mathbf{\lambda}$', labelpad=3)
    ax.legend([r'$\mathregular{\lambda_{acl}}$', r'$\mathregular{\lambda_{mem}}$',
               r'$\mathregular{\lambda_{ccl}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)


def plot_s(variables, operating_inputs, parameters, ax):
    """This function plots the liquid water saturation at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the liquid water saturation will be plotted.
    """

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    n_gdl, n_mpl, pola_current_parameters = parameters['n_gdl'], parameters['n_mpl'], parameters['pola_current_parameters']
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    s_agdl_t = np.array(variables[f's_agdl_{int(np.ceil(n_gdl / 2))}'])[mask]
    s_ampl_t = np.array(variables[f's_ampl_{int(np.ceil(n_mpl / 2))}'])[mask]
    s_acl_t = np.array(variables['s_acl'])[mask]
    s_ccl_t = np.array(variables['s_ccl'])[mask]
    s_cmpl_t = np.array(variables[f's_cmpl_{int(np.ceil(n_mpl / 2))}'])[mask]
    s_cgdl_t = np.array(variables[f's_cgdl_{int(np.ceil(n_gdl / 2))}'])[mask]

    # Plot the liquid water saturation at different spatial localisations: s
    if type_current == "polarization":
        n = len(t)
        ifc_t = np.zeros(n)
        for i in range(n):  # Creation of i_fc
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        # Recovery of the internal states from the model after each stack stabilisation
        delta_t_ini_pola = pola_current_parameters['delta_t_ini_pola']
        delta_t_load_pola = pola_current_parameters['delta_t_load_pola']
        delta_t_break_pola = pola_current_parameters['delta_t_break_pola']
        nb_loads = int(pola_current_parameters['i_max_pola'] / pola_current_parameters['delta_i_pola']) # Number of loads
        ifc_discretized_t = np.zeros(nb_loads)
        s_agdl_discretized_t, s_ampl_discretized_t, s_acl_discretized_t = [np.zeros(nb_loads)] * 3
        s_ccl_discretized_t, s_cmpl_discretized_t, s_cgdl_discretized_t = [np.zeros(nb_loads)] * 3
        for i in range(nb_loads):
            t_load = delta_t_ini_pola + (i + 1) * (delta_t_load_pola + delta_t_break_pola)  # time for measurement
            idx = (np.abs(t - t_load)).argmin()  # the corresponding index
            ifc_discretized_t[i] = ifc_t[idx]  # the last value at the end of each load
            s_agdl_discretized_t[i] = s_agdl_t[idx]  # the last value at the end of each load
            s_ampl_discretized_t[i] = s_ampl_t[idx]  # the last value at the end of each load
            s_acl_discretized_t[i] = s_acl_t[idx]  # the last value at the end of each load
            s_ccl_discretized_t[i] = s_ccl_t[idx]  # the last value at the end of each load
            s_cmpl_discretized_t[i] = s_cmpl_t[idx]  # the last value at the end of each load
            s_cgdl_discretized_t[i] =s_cgdl_t[idx]  # the last value at the end of each load
        ax.scatter(ifc_discretized_t, s_agdl_discretized_t, marker='o', color=colors(1))
        ax.scatter(ifc_discretized_t, s_ampl_discretized_t, marker='o', color=colors(2))
        ax.scatter(ifc_discretized_t, s_acl_discretized_t, marker='o', color=colors(3))
        ax.scatter(ifc_discretized_t, s_ccl_discretized_t, marker='o', color=colors(5))
        ax.scatter(ifc_discretized_t, s_cmpl_discretized_t, marker='o', color=colors(6))
        ax.scatter(ifc_discretized_t, s_cgdl_discretized_t, marker='o', color=colors(7))
        ax.set_xlabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                      labelpad=3)
    else:
        ax.plot(t, s_agdl_t, color=colors(1))
        ax.plot(t, s_ampl_t, color=colors(2))
        ax.plot(t, s_acl_t, color=colors(3))
        ax.plot(t, s_ccl_t, color=colors(5))
        ax.plot(t, s_cmpl_t, color=colors(6))
        ax.plot(t, s_cgdl_t, color=colors(7))
        ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Liquid}$ $\mathbf{water}$ $\mathbf{saturation}$ $\mathbf{s}$', labelpad=3)
    ax.legend([r'$\mathregular{s_{agdl}}$', r'$\mathregular{s_{ampl}}$', r'$\mathregular{s_{acl}}$',
               r'$\mathregular{s_{ccl}}$', r'$\mathregular{s_{cmpl}}$', r'$\mathregular{s_{cgdl}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)


def plot_C_H2(variables, parameters, ax):
    """This function plots the hydrogen concentration at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the hydrogen concentration will be plotted.
    """

    # Extraction of the parameters
    n_gdl, n_mpl, type_current, type_plot = parameters['n_gdl'], parameters['n_mpl'], parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    C_H2_agc_t = np.array(variables['C_H2_agc'])[mask]
    C_H2_agdl_t = np.array(variables[f'C_H2_agdl_{int(np.ceil(n_gdl / 2))}'])[mask]
    C_H2_ampl_t = np.array(variables[f'C_H2_ampl_{int(np.ceil(n_mpl / 2))}'])[mask]
    C_H2_acl_t = np.array(variables['C_H2_acl'])[mask]

    # Plot the hydrogen concentration at different spatial localisations: C_H2
    ax.plot(t, C_H2_agc_t, color=colors(0))
    ax.plot(t, C_H2_agdl_t, color=colors(1))
    ax.plot(t, C_H2_ampl_t, color=colors(2))
    ax.plot(t, C_H2_acl_t, color=colors(3))
    ax.legend([r'$\mathregular{C_{H_{2},agc}}$', r'$\mathregular{C_{H_{2},agdl}}$', r'$\mathregular{C_{H_{2},ampl}}$',
               r'$\mathregular{C_{H_{2},acl}}$'], loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Hydrogen}$ $\mathbf{concentration}$ $\mathbf{C_{H_{2}}}$ $\mathbf{\left( mol.m^{-3} \right)}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_C_O2(variables, parameters, ax):
    """This function plots the oxygen concentration at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the oxygen concentration will be plotted.
    """

    # Extraction of the parameters
    n_gdl, n_mpl, type_current, type_plot = parameters['n_gdl'], parameters['n_mpl'], parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    C_O2_ccl_t = np.array(variables['C_O2_ccl'])[mask]
    C_O2_cmpl_t = np.array(variables[f'C_O2_cmpl_{int(np.ceil(n_mpl / 2))}'])[mask]
    C_O2_cgdl_t = np.array(variables[f'C_O2_cgdl_{int(np.ceil(n_gdl / 2))}'])[mask]
    C_O2_cgc_t = np.array(variables['C_O2_cgc'])[mask]

    # Plot the oxygen concentration at different spatial localisations: C_O2
    ax.plot(t, C_O2_ccl_t, color=colors(5))
    ax.plot(t, C_O2_cmpl_t, color=colors(6))
    ax.plot(t, C_O2_cgdl_t, color=colors(7))
    ax.plot(t, C_O2_cgc_t, color=colors(8))
    ax.legend([r'$\mathregular{C_{O_{2},ccl}}$', r'$\mathregular{C_{O_{2},cmpl}}$', r'$\mathregular{C_{O_{2},cgdl}}$',
               r'$\mathregular{C_{O_{2},cgc}}$'], loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Oxygen}$ $\mathbf{concentration}$ $\mathbf{C_{O_{2}}}$ $\mathbf{\left( mol.m^{-3} \right)}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_C_N2(variables, parameters, ax):
    """This function plots the nitrogen concentration as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    ax : matplotlib.axes.Axes
        Axes on which the nitrogen concentration will be plotted.
    """

    # Extraction of the parameters
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    C_N2_a_t = np.array(variables['C_N2_a'])[mask]
    C_N2_c_t = np.array(variables['C_N2_c'])[mask]

    # Plot C_N2
    ax.plot(t, C_N2_a_t, color=colors(6))
    ax.plot(t, C_N2_c_t, color=colors(6))
    ax.legend([r'$\mathregular{C_{N_{2},a}}$', r'$\mathregular{C_{N_{2},c}}$'], loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Nitrogen}$ $\mathbf{concentration}$ $\mathbf{C_{N_{2}}}$ $\mathbf{\left( mol.m^{-3} \right)}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_T(variables, operating_inputs, parameters, ax):
    """This function plots the vapor concentrations at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the vapor concentration will be plotted.
    """

    # Extraction of the operating inputs and parameters
    T_des = operating_inputs['T_des']
    n_gdl, n_mpl, type_current, type_plot = parameters['n_gdl'], parameters['n_mpl'], parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables and the operating inputs
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    T_agc_t = np.array(variables['T_agc'])[mask] - 273.15 # Conversion in °C.
    T_agdl_t = np.array(variables[f'T_agdl_{int(np.ceil(n_gdl / 2))}'])[mask] - 273.15 # Conversion in °C.
    T_ampl_t = np.array(variables[f'T_ampl_{int(np.ceil(n_mpl / 2))}'])[mask] - 273.15 # Conversion in °C.
    T_acl_t = np.array(variables['T_acl'])[mask] - 273.15  # Conversion in °C.
    T_mem_t = np.array(variables['T_mem'])[mask] - 273.15 # Conversion in °C.
    T_ccl_t = np.array(variables['T_ccl'])[mask] - 273.15 # Conversion in °C.
    T_cmpl_t = np.array(variables[f'T_cmpl_{int(np.ceil(n_mpl / 2))}'])[mask] - 273.15  # Conversion in °C.
    T_cgdl_t = np.array(variables[f'T_cgdl_{int(np.ceil(n_gdl / 2))}'])[mask] - 273.15 # Conversion in °C.
    T_cgc_t = np.array(variables['T_cgc'])[mask] - 273.15 # Conversion in °C.

    # Plot the temperature at different spatial localisations
    T_des_t = np.array([T_des - 273.15] * len(t))
    ax.plot(t, T_agc_t, color=colors(0))
    ax.plot(t, T_agdl_t, color=colors(1))
    ax.plot(t, T_ampl_t, color=colors(2))
    ax.plot(t, T_acl_t, color=colors(3))
    ax.plot(t, T_mem_t, color=colors(4))
    ax.plot(t, T_ccl_t, color=colors(5))
    ax.plot(t, T_cmpl_t, color=colors(6))
    ax.plot(t, T_cgdl_t, color=colors(7))
    ax.plot(t, T_cgc_t, color=colors(8))
    ax.plot(t, T_des_t, color='k')
    ax.legend([r'$\mathregular{T_{agc}}$', r'$\mathregular{T_{agdl}}$', r'$\mathregular{T_{ampl}}$',
               r'$\mathregular{T_{acl}}$', r'$\mathregular{T_{mem}}$', r'$\mathregular{T_{ccl}}$',
               r'$\mathregular{T_{cmpl}}$', r'$\mathregular{T_{cgdl}}$', r'$\mathregular{T_{cgc}}$',
               r'$\mathregular{T_{des}}$'], loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r"$\mathbf{Temperature}$ $\mathbf{T}$ $\mathbf{\left( °C \right)}$", labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_Ucell(variables, parameters, ax):
    """This function plots the cell voltage as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    ax : matplotlib.axes.Axes
        Axes on which the cell voltage will be plotted.
    """

    # Extraction of the parameters
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    Ucell_t = np.array(variables['Ucell'])[mask]

    # Plot the cell voltage: Ucell
    ax.plot(t, Ucell_t, color=colors(0), label=r'$\mathregular{U_{cell}}$')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Cell}$ $\mathbf{voltage}$ $\mathbf{U_{cell}}$ $\mathbf{\left( V \right)}$', labelpad=3)
    ax.legend([r'$\mathregular{U_{cell}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)


def plot_P(variables, parameters, ax):
    """This function plots the pressure at different spatial localisations as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    ax : matplotlib.axes.Axes
        Axes on which the pressure will be plotted.
    """

    # Extraction of the parameters
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    Pagc_t = np.array(variables['Pagc'])[mask] / 1e5 # Conversion in atm
    Pcgc_t = np.array(variables['Pcgc'])[mask] / 1e5 # Conversion in atm
    Pasm_t = np.array(variables['Pasm'])[mask] / 1e5 # Conversion in atm
    Paem_t = np.array(variables['Paem'])[mask] / 1e5 # Conversion in atm
    Pcsm_t = np.array(variables['Pcsm'])[mask] / 1e5 # Conversion in atm
    Pcem_t = np.array(variables['Pcem'])[mask] / 1e5 # Conversion in atm

    # Plot the pressure at different spatial localisations: P
    ax.plot(t, Pagc_t, color=colors(0))
    ax.plot(t, Pcgc_t, color=colors(6))
    ax.plot(t, Pasm_t, color=colors(7))
    ax.plot(t, Paem_t, color=colors(8))
    ax.plot(t, Pcsm_t, color=colors(9))
    ax.plot(t, Pcem_t, color=colors(3))
    ax.legend([r'$\mathregular{P_{agc}}$', r'$\mathregular{P_{cgc}}$', r'$\mathregular{P_{asm}}$',
               r'$\mathregular{P_{aem}}$', r'$\mathregular{P_{csm}}$', r'$\mathregular{P_{cem}}$'], loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Pressure}$ $\mathbf{P}$ $\mathbf{\left( bar \right)}$', labelpad=3)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # Plot instructions
    plot_general_instructions(ax)


def plot_Phi_a(variables, operating_inputs, parameters, ax):
    """This function plots the humidity at the anode side, at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    ax : matplotlib.axes.Axes
        Axes on which the humidity will be plotted.
    """

    # Extraction of the operating inputs and parameters
    Phi_a_des = operating_inputs['Phi_a_des']
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    C_v_agc_t = np.array(variables['C_v_agc'])[mask]
    T_agc_t = np.array(variables['T_agc'])[mask]
    Phi_asm_t = np.array(variables['Phi_asm'])[mask]
    Phi_aem_t = np.array(variables['Phi_aem'])[mask]

    # Calculate the humidity Phi
    Phi_agc_t = C_v_agc_t * R * T_agc_t / Psat(T_agc_t)

    # Plot the humidity at different spatial localisations: Phi
    ax.plot(t, Phi_agc_t, color=colors(0), label=r'$\mathregular{\Phi_{agc}}$')
    ax.plot(t, Phi_asm_t, color=colors(1), label=r'$\mathregular{\Phi_{asm}}$')
    ax.plot(t, Phi_aem_t, color=colors(2), label=r'$\mathregular{\Phi_{aem}}$')
    ax.plot(t, np.array([Phi_a_des]*len(t)), color='black', label=r'$\mathregular{\Phi_{a,des}}$')
    ax.legend(loc='center right', bbox_to_anchor=(1, 0.67))
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Humidity}$ $\mathbf{at}$ $\mathbf{the}$ $\mathbf{anode}$ $\mathbf{side}$ $\mathbf{\Phi}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_Phi_c(variables, operating_inputs, parameters, ax):
    """This function plots the humidity, at the cathode side, at different spatial localisations as a function of time.

    Parameters
    ----------
    ax.plot(t, np.array([Phi_a_des]*len(t)), color='black', label=r'$\mathregular{\Phi_{a,des}}$')
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    ax : matplotlib.axes.Axes
        Axes on which the humidity will be plotted.
    """

    # Extraction of the operating inputs and parameters
    Phi_c_des = operating_inputs['Phi_c_des']
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    C_v_cgc_t = np.array(variables['C_v_cgc'])[mask]
    T_cgc_t = np.array(variables['T_cgc'])[mask]
    Phi_csm_t = np.array(variables['Phi_csm'])[mask]
    Phi_cem_t = np.array(variables['Phi_cem'])[mask]

    # Calculate the humidity Phi
    Phi_cgc_t = C_v_cgc_t * R * T_cgc_t / Psat(T_cgc_t)

    # Plot the humidity at different spatial localisations: Phi
    ax.plot(t, Phi_cgc_t, color=colors(0), label=r'$\mathregular{\Phi_{cgc}}$')
    ax.plot(t, Phi_csm_t, color=colors(1), label=r'$\mathregular{\Phi_{csm}}$')
    ax.plot(t, Phi_cem_t, color=colors(2), label=r'$\mathregular{\Phi_{cem}}$')
    ax.plot(t, np.array([Phi_c_des]*len(t)), color='black', label=r'$\mathregular{\Phi_{c,des}}$')
    ax.legend(loc='best')
    ax.set_xlabel(r'$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'$\mathbf{Humidity}$ $\mathbf{at}$ $\mathbf{the}$ $\mathbf{cathode}$ $\mathbf{side}$ $\mathbf{\Phi}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_Phi_des(variables, operating_inputs, parameters, ax):
    """This function plots the controlled or uncontrolled desired humidity at the anode and cathode as a function of the
    current density.

    Parameters
    ax.plot(t, np.array([Phi_c_des]*len(t)), color='black', label=r'$\mathregular{\Phi_{c,des}}$')
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the humidity will be plotted.
    """

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    pola_current_parameters = parameters['pola_current_parameters']
    type_current, type_plot = parameters['type_current'], parameters['type_plot']
    if type_current == 'step':
        delta_t_ini = parameters['step_current_parameters']['delta_t_ini_step']
    elif type_current == 'polarization':
        delta_t_ini = parameters['pola_current_parameters']['delta_t_ini_pola']
    elif type_current == 'polarization_for_cali':
        delta_t_ini = parameters['pola_current_for_cali_parameters']['delta_t_ini_pola_cali']
    else:
        delta_t_ini = 0
    # Extraction of the variables
    if type_plot == "fixed":
        mask = np.array(variables['t']) >= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else: # type_plot == "dynamic"
        mask = np.ones_like(variables['t'], dtype=bool)
    t = np.array(variables['t'])[mask]
    if parameters['type_control'] == "Phi_des":
        Phi_a_des_t = variables['Phi_a_des'][mask]
        Phi_c_des_t = variables['Phi_c_des'][mask]
        ax.set_ylabel(r'$\mathbf{Controlled}$ $\mathbf{inlet}$ $\mathbf{humidity}$  $\mathbf{\Phi_{des}}$', labelpad=3)
    else:
        Phi_a_des_t = np.array([operating_inputs['Phi_a_des']] * len(t))
        Phi_c_des_t = np.array([operating_inputs['Phi_c_des']] * len(t))
        ax.set_ylabel(r'$\mathbf{Uncontrolled}$ $\mathbf{inlet}$ $\mathbf{humidity}$ $\mathbf{\Phi_{des}}$', labelpad=3)

    # Plot Phi_des
    n = len(t)
    ifc_t = np.zeros(n)
    for i in range(n):  # Creation of ifc_t
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²

    # Recovery of the internal states from the model after each stack stabilisation
    delta_t_ini_pola = pola_current_parameters['delta_t_ini_pola']
    delta_t_load_pola = pola_current_parameters['delta_t_load_pola']
    delta_t_break_pola = pola_current_parameters['delta_t_break_pola']
    nb_loads = int(pola_current_parameters['i_max_pola'] / pola_current_parameters['delta_i_pola'])  # Number of loads
    ifc_discretized_t = np.zeros(nb_loads)
    Phi_a_des_discretized_t, Phi_c_des_discretized_t = np.zeros(nb_loads), np.zeros(nb_loads)
    for i in range(nb_loads):
        t_load = delta_t_ini_pola + (i + 1) * (delta_t_load_pola + delta_t_break_pola)  # time for measurement
        idx = (np.abs(t - t_load)).argmin()  # the corresponding index
        ifc_discretized_t[i] = ifc_t[idx]  # the last value at the end of each load
        Phi_a_des_discretized_t[i] = Phi_a_des_t[idx]  # the last value at the end of each load
        Phi_c_des_discretized_t[i] = Phi_c_des_t[idx]  # the last value at the end of each load

    ax.scatter(ifc_discretized_t, Phi_c_des_discretized_t, color=colors(6), label=r'$\mathregular{\Phi_{c,des}}$')
    ax.set_xlabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                  labelpad=3)
    if parameters['type_auxiliary'] == "forced-convective_cathode_with_flow-through_anode" or \
       parameters['type_auxiliary'] == "no_auxiliary":
        ax.scatter(ifc_discretized_t, Phi_a_des_discretized_t, color=colors(0), label=r'$\mathregular{\Phi_{a,des}}$')
        ax.legend([r'$\mathregular{\Phi_{a,des}}$', r'$\mathregular{\Phi_{c,des}}$'], loc='best')
    else:
        ax.legend([r'$\mathregular{\Phi_{c,des}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)


# ___________________________________________________Global indicators__________________________________________________

def plot_power_density_curve(variables, operating_inputs, parameters, n, ax):
    """This function plots the power density curve Pfc, produced by a cell, as a function of the current density.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    n : int
        Number of points used to plot the power density curve.
    ax : matplotlib.axes.Axes
        Axes on which the power density curve will be plotted.
    """

    # Extraction of the variables
    t, Ucell_t = variables['t'], variables['Ucell']
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    type_fuel_cell, type_current = parameters['type_fuel_cell'], parameters['type_current']
    type_auxiliary, type_control = parameters['type_auxiliary'], parameters['type_control']

    # Creation of the power density function: Pfc
    ifc_t, Pfc_t = np.zeros(n), np.zeros(n)
    for i in range(n):
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        Pfc_t[i] = Ucell_t[i] * ifc_t[i]

    # Plot of the power density function: Pfc
    plot_specific_line(ifc_t, Pfc_t, type_fuel_cell, type_current, type_auxiliary, type_control, None, ax)
    ax.set_xlabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                  labelpad=0)
    ax.set_ylabel(r'$\mathbf{Fuel}$ $\mathbf{cell}$ $\mathbf{power}$ $\mathbf{density}$ $\mathbf{P_{fc}}$ $\mathbf{\left( W.cm^{-2} \right)}$',
                  labelpad=0)
    ax.legend(loc='best')

    # Plot instructions
    plot_general_instructions(ax)


def plot_cell_efficiency(variables, operating_inputs, parameters, n, ax):
    """This function plots the fuel cell efficiency eta_fc as a function of the current density.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    n : int
        Number of points used to plot the fuel cell efficiency.
    ax : matplotlib.axes.Axes
        Axes on which the fuel cell efficiency will be plotted.
    """

    # Extraction of the variables
    t, Ucell_t, lambda_mem_t = variables['t'], variables['Ucell'], variables['lambda_mem']
    C_H2_acl_t, C_O2_ccl_t = variables['C_H2_acl'], variables['C_O2_ccl']
    T_acl_t, T_mem_t, T_ccl_t = variables['T_acl'], variables['T_mem'], variables['T_ccl']
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    Hmem, Hacl, Hccl, kappa_co = parameters['Hmem'], parameters['Hacl'], parameters['Hccl'], parameters['kappa_co']
    type_fuel_cell, type_current = parameters['type_fuel_cell'], parameters['type_current']
    type_auxiliary, type_control = parameters['type_auxiliary'], parameters['type_control']

    # Creation of the fuel cell efficiency: eta_fc
    ifc_t, Pfc_t, eta_fc_t = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        Pfc_t[i] = Ucell_t[i] * ifc_t[i]
        Ueq = E0 - 8.5e-4 * (T_ccl_t[i] - 298.15) + \
              R * T_ccl_t[i] / (2 * F) * (np.log(R * T_acl_t[i] * C_H2_acl_t[i] / Pref) +
                                          0.5 * np.log(R * T_ccl_t[i] * C_O2_ccl_t[i] / Pref))
        T_acl_mem_ccl = average([T_acl_t[i], T_mem_t[i], T_ccl_t[i]],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
        i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl_t[i] * k_H2(lambda_mem_t[i], T_mem_t[i], kappa_co)
        i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl_t[i] * k_O2(lambda_mem_t[i], T_mem_t[i], kappa_co)
        i_n = (i_H2 + i_O2) / 1e4  # Conversion in A/cm²
        eta_fc_t[i] = Pfc_t[i] / (Ueq * (ifc_t[i] + i_n))

    # Plot of the fuel cell efficiency: eta_fc
    plot_specific_line(ifc_t, eta_fc_t, type_fuel_cell, type_current, type_auxiliary, type_control, None, ax)
    ax.set_xlabel(r'$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$',
                  labelpad=0)
    ax.set_ylabel(r'$\mathbf{Fuel}$ $\mathbf{cell}$ $\mathbf{efficiency}$ $\mathbf{\eta_{fc}}$', labelpad=0)
    ax.legend(loc='best')

    # Plot instructions
    plot_general_instructions(ax)


# ______________________________________________________Instructions____________________________________________________

def calculate_simulation_error(Ucell, U_exp_t):
    """This function calculates the simulation error between the simulated cell voltage and the experimental cell
    voltage. It is calculated as the maximum relative difference between the two voltages (in %).

    Parameters
    ----------
    Ucell : numpy.ndarray
        Simulated cell voltage.
    U_exp_t : numpy.ndarray
        Experimental cell voltage.

    Returns
    -------
    float
        Simulation error between the simulated cell voltage and the experimental cell voltage (in %).
    """
    return np.round(np.max(np.abs(Ucell - U_exp_t) / U_exp_t * 100), 2)  # in %.


def plot_specific_line(x, y, type_fuel_cell, type_current, type_auxiliary, type_control, sim_error, ax):
    """ This function adds the appropriate plot configuration according to the type_input to the ax object.

    Parameters
    ----------
    x : numpy.ndarray
        x-axis values.
    y : numpy.ndarray
        y-axis values.
    type_fuel_cell : str
        Type of fuel cell configuration.
    type_current : str
        Type of current density.
    type_auxiliary : str
        Type of auxiliary system.
    type_control : str
        Type of control system.
    sim_error : float
        Simulation error between the simulated cell voltage and the experimental cell voltage (in %).
    ax : matplotlib.axes.Axes
        Axes on which the line will be plotted.
    """

    # For EH-31 fuel cell
    if type_current == "polarization":
        # ZSW fuel cell
        if type_fuel_cell == "ZSW-GenStack" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(0), label='Sim. - nominal operating conditions' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(0), label='Sim. - nominal operating conditions')

        elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(1), label='Sim. - P$_a$ = 1.61 bar - P$_c$ = 1.41 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(1), label='Sim. - P$_a$ = 1.61 bar - P$_c$ = 1.41 bar')

        elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(2), label='Sim. - P$_a$ = 2.01 bar - P$_c$ = 1.81 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(2), label='Sim. - P$_a$ = 2.01 bar - P$_c$ = 1.81 bar')

        elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(3), label='Sim. - P$_a$ = 2.4 bar - P$_c$ = 2.2 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(3), label='Sim. - P$_a$ = 2.4 bar - P$_c$ = 2.2 bar')

        elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(4), label='Sim. - P$_a$ = 2.8 bar - P$_c$ = 2.6 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(4), label='Sim. - P$_a$ = 2.8 bar - P$_c$ = 2.6 bar')

        elif type_fuel_cell == "ZSW-GenStack_T_62" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(5), label='Sim. - T = 62 °C' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_T_62" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(5), label='Sim. - T = 62 °C')

        elif type_fuel_cell == "ZSW-GenStack_T_76" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(6), label='Sim. - T = 76 °C' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_T_76" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(6), label='Sim. - T = 76 °C')

        elif type_fuel_cell == "ZSW-GenStack_T_84" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(7), label='Sim. - T = 84 °C' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_T_84" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(7), label='Sim. - T = 84 °C')

        # EH-31 fuel cell
        elif type_fuel_cell == "EH-31_1.5" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(0), label='Sim. - P = 1.5 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_1.5" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(0), label='Sim. - P = 1.5 bar')

        elif type_fuel_cell == "EH-31_2.0" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(1),
                    label='Sim. - P = 2.0 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.0" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            if type_control == "Phi_des":
                ax.plot(x, y, color=colors(5),
                        label=r'Sim. - P = 2.0 bar - controlled $\mathregular{\Phi_{des}}$')
            else:
                ax.plot(x, y, color=colors(1),
                        label=r'Sim. - P = 2.0 bar - uncontrolled $\mathregular{\Phi_{des}}$')

        elif type_fuel_cell == "EH-31_2.25" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, '--', color=colors(2),
                    label='Sim. - P = 2.25 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.25" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(2), label='Sim. - P = 2.25 bar')

        elif type_fuel_cell == "EH-31_2.5" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(3), label='Sim - P = 2.5 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.5" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.plot(x, y, color=colors(3), label='Sim - P = 2.5 bar')

        # For other fuel cell
        else:
            ax.plot(x, y, color=colors(0), label='Simulation')

    elif type_current == "polarization_for_cali":
        # ZSW fuel cell
        if type_fuel_cell == "ZSW-GenStack" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(0),
                       label='Sim. - nominal operating conditions' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(0), label='Sim. - nominal operating conditions')

        elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(1),
                       label='Sim. - P$_a$ = 1.61 bar - P$_c$ = 1.41 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(1), label='Sim. - P$_a$ = 1.61 bar - P$_c$ = 1.41 bar')

        elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(2),
                       label='Sim. - P$_a$ = 2.01 bar - P$_c$ = 1.81 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(2), label='Sim. - P$_a$ = 2.01 bar - P$_c$ = 1.81 bar')

        elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(3),
                       label='Sim. - P$_a$ = 2.4 bar - P$_c$ = 2.2 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(3), label='Sim. - P$_a$ = 2.4 bar - P$_c$ = 2.2 bar')

        elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(4),
                       label='Sim. - P$_a$ = 2.8 bar - P$_c$ = 2.6 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(4), label='Sim. - P$_a$ = 2.8 bar - P$_c$ = 2.6 bar')

        elif type_fuel_cell == "ZSW-GenStack_T_62" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(5),
                       label='Sim. - T = 62 °C' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_T_62" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(5), label='Sim. - T = 62 °C')

        elif type_fuel_cell == "ZSW-GenStack_T_76" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(6),
                       label='Sim. - T = 76 °C' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_T_76" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(6), label='Sim. - T = 76 °C')

        elif type_fuel_cell == "ZSW-GenStack_T_84" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(7),
                       label='Sim. - T = 84 °C' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "ZSW-GenStack_T_84" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(7), label='Sim. - T = 84 °C')

        # EH-31 fuel cell
        if type_fuel_cell == "EH-31_1.5" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(0),
                       label='Sim. - P = 1.5 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_1.5" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(0), label='Sim. - P = 1.5 bar')

        elif type_fuel_cell == "EH-31_2.0" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(1),
                       label='Sim. - P = 2.0 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.0" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            if type_control == "Phi_des":
                ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(5),
                           label=r'Sim. - P = 2.0 bar - controlled $\mathregular{\Phi_{des}}$')
            else:
                ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(1),
                           label=r'Sim. - P = 2.0 bar - uncontrolled $\mathregular{\Phi_{des}}$')

        elif type_fuel_cell == "EH-31_2.25" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(2),
                       label='Sim. - P = 2.25 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.25" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(2), label='Sim. - P = 2.25 bar')

        elif type_fuel_cell == "EH-31_2.5" and type_auxiliary == "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(3),
                       label='Sim - P = 2.5 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.5" and type_auxiliary != "forced-convective_cathode_with_flow-through_anode":
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(3), label='Sim - P = 2.5 bar')

        # For other fuel cell
        else:
            ax.scatter(x, y, marker='o', linewidths=1.5, color=colors(0), label='Simulation')

    else:
        raise ValueError('Only "polarization_current" and "polarization_current_for_cali" are considered here.')


def plot_general_instructions(ax, set_y = True):
    """This function adds the common instructions for all the plots displayed by AlphaPEM to the ax object.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    """
    # Get the current x-axis and y-axis limits
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    # Calculate the major step for the x-axis and y-axis ticks
    major_step_x = (x_max - x_min) / 5
    major_step_y = (y_max - y_min) / 5
    major_step_x_rounded = round_nice(major_step_x)
    major_step_y_rounded = round_nice(major_step_y)
    # Set the major and minor locators for the x-axis and y-axis
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(major_step_x_rounded))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(major_step_x_rounded / 5))
    if set_y:
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(major_step_y_rounded))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(major_step_y_rounded / 5))
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis='both', which='major', size=10, width=1.5, direction='out')
    ax.tick_params(axis='both', which='minor', size=5, width=1.5, direction='out')
    plt.tight_layout() # Adjust layout to prevent overlap between labels and the figure
    plt.show() # Show the figure

def round_nice(x):
    """Round the main step to a "nice" number.

    Parameters
    ----------
    x : float
        The value to be rounded.

    Returns
    -------
    float
        The value rounded to a "nice" number.
    """
    exp = np.floor(np.log10(x))
    f = x / 10**exp
    if f < 1.5:
        nice = 1
    elif f < 3:
        nice = 2
    elif f < 7:
        nice = 5
    else:
        nice = 10
    return nice * 10**exp

def plot_pola_instructions(type_fuel_cell, ax, show = True):
    """This function adds the specific instructions for polarisation plots according to the type_input to the ax object.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    show : bool, optional
        If True, the figure will be displayed. Default is True.
    """

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or \
            type_fuel_cell == "EH-31_2.25" or type_fuel_cell == "EH-31_2.5":
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1 / 5))
        ax.set_xlim(0, 3.0)
        ax.set_ylim(0.4, 1.04)

    # For other fuel cell
    else:
        pass

    # Configure the appearance of major and minor ticks
    ax.tick_params(axis='both', which='major', size=10, width=1.5, direction='out')
    ax.tick_params(axis='both', which='minor', size=5, width=1.5, direction='out')
    plt.tight_layout()  # Adjust layout to prevent overlap between labels and the figure
    if show:
        plt.show()  # Show the figure

def plot_EIS_Nyquist_instructions(type_fuel_cell, f_Fourier, x, y, ax):
    """This function adds the instructions for EIS plots according to the type_input to the ax object.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.
    f_Fourier : numpy.ndarray
        Frequency at which the EIS is simulated.
    x : numpy.ndarray
        x-axis values for plotting the annotation.
    y : numpy.ndarray
        y-axis values for plotting the annotation.
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    """

    # Commun instructions
    ax.set_aspect('equal', adjustable='box')  # Set orthonormal axis.
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis='both', which='major', size=10, width=1.5, direction='out')
    ax.tick_params(axis='both', which='minor', size=5, width=1.5, direction='out')
    plt.tight_layout()  # Adjust layout to prevent overlap between labels and the figure
    plt.show()  # Show the figure

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or \
            type_fuel_cell == "EH-31_2.25" or type_fuel_cell == "EH-31_2.5":
        # Double charge transfer
        if (f_Fourier >= 70 and f_Fourier <= 80):
            freq_str = str(int(f_Fourier)) + ' Hz'  # Frequency annotation.
            ax.annotate(freq_str, (x, y), textcoords="offset points", xytext=(0, -40), ha='center', fontsize=14,
                        rotation=90, weight='bold')
        # Auxiliary system
        if (f_Fourier >= 0.14 and f_Fourier <= 0.16):
            freq_str = f'{f_Fourier:.2g} Hz'  # Frequency annotation.
            ax.annotate(freq_str, (x, y), textcoords="offset points", xytext=(0, 7), ha='center', fontsize=14,
                        rotation=90, weight='bold')
        if (f_Fourier >= 1.2 and f_Fourier <= 1.4):
            freq_str = f'{f_Fourier:.2g} Hz'  # Frequency annotation.
            ax.annotate(freq_str, (x, y), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=14,
                        rotation=90, weight='bold')
        # Diffusion
        if (f_Fourier >= 0.015 and f_Fourier <= 0.020):
            freq_str = f'{f_Fourier:.2g} Hz'  # Frequency annotation.
            ax.annotate(freq_str, (x, y), textcoords="offset points", xytext=(30, 0), ha='center', fontsize=14,
                        rotation=0, weight='bold')
        if (f_Fourier >= 0.9 and f_Fourier <= 1.1):
            freq_str = f'{f_Fourier:.2g} Hz'  # Frequency annotation.
            ax.annotate(freq_str, (x, y), textcoords="offset points", xytext=(0, 10), ha='center', fontsize=14,
                        rotation=90, weight='bold')
        if (f_Fourier >= 70 and f_Fourier <= 90):
            freq_str = str(int(f_Fourier)) + ' Hz'  # Frequency annotation.
            ax.annotate(freq_str, (x, y), textcoords="offset points", xytext=(0, -40), ha='center', fontsize=14,
                        rotation=90, weight='bold')
        if (f_Fourier >= 10000 and f_Fourier <= 12000):
            freq_str = str(int(f_Fourier)) + ' Hz'  # Frequency annotation.
            ax.annotate(freq_str, (x, y), textcoords="offset points", xytext=(35, 0), ha='center', fontsize=14,
                        rotation=0, weight='bold')
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(20 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10 / 5))
        ax.set_xlim(30, 200)
        ax.set_ylim(-25, 55)

def plot_Bode_amplitude_instructions(f_EIS, type_fuel_cell, ax):
    """This function adds the instructions for amplitude Bode plots according to the type_input to the ax object.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    """

    # Commun instructions
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS  # They are the frequency parameters for the EIS
    #                                                                    simulation.
    ax.set_xscale('log') # set logarithmic scale for the x-axis
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis='both', which='major', size=10, width=1.5, direction='out')
    ax.tick_params(axis='both', which='minor', size=5, width=1.5, direction='out')
    plt.tight_layout()  # Adjust layout to prevent overlap between labels and the figure
    plt.show()  # Show the figure

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or \
            type_fuel_cell == "EH-31_2.25" or type_fuel_cell == "EH-31_2.5":
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=f_power_max_EIS - f_power_min_EIS + 1))
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
                                              numticks=(f_power_max_EIS - f_power_min_EIS + 1) * len(np.arange(2, 10))))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(30))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(30 / 5))
        ax.set_xlim([10**f_power_min_EIS, 10**f_power_max_EIS])
        # ax.set_ylim(0, 200)

def plot_Bode_phase_instructions(f_EIS, type_fuel_cell, ax):
    """This function adds the instructions for phase Bode plots according to the type_input to the ax object.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    """

    # Commun instructions
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS  # They are the frequency parameters for the EIS
    #                                                                    simulation.
    ax.set_xscale('log')  # set logarithmic scale for the x-axis
    if not ax.yaxis_inverted():
        ax.invert_yaxis()  # Invert the y-axis
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis='both', which='major', size=10, width=1.5, direction='out')
    ax.tick_params(axis='both', which='minor', size=5, width=1.5, direction='out')
    plt.tight_layout()  # Adjust layout to prevent overlap between labels and the figure
    plt.show()  # Show the figure

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or \
            type_fuel_cell == "EH-31_2.25" or type_fuel_cell == "EH-31_2.5":
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks = f_power_max_EIS-f_power_min_EIS+1))
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
                                              numticks = (f_power_max_EIS-f_power_min_EIS+1)*len(np.arange(2, 10))))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(5))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5 / 5))
        ax.set_xlim([10**f_power_min_EIS, 10**f_power_max_EIS])
        # ax.set_ylim(0, 360)
