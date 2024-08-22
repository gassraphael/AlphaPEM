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
from modules.transitory_functions import Psat, k_H2, k_O2
from calibration.experimental_values import pola_exp_values, plot_experimental_polarisation_curve

# General edition
mpl.rcParams['font.family'] = 'cmr10'  # 'cmr10' for English characters and 'DejaVu Serif' for French ones
plt.rcParams['axes.formatter.use_mathtext'] = True
plt.rcParams['axes.linewidth'] = 2.0
colors = mpl.colormaps['tab10']

# PyCharm requirement for dynamic plot display
mpl.use("Qt5Agg")


# __________________________________________________Polarisation curve__________________________________________________

def plot_polarisation_curve(variables, operating_inputs, parameters, ax):
    """
    This function plots the model polarisation curve, and compare it to the experimental one (if it exists). The
    polarisation curve is a classical representation of the cell performances, showing the cell voltage as a function
    of the current density.
    To generate it, the current density is increased step by step, and the cell voltage is recorded at each step.
    The time for which this point is captured is determined using the following approach: at the beginning of each load,
    a delta_t_load_pola time is needed to raise the current density to its next value. Subsequently, a delta_t_break_pola
    time is observed to ensure the dynamic stability of the stack's variables before initiating a new load. Ideally,
    each polarisation point should be recorded at the end of each delta_t_break_pola time. However, due to the design of the
    increments to minimize program instability (as observed in step_current function), the end of each delta_t_break_pola
    time corresponds to the beginning of a new load. To ensure a stationary operation and accurate polarisation point
    measurements, it is recommended to take the polarisation point just before by subtracting a delta_t value from it.
    This adjustment allows for stable and consistent measurements during the stationary period.

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
    t_step, i_step, i_max_pola = parameters['t_step'], parameters['i_step'], parameters['i_max_pola']
    delta_pola = parameters['delta_pola']
    i_EIS, t_EIS, f_EIS = parameters['i_EIS'], parameters['t_EIS'], parameters['f_EIS']
    type_fuel_cell, type_auxiliary = parameters['type_fuel_cell'], parameters['type_auxiliary']
    type_control, type_plot = parameters['type_control'], parameters['type_plot']

    if type_plot == "fixed":
        # Creation of ifc_t
        n = len(t)
        ifc_t = np.zeros(n)
        for i in range(n):
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²

        # Recovery of ifc and Ucell from the model after each stack stabilisation
        delta_t_load_pola, delta_t_break_pola, delta_i_pola, delta_t_ini_pola = delta_pola
        nb_loads = int(i_max_pola / delta_i_pola + 1)  # Number of loads which are made
        ifc_discretized = np.zeros(nb_loads)
        Ucell_discretized = np.zeros(nb_loads)
        for i in range(nb_loads):
            t_load = delta_t_ini_pola + (i + 1) * (delta_t_load_pola + delta_t_break_pola) - delta_t_break_pola / 10
            #                                                                                    # time for measurement
            idx = (np.abs(t - t_load)).argmin()  # the corresponding index
            ifc_discretized[i] = ifc_t[idx]  # the last value at the end of each load
            Ucell_discretized[i] = Ucell_t[idx]  # the last value at the end of each load

        # Plot the experimental polarization curve and calculate the simulation error compared with experimental data
        if type_fuel_cell != "manual_setup" and type_auxiliary == "opened_anode":  # Experimental points are accessible
            # Plot of the experimental polarization curve
            i_exp_t, U_exp_t = pola_exp_values(type_fuel_cell)
            plot_experimental_polarisation_curve(type_fuel_cell, i_exp_t, U_exp_t, ax)
            # Calculate the simulation error compared with experimental data
            #       i_fc and Ucell are reduced to remain within experimental limits for comparison
            i_fc_reduced = ifc_discretized[(ifc_discretized >= i_exp_t[0]) & (ifc_discretized <= i_exp_t[-1])]
            Ucell_reduced = Ucell_discretized[(ifc_discretized >= i_exp_t[0]) & (ifc_discretized <= i_exp_t[-1])]
            #       Experimental points are interpolated to correspond to the model points
            U_exp_interpolated = interp1d(i_exp_t, U_exp_t, kind='linear')(i_fc_reduced)
            sim_error = calculate_simulation_error(Ucell_reduced, U_exp_interpolated)
        else:
            sim_error = None

        # Plot the model polarisation curve
        plot_specific_line(ifc_discretized, Ucell_discretized, type_fuel_cell, type_auxiliary, type_control, sim_error,
                           ax)
        plot_pola_instructions(type_fuel_cell, ax)

    else:  # type_plot == "dynamic"
        # Plot of the polarisation curve produced by the model
        delta_t_load_pola, delta_t_break_pola, delta_i_pola, delta_t_ini_pola = delta_pola
        idx = (np.abs(t - t[-1] + delta_t_break_pola / 10)).argmin()  # index for polarisation measurement
        ifc = np.array(current_density(t[idx], parameters) / 1e4)  # time for polarisation measurement
        Ucell = np.array(Ucell_t[idx])  # voltage measurement
        ax.plot(ifc, Ucell, 'og', markersize=2)

    # Add the common instructions for the plot
    ax.set_xlabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$',
                  labelpad=3, fontsize=18)
    ax.set_ylabel(r'Cell voltage $\mathregular{U_{cell}}$ $\mathregular{\left( V \right)}$',
                  labelpad=3, fontsize=18)
    plot_general_instructions(ax)
    if type_plot == "fixed":
        ax.legend(loc='best', frameon=False)


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
    t_EIS, max_step = parameters['t_EIS'], parameters['max_step']

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
    freq_t = fftfreq(N, max_step)[:N // 2]  # Recovery of all frequency values used by fft
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
    ax.set_xlabel(r'$Z_{real}$ ($m\Omega.cm^{2}$)', labelpad=3)
    ax.set_ylabel(r'$-Z_{imag}$ ($m\Omega.cm^{2}$)', labelpad=3)
    #       Plot instructions
    plot_general_instructions(ax)
    plot_EIS_Nyquist_instructions(type_fuel_cell, ax)


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
    ax.set_xlabel('Frequency (Hz, logarithmic scale)', labelpad=3)
    ax.set_ylabel(r'Impedance amplitude ($m\Omega.cm^{2}$)', labelpad=3)
    #   Plot instructions
    plot_general_instructions(ax)
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
    theta_i_t = np.angle(ifc_Fourier[0:N // 2])  # Recovery of all dephasing values calculated by fft
    theta_U = theta_U_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation
    theta_i = theta_i_t[np.argmax(A_period_t == A)]  # Dephasing at the frequency of the perturbation

    # Plot the angle Bode diagram
    ax.plot(f, ((theta_U - theta_i) * 180 / np.pi) % 360, 'o', color=colors(2),
            label='Angle Bode diagram')
    ax.set_xlabel('Frequency (Hz, logarithmic scale)', labelpad=3)
    ax.set_ylabel(r'Phase (°)', labelpad=3)
    #   Plot instructions
    plot_general_instructions(ax)
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
    max_step = parameters['max_step']
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

def plot_ifc(variables, operating_inputs, parameters, n, ax):
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

    # Extraction of the variables
    t = variables['t']
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']

    # Plot the current density: ifc
    ifc_t = np.zeros(n)
    for i in range(n):  # Creation of ifc_t
        ifc_t[i] = current_density(t[i], parameters) / 10000  # Conversion in A/cm²
    ax.plot(t, ifc_t, color=colors(0), label=r'$\mathregular{i_{fc}}$')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax.legend([r'$\mathregular{i_{fc}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))


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
    # Extraction of the variables
    t, S_sorp_acl_t, S_sorp_ccl_t = variables['t'], variables['S_sorp_acl'], variables['S_sorp_ccl'],
    J_lambda_mem_acl_t, J_lambda_mem_ccl_t = variables['J_lambda_mem_acl'], variables['J_lambda_mem_ccl']
    # Extraction of the operating inputs and the parameters
    Hcl = parameters['Hcl']

    # Plot the sorption and dissolved water flows: J
    J_sorp_acl, J_sorp_ccl = [x * Hcl for x in S_sorp_acl_t], [x * Hcl for x in S_sorp_ccl_t]  # Conversion in
    #                                                                                         mol.m⁻².s⁻¹ for comparison
    ax.plot(t, J_sorp_acl, color=colors(2))
    ax.plot(t, J_lambda_mem_acl_t, color=colors(3))
    ax.plot(t, J_sorp_ccl, color=colors(4))
    ax.plot(t, J_lambda_mem_ccl_t, color=colors(7))
    ax.legend([r'$\mathregular{J_{sorp,acl}}$', r'$\mathregular{J_{\lambda,mem,acl}}$', r'$\mathregular{J_{sorp,ccl}}$',
               r'$\mathregular{J_{\lambda,mem,ccl}}$'], loc='best')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Flows $\mathregular{J}$ $\mathregular{\left( mol.m^{-2}.s^{-1} \right)}$', labelpad=0)

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.02 / 5))
    plt.show()


def plot_C_v(variables, n_gdl, C_v_sat, n, ax):
    """This function plots the vapor concentrations at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    C_v_sat : float
        Saturation vapor concentration.
    n : int
        Number of points used to plot the vapor concentration.
    ax : matplotlib.axes.Axes
        Axes on which the vapor concentration will be plotted.
    """

    # Extraction of the variables
    t, C_v_agc_t, C_v_agdl_t = variables['t'], variables['C_v_agc'], variables[f'C_v_agdl_{n_gdl // 2}']
    C_v_acl_t, C_v_ccl_t = variables['C_v_acl'], variables['C_v_ccl']
    C_v_cgdl_t, C_v_cgc_t = variables[f'C_v_cgdl_{n_gdl // 2}'], variables['C_v_cgc']

    # Plot the vapor concentrations at different spatial localisations Cv
    C_v_sat_t = np.ones(n) * C_v_sat
    ax.plot(t, C_v_agc_t, color=colors(0))
    ax.plot(t, C_v_agdl_t, color=colors(1))
    ax.plot(t, C_v_acl_t, color=colors(2))
    ax.plot(t, C_v_ccl_t, color=colors(4))
    ax.plot(t, C_v_cgdl_t, color=colors(5))
    ax.plot(t, C_v_cgc_t, color=colors(6))
    ax.plot(t, C_v_sat_t, color='k')
    ax.legend([r'$\mathregular{C_{v,agc}}$', r'$\mathregular{C_{v,agdl}}$', r'$\mathregular{C_{v,acl}}$',
               r'$\mathregular{C_{v,ccl}}$', r'$\mathregular{C_{v,cgdl}}$', r'$\mathregular{C_{v,cgc}}$',
               r'$\mathregular{C_{v,sat}}$'], loc='best')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r"Vapor concentration $\mathregular{C_{v}}$ $\mathregular{\left( mol.m^{-3} \right)}$",
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1 / 5))
    ax.set_ylim(11, 15.5)


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

    # Extraction of the variables
    t, lambda_acl_t = variables['t'], variables['lambda_acl']
    lambda_mem_t, lambda_ccl_t = variables['lambda_mem'], variables['lambda_ccl']
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    type_current = parameters['type_current']

    # Plot the water content at different spatial localisations: lambda
    if type_current == "polarization":
        n = len(t)
        ifc_t = np.zeros(n)
        for i in range(n):  # Creation of i_fc
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        ax.plot(ifc_t, lambda_acl_t, color=colors(2))
        ax.plot(ifc_t, lambda_mem_t, color=colors(3))
        ax.plot(ifc_t, lambda_ccl_t, color=colors(4))
        ax.set_xlabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    else:
        ax.plot(t, lambda_acl_t, color=colors(2))
        ax.plot(t, lambda_mem_t, color=colors(3))
        ax.plot(t, lambda_ccl_t, color=colors(4))
        ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Water content $\mathregular{\lambda}$', labelpad=3)
    ax.legend([r'$\mathregular{\lambda_{acl}}$', r'$\mathregular{\lambda_{mem}}$',
               r'$\mathregular{\lambda_{ccl}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(3))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(3 / 5))


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
    n_gdl, type_current = parameters['n_gdl'], parameters['type_current']
    # Extraction of the variables
    t, s_agdl_t, s_acl_t = variables['t'], variables[f's_agdl_{n_gdl // 2}'], variables['s_acl']
    s_ccl_t, s_cgdl_t = variables['s_ccl'], variables[f's_cgdl_{n_gdl // 2}']

    # Plot the liquid water saturation at different spatial localisations: s
    if type_current == "polarization":
        n = len(t)
        ifc_t = np.zeros(n)
        for i in range(n):  # Creation of i_fc
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        ax.plot(ifc_t, s_agdl_t, color=colors(1))
        ax.plot(ifc_t, s_acl_t, color=colors(2))
        ax.plot(ifc_t, s_ccl_t, color=colors(4))
        ax.plot(ifc_t, s_cgdl_t, color=colors(5))
        ax.set_xlabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    else:
        ax.plot(t, s_agdl_t, color=colors(1))
        ax.plot(t, s_acl_t, color=colors(2))
        ax.plot(t, s_ccl_t, color=colors(4))
        ax.plot(t, s_cgdl_t, color=colors(5))
        ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Liquid water saturation $\mathregular{s}$', labelpad=3)
    ax.legend([r'$\mathregular{s_{agdl}}$', r'$\mathregular{s_{acl}}$',
               r'$\mathregular{s_{ccl}}$', r'$\mathregular{s_{cgdl}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.04))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.04 / 5))


def plot_C_H2(variables, n_gdl, ax):
    """This function plots the hydrogen concentration at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    ax : matplotlib.axes.Axes
        Axes on which the hydrogen concentration will be plotted.
    """

    # Extraction of the variables
    t, C_H2_agc_t = variables['t'], variables['C_H2_agc']
    C_H2_agdl_t, C_H2_acl_t = variables[f'C_H2_agdl_{n_gdl // 2}'], variables['C_H2_acl']

    # Plot the hydrogen concentration at different spatial localisations: C_H2
    ax.plot(t, C_H2_agc_t, color=colors(0))
    ax.plot(t, C_H2_agdl_t, color=colors(1))
    ax.plot(t, C_H2_acl_t, color=colors(2))
    ax.legend([r'$\mathregular{C_{H_{2},agc}}$', r'$\mathregular{C_{H_{2},agdl}}$',
               r'$\mathregular{C_{H_{2},acl}}$'], loc='best')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Hydrogen concentration $\mathregular{C_{H_{2}}}$ $\mathregular{\left( mol.m^{-3} \right)}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1 / 5))
    ax.set_ylim(55, 58)


def plot_C_O2(variables, n_gdl, ax):
    """This function plots the oxygen concentration at different spatial localisations, as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    n_gdl : int
        Number of model nodes placed inside each GDL.
    ax : matplotlib.axes.Axes
        Axes on which the oxygen concentration will be plotted.
    """

    # Extraction of the variables
    t, C_O2_ccl_t = variables['t'], variables['C_O2_ccl']
    C_O2_cgdl_t, C_O2_cgc_t = variables[f'C_O2_cgdl_{n_gdl // 2}'], variables['C_O2_cgc']

    # Plot the oxygen concentration at different spatial localisations: C_O2
    ax.plot(t, C_O2_ccl_t, color=colors(4))
    ax.plot(t, C_O2_cgdl_t, color=colors(5))
    ax.plot(t, C_O2_cgc_t, color=colors(6))
    ax.legend([r'$\mathregular{C_{O_{2},ccl}}$', r'$\mathregular{C_{O_{2},cgdl}}$',
               r'$\mathregular{C_{O_{2},cgc}}$'], loc='best')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Oxygen concentration $\mathregular{C_{O_{2}}}$ $\mathregular{\left( mol.m^{-3} \right)}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1 / 5))
    ax.set_ylim(6, 11)


def plot_C_N2(variables, ax):
    """This function plots the nitrogen concentration as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    ax : matplotlib.axes.Axes
        Axes on which the nitrogen concentration will be plotted.
    """

    # Extraction of the variables
    t, C_N2_t = variables['t'], variables['C_N2']

    # Plot C_N2
    ax.plot(t, C_N2_t, color=colors(6))
    ax.legend([r'$\mathregular{C_{N_{2}}}$'], loc='best')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Nitrogen concentration $\mathregular{C_{N_{2}}}$ $\mathregular{\left( mol.m^{-3} \right)}$',
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_Ucell(variables, ax):
    """This function plots the cell voltage as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    ax : matplotlib.axes.Axes
        Axes on which the cell voltage will be plotted.
    """

    # Extraction of the variables
    t, Ucell_t = variables['t'], variables['Ucell']

    # Plot the cell voltage: Ucell
    ax.plot(t, Ucell_t, color=colors(0), label=r'$\mathregular{U_{cell}}$')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Cell voltage $\mathregular{U_{cell}}$ $\mathregular{\left( V \right)}$', labelpad=3)
    ax.legend([r'$\mathregular{U_{cell}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.05 / 5))


def plot_P(variables, ax):
    """This function plots the pressure at different spatial localisations as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    ax : matplotlib.axes.Axes
        Axes on which the pressure will be plotted.
    """

    # Extraction of the variables
    t, Pagc_t, Pcgc_t = variables['t'], variables['Pagc'], variables['Pcgc']
    Pasm_t, Paem_t, Pcsm_t, Pcem_t = variables['Pasm'], variables['Paem'], variables['Pcsm'], variables['Pcem']

    # Conversion in atm
    Pagc_t, Pcgc_t, Pasm_t = [x / 1e5 for x in Pagc_t], [x / 1e5 for x in Pcgc_t], [x / 1e5 for x in Pasm_t]
    Paem_t, Pcsm_t, Pcem_t = [x / 1e5 for x in Paem_t], [x / 1e5 for x in Pcsm_t], [x / 1e5 for x in Pcem_t]

    # Plot the pressure at different spatial localisations: P
    ax.plot(t, Pagc_t, color=colors(0))
    ax.plot(t, Pcgc_t, color=colors(6))
    ax.plot(t, Pasm_t, color=colors(7))
    ax.plot(t, Paem_t, color=colors(8))
    ax.plot(t, Pcsm_t, color=colors(9))
    ax.plot(t, Pcem_t, color=colors(3))
    ax.legend([r'$\mathregular{P_{agc}}$', r'$\mathregular{P_{cgc}}$', r'$\mathregular{P_{asm}}$',
               r'$\mathregular{P_{aem}}$', r'$\mathregular{P_{csm}}$', r'$\mathregular{P_{cem}}$'], loc='best')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Pressure $\mathregular{P}$ $\mathregular{\left( bar \right)}$', labelpad=3)
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(200))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(200 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5e-4))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5e-4 / 5))


def plot_Phi(variables, operating_inputs, ax):
    """This function plots the humidity at different spatial localisations as a function of time.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    ax : matplotlib.axes.Axes
        Axes on which the humidity will be plotted.
    """

    # Extraction of the variables
    t, C_agc_t, C_cgc_t = variables['t'], variables['C_agc'], variables['C_cgc']
    Phi_asm_t, Phi_aem_t = variables['Phi_asm'], variables['Phi_aem']
    Phi_csm_t, Phi_cem_t = variables['Phi_csm'], variables['Phi_cem']
    # Extraction of the operating inputs
    Tfc = operating_inputs['Tfc']
    Phi_a_des, Phi_c_des = operating_inputs['Phi_a_des'], operating_inputs['Phi_c_des']

    # Calculate the humidity Phi
    Phi_agc_t, Phi_cgc_t = [0] * len(t), [0] * len(t)
    for i in range(len(t)):
        Phi_agc_t[i] = C_agc_t[i] * R * Tfc / Psat(Tfc)
        Phi_cgc_t[i] = C_cgc_t[i] * R * Tfc / Psat(Tfc)

    # Plot the humidity at different spatial localisations: Phi
    ax.plot(t, Phi_agc_t, color=colors(0))
    ax.plot(t, Phi_cgc_t, color=colors(6))
    ax.plot(t, Phi_asm_t, color=colors(7))
    ax.plot(t, Phi_aem_t, color=colors(8))
    ax.plot(t, Phi_csm_t, color=colors(9))
    ax.plot(t, Phi_cem_t, color=colors(3))
    ax.plot(t, [Phi_a_des] * len(t), color='black')
    ax.plot(t, [Phi_c_des] * len(t), color='dimgrey')
    ax.legend([r'$\mathregular{\Phi_{agc}}$', r'$\mathregular{\Phi_{cgc}}$',
               r'$\mathregular{\Phi_{asm}}$', r'$\mathregular{\Phi_{aem}}$',
               r'$\mathregular{\Phi_{csm}}$', r'$\mathregular{\Phi_{cem}}$',
               r'$\mathregular{\Phi_{a,des}}$', r'$\mathregular{\Phi_{c,des}}$'], loc='best')
    ax.set_xlabel(r'Time $\mathregular{t}$ $\mathregular{\left( s \right)}$', labelpad=3)
    ax.set_ylabel(r'Humidity $\mathregular{\Phi}$', labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)


def plot_Phi_des(variables, operating_inputs, parameters, ax):
    """This function plots the controlled or uncontrolled desired humidity at the anode and cathode as a function of the
    current density.

    Parameters
    ----------
    variables : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
    operating_inputs : dict
        Operating inputs of the fuel cell.
    parameters : dict
        Parameters of the fuel cell model.
    ax : matplotlib.axes.Axes
        Axes on which the humidity will be plotted.
    """

    # Extraction of the variables
    t = variables['t']
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs['current_density']
    if parameters['type_control'] == "Phi_des":
        Phi_a_des, Phi_c_des = variables['Phi_a_des'], variables['Phi_c_des']
        ax.set_ylabel(r'Controlled inlet humidity  $\mathregular{\Phi_{des}}$', labelpad=3)
    else:
        Phi_a_des, Phi_c_des = [operating_inputs['Phi_a_des']] * len(t), [operating_inputs['Phi_c_des']] * len(t)
        ax.set_ylabel(r'Uncontrolled inlet humidity $\mathregular{\Phi_{des}}$', labelpad=3)

    # Plot Phi_des
    n = len(t)
    ifc_t = np.zeros(n)
    for i in range(n):  # Creation of ifc_t
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
    ax.plot(ifc_t, Phi_c_des, color=colors(6), label=r'$\mathregular{\Phi_{c,des}}$')
    ax.set_xlabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    if parameters['type_auxiliary'] == "opened_anode" or parameters['type_auxiliary'] == "no_auxiliary":
        ax.plot(t, Phi_a_des, color=colors(0), label=r'$\mathregular{\Phi_{a,des}}$')
        ax.legend([r'$\mathregular{\Phi_{a,des}}$', r'$\mathregular{\Phi_{c,des}}$'], loc='best')
    else:
        ax.legend([r'$\mathregular{\Phi_{c,des}}$'], loc='best')

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2 / 5))
    ax.set_xlim(0, 4.1)


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
    type_fuel_cell, type_auxiliary = parameters['type_fuel_cell'], parameters['type_auxiliary']
    type_control = parameters['type_control']

    # Creation of the power density function: Pfc
    ifc_t, Pfc_t = np.zeros(n), np.zeros(n)
    for i in range(n):
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        Pfc_t[i] = Ucell_t[i] * ifc_t[i]

    # Plot of the power density function: Pfc
    plot_specific_line(ifc_t, Pfc_t, type_fuel_cell, type_auxiliary, type_control, None, ax)
    ax.set_xlabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax.set_ylabel(r'Fuel cell power density $\mathregular{P_{fc}}$ $\mathregular{\left( W.cm^{-2} \right)}$',
                  labelpad=3)
    ax.legend(loc='best', frameon=False)

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.3))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.3 / 5))
    ax.set_xlim(0, 4.1)
    ax.set_ylim(0, 2.1)


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
    # Extraction of the operating inputs and the parameters
    current_density, Tfc = operating_inputs['current_density'], operating_inputs['Tfc']
    Hmem, kappa_co = parameters['Hmem'], parameters['kappa_co']
    type_fuel_cell, type_auxiliary = parameters['type_fuel_cell'], parameters['type_auxiliary']
    type_control = parameters['type_control']

    # Creation of the fuel cell efficiency: eta_fc
    ifc_t, Pfc_t, eta_fc_t = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(n):
        ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        Pfc_t[i] = Ucell_t[i] * ifc_t[i]
        Ueq = E0 - 8.5e-4 * (Tfc - 298.15) + R * Tfc / (2 * F) * (np.log(R * Tfc * C_H2_acl_t[i] / Pref) +
                                                                  0.5 * np.log(R * Tfc * C_O2_ccl_t[i] / Pref))
        i_H2 = 2 * F * R * Tfc / Hmem * C_H2_acl_t[i] * k_H2(lambda_mem_t[i], Tfc, kappa_co)
        i_O2 = 4 * F * R * Tfc / Hmem * C_O2_ccl_t[i] * k_O2(lambda_mem_t[i], Tfc, kappa_co)
        i_n = (i_H2 + i_O2) / 1e4  # Conversion in A/cm²
        eta_fc_t[i] = Pfc_t[i] / (Ueq * (ifc_t[i] + i_n))

    # Plot of the fuel cell efficiency: eta_fc
    plot_specific_line(ifc_t, eta_fc_t, type_fuel_cell, type_auxiliary, type_control, None, ax)
    ax.set_xlabel(r'Current density $\mathregular{i_{fc}}$ $\mathregular{\left( A.cm^{-2} \right)}$', labelpad=3)
    ax.set_ylabel(r'Fuel cell efficiency $\mathregular{\eta_{fc}}$', labelpad=3)
    ax.legend(loc='best', frameon=False)

    # Plot instructions
    plot_general_instructions(ax)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1 / 5))
    ax.set_xlim(0, 4.1)
    ax.set_ylim(0, 0.7)


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


def plot_specific_line(x, y, type_fuel_cell, type_auxiliary, type_control, sim_error, ax):
    """ This function adds the appropriate plot configuration according to the type_input to the ax object.

    Parameters
    ----------
    x : numpy.ndarray
        x-axis values.
    y : numpy.ndarray
        y-axis values.
    type_fuel_cell : str
        Type of fuel cell configuration.
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
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or type_fuel_cell == "EH-31_2.25" or \
            type_fuel_cell == "EH-31_2.5":
        if type_fuel_cell == "EH-31_1.5" and type_auxiliary == "opened_anode":
            ax.plot(x, y, color=colors(0), label='Sim. - P = 1.5 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_1.5" and type_auxiliary != "opened_anode":  # no auxiliaries
            ax.plot(x, y, color=colors(0), label='Sim. - P = 1.5 bar')
        elif type_fuel_cell == "EH-31_2.0" and type_auxiliary == "opened_anode":
            ax.plot(x, y, '--', color=colors(1),
                    label='Sim. - P = 2.0 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.0" and type_auxiliary != "opened_anode":  # no auxiliaries
            if type_control == "Phi_des":
                ax.plot(x, y, color=colors(5),
                        label=r'Sim. - P = 2.0 bar - controlled $\mathregular{\Phi_{des}}$')
            else:
                ax.plot(x, y, color=colors(1),
                        label=r'Sim. - P = 2.0 bar - uncontrolled $\mathregular{\Phi_{des}}$')
        elif type_fuel_cell == "EH-31_2.25" and type_auxiliary == "opened_anode":
            ax.plot(x, y, '--', color=colors(2),
                    label='Sim. - P = 2.25 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.25" and type_auxiliary != "opened_anode":  # no auxiliaries
            ax.plot(x, y, color=colors(2), label='Sim. - P = 2.25 bar')
        elif type_fuel_cell == "EH-31_2.5" and type_auxiliary == "opened_anode":
            ax.plot(x, y, color=colors(3), label='Sim - P = 2.5 bar' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "EH-31_2.5" and type_auxiliary != "opened_anode":  # no auxiliaries
            ax.plot(x, y, color=colors(3), label='Sim - P = 2.5 bar')

    # For BX fuel cell
    elif type_fuel_cell == "BX_1.0" or type_fuel_cell == "BX_1.35":
        if type_fuel_cell == "BX_1.0" and type_auxiliary == "opened_anode":
            ax.plot(x, y, '--', color=colors(0),
                    label='Sim. - P = 1.35 atm' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "BX_1.0" and type_auxiliary != "opened_anode":  # no auxiliaries
            ax.plot(x, y, '--', color=colors(0), label='Sim. - P = 1.35 atm')
        elif type_fuel_cell == "BX_1.35" and type_auxiliary == "opened_anode":
            ax.plot(x, y, color=colors(1), label='Sim. - P = 1.0 atm' + r' - $ΔU_{max}$ =' f' {sim_error} %')
        elif type_fuel_cell == "BX_1.35" and type_auxiliary != "opened_anode":  # no auxiliaries
            ax.plot(x, y, color=colors(1), label='Sim. - P = 1.0 atm')

    # For LF fuel cell
    elif type_fuel_cell == "LF":
        ax.plot(x, y, color=colors(0), label='Simulation')

    # For other fuel cell
    else:
        ax.plot(x, y, color=colors(0), label='Simulation')


def plot_general_instructions(ax):
    """This function adds the common instructions for all the plots displayed by AlphaPEM to the ax object.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    """

    ax.tick_params(axis='both', which='major', size=10, width=1.5, direction='out', labelsize=15)
    ax.tick_params(axis='both', which='minor', size=5, width=1.5, direction='out', labelsize=15)
    plt.show()


def plot_pola_instructions(type_fuel_cell, ax):
    """This function adds the specific instructions for polarisation plots according to the type_input to the ax object.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    """

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or \
            type_fuel_cell == "EH-31_2.25" or type_fuel_cell == "EH-31_2.5":
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1 / 5))
        ax.set_xlim(0, 3.0)
        ax.set_ylim(0.4, 1.0)

    # For BX fuel cell
    elif type_fuel_cell == "BX_1.0" or type_fuel_cell == "BX_1.35":
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.3))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.3 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2 / 5))
        ax.set_xlim(-0.1, 4.5)
        ax.set_ylim(0.1, 1.0)

    # For LF fuel cell
    elif type_fuel_cell == "LF":
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.4))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.4 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.2 / 5))
        ax.set_xlim(-0.05, 1.6)
        ax.set_ylim(0, 1.0)

    # For other fuel cell
    else:
        pass

def plot_EIS_Nyquist_instructions(type_fuel_cell, ax):
    """This function adds the instructions for EIS plots according to the type_input to the ax object.

    Parameters
    ----------
    type_fuel_cell : str
        Type of fuel cell configuration.
    ax : matplotlib.axes.Axes
        Axes on which the instructions will be added.
    """

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" or type_fuel_cell == "EH-31_2.0" or \
            type_fuel_cell == "EH-31_2.25" or type_fuel_cell == "EH-31_2.5":
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(40))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(40 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10 / 5))
        # ax.set_xlim(0, 200)
        # ax.set_ylim(-30, 60)

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