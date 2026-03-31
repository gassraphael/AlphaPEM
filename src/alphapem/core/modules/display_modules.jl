# -*- coding: utf-8 -*-

"""This module is used to accurately plot the figures.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using PyCall
using Interpolations

const np = pyimport("numpy")
const mpl = pyimport("matplotlib")
const plt = pyimport("matplotlib.pyplot")
const LogLocator = pyimport("matplotlib.ticker").LogLocator

# Importing constants' value and functions
include(joinpath(@__DIR__, "../../utils/physics_constants.jl"))
include(joinpath(@__DIR__, "../../utils/maths_functions.jl"))
include(joinpath(@__DIR__, "flows_1D_MEA_modules.jl"))
include(joinpath(@__DIR__, "../../utils/physics_functions.jl"))
include(joinpath(@__DIR__, "../../config/pola_exp_values.jl"))
include(joinpath(@__DIR__, "display_calc_modules.jl"))

# General edition
const colors = pyimport("matplotlib.pyplot").get_cmap("tab20")


# __________________________________________________Polarisation curve__________________________________________________

"""
    plot_polarisation_curve(variables, operating_inputs, parameters, ax, show=true)

This function plots the model polarisation curve, and compare it to the experimental one (if it exists). The
polarisation curve is a classical representation of the cell performances, showing the cell voltage as a function
of the current density.
To generate it, the current density is increased step by step, and the cell voltage is recorded at each step.
The time for which this point is captured is determined using the following approach: at the beginning of each load,
a delta_t_load_pola time is needed to raise the current density to its next value. Subsequently, a delta_t_break_pola
time is observed to ensure the dynamic stability of the stack's variables before initiating a new load. Finally,
each polarisation point is recorded at the end of each delta_t_break_pola time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver. They correspond to the fuel
  cell internal states.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the polarisation curve will be plotted.
- `show::Bool=true`: If true, the polarisation curve will be displayed. If false, it will not be displayed.
"""
function plot_polarisation_curve(variables::Dict{String, Any},
                                 operating_inputs::Dict{String, Any},
                                 parameters::Dict{String, Any},
                                 ax,
                                 show::Bool=true)

    # Extraction of the variables
    t       = collect(variables["t"])
    Ucell_t = collect(variables["Ucell"])
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    pola_current_parameters = parameters["pola_current_parameters"]
    delta_t_ini_pola = pola_current_parameters["delta_t_ini_pola"]
    delta_t_load_pola = pola_current_parameters["delta_t_load_pola"]
    delta_t_break_pola = pola_current_parameters["delta_t_break_pola"]
    delta_i_pola = pola_current_parameters["delta_i_pola"]
    i_max_pola = pola_current_parameters["i_max_pola"]
    type_fuel_cell = parameters["type_fuel_cell"]
    type_current = parameters["type_current"]
    voltage_zone = parameters["voltage_zone"]
    type_auxiliary = parameters["type_auxiliary"]
    type_plot = parameters["type_plot"]

    if type_plot == "fixed"
        # Creation of ifc_t
        n = length(t)
        ifc_t = zeros(n)
        for i in 1:n
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        end

        # Recovery of ifc and Ucell from the model after each stack stabilisation
            nb_loads = floor(Int, i_max_pola / delta_i_pola)  # Number of loads which are made
        ifc_discretized = zeros(nb_loads + 1)      # One point is taken at ifc = 0, before the first load.
        Ucell_discretized = zeros(nb_loads + 1)    # One point is taken at ifc = 0, before the first load.
        for i in 0:nb_loads
            t_load = delta_t_ini_pola + i * (delta_t_load_pola + delta_t_break_pola)  # time for measurement
            idx = argmin(abs.(t .- t_load))  # the corresponding index
            ifc_discretized[i + 1] = ifc_t[idx]  # the last value at the end of each load
            Ucell_discretized[i + 1] = Ucell_t[idx]  # the last value at the end of each load
        end

        # Plot the experimental polarization curve and calculate the simulation error compared with experimental data
        if type_fuel_cell != "manual_setup" &&
           (type_auxiliary == "forced-convective_cathode_with_flow-through_anode" || type_auxiliary == "no_auxiliary")
            # Extraction of the experimental current density and voltage values.
            i_exp_t, U_exp_t = pola_exp_values(type_fuel_cell, voltage_zone)  # (A.m-2, V).
            # Plot of the experimental polarization curve
            i_exp_t = i_exp_t ./ 1e4  # Conversion in A/cm²
            plot_experimental_polarisation_curve(type_fuel_cell, i_exp_t, U_exp_t, ax)
            # Calculate the simulation error compared with experimental data
            # Experimental points are interpolated to correspond to the model points
            itp = linear_interpolation(ifc_discretized, Ucell_discretized; extrapolation_bc=Line())
            Ucell_interpolated = itp.(i_exp_t)  # Extrapolates linearly beyond the input x-range.
            sim_error = calculate_simulation_error(Ucell_interpolated, U_exp_t)
        else
            sim_error = nothing
        end

        # Plot the model polarisation curve
        plot_specific_line(ifc_discretized, Ucell_discretized, type_fuel_cell, type_current, type_auxiliary,
                           sim_error, ax)
        plot_pola_instructions(type_fuel_cell, ax, show)

    else  # type_plot == "dynamic"
        # Plot of the polarisation curve produced by the model
        idx = argmin(abs.(t .- t[end]))  # index for polarisation measurement
        ifc = current_density(t[idx], parameters) / 1e4  # time for polarisation measurement
        Ucell = Ucell_t[idx]  # voltage measurement
        ax.plot([ifc], [Ucell], "og"; markersize=2)
    end

    # Add the common instructions for the plot
    ax.set_xlabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$";
                  labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Cell}$ $\mathbf{voltage}$ $\mathbf{U_{cell}}$ $\mathbf{\left( V \right)}$";
                  labelpad=3)
    if type_plot == "fixed"
        ax.legend(loc="best")
    end

    return nothing
end


"""
    plot_polarisation_curve_for_cali(variables, operating_inputs, parameters, ax)

This function plots the model polarisation curve, and compare it to the experimental one. The
polarisation curve is a classical representation of the cell performances, showing the cell voltage as a function
of the current density.
To generate it, the current density is increased step by step, and the cell voltage is recorded at each step.
The time for which this point is captured is determined using the following approach: at the beginning of each load,
a delta_t_load_pola time is needed to raise the current density to its next value. Subsequently, a delta_t_break_pola
time is observed to ensure the dynamic stability of the stack's variables before initiating a new load. Finally,
each polarisation point is recorded at the end of each delta_t_break_pola time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver. They correspond to the fuel
  cell internal states.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the polarisation curve will be plotted.
"""
function plot_polarisation_curve_for_cali(variables::Dict{String, Any},
                                          operating_inputs::Dict{String, Any},
                                          parameters::Dict{String, Any},
                                          ax)

    # Extraction of the variables
    t       = collect(variables["t"])
    Ucell_t = collect(variables["Ucell"])
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    pola_current_for_cali_parameters = parameters["pola_current_for_cali_parameters"]
    delta_t_ini_pola_cali = pola_current_for_cali_parameters["delta_t_ini_pola_cali"]
    delta_t_load_pola_cali = pola_current_for_cali_parameters["delta_t_load_pola_cali"]
    delta_t_break_pola_cali = pola_current_for_cali_parameters["delta_t_break_pola_cali"]
    type_fuel_cell = parameters["type_fuel_cell"]
    type_current = parameters["type_current"]
    voltage_zone = parameters["voltage_zone"]
    type_auxiliary = parameters["type_auxiliary"]
    type_plot = parameters["type_plot"]
    # Extraction of the experimental current density and voltage values for the calibration.
    i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(type_fuel_cell, voltage_zone)  # (A.m-2, V).

    if type_plot == "fixed"
        # Creation of ifc_t
        n = length(t)
        ifc_t = zeros(n)
        for i in 1:n
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        end

        # Recovery of ifc and Ucell from the model after each stack stabilisation
        nb_loads = length(i_exp_cali_t)  # Number of loads which are made
        delta_t_cali = delta_t_load_pola_cali + delta_t_break_pola_cali  # s. It is the time of one load.
        ifc_discretized = zeros(nb_loads)
        Ucell_discretized = zeros(nb_loads)
        for i in 1:nb_loads
            t_load = delta_t_ini_pola_cali + i * delta_t_cali  # time for measurement
            idx = argmin(abs.(t .- t_load))  # the corresponding index
            ifc_discretized[i] = ifc_t[idx]  # the last value at the end of each load
            Ucell_discretized[i] = Ucell_t[idx]  # the last value at the end of each load
        end

        # Plot the experimental polarization curve
        i_exp_cali_t = i_exp_cali_t ./ 1e4  # Conversion in A/cm²
        plot_experimental_polarisation_curve(type_fuel_cell, i_exp_cali_t, U_exp_cali_t, ax)

        # Plot the model polarisation curve
        sim_error = calculate_simulation_error(Ucell_discretized, U_exp_cali_t)  # Calculate the simulation error
        plot_specific_line(ifc_discretized, Ucell_discretized, type_fuel_cell, type_current, type_auxiliary,
                           sim_error, ax)
        plot_pola_instructions(type_fuel_cell, ax)

    else  # type_plot == "dynamic"
        # Plot of the polarisation curve produced by the model
        idx = argmin(abs.(t .- t[end]))  # index for polarisation measurement
        ifc = current_density(t[idx], parameters) / 1e4  # time for polarisation measurement
        Ucell = Ucell_t[idx]  # voltage measurement
        ax.plot([ifc], [Ucell], "og"; markersize=2)
    end

    # Add the common instructions for the plot
    ax.set_xlabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$";
                  labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Cell}$ $\mathbf{voltage}$ $\mathbf{U_{cell}}$ $\mathbf{\left( V \right)}$";
                  labelpad=3)
    if type_plot == "fixed"
        ax.legend(loc="best")
    end

    return nothing
end


# _______________________________________________________EIS curves_____________________________________________________

"""
    plot_EIS_curve_Nyquist(parameters, Fourier_results, ax)

This function is used to plot the Nyquist diagram of the EIS curves.

# Arguments
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `Fourier_results::Dict{String, Any}`: Dictionary containing the Fourier transformation (FT)
  post-processing outputs.
- `ax`: Axes on which the Nyquist diagram will be plotted.
"""
function plot_EIS_curve_Nyquist(parameters::Dict{String, Any},
                                Fourier_results::Dict{String, Any},
                                ax)

    # Extraction of the parameters
    i_EIS, ratio_EIS, type_fuel_cell = parameters["i_EIS"], parameters["ratio_EIS"], parameters["type_fuel_cell"]
    # Extraction of the Fourier results
    Ucell_Fourier, ifc_Fourier = Fourier_results["Ucell_Fourier"], Fourier_results["ifc_Fourier"]
    f_Fourier = Fourier_results["f"]
    A_period_t, A, N = Fourier_results["A_period_t"], Fourier_results["A"], Fourier_results["N"]

    # Calculation of the real and imaginary component of the impedance for each period
    Z0 = A / (ratio_EIS * (-i_EIS)) * 1e7  # Impedance of the perturbation in mΩ.cm². The sign of i is inverted to
    # comply with the standards of EIS, which measure a device under load rather than a current source.
    theta_U_t = angle.(Ucell_Fourier[1:N÷2])  # Recovery of all dephasing values calculated by fft
    theta_i_t = angle.(ifc_Fourier[1:N÷2])  # Recovery of all dephasing values calculated by fft
    idx_A = findfirst(A_period_t .== A)
    theta_U = theta_U_t[idx_A]  # Dephasing at the frequency of the perturbation
    theta_i = theta_i_t[idx_A]  # Dephasing at the frequency of the perturbation
    Z_real = Z0 * cos(theta_U - theta_i)  # Real component of the impedance for each period
    Z_imag = Z0 * sin(theta_U - theta_i)  # Imaginary component of the impedance for each period

    # Plot the Nyquist diagram
    ax.plot([Z_real], [-Z_imag], "o"; color=colors(0), label="Nyquist diagram")
    ax.set_xlabel(raw"$\mathbf{Z_{real}}$ $\mathbf{(m\Omega.cm^{2})}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{-Z_{imag}}$ $\mathbf{(m\Omega.cm^{2})}$"; labelpad=3)
    # Plot instructions
    plot_EIS_Nyquist_instructions(type_fuel_cell, f_Fourier, Z_real, -Z_imag, ax)

    return nothing
end


"""
    plot_EIS_curve_Bode_amplitude(parameters, Fourier_results, ax)

This function is used to plot the amplitude Bode diagram of the EIS curves.

# Arguments
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `Fourier_results::Dict{String, Any}`: Dictionary containing Fourier post-processing outputs.
- `ax`: Axes on which the amplitude Bode diagram will be plotted.
"""
function plot_EIS_curve_Bode_amplitude(parameters::Dict{String, Any},
                                       Fourier_results::Dict{String, Any},
                                       ax)

    # Extraction of the parameters
    i_EIS = parameters["i_EIS"]
    ratio_EIS = parameters["ratio_EIS"]
    f_EIS = parameters["f_EIS"]
    type_fuel_cell = parameters["type_fuel_cell"]
    # Extraction of the Fourier results
    A = Fourier_results["A"]
    f = Fourier_results["f"]

    # Calculation of the impedance of the perturbation
    Z0 = A / (ratio_EIS * (-i_EIS)) * 1e7  # in mΩ.cm². The sign of i is inverted to comply with EIS standards.

    # Plot the amplitude Bode diagram
    ax.plot([f], [abs(Z0)], "o"; color=colors(1), label="Amplitude Bode diagram")
    ax.set_xlabel(raw"$\mathbf{Frequency}$ $\mathbf{(Hz,}$ $\mathbf{logarithmic}$ $\mathbf{scale)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Impedance}$ $\mathbf{amplitude}$ $\mathbf{(m\Omega.cm^{2})}$"; labelpad=3)
    # Plot instructions
    plot_Bode_amplitude_instructions(f_EIS, type_fuel_cell, ax)

    return nothing
end


"""
    plot_EIS_curve_Bode_angle(parameters, Fourier_results, ax)

This function is used to plot the angle Bode diagram. It only works with an entry signal made with a cosinus
(not a sinus).

# Arguments
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `Fourier_results::Dict{String, Any}`: Dictionary containing Fourier post-processing outputs.
- `ax`: Axes on which the angle Bode diagram will be plotted.
"""
function plot_EIS_curve_Bode_angle(parameters::Dict{String, Any},
                                   Fourier_results::Dict{String, Any},
                                   ax)

    # Extraction of the parameters
    f_EIS = parameters["f_EIS"]
    type_fuel_cell = parameters["type_fuel_cell"]
    # Extraction of the Fourier results
    Ucell_Fourier = Fourier_results["Ucell_Fourier"]
    ifc_Fourier = Fourier_results["ifc_Fourier"]
    A_period_t = Fourier_results["A_period_t"]
    A = Fourier_results["A"]
    f = Fourier_results["f"]
    N = Fourier_results["N"]

    # Calculation of the dephasing values at the frequency of the perturbation
    theta_U_t = angle.(Ucell_Fourier[1:N÷2])  # Recovery of all dephasing values calculated by fft
    theta_i_t = angle.(ifc_Fourier[1:N÷2]) .+ π  # Recovery of all dephasing values calculated by fft.
    # An angle of pi is added to comply with EIS standards (device under load rather than current source).
    idx_A = findfirst(A_period_t .== A)
    theta_U = theta_U_t[idx_A]  # Dephasing at the frequency of the perturbation
    theta_i = theta_i_t[idx_A]  # Dephasing at the frequency of the perturbation
    phi_U_i = mod((theta_U - theta_i) * 180 / π, 360)  # Dephasing between Ucell and ifc with a value between 0 and 360
    if phi_U_i > 180
        phi_U_i -= 360  # To have a value between -180 and 180
    end

    # Plot the angle Bode diagram
    ax.plot([f], [phi_U_i], "o"; color=colors(2), label="Angle Bode diagram")
    ax.set_xlabel(raw"$\mathbf{Frequency}$ $\mathbf{(Hz,}$ $\mathbf{logarithmic}$ $\mathbf{scale)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Phase}$ $\mathbf{(^\circ)}$"; labelpad=3)
    # Plot instructions
    plot_Bode_phase_instructions(f_EIS, type_fuel_cell, ax)

    return nothing
end


"""
    plot_EIS_curve_tests(variables, operating_inputs, parameters, Fourier_results)

This function is used to test the accuracy of the EIS results. It compares the reconstructed Ucell_Fourier(t)
from the Fourier transformation with the current density ifc(t), and displays Ucell(t) given by the model with the
reconstructed Ucell_Fourier(t).

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver. They correspond to the fuel
  cell internal states.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `Fourier_results::Dict{String, Any}`: Dictionary containing Fourier post-processing outputs.
"""
function plot_EIS_curve_tests(variables::Dict{String, Any},
                              operating_inputs::Dict{String, Any},
                              parameters::Dict{String, Any},
                              Fourier_results::Dict{String, Any})

    # Extraction of the variables
    t = collect(variables["t"])
    Ucell_t = variables["Ucell"]
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    i_EIS = parameters["i_EIS"]
    ratio_EIS = parameters["ratio_EIS"]
    t_EIS = parameters["t_EIS"]
    f_EIS = parameters["f_EIS"]
    # Extraction of the Fourier results
    Ucell_Fourier = Fourier_results["Ucell_Fourier"]
    ifc_Fourier = Fourier_results["ifc_Fourier"]
    A_period_t = Fourier_results["A_period_t"]
    A = Fourier_results["A"]
    f = Fourier_results["f"]
    N = Fourier_results["N"]

    # Reconstructed Ucell with a cosinus form, and comparison of its form with the current density one.
    t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS
    n_inf = findlast(t_new_start_EIS .<= t[1])  # The number of frequency changes which has been made.
    f_current = 10 .^ range(f_power_min_EIS, f_power_max_EIS; length=nb_f_EIS)
    theta_U_t = angle.(Ucell_Fourier[1:N÷2])  # Recovery of all dephasing values calculated by fft
    theta_i_t = angle.(ifc_Fourier[1:N÷2])  # Recovery of all dephasing values calculated by fft
    idx_A = findfirst(A_period_t .== A)
    theta_U = theta_U_t[idx_A]  # Dephasing at the frequency of the perturbation
    theta_i = theta_i_t[idx_A]  # Dephasing at the frequency of the perturbation

    println("Ucell: ", round(A_period_t[1], digits=4), " + ", round(A, digits=6),
            " * cos(2*pi*", round(f, digits=4), "*t + ", round(theta_U, digits=4), ").")
    println("Current: ", i_EIS, " + ", ratio_EIS * i_EIS,
            " * cos(2*pi*", round(f_current[n_inf], digits=4), "*t + ", round(theta_i, digits=4), ").\n")

    # Display ifc(t)
    plt.figure(3)
    plt.subplot(2, 1, 1)
    # Creation of ifc_t
    n = length(t)
    ifc_t = zeros(n)
    for i in 1:n  # Conversion in A/cm²
        ifc_t[i] = current_density(t[i], parameters) / 1e4
    end
    # Plot of ifc_t
    plt.plot(t, ifc_t; color="blue", label="ifc")
    plt.xlabel("Time (s)")
    plt.ylabel("Current density (A/cm²)")
    plt.title("The current density\nbehaviour over time")

    # Display Ucell(t) and compare it with the reconstructed Ucell_Fourier(t) from the Fourier transformation
    plt.subplot(2, 1, 2)
    Ucell_Fourier_reconstructed = A_period_t[1] .+ A .* cos.(2 .* π .* f .* t .+ theta_U)
    plt.plot(t, Ucell_t; color="blue", label="Ucell")
    plt.plot(t, Ucell_Fourier_reconstructed; color="black", label="Ucell_Fourier")
    plt.xlabel("Time (s)")
    plt.ylabel("Cell voltage (V)")
    plt.title("The cell voltage\nbehaviour over time")

    return nothing
end


# ____________________________________________Internal variables - 1D temporal__________________________________________

"""
    plot_ifc_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the current density as a function of time. The different current density values at different
spatial localisation through the gas channel are plotted on the same graph, to compare their behaviour over time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the current density will be plotted.
"""
function plot_ifc_1D_temporal(variables::Dict{String, Any},
                              operating_inputs::Dict{String, Any},
                              parameters::Dict{String, Any},
                              ax)

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    nb_gc = parameters["nb_gc"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]

    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini  # select the time after 0.9*delta_t_ini
    else
        trues(length(t_full))
    end
    t = t_full[mask]
    # Python previously used [None] to force 1-based indexing. Julia is natively 1-based.
    i_fc_t = [collect(variables["i_fc"][i])[mask] ./ 1e4 for i in 1:nb_gc]  # Conversion in A/cm²

    # Collect handles and labels then build manually the legend
    handles = Any[]
    labels = String[]

    # Plot the current density: ifc
    n = length(t)
    i_fc_cell_t = zeros(n)
    for i in 1:n  # Creation of i_fc_cell_t
        i_fc_cell_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
    end
    h = ax.plot(t, i_fc_cell_t; color=colors(0))[1]
    push!(handles, h)
    push!(labels, raw"$\mathregular{i_{fc,cell}}$")
    for i in 1:nb_gc
        h = ax.plot(t, i_fc_t[i]; color=colors(i))[1]
        push!(handles, h)
        push!(labels, "\$\\mathregular{i_{fc,$(i)}}\$")
    end
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$";
                  labelpad=3)
    ax.legend(handles=handles, labels=labels, loc="best")

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_C_v_1D_temporal(variables, parameters, ax)

This function plots the vapor concentrations at different spatial localisations through the thickness of the
cell, as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the vapor concentration will be plotted.
"""
function plot_C_v_1D_temporal(variables::Dict{String, Any},
                              parameters::Dict{String, Any},
                              ax)

    # Extraction of the parameter
    nb_gc, nb_gdl, nb_mpl = parameters["nb_gc"], parameters["nb_gdl"], parameters["nb_mpl"]
    type_current, type_plot = parameters["type_current"], parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))  # Middle of the gas channel
    t = t_full[mask]
    C_v_agc_t = collect(variables["C_v_agc"][nb_gc_mid])[mask]
    C_v_agdl_t = collect(variables["C_v_agdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask]
    C_v_ampl_t = collect(variables["C_v_ampl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask]
    C_v_acl_t = collect(variables["C_v_acl"][nb_gc_mid])[mask]
    C_v_ccl_t = collect(variables["C_v_ccl"][nb_gc_mid])[mask]
    C_v_cmpl_t = collect(variables["C_v_cmpl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask]
    C_v_cgdl_t = collect(variables["C_v_cgdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask]
    C_v_cgc_t = collect(variables["C_v_cgc"][nb_gc_mid])[mask]
    T_ccl = collect(variables["T_ccl"][nb_gc_mid])[mask]

    # Plot the vapor concentrations at different spatial localisations Cv
    C_v_sat_ccl_t = [C_v_sat(T) for T in T_ccl]
    ax.plot(t, C_v_agc_t; color=colors(0))
    ax.plot(t, C_v_agdl_t; color=colors(1))
    ax.plot(t, C_v_ampl_t; color=colors(2))
    ax.plot(t, C_v_acl_t; color=colors(3))
    ax.plot(t, C_v_ccl_t; color=colors(5))
    ax.plot(t, C_v_cmpl_t; color=colors(6))
    ax.plot(t, C_v_cgdl_t; color=colors(7))
    ax.plot(t, C_v_cgc_t; color=colors(8))
    ax.plot(t, C_v_sat_ccl_t; color="k", linewidth=3)
    ax.legend([raw"$\mathregular{C_{v,agc}}$", raw"$\mathregular{C_{v,agdl}}$", raw"$\mathregular{C_{v,ampl}}$",
               raw"$\mathregular{C_{v,acl}}$", raw"$\mathregular{C_{v,ccl}}$", raw"$\mathregular{C_{v,cmpl}}$",
               raw"$\mathregular{C_{v,cgdl}}$", raw"$\mathregular{C_{v,cgc}}$", raw"$\mathregular{C_{v,sat,ccl}}$"];
              loc="best")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Vapor}$ $\mathbf{concentration}$ $\mathbf{C_{v}}$ $\mathbf{\left( mol.m^{-3} \right)}$";
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_lambda_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the water content at different spatial localisations through the thickness of the cell,
as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the water content will be plotted.
"""
function plot_lambda_1D_temporal(variables::Dict{String, Any},
                                 operating_inputs::Dict{String, Any},
                                 parameters::Dict{String, Any},
                                 ax)

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    nb_gc = parameters["nb_gc"]
    pola_current_parameters = parameters["pola_current_parameters"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))  # Middle of the gas channel
    t = t_full[mask]
    lambda_acl_t = collect(variables["lambda_acl"][nb_gc_mid])[mask]
    lambda_mem_t = collect(variables["lambda_mem"][nb_gc_mid])[mask]
    lambda_ccl_t = collect(variables["lambda_ccl"][nb_gc_mid])[mask]

    # Plot the water content at different spatial localisations: lambda
    if type_current == "polarization"
        n = length(t)
        ifc_t = zeros(n)
        for i in 1:n
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        end

        # Recovery of the internal states from the model after each stack stabilisation
        delta_t_ini_pola = pola_current_parameters["delta_t_ini_pola"]
        delta_t_load_pola = pola_current_parameters["delta_t_load_pola"]
        delta_t_break_pola = pola_current_parameters["delta_t_break_pola"]
            nb_loads = floor(Int, pola_current_parameters["i_max_pola"] / pola_current_parameters["delta_i_pola"])
        ifc_discretized_t = zeros(nb_loads)
        lambda_acl_discretized_t = zeros(nb_loads)
        lambda_mem_discretized_t = zeros(nb_loads)
        lambda_ccl_discretized_t = zeros(nb_loads)
        for i in 1:nb_loads
            t_load = delta_t_ini_pola + i * (delta_t_load_pola + delta_t_break_pola)
            idx = argmin(abs.(t .- t_load))
            ifc_discretized_t[i] = ifc_t[idx]
            lambda_acl_discretized_t[i] = lambda_acl_t[idx]
            lambda_mem_discretized_t[i] = lambda_mem_t[idx]
            lambda_ccl_discretized_t[i] = lambda_ccl_t[idx]
        end
        ax.scatter(ifc_discretized_t, lambda_acl_discretized_t; marker="o", color=colors(2))
        ax.scatter(ifc_discretized_t, lambda_mem_discretized_t; marker="o", color=colors(3))
        ax.scatter(ifc_discretized_t, lambda_ccl_discretized_t; marker="o", color=colors(4))
        ax.set_xlabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$";
                      labelpad=3)
    else
        ax.plot(t, lambda_acl_t; color=colors(3))
        ax.plot(t, lambda_mem_t; color=colors(4))
        ax.plot(t, lambda_ccl_t; color=colors(5))
        ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    end
    ax.set_ylabel(raw"$\mathbf{Water}$ $\mathbf{content}$ $\mathbf{\lambda}$"; labelpad=3)
    ax.legend([raw"$\mathregular{\lambda_{acl}}$", raw"$\mathregular{\lambda_{mem}}$",
               raw"$\mathregular{\lambda_{ccl}}$"]; loc="best")

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_s_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the liquid water saturation at different spatial localisations through the thickness of the
cell, as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the liquid water saturation will be plotted.
"""
function plot_s_1D_temporal(variables::Dict{String, Any},
                            operating_inputs::Dict{String, Any},
                            parameters::Dict{String, Any},
                            ax)

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    nb_gc = parameters["nb_gc"]
    nb_gdl = parameters["nb_gdl"]
    nb_mpl = parameters["nb_mpl"]
    pola_current_parameters = parameters["pola_current_parameters"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))
    t = t_full[mask]
    s_agc_t = collect(variables["s_agc"][nb_gc_mid])[mask]
    s_agdl_t = collect(variables["s_agdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask]
    s_ampl_t = collect(variables["s_ampl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask]
    s_acl_t = collect(variables["s_acl"][nb_gc_mid])[mask]
    s_ccl_t = collect(variables["s_ccl"][nb_gc_mid])[mask]
    s_cmpl_t = collect(variables["s_cmpl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask]
    s_cgdl_t = collect(variables["s_cgdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask]
    s_cgc_t = collect(variables["s_cgc"][nb_gc_mid])[mask]

    # Plot the liquid water saturation at different spatial localisations: s
    if type_current == "polarization"
        n = length(t)
        ifc_t = zeros(n)
        for i in 1:n
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        end
        # Recovery of the internal states from the model after each stack stabilisation
        delta_t_ini_pola = pola_current_parameters["delta_t_ini_pola"]
        delta_t_load_pola = pola_current_parameters["delta_t_load_pola"]
        delta_t_break_pola = pola_current_parameters["delta_t_break_pola"]
            nb_loads = floor(Int, pola_current_parameters["i_max_pola"] / pola_current_parameters["delta_i_pola"])
        ifc_discretized_t = zeros(nb_loads)
        s_agc_discretized_t = zeros(nb_loads)
        s_agdl_discretized_t = zeros(nb_loads)
        s_ampl_discretized_t = zeros(nb_loads)
        s_acl_discretized_t = zeros(nb_loads)
        s_ccl_discretized_t = zeros(nb_loads)
        s_cmpl_discretized_t = zeros(nb_loads)
        s_cgdl_discretized_t = zeros(nb_loads)
        s_cgc_discretized_t = zeros(nb_loads)
        for i in 1:nb_loads
            t_load = delta_t_ini_pola + i * (delta_t_load_pola + delta_t_break_pola)
            idx = argmin(abs.(t .- t_load))
            ifc_discretized_t[i] = ifc_t[idx]
            s_agc_discretized_t[i] = s_agc_t[idx]
            s_agdl_discretized_t[i] = s_agdl_t[idx]
            s_ampl_discretized_t[i] = s_ampl_t[idx]
            s_acl_discretized_t[i] = s_acl_t[idx]
            s_ccl_discretized_t[i] = s_ccl_t[idx]
            s_cmpl_discretized_t[i] = s_cmpl_t[idx]
            s_cgdl_discretized_t[i] = s_cgdl_t[idx]
            s_cgc_discretized_t[i] = s_cgc_t[idx]
        end
        ax.scatter(ifc_discretized_t, s_agc_discretized_t; marker="o", color=colors(0))
        ax.scatter(ifc_discretized_t, s_agdl_discretized_t; marker="o", color=colors(1))
        ax.scatter(ifc_discretized_t, s_ampl_discretized_t; marker="o", color=colors(2))
        ax.scatter(ifc_discretized_t, s_acl_discretized_t; marker="o", color=colors(3))
        ax.scatter(ifc_discretized_t, s_ccl_discretized_t; marker="o", color=colors(5))
        ax.scatter(ifc_discretized_t, s_cmpl_discretized_t; marker="o", color=colors(6))
        ax.scatter(ifc_discretized_t, s_cgdl_discretized_t; marker="o", color=colors(7))
        ax.scatter(ifc_discretized_t, s_cgc_discretized_t; marker="o", color=colors(8))
        ax.set_xlabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$";
                      labelpad=3)
    else
        ax.plot(t, s_agc_t; color=colors(0))
        ax.plot(t, s_agdl_t; color=colors(1))
        ax.plot(t, s_ampl_t; color=colors(2))
        ax.plot(t, s_acl_t; color=colors(3))
        ax.plot(t, s_ccl_t; color=colors(5))
        ax.plot(t, s_cmpl_t; color=colors(6))
        ax.plot(t, s_cgdl_t; color=colors(7))
        ax.plot(t, s_cgc_t; color=colors(8))
        ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    end
    ax.set_ylabel(raw"$\mathbf{Liquid}$ $\mathbf{water}$ $\mathbf{saturation}$ $\mathbf{s}$"; labelpad=3)
    ax.legend([raw"$\mathregular{s_{agc}}$", raw"$\mathregular{s_{agdl}}$", raw"$\mathregular{s_{ampl}}$",
               raw"$\mathregular{s_{acl}}$", raw"$\mathregular{s_{ccl}}$", raw"$\mathregular{s_{cmpl}}$",
               raw"$\mathregular{s_{cgdl}}$", raw"$\mathregular{s_{cgc}}$"];
              loc="best")

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_C_H2_1D_temporal(variables, parameters, ax)

This function plots the hydrogen concentration at different spatial localisations through the thickness of the
cell, as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the hydrogen concentration will be plotted.
"""
function plot_C_H2_1D_temporal(variables::Dict{String, Any},
                               parameters::Dict{String, Any},
                               ax)

    # Extraction of the parameters
    nb_gc = parameters["nb_gc"]
    nb_gdl = parameters["nb_gdl"]
    nb_mpl = parameters["nb_mpl"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))  # Middle of the gas channel
    t = t_full[mask]
    C_H2_agc_t = collect(variables["C_H2_agc"][nb_gc_mid])[mask]
    C_H2_agdl_t = collect(variables["C_H2_agdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask]
    C_H2_ampl_t = collect(variables["C_H2_ampl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask]
    C_H2_acl_t = collect(variables["C_H2_acl"][nb_gc_mid])[mask]

    # Plot the hydrogen concentration at different spatial localisations: C_H2
    ax.plot(t, C_H2_agc_t; color=colors(0))
    ax.plot(t, C_H2_agdl_t; color=colors(1))
    ax.plot(t, C_H2_ampl_t; color=colors(2))
    ax.plot(t, C_H2_acl_t; color=colors(3))
    ax.legend([raw"$\mathregular{C_{H_{2},agc}}$", raw"$\mathregular{C_{H_{2},agdl}}$", raw"$\mathregular{C_{H_{2},ampl}}$",
               raw"$\mathregular{C_{H_{2},acl}}$"];
              loc="best")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Hydrogen}$ $\mathbf{concentration}$ $\mathbf{C_{H_{2}}}$ $\mathbf{\left( mol.m^{-3} \right)}$";
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_C_O2_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the oxygen concentration at different spatial localisations through the thickness of the
cell, as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the oxygen concentration will be plotted.
"""
function plot_C_O2_1D_temporal(variables::Dict{String, Any},
                               operating_inputs::Dict{String, Any},
                               parameters::Dict{String, Any},
                               ax)

    # Extraction of the parameters
    current_density = operating_inputs["current_density"]
    Hccl = parameters["Hccl"]
    nb_gc = parameters["nb_gc"]
    nb_gdl = parameters["nb_gdl"]
    nb_mpl = parameters["nb_mpl"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))
    t = t_full[mask]
    C_O2_ccl_t = collect(variables["C_O2_ccl"][nb_gc_mid])[mask]
    C_O2_cmpl_t = collect(variables["C_O2_cmpl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask]
    C_O2_cgdl_t = collect(variables["C_O2_cgdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask]
    C_O2_cgc_t = collect(variables["C_O2_cgc"][nb_gc_mid])[mask]
    C_O2_Pt_t = collect(variables["C_O2_Pt"][nb_gc_mid])[mask]

    # Plot the oxygen concentration at different spatial localisations: C_O2
    ax.plot(t, C_O2_Pt_t; color=colors(10))
    ax.plot(t, C_O2_ccl_t; color=colors(5))
    ax.plot(t, C_O2_cmpl_t; color=colors(6))
    ax.plot(t, C_O2_cgdl_t; color=colors(7))
    ax.plot(t, C_O2_cgc_t; color=colors(8))
    ax.legend([raw"$\mathregular{C_{O_{2},P_t}}$", raw"$\mathregular{C_{O_{2},ccl}}$", raw"$\mathregular{C_{O_{2},cmpl}}$",
               raw"$\mathregular{C_{O_{2},cgdl}}$", raw"$\mathregular{C_{O_{2},cgc}}$"];
              loc="best")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Oxygen}$ $\mathbf{concentration}$ $\mathbf{C_{O_{2}}}$ $\mathbf{\left( mol.m^{-3} \right)}$";
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_C_N2_1D_temporal(variables, parameters, ax)

This function plots the nitrogen concentration as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the nitrogen concentration will be plotted.
"""
function plot_C_N2_1D_temporal(variables::Dict{String, Any},
                               parameters::Dict{String, Any},
                               ax)

    # Extraction of the parameters
    nb_gc = parameters["nb_gc"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))
    t = t_full[mask]
    C_N2_agc_t = collect(variables["C_N2_agc"][nb_gc_mid])[mask]
    C_N2_cgc_t = collect(variables["C_N2_cgc"][nb_gc_mid])[mask]

    # Plot C_N2
    ax.plot(t, C_N2_agc_t; color=colors(6))
    ax.plot(t, C_N2_cgc_t; color=colors(6))
    ax.legend([raw"$\mathregular{C_{N_{2},agc}}$", raw"$\mathregular{C_{N_{2},cgc}}$"];
              loc="best")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Nitrogen}$ $\mathbf{concentration}$ $\mathbf{C_{N_{2}}}$ $\mathbf{\left( mol.m^{-3} \right)}$";
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_P_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the pressure at different spatial localisations as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the pressure will be plotted.
"""
function plot_P_1D_temporal(variables::Dict{String, Any},
                            operating_inputs::Dict{String, Any},
                            parameters::Dict{String, Any},
                            ax)

    # Extraction of the parameters
    Pa_des = operating_inputs["Pa_des"]
    Pc_des = operating_inputs["Pc_des"]
    nb_gc = parameters["nb_gc"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    t = t_full[mask]
    nb_gc_mid = Int(ceil(nb_gc / 2))
    C_v_agc = collect(variables["C_v_agc"][nb_gc_mid])[mask]
    C_H2_agc = collect(variables["C_H2_agc"][nb_gc_mid])[mask]
    C_N2_agc = collect(variables["C_N2_agc"][nb_gc_mid])[mask]
    T_agc = collect(variables["T_agc"][nb_gc_mid])[mask]
    C_v_cgc = collect(variables["C_v_cgc"][nb_gc_mid])[mask]
    C_O2_cgc = collect(variables["C_O2_cgc"][nb_gc_mid])[mask]
    C_N2_cgc = collect(variables["C_N2_cgc"][nb_gc_mid])[mask]
    T_cgc = collect(variables["T_cgc"][nb_gc_mid])[mask]
    Pagc_t = (C_v_agc .+ C_H2_agc .+ C_N2_agc) .* R .* T_agc ./ 1e5  # Conversion in bar
    Pcgc_t = (C_v_cgc .+ C_O2_cgc .+ C_N2_cgc) .* R .* T_cgc ./ 1e5  # Conversion in bar
    if parameters["type_auxiliary"] != "no_auxiliary"
        Pasm_t = collect(variables["Pasm"])[mask] ./ 1e5
        Paem_t = collect(variables["Paem"])[mask] ./ 1e5
        Pcsm_t = collect(variables["Pcsm"])[mask] ./ 1e5
        Pcem_t = collect(variables["Pcem"])[mask] ./ 1e5
    else
        Pa_in_t = collect(variables["Pa_in"])[mask] ./ 1e5
        Pa_out_t = fill(Pa_des, length(t)) ./ 1e5
        Pc_in_t = collect(variables["Pc_in"])[mask] ./ 1e5
        Pc_out_t = fill(Pc_des, length(t)) ./ 1e5
    end

    # Plot the pressure at different spatial localisations: P
    ax.plot(t, Pagc_t; color=colors(0))
    ax.plot(t, Pcgc_t; color=colors(6))
    if parameters["type_auxiliary"] != "no_auxiliary"
        ax.plot(t, Pasm_t; color=colors(7))
        ax.plot(t, Paem_t; color=colors(8))
        ax.plot(t, Pcsm_t; color=colors(9))
        ax.plot(t, Pcem_t; color=colors(10))
    else
        ax.plot(t, Pa_in_t; color=colors(7))
        ax.plot(t, Pa_out_t; color=colors(8))
        ax.plot(t, Pc_in_t; color=colors(9))
        ax.plot(t, Pc_out_t; color=colors(10))
    end
    if parameters["type_auxiliary"] != "no_auxiliary"
        ax.legend([raw"$\mathregular{P_{agc}}$", raw"$\mathregular{P_{cgc}}$", raw"$\mathregular{P_{asm}}$",
                   raw"$\mathregular{P_{aem}}$", raw"$\mathregular{P_{csm}}$", raw"$\mathregular{P_{cem}}$"];
                  loc="best")
    else
        ax.legend([raw"$\mathregular{P_{agc}}$", raw"$\mathregular{P_{cgc}}$", raw"$\mathregular{P_{a,in}}$",
                   raw"$\mathregular{P_{a,out}}$", raw"$\mathregular{P_{c,in}}$", raw"$\mathregular{P_{c,out}}$"];
                  loc="best")
    end
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Pressure}$ $\mathbf{P}$ $\mathbf{\left( bar \right)}$"; labelpad=3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_T_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the temperature at different spatial localisations through the thickness of the
cell, as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the temperature will be plotted.
"""
function plot_T_1D_temporal(variables::Dict{String, Any},
                            operating_inputs::Dict{String, Any},
                            parameters::Dict{String, Any},
                            ax)

    # Extraction of the operating inputs and parameters
    current_density = operating_inputs["current_density"]
    T_des = operating_inputs["T_des"]
    nb_gc = parameters["nb_gc"]
    nb_gdl = parameters["nb_gdl"]
    nb_mpl = parameters["nb_mpl"]
    pola_current_parameters = parameters["pola_current_parameters"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables and the operating inputs
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))  # Middle of the gas channel
    t = t_full[mask]
    T_agc_t = collect(variables["T_agc"][nb_gc_mid])[mask] .- 273.15
    T_agdl_t = collect(variables["T_agdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask] .- 273.15
    T_ampl_t = collect(variables["T_ampl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask] .- 273.15
    T_acl_t = collect(variables["T_acl"][nb_gc_mid])[mask] .- 273.15
    T_mem_t = collect(variables["T_mem"][nb_gc_mid])[mask] .- 273.15
    T_ccl_t = collect(variables["T_ccl"][nb_gc_mid])[mask] .- 273.15
    T_cmpl_t = collect(variables["T_cmpl_$(Int(ceil(nb_mpl / 2)))"][nb_gc_mid])[mask] .- 273.15
    T_cgdl_t = collect(variables["T_cgdl_$(Int(ceil(nb_gdl / 2)))"][nb_gc_mid])[mask] .- 273.15
    T_cgc_t = collect(variables["T_cgc"][nb_gc_mid])[mask] .- 273.15

    # Plot the temperature at different spatial localisations
    if type_current == "polarization"
        n = length(t)
        ifc_t = zeros(n)
        for i in 1:n
            ifc_t[i] = current_density(t[i], parameters) / 1e4  # Conversion in A/cm²
        end
        # Recovery of the internal states from the model after each stack stabilisation
        delta_t_ini_pola = pola_current_parameters["delta_t_ini_pola"]
        delta_t_load_pola = pola_current_parameters["delta_t_load_pola"]
        delta_t_break_pola = pola_current_parameters["delta_t_break_pola"]
            nb_loads = floor(Int, pola_current_parameters["i_max_pola"] / pola_current_parameters["delta_i_pola"])
        ifc_discretized_t = zeros(nb_loads)
        T_agc_discretized_t = zeros(nb_loads)
        T_agdl_discretized_t = zeros(nb_loads)
        T_ampl_discretized_t = zeros(nb_loads)
        T_acl_discretized_t = zeros(nb_loads)
        T_mem_discretized_t = zeros(nb_loads)
        T_ccl_discretized_t = zeros(nb_loads)
        T_cmpl_discretized_t = zeros(nb_loads)
        T_cgdl_discretized_t = zeros(nb_loads)
        T_cgc_discretized_t = zeros(nb_loads)
        for i in 1:nb_loads
            t_load = delta_t_ini_pola + i * (delta_t_load_pola + delta_t_break_pola)
            idx = argmin(abs.(t .- t_load))
            ifc_discretized_t[i] = ifc_t[idx]
            T_agc_discretized_t[i] = T_agc_t[idx]
            T_agdl_discretized_t[i] = T_agdl_t[idx]
            T_ampl_discretized_t[i] = T_ampl_t[idx]
            T_acl_discretized_t[i] = T_acl_t[idx]
            T_mem_discretized_t[i] = T_mem_t[idx]
            T_ccl_discretized_t[i] = T_ccl_t[idx]
            T_cmpl_discretized_t[i] = T_cmpl_t[idx]
            T_cgdl_discretized_t[i] = T_cgdl_t[idx]
            T_cgc_discretized_t[i] = T_cgc_t[idx]
        end
        T_des_discretized_t = fill(T_des - 273.15, length(ifc_discretized_t))
        ax.scatter(ifc_discretized_t, T_agc_discretized_t; marker="o", color=colors(0))
        ax.scatter(ifc_discretized_t, T_agdl_discretized_t; marker="o", color=colors(1))
        ax.scatter(ifc_discretized_t, T_ampl_discretized_t; marker="o", color=colors(2))
        ax.scatter(ifc_discretized_t, T_acl_discretized_t; marker="o", color=colors(3))
        ax.scatter(ifc_discretized_t, T_mem_discretized_t; marker="o", color=colors(4))
        ax.scatter(ifc_discretized_t, T_ccl_discretized_t; marker="o", color=colors(5))
        ax.scatter(ifc_discretized_t, T_cmpl_discretized_t; marker="o", color=colors(6))
        ax.scatter(ifc_discretized_t, T_cgdl_discretized_t; marker="o", color=colors(7))
        ax.scatter(ifc_discretized_t, T_cgc_discretized_t; marker="o", color=colors(8))
        ax.scatter(ifc_discretized_t, T_des_discretized_t; marker="o", color="k")
        ax.set_xlabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$";
                      labelpad=3)
    else
        T_des_t = fill(T_des - 273.15, length(t))
        ax.plot(t, T_agc_t; color=colors(0))
        ax.plot(t, T_agdl_t; color=colors(1))
        ax.plot(t, T_ampl_t; color=colors(2))
        ax.plot(t, T_acl_t; color=colors(3))
        ax.plot(t, T_mem_t; color=colors(4))
        ax.plot(t, T_ccl_t; color=colors(5))
        ax.plot(t, T_cmpl_t; color=colors(6))
        ax.plot(t, T_cgdl_t; color=colors(7))
        ax.plot(t, T_cgc_t; color=colors(8))
        ax.plot(t, T_des_t; color="k", linewidth=3)
        ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    end
    ax.legend([raw"$\mathregular{T_{agc}}$", raw"$\mathregular{T_{agdl}}$", raw"$\mathregular{T_{ampl}}$",
               raw"$\mathregular{T_{acl}}$", raw"$\mathregular{T_{mem}}$", raw"$\mathregular{T_{ccl}}$",
               raw"$\mathregular{T_{cmpl}}$", raw"$\mathregular{T_{cgdl}}$", raw"$\mathregular{T_{cgc}}$",
               raw"$\mathregular{T_{des}}$"]; loc="best")
    ax.set_ylabel(raw"$\mathbf{Temperature}$ $\mathbf{T}$ $\mathbf{\left( ^\circ C \right)}$"; labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_Ucell(variables, parameters, ax)

This function plots the cell voltage as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the cell voltage will be plotted.
"""
function plot_Ucell(variables::Dict{String, Any},
                    parameters::Dict{String, Any},
                    ax)

    # Extraction of the parameters
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    t = t_full[mask]
    Ucell_t = collect(variables["Ucell"])[mask]

    # Plot the cell voltage: Ucell
    ax.plot(t, Ucell_t; color=colors(0), label=raw"$\mathregular{U_{cell}}$")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Cell}$ $\mathbf{voltage}$ $\mathbf{U_{cell}}$ $\mathbf{\left( V \right)}$"; labelpad=3)
    ax.legend([raw"$\mathregular{U_{cell}}$"]; loc="best")

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_Phi_a_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the humidity at the anode side, at different spatial localisations, as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the humidity will be plotted.
"""
function plot_Phi_a_1D_temporal(variables::Dict{String, Any},
                                operating_inputs::Dict{String, Any},
                                parameters::Dict{String, Any},
                                ax)
    # Extraction of the operating inputs and parameters
    Phi_a_des = operating_inputs["Phi_a_des"]
    nb_gc = parameters["nb_gc"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))
    t = t_full[mask]
    C_v_agc_t = collect(variables["C_v_agc"][nb_gc_mid])[mask]
    T_agc_t = collect(variables["T_agc"][nb_gc_mid])[mask]
    Phi_asm_t = collect(variables["Phi_asm"])[mask]
    Phi_aem_t = collect(variables["Phi_aem"])[mask]

    # Calculate the humidity Phi
    Phi_agc_t = C_v_agc_t .* R .* T_agc_t ./ Psat.(T_agc_t)

    # Plot the humidity at different spatial localisations: Phi
    ax.plot(t, Phi_agc_t; color=colors(0), label=raw"$\mathregular{\Phi_{agc}}$")
    ax.plot(t, Phi_asm_t; color=colors(1), label=raw"$\mathregular{\Phi_{asm}}$")
    ax.plot(t, Phi_aem_t; color=colors(2), label=raw"$\mathregular{\Phi_{aem}}$")
    ax.plot(t, fill(Phi_a_des, length(t)); color="black", label=raw"$\mathregular{\Phi_{a,des}}$")
    ax.legend(loc="center right", bbox_to_anchor=(1, 0.67))
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Humidity}$ $\mathbf{at}$ $\mathbf{the}$ $\mathbf{anode}$ $\mathbf{side}$ $\mathbf{\Phi}$";
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_Phi_c_1D_temporal(variables, operating_inputs, parameters, ax)

This function plots the humidity at the cathode side, at different spatial localisations, as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `operating_inputs::Dict{String, Any}`: Operating inputs of the fuel cell.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the humidity will be plotted.
"""
function plot_Phi_c_1D_temporal(variables::Dict{String, Any},
                                operating_inputs::Dict{String, Any},
                                parameters::Dict{String, Any},
                                ax)
    # Extraction of the operating inputs and parameters
    Phi_c_des = operating_inputs["Phi_c_des"]
    nb_gc = parameters["nb_gc"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))
    t = t_full[mask]
    C_v_cgc_t = collect(variables["C_v_cgc"][nb_gc_mid])[mask]
    T_cgc_t = collect(variables["T_cgc"][nb_gc_mid])[mask]
    Phi_csm_t = collect(variables["Phi_csm"])[mask]
    Phi_cem_t = collect(variables["Phi_cem"])[mask]

    # Calculate the humidity Phi
    Phi_cgc_t = C_v_cgc_t .* R .* T_cgc_t ./ Psat.(T_cgc_t)

    # Plot the humidity at different spatial localisations: Phi
    ax.plot(t, Phi_cgc_t; color=colors(0), label=raw"$\mathregular{\Phi_{cgc}}$")
    ax.plot(t, Phi_csm_t; color=colors(1), label=raw"$\mathregular{\Phi_{csm}}$")
    ax.plot(t, Phi_cem_t; color=colors(2), label=raw"$\mathregular{\Phi_{cem}}$")
    ax.plot(t, fill(Phi_c_des, length(t)); color="black", label=raw"$\mathregular{\Phi_{c,des}}$")
    ax.legend(loc="best")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Humidity}$ $\mathbf{at}$ $\mathbf{the}$ $\mathbf{cathode}$ $\mathbf{side}$ $\mathbf{\Phi}$";
                  labelpad=3)

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_v_1D_temporal(variables, parameters, ax)

This function plots the velocity at the anode and the cathode as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the pressure will be plotted.
"""
function plot_v_1D_temporal(variables::Dict{String, Any},
                            parameters::Dict{String, Any},
                            ax)

    # Extraction of the parameters
    nb_gc = parameters["nb_gc"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2)) # Middle of the gas channel
    t = t_full[mask]
    v_a_t = collect(variables["v_a"][nb_gc_mid])[mask]
    v_c_t = collect(variables["v_c"][nb_gc_mid])[mask]

    # Plot the pressure at different spatial localisations: P
    ax.plot(t, v_a_t; color=colors(0))
    ax.plot(t, v_c_t; color=colors(6))
    ax.legend([raw"$\mathregular{v_{a}}$", raw"$\mathregular{v_{c}}$"]; loc="best")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Velocities at the middle of the GC}$ $\mathbf{v}$ $\mathbf{\left( m.s^{-1} \right)}$";
                  labelpad=3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_Re_nb_1D_temporal(variables, parameters, ax)

This function plots the Reynolds number at the center of the AGC and the CGC as a function of time.

# Arguments
- `variables::Dict{String, Any}`: Variables calculated by the solver.
- `parameters::Dict{String, Any}`: Parameters of the fuel cell model.
- `ax`: Axes on which the pressure will be plotted.
"""
function plot_Re_nb_1D_temporal(variables::Dict{String, Any},
                                parameters::Dict{String, Any},
                                ax)

    # Extraction of the parameters
    Hcgc = parameters["Hcgc"]
    Wcgc = parameters["Wcgc"]
    nb_gc = parameters["nb_gc"]
    type_current = parameters["type_current"]
    type_plot = parameters["type_plot"]
    if type_current == "step"
        delta_t_ini = parameters["step_current_parameters"]["delta_t_ini_step"]
    elseif type_current == "polarization"
        delta_t_ini = parameters["pola_current_parameters"]["delta_t_ini_pola"]
    elseif type_current == "polarization_for_cali"
        delta_t_ini = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]
    else
        delta_t_ini = 0
    end

    # Extraction of the variables
    t_full = collect(variables["t"])
    mask = if type_plot == "fixed"
        t_full .>= 0.9 * delta_t_ini
    else
        trues(length(t_full))
    end
    nb_gc_mid = Int(ceil(nb_gc / 2))
    t = t_full[mask]
    v_a_t = collect(variables["v_a"][nb_gc_mid])[mask]
    v_c_t = collect(variables["v_c"][nb_gc_mid])[mask]
    C_v_agc_t = collect(variables["C_v_agc"][nb_gc_mid])[mask]
    C_v_cgc_t = collect(variables["C_v_cgc"][nb_gc_mid])[mask]
    C_H2_agc_t = collect(variables["C_H2_agc"][nb_gc_mid])[mask]
    C_O2_cgc_t = collect(variables["C_O2_cgc"][nb_gc_mid])[mask]
    C_N2_agc_t = collect(variables["C_N2_agc"][nb_gc_mid])[mask]
    C_N2_cgc_t = collect(variables["C_N2_cgc"][nb_gc_mid])[mask]
    T_agc_t = collect(variables["T_agc"][nb_gc_mid])[mask]
    T_cgc_t = collect(variables["T_cgc"][nb_gc_mid])[mask]

    # Calculation of the Reynolds Number
    d_pipe = sqrt(4 * Hcgc * Wcgc / π)
    P_agc = (C_v_agc_t .+ C_H2_agc_t .+ C_N2_agc_t) .* R .* T_agc_t
    P_cgc = (C_v_cgc_t .+ C_O2_cgc_t .+ C_N2_cgc_t) .* R .* T_cgc_t
    M_agc = C_v_agc_t .* R .* T_agc_t ./ P_agc .* M_H2O .+
            C_H2_agc_t .* R .* T_agc_t ./ P_agc .* M_H2 .+
            C_N2_agc_t .* R .* T_agc_t ./ P_agc .* M_N2
    M_cgc = C_v_cgc_t .* R .* T_cgc_t ./ P_cgc .* M_H2O .+
            C_O2_cgc_t .* R .* T_cgc_t ./ P_cgc .* M_O2 .+
            C_N2_cgc_t .* R .* T_cgc_t ./ P_cgc .* M_N2
    rho_agc = P_agc ./ (R .* T_agc_t) .* M_agc
    rho_cgc = P_cgc ./ (R .* T_cgc_t) .* M_cgc
    x_H2O_v_agc = C_v_agc_t ./ (C_v_agc_t .+ C_H2_agc_t .+ C_N2_agc_t)
    x_H2O_v_cgc = C_v_cgc_t ./ (C_v_cgc_t .+ C_O2_cgc_t .+ C_N2_cgc_t)
    y_H2_agc = C_H2_agc_t ./ (C_H2_agc_t .+ C_N2_agc_t)
    y_O2_cgc = C_O2_cgc_t ./ (C_O2_cgc_t .+ C_N2_cgc_t)
    mu_agc = [mu_mixture_gases(["H2O_v", "H2", "N2"],
                               [x_H2O_v_agc[i], y_H2_agc[i] * (1 - x_H2O_v_agc[i]), (1 - y_H2_agc[i]) * (1 - x_H2O_v_agc[i])],
                               T_agc_t[i]) for i in eachindex(T_agc_t)]
    mu_cgc = [mu_mixture_gases(["H2O_v", "O2", "N2"],
                               [x_H2O_v_cgc[i], y_O2_cgc[i] * (1 - x_H2O_v_cgc[i]), (1 - y_O2_cgc[i]) * (1 - x_H2O_v_cgc[i])],
                               T_cgc_t[i]) for i in eachindex(T_cgc_t)]
    Re_nb_a_t = (rho_agc .* v_a_t .* d_pipe) ./ mu_agc
    Re_nb_c_t = (rho_cgc .* v_c_t .* d_pipe) ./ mu_cgc

    # Plot the pressure at different spatial localisations: P
    ax.plot(t, Re_nb_a_t; color=colors(0))
    ax.plot(t, Re_nb_c_t; color=colors(6))
    ax.legend([raw"$\mathregular{Re_{a}}$", raw"$\mathregular{Re_{c}}$"]; loc="best")
    ax.set_xlabel(raw"$\mathbf{Time}$ $\mathbf{t}$ $\mathbf{\left( s \right)}$"; labelpad=3)
    ax.set_ylabel(raw"$\mathbf{Reynold}$ $\mathbf{number}$ $\mathbf{at}$ $\mathbf{the}$ $\mathbf{inlet,}$ $\mathbf{Re}$";
                  labelpad=3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


# ____________________________________________Internal variables - 1D+1D final__________________________________________

"""
    plot_T_pseudo_2D_final(variables, operating_inputs, parameters, ax)

This function plots the temperature at different spatial localisations inside the fuel cell in pseudo 2D,
for the last time step of the simulation.
"""
function plot_T_pseudo_2D_final(variables::Dict{String, Any},
                                operating_inputs::Dict{String, Any},
                                parameters::Dict{String, Any},
                                ax)

    # Extraction of the operating inputs and parameters
    T_des = operating_inputs["T_des"] - 273.15  # Conversion in °C.
    nb_gc = parameters["nb_gc"]
    nb_gdl = parameters["nb_gdl"]
    nb_mpl = parameters["nb_mpl"]

    # Construction of the temperature matrix for the pseudo-2D plot
    var_order = vcat(["agc"],
                     ["agdl$(j)" for j in 1:nb_gdl],
                     ["ampl$(j)" for j in 1:nb_mpl],
                     ["acl", "mem", "ccl"],
                     ["cmpl$(j)" for j in 1:nb_mpl],
                     ["cgdl$(j)" for j in 1:nb_gdl],
                     ["cgc"])
    n_cols = length(var_order)
    temp_matrix = fill(NaN, nb_gc, n_cols)

    # Helper to fetch the last element of a variable and convert to °C
    last_celsius(arr) = collect(arr)[end] - 273.15

    # Extraction of the last temperature values and insertion into the matrix
    col = 1
    for name in var_order
        if name == "agc"
            temp_matrix[:, col] = [last_celsius(variables["T_agc"][i]) for i in 1:nb_gc]
        elseif startswith(name, "agdl")
            j = parse(Int, replace(name, "agdl" => ""))
            temp_matrix[:, col] = [last_celsius(variables["T_agdl_$(j)"][i]) for i in 1:nb_gc]
        elseif startswith(name, "ampl")
            j = parse(Int, replace(name, "ampl" => ""))
            temp_matrix[:, col] = [last_celsius(variables["T_ampl_$(j)"][i]) for i in 1:nb_gc]
        elseif name in ("acl", "mem", "ccl", "cgc")
            temp_matrix[:, col] = [last_celsius(variables["T_$(name)"][i]) for i in 1:nb_gc]
        elseif startswith(name, "cmpl")
            j = parse(Int, replace(name, "cmpl" => ""))
            temp_matrix[:, col] = [last_celsius(variables["T_cmpl_$(j)"][i]) for i in 1:nb_gc]
        elseif startswith(name, "cgdl")
            j = parse(Int, replace(name, "cgdl" => ""))
            temp_matrix[:, col] = [last_celsius(variables["T_cgdl_$(j)"][i]) for i in 1:nb_gc]
        end
        col += 1
    end

    # Plot the figure
    vmin = np.nanmin(temp_matrix)
    vmax = np.nanmax(temp_matrix)
    cmap = plt.get_cmap("Reds")
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(temp_matrix; aspect="auto", origin="lower", cmap=cmap, norm=norm)

    # X-axis labels
    x_labels = vcat(["agc"],
                    ["\$\\mathregular{agdl_{$j}}\$" for j in 1:nb_gdl],
                    ["\$\\mathregular{ampl_{$j}}\$" for j in 1:nb_mpl],
                    ["acl", "mem", "ccl"],
                    ["\$\\mathregular{cmpl_{$j}}\$" for j in 1:nb_mpl],
                    ["\$\\mathregular{cgdl_{$j}}\$" for j in 1:nb_gdl],
                    ["cgc"])
    ax.set_xticks(collect(0:n_cols-1))
    ax.set_xticklabels(x_labels; rotation=45, ha="right", fontsize=9)

    # Y-axis labels
    y_labels = [string(i) for i in 1:nb_gc]
    if nb_gc > 1
        y_labels[1] = "inlet - 1"
        y_labels[end] = "outlet - $(nb_gc)"
    end
    ax.set_yticks(collect(0:nb_gc-1))
    ax.set_yticklabels(y_labels; fontsize=9)
    ax.invert_yaxis()  # To have the first gas channel at the top

    ax.set_xlabel("Through the thickness of the cell"; fontsize=11)
    ax.set_ylabel("Through the gas channel"; fontsize=11)
    ax.set_title(raw"$\mathbf{Final\ temperature\ distribution\ \left(^\circ\!C\right)}$"; fontsize=13, pad=15)

    cbar = plt.colorbar(im; ax=ax, orientation="vertical")
    cbar.set_label(raw"$\mathbf{Temperature\ \left(^\circ\!C\right)}$"; rotation=270, labelpad=15)
    ax.text(0.95, 0.05,
            "\$T_{des} = $(round(T_des, digits=1))\\,^\\circ C\$";
            transform=ax.transAxes, ha="right", va="bottom", fontsize=10,
            bbox=Dict{String, Any}("boxstyle"=>"round,pad=0.3", "fc"=>"white", "ec"=>"gray", "lw"=>0.5))
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

    return nothing
end


# ___________________________________________________Global indicators__________________________________________________

"""
    plot_power_density_curve(variables, operating_inputs, parameters, n, ax)

This function plots the power density curve Pfc, produced by a cell, as a function of the current density.
"""
function plot_power_density_curve(variables::Dict{String, Any},
                                  operating_inputs::Dict{String, Any},
                                  parameters::Dict{String, Any},
                                  n::Integer,
                                  ax)

    # Extraction of the variables
    t = variables["t"]
    Ucell_t = variables["Ucell"]
    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    type_fuel_cell = parameters["type_fuel_cell"]
    type_current = parameters["type_current"]
    type_auxiliary = parameters["type_auxiliary"]

    # Creation of the power density function: Pfc
    ifc_t = zeros(n)
    Pfc_t = zeros(n)
    for i in 1:n
        ifc_t[i] = current_density(t[i], parameters) / 1e4
        Pfc_t[i] = Ucell_t[i] * ifc_t[i]
    end

    # Plot of the power density function: Pfc
    plot_specific_line(ifc_t, Pfc_t, type_fuel_cell, type_current, type_auxiliary, nothing, ax)
    ax.set_xlabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$"; labelpad=0)
    ax.set_ylabel(raw"$\mathbf{Fuel}$ $\mathbf{cell}$ $\mathbf{power}$ $\mathbf{density}$ $\mathbf{P_{fc}}$ $\mathbf{\left( W.cm^{-2} \right)}$";
                  labelpad=0)
    ax.legend(loc="best")

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


"""
    plot_cell_efficiency(variables, operating_inputs, parameters, n, ax)

This function plots the fuel cell efficiency eta_fc as a function of the current density.
"""
function plot_cell_efficiency(variables::Dict{String, Any},
                              operating_inputs::Dict{String, Any},
                              parameters::Dict{String, Any},
                              n::Integer,
                              ax)

    # Extraction of the operating inputs and the parameters
    current_density = operating_inputs["current_density"]
    Hmem = parameters["Hmem"]
    Hacl = parameters["Hacl"]
    Hccl = parameters["Hccl"]
    kappa_co = parameters["kappa_co"]
    nb_gc = parameters["nb_gc"]
    type_fuel_cell = parameters["type_fuel_cell"]
    type_current = parameters["type_current"]
    type_auxiliary = parameters["type_auxiliary"]

    # Extraction of the variables
    nb_gc_mid = Int(ceil(nb_gc / 2))
    t = variables["t"]
    Ucell_t = variables["Ucell"]
    lambda_mem_t = variables["lambda_mem"][nb_gc_mid]
    C_H2_acl_t = variables["C_H2_acl"][nb_gc_mid]
    C_O2_ccl_t = variables["C_O2_ccl"][nb_gc_mid]
    T_acl_t = variables["T_acl"][nb_gc_mid]
    T_mem_t = variables["T_mem"][nb_gc_mid]
    T_ccl_t = variables["T_ccl"][nb_gc_mid]

    # Creation of the fuel cell efficiency: eta_fc
    ifc_t = zeros(n)
    Pfc_t = zeros(n)
    eta_fc_t = zeros(n)
    for i in 1:n
        ifc_t[i] = current_density(t[i], parameters) / 1e4
        Pfc_t[i] = Ucell_t[i] * ifc_t[i]
        Ueq = E0 - 8.5e-4 * (T_ccl_t[i] - 298.15) +
              R * T_ccl_t[i] / (2 * F) * (log(R * T_acl_t[i] * C_H2_acl_t[i] / Pref_eq) +
                                          0.5 * log(R * T_ccl_t[i] * C_O2_ccl_t[i] / Pref_eq))
        T_acl_mem_ccl = average([T_acl_t[i], T_mem_t[i], T_ccl_t[i]],
                                [Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
        i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl_t[i] * k_H2(lambda_mem_t[i], T_mem_t[i], kappa_co)
        i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl_t[i] * k_O2(lambda_mem_t[i], T_mem_t[i], kappa_co)
        i_n = (i_H2 + i_O2) / 1e4
        eta_fc_t[i] = Pfc_t[i] / (Ueq * (ifc_t[i] + i_n))
    end

    # Plot of the fuel cell efficiency: eta_fc
    plot_specific_line(ifc_t, eta_fc_t, type_fuel_cell, type_current, type_auxiliary, nothing, ax)
    ax.set_xlabel(raw"$\mathbf{Current}$ $\mathbf{density}$ $\mathbf{i_{fc}}$ $\mathbf{\left( A.cm^{-2} \right)}$"; labelpad=0)
    ax.set_ylabel(raw"$\mathbf{Fuel}$ $\mathbf{cell}$ $\mathbf{efficiency}$ $\mathbf{\eta_{fc}}$"; labelpad=0)
    ax.legend(loc="best")

    # Plot instructions
    plot_general_instructions(ax)
    return nothing
end


# _________________________________________________Plot instructions____________________________________________________

"""
    plot_general_instructions(ax, true)

This function adds the common instructions for all the plots displayed by AlphaPEM to the ax object.

# Arguments
- `ax`: Axes on which the instructions will be added.
- `set_y::Bool=true`: If true, set y-axis major/minor locators.
"""
function plot_general_instructions(ax, set_y::Bool=true)
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
    if set_y
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(major_step_y_rounded))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(major_step_y_rounded / 5))
    end
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis="both", which="major", size=10, width=1.5, direction="out")
    ax.tick_params(axis="both", which="minor", size=5, width=1.5, direction="out")
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
    # Adjust layout to prevent overlap between labels and the figure
    plt.show()  # Show the figure

    return nothing
end


"""
    plot_pola_instructions(type_fuel_cell, ax, true)

This function adds the specific instructions for polarisation plots according to the type_input to the ax object.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `ax`: Axes on which the instructions will be added.
- `show::Bool=true`: If true, the figure will be displayed.
"""
function plot_pola_instructions(type_fuel_cell::String, ax, show::Bool=true)

    # For ZSW-GenStack fuel cell
    if type_fuel_cell == "ZSW-GenStack" || type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" ||
       type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" || type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" ||
       type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" || type_fuel_cell == "ZSW-GenStack_T_62" ||
       type_fuel_cell == "ZSW-GenStack_T_76" || type_fuel_cell == "ZSW-GenStack_T_84"
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1 / 5))
        ax.set_xlim(0, 3.0)
        ax.set_ylim(0.4, 1.24)

    # For EH-31 fuel cell
    elseif type_fuel_cell == "EH-31_1.5" || type_fuel_cell == "EH-31_2.0" ||
           type_fuel_cell == "EH-31_2.25" || type_fuel_cell == "EH-31_2.5"
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1 / 5))
        ax.set_xlim(0, 3.0)
        ax.set_ylim(0.4, 1.04)
    else
        nothing
    end

    # Configure the appearance of major and minor ticks
    ax.tick_params(axis="both", which="major", size=10, width=1.5, direction="out")
    ax.tick_params(axis="both", which="minor", size=5, width=1.5, direction="out")
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
    # Adjust layout to prevent overlap between labels and the figure
    if show
        plt.show()  # Show the figure
    end

    return nothing
end


"""
    plot_specific_line(x, y, type_fuel_cell, type_current, type_auxiliary, sim_error, ax)

This function adds the appropriate plot configuration according to the type_input to the ax object.
"""
function plot_specific_line(x,
                            y,
                            type_fuel_cell::String,
                            type_current::String,
                            type_auxiliary::String,
                            sim_error,
                            ax)

    aux_ok = (type_auxiliary == "forced-convective_cathode_with_flow-through_anode" || type_auxiliary == "no_auxiliary")

    if type_current == "polarization"
        if type_fuel_cell == "ZSW-GenStack" && aux_ok
            ax.plot(x, y, "--"; color=colors(0), label="Sim. - nominal - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack"
            ax.plot(x, y; color=colors(0), label="Sim. - nominal")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" && aux_ok
            ax.plot(x, y, "--"; color=colors(1), label="Sim. - P\$_a\$/P\$_c\$ = 1.61/1.41 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
            ax.plot(x, y; color=colors(1), label="Sim. - P\$_a\$/P\$_c\$ = 1.61/1.41 bar")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" && aux_ok
            ax.plot(x, y; color=colors(2), label="Sim. - P\$_a\$/P\$_c\$ = 2.01/1.81 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
            ax.plot(x, y; color=colors(2), label="Sim. - P\$_a\$/P\$_c\$ = 2.01/1.81 bar")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" && aux_ok
            ax.plot(x, y; color=colors(3), label="Sim. - P\$_a\$/P\$_c\$ = 2.4/2.2 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
            ax.plot(x, y; color=colors(3), label="Sim. - P\$_a\$/P\$_c\$ = 2.4/2.2 bar")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" && aux_ok
            ax.plot(x, y; color=colors(4), label="Sim. - P\$_a\$/P\$_c\$ = 2.8/2.6 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
            ax.plot(x, y; color=colors(4), label="Sim. - P\$_a\$/P\$_c\$ = 2.8/2.6 bar")

        elseif type_fuel_cell == "ZSW-GenStack_T_62" && aux_ok
            ax.plot(x, y, "--"; color=colors(5), label="Sim. - T = 62 \$^\\circ\$C - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_T_62"
            ax.plot(x, y; color=colors(5), label="Sim. - T = 62 \$^\\circ\$C")

        elseif type_fuel_cell == "ZSW-GenStack_T_76" && aux_ok
            ax.plot(x, y; color=colors(6), label="Sim. - T = 76 \$^\\circ\$C - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_T_76"
            ax.plot(x, y; color=colors(6), label="Sim. - T = 76 \$^\\circ\$C")

        elseif type_fuel_cell == "ZSW-GenStack_T_84" && aux_ok
            ax.plot(x, y; color=colors(7), label="Sim. - T = 84 \$^\\circ\$C - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_T_84"
            ax.plot(x, y; color=colors(7), label="Sim. - T = 84 \$^\\circ\$C")

        elseif type_fuel_cell == "EH-31_1.5" && aux_ok
            ax.plot(x, y; color=colors(0), label="Sim. - P = 1.5 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_1.5"
            ax.plot(x, y; color=colors(0), label="Sim. - P = 1.5 bar")

        elseif type_fuel_cell == "EH-31_2.0" && aux_ok
            ax.plot(x, y, "--"; color=colors(1), label="Sim. - P = 2.0 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_2.0"
            ax.plot(x, y; color=colors(1), label="Sim. - P = 2.0 bar")

        elseif type_fuel_cell == "EH-31_2.25" && aux_ok
            ax.plot(x, y, "--"; color=colors(2), label="Sim. - P = 2.25 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_2.25"
            ax.plot(x, y; color=colors(2), label="Sim. - P = 2.25 bar")

        elseif type_fuel_cell == "EH-31_2.5" && aux_ok
            ax.plot(x, y; color=colors(3), label="Sim - P = 2.5 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_2.5"
            ax.plot(x, y; color=colors(3), label="Sim - P = 2.5 bar")
        else
            ax.plot(x, y; color=colors(0), label="Simulation")
        end

    elseif type_current == "polarization_for_cali"
        if type_fuel_cell == "ZSW-GenStack" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(0), label="Sim. - nominal operating conditions - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(0), label="Sim. - nominal operating conditions")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(1), label="Sim. - P\$_a\$/P\$_c\$ = 1.61/1.41 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(1), label="Sim. - P\$_a\$/P\$_c\$ = 1.61/1.41 bar")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(2), label="Sim. - P\$_a\$/P\$_c\$ = 2.01/1.81 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(2), label="Sim. - P\$_a\$/P\$_c\$ = 2.01/1.81 bar")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(3), label="Sim. - P\$_a\$/P\$_c\$ = 2.4/2.2 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(3), label="Sim. - P\$_a\$/P\$_c\$ = 2.4/2.2 bar")

        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(4), label="Sim. - P\$_a\$/P\$_c\$ = 2.8/2.6 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(4), label="Sim. - P\$_a\$/P\$_c\$ = 2.8/2.6 bar")

        elseif type_fuel_cell == "ZSW-GenStack_T_62" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(5), label="Sim. - T = 62 \$^\\circ\$C - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_T_62"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(5), label="Sim. - T = 62 \$^\\circ\$C")

        elseif type_fuel_cell == "ZSW-GenStack_T_76" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(6), label="Sim. - T = 76 \$^\\circ\$C - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_T_76"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(6), label="Sim. - T = 76 \$^\\circ\$C")

        elseif type_fuel_cell == "ZSW-GenStack_T_84" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(7), label="Sim. - T = 84 \$^\\circ\$C - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "ZSW-GenStack_T_84"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(7), label="Sim. - T = 84 \$^\\circ\$C")

        elseif type_fuel_cell == "EH-31_1.5" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(0), label="Sim. - P = 1.5 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_1.5"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(0), label="Sim. - P = 1.5 bar")

        elseif type_fuel_cell == "EH-31_2.0" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(1), label="Sim. - P = 2.0 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_2.0"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(1), label="Sim. - P = 2.0 bar")

        elseif type_fuel_cell == "EH-31_2.25" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(2), label="Sim. - P = 2.25 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_2.25"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(2), label="Sim. - P = 2.25 bar")

        elseif type_fuel_cell == "EH-31_2.5" && aux_ok
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(3), label="Sim - P = 2.5 bar - \$ΔU_{RMSE}\$ = $(sim_error) %")
        elseif type_fuel_cell == "EH-31_2.5"
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(3), label="Sim - P = 2.5 bar")
        else
            ax.scatter(x, y; marker="o", linewidths=1.5, color=colors(0), label="Simulation")
        end
    else
        throw(ArgumentError("Only 'polarization' and 'polarization_for_cali' are considered here."))
    end

    return nothing
end






"""
    plot_experimental_polarisation_curve(type_fuel_cell, i_fc_t, U_exp_t, ax)

This function plots the experimental polarisation curve on the same graph as the model results.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell used in the model. This parameter includes the fuel cell used
  in the model and the corresponding operating conditions.
- `i_fc_t::Vector{<:Real}`: Current density values.
- `U_exp_t::Vector{<:Real}`: Experimental values of the voltage.
- `ax`: Axes-like object exposing `scatter` and `legend` methods.
"""
function plot_experimental_polarisation_curve(type_fuel_cell::String,
                                              i_fc_t::Vector{<:Real},
                                              U_exp_t::Vector{<:Real},
                                              ax)
    # ZSW-GenStack
    if type_fuel_cell == "ZSW-GenStack"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="s", color="black", label="Exp. - nominal")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="v", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 1.61/1.41 bar")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="^", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 2.01/1.81 bar")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="p", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 2.4/2.2 bar")
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="D", color="black", label="Exp. - P\$_a\$/P\$_c\$ = 2.8/2.6 bar")
    elseif type_fuel_cell == "ZSW-GenStack_T_62"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="P", color="black", label="Exp. - T = 62 \$^\\circ\$C")
    elseif type_fuel_cell == "ZSW-GenStack_T_76"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="X", color="black", label="Exp. - T = 76 \$^\\circ\$C")
    elseif type_fuel_cell == "ZSW-GenStack_T_84"
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="*", color="black", label="Exp. - T = 84 \$^\\circ\$C")

    # EH-31
    elseif type_fuel_cell == "EH-31_1.5"  # at 1.5 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="s", color="black", label="Exp. - P = 1.5 bar")
    elseif type_fuel_cell == "EH-31_2.0"  # at 2.0 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="v", color="black", label="Exp. - P = 2.0 bar")
    elseif type_fuel_cell == "EH-31_2.25"  # at 2.25 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="^", color="black", label="Exp. - P = 2.25 bar")
    elseif type_fuel_cell == "EH-31_2.5"  # at 2.5 bar
        ax.scatter(i_fc_t, U_exp_t; linewidths=1.5, marker="p", color="black", label="Exp. - P = 2.5 bar")
    else
        throw(ArgumentError("Unknown type_fuel_cell: $type_fuel_cell"))
    end

    ax.legend(loc="best", markerscale=0.5)
    return nothing
end



"""
    plot_EIS_Nyquist_instructions(type_fuel_cell, f_Fourier, x, y, ax)

This function adds the instructions for EIS plots according to the type_input to the ax object.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `f_Fourier::Number`: Frequency at which the EIS is simulated.
- `x::Number`: x-axis value for plotting the annotation.
- `y::Number`: y-axis value for plotting the annotation.
- `ax`: Axes on which the instructions will be added.
"""
function plot_EIS_Nyquist_instructions(type_fuel_cell::String,
                                       f_Fourier::Number,
                                       x::Number,
                                       y::Number,
                                       ax)

    # Common instructions
    ax.set_aspect("equal"; adjustable="box")  # Set orthonormal axis.
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis="both", which="major", size=10, width=1.5, direction="out")
    ax.tick_params(axis="both", which="minor", size=5, width=1.5, direction="out")
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
    # Adjust layout to prevent overlap between labels and the figure
    plt.show()  # Show the figure

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" || type_fuel_cell == "EH-31_2.0" ||
       type_fuel_cell == "EH-31_2.25" || type_fuel_cell == "EH-31_2.5"

        # Double charge transfer
        if f_Fourier >= 70 && f_Fourier <= 80
            freq_str = string(Int(f_Fourier)) * " Hz"  # Frequency annotation.
            ax.annotate(freq_str, (x, y); textcoords="offset points", xytext=(0, -40), ha="center", fontsize=14,
                        rotation=90, weight="bold")
        end
        # Auxiliary system
        if f_Fourier >= 0.14 && f_Fourier <= 0.16
            freq_str = string(round(f_Fourier, sigdigits=2)) * " Hz"  # Frequency annotation.
            ax.annotate(freq_str, (x, y); textcoords="offset points", xytext=(0, 7), ha="center", fontsize=14,
                        rotation=90, weight="bold")
        end
        if f_Fourier >= 1.2 && f_Fourier <= 1.4
            freq_str = string(round(f_Fourier, sigdigits=2)) * " Hz"  # Frequency annotation.
            ax.annotate(freq_str, (x, y); textcoords="offset points", xytext=(0, 10), ha="center", fontsize=14,
                        rotation=90, weight="bold")
        end
        # Diffusion
        if f_Fourier >= 0.015 && f_Fourier <= 0.020
            freq_str = string(round(f_Fourier, sigdigits=2)) * " Hz"  # Frequency annotation.
            ax.annotate(freq_str, (x, y); textcoords="offset points", xytext=(30, 0), ha="center", fontsize=14,
                        rotation=0, weight="bold")
        end
        if f_Fourier >= 0.9 && f_Fourier <= 1.1
            freq_str = string(round(f_Fourier, sigdigits=2)) * " Hz"  # Frequency annotation.
            ax.annotate(freq_str, (x, y); textcoords="offset points", xytext=(0, 10), ha="center", fontsize=14,
                        rotation=90, weight="bold")
        end
        if f_Fourier >= 70 && f_Fourier <= 90
            freq_str = string(Int(f_Fourier)) * " Hz"  # Frequency annotation.
            ax.annotate(freq_str, (x, y); textcoords="offset points", xytext=(0, -40), ha="center", fontsize=14,
                        rotation=90, weight="bold")
        end
        if f_Fourier >= 10000 && f_Fourier <= 12000
            freq_str = string(Int(f_Fourier)) * " Hz"  # Frequency annotation.
            ax.annotate(freq_str, (x, y); textcoords="offset points", xytext=(35, 0), ha="center", fontsize=14,
                        rotation=0, weight="bold")
        end
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(20))
        ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(20 / 5))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(10 / 5))
        ax.set_xlim(30, 200)
        ax.set_ylim(-25, 55)
    end

    return nothing
end


"""
    plot_Bode_amplitude_instructions(f_EIS, type_fuel_cell, ax)

This function adds the instructions for amplitude Bode plots according to the type_input to the ax object.

# Arguments
- `f_EIS`: Frequency parameters for EIS simulation.
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `ax`: Axes on which the instructions will be added.
"""
function plot_Bode_amplitude_instructions(f_EIS, type_fuel_cell::String, ax)

    # Common instructions
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS
    ax.set_xscale("log")  # Set logarithmic scale for the x-axis
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis="both", which="major", size=10, width=1.5, direction="out")
    ax.tick_params(axis="both", which="minor", size=5, width=1.5, direction="out")
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
    # Adjust layout to prevent overlap between labels and the figure
    plt.show()  # Show the figure

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" || type_fuel_cell == "EH-31_2.0" ||
       type_fuel_cell == "EH-31_2.25" || type_fuel_cell == "EH-31_2.5"
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=f_power_max_EIS - f_power_min_EIS + 1))
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) .* 0.1,
                                              numticks=(f_power_max_EIS - f_power_min_EIS + 1) * length(np.arange(2, 10))))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(30))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(30 / 5))
        ax.set_xlim([10^f_power_min_EIS, 10^f_power_max_EIS])
        # ax.set_ylim(0, 200)
    end

    return nothing
end


"""
    plot_Bode_phase_instructions(f_EIS, type_fuel_cell, ax)

This function adds the instructions for phase Bode plots according to the type_input to the ax object.

# Arguments
- `f_EIS`: Frequency parameters for EIS simulation.
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `ax`: Axes on which the instructions will be added.
"""
function plot_Bode_phase_instructions(f_EIS, type_fuel_cell::String, ax)

    # Common instructions
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS
    ax.set_xscale("log")  # Set logarithmic scale for the x-axis
    if !ax.yaxis_inverted()
        ax.invert_yaxis()  # Invert the y-axis
    end
    # Configure the appearance of major and minor ticks
    ax.tick_params(axis="both", which="major", size=10, width=1.5, direction="out")
    ax.tick_params(axis="both", which="minor", size=5, width=1.5, direction="out")
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
    # Adjust layout to prevent overlap between labels and the figure
    plt.show()  # Show the figure

    # For EH-31 fuel cell
    if type_fuel_cell == "EH-31_1.5" || type_fuel_cell == "EH-31_2.0" ||
       type_fuel_cell == "EH-31_2.25" || type_fuel_cell == "EH-31_2.5"
        ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=f_power_max_EIS - f_power_min_EIS + 1))
        ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) .* 0.1,
                                              numticks=(f_power_max_EIS - f_power_min_EIS + 1) * length(np.arange(2, 10))))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(5))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5 / 5))
        ax.set_xlim([10^f_power_min_EIS, 10^f_power_max_EIS])
        # ax.set_ylim(0, 360)
    end

    return nothing
end