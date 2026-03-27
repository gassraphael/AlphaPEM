# -*- coding: utf-8 -*-

"""This file contains the functions that generate the current densities for the simulation."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries

# Importing functions
if !isdefined(@__MODULE__, :pola_exp_values_calibration)
    include(joinpath(@__DIR__, "pola_exp_values.jl"))
end


# __________________________________________________Current densities___________________________________________________

"""
    step_current(t, parameters)

This function represents a step current density experiment. For the first delta_t_ini_step seconds, the current
 density is set to i_ini A.m-2 to allow the internal states of the fuel cell to stabilise. Then, the current density
 increases from i_ini to i_step A.m-2 in a step change over delta_t_load seconds. Finally, the current density
 remains at i_step A.m-2.
This is a C∞ function, which is advantageous for enhancing the overall stability of the results.

# Arguments
- `t: Time in seconds.
- `parameters::Dict`: A dictionary containing the parameters for the current density
  function.

# Returns
- The step current density at time t.
"""
function step_current(t, parameters::Dict)

    # Extraction of the parameters
    #   Initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_step = parameters["step_current_parameters"]["delta_t_ini_step"]  # (s).
    #   Initial current density used for the stabilisation of the internal states.
    i_ini = 1.0e4  # (A.m-2). This is the standard value for the initialisation.
    #   Loading time for the step current density function, from 0 to i_step.
    delta_t_load_step = parameters["step_current_parameters"]["delta_t_load_step"]  # (s).
    #   Current density for the step current density function.
    i_step = parameters["step_current_parameters"]["i_step"]  # (A.m-2).

    # Step current density
    return i_ini * (1.0 + tanh(4 * (t - (delta_t_load_step / 2)) / (delta_t_load_step / 2))) / 2 +
           (i_step - i_ini) * (1.0 + tanh(4 * (t - delta_t_ini_step - (delta_t_load_step / 2)) /
                                          (delta_t_load_step / 2))) / 2
end


"""
    polarization_current(t, parameters)

This function represents a current density used for creating a polarization curve. For the first
 delta_t_ini_step seconds, the current density is set to i_ini A.m-2 to allow the internal states of the fuel cell to
 stabilise. Then, the current density increases by the value of delta_i_pola every delta_t, following C∞ step current
 increments, until it reaches i_max_pola. Each increment lasts for delta_t_load_pola seconds. After each increment,
 there is a pause of delta_t_break_pola seconds to allow the stack to reach equilibrium.

# Arguments
- `t: Time in seconds.
- `parameters::Dict`: A dictionary containing the parameters for the current density
  function.

# Returns
- The polarization current density at time t.
"""
function polarization_current(t, parameters::Dict)

    # Extraction of the parameters
    #   Initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola = parameters["pola_current_parameters"]["delta_t_ini_pola"]  # (s).
    #   Loading time for one step current of the polarisation current density function.
    delta_t_load_pola = parameters["pola_current_parameters"]["delta_t_load_pola"]  # (s).
    #   Breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola = parameters["pola_current_parameters"]["delta_t_break_pola"]  # (s).
    #   Current density step for the polarisation current density function.
    delta_i_pola = parameters["pola_current_parameters"]["delta_i_pola"]  # (A.m-2).
    #   Maximum current density for the polarization curve.
    i_max_pola = parameters["pola_current_parameters"]["i_max_pola"]  # (A.m-2).

    # Calculation of the time parameters
    #   Time of one step.
    delta_t = delta_t_load_pola + delta_t_break_pola  # (s).
    #   Duration of this polarization curve.
    tf = delta_t_ini_pola + Int(floor(i_max_pola / delta_i_pola + 1)) * delta_t  # (s).
    #   Number of loads made for this polarization curve.
    n = Int(floor(tf / delta_t))

    # Current density for the polarization curve
    i_fc = 0.0  # A.m-2. Initialisation of the current density.
    for i in 0:(n - 1)
        i_fc += delta_i_pola * (1.0 + tanh(4 * (t - delta_t_ini_pola - i * delta_t - (delta_t_load_pola / 2)) /
                                            (delta_t_load_pola / 2))) / 2
    end
    return i_fc
end


"""
    polarization_current_for_calibration(t, parameters)

This function represents a current density used for creating a polarisation curve dedicated to the calibration of
 a specific fuel cell. The principle is similar to the polarization_current function, but it uses experimental values
 for the current density load.

# Arguments
- `t`: Time in seconds.
- `parameters::Dict`: A dictionary containing the parameters for the current density
  function.

# Returns
- The polarisation current density at time t.
"""
function polarization_current_for_calibration(t, parameters::Dict)

    # Extraction of the parameters
    #   Initial time at zero current density for the stabilisation of the internal states.
    delta_t_ini_pola_cali = parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]  # (s).
    #   Loading time for one step current of the polarisation current density function.
    delta_t_load_pola_cali = parameters["pola_current_for_cali_parameters"]["delta_t_load_pola_cali"]  # (s).
    #   Breaking time for one step current, for the stabilisation of the internal states.
    delta_t_break_pola_cali = parameters["pola_current_for_cali_parameters"]["delta_t_break_pola_cali"]  # (s).
    type_fuel_cell = parameters["type_fuel_cell"]
    voltage_zone = parameters["voltage_zone"]  # The fuel cell for which the calibration is performed.
    i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(type_fuel_cell, voltage_zone)  # (A.m-2, V). It is the experimental
    #                                                           current density and voltage values for the calibration.

    # Calculation of the time parameters
    #   Time of one step.
    delta_t = delta_t_load_pola_cali + delta_t_break_pola_cali  # (s).

    # Current density for the polarization curve used for calibration
    i_fc = 0.0  # A.m-2. Initialisation of the current density.
    for e in eachindex(i_exp_cali_t)
        if e == 1
            i_fc += i_exp_cali_t[1] * (1.0 + tanh(4 * (t - delta_t_ini_pola_cali - (delta_t_load_pola_cali / 2)) /
                                                   (delta_t_load_pola_cali / 2))) / 2
        else
            delta_i_exp_cali = (i_exp_cali_t[e] - i_exp_cali_t[e - 1])  # (A.m-2). It is the difference between the
            #                                                  current density at the current step and the previous one.
            i_fc += delta_i_exp_cali * (1.0 + tanh(4 * (t - delta_t_ini_pola_cali - (e - 1) * delta_t - (delta_t_load_pola_cali / 2)) /
                                                    (delta_t_load_pola_cali / 2))) / 2
        end
    end
    return i_fc
end


"""
    EIS_current(t, parameters)

Represents a current density used for creating an EIS curve and Bode diagrams.
The current density is first equilibrated at i_EIS A.m-2 from 0 to t0_EIS seconds using a step increase.
Then, a sinusoidal perturbation is added to the current density. This perturbation has an amplitude of
(ratio_EIS * i_EIS) A.m-2 and a frequency of f[n_inf] Hz.

# Arguments
- `t`: Time in seconds.
- `parameters::Dict`: A dictionary containing the parameters for the current density
  function.

# Returns
- The polarization current density at time t.
"""
function EIS_current(t, parameters::Dict)

    # Initialisation
    i_EIS = parameters["i_EIS"]
    ratio_EIS = parameters["ratio_EIS"]  # (A/m², ). i_EIS is the current for which a
    #                                      ratio_EIS perturbation is added.
    t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = parameters["t_EIS"]  # It is the initial
    #         EIS time after stack equilibrium, a list of time parameters which gives the beginning of each frequency
    #         change, the final time, a list of time parameters which gives the estimated time for reaching equilibrium
    #         at each frequency, and a list of time parameters which gives the estimated time for measuring the voltage
    #         response at each frequency.
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = parameters["f_EIS"]  # It is the power of the initial
    #         frequency: f_min_EIS = 10^f_power_min_EIS, the power of the final frequency, the number of frequencies
    #         tested and the number of points calculated per specific period.
    f = 10.0 .^ range(f_power_min_EIS, f_power_max_EIS; length=Int(nb_f_EIS))  # It is a list of all the frequency tested,
    #                                                                             ranged logarithmically.

    # Current density for the EIS curve
    if t < t0_EIS
        delta_t_ini = 3 * 60  # s. It is the required time for elevating i_fc from 0 to i_EIS without starving the cell.
        i_fc = i_EIS * (1.0 + tanh(4 * (t - 2 * (delta_t_ini / 2)) / delta_t_ini)) / 2
    else
        n_inf = searchsortedlast(t_new_start_EIS, t)  # It is the number of frequency changes which has been made so far.
        n_inf = max(1, n_inf)
        i_disruption = (ratio_EIS * i_EIS) * cos(2 * pi * f[n_inf] * t)
        i_fc = i_EIS + i_disruption
    end

    return i_fc
end


"""
    EIS_parameters(f_EIS)

This function gives the time parameters for the EIS_current density function.

# Arguments
- f_power_min_EIS::Float64
        Power of the initial frequency for the EIS curve, such that f_min_EIS = 10^f_power_min_EIS.
- f_power_max_EIS::Float64
        Power of the final frequency for the EIS curve, such that f_max_EIS = 10^f_power_max_EIS.
- nb_f_EIS::Float64
        Number of frequencies tested for the EIS curve.
- nb_points_EIS::Float64
        Number of points calculated per specific period for the EIS curve.

# Returns
- `Tuple{Float64, Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}}`: EIS parameters. It is a tuple
  containing the initial EIS time after stack equilibrium `t0_EIS`, a list of time parameters which gives the
  beginning of each frequency change `t_new_start_EIS`, the final time `tf_EIS`, a list of time parameters which gives
  the estimated time for reaching equilibrium at each frequency `delta_t_break_EIS`, and a list of time parameters
  which gives the estimated time for measuring the voltage response at each frequency `delta_t_measurement_EIS`.
"""
function EIS_parameters(f_power_min_EIS::Float64, f_power_max_EIS::Float64, nb_f_EIS::Float64,
                        nb_points_EIS::Float64)::Tuple{Float64, Vector{Float64}, Float64, Vector{Float64}, Vector{Float64}}

    # Initialisation
    #       Frequencies
    f = 10.0 .^ range(f_power_min_EIS, f_power_max_EIS; length=nb_f_EIS)  # It is the tested frequencies
    nb_period_break_EIS, nb_period_measurement_EIS = 50, 50  # They are the nu    f_EIS = (-3.0, 5.0, 90.0, 50.0)   # → NTuple{4, Float64} ✓    f_EIS = (-3.0, 5.0, 90.0, 50.0)   # → NTuple{4, Float64} ✓    f_EIS = (-3.0, 5.0, 90.0, 50.0)   # → NTuple{4, Float64} ✓mber of temporal periods which are used
    #                                                          for break and for measurement. It is more accurate to use
    #                                                          periods than time as the frequency range is big.
    #       Time parameters
    delta_t_break_EIS = Float64[]  # It is the estimated time for reaching equilibrium at each frequency.
    delta_t_measurement_EIS = Float64[]  # It is the estimated time for measuring the voltage response.

    # Time parameters calculation
    t0_EIS = 120 * 60.0  # s. It is the simulation starting time. [0, t0_EIS] is used to let the stack equilibrate to i_EIS.
    t_new_start_EIS = [t0_EIS]  # It is a list of time parameters which gives the beginning of each frequency
    #                            change.
    tf_EIS = t0_EIS  # s. Default value; will be updated in the loop below.
    for i in 1:nb_f_EIS  # The goal is to measure nb_f_EIS periods of the signal in order to have precise enough values.
        T_i = 1 / f[i]  # s. It is the period of the signal.
        push!(delta_t_break_EIS, nb_period_break_EIS * T_i)
        push!(delta_t_measurement_EIS, nb_period_measurement_EIS * T_i)
        if i < nb_f_EIS
            next_start_EIS = t_new_start_EIS[i] + delta_t_break_EIS[i] + delta_t_measurement_EIS[i]
            push!(t_new_start_EIS, next_start_EIS)
        else
            tf_EIS = t_new_start_EIS[end] + delta_t_break_EIS[end] + delta_t_measurement_EIS[end]  # s. It is the
            #                                                                                         simulation ending time
        end
    end

    t_EIS = t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS
    return t_EIS
end


# _________________________________________________Test of the program__________________________________________________


# NOTE:
# The plotting test block from Python (__main__) is intentionally left out here to keep this
# file free of plotting-package dependencies at runtime. The numerical functions above are fully translated.

# Translated test block (kept commented on purpose)
# if abspath(PROGRAM_FILE) == @__FILE__
#     using Plots
#     colors = palette(:tab10)
#     p = plot(layout=(1, 4), size=(2400, 600))
#
#     # Tests for step_current curve:
#     #   Step current parameters
#     delta_t_ini_step = 30 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
#     delta_t_load_step = 30  # (s). Loading time for the step current density function, from 0 to i_step.
#     delta_t_break_step = 15 * 60  # (s). Time at i_step current density for the stabilization of the internal states.
#     i_step = 1.5e4  # (A.m-2). Current density for the step current density function.
#     step_current_parameters = Dict("delta_t_ini_step" => delta_t_ini_step, "delta_t_load_step" => delta_t_load_step,
#                                    "delta_t_break_step" => delta_t_break_step, "i_step" => i_step)
#     parameters = Dict("step_current_parameters" => step_current_parameters)
#     #   Display
#     n = 10000
#     t = range(0, delta_t_ini_step + delta_t_load_step + delta_t_break_step; length=n)
#     i_fc_step = [step_current(ti, parameters) / 10000 for ti in t]  # Conversion in A/cm²
#     plot!(p[1], t, i_fc_step, color=colors[1], label="i_fc (step)")
#
#     # Tests for polarization curves:
#     #   Polarization current parameters
#     delta_t_ini_pola = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
#     delta_t_load_pola = 30  # (s). Loading time for one step current of the polarisation current density function.
#     delta_t_break_pola = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
#     delta_i_pola = 0.1e4  # (A.m-2). Current density step for the polarisation current density function.
#     i_max_pola = 2.5e4  # (A.m-2). Maximum current density for the polarization curve.
#     pola_current_parameters = Dict("delta_t_ini_pola" => delta_t_ini_pola, "delta_t_load_pola" => delta_t_load_pola,
#                                    "delta_t_break_pola" => delta_t_break_pola, "delta_i_pola" => delta_i_pola,
#                                    "i_max_pola" => i_max_pola)
#     parameters = Dict("pola_current_parameters" => pola_current_parameters)
#     #   Display
#     n = 10000
#     tf = delta_t_ini_pola + Int(floor(i_max_pola / delta_i_pola)) * (delta_t_load_pola + delta_t_break_pola)
#     t = range(0, tf; length=n)
#     i_fc_pola = [polarization_current(ti, parameters) / 1e4 for ti in t]  # Conversion in A/cm²
#     plot!(p[2], t, i_fc_pola, color=colors[2], label="i_fc (pola)")
#
#     # Tests for calibration curves:
#     #   Polarization current for calibration parameters
#     delta_t_ini_pola_cali = 120 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
#     delta_t_load_pola_cali = 30  # (s). Loading time for one step current of the polarisation current density function.
#     delta_t_break_pola_cali = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
#     pola_current_for_cali_parameters = Dict("delta_t_ini_pola_cali" => delta_t_ini_pola_cali,
#                                             "delta_t_load_pola_cali" => delta_t_load_pola_cali,
#                                             "delta_t_break_pola_cali" => delta_t_break_pola_cali)
#     parameters = Dict("pola_current_for_cali_parameters" => pola_current_for_cali_parameters,
#                       "type_fuel_cell" => "EH-31_2.0", "voltage_zone" => "full")
#     i_exp_cali_t, U_exp_cali_t = pola_exp_values_calibration(parameters["type_fuel_cell"], parameters["voltage_zone"])
#     #   Display
#     n = 10000
#     tf = delta_t_ini_pola_cali + length(i_exp_cali_t) * (delta_t_load_pola_cali + delta_t_break_pola_cali)
#     t = range(0, tf; length=n)
#     i_fc_cali = [polarization_current_for_calibration(ti, parameters) / 1e4 for ti in t]  # Conversion in A/cm²
#     plot!(p[3], t, i_fc_cali, color=colors[3], label="i_fc (pola)")
#
#     # Tests for EIS curve:
#     #   EIS_current parameters
#     i_EIS, ratio_EIS = 1.0e4, 5 / 100  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
#     f_EIS = (-3, 5, 90, 50)  # Frequency parameters for the EIS_current density function.
#     t_EIS = EIS_parameters(f_EIS)  # It is the EIS parameters.
#     parameters = Dict("i_EIS" => i_EIS, "ratio_EIS" => ratio_EIS, "f_EIS" => f_EIS, "t_EIS" => t_EIS)
#
#     f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = f_EIS
#     t0_EIS, t_new_start_EIS, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = t_EIS
#     f = 10.0 .^ range(f_power_min_EIS, f_power_max_EIS; length=nb_f_EIS)  # It is a list of all the frequency tested.
#
#     # Step current for reaching i_EIS
#     n = 1000
#     t = range(0, t0_EIS; length=n)
#     i_fc_EIS = [EIS_current(ti, parameters) / 10000 for ti in t]  # Conversion in A/cm²
#     plot!(p[4], t, i_fc_EIS, color=colors[4], label="")
#
#     # EIS current density
#     for i in eachindex(t_new_start_EIS)  # The EIS curve is displayed for each frequency change to reduce points.
#         t0_EIS_temp = t_new_start_EIS[i]
#         tf_EIS_temp = t_new_start_EIS[i] + delta_t_break_EIS[i] + delta_t_measurement_EIS[i]
#         n_inf = findlast(x -> x <= t0_EIS_temp, t_new_start_EIS)  # It is the number of frequency changes so far.
#         n = Int(floor(f[n_inf] * (tf_EIS_temp - t0_EIS_temp) * nb_points_EIS))  # It is the number of points to calculate.
#         t = range(t0_EIS_temp, tf_EIS_temp; length=max(n, 2))  # It is the time interval for this frequency.
#         i_fc_EIS = [EIS_current(tj, parameters) / 10000 for tj in t]  # Conversion in A/cm²
#         plot!(p[4], t, i_fc_EIS, color=colors[4], label="")
#     end
#
#     plot!(p[4], [], [], label="i_fc (EIS)")
#     display(p)
# end

