# -*- coding: utf-8 -*-

"""This module contains purely computational functions used for display purposes: Fourier
transformations, simulation error calculation, and axis-tick rounding helpers.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using FFTW
using Statistics


# ___________________________________________Computational display functions____________________________________________

"""
    make_Fourier_transformation(outputs, cd, cfg)

This function calculates the Fourier transformation of both cell voltage and current density. It will be used to
display the Nyquist and Bode diagrams.
To generate it at each frequency change, the cell voltage and the current density are recorded. The time for which
these points are captured is determined using the following approach: at the beginning of each frequency change, a
delta_t_break_EIS time is observed to ensure the dynamic stability of the stack's variables. Subsequently, a
delta_t_measurement_EIS time is needed to record the cell voltage and the current density.

# Arguments
- `outputs::SimulationOutputs`: Typed simulation outputs used by post-processing.
- `cd::AbstractCurrent`: Current profile used by the simulation.
- `cfg::SimulationConfig`: Simulation configuration.

# Returns
- `FourierOutputs`: Structured Fourier post-processing outputs.
"""
function make_Fourier_transformation(outputs::SimulationOutputs,
                                     cd::AbstractCurrent,
                                     cfg::SimulationConfig)::FourierOutputs

    # Extraction of the variables
    t = time_history(outputs)
    Ucell_t = derived_outputs(outputs).Ucell
    # EIS timing is only available for EIS current profiles.
    cfg.type_current isa EISParams ||
        throw(ArgumentError("make_Fourier_transformation requires type_current isa EISParams."))

    # Creation of the current density vector at the same time points as the cell voltage.
    ifc_t = current(cd, t)

    # Identify the areas where Ucell and ifc can be measured for the EIS: after equilibrium and at each frequency change
    t_new_start_EIS = cd.t_new_start
    delta_t_break_EIS = cd.delta_t_break
    delta_t_measurement_EIS = cd.delta_t_measurement
    n_inf = findlast(t_new_start_EIS .<= t[1])  # Number of frequency changes already applied.
    mask_EIS = (t .> (t[1] + delta_t_break_EIS[n_inf])) .& (t .< (t[1] + delta_t_break_EIS[n_inf] + delta_t_measurement_EIS[n_inf]))
    Ucell_EIS_measured = Ucell_t[mask_EIS]
    ifc_EIS_measured   = ifc_t[mask_EIS]

    # Determination of the Fourier transformation
    N             = length(Ucell_EIS_measured)              # Number of points used for the Fourier transformation
    Ucell_Fourier = fft(Ucell_EIS_measured)                  # Ucell Fourier transformation
    ifc_Fourier   = fft(ifc_EIS_measured)                    # ifc Fourier transformation
    A_period_t    = vcat([abs(Ucell_Fourier[1]) / N],        # Recovery of all amplitude values calculated by fft
                          abs.(Ucell_Fourier[2:N÷2]) .* 2 ./ N)
    A      = maximum(A_period_t[2:end])                      # Amplitude at the frequency of the perturbation
    freq_t = fftfreq(N)[1:N÷2]                              # Recovery of all frequency values used by fft
    f      = freq_t[findfirst(A_period_t .== A)]             # Recovery of the studied frequency

    return FourierOutputs(
        ComplexF64.(Ucell_Fourier),
        ComplexF64.(ifc_Fourier),
        Float64.(A_period_t),
        Float64(A),
        Float64.(freq_t),
        Float64(f),
        N,
    )
end


"""
    calculate_simulation_error(Ucell, U_exp_t)

This function calculates the simulation error between the simulated cell voltage and the experimental cell
voltage. It is calculated as the RMSE (root-mean-square error) of the relative differences (in %).

# Arguments
- `Ucell::Vector`: Simulated cell voltage, interpolated at the experimental measurement points.
- `U_exp_t::Vector`: Experimental cell voltage.

# Returns
- RMSE between the simulated cell voltage and the experimental cell voltage (in %).
"""
function calculate_simulation_error(Ucell::Vector,
                                    U_exp_t::Vector)

    # Distance between the simulated and the experimental polarization curves (RMSE: root-mean-square error).
    res1 = (Ucell .- U_exp_t) ./ U_exp_t .* 100  # in %.
    return round(sqrt(mean(res1 .^ 2)), digits=2)
end


"""
    round_nice(x)

Round the main step to a "nice" number.

# Arguments
- `x`: The value to be rounded.

# Returns
- The value rounded to a "nice" number.
"""
function round_nice(x)
    exp  = floor(log10(x))
    f    = x / 10^exp
    nice = if f < 1.5
        1
    elseif f < 3
        2
    elseif f < 7
        5
    else
        10
    end
    return nice * 10^exp
end



