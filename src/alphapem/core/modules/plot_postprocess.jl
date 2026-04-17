# -*- coding: utf-8 -*-

"""This module contains purely computational functions used for display purposes: Fourier
transformations, simulation error calculation, and axis-tick rounding helpers.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using FFTW
using Interpolations: linear_interpolation, Line
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

    isempty(t) && return FourierOutputs(ComplexF64[], ComplexF64[], Float64[], NaN, Float64[], NaN, 0)

    # Identify the active EIS segment for the current live run and keep only its measurement window.
    n_inf = searchsortedlast(cd.t_new_start, t[1])
    n_inf = clamp(n_inf, 1, length(cd.f))
    t_start = cd.t_new_start[n_inf]
    t_measure_start = t_start + cd.delta_t_break[n_inf]
    t_measure_end = t_measure_start + cd.delta_t_measurement[n_inf]

    mask_EIS = (t .>= t_measure_start) .& (t .<= t_measure_end)
    t_measured = t[mask_EIS]
    Ucell_measured = Ucell_t[mask_EIS]
    ifc_measured = ifc_t[mask_EIS]

    # Not enough raw samples yet to reconstruct one period at the target EIS resolution:
    # return NaNs so dynamic plotting can safely skip this update.
    min_raw_points = max(2, cd.nb_points)
    if length(t_measured) < min_raw_points
        return FourierOutputs(ComplexF64[], ComplexF64[], Float64[], NaN, Float64[], NaN, 0)
    end

    # FFT requires uniformly sampled data. Re-sample each segment using nb_points per period.
    dt = 1.0 / (cd.f[n_inf] * cd.nb_points)
    n_uniform = floor(Int, (t_measure_end - t_measure_start) / dt) + 1
    n_uniform = max(n_uniform, cd.nb_points + 1)
    t_uniform = collect(range(t_measure_start, stop=t_measure_end, length=n_uniform))

    itp_U = linear_interpolation(t_measured, Ucell_measured; extrapolation_bc=Line())
    itp_i = linear_interpolation(t_measured, ifc_measured; extrapolation_bc=Line())
    Ucell_EIS_measured = itp_U.(t_uniform)
    ifc_EIS_measured = itp_i.(t_uniform)

    # Determination of the Fourier transformation
    N             = length(Ucell_EIS_measured)              # Number of points used for the Fourier transformation
    Ucell_Fourier = fft(Ucell_EIS_measured)                 # Ucell Fourier transformation
    ifc_Fourier   = fft(ifc_EIS_measured)                   # ifc Fourier transformation
    A_period_t    = vcat([abs(Ucell_Fourier[1]) / N],       # Recovery of all amplitude values calculated by fft
                          abs.(Ucell_Fourier[2:N÷2]) .* 2 ./ N)

    # Ignore the DC component when searching the perturbation amplitude.
    if length(A_period_t) <= 1
        return FourierOutputs(ComplexF64.(Ucell_Fourier), ComplexF64.(ifc_Fourier), Float64.(A_period_t), NaN, Float64[], NaN, N)
    end

    idx_A_local = argmax(@view A_period_t[2:end])
    idx_A = idx_A_local + 1
    A = A_period_t[idx_A]

    dt_uniform = t_uniform[2] - t_uniform[1]
    freq_t = collect(0:(N ÷ 2 - 1)) ./ (N * dt_uniform)    # Frequencies in Hz associated with FFT bins
    f = freq_t[idx_A]

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


"""Return polarization points sampled at stabilization times (fixed mode).

The returned current density is in A.cm^-2 to match plotting conventions."""
function _polarization_points(outputs::SimulationOutputs,
                              cd::AbstractCurrent)
    t_hist = time_history(outputs)
    Ucell_t = derived_outputs(outputs).Ucell
    ifc_t = [current(cd, t) / 1e4 for t in t_hist]
    sample_indices = polarisation_sampling_indices(outputs, cd)
    return ifc_t[sample_indices], Ucell_t[sample_indices]
end

"""Return true when legacy RMSE comparison against experiments is enabled."""
function _pola_rmse_enabled(cfg::SimulationConfig)::Bool
    return cfg.type_fuel_cell != :manual_setup &&
           cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode, :no_auxiliary)
end

"""Compute polarization RMSE using legacy logic (interpolate model voltage on experimental currents)."""
function _polarization_rmse(ifc_discretized::AbstractVector{<:Real},
                            Ucell_discretized::AbstractVector{<:Real},
                            i_exp::AbstractVector{<:Real},
                            U_exp::AbstractVector{<:Real})
    itp = linear_interpolation(collect(ifc_discretized), collect(Ucell_discretized); extrapolation_bc=Line())
    Ucell_interpolated = itp.(collect(i_exp))
    return calculate_simulation_error(Ucell_interpolated, collect(U_exp))
end

"""Compute one EIS impedance point from Fourier outputs."""
function _eis_point(cd::AbstractCurrent,
                    Fourier_results::FourierOutputs)
    i_EIS = cd.i_EIS
    ratio_EIS = cd.ratio

    Z0 = Fourier_results.A / (ratio_EIS * (-i_EIS)) * 1e7

    theta_U_t = angle.(Fourier_results.Ucell_Fourier[1:Fourier_results.N÷2])
    theta_i_t = angle.(Fourier_results.ifc_Fourier[1:Fourier_results.N÷2])
    idx_A = findfirst(Fourier_results.A_period_t .== Fourier_results.A)
    idx_A === nothing && return (NaN, NaN, NaN, NaN, NaN)

    theta_U = theta_U_t[idx_A]
    theta_i = theta_i_t[idx_A]

    Z_real = Z0 * cos(theta_U - theta_i)
    Z_imag = Z0 * sin(theta_U - theta_i)

    # Same convention as the legacy implementation for Bode phase.
    phi_deg = mod((theta_U - (theta_i + π)) * 180 / π, 360)
    phi_deg > 180 && (phi_deg -= 360)

    return (Z_real, -Z_imag, Fourier_results.f, abs(Z0), phi_deg)
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
