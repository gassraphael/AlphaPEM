# -*- coding: utf-8 -*-

"""This module contains purely computational functions used for display purposes: Fourier
transformations, simulation error calculation, and axis-tick rounding helpers.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using FFTW
using Interpolations: linear_interpolation, Line, deduplicate_knots!
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
- `Vector{FourierOutputs}`: List of structured Fourier post-processing outputs.
"""
function make_Fourier_transformation(outputs::SimulationOutputs,
                                     cd::AbstractCurrent,
                                     cfg::SimulationConfig)::Vector{FourierOutputs}

    # Extraction of the variables
    t = time_history(outputs)
    Ucell_t = derived_outputs(outputs).Ucell

    # EIS timing is only available for EIS current profiles.
    cfg.type_current isa EISParams ||
        throw(ArgumentError("make_Fourier_transformation requires type_current isa EISParams."))

    if isempty(t)
        return [FourierOutputs(ComplexF64[], ComplexF64[], Float64[], NaN, Float64[], NaN, 0)]
    end

    if cfg.display_timing == :live
        # Identify the active EIS segment for the current live run (latest point in history).
        n_inf = searchsortedlast(cd.t_new_start, t[end])
        n_inf = clamp(n_inf, 1, length(cd.f))
        return [_compute_fourier_for_segment(t, Ucell_t, cd, n_inf)]
    else
        # :postrun mode: Iterate through all frequency segments and extract all points.
        results = FourierOutputs[]
        for i in 1:length(cd.f)
            res = _compute_fourier_for_segment(t, Ucell_t, cd, i)
            # Only include segments that have been fully simulated and processed.
            if !isnan(res.f)
                push!(results, res)
            end
        end
        return results
    end
end


"""
    _compute_fourier_for_segment(t, Ucell_t, cd, n_inf)

Helper function to calculate the Fourier transform for a specific EIS frequency segment.
"""
function _compute_fourier_for_segment(t::AbstractVector{Float64},
                                      Ucell_t::AbstractVector{Float64},
                                      cd::AbstractCurrent,
                                      n_inf::Int)::FourierOutputs

    # Identify measurement window for this frequency segment.
    t_measure_start = cd.t_new_start[n_inf] + cd.delta_t_break[n_inf]
    t_measure_end = t_measure_start + cd.delta_t_measurement[n_inf]

    mask_EIS = (t .>= t_measure_start) .& (t .<= t_measure_end)
    t_measured = t[mask_EIS]
    Ucell_measured = Ucell_t[mask_EIS]

    # Check if we have enough points for reconstruction.
    min_raw_points = max(2, cd.nb_points)
    if length(t_measured) < min_raw_points
        return FourierOutputs(ComplexF64[], ComplexF64[], Float64[], NaN, Float64[], NaN, 0)
    end

    # Current density computed only on the measurement window.
    ifc_measured = current(cd, t_measured)

    # FFT requires uniformly sampled data. Re-sample each segment using nb_points per period.
    dt = 1.0 / (cd.f[n_inf] * cd.nb_points)
    n_uniform = floor(Int, (t_measure_end - t_measure_start) / dt) + 1
    n_uniform = max(n_uniform, cd.nb_points + 1)
    t_uniform = collect(range(t_measure_start, stop=t_measure_end, length=n_uniform))

    # Interpolation to uniform grid.
    t_meas = copy(t_measured)
    deduplicate_knots!(t_meas)
    itp_U = linear_interpolation(t_meas, Ucell_measured; extrapolation_bc=Line())
    itp_i = linear_interpolation(t_meas, ifc_measured; extrapolation_bc=Line())
    Ucell_EIS_measured = itp_U.(t_uniform)
    ifc_EIS_measured = itp_i.(t_uniform)

    # Determination of the Fourier transformation
    N             = length(Ucell_EIS_measured)              # Number of points used for the Fourier transformation
    Ucell_Fourier = fft(Ucell_EIS_measured)                 # Ucell Fourier transformation
    ifc_Fourier   = fft(ifc_EIS_measured)                   # ifc Fourier transformation
    A_period_t    = vcat([abs(Ucell_Fourier[1]) / N],       # Recovery of all amplitude values calculated by fft
                          abs.(Ucell_Fourier[2:N÷2]) .* 2 ./ N)

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


"""Return polarization points sampled at stabilization times (fixed mode).

The returned current density is in A.cm^-2 to match plotting conventions.
If `average` is true, multiple points at the same current are averaged.
"""
function _polarization_points(outputs::SimulationOutputs,
                              cd::AbstractCurrent;
                              average::Bool=true)
    # Extract time history and cell voltage
    t_hist = time_history(outputs)
    Ucell_t = derived_outputs(outputs).Ucell

    # Calculate current density and convert to A/cm² for display
    ifc_t = [current(cd, t) / 1e4 for t in t_hist]

    # Identify indices corresponding to stabilized points (end of step)
    sample_indices = polarisation_sampling_indices(outputs, cd)
    
    ifc_samples = ifc_t[sample_indices]
    Ucell_samples = Ucell_t[sample_indices]

    # If requested, group and average voltages for each current level (e.g., multiple cycles)
    if average
        seen_i = Dict{Float64, Vector{Float64}}()
        for (i, u) in zip(ifc_samples, Ucell_samples)
            # Round to group points despite potential numerical micro-variations
            i_key = round(i, digits=3)
            if haskey(seen_i, i_key)
                push!(seen_i[i_key], u)
            else
                seen_i[i_key] = [u]
            end
        end
        
        # Sort by current density to ensure a monotonic curve for display
        sorted_i = sort(collect(keys(seen_i)))
        return sorted_i, [mean(seen_i[i]) for i in sorted_i]
    end

    return ifc_samples, Ucell_samples
end

"""Return true when legacy RMSE comparison against experiments is enabled."""
function _pola_rmse_enabled(cfg::SimulationConfig)::Bool
    return cfg.type_fuel_cell != :manual_setup &&
           cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode, :no_auxiliary)
end

"""Compute polarization RMSE with interpolation (interpolate model voltage on experimental currents)."""
function _interpolated_rmse(ifc_discretized::AbstractVector{<:Real},
                            Ucell_discretized::AbstractVector{<:Real},
                            i_exp::AbstractVector{<:Real},
                            U_exp::AbstractVector{<:Real})
    ifc_vec = collect(ifc_discretized)
    ucell_vec = collect(Ucell_discretized)
    deduplicate_knots!(ifc_vec)
    itp = linear_interpolation(ifc_vec, ucell_vec; extrapolation_bc=Line())
    Ucell_interpolated = itp.(collect(i_exp))
    return _calculate_rmse(Ucell_interpolated, collect(U_exp))
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


"""
    calculate_reynolds_numbers(outputs, fc)

Compute the Reynolds number at each gas-channel node for both anode (AGC) and cathode (CGC).

The hydraulic diameter is computed as Dh = 2HW/(H+W) for a rectangular cross-section.
Gas-mixture density is derived from concentrations; dynamic viscosity uses the
mixture rule from `mu_mixture_gases`.

# Returns
- `Re_a`: Vector of vectors (per node) of anode Reynolds numbers over time.
- `Re_c`: Vector of vectors (per node) of cathode Reynolds numbers over time.
"""
function calculate_reynolds_numbers(outputs::SimulationOutputs,
                                    fc::AbstractFuelCell)
    nb_gc = gas_channel_count(outputs)
    pp = fc.physical_parameters
    t = time_history(outputs)
    n_t = length(t)

    # Hydraulic diameters of the rectangular gas channels.
    Dh_a = 2 * pp.Hagc * pp.Wagc / (pp.Hagc + pp.Wagc)
    Dh_c = 2 * pp.Hcgc * pp.Wcgc / (pp.Hcgc + pp.Wcgc)

    Re_a = [Vector{Float64}(undef, n_t) for _ in 1:nb_gc]
    Re_c = [Vector{Float64}(undef, n_t) for _ in 1:nb_gc]

    for i in 1:nb_gc
        # Gas velocities.
        v_a_i = extract_derived_gc_series(outputs, i, x -> x.v_a)
        v_c_i = extract_derived_gc_series(outputs, i, x -> x.v_c)
        # Anode GC state.
        C_v_agc  = extract_mea_series(outputs, i, mea -> mea.agc.C_v)
        C_H2_agc = extract_mea_series(outputs, i, mea -> mea.agc.C_H2)
        T_agc    = extract_mea_series(outputs, i, mea -> mea.agc.T)
        # Cathode GC state.
        C_v_cgc  = extract_mea_series(outputs, i, mea -> mea.cgc.C_v)
        C_O2_cgc = extract_mea_series(outputs, i, mea -> mea.cgc.C_O2)
        C_N2_cgc = extract_mea_series(outputs, i, mea -> mea.cgc.C_N2)
        T_cgc    = extract_mea_series(outputs, i, mea -> mea.cgc.T)

        for j in 1:n_t
            # Anode: H₂O + H₂ mixture.
            x_H2O_a = C_v_agc[j] + C_H2_agc[j] > 0 ?
                      C_v_agc[j] / (C_v_agc[j] + C_H2_agc[j]) : 0.0
            rho_a   = C_v_agc[j] * M_H2O + C_H2_agc[j] * M_H2
            mu_a    = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_a, 1 - x_H2O_a], T_agc[j])
            Re_a[i][j] = mu_a > 0 ? rho_a * abs(v_a_i[j]) * Dh_a / mu_a : 0.0

            # Cathode: H₂O + O₂ + N₂ mixture.
            C_dry_c  = max(C_O2_cgc[j] + C_N2_cgc[j], eps(Float64))
            y_O2_c   = C_O2_cgc[j] / C_dry_c
            x_H2O_c  = C_v_cgc[j] + C_O2_cgc[j] + C_N2_cgc[j] > 0 ?
                       C_v_cgc[j] / (C_v_cgc[j] + C_O2_cgc[j] + C_N2_cgc[j]) : 0.0
            x_O2_c   = y_O2_c * (1 - x_H2O_c)
            x_N2_c   = (1 - y_O2_c) * (1 - x_H2O_c)
            rho_c    = C_v_cgc[j] * M_H2O + C_O2_cgc[j] * M_O2 + C_N2_cgc[j] * M_N2
            mu_c     = mu_mixture_gases(["H2O_v", "O2", "N2"], [x_H2O_c, x_O2_c, x_N2_c], T_cgc[j])
            Re_c[i][j] = mu_c > 0 ? rho_c * abs(v_c_i[j]) * Dh_c / mu_c : 0.0
        end
    end
    return Re_a, Re_c
end
