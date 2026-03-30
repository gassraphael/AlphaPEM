# -*- coding: utf-8 -*-

"""This file is used to construct the EIS current density."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing functions
include(joinpath(@__DIR__, "abstract.jl"))


# __________________________________________________EIS current density_________________________________________________

"""
    EISCurrent

Current density model used for Electrochemical Impedance Spectroscopy (EIS).

After an initial stabilization phase at constant current, a sinusoidal
perturbation is applied with varying frequencies.

# Fields
- `i_EIS::Float64`: Base current density (A/m²)
- `ratio::Float64`: Relative amplitude of the perturbation
- `f::Vector{Float64}`: Frequencies (Hz)
- `t0::Float64`: Initial stabilization time (s)
- `t_new_start::Vector{Float64}`: Start times for each frequency segment (s)
- `tf::Float64`: Final simulation time (s)
- `delta_t_break::Vector{Float64}`: Stabilization durations (s)
- `delta_t_measurement::Vector{Float64}`: Measurement durations (s)
"""
struct EISCurrent <: AbstractCurrent
    i_EIS::Float64
    ratio::Float64
    f::Vector{Float64}
    t0::Float64
    t_new_start::Vector{Float64}
    tf::Float64
    delta_t_break::Vector{Float64}
    delta_t_measurement::Vector{Float64}
end

# --- Constructor with precomputation --------------------------------------

# This constructor computes all the time intervals and frequency steps for the EIS protocol.
# It generates logarithmically spaced frequencies, and for each frequency, it computes the break and measurement durations.
# The start time for each frequency segment is also precomputed for efficient lookup during simulation.
function EISCurrent(
    i_EIS::Real,
    ratio::Real,
    f_power_min::Real,
    f_power_max::Real,
    nb_f::Integer,
    nb_points::Integer
)
    # Check input validity
    nb_f > 0 || throw(ArgumentError("nb_f must be > 0"))
    nb_points > 0 || throw(ArgumentError("nb_points must be > 0"))
    ratio ≥ 0 || throw(ArgumentError("ratio must be ≥ 0"))

    # Generate logarithmically spaced frequencies
    f = 10.0 .^ range(f_power_min, f_power_max; length=nb_f)

    nb_period_break = 50  # Number of periods for stabilization at each frequency
    nb_period_measurement = 50  # Number of periods for measurement at each frequency

    delta_t_break = Float64[]  # List of stabilization durations for each frequency
    delta_t_measurement = Float64[]  # List of measurement durations for each frequency

    t0 = 120.0 * 60.0  # Initial stabilization time (s)
    t_new_start = [t0]  # Start times for each frequency segment
    tf = t0  # Final simulation time (will be updated in the loop below)

    # For each frequency, compute the durations and update the start times
    for i in 1:nb_f
        T = 1 / f[i]  # Period for current frequency

        push!(delta_t_break, nb_period_break * T)
        push!(delta_t_measurement, nb_period_measurement * T)

        if i < nb_f
            # Compute the start time for the next frequency segment
            next_start = t_new_start[i] + delta_t_break[i] + delta_t_measurement[i]
            push!(t_new_start, next_start)
        else
            # For the last frequency, set the final simulation time
            tf = t_new_start[end] + delta_t_break[end] + delta_t_measurement[end]
        end
    end

    return EISCurrent(
        Float64(i_EIS),
        Float64(ratio),
        f,
        t0,
        t_new_start,
        tf,
        delta_t_break,
        delta_t_measurement
    )
end

# --- Interface implementation ---------------------------------------------

"""
    current(c::EISCurrent, t)

Compute the EIS current density at time `t`.
- For t < t0, a smooth ramp is applied to reach the base current.
- For t >= t0, the function determines the current frequency segment and applies a sinusoidal perturbation.
"""
function current(c::EISCurrent, t::Real)
    if t < c.t0
        # Smooth ramp to i_EIS (using tanh for a gradual increase)
        delta_t_ini = 3 * 60.0  # Duration of the initial ramp (s)
        return c.i_EIS *
            (1.0 + tanh(4 * (t - 2 * (delta_t_ini / 2)) / delta_t_ini)) / 2
    end

    # Find current frequency index (which frequency segment are we in?)
    n_inf = searchsortedlast(c.t_new_start, t)
    n_inf = max(1, n_inf)  # Ensure index is at least 1

    # Compute the sinusoidal perturbation for the current frequency
    i_disruption = (c.ratio * c.i_EIS) * cos(2 * π * c.f[n_inf] * t)

    return c.i_EIS + i_disruption
end

"""
    time_interval(c::EISCurrent)

Return the default simulation time interval `(t0, tf)`.
- t0 is always 0.0 (simulation start)
- tf is the end of the last frequency segment
"""
function time_interval(c::EISCurrent)
    return (0.0, c.tf)
end