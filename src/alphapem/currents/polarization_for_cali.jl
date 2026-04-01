# -*- coding: utf-8 -*-

"""This file is used to construct the polarization current density used for the calibration algorithm."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing functions
include(joinpath(@__DIR__, "abstract.jl"))


# _____________________________________Polarization current density for calibration_____________________________________

"""
    PolarizationCalibrationCurrent

Current density model for polarization curve used in the calibration algorithm.
Based on experimental data.

The current follows experimental current density values with smooth transitions between each step.

# Fields
- `delta_t_ini::Float64`: Initial stabilization time (s)
- `v_load::Float64`: Loading rate (A/m²/s)
- `delta_t_break::Float64`: Stabilization time per step (s)
- `i_exp::Vector{Float64}`: Experimental current density values (A/m²)
"""
struct PolarizationCalibrationCurrent <: AbstractCurrent
    delta_t_ini::Float64
    v_load::Float64
    delta_t_break::Float64
    i_exp::Vector{Float64}
    time_interval::Tuple{Float64, Float64}

    function PolarizationCalibrationCurrent(p::PolarizationCalibrationParams)
        p.delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        p.v_load > 0 || throw(ArgumentError("v_load must be > 0"))
        p.delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))
        length(p.i_exp) > 0 || throw(ArgumentError("i_exp must not be empty"))

        return new(
            Float64(p.delta_t_ini),
            Float64(p.v_load),
            Float64(p.delta_t_break),
            Float64.(p.i_exp),
            time_interval(p)
        )
    end
end

# --- Internal utilities ----------------------------------------------------

"""
Returns a vector containing the loading duration for each current transition, calculated as abs(Δi) / v_load.
Δi = difference between two successive values of i_exp.
"""

function delta_t_load(c::PolarizationCalibrationCurrent)
    n = length(c.i_exp)
    dt_load = Vector{Float64}(undef, n)
    dt_load[1] = abs(c.i_exp[1]) / c.v_load
    for k in 2:n
        dt_load[k] = abs(c.i_exp[k] - c.i_exp[k-1]) / c.v_load
    end
    return dt_load
end
function delta_t_load(p::PolarizationCalibrationParams)
    n = length(p.i_exp)
    dt_load = Vector{Float64}(undef, n)
    dt_load[1] = abs(p.i_exp[1]) / p.v_load
    for k in 2:n
        dt_load[k] = abs(p.i_exp[k] - p.i_exp[k-1]) / p.v_load
    end
    return dt_load
end

"""
Returns a vector containing the total duration of each step (loading + stabilization).
"""

function step_duration(c::PolarizationCalibrationCurrent)
    dt_load = delta_t_load(c)
    return dt_load .+ c.delta_t_break
end
function step_duration(p::PolarizationCalibrationParams)
    dt_load = delta_t_load(p)
    return dt_load .+ p.delta_t_break
end

# --- Interface implementation ---------------------------------------------

"""
    current(c::PolarizationCalibrationCurrent, t)

Compute the calibration current density at time t.
"""
function current(c::PolarizationCalibrationCurrent, t::Real)
    dt_load = delta_t_load(c)
    dt_step = step_duration(c)
    i_fc = 0.0

    for e in eachindex(c.i_exp)
        if e == 1
            i_fc += c.i_exp[1] *
                    (1.0 + tanh(4 * (t - c.delta_t_ini - (dt_load[1] / 2)) /
                    (dt_load[1] / 2))) / 2
        else
            Δi = c.i_exp[e] - c.i_exp[e - 1]
            i_fc += Δi *
                   (1.0 + tanh(4 * (t - c.delta_t_ini - sum(dt_step[1:(e-1)]) - (dt_load[e] / 2)) /
                   (dt_load[e] / 2))) / 2
        end
    end
    return i_fc
end

"""
    time_interval(c::PolarizationCalibrationParams)

Return the default simulation time interval `(t0, tf)`.

# Returns
- `(0.0, tf)` where:
    tf = delta_t_ini + sum of the durations of each step (loading + stabilization)
"""
function time_interval(p::PolarizationCalibrationParams)
    return (0.0, p.delta_t_ini + sum(step_duration(p)))
end