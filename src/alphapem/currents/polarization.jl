# -*- coding: utf-8 -*-

"""This file is used to construct the polarization current density."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing functions
include(joinpath(@__DIR__, "abstract.jl"))


# _____________________________________________Polarization current density_____________________________________________

"""
    PolarizationCurrent

Current density model used to generate polarization curves.

The current increases stepwise from 0 up to `i_max` using smooth transitions
based on hyperbolic tangent functions. Each increment is followed by a
stabilization period.

# Fields
- `delta_t_ini::Float64`: Initial stabilization time (s)
- `v_load::Float64`: Loading rate for one step current (A/m²/s)
- `delta_t_break::Float64`: Stabilization time after each increment (s)
- `delta_i::Float64`: Current increment (A/m²)
- `i_max::Float64`: Maximum current density (A/m²)
"""
struct PolarizationCurrent <: AbstractCurrent
    delta_t_ini::Float64
    v_load::Float64
    delta_t_break::Float64
    delta_i::Float64
    i_max::Float64

    function PolarizationCurrent(p::PolarizationParams)
        p.delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        p.v_load > 0 || throw(ArgumentError("v_load must be > 0"))
        p.delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))
        p.delta_i > 0 || throw(ArgumentError("delta_i must be > 0"))
        p.i_max ≥ 0 || throw(ArgumentError("i_max must be ≥ 0"))

        return new(
            Float64(p.delta_t_ini),
            Float64(p.v_load),
            Float64(p.delta_t_break),
            Float64(p.delta_i),
            Float64(p.i_max)
        )
    end
end



# --- Internal utilities ----------------------------------------------------

"""
    n_steps(c::PolarizationCurrent)

Number of current increments required to reach `i_max`.
"""
n_steps(c::PolarizationCurrent) = Int(floor(c.i_max / c.delta_i))

"""
delta_t_load(c::PolarizationCurrent)

Loading time for one step current.
"""
delta_t_load(c::PolarizationCurrent) = c.delta_i / c.v_load


"""
step_duration(c::PolarizationCurrent)

Total duration of one increment cycle (load + break).
"""
step_duration(c::PolarizationCurrent) = delta_t_load(c) + c.delta_t_break


# --- Interface implementation ---------------------------------------------

"""
    current(c::PolarizationCurrent, t)

Compute the polarization current density at time `t`.
"""
function current(c::PolarizationCurrent, t::Real)

    dt_step = step_duration(c)
    dt_load = delta_t_load(c)
    i_fc = 0.0

    for k in 0:(n_steps(c) - 1)
        i_fc += c.delta_i *
                (1.0 + tanh(4 * (t - c.delta_t_ini - k * dt_step - (dt_load / 2)) /
                         (dt_load / 2))) / 2
    end

    return i_fc
end


"""
    time_interval(c::PolarizationCurrent)

Return the default simulation time interval `(t0, tf)`.

# Returns
- `(0.0, tf)` where:
    tf = delta_t_ini + n_steps * (delta_t_load + delta_t_break)
"""
function time_interval(c::PolarizationCurrent)
    t0 = 0.0
    tf = c.delta_t_ini + n_steps(c) * step_duration(c)
    return (t0, tf)
end