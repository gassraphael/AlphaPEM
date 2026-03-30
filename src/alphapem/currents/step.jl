# -*- coding: utf-8 -*-

"""This file is used to construct the step current density."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing functions
include(joinpath(@__DIR__, "abstract.jl"))


# _________________________________________________Step current density_________________________________________________


"""
    StepCurrent

Step current density model with a smooth transition using hyperbolic tangent functions.
For the first delta_t_ini_step seconds, the current density is set to i_ini A.m-2 to allow the internal states of the
fuel cell to stabilise.
Then, the current density increases from i_ini to i_step A.m-2 in a step change over delta_t_load seconds.
Finally, the current density remains at i_step A.m-2.
"""
struct StepCurrent <: AbstractCurrent
    delta_t_ini::Float64
    delta_t_load::Float64
    delta_t_break::Float64
    i_ini::Float64
    i_step::Float64

    function StepCurrent(
        delta_t_ini::Real,
        delta_t_load::Real,
        delta_t_break::Real,
        i_ini::Real,
        i_step::Real
    )
        delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        delta_t_load > 0 || throw(ArgumentError("delta_t_load must be > 0"))
        delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))

        return new(
            Float64(delta_t_ini),
            Float64(delta_t_load),
            Float64(delta_t_break),
            Float64(i_ini),
            Float64(i_step)
        )
    end
end


"""
    current(c::StepCurrent, t)

Compute the current density at time `t`.
"""
function current(c::StepCurrent, t::Real)

    term1 = c.i_ini *
        (1.0 + tanh(4 * (t - (c.delta_t_load / 2)) / (c.delta_t_load / 2))) / 2

    term2 = (c.i_step - c.i_ini) *
        (1.0 + tanh(4 * (t - c.delta_t_ini - (c.delta_t_load / 2)) /
                     (c.delta_t_load / 2))) / 2

    return term1 + term2
end


"""
    time_interval(c::StepCurrent)

Return the default simulation time interval `(t0, tf)`.

# Returns
- `(0.0, tf)` where:
    tf = delta_t_ini + delta_t_load + delta_t_break
"""
function time_interval(c::StepCurrent)
    t0 = 0.0
    tf = c.delta_t_ini + c.delta_t_load + c.delta_t_break
    return (t0, tf)
end