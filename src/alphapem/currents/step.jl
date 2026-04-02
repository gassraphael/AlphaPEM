# -*- coding: utf-8 -*-

"""This file is used to construct the step current density."""

# _____________________________________________________Preliminaries____________________________________________________


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
    time_interval::Tuple{Float64, Float64}

    function StepCurrent(p::StepParams)
        p.delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        p.delta_t_load > 0 || throw(ArgumentError("delta_t_load must be > 0"))
        p.delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))

        return new(
            Float64(p.delta_t_ini),
            Float64(p.delta_t_load),
            Float64(p.delta_t_break),
            Float64(p.i_ini),
            Float64(p.i_step),
            time_interval(p)
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
    time_interval(p::StepParams)

Compute the simulation time interval from the raw parameters, before the struct is built.
Used internally by the `StepCurrent` constructor.
"""
function time_interval(p::StepParams)
    return (0.0, p.delta_t_ini + p.delta_t_load + p.delta_t_break)
end