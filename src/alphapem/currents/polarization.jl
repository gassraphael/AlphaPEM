# -*- coding: utf-8 -*-

"""This file is used to construct the polarization current density."""

# _____________________________________________________Preliminaries____________________________________________________


# _____________________________________________Polarization current density_____________________________________________

"""
    PolarizationCurrent

Current density model used to generate polarization curves.

The current follows a complete cycle:
1. Initial jump from 0 to 1.0 A/cm² followed by a long stabilization period (`delta_t_ini`).
2. Stepwise increase from 1.0 A/cm² up to `i_max`.
3. Stepwise decrease from `i_max` down to 0.
4. Stepwise increase from 0 back up to 1.0 A/cm².

Transitions are based on smooth hyperbolic tangent functions. Each increment
is followed by a stabilization period (`delta_t_break`), which varies based on current density.

# Fields
- `delta_t_ini::Float64`: Initial stabilization time (s)
- `v_load::Float64`: Loading rate for one step current (A/m²/s)
- `delta_t_break::Float64`: Standard stabilization time after each increment (s)
- `delta_t_break_lower_0_3::Float64`: Stabilization time when i_fc <= 0.3 A/cm² (s)
- `delta_t_break_OCV::Float64`: Stabilization time when i_fc = 0 A/cm² (s)
- `delta_t_measurement::Float64`: Measurement time at the end of each increment (s)
- `di_step::Float64`: Nominal current increment per step (A/m²)
- `i_max::Float64`: Maximum current density (A/m²)
- `di_transitions::Vector{Float64}`: Sequence of actual current variations (A/m²)
- `t_starts::Vector{Float64}`: Start times for each increment (s)
- `dt_loads::Vector{Float64}`: Duration of each loading phase (s)
- `delta_t_breaks::Vector{Float64}`: Actual stabilization times used for each step (s)
- `time_interval::Tuple{Float64, Float64}`: Simulation time interval (s)
"""
struct PolarizationCurrent <: AbstractCurrent
    delta_t_ini::Float64
    i_ini::Float64
    v_load::Float64
    delta_t_break::Float64
    delta_t_break_lower_0_3::Float64
    delta_t_break_OCV::Float64
    delta_t_measurement::Float64
    di_step::Float64
    i_max::Float64
    di_transitions::Vector{Float64}
    t_starts::Vector{Float64}
    dt_loads::Vector{Float64}
    delta_t_breaks::Vector{Float64}
    time_interval::Tuple{Float64, Float64}

    function PolarizationCurrent(p::PolarizationParams)
        p.delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        p.i_ini ≥ 0 || throw(ArgumentError("i_ini must be ≥ 0"))
        p.v_load > 0 || throw(ArgumentError("v_load must be > 0"))
        p.delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))
        p.delta_t_break_lower_0_3 ≥ 0 || throw(ArgumentError("delta_t_break_lower_0_3 must be ≥ 0"))
        p.delta_t_break_OCV ≥ 0 || throw(ArgumentError("delta_t_break_OCV must be ≥ 0"))
        p.delta_t_measurement ≥ 0 || throw(ArgumentError("delta_t_measurement must be ≥ 0"))
        p.di_step > 0 || throw(ArgumentError("di_step must be > 0"))
        p.i_max ≥ 0 || throw(ArgumentError("i_max must be ≥ 0"))

        di_transitions, t_starts, dt_loads, delta_t_breaks, tf, i_max_rounded = _polarization_transitions(p)

        return new(
            Float64(p.delta_t_ini),
            Float64(p.i_ini),
            Float64(p.v_load),
            Float64(p.delta_t_break),
            Float64(p.delta_t_break_lower_0_3),
            Float64(p.delta_t_break_OCV),
            Float64(p.delta_t_measurement),
            Float64(p.di_step),
            Float64(i_max_rounded),
            di_transitions,
            t_starts,
            dt_loads,
            delta_t_breaks,
            (0.0, tf)
        )
    end
end



# --- Internal utilities ----------------------------------------------------

"""
    _round_i_max_to_di_step(i_max, di_step)

Round i_max up to the nearest multiple of di_step to ensure symmetric ramp steps.
"""
function _round_i_max_to_di_step(i_max::Real, di_step::Real)::Float64
    n_steps = ceil(i_max / di_step)
    return n_steps * di_step
end

"""
    _select_break_time(i_fc, p)

Select the appropriate break time based on current density.
"""
function _select_break_time(i_fc::Real, p::PolarizationParams)::Float64
    i_fc_cm2 = i_fc / 1e4  # Convert from A/m² to A/cm²
    if abs(i_fc_cm2) < 1e-6  # At OCV (i_fc = 0)
        return p.delta_t_break_OCV
    elseif i_fc_cm2 <= 0.3
        return p.delta_t_break_lower_0_3
    else
        return p.delta_t_break
    end
end

"""
Generate the sequence of current transitions for the polarization cycle:
0 -> i_start (stabilization) -> i_max -> 0 -> i_start.
Returns (transitions, t_starts, dt_loads, delta_t_breaks, tf, i_max_rounded) to ensure symmetric stepping.
"""
function _polarization_transitions(p::PolarizationParams)
    i_max = _round_i_max_to_di_step(p.i_max, p.di_step)

    di_transitions = Float64[]
    t_starts = Float64[]
    dt_loads = Float64[]
    delta_t_breaks = Float64[]

    # 1. Initial loading: 0 -> i_ini (no measurement)
    push!(di_transitions, p.i_ini)
    push!(t_starts, 0.0)
    dt_load_ini = p.i_ini / p.v_load
    push!(dt_loads, dt_load_ini)
    push!(delta_t_breaks, p.delta_t_ini)  # Stabilization time after initial loading
    t_curr = dt_load_ini + p.delta_t_ini
    i_curr = p.i_ini

    # 2. Single large ramp: i_ini -> i_max (no measurement during ramp)
    # Note: i_max can be lower than i_ini (e.g. :before_voltage_drop zones), in which
    # case this is a ramp DOWN. The duration must stay positive in both directions.
    di_to_max = i_max - p.i_ini
    push!(di_transitions, di_to_max)
    push!(t_starts, t_curr)
    dt = abs(di_to_max) / p.v_load
    push!(dt_loads, dt)
    dt_break = _select_break_time(i_max, p)
    push!(delta_t_breaks, dt_break)
    t_curr += dt + dt_break
    i_curr = i_max

    # 3. Ramp down: i_max -> 0 (measurement steps)
    while i_curr > 1e-6
        di = min(p.di_step, i_curr)
        push!(di_transitions, -di)
        push!(t_starts, t_curr)
        dt = di / p.v_load
        push!(dt_loads, dt)
        i_next = i_curr - di
        dt_break = _select_break_time(i_next, p)
        push!(delta_t_breaks, dt_break)
        t_curr += dt + dt_break
        i_curr = i_next
    end

    # 4. Ramp back: 0 -> i_max (measurement steps)
    while i_curr < i_max - 1e-6
        di = min(p.di_step, i_max - i_curr)
        push!(di_transitions, di)
        push!(t_starts, t_curr)
        dt = di / p.v_load
        push!(dt_loads, dt)
        i_next = i_curr + di
        dt_break = _select_break_time(i_next, p)
        push!(delta_t_breaks, dt_break)
        t_curr += dt + dt_break
        i_curr = i_next
    end

    return di_transitions, t_starts, dt_loads, delta_t_breaks, t_curr, i_max
end


"""
Number of current increments required to reach `i_max`.
"""
n_steps(p::PolarizationParams) = length(_polarization_transitions(p)[1])  # Returns transitions as first element
n_steps(c::PolarizationCurrent) = length(c.di_transitions)



"""
Total duration of one increment cycle (load + break).
Returns a vector of durations for each transition in the cycle.
"""
function step_duration(c::PolarizationCurrent)
    dts = Vector{Float64}(undef, length(c.di_transitions))
    for k in eachindex(c.di_transitions)
        dts[k] = c.dt_loads[k] + (k == 1 ? c.delta_t_ini : c.delta_t_breaks[k])
    end
    return dts
end
function step_duration(p::PolarizationParams)
    di_transitions, t_starts, dt_loads, delta_t_breaks, tf, _ = _polarization_transitions(p)
    dts = Vector{Float64}(undef, length(di_transitions))
    for k in eachindex(di_transitions)
        dts[k] = dt_loads[k] + (k == 1 ? p.delta_t_ini : delta_t_breaks[k])
    end
    return dts
end


# --- Interface implementation ---------------------------------------------

"""
    current(c::PolarizationCurrent, t)

Compute the polarization current density at time `t`.
"""
function current(c::PolarizationCurrent, t::Real)
    i_fc = 0.0
    for k in eachindex(c.di_transitions)
        i_fc += c.di_transitions[k] *
                (1.0 + tanh(4 * (t - c.t_starts[k] - (c.dt_loads[k] / 2)) /
                         (c.dt_loads[k] / 2))) / 2
    end
    return i_fc
end


"""
    solver_tstops(c::PolarizationCurrent, tspan)

Return the stop times associated with each smooth current increment.

- `times` are the candidate absolute instants generated by the profile itself
  (here: the beginning of each ramp).
- `tspan` is the global simulation interval `(t0, tf)` used to keep only the
  candidates that fall inside the solve window.

For the polarization profile, only the beginning of each ramp is retained as a
forced stop point.
"""
function solver_tstops(c::PolarizationCurrent, tspan::Tuple{<:Real, <:Real})::Vector{Float64}
    return _solver_tstops_in_range(c.t_starts, tspan)
end


"""
    time_interval(c::PolarizationParams)

Return the default simulation time interval `(t0, tf)`.
"""
function time_interval(p::PolarizationParams)
    _, _, _, _, tf, _ = _polarization_transitions(p)
    return (0.0, tf)
end