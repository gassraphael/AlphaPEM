# -*- coding: utf-8 -*-

"""This file is used to construct the polarization current density used for the calibration algorithm."""

# _____________________________________________________Preliminaries____________________________________________________


# _____________________________________Polarization current density for calibration_____________________________________

"""
    PolarizationCalibrationCurrent

Current density model for polarization curve used in the calibration algorithm.
Based on experimental data.

The current follows a complete cycle using experimental current density values:
1. Initial jump from 0 to 1.0 A/cm² followed by a long stabilization period (`delta_t_ini`).
2. Stepwise increase from 1.0 A/cm² up to `i_max`.
3. Stepwise decrease from `i_max` down to 0.
4. Stepwise increase from 0 back up to 1.0 A/cm².

Transitions are smooth between each step.

# Fields
- `delta_t_ini::Float64`: Initial stabilization time (s)
- `v_load::Float64`: Loading rate (A/m²/s)
- `delta_t_break::Float64`: Stabilization time per step (s)
- `i_exp::Vector{Float64}`: Experimental current density values (A/m²)
- `di_transitions::Vector{Float64}`: Sequence of actual current variations (A/m²)
"""
struct PolarizationCalibrationCurrent <: AbstractCurrent
    delta_t_ini::Float64
    v_load::Float64
    delta_t_break::Float64
    i_exp::Vector{Float64}
    di_transitions::Vector{Float64}
    t_starts::Vector{Float64}
    dt_loads::Vector{Float64}
    time_interval::Tuple{Float64, Float64}

    function PolarizationCalibrationCurrent(p::PolarizationCalibrationParams)
        p.delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        p.v_load > 0 || throw(ArgumentError("v_load must be > 0"))
        p.delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))
        length(p.i_exp) > 0 || throw(ArgumentError("i_exp must not be empty"))

        di_transitions, t_starts, dt_loads, tf = _polarization_cali_transitions(p)

        return new(
            Float64(p.delta_t_ini),
            Float64(p.v_load),
            Float64(p.delta_t_break),
            Float64.(p.i_exp),
            di_transitions,
            t_starts,
            dt_loads,
            (0.0, tf)
        )
    end
end

# --- Internal utilities ----------------------------------------------------

"""
Generate the sequence of current transitions for the calibration polarization cycle:
0 -> 1.0 A/cm² (stabilization) -> i_max -> 0 -> 1.0 A/cm².
Points from `i_exp` are used as steps.
"""
function _polarization_cali_transitions(p::PolarizationCalibrationParams)
    i_start = 1.0e4 # 1.0 A/cm²
    i_ref = sort(unique(p.i_exp))
    i_max = isempty(i_ref) ? i_start : i_ref[end]

    di_transitions = Float64[]
    t_starts = Float64[]
    dt_loads = Float64[]

    # 1. Initial Step: 0 -> i_start
    push!(di_transitions, i_start)
    push!(t_starts, 0.0)
    dt_load_ini = i_start / p.v_load
    push!(dt_loads, dt_load_ini)

    t_curr = dt_load_ini + p.delta_t_ini
    i_curr = i_start

    # 2. Ramp up: i_start -> i_max
    for i_target in i_ref
        if i_target > i_curr + 1e-6
            di = i_target - i_curr
            push!(di_transitions, di)
            push!(t_starts, t_curr)
            dt = di / p.v_load
            push!(dt_loads, dt)
            t_curr += dt + p.delta_t_break
            i_curr = i_target
        end
    end

    # 3. Ramp down: i_max -> 0
    i_ref_desc = reverse(i_ref)
    for i_target in i_ref_desc
        if i_target < i_curr - 1e-6
            di = i_curr - i_target
            push!(di_transitions, -di)
            push!(t_starts, t_curr)
            dt = di / p.v_load
            push!(dt_loads, dt)
            t_curr += dt + p.delta_t_break
            i_curr = i_target
        end
    end
    if i_curr > 1e-6
        di = i_curr
        push!(di_transitions, -di)
        push!(t_starts, t_curr)
        dt = di / p.v_load
        push!(dt_loads, dt)
        t_curr += dt + p.delta_t_break
        i_curr = 0.0
    end

    # 4. Ramp back: 0 -> i_start
    for i_target in i_ref
        if i_target > i_curr + 1e-6 && i_target <= i_start + 1e-6
            di = i_target - i_curr
            push!(di_transitions, di)
            push!(t_starts, t_curr)
            dt = di / p.v_load
            push!(dt_loads, dt)
            t_curr += dt + p.delta_t_break
            i_curr = i_target
        end
    end
    if i_curr < i_start - 1e-6
        di = i_start - i_curr
        push!(di_transitions, di)
        push!(t_starts, t_curr)
        dt = di / p.v_load
        push!(dt_loads, dt)
        t_curr += dt + p.delta_t_break
        i_curr = i_start
    end

    return di_transitions, t_starts, dt_loads, t_curr
end

"""
Returns a vector containing the total duration of each step (loading + stabilization).
"""
function step_duration(c::PolarizationCalibrationCurrent)
    dts = Vector{Float64}(undef, length(c.di_transitions))
    for k in eachindex(c.di_transitions)
        dts[k] = c.dt_loads[k] + (k == 1 ? c.delta_t_ini : c.delta_t_break)
    end
    return dts
end
function step_duration(p::PolarizationCalibrationParams)
    di_transitions, t_starts, dt_loads, tf = _polarization_cali_transitions(p)
    dts = Vector{Float64}(undef, length(di_transitions))
    for k in eachindex(di_transitions)
        dts[k] = dt_loads[k] + (k == 1 ? p.delta_t_ini : p.delta_t_break)
    end
    return dts
end

# --- Interface implementation ---------------------------------------------

"""
    current(c::PolarizationCalibrationCurrent, t)

Compute the calibration current density at time t.
"""
function current(c::PolarizationCalibrationCurrent, t::Real)
    i_fc = 0.0
    for k in eachindex(c.di_transitions)
        i_fc += c.di_transitions[k] *
                (1.0 + tanh(4 * (t - c.t_starts[k] - (c.dt_loads[k] / 2)) /
                (c.dt_loads[k] / 2))) / 2
    end
    return i_fc
end


"""
    solver_tstops(c::PolarizationCalibrationCurrent, tspan)

Return the stop times associated with each experimental loading transition.
"""
function solver_tstops(c::PolarizationCalibrationCurrent, tspan::Tuple{<:Real, <:Real})::Vector{Float64}
    return _solver_tstops_in_range(c.t_starts, tspan)
end

"""
    time_interval(c::PolarizationCalibrationParams)

Return the default simulation time interval `(t0, tf)`.
"""
function time_interval(p::PolarizationCalibrationParams)
    _, _, _, tf = _polarization_cali_transitions(p)
    return (0.0, tf)
end