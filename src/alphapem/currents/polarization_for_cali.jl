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

Transitions are smooth between each step. Stabilization times vary based on current density.

# Fields
- `delta_t_ini::Float64`: Initial stabilization time (s)
- `v_load::Float64`: Loading rate (A/m²/s)
- `delta_t_break::Float64`: Standard stabilization time per step (s)
- `delta_t_break_lower_0_3::Float64`: Stabilization time when i_fc <= 0.3 A/cm² (s)
- `delta_t_break_OCV::Float64`: Stabilization time when i_fc = 0 A/cm² (s)
- `delta_t_measurement::Float64`: Measurement time at the end of each step (s)
- `i_exp::Vector{Float64}`: Experimental current density values (A/m²)
- `di_transitions::Vector{Float64}`: Sequence of actual current variations (A/m²)
- `t_starts::Vector{Float64}`: Start times for each increment (s)
- `dt_loads::Vector{Float64}`: Duration of each loading phase (s)
- `delta_t_breaks::Vector{Float64}`: Actual stabilization times used for each step (s)
- `time_interval::Tuple{Float64, Float64}`: Simulation time interval (s)
"""
struct PolarizationCalibrationCurrent <: AbstractCurrent
    delta_t_ini::Float64
    i_ini::Float64
    v_load::Float64
    delta_t_break::Float64
    delta_t_break_lower_0_3::Float64
    delta_t_break_OCV::Float64
    delta_t_measurement::Float64
    i_exp::Vector{Float64}
    di_transitions::Vector{Float64}
    t_starts::Vector{Float64}
    dt_loads::Vector{Float64}
    delta_t_breaks::Vector{Float64}
    time_interval::Tuple{Float64, Float64}

    function PolarizationCalibrationCurrent(p::PolarizationCalibrationParams)
        p.delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        p.i_ini ≥ 0 || throw(ArgumentError("i_ini must be ≥ 0"))
        p.v_load > 0 || throw(ArgumentError("v_load must be > 0"))
        p.delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))
        p.delta_t_break_lower_0_3 ≥ 0 || throw(ArgumentError("delta_t_break_lower_0_3 must be ≥ 0"))
        p.delta_t_break_OCV ≥ 0 || throw(ArgumentError("delta_t_break_OCV must be ≥ 0"))
        p.delta_t_measurement ≥ 0 || throw(ArgumentError("delta_t_measurement must be ≥ 0"))
        length(p.i_exp) > 0 || throw(ArgumentError("i_exp must not be empty"))

        di_transitions, t_starts, dt_loads, delta_t_breaks, tf = _polarization_cali_transitions(p)

        return new(
            Float64(p.delta_t_ini),
            Float64(p.i_ini),
            Float64(p.v_load),
            Float64(p.delta_t_break),
            Float64(p.delta_t_break_lower_0_3),
            Float64(p.delta_t_break_OCV),
            Float64(p.delta_t_measurement),
            Float64.(p.i_exp),
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
    _select_break_time(i_fc, p)

Select the appropriate break time based on current density.
"""
function _select_break_time(i_fc::Real, p::PolarizationCalibrationParams)::Float64
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
Generate the sequence of current transitions for the calibration polarization cycle:
0 -> i_ini (initialization, no measurement) -> i_max -> 0 -> i_max (measurement steps using i_exp).
Points from `i_exp` are used as steps in descent and ascent.
"""
function _polarization_cali_transitions(p::PolarizationCalibrationParams)
    i_ref = sort(unique(p.i_exp))
    i_max = isempty(i_ref) ? p.i_ini : i_ref[end]

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
    push!(delta_t_breaks, _select_break_time(i_max, p))
    t_curr += dt + delta_t_breaks[2]
    i_curr = i_max

    # 3. Ramp down: i_max -> 0 (measurement steps)
    i_ref_desc = reverse(i_ref)
    for i_target in i_ref_desc
        if i_target < i_curr - 1e-6
            di = i_curr - i_target
            push!(di_transitions, -di)
            push!(t_starts, t_curr)
            dt = di / p.v_load
            push!(dt_loads, dt)
            dt_break = _select_break_time(i_target, p)
            push!(delta_t_breaks, dt_break)
            t_curr += dt + dt_break
            i_curr = i_target
        end
    end
    if i_curr > 1e-6
        di = i_curr
        push!(di_transitions, -di)
        push!(t_starts, t_curr)
        dt = di / p.v_load
        push!(dt_loads, dt)
        dt_break = _select_break_time(0.0, p)
        push!(delta_t_breaks, dt_break)
        t_curr += dt + dt_break
        i_curr = 0.0
    end

    # 4. Ramp back: 0 -> i_max (measurement steps)
    for i_target in i_ref
        if i_target > i_curr + 1e-6
            di = i_target - i_curr
            push!(di_transitions, di)
            push!(t_starts, t_curr)
            dt = di / p.v_load
            push!(dt_loads, dt)
            dt_break = _select_break_time(i_target, p)
            push!(delta_t_breaks, dt_break)
            t_curr += dt + dt_break
            i_curr = i_target
        end
    end

    return di_transitions, t_starts, dt_loads, delta_t_breaks, t_curr
end

"""
Returns a vector containing the total duration of each step (loading + stabilization).
"""
function step_duration(c::PolarizationCalibrationCurrent)
    dts = Vector{Float64}(undef, length(c.di_transitions))
    for k in eachindex(c.di_transitions)
        dts[k] = c.dt_loads[k] + c.delta_t_breaks[k]
    end
    return dts
end
function step_duration(p::PolarizationCalibrationParams)
    di_transitions, t_starts, dt_loads, delta_t_breaks, tf = _polarization_cali_transitions(p)
    dts = Vector{Float64}(undef, length(di_transitions))
    for k in eachindex(di_transitions)
        dts[k] = dt_loads[k] + delta_t_breaks[k]
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
    _, _, _, _, tf = _polarization_cali_transitions(p)
    return (0.0, tf)
end