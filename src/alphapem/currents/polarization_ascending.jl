# -*- coding: utf-8 -*-

"""This file is used to construct the ascending current sweep used by the extended local sensitivity analysis."""

# _____________________________________________Ascending polarization current___________________________________________

"""
    AscendingPolarizationCurrent

Current density model for one ascending sweep, without a return pass or hysteresis.

The profile first ramps from OCV to 1 A/cm² and stabilises there for `delta_t_ini`. It then ramps
back down to OCV and begins the requested ascending sweep. Every grid density is reached once, in
increasing order, and measured at the end of its stabilisation period. Transitions use the same
smooth hyperbolic tangent shape as the other polarization profiles.

# Fields
- `delta_t_ini::Float64`: Stabilisation time at 1 A/cm² before returning to OCV (s)
- `v_load::Float64`: Loading rate for initialization and sweep ramps (A/m²/s)
- `delta_t_break::Float64`: Stabilisation time above 0.3 A/cm² (s)
- `delta_t_break_lower_0_3::Float64`: Stabilisation time up to 0.3 A/cm², excluding OCV (s)
- `delta_t_break_OCV::Float64`: Stabilisation time at OCV after returning from 1 A/cm² (s)
- `delta_t_measurement::Float64`: Measurement time at the end of each stabilisation (s)
- `i_targets::Vector{Float64}`: Current densities to reach and measure, in increasing order (A/m²)
- `t_meas_starts::Vector{Float64}`: Start of the measurement window of each target (s)
- `t_meas_ends::Vector{Float64}`: End of the measurement window of each target (s)
- `di_transitions::Vector{Float64}`: Current variation of each ramp (A/m²)
- `t_starts::Vector{Float64}`: Start time of each ramp (s)
- `dt_loads::Vector{Float64}`: Duration of each ramp (s)
- `time_interval::Tuple{Float64, Float64}`: Simulation time interval (s)

The ramp to 1 A/cm² and the return ramp to OCV are also included in `di_transitions`. When the
first target is 0 A/m², it is measured during the OCV stabilisation before the ascending sweep.
"""
struct AscendingPolarizationCurrent <: AbstractCurrent
    delta_t_ini::Float64
    v_load::Float64
    delta_t_break::Float64
    delta_t_break_lower_0_3::Float64
    delta_t_break_OCV::Float64
    delta_t_measurement::Float64
    i_targets::Vector{Float64}
    t_meas_starts::Vector{Float64}
    t_meas_ends::Vector{Float64}
    di_transitions::Vector{Float64}
    t_starts::Vector{Float64}
    dt_loads::Vector{Float64}
    time_interval::Tuple{Float64, Float64}

    function AscendingPolarizationCurrent(p::AscendingPolarizationParams)
        # === Preliminary checks ===
        p.delta_t_ini ≥ 0 || throw(ArgumentError("delta_t_ini must be ≥ 0"))
        p.v_load > 0 || throw(ArgumentError("v_load must be > 0"))
        p.delta_t_break ≥ 0 || throw(ArgumentError("delta_t_break must be ≥ 0"))
        p.delta_t_break_lower_0_3 ≥ 0 || throw(ArgumentError("delta_t_break_lower_0_3 must be ≥ 0"))
        p.delta_t_break_OCV ≥ 0 || throw(ArgumentError("delta_t_break_OCV must be ≥ 0"))
        p.delta_t_measurement ≥ 0 || throw(ArgumentError("delta_t_measurement must be ≥ 0"))
        !isempty(p.i_grid) || throw(ArgumentError("i_grid must not be empty"))
        all(isfinite, p.i_grid) && all(i -> i ≥ 0.0, p.i_grid) ||
            throw(ArgumentError("i_grid must contain finite, non-negative values"))
        all(diff(p.i_grid) .> 0.0) ||
            throw(ArgumentError("i_grid must be strictly increasing"))
        iszero(first(p.i_grid)) && p.delta_t_break_OCV < p.delta_t_measurement &&
            throw(ArgumentError("delta_t_measurement cannot exceed delta_t_break_OCV for an OCV grid point"))

        # === Ramps, stabilisations and measurement windows ===
        t_meas_starts, t_meas_ends, di_transitions, t_starts, dt_loads, tf =
            _ascending_polarization_timeline(p)

        return new(
            Float64(p.delta_t_ini),
            Float64(p.v_load),
            Float64(p.delta_t_break),
            Float64(p.delta_t_break_lower_0_3),
            Float64(p.delta_t_break_OCV),
            Float64(p.delta_t_measurement),
            Float64.(p.i_grid),
            t_meas_starts,
            t_meas_ends,
            di_transitions,
            t_starts,
            dt_loads,
            (0.0, tf),
        )
    end
end


# --- Internal utilities ----------------------------------------------------

"""
    _ascending_break_time(i_fc, p)

Return the stabilisation time of a current density `i_fc` (A/m²) other than the first one.
"""
function _ascending_break_time(i_fc::Real, p::AscendingPolarizationParams)::Float64
    return i_fc / 1e4 <= 0.3 ? p.delta_t_break_lower_0_3 : p.delta_t_break
end

"""
Build the timeline of the initialization and ascending sweep.

The sequence is: ramp to 1 A/cm² and stabilize, ramp down to OCV and stabilize, then reach each
requested current density in ascending order. The voltage is measured over the last
`delta_t_measurement` seconds of each requested point's stabilisation. Returns
`(t_meas_starts, t_meas_ends, di_transitions, t_starts, dt_loads, tf)`.
"""
function _ascending_polarization_timeline(p::AscendingPolarizationParams)
    t_meas_starts = Float64[]
    t_meas_ends = Float64[]
    di_transitions = Float64[]
    t_starts = Float64[]
    dt_loads = Float64[]

    # Initialize the fuel cell at 1 A/cm² before the requested ascending sweep.
    i_initial = 1.0e4
    dt_initial_load = i_initial / p.v_load
    push!(di_transitions, i_initial)
    push!(t_starts, 0.0)
    push!(dt_loads, dt_initial_load)
    t_curr = dt_initial_load + p.delta_t_ini

    # Return to OCV and stabilize there before starting the ascending sweep.
    dt_return_to_ocv = i_initial / p.v_load
    push!(di_transitions, -i_initial)
    push!(t_starts, t_curr)
    push!(dt_loads, dt_return_to_ocv)
    t_curr += dt_return_to_ocv
    t_curr += p.delta_t_break_OCV

    # If OCV is requested, measure it during the hold immediately before the ascent.
    if iszero(first(p.i_grid))
        push!(t_meas_starts, t_curr - p.delta_t_measurement)
        push!(t_meas_ends, t_curr)
    end

    # Start the measured sweep at OCV.
    i_curr = 0.0
    for i_target in p.i_grid
        # OCV was already stabilized and measured before the ascending sweep.
        iszero(i_target) && continue

        # Ramp from the previous requested density to this ascending target.
        t_ramp_end = t_curr
        di = i_target - i_curr
        dt_load = di / p.v_load
        push!(di_transitions, di)
        push!(t_starts, t_curr)
        push!(dt_loads, dt_load)
        t_ramp_end = t_curr + dt_load

        # Stabilisation, whose last seconds form the measurement window.
        stabilisation = _ascending_break_time(i_target, p)
        stabilisation >= p.delta_t_measurement ||
            throw(ArgumentError(
                "delta_t_measurement ($(p.delta_t_measurement) s) cannot exceed the " *
                "stabilisation time ($(stabilisation) s) at $(i_target / 1e4) A/cm²"))
        t_stabilised = t_ramp_end + stabilisation
        push!(t_meas_starts, t_stabilised - p.delta_t_measurement)
        push!(t_meas_ends, t_stabilised)

        t_curr = t_stabilised
        i_curr = i_target
    end

    return t_meas_starts, t_meas_ends, di_transitions, t_starts, dt_loads, t_curr
end


# --- Interface implementation ---------------------------------------------

"""
    current(c::AscendingPolarizationCurrent, t)

Compute the ascending-sweep current density at time `t`.
"""
function current(c::AscendingPolarizationCurrent, t::Real)
    i_fc = 0.0
    for k in eachindex(c.di_transitions)
        i_fc += c.di_transitions[k] *
                (1.0 + tanh(4 * (t - c.t_starts[k] - (c.dt_loads[k] / 2)) /
                         (c.dt_loads[k] / 2))) / 2
    end
    return i_fc
end


"""
    solver_tstops(c::AscendingPolarizationCurrent, tspan)

Return the start time of each ramp as forced solver stop point.
"""
function solver_tstops(c::AscendingPolarizationCurrent, tspan::Tuple{<:Real, <:Real})::Vector{Float64}
    return _solver_tstops_in_range(c.t_starts, tspan)
end


"""
    time_interval(p::AscendingPolarizationParams)

Return the default simulation time interval `(t0, tf)`.
"""
function time_interval(p::AscendingPolarizationParams)
    return (0.0, _ascending_polarization_timeline(p)[end])
end
