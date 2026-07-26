# src/alphapem/currents/abstract.jl

"""
    AbstractCurrent

Abstract type representing a current density law.

All concrete implementations must subtype `AbstractCurrent` and implement:

    current(c::AbstractCurrent, t::Real)

which returns the current density at time `t`.

They should also implement:

    time_interval(c::AbstractCurrent)

which provides a default simulation time interval `(t0, tf)`.
"""
abstract type AbstractCurrent end


"""
    current(c::AbstractCurrent, t)

Compute the current density at time `t`.

# Arguments
- `c::AbstractCurrent`: Current density model
- `t::Real`: Time (s)

# Returns
- Current density (A/m²)
"""
function current(c::AbstractCurrent, t::Real)
    throw(MethodError(current, (c, t)))
end


"""
    current(c::AbstractCurrent, t::AbstractArray)

Vectorized evaluation of the current density.
"""
current(c::AbstractCurrent, t::AbstractArray) = current.(Ref(c), t)


"""
    solver_tstops(c::AbstractCurrent, tspan)

Return the solver stop times associated with a current profile.

Arguments
---------
- `tspan::Tuple{<:Real, <:Real}`: global simulation interval `(t0, tf)` used
  to filter candidate stop times.

Notes
-----
- `times` denotes a list of *candidate instants* generated from the current
  profile itself (for example, the beginnings of ramps).
- `tspan` denotes the *global time window* in which the solver is allowed to
  stop.  Only candidate times strictly inside this interval are kept.

The default implementation returns an empty vector, which means the solver
will use its own adaptive stepping with no forced stop points.
"""
solver_tstops(::AbstractCurrent, ::Tuple{<:Real, <:Real}) = Float64[]


"""
    solver_dtmax(c::AbstractCurrent, t::Real)

Return the maximum time step allowed for the solver at time `t`.
Returns `Inf` by default (no constraint).
"""
solver_dtmax(::AbstractCurrent, ::Real) = Inf


"""
    saveat_times(c::AbstractCurrent, tspan, save_freq::Float64) -> Vector{Float64}

Mandatory output-save instants for one IDA solve covering `tspan = (t0, tf)`,
passed as `saveat` to `solve()` and thus forcing the solver to stop exactly
there. Use this only when a current profile needs guaranteed samples
regardless of the solver's natural step size.

General density control (bounding the size of `sol.t`/`sol.u` to `save_freq`
points per second) is handled separately and automatically by a rate-limited
saving callback (see `_build_rate_limited_saving_callback` in `AlphaPEM.jl`),
which rides along on whatever steps IDA naturally takes instead of forcing
extra ones. `saveat_times` only adds points on top of that.

Default: no mandatory points. Override when dense, exact sampling is required
independent of the solver's step size — for example `EISCurrent`, which needs
`nb_points` per period during its stabilization/measurement windows for
accurate Fourier reconstruction.
"""
function saveat_times(::AbstractCurrent, ::Tuple{<:Real,<:Real}, ::Float64)::Vector{Float64}
    return Float64[]
end


"""
    _solver_tstops_in_range(times, tspan)

Filter a collection of candidate stop times.

- `times`: candidate event times produced by the current profile.
- `tspan`: global solver interval `(t0, tf)`.

Only times strictly inside `tspan` are kept, then they are sorted and
deduplicated. This is particularly useful for the :live display mode,
where the solver may be called with a time interval that is shorter
than the full simulation time interval.
"""
function _solver_tstops_in_range(times::AbstractVector{<:Real}, tspan::Tuple{<:Real, <:Real})::Vector{Float64}
    t0 = Float64(tspan[1])
    tf = Float64(tspan[2])
    stops = Float64[]

    # We use a small absolute tolerance to avoid stop times that are too close
    # to the interval boundaries. SUNDIALS IDA may fail if a tstop is
    # practically identical to the initial time t0.
    tol = 1e-9

    for t in times
        t_ft = Float64(t)
        if t0 + tol < t_ft < tf - tol
            push!(stops, t_ft)
        end
    end

    sort!(unique!(stops))
    return stops
end
