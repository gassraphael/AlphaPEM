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

Return the output-save instants for one IDA solve covering `tspan = (t0, tf)`.

These times are passed as `saveat` to `solve()` along with `save_everystep=false`
and `dense=false`, so that `sol.u` only grows with the number of returned points
regardless of how many internal IDA steps were taken.  This has two benefits:

1. **Memory** — `sol.u` stays small even for stiff trajectories that require
   hundreds of thousands of tiny steps.  Each entry of `sol.u` is a full state
   copy (~400 Float64 for a typical AlphaPEM system); without this bound, a run
   with 5×10⁵ internal steps accumulates ~1.6 GB and makes Julia's GC quadratically
   slower as the heap grows.

2. **GC pressure** — every accepted IDA step normally appends a new Julia object to
   `sol.u`.  With 5×10⁵ such objects the GC must scan them on every minor collection
   triggered by allocations inside `residual!`.  Limiting `sol.u` to a few thousand
   points eliminates this feedback loop.

Default implementation: sampling at the rate specified by `save_freq` (Hz).
This is sufficient for smooth profiles (step, polarization, calibration) and gives
a predictable, bounded output size independent of stiffness or `nb_gc`.

Concrete subtypes should override this method when the dynamics require a different
sampling rate — for example `EISCurrent`, which needs dense sampling during
measurement windows to accurately reconstruct the frequency response.

# Arguments
- `c::AbstractCurrent`: Current density model
- `tspan::Tuple{<:Real, <:Real}`: Time interval (t0, tf)
- `save_freq::Float64`: Output save frequency in Hz (points per second)
"""
function saveat_times(::AbstractCurrent, tspan::Tuple{<:Real,<:Real}, save_freq::Float64)::Vector{Float64}
    t0 = Float64(tspan[1])
    tf = Float64(tspan[2])
    dt = 1.0 / save_freq
    # Always include at least the boundary points so recovery! has t0 and tf.
    tf <= t0 && return Float64[t0]
    return collect(t0 : dt : tf)
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
