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
    _solver_tstops_in_range(times, tspan)

Filter a collection of candidate stop times.

- `times`: candidate event times produced by the current profile.
- `tspan`: global solver interval `(t0, tf)`.

Only times strictly inside `tspan` are kept, then they are sorted and
deduplicated.
"""
function _solver_tstops_in_range(times::AbstractVector{<:Real}, tspan::Tuple{<:Real, <:Real})::Vector{Float64}
    t0 = Float64(tspan[1])
    tf = Float64(tspan[2])
    stops = Float64[]

    for t in times
        t_ft = Float64(t)
        if t0 < t_ft < tf
            push!(stops, t_ft)
        end
    end

    sort!(unique!(stops))
    return stops
end
