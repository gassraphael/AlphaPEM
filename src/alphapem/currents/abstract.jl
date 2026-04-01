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