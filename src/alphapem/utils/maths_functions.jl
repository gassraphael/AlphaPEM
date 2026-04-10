# -*- coding: utf-8 -*-

"""This module contains mathematical functions which are used for modeling the PEM fuel cell."""

# ________________________________________________Mathematical functions________________________________________________

"""
    hmean(terms, weights=nothing)

Calculate the weighted harmonic mean of a list of terms with corresponding weights.
It is more efficient to express this function in the code than calling hmean from a library.

# Arguments
- `terms::AbstractVector`: The terms to calculate the harmonic mean for.
- `weights::Union{AbstractVector, Nothing}`: The weights corresponding to each term. If nothing, uniform weights are assumed.

# Returns
- The weighted harmonic mean.
"""
function hmean(terms::AbstractVector, weights::Union{AbstractVector, Nothing}=nothing)

    n = length(terms)
    # Calculate the weighted harmonic mean.
    weighted_sum = 0.0
    total_weight = 0.0
    if weights === nothing
        @inbounds for i in 1:n
            t = terms[i]
            if t != 0
                weighted_sum += 1.0 / t
            end
        end
        total_weight = Float64(n)
    else
        if length(weights) != n
            throw(ArgumentError("The length of terms and weights must be the same."))
        end
        @inbounds for i in 1:n
            w, t = weights[i], terms[i]
            if t != 0
                weighted_sum += w / t
            end
            total_weight += w
        end
    end

    if weighted_sum == 0
        throw(ArgumentError("All weights are zero in hmean calculation"))
    end

    return total_weight / weighted_sum
end


"""
    average(terms, weights=nothing)

Calculate the weighted arithmetic mean of a list of terms with corresponding weights.
It is more efficient to express this function in the code than calling average from a library.

# Arguments
- `terms::AbstractVector`: The terms to calculate the average for.
- `weights::Union{AbstractVectorVector, Nothing}`: The weights corresponding to each term. If nothing, uniform weights are assumed.

# Returns
- The weighted arithmetic mean.
"""
function average(terms::AbstractVector, weights::Union{AbstractVector, Nothing}=nothing)

    n = length(terms)
    if weights === nothing
        total_weight = Float64(n)
        weighted_sum = 0.0
        @inbounds for t in terms
            weighted_sum += t
        end
    else
        if n != length(weights)
            throw(ArgumentError("The length of terms and weights must be the same."))
        end
        total_weight = 0.0
        weighted_sum = 0.0
        @inbounds for i in 1:n
            w = weights[i]
            total_weight += w
            weighted_sum += w * terms[i]
        end
    end

    if total_weight == 0
        return NaN
    end
    return weighted_sum / total_weight
end


"""
    interpolate(terms, distances)

Fast inverse distance interpolation for exactly 2 points.

# Arguments
- `terms::AbstractVector`: The values at each node ([y1, y2]).
- `distances::AbstractVector`: The distances from each node to the interpolation point ([d1, d2]).

# Returns
- The interpolated value.
"""
function interpolate(terms::AbstractVector,
                     distances::AbstractVector)
    if length(terms) != 2 || length(distances) != 2
        throw(ArgumentError("This function only supports interpolation with 2 points."))
    end
    y1, y2 = terms
    d1, d2 = distances
    if d1 == 0 return y1 end
    if d2 == 0 return y2 end
    return (d2 * y1 + d1 * y2) / (d1 + d2)
end


"""
    d_dx(y_minus, y_plus, dx)

    d_dx(y_minus, y_plus, dx_minus, dx_plus)

Computes the centered first derivative (second order) with different steps to the left and right.

# Arguments
- `y_minus`: Value at the left point (i-1).
- `y_plus`: Value at the right point (i+1).
- `dx`: Half step between (i-1) and (i+1) when the spacing is uniform.
- `dx_minus`: Step between (i-1) and i.
- `dx_plus`: Step between i and (i+1).

# Returns
- `d_dx`: Approximation of the first derivative at i.
"""
function d_dx(y_minus, y_plus, dx)
    if dx == 0
        throw(ArgumentError("dx must be non-zero."))
    end
    return (y_plus - y_minus) / (2.0 * dx)
end


function d_dx(y_minus,
              y_plus,
              dx_minus,
              dx_plus)
    if dx_minus <= 0 || dx_plus <= 0
        throw(ArgumentError("dx_minus and dx_plus must be positive non-zero values."))
    end
    y_0 = interpolate([y_minus, y_plus], [dx_minus, dx_plus])
    return (y_plus * dx_minus^2 + y_0 * (dx_plus^2 - dx_minus^2) - y_minus * dx_plus^2) /
           (dx_minus * dx_plus * (dx_minus + dx_plus))
end
