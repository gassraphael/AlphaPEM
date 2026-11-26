# -*- coding: utf-8 -*-

"""
    _enforce_integer_params(X, params)

Round and clamp the rows of sample matrix `X` that belong to integer parameters.

Parameters
----------
X : Matrix{Float64}
    Samples, one parameter per row and one sample per column.
params : AbstractVector
    Parameter descriptions in row order, each with the fields `type`, `min` and `max`
    (`InputParameter` or `ParameterBound`).

Returns
-------
Matrix{Float64}
    Copy of `X` with integer parameters rounded and clamped to their bounds.
"""
function _enforce_integer_params(X::Matrix{Float64}, params::AbstractVector)::Matrix{Float64}
    X_out = copy(X)
    for (j, p) in enumerate(params)
        if p.type == :int
            for i in 1:size(X_out, 2)
                X_out[j, i] = clamp(round(X_out[j, i]), p.min, p.max)
            end
        end
    end
    return X_out
end


"""
    _physical_params_with_overrides(base, params, values) -> PhysicalParams

Return `base` with the analysed parameters replaced by `values`, all other fields unchanged.

Parameters
----------
base : PhysicalParams
    Nominal physical parameters.
params : AbstractVector
    Analysed parameters, each with the fields `name` and `type`.
values : AbstractVector{<:Real}
    Physical values in the order of `params`; integer parameters are rounded.

Returns
-------
PhysicalParams
    New parameter set.
"""
function _physical_params_with_overrides(base::PhysicalParams, params::AbstractVector,
                                         values::AbstractVector{<:Real})::PhysicalParams
    base_fields = (; (f => getfield(base, f) for f in fieldnames(PhysicalParams))...)
    overrides = (; (p.name => (p.type == :int ? Int(round(v)) : Float64(v))
                    for (p, v) in zip(params, values))...)
    return PhysicalParams(; merge(base_fields, overrides)...)
end
