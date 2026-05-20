# -*- coding: utf-8 -*-

"""Small helper functions used to protect physical and numerical calculations.

These helpers only sanitize or bound inputs before they are passed to the actual
models, so they live in the shared utility layer.
"""

@inline _positive_temperature_value(T::Real) = max(Float64(T), 1.0)
@inline _liquid_water_temperature_value(T::Real) = clamp(Float64(T), 140.0 + 1e-6, 647.15 - 1e-6)
@inline _bounded_saturation_value(s::Real) = clamp(Float64(s), 1e-9, 1.0 - 1e-9)
@inline _positive_pressure_value(P::Real) = max(Float64(P), 1.0)
@inline _nonnegative_value(x::Real) = max(Float64(x), eps(Float64))

"""
    _positive_concentration_value(x)

Lower-bound a species concentration to 1e-4 mol/m³.
Use this when the concentration appears inside a logarithm to avoid Newton/Jacobian blow-up.
"""
@inline _positive_concentration_value(x::Real) = max(Float64(x), 1e-4)

function _safe_porous_phase_weights(epsilon::Float64, s)
    s_eff = _bounded_saturation_value(s)
    return [max(1 - epsilon, 0.0), max(epsilon * s_eff, 0.0), max(epsilon * (1 - s_eff), 0.0)]
end

function _safe_cl_phase_weights(epsilon_cl_val::Float64, epsilon_mc_val::Float64, s)
    s_eff = _bounded_saturation_value(s)
    return [
        max(1 - epsilon_cl_val - epsilon_mc_val, 0.0),
        max(epsilon_mc_val, 0.0),
        max(epsilon_cl_val * s_eff, 0.0),
        max(epsilon_cl_val * (1 - s_eff), 0.0),
    ]
end

