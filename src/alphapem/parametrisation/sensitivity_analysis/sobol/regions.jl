# -*- coding: utf-8 -*-

"""
    build_regions(thresholds)

Build the three characteristic polarization regions from current-density thresholds.

`thresholds` is a tuple `(t1, t2)` in A/cm². The regions are closed on both sides:
- activation:      [0, t1]
- ohmic:           [t1, t2]
- mass_transport:  [t2, ∞)

Returns a vector of `PolarizationRegion`.
"""
function build_regions(thresholds::Tuple{Float64, Float64})::Vector{PolarizationRegion}
    t1, t2 = thresholds
    0.0 < t1 < t2 || throw(ArgumentError("Region thresholds must satisfy 0 < t1 < t2 (got $thresholds)"))
    return [
        PolarizationRegion(:activation,     0.0, t1),
        PolarizationRegion(:ohmic,          t1,  t2),
        PolarizationRegion(:mass_transport, t2,  Inf),
    ]
end


"""
    region_auc(ifc, Ucell, region; clip_negative=true)

Compute the area under the polarization curve `U = f(i)` inside a `PolarizationRegion`.

Arguments:
- `ifc::Vector{Float64}`: Current densities in A/cm².
- `Ucell::Vector{Float64}`: Cell voltages in V.
- `region::PolarizationRegion`: Region of interest.
- `clip_negative::Bool`: If `true`, negative voltages are clipped to 0 before integration.

Returns the AUC in V·A/cm². If no point falls inside the region, returns `NaN`.
"""
function region_auc(ifc::Vector{Float64},
                    Ucell::Vector{Float64},
                    region::PolarizationRegion;
                    clip_negative::Bool = true)::Float64
    length(ifc) == length(Ucell) || throw(ArgumentError("ifc and Ucell must have the same length"))
    isempty(ifc) && return NaN

    U = clip_negative ? max.(Ucell, 0.0) : Ucell

    # Sort by current density for integration
    idx = sortperm(ifc)
    i_sorted = ifc[idx]
    U_sorted = U[idx]

    i_min = region.i_min
    i_max = region.i_max

    # Select points strictly inside the closed region (no boundary interpolation)
    mask = (i_sorted .>= i_min) .& (i_sorted .<= i_max)
    if !any(mask)
        return NaN
    end

    i_in = i_sorted[mask]
    U_in = U_sorted[mask]

    # Trapezoidal integration
    auc = 0.0
    for k in 1:(length(i_in) - 1)
        auc += 0.5 * (U_in[k] + U_in[k+1]) * (i_in[k+1] - i_in[k])
    end

    return max(auc, 0.0)
end


"""
    _linear_interp(x, x1, x2, y1, y2)

Linear interpolation of `y` at `x` from two points `(x1, y1)` and `(x2, y2)`.
"""
function _linear_interp(x::Float64, x1::Float64, x2::Float64, y1::Float64, y2::Float64)::Float64
    if x2 ≈ x1
        return 0.5 * (y1 + y2)
    end
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)
end


"""
    compute_regional_aucs(ifc, Ucell, regions)

Compute the AUC for each region in `regions`.

Returns a `Vector{Float64}` in the same order as `regions`.
"""
function compute_regional_aucs(ifc::Vector{Float64},
                               Ucell::Vector{Float64},
                               regions::Vector{PolarizationRegion})::Vector{Float64}
    return [region_auc(ifc, Ucell, r) for r in regions]
end
