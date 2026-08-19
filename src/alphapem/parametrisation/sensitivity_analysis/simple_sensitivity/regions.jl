# -*- coding: utf-8 -*-

"""
Polarization regions for simple sensitivity analysis.

Defines the `PolarizationRegion` container and the `build_regions` constructor used to
split the polarization curve into activation, ohmic and mass-transport zones.
"""

"""
    PolarizationRegion

Named container for a characteristic region of the polarization curve.
`i_min` and `i_max` are stored in A/m².
"""
struct PolarizationRegion
    name::Symbol
    i_min::Float64
    i_max::Float64
end

"""
    build_regions(thresholds)

Build the three characteristic polarization regions from current-density thresholds.

`thresholds` is a tuple `(t1, t2)` in A/cm². The regions are closed on both sides
and converted internally to A/m²:
- activation:      [0, t1]
- ohmic:           [t1, t2]
- mass_transport:  [t2, ∞)

Returns a vector of `PolarizationRegion`.
"""
function build_regions(thresholds::Tuple{Float64, Float64})::Vector{PolarizationRegion}
    t1, t2 = thresholds
    0.0 < t1 < t2 || throw(ArgumentError("Region thresholds must satisfy 0 < t1 < t2 (got $thresholds A/cm²)"))
    t1_m2 = t1 * 1e4
    t2_m2 = t2 * 1e4
    return [
        PolarizationRegion(:activation,     0.0, t1_m2),
        PolarizationRegion(:ohmic,          t1_m2, t2_m2),
        PolarizationRegion(:mass_transport, t2_m2, Inf),
    ]
end
