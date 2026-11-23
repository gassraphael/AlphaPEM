# -*- coding: utf-8 -*-

# ============================================================
# CURRENT-DENSITY SENSITIVITY SUMMARIES
# ============================================================

"""
    summarize_current_density_sensitivities(S, params, cfg) -> DataFrame

Summarize signed local sensitivities over the full current grid and configured current-density zones.

For each parameter and region, the function reports the median absolute sensitivity as the measure
of influence, while also reporting the signed median to retain a concise indication of direction.
Here, "signed" means that the sensitivity's `+` or `-` sign is retained: with a positive voltage
output, a positive sensitivity means voltage increases as the parameter increases, and a negative
sensitivity means voltage decreases. The sign comes from the output slope; it is not the sign of
the current density. The input matrix remains unchanged, so each signed sensitivity at every
current density stays available to the caller. Only finite sensitivity values are included; no
missing value is interpolated or imputed. A summary is published when at least one current point is
valid.

Parameters
----------
S : AbstractMatrix{<:Real}
    Signed sensitivities with shape `n_params × n_current_densities`.
params : AbstractVector{ParameterBound}
    Analysed parameters in the same order as the rows of `S`.
cfg : ExtendedLocalSensitivityConfig
    Configuration providing the current-density grid and named, non-overlapping zones.

Returns
-------
DataFrame
    One row per parameter and region, with columns `parameter`, `region`, `median_magnitude`,
    `median_signed`, `valid_points`, `total_points`, and `coverage`. `median_signed` is the median
    of the algebraically signed sensitivities and indicates the typical direction of the voltage
    effect in that region. A region with no valid points has `NaN` medians and zero coverage.

Raises
------
ArgumentError
    If matrix dimensions do not match the parameters and grid.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function summarize_current_density_sensitivities(
    S::AbstractMatrix{<:Real},
    params::AbstractVector{ParameterBound},
    cfg::ExtendedLocalSensitivityConfig,
)::DataFrame
    # === Validate the configuration and matrix dimensions ===
    validate_config(cfg)
    size(S, 1) == length(params) ||
        throw(ArgumentError("S must have one row per analysed parameter"))
    size(S, 2) == length(cfg.current_density_grid) ||
        throw(ArgumentError("S must have one column per configured current density"))

    # === Prepare the full-range and named-zone selections ===
    current_densities = cfg.current_density_grid
    region_masks = Pair{Symbol,BitVector}[
        :all_current_densities => trues(length(current_densities))
    ]
    for zone in cfg.current_density_zones
        # Include or exclude each boundary exactly as specified by the zone.
        above_lower = zone.include_lower ?
            current_densities .>= zone.lower : current_densities .> zone.lower
        below_upper = zone.include_upper ?
            current_densities .<= zone.upper : current_densities .< zone.upper
        push!(region_masks, zone.name => BitVector(above_lower .& below_upper))
    end

    # === Summarize finite sensitivities without altering the input matrix ===
    summary = DataFrame(
        parameter = Symbol[],
        region = Symbol[],
        median_magnitude = Float64[],
        median_signed = Float64[],
        valid_points = Int[],
        total_points = Int[],
        coverage = Float64[],
    )
    for (parameter_index, parameter) in enumerate(params)
        # Process the full grid first, then each configured zone in declaration order.
        for (region_name, region_mask) in region_masks
            selected = @view S[parameter_index, region_mask]
            finite_values = Float64[value for value in selected if isfinite(value)]
            valid_count = length(finite_values)
            total_count = count(region_mask)

            # Keep a row for empty zones and failed parameters so coverage remains explicit.
            magnitude_median = valid_count == 0 ? NaN : median(abs.(finite_values))
            signed_median = valid_count == 0 ? NaN : median(finite_values)
            coverage = total_count == 0 ? NaN : valid_count / total_count
            push!(summary, (
                parameter.name,
                region_name,
                magnitude_median,
                signed_median,
                valid_count,
                total_count,
                coverage,
            ))
        end
    end

    return summary
end
