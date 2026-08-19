# -*- coding: utf-8 -*-

"""
Polarization point extraction and RMSE computation for simple sensitivity analysis.
"""

using AlphaPEM.Core.Models: _polarization_points_cali, _calculate_rmse

"""
    extract_polarization_points_cali(simu, i_exp)

Extract polarization points sampled at the known calibration current densities `i_exp`
(A/m²). Every sensitivity-analysis run uses the same `PolarizationCalibrationParams`
current profile (independent of the varied PhysicalParams), so the requested current
densities are identical across the nominal and every modified run.
Returns tuple of (ifc_samples, Ucell_samples), in the same order for every run.
"""
function extract_polarization_points_cali(simu, i_exp::AbstractVector{<:Real})
    if isnothing(simu.outputs)
        return nothing
    end

    try
        ifc_samples, Ucell_samples = _polarization_points_cali(simu.outputs, simu.current_density, i_exp)
        return (ifc_samples, Ucell_samples)
    catch e
        @warn "Error extracting polarization points: $e"
        return nothing
    end
end

"""
    _rmse_for_indices(U_modified, U_nominal, idx)

Compute RMSE (relative error, %) between modified and nominal voltages for the given indices.
"""
function _rmse_for_indices(U_modified::Vector{Float64}, U_nominal::Vector{Float64}, idx::Vector{Int})
    isempty(idx) && return NaN
    U_mod = U_modified[idx]
    U_nom = U_nominal[idx]
    if any(isnan.(U_mod)) || any(isinf.(U_mod)) || any(isnan.(U_nom)) || any(isinf.(U_nom))
        return NaN
    end
    return _calculate_rmse(U_mod, U_nom; digits = nothing)
end

"""
    compute_rmse_from_nominal(simu_modified, i_exp, Ucell_nominal)

Compute RMSE between the modified and nominal polarization curves. Because both runs are
sampled at the exact same current densities `i_exp` and in the same order, no point-matching
is needed. The only case where the two curves can differ in length is when the modified run
stops early (safety stop at a lower current density before completing the sweep) — in that
case only the common leading points are compared.

Returns `Inf` if the modified simulation failed; otherwise returns the global RMSE in %.
"""
function compute_rmse_from_nominal(simu_modified, i_exp::Vector{Float64}, Ucell_nominal::Vector{Float64})
    rmse_dict = compute_regional_rmse_from_nominal(simu_modified, i_exp, Ucell_nominal, PolarizationRegion[])
    return rmse_dict[:global]
end

"""
    compute_regional_rmse_from_nominal(simu_modified, i_exp, Ucell_nominal, regions)

Compute the global RMSE and the RMSE inside each `PolarizationRegion`. The current densities
`i_exp` are expected in A/m², consistent with `PolarizationCalibrationParams`.

Returns a `Dict{Symbol, Float64}` with keys `:global`, `:activation`, `:ohmic`, and
`:mass_transport`. If `regions` is empty, only `:global` is returned. Failed simulations
return `Inf` for all entries.
"""
function compute_regional_rmse_from_nominal(simu_modified, i_exp::Vector{Float64},
                                            Ucell_nominal::Vector{Float64},
                                            regions::Vector{PolarizationRegion})
    rmse_dict = Dict{Symbol, Float64}(:global => Inf)
    for r in regions
        rmse_dict[r.name] = Inf
    end

    pola_points = extract_polarization_points_cali(simu_modified, i_exp)
    if isnothing(pola_points)
        return rmse_dict
    end

    _, Ucell_modified = pola_points
    if isempty(Ucell_modified) || any(isnan.(Ucell_modified)) || any(isinf.(Ucell_modified))
        return rmse_dict
    end

    n = min(length(Ucell_modified), length(Ucell_nominal))
    if n == 0
        return rmse_dict
    end

    i_common = i_exp[1:n]
    U_mod_common = Ucell_modified[1:n]
    U_nom_common = Ucell_nominal[1:n]

    try
        rmse_dict[:global] = _calculate_rmse(U_mod_common, U_nom_common; digits = nothing)
    catch e
        @warn "Error computing global RMSE: $e"
        return rmse_dict
    end

    for r in regions
        mask = (i_common .>= r.i_min) .& (i_common .<= r.i_max)
        idx = findall(mask)
        if !isempty(idx)
            regional_rmse = _rmse_for_indices(U_mod_common, U_nom_common, idx)
            rmse_dict[r.name] = isfinite(regional_rmse) ? regional_rmse : Inf
        else
            rmse_dict[r.name] = Inf
        end
    end

    return rmse_dict
end
