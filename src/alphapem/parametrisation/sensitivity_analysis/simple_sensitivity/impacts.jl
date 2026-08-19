# -*- coding: utf-8 -*-

"""
Impact computation and parameter ranking for simple sensitivity analysis.
"""

"""
    compute_parameter_impacts(params_dict, results; variation_pct=5.0)

Compute impact of each parameter modification on RMSE vs baseline, both globally and in
each polarization region. `results` is expected to map `param_minusN`/`param_plusN` keys to
a `Dict{Symbol, Float64}` containing `:global`, `:activation`, `:ohmic`, and
`:mass_transport` RMSE values.

Returns a vector of named tuples sorted by the global average impact. Each tuple also
contains a `regional_impacts` dictionary mapping region name to average impact in that region.
"""
function compute_parameter_impacts(params_dict::Dict, results::Dict; variation_pct::Float64=5.0)
    impacts = []
    pct_int = round(Int, variation_pct)
    region_names = [:activation, :ohmic, :mass_transport]

    for (param, nominal_value) in params_dict
        rmse_minus_key = "$(param)_minus$(pct_int)"
        rmse_plus_key = "$(param)_plus$(pct_int)"

        if haskey(results, rmse_minus_key) && haskey(results, rmse_plus_key)
            rmse_minus = results[rmse_minus_key]
            rmse_plus = results[rmse_plus_key]

            # Only include if both simulations succeeded globally (RMSE < 1000%)
            global_minus = get(rmse_minus, :global, Inf)
            global_plus = get(rmse_plus, :global, Inf)
            if isinf(global_minus) || isinf(global_plus) || global_minus >= 1000 || global_plus >= 1000
                continue
            end

            regional_impacts = Dict{Symbol, Float64}()
            for r in region_names
                r_minus = get(rmse_minus, r, Inf)
                r_plus = get(rmse_plus, r, Inf)
                if isfinite(r_minus) && isfinite(r_plus) && r_minus < 1000 && r_plus < 1000
                    regional_impacts[r] = (r_minus + r_plus) / 2
                else
                    regional_impacts[r] = Inf
                end
            end

            avg_impact = (global_minus + global_plus) / 2

            push!(impacts, (
                parameter = param,
                nominal_value = nominal_value,
                rmse_minus = global_minus,
                rmse_plus = global_plus,
                avg_impact_pct = avg_impact,
                regional_rmse_minus = rmse_minus,
                regional_rmse_plus = rmse_plus,
                regional_impacts = regional_impacts
            ))
        end
    end

    sort!(impacts, by=x -> x.avg_impact_pct, rev=true)
    return impacts
end
