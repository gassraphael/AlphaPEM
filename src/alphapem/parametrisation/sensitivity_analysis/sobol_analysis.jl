# -*- coding: utf-8 -*-

"""
    compute_sobol_indices(cfg, params, A, B, df, regions)

Compute Sobol sensitivity indices for each polarization region.

Returns a `Dict{Symbol, DataFrame}` mapping region name to a DataFrame with columns:
- `parameter`: input parameter name
- `S1`: first-order index
- `S1_conf`: confidence interval for S1
- `ST`: total-order index
- `ST_conf`: confidence interval for ST
- `S2`: second-order index (if requested)
- `S2_conf`: confidence interval for S2
"""
function compute_sobol_indices(cfg::SobolAnalysisConfig,
                               params::Vector{InputParameter},
                               A::Matrix{Float64},
                               B::Matrix{Float64},
                               df::DataFrame,
                               regions::Vector{PolarizationRegion})::Dict{Symbol, DataFrame}
    n_A = size(A, 2)
    n_B = size(B, 2)

    # Fuse design matrices in the exact order GlobalSensitivity expects
    all_points = GlobalSensitivity.fuse_designs(A, B; second_order = cfg.second_order)

    # Build output matrix Y aligned with all_points columns
    Y = _build_output_matrix(df, regions, all_points)

    param_names = [p.name for p in params]

    # Run GlobalSensitivity Sobol analysis using a dummy batch model
    order_vec = cfg.second_order ? [0, 1, 2] : [0, 1]
    method = Sobol(order = order_vec, nboot = 100, conf_level = 0.95)
    sobol_result = gsa((X) -> _model_for_gsa(X, Y), method, A, B; batch = true)

    # Assemble results per region
    results = Dict{Symbol, DataFrame}()
    for (r_idx, region) in enumerate(regions)
        df_region = DataFrame(parameter = param_names)
        df_region.S1 = sobol_result.S1[r_idx, :]
        df_region.S1_conf = sobol_result.S1_Conf_Int[r_idx, :]
        df_region.ST = sobol_result.ST[r_idx, :]
        df_region.ST_conf = sobol_result.ST_Conf_Int[r_idx, :]

        results[region.name] = df_region

        if cfg.second_order
            S2 = sobol_result.S2[r_idx, :, :]
            S2_conf = sobol_result.S2_Conf_Int[r_idx, :, :]
            results[Symbol(:S2_, region.name)] = _build_s2_dataframe(param_names, S2, S2_conf)
        end
    end

    return results
end


"""
    _build_output_matrix(df, regions, all_points)

Build the output matrix `Y` of size `(n_regions, n_samples)` aligned with `all_points`.
"""
function _build_output_matrix(df::DataFrame,
                              regions::Vector{PolarizationRegion},
                              all_points::Matrix{Float64})::Matrix{Float64}
    n_regions = length(regions)
    n_total = size(all_points, 2)
    Y = Matrix{Float64}(undef, n_regions, n_total)

    # Build a lookup from sample_id to row in df
    row_lookup = Dict(df.sample_id[i] => i for i in 1:nrow(df))

    for col in 1:n_total
        row_idx = get(row_lookup, col, nothing)
        if row_idx === nothing
            Y[:, col] .= NaN
            continue
        end
        ifc = df.ifc[row_idx]
        Ucell = df.Ucell[row_idx]
        if ifc === nothing || Ucell === nothing
            Y[:, col] .= NaN
        else
            for (r_idx, region) in enumerate(regions)
                Y[r_idx, col] = region_auc(ifc, Ucell, region)
            end
        end
    end

    if any(isnan, Y)
        error("Missing outputs remain after KNN imputation; cannot compute Sobol indices.")
    end

    return Y
end


"""
    _model_for_gsa(X, Y_ref)

Dummy batch model used to feed pre-computed outputs into `GlobalSensitivity.gsa`.

`X` is ignored; the function returns the corresponding columns of `Y_ref`.
"""
function _model_for_gsa(X::AbstractMatrix{<:Real}, Y_ref::Matrix{Float64})
    n_cols = size(X, 2)
    return Y_ref[:, 1:n_cols]
end


"""
    _build_s2_dataframe(param_names, S2, S2_conf)

Build a DataFrame of second-order indices.
"""
function _build_s2_dataframe(param_names::Vector{Symbol},
                             S2::AbstractMatrix{<:Real},
                             S2_conf::AbstractMatrix{<:Real})::DataFrame
    n = length(param_names)
    rows = []
    for i in 1:n
        for j in (i+1):n
            push!(rows, (
                param_i = param_names[i],
                param_j = param_names[j],
                S2 = S2[i, j],
                S2_conf = S2_conf[i, j]
            ))
        end
    end
    return DataFrame(rows)
end
