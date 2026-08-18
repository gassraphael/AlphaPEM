# -*- coding: utf-8 -*-

"""
    _sobol_nboot(N; max_boot=100, min_block_size=2)

Return the largest bootstrap block count `nboot <= max_boot` that divides `N`
and gives at least `min_block_size` columns per block.

`GlobalSensitivity.gsa` splits the base matrices into `nboot` contiguous blocks;
choosing a divisor of `N` guarantees that every column of the fused design matrix
is used and that confidence intervals are available when `nboot > 1`.
"""
function _sobol_nboot(N::Int; max_boot::Int = 100, min_block_size::Int = 2)::Int
    N <= 1 && return 1
    for nboot in min(N, max_boot):-1:2
        (N % nboot == 0) || continue
        (N ÷ nboot) >= min_block_size && return nboot
    end
    return 1
end


"""
    _fuse_designs_bootstrap(A, B; second_order, nboot)

Fuse Sobol design matrices in the same block-wise order used internally by
`GlobalSensitivity.gsa` when `nboot > 1`.

Each bootstrap block is fused with `GlobalSensitivity.fuse_designs`, and the
resulting blocks are concatenated horizontally. This ordering must be used both
for running the simulations and for building the output matrix `Y`.
"""
function _fuse_designs_bootstrap(A::Matrix{Float64},
                                 B::Matrix{Float64};
                                 second_order::Bool,
                                 nboot::Int)::Matrix{Float64}
    N = size(A, 2)
    @assert N == size(B, 2) "A and B must have the same number of columns"
    @assert N % nboot == 0 "nboot must divide the base sample size N"
    n = N ÷ nboot
    blocks = Matrix{Float64}[]
    for b in 1:nboot
        A_b = A[:, ((b - 1) * n + 1):(b * n)]
        B_b = B[:, ((b - 1) * n + 1):(b * n)]
        push!(blocks, GlobalSensitivity.fuse_designs(A_b, B_b; second_order = second_order))
    end
    return reduce(hcat, blocks)
end


"""
    compute_sobol_indices(cfg, params, A, B, Y, regions)

Compute Sobol sensitivity indices for each polarization region.

`Y` must be aligned with the block-fused design matrix used internally by
`GlobalSensitivity.gsa` (see `_fuse_designs_bootstrap`).

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
                               Y::Matrix{Float64},
                               regions::Vector{PolarizationRegion})::Dict{Symbol, DataFrame}
    param_names = [p.name for p in params]

    # Choose a bootstrap count that divides N so every column is used
    nboot = _sobol_nboot(cfg.N)
    order_vec = cfg.second_order ? [0, 1, 2] : [0, 1]
    method = Sobol(order = order_vec, nboot = nboot, conf_level = 0.95)

    expected_cols = size(Y, 2)
    sobol_result = gsa((X) -> _model_for_gsa(X, Y, expected_cols), method, A, B; batch = true)

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
            # GlobalSensitivity stores S2 as (param_i, param_j, output);
            # select the slice for the current output (region).
            S2 = sobol_result.S2[:, :, r_idx]
            S2_conf = sobol_result.S2_Conf_Int[:, :, r_idx]
            results[Symbol(:S2_, region.name)] = _build_s2_dataframe(param_names, S2, S2_conf)
        end
    end

    return results
end


"""
    build_output_matrix(df, regions, all_points)

Build the output matrix `Y` of size `(n_regions, n_samples)` aligned with `all_points`.

Rows correspond to regions and columns to samples. Entries are `NaN` when a region
has no point in the polarization curve; such entries can be imputed with
`impute_missing_aucs!` before computing Sobol indices.
"""
function build_output_matrix(df::DataFrame,
                              regions::Vector{PolarizationRegion},
                              all_points::Matrix{Float64})::Matrix{Float64}
    n_regions = length(regions)
    n_total = size(all_points, 2)
    Y = fill(NaN, n_regions, n_total)

    # Build a lookup from sample_id to row in df
    row_lookup = Dict(df.sample_id[i] => i for i in 1:nrow(df))

    for col in 1:n_total
        row_idx = get(row_lookup, col, nothing)
        if row_idx === nothing
            continue
        end
        ifc = df.ifc[row_idx]
        Ucell = df.Ucell[row_idx]
        if ifc === nothing || Ucell === nothing
            continue
        end
        for (r_idx, region) in enumerate(regions)
            Y[r_idx, col] = region_auc(ifc, Ucell, region)
        end
    end

    return Y
end


"""
    _model_for_gsa(X, Y_ref, expected_cols)

Dummy batch model used to feed pre-computed outputs into `GlobalSensitivity.gsa`.

`X` is ignored, but its number of columns is checked against `expected_cols`
(the number of columns in the fused design matrix) to detect any mismatch in
the batching strategy used by `GlobalSensitivity.gsa`.
"""
function _model_for_gsa(X::AbstractMatrix{<:Real},
                         Y_ref::Matrix{Float64},
                         expected_cols::Int)
    n_cols = size(X, 2)
    if n_cols != expected_cols
        error("Batch model received $(n_cols) columns but expected $(expected_cols). " *
              "The output matrix Y_ref is no longer aligned with the design matrix; " *
              "Sobol indices would be silently wrong.")
    end
    if n_cols > size(Y_ref, 2)
        error("Batch model requested $(n_cols) columns but Y_ref only has $(size(Y_ref, 2)).")
    end
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
