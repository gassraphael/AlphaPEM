# -*- coding: utf-8 -*-

"""
    impute_missing_curves!(df, params; k=10)

Impute missing or failed polarization curves using K-Nearest Neighbors in parameter space.

For each row with missing `ifc`/`Ucell`, the function finds the `k` valid rows with the
smallest Euclidean distance in the normalized parameter space and averages their voltages
at each common current density.

Modifies `df` in place: missing curves are replaced by imputed ones and `status` becomes
`:imputed`.

Returns the number of imputed rows.
"""
function impute_missing_curves!(df::DataFrame,
                                params::Vector{InputParameter};
                                k::Int = 10)::Int
    param_names = [p.name for p in params]

    # Identify valid and missing rows
    valid_mask = [df.status[i] == :ok && df.ifc[i] !== nothing && df.Ucell[i] !== nothing for i in 1:nrow(df)]
    missing_mask = .!valid_mask

    n_valid = count(valid_mask)
    n_missing = count(missing_mask)

    if n_missing == 0
        return 0
    end
    if n_valid == 0
        @warn "No valid curves available for KNN imputation."
        return 0
    end

    k_eff = min(k, n_valid)

    valid_idx = findall(valid_mask)
    missing_idx = findall(missing_mask)

    # Normalize parameter matrix
    X_valid = _normalize_params(df[valid_idx, param_names], params)
    X_missing = _normalize_params(df[missing_idx, param_names], params)

    @info "Imputing $(n_missing) missing curves with KNN(k=$(k_eff))..."

    for (m_idx, row_idx) in enumerate(missing_idx)
        x = X_missing[m_idx, :]
        distances = [sqrt(sum((x .- X_valid[v, :]).^2)) for v in 1:n_valid]
        neighbor_order = sortperm(distances)
        neighbors = valid_idx[neighbor_order[1:k_eff]]

        # Aggregate neighbor curves on a common current grid
        ifc_imp, Ucell_imp = _aggregate_neighbor_curves(df, neighbors)

        if !isempty(ifc_imp) && !isempty(Ucell_imp)
            df.ifc[row_idx] = ifc_imp
            df.Ucell[row_idx] = Ucell_imp
            df.status[row_idx] = :imputed
        end
    end

    return n_missing
end


"""
    _normalize_params(df_cols, params)

Normalize parameter columns to [0, 1] for distance computation.
"""
function _normalize_params(df_cols::DataFrame,
                           params::Vector{InputParameter})::Matrix{Float64}
    n = nrow(df_cols)
    X = Matrix{Float64}(undef, n, length(params))
    for (j, p) in enumerate(params)
        span = p.max - p.min
        for i in 1:n
            v = df_cols[i, p.name]
            X[i, j] = span > 0 ? (Float64(v) - p.min) / span : 0.0
        end
    end
    return X
end


"""
    _aggregate_neighbor_curves(df, neighbors)

Average the polarization curves of `neighbors` on a common current-density grid.

Uses linear interpolation to align curves, then takes the mean voltage at each grid point.
"""
function _aggregate_neighbor_curves(df::DataFrame,
                                    neighbors::Vector{Int})
    # Collect all current densities from neighbors
    all_ifc = Float64[]
    for idx in neighbors
        append!(all_ifc, df.ifc[idx])
    end
    if isempty(all_ifc)
        return Float64[], Float64[]
    end

    # Common grid: unique sorted currents
    grid = sort(unique(all_ifc))
    if length(grid) < 2
        return Float64[], Float64[]
    end

    # Interpolate each neighbor onto the grid and accumulate
    U_sum = zeros(length(grid))
    U_count = zeros(Int, length(grid))

    for idx in neighbors
        ifc = df.ifc[idx]
        U = df.Ucell[idx]
        if length(ifc) != length(U) || isempty(ifc)
            continue
        end

        # Sort neighbor curve
        order = sortperm(ifc)
        ifc_s = ifc[order]
        U_s = U[order]

        for (g_idx, i_g) in enumerate(grid)
            if i_g < ifc_s[1] || i_g > ifc_s[end]
                continue
            end
            k = searchsortedlast(ifc_s, i_g)
            if k < 1
                k = 1
            end
            if k >= length(ifc_s)
                k = length(ifc_s) - 1
            end
            U_interp = _linear_interp(i_g, ifc_s[k], ifc_s[k+1], U_s[k], U_s[k+1])
            if isfinite(U_interp)
                U_sum[g_idx] += U_interp
                U_count[g_idx] += 1
            end
        end
    end

    mask = U_count .> 0
    if !any(mask)
        return Float64[], Float64[]
    end

    return grid[mask], U_sum[mask] ./ U_count[mask]
end
