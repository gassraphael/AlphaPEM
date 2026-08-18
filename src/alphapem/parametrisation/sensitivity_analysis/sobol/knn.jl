# -*- coding: utf-8 -*-

"""
    impute_missing_curves!(df, params; k=10)

Impute missing or failed polarization curves using K-Nearest Neighbors in parameter space.

Only rows whose curve is classified as `:valid` by `ValidityCriteria` are used as
neighbours. Parameters are normalized with an empirical z-score (mean and standard
deviation computed on the combined pool). Voltages are averaged point-to-point,
assuming all successful simulations share the same current-density grid.

Modifies `df` in place: missing curves are replaced by imputed ones and `status`
becomes `:imputed`.

Returns the number of curves actually imputed with success.
"""
function impute_missing_curves!(df::DataFrame,
                                params::Vector{InputParameter};
                                k::Int = 10)::Int
    # Identify rows to impute
    missing_mask = [df.ifc[i] === nothing || df.Ucell[i] === nothing for i in 1:nrow(df)]
    n_missing = count(missing_mask)

    if n_missing == 0
        return 0
    end

    # Build the pool of valid neighbour curves
    valid_idx = _find_valid_neighbor_rows(df)
    n_valid = length(valid_idx)

    if n_valid == 0
        @warn "No valid curves available for KNN imputation."
        return 0
    end

    k_eff = min(k, n_valid)
    missing_idx = findall(missing_mask)

    # Empirical z-score normalization on the combined population
    X_valid, X_missing = _normalize_params_zscore(df, valid_idx, missing_idx, params)

    @info "Imputing $(n_missing) missing curves with KNN(k=$(k_eff))..."

    n_imputed = 0
    prog = Progress(n_missing; desc = "KNN imputation: ", barlen = 40, color = :magenta)
    for (m_idx, row_idx) in enumerate(missing_idx)
        x = X_missing[m_idx, :]
        distances = [sqrt(sum((x .- X_valid[v, :]).^2)) for v in 1:n_valid]
        neighbor_order = sortperm(distances)
        neighbors = valid_idx[neighbor_order[1:k_eff]]

        # Aggregate neighbor curves point-to-point
        ifc_imp, Ucell_imp = _aggregate_neighbor_curves(df, neighbors)

        if !isempty(ifc_imp) && !isempty(Ucell_imp)
            df.ifc[row_idx] = ifc_imp
            df.Ucell[row_idx] = Ucell_imp
            df.status[row_idx] = :imputed
            n_imputed += 1
        end
        next!(prog)
    end
    finish!(prog)

    if n_imputed < n_missing
        failed_rows = setdiff(missing_idx, findall(df.status .== :imputed))
        @warn "KNN imputation succeeded for $(n_imputed)/$(n_missing) rows; rows $(failed_rows) could not be imputed."
    end

    return n_imputed
end


"""
    _find_valid_neighbor_rows(df)

Return the row indices of curves that can be used as KNN neighbours.

A row is eligible when its simulation succeeded (`status == :ok`) and its
polarization curve is classified as `:valid` by `ValidityCriteria`.
"""
function _find_valid_neighbor_rows(df::DataFrame)::Vector{Int}
    valid_idx = Int[]
    for i in 1:nrow(df)
        df.status[i] == :ok || continue
        ifc = df.ifc[i]
        Ucell = df.Ucell[i]
        ifc === nothing && continue
        Ucell === nothing && continue
        isempty(ifc) && continue
        isempty(Ucell) && continue
        result = ValidityCriteria.classify_polarization_curve(Ucell, ifc)
        result.classification == :valid && push!(valid_idx, i)
    end
    return valid_idx
end


"""
    _normalize_params_zscore(df, valid_idx, missing_idx, params)

Normalize parameter columns with an empirical z-score (zero mean, unit variance).

Mean and standard deviation are computed on the combined valid + missing population,
mimicking the behaviour of `sklearn.preprocessing.StandardScaler`.
"""
function _normalize_params_zscore(df::DataFrame,
                                  valid_idx::Vector{Int},
                                  missing_idx::Vector{Int},
                                  params::Vector{InputParameter})::Tuple{Matrix{Float64}, Matrix{Float64}}
    all_idx = vcat(valid_idx, missing_idx)
    n_params = length(params)

    X_all = Matrix{Float64}(undef, length(all_idx), n_params)
    for (j, p) in enumerate(params)
        col = Float64[df[i, p.name] for i in all_idx]
        μ = mean(col)
        σ = std(col; corrected = false)
        for (i, v) in enumerate(col)
            X_all[i, j] = σ > 0 ? (v - μ) / σ : 0.0
        end
    end

    n_valid = length(valid_idx)
    X_valid = X_all[1:n_valid, :]
    X_missing = X_all[(n_valid + 1):end, :]
    return X_valid, X_missing
end


"""
    _aggregate_neighbor_curves(df, neighbors)

Average the polarization curves of `neighbors` point-to-point.

All neighbour curves are assumed to share the same current-density grid (the same
`PolarizationParams` profile is used for every simulation). The returned curve has
the same grid as the neighbours.
"""
function _aggregate_neighbor_curves(df::DataFrame,
                                    neighbors::Vector{Int})
    isempty(neighbors) && return Float64[], Float64[]

    # Use the first neighbour as the reference grid
    ref_ifc = df.ifc[neighbors[1]]
    ref_n = length(ref_ifc)
    if ref_n == 0
        return Float64[], Float64[]
    end

    # Check that all neighbours share the same grid
    for idx in neighbors
        ifc = df.ifc[idx]
        if length(ifc) != ref_n || ifc != ref_ifc
            @warn "Neighbour curves do not share the same current-density grid; cannot aggregate point-to-point."
            return Float64[], Float64[]
        end
    end

    U_sum = zeros(ref_n)
    for idx in neighbors
        U = df.Ucell[idx]
        if length(U) != ref_n
            @warn "Inconsistent curve length for neighbour row $idx; skipping in aggregation."
            continue
        end
        U_sum .+= U
    end

    return ref_ifc, U_sum ./ length(neighbors)
end


"""
    impute_missing_aucs!(Y, df, params, regions; k=10)

Impute missing regional AUC entries in `Y` using KNN in parameter space.

Only samples that have finite AUCs for all regions are used as neighbors. For each
missing entry, the imputed value is the mean of the corresponding regional AUCs of
the `k` closest valid neighbors.

Returns a report dictionary with:
- `:n_auc_missing`: number of samples with at least one missing regional AUC
- `:n_auc_imputed`: number of (sample, region) entries actually imputed
- `:auc_imputed_sample_ids`: sample IDs with at least one imputed AUC
- `:auc_imputed_details`: vector of `(sample_id, region_name)` tuples
"""
function impute_missing_aucs!(Y::Matrix{Float64},
                               df::DataFrame,
                               params::Vector{InputParameter},
                               regions::Vector{PolarizationRegion};
                               k::Int = 10)::Dict{Symbol, Any}
    n_regions, n_total = size(Y)

    # Identify samples with any missing regional AUC
    missing_cols = Int[]
    missing_regions_per_col = Vector{Symbol}[]
    for col in 1:n_total
        missing_regs = Symbol[]
        for (r_idx, region) in enumerate(regions)
            if isnan(Y[r_idx, col])
                push!(missing_regs, region.name)
            end
        end
        if !isempty(missing_regs)
            push!(missing_cols, col)
            push!(missing_regions_per_col, missing_regs)
        end
    end

    n_missing = length(missing_cols)
    if n_missing == 0
        return Dict{Symbol, Any}(
            :n_auc_missing => 0,
            :n_auc_imputed => 0,
            :auc_imputed_sample_ids => Int[],
            :auc_imputed_details => Vector{Tuple{Int, Symbol}}(),
        )
    end

    # Valid neighbors have finite AUCs for every region
    valid_cols = Int[]
    for col in 1:n_total
        if all(isfinite, Y[:, col])
            push!(valid_cols, col)
        end
    end

    if isempty(valid_cols)
        error("No sample has finite AUCs for all regions; cannot impute missing regional AUCs.")
    end

    k_eff = min(k, length(valid_cols))
    X_valid, X_missing = _normalize_params_zscore(df, valid_cols, missing_cols, params)

    imputed_details = Tuple{Int, Symbol}[]

    @info "Imputing missing regional AUCs for $(n_missing) sample(s) with KNN(k=$(k_eff))..."
    for (m_idx, col) in enumerate(missing_cols)
        x = X_missing[m_idx, :]
        distances = [sqrt(sum((x .- X_valid[v, :]).^2)) for v in 1:length(valid_cols)]
        neighbor_order = sortperm(distances)
        neighbors = valid_cols[neighbor_order[1:k_eff]]

        for (r_idx, region) in enumerate(regions)
            if isnan(Y[r_idx, col])
                neighbor_vals = [Y[r_idx, nb] for nb in neighbors if isfinite(Y[r_idx, nb])]
                if !isempty(neighbor_vals)
                    Y[r_idx, col] = mean(neighbor_vals)
                    push!(imputed_details, (col, region.name))
                end
            end
        end
    end

    imputed_sample_ids = sort!(collect(Set([sid for (sid, _) in imputed_details])))

    return Dict{Symbol, Any}(
        :n_auc_missing => n_missing,
        :n_auc_imputed => length(imputed_details),
        :auc_imputed_sample_ids => imputed_sample_ids,
        :auc_imputed_details => imputed_details,
    )
end
