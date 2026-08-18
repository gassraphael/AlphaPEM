# -*- coding: utf-8 -*-

"""
    add_confidence_intervals(df)

Add derived confidence-interval columns to a DataFrame of Sobol indices.

For each index column ending in `_conf`, adds:
- `CI_lower`  = index - 0.5 * conf
- `CI_upper`  = index + 0.5 * conf
- `CI_contains_0` = true if the interval contains 0
"""
function add_confidence_intervals(df::DataFrame)::DataFrame
    out = copy(df)
    for col in names(df)
        endswith(col, "_conf") || continue
        base = Symbol(chop(col; tail = 5))  # remove "_conf"
        hasproperty(out, base) || continue
        lower = out[!, base] .- 0.5 .* out[!, col]
        upper = out[!, base] .+ 0.5 .* out[!, col]
        out[!, Symbol(base, :_CI_lower)] = lower
        out[!, Symbol(base, :_CI_upper)] = upper
        out[!, Symbol(base, :_CI_contains_0)] = (lower .<= 0.0) .& (upper .>= 0.0)
    end
    return out
end


"""
    build_sobol_region_summary(result, region_name)

Build a single-region summary DataFrame with S1/ST/(S2) sums and significance checks.

Returns a one-row DataFrame with columns:
- `region`
- `sum_S1`, `sum_S1_positive`
- `sum_ST`
- `sum_S2`, `sum_S2_significant_positive` (if S2 is available)
- `negative_S1_warning`, `negative_S2_warning`
"""
function build_sobol_region_summary(result::SobolAnalysisResult,
                                     region_name::Symbol)::DataFrame
    df = result.sobol_indices[region_name]
    df = add_confidence_intervals(df)

    sum_S1 = sum(df.S1)
    sum_S1_positive = sum(df.S1[df.S1 .> 0])
    sum_ST = sum(df.ST)
    negative_S1_warning = any(df.S1 .< 0)

    s2_key = Symbol(:S2_, region_name)
    sum_S2 = NaN
    sum_S2_significant_positive = NaN
    negative_S2_warning = false

    if haskey(result.sobol_indices, s2_key)
        s2_df = result.sobol_indices[s2_key]
        s2_df = add_confidence_intervals(s2_df)
        sum_S2 = sum(s2_df.S2)
        sig_positive = (s2_df.S2 .> 0) .& .!s2_df.S2_CI_contains_0
        sum_S2_significant_positive = sum(s2_df.S2[sig_positive])
        negative_S2_warning = any(s2_df.S2 .< 0)
    end

    return DataFrame(
        region = [region_name],
        sum_S1 = [sum_S1],
        sum_S1_positive = [sum_S1_positive],
        sum_ST = [sum_ST],
        sum_S2 = [sum_S2],
        sum_S2_significant_positive = [sum_S2_significant_positive],
        negative_S1_warning = [negative_S1_warning],
        negative_S2_warning = [negative_S2_warning],
    )
end


"""
    build_sobol_summary_table(result; index_type=:ST)

Build a multi-region summary table with indices, confidence intervals and ranks.

Rows are parameters. Columns are grouped by region:
- `activation_value`, `activation_CI`, `activation_rank`
- `ohmic_value`, `ohmic_CI`, `ohmic_rank`
- `mass_transport_value`, `mass_transport_CI`, `mass_transport_rank`

The last rows contain `Sum` and `Avg_Conf`.
"""
function build_sobol_summary_table(result::SobolAnalysisResult;
                                   index_type::Symbol = :ST)::DataFrame
    params = result.param_names
    n_params = length(params)

    col_names = Symbol[]
    for r in result.regions
        push!(col_names, Symbol(r.name, :_value))
        push!(col_names, Symbol(r.name, :_CI))
        push!(col_names, Symbol(r.name, :_rank))
    end

    data = Dict{Symbol, Vector{Any}}()
    for c in col_names
        data[c] = Vector{Any}(undef, n_params)
    end

    for (r_idx, region) in enumerate(result.regions)
        df = add_confidence_intervals(result.sobol_indices[region.name])
        value_col = Symbol(region.name, :_value)
        ci_col = Symbol(region.name, :_CI)
        rank_col = Symbol(region.name, :_rank)

        # Sort to compute ranks
        perm = sortperm(df[!, index_type], rev = true)
        ranks = zeros(Int, nrow(df))
        for (rank, idx) in enumerate(perm)
            ranks[idx] = rank
        end

        for (i, p) in enumerate(params)
            row = findfirst(df.parameter .== p)
            if row === nothing
                data[value_col][i] = NaN
                data[ci_col][i] = NaN
                data[rank_col][i] = NaN
            else
                val = df[row, index_type]
                conf = df[row, Symbol(index_type, :_conf)]
                data[value_col][i] = val
                data[ci_col][i] = conf
                data[rank_col][i] = ranks[row]
            end
        end
    end

    summary = DataFrame(parameter = params)
    for c in col_names
        summary[!, c] = data[c]
    end

    # Add Sum and Avg_Conf rows
    sum_row = Dict{Symbol, Any}(:parameter => :Sum)
    avg_conf_row = Dict{Symbol, Any}(:parameter => :Avg_Conf)
    for r in result.regions
        value_col = Symbol(r.name, :_value)
        ci_col = Symbol(r.name, :_CI)
        rank_col = Symbol(r.name, :_rank)
        sum_row[value_col] = sum(skipmissing(replace(summary[!, value_col], NaN => missing)))
        sum_row[ci_col] = NaN
        sum_row[rank_col] = NaN
        avg_conf_row[value_col] = NaN
        avg_conf_row[ci_col] = mean(skipmissing(replace(summary[!, ci_col], NaN => missing)))
        avg_conf_row[rank_col] = NaN
    end
    push!(summary, sum_row)
    push!(summary, avg_conf_row)

    return summary
end


"""
    select_top_features(result; threshold=0.85, index_type=:ST)

Select the smallest set of parameters explaining at least `threshold` of the
cumulative sum of `index_type` per region.

Returns a `Dict{Symbol, Vector{Symbol}}` mapping region name to the selected
parameter names in order of decreasing importance.
"""
function select_top_features(result::SobolAnalysisResult;
                             threshold::Float64 = 0.85,
                             index_type::Symbol = :ST)::Dict{Symbol, Vector{Symbol}}
    selected = Dict{Symbol, Vector{Symbol}}()
    for region in result.regions
        df = result.sobol_indices[region.name]
        df_sorted = sort(df, index_type, rev = true)
        vals = df_sorted[!, index_type]
        total = sum(vals)
        if total <= 0
            selected[region.name] = Symbol[]
            continue
        end
        cumulative = cumsum(vals) ./ total
        n = findfirst(cumulative .>= threshold)
        n = n === nothing ? length(vals) : n
        selected[region.name] = df_sorted.parameter[1:n]
    end
    return selected
end
