# -*- coding: utf-8 -*-

"""
    plot_sobol_indices(result; index_type=:ST, top_k=13, figsize=(1000, 600))

Create a grouped bar plot of Sobol indices across all regions.

Saves the figure to `output_dir/sobol_indices_<index_type>.png`.
"""
function plot_sobol_indices(result::SobolAnalysisResult;
                            index_type::Symbol = :ST,
                            top_k::Int = 13,
                            figsize::Tuple{Int, Int} = (1000, 600))
    regions = result.regions
    n_regions = length(regions)

    # Collect top-k parameters per region
    top_params = Set{Symbol}()
    region_values = Dict{Symbol, Dict{Symbol, Float64}}()
    for r in regions
        df = result.sobol_indices[r.name]
        df_sorted = sort(df, index_type, rev = true)
        for p in df_sorted.parameter[1:min(top_k, nrow(df_sorted))]
            push!(top_params, p)
        end
        region_values[r.name] = Dict(df.parameter .=> df[!, index_type])
    end

    params = sort(collect(top_params))
    n_params = length(params)

    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              title = "Sobol $(index_type) indices by region",
              xlabel = "Parameter",
              ylabel = "$(index_type)",
              xticks = (1:n_params, string.(params)),
              xticklabelrotation = π / 4,
              xticklabelalign = (:right, :center))

    bar_width = 0.25
    offsets = range(-bar_width * (n_regions - 1) / 2, bar_width * (n_regions - 1) / 2; length = n_regions)
    colors = Makie.wong_colors()

    for (r_idx, r) in enumerate(regions)
        vals = [get(region_values[r.name], p, 0.0) for p in params]
        xs = 1:n_params
        barplot!(ax, xs .+ offsets[r_idx], vals;
                 width = bar_width,
                 color = colors[r_idx],
                 label = string(r.name))
    end

    axislegend(ax; position = :rt)

    outpath = joinpath(result.config.output_dir, "sobol_indices_$(index_type).png")
    save(outpath, fig)
    return outpath
end


"""
    plot_sobol_ranking(result; index_type=:ST, top_k=13, figsize=(1000, 600))

Create a horizontal lollipop plot ranking parameters by Sobol index per region.
"""
function plot_sobol_ranking(result::SobolAnalysisResult;
                            index_type::Symbol = :ST,
                            top_k::Int = 13,
                            figsize::Tuple{Int, Int} = (1000, 600))
    regions = result.regions
    n_regions = length(regions)

    fig = Figure(size = figsize)

    for (r_idx, r) in enumerate(regions)
        df = result.sobol_indices[r.name]
        df_sorted = sort(df, index_type, rev = true)
        df_top = df_sorted[1:min(top_k, nrow(df_sorted)), :]

        ax = Axis(fig[r_idx, 1],
                  title = string(r.name),
                  ylabel = "Parameter",
                  xlabel = "$(index_type)",
                  yticks = (1:nrow(df_top), string.(df_top.parameter)),
                  yreversed = true,
                  yticklabelalign = (:right, :center))

        for i in 1:nrow(df_top)
            val = df_top[i, index_type]
            lines!(ax, [0.0, val], [i, i])
            scatter!(ax, [val], [i])
        end
    end

    outpath = joinpath(result.config.output_dir, "sobol_ranking_$(index_type).png")
    save(outpath, fig)
    return outpath
end


"""
    plot_top_k_rankings_across_regions(result; index_type=:ST, top_k=13, figsize=(1000, 600))

Plot the ranking of the top-k parameters across all regions as a heatmap.

Rows are parameters, columns are regions, and cell color encodes the rank.
Saves the figure to `output_dir/sobol_top_k_rankings_<index_type>.png`.
"""
function plot_top_k_rankings_across_regions(result::SobolAnalysisResult;
                                             index_type::Symbol = :ST,
                                             top_k::Int = 13,
                                             figsize::Tuple{Int, Int} = (1000, 600))
    regions = result.regions
    n_regions = length(regions)

    # Collect top-k parameters per region
    top_per_region = Dict{Symbol, Vector{Symbol}}()
    for r in regions
        df = result.sobol_indices[r.name]
        df_sorted = sort(df, index_type, rev = true)
        top_per_region[r.name] = df_sorted.parameter[1:min(top_k, nrow(df_sorted))]
    end

    all_top = sort(collect(union(values(top_per_region)...)))
    n_params = length(all_top)

    # Build rank matrix (NaN if parameter not in top-k of region)
    rank_matrix = Matrix{Float64}(undef, n_params, n_regions)
    for (j, r) in enumerate(regions)
        df = result.sobol_indices[r.name]
        df_sorted = sort(df, index_type, rev = true)
        rank_lookup = Dict(df_sorted.parameter[i] => i for i in 1:nrow(df_sorted))
        for (i, p) in enumerate(all_top)
            rank_matrix[i, j] = get(rank_lookup, p, NaN)
        end
    end

    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              title = "Top-$(top_k) parameter rankings across regions",
              xticks = (1:n_regions, string.([r.name for r in regions])),
              yticks = (1:n_params, string.(all_top)),
              yreversed = true)

    hm = heatmap!(ax, rank_matrix; colorrange = (1, top_k), colormap = :viridis)
    Colorbar(fig[1, 2], hm; label = "Rank")

    outpath = joinpath(result.config.output_dir, "sobol_top_k_rankings_$(index_type).png")
    save(outpath, fig)
    return outpath
end


"""
    plot_sobol_heatmap(result; region_name=:activation, figsize=(900, 700))

Plot a heatmap of second-order Sobol indices S2 for one region.

Saves the figure to `output_dir/sobol_second_order_heatmap_<region_name>.png`.
"""
function plot_sobol_heatmap(result::SobolAnalysisResult;
                             region_name::Symbol = :activation,
                             figsize::Tuple{Int, Int} = (900, 700))
    s2_key = Symbol(:S2_, region_name)
    haskey(result.sobol_indices, s2_key) || error("No second-order indices for region $(region_name)")

    s2_df = result.sobol_indices[s2_key]
    params = unique(vcat(s2_df.param_i, s2_df.param_j))
    n = length(params)
    param_idx = Dict(params[i] => i for i in 1:n)

    M = Matrix{Float64}(undef, n, n)
    fill!(M, NaN)
    for row in eachrow(s2_df)
        i = param_idx[row.param_i]
        j = param_idx[row.param_j]
        M[i, j] = row.S2
        M[j, i] = row.S2
    end

    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              title = "Second-order Sobol indices S2 ($(region_name))",
              xticks = (1:n, string.(params)),
              yticks = (1:n, string.(params)),
              xticklabelrotation = π / 4)

    hm = heatmap!(ax, M; colorrange = (0, maximum(filter(isfinite, M))), colormap = :viridis)
    Colorbar(fig[1, 2], hm; label = "S2")

    outpath = joinpath(result.config.output_dir, "sobol_second_order_heatmap_$(region_name).png")
    save(outpath, fig)
    return outpath
end
