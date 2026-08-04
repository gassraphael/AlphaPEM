# -*- coding: utf-8 -*-

"""
    plot_sobol_indices(result; index_type=:ST, top_k=10, figsize=(1000, 600))

Create a grouped bar plot of Sobol indices across all regions.

Saves the figure to `output_dir/sobol_indices_<index_type>.png`.
"""
function plot_sobol_indices(result::SobolAnalysisResult;
                            index_type::Symbol = :ST,
                            top_k::Int = 10,
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
    plot_sobol_ranking(result; index_type=:ST, top_k=10, figsize=(1000, 600))

Create a horizontal lollipop plot ranking parameters by Sobol index per region.
"""
function plot_sobol_ranking(result::SobolAnalysisResult;
                            index_type::Symbol = :ST,
                            top_k::Int = 10,
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
