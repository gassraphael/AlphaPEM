# -*- coding: utf-8 -*-

"""
    run_sobol_convergence_analysis(cfg, params, df_curves, regions; Ns=[128, 256, 512, 1024])

Recompute Sobol indices for increasing base sample sizes `Ns` using the same
polarization curves stored in `df_curves`.

For each `N` in `Ns`, a Saltelli plan of size `N` is generated, the corresponding
columns of `df_curves` are selected, and the Sobol indices are recomputed.

Returns a vector of `(N, sobol_indices_dict)` tuples.
"""
function run_sobol_convergence_analysis(cfg::SobolAnalysisConfig,
                                         params::Vector{InputParameter},
                                         df_curves::DataFrame,
                                         regions::Vector{PolarizationRegion};
                                         Ns::Vector{Int} = [128, 256, 512, 1024])::Vector{Tuple{Int, Dict{Symbol, DataFrame}}}
    results = Tuple{Int, Dict{Symbol, DataFrame}}[]
    max_available = nrow(df_curves)

    d = length(params)
    for N in sort(unique(Ns))
        n_total = cfg.second_order ? N * (2 * d + 2) : N * (d + 2)

        if n_total > max_available
            @warn "Convergence N=$N requires $(n_total) evaluations, but only $(max_available) are available; skipping."
            continue
        end

        # Generate Saltelli plan for this N
        cfg_N = SobolAnalysisConfig(
            fuel_cell_type = cfg.fuel_cell_type,
            year = cfg.year,
            voltage_zone = cfg.voltage_zone,
            N = N,
            second_order = cfg.second_order,
            seed = cfg.seed,
            region_thresholds = cfg.region_thresholds,
            include_operating_conditions = cfg.include_operating_conditions,
            parameter_bounds = cfg.parameter_bounds,
            operating_condition_constraints = cfg.operating_condition_constraints,
            excluded_operating_conditions = cfg.excluded_operating_conditions,
            polarization_params = cfg.polarization_params,
            numerical_params = cfg.numerical_params,
            nb_gc = cfg.nb_gc,
            parallel = cfg.parallel,
            max_run_time_s = cfg.max_run_time_s,
            knn_k = cfg.knn_k,
            output_dir = cfg.output_dir,
            save_curves = cfg.save_curves,
        )
        A_N, B_N = generate_sobol_design_matrices(cfg_N, params)

        # KNN imputation on the truncated dataframe
        df_N = df_curves[df_curves.sample_id .<= n_total, :]
        impute_missing_curves!(df_N, params; k = cfg.knn_k)

        # Build and impute regional AUC output matrix
        nboot_N = _sobol_nboot(N)
        all_points_N = _fuse_designs_bootstrap(A_N, B_N; second_order = cfg_N.second_order, nboot = nboot_N)
        Y_N = build_output_matrix(df_N, regions, all_points_N)
        auc_report_N = impute_missing_aucs!(Y_N, df_N, params, regions; k = cfg.knn_k)
        if any(isnan, Y_N)
            @warn "Convergence N=$(N) skipped: missing regional AUCs remain after imputation."
            continue
        end

        # Compute indices
        sobol_indices_N = compute_sobol_indices(cfg_N, params, A_N, B_N, Y_N, regions)
        push!(results, (N, sobol_indices_N))
    end

    return results
end


"""
    plot_sobol_index_convergence(convergence_results, param_name, region_name, index_type=:ST)

Plot the convergence of a single Sobol index for one parameter and one region
across the sample sizes stored in `convergence_results`.

`convergence_results` is the output of `run_sobol_convergence_analysis`.
Saves the figure to `output_dir/sobol_convergence_<param>_<region>_<index_type>.png`.
"""
function plot_sobol_index_convergence(convergence_results::Vector{Tuple{Int, Dict{Symbol, DataFrame}}},
                                       param_name::Symbol,
                                       region_name::Symbol,
                                       index_type::Symbol = :ST;
                                       output_dir::String = ".",
                                       figsize::Tuple{Int, Int} = (800, 500))
    Ns = Int[N for (N, _) in convergence_results]
    vals = Float64[]
    confs = Float64[]

    for (_, sobol_indices) in convergence_results
        df = sobol_indices[region_name]
        row = findfirst(df.parameter .== param_name)
        if row === nothing
            push!(vals, NaN)
            push!(confs, NaN)
        else
            push!(vals, df[row, index_type])
            push!(confs, df[row, Symbol(index_type, :_conf)])
        end
    end

    fig = Figure(size = figsize)
    ax = Axis(fig[1, 1],
              title = "Convergence of $(index_type) for $(param_name) ($(region_name))",
              xlabel = "N",
              ylabel = "$(index_type)")

    scatter!(ax, Ns, vals)
    errorbars!(ax, Ns, vals, confs)

    outpath = joinpath(output_dir, "sobol_convergence_$(param_name)_$(region_name)_$(index_type).png")
    save(outpath, fig)
    return outpath
end
