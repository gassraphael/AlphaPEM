# -*- coding: utf-8 -*-

"""
    export_sobol_results(result, output_dir; curves_df=nothing)

Export Sobol sensitivity analysis results to CSV and YAML files.

If `curves_df` is provided, an additional CSV file listing the imputed samples and
their parameters is written.

Returns a dictionary of output file paths.
"""
function export_sobol_results(result::SobolAnalysisResult,
                               output_dir::String;
                               curves_df::Union{DataFrame, Nothing} = nothing)::Dict{Symbol, String}
    mkpath(output_dir)
    output_files = Dict{Symbol, String}()

    # Main indices per region
    for (region_name, df) in result.sobol_indices
        # Skip second-order composite keys in the main loop
        startswith(string(region_name), "S2_") && continue

        fname = "sobol_indices_$(region_name).csv"
        fpath = abspath(joinpath(output_dir, fname))
        CSV.write(fpath, df)
        output_files[Symbol(:sobol_indices_, region_name)] = fpath
    end

    # Second-order indices if present
    for (key, df) in result.sobol_indices
        if startswith(string(key), "S2_")
            region = replace(string(key), "S2_" => "")
            fname = "sobol_second_order_$(region).csv"
            fpath = abspath(joinpath(output_dir, fname))
            CSV.write(fpath, df)
            output_files[Symbol(:sobol_second_order_, Symbol(region))] = fpath
        end
    end

    # Imputed samples CSVs
    if curves_df !== nothing
        imputed_path = _write_imputed_samples_csv(result, curves_df, output_dir)
        if imputed_path !== nothing
            output_files[:sobol_imputed_samples_csv] = imputed_path
        end

        auc_imputed_path = _write_auc_imputed_samples_csv(result, output_dir)
        if auc_imputed_path !== nothing
            output_files[:sobol_auc_imputed_csv] = auc_imputed_path
        end
    end

    # Summary tables
    summary_path = abspath(joinpath(output_dir, "sobol_summary_table.csv"))
    CSV.write(summary_path, build_sobol_summary_table(result; index_type = :ST))
    output_files[:sobol_summary_table_csv] = summary_path

    region_summary_path = abspath(joinpath(output_dir, "sobol_region_summary.csv"))
    region_summary = vcat([build_sobol_region_summary(result, r.name) for r in result.regions]...)
    CSV.write(region_summary_path, region_summary)
    output_files[:sobol_region_summary_csv] = region_summary_path

    # Top features per region
    top_features = select_top_features(result; threshold = 0.85, index_type = :ST)
    top_features_path = abspath(joinpath(output_dir, "sobol_top_features.yaml"))
    _write_top_features_yaml(top_features, top_features_path)
    output_files[:sobol_top_features_yaml] = top_features_path

    # Indices with confidence intervals per region
    for (region_name, df) in result.sobol_indices
        startswith(string(region_name), "S2_") && continue
        fname = "sobol_indices_$(region_name)_with_ci.csv"
        fpath = abspath(joinpath(output_dir, fname))
        CSV.write(fpath, add_confidence_intervals(df))
        output_files[Symbol(:sobol_indices_, region_name, :_with_ci)] = fpath
    end

    # Summary YAML
    summary_path = abspath(joinpath(output_dir, "sobol_summary.yaml"))
    _write_sobol_summary_yaml(result, summary_path)
    output_files[:sobol_summary_yaml] = summary_path

    return output_files
end


"""
    _write_imputed_samples_csv(result, curves_df, output_dir)

Write a CSV file containing the parameters of all imputed samples.
Returns the file path or `nothing` if no sample was imputed.
"""
function _write_imputed_samples_csv(result::SobolAnalysisResult,
                                     curves_df::DataFrame,
                                     output_dir::String)::Union{String, Nothing}
    imputed_mask = curves_df.status .== :imputed
    count(imputed_mask) == 0 && return nothing

    cols = [:sample_id; result.param_names; :status]
    df_out = curves_df[imputed_mask, cols]
    fpath = abspath(joinpath(output_dir, "sobol_imputed_samples.csv"))
    CSV.write(fpath, df_out)
    return fpath
end


"""
    _write_auc_imputed_samples_csv(result, output_dir)

Write a CSV file listing samples whose regional AUCs were imputed.
Returns the file path or `nothing` if no AUC imputation occurred.
"""
function _write_auc_imputed_samples_csv(result::SobolAnalysisResult,
                                         output_dir::String)::Union{String, Nothing}
    report = result.imputation_report
    details = get(report, :auc_imputed_details, Vector{Tuple{Int, Symbol}}())
    isempty(details) && return nothing

    rows = []
    for (sample_id, region_name) in details
        push!(rows, (sample_id = sample_id, region = region_name))
    end
    df_out = DataFrame(rows)
    fpath = abspath(joinpath(output_dir, "sobol_auc_imputed_samples.csv"))
    CSV.write(fpath, df_out)
    return fpath
end


"""
    _write_top_features_yaml(top_features, filepath)

Write the top-feature selection per region to a YAML file.
"""
function _write_top_features_yaml(top_features::Dict{Symbol, Vector{Symbol}},
                                   filepath::String)
    open(filepath, "w") do io
        println(io, "top_features:")
        for (region, params) in top_features
            println(io, "  $(region):")
            for p in params
                println(io, "    - $(p)")
            end
        end
    end
    return nothing
end


"""
    _write_sobol_summary_yaml(result, filepath)

Write a YAML summary of the Sobol analysis run.
"""
function _write_sobol_summary_yaml(result::SobolAnalysisResult, filepath::String)
    cfg = result.config
    open(filepath, "w") do io
        println(io, "metadata:")
        println(io, "  fuel_cell_type: $(cfg.fuel_cell_type)")
        println(io, "  voltage_zone: $(cfg.voltage_zone)")
        println(io, "  N: $(cfg.N)")
        println(io, "  second_order: $(cfg.second_order)")
        println(io, "  seed: $(cfg.seed)")
        println(io, "  region_thresholds: $(cfg.region_thresholds)")
        println(io, "  include_operating_conditions: $(cfg.include_operating_conditions)")
        println(io, "  execution_time_s: $(round(result.execution_time; digits=2))")
        ts = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
        println(io, "  timestamp: $(ts)")

        println(io, "regions:")
        for r in result.regions
            println(io, "  $(r.name):")
            println(io, "    i_min: $(r.i_min)")
            println(io, "    i_max: $(r.i_max)")
        end

        println(io, "parameters:")
        for p in result.param_names
            println(io, "  - $(p)")
        end

        println(io, "imputation_report:")
        report = result.imputation_report
        println(io, "  n_total: $(report[:n_total])")
        println(io, "  n_failed_before_imputation: $(report[:n_failed_before_imputation])")
        println(io, "  n_curve_imputed: $(report[:n_imputed])")
        println(io, "  curve_imputation_rate: $(report[:imputation_rate])")
        println(io, "  sample_id_curve_imputed: $(report[:sample_id_imputed])")
        println(io, "  sample_id_curve_failed: $(report[:sample_id_failed])")
        println(io, "  n_auc_missing: $(get(report, :n_auc_missing, 0))")
        println(io, "  n_auc_imputed: $(get(report, :n_auc_imputed, 0))")
        println(io, "  sample_id_auc_imputed: $(get(report, :auc_imputed_sample_ids, Int[]))")

        println(io, "top_parameters_per_region:")
        for (region_name, df) in result.sobol_indices
            startswith(string(region_name), "S2_") && continue
            top = sort(df, :ST, rev=true)
            top_names = top.parameter[1:min(5, nrow(top))]
            top_str = join(string.(top_names), ", ")
            println(io, "  $(region_name): [$(top_str)]")
        end
    end
    return nothing
end
