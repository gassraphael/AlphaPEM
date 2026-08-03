# -*- coding: utf-8 -*-

"""
    export_sobol_results(result, output_dir)

Export Sobol sensitivity analysis results to CSV and YAML files.

Returns a dictionary of output file paths.
"""
function export_sobol_results(result::SobolAnalysisResult,
                               output_dir::String)::Dict{Symbol, String}
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

    # Summary YAML
    summary_path = abspath(joinpath(output_dir, "sobol_summary.yaml"))
    _write_sobol_summary_yaml(result, summary_path)
    output_files[:sobol_summary_yaml] = summary_path

    return output_files
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

        println(io, "top_parameters_per_region:")
        for (region_name, df) in result.sobol_indices
            startswith(string(region_name), "S2_") && continue
            top = sort(df, :ST, rev=true)
            top_names = top.parameter[1:min(5, nrow(top))]
            println(io, "  $(region_name): [$(join(top_names, ", "))]")
        end
    end
    return nothing
end
