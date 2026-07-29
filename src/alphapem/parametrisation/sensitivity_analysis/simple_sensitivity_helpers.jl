# -*- coding: utf-8 -*-

"""
Helper functions for simple sensitivity analysis.

Contains utility functions for:
- Parameter extraction and modification
- Polarization point extraction and comparison
- Results processing and reporting
"""

using Dates
using Printf
using Statistics: mean, median
using AlphaPEM.Core.Models: _polarization_points_cali, _calculate_rmse

# ═══════════════════════════════════════════════════════════════════════════════
#  Polarization Point Extraction & RMSE Calculation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    extract_polarization_points_cali(simu, i_exp)

Extract polarization points sampled at the known calibration current densities `i_exp`
(A/m²). Every sensitivity-analysis run uses the same `PolarizationCalibrationParams`
current profile (independent of the varied PhysicalParams), so the requested current
densities are identical across the nominal and every modified run.
Returns tuple of (ifc_samples, Ucell_samples), in the same order for every run.
"""
function extract_polarization_points_cali(simu, i_exp::AbstractVector{<:Real})
    if isnothing(simu.outputs)
        return nothing
    end

    try
        ifc_samples, Ucell_samples = _polarization_points_cali(simu.outputs, simu.current_density, i_exp)
        return (ifc_samples, Ucell_samples)
    catch e
        @warn "Error extracting polarization points: $e"
        return nothing
    end
end

"""
    compute_rmse_from_nominal(simu_modified, i_exp, Ucell_nominal)

Compute RMSE between the modified and nominal polarization curves. Because both runs are
sampled at the exact same current densities `i_exp` and in the same order, no point-matching
is needed. The only case where the two curves can differ in length is when the modified run
stops early (safety stop at a lower current density before completing the sweep) — in that
case only the common leading points are compared.
"""
function compute_rmse_from_nominal(simu_modified, i_exp::Vector{Float64}, Ucell_nominal::Vector{Float64})
    pola_points = extract_polarization_points_cali(simu_modified, i_exp)
    if isnothing(pola_points)
        return Inf
    end

    _, Ucell_modified = pola_points
    if isempty(Ucell_modified) || any(isnan.(Ucell_modified)) || any(isinf.(Ucell_modified))
        return Inf
    end

    n = min(length(Ucell_modified), length(Ucell_nominal))
    if n == 0
        return Inf
    end

    try
        return _calculate_rmse(Ucell_modified[1:n], Ucell_nominal[1:n])
    catch e
        @warn "Error computing RMSE: $e"
        return Inf
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Results Processing & Output
# ═══════════════════════════════════════════════════════════════════════════════

"""
    compute_parameter_impacts(params_dict, results)

Compute impact of each parameter modification on RMSE vs baseline.
"""
function compute_parameter_impacts(params_dict::Dict, results::Dict)
    impacts = []

    for (param, nominal_value) in params_dict
        rmse_minus_key = "$(param)_minus20"
        rmse_plus_key = "$(param)_plus20"

        if haskey(results, rmse_minus_key) && haskey(results, rmse_plus_key)
            rmse_minus = results[rmse_minus_key]
            rmse_plus = results[rmse_plus_key]

            # Only include if both simulations succeeded (RMSE < 1000%)
            if !isinf(rmse_minus) && !isinf(rmse_plus) && rmse_minus < 1000 && rmse_plus < 1000
                avg_impact = (rmse_minus + rmse_plus) / 2

                push!(impacts, (
                    parameter = param,
                    nominal_value = nominal_value,
                    rmse_minus20 = rmse_minus,
                    rmse_plus20 = rmse_plus,
                    avg_impact_pct = avg_impact
                ))
            end
        end
    end

    sort!(impacts, by=x -> x.avg_impact_pct, rev=true)
    return impacts
end

"""
    generate_output_path(base_dir, base_filename)

Generate output path with date and optional counter to avoid overwrites.
"""
function generate_output_path(base_dir::String, base_filename::String)
    date_str = Dates.format(Dates.today(), dateformat"yyyy-mm-dd")
    base_name = split(base_filename, '.')[1]
    ext = "csv"

    path = joinpath(base_dir, "$(base_name)_$(date_str).$(ext)")
    if !isfile(path)
        return path
    end

    counter = 1
    while isfile(joinpath(base_dir, "$(base_name)_$(date_str)_$(counter).$(ext)"))
        counter += 1
    end

    return joinpath(base_dir, "$(base_name)_$(date_str)_$(counter).$(ext)")
end

"""
    write_sensitivity_csv(path, impacts)

Write sensitivity analysis results to CSV.
"""
function write_sensitivity_csv(path, impacts)
    open(path, "w") do io
        date_str = Dates.format(Dates.now(), dateformat"yyyy-mm-dd HH:MM:SS")
        println(io, "# Generated: $(date_str)")
        println(io, "# RMSE values are relative errors (%) between modified and nominal polarization curves")
        println(io, "Rank,Parameter,Nominal_Value,RMSE_-20%,RMSE_+20%,Avg_Impact_%")

        for (rank, impact) in enumerate(impacts)
            println(io, string(
                rank, ",",
                impact.parameter, ",",
                impact.nominal_value, ",",
                round(impact.rmse_minus20, digits=4), ",",
                round(impact.rmse_plus20, digits=4), ",",
                round(impact.avg_impact_pct, digits=4)
            ))
        end
    end
end

"""
    print_sensitivity_report(impacts, out_csv)

Print formatted sensitivity analysis report to console.
"""
function print_sensitivity_report(impacts, out_csv)
    println("\n" * "=" ^ 110)
    println("SENSITIVITY ANALYSIS REPORT")
    println("=" ^ 110)
    println("RMSE metric: Relative error (%) between modified and nominal polarization curves")
    println("Successful parameter variations: $(length(impacts))\n")

    if isempty(impacts)
        println("No successful parameter variations to report.")
        println("=" ^ 110)
        return
    end

    println(lpad("Rank", 5) * " | " * rpad("Parameter", 20) * " | " *
            rpad("Nominal Value", 15) * " | " * rpad("RMSE -20%", 14) * " | " *
            rpad("RMSE +20%", 14) * " | " * rpad("Avg Impact %", 14))
    println("-" ^ 110)

    # Parameters whose average impact exceeds this threshold are highlighted in red.
    highlight_threshold_pct = 0.05
    red = "\e[31m"
    reset = "\e[0m"

    num_display = min(30, length(impacts))
    for (rank, impact) in enumerate(impacts[1:num_display])
        nom_val = string(round(Float64(impact.nominal_value), sigdigits=4))
        rmse_m = string(round(impact.rmse_minus20, digits=4))
        rmse_p = string(round(impact.rmse_plus20, digits=4))
        avg_imp = string(round(impact.avg_impact_pct, digits=4))

        line = lpad(string(rank), 5) * " | " * rpad(impact.parameter, 20) * " | " *
               lpad(nom_val, 15) * " | " * lpad(rmse_m, 14) * " | " *
               lpad(rmse_p, 14) * " | " * lpad(avg_imp, 14)

        if impact.avg_impact_pct > highlight_threshold_pct
            println(red * line * reset)
        else
            println(line)
        end
    end

    if length(impacts) > num_display
        println("... ($(length(impacts) - num_display) more parameters)")
    end
    println("-" ^ 110)

    if !isempty(impacts)
        avg_impacts = [i.avg_impact_pct for i in impacts]
        println("\nStatistics (Avg Impact %):")
        @printf("  Mean:          %.4f %%\n", mean(avg_impacts))
        @printf("  Median:        %.4f %%\n", median(avg_impacts))
        @printf("  Max:           %.4f %% (parameter: %s)\n", maximum(avg_impacts), impacts[1].parameter)
        @printf("  Min:           %.4f %%\n", minimum(avg_impacts))
    end

    println("\nCSV written: $(out_csv)")
    println("=" ^ 110)
end
