# -*- coding: utf-8 -*-

"""
Output generation and console reporting for simple sensitivity analysis.
"""

using Dates
using Printf

"""
    generate_output_path(base_dir, base_filename)

Generate output path with date and optional counter to avoid overwrites.
The extension of `base_filename` is preserved if it is `.txt`; otherwise `.txt` is used.
"""
function generate_output_path(base_dir::String, base_filename::String)
    date_str = Dates.format(Dates.today(), dateformat"yyyy-mm-dd")
    base_name = splitext(base_filename)[1]
    ext = lowercase(splitext(base_filename)[2]) == ".txt" ? "txt" : "txt"

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
    write_sensitivity_txt(path, impacts, regions; variation_pct=5.0)

Write sensitivity analysis results to a formatted text file. The document contains four
separate, clearly delimited tables: one global ranking and one ranking per polarization
region (activation, ohmic, mass transport).
"""
function write_sensitivity_txt(path, impacts, regions::Vector{PolarizationRegion}; variation_pct::Float64=5.0)
    pct_int = round(Int, variation_pct)
    date_str = Dates.format(Dates.now(), dateformat"yyyy-mm-dd HH:MM:SS")
    region_order = [:activation, :ohmic, :mass_transport]
    region_lookup = Dict(r.name => r for r in regions)

    open(path, "w") do io
        println(io, "=" ^ 80)
        println(io, "  AlphaPEM — Simple Sensitivity Analysis")
        println(io, "=" ^ 80)
        println(io, "Generated:        $(date_str)")
        println(io, "Variation:        ±$(variation_pct)%")
        println(io, "RMSE metric:      Relative error (%) between modified and nominal curves")
        println(io, "Successful params: $(length(impacts))")
        println(io, "=" ^ 80)

        # Global table
        println(io, "\n" * "─" ^ 80)
        println(io, "  TABLE 1 — GLOBAL IMPACT RANKING")
        println(io, "─" ^ 80)
        _write_impact_table(io, impacts, :global; variation_pct = variation_pct)

        # Regional tables
        table_num = 2
        for rname in region_order
            if !haskey(region_lookup, rname)
                continue
            end
            r = region_lookup[rname]
            println(io, "\n" * "─" ^ 80)
            println(io, "  TABLE $(table_num) — REGION: $(uppercase(String(rname))) ($(r.i_min/1e4) – $(isfinite(r.i_max) ? r.i_max/1e4 : "∞") A/cm²)")
            println(io, "─" ^ 80)
            _write_impact_table(io, impacts, rname; variation_pct = variation_pct)
            table_num += 1
        end

        println(io, "\n" * "=" ^ 80)
        println(io, "End of report")
        println(io, "=" ^ 80)
    end
end

"""
    _write_impact_table(io, impacts, region; variation_pct=5.0)

Write a single formatted ranking table to `io` for the requested `region` (`:global` or a
region name such as `:activation`). Parameters are sorted by decreasing impact in that region.
"""
function _write_impact_table(io::IO, impacts, region::Symbol; variation_pct::Float64=5.0)
    pct_int = round(Int, variation_pct)
    col_rank = 6
    col_param = 20
    col_nom = 16
    col_minus = 14
    col_plus = 14
    col_avg = 14

    header = (
        lpad("Rank", col_rank) * " │ " *
        rpad("Parameter", col_param) * " │ " *
        lpad("Nominal", col_nom) * " │ " *
        lpad("RMSE -$(pct_int)%", col_minus) * " │ " *
        lpad("RMSE +$(pct_int)%", col_plus) * " │ " *
        lpad("Avg Impact %", col_avg)
    )
    sep = "─" ^ col_rank * "─┼─" * "─" ^ col_param * "─┼─" * "─" ^ col_nom *
          "─┼─" * "─" ^ col_minus * "─┼─" * "─" ^ col_plus * "─┼─" * "─" ^ col_avg

    println(io, header)
    println(io, sep)

    if isempty(impacts)
        println(io, "No successful parameter variations to report.")
        return
    end

    # Sort by the requested region
    sorted = if region === :global
        impacts
    else
        sort(impacts, by = x -> get(x.regional_impacts, region, -Inf), rev = true)
    end

    for (rank, impact) in enumerate(sorted)
        if region === :global
            rmse_m = impact.rmse_minus
            rmse_p = impact.rmse_plus
            avg_imp = impact.avg_impact_pct
        else
            rmse_m = get(impact.regional_rmse_minus, region, NaN)
            rmse_p = get(impact.regional_rmse_plus, region, NaN)
            avg_imp = get(impact.regional_impacts, region, NaN)
        end

        line = (
            lpad(string(rank), col_rank) * " │ " *
            rpad(impact.parameter, col_param) * " │ " *
            lpad(_fmt_value(impact.nominal_value), col_nom) * " │ " *
            lpad(_fmt_float(rmse_m), col_minus) * " │ " *
            lpad(_fmt_float(rmse_p), col_plus) * " │ " *
            lpad(_fmt_float(avg_imp), col_avg)
        )
        println(io, line)
    end
end

"""
    _fmt_value(x)

Format a nominal parameter value for tabular output.
"""
_fmt_value(x) = string(x)
function _fmt_value(x::AbstractFloat)
    isfinite(x) || return "—"
    abs(x) == 0 && return "0"
    abs(x) >= 1e4 || abs(x) <= 1e-3 ? @sprintf("%.4e", x) : @sprintf("%.4f", x)
end

"""
    _fmt_float(x)

Format a floating-point RMSE for tabular output.
"""
_fmt_float(x) = isfinite(x) ? @sprintf("%.4f", x) : "—"

"""
    print_sensitivity_report(impacts, out_txt, regions; variation_pct=5.0)

Print formatted sensitivity analysis report to console, including global and regional rankings.
"""
function print_sensitivity_report(impacts, out_txt, regions::Vector{PolarizationRegion}; variation_pct::Float64=5.0)
    pct_int = round(Int, variation_pct)
    region_order = [:activation, :ohmic, :mass_transport]
    region_lookup = Dict(r.name => r for r in regions)

    println("\n" * "=" ^ 110)
    println("SENSITIVITY ANALYSIS REPORT")
    println("=" ^ 110)
    println("RMSE metric: Relative error (%) between modified and nominal polarization curves")
    println("Variation amplitude: ±$(variation_pct)%")
    println("Successful parameter variations: $(length(impacts))\n")

    if isempty(impacts)
        println("No successful parameter variations to report.")
        println("=" ^ 110)
        return
    end

    # Global table
    println("GLOBAL RANKING")
    println("-" ^ 110)
    _print_impact_table(impacts, :global; variation_pct = variation_pct)

    # Regional tables
    for rname in region_order
        if !haskey(region_lookup, rname)
            continue
        end
        r = region_lookup[rname]
        println("\nREGION: $(uppercase(String(rname))) ($(r.i_min/1e4) – $(isfinite(r.i_max) ? r.i_max/1e4 : "∞") A/cm²)")
        println("-" ^ 110)
        sorted = sort(impacts, by = x -> get(x.regional_impacts, rname, -Inf), rev = true)
        _print_impact_table(sorted, rname; variation_pct = variation_pct)
    end

    println("\nTXT written: $(out_txt)")
    println("=" ^ 110)
end

"""
    _print_impact_table(impacts, region; variation_pct=5.0)

Print a single ranking table to the console for the requested region.
"""
function _print_impact_table(impacts, region::Symbol; variation_pct::Float64=5.0)
    pct_int = round(Int, variation_pct)
    println(lpad("Rank", 5) * " | " * rpad("Parameter", 20) * " | " *
            rpad("Nominal Value", 15) * " | " * rpad("RMSE -$(pct_int)%", 14) * " | " *
            rpad("RMSE +$(pct_int)%", 14) * " | " * rpad("Avg Impact %", 14))
    println("-" ^ 110)

    highlight_threshold_pct = 1.0
    red = "\e[31m"
    reset = "\e[0m"

    num_display = min(30, length(impacts))
    for (rank, impact) in enumerate(impacts[1:num_display])
        if region === :global
            rmse_m = impact.rmse_minus
            rmse_p = impact.rmse_plus
            avg_imp = impact.avg_impact_pct
        else
            rmse_m = get(impact.regional_rmse_minus, region, NaN)
            rmse_p = get(impact.regional_rmse_plus, region, NaN)
            avg_imp = get(impact.regional_impacts, region, NaN)
        end

        nom_val = _fmt_value(impact.nominal_value)
        rmse_m_str = _fmt_float(rmse_m)
        rmse_p_str = _fmt_float(rmse_p)
        avg_imp_str = _fmt_float(avg_imp)

        line = lpad(string(rank), 5) * " | " * rpad(impact.parameter, 20) * " | " *
               lpad(nom_val, 15) * " | " * lpad(rmse_m_str, 14) * " | " *
               lpad(rmse_p_str, 14) * " | " * lpad(avg_imp_str, 14)

        if isfinite(avg_imp) && avg_imp > highlight_threshold_pct
            println(red * line * reset)
        else
            println(line)
        end
    end

    if length(impacts) > num_display
        println("... ($(length(impacts) - num_display) more parameters)")
    end
end
