# -*- coding: utf-8 -*-

"""
CPU profiling for AlphaPEM step simulation.
- N warm-up runs (compilation/caches)
- N profiled runs aggregated in one profile buffer
- Flat profile report export in results/profiling/profile_step.txt
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Profile
using Printf
using AlphaPEM.Config: SimulationConfig, StepParams
using AlphaPEM.Application: run_simulation

function collect_flat_profile_text()
    buffer = IOBuffer()
    Profile.print(buffer; format = :flat, sortedby = :count)
    return String(take!(buffer))
end

function is_actionable_hotspot(file, function_name)
    orchestration_files = (
        "application/run_simulation.jl",
        "application/run_simulation_modules.jl",
        "examples/profile_step.jl",
        "core/models/AlphaPEM.jl"
    )
    orchestration_functions = (
        "run_simulation",
        "launch_AlphaPEM",
        "simulate_model!",
        "#simulate_model!##0",
        "recovery!",
        "main()",
        "macro expansion"
    )

    any(occursin(fragment, file) for fragment in orchestration_files) && return false
    any(occursin(fragment, function_name) for fragment in orchestration_functions) && return false
    return true
end

function extract_top_hotspots(flat_profile_text; n=5)
    hotspots = NamedTuple[]

    for line in split(flat_profile_text, '\n')
        stripped = strip(line)
        isempty(stripped) && continue
        occursin("@AlphaPEM/", line) || continue

        m = match(r"^\s*(\d+)\s+(\d+)\s+(.+?)\s+(\S+)\s+(.+)$", line)
        m === nothing && continue

        file = strip(m.captures[3])
        function_name = strip(m.captures[5])
        is_actionable_hotspot(file, function_name) || continue

        push!(hotspots, (
            count = parse(Int, m.captures[1]),
            overhead = parse(Int, m.captures[2]),
            file = file,
            line = m.captures[4],
            function_name = function_name
        ))
    end

    aggregated = Dict{Tuple{String, String}, NamedTuple}()
    for hotspot in hotspots
        key = (hotspot.file, hotspot.function_name)
        if haskey(aggregated, key)
            previous = aggregated[key]
            aggregated[key] = (
                count = previous.count + hotspot.count,
                overhead = previous.overhead + hotspot.overhead,
                file = hotspot.file,
                line = previous.line == hotspot.line ? previous.line : string(previous.line, ",", hotspot.line),
                function_name = hotspot.function_name
            )
        else
            aggregated[key] = (
                count = hotspot.count,
                overhead = hotspot.overhead,
                file = hotspot.file,
                line = hotspot.line,
                function_name = hotspot.function_name
            )
        end
    end

    aggregated_hotspots = collect(values(aggregated))
    sort!(aggregated_hotspots; by = h -> (-h.count, -h.overhead, h.file, h.function_name))
    return first(aggregated_hotspots, min(n, length(aggregated_hotspots)))
end

function format_line_span(line_text)
    line_parts = split(line_text, ',')
    unique_parts = unique(strip.(line_parts))
    isempty(unique_parts) && return "?"
    length(unique_parts) == 1 && return unique_parts[1]
    return string(unique_parts[1], " (+", length(unique_parts) - 1, " more)")
end

function compact_path(file; maxlen=46)
    normalized = replace(file, "\\" => "/")
    trimmed = replace(normalized, "@AlphaPEM/" => "")
    length(trimmed) <= maxlen && return trimmed
    return string("...", last(trimmed, maxlen - 3))
end

function print_hotspots(io, hotspots, total_profiled_time)
    total_samples = sum(h.count for h in hotspots)
    println(io, rpad("Rank", 6), lpad("Samples", 10), lpad("Share", 9),
            lpad("Est. s", 10), "  ", rpad("Location", 56), "Function")
    println(io, repeat("-", 122))

    isempty(hotspots) && begin
        println(io, "No actionable AlphaPEM hotspots found.")
        return
    end

    for (i, hotspot) in enumerate(hotspots)
        pct = total_samples == 0 ? 0.0 : 100 * hotspot.count / total_samples
        est_seconds = total_profiled_time * (pct / 100)
        location = string(compact_path(hotspot.file), ":", format_line_span(hotspot.line))
        sample_txt = string(hotspot.count)
        pct_txt = @sprintf("%.1f%%", pct)
        est_txt = @sprintf("%.4f", est_seconds)
        println(io,
                rpad(string(i), 6),
                lpad(sample_txt, 10),
                lpad(pct_txt, 9),
                lpad(est_txt, 10),
                "  ",
                rpad(location, 56),
                hotspot.function_name)
    end
end

function make_step_config()
    current_params = StepParams(
        delta_t_ini = 30.0 * 60.0,
        delta_t_load = 30.0,
        delta_t_break = 2.0 * 60.0,
        i_ini = 1.0e4,
        i_step = 2.0e4
    )

    return SimulationConfig(
        type_fuel_cell = :ZSW_GenStack,
        type_current = current_params,
        voltage_zone = :full,
        type_auxiliary = :no_auxiliary,
        type_purge = :no_purge,
        type_display = :no_display,
        type_plot = :fixed
    )
end

function main()
    warmup_runs = parse(Int, get(ENV, "PROFILE_WARMUP_RUNS", "1"))
    profiled_runs = parse(Int, get(ENV, "PROFILE_RUNS", "3"))

    out_dir = joinpath(@__DIR__, "..", "results", "profiling")
    mkpath(out_dir)
    out_txt = joinpath(out_dir, "profile_step.txt")
    cfg = make_step_config()

    for i in 1:warmup_runs
        println("Warm-up run ", i, "/", warmup_runs, "...")
        run_simulation(cfg)
    end

    println("Profiled runs (aggregated samples)...")
    Profile.clear()
    elapsed_runs = Float64[]
    @profile begin
        for i in 1:profiled_runs
            println("Profiled run ", i, "/", profiled_runs, "...")
            t = @elapsed run_simulation(cfg)
            push!(elapsed_runs, t)
        end
    end

    total_profiled_time = sum(elapsed_runs)
    mean_profiled_time = total_profiled_time / length(elapsed_runs)
    min_profiled_time = minimum(elapsed_runs)
    max_profiled_time = maximum(elapsed_runs)
    flat_profile_text = collect_flat_profile_text()
    top_hotspots = extract_top_hotspots(flat_profile_text; n=5)

    open(out_txt, "w") do io
        println(io, "AlphaPEM step profile")
        println(io, "timestamp: ", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
        println(io, "warmup_runs: ", warmup_runs)
        println(io, "profiled_runs: ", profiled_runs)
        println(io, "total_profiled_time_s: ", total_profiled_time)
        println(io, "mean_profiled_time_s: ", mean_profiled_time)
        println(io, "min_profiled_time_s: ", min_profiled_time)
        println(io, "max_profiled_time_s: ", max_profiled_time)
        println(io, "simulation_config: ", cfg)
        println(io, "")
        println(io, "Top 5 actionable AlphaPEM hotspots (aggregated by function, orchestration excluded):")
        println(io, "Columns: Rank | Samples | Share (within top 5) | Est. s (Share * total_profiled_time_s) | Location | Function")
        print_hotspots(io, top_hotspots, total_profiled_time)
        println(io, "")
        println(io, "Top stack samples (flat view, sorted by count):")
        print(io, flat_profile_text)
    end

    println("Profile report written: ", out_txt)
    println("Profiled runs: ", profiled_runs)
    println("Total profiled time (s): ", total_profiled_time)
    println("Mean profiled run time (s): ", mean_profiled_time)
    println("Min/Max profiled run time (s): ", min_profiled_time, " / ", max_profiled_time)
    println("Top 5 actionable hotspots:")
    println("Columns: Rank | Samples | Share (within top 5) | Est. s (Share * total_profiled_time_s) | Location | Function")
    print_hotspots(stdout, top_hotspots, total_profiled_time)
end

main()

