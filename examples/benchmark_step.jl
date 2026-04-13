# -*- coding: utf-8 -*-

"""
Performance benchmark for AlphaPEM.
- 1 first run including compilation
- N measured runs (N=1 by default)
- Simple CSV export in results/benchmark/benchmark_step.csv
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Printf
using AlphaPEM.Config: SimulationConfig, StepParams
using AlphaPEM.Application: run_simulation

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

function write_csv(path, rows)
    open(path, "w") do io
        println(io, "timestamp,phase,run_index,time_s,memory_gb,gc_time_s")
        for r in rows
            println(io, string(
                r.timestamp, ",",
                r.phase, ",",
                r.run_index, ",",
                r.time_s, ",",
                r.memory_gb, ",",
                r.gc_time_s
            ))
        end
    end
end

function print_benchmark_report(rows, out_csv)
    measured_rows = filter(r -> r.phase == "measured", rows)
    first_row = only(filter(r -> r.phase == "first_run_with_compilation", rows))

    times = [r.time_s for r in measured_rows]
    memory_gb = [r.memory_gb for r in measured_rows]
    gc_times = [r.gc_time_s for r in measured_rows]

    println("\nBenchmark results")
    println("-----------------")
    println("GC (s) = time spent in garbage collection.")
    @printf("First launch (with compilation): time = %.6f s | memory = %.6f GB | GC = %.6f s\n", first_row.time_s, first_row.memory_gb, first_row.gc_time_s)

    println("\nMeasured runs (without first-launch compilation)")
    @printf("%-6s | %-10s | %-12s | %-10s\n", "Run", "Time (s)", "Memory (GB)", "GC (s)")
    println("---------------------------------------------------")
    for r in measured_rows
        @printf("%-6d | %-10.6f | %-12.6f | %-10.6f\n", r.run_index, r.time_s, r.memory_gb, r.gc_time_s)
    end
    println("---------------------------------------------------")
    @printf("Avg    | %-10.6f | %-12.6f | %-10.6f\n", sum(times) / length(times), sum(memory_gb) / length(memory_gb), sum(gc_times) / length(gc_times))
    @printf("Min    | %-10.6f | %-12.6f | %-10.6f\n", minimum(times), minimum(memory_gb), minimum(gc_times))
    @printf("Max    | %-10.6f | %-12.6f | %-10.6f\n", maximum(times), maximum(memory_gb), maximum(gc_times))
    println("CSV written: ", out_csv)
end

function main()
    runs = parse(Int, get(ENV, "BENCHMARK_RUNS", "1"))
    out_dir = joinpath(@__DIR__, "..", "results", "benchmark")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "benchmark_step.csv")

    rows = NamedTuple[]
    println("First launch (includes compilation)...")
    m0 = @timed run_simulation(make_step_config())
    push!(rows, (
        timestamp = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
        phase = "first_run_with_compilation",
        run_index = 0,
        time_s = m0.time,
        memory_gb = m0.bytes / 1.0e9,
        gc_time_s = m0.gctime
    ))

    for i in 1:runs
        println("Measured run ", i, "/", runs, "...")
        m = @timed run_simulation(make_step_config())
        push!(rows, (
            timestamp = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            phase = "measured",
            run_index = i,
            time_s = m.time,
            memory_gb = m.bytes / 1.0e9,
            gc_time_s = m.gctime
        ))
    end

    write_csv(out_csv, rows)
    print_benchmark_report(rows, out_csv)
end

main()


