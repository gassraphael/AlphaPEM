# -*- coding: utf-8 -*-

"""
Benchmark `run_step` and `run_pola` across configurable `nb_gc` values.

Defaults:
- nb_gc values: 1, 5, 10
- one measured run per scenario/nb_gc (BENCHMARK_RUNS=1)

Mandatory warm-up sequence (not measured):
1) run_step with nb_gc = 1
2) run_pola with nb_gc = 1

Outputs:
- Console summary
- CSV file: results/benchmark/benchmark_nb_gc.csv
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Dates
using Printf
using AlphaPEM.Config: SimulationConfig, StepParams, PolarizationParams, NumericalParams
using AlphaPEM.Application: run_simulation

const DEFAULT_NB_GC_VALUES = [1, 5, 10]

function parse_nb_gc_values()
    raw = strip(get(ENV, "BENCHMARK_NB_GC", ""))
    isempty(raw) && return DEFAULT_NB_GC_VALUES

    vals = Int[]
    for token in split(raw, ',')
        parsed = tryparse(Int, strip(token))
        parsed === nothing && throw(ArgumentError("Invalid BENCHMARK_NB_GC token: $(token)"))
        parsed >= 1 || throw(ArgumentError("BENCHMARK_NB_GC values must be >= 1."))
        push!(vals, parsed)
    end
    return unique(vals)
end

function make_step_cfg(nb_gc::Int)
    step = StepParams(
        delta_t_ini = 30.0 * 60.0,
        delta_t_load = 30.0,
        delta_t_break = 2.0 * 60.0,
        i_ini = 1.0e4,
        i_step = 1.5e4,
    )
    return SimulationConfig(
        type_fuel_cell = :ZSW_GenStack,
        type_current = step,
        numerical_parameters = NumericalParams(nb_gc = nb_gc),
        voltage_zone = :full,
        type_auxiliary = :no_auxiliary,
        type_purge = :no_purge,
        type_display = :no_display,
        display_timing = :postrun,
    )
end

function make_pola_cfg(nb_gc::Int)
    pola = PolarizationParams(
        delta_t_ini = 30.0 * 60.0,
        delta_i = 0.05e4,
        v_load = 0.01e4,
        delta_t_break = 15.0 * 60.0,
        i_max = 3.0e4,
    )
    return SimulationConfig(
        type_fuel_cell = :ZSW_GenStack,
        type_current = pola,
        numerical_parameters = NumericalParams(nb_gc = nb_gc),
        voltage_zone = :full,
        type_auxiliary = :no_auxiliary,
        type_purge = :no_purge,
        type_display = :no_display,
        display_timing = :postrun,
    )
end

function timed_run(cfg::SimulationConfig)
    GC.gc()
    m = @timed try
        run_simulation(cfg)
        :ok
    catch err
        err
    end
    status = m.value == :ok ? "ok" : "fail"
    err_msg = status == "ok" ? "" : sprint(showerror, m.value)
    return (status = status, err = err_msg, time_s = m.time, alloc_gb = m.bytes / 1e9, gc_s = m.gctime)
end

function write_csv(path, rows)
    open(path, "w") do io
        println(io, "timestamp,phase,scenario,nb_gc,run_index,status,time_s,alloc_gb,gc_s,error")
        for r in rows
            println(io, string(
                r.timestamp, ",",
                r.phase, ",",
                r.scenario, ",",
                r.nb_gc, ",",
                r.run_index, ",",
                r.status, ",",
                r.time_s, ",",
                r.alloc_gb, ",",
                r.gc_s, ",",
                replace(r.err, "," => ";"),
            ))
        end
    end
end

function print_summary(rows)
    measured = filter(r -> r.phase == "measured", rows)

    println("\nBenchmark results")
    println("-----------------")
    for scenario in ("step", "pola")
        sub = filter(r -> r.scenario == scenario, measured)
        isempty(sub) && continue
        println("\nScenario: ", scenario)
        @printf("%-6s | %-6s | %-8s | %-10s | %-10s | %-10s\n", "nb_gc", "Run", "Status", "Time (s)", "Alloc (GB)", "GC (s)")
        println("----------------------------------------------------------------")
        for r in sub
            @printf("%-6d | %-6d | %-8s | %-10.4f | %-10.4f | %-10.4f\n", r.nb_gc, r.run_index, r.status, r.time_s, r.alloc_gb, r.gc_s)
        end
    end
end

function main()
    runs = parse(Int, get(ENV, "BENCHMARK_RUNS", "1"))
    nb_gc_values = parse_nb_gc_values()

    out_dir = joinpath(@__DIR__, "..", "results", "benchmark")
    mkpath(out_dir)
    out_csv = joinpath(out_dir, "benchmark_nb_gc.csv")

    rows = NamedTuple[]

    println("Warm-up #1: run_step with nb_gc = 1")
    warm_step = timed_run(make_step_cfg(1))
    push!(rows, (
        timestamp = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
        phase = "warmup",
        scenario = "step",
        nb_gc = 1,
        run_index = 0,
        status = warm_step.status,
        time_s = warm_step.time_s,
        alloc_gb = warm_step.alloc_gb,
        gc_s = warm_step.gc_s,
        err = warm_step.err,
    ))

    println("Warm-up #2: run_pola with nb_gc = 1")
    warm_pola = timed_run(make_pola_cfg(1))
    push!(rows, (
        timestamp = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
        phase = "warmup",
        scenario = "pola",
        nb_gc = 1,
        run_index = 0,
        status = warm_pola.status,
        time_s = warm_pola.time_s,
        alloc_gb = warm_pola.alloc_gb,
        gc_s = warm_pola.gc_s,
        err = warm_pola.err,
    ))

    for scenario in ("step", "pola")
        for nb_gc in nb_gc_values
            for i in 1:runs
                println("Measured run ", i, "/", runs, " | scenario=", scenario, " nb_gc=", nb_gc)
                result = scenario == "step" ? timed_run(make_step_cfg(nb_gc)) : timed_run(make_pola_cfg(nb_gc))
                push!(rows, (
                    timestamp = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
                    phase = "measured",
                    scenario = scenario,
                    nb_gc = nb_gc,
                    run_index = i,
                    status = result.status,
                    time_s = result.time_s,
                    alloc_gb = result.alloc_gb,
                    gc_s = result.gc_s,
                    err = result.err,
                ))
            end
        end
    end

    write_csv(out_csv, rows)
    print_summary(rows)
    println("\nCSV written: ", out_csv)

    failures = filter(r -> r.phase == "measured" && r.status != "ok", rows)
    isempty(failures) || error("Some benchmark runs failed. See CSV for details: $(out_csv)")
end

main()

