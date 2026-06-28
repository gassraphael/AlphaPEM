# -*- coding: utf-8 -*-

"""
Benchmark `run_step`, `run_pola` and `run_eis` across configurable `nb_gc` values.

Defaults:
- nb_gc values: 1, 5, 10
- 5 measured run per scenario/nb_gc (BENCHMARK_RUNS=5)

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
using Statistics
using AlphaPEM.Config: SimulationConfig, StepParams, PolarizationParams, EISParams, NumericalParams
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
        type_flow = :co_flow, # :co_flow, :counter_flow.
        type_purge = :no_purge,
        type_display = :no_display,
        display_timing = :postrun,
    )
end

function make_pola_cfg(nb_gc::Int)
    pola = PolarizationParams(
        delta_t_ini = 30.0 * 60.0,
        di_step = 0.05e4,
        v_load = 0.01e4,
        delta_t_break = 5 * 60.0,
        i_max = 3.0e4,
    )
    return SimulationConfig(
        type_fuel_cell = :ZSW_GenStack,
        type_current = pola,
        numerical_parameters = NumericalParams(nb_gc = nb_gc),
        voltage_zone = :full,
        type_auxiliary = :no_auxiliary,
        type_flow = :co_flow, # :co_flow, :counter_flow.
        type_purge = :no_purge,
        type_display = :no_display,
        display_timing = :postrun,
    )
end

function make_eis_cfg(nb_gc::Int)
    current_params = EISParams(
        i_EIS = 1.0e4,        # (A/m²). Parameters for the EIS curve.
        ratio = 5.0 / 100.0,  # (.). Parameters for the EIS curve.
        f_power_min = -3.0,   # (.). Power of the minimum frequency for the EIS current density function.
        f_power_max = 5.0,    # (.). Power of the maximum frequency for the EIS current density function.
        nb_f = 90,            # (.). Number of frequencies tested for the EIS current density function.
        nb_points = 50,       # (.). Number of points calculated per specific period for the EIS current density function.
    )

    return SimulationConfig(
        type_fuel_cell = :ZSW_GenStack,
        type_current = current_params,
        numerical_parameters = NumericalParams(nb_gc = nb_gc),
        voltage_zone = :full,
        type_auxiliary = :no_auxiliary,
        type_flow = :co_flow, # :co_flow, :counter_flow.
        type_purge = :no_purge,
        type_display = :no_display,
        display_timing = :live
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

function compute_statistics(rows)
    measured = filter(r -> r.phase == "measured", rows)
    stats = Dict()

    for scenario in ("step", "pola", "eis")
        scenario_rows = filter(r -> r.scenario == scenario, measured)
        isempty(scenario_rows) && continue
        for nb_gc in sort(unique(r.nb_gc for r in scenario_rows))
            subset = filter(r -> r.scenario == scenario && r.nb_gc == nb_gc, measured)
            if !isempty(subset)
                stats[(scenario, nb_gc)] = (
                    avg_time = mean(r.time_s for r in subset),
                    avg_alloc = mean(r.alloc_gb for r in subset),
                    avg_gc = mean(r.gc_s for r in subset),
                    runs = length(subset),
                )
            end
        end
    end
    return stats
end

function print_statistics_table(rows)
    stats = compute_statistics(rows)
    isempty(stats) && return

    println("\n=== Summary Statistics (Average across all runs) ===")
    @printf("%-8s | %-6s | %-8s | %-10s | %-10s | %-10s\n", "Scenario", "nb_gc", "Runs", "Time (s)", "Alloc (GB)", "GC (s)")
    println("-----------------------------------------------------------")

    for scenario in ("step", "pola", "eis")
        for nb_gc in sort(collect(Set(r.nb_gc for r in rows if r.scenario == scenario && r.phase == "measured")))
            if haskey(stats, (scenario, nb_gc))
                s = stats[(scenario, nb_gc)]
                @printf("%-8s | %-6d | %-8d | %-10.4f | %-10.4f | %-10.4f\n",
                    scenario, nb_gc, s.runs, s.avg_time, s.avg_alloc, s.avg_gc)
            end
        end
    end
end

function write_csv(path, rows)
    stats = compute_statistics(rows)

    open(path, "w") do io
        if !isempty(stats)
            println(io, "# SUMMARY STATISTICS (Average across all runs)")
            println(io, "scenario,nb_gc,avg_time_s,avg_alloc_gb,avg_gc_s")
            for scenario in ("step", "pola", "eis")
                for nb_gc in sort(collect(Set(r.nb_gc for r in rows if r.scenario == scenario && r.phase == "measured")))
                    if haskey(stats, (scenario, nb_gc))
                        s = stats[(scenario, nb_gc)]
                        println(io, "$(scenario),$(nb_gc),$(s.avg_time),$(s.avg_alloc),$(s.avg_gc)")
                    end
                end
            end
            println(io, "")
        end

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
#    for scenario in ("step", "pola", "eis")
        sub = filter(r -> r.scenario == scenario, measured)
        isempty(sub) && continue
        println("\nScenario: ", scenario)
        @printf("%-6s | %-6s | %-8s | %-10s | %-10s | %-10s\n", "nb_gc", "Run", "Status", "Time (s)", "Alloc (GB)", "GC (s)")
        println("----------------------------------------------------------------")
        for r in sub
            @printf("%-6d | %-6d | %-8s | %-10.4f | %-10.4f | %-10.4f\n", r.nb_gc, r.run_index, r.status, r.time_s, r.alloc_gb, r.gc_s)
        end
    end

    print_statistics_table(rows)
end

function main()
    runs = parse(Int, get(ENV, "BENCHMARK_RUNS", "5"))
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

    println("Warm-up #3: run_eis with nb_gc = 1")
    warm_eis = timed_run(make_eis_cfg(1))
    push!(rows, (
        timestamp = Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
        phase = "warmup",
        scenario = "eis",
        nb_gc = 1,
        run_index = 0,
        status = warm_eis.status,
        time_s = warm_eis.time_s,
        alloc_gb = warm_eis.alloc_gb,
        gc_s = warm_eis.gc_s,
        err = warm_eis.err,
    ))

    for scenario in ("step", "pola", "eis")
        for nb_gc in nb_gc_values
            for i in 1:runs
                println("Measured run ", i, "/", runs, " | scenario=", scenario, " nb_gc=", nb_gc)
                result = if scenario == "step"
                    timed_run(make_step_cfg(nb_gc))
                elseif scenario == "pola"
                    timed_run(make_pola_cfg(nb_gc))
                else
                    timed_run(make_eis_cfg(nb_gc))
                end
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
