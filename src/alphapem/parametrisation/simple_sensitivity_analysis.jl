# -*- coding: utf-8 -*-

"""
Simple sensitivity analysis for AlphaPEM's PhysicalParams

This performs a sensitivity analysis on the fuel cell's physical parameters by:
1. Running a baseline polarization simulation with the nominal PhysicalParams of :ZSW_nominal
2. For each field of PhysicalParams:
   - Modifying it by ±20% (respecting integer/domain constraints, see `modify_physical_param`)
   - Running a polarization simulation in parallel, passing the modified PhysicalParams
     directly to SimulationConfig
   - Computing the RMSE vs. the nominal polarization curve
3. Ranking parameters by impact on the polarization curve
4. Saving results to results/simple_sensitivity/

Parallelization: Uses Distributed.jl with pmap for concurrent simulations. Since parameter
variations are passed as plain data (PhysicalParams instances) instead of patched source
files, workers never need to recompile or reload AlphaPEM.

This file is meant to be broadcast to every worker process (via `@everywhere include(...)`)
by the launcher script in examples/run_simple_sensitivity_analysis.jl, which is responsible
for activating the project and spawning the Distributed workers.
"""

using Distributed
using Printf
using AlphaPEM.Config: SimulationConfig, PolarizationCalibrationParams, NumericalParams, PhysicalParams
using AlphaPEM.Application: run_simulation
using AlphaPEM.Core.Models: _polarization_points_cali, _calculate_rmse
using AlphaPEM.Fuelcell: create_fuelcell

include(joinpath(@__DIR__, "sensitivity_analysis", "simple_sensitivity_helpers.jl"))

# Integer fields of PhysicalParams that can never go below 1.
const INTEGER_MIN1_FIELDS = (:nb_cell, :nb_channel_in_gc, :Kshape)

"""
   make_pola_config([base_config]; [physical_parameters])

Create the polarization SimulationConfig used throughout the sensitivity analysis.
If a `base_config` is provided, all its fields are preserved except `physical_parameters`
(which is overridden by the sensitivity variation) and `type_current` (which is forced to
`PolarizationCalibrationParams` so that every run uses the same calibration current points).

Uses `PolarizationCalibrationParams` (an empty `i_exp` gets auto-filled from the fuel
cell's calibration data by `create_current`) so that every run — nominal and modified —
is stepped through the exact same, known current densities. This removes the need to
match points between curves that don't otherwise land on identical current values.
"""
function make_pola_config(base_config::SimulationConfig = SimulationConfig();
                          physical_parameters::Union{Nothing, PhysicalParams} = nothing)
    current_params = PolarizationCalibrationParams(
        delta_t_ini = 30 * 60.0,      # (s). Initial stabilisation time.
        v_load = 0.01e4,              # (A.m-2.s-1). Loading rate.
        delta_t_break = 5 * 60.0,     # (s). Breaking time.
        i_exp = Float64[]             # Auto-filled from the fuel cell's calibration data.
    )

    return SimulationConfig(
        type_fuel_cell = base_config.type_fuel_cell,
        year = base_config.year,
        type_current = current_params,
        numerical_parameters = base_config.numerical_parameters,
        voltage_zone = base_config.voltage_zone,
        type_auxiliary = base_config.type_auxiliary,
        type_flow = base_config.type_flow,
        type_purge = base_config.type_purge,
        type_display = :no_display,
        display_timing = base_config.display_timing,
        state_scaling = base_config.state_scaling,
        physical_parameters = physical_parameters,
        operating_conditions = base_config.operating_conditions
    )
end

"""
    modify_physical_param(params, field, factor)

Return a copy of `params` with `field` scaled by `factor`, respecting the constraints of
PhysicalParams: integer fields stay integers, `nb_cell`/`nb_channel_in_gc`/`Kshape` can
never go below 1, and the capillary exponent `e` can only take the values 3, 4 or 5.
"""
function modify_physical_param(params::PhysicalParams, field::Symbol, factor::Float64)::PhysicalParams
    original_value = getfield(params, field)
    modified_value = original_value * factor

    if original_value isa Integer
        modified_value = round(Int64, modified_value)
        if field in INTEGER_MIN1_FIELDS
            modified_value = max(1, modified_value)
        elseif field === :e
            modified_value = clamp(modified_value, 3, 5)
        end
    end

    return PhysicalParams(;
        (f => (f === field ? modified_value : getfield(params, f)) for f in fieldnames(PhysicalParams))...
    )
end

"""Worker task: run one parameter variation and return its RMSE vs. the nominal polarization curve."""
function run_parameter_variation_task(nominal_params::PhysicalParams, field::Symbol, factor::Float64,
                                       i_exp::Vector{Float64}, Ucell_nominal::Vector{Float64},
                                       base_config::SimulationConfig)
    try
        modified_params = modify_physical_param(nominal_params, field, factor)
        cfg = make_pola_config(base_config, physical_parameters = modified_params)
        simu = run_simulation(cfg)
        return compute_rmse_from_nominal(simu, i_exp, Ucell_nominal)
    catch e
        return Inf
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Baseline Simulation
# ═══════════════════════════════════════════════════════════════════════════════

function run_baseline_simulation(i_exp::Vector{Float64}, base_config::SimulationConfig)
    println("=" ^ 80)
    println("WARM-UP: Running baseline polarization simulation...")
    println("=" ^ 80)

    cfg = make_pola_config(base_config)
    simu = run_simulation(cfg)

    pola_points = extract_polarization_points_cali(simu, i_exp)
    if isnothing(pola_points)
        error("Failed to extract baseline polarization points")
    end

    _, Ucell_nominal = pola_points
    println("\n✓ Baseline polarization extracted: $(length(Ucell_nominal)) points")
    return Ucell_nominal
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Project root resolution
# ═══════════════════════════════════════════════════════════════════════════════

"""Locate the project root by walking up from this file until a repo marker is found."""
function find_project_root()
    markers = [".git", "Project.toml", "Manifest.toml", "README.md"]
    parent = dirname(abspath(@__FILE__))
    while true
        any(ispath(joinpath(parent, m)) for m in markers) && return parent
        new_parent = dirname(parent)
        new_parent == parent && return pwd()
        parent = new_parent
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Main Execution
# ═══════════════════════════════════════════════════════════════════════════════

function run_simple_sensitivity_analysis(base_config::SimulationConfig = make_pola_config())
    println("\n")
    println("╔" * "═" ^ 78 * "╗")
    println("║" * " " ^ 20 * "ALPHAPEM SENSITIVITY ANALYSIS" * " " ^ 28 * "║")
    println("╚" * "═" ^ 78 * "╝")

    out_dir = joinpath(find_project_root(), "results", "simple_sensitivity_analysis")
    mkpath(out_dir)
    out_csv = generate_output_path(out_dir, "sensitivity_analysis.csv")

    # Step 1: Nominal PhysicalParams for the requested fuel cell. i_exp is fixed regardless
    # of PhysicalParams, so it can be shared as-is across the nominal run and every modified run.
    nominal_fc = create_fuelcell(base_config.type_fuel_cell, base_config.voltage_zone; year=base_config.year)
    nominal_params = nominal_fc.physical_parameters
    i_exp = Float64.(nominal_fc.pola_exp_data_cali.i_exp)

    # Step 2: Baseline simulation (warm-up)
    Ucell_nominal = run_baseline_simulation(i_exp, base_config)

    # Step 3: Build the (field, factor) task list over every PhysicalParams field
    tunable_fields = fieldnames(PhysicalParams)
    println("\n" * "=" ^ 80)
    println("PARAMETER VARIATIONS: Running $(length(tunable_fields)) parameters × 2 variations...")
    println("=" ^ 80)

    tasks = []
    for field in tunable_fields
        push!(tasks, (field, 0.95))  # -5%
        push!(tasks, (field, 1.05))  # +5%
    end

    # Step 4: Run variations in parallel using pmap
    println("Initializing $(nprocs() - 1) workers for parallel processing...")

    @everywhere begin
        const NOMINAL_PARAMS = $nominal_params
        const I_EXP = $i_exp
        const UCELL_NOMINAL = $Ucell_nominal
        const BASE_CONFIG = $base_config
    end

    results_list = pmap(
        (task) -> run_parameter_variation_task(NOMINAL_PARAMS, task[1], task[2], I_EXP, UCELL_NOMINAL, BASE_CONFIG),
        tasks,
        batch_size = 1,
        on_error = (task, e) -> (
            @warn "Task failed for $(task[1]) with factor $(task[2]): $e";
            Inf
        )
    )

    # Convert results to dictionary format
    results = Dict{String, Float64}()
    for (i, (field, factor)) in enumerate(tasks)
        key = factor == 0.80 ? "$(field)_minus20" : "$(field)_plus20"
        results[key] = results_list[i]

        @printf("\r[%3d/%3d] %-25s (%+3.0f%%)  RMSE: %.6f",
            i, length(tasks), string(field), (factor - 1) * 100, results_list[i])
        flush(stdout)
    end
    println("\r✓ All variations completed" * " " ^ 50)

    # Step 5: Compute impacts and rank
    params_dict = Dict{String, Any}(string(f) => getfield(nominal_params, f) for f in tunable_fields)
    impacts = compute_parameter_impacts(params_dict, results)

    # Step 6: Save results
    write_sensitivity_csv(out_csv, impacts)
    print_sensitivity_report(impacts, out_csv)

    println("\n✓ Sensitivity analysis completed successfully!\n")
end
