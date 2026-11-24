# -*- coding: utf-8 -*-

# Usage: julia --project=. examples/run_extended_local_sensitivity_analysis.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

# ============================================================
# CONFIGURATION
# ============================================================

# Sobol balance is guaranteed for a power-of-two sample size. A non-power-of-two edit is rounded up automatically.
const N_CANDIDATES = 2^9
const CURRENT_DENSITY_GRID = Float64[
    0.0, 0.03, 0.05, 0.1, 0.2, 0.6, 1.0, 1.2, 1.5, 1.7, 1.9, 2.1, 2.4, 2.7, 3.0
]
const OPERATING_BOUNDS_OVERRIDES = Dict{Symbol, Tuple{Float64, Float64, Symbol}}()

cfg = ExtendedLocalSensitivityConfig(
    n_candidates = N_CANDIDATES,                                 # Sobol candidates generated before experimental screening.
    current_density_grid = CURRENT_DENSITY_GRID,                 # Shared current-density points for every operating case (A/cm²).
    current_density_zones = copy(DEFAULT_CURRENT_DENSITY_ZONES), # Magnitude summaries by current range; copy before customizing zone boundaries.
    pressure_gap = 0.2e5,                                        # Exact pressure difference: Pa_des - Pc_des (Pa).
    operating_bounds_overrides = OPERATING_BOUNDS_OVERRIDES,     # Optional partial bound overrides; empty keeps framework defaults.
    base_operating_conditions = nothing,                         # Use selected fuel-cell nominal conditions as the base; swept fields vary, others stay nominal.
    screening_datasets = [
        ScreeningDataset(:ZSW_nominal, 2024, :full),             # Candidate must match this experimental polarization curve.
        ScreeningDataset(:ZSW_Pa_1_61_Pc_1_41, 2024, :full),     # Candidate must also match this curve.
    ],                                                           # A candidate is retained only if it passes screening on every listed dataset.
)
validate_config(cfg)

# ============================================================
# PARAMETER SELECTION AND SOBOL DESIGN
# ============================================================

# Get the analysed physical parameters; integer parameters excluded by policy are reported separately.
params, fixed_parameters = select_parameters(cfg)

# Use the selected fuel cell's nominal parameters when an explicit nominal candidate is requested.
base_fc = AlphaPEM.Fuelcell.create_fuelcell(
    cfg.fuel_cell_type, cfg.voltage_zone; year = cfg.year, nb_gc = cfg.nb_gc)
nominal_unit = cfg.include_nominal ?
    nominal_unit_point(params, base_fc.physical_parameters) : nothing

# Generate candidates in the unit hypercube, then map them to the physical parameter bounds.
unit_design = generate_unit_design(N_CANDIDATES, length(params);
                                   seed = cfg.seed, nominal_unit = nominal_unit)
physical_design = unit_design_to_physical(unit_design, params)

# ============================================================
# CANDIDATE SCREENING
# ============================================================

# Bind each candidate to real AlphaPEM simulations and classify it against every screening dataset.
simulate_fn = alphapem_screening_simulator(params, cfg)
records = screen_candidates(physical_design, cfg.screening_datasets, cfg;
                            simulate_fn = simulate_fn)

# ============================================================
# FACTORIAL OPERATING-CASE LIBRARY AND BUDGET
# ============================================================

# Build all 729 operating cases; the configured common current grid is shared by the later simulations.
accepted_count = count(r -> r.status == :accepted, records)
operating_cases = generate_operating_case_library(cfg)
# Estimate future local-sensitivity evaluations; this is a warning only and never reduces the library.
projected_local_evaluations =
    accepted_count * length(params) * cfg.n_levels * length(operating_cases)

@info "Projected local-sensitivity run count; all 729 operating cases will be used." accepted_parameter_points = accepted_count parameter_count = length(params) finite_difference_levels = cfg.n_levels operating_case_count = length(operating_cases) projected_simulations = projected_local_evaluations

# ============================================================
# LOCAL-SENSITIVITY SIMULATION ENGINE
# ============================================================

# Open the live convergence figure before simulations start; callbacks add up to six cumulative
# states only after the corresponding candidate prefix is complete.
convergence_plot = create_sensitivity_convergence_plot(params, accepted_count)
# Open the voltage-sensitivity figure too; each finished candidate adds its gray points.
voltage_plot = create_voltage_sensitivity_plot(params)
# Run every accepted candidate against every operating case; the engine checkpoints completed blocks.
local_result = run_local_sensitivity_simulations(
    records, params, operating_cases, cfg;
    on_convergence_update = snapshot ->
        update_sensitivity_convergence_plot!(convergence_plot, snapshot),
    on_candidate_update = (_, magnitudes) ->
        add_voltage_sensitivity_candidate!(voltage_plot, magnitudes),
)
# Medians are drawn once, from every accepted candidate, after the run is complete.
final_snapshot = convergence_snapshot(
    local_result.sensitivities, params, length(local_result.candidates), 1)
finish_voltage_sensitivity_plot!(voltage_plot, final_snapshot.median_magnitude)
@info "Local sensitivity simulation complete" candidate_count = length(local_result.candidates) case_count = length(local_result.operating_cases) block_status_counts = Dict(status => count(==(status), local_result.statuses) for status in (:complete, :partial, :invalid)) execution_time_s = round(local_result.execution_time; digits = 2)

# ============================================================
# RUN SUMMARY
# ============================================================

println("Accepted $accepted_count / $(length(records)) candidates.")
println("Fixed parameters: $(join(string.(fixed_parameters), ", "))")

# In scripts, keep the browser page alive after all live figures are complete.
finish_sensitivity_convergence_plot!(convergence_plot)
