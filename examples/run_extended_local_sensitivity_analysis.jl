# -*- coding: utf-8 -*-

# Steps 4-5: generate and screen Sobol candidates, then build the factorial case library.
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
    n_candidates = N_CANDIDATES,
    current_density_grid = CURRENT_DENSITY_GRID,
    pressure_gap = 0.2e5,
    operating_bounds_overrides = OPERATING_BOUNDS_OVERRIDES,
    # Set this to an OperatingConditions value to override the selected fuel cell's nominal base.
    base_operating_conditions = nothing,
    screening_datasets = [
        ScreeningDataset(:ZSW_nominal, 2024, :full),
        ScreeningDataset(:ZSW_Pa_1_61_Pc_1_41, 2024, :full),
    ],
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

@info "Factorial operating-case library generated" case_count = length(operating_cases) current_points = length(cfg.current_density_grid)
@warn "Projected local-sensitivity run count; all 729 operating cases will be used." accepted_parameter_points = accepted_count parameter_count = length(params) finite_difference_levels = cfg.n_levels operating_case_count = length(operating_cases) projected_simulations = projected_local_evaluations

# ============================================================
# RUN SUMMARY
# ============================================================

println("Accepted $accepted_count / $(length(records)) candidates.")
println("Fixed parameters: $(join(string.(fixed_parameters), ", "))")
