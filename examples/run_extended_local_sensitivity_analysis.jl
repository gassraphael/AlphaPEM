# -*- coding: utf-8 -*-

# Step 4 runner: generate and screen the scrambled Sobol candidates.
# Usage: julia --project=. examples/run_extended_local_sensitivity_analysis.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

# Sobol balance is guaranteed for a power-of-two sample size. A non-power-of-two edit is rounded up automatically.
const N_CANDIDATES = 2^8

cfg = ExtendedLocalSensitivityConfig(
    n_candidates = N_CANDIDATES,
    screening_datasets = [
        ScreeningDataset(:ZSW_nominal, 2024, :full),
        ScreeningDataset(:ZSW_Pa_1_61_Pc_1_41, 2024, :full),
    ],
)
validate_config(cfg)

params, fixed_parameters = select_parameters(cfg)
base_fc = AlphaPEM.Fuelcell.create_fuelcell(
    cfg.fuel_cell_type, cfg.voltage_zone; year = cfg.year, nb_gc = cfg.nb_gc)
nominal_unit = cfg.include_nominal ?
    nominal_unit_point(params, base_fc.physical_parameters) : nothing
unit_design = generate_unit_design(N_CANDIDATES, length(params);
                                   seed = cfg.seed, nominal_unit = nominal_unit)
physical_design = unit_design_to_physical(unit_design, params)
simulate_fn = alphapem_screening_simulator(params, cfg)
records = screen_candidates(physical_design, cfg.screening_datasets, cfg;
                            simulate_fn = simulate_fn)

println("Accepted $(count(r -> r.status == :accepted, records)) / $(length(records)) candidates.")
println("Fixed parameters: $(join(string.(fixed_parameters), ", "))")
