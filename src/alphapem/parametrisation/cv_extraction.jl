# -*- coding: utf-8 -*-

"""
    AlphaPEM.Parametrisation.CVExtraction

Cyclic voltammetry (CV) parameter extraction for PEM fuel cells.

This module provides a Julia port of the CV analysis routines originally
developed by ZSW in Matlab (`external/Parameter extraction - from ZSW/CV/`).
Given an experimental CV file, it corrects the curve for ohmic drop and H₂
crossover and extracts:

- `ECSA` — electrochemical surface area from H₂ adsorption/desorption peaks
- `H₂ crossover current density`
- `double-layer capacitance (DLC)`
- `ohmic-drop slope`

## Example

```julia
using AlphaPEM.Parametrisation.CVExtraction

cfg = CVExtractionConfig(
    fuel_cell_type = :ZSW_GenStack_Pa_1_61_Pc_1_41,
    area_cm2 = 280.0,
)

result = extract_cv_parameters("data/experimental/ZSW/cv/genstack1516-010c_2025-01-10_cv_cathode_cell01.txt", cfg)

@show result.ecsa_adsorption_cm2
@show result.ecsa_desorption_cm2
@show result.crossover_a_cm2
@show result.dlc_f_cm2
```
"""
module CVExtraction

using Printf

include("cv_extraction/types.jl")
include("cv_extraction/io.jl")
include("cv_extraction/processing.jl")
include("cv_extraction/analysis.jl")
include("cv_extraction/plotting.jl")

export CVData,
       CVCycle,
       CVExtractionConfig,
       CVExtractionResult,
       read_cv_file,
       clean_cv_data,
       resample_cv_data,
       split_cycles,
       mean_cycle,
       interpolate_cycle,
       smallest_distance,
       ohmic_drop_fit,
       area_integral,
       extract_parameters,
       extract_cv_parameters,
       plot_cv_analysis

"""
    extract_cv_parameters(path::String, cfg::CVExtractionConfig) -> CVExtractionResult

High-level entry point: read a CV file, process cycles, compute the mean cycle,
and extract physical parameters.
"""
function extract_cv_parameters(path::String, cfg::CVExtractionConfig)::CVExtractionResult
    raw = read_cv_file(path)
    cleaned = clean_cv_data(raw)
    resampled, scan_rate = resample_cv_data(cleaned; scan_rate = cfg.scan_rate_vs)

    if cfg.scan_rate_vs <= 0.0
        cfg = CVExtractionConfig(
            fuel_cell_type = cfg.fuel_cell_type,
            area_cm2 = cfg.area_cm2,
            scan_rate_vs = scan_rate,
            double_layer_limit_min = cfg.double_layer_limit_min,
            double_layer_limit_max = cfg.double_layer_limit_max,
            ohmic_drop_limit_min = cfg.ohmic_drop_limit_min,
            ohmic_drop_limit_max = cfg.ohmic_drop_limit_max,
            ecsa_limit_min = cfg.ecsa_limit_min,
            ecsa_limit_max = cfg.ecsa_limit_max,
            conversion_factor_uc_cm2 = cfg.conversion_factor_uc_cm2,
            covering_degree = cfg.covering_degree,
            compensate_ohmic_drop = cfg.compensate_ohmic_drop,
            compensate_crossover = cfg.compensate_crossover,
            cycle_for_extraction = cfg.cycle_for_extraction,
            ignore_cycles_for_mean = cfg.ignore_cycles_for_mean,
            interpolation_factor = cfg.interpolation_factor,
        )
    end

    cycles = split_cycles(resampled)
    if isempty(cycles)
        error("No CV cycles could be detected in $path")
    end

    interpolated = [interpolate_cycle(c, scan_rate, cfg.interpolation_factor) for c in cycles]

    # Select the cycle to analyse
    if cfg.cycle_for_extraction == 0
        selected = mean_cycle(interpolated, cfg.ignore_cycles_for_mean)
    else
        if cfg.cycle_for_extraction > length(interpolated)
            error("Requested cycle $(cfg.cycle_for_extraction) but only $(length(interpolated)) cycle(s) available")
        end
        selected = interpolated[cfg.cycle_for_extraction]
    end

    # Convert to current density
    selected_density = CVCycle(selected.t, selected.U, selected.I, selected.I ./ cfg.area_cm2)

    return extract_parameters(selected_density, cfg)
end

end  # module CVExtraction
