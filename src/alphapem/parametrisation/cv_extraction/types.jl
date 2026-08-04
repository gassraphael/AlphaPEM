# -*- coding: utf-8 -*-

"""
    CVExtraction types

Data structures for cyclic voltammetry (CV) parameter extraction.
"""

"""
    CVData

Raw cyclic voltammetry data.

# Fields
- `t::Vector{Float64}`: time in seconds
- `U::Vector{Float64}`: potential in V
- `I::Vector{Float64}`: current in A
"""
struct CVData
    t::Vector{Float64}
    U::Vector{Float64}
    I::Vector{Float64}
end

"""
    CVCycle

A single CV cycle with computed current density.

# Fields
- `t::Vector{Float64}`: time in seconds
- `U::Vector{Float64}`: potential in V
- `I::Vector{Float64}`: current in A
- `j::Vector{Float64}`: current density in A/cm²
"""
struct CVCycle
    t::Vector{Float64}
    U::Vector{Float64}
    I::Vector{Float64}
    j::Vector{Float64}
end

"""
    CVExtractionConfig

Configuration for CV parameter extraction.

# Fields
- `fuel_cell_type::Symbol`: fuel cell identifier (e.g. `:ZSW_GenStack_Pa_1_61_Pc_1_41`)
- `area_cm2::Float64`: geometric electrode area in cm²
- `scan_rate_vs::Float64`: scan rate in V·s⁻¹ (optional; if ≤ 0, inferred from data)
- `double_layer_limit_min::Float64`: lower potential bound for double-layer analysis (V vs reference electrode)
- `double_layer_limit_max::Float64`: upper potential bound for double-layer analysis (V vs reference electrode)
- `ohmic_drop_limit_min::Float64`: lower potential bound for ohmic-drop fit (V vs reference electrode)
- `ohmic_drop_limit_max::Float64`: upper potential bound for ohmic-drop fit (V vs reference electrode)
- `ecsa_limit_min::Float64`: lower potential bound for ECSA integration (V vs reference electrode)
- `ecsa_limit_max::Float64`: upper potential bound for ECSA integration (V vs reference electrode)
- `conversion_factor_uc_cm2::Float64`: charge for a monolayer of H on polycrystalline Pt (μC·cm⁻²)
- `covering_degree::Float64`: empirical H coverage on Pt (dimensionless, 0–1)
- `compensate_ohmic_drop::Bool`: subtract ohmic-drop slope from the curve
- `compensate_crossover::Bool`: subtract H₂ crossover current from the curve
- `ignore_cycles_for_mean::Vector{Int}`: cycle indices to exclude from the mean cycle
- `interpolation_factor::Int`: upsampling factor for cycle interpolation
"""
Base.@kwdef struct CVExtractionConfig
    fuel_cell_type::Symbol
    area_cm2::Float64
    scan_rate_vs::Float64 = -1.0
    double_layer_limit_min::Float64 = 0.35
    double_layer_limit_max::Float64 = 0.55
    ohmic_drop_limit_min::Float64 = 0.35
    ohmic_drop_limit_max::Float64 = 0.55
    ecsa_limit_min::Float64 = 0.05
    ecsa_limit_max::Float64 = 0.40
    conversion_factor_uc_cm2::Float64 = 210.0
    covering_degree::Float64 = 0.77
    compensate_ohmic_drop::Bool = true
    compensate_crossover::Bool = true
    ignore_cycles_for_mean::Vector{Int} = Int[]
    interpolation_factor::Int = 20
end

"""
    CVExtractionResult

Results of CV parameter extraction.

# Fields
- `fuel_cell_type::Symbol`: fuel cell identifier
- `ecsa_adsorption_cm2::Float64`: ECSA from H₂ adsorption peak (cm²_Pt·cm⁻²_electrode)
- `ecsa_desorption_cm2::Float64`: ECSA from H₂ desorption peak (cm²_Pt·cm⁻²_electrode)
- `crossover_a_cm2::Float64`: H₂ crossover current density (A·cm⁻²)
- `dlc_f_cm2::Float64`: double-layer capacitance (F·cm⁻²)
- `ohmic_drop_slope::Float64`: slope of the ohmic-drop line (A·V⁻¹·cm⁻²)
- `scan_rate_vs::Float64`: scan rate used for the extraction (V·s⁻¹)
- `mean_cycle::CVCycle`: corrected mean cycle used for the extraction
- `raw_mean_cycle::CVCycle`: original (uncorrected) mean cycle
- `anodic::CVCycle`: anodic branch of the corrected mean cycle
- `cathodic::CVCycle`: cathodic branch of the corrected mean cycle
"""
struct CVExtractionResult
    fuel_cell_type::Symbol
    ecsa_adsorption_cm2::Float64
    ecsa_desorption_cm2::Float64
    crossover_a_cm2::Float64
    dlc_f_cm2::Float64
    ohmic_drop_slope::Float64
    scan_rate_vs::Float64
    mean_cycle::CVCycle
    raw_mean_cycle::CVCycle
    anodic::CVCycle
    cathodic::CVCycle
end
