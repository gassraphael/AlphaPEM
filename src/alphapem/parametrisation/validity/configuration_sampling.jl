# -*- coding: utf-8 -*-

"""
    ConfigurationSampling

Module for generating parameter samples that cover the undetermined-parameter space
of AlphaPEM.

## Undetermined parameters

The parameters that cannot be directly measured and must be inferred are:

| Parameter     | Description                                           | Unit   |
|---------------|-------------------------------------------------------|--------|
| `Hacl`        | Anode catalyst-layer thickness                        | m      |
| `Hccl`        | Cathode catalyst-layer thickness                      | m      |
| `Hmem`        | Membrane thickness                                    | m      |
| `Hgdl`        | Gas-diffusion-layer thickness                         | m      |
| `Hmpl`        | Microporous-layer thickness                           | m      |
| `epsilon_gdl` | GDL porosity                                          | —      |
| `e`           | Capillary exponent (integer)                          | —      |
| `Re`          | Electron-conduction resistance                        | Ω m²   |
| `i0_c_ref`    | Reference cathode exchange current density            | A m⁻²  |
| `kappa_co`    | Crossover correction coefficient                      | —      |
| `kappa_c`     | Overpotential correction exponent                     | —      |
| `K_l_ads`     | Liquid/vapor water-sorption rate ratio (`:full` only) | —      |
| `K_O2_ad_Pt`  | O₂ adsorption resistance coefficient (`:full` only)   | —      |


## Sampling strategy

**Latin Hypercube Sampling (LHS)** is used by default.  LHS divides each parameter
dimension into equally-probable strata and draws exactly one sample per stratum,
giving better space-filling properties than simple random sampling for the same budget.

# Exports
- `ParameterBound`       — Bounds and metadata for one parameter
- `ParameterBounds`      — Full set of bounds for a fuel-cell type / voltage zone
- `SamplingConfig`       — Sampling options (n, method, seed)
- `bounds_for_fuel_cell` — Build `ParameterBounds` for a supported fuel-cell type
- `generate_lhs_samples` — Draw samples within given bounds
- `apply_bounds_to_params` — Map a sample vector onto a `PhysicalParams` struct
- `get_reference_config` — Return the nominal reference configuration
"""
module ConfigurationSampling

# `LatinHypercubeSampling` provides the LHS engine; `Random` is kept only for
# RNG seeding / the optional Monte-Carlo branch.
using Random
using LatinHypercubeSampling: randomLHC, scaleLHC

using AlphaPEM.Config: PhysicalParams
using AlphaPEM.Fuelcell: create_fuelcell

export ParameterBound,
       ParameterBounds,
       SamplingConfig,
       bounds_for_fuel_cell,
       generate_lhs_samples,
       apply_bounds_to_params,
       get_reference_config

# ─────────────────────────────────────────────────────────────────────────────
# DATA STRUCTURES
# ─────────────────────────────────────────────────────────────────────────────

"""
    ParameterBound

Bounds and metadata for a single undetermined parameter.

# Fields
- `name::Symbol`: Parameter identifier (e.g. `:Hacl`, `:Re`, `:epsilon_gdl`).
- `min::Float64`: Lower bound (inclusive).
- `max::Float64`: Upper bound (inclusive).
- `type::Symbol`: `:real` for continuous parameters, `:int` for integer ones (e.g. `:e`).
- `unit::String`: Physical unit string used for display and export.
- `description::String`: Short human-readable description.
"""
struct ParameterBound
    name::Symbol
    min::Float64
    max::Float64
    type::Symbol          # :real or :int
    unit::String
    description::String
end


"""
    ParameterBounds

Complete set of bounds for all undetermined parameters of a given fuel-cell type.

Built by `bounds_for_fuel_cell` and consumed by `generate_lhs_samples` and
`apply_bounds_to_params`.

# Fields
- `bounds::Vector{ParameterBound}`: Bound specifications in calibration order.
- `fuel_cell_type::Symbol`: Fuel-cell type (e.g. `:ZSW_GenStack`).
- `voltage_zone::Symbol`: `:before_voltage_drop` or `:full`.
- `n_params::Int`: Number of undetermined parameters (= `length(bounds)`).

# Example
```julia
pb = bounds_for_fuel_cell(:ZSW_GenStack, :full)
println(pb.n_params)        # 13
println(pb.bounds[1].name)  # :Hacl
```
"""
struct ParameterBounds
    bounds::Vector{ParameterBound}
    fuel_cell_type::Symbol
    voltage_zone::Symbol
    n_params::Int
end


"""
    SamplingConfig

Options controlling how parameter samples are generated.

# Fields
- `n_samples::Int`: Number of configurations to draw. Default `10000`.
- `method::Symbol`: `:lhs` (Latin Hypercube, recommended), `:random` (Monte Carlo),
  or `:grid` (full factorial). Default `:lhs`.
- `seed::Union{Int,Nothing}`: Random seed for reproducibility. Default `42`.
- `include_reference::Bool`: Prepend the reference configuration as sample row 1,
  ensuring the known-valid point is always present. Default `true`.

# Example
```julia
cfg = SamplingConfig(n_samples = 500, seed = 0)
```
"""
Base.@kwdef struct SamplingConfig
    n_samples::Int            = 10000
    method::Symbol            = :lhs
    seed::Union{Int, Nothing} = 42
    include_reference::Bool   = true
end


# ─────────────────────────────────────────────────────────────────────────────
# BOUNDS EXTRACTION
# ─────────────────────────────────────────────────────────────────────────────

"""
    bounds_for_fuel_cell(fuel_cell_type, voltage_zone = :full) -> ParameterBounds

Return the undetermined-parameter bounds for a given fuel-cell type and voltage zone.

Bounds are adapted from `parameter_bounds_for_calibration()` in `calibration_modules.jl`,
restricted to the undetermined physical parameters (operating inputs are excluded).

# Supported `fuel_cell_type` values
- `:ZSW_GenStack` (and pressure/temperature variants)
- `:EH_31_1_5`, `:EH_31_2_0`, `:EH_31_2_25`, `:EH_31_2_5`
"""
function bounds_for_fuel_cell(fuel_cell_type::Symbol,
                               voltage_zone::Symbol = :full)::ParameterBounds
    voltage_zone in (:full, :before_voltage_drop) ||
        throw(ArgumentError("voltage_zone must be :full or :before_voltage_drop (got $voltage_zone)"))

    bounds = ParameterBound[]

    # Helper for concise insertion.
    _push!(name::Symbol, mn::Real, mx::Real;
           type::Symbol = :real,
           unit::String = "",
           description::String = "") = push!(bounds,
                                            ParameterBound(name, Float64(mn), Float64(mx), type, unit, description))

    # Bounds are adapted from `parameter_bounds_for_calibration()`.
    if fuel_cell_type in (
        :ZSW_GenStack,
        :ZSW_GenStack_Pa_1_61_Pc_1_41,
        :ZSW_GenStack_Pa_2_01_Pc_1_81,
        :ZSW_GenStack_Pa_2_4_Pc_2_2,
        :ZSW_GenStack_Pa_2_8_Pc_2_6,
        :ZSW_GenStack_T_62,
        :ZSW_GenStack_T_76,
        :ZSW_GenStack_T_84,
    )
        _push!(:Hacl, 5e-6, 15e-6; unit = "m", description = "Anode catalyst-layer thickness")
        _push!(:Hccl, 5e-6, 20e-6; unit = "m", description = "Cathode catalyst-layer thickness")
        _push!(:Hmem, 5e-6, 30e-6; unit = "m", description = "Membrane thickness")
        _push!(:Hgdl, 100e-6, 150e-6; unit = "m", description = "Gas-diffusion-layer thickness")
        _push!(:Hmpl, 40e-6, 100e-6; unit = "m", description = "Microporous-layer thickness")
        _push!(:epsilon_gdl, 0.5, 0.9; unit = "—", description = "GDL porosity")
        _push!(:e, 3, 5; type = :int, unit = "—", description = "Capillary exponent")
        _push!(:Re, 5e-8, 5e-6; unit = "Ω·m²", description = "Electron-conduction resistance")
        _push!(:i0_c_ref, 1e-1, 100.0; unit = "A·m⁻²", description = "Reference cathode exchange current density")
        _push!(:kappa_co, 0.01, 40.0; unit = "—", description = "Crossover correction coefficient")
        _push!(:kappa_c, 0.25, 4.0; unit = "—", description = "Overpotential correction exponent")

        if voltage_zone == :full
            _push!(:K_l_ads, 1.0, 100.0; unit = "—", description = "Liquid/vapor water-sorption rate ratio")
            _push!(:K_O2_ad_Pt, 1.0, 10.0; unit = "—", description = "O₂ adsorption resistance coefficient")
        end

    elseif fuel_cell_type in (:EH_31_1_5, :EH_31_2_0, :EH_31_2_25, :EH_31_2_5)
        # The EH-31 calibration assumes Hccl == Hacl. We therefore only sample Hacl and
        # will mirror Hccl in `apply_bounds_to_params`.
        _push!(:Hacl, 8e-6, 20e-6; unit = "m", description = "Anode/cathode catalyst-layer thickness (Hccl = Hacl)")
        _push!(:Hmem, 15e-6, 50e-6; unit = "m", description = "Membrane thickness")
        _push!(:epsilon_gdl, 0.40, 0.95; unit = "—", description = "GDL porosity")
        _push!(:e, 3, 5; type = :int, unit = "—", description = "Capillary exponent")
        _push!(:Re, 5e-7, 5e-6; unit = "Ω·m²", description = "Electron-conduction resistance")
        _push!(:i0_c_ref, 1e-1, 100.0; unit = "A·m⁻²", description = "Reference cathode exchange current density")
        _push!(:kappa_co, 0.01, 40.0; unit = "—", description = "Crossover correction coefficient")
        _push!(:kappa_c, 0.25, 4.0; unit = "—", description = "Overpotential correction exponent")
        if voltage_zone == :full
            _push!(:K_O2_ad_Pt, 1.0, 10.0; unit = "—", description = "O₂ adsorption resistance coefficient")
        end

    else
        throw(ArgumentError("Unsupported fuel_cell_type: $fuel_cell_type"))
    end

    # Sanity checks.
    for b in bounds
        b.min <= b.max || throw(ArgumentError("Invalid bounds for $(b.name): min=$(b.min) > max=$(b.max)"))
        b.type in (:real, :int) || throw(ArgumentError("Invalid type for $(b.name): $(b.type)"))
    end

    return ParameterBounds(bounds, fuel_cell_type, voltage_zone, length(bounds))
end


# ─────────────────────────────────────────────────────────────────────────────
# SAMPLE GENERATION
# ─────────────────────────────────────────────────────────────────────────────

"""
    generate_lhs_samples(pb, cfg = SamplingConfig()) -> Matrix{Float64}

Draw a sample matrix using the strategy specified by `cfg`.

Each row is one configuration; columns correspond to parameters in the order defined
by `pb.bounds`. Integer parameters (type `:int`, e.g. `e`) are rounded to the nearest
integer after sampling.

# Arguments
- `pb::ParameterBounds`: Bounds and metadata.
- `cfg::SamplingConfig`: Sampling options.

# Returns
`Matrix{Float64}` of size `(cfg.n_samples, pb.n_params)`.

# Example
```julia
pb  = bounds_for_fuel_cell(:ZSW_GenStack, :full)
cfg = SamplingConfig(n_samples = 500, seed = 1)
X   = generate_lhs_samples(pb, cfg)
size(X)   # (500, 13)
```
"""
function generate_lhs_samples(pb::ParameterBounds,
                               cfg::SamplingConfig = SamplingConfig())::Matrix{Float64}
    n_samples = cfg.n_samples
    n_params = pb.n_params
    n_samples > 0 || throw(ArgumentError("n_samples must be > 0 (got $n_samples)"))
    n_params > 0 || throw(ArgumentError("No parameters to sample (pb.n_params == 0)"))

    # RNG for reproducibility (LHS + Monte Carlo).
    rng = cfg.seed === nothing ? Random.default_rng() : MersenneTwister(cfg.seed)

    # Build the LHS or Monte Carlo plan, then scale it to physical parameter bounds.
    X = if cfg.method == :lhs
        raw = randomLHC(rng, n_samples, n_params)
        scaleLHC(raw, [(b.min, b.max) for b in pb.bounds])
    elseif cfg.method == :random
        # Map uniform [0,1] draws to each parameter's physical bounds for Monte Carlo sampling.
        raw = rand(rng, n_samples, n_params)
        for (j, b) in enumerate(pb.bounds)
            span = b.max - b.min
            @inbounds for i in 1:n_samples
                raw[i, j] = b.min + raw[i, j] * span
            end
        end
        raw
    else
        throw(ArgumentError("Unsupported sampling method: $(cfg.method). Supported: :lhs, :random"))
    end
    X = Matrix{Float64}(X)

    # Enforce integer parameters (e.g. e) within bounds.
    for (j, b) in enumerate(pb.bounds)
        if b.type == :int
            @inbounds for i in 1:n_samples
                X[i, j] = clamp(round(X[i, j]), b.min, b.max)
            end
        end
    end

    # Prepend the reference configuration as the first sample row when requested.
    if cfg.include_reference
        # Force row 1 to the nominal reference configuration.
        ref = get_reference_config(pb.fuel_cell_type)
        ref_vec = Float64[]
        sizehint!(ref_vec, n_params)
        for b in pb.bounds
            v = getfield(ref, b.name)
            push!(ref_vec, Float64(v))
        end
        # Ensure integer parameters are integral in the reference row.
        for (j, b) in enumerate(pb.bounds)
            if b.type == :int
                ref_vec[j] = clamp(round(ref_vec[j]), b.min, b.max)
            end
        end
        @inbounds X[1, :] .= ref_vec
    end

    return X
end


"""
    apply_bounds_to_params(sample, pb, base_params) -> PhysicalParams

Return a new `PhysicalParams` with undetermined parameters replaced by `sample`.

All parameters not listed in `pb` (e.g. `Aact`, `nb_cell`, channel geometry) are
copied unchanged from `base_params`.

# Arguments
- `sample::Vector{Float64}`: One row from `generate_lhs_samples`.
- `pb::ParameterBounds`: Column-to-name mapping.
- `base_params::PhysicalParams`: Reference struct supplying all fixed parameters.
"""
function apply_bounds_to_params(sample::Vector{Float64},
                                pb::ParameterBounds,
                                base_params)::Any    # ::PhysicalParams
    length(sample) == pb.n_params ||
        throw(ArgumentError("Sample length $(length(sample)) does not match pb.n_params=$(pb.n_params)"))

    base_params isa PhysicalParams ||
        throw(ArgumentError("base_params must be a PhysicalParams (got $(typeof(base_params)))"))

    # Copy all fields, then override only sampled parameters.
    all_fields = fieldnames(PhysicalParams)
    base_nt = (; (f => getfield(base_params, f) for f in all_fields)...)

    overrides = Dict{Symbol, Any}()
    for (j, b) in enumerate(pb.bounds)
        # Convert / clamp to enforce type and bounds.
        if b.type == :int
            v = Int(clamp(round(sample[j]), b.min, b.max))
        else
            v = Float64(clamp(sample[j], b.min, b.max))
        end
        overrides[b.name] = v
    end

    # EH-31 constraint: Hccl mirrors Hacl when Hccl is not sampled.
    if pb.fuel_cell_type in (:EH_31_1_5, :EH_31_2_0, :EH_31_2_25, :EH_31_2_5)
        if haskey(overrides, :Hacl) && !haskey(overrides, :Hccl)
            overrides[:Hccl] = overrides[:Hacl]
        end
    end

    new_nt = merge(base_nt, (; (k => overrides[k] for k in keys(overrides))...))
    return PhysicalParams(; new_nt...)
end


# ─────────────────────────────────────────────────────────────────────────────
# REFERENCE CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

"""
    get_reference_config(fuel_cell_type) -> PhysicalParams

Return the nominal `PhysicalParams` for `fuel_cell_type`.

This configuration:
- Serves as the starting point for calibration.
- Is passed to the PRIM analysis as the *x_interest* point — the known-valid
  configuration that the found hyperbox must contain.
"""
function get_reference_config(fuel_cell_type::Symbol)::Any   # ::PhysicalParams
    # Physical parameters do not depend on the voltage zone in the current AlphaPEM
    # implementation, but `create_fuelcell` requires it.
    fc = create_fuelcell(fuel_cell_type, :full)
    return getfield(fc, :physical_parameters)
end

end # module ConfigurationSampling
