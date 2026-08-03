# -*- coding: utf-8 -*-

"""
    InputParameter

Descriptor for one input parameter of the Sobol analysis.

# Fields
- `name::Symbol`: Parameter name.
- `min::Float64`: Lower bound.
- `max::Float64`: Upper bound.
- `type::Symbol`: `:real` or `:int`.
- `source::Symbol`: `:physical` or `:operating`.
"""
struct InputParameter
    name::Symbol
    min::Float64
    max::Float64
    type::Symbol
    source::Symbol
end


"""
    build_input_parameters(cfg)

Build the list of input parameters to vary in the Sobol analysis.

Includes undetermined physical parameters and, if requested, operating conditions.
Returns a vector of `InputParameter`.
"""
function build_input_parameters(cfg::SobolAnalysisConfig)::Vector{InputParameter}
    params = InputParameter[]

    # Physical undetermined parameters
    pb = cfg.parameter_bounds !== nothing ? cfg.parameter_bounds : bounds_for_fuel_cell(cfg.fuel_cell_type, cfg.voltage_zone)
    for b in pb.bounds
        push!(params, InputParameter(b.name, b.min, b.max, b.type, :physical))
    end

    if cfg.include_operating_conditions
        # Operating conditions from the reference fuel cell
        fc = create_fuelcell(cfg.fuel_cell_type, cfg.voltage_zone)
        oc = fc.operating_conditions
        for f in fieldnames(typeof(oc))
            v = getfield(oc, f)
            if v isa Real
                # Use ±20% around the nominal value as a sensible range
                val = Float64(v)
                delta = abs(val) * 0.2
                if delta == 0.0
                    delta = 1.0
                end
                lo = val - delta
                hi = val + delta
                push!(params, InputParameter(f, lo, hi, :real, :operating))
            end
        end
    end

    return params
end


"""
    parameter_ranges(params)

Return a vector of `(min, max)` tuples compatible with `GlobalSensitivity.gsa`.
"""
function parameter_ranges(params::Vector{InputParameter})::Vector{Tuple{Float64, Float64}}
    return [(p.min, p.max) for p in params]
end


"""
    generate_sobol_design_matrices(cfg, params)

Generate the Sobol design matrices required by `GlobalSensitivity.Sobol`.

Returns `(A, B)` where each column is one sample and rows correspond to parameters.
The matrices are scaled to the parameter ranges and integer parameters are rounded.
"""
function generate_sobol_design_matrices(cfg::SobolAnalysisConfig,
                                        params::Vector{InputParameter})
    ranges = parameter_ranges(params)
    lb = [r[1] for r in ranges]
    ub = [r[2] for r in ranges]

    sampler = QuasiMonteCarlo.SobolSample()
    R = QuasiMonteCarlo.Shift()

    # Generate A and B matrices: each column is a sample
    n_params = length(lb)
    A, B = QuasiMonteCarlo.generate_design_matrices(cfg.N, n_params, sampler, R, 2)

    # Scale to parameter ranges
    for j in 1:n_params
        span = ub[j] - lb[j]
        for i in 1:cfg.N
            A[j, i] = lb[j] + A[j, i] * span
            B[j, i] = lb[j] + B[j, i] * span
        end
    end

    # Enforce integer parameters
    A = _enforce_integer_params(A, params)
    B = _enforce_integer_params(B, params)

    return A, B
end


"""
    _enforce_integer_params(X, params)

Round and clamp integer parameters in sample matrix `X`.
"""
function _enforce_integer_params(X::Matrix{Float64},
                                  params::Vector{InputParameter})::Matrix{Float64}
    X_out = copy(X)
    for (j, p) in enumerate(params)
        if p.type == :int
            for i in 1:size(X_out, 2)
                X_out[j, i] = clamp(round(X_out[j, i]), p.min, p.max)
            end
        end
    end
    return X_out
end


"""
    sample_to_physical_params(sample, params, base_params)

Map a Sobol sample vector to a `PhysicalParams` object.

Only physical parameters are overridden; operating conditions are handled separately.
"""
function sample_to_physical_params(sample::Vector{Float64},
                                   params::Vector{InputParameter},
                                   base_params)::PhysicalParams
    physical_overrides = Dict{Symbol, Any}()
    for (j, p) in enumerate(params)
        if p.source == :physical
            if p.type == :int
                physical_overrides[p.name] = Int(clamp(round(sample[j]), p.min, p.max))
            else
                physical_overrides[p.name] = Float64(clamp(sample[j], p.min, p.max))
            end
        end
    end

    # Use ParametrisationCommon mapping for consistency (handles EH-31 constraints)
    pb = _params_to_bounds(params)
    if isempty(pb.bounds)
        return base_params
    end

    sample_physical = Float64[sample[j] for j in 1:length(params) if params[j].source == :physical]
    return new_PhysicalParams_from_sample(sample_physical, pb, base_params)
end


"""
    sample_to_operating_conditions(sample, params, base_oc)

Map a Sobol sample vector to an `OperatingConditions` object.
"""
function sample_to_operating_conditions(sample::Vector{Float64},
                                        params::Vector{InputParameter},
                                        base_oc)::Any
    overrides = Dict{Symbol, Any}()
    for (j, p) in enumerate(params)
        if p.source == :operating
            if p.type == :int
                overrides[p.name] = Int(clamp(round(sample[j]), p.min, p.max))
            else
                overrides[p.name] = Float64(clamp(sample[j], p.min, p.max))
            end
        end
    end

    if isempty(overrides)
        return base_oc
    end

    all_fields = fieldnames(typeof(base_oc))
    nt = (; (f => getfield(base_oc, f) for f in all_fields)...)
    new_nt = merge(nt, (; (k => overrides[k] for k in keys(overrides))...))
    return typeof(base_oc)(; new_nt...)
end


"""
    _params_to_bounds(params)

Convert a vector of `InputParameter` to a `ParameterBounds` containing only physical params.
"""
function _params_to_bounds(params::Vector{InputParameter})::ParameterBounds
    bounds = ParameterBound[]
    for p in params
        if p.source == :physical
            unit, description = get(PARAMETER_METADATA, p.name, ("", ""))
            push!(bounds, ParameterBound(p.name, p.min, p.max, p.type, unit, description))
        end
    end
    return ParameterBounds(bounds, :unknown, :full, length(bounds))
end
