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
    pb = cfg.parameter_bounds !== nothing ? cfg.parameter_bounds : bounds_for_fuel_cell(cfg.fuel_cell_type, cfg.voltage_zone; year=cfg.year)
    for b in pb.bounds
        push!(params, InputParameter(b.name, b.min, b.max, b.type, :physical))
    end

    if cfg.include_operating_conditions
        # Operating conditions from OPERATING_CONDITIONS_BOUNDS
        excluded = Set(cfg.excluded_operating_conditions)
        for (name, (lo, hi, ptype)) in OPERATING_CONDITIONS_BOUNDS
            name in excluded && continue
            push!(params, InputParameter(name, lo, hi, ptype, :operating))
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

This implementation follows the Saltelli scheme used by SALib: a single scrambled
Sobol sequence of length `2N` is generated, then split into matrix `A` (first `N`
points) and matrix `B` (next `N` points). The random seed is taken from `cfg.seed`
to ensure reproducibility.
"""
function generate_sobol_design_matrices(cfg::SobolAnalysisConfig,
                                        params::Vector{InputParameter})
    ranges = parameter_ranges(params)
    lb = [r[1] for r in ranges]
    ub = [r[2] for r in ranges]
    n_params = length(lb)

    # Scrambled Sobol sequence with a reproducible RNG.
    rng = Random.MersenneTwister(cfg.seed)
    scramble = QuasiMonteCarlo.OwenScramble(; base = 2, rng = rng)
    sampler = QuasiMonteCarlo.SobolSample(scramble)

    # Generate 2N points and split into the Saltelli A/B matrices.
    all_points = QuasiMonteCarlo.sample(2 * cfg.N, n_params, sampler)
    A = all_points[:, 1:cfg.N]
    B = all_points[:, (cfg.N + 1):(2 * cfg.N)]

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
    # Use ParametrisationCommon mapping for consistency (handles EH constraints)
    pb = _params_to_bounds(params)
    if isempty(pb.bounds)
        return base_params
    end

    sample_physical = Float64[sample[j] for j in 1:length(params) if params[j].source == :physical]
    return new_PhysicalParams_from_sample(sample_physical, pb, base_params)
end


"""
    sample_to_operating_conditions(sample, params, base_oc, constraints, excluded)

Map a Sobol sample vector to an `OperatingConditions` object, applying constraints.

Operating conditions listed in `excluded` are not sampled; they keep their nominal
value from `base_oc`.

Constraints are applied according to their `kind`:
- `:(=)`: `target = fn(overrides)`.
- `:(>=)`: `target = max(target_sampled, fn(overrides))`, then clamped to the
  declared bounds of `target`.
- `:(<=)`: `target = min(target_sampled, fn(overrides))`, then clamped to the
  declared bounds of `target`.
"""
function sample_to_operating_conditions(sample::Vector{Float64},
                                        params::Vector{InputParameter},
                                        base_oc,
                                        constraints::Vector{OperatingConditionConstraint} = OperatingConditionConstraint[],
                                        excluded::Vector{Symbol} = Symbol[])::Any
    overrides = Dict{Symbol, Any}()

    # Sampled operating conditions
    for (j, p) in enumerate(params)
        if p.source == :operating
            if p.type == :int
                overrides[p.name] = Int(clamp(round(sample[j]), p.min, p.max))
            else
                overrides[p.name] = Float64(clamp(sample[j], p.min, p.max))
            end
        end
    end

    # Excluded operating conditions keep their nominal value
    for name in excluded
        if hasproperty(base_oc, name)
            overrides[name] = getfield(base_oc, name)
        end
    end

    # Apply active constraints
    for c in constraints
        if c.active && haskey(overrides, c.target)
            # Ensure all source values are present in overrides
            if all(haskey(overrides, s) for s in c.sources)
                target_bounds = get(OPERATING_CONDITIONS_BOUNDS, c.target, nothing)
                current = Float64(overrides[c.target])
                bound_value = Float64(c.fn(overrides))

                if c.kind == :(=)
                    current = bound_value
                elseif c.kind == :(>=)
                    current = max(current, bound_value)
                elseif c.kind == :(<=)
                    current = min(current, bound_value)
                end

                # Clamp to declared bounds if they exist
                if target_bounds !== nothing
                    lo, hi, _ = target_bounds
                    current = clamp(current, lo, hi)
                end

                overrides[c.target] = current
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
    return ParameterBounds(bounds, :unknown, nothing, :full, length(bounds))
end


"""
    is_valid_operating_conditions(oc, constraints=OperatingConditionConstraint[])

Check that operating conditions are physically consistent and that all active
operating-condition constraints are satisfied.

Current checks:
- Pressures are positive.
- Humidities are in [0, 1].
- Stoichiometric ratios are >= 1.
- Temperature is positive in °C.
- Active operating-condition constraints are respected (with tolerance).
"""
function is_valid_operating_conditions(oc,
                                       constraints::Vector{OperatingConditionConstraint} = OperatingConditionConstraint[])::Bool
    try
        # Pressures
        Pa_des = oc.Pa_des
        Pc_des = oc.Pc_des
        Pa_des >= 1.0 && Pc_des >= 1.0 || return false

        # Humidities
        0.0 <= oc.Phi_a_des <= 1.0 || return false
        0.0 <= oc.Phi_c_des <= 1.0 || return false

        # Stoichiometric ratios
        oc.Sa >= 1.0 || return false
        oc.Sc >= 1.0 || return false

        # Temperature
        oc.T_des > 273.15 || return false

        # Hydrogen molar fraction
        0.0 <= oc.y_H2_in <= 1.0 || return false

        # Active operating-condition constraints
        oc_dict = Dict{Symbol, Float64}()
        for f in fieldnames(typeof(oc))
            v = getfield(oc, f)
            if v isa Real
                oc_dict[f] = Float64(v)
            end
        end
        tol = 1e-9
        for c in constraints
            c.active || continue
            if all(haskey(oc_dict, s) for s in c.sources) && haskey(oc_dict, c.target)
                bound_value = Float64(c.fn(oc_dict))
                target_value = oc_dict[c.target]
                if c.kind == :(=)
                    abs(target_value - bound_value) <= tol || return false
                elseif c.kind == :(>=)
                    target_value + tol >= bound_value || return false
                elseif c.kind == :(<=)
                    target_value - tol <= bound_value || return false
                end
            end
        end

        return true
    catch
        return false
    end
end
