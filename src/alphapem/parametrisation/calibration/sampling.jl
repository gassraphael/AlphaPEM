"""
    _scale_unit_interval(u, b)

Map a normalized value `u ∈ [0, 1]` to the physical range of `b`, respecting
`b.scale` (`:linear` or `:log`).  For `:log` bounds the mapping is
`exp(log(min) + u * (log(max) - log(min)))` so that each decade is equally
sampled.
"""
function _scale_unit_interval(u::Real, b::ParameterBound)::Float64
    if b.scale == :log
        return clamp(exp(log(b.min) + u * (log(b.max) - log(b.min))), b.min, b.max)
    else
        return clamp(b.min + u * (b.max - b.min), b.min, b.max)
    end
end

"""
    generate_lhs_samples(pb, cfg = SamplingConfig()) -> Matrix{Float64}

Draw a sample matrix using the strategy specified by `cfg`.  Parameters whose
bounds span at least two orders of magnitude are sampled in log-space.
"""
function generate_lhs_samples(pb::ParameterBounds,
                               cfg::SamplingConfig = SamplingConfig())::Matrix{Float64}
    n_samples = cfg.n_samples
    n_params = pb.n_params
    n_samples > 0 || throw(ArgumentError("n_samples must be > 0 (got $n_samples)"))

    rng = cfg.seed === nothing ? Random.default_rng() : MersenneTwister(cfg.seed)

    X = if cfg.method == :lhs
        raw = randomLHC(rng, n_samples, n_params)
        # First scale linearly to the declared bounds; this preserves the LHS
        # stratification because the subsequent log transform is monotonic.
        X_lin = scaleLHC(raw, [(b.min, b.max) for b in pb.bounds])
        X_lin = Matrix{Float64}(X_lin)
        for (j, b) in enumerate(pb.bounds)
            b.scale == :log || continue
            span = b.max - b.min
            log_min = log(b.min)
            log_max = log(b.max)
            @inbounds for i in 1:n_samples
                u = (X_lin[i, j] - b.min) / span
                X_lin[i, j] = clamp(exp(log_min + u * (log_max - log_min)), b.min, b.max)
            end
        end
        X_lin
    elseif cfg.method == :random
        raw = rand(rng, n_samples, n_params)
        for (j, b) in enumerate(pb.bounds)
            @inbounds for i in 1:n_samples
                raw[i, j] = _scale_unit_interval(raw[i, j], b)
            end
        end
        raw
    else
        throw(ArgumentError("Unsupported sampling method: $(cfg.method)"))
    end
    X = Matrix{Float64}(X)

    for (j, b) in enumerate(pb.bounds)
        if b.type == :int
            @inbounds for i in 1:n_samples
                X[i, j] = clamp(round(X[i, j]), b.min, b.max)
            end
        end
    end

    if cfg.include_reference
        ref = get_reference_config(pb.fuel_cell_type; year=pb.year)
        ref_vec = Float64[]
        for b in pb.bounds
            v = getfield(ref, b.name)
            push!(ref_vec, Float64(v))
        end
        for (j, b) in enumerate(pb.bounds)
            if b.type == :int
                ref_vec[j] = clamp(round(ref_vec[j]), b.min, b.max)
            end
        end
        @inbounds X[1, :] .= ref_vec
    end

    return X
end
