"""
    generate_lhs_samples(pb, cfg = SamplingConfig()) -> Matrix{Float64}

Draw a sample matrix using the strategy specified by `cfg`.
"""
function generate_lhs_samples(pb::ParameterBounds,
                               cfg::SamplingConfig = SamplingConfig())::Matrix{Float64}
    n_samples = cfg.n_samples
    n_params = pb.n_params
    n_samples > 0 || throw(ArgumentError("n_samples must be > 0 (got $n_samples)"))

    rng = cfg.seed === nothing ? Random.default_rng() : MersenneTwister(cfg.seed)

    X = if cfg.method == :lhs
        raw = randomLHC(rng, n_samples, n_params)
        scaleLHC(raw, [(b.min, b.max) for b in pb.bounds])
    elseif cfg.method == :random
        raw = rand(rng, n_samples, n_params)
        for (j, b) in enumerate(pb.bounds)
            span = b.max - b.min
            @inbounds for i in 1:n_samples
                raw[i, j] = b.min + raw[i, j] * span
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
        ref = get_reference_config(pb.fuel_cell_type)
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
