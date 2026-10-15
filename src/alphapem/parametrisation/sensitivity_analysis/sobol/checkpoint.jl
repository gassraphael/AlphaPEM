# -*- coding: utf-8 -*-

"""
    SOBOL_CHECKPOINT_FILENAME

Name of the checkpoint file written in `SobolAnalysisConfig.output_dir`.
"""
const SOBOL_CHECKPOINT_FILENAME = "sobol_checkpoint.jld2"


"""
    _sobol_checkpoint_path(output_dir)

Return the full path to the Sobol checkpoint file.
"""
_sobol_checkpoint_path(output_dir::String)::String = joinpath(output_dir, SOBOL_CHECKPOINT_FILENAME)


"""
    _config_hash_for_checkpoint(cfg)

Return a stable SHA-256 digest of the configuration fields that affect the
Sobol design and the AlphaPEM simulations.

Fields that only affect post-processing (`output_dir`, `top_k`, `save_curves`,
`resume`, `checkpoint_frequency`, `parallel`) are intentionally excluded so that
adjusting the output or the checkpoint policy does not invalidate a partial run.
"""
function _config_hash_for_checkpoint(cfg::SobolAnalysisConfig)::String
    buf = IOBuffer()

    print(buf, cfg.fuel_cell_type)
    print(buf, cfg.year)
    print(buf, cfg.voltage_zone)
    print(buf, cfg.N)
    print(buf, cfg.second_order)
    print(buf, cfg.seed)
    print(buf, join(string.(cfg.region_thresholds), ","))
    print(buf, cfg.include_operating_conditions)
    print(buf, cfg.max_run_time_s)
    print(buf, cfg.knn_k)

    if cfg.parameter_bounds === nothing
        print(buf, "nothing")
    else
        print(buf, _struct_to_string(cfg.parameter_bounds))
    end

    for c in cfg.operating_condition_constraints
        print(buf, c.name)
        print(buf, c.target)
        print(buf, join(string.(c.sources), ","))
        print(buf, c.kind)
        print(buf, c.active)
    end

    print(buf, join(string.(cfg.excluded_operating_conditions), ","))
    print(buf, _struct_to_string(cfg.polarization_params))
    print(buf, cfg.nb_gc)
    print(buf, _struct_to_string(cfg.numerical_params))

    return bytes2hex(sha256(take!(buf)))
end


"""
    _save_checkpoint(output_dir, cfg, A, B, params, regions, df_curves, elapsed_time)

Atomically persist the partial state of a Sobol analysis.

The checkpoint is first written to a temporary file and then renamed so that an
interrupted write never leaves a corrupt checkpoint behind.
"""
function _save_checkpoint(output_dir::String,
                          cfg::SobolAnalysisConfig,
                          A::Matrix{Float64},
                          B::Matrix{Float64},
                          params::Vector{InputParameter},
                          regions::Vector{PolarizationRegion},
                          df_curves::DataFrame,
                          elapsed_time::Float64)
    _save_jld2_atomic(
        _sobol_checkpoint_path(output_dir), "Sobol";
        config_hash = _config_hash_for_checkpoint(cfg),
        A = A,
        B = B,
        params = params,
        regions = regions,
        df_curves = df_curves,
        elapsed_time = elapsed_time,
    )
    return nothing
end


"""
    _load_checkpoint(output_dir, cfg)

Load a Sobol checkpoint if it exists and matches the current configuration.

Returns a `NamedTuple` `(A, B, params, regions, df_curves, elapsed_time)` on
success, or `nothing` if no usable checkpoint is found. If a checkpoint exists
but is incompatible, a warning is emitted and `nothing` is returned.
"""
function _load_checkpoint(output_dir::String, cfg::SobolAnalysisConfig)
    try
        data = _load_jld2_checked(
            _sobol_checkpoint_path(output_dir), _config_hash_for_checkpoint(cfg), "Sobol"
        )
        data === nothing && return nothing

        A = data["A"]::Matrix{Float64}
        B = data["B"]::Matrix{Float64}
        params = data["params"]::Vector{InputParameter}
        regions = data["regions"]::Vector{PolarizationRegion}
        df_curves = data["df_curves"]::DataFrame
        elapsed_time = data["elapsed_time"]::Float64

        # Sanity checks
        n_evals_expected = size(_fuse_designs_bootstrap(
            A, B; second_order = cfg.second_order, nboot = _sobol_nboot(cfg.N)
        ), 2)
        if nrow(df_curves) > 0 && maximum(df_curves.sample_id) > n_evals_expected
            @warn "Sobol checkpoint contains sample_ids beyond the expected design size. Ignoring checkpoint."
            return nothing
        end

        @info "Resumed Sobol analysis from checkpoint: $(nrow(df_curves)) sample(s) already completed."
        return (; A, B, params, regions, df_curves, elapsed_time)
    catch e
        @warn "Failed to load Sobol checkpoint: $e. Starting fresh."
        return nothing
    end
end
