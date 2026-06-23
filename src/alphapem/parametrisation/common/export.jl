"""
    export_parameter_bounds(bounds, filepath; method = :PRIM, metadata = Dict())::Nothing

Write parameter bounds to a YAML file.
"""
function export_parameter_bounds(bounds::Dict{Symbol, Tuple{Float64, Float64}},
                                 filepath::String;
                                 method::Symbol = :PRIM,
                                 metadata::Dict = Dict())::Nothing
    mkpath(dirname(filepath))
    ts = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    meta_block = merge(
        Dict("method" => string(method), "timestamp" => ts),
        Dict(string(k) => v for (k, v) in metadata)
    )

     open(filepath, "w") do io
         println(io, "metadata:")
         for k in sort(collect(keys(meta_block)))
             v = meta_block[k]
             println(io, "  $(k): $(v)")
         end
        println(io, "parameters:")
        for name in sort(collect(keys(bounds)))
            lo, hi = bounds[Symbol(name)]
            println(io, "  $(name):")
            println(io, "    min: $(lo)")
            println(io, "    max: $(hi)")
        end
    end
    return nothing
end

"""
    export_calibrated_params(params, filepath; method = :calibration, metadata = Dict())::Nothing

Write calibrated parameter values (single value per parameter) to a YAML file.
"""
function export_calibrated_params(params::Dict{Symbol, Float64},
                                  filepath::String;
                                  method::Symbol = :calibration,
                                  metadata::Dict = Dict())::Nothing
    mkpath(dirname(filepath))
    ts = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")
    meta_block = merge(
        Dict("method" => string(method), "timestamp" => ts),
        Dict(string(k) => v for (k, v) in metadata)
    )

    open(filepath, "w") do io
        println(io, "metadata:")
        for k in sort(collect(keys(meta_block)))
            v = meta_block[k]
            println(io, "  $(k): $(v)")
        end
        println(io, "parameters:")
        for name in sort(collect(keys(params)))
            println(io, "  $(name): $(params[name])")
        end
    end
    return nothing
end
