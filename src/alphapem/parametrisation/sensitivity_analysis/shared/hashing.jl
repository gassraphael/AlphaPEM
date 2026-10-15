# -*- coding: utf-8 -*-

"""
    _struct_to_string(x)

Return a deterministic string representation of a composite type for hashing.

Only scalar fields are serialised; this is sufficient for the configuration
objects used by the sensitivity analyses (`PolarizationParams`, `NumericalParams`,
`ParameterBounds`, `OperatingConditionConstraint`).
"""
function _struct_to_string(x)::String
    buf = IOBuffer()
    print(buf, typeof(x), "(")
    for (i, fn) in enumerate(fieldnames(typeof(x)))
        i > 1 && print(buf, ",")
        val = getfield(x, fn)
        if isa(val, AbstractString) || isa(val, Symbol)
            print(buf, val)
        elseif isa(val, Number) || isa(val, Bool)
            print(buf, Float64(val))
        elseif isa(val, Tuple)
            print(buf, join(string.(val), ","))
        elseif isa(val, AbstractVector)
            print(buf, join(string.(val), ","))
        else
            print(buf, string(val))
        end
    end
    print(buf, ")")
    return String(take!(buf))
end


"""
    _hash_sample(sample)

Return a SHA-256 hex digest of a sample vector.

This provides a reproducible, content-based identifier for each sample of a
design, enabling traceability and deduplication across plan, simulations and outputs.
"""
function _hash_sample(sample::AbstractVector{<:Real})::String
    bytes = Vector{UInt8}(join(string.(sample), ","))
    return bytes2hex(sha256(bytes))
end
