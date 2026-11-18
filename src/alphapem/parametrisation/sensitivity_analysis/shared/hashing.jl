# -*- coding: utf-8 -*-

"""
    _struct_to_string(x)

Serialize a composite value to a deterministic string representation.

Only scalar fields are serialized; this is sufficient for the configuration objects used by the
sensitivity analyses (`PolarizationParams`, `NumericalParams`, `ParameterBounds`, and
`OperatingConditionConstraint`).

Parameters
----------
x : Any
    Composite value whose fields are serialized.

Returns
-------
String
    Deterministic representation of the type and its field values.
"""
function _struct_to_string(x)::String
    # Initialize the serialized type representation
    buf = IOBuffer()
    print(buf, typeof(x), "(")

    # Append each field using a stable scalar representation
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

    # Close and return the serialized value
    print(buf, ")")
    return String(take!(buf))
end


"""
    _hash_sample(sample)

Compute a SHA-256 digest of a sample vector.

Parameters
----------
sample : AbstractVector{<:Real}
    Sample values to identify.

Returns
-------
String
    Lowercase hexadecimal SHA-256 digest of the sample representation.
"""
function _hash_sample(sample::AbstractVector{<:Real})::String
    bytes = Vector{UInt8}(join(string.(sample), ","))
    return bytes2hex(sha256(bytes))
end
