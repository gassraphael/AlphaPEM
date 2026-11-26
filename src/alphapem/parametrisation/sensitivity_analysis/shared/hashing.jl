# -*- coding: utf-8 -*-

"""
    _struct_to_string(x)

Serialize a value and its nested fields/collections to a deterministic string representation.

Parameters
----------
x : Any
    Composite value whose fields are serialized.

Returns
-------
String
    Deterministic representation of the type, dimensions, elements, and field values.
"""
function _struct_to_string(x)::String
    # Recursively serialize collections instead of relying on display strings that may include
    # process-specific memory addresses on newer Julia versions.
    buf = IOBuffer()
    _write_stable_value(buf, x)
    return String(take!(buf))
end


"""
    _write_stable_value(io, value)

Recursively write stable representations of scalar, collection, and struct values.

Parameters
----------
io : IO
    Destination buffer.
value : Any
    Value to serialize.

Returns
-------
Nothing
    Writes the representation to `io`.
"""
function _write_stable_value(io::IO, value)
    if value === nothing
        print(io, "nothing")
    elseif value isa AbstractString
        print(io, "String(", repr(String(value)), ")")
    elseif value isa Symbol
        print(io, "Symbol(", repr(String(value)), ")")
    elseif value isa Number || value isa Bool
        print(io, typeof(value), ":", repr(value))
    elseif value isa Tuple
        print(io, typeof(value), "(")
        for (index, item) in enumerate(value)
            index > 1 && print(io, ",")
            _write_stable_value(io, item)
        end
        print(io, ")")
    elseif value isa AbstractDict
        # Sort by serialized key: dictionary iteration order and internal storage are not stable
        entries = sort!([(_struct_to_string(k), v) for (k, v) in value]; by = first)
        print(io, typeof(value), "{")
        for (index, (key, item)) in enumerate(entries)
            index > 1 && print(io, ",")
            print(io, key, "=>")
            _write_stable_value(io, item)
        end
        print(io, "}")
    elseif value isa AbstractArray
        print(io, typeof(value), "size=", size(value), "[")
        for (index, item) in enumerate(value)
            index > 1 && print(io, ",")
            _write_stable_value(io, item)
        end
        print(io, "]")
    elseif isstructtype(typeof(value))
        print(io, typeof(value), "(")
        for (index, field) in enumerate(fieldnames(typeof(value)))
            index > 1 && print(io, ",")
            print(io, field, "=")
            _write_stable_value(io, getfield(value, field))
        end
        print(io, ")")
    else
        print(io, typeof(value), ":", repr(value))
    end
    return nothing
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
