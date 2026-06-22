"""
    ParameterBound

Bounds and metadata for a single undetermined parameter.
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
"""
Base.@kwdef struct SamplingConfig
    n_samples::Int            = 10000
    method::Symbol            = :lhs
    seed::Union{Int, Nothing} = 42
    include_reference::Bool   = true
end
