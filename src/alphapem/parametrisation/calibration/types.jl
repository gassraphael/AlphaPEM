"""
    ParameterBound

Bounds and metadata for a single parameter. `latex` is the display symbol (empty if unknown).
"""
struct ParameterBound
    name::Symbol
    min::Float64
    max::Float64
    type::Symbol          # :real or :int
    scale::Symbol         # :linear or :log
    unit::String
    description::String
    latex::String
end

# Backward-compatible constructors: no explicit scale means linear, no symbol means empty.
ParameterBound(name, min, max, type, scale::Symbol, unit, description) =
    ParameterBound(name, min, max, type, scale, unit, description, "")
ParameterBound(name, min, max, type, unit::AbstractString, description) =
    ParameterBound(name, min, max, type, :linear, unit, description, "")

"""
    ParameterBounds

Complete set of bounds for all undetermined parameters of a given fuel-cell type.
"""
struct ParameterBounds
    bounds::Vector{ParameterBound}
    fuel_cell_type::Symbol
    year::Union{Int,Nothing}
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
