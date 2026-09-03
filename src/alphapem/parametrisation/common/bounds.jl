"""
    bounds_for_fuel_cell(fuel_cell_type, voltage_zone = :full; year = nothing, nb_gc = 5) -> ParameterBounds

Return the undetermined-parameter bounds for a given fuel-cell type, voltage zone,
and GC spatial resolution (`nb_gc`).
"""
function bounds_for_fuel_cell(fuel_cell_type::Symbol,
                               voltage_zone::Symbol = :full;
                               year::Union{Int,Nothing} = nothing,
                               nb_gc::Int = 5)::ParameterBounds
    voltage_zone in (:full, :before_voltage_drop) ||
        throw(ArgumentError("voltage_zone must be :full or :before_voltage_drop (got $voltage_zone)"))

    bounds = ParameterBound[]
    fc = create_fuelcell(fuel_cell_type, voltage_zone; year=year, nb_gc=nb_gc)

    undetermined_params = undetermined_parameter_bounds(fc, voltage_zone)

    for (param_name, min_val, max_val) in undetermined_params
        if !haskey(PARAMETER_METADATA, param_name)
            throw(ArgumentError("No metadata found for parameter $param_name"))
        end
        unit, description = PARAMETER_METADATA[param_name]
        param_type = param_name == :e ? :int : :real
        push!(bounds, ParameterBound(param_name, Float64(min_val), Float64(max_val),
                                     param_type, unit, description))
    end

    for b in bounds
        b.min <= b.max || throw(ArgumentError("Invalid bounds for $(b.name): min=$(b.min) > max=$(b.max)"))
    end

    return ParameterBounds(bounds, fuel_cell_type, year, voltage_zone, length(bounds))
end
