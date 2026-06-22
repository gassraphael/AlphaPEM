"""
    new_PhysicalParams_from_sample(sample, pb, base_params) -> PhysicalParams

Return a new `PhysicalParams` object with undetermined parameters replaced by `sample`.
"""
function new_PhysicalParams_from_sample(sample::Vector{Float64},
                                pb::ParameterBounds,
                                base_params)::Any
    length(sample) == pb.n_params || throw(ArgumentError("Sample length mismatch"))
    base_params isa PhysicalParams || throw(ArgumentError("base_params must be PhysicalParams"))

    all_fields = fieldnames(PhysicalParams)
    base_nt = (; (f => getfield(base_params, f) for f in all_fields)...)

    overrides = Dict{Symbol, Any}()
    for (j, b) in enumerate(pb.bounds)
        if b.type == :int
            v = Int(clamp(round(sample[j]), b.min, b.max))
        else
            v = Float64(clamp(sample[j], b.min, b.max))
        end
        overrides[b.name] = v
    end

    if pb.fuel_cell_type in (:EH_31_1_5, :EH_31_2_0, :EH_31_2_25, :EH_31_2_5)
        if haskey(overrides, :Hacl) && !haskey(overrides, :Hccl)
            overrides[:Hccl] = overrides[:Hacl]
        end
    end

    new_nt = merge(base_nt, (; (k => overrides[k] for k in keys(overrides))...))
    return PhysicalParams(; new_nt...)
end

"""
    get_reference_config(fuel_cell_type) -> PhysicalParams

Return the nominal `PhysicalParams` for `fuel_cell_type`.
"""
function get_reference_config(fuel_cell_type::Symbol)::Any
    fc = create_fuelcell(fuel_cell_type, :full)
    return getfield(fc, :physical_parameters)
end
