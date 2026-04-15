# outputs_accessors.jl
#
# Typed accessors shared by display and post-processing functions.
# They centralise navigation through `SimulationOutputs` so that plotting code
# can work directly with specialised runtime structures.

"""Return the full simulation time history."""
time_history(outputs::SimulationOutputs) = outputs.solver.t


"""Return the typed solver-state history."""
solver_state_history(outputs::SimulationOutputs) = outputs.solver.states


"""Return the bundle of derived post-processing outputs."""
derived_outputs(outputs::SimulationOutputs) = outputs.derived


"""Collect a derived time series stored once per time step."""
extract_derived_series(outputs::SimulationOutputs,
                       accessor::Function) = accessor(derived_outputs(outputs))


"""Return the number of gas-channel nodes carried by the typed outputs."""
gas_channel_count(::SimulationOutputs{nb_gdl, nb_mpl, nb_gc}) where {nb_gdl, nb_mpl, nb_gc} = nb_gc


"""Return the middle gas-channel index used by legacy display routines."""
middle_gas_channel_index(fc::AbstractFuelCell) = cld(fc.numerical_parameters.nb_gc, 2)


"""Return the middle GDL index used by legacy display routines."""
middle_gdl_index(fc::AbstractFuelCell) = cld(fc.numerical_parameters.nb_gdl, 2)


"""Return the middle MPL index used by legacy display routines."""
middle_mpl_index(fc::AbstractFuelCell) = cld(fc.numerical_parameters.nb_mpl, 2)


"""Return the display start time when a current profile exposes one.

When no dedicated initialisation time exists, the full history is kept.
"""
display_start_time(cd::AbstractCurrent) = hasproperty(cd, :delta_t_ini) ? getproperty(cd, :delta_t_ini) : 0.0


"""Build the temporal mask used by display functions.

The historical display logic keeps only the tail of fixed simulations and all
points in dynamic mode.
"""
function display_time_mask(outputs::SimulationOutputs,
                           cd::AbstractCurrent,
                           cfg::SimulationConfig)::BitVector
    t_hist = time_history(outputs)
    if cfg.type_plot == :fixed
        return BitVector(t_hist .>= 0.9 * display_start_time(cd))
    end
    return trues(length(t_hist))
end


"""Return the masked simulation time history used for display."""
masked_time_history(outputs::SimulationOutputs,
                    cd::AbstractCurrent,
                    cfg::SimulationConfig) = time_history(outputs)[display_time_mask(outputs, cd, cfg)]


"""Collect and mask a derived time series stored once per time step."""
function extract_masked_derived_series(outputs::SimulationOutputs,
                                       cd::AbstractCurrent,
                                       cfg::SimulationConfig,
                                       accessor::Function)
    mask = display_time_mask(outputs, cd, cfg)
    return collect(extract_derived_series(outputs, accessor))[mask]
end


"""Return one MEA state at a given gas-channel position."""
mea_state_at(state::FuelCellStateP2D, gc_index::Integer) = state.nodes[gc_index]


"""Return the MEA state located at the middle of the gas channel."""
middle_mea_state(state::FuelCellStateP2D, fc::AbstractFuelCell) = mea_state_at(state, middle_gas_channel_index(fc))


"""Collect a time series from the typed MEA-state history for one GC position."""
function extract_mea_series(outputs::SimulationOutputs,
                            gc_index::Integer,
                            extractor::Function)
    return [extractor(mea_state_at(state, gc_index)) for state in solver_state_history(outputs)]
end


"""Collect a time series from the typed MEA-state history at the middle GC position."""
extract_mid_mea_series(outputs::SimulationOutputs,
                       fc::AbstractFuelCell,
                       extractor::Function) = extract_mea_series(outputs, middle_gas_channel_index(fc), extractor)


"""Collect and mask a time series from the typed MEA-state history for one GC position."""
function extract_masked_mea_series(outputs::SimulationOutputs,
                                   gc_index::Integer,
                                   cd::AbstractCurrent,
                                   cfg::SimulationConfig,
                                   extractor::Function)
    mask = display_time_mask(outputs, cd, cfg)
    return collect(extract_mea_series(outputs, gc_index, extractor))[mask]
end


"""Collect and mask a time series from the typed MEA-state history at the middle GC position."""
extract_masked_mid_mea_series(outputs::SimulationOutputs,
                              fc::AbstractFuelCell,
                              cd::AbstractCurrent,
                              cfg::SimulationConfig,
                              extractor::Function) =
    extract_masked_mea_series(outputs, middle_gas_channel_index(fc), cd, cfg, extractor)


"""Collect a derived time series stored per gas-channel node."""
function extract_derived_gc_series(outputs::SimulationOutputs,
                                   gc_index::Integer,
                                   accessor::Function)
    return accessor(derived_outputs(outputs))[gc_index]
end


"""Collect and mask a derived time series stored per gas-channel node."""
function extract_masked_derived_gc_series(outputs::SimulationOutputs,
                                          gc_index::Integer,
                                          cd::AbstractCurrent,
                                          cfg::SimulationConfig,
                                          accessor::Function)
    mask = display_time_mask(outputs, cd, cfg)
    return collect(extract_derived_gc_series(outputs, gc_index, accessor))[mask]
end


"""Collect a derived time series stored per gas-channel node at the middle GC position."""
extract_mid_derived_gc_series(outputs::SimulationOutputs,
                              fc::AbstractFuelCell,
                              accessor::Function) =
    extract_derived_gc_series(outputs, middle_gas_channel_index(fc), accessor)


"""Collect and mask a derived time series stored per gas-channel node at the middle GC position."""
extract_masked_mid_derived_gc_series(outputs::SimulationOutputs,
                                     fc::AbstractFuelCell,
                                     cd::AbstractCurrent,
                                     cfg::SimulationConfig,
                                     accessor::Function) =
    extract_masked_derived_gc_series(outputs, middle_gas_channel_index(fc), cd, cfg, accessor)


"""Return the nearest indices in `t_hist` associated with sampling times."""
function nearest_time_indices(t_hist::AbstractVector{<:Real},
                              sample_times::AbstractVector{<:Real})::Vector{Int}
    return [argmin(abs.(t_hist .- t_sample)) for t_sample in sample_times]
end


"""Return the characteristic sampling indices used for polarisation-like displays."""
polarisation_sampling_indices(outputs::SimulationOutputs,
                              cd::AbstractCurrent) =
    nearest_time_indices(time_history(outputs), polarisation_sampling_times(cd))


"""Return the characteristic sampling times used for polarisation-like displays."""
function polarisation_sampling_times(cd::AbstractCurrent)::Vector{Float64}
    hasproperty(cd, :delta_t_ini) ||
        throw(ArgumentError("The current profile does not expose `delta_t_ini` for polarisation sampling."))
    hasproperty(cd, :delta_t_break) ||
        throw(ArgumentError("The current profile does not expose `delta_t_break` for polarisation sampling."))

    delta_t_ini = getproperty(cd, :delta_t_ini)
    delta_t_break = getproperty(cd, :delta_t_break)

    if hasproperty(cd, :delta_i) && hasproperty(cd, :i_max) && hasproperty(cd, :v_load)
        delta_i = getproperty(cd, :delta_i)
        delta_t_load = delta_i / getproperty(cd, :v_load)
        nb_loads = floor(Int, getproperty(cd, :i_max) / delta_i)
        return [delta_t_ini + i * (delta_t_load + delta_t_break) for i in 1:nb_loads]
    elseif hasproperty(cd, :i_exp) && hasproperty(cd, :v_load)
        i_exp = getproperty(cd, :i_exp)
        delta_t_load = abs(i_exp[1]) / getproperty(cd, :v_load)
        delta_t_cali = delta_t_load + delta_t_break
        return [delta_t_ini + i * delta_t_cali for i in 1:length(i_exp)]
    end

    throw(ArgumentError("Unsupported current profile for polarisation sampling."))
end


"""Build the final pseudo-2D temperature matrix in °C.

Rows follow the gas-channel direction, while columns follow the through-plane
MEA ordering used in the display layer.
"""
function final_temperature_matrix_celsius(outputs::SimulationOutputs{nb_gdl, nb_mpl, nb_gc}) where {nb_gdl, nb_mpl, nb_gc}
    last_state = solver_state_history(outputs)[end]
    n_cols = 2 * nb_gdl + 2 * nb_mpl + 5
    temp_matrix = Matrix{Float64}(undef, nb_gc, n_cols)

    for k in 1:nb_gc
        mea = mea_state_at(last_state, k)
        temp_matrix[k, :] = vcat(
            mea.agc.T,
            [node.T for node in mea.agdl],
            [node.T for node in mea.ampl],
            mea.acl.T,
            mea.mem.T,
            mea.ccl.T,
            [node.T for node in mea.cmpl],
            [node.T for node in mea.cgdl],
            mea.cgc.T,
        ) .- 273.15
    end

    return temp_matrix
end


