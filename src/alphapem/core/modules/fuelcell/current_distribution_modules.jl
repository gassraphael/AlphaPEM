# -*- coding: utf-8 -*-

"""Helpers for the GC current-distribution algebraic block.

The algebraic unknown layout is:
    [U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc]]
"""

"""Build unknown and residual scaling vectors for the GC current solve."""
function _build_gc_current_density_scaling(cfg::SimulationConfig)::Tuple{Vector{Float64}, Vector{Float64}}
    nb_gc = cfg.numerical_parameters.nb_gc
    cd_scaling = StateScaling().current_distribution

    x_scales = Vector{Float64}(undef, 2 * nb_gc + 1)
    x_scales[1] = cd_scaling.U
    x_scales[2:nb_gc+1] .= cd_scaling.i_fc
    x_scales[nb_gc+2:2*nb_gc+1] .= cd_scaling.C_O2_Pt

    # Residual blocks follow the same physical-family scaling.
    res_scales = Vector{Float64}(undef, 2 * nb_gc + 1)
    res_scales[1:nb_gc] .= cd_scaling.U
    res_scales[nb_gc+1] = cd_scaling.i_fc
    res_scales[nb_gc+2:2*nb_gc+1] .= cd_scaling.C_O2_Pt

    return x_scales, res_scales
end


