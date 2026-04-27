# -*- coding: utf-8 -*-

"""Helpers for the GC current-distribution nonlinear solve scaling.

The nonlinear unknown vector layout is:
    [U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc]]

Helpers here keep scaling definitions out of `core/models` while preserving a
physical API for `calculate_1D_GC_current_density`.
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

"""Fill algebraic residuals for GC current distribution in physical units.

Unknown layout is fixed as:
    [U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc]]
"""
function gc_current_distribution_residuals!(res::AbstractVector,
                                            U_cell::Real,
                                            i_fc::AbstractVector{<:Real},
                                            C_O2_Pt::AbstractVector{<:Real},
                                            i_fc_cell::Real,
                                            sv::AbstractVector{<:CellState1D},
                                            fc::AbstractFuelCell,
                                            cfg::SimulationConfig)
    nb_gc = cfg.numerical_parameters.nb_gc
    length(sv) == nb_gc || throw(ArgumentError("sv size mismatch with cfg.numerical_parameters.nb_gc in gc_current_distribution_residuals!."))
    length(i_fc) == nb_gc || throw(ArgumentError("i_fc size mismatch in gc_current_distribution_residuals!."))
    length(C_O2_Pt) == nb_gc || throw(ArgumentError("C_O2_Pt size mismatch in gc_current_distribution_residuals!."))
    length(res) == 2 * nb_gc + 1 || throw(ArgumentError("res size mismatch in gc_current_distribution_residuals!."))

    @inbounds for i in 1:nb_gc
        res[i] = calculate_cell_voltage(i_fc[i], C_O2_Pt[i], sv[i], fc) - U_cell
    end
    res[nb_gc + 1] = i_fc_cell - average(i_fc)
    @inbounds for i in 1:nb_gc
        res[nb_gc + 1 + i] = calculate_C_O2_Pt(i_fc[i], sv[i], fc) - C_O2_Pt[i]
    end
    return nothing
end

