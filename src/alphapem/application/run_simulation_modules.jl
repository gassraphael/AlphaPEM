# -*- coding: utf-8 -*-

"""
Utility helpers used by `run_simulation.jl`.

This file intentionally contains only secondary support functions:
- matplotlib/figure preparation
- internal-state extraction for dynamic segmented runs
"""

# ═══════════════════════════════════════════════════════════════════════════════
#  Matplotlib helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""Apply global matplotlib style settings shared by all simulations."""
function _setup_matplotlib_style!()
    mpl.rcParams["font.family"] = "cmr10" # "cmr10" for English characters and "DejaVu Serif" for French ones
    mpl.rcParams["axes.formatter.use_mathtext"] = true # For scientific notation
    mpl.rcParams["lines.linewidth"] = 2.0
    mpl.rcParams["lines.markersize"] = 5.0
end

"""Return figure/axes layout for step-current runs."""
function _create_figures(cfg::SimulationConfig{<:StepParams})
    if cfg.type_display == :multiple
        mpl.rcParams["font.size"] = 18 # Font size for all text
        # Here, additional plots are unnecessary: saving is handled in `Display` for this mode.
        return (nothing, nothing, nothing, nothing, nothing, nothing)
    elseif cfg.type_display == :synthetic
        mpl.rcParams["font.size"] = 13 # Font size for all text
        fig1, ax1 = plt.subplots(3, 3; figsize=(14, 14))
        fig2, ax2 = cfg.type_plot == :fixed ? plt.subplots(; figsize=(8, 8)) : (nothing, nothing)
        plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        return (fig1, ax1, fig2, ax2, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end

"""Return figure/axes layout for polarization runs."""
function _create_figures(cfg::SimulationConfig{<:PolarizationParams})
    if cfg.type_display == :multiple
        mpl.rcParams["font.size"] = 11 # Font size for all text
        fig1, ax1 = plt.subplots(1, 3; figsize=(14, 4.7))
        fig2, ax2 = plt.subplots(1, 4; figsize=(18.7, 4.7))
        plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        return (fig1, ax1, fig2, ax2, nothing, nothing)
    elseif cfg.type_display == :synthetic
        mpl.rcParams["font.size"] = 18 # Font size for all text
        mpl.rcParams["legend.fontsize"] = 15 # Legend font size only
        fig1, ax1 = plt.subplots(; figsize=(8, 8))
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end

"""Return figure/axes layout for calibration-polarization runs."""
function _create_figures(cfg::SimulationConfig{<:PolarizationCalibrationParams})
    if cfg.type_display == :multiple
        mpl.rcParams["font.size"] = 11 # Font size for all text
        fig1, ax1 = plt.subplots(1, 3; figsize=(14, 4.7))
        plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    elseif cfg.type_display == :synthetic
        mpl.rcParams["font.size"] = 18 # Font size for all text
        fig1, ax1 = plt.subplots(; figsize=(8, 8))
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end

"""Return figure/axes layout for EIS runs."""
function _create_figures(cfg::SimulationConfig{<:EISParams})
    if cfg.type_display == :multiple
        mpl.rcParams["font.size"] = 18 # Font size for all text
        fig1, ax1 = plt.subplots(; figsize=(8, 8))
        fig2, ax2 = plt.subplots(; figsize=(8, 8))
        fig3, ax3 = plt.subplots(; figsize=(8, 8))
        return (fig1, ax1, fig2, ax2, fig3, ax3)
    elseif cfg.type_display == :synthetic
        mpl.rcParams["font.size"] = 13 # Font size for all text
        fig1, ax1 = plt.subplots(1, 3; figsize=(14, 4.7))
        plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end

"""
Configure matplotlib and return the figure/axes layout required for `cfg`.
Tuple elements are `nothing` when unused.
"""
function figures_preparation(cfg::SimulationConfig)
    _setup_matplotlib_style!()
    # Enable interactive mode for non-blocking display.
    plt.ion()
    cfg.type_display == :no_display &&
        return (nothing, nothing, nothing, nothing, nothing, nothing)
    return _create_figures(cfg)
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Internal-state helpers
# ═══════════════════════════════════════════════════════════════════════════════

"""Append one typed 1-D MEA state to canonical solver ordering."""
function _append_mea_state_1D_to_solver_order!(dest::Vector{Float64}, state)
    # C_v block
    push!(dest, state.agc.C_v)
    for node in state.agdl; push!(dest, node.C_v); end
    for node in state.ampl; push!(dest, node.C_v); end
    push!(dest, state.acl.C_v, state.ccl.C_v)
    for node in state.cmpl; push!(dest, node.C_v); end
    for node in state.cgdl; push!(dest, node.C_v); end
    push!(dest, state.cgc.C_v)

    # s block
    push!(dest, state.agc.s)
    for node in state.agdl; push!(dest, node.s); end
    for node in state.ampl; push!(dest, node.s); end
    push!(dest, state.acl.s, state.ccl.s)
    for node in state.cmpl; push!(dest, node.s); end
    for node in state.cgdl; push!(dest, node.s); end
    push!(dest, state.cgc.s)

    # lambda block
    push!(dest, state.acl.lambda, state.mem.lambda, state.ccl.lambda)

    # species block
    push!(dest, state.agc.C_H2)
    for node in state.agdl; push!(dest, node.C_H2); end
    for node in state.ampl; push!(dest, node.C_H2); end
    push!(dest, state.acl.C_H2, state.ccl.C_O2)
    for node in state.cmpl; push!(dest, node.C_O2); end
    for node in state.cgdl; push!(dest, node.C_O2); end
    push!(dest, state.cgc.C_O2, state.agc.C_N2, state.cgc.C_N2)

    # temperature block
    push!(dest, state.agc.T)
    for node in state.agdl; push!(dest, node.T); end
    for node in state.ampl; push!(dest, node.T); end
    push!(dest, state.acl.T, state.mem.T, state.ccl.T)
    for node in state.cmpl; push!(dest, node.T); end
    for node in state.cgdl; push!(dest, node.T); end
    push!(dest, state.cgc.T)

    # voltage block
    push!(dest, state.ccl.eta_c)
    return dest
end

"""
Extract final internal states from a simulation to initialize the next segment.
"""
function _extract_last_internal_state(simu::AlphaPEM, cfg::SimulationConfig)::Vector{Float64}
    dest = Float64[]

    # Fast typed path: read directly from typed simulation outputs.
    last_state = simu.outputs.solver.states[end]

    for k in 1:simu.fuel_cell.numerical_parameters.nb_gc
        _append_mea_state_1D_to_solver_order!(dest, last_state.nodes[k])
    end

    # Auxiliary typed state extraction is intentionally deferred until subsystem reactivation.
    if cfg.type_auxiliary in (
        :forced_convective_cathode_with_flow_through_anode,
        :forced_convective_cathode_with_anodic_recirculation,
    )
        error("Auxiliary state extraction is not yet implemented in the typed API.")
    end

    return dest
end



