# -*- coding: utf-8 -*-

"""
Utility helpers used by `run_simulation.jl`.

This file intentionally contains only secondary support functions:
- native CairoMakie figure preparation
- internal-state extraction for dynamic segmented runs
"""

# ═══════════════════════════════════════════════════════════════════════════════
#  Figure-preparation helpers (Makie backends)
# ═══════════════════════════════════════════════════════════════════════════════

const _interactive_screens = Any[]
const _interactive_fig_to_screen = IdDict{Any, Any}()
const _active_display_backend = Ref{Symbol}(:cairo)
const _prepared_interactive_figures = Any[]

"""Return true when the simulation should be shown in a native interactive window."""
_use_interactive_display(cfg::SimulationConfig) =
    cfg.type_display != :no_display

"""Apply a publication-oriented Makie theme and activate the required backend."""
function _setup_makie_theme!(cfg::SimulationConfig)
    publication_palette = [
        "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
        "#E69F00", "#882255", "#44AA99", "#999933", "#332288", "#DDCC77",
    ]

    if _use_interactive_display(cfg)
        try
            GLMakie.activate!()
            _active_display_backend[] = :gl
        catch err
            @warn "GLMakie could not be activated; falling back to CairoMakie non-interactive display." exception=(err, catch_backtrace())
            CairoMakie.activate!()
            _active_display_backend[] = :cairo
        end
    else
        CairoMakie.activate!()
        _active_display_backend[] = :cairo
    end

    set_theme!(Theme(
        fontsize=15,
        linewidth=2.2,
        markersize=9,
        palette=(color=publication_palette,),
        Axis=(
            xgridvisible=true,
            ygridvisible=true,
            xminorticksvisible=true,
            yminorticksvisible=true,
            xminorgridvisible=false,
            yminorgridvisible=false,
        ),
    ))
    return nothing
end


"""Display figures in native GLMakie windows and enable pointer-based inspection."""
function _open_interactive_figures!(figs...)
    for fig in figs
        fig === nothing && continue
        fig isa AbstractVector && (foreach(f -> _open_interactive_figures!(f), fig); continue)

        if haskey(_interactive_fig_to_screen, fig)
            screen = _interactive_fig_to_screen[fig]
            if isopen(screen)
                continue
            else
                delete!(_interactive_fig_to_screen, fig)
            end
        end

        DataInspector(fig)
        screen = GLMakie.Screen()
        display(screen, fig)
        _interactive_fig_to_screen[fig] = screen
        screen in _interactive_screens || push!(_interactive_screens, screen)
    end
    return nothing
end


"""Store the latest prepared figure handles for deferred interactive display."""
function _store_prepared_figures!(figs...)
    empty!(_prepared_interactive_figures)
    for fig in figs
        fig === nothing && continue
        push!(_prepared_interactive_figures, fig)
    end
    return nothing
end


"""Create a figure with a regular 2D grid of axes."""
function _make_grid_axes(n_rows::Integer,
                         n_cols::Integer;
                         size::Tuple{Int, Int})
    fig = Figure(; size=size)
    axes = [Axis(fig[r, c]) for r in 1:n_rows, c in 1:n_cols]
    return fig, axes
end


"""Create a figure with one horizontal row of axes."""
function _make_row_axes(n_cols::Integer;
                        size::Tuple{Int, Int})
    fig = Figure(; size=size)
    axes = [Axis(fig[1, c]) for c in 1:n_cols]
    return fig, axes
end


"""Create a figure with a single axis."""
function _make_single_axis(; size::Tuple{Int, Int})
    fig = Figure(; size=size)
    ax = Axis(fig[1, 1])
    return fig, ax
end


"""Return figure/axes layout for step-current runs."""
function _create_figures(cfg::SimulationConfig{<:StepParams}, nb_gc::Union{Nothing, Integer}=nothing)
    has_extended_gc_profiles = isnothing(nb_gc) || nb_gc >= 3
    if cfg.type_display == :multiple
        # Multiple mode: one figure per internal-state plot.
        fig1 = Figure[]
        ax1 = Axis[]
        for _ in 1:14
            f, a = _make_single_axis(size=(760, 600))
            push!(fig1, f)
            push!(ax1, a)
        end
        fig2, ax2 = (cfg.display_timing == :postrun && has_extended_gc_profiles) ?
                    _make_single_axis(size=(740, 620)) : (nothing, nothing)
        if cfg.display_timing == :postrun && has_extended_gc_profiles
            fig3 = Figure[]
            ax3 = Axis[]
            for _ in 1:3
                f, a = _make_single_axis(size=(740, 620))
                push!(fig3, f)
                push!(ax3, a)
            end
        else
            fig3, ax3 = nothing, nothing
        end
        return (fig1, ax1, fig2, ax2, fig3, ax3)
    elseif cfg.type_display == :synthetic
        fig1, ax1 = _make_grid_axes(3, 3; size=(1360, 1120))
        if cfg.display_timing == :postrun && has_extended_gc_profiles
            fig2 = Figure(; size=(1260, 1120))
            # Keep the middle column free for the heatmap colorbar.
            ax2 = [Axis(fig2[1, 1]) Axis(fig2[1, 3]);
                   Axis(fig2[2, 1]) Axis(fig2[2, 3])]
            return (fig1, ax1, fig2, ax2, nothing, nothing)
        end
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end


"""Return figure/axes layout for polarization runs."""
function _create_figures(cfg::SimulationConfig{<:PolarizationParams})
    if cfg.type_display == :multiple
        # Multiple mode: internal states + derived curves in individual figures,
        # plus a dedicated polarization curve figure.
        fig1 = Figure[]
        ax1 = Axis[]
        for _ in 1:5
            f, a = _make_single_axis(size=(760, 600))
            push!(fig1, f)
            push!(ax1, a)
        end
        fig2, ax2 = _make_single_axis(size=(760, 600))
        return (fig1, ax1, fig2, ax2, nothing, nothing)
    elseif cfg.type_display == :synthetic
        fig1, ax1 = _make_single_axis(size=(760, 600))
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end


"""Return figure/axes layout for calibration-polarization runs."""
function _create_figures(cfg::SimulationConfig{<:PolarizationCalibrationParams})
    if cfg.type_display == :multiple
        # Multiple mode: internal states + derived curves in individual figures,
        # plus a dedicated calibration polarization figure.
        fig1 = Figure[]
        ax1 = Axis[]
        for _ in 1:5
            f, a = _make_single_axis(size=(760, 600))
            push!(fig1, f)
            push!(ax1, a)
        end
        fig2, ax2 = _make_single_axis(size=(760, 600))
        return (fig1, ax1, fig2, ax2, nothing, nothing)
    elseif cfg.type_display == :synthetic
        fig1, ax1 = _make_single_axis(size=(760, 600))
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end


"""Return figure/axes layout for EIS runs."""
function _create_figures(cfg::SimulationConfig{<:EISParams})
    if cfg.type_display == :multiple
        fig1, ax1 = _make_single_axis(size=(760, 620))
        fig2, ax2 = _make_single_axis(size=(760, 620))
        fig3, ax3 = _make_single_axis(size=(760, 620))
        return (fig1, ax1, fig2, ax2, fig3, ax3)
    elseif cfg.type_display == :synthetic
        fig1, ax1 = _make_row_axes(3; size=(1420, 480))
        return (fig1, ax1, nothing, nothing, nothing, nothing)
    end
    return (nothing, nothing, nothing, nothing, nothing, nothing)
end


"""Configure the appropriate Makie backend and return the figure/axes layout required for `cfg`."""
function figures_preparation(cfg::SimulationConfig, nb_gc::Union{Nothing, Integer}=nothing)
    cfg.type_display == :no_display &&
        return (nothing, nothing, nothing, nothing, nothing, nothing)
    _setup_makie_theme!(cfg)
    figs = cfg.type_current isa StepParams ? _create_figures(cfg, nb_gc) : _create_figures(cfg)
    _store_prepared_figures!(figs[1], figs[3], figs[5])
    # Dynamic mode opens windows immediately so plots can be refreshed while solving.
    _active_display_backend[] == :gl && cfg.display_timing == :live &&
        _open_interactive_figures!(figs[1], figs[3], figs[5])
    return figs
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Internal-state helpers
# ═══════════════════════════════════════════════════════════════════════════════


"""
Extract final internal states from a simulation to initialize the next segment.

For DAE segmented runs, the restart vector must keep the full canonical solver
layout (differential + algebraic) so that IDA resumes from the exact end state
without external algebraic re-solves.
"""
function _extract_last_internal_state(simu::AlphaPEM)::Vector{Float64}
    simu.sol === nothing && throw(ArgumentError("Cannot extract internal state before a successful solve."))

    solver_state_scaling = build_internal_solver_state_scaling(simu.fuel_cell, simu.cfg; include_algebraic=true)
    length(solver_state_scaling) == length(simu.sol.u[end]) ||
        throw(ArgumentError("Internal solver scaling size mismatch in _extract_last_internal_state."))

    # `simu.sol.u[end]` is stored in scaled solver coordinates.
    # Convert back to physical units because `simulate_model!` expects physical
    # initial values before applying its internal scaling pipeline.
    return unscale_values(simu.sol.u[end], solver_state_scaling)
end

