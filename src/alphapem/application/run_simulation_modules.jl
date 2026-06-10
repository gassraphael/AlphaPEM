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

using CairoMakie, GLMakie, WGLMakie

const _interactive_screens = Any[]
const _interactive_fig_to_screen = IdDict{Any, Any}()
const _active_display_backend = Ref{Symbol}(:cairo)
const _prepared_interactive_figures = Any[]

"""Return true when the simulation should be shown in a native interactive window."""
_use_interactive_display(cfg::SimulationConfig) =
    cfg.type_display != :no_display

"""Apply a publication-oriented Makie theme and activate the required backend."""
function _setup_makie_theme!(cfg::SimulationConfig; backend::Symbol=:auto)
    publication_palette = [
        "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9",
        "#E69F00", "#882255", "#44AA99", "#999933", "#332288", "#DDCC77",
    ]

    if backend == :wgl
        WGLMakie.activate!()
        _active_display_backend[] = :wgl
    elseif backend == :cairo
        CairoMakie.activate!()
        _active_display_backend[] = :cairo
    elseif backend == :gl
        try
            GLMakie.activate!()
            _active_display_backend[] = :gl
        catch err
            @warn "GLMakie could not be activated; falling back to CairoMakie non-interactive display." exception=(err, catch_backtrace())
            CairoMakie.activate!()
            _active_display_backend[] = :cairo
        end
    elseif _use_interactive_display(cfg)
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
function figures_preparation(cfg::SimulationConfig,
                             nb_gc::Union{Nothing, Integer}=nothing;
                             backend::Symbol=:auto)
    cfg.type_display == :no_display &&
        return (nothing, nothing, nothing, nothing, nothing, nothing)
    _setup_makie_theme!(cfg; backend=backend)
    figs = cfg.type_current isa StepParams ? _create_figures(cfg, nb_gc) : _create_figures(cfg)
    _store_prepared_figures!(figs[1], figs[3], figs[5])
    return figs
end


"""Return web-plot metadata in the same order as figures are populated by `display!`."""
function _web_plot_specs(cfg::SimulationConfig, nb_gc::Integer)::Vector{Dict{String, String}}
    has_extended_gc_profiles = nb_gc >= 3 && cfg.display_timing == :postrun
    specs = Dict{String, String}[]

    add_spec!(filename::String, title::String, group::String, key::String) =
        push!(specs, Dict(
            "filename" => filename,
            "title"    => title,
            "group"    => group,
            "key"      => key,
        ))

    if cfg.type_current isa StepParams
        add_spec!("01_ifc_1d_temporal.html",     "Current density evolution along the MEA",   "Internal states",      "ifc_1d_temporal")
        add_spec!("02_ucell_temporal.html",      "Cell voltage",                               "Internal states",      "Ucell_temporal")
        add_spec!("03_t_1d_temporal.html",       "Temperature evolution",                      "Internal states",      "T_1D_temporal")
        add_spec!("04_cv_1d_temporal.html",      "Water vapour concentration evolution",       "Internal states",      "Cv_1D_temporal")
        add_spec!("05_s_1d_temporal.html",       "Liquid saturation evolution",                "Internal states",      "s_1D_temporal")
        add_spec!("06_lambda_1d_temporal.html",  "Membrane water content evolution",           "Internal states",      "lambda_1D_temporal")
        add_spec!("07_ch2_1d_temporal.html",     "Hydrogen concentration evolution",           "Internal states",      "CH2_1D_temporal")
        add_spec!("08_co2_1d_temporal.html",     "Oxygen concentration evolution",             "Internal states",      "CO2_1D_temporal")
        add_spec!("09_p_1d_temporal.html",       "Pressure evolution",                         "Internal states",      "P_1D_temporal")
        add_spec!("10_cn2_1d_temporal.html",     "Nitrogen concentration evolution",           "Internal states",      "CN2_1D_temporal")
        add_spec!("11_phi_a_1d_temporal.html",   "Anode relative humidity evolution",          "Internal states",      "Phi_a_1D_temporal")
        add_spec!("12_phi_c_1d_temporal.html",   "Cathode relative humidity evolution",        "Internal states",      "Phi_c_1D_temporal")
        add_spec!("13_v_1d_temporal.html",       "Gas velocity evolution",                     "Internal states",      "v_1D_temporal")
        add_spec!("14_re_1d_temporal.html",      "Reynolds number evolution",                  "Internal states",      "Re_1D_temporal")
        if has_extended_gc_profiles
            add_spec!("15_t_pseudo_2d_final.html",      "Final pseudo-2D temperature map",             "Final gas-channel profiles", "T_pseudo_2D_final")
            add_spec!("16_ifc_gc_final.html",           "Final gas-channel current-density profile",    "Final gas-channel profiles", "ifc_GC_final")
            add_spec!("17_co2_pt_gc_final.html",        "Final Pt oxygen concentration profile",        "Final gas-channel profiles", "CO2_Pt_GC_final")
            add_spec!("18_lambda_mem_gc_final.html",    "Final membrane water-content profile",         "Final gas-channel profiles", "lambda_mem_GC_final")
        end
    elseif cfg.type_current isa PolarizationParams
        add_spec!("01_polarization_curve.html",  "Polarization curve",                         "Performance curves",   "polarization_curve")
        add_spec!("02_power_density_curve.html", "Power-density curve",                        "Performance curves",   "power_density_curve")
        add_spec!("03_efficiency_curve.html",    "Cell efficiency",                            "Performance curves",   "efficiency_curve")
        add_spec!("04_ifc_1d_temporal.html",     "Current density evolution along the MEA",   "Internal states",      "ifc_1d_temporal")
        add_spec!("05_ucell_temporal.html",      "Cell voltage",                               "Internal states",      "Ucell_temporal")
        add_spec!("06_t_1d_temporal.html",       "Temperature evolution",                      "Internal states",      "T_1D_temporal")
    elseif cfg.type_current isa PolarizationCalibrationParams
        add_spec!("01_polarization_curve_cali.html",  "Calibration polarization curve",             "Performance curves",   "polarization_curve_cali")
        add_spec!("02_power_density_curve_cali.html", "Power-density curve",                        "Performance curves",   "power_density_curve_cali")
        add_spec!("03_efficiency_curve_cali.html",    "Cell efficiency",                            "Performance curves",   "efficiency_curve_cali")
        add_spec!("04_ifc_1d_temporal_cali.html",     "Current density evolution along the MEA",   "Internal states",      "ifc_1d_temporal_cali")
        add_spec!("05_ucell_temporal_cali.html",      "Cell voltage",                               "Internal states",      "Ucell_temporal_cali")
        add_spec!("06_t_1d_temporal_cali.html",       "Temperature evolution",                      "Internal states",      "T_1D_temporal_cali")
    elseif cfg.type_current isa EISParams
        add_spec!("01_nyquist_plot.html",        "Nyquist plot",                                "EIS",                  "Nyquist_plot")
        add_spec!("02_bode_amplitude_plot.html", "Bode amplitude plot",                         "EIS",                  "Bode_amplitude_curve")
        add_spec!("03_bode_angle_plot.html",     "Bode angle plot",                             "EIS",                  "Bode_angle_curve")
    end

    return specs
end


"""Flatten the `(fig1, fig2, fig3)` layout returned by `figures_preparation` in display order."""
function _collect_web_figures(fig1, fig2, fig3)::Vector{Any}
    figures = Any[]
    for fig in (fig1, fig2, fig3)
        fig === nothing && continue
        if fig isa AbstractVector
            append!(figures, fig)
        else
            push!(figures, fig)
        end
    end
    return figures
end


"""Prepare interactive WGLMakie figures for a completed simulation, keeping them in memory.

This function calls the same `plot.jl` routines used by the native display, so the
web interface renders identical figures. The figures are returned alongside their
metadata so a live Bonito server can host them as dynamic routes — this is what
guarantees real interactivity (zoom, pan, hover, ...) and avoids the limitations
of any HTML pre-export (which can only ever serialize a frozen snapshot).

# Arguments
- `simu::AlphaPEM`: Completed simulation instance (must have `outputs != nothing`).

# Returns
- `Tuple{Vector{Dict{String,String}}, Vector{Any}}`: `(specs, figures)`, with one
  spec per figure, in the same order. Returns `(Dict[], Any[])` if no outputs
  are available.
"""
function prepare_web_figures(simu::AlphaPEM)::Tuple{Vector{Dict{String, String}}, Vector{Any}}
    simu.outputs === nothing && return (Dict{String, String}[], Any[])

    # Temporarily switch to :multiple / :postrun for web rendering.
    orig_display = simu.cfg.type_display
    orig_timing  = simu.cfg.display_timing
    simu.cfg.type_display   = :multiple
    simu.cfg.display_timing = :postrun

    plot_specs = Dict{String, String}[]
    rendered_figures = Any[]

    try
        @debug "Activating WGLMakie for web plot generation..."
        WGLMakie.activate!()
        nb_gc  = simu.cfg.numerical_parameters.nb_gc
        @debug "Preparing figures with nb_gc=$nb_gc..."
        fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg, nb_gc; backend=:wgl)

        if fig1 === nothing
            @warn "No figures prepared; skipping web plot generation"
            return (Dict{String, String}[], Any[])
        end

        @debug "Populating figures with simulation data..."
        # Populate figures using the same routines as the native display.
        display!(simu, ax1, ax2, ax3)

        rendered_figures = _collect_web_figures(fig1, fig2, fig3)
        plot_specs = _web_plot_specs(simu.cfg, nb_gc)

        length(rendered_figures) == length(plot_specs) ||
            throw(ArgumentError("Web plot mismatch: $(length(rendered_figures)) figures rendered for $(length(plot_specs)) plot specifications."))

    catch e
        @error "Error during web plot preparation: $e" exception=e
        rethrow(e)
    finally
        simu.cfg.type_display   = orig_display
        simu.cfg.display_timing = orig_timing
    end

    @debug "Web plot preparation completed: $(length(plot_specs)) figures ready"
    return (plot_specs, rendered_figures)
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
function _extract_last_internal_state(simu::AlphaPEM)::Tuple{Vector{Float64}, Vector{Float64}}
    simu.sol === nothing && throw(ArgumentError("Cannot extract internal state before a successful solve."))

    solver_state_scaling = build_solver_state_scaling(simu.cfg; include_algebraic=true)
    length(solver_state_scaling) == length(simu.sol.u[end]) ||
        throw(ArgumentError("Internal solver scaling size mismatch in _extract_last_internal_state."))

    # `simu.sol.u[end]` is stored in scaled solver coordinates.
    # Convert back to physical units because `simulate_model!` expects physical
    # initial values before applying its internal scaling pipeline.
    return unscale_values(simu.sol.u[end], solver_state_scaling), copy(simu.sol.du[end])
end
