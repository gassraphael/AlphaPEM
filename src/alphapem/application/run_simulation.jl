# -*- coding: utf-8 -*-

"""
    AlphaPEM.Application.run_simulation

Main execution module for AlphaPEM.
Core orchestration logic lives here (launch, dispatch, API), while the rest is kept in:
- `run_simulation_modules.jl`: figure preparation and internal-state helpers
- `interactive_display.jl`: serving the figures in a web browser
- `run_simulation_live.jl`: segmented `:live` runs (a debugging aid, read it apart)

Public API
----------
- `run_simulation(cfg::SimulationConfig)  -> AlphaPEM`
- `run_simulation(cfgs::AbstractVector{<:SimulationConfig}) -> Vector{AlphaPEM}`
"""

# ─────────────────────────────────────── Imports ─────────────────────────────

using ..Config:   AbstractCurrentParams,
                  SimulationConfig,
                  StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using ..Fuelcell: create_fuelcell
using ..Currents: create_current, step_duration
using SciMLBase: ReturnCode
using ..Core.Models: AlphaPEM, simulate_model!, display!, save_plot!,
                     build_solver_state_scaling, unscale_values, EIS_STATE_RESPONSES
using ..Core.Types: SimulationOutputs

include(joinpath(@__DIR__, "run_simulation_modules.jl"))
include(joinpath(@__DIR__, "interactive_display.jl"))
include(joinpath(@__DIR__, "run_simulation_live.jl"))

# ═══════════════════════════════════════════════════════════════════════════════
#  Section 1 – Launch functions (single simulator)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    launch_AlphaPEM_for_step_current(simu) -> AlphaPEM

Run a step-current simulation.

- `:postrun` display timing : one full ODE solve followed by a display update.
- `:live` display timing    : incremental solves separated by `_LIVE_STEP_SEGMENT_DURATION`, with a
  display refresh after each (see `run_simulation_live.jl`).
"""
function launch_AlphaPEM_for_step_current(simu::AlphaPEM)::AlphaPEM
    # Figures preparation
    nb_gc = simu.cfg.numerical_parameters.nb_gc
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg, nb_gc)

    # Dynamic display requires a dedicated segmented simulation flow.
    if simu.cfg.display_timing == :live
        _simulate_live_step!(simu, ax1, ax2, ax3)
    else  # :postrun
        # Simulation
        simulate_model!(simu)
        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)
    end

    # Plot saving
    save_plot!(simu, fig1, fig2, fig3)
    return simu
end


"""
    launch_AlphaPEM_for_polarization_current(simu) -> AlphaPEM

Run a polarization-curve simulation.

- `:postrun` display timing : one full solve.
- `:live` display timing    : one solve per current step, with display refresh after each.
"""
function launch_AlphaPEM_for_polarization_current(simu::AlphaPEM)::AlphaPEM
    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg)

    # Dynamic display requires a dedicated segmented simulation flow.
    if simu.cfg.display_timing == :live
        _simulate_live_polarization!(simu, ax1, ax2, ax3)
    else  # :postrun
        # Simulation
        simulate_model!(simu)
        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)
    end

    # Plot saving
    save_plot!(simu, fig1, fig2, fig3)
    return simu
end


"""
    launch_AlphaPEM_for_polarization_current_for_calibration(simu) -> AlphaPEM

Run a calibration polarization simulation.

- `:postrun` display timing : one full solve.
- `:live` display timing    : one solve per experimental current step (variable duration).
"""
function launch_AlphaPEM_for_polarization_current_for_calibration(simu::AlphaPEM)::AlphaPEM
    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg)

    # Dynamic display requires a dedicated segmented simulation flow.
    if simu.cfg.display_timing == :live
        _simulate_live_calibration!(simu, ax1, ax2, ax3)
    else  # :postrun
        # Simulation
        simulate_model!(simu)
        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)
    end

    # Plot saving
    save_plot!(simu, fig1, fig2, fig3)
    return simu
end


"""
    launch_AlphaPEM_for_EIS_current(simu) -> AlphaPEM

Run an EIS simulation.
- In `:live` mode: one ODE segment is solved per frequency point with intermediate displays.
- In `:postrun` mode: the sweep is solved frequency by frequency, followed by a single display.
"""
function launch_AlphaPEM_for_EIS_current(simu::AlphaPEM)::AlphaPEM
    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg)

    if simu.cfg.display_timing == :live
        _simulate_live_EIS!(simu, ax1, ax2, ax3)
    else
        # :postrun mode: segment the solve by frequency, exactly like :live, but without
        # displaying between segments. A single continuous solve over the whole EIS sweep
        # would need one saveat point every ~1/(nb_points*f) seconds for every frequency,
        # accumulating sol.t/sol.u to hundreds of thousands of entries in one live array —
        # risking a max_run_time_s timeout long before the sweep completes. Segmenting keeps
        # each frequency's sol.t/sol.u small and gives it a fresh maxiters budget, then the
        # segments are concatenated for a single final display.

        # Initialization
        cd = simu.current_density   # ::EISCurrent
        n  = length(cd.t_new_start)
        initial_variable_values = nothing
        initial_derivative_values = nothing


        # A preliminary simulation run is necessary to equilibrate internal variables at i_EIS.
        simulate_model!(simu, initial_variable_values, initial_derivative_values, (0.0, cd.t0))
        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values, initial_derivative_values = _extract_last_internal_state(simu)
        # Conserving the outputs of the preliminary run for later concatenation.
        accumulated = simu.outputs

        # Dynamic simulation: one segment per frequency.
        for i in 1:n
            # Time interval actualization
            t0 = accumulated.solver.t[end]
            tf = cd.t_new_start[i] + cd.delta_t_break[i] + cd.delta_t_measurement[i]

            simulate_model!(simu, initial_variable_values, initial_derivative_values, (t0, tf))

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values, initial_derivative_values = _extract_last_internal_state(simu)

            # Concatenate the outputs of this segment with the accumulated outputs of previous segments.
            _append_simulation_segment!(accumulated, simu.outputs)
        end

        # Replace the simulator's outputs with the accumulated outputs of all segments.
        simu.outputs = accumulated

        # Trigger a single display! call at the end.
        if simu.cfg.type_display != :no_display
            display!(simu, ax1, ax2, ax3)
        end
    end

    # Plot saving
    save_plot!(simu, fig1, fig2, fig3)
    return simu
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 2 – Launch functions (multi-simulator batch runs)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    launch_AlphaPEM_for_polarization_current(simulators) -> Vector{AlphaPEM}

Run a batch of polarization-curve simulations sequentially, sharing one set of axes.
"""
function launch_AlphaPEM_for_polarization_current(
    simulators::AbstractVector{<:AlphaPEM},
)::Vector{<:AlphaPEM}
    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(first(simulators).cfg)
    for (i, simu) in enumerate(simulators)
        # Simulation
        simulate_model!(simu)
        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3; series_index=i)
        # Release the solver solution before starting the next simulation.
        simu.sol = nothing
        GC.gc(true)
    end
    # Plot saving
    save_plot!(first(simulators), fig1, fig2, fig3)
    return simulators
end


"""
    launch_AlphaPEM_for_polarization_current_for_calibration(simulators) -> Vector{AlphaPEM}

Run a batch of calibration polarization simulations sequentially, sharing one set of axes.
"""
function launch_AlphaPEM_for_polarization_current_for_calibration(
    simulators::AbstractVector{<:AlphaPEM},
)::Vector{<:AlphaPEM}
    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(first(simulators).cfg)
    for (i, simu) in enumerate(simulators)
        # Simulation
        simulate_model!(simu)
        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3; series_index=i)
        # Release the solver solution before starting the next simulation.
        simu.sol = nothing
        GC.gc(true)
    end
    # Plot saving
    save_plot!(first(simulators), fig1, fig2, fig3)
    return simulators
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 3 – Dispatch helpers
#
#  These thin wrappers let `run_simulation` stay free of if/elseif chains.
#  The compiler resolves the correct launch function at compile time from the
#  concrete type of `cfg.type_current`.
# ═══════════════════════════════════════════════════════════════════════════════

_dispatch_launch!(simu::AlphaPEM) =
    _dispatch_launch!(simu, simu.cfg.type_current)

_dispatch_launch!(simu::AlphaPEM, ::StepParams) =
    launch_AlphaPEM_for_step_current(simu)

_dispatch_launch!(simu::AlphaPEM, ::PolarizationParams) =
    launch_AlphaPEM_for_polarization_current(simu)

_dispatch_launch!(simu::AlphaPEM, ::PolarizationCalibrationParams) =
    launch_AlphaPEM_for_polarization_current_for_calibration(simu)

_dispatch_launch!(simu::AlphaPEM, ::EISParams) =
    launch_AlphaPEM_for_EIS_current(simu)


function _dispatch_launch!(simulators::AbstractVector{<:AlphaPEM})
    isempty(simulators) && throw(ArgumentError("At least one simulator is required for a multi-simulator run."))

    # Homogeneity check: all current-profile types must match in multi-simulator mode.
    T = typeof(first(simulators).cfg.type_current)
    all(s -> s.cfg.type_current isa T, simulators) ||
        throw(ArgumentError(
            "All configurations in a multi-simulator run must share the same " *
            "current profile type (got mixed types)."))
    return _dispatch_launch!(simulators, first(simulators).cfg.type_current)
end

_dispatch_launch!(sims::AbstractVector{<:AlphaPEM}, ::PolarizationParams) =
    launch_AlphaPEM_for_polarization_current(sims)

_dispatch_launch!(sims::AbstractVector{<:AlphaPEM}, ::PolarizationCalibrationParams) =
    launch_AlphaPEM_for_polarization_current_for_calibration(sims)

function _dispatch_launch!(sims::AbstractVector{<:AlphaPEM}, ::AbstractCurrentParams)
    throw(ArgumentError(
        "Multi-simulator runs only support `PolarizationParams` and " *
        "`PolarizationCalibrationParams` current profiles."))
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 4 – Display finalisation
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _finalize_display!(cfg)
    _finalize_display!(cfgs)

Serve the final figures and, in a script, keep them served until Ctrl+C.

In an interactive Julia session (REPL), the call returns immediately and the figures stay
served until the next run replaces them (see `_plot_block_requested`).

Parameters
----------
cfg : SimulationConfig
    Configuration of the run (or `cfgs`, the configurations of a batch, of which the first
    one decides the display).

Returns
-------
Nothing
"""
function _finalize_display!(cfg::SimulationConfig)
    if _active_display_backend[] == :wgl && _use_interactive_display(cfg)
        # Fixed mode: show the fully populated figures only once computations are done.
        if cfg.display_timing == :postrun
            _open_interactive_figures!()
        end

        if !isempty(_served_figure_routes) && _plot_block_requested()
            _serve_until_interrupted!("Press Ctrl+C to stop serving the interactive figures.")
        end
    end
    empty!(_prepared_interactive_figures)
    return nothing
end

function _finalize_display!(cfgs::AbstractVector{<:SimulationConfig})
    if !isempty(cfgs)
        _finalize_display!(first(cfgs))
    end
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Section 5 – Public API
# ═══════════════════════════════════════════════════════════════════════════════

"""
    run_simulation(cfg::SimulationConfig) -> AlphaPEM

Run a single PEM fuel cell simulation according to `cfg`.

# Arguments
- `cfg::SimulationConfig`: Configuration object specifying the fuel-cell type,
  current profile, display mode, and plot style.

# Returns
- `AlphaPEM`: The completed simulator instance with all outputs populated.

# Example
```julia
using AlphaPEM.Config: SimulationConfig, StepParams
using AlphaPEM.Application: run_simulation

cfg = SimulationConfig(
    type_fuel_cell = :ZSW_nominal,
    year           = 2024,
    type_current   = StepParams(),
    type_display   = :synthetic,
    display_timing = :postrun,
)
simu = run_simulation(cfg)
```
"""
function run_simulation(cfg::SimulationConfig)::AlphaPEM
    # Build a Fuelcell object with the given configuration.
    fuel_cell       = create_fuelcell(cfg.type_fuel_cell, cfg.voltage_zone; year=cfg.year, nb_gc=cfg.numerical_parameters.nb_gc)

    # Override with custom parameters if provided
    if cfg.physical_parameters !== nothing
        fuel_cell.physical_parameters = cfg.physical_parameters
    end
    if cfg.operating_conditions !== nothing
        fuel_cell.operating_conditions = cfg.operating_conditions
    end

    # Build a Current object with the given configuration.
    current_density = create_current(cfg.type_current, fuel_cell)
    # Create a simulator.
    simu            = AlphaPEM(fuel_cell, current_density, cfg)

    # Launch the simulation.
    _dispatch_launch!(simu)
    _finalize_display!(cfg)

    return simu
end


"""
    run_simulation(cfgs::AbstractVector{<:SimulationConfig}) -> Vector{AlphaPEM}

Run a batch of PEM fuel cell simulations.

All configurations must share the same current profile type
(`PolarizationParams` or `PolarizationCalibrationParams`).  The simulators share
a single set of figure axes so that their results can be compared on the same plot.

# Arguments
- `cfgs::AbstractVector{<:SimulationConfig}`: Vector of configuration objects.

# Returns
- `Vector{AlphaPEM}`: Completed simulator instances in the same order as `cfgs`.
"""
function run_simulation(cfgs::AbstractVector{<:SimulationConfig})::Vector{AlphaPEM}
    isempty(cfgs) && throw(ArgumentError("`cfgs` must contain at least one SimulationConfig."))

    # Determine the number of simulators to create.
    n                 = length(cfgs)
    # Build Fuelcell objects for each configuration.
    fuel_cells = []
    for i in 1:n
        fc = create_fuelcell(cfgs[i].type_fuel_cell, cfgs[i].voltage_zone; year=cfgs[i].year, nb_gc=cfgs[i].numerical_parameters.nb_gc)
        if cfgs[i].physical_parameters !== nothing
            fc.physical_parameters = cfgs[i].physical_parameters
        end
        if cfgs[i].operating_conditions !== nothing
            fc.operating_conditions = cfgs[i].operating_conditions
        end
        push!(fuel_cells, fc)
    end
    # Build Current objects for each configuration.
    current_densities = [create_current(cfgs[i].type_current, fuel_cells[i])           for i in 1:n]
    # Create one simulator per selected fuel cell / current configuration.
    simulators        = [AlphaPEM(fuel_cells[i], current_densities[i], cfgs[i])        for i in 1:n]

    # Launch the simulations.
    _dispatch_launch!(simulators)
    _finalize_display!(cfgs)

    return simulators
end
