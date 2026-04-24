# -*- coding: utf-8 -*-

"""
    AlphaPEM.Application.run_simulation

Main execution module for AlphaPEM.
Core orchestration logic lives here (launch, dispatch, API), while utility
helpers are kept in `run_simulation_modules.jl`.

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
using ..Core.Models: AlphaPEM, simulate_model!, display!, save_plot!,
                     build_internal_solver_state_scaling, unscale_values

include(joinpath(@__DIR__, "run_simulation_modules.jl"))

# ═══════════════════════════════════════════════════════════════════════════════
#  Section 1 – Launch functions (single simulator)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    launch_AlphaPEM_for_step_current(simu) -> AlphaPEM

Run a step-current simulation.

- `:postrun` display timing : one full ODE solve followed by a display update.
- `:live` display timing    : incremental solves separated by `delta_t_dyn_step`
  (from `simu.fuel_cell.numerical_parameters`).
"""
function launch_AlphaPEM_for_step_current(simu::AlphaPEM)::AlphaPEM
    # Figures preparation
    nb_gc = simu.fuel_cell.numerical_parameters.nb_gc
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg, nb_gc)

    # Dynamic display requires a dedicated segmented simulation flow.
    if simu.cfg.display_timing == :live
        # Certain conditions must be met.
        if simu.cfg.type_display == :multiple
            throw(ArgumentError("step current cannot be run with `display_timing = :live` and `type_display = :multiple` (too many plots to handle)."))
        end

        # Extraction of parameters
        p       = simu.cfg.type_current       # ::StepParams
        tf_step = p.delta_t_ini + p.delta_t_load + p.delta_t_break
        delta_t_dyn_step = simu.fuel_cell.numerical_parameters.delta_t_dyn_step
        delta_t_dyn_step > 0 || throw(ArgumentError("delta_t_dyn_step must be > 0 for dynamic step runs."))

        # Initialization
        n             = floor(Int, tf_step / delta_t_dyn_step)
        t0            = 0.0
        initial_state = nothing

        # Dynamic simulation
        for i in 1:n
            tf = i * delta_t_dyn_step
            simulate_model!(simu, initial_state, (t0, tf))

            # Display
            simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)

            # Recovery of the internal states from the end of the preceding simulation.
            initial_state = _extract_last_internal_state(simu)
            # Time interval actualization.
            t0            = simu.outputs.solver.t[end]
        end

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
        # Initialization
        p            = simu.current_density            # ::PolarizationCurrent
        delta_t_step = step_duration(p)
        _, tf_full   = p.time_interval
        n            = round(Int, (tf_full - p.delta_t_ini) / delta_t_step)
        t0           = 0.0
        tf           = p.delta_t_ini + delta_t_step
        initial_state = nothing

        # Dynamic simulation
        for i in 1:n
            simulate_model!(simu, initial_state, (t0, tf))

            # Display
            simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)

            # Recovery of the internal states from the end of the preceding simulation.
            initial_state = _extract_last_internal_state(simu)
            # Time interval actualization.
            t0 = simu.outputs.solver.t[end]
            tf = p.delta_t_ini + (i + 1) * delta_t_step
        end

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
        # Initialization
        p          = simu.current_density           # ::PolarizationCalibrationCurrent
        step_dts   = step_duration(p)               # Vector{Float64}: one duration per exp. step
        n          = length(step_dts)
        boundaries = cumsum(step_dts)               # absolute end-time of each step (after t_ini)
        t0         = 0.0
        initial_state = nothing

        # Dynamic simulation
        for i in 1:n
            tf = p.delta_t_ini + boundaries[i]
            simulate_model!(simu, initial_state, (t0, tf))

            # Display
            simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)

            # Recovery of the internal states from the end of the preceding simulation.
            initial_state = _extract_last_internal_state(simu)
            # Time interval actualization.
            t0 = simu.outputs.solver.t[end]
        end

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

Run an EIS simulation.  One ODE segment is solved per frequency point.
"""
function launch_AlphaPEM_for_EIS_current(simu::AlphaPEM)::AlphaPEM
    # Check if `display_timing` is valid for EIS.
    simu.cfg.display_timing == :live ||
        throw(ArgumentError(
            "EIS requires `display_timing = :live` so that `max_step` can be " *
            "adjusted at each frequency."))

    # Warnings
    @warn "EIS simulation is currently not maintained.  Some unexpected bugs or " *
          "incorrect results may occur.  This will be resolved in a future release."

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg)

    # Initialization
    cd = simu.current_density   # ::EISCurrent
    n  = length(cd.t_new_start)

    # A preliminary simulation run is necessary to equilibrate internal variables at i_EIS.
    simulate_model!(simu, nothing, (0.0, cd.t0))
    # Recovery of the internal states from the end of the preceding simulation.
    initial_state = _extract_last_internal_state(simu)

    if simu.cfg.type_display == :multiple
        @warn "Dynamic graph refresh is not available with `:multiple` display " *
              "due to the volume of data involved.  Data are computed correctly " *
              "and saved to the 'results' folder.  Use `:synthetic` to avoid this."
    end

    # Dynamic simulation: one segment per frequency.
    for i in 1:n
        # Time interval actualization
        t0 = simu.outputs.solver.t[end]
        tf = cd.t_new_start[i] + cd.delta_t_break[i] + cd.delta_t_measurement[i]

        simulate_model!(simu, initial_state, (t0, tf))

        # Recovery of the internal states from the end of the preceding simulation.
        initial_state = _extract_last_internal_state(simu)

        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)
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
    for simu in simulators
        # Simulation
        simulate_model!(simu)
        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)
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
    for simu in simulators
        # Simulation
        simulate_model!(simu)
        # Display
        simu.cfg.type_display != :no_display && display!(simu, ax1, ax2, ax3)
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

"""Disable interactive mode and block until all figure windows are closed."""
function _finalize_display!(cfg::SimulationConfig)
    if _active_display_backend[] == :gl && _use_interactive_display(cfg)
        # Fixed mode: show the fully populated figures only once computations are done.
        if cfg.display_timing == :postrun && isempty(_interactive_screens)
            _open_interactive_figures!(_prepared_interactive_figures...)
        end

        for screen in copy(_interactive_screens)
            try
                wait(screen)
            catch err
                @warn "Interactive GLMakie screen wait failed." exception=(err, catch_backtrace())
            end
        end
        empty!(_interactive_screens)
        empty!(_interactive_fig_to_screen)
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
    type_fuel_cell = :ZSW_GenStack,
    type_current   = StepParams(),
    type_display   = :synthetic,
    display_timing = :postrun,
)
simu = run_simulation(cfg)
```
"""
function run_simulation(cfg::SimulationConfig)::AlphaPEM
    # Build a Fuelcell object with the given configuration.
    fuel_cell       = create_fuelcell(cfg.type_fuel_cell, cfg.voltage_zone)
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
    fuel_cells        = [create_fuelcell(cfgs[i].type_fuel_cell, cfgs[i].voltage_zone) for i in 1:n]
    # Build Current objects for each configuration.
    current_densities = [create_current(cfgs[i].type_current, fuel_cells[i])           for i in 1:n]
    # Create one simulator per selected fuel cell / current configuration.
    simulators        = [AlphaPEM(fuel_cells[i], current_densities[i], cfgs[i])        for i in 1:n]

    # Launch the simulations.
    _dispatch_launch!(simulators)
    _finalize_display!(cfgs)

    return simulators
end
