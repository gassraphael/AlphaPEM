# -*- coding: utf-8 -*-

"""
Segmented `:live` runs and their display.

With `display_timing = :live`, the simulation is solved in short segments and the figures
are refreshed between the segments, so that the internal states can be watched while the
run progresses. It is a debugging aid, so everything that is specific to it lives here (and
in `core/modules/fuelcell/plot_live.jl`), apart from the launch functions of
`run_simulation.jl`, which only call the functions of the first section.

Sections:
1. Segmented live runs (one function per current profile)
2. Live display (history, refresh policy, failure display)
"""

using ..Core.Models: recovery!

# Live display timing (kept out of the general numerical-solver parameters).
const _LIVE_REFRESH_INTERVAL = 0.2  # Minimum wall-clock time between refreshes (s); 0 disables throttling.
const _LIVE_FULL_RATE_TIME = 10.0   # Simulated time during which every segment is displayed (s).
const _LIVE_STEP_SEGMENT_DURATION = 0.1  # Simulated duration of each step-current segment (s).

# ═══════════════════════════════════════════════════════════════════════════════
#  1. Segmented live runs
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _simulate_live_step!(simu, ax1, ax2, ax3)

Solve a step-current run in segments of `_LIVE_STEP_SEGMENT_DURATION`, refreshing the figures after
each segment.

The figures show the whole trajectory since the beginning of the run. Every segment is
displayed during the first `_LIVE_FULL_RATE_TIME` simulated seconds, then the refreshes
are spaced by at least `_LIVE_REFRESH_INTERVAL` seconds (0 = never throttled). The last
segment, and the trajectory reached before a failure, are always displayed. Both
`:synthetic` and `:multiple` displays are supported.

Parameters
----------
simu : AlphaPEM
    Simulator to run.
ax1, ax2, ax3 : Any
    Axes returned by `figures_preparation` (`nothing` when nothing is displayed).

Returns
-------
Nothing
    On return, `simu.outputs` holds the whole trajectory.
"""
function _simulate_live_step!(simu::AlphaPEM, ax1, ax2, ax3)
    # Extraction of parameters
    p       = simu.cfg.type_current       # ::StepParams
    tf_step = p.delta_t_ini + p.delta_t_load + p.delta_t_break
    _LIVE_STEP_SEGMENT_DURATION > 0 ||
        throw(ArgumentError("_LIVE_STEP_SEGMENT_DURATION must be > 0 for dynamic step runs."))

    # Initialization
    n             = floor(Int, tf_step / _LIVE_STEP_SEGMENT_DURATION)
    t0            = 0.0
    initial_variable_values = nothing
    initial_derivative_values = nothing
    history  = nothing                 # outputs accumulated since the beginning of the run.
    throttle = _LiveThrottle()

    # Segmented simulation
    try
        for i in 1:n
            tf = i * _LIVE_STEP_SEGMENT_DURATION
            simulate_model!(simu, initial_variable_values, initial_derivative_values, (t0, tf))
            history = _accumulate_live_outputs!(simu, history)

            # Display: the figures are served after the first segment, then refreshed with the
            # whole trajectory; the last segment is always displayed.
            if simu.cfg.type_display != :no_display
                _refresh_live_display!(simu, ax1, ax2, ax3; throttle=throttle, force=(i == n))
            end

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values, initial_derivative_values = _extract_last_internal_state(simu)
            # Time interval actualization.
            t0            = simu.outputs.solver.t[end]
        end
    catch err
        # On failure, show the states reached so far (the ones that matter to debug).
        err isa InterruptException || _display_live_failure!(simu, history, ax1, ax2, ax3)
        rethrow()
    end
    return nothing
end

"""
    _simulate_live_polarization!(simu, ax1, ax2, ax3)

Solve a polarization-curve run with one segment per current step, refreshing the figures
after each segment. The run stops early when the solver stops on a non-physical state.

Parameters
----------
simu : AlphaPEM
    Simulator to run.
ax1, ax2, ax3 : Any
    Axes returned by `figures_preparation` (`nothing` when nothing is displayed).

Returns
-------
Nothing
"""
function _simulate_live_polarization!(simu::AlphaPEM, ax1, ax2, ax3)
    # Initialization
    p          = simu.current_density            # ::PolarizationCurrent
    dt_steps = step_duration(p)                  # Vector{Float64}: one duration per transition
    n          = length(dt_steps)
    tf_steps = cumsum(dt_steps)                # absolute end-time of each step
    t0         = 0.0
    initial_variable_values = nothing
    initial_derivative_values = nothing

    # Segmented simulation
    for i in 1:n
        tf = tf_steps[i]
        simulate_model!(simu, initial_variable_values, initial_derivative_values, (t0, tf))

        # Safety stop: cell voltage or O2 concentration became non-positive.
        simu.sol.retcode === ReturnCode.Terminated && break

        # Display: the first call serves the figures now that they have content;
        # later calls only refresh them.
        simu.cfg.type_display != :no_display && _refresh_live_display!(simu, ax1, ax2, ax3)

        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values, initial_derivative_values = _extract_last_internal_state(simu)
        # Time interval actualization.
        t0 = simu.outputs.solver.t[end]
    end
    return nothing
end

"""
    _simulate_live_calibration!(simu, ax1, ax2, ax3)

Solve a calibration polarization run with one segment per experimental current step
(variable duration), refreshing the figures after each segment. The run stops early when
the solver stops on a non-physical state.

Parameters
----------
simu : AlphaPEM
    Simulator to run.
ax1, ax2, ax3 : Any
    Axes returned by `figures_preparation` (`nothing` when nothing is displayed).

Returns
-------
Nothing
"""
function _simulate_live_calibration!(simu::AlphaPEM, ax1, ax2, ax3)
    # Initialization
    p          = simu.current_density           # ::PolarizationCalibrationCurrent
    step_dts   = step_duration(p)               # Vector{Float64}: one duration per exp. step
    n          = length(step_dts)
    boundaries = cumsum(step_dts)               # absolute end-time of each step (after t_ini)
    t0         = 0.0
    initial_variable_values = nothing
    initial_derivative_values = nothing

    # Segmented simulation
    for i in 1:n
        tf = boundaries[i]
        simulate_model!(simu, initial_variable_values, initial_derivative_values, (t0, tf))

        # Safety stop: cell voltage or O2 concentration became non-positive.
        simu.sol.retcode === ReturnCode.Terminated && break

        # Display
        simu.cfg.type_display != :no_display && _refresh_live_display!(simu, ax1, ax2, ax3)

        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values, initial_derivative_values = _extract_last_internal_state(simu)
        # Time interval actualization.
        t0 = simu.outputs.solver.t[end]
    end
    return nothing
end

"""
    _simulate_live_EIS!(simu, ax1, ax2, ax3)

Solve an EIS run with one segment per frequency point, refreshing the figures after each
segment. The `:multiple` display is not refreshed dynamically (too much data): a warning is
emitted and the figures are only saved at the end.

Parameters
----------
simu : AlphaPEM
    Simulator to run.
ax1, ax2, ax3 : Any
    Axes returned by `figures_preparation` (`nothing` when nothing is displayed).

Returns
-------
Nothing
"""
function _simulate_live_EIS!(simu::AlphaPEM, ax1, ax2, ax3)
    # Initialization
    cd = simu.current_density   # ::EISCurrent
    n  = length(cd.t_new_start)
    initial_variable_values = nothing
    initial_derivative_values = nothing

    # A preliminary simulation run is necessary to equilibrate internal variables at i_EIS.
    simulate_model!(simu, initial_variable_values, initial_derivative_values, (0.0, cd.t0))
    # Recovery of the internal states from the end of the preceding simulation.
    initial_variable_values, initial_derivative_values = _extract_last_internal_state(simu)

    if simu.cfg.type_display == :multiple
        @warn "Dynamic graph refresh is not available with `:multiple` display " *
              "due to the volume of data involved.  Data are computed correctly " *
              "and saved to the 'results' folder.  Use `:synthetic` to avoid this."
    end

    # Segmented simulation: one segment per frequency.
    for i in 1:n
        # Time interval actualization
        t0 = simu.outputs.solver.t[end]
        tf = cd.t_new_start[i] + cd.delta_t_break[i] + cd.delta_t_measurement[i]

        simulate_model!(simu, initial_variable_values, initial_derivative_values, (t0, tf))

        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values, initial_derivative_values = _extract_last_internal_state(simu)

        # Display
        simu.cfg.type_display != :no_display && _refresh_live_display!(simu, ax1, ax2, ax3)
    end
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  2. Live display
# ═══════════════════════════════════════════════════════════════════════════════

# Each refresh sends the new data to the browser, which can take far longer than solving one
# short segment, so long runs throttle the refreshes by wall-clock time. The throttle never
# hides what matters for debugging: every segment is displayed during the first
# `_LIVE_FULL_RATE_TIME` simulated seconds, the last segment is always displayed, and the
# trajectory reached before a failure is displayed too.

const _LIVE_COST_RATIO = 3.0   # a throttled refresh may use at most ~1/(1+ratio) of the wall time.

"""
    _LiveThrottle

Refresh policy and wall-clock bookkeeping of a `:live` run. `_LIVE_REFRESH_INTERVAL = 0`
disables the throttle so that every segment is displayed.

Fields
------
min_interval : Float64
    Minimum wall-clock time between two refreshes (s).
full_rate_time : Float64
    Simulated time during which every segment is displayed (s).
last_end : Float64
    Wall-clock time at which the last refresh ended.
last_cost : Float64
    Duration of the last refresh (the first one, dominated by compilation, is ignored).
refreshes : Int
    Number of refreshes done so far.
"""
mutable struct _LiveThrottle
    min_interval::Float64
    full_rate_time::Float64
    last_end::Float64
    last_cost::Float64
    refreshes::Int
end
function _LiveThrottle()
    _LIVE_REFRESH_INTERVAL >= 0.0 ||
        throw(ArgumentError("_LIVE_REFRESH_INTERVAL must be >= 0."))
    _LIVE_FULL_RATE_TIME >= 0.0 ||
        throw(ArgumentError("_LIVE_FULL_RATE_TIME must be >= 0."))
    return _LiveThrottle(_LIVE_REFRESH_INTERVAL, _LIVE_FULL_RATE_TIME, -Inf, 0.0, 0)
end

"""
    _live_refresh_due(th, t_sim) -> Bool

Tell whether the refresh policy asks for a display at simulated time `t_sim`.

Parameters
----------
th : _LiveThrottle
    Refresh policy and bookkeeping.
t_sim : Real
    Simulated time reached (s).

Returns
-------
Bool
    True during the full-rate window, when the throttle is disabled, or when the minimum
    interval (and `_LIVE_COST_RATIO` times the cost of the last refresh) has elapsed.
"""
function _live_refresh_due(th::_LiveThrottle, t_sim::Real)
    th.min_interval <= 0.0 && return true
    t_sim <= th.full_rate_time && return true
    return time() - th.last_end >= max(th.min_interval, _LIVE_COST_RATIO * th.last_cost)
end

"""
    _is_plot_session_error(err) -> Bool

Tell whether `err` was raised because the browser tab holding a figure is closed, still
loading or unresponsive. Such errors must not stop the simulation.

Parameters
----------
err : Any
    Caught exception.

Returns
-------
Bool
    True for the WGLMakie/Bonito session errors.
"""
_is_plot_session_error(err) =
    err isa ErrorException &&
    (occursin("Can't get three", err.msg) || occursin("Session", err.msg) || occursin("Timed out", err.msg))

"""
    _prune_closed_screens!()

Detach the WGLMakie screens whose browser session is closed from every scene of the
served figures.

WGLMakie keeps such screens attached to the scenes, and inserting a plot then raises an
error for the closed session, which would freeze the live view for good after a tab is
closed or reloaded.

References
----------
- WGLMakie, https://docs.makie.org/stable/explanations/backends/wglmakie
"""
function _prune_closed_screens!()
    Bonito = WGLMakie.Bonito
    closed(s) = s isa WGLMakie.Screen && s.session !== nothing &&
                s.session.status in (Bonito.CLOSED, Bonito.SOFT_CLOSED)

    function prune!(scene)
        filter!(!closed, scene.current_screens)
        foreach(prune!, scene.children)
    end

    for fig in _collect_web_figures(_prepared_interactive_figures...)
        prune!(fig.scene)
    end
    return nothing
end

"""
    _refresh_live_display!(simu, ax1, ax2, ax3; throttle=nothing, force=false)

Refresh the interactive figures of a `:live` run with the current simulation history.

Without `throttle` (or with `force=true`) the figures are always refreshed. With a
`throttle`, the refresh is skipped unless the policy asks for it (see `_LiveThrottle`), so
that the display does not slow long runs down much. A closed, still-loading or unresponsive
browser tab does not stop the simulation: the refresh is skipped, and the view resumes at
the next refresh once the page is (re)loaded.

Parameters
----------
simu : AlphaPEM
    Simulator whose `outputs` hold the trajectory to display.
ax1, ax2, ax3 : Any
    Axes returned by `figures_preparation`.
throttle : Union{Nothing, _LiveThrottle}
    Refresh policy; `nothing` refreshes every time.
force : Bool
    Refresh even if the policy does not ask for it (last segment, failure).

Returns
-------
Nothing
"""
function _refresh_live_display!(simu::AlphaPEM, ax1, ax2, ax3;
                                throttle::Union{Nothing, _LiveThrottle}=nothing,
                                force::Bool=false)
    yield()   # lets the plot server answer pending browser requests.
    if throttle !== nothing && !force && !_live_refresh_due(throttle, simu.outputs.solver.t[end])
        return nothing
    end

    # Refresh. Bonito flags a closed tab late: the first attempt can still fail on it, after
    # which the stale screen is known and can be detached for the second attempt.
    t_begin = time()
    for attempt in 1:2
        try
            _active_display_backend[] == :wgl && _prune_closed_screens!()
            display!(simu, ax1, ax2, ax3)
            _active_display_backend[] == :wgl && _open_interactive_figures!()
            break
        catch err
            _is_plot_session_error(err) || rethrow()
            attempt == 2 &&
                @warn "Live display update skipped (browser tab closed, still loading or unresponsive); the simulation continues. Reload the page to resume the live view." exception=err maxlog=3
        end
    end

    # Bookkeeping of the refresh policy
    if throttle !== nothing
        throttle.refreshes > 0 && (throttle.last_cost = time() - t_begin)
        throttle.refreshes += 1
        throttle.last_end = time()
    end
    yield()
    return nothing
end

"""
    _accumulate_live_outputs!(simu, history) -> SimulationOutputs

Accumulate the outputs of the segment just solved into the history of a `:live` run.

The first segment becomes the history; later segments are appended to it (without their
duplicated restart point). `simu.outputs` is then set to the full history, so that
`display!` and `save_plot!` see the whole trajectory since the beginning of the run.

Parameters
----------
simu : AlphaPEM
    Simulator whose `outputs` hold the segment that has just been solved.
history : Union{Nothing, SimulationOutputs}
    History accumulated so far (`nothing` before the first segment).

Returns
-------
SimulationOutputs
    The updated history (also stored in `simu.outputs`).
"""
function _accumulate_live_outputs!(simu::AlphaPEM, history)
    segment = simu.outputs
    segment === nothing && throw(ArgumentError("No segment outputs to accumulate."))
    if history === nothing
        history = segment
    elseif segment !== history
        _append_simulation_segment!(history, segment; skip_first=true)
    end
    simu.outputs = history
    return history
end

"""
    _display_live_failure!(simu, history, ax1, ax2, ax3)

Display the trajectory reached when a `:live` run fails, so that the states preceding the
failure stay visible.

The segment that failed is included when IDA kept its states up to the failure (this also
covers a failure in the very first segment). In a script, the figures are then served until
Ctrl+C, before the original error is reported. Best effort: it never raises itself.

Parameters
----------
simu : AlphaPEM
    Simulator whose last solve failed.
history : Union{Nothing, SimulationOutputs}
    Trajectory accumulated before the failure (`nothing` if the first segment failed).
ax1, ax2, ax3 : Any
    Axes returned by `figures_preparation`.

Returns
-------
Nothing
"""
function _display_live_failure!(simu::AlphaPEM, history, ax1, ax2, ax3)
    simu.cfg.type_display == :no_display && return nothing

    # Recovery of the failed segment, then display of the whole trajectory
    shown = false
    try
        sol = simu.sol
        # `sol` is the failed segment unless it is the previous one (nothing new was solved).
        if sol !== nothing && length(sol.t) > 1 &&
           (history === nothing || sol.t[1] ≈ history.solver.t[end])
            recovery!(simu)
            history = _accumulate_live_outputs!(simu, history)
        end
        if history !== nothing
            simu.outputs = history
            _refresh_live_display!(simu, ax1, ax2, ax3; force=true)
            shown = true
        end
    catch err
        @warn "Could not display the trajectory reached before the failure." exception=err maxlog=1
    end

    # Keep the figures reachable until the user has looked at them (scripts only)
    if shown && !isempty(_served_figure_routes) && _plot_block_requested()
        _serve_until_interrupted!("The simulation failed: the figures show the trajectory reached before the failure. Press Ctrl+C to stop serving and see the error.")
    end
    return nothing
end
