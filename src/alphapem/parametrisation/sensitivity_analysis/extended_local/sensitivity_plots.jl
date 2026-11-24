# -*- coding: utf-8 -*-

# ============================================================
# INTERPRETATION WARNING
# ============================================================

"""
Caution shown with every sensitivity figure.

A sensitivity magnitude is conditional on the model structure, its assumptions, the parameter bounds
and the operating cases. A large value does not by itself prove that a parameter is physically
important, and a small value does not prove that it is unimportant.
"""
const SENSITIVITY_INTERPRETATION_WARNING =
    "Sensitivity depends on the model, its assumptions and the parameter bounds; " *
    "it does not by itself show that a parameter is physically important."


# ============================================================
# SHARED PLOTTING HELPERS
# ============================================================

# Build a log-scale axis whose x ticks are the analysed parameter names.
function _parameter_axis(
    cell,
    parameter_names::Vector{Symbol};
    title::String,
    ylabel::String,
    xlabel::String = "",
    logarithmic::Bool = true,
    limits = (nothing, nothing),
)
    positions = collect(eachindex(parameter_names))
    return WGLMakie.Axis(
        cell;
        title,
        xlabel,
        ylabel,
        limits,
        yscale = logarithmic ? log10 : identity,
        xticks = (positions, string.(parameter_names)),
        xticklabelrotation = pi / 2,
    )
end


# Color of one cumulative state: increasingly opaque blue, with the final state in red.
function _state_color(state_index::Int, n_states::Int)
    state_index == n_states && return :red
    return (:blue, 0.18 + 0.67 * state_index / max(n_states - 1, 1))
end


# Keep only strictly positive finite values, as required by logarithmic axes.
_positive_or_nan(values) = [isfinite(value) && value > 0.0 ? Float64(value) : NaN
                            for value in values]


# Serve a figure in the browser page; return whether it was served.
function _serve_live_figure(figure; name::String, title::String, display::Bool)::Bool
    display || return false
    serve_interactive_figure!(figure; name, title)
    return true
end


# Keep a served figure alive when running as a script.
function _keep_figure_served(displayed::Bool, message::String)
    displayed || return nothing
    wait_for_interactive_figures!(message)
    return nothing
end


# ============================================================
# LIVE CONVERGENCE FIGURE
# ============================================================

"""
    SensitivityConvergencePlot

Mutable handles and progress state for the live three-panel convergence figure.
"""
mutable struct SensitivityConvergencePlot
    figure
    median_axis
    deviation_axis
    coverage_axis
    parameter_names::Vector{Symbol}
    sample_counts::Vector{Int}
    median_values::Vector{Any}
    deviation_values::Vector{Any}
    coverage_values::Vector{Any}
    next_state::Int
    displayed::Bool
end


"""
    create_sensitivity_convergence_plot(params, n_samples; display=true) -> SensitivityConvergencePlot

Create the live median/standard-deviation convergence figure with measurement coverage.

Each panel has one curve per cumulative sample prefix. Intermediate states are blue with increasing
opacity; the final state is red. The coverage panel shows the median fraction of available
measurements while partial curves contribute to the RMS. When `display` is true, the figure is
served in the AlphaPEM browser plot page before simulations start.

$(SENSITIVITY_INTERPRETATION_WARNING)

Parameters
----------
params : Vector{ParameterBound}
    Analysed parameters in the order used by the sensitivity matrices.
n_samples : Int
    Number of accepted samples expected from the screening stage.
display : Bool, optional
    Serve the figure in the browser as it is created. Default is `true`.

Returns
-------
SensitivityConvergencePlot
    Figure and state needed by [`update_sensitivity_convergence_plot!`](@ref).
"""
function create_sensitivity_convergence_plot(
    params::Vector{ParameterBound},
    n_samples::Int;
    display::Bool = true,
)::SensitivityConvergencePlot
    isempty(params) && throw(ArgumentError("params must not be empty"))
    sample_counts = convergence_sample_counts(n_samples)
    parameter_names = [parameter.name for parameter in params]
    positions = Float64.(eachindex(parameter_names))

    # The magnitude panels use a logarithmic y-axis, as the sensitivities span several orders
    # of magnitude.
    figure = WGLMakie.Figure(size = (1500, 1250))
    median_axis = _parameter_axis(
        figure[1, 1], parameter_names;
        title = "Convergence of voltage-sensitivity median",
        ylabel = "Median sensitivity magnitude")
    deviation_axis = _parameter_axis(
        figure[2, 1], parameter_names;
        title = "Convergence of voltage-sensitivity standard deviation",
        ylabel = "Sensitivity standard deviation")
    coverage_axis = _parameter_axis(
        figure[3, 1], parameter_names;
        title = "Convergence of median measurement coverage",
        xlabel = "Parameter", ylabel = "Median coverage",
        logarithmic = false, limits = (nothing, (0.0, 1.05)))
    WGLMakie.Label(figure[4, 1], SENSITIVITY_INTERPRETATION_WARNING;
                   fontsize = 12, halign = :left, tellwidth = false)

    # Reserve one observable line per state before opening the browser page. Updating these
    # observables later changes the existing WGL figure without adding plots from worker tasks.
    n_states = length(sample_counts)
    legend_handles = Any[]
    legend_labels = ["n = $count" for count in sample_counts]
    median_values = Any[]
    deviation_values = Any[]
    coverage_values = Any[]
    for state_index in 1:n_states
        color = _state_color(state_index, n_states)
        push!(legend_handles, WGLMakie.LineElement(color = color, linewidth = 2))
        median_state = WGLMakie.Observable(fill(NaN, length(params)))
        deviation_state = WGLMakie.Observable(fill(NaN, length(params)))
        coverage_state = WGLMakie.Observable(fill(0.0, length(params)))
        push!(median_values, median_state)
        push!(deviation_values, deviation_state)
        push!(coverage_values, coverage_state)
        for (axis, values) in ((median_axis, median_state),
                               (deviation_axis, deviation_state),
                               (coverage_axis, coverage_state))
            WGLMakie.lines!(axis, positions, values; color, linewidth = 2)
            WGLMakie.scatter!(axis, positions, values; color, markersize = 6)
        end
    end
    WGLMakie.Legend(figure[1:2, 2], legend_handles, legend_labels;
                    title = "Cumulative samples")

    displayed = _serve_live_figure(
        figure; name = "extended_local_convergence",
        title = "Extended local sensitivity convergence", display)
    return SensitivityConvergencePlot(
        figure, median_axis, deviation_axis, coverage_axis, parameter_names, sample_counts,
        median_values, deviation_values, coverage_values, 1, displayed)
end


"""
    update_sensitivity_convergence_plot!(plot, snapshot)

Update the pre-created median and standard-deviation curves for one completed state.

Parameters
----------
plot : SensitivityConvergencePlot
    Live figure returned by [`create_sensitivity_convergence_plot`](@ref).
snapshot : SensitivityConvergenceSnapshot
    Statistics for the next expected cumulative prefix.

Raises
------
ArgumentError
    If states arrive out of order or the parameter dimensions do not match the figure.
"""
function update_sensitivity_convergence_plot!(
    plot::SensitivityConvergencePlot,
    snapshot::SensitivityConvergenceSnapshot,
)
    snapshot.state_index == plot.next_state ||
        throw(ArgumentError("Expected convergence state $(plot.next_state), got $(snapshot.state_index)"))
    length(snapshot.median_magnitude) == length(plot.parameter_names) ||
        throw(ArgumentError("Convergence snapshot does not match the plotted parameter count"))
    length(snapshot.standard_deviation) == length(plot.parameter_names) ||
        throw(ArgumentError("Convergence snapshot standard deviations have the wrong length"))
    length(snapshot.median_coverage) == length(plot.parameter_names) ||
        throw(ArgumentError("Convergence snapshot coverage has the wrong length"))

    coverages = [isfinite(value) ? clamp(value, 0.0, 1.0) : 0.0
                 for value in snapshot.median_coverage]

    # Setting the state observables refreshes all panels in the already-served browser figure.
    plot.median_values[snapshot.state_index][] = _positive_or_nan(snapshot.median_magnitude)
    plot.deviation_values[snapshot.state_index][] = _positive_or_nan(snapshot.standard_deviation)
    plot.coverage_values[snapshot.state_index][] = coverages
    plot.next_state += 1
    return nothing
end


"""
    finish_sensitivity_convergence_plot!(plot)

Keep the live browser figures available after the simulation finishes when running as a script.

Parameters
----------
plot : SensitivityConvergencePlot
    Figure whose final convergence state has been drawn. Call this last: in a script it blocks
    while serving every live figure.

Returns
-------
Nothing
    In a script, waits until the user interrupts the plot server; in a REPL, returns immediately.
"""
function finish_sensitivity_convergence_plot!(plot::SensitivityConvergencePlot)
    _keep_figure_served(
        plot.displayed,
        "Sensitivity figures are complete. Press Ctrl+C to close the browser plot server.")
end


# ============================================================
# LIVE VOLTAGE-SENSITIVITY FIGURE
# ============================================================

"""
    VoltageSensitivityPlot

Mutable handles for the live voltage-sensitivity figure.

# Fields
- `figure`, `axis`: Makie figure and axis.
- `parameter_names::Vector{Symbol}`: Analysed parameters, one x position each.
- `points`: Observable vector of gray per-candidate points.
- `median_values`, `median_colors`: Observables for the median markers and their color values.
- `color_range`: Observable color limits, set when the medians are known.
- `displayed::Bool`: Whether the figure was served in the browser.
"""
mutable struct VoltageSensitivityPlot
    figure
    axis
    parameter_names::Vector{Symbol}
    points
    median_values
    median_colors
    color_range
    displayed::Bool
end


"""
    create_voltage_sensitivity_plot(params; display=true) -> VoltageSensitivityPlot

Create the live voltage-sensitivity figure.

Each accepted candidate contributes one gray point per parameter: the RMS of its signed voltage
sensitivities over all operating cases and current densities. "Signed" means that each underlying
sensitivity keeps the sign of the voltage slope with respect to the normalized parameter; the RMS
discards that sign to give a magnitude. Points are appended as candidates finish. The per-parameter
median is shown only once, at the end, as markers colored by a logarithmic color scale.

$(SENSITIVITY_INTERPRETATION_WARNING)

Parameters
----------
params : Vector{ParameterBound}
    Analysed parameters in the order used by the sensitivity matrices.
display : Bool, optional
    Serve the figure in the browser as it is created. Default is `true`.

Returns
-------
VoltageSensitivityPlot
    Figure and observables needed by [`add_voltage_sensitivity_candidate!`](@ref) and
    [`finish_voltage_sensitivity_plot!`](@ref).

Raises
------
ArgumentError
    If `params` is empty.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function create_voltage_sensitivity_plot(
    params::Vector{ParameterBound};
    display::Bool = true,
)::VoltageSensitivityPlot
    isempty(params) && throw(ArgumentError("params must not be empty"))
    parameter_names = [parameter.name for parameter in params]
    positions = Float64.(eachindex(parameter_names))

    figure = WGLMakie.Figure(size = (1500, 800))
    axis = _parameter_axis(
        figure[1, 1], parameter_names;
        title = "Voltage sensitivity",
        xlabel = "Parameter", ylabel = "Sensitivity magnitude")
    WGLMakie.Label(figure[2, 1:2], SENSITIVITY_INTERPRETATION_WARNING;
                   fontsize = 12, halign = :left, tellwidth = false)

    # Pre-create every plot so later updates only assign observables; the median layer stays
    # empty (NaN) until the whole run is complete.
    points = WGLMakie.Observable(WGLMakie.Point2f[])
    median_values = WGLMakie.Observable(fill(NaN, length(params)))
    median_colors = WGLMakie.Observable(fill(0.0, length(params)))
    color_range = WGLMakie.Observable((0.0, 1.0))
    WGLMakie.scatter!(axis, points; color = (:gray, 0.45), markersize = 6)
    median_plot = WGLMakie.scatter!(
        axis, positions, median_values;
        color = median_colors, colormap = :viridis, colorrange = color_range,
        markersize = 16, strokecolor = :black, strokewidth = 1)
    WGLMakie.Colorbar(figure[1, 2], median_plot; label = "log10 median sensitivity magnitude")

    displayed = _serve_live_figure(
        figure; name = "extended_local_voltage_sensitivity",
        title = "Extended local voltage sensitivity", display)
    return VoltageSensitivityPlot(
        figure, axis, parameter_names, points, median_values, median_colors, color_range,
        displayed)
end


"""
    add_voltage_sensitivity_candidate!(plot, magnitudes)

Add one completed candidate's gray points, one per parameter.

Parameters
----------
plot : VoltageSensitivityPlot
    Figure returned by [`create_voltage_sensitivity_plot`](@ref).
magnitudes : AbstractVector{<:Real}
    One RMS sensitivity magnitude per parameter. Non-positive or non-finite values are skipped
    because the y-axis is logarithmic.

Raises
------
ArgumentError
    If `magnitudes` does not have one entry per plotted parameter.
"""
function add_voltage_sensitivity_candidate!(
    plot::VoltageSensitivityPlot,
    magnitudes::AbstractVector{<:Real},
)
    length(magnitudes) == length(plot.parameter_names) ||
        throw(ArgumentError("magnitudes must contain one value per plotted parameter"))
    new_points = [WGLMakie.Point2f(index, value)
                  for (index, value) in enumerate(magnitudes) if isfinite(value) && value > 0]
    isempty(new_points) && return nothing
    append!(plot.points[], new_points)
    WGLMakie.notify(plot.points)
    return nothing
end


"""
    finish_voltage_sensitivity_plot!(plot, median_magnitudes)

Draw the final per-parameter medians on the color scale.

Parameters
----------
plot : VoltageSensitivityPlot
    Figure whose candidate points have all been added.
median_magnitudes : AbstractVector{<:Real}
    Median RMS sensitivity magnitude per parameter; non-positive or non-finite values are hidden.

Raises
------
ArgumentError
    If `median_magnitudes` does not have one entry per plotted parameter.
"""
function finish_voltage_sensitivity_plot!(
    plot::VoltageSensitivityPlot,
    median_magnitudes::AbstractVector{<:Real},
)
    length(median_magnitudes) == length(plot.parameter_names) ||
        throw(ArgumentError("median_magnitudes must contain one value per plotted parameter"))
    medians = _positive_or_nan(median_magnitudes)
    colors = log10.(medians)
    finite_colors = filter(isfinite, colors)

    # Equal limits would make the color scale degenerate, so widen a single-valued range.
    lower, upper = isempty(finite_colors) ? (0.0, 1.0) : extrema(finite_colors)
    lower == upper && (lower -= 0.5; upper += 0.5)
    plot.color_range[] = (lower, upper)
    plot.median_colors[] = [isfinite(value) ? value : lower for value in colors]
    plot.median_values[] = medians
    return nothing
end
