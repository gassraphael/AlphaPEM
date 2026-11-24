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


# ============================================================
# FIGURE AUTOSAVE
# ============================================================

"""
    FigureAutosave(directory, name)

Overwritable image file that mirrors a live figure while the analysis runs.

Every save writes the same PDF file, so the file always shows the latest drawn state and ends as the
final figure. A save is attempted each time the live plot is updated.

Parameters
----------
directory : AbstractString
    Output directory, created on the first save.
name : AbstractString
    File name without extension.
"""
mutable struct FigureAutosave
    path::String

    function FigureAutosave(directory::AbstractString, name::AbstractString)
        return new(joinpath(directory, "$(name).pdf"))
    end
end


"""
    autosave_figure!(autosave, figure) -> Bool

Overwrite the autosave file with the current state of `figure`.

A rendering or disk failure only logs a warning, so it cannot interrupt a long analysis. `nothing`
as `autosave` does nothing.

Parameters
----------
autosave : FigureAutosave or nothing
    Target PDF file.
figure : Any
    Makie figure to render.

Returns
-------
Bool
    Whether the file was written.
"""
function autosave_figure!(::Nothing, figure)::Bool
    return false
end

function autosave_figure!(autosave::FigureAutosave, figure)::Bool
    # Render to a temporary file first so a crash never leaves a truncated image.
    try
        mkpath(dirname(autosave.path))
        temporary_path = autosave.path * ".tmp"
        write(temporary_path, export_figure_bytes(figure, "pdf"))
        mv(temporary_path, autosave.path; force = true)
        return true
    catch err
        err isa InterruptException && rethrow()
        @warn "Could not save sensitivity figure" path = autosave.path exception = (err, catch_backtrace())
        return false
    end
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
    sensitivity_threshold::Float64
    autosave::Union{Nothing, FigureAutosave}
    next_state::Int
end


"""
    create_sensitivity_convergence_plot(params, n_samples; sensitivity_threshold=0.1,
                                        autosave=nothing) -> SensitivityConvergencePlot

Create the live median/standard-deviation convergence figure with measurement coverage.

Each panel has one curve per cumulative sample prefix. Intermediate states are blue with increasing
opacity; the final state is red. The coverage panel shows the median fraction of available
measurements while partial curves contribute to the RMS. The figure is
served with [`serve_sensitivity_plots!`](@ref) before simulations start.

$(SENSITIVITY_INTERPRETATION_WARNING)

Parameters
----------
params : Vector{ParameterBound}
    Analysed parameters in the order used by the sensitivity matrices.
n_samples : Int
    Number of accepted samples expected from the screening stage.
sensitivity_threshold : Real, optional
    Horizontal reference line for non-negligible median sensitivity. Default is `0.1`.
autosave : FigureAutosave or nothing, optional
    File overwritten after every updated state; the last state gives the final figure.

Returns
-------
SensitivityConvergencePlot
    Figure and state needed by [`update_sensitivity_convergence_plot!`](@ref).
"""
function create_sensitivity_convergence_plot(
    params::Vector{ParameterBound},
    n_samples::Int,
    ;
    sensitivity_threshold::Real = 0.1,
    autosave::Union{Nothing, FigureAutosave} = nothing,
)::SensitivityConvergencePlot
    isempty(params) && throw(ArgumentError("params must not be empty"))
    isfinite(sensitivity_threshold) && sensitivity_threshold >= 0 ||
        throw(ArgumentError("sensitivity_threshold must be finite and non-negative"))
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
    # A zero threshold cannot be represented on this logarithmic axis.
    sensitivity_threshold > 0 &&
        WGLMakie.hlines!(median_axis, [sensitivity_threshold];
                         color = :black, linestyle = :dot, linewidth = 1.5)
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

    return SensitivityConvergencePlot(
        figure, median_axis, deviation_axis, coverage_axis, parameter_names, sample_counts,
        median_values, deviation_values, coverage_values, Float64(sensitivity_threshold), autosave, 1)
end


"""
    update_sensitivity_convergence_plot!(plot, snapshot)

Update the pre-created median and standard-deviation curves for one completed state, then
refresh the autosave file.

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

    # Overwrite the saved PDF after every new convergence curve.
    autosave_figure!(plot.autosave, plot.figure)
    return nothing
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
"""
mutable struct VoltageSensitivityPlot
    figure
    axis
    parameter_names::Vector{Symbol}
    points
    median_values
    median_colors
    color_range
    sensitivity_threshold::Float64
    autosave::Union{Nothing, FigureAutosave}
end


"""
    create_voltage_sensitivity_plot(params; sensitivity_threshold=0.1, autosave=nothing)
        -> VoltageSensitivityPlot

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
sensitivity_threshold : Real, optional
    Horizontal reference line for non-negligible median sensitivity. Default is `0.1`.
autosave : FigureAutosave or nothing, optional
    File overwritten as candidate points are added and when the medians are drawn.

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
    params::Vector{ParameterBound},
    ;
    sensitivity_threshold::Real = 0.1,
    autosave::Union{Nothing, FigureAutosave} = nothing,
)::VoltageSensitivityPlot
    isempty(params) && throw(ArgumentError("params must not be empty"))
    isfinite(sensitivity_threshold) && sensitivity_threshold >= 0 ||
        throw(ArgumentError("sensitivity_threshold must be finite and non-negative"))
    parameter_names = [parameter.name for parameter in params]
    positions = Float64.(eachindex(parameter_names))

    figure = WGLMakie.Figure(size = (1500, 800))
    axis = _parameter_axis(
        figure[1, 1], parameter_names;
        title = "Voltage sensitivity",
        xlabel = "Parameter", ylabel = "Sensitivity magnitude")
    # A zero threshold cannot be represented on this logarithmic axis.
    sensitivity_threshold > 0 &&
        WGLMakie.hlines!(axis, [sensitivity_threshold];
                         color = :black, linestyle = :dot, linewidth = 1.5)
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

    return VoltageSensitivityPlot(
        figure, axis, parameter_names, points, median_values, median_colors, color_range,
        Float64(sensitivity_threshold), autosave)
end


"""
    add_voltage_sensitivity_candidate!(plot, magnitudes)

Add one completed candidate's gray points, one per parameter, and refresh the autosave file
when its rate limit allows.

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

    # Overwrite the saved PDF after this candidate's batch of parameter points is plotted.
    autosave_figure!(plot.autosave, plot.figure)
    return nothing
end


"""
    finish_voltage_sensitivity_plot!(plot, median_magnitudes)

Draw the final per-parameter medians on the color scale and write the final autosave file.

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

    # The final median update overwrites the last candidate-only version.
    autosave_figure!(plot.autosave, plot.figure)
    return nothing
end


# ============================================================
# COLLINEARITY FIGURE
# ============================================================

"""
Caution shown with every collinearity figure.

High collinearity is a statement about the model structure and the measurements used, not about
the physics of the parameters.
"""
const COLLINEARITY_INTERPRETATION_WARNING =
    "High collinearity does not mean the parameters are physically related: the model structure " *
    "cannot distinguish their effects from these measurements alone, so the same curve can be " *
    "obtained by tuning either. If experiments clearly separate them, the model needs refinement."


"""
    CollinearityPlot

Mutable handles for the voltage-sensitivity collinearity figure.

# Fields
- `figure`, `axis`: Makie figure and pair axis.
- `parameter_names::Vector{Symbol}`: Analysed parameters, in matrix order.
- `faded_values`, `highlighted_values`: Observables for all pairs (faded) and for pairs above the
  collinearity threshold (full color).
"""
mutable struct CollinearityPlot
    figure
    axis
    parameter_names::Vector{Symbol}
    faded_values
    highlighted_values
    autosave::Union{Nothing, FigureAutosave}
end


"""
    create_collinearity_plot(params; autosave=nothing) -> CollinearityPlot

Create the empty voltage-sensitivity collinearity figure, to be filled by
[`finish_collinearity_plot!`](@ref).

The layout follows the article's collinearity figures: lower triangle of pairwise collinearity
indices on a fixed `[0, 1]` jet color scale, pairs above the threshold in full color and the others
faded, with a color bar from 0 (blue) to 1 (red).

$(COLLINEARITY_INTERPRETATION_WARNING)

Parameters
----------
params : Vector{ParameterBound}
    Analysed parameters in the order used by the sensitivity matrices.
autosave : FigureAutosave or nothing, optional
    File written when the collinearity matrix is drawn.

Returns
-------
CollinearityPlot
    Figure and observables updated at the end of the analysis.

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
function create_collinearity_plot(
    params::Vector{ParameterBound};
    autosave::Union{Nothing, FigureAutosave} = nothing,
)::CollinearityPlot
    isempty(params) && throw(ArgumentError("params must not be empty"))
    parameter_names = [parameter.name for parameter in params]
    n_params = length(parameter_names)
    positions = collect(1:n_params)
    labels = string.(parameter_names)

    figure = WGLMakie.Figure(size = (1100, 1250))
    axis = WGLMakie.Axis(
        figure[1, 1];
        title = "Voltage Sensitivity Collinearity",
        aspect = WGLMakie.DataAspect(),
        yreversed = true,
        xticks = (positions, labels), yticks = (positions, labels),
        xticklabelrotation = pi / 2,
        xticklabelsize = 9, yticklabelsize = 9,
    )

    # Both layers use the same colormap; the faded layer sits below the highlighted one.
    faded_values = WGLMakie.Observable(fill(NaN, n_params, n_params))
    highlighted_values = WGLMakie.Observable(fill(NaN, n_params, n_params))
    WGLMakie.heatmap!(axis, positions, positions, faded_values;
                      colormap = :jet, colorrange = (0, 1), alpha = 0.25, nan_color = :transparent)
    highlight = WGLMakie.heatmap!(axis, positions, positions, highlighted_values;
                                  colormap = :jet, colorrange = (0, 1), nan_color = :transparent)
    WGLMakie.Colorbar(figure[1, 2], highlight; label = "Collinearity index")
    WGLMakie.Label(figure[2, 1:2], COLLINEARITY_INTERPRETATION_WARNING;
                   fontsize = 12, halign = :left, tellwidth = false, word_wrap = true)

    return CollinearityPlot(
        figure, axis, parameter_names, faded_values, highlighted_values, autosave)
end


"""
    finish_collinearity_plot!(plot, collinearity, sensitivity_median; collinearity_threshold,
                              sensitivity_threshold)

Draw the median pairwise collinearity and write the final autosave file.

Only pairs where both parameters have a non-negligible median sensitivity are drawn, and only pairs above `collinearity_threshold` use the full color.

Parameters
----------
plot : CollinearityPlot
    Figure returned by [`create_collinearity_plot`](@ref).
collinearity : AbstractMatrix{<:Real}
    Symmetric `n_params × n_params` median collinearity indices; `NaN` marks unavailable pairs.
sensitivity_median : AbstractVector{<:Real}
    Median voltage-sensitivity magnitude per parameter, used for the non-negligible criterion.
collinearity_threshold : Real, optional
    Index above which a pair is highlighted. Default is `0.707` (angle below 45°).
sensitivity_threshold : Real, optional
    Minimum median magnitude for a parameter to be non-negligible. Default is `0.1`.

Raises
------
ArgumentError
    If the input sizes do not match the plotted parameters.
"""
function finish_collinearity_plot!(
    plot::CollinearityPlot,
    collinearity::AbstractMatrix{<:Real},
    sensitivity_median::AbstractVector{<:Real};
    collinearity_threshold::Real = 0.707,
    sensitivity_threshold::Real = 0.1,
)
    n_params = length(plot.parameter_names)
    size(collinearity) == (n_params, n_params) ||
        throw(ArgumentError("collinearity must be n_params × n_params"))
    length(sensitivity_median) == n_params ||
        throw(ArgumentError("sensitivity_median must contain one value per parameter"))

    # Lower triangle only; the heatmap is indexed (column parameter, row parameter).
    influential = [isfinite(value) && value >= sensitivity_threshold
                   for value in sensitivity_median]
    faded = fill(NaN, n_params, n_params)
    highlighted = fill(NaN, n_params, n_params)
    for row in 1:n_params, column in 1:(row - 1)
        value = collinearity[row, column]
        (isfinite(value) && influential[row] && influential[column]) || continue
        faded[column, row] = value
        value >= collinearity_threshold && (highlighted[column, row] = value)
    end

    plot.faded_values[] = faded
    plot.highlighted_values[] = highlighted
    autosave_figure!(plot.autosave, plot.figure)
    return nothing
end


# ============================================================
# BROWSER SERVICE
# ============================================================

"""
    serve_sensitivity_plots!(plots...) -> String

Serve every given sensitivity figure together in the AlphaPEM browser plot page.

Call it once before the simulations start so the figures update as results arrive.

Parameters
----------
plots : SensitivityConvergencePlot, VoltageSensitivityPlot, CollinearityPlot
    Figures to publish; each type has a fixed route name and title.

Returns
-------
String
    Entry URL of the plot page.
"""
function serve_sensitivity_plots!(plots...)::String
    descriptions = Dict(
        SensitivityConvergencePlot =>
            ("extended_local_convergence", "Extended local sensitivity convergence"),
        VoltageSensitivityPlot =>
            ("extended_local_voltage_sensitivity", "Extended local voltage sensitivity"),
        CollinearityPlot =>
            ("extended_local_collinearity", "Extended local voltage collinearity"),
    )
    return serve_interactive_figures!(
        [plot.figure for plot in plots];
        names = [descriptions[typeof(plot)][1] for plot in plots],
        titles = [descriptions[typeof(plot)][2] for plot in plots])
end
