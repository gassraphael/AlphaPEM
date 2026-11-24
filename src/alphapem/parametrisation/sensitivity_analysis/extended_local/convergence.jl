# -*- coding: utf-8 -*-

# ============================================================
# CONVERGENCE SAMPLE COUNTS AND SUMMARY DATA
# ============================================================

"""
    convergence_sample_counts(n_samples; n_states=6) -> Vector{Int}

Choose cumulative sample counts for the visual convergence comparison.

The counts divide the available sample range into approximately equal intervals, are rounded to
distinct integers, and always include the final sample count. If fewer than `n_states` samples are
available, every sample forms one state.

Parameters
----------
n_samples : Int
    Number of accepted parameter-space samples.
n_states : Int, optional
    Maximum number of plotted states, including the final state. Default is `6`.

Returns
-------
Vector{Int}
    Strictly increasing cumulative sample counts ending at `n_samples`.

Raises
------
ArgumentError
    If either input is not positive.
"""
function convergence_sample_counts(n_samples::Int; n_states::Int = 6)::Vector{Int}
    n_samples > 0 || throw(ArgumentError("n_samples must be positive"))
    n_states > 0 || throw(ArgumentError("n_states must be positive"))

    # Use no more states than distinct prefixes can exist.
    n_states = min(n_states, n_samples)
    counts = [round(Int, state * n_samples / n_states) for state in 1:n_states]
    counts[end] = n_samples
    return unique!(counts)
end


"""
    SensitivityConvergenceSnapshot

Statistics for one completed cumulative prefix of parameter-space samples.

# Fields
- `state_index::Int`: One-based position of this convergence state.
- `sample_count::Int`: Number of accepted samples included in this prefix.
- `median_magnitude::Vector{Float64}`: Per-parameter median RMS sensitivity over available measurements.
- `standard_deviation::Vector{Float64}`: Sample standard deviation of per-candidate RMS sensitivities.
- `valid_sample_counts::Vector{Int}`: Number of prefix samples with at least one valid measurement.
- `median_coverage::Vector{Float64}`: Median fraction of expected measurements available per parameter.
"""
struct SensitivityConvergenceSnapshot
    state_index::Int
    sample_count::Int
    median_magnitude::Vector{Float64}
    standard_deviation::Vector{Float64}
    valid_sample_counts::Vector{Int}
    median_coverage::Vector{Float64}
end


"""
    convergence_snapshot(sensitivities, params, n_samples, state_index)
        -> SensitivityConvergenceSnapshot

Compute median and standard deviation of RMS voltage sensitivities for one cumulative sample prefix.

For each candidate and parameter, take the root-mean-square of all finite sensitivities across the
available current-density points and operating cases. A candidate with a partial curve still
contributes its valid measurements; missing values are never treated as zero. Per-parameter median
coverage reports how much of the expected measurement vector was available across the prefix.

Parameters
----------
sensitivities : AbstractMatrix
    Candidate-by-case matrix whose completed entries are sensitivity matrices with shape
    `n_params × n_current_densities`; unfinished entries are `nothing`.
params : Vector{ParameterBound}
    Analysed parameters in sensitivity-matrix row order.
n_samples : Int
    Number of initial candidate rows included in this cumulative state.
state_index : Int
    One-based state position, retained in the returned snapshot.

Returns
-------
SensitivityConvergenceSnapshot
    Per-parameter median RMS, standard deviation, number of contributing samples and median
    measurement coverage.

Raises
------
ArgumentError
    If the prefix is empty, exceeds the available candidates, or completed blocks have
    inconsistent shapes.
"""
function convergence_snapshot(
    sensitivities::AbstractMatrix,
    params::Vector{ParameterBound},
    n_samples::Int,
    state_index::Int,
)::SensitivityConvergenceSnapshot
    1 <= n_samples <= size(sensitivities, 1) ||
        throw(ArgumentError("n_samples must select a non-empty prefix of candidate rows"))
    state_index >= 1 || throw(ArgumentError("state_index must be positive"))

    # Parameter metadata defines the expected row count and total measurement width.
    n_params = length(params)
    n_cases = size(sensitivities, 2)
    n_current_densities = 0
    for candidate_index in 1:n_samples, case_index in 1:n_cases
        block = sensitivities[candidate_index, case_index]
        block === nothing && continue
        size(block, 1) == n_params ||
            throw(ArgumentError("Sensitivity blocks have inconsistent parameter dimensions"))
        n_current_densities = size(block, 2)
        break
    end
    expected_measurements = n_cases * n_current_densities

    # Store each candidate's RMS over its available measurements and its corresponding coverage.
    sample_magnitudes = fill(NaN, n_samples, n_params)
    sample_coverages = zeros(n_samples, n_params)
    for candidate_index in 1:n_samples, parameter_index in 1:n_params
        measurement_values = Float64[]

        # Collect finite entries from every case; partial high-current curves remain informative.
        for case_index in axes(sensitivities, 2)
            block = sensitivities[candidate_index, case_index]
            block === nothing && continue
            size(block, 1) == n_params && size(block, 2) == n_current_densities ||
                throw(ArgumentError("Sensitivity blocks have inconsistent dimensions"))
            parameter_values = @view block[parameter_index, :]
            append!(measurement_values, filter(isfinite, parameter_values))
        end

        sample_coverages[candidate_index, parameter_index] =
            expected_measurements == 0 ? 0.0 : length(measurement_values) / expected_measurements
        isempty(measurement_values) && continue
        sample_magnitudes[candidate_index, parameter_index] =
            sqrt(sum(abs2, measurement_values) / length(measurement_values))
    end

    # Summarize RMS values from candidates with at least one valid measurement.
    median_magnitude = fill(NaN, n_params)
    standard_deviation = fill(NaN, n_params)
    valid_sample_counts = zeros(Int, n_params)
    median_coverage = zeros(n_params)
    for parameter_index in 1:n_params
        valid_values = filter(isfinite, @view sample_magnitudes[:, parameter_index])
        valid_sample_counts[parameter_index] = length(valid_values)
        median_coverage[parameter_index] =
            median(@view sample_coverages[:, parameter_index])
        if !isempty(valid_values)
            median_magnitude[parameter_index] = median(valid_values)
            length(valid_values) >= 2 &&
                (standard_deviation[parameter_index] = std(valid_values))
        end
    end

    return SensitivityConvergenceSnapshot(
        state_index, n_samples, median_magnitude, standard_deviation,
        valid_sample_counts, median_coverage)
end


# ============================================================
# LIVE CONVERGENCE FIGURE
# ============================================================

"""
    SensitivityConvergencePlot

Mutable handles and progress state for the live two-panel convergence figure.
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
    positions = collect(eachindex(parameter_names))
    labels = string.(parameter_names)

    # Both panels use a logarithmic y-axis, as the article's convergence figure spans
    # several orders of magnitude in sensitivity.
    figure = WGLMakie.Figure(size = (1500, 1250))
    median_axis = WGLMakie.Axis(
        figure[1, 1],
        title = "Convergence of voltage-sensitivity median",
        ylabel = "Median sensitivity magnitude",
        yscale = log10,
        xticks = (positions, labels),
        xticklabelrotation = pi / 2,
    )
    deviation_axis = WGLMakie.Axis(
        figure[2, 1],
        title = "Convergence of voltage-sensitivity standard deviation",
        ylabel = "Sensitivity standard deviation",
        yscale = log10,
        xticks = (positions, labels),
        xticklabelrotation = pi / 2,
    )
    coverage_axis = WGLMakie.Axis(
        figure[3, 1],
        title = "Convergence of median measurement coverage",
        xlabel = "Parameter",
        ylabel = "Median coverage",
        limits = (nothing, (0.0, 1.05)),
        xticks = (positions, labels),
        xticklabelrotation = pi / 2,
    )

    # Reserve one observable line per state before opening the browser page. Updating these
    # observables later changes the existing WGL figure without adding plots from worker tasks.
    n_states = length(sample_counts)
    legend_handles = Any[]
    legend_labels = ["n = $count" for count in sample_counts]
    median_values = Any[]
    deviation_values = Any[]
    coverage_values = Any[]
    positions_float = Float64.(positions)
    for state_index in 1:n_states
        color = if state_index == n_states
            :red
        else
            alpha = 0.18 + 0.67 * state_index / max(n_states - 1, 1)
            (:blue, alpha)
        end
        push!(legend_handles, WGLMakie.LineElement(color = color, linewidth = 2))
        median_state = WGLMakie.Observable(fill(NaN, length(params)))
        deviation_state = WGLMakie.Observable(fill(NaN, length(params)))
        coverage_state = WGLMakie.Observable(fill(0.0, length(params)))
        push!(median_values, median_state)
        push!(deviation_values, deviation_state)
        push!(coverage_values, coverage_state)
        WGLMakie.lines!(median_axis, positions_float, median_state; color, linewidth = 2)
        WGLMakie.scatter!(median_axis, positions_float, median_state; color, markersize = 6)
        WGLMakie.lines!(deviation_axis, positions_float, deviation_state; color, linewidth = 2)
        WGLMakie.scatter!(deviation_axis, positions_float, deviation_state; color, markersize = 6)
        WGLMakie.lines!(coverage_axis, positions_float, coverage_state; color, linewidth = 2)
        WGLMakie.scatter!(coverage_axis, positions_float, coverage_state; color, markersize = 6)
    end
    WGLMakie.Legend(figure[1:2, 2], legend_handles, legend_labels;
                    title = "Cumulative samples")

    plot = SensitivityConvergencePlot(
        figure, median_axis, deviation_axis, coverage_axis, parameter_names, sample_counts,
        median_values, deviation_values, coverage_values, 1, false)
    if display
        serve_interactive_figure!(
            figure;
            name = "extended_local_convergence",
            title = "Extended local sensitivity convergence",
        )
        plot.displayed = true
    end
    return plot
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

    # Keep only positive finite values because both convergence axes are logarithmic.
    medians = [isfinite(value) && value > 0.0 ? value : NaN
               for value in snapshot.median_magnitude]
    deviations = [isfinite(value) && value > 0.0 ? value : NaN
                  for value in snapshot.standard_deviation]
    coverages = [isfinite(value) ? clamp(value, 0.0, 1.0) : 0.0
                 for value in snapshot.median_coverage]

    # Setting the state observables refreshes all panels in the already-served browser figure.
    plot.median_values[snapshot.state_index][] = medians
    plot.deviation_values[snapshot.state_index][] = deviations
    plot.coverage_values[snapshot.state_index][] = coverages
    plot.next_state += 1
    return nothing
end


"""
    finish_sensitivity_convergence_plot!(plot)

Keep the live browser figure available after the simulation finishes when running as a script.

Parameters
----------
plot : SensitivityConvergencePlot
    Figure whose final convergence state has been drawn.

Returns
-------
Nothing
    In a script, waits until the user interrupts the plot server; in a REPL, returns immediately.
"""
function finish_sensitivity_convergence_plot!(plot::SensitivityConvergencePlot)
    plot.displayed || return nothing
    wait_for_interactive_figures!(
        "Convergence figure is complete. Press Ctrl+C to close the browser plot server.")
    return nothing
end
