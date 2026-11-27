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
    candidate_rms_magnitudes(sensitivities, params, candidate_index; n_current_densities=nothing)
        -> Tuple{Vector{Float64}, Vector{Float64}}

Summarize one candidate's voltage sensitivities as one RMS magnitude per parameter.

The RMS is taken over every finite sensitivity of the parameter across all operating cases and
current densities, so partially completed high-current curves still contribute. Missing values are
never treated as zero. Coverage is the fraction of expected measurements that were finite.

Parameters
----------
sensitivities : AbstractMatrix
    Candidate-by-case matrix of `n_params × n_current_densities` blocks, or `nothing` if unfinished.
params : Vector{ParameterBound}
    Analysed parameters in sensitivity-matrix row order.
candidate_index : Int
    Row of `sensitivities` to summarize.
n_current_densities : Int or nothing, optional
    Expected block width. If `nothing`, it is read from the candidate's first completed block.

Returns
-------
Tuple{Vector{Float64}, Vector{Float64}}
    RMS magnitudes (`NaN` without finite measurements) and measurement coverages.

Raises
------
ArgumentError
    If a completed block has inconsistent dimensions.
"""
function candidate_rms_magnitudes(
    sensitivities::AbstractMatrix,
    params::Vector{ParameterBound},
    candidate_index::Int;
    n_current_densities::Union{Nothing, Int} = nothing,
)::Tuple{Vector{Float64}, Vector{Float64}}
    n_params = length(params)
    n_cases = size(sensitivities, 2)
    if n_current_densities === nothing
        first_block = findfirst(!isnothing, @view sensitivities[candidate_index, :])
        n_current_densities = first_block === nothing ? 0 :
            size(sensitivities[candidate_index, first_block], 2)
    end
    expected_measurements = n_cases * n_current_densities

    magnitudes = fill(NaN, n_params)
    coverages = zeros(n_params)
    for parameter_index in 1:n_params
        measurement_values = Float64[]

        # Collect finite entries from every case; partial high-current curves remain informative.
        for case_index in 1:n_cases
            block = sensitivities[candidate_index, case_index]
            block === nothing && continue
            size(block, 1) == n_params && size(block, 2) == n_current_densities ||
                throw(ArgumentError("Sensitivity blocks have inconsistent dimensions"))
            append!(measurement_values, filter(isfinite, @view block[parameter_index, :]))
        end

        coverages[parameter_index] =
            expected_measurements == 0 ? 0.0 : length(measurement_values) / expected_measurements
        isempty(measurement_values) && continue
        magnitudes[parameter_index] =
            sqrt(sum(abs2, measurement_values) / length(measurement_values))
    end
    return magnitudes, coverages
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
    for candidate_index in 1:n_samples
        sample_magnitudes[candidate_index, :], sample_coverages[candidate_index, :] =
            candidate_rms_magnitudes(
                sensitivities, params, candidate_index; n_current_densities)
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
