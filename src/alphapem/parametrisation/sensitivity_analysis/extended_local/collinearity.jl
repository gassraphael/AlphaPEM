# -*- coding: utf-8 -*-

# ============================================================
# COLLINEARITY OF VOLTAGE-SENSITIVITY VECTORS
# ============================================================

"""
    candidate_collinearity(sensitivities, params, candidate_index; floor=1e-12) -> Matrix{Float64}

Compute the pairwise collinearity of one candidate's voltage-sensitivity vectors.

The vector of a parameter concatenates its sensitivities over all operating cases and current
densities. Each pair uses only the measurements where both sensitivities are finite, so a partial
high-current curve still contributes and missing values are never treated as zero. Pairs without
a common measurement are `NaN`.

Parameters
----------
sensitivities : AbstractMatrix
    Candidate-by-case matrix of `n_params × n_current_densities` blocks, or `nothing` if unfinished.
params : Vector{ParameterBound}
    Analysed parameters in sensitivity-matrix row order.
candidate_index : Int
    Row of `sensitivities` to analyse.
floor : Real, optional
    Minimum denominator, as in [`collinearity_index`](@ref). Default is `1e-12`.

Returns
-------
Matrix{Float64}
    Symmetric `n_params × n_params` matrix with values in `[0, 1]`.

Raises
------
ArgumentError
    If a completed block has the wrong number of parameter rows.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function candidate_collinearity(
    sensitivities::AbstractMatrix,
    params::Vector{ParameterBound},
    candidate_index::Int;
    floor::Real = 1e-12,
)::Matrix{Float64}
    n_params = length(params)

    # Concatenate every completed block into one parameter-by-measurement matrix.
    blocks = Matrix{Float64}[]
    for case_index in axes(sensitivities, 2)
        block = sensitivities[candidate_index, case_index]
        block === nothing && continue
        size(block, 1) == n_params ||
            throw(ArgumentError("Sensitivity blocks have inconsistent parameter dimensions"))
        push!(blocks, block)
    end
    S = isempty(blocks) ? zeros(n_params, 0) : reduce(hcat, blocks)

    # Pairwise computation on shared finite measurements.
    psi = fill(NaN, n_params, n_params)
    for i in 1:n_params, j in i:n_params
        shared = [isfinite(S[i, k]) && isfinite(S[j, k]) for k in axes(S, 2)]
        any(shared) || continue
        si = @view S[i, shared]
        sj = @view S[j, shared]
        value = abs(sum(si .* sj)) / max(sqrt(sum(abs2, si)) * sqrt(sum(abs2, sj)), floor)
        psi[i, j] = psi[j, i] = min(value, 1.0)
    end
    return psi
end


"""
    CollinearitySummary

Median collinearity and sensitivity magnitudes over the accepted candidates.

# Fields
- `median_collinearity::Matrix{Float64}`: Median pairwise index over candidates.
- `maximum_collinearity::Matrix{Float64}`: Maximum pairwise index over candidates.
- `sensitivity_median::Vector{Float64}`: Median RMS voltage-sensitivity magnitude per parameter.
"""
struct CollinearitySummary
    median_collinearity::Matrix{Float64}
    maximum_collinearity::Matrix{Float64}
    sensitivity_median::Vector{Float64}
end


"""
    collinearity_summary(result; floor=result.config.collinearity_floor) -> CollinearitySummary

Summarize collinearity and sensitivity magnitude over every candidate of a simulation result.

Parameters
----------
result : ExtendedLocalSimulationResult
    Completed sensitivity simulations.
floor : Real, optional
    Minimum denominator of the collinearity index.

Returns
-------
CollinearitySummary
    Candidate medians and maxima; `NaN` where no candidate has a valid value.
"""
function collinearity_summary(
    result::ExtendedLocalSimulationResult;
    floor::Real = result.config.collinearity_floor,
)::CollinearitySummary
    params = result.params
    n_params = length(params)
    n_candidates = size(result.sensitivities, 1)
    psi = [candidate_collinearity(result.sensitivities, params, index; floor)
           for index in 1:n_candidates]
    magnitudes = [candidate_rms_magnitudes(result.sensitivities, params, index)[1]
                  for index in 1:n_candidates]

    # Aggregate finite candidate values independently for each pair or parameter.
    summarize(values, reducer) = isempty(values) ? NaN : reducer(values)
    median_psi = [summarize(filter(isfinite, [matrix[i, j] for matrix in psi]), median)
                  for i in 1:n_params, j in 1:n_params]
    maximum_psi = [summarize(filter(isfinite, [matrix[i, j] for matrix in psi]), maximum)
                   for i in 1:n_params, j in 1:n_params]
    sensitivity_median = [summarize(filter(isfinite, [values[i] for values in magnitudes]), median)
                          for i in 1:n_params]
    return CollinearitySummary(median_psi, maximum_psi, sensitivity_median)
end


"""
    collinear_pairs(summary, params; collinearity_threshold=0.707, sensitivity_threshold=0.1)
        -> DataFrame

List collinear parameter pairs for which both parameters have a non-negligible sensitivity.

Parameters
----------
summary : CollinearitySummary
    Result of [`collinearity_summary`](@ref).
params : Vector{ParameterBound}
    Analysed parameters in summary order.
collinearity_threshold : Real, optional
    Minimum median collinearity. Default is `0.707` (angle below 45°).
sensitivity_threshold : Real, optional
    Minimum median sensitivity magnitude of both parameters of the pair. Default is `0.1`.

Returns
-------
DataFrame
    Columns `i`, `j`, `psi_median`, `psi_max`, `s_i`, `s_j`, sorted by decreasing `psi_median`.
"""
function collinear_pairs(
    summary::CollinearitySummary,
    params::Vector{ParameterBound};
    collinearity_threshold::Real = 0.707,
    sensitivity_threshold::Real = 0.1,
)::DataFrame
    n_params = length(params)
    size(summary.median_collinearity) == (n_params, n_params) ||
        throw(ArgumentError("summary does not match params"))
    rows = NamedTuple[]
    for i in 1:n_params, j in (i + 1):n_params
        psi = summary.median_collinearity[i, j]
        s_i, s_j = summary.sensitivity_median[i], summary.sensitivity_median[j]
        isfinite(psi) && psi >= collinearity_threshold || continue
        isfinite(s_i) && s_i >= sensitivity_threshold &&
            isfinite(s_j) && s_j >= sensitivity_threshold || continue
        push!(rows, (i = params[i].name, j = params[j].name, psi_median = psi,
                     psi_max = summary.maximum_collinearity[i, j], s_i = s_i, s_j = s_j))
    end
    pairs = DataFrame(
        i = Symbol[row.i for row in rows], j = Symbol[row.j for row in rows],
        psi_median = Float64[row.psi_median for row in rows],
        psi_max = Float64[row.psi_max for row in rows],
        s_i = Float64[row.s_i for row in rows], s_j = Float64[row.s_j for row in rows])
    return sort(pairs, :psi_median; rev = true)
end
