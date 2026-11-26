# -*- coding: utf-8 -*-

# ============================================================
# PARAMETER-SUBSET SELECTION: COLLINEARITY FILTER
# ============================================================

"""
    CollinearityFilterResult

Outcome of [`filter_collinear_parameters`](@ref).

# Fields
- `retained::Vector{Symbol}`: Names of the parameters kept in the candidate list, in `params` order.
- `removed::DataFrame`: One row per removed parameter, with columns `removed`, `kept` (the partner
  it was collinear with), `psi_median`, `s_removed` and `s_kept` (median sensitivity magnitudes).
"""
struct CollinearityFilterResult
    retained::Vector{Symbol}
    removed::DataFrame
end


"""
    filter_collinear_parameters(summary, params; psi_th=0.9)
        -> CollinearityFilterResult

Remove strongly collinear parameters from the candidate list (eq. 1-2 of the reference below).

Pairs are examined by decreasing median collinearity `ψ_ij`. Whenever `ψ_ij` exceeds the
threshold and both parameters are still in the list, the one with the smaller median sensitivity
magnitude is removed. A parameter that has already been removed no longer takes part in any pair.
A `NaN` magnitude is treated as the weakest possible, so the parameter is removed first.

This threshold is a selection threshold, distinct from the graphical `collinearity_threshold`
used to highlight pairs in the collinearity figure.

Parameters
----------
summary : CollinearitySummary
    Result of [`collinearity_summary`](@ref).
params : Vector{ParameterBound}
    Analysed parameters in summary order.
psi_th : Real, optional
    Median collinearity above which the weaker parameter of a pair is removed. Default is `0.9`
    (Goshtasbi et al.).

Returns
-------
CollinearityFilterResult
    Retained parameter names and the table of removed parameters.

Raises
------
ArgumentError
    If the threshold is not in `(0, 1)` or `summary` does not match `params`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part ii: Robust parameter subset selection, robust optimal experimental design,
and multi-step parameter identification algorithm," Journal of The Electrochemical Society,
vol. 167, p. 044505, Mar. 2020.
"""
function filter_collinear_parameters(
    summary::CollinearitySummary,
    params::Vector{ParameterBound};
    psi_th::Real = 0.9,
)::CollinearityFilterResult
    # Preliminary checks
    n_params = length(params)
    0.0 < psi_th < 1.0 ||
        throw(ArgumentError("psi_th must lie in (0, 1)"))
    size(summary.median_collinearity) == (n_params, n_params) &&
        length(summary.sensitivity_median) == n_params ||
        throw(ArgumentError("summary does not match params"))

    # Strongly collinear pairs, strongest first
    strength(s) = isfinite(s) ? s : -Inf
    pairs = [(i, j) for i in 1:n_params for j in (i + 1):n_params
             if isfinite(summary.median_collinearity[i, j]) &&
                summary.median_collinearity[i, j] > psi_th]
    sort!(pairs; by = p -> summary.median_collinearity[p...], rev = true)

    # Remove the weaker parameter of each pair still fully present
    active = trues(n_params)
    rows = NamedTuple[]
    for (i, j) in pairs
        active[i] && active[j] || continue
        weak, strong = strength(summary.sensitivity_median[i]) <= strength(summary.sensitivity_median[j]) ?
                       (i, j) : (j, i)
        active[weak] = false
        push!(rows, (removed = params[weak].name, kept = params[strong].name,
                     psi_median = summary.median_collinearity[i, j],
                     s_removed = summary.sensitivity_median[weak],
                     s_kept = summary.sensitivity_median[strong]))
    end

    removed = DataFrame(
        removed = Symbol[r.removed for r in rows], kept = Symbol[r.kept for r in rows],
        psi_median = Float64[r.psi_median for r in rows],
        s_removed = Float64[r.s_removed for r in rows],
        s_kept = Float64[r.s_kept for r in rows])
    return CollinearityFilterResult(Symbol[p.name for p in params[active]], removed)
end


# ============================================================
# PARAMETER-SUBSET SELECTION: SENSITIVITY THRESHOLD
# ============================================================

"""
    select_parameters_by_sensitivity(summary, params, retained; sensitivity_threshold=0.1)
        -> NamedTuple{(:selected, :below_threshold)}

Select, among the retained parameters, those whose median sensitivity reaches a threshold.

In plain terms: once the redundant parameters have been removed by filter_collinear_parameters,
every remaining parameter is kept if the voltage reacts strongly enough to it
(median RMS sensitivity magnitude over the sampled points), and discarded
otherwise. A parameter with a non-finite median magnitude is discarded.

This is a deliberate simplification and **not the method proposed by Goshtasbi et al.** Their
second paper chooses the subset size from the singular values of the sensitivity matrix, then
selects the subset maximising the D-optimality criterion `log det(S_I S_Iᵀ)` with a genetic
search, robustly over all sampled points. Unlike that method, the threshold selection judges each
parameter on its own: it does not detect redundancy among three or more parameters that no single
pair reveals, and the number of selected parameters follows from the threshold instead of from
the information content of the data.

Parameters
----------
summary : CollinearitySummary
    Result of [`collinearity_summary`](@ref), providing the median sensitivity magnitudes.
params : Vector{ParameterBound}
    Analysed parameters in summary order.
retained : Vector{Symbol}
    Parameters to choose from, typically `filter_collinear_parameters(...).retained`.
sensitivity_threshold : Real, optional
    Minimum median sensitivity magnitude to be selected. Default is `0.1`.

Returns
-------
NamedTuple
    `selected` and `below_threshold`, both vectors of parameter names in `retained` order.

Raises
------
ArgumentError
    If the threshold is negative or non-finite, `summary` does not match `params`, or `retained`
    contains an unknown parameter.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part ii: Robust parameter subset selection, robust optimal experimental design,
and multi-step parameter identification algorithm," Journal of The Electrochemical Society,
vol. 167, p. 044505, Mar. 2020. (Method not followed here, see above.)
"""
function select_parameters_by_sensitivity(
    summary::CollinearitySummary,
    params::Vector{ParameterBound},
    retained::Vector{Symbol};
    sensitivity_threshold::Real = 0.1,
)
    # Preliminary checks
    length(summary.sensitivity_median) == length(params) ||
        throw(ArgumentError("summary does not match params"))
    names = Symbol[p.name for p in params]
    indices = Int[something(findfirst(==(name), names), 0) for name in retained]
    all(>(0), indices) ||
        throw(ArgumentError("retained contains a parameter that was not analysed"))

    # Threshold on the median sensitivity magnitude
    influential = influential_mask(summary.sensitivity_median[indices], sensitivity_threshold)
    return (selected = retained[influential], below_threshold = retained[.!influential])
end
