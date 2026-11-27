# -*- coding: utf-8 -*-

"""
    output_sensitivity_magnitude(S_n) -> Matrix{Float64}

Combine sensitivities across all measurable model outputs into the Euclidean norm `‖s_i^k‖₂` (eq. 8 of the reference bellow):

    M[k, i] = √(Σ_j S_n[k, i, j]²)

`S_n` has shape `n_points × n_params × n_measurements`; each entry is the normalised sensitivity
of one measurable output to one parameter at one sampled point. The measurement axis includes
individual output values, such as voltages at different current densities and operating cases.

Parameters
----------
S_n : AbstractArray{<:Real,3}
    Sensitivity array with axes for sampled points, parameters, and measurable outputs.

Returns
-------
Matrix{Float64}
    Magnitudes with shape `n_points × n_params`. A point-parameter magnitude is `NaN` if any
    sensitivity in its output vector is `NaN`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function output_sensitivity_magnitude(S_n::AbstractArray{<:Real,3})::Matrix{Float64}
    return dropdims(sqrt.(sum(abs2, S_n; dims = 3)); dims = 3)
end


"""
    sensitivity_statistics(S_n) -> NamedTuple{(:median, :std)}

Compute summary statistics of sensitivity magnitudes over sampled points.

- `median[i] = s̄_i = med_k ‖s_i^k‖₂`: robust measure of the influence of parameter `i`.
- `std[i] = s*_i = √var_k ‖s_i^k‖₂` (corrected sample variance, `n − 1`; a single valid point gives `NaN`): dependence of that influence on the other parameters.

Points with a `NaN` magnitude (failed runs) are ignored.

Parameters
----------
S_n : AbstractArray{<:Real,3}
    Sensitivity array with shape `n_points × n_params × n_measurements`.

Returns
-------
NamedTuple
    Named tuple containing `median` and corrected sample `std` vectors, one value per parameter.
    A parameter with no valid points has `NaN` statistics; one valid point has a `NaN` standard deviation.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function sensitivity_statistics(S_n::AbstractArray{<:Real,3})
    # Compute one magnitude per point and parameter
    M = output_sensitivity_magnitude(S_n)

    # Exclude failed points and summarize each parameter across the remaining points
    valid(i) = filter(!isnan, @view M[:, i])
    med = [(v = valid(i); isempty(v) ? NaN : median(v)) for i in axes(M, 2)]
    sd  = [(v = valid(i); isempty(v) ? NaN : std(v)) for i in axes(M, 2)]
    return (median = med, std = sd)
end


"""
    collinearity_index(S; floor = 1e-12) -> Matrix{Float64}

Compute pairwise collinearity of parameter sensitivity vectors at one sampled point.

For rows `s_i` and `s_j` of `S`, the index is `ψ_ij = |s_i · s_j| / (‖s_i‖ ‖s_j‖)`. Values near one
indicate effects in the same direction; zero indicates orthogonal effects. The norm product is bounded
below by `floor` to avoid division by zero. A row containing `NaN` propagates `NaN` to the corresponding
row and column.

Parameters
----------
S : AbstractMatrix{<:Real}
    Sensitivity matrix of size `n_params × n_measurements`.
floor : Real, optional
    Minimum denominator used for zero or very small norm products. Default is `1e-12`.

Returns
-------
Matrix{Float64}
    Pairwise collinearity matrix of size `n_params × n_params`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function collinearity_index(S::AbstractMatrix{<:Real}; floor::Real = 1e-12)::Matrix{Float64}
    # Compute one sensitivity-vector norm per parameter
    norms = vec(sqrt.(sum(abs2, S; dims = 2)))              # ‖s_i‖ for each parameter i

    # Normalize pairwise dot products, protecting the denominator with floor
    return abs.(S * S') ./ max.(norms * norms', floor)      # |s_i · s_j| / (‖s_i‖ ‖s_j‖)
end


"""
    collinearity_indices(S_n; floor = 1e-12) -> Array{Float64,3}

Compute pairwise parameter collinearity at every sampled point.

Parameters
----------
S_n : AbstractArray{<:Real,3}
    Sensitivity array with shape `n_points × n_params × n_measurements`.
floor : Real, optional
    Minimum denominator passed to [`collinearity_index`](@ref). Default is `1e-12`.

Returns
-------
Array{Float64,3}
    Collinearity array with shape `n_points × n_params × n_params`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function collinearity_indices(S_n::AbstractArray{<:Real,3}; floor::Real = 1e-12)::Array{Float64,3}
    # Allocate an output matrix for each sampled point
    n_points, n_params, _ = size(S_n)
    Ψ = Array{Float64,3}(undef, n_points, n_params, n_params)

    # Compute the parameter-pair indices independently at each point
    for k in 1:n_points
        Ψ[k, :, :] = collinearity_index(@view S_n[k, :, :]; floor)
    end
    return Ψ
end


"""
    median_collinearity(S_n; floor = 1e-12) -> Matrix{Float64}

Compute the median collinearity of every parameter pair over sampled points.

Parameters
----------
S_n : AbstractArray{<:Real,3}
    Sensitivity array with shape `n_points × n_params × n_measurements`.
floor : Real, optional
    Minimum denominator passed to [`collinearity_index`](@ref). Default is `1e-12`.

Returns
-------
Matrix{Float64}
    Median collinearity matrix of size `n_params × n_params`. NaN values from failed points are
    ignored; pairs with no valid points have `NaN`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function median_collinearity(S_n::AbstractArray{<:Real,3}; floor::Real = 1e-12)::Matrix{Float64}
    # Compute pointwise collinearity before aggregating across points
    Ψ = collinearity_indices(S_n; floor)
    n_params = size(Ψ, 2)

    # Ignore failed points when taking each pair's median
    return [(v = filter(!isnan, @view Ψ[:, i, j]); isempty(v) ? NaN : median(v))
            for i in 1:n_params, j in 1:n_params]
end


"""
    influential_mask(sensitivity_median, threshold) -> BitVector

Identify parameters whose median sensitivity magnitude reaches a threshold.

Non-finite magnitudes (`NaN` or `Inf`) are excluded. Collinearity is interpreted only for influential
parameters.

Parameters
----------
sensitivity_median : AbstractVector{<:Real}
    Median sensitivity magnitude for each parameter.
threshold : Real
    Minimum magnitude for a parameter to be considered influential.

Returns
-------
BitVector
    Boolean mask indicating influential parameters.

Raises
------
ArgumentError
    If `threshold` is negative or non-finite.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function influential_mask(sensitivity_median::AbstractVector{<:Real}, threshold::Real)::BitVector
    # Validate the cutoff before classifying parameters
    isfinite(threshold) && threshold >= 0 ||
        throw(ArgumentError("threshold must be finite and non-negative"))

    # Keep only finite sensitivities that meet the influence threshold
    return BitVector(isfinite.(sensitivity_median) .& (sensitivity_median .>= threshold))
end
