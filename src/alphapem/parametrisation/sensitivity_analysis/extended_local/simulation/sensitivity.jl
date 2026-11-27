# -*- coding: utf-8 -*-

"""
    equidistant_levels(n_levels::Int) -> Vector{Float64}

Generate the equidistant unit-interval levels used to sweep each parameter.

Parameters
----------
n_levels : Int
    Number of levels to generate.

Returns
-------
Vector{Float64}
    `n_levels` values spanning `[0, 1]`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
equidistant_levels(n_levels::Int)::Vector{Float64} = collect(range(0.0, 1.0; length = n_levels))


"""
    extended_local_sensitivity_matrix(simulate, θ_k, n_levels) -> Matrix{Float64}

Compute normalised local sensitivities at one point of the unit hypercube.

Each parameter is swept independently over `n_levels` equidistant values while all other coordinates
remain fixed. For measurement `j`, with `ỹ_m` the outputs at the successful levels `θ̄_m`:

    s_{j,i} = mean_m(ỹ_{j,m})⁻¹ · mean_m[(ỹ_{j,m+1} − ỹ_{j,m}) / (θ̄_{m+1} − θ̄_m)]

The sensitivity is signed: its sign is the sign of the output slope as the normalized parameter
increases. For the positive voltage outputs analysed here, a positive value means voltage increases
with the parameter, while a negative value means voltage decreases.

Parameters
----------
simulate : Function
    Maps a unit-hypercube point to a measurement vector, or returns `nothing` if simulation fails.
θ_k : AbstractVector{<:Real}
    Evaluation point in the unit hypercube.
n_levels : Int
    Number of equidistant sweep levels.

Returns
-------
Matrix{Float64}
    Matrix of size `n_params × n_measurements`. Failed levels and non-finite outputs are omitted;
    if fewer than two valid levels remain for a parameter, its row is `NaN`.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function extended_local_sensitivity_matrix(simulate::Function,
                                           θ_k::AbstractVector{<:Real},
                                           n_levels::Int)::Matrix{Float64}
    # Step 1: define the equidistant sweep levels and the number of parameters in θ_k.
    n_levels >= 2 || throw(ArgumentError("n_levels must be at least 2"))
    levels = equidistant_levels(n_levels)
    n_params = length(θ_k)  # Number of model parameters (dimension of θ_k).
    S = nothing # S[i,j] stores the sensitivity of output j to parameter i.

    # Step 2: sweep one parameter i at a time, keeping all other coordinates at θ_k.
    for i in 1:n_params # i is the parameter index varied in this sweep; every other coordinate stays at θ_k.
        ys = Vector{Vector{Float64}}()   # Successful output vectors ỹ_m for parameter i.
        valid_levels = Float64[]         # Matching values θ̄_{i,m} of the swept coordinate.

        # Step 3: simulate the model at each perturbed vector θ̄ and retain successful evaluations.
        for level in levels
            θ = collect(Float64, θ_k)    # Perturbed parameter vector θ̄; starts from the sample θ_k.
            θ[i] = level                 # Set θ̄_i to the current equidistant level.
            y = simulate(θ)              # Model outputs ỹ at this θ̄ (one value per measurement j).
            (y === nothing || !all(isfinite, y)) && continue  # Skip failed runs or outputs containing NaN/Inf.
            push!(ys, y)
            push!(valid_levels, level)
        end
        n_meas = isempty(ys) ? 0 : length(ys[1])  # Number of output measurements j.
        if S === nothing
            isempty(ys) && continue
            S = fill(NaN, n_params, n_meas)
        end

        # Leave the row as NaN when fewer than two valid levels remain
        length(ys) >= 2 || continue

        # Step 4: assemble outputs by measurement and compute the normalised finite-difference sensitivity for parameter i.
        Y = reduce(hcat, ys)  # Y[j,m] = ỹ_{j,m}, with columns ordered by valid θ̄_{i,m}.
        slopes = diff(Y; dims = 2) ./ reshape(diff(valid_levels), 1, :) # Adjacent finite-difference slopes: (ỹ_{j,m+1} - ỹ_{j,m}) / (θ̄_{i,m+1} - θ̄_{i,m}).
        S[i, :] = vec(sum(slopes; dims = 2)) ./ (size(slopes, 2) .* vec(sum(Y; dims = 2)) ./ size(Y, 2)) # Divide the mean slope by mean(ỹ_j), equivalent to (n_s / Σ_m ỹ_{j,m}) times the mean slope. Averages use only successful levels.
    end

    # Step 5: return the complete parameter-by-measurement sensitivity matrix.
    S === nothing && throw(ArgumentError("All simulations failed at θ_k: cannot size the sensitivity matrix"))
    return S
end
