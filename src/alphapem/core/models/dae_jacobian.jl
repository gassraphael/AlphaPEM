# -*- coding: utf-8 -*-

"""Sparse Jacobian helpers for the DAE/IDA integration path.
"""

using SparseArrays: sparse, SparseMatrixCSC, rowvals, nzrange, nonzeros

"""Build a sparse Jacobian prototype for the DAE residual.

 The prototype is built once at solve setup time using a conservative strategy:
 - always include the diagonal,
 - add numerically detected couplings around the initial state.

This keeps the setup robust for `IDA(linear_solver=:KLU)` while avoiding a fully
manual Jacobian implementation.
"""
function _build_dae_jacobian_prototype(residual!,
                                       packed,
                                       initial_solver_derivatives::Vector{Float64},
                                       initial_solver_values::Vector{Float64},
                                       t0::Float64,
                                       differential_vars::BitVector)::SparseMatrixCSC{Float64, Int}
    # Initialisation
    n = length(initial_solver_values)
    n == length(initial_solver_derivatives) ||
        throw(ArgumentError("State/derivative size mismatch in _build_dae_jacobian_prototype."))
    n == length(differential_vars) ||
        throw(ArgumentError("differential_vars size mismatch in _build_dae_jacobian_prototype."))
    # Baseline residual at the initial scaled state.
    res0 = zeros(Float64, n)
    residual!(res0, initial_solver_derivatives, initial_solver_values, packed, t0)

    # `rows[k], cols[k]` is the position of the k-th stored nonzero entry.
    rows = Int[] # Row indices of nonzero entries in the Jacobian prototype.
    cols = Int[] # Column indices of nonzero entries in the Jacobian prototype.

    # 1) Always include the diagonal.
    for i in 1:n
        push!(rows, i)
        push!(cols, i)
    end

    # 2) Numerically probe local couplings around (t0, y0, ydot0):
    # perturb each state one at a time, recompute the residual, estimate sensitivities
    # by central finite differences, and keep only entries above an adaptive noise threshold.
    fd_eps = cbrt(eps(Float64)) # Robust FD perturbation step for structural detection.
    sensitivity_atol = 1e-14 # Ignore finite difference noise
    sensitivity_rtol = 1e-8
    for j in 1:n
        yj = initial_solver_values[j]
        delta = fd_eps * max(abs(yj), 1.0)

        y_perturbed_plus = copy(initial_solver_values)
        y_perturbed_plus[j] += delta
        res_perturbed_plus = similar(res0) # allocates same type/size uninitialized;
        residual!(res_perturbed_plus, initial_solver_derivatives, y_perturbed_plus, packed, t0) # computes the residual at the perturbed state

        y_perturbed_minus = copy(initial_solver_values)
        y_perturbed_minus[j] -= delta
        res_perturbed_minus = similar(res0)
        residual!(res_perturbed_minus, initial_solver_derivatives, y_perturbed_minus, packed, t0) # computes the residual at the perturbed state

        inv_2delta = 0.5 / delta # Precompute the inverse of 2*delta for central difference sensitivity estimation.
        @inbounds for i in 1:n
            sensitivity = (res_perturbed_plus[i] - res_perturbed_minus[i]) * inv_2delta # Central finite difference sensitivity estimate for dF[i]/dy[j]
            local_scale = max(abs(res0[i]), abs(res_perturbed_plus[i]), abs(res_perturbed_minus[i]), 1.0)
            abs(sensitivity) > (sensitivity_atol + sensitivity_rtol * local_scale) || continue # Skip entries that are below the noise threshold.
            push!(rows, i)
            push!(cols, j)
        end
    end

    # Sparse structural pattern. Values are placeholders; only the pattern is used.
    return sparse(rows, cols, ones(Float64, length(rows)), n, n)
end

"""Fill the DAE Jacobian matrix in place for IDA.

The matrix expected by IDA is:
    J = dF/dy + gamma * dF/d(dy/dt)

`dF/dy` is approximated by finite differences on `y` at fixed `dy/dt`.
The second term is exact for this model because differential residuals are
written as `F_diff = dydt_IDA - f(y)`, hence `dF/d(dy/dt)` is the identity on
`differential_vars` and zero on algebraic rows.
"""
function _dae_jacobian_fd!(J,
                           dydt_IDA::Vector{Float64},
                           y::Vector{Float64},
                           packed,
                           gamma::Float64,
                           t::Float64,
                           residual!,
                           differential_vars::BitVector)
    n = length(y)
    n == length(dydt_IDA) ||
        throw(ArgumentError("State/derivative size mismatch in _dae_jacobian_fd!."))
    n == length(differential_vars) ||
        throw(ArgumentError("differential_vars size mismatch in _dae_jacobian_fd!."))

    # Baseline residual.
    res0 = zeros(Float64, n)
    residual!(res0, dydt_IDA, y, packed, t)

    # Finite-difference Jacobian of dF/dy.
    # Important: write only to entries that already exist in the sparse matrix
    # to keep a fixed sparsity pattern compatible with Sundials sparse matrix handles.
    fd_eps = cbrt(eps(Float64))  # Robust FD step for Jacobian values.
    fill!(J, 0.0)      # Reset nonzero values, keep sparse structure
    rows = rowvals(J)  # for each stored position `idx` in column `j`, `rows[idx]` is the row `i` of that exact nonzero entry.
    vals = nonzeros(J) # `vals[idx]` is `J[rows[idx], j]`, so `vals[idx] = ...` updates that entry in place, with no sparsity change and no extra allocation.

    for j in 1:n
        yj = y[j]
        delta = fd_eps * max(abs(yj), 1.0)

        y_perturbed = copy(y)
        y_perturbed[j] += delta

        res_perturbed = similar(res0) # allocates same type/size uninitialized;
        residual!(res_perturbed, dydt_IDA, y_perturbed, packed, t) # computes the residual at the perturbed state

        inv_delta = 1.0 / delta
        @inbounds for idx in nzrange(J, j) # Iterate over the stored nonzero entries of column `j` in `J`
            i = rows[idx]
            vals[idx] = (res_perturbed[i] - res0[i]) * inv_delta # Finite-difference sensitivity for dF[i]/dy[j]
        end
    end

    # Add gamma * dF/d(dy/dt): identity on differential rows only.
    @inbounds for i in eachindex(differential_vars)
        differential_vars[i] || continue
        J[i, i] += gamma
    end

    return nothing
end

