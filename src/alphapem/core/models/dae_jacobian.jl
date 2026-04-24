# -*- coding: utf-8 -*-

"""Sparse Jacobian helpers for the DAE/IDA integration path.

This file keeps Jacobian-related setup outside `AlphaPEM.jl` so that
`simulate_model!` remains focused on orchestration while the sparse linear-algebra
plumbing required by `IDA(linear_solver=:KLU)` stays isolated and testable.
"""

using SparseArrays: sparse, SparseMatrixCSC, rowvals, nzrange, nonzeros

"""Build a sparse Jacobian prototype for the DAE residual.

The prototype is built once at solve setup time using a conservative strategy:
- always include the diagonal,
- always include couplings involving algebraic rows/columns,
- add numerically detected differential couplings around the initial state.

This keeps the setup robust for `IDA(linear_solver=:KLU)` while avoiding a fully
manual Jacobian implementation.
"""
function _build_dae_jacobian_prototype(residual!,
                                       packed,
                                       initial_solver_derivatives::Vector{Float64},
                                       initial_solver_values::Vector{Float64},
                                       t0::Float64,
                                       differential_vars::BitVector)::SparseMatrixCSC{Float64, Int}
    n = length(initial_solver_values)
    n == length(initial_solver_derivatives) ||
        throw(ArgumentError("State/derivative size mismatch in _build_dae_jacobian_prototype."))
    n == length(differential_vars) ||
        throw(ArgumentError("differential_vars size mismatch in _build_dae_jacobian_prototype."))

    n_diff = count(differential_vars)
    n_alg = n - n_diff

    # Baseline residual at the initial scaled state.
    res0 = zeros(Float64, n)
    residual!(res0, initial_solver_derivatives, initial_solver_values, packed, t0)

    rows = Int[]
    cols = Int[]

    # 1) Always include the diagonal.
    for i in 1:n
        push!(rows, i)
        push!(cols, i)
    end

    # 2) Conservative dense couplings for algebraic rows/columns.
    #    The algebraic block is small (`2*nb_gc + 3`) and can couple broadly.
    if n_alg > 0
        alg_first = n_diff + 1
        for j in alg_first:n
            for i in 1:n
                push!(rows, i)
                push!(cols, j)
            end
        end
        for i in alg_first:n
            for j in 1:n
                push!(rows, i)
                push!(cols, j)
            end
        end
    end

    # 3) Numerical probing for differential-to-differential couplings near (t0, y0, ydot0).
    fd_eps = sqrt(eps(Float64))
    sensitivity_threshold = 1e-12
    for j in 1:n_diff
        yj = initial_solver_values[j]
        delta = fd_eps * max(abs(yj), 1.0)

        y_perturbed = copy(initial_solver_values)
        y_perturbed[j] += delta
        res_perturbed = similar(res0)
        residual!(res_perturbed, initial_solver_derivatives, y_perturbed, packed, t0)

        inv_delta = 1.0 / delta
        @inbounds for i in 1:n_diff
            abs((res_perturbed[i] - res0[i]) * inv_delta) > sensitivity_threshold || continue
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
    # Important: write only existing CSC nonzeros to keep a fixed sparsity
    # pattern compatible with Sundials sparse matrix handles.
    fd_eps = sqrt(eps(Float64))
    fill!(J, 0.0)
    rows = rowvals(J)
    vals = nonzeros(J)
    for j in 1:n
        yj = y[j]
        delta = fd_eps * max(abs(yj), 1.0)

        y_perturbed = copy(y)
        y_perturbed[j] += delta

        res_perturbed = similar(res0)
        residual!(res_perturbed, dydt_IDA, y_perturbed, packed, t)

        inv_delta = 1.0 / delta
        @inbounds for idx in nzrange(J, j)
            i = rows[idx]
            vals[idx] = (res_perturbed[i] - res0[i]) * inv_delta
        end
    end

    # Add gamma * dF/d(dy/dt): identity on differential rows only.
    @inbounds for i in eachindex(differential_vars)
        differential_vars[i] || continue
        J[i, i] += gamma
    end

    return nothing
end

