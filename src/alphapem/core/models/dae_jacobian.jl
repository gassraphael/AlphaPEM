# -*- coding: utf-8 -*-

"""Sparse Jacobian helpers for the DAE/IDA integration path.

The coloring-based approach reduces the cost of finite-difference Jacobian
evaluation from 2n residual calls (one per column) to 2·ncolors calls, where
ncolors is the chromatic number of the column-conflict graph of the Jacobian
pattern.  For the AlphaPEM DAE, whose Jacobian is block-sparse with local
1D couplings, ncolors ≪ n, leading to a significant speedup in every Newton
step.
"""

using SparseArrays: sparse, SparseMatrixCSC, rowvals, nzrange, nonzeros
using SparseMatrixColorings: ColoringProblem, GreedyColoringAlgorithm,
                              coloring, ncolors, column_colors


# ──────────────────────────────────────────────────────────────────────────────
# Jacobian coloring cache
# ──────────────────────────────────────────────────────────────────────────────

"""Pre-computed coloring data and reusable work buffers for the FD Jacobian.

Built once at solve-setup time by `_build_jacobian_coloring_cache` and passed
to every `_dae_jacobian_fd!` call, eliminating all per-call allocations.

Fields
------
- `coloring_result` : column coloring of the Jacobian sparsity pattern.
- `colors` : `column_colors(coloring_result)` copied out for tight-loop access.
- `ncolors_count` : number of color groups = number of FD perturbations needed.
- `J_compressed` : dense `n × ncolors_count` buffer; column c stores the raw
  two-point differences `F⁺[i] − F⁻[i]` for color group c (before δ-division).
- `y_work` : length-n copy of the current state used during perturbations.
- `res_plus` : residual buffer at `y + Δ`.
- `res_minus` : residual buffer at `y − Δ`.
- `deltas` : per-column FD step sizes `δⱼ = cbrt(ε)·max(|yⱼ|,1)`, recomputed
  each call; pre-allocated to avoid per-call allocation.
- `diag_nz_indices` : for each state index `i`, the CSC storage index of `J[i,i]`,
  pre-cached to add the `γ·I` term on differential rows in O(n).
"""
struct JacobianColoringCache
    coloring_result
    colors::Vector{Int}
    ncolors_count::Int
    J_compressed::Matrix{Float64}
    y_work::Vector{Float64}
    res_plus::Vector{Float64}
    res_minus::Vector{Float64}
    deltas::Vector{Float64}
    diag_nz_indices::Vector{Int}
end

"""Build a `JacobianColoringCache` from the Jacobian sparsity prototype.

Runs the greedy column coloring once at setup, then pre-allocates all work
buffers so that subsequent `_dae_jacobian_fd!` calls are allocation-free.

The diagonal-index cache is computed here rather than in `simulate_model!`,
centralising all Jacobian-related setup.
"""
function _build_jacobian_coloring_cache(jac_prototype::SparseMatrixCSC,
                                        differential_vars::BitVector)::JacobianColoringCache
    n = size(jac_prototype, 1)
    n == size(jac_prototype, 2) ||
        throw(ArgumentError("jac_prototype must be square."))
    n == length(differential_vars) ||
        throw(ArgumentError("differential_vars size mismatch in _build_jacobian_coloring_cache."))

    # Greedy column coloring of the Jacobian sparsity pattern.
    # :nonsymmetric/:column coloring guarantees that no two same-colored columns
    # share a nonzero row — the precondition for compressed FD evaluation.
    result  = coloring(jac_prototype,
                       ColoringProblem{:nonsymmetric, :column}(),
                       GreedyColoringAlgorithm())
    nc      = ncolors(result)
    colors  = column_colors(result)

    # Pre-allocate work buffers — reused at every Newton step, zero runtime allocation.
    J_compressed = zeros(Float64, n, nc)
    y_work       = zeros(Float64, n)
    res_plus     = Vector{Float64}(undef, n)
    res_minus    = Vector{Float64}(undef, n)
    deltas       = Vector{Float64}(undef, n)

    # Cache the CSC storage index of each diagonal entry J[i,i].
    jac_rows        = rowvals(jac_prototype)
    diag_nz_indices = Vector{Int}(undef, n)
    @inbounds for j in 1:n
        diag_idx = 0
        for idx in nzrange(jac_prototype, j)
            if jac_rows[idx] == j
                diag_idx = idx
                break
            end
        end
        diag_idx == 0 &&
            throw(ArgumentError("Jacobian prototype must contain all diagonal entries (missing J[$j,$j])."))
        diag_nz_indices[j] = diag_idx
    end

    return JacobianColoringCache(result, colors, nc, J_compressed,
                                 y_work, res_plus, res_minus, deltas, diag_nz_indices)
end


# ──────────────────────────────────────────────────────────────────────────────
# Jacobian sparsity prototype
# ──────────────────────────────────────────────────────────────────────────────

"""Build a sparse Jacobian prototype for the DAE residual.

The prototype is built once at solve setup time using a conservative strategy:
- always include the diagonal,
- add numerically detected couplings around the initial state.

This keeps the setup robust for `IDA(linear_solver=:KLU)` while avoiding a
fully manual Jacobian implementation.
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

    # 2) Numerically probe local couplings by central finite differences.
    # Two probe points are used to capture the complete structural pattern:
    #
    # Pass A — initial state `y0`:  detects first-order couplings active at
    #   the start of the simulation.  Uses a relative threshold scaled by the
    #   local residual magnitude so that small sensitivities in large-residual
    #   rows are not falsely included.
    #
    # Pass B — activated state `max.(y0, 0.1)`:  forces all state variables
    #   to at least 0.1 in scaled space.  This reveals dormant higher-order
    #   couplings (e.g. liquid-water transport) that vanish at y0 but
    #   become significant once the state evolves away from zero.
    #   The activated state is far from equilibrium, so residuals are large and
    #   `local_scale` can reach O(10²–10³).  Using a relative threshold would
    #   mask genuine couplings of O(0.01–0.1), so Pass B uses an absolute-only
    #   threshold (local_rtol = 0.0) to avoid missing those entries.
    #
    # The union of both passes gives a conservative structural pattern that
    #   remains valid throughout the entire simulation.
    fd_eps           = cbrt(eps(Float64))
    sensitivity_atol = 1e-14
    sensitivity_rtol = 1e-8

    y_work              = copy(initial_solver_values)
    res_perturbed_plus  = Vector{Float64}(undef, n)
    res_perturbed_minus = Vector{Float64}(undef, n)

    # Activated probe state: lift every variable to at least 0.1 in scaled
    # space so that higher-order terms (e.g. liquid-water transport ∝ s^3)
    # produce detectable FD sensitivities despite vanishing at y0.
    y_activated    = max.(initial_solver_values, 0.1)
    res0_activated = zeros(Float64, n)
    residual!(res0_activated, initial_solver_derivatives, y_activated, packed, t0)

    for (y_probe, res0_probe, local_rtol) in (
            (initial_solver_values, res0,           sensitivity_rtol),  # Pass A: relative threshold
            (y_activated,           res0_activated, 0.0))               # Pass B: absolute threshold only
        for j in 1:n
            yj    = y_probe[j]
            delta = fd_eps * max(abs(yj), 1.0)

            copyto!(y_work, y_probe)
            y_work[j] = yj + delta
            residual!(res_perturbed_plus, initial_solver_derivatives, y_work, packed, t0)

            y_work[j] = yj - delta
            residual!(res_perturbed_minus, initial_solver_derivatives, y_work, packed, t0)

            inv_2delta = 0.5 / delta
            @inbounds for i in 1:n
                i == j && continue  # Diagonal already added; skip.
                sensitivity  = (res_perturbed_plus[i] - res_perturbed_minus[i]) * inv_2delta
                local_scale  = max(abs(res0_probe[i]), abs(res_perturbed_plus[i]),
                                   abs(res_perturbed_minus[i]), 1.0)
                abs(sensitivity) > (sensitivity_atol + local_rtol * local_scale) || continue
                push!(rows, i)
                push!(cols, j)
            end
        end
    end

    proto = sparse(rows, cols, ones(Float64, length(rows)), n, n)
    return proto
end


# ──────────────────────────────────────────────────────────────────────────────
# Coloring-based FD Jacobian evaluator
# ──────────────────────────────────────────────────────────────────────────────

"""Fill the DAE Jacobian in place for IDA using coloring-compressed finite differences.

The matrix expected by IDA is:
    J = dF/dy + γ · dF/d(dy/dt)

`dF/dy` is approximated by coloring-compressed central finite differences.
Same-colored columns are independent (no shared nonzero row by construction of
the column coloring) and can be perturbed simultaneously, reducing the residual
call count from 2n to 2·ncolors.

`dF/d(dy/dt)` equals the identity on differential rows and zero on algebraic
rows (exact for the residual convention `F_diff = dydt_IDA − f(y)`), and is
added analytically using the pre-cached `diag_nz_indices`.

FD step choice
--------------
Per-column adaptive steps `δⱼ = cbrt(ε)·max(|yⱼ|,1)` match the original
column-by-column implementation exactly.  When multiple columns of the same
color group are perturbed simultaneously, each column j uses its own δⱼ.
Due to color independence (no two same-colored columns share a nonzero row),
row i only sees the contribution of one column j(i,c) per color group c, so:
    (F⁺[i] − F⁻[i]) = 2·δ_{j(i,c)}·J[i, j(i,c)]
The raw difference (without δ-division) is stored in `J_compressed[:, c]`, and
the per-column division 1/(2·δⱼ) is applied during decompression.
"""
function _dae_jacobian_fd!(J,
                           dydt_IDA::Vector{Float64},
                           y::Vector{Float64},
                           packed,
                           gamma::Float64,
                           t::Float64,
                           residual!,
                           differential_vars::BitVector,
                           cache::JacobianColoringCache)
    n = length(y)
    n == length(dydt_IDA) ||
        throw(ArgumentError("State/derivative size mismatch in _dae_jacobian_fd!."))
    n == length(differential_vars) ||
        throw(ArgumentError("differential_vars size mismatch in _dae_jacobian_fd!."))

    fd_eps = cbrt(eps(Float64))

    # Compute per-column adaptive FD step sizes (mirrors the original formula).
    @inbounds for j in 1:n
        cache.deltas[j] = fd_eps * max(abs(y[j]), 1.0)
    end

    # Reset the compressed buffer and initialise the working state copy.
    fill!(cache.J_compressed, 0.0)
    copyto!(cache.y_work, y)

    # ── Compressed finite-difference sweep ──────────────────────────────────
    # For each color group c, columns with `colors[j] == c` share no nonzero
    # row, so perturbing them all simultaneously gives one independent FD
    # estimate per nonzero row from a single pair of residual evaluations.
    # Each column j is perturbed by its own δⱼ; the raw two-point differences
    # are stored in J_compressed[:, c] for division during decompression.
    for c in 1:cache.ncolors_count

        # Apply +δⱼ to every column in group c.
        @inbounds for j in 1:n
            cache.colors[j] == c || continue
            cache.y_work[j] = y[j] + cache.deltas[j]
        end
        residual!(cache.res_plus, dydt_IDA, cache.y_work, packed, t)

        # Apply −δⱼ to every column in group c.
        @inbounds for j in 1:n
            cache.colors[j] == c || continue
            cache.y_work[j] = y[j] - cache.deltas[j]
        end
        residual!(cache.res_minus, dydt_IDA, cache.y_work, packed, t)

        # Restore working state for group c and store raw two-point differences.
        @inbounds for j in 1:n
            cache.colors[j] == c || continue
            cache.y_work[j] = y[j]
        end
        @inbounds for i in 1:n
            cache.J_compressed[i, c] = cache.res_plus[i] - cache.res_minus[i]
        end
    end

    # ── Decompression with per-column δ-division ─────────────────────────────
    # J[i, j] = J_compressed[i, colors[j]] / (2·δⱼ) for each stored nonzero.
    # Due to color independence, J_compressed[i, colors[j]] = F⁺[i] − F⁻[i]
    # = 2·δⱼ·J[i,j], so dividing by 2·δⱼ recovers the correct entry.
    fill!(J, 0.0)
    vals     = nonzeros(J)
    jac_rows = rowvals(J)
    @inbounds for j in 1:n
        inv_2delta_j = 0.5 / cache.deltas[j]
        c = cache.colors[j]
        for idx in nzrange(J, j)
            vals[idx] = cache.J_compressed[jac_rows[idx], c] * inv_2delta_j
        end
    end

    # ── Add γ · dF/d(dy/dt) on differential rows ────────────────────────────
    @inbounds for i in eachindex(differential_vars)
        differential_vars[i] || continue
        vals[cache.diag_nz_indices[i]] += gamma
    end

    return nothing
end
