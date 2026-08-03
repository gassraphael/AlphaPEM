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

"""Build the set of probe states used for numerical sparsity detection.

Returns a vector of `(y_probe, local_rtol)` pairs.

Why several probes, and why non-uniform ones
--------------------------------------------
The pattern is detected by finite differences, so a coupling is only visible at
a probe state where its partial derivative is *numerically* non-zero. The
initial state `y0` is a no-flow thermodynamic equilibrium: every field
(`C_v`, `C_H2`, `C_O2`, `C_N2`, `T`, `s`, `lambda`) is spatially **uniform**
across the through-plane nodes, so every diffusive/convective flux is exactly
zero.

That matters because most transport terms have the form `J = -D(s, T, C) · ∇C`.
At a uniform state `∇C = 0`, hence

    ∂J/∂s = -(∂D/∂s)·∇C = 0,   ∂J/∂T = 0,   ∂J/∂C_other = 0

*identically*, even though these couplings are large as soon as the simulation
develops gradients. A uniform probe therefore cannot see any
diffusivity-modulating coupling, no matter how permissive the threshold is.
Lifting all variables by a constant (`max.(y0, 0.1)`) does not help: the state
stays uniform.

Missing entries are not benign here. The Jacobian is evaluated with
*coloring-compressed* finite differences, which assume that same-coloured
columns share no non-zero row. An entry absent from the pattern lets a second
column of the same colour write into a row it was assumed not to touch, so the
recovered value for the *legitimate* entry becomes the sum of two derivatives.
The Jacobian is then not merely incomplete but wrong, IDA's modified Newton
stops converging, and the solver collapses to tiny steps with a near-100%
non-linear-solver failure rate.

The probe set is therefore:
- Pass A  — `y0`, relative threshold: couplings active at the start.
- Pass B  — `max.(y0, 0.1)`, absolute threshold: terms that vanish at `s = 0`.
- Pass C  — several **non-uniform** states obtained by scattering `max.(|y0|, 0.1)`
  with a deterministic low-discrepancy (golden-ratio) per-index factor. Every
  node then differs from its neighbours, so all gradients are non-zero and the
  flux-modulating couplings become visible.

Pass B and Pass C fix two *different* degeneracies
------------------------------------------------------------------------
- Pass B stays spatially **uniform** (every node still equals every other node);
  it only moves state variables away from `0`. It reveals *local*, same-node
  couplings whose derivative vanishes at `y0` because a variable itself is zero
  there (e.g. liquid saturation `s = 0` at rest kills `∂(term ∝ s^3)/∂s`, which
  is present but dormant at equilibrium).
- Pass C breaks spatial **uniformity**. It reveals *transport*, cross-node
  couplings whose derivative is proportional to a gradient `∇C`, `∇T`, ...  —
  these are zero at *any* uniform state, however far from `0` each variable
  individually sits, so Pass B alone cannot expose them. In practice this is
  the majority of the entries missed by Pass A/B alone (diffusivity-modulating
  terms of the form `∂(D(s,T,C)·∇C)/∂s`, `.../∂T`, etc., see above).
An entry may need either mechanism, both, or neither; hence probing with all
of Pass A, B and C rather than trying to guess which applies where.

Golden-ratio scattering (Pass C)
---------------------------------
Pass C needs a deterministic, reproducible spatial pattern that puts every
solver index at a *different* value from its neighbours, for any mesh size
(`nb_gdl`, `nb_mpl` here are small, typically 2-5). This is exactly what a Weyl
sequence provides:

    spread[k] = 2 * frac(k·α + phase) - 1,     α irrational

By Weyl's equidistribution theorem, `frac(k·α)` never repeats and spreads
uniformly over `[0,1)` whenever `α` is irrational. The golden ratio
`α = φ-1 = (√5-1)/2` is the irrational number *hardest to approximate by a
rational `p/q`* (every term of its continued-fraction expansion is `1`), so
`frac(k·α)` avoids near-periodicities even for the small, regular mesh sizes
used here — a generic irrational could resonate with a mesh size such as 3 or
5 and, by accident, put two neighbouring nodes back at nearly equal values.
A plain linear ramp does not offer this guarantee (it can create spurious
cancellations between two physically-coupled variables that happen to vary at
the same rate along the index), and `rand()` would make the detected pattern —
hence the coloring, hence the solver's step-by-step behaviour — different on
every run. The per-pass phase shift `0.37*m` only decorrelates the three
amplitude levels (`m = 1, 2, 3`) from one another; it carries no special
meaning beyond "not the same offset twice".

The Pass C states are deterministic, so the pattern is reproducible run to run.
"""
function _dae_jacobian_probe_states(initial_solver_values::Vector{Float64};
                                    amplitudes = (0.10, 0.25, 0.45),
                                    sensitivity_rtol::Float64 = 1e-8)
    n = length(initial_solver_values)
    probes = Vector{Tuple{Vector{Float64}, Float64}}()

    # Pass A — equilibrium state, relative threshold.
    push!(probes, (copy(initial_solver_values), sensitivity_rtol))
    # Pass B — activated state, absolute threshold only.
    push!(probes, (max.(initial_solver_values, 0.1), 0.0))

    # Pass C — non-uniform states. The golden-ratio sequence gives a spread that
    # never repeats between neighbouring indices, so no gradient stays zero.
    golden = 0.6180339887498949     # golden ratio:  φ - 1 = (√5 - 1)/2
    base = max.(abs.(initial_solver_values), 0.1)
    for (m, amplitude) in enumerate(amplitudes)
        spread = [2.0 * ((k * golden + 0.37 * m) % 1.0) - 1.0 for k in 1:n] # Weyl sequence
        push!(probes, (base .* (1.0 .+ amplitude .* spread), 0.0))
    end
    return probes
end


"""Build a sparse Jacobian prototype for the DAE residual.

The prototype is built once at solve setup time using a conservative strategy:
- always include the diagonal,
- add couplings numerically detected over the probe states returned by
  `_dae_jacobian_probe_states` (see there for why several, non-uniform states
  are required).

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

    rows = Int[]
    cols = Int[]

    # 1) Always include the diagonal.
    for i in 1:n
        push!(rows, i)
        push!(cols, i)
    end

    # 2) Numerically probe couplings by central finite differences at every probe state.
    # `fd_eps = cbrt(eps)` matches the step formula used later by `_dae_jacobian_fd!`
    # (the actual per-solve evaluator), so a coupling detected here is evaluated
    # with the same sensitivity there — no entry can be "seen" by the prototype
    # and then silently fall under the runtime FD noise floor, or vice versa.
    fd_eps           = cbrt(eps(Float64))
    sensitivity_atol = 1e-14

    y_work              = Vector{Float64}(undef, n)
    res_probe           = Vector{Float64}(undef, n)
    res_perturbed_plus  = Vector{Float64}(undef, n)
    res_perturbed_minus = Vector{Float64}(undef, n)

    probe_states = _dae_jacobian_probe_states(initial_solver_values)
    for (probe_index, (y_probe, local_rtol)) in enumerate(probe_states)
        # Pass A (probe_index == 1) is y0 itself: if the residual cannot be
        # evaluated there, the model is broken and must fail loudly rather than
        # silently degrade to a diagonal-only (hence useless) prototype.
        # Pass B/C states are synthetic and can wander outside the domain of
        # validity of a physical correlation (they throw on non-physical
        # inputs, e.g. a log of a negative concentration); such a probe is
        # simply skipped — it contributes nothing, the other probes still run.
        evaluable = if probe_index == 1
            residual!(res_probe, initial_solver_derivatives, y_probe, packed, t0)
            all(isfinite, res_probe) ||
                throw(ArgumentError("Non-finite DAE residual at the initial state in " *
                                    "_build_dae_jacobian_prototype."))
            true
        else
            try
                residual!(res_probe, initial_solver_derivatives, y_probe, packed, t0)
                all(isfinite, res_probe)
            catch
                false
            end
        end
        evaluable || continue

        # Central-FD sweep: perturb one solver variable j at a time and record
        # every row i with a non-negligible response, for this probe state.
        for j in 1:n
            yj    = y_probe[j]
            delta = fd_eps * max(abs(yj), 1.0)  # adaptive step, scale-invariant near yj = 0.

            copyto!(y_work, y_probe)
            y_work[j] = yj + delta
            try
                residual!(res_perturbed_plus, initial_solver_derivatives, y_work, packed, t0)
            catch
                continue  # column j non-evaluable at this probe: skip it, keep scanning.
            end

            y_work[j] = yj - delta
            try
                residual!(res_perturbed_minus, initial_solver_derivatives, y_work, packed, t0)
            catch
                continue
            end

            inv_2delta = 0.5 / delta
            @inbounds for i in 1:n
                i == j && continue  # Diagonal already added; skip.
                sensitivity = (res_perturbed_plus[i] - res_perturbed_minus[i]) * inv_2delta
                isfinite(sensitivity) || continue
                # Mixed threshold: `local_rtol` is 0 for Pass B/C (absolute-only —
                # see `_dae_jacobian_probe_states`), non-zero for Pass A. `local_scale`
                # normalises by the residual's own magnitude so a genuinely small
                # sensitivity in a large-residual row isn't mistaken for noise.
                local_scale = max(abs(res_probe[i]), abs(res_perturbed_plus[i]),
                                  abs(res_perturbed_minus[i]), 1.0)
                abs(sensitivity) > (sensitivity_atol + local_rtol * local_scale) || continue
                push!(rows, i)
                push!(cols, j)
            end
        end
    end

    # 3) Union of every probe's detections (entries repeated across probes collapse
    # via `sparse`'s default summing combiner; values themselves are unused, only
    # the nonzero structure matters downstream).
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
