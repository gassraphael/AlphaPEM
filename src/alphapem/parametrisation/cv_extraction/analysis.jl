# -*- coding: utf-8 -*-

"""
    CVExtraction analysis

Extract physical parameters (ECSA, H₂ crossover, double-layer capacitance,
ohmic-drop slope) from a corrected cyclic voltammetry cycle.
"""

using LinearAlgebra
using Statistics: mean

"""
    _polyfit(x::Vector{Float64}, y::Vector{Float64}, degree::Int) -> Vector{Float64}

Fit a polynomial of the given degree to `(x, y)` using least squares.
Coefficients are returned in descending order (highest degree first).
"""
function _polyfit(x::Vector{Float64}, y::Vector{Float64}, degree::Int)::Vector{Float64}
    n = length(x)
    A = [x[i]^j for i in 1:n, j in degree:-1:0]
    return A \ y
end

"""
    _polyval(p::Vector{Float64}, x) -> Union{Float64, Vector{Float64}}

Evaluate a polynomial with coefficients `p` (descending order) at `x`.
"""
function _polyval(p::Vector{Float64}, x::Real)::Float64
    y = zero(Float64)
    for c in p
        y = y * x + c
    end
    return y
end

function _polyval(p::Vector{Float64}, x::AbstractVector{<:Real})::Vector{Float64}
    return [_polyval(p, xi) for xi in x]
end

"""
    _polyder(p::Vector{Float64}) -> Vector{Float64}

Return the derivative of polynomial `p` (descending order).
"""
function _polyder(p::Vector{Float64})::Vector{Float64}
    n = length(p)
    if n <= 1
        return [0.0]
    end
    deg = n - 1
    return [p[i] * (deg - i + 1) for i in 1:(n - 1)]
end

"""
    smallest_distance(U::Vector{Float64}, j::Vector{Float64}, U_min::Float64, U_max::Float64) -> Tuple{Vector{Float64}, Vector{Float64}}

Find the smallest distance between the anodic and cathodic branches inside the
potential window `[U_min, U_max]`. Returns two points `P1` (anodic) and `P2`
(cathodic) in `(U, j)` coordinates.
"""
function smallest_distance(
    U::Vector{Float64},
    j::Vector{Float64},
    U_min::Float64,
    U_max::Float64,
)::Tuple{Vector{Float64}, Vector{Float64}}
    a = min(U_min, U_max)
    b = max(U_min, U_max)

    # Split at maximum potential
    split_idx = argmax(U)
    Uf = U[1:split_idx]
    jf = j[1:split_idx]
    Ug = reverse(U[(split_idx + 1):end])
    jg = reverse(j[(split_idx + 1):end])

    # Trim to window
    idx_f = [findfirst(>=(a), Uf), findfirst(>=(b), Uf)]
    idx_g = [findfirst(>=(a), Ug), findfirst(>=(b), Ug)]

    if any(isnothing, idx_f) || any(isnothing, idx_g)
        error("Potential window [$a, $b] is outside the CV data range")
    end

    xf = Uf[idx_f[1]:idx_f[2]]
    yf = jf[idx_f[1]:idx_f[2]]
    xg = Ug[idx_g[1]:idx_g[2]]
    yg = jg[idx_g[1]:idx_g[2]]

    pf = _polyfit(xf, yf, 4)
    pg = _polyfit(xg, yg, 4)

    f(x) = _polyval(pf, x)
    g(x) = _polyval(pg, x)
    df(x) = _polyval(_polyder(pf), x)
    dg(x) = _polyval(_polyder(pg), x)

    # Minimise F(x1, x2) = (x1 - x2)^2 + (f(x1) - g(x2))^2
    function F(x::Vector{Float64})::Float64
        x1, x2 = x
        return (x1 - x2)^2 + (f(x1) - g(x2))^2
    end

    function gradF(x::Vector{Float64})::Vector{Float64}
        x1, x2 = x
        dx1 = 2 * (x1 - x2) + 2 * (f(x1) - g(x2)) * df(x1)
        dx2 = -2 * (x1 - x2) - 2 * (f(x1) - g(x2)) * dg(x2)
        return [dx1, dx2]
    end

    # Gradient descent with Armijo backtracking
    x = [0.4, 0.4]
    tol = 1e-7
    sigma = 0.01
    beta = 0.2

    while norm(gradF(x)) >= tol
        d = -gradF(x)
        alpha = 1.0
        while F(x + alpha * d) > F(x) + sigma * alpha * dot(gradF(x), d)
            alpha *= beta
            if alpha < 1e-12
                break
            end
        end
        x .= x + alpha * d
    end

    U_dlc = mean(x)
    P1 = [U_dlc, f(U_dlc)]
    P2 = [U_dlc, g(U_dlc)]

    return P1, P2
end

"""
    ohmic_drop_fit(U::Vector{Float64}, j::Vector{Float64}, U_min::Float64, U_max::Float64) -> Tuple{Float64, Float64}

Linear fit `j = m * U + c` inside the potential window. Returns `(m, c)`.
`c` is interpreted as the H₂ crossover current density.
"""
function ohmic_drop_fit(
    U::Vector{Float64},
    j::Vector{Float64},
    U_min::Float64,
    U_max::Float64,
)::Tuple{Float64, Float64}
    a = min(U_min, U_max)
    b = max(U_min, U_max)
    idx = findall(u -> a < u < b, U)
    if length(idx) < 2
        error("Not enough data points in ohmic-drop window [$a, $b]")
    end
    p = _polyfit(U[idx], j[idx], 1)
    return p[1], p[2]
end

"""
    area_integral(x::Vector{Float64}, y::Vector{Float64}) -> Float64

Compute the area enclosed between `y` and its baseline, using the trapezoidal
rule.

The baseline is the horizontal line at `y[1]`. This matches Matlab's
`AreaIntegral`, where `polyfit(x(1, end), y(1, end), 1)` degenerates to a single
sample (the comma indexes a column vector as a 2-D array, so `x(1, end)` is
`x(1)`) and the rank-deficient fit collapses to the constant `y(1)`. Callers
prepend the peak's closing sample so that `y[1]` carries the double-layer level.
"""
function area_integral(x::Vector{Float64}, y::Vector{Float64})::Float64
    if length(x) < 2
        return 0.0
    end
    y_abs = abs.(y .- first(y))
    return _trapz(x, y_abs)
end

function _trapz(x::Vector{Float64}, y::Vector{Float64})::Float64
    n = length(x)
    if n != length(y) || n < 2
        return 0.0
    end
    s = 0.0
    for i in 2:n
        s += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1])
    end
    return s
end

"""
    extract_parameters(cycle::CVCycle, cfg::CVExtractionConfig, cv_file_name::String = "") -> CVExtractionResult

Extract ECSA, H₂ crossover, double-layer capacitance and ohmic-drop slope from
a single CV cycle. The optional `cv_file_name` is stored in the result for
identification.
"""
function extract_parameters(cycle::CVCycle, cfg::CVExtractionConfig, cv_file_name::String = "")::CVExtractionResult
    U = cycle.U
    j = cycle.j

    # 1. Smallest distance in double-layer region
    P1, P2 = smallest_distance(
        U,
        j,
        cfg.double_layer_limit_min,
        cfg.double_layer_limit_max,
    )

    # 2. Ohmic-drop fit
    m, c = ohmic_drop_fit(
        U,
        j,
        cfg.ohmic_drop_limit_min,
        cfg.ohmic_drop_limit_max,
    )

    crossover = c
    slope = m .* U
    j_baseline = j .- slope .- crossover

    # Optional compensations
    if cfg.compensate_ohmic_drop
        j = j .- slope
        P1[2] -= P1[1] * m
        P2[2] -= P2[1] * m
    end
    if cfg.compensate_crossover
        j = j .- crossover
        P1[2] -= crossover
        P2[2] -= crossover
    end

    # Rebuild corrected cycle
    corrected = CVCycle(cycle.t, U, cycle.I, j)
    raw = CVCycle(cycle.t, U, cycle.I, cycle.j)

    # 3. Split anodic / cathodic branches
    dU = diff(U)
    anodic_idx = findall(>=(0), dU)
    cathodic_idx = findall(<=(0), dU)

    anodic = _sort_cycle(CVCycle(cycle.t[anodic_idx], U[anodic_idx], cycle.I[anodic_idx], j[anodic_idx]))
    cathodic = _sort_cycle(CVCycle(cycle.t[cathodic_idx], U[cathodic_idx], cycle.I[cathodic_idx], j[cathodic_idx]))

    # 4. Double-layer capacitance
    id_anodic_dlc = argmin(abs.(anodic.U .- P1[1]))
    id_cathodic_dlc = argmin(abs.(cathodic.U .- P2[1]))

    dlc_current = 0.5 * abs(anodic.I[id_anodic_dlc] - cathodic.I[id_cathodic_dlc]) / cfg.area_cm2
    scan_rate = cfg.scan_rate_vs > 0.0 ? cfg.scan_rate_vs : _infer_scan_rate(U, cycle.t)
    dlc = dlc_current / scan_rate

    # 5. ECSA from H₂ adsorption / desorption
    dlc_min_j = min(anodic.j[id_anodic_dlc], cathodic.j[id_cathodic_dlc])
    dlc_max_j = max(anodic.j[id_anodic_dlc], cathodic.j[id_cathodic_dlc])
    dlc_min_U = min(anodic.U[id_anodic_dlc], cathodic.U[id_cathodic_dlc])

    ecsa_max_U = min(cfg.ecsa_limit_max, dlc_min_U)
    ecsa_limits = (min(cfg.ecsa_limit_min, ecsa_max_U), ecsa_max_U)

    conversion_factor = cfg.conversion_factor_uc_cm2 * 1e-6  # μC/cm² -> C/cm²

    ecsa_adsorption = _ecsa_from_branch(
        cathodic.U,
        cathodic.j,
        ecsa_limits,
        dlc_min_j,
        scan_rate,
        conversion_factor,
        cfg.covering_degree,
        <=,
    )

    ecsa_desorption = _ecsa_from_branch(
        anodic.U,
        anodic.j,
        ecsa_limits,
        dlc_max_j,
        scan_rate,
        conversion_factor,
        cfg.covering_degree,
        >=,
    )

    return CVExtractionResult(
        cv_file_name,
        ecsa_adsorption,
        ecsa_desorption,
        crossover,
        dlc,
        m,
        scan_rate,
        corrected,
        raw,
        anodic,
        cathodic,
    )
end

function _sort_cycle(cycle::CVCycle)::CVCycle
    order = sortperm(cycle.U)
    return CVCycle(
        cycle.t[order],
        cycle.U[order],
        cycle.I[order],
        cycle.j[order],
    )
end

function _infer_scan_rate(U::Vector{Float64}, t::Vector{Float64})::Float64
    return mean(abs.(diff(U))) / mean(diff(t))
end

function _ecsa_from_branch(
    U::Vector{Float64},
    j::Vector{Float64},
    limits::Tuple{Float64, Float64},
    dlc_j::Float64,
    scan_rate::Float64,
    conversion_factor::Float64,
    covering_degree::Float64,
    cmp_op,
)::Float64
    a, b = limits
    lo = min(a, b)
    hi = max(a, b)
    idx = Int[]
    for i in eachindex(U)
        u = U[i]
        if lo < u < hi && cmp_op(j[i], dlc_j)
            push!(idx, i)
        end
    end
    if length(idx) < 2
        return 0.0
    end
    U_peak = U[idx]
    j_peak = j[idx]

    # Close the peak with the sample at its high-potential end, which sets the
    # baseline level in `area_integral` (Matlab's `AddLine`). The prepended
    # sample repeats the first potential, so it adds no area itself.
    U_closed = [first(U_peak); U_peak]
    j_closed = [last(j_peak); j_peak]

    area = area_integral(U_closed, j_closed)
    q = area / scan_rate
    return q / conversion_factor * covering_degree
end
