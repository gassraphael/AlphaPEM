# -*- coding: utf-8 -*-

"""
    CVExtraction processing

Clean, resample, split, and average cyclic voltammetry cycles.
"""

using Statistics: mean, median
using Interpolations: linear_interpolation, interpolate, scale, BSpline, Cubic, Free, OnGrid, Flat

"""
    clean_cv_data(data::CVData) -> CVData

Remove duplicate time stamps and duplicate potential values.
"""
function clean_cv_data(data::CVData)::CVData
    t = data.t
    U = data.U
    I = data.I

    # Unique time stamps (keep first occurrence)
    seen = Set{Float64}()
    keep = Int[]
    for i in eachindex(t)
        if !(t[i] in seen)
            push!(seen, t[i])
            push!(keep, i)
        end
    end
    t = t[keep]
    U = U[keep]
    I = I[keep]

    # Remove consecutive duplicate potentials. Matlab drops the *first* element
    # of each equal pair (`DataTable(find(diff(U) == 0), :) = []`), which is what
    # is reproduced here: dropping the second one instead shifts the retained
    # samples by one and perturbs the inferred scan rate.
    dup = findall(diff(U) .== 0)
    keep = setdiff(1:length(U), dup)
    t = t[keep]
    U = U[keep]
    I = I[keep]

    return CVData(t, U, I)
end

"""
    resample_cv_data(data::CVData; scan_rate::Float64 = -1.0) -> Tuple{CVData, Float64}

Resample the CV data on an equidistant time grid and return the resampled data
plus the inferred scan rate in V/s.
"""
function resample_cv_data(data::CVData; scan_rate::Float64 = -1.0)::Tuple{CVData, Float64}
    t = data.t
    U = data.U
    I = data.I

    if scan_rate <= 0.0
        scan_rate = mean(abs.(diff(U))) / mean(diff(t))
    end

    dt = scan_rate
    t_new = range(t[1], t[end]; step = dt)

    itp_U = linear_interpolation(t, U; extrapolation_bc = Flat())
    itp_I = linear_interpolation(t, I; extrapolation_bc = Flat())

    U_new = itp_U(t_new)
    I_new = itp_I(t_new)

    return CVData(collect(t_new), U_new, I_new), scan_rate
end

"""
    _strict_local_extrema(x::Vector{Float64}, mode::Symbol) -> Vector{Int}

Return the indices of the strict local extrema of `x`, with `mode` either
`:min` or `:max`. This mirrors Matlab's `islocalmin` / `islocalmax` with default
options: a point qualifies when it is strictly beyond both of its neighbours,
and a flat plateau is reported at its centre index.
"""
function _strict_local_extrema(x::Vector{Float64}, mode::Symbol)::Vector{Int}
    s = mode === :max ? 1.0 : -1.0
    y = s .* x
    n = length(y)
    idx = Int[]

    i = 2
    while i <= n - 1
        if y[i] > y[i - 1]
            # Walk to the end of a possible plateau at this level
            j = i
            while j + 1 <= n && y[j + 1] == y[i]
                j += 1
            end
            if j + 1 <= n && y[j + 1] < y[i]
                push!(idx, (i + j) ÷ 2)
            end
            i = j + 1
        else
            i += 1
        end
    end

    return idx
end

"""
    _local_extrema(x::Vector{Float64}) -> Tuple{Vector{Int}, Vector{Int}}

Return the indices of the local minima and local maxima of `x`.
"""
function _local_extrema(x::Vector{Float64})::Tuple{Vector{Int}, Vector{Int}}
    return _strict_local_extrema(x, :min), _strict_local_extrema(x, :max)
end

"""
    split_cycles(data::CVData) -> Vector{CVData}

Split a CV dataset into individual cycles based on local minima of the potential.
"""
function split_cycles(data::CVData)::Vector{CVData}
    minima, maxima = _local_extrema(data.U)

    # Boundaries include the first and last data points
    boundaries = [1; minima; length(data.U)]
    sort!(boundaries)
    unique!(boundaries)

    if length(boundaries) < 3
        return [data]
    end

    cycle_lengths = diff(boundaries)
    min_cycle_length = (1 / 3) * median(cycle_lengths)

    cycles = CVData[]
    for i in 1:(length(boundaries) - 1)
        range = boundaries[i]:boundaries[i + 1]
        if length(range) < min_cycle_length
            continue
        end
        # Keep only cycles that contain a local maximum
        if any(m -> boundaries[i] <= m <= boundaries[i + 1], maxima)
            push!(cycles, CVData(data.t[range], data.U[range], data.I[range]))
        end
    end

    return cycles
end

"""
    interpolate_cycle(data::CVData, scan_rate::Float64, factor::Int) -> CVData

Upsample a single cycle onto a time grid `factor` times finer, using a cubic
spline. Matlab calls `interp1(..., 'spline')` here, i.e. a not-a-knot cubic
spline; a linear interpolation flattens the ripple of noisy CV data and biases
the peak integrals.
"""
function interpolate_cycle(data::CVData, scan_rate::Float64, factor::Int)::CVData
    dt = scan_rate / factor
    t_new = range(data.t[1], data.t[end]; step = dt)

    U_new = _spline_resample(data.t, data.U, t_new)
    I_new = _spline_resample(data.t, data.I, t_new)

    return CVData(collect(t_new), U_new, I_new)
end

"""
    _spline_resample(t::Vector{Float64}, y::Vector{Float64}, t_new) -> Vector{Float64}

Evaluate a not-a-knot cubic spline through `(t, y)` at `t_new`. `t` is assumed
to be an equidistant grid, which is the case for cycles cut out of the
resampled dataset. Falls back to linear interpolation when there are too few
points for a cubic spline.
"""
function _spline_resample(t::Vector{Float64}, y::Vector{Float64}, t_new)::Vector{Float64}
    if length(t) < 4
        itp = linear_interpolation(t, y; extrapolation_bc = Flat())
        return [itp(x) for x in t_new]
    end
    grid = range(t[1], t[end]; length = length(t))
    itp = scale(interpolate(y, BSpline(Cubic(Free(OnGrid())))), grid)
    lo, hi = first(grid), last(grid)
    return [itp(clamp(x, lo, hi)) for x in t_new]
end

"""
    mean_cycle(cycles::Vector{CVData}, ignore::Vector{Int} = Int[]) -> CVData

Build the mean cycle from a list of cycles, optionally ignoring some cycle
indices (1-based, as displayed to the user).

The selected cycles are concatenated, split into a rising (anodic) and a falling
(cathodic) branch, each branch is sorted by potential and smoothed, and the two
branches are joined back into a single anodic-then-cathodic sweep. This is the
`Create Mean Cycle` block of Matlab's `ImportCV.m`; note that the cycles fed in
must be the *resampled* ones, not the finely interpolated ones, because Matlab
builds the mean cycle before its per-cycle interpolation loop.
"""
function mean_cycle(cycles::Vector{CVData}, ignore::Vector{Int} = Int[])::CVData
    n = length(cycles)
    use = setdiff(1:n, ignore)
    if isempty(use)
        use = collect(1:n)
    end

    # Concatenate selected cycles
    t_all = vcat([cycles[i].t for i in use]...)
    U_all = vcat([cycles[i].U for i in use]...)
    I_all = vcat([cycles[i].I for i in use]...)

    # Split into rising and falling potential branches
    dU = diff(U_all)
    rising_idx = findall(>(0), dU)
    falling_idx = findall(<(0), dU)

    # Sort by potential only. Matlab's `sortrows(..., {'U'})` is stable, so
    # samples sharing a potential keep their acquisition order; sorting
    # (U, I) pairs instead would reorder them by current and distort the branch.
    rising_idx = rising_idx[sortperm(U_all[rising_idx]; alg = MergeSort)]
    falling_idx = falling_idx[sortperm(U_all[falling_idx]; alg = MergeSort, rev = true)]

    Ir = _matlab_smooth_interior(I_all[rising_idx], 50)
    If = _matlab_smooth_interior(I_all[falling_idx], 50)

    U_mean = vcat(U_all[rising_idx], U_all[falling_idx])
    I_mean = vcat(Ir, If)
    t_mean = vcat(t_all[rising_idx], t_all[falling_idx])

    return CVData(t_mean, U_mean, I_mean)
end

"""
    _matlab_smooth_interior(x::Vector{Float64}, span::Int) -> Vector{Float64}

Apply [`_matlab_smooth`](@ref) to `x[2:end-2]` and leave the remaining samples
untouched, as Matlab does with `A.I(2:end-2) = smooth(A.I(2:end-2), 50)`.
"""
function _matlab_smooth_interior(x::Vector{Float64}, span::Int)::Vector{Float64}
    y = copy(x)
    if length(y) >= 4
        y[2:(end - 2)] = _matlab_smooth(y[2:(end - 2)], span)
    end
    return y
end

"""
    _matlab_smooth(x::Vector{Float64}, span::Int) -> Vector{Float64}

Moving average reproducing Matlab's `smooth(x, span)`: the span is reduced to
the next odd value, and the window shrinks symmetrically near both ends instead
of being truncated on one side only.
"""
function _matlab_smooth(x::Vector{Float64}, span::Int)::Vector{Float64}
    n = length(x)
    w = min(span, n)
    if iseven(w)
        w -= 1
    end
    if w < 3
        return copy(x)
    end

    cumulative = cumsum(vcat(0.0, x))
    y = similar(x)
    for i in 1:n
        width = min(w, 2 * i - 1, 2 * (n - i) + 1)
        half = width ÷ 2
        y[i] = (cumulative[i + half + 1] - cumulative[i - half]) / width
    end
    return y
end
