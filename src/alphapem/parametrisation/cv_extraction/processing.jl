# -*- coding: utf-8 -*-

"""
    CVExtraction processing

Clean, resample, split, and average cyclic voltammetry cycles.
"""

using Statistics: mean, median
using Interpolations: linear_interpolation, Flat

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

    # Remove consecutive duplicate potentials
    dup = findall(diff(U) .== 0)
    keep = setdiff(1:length(U), dup .+ 1)
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
    _local_extrema(x::Vector{Float64}; window::Int = 5) -> Tuple{Vector{Int}, Vector{Int}}

Return indices of local minima and local maxima using a sliding window.
"""
function _local_extrema(x::Vector{Float64}; window::Int = 5)::Tuple{Vector{Int}, Vector{Int}}
    n = length(x)
    minima = Int[]
    maxima = Int[]
    half = max(1, window ÷ 2)

    for i in (half + 1):(n - half)
        neighborhood = x[(i - half):(i + half)]
        if x[i] == minimum(neighborhood) && x[i] < x[i - 1] && x[i] < x[i + 1]
            push!(minima, i)
        elseif x[i] == maximum(neighborhood) && x[i] > x[i - 1] && x[i] > x[i + 1]
            push!(maxima, i)
        end
    end

    return minima, maxima
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

Upsample a single cycle using spline-like linear interpolation.
"""
function interpolate_cycle(data::CVData, scan_rate::Float64, factor::Int)::CVData
    dt = scan_rate / factor
    t_new = range(data.t[1], data.t[end]; step = dt)

    itp_U = linear_interpolation(data.t, data.U; extrapolation_bc = Flat())
    itp_I = linear_interpolation(data.t, data.I; extrapolation_bc = Flat())

    U_new = itp_U(t_new)
    I_new = itp_I(t_new)

    return CVData(collect(t_new), U_new, I_new)
end

"""
    mean_cycle(cycles::Vector{CVData}, ignore::Vector{Int} = Int[]) -> CVData

Compute a mean cycle from a list of cycles, optionally ignoring some cycle indices.
"""
function mean_cycle(cycles::Vector{CVData}, ignore::Vector{Int} = Int[])::CVData
    n = length(cycles)
    use = setdiff(1:n, ignore)
    if isempty(use)
        use = 1:n
    end

    # Concatenate selected cycles
    t_all = vcat([cycles[i].t for i in use]...)
    U_all = vcat([cycles[i].U for i in use]...)
    I_all = vcat([cycles[i].I for i in use]...)

    # Split into rising and falling potential branches
    dU = diff(U_all)
    rising_idx = findall(>(0), dU)
    falling_idx = findall(<(0), dU)

    rising = sort(collect(zip(U_all[rising_idx], I_all[rising_idx])))
    falling = sort(collect(zip(U_all[falling_idx], I_all[falling_idx])); rev = true)

    Ur = first.(rising)
    Ir = last.(rising)
    Uf = first.(falling)
    If = last.(falling)

    # Light smoothing (moving average, window = 50)
    Ir = _moving_average(Ir, 50)
    If = _moving_average(If, 50)

    U_mean = vcat(Ur, Uf)
    I_mean = vcat(Ir, If)
    t_mean = range(0.0, 1.0; length = length(U_mean))

    return CVData(collect(t_mean), U_mean, I_mean)
end

function _moving_average(x::Vector{Float64}, window::Int)::Vector{Float64}
    w = max(1, window)
    n = length(x)
    y = similar(x)
    half = w ÷ 2
    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        y[i] = mean(x[lo:hi])
    end
    return y
end
