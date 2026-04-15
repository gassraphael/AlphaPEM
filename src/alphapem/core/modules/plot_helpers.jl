# plot_helpers.jl
#
# Helper utilities shared by CairoMakie plotting functions.

module PlotHelpers
using CairoMakie

export _publication_colors,
       _finalize_axis!,
       _format_fixed,
       _compact_tick_labels,
       _nice_tick_step,
       _colorbar_ticks_auto,
       _rounded_major_ticks,
       _set_dense_ticks!,
       lsub

"""Return a reproducible publication-oriented color set."""
function _publication_colors()
    return [
        RGBf(0.00, 0.45, 0.70),
        RGBf(0.84, 0.37, 0.00),
        RGBf(0.00, 0.62, 0.45),
        RGBf(0.80, 0.47, 0.65),
        RGBf(0.34, 0.71, 0.91),
        RGBf(0.90, 0.62, 0.00),
        RGBf(0.53, 0.13, 0.33),
        RGBf(0.27, 0.67, 0.60),
        RGBf(0.60, 0.60, 0.20),
        RGBf(0.20, 0.13, 0.53),
        RGBf(0.86, 0.80, 0.47),
    ]
end

"""Apply common axis formatting for scientific figures."""
function _finalize_axis!(ax;
                         xlabel="",
                         ylabel="",
                         title="",
                         legend::Bool=false,
                         legend_position=:rb)
    # Always convert labels to rich text for consistent bold rendering.
    ax.xlabel = rich(xlabel; font=:bold)
    ax.ylabel = rich(ylabel; font=:bold)
    isempty(title) || (ax.title = title)

    # Dense major/minor ticks with differentiated sizes for publication readability.
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xminorgridvisible = false
    ax.yminorgridvisible = false
    ax.xminorticks = IntervalsBetween(5)
    ax.yminorticks = IntervalsBetween(5)
    ax.xminorticksvisible = true
    ax.yminorticksvisible = true
    ax.xticksize = 13
    ax.yticksize = 13
    ax.xminorticksize = 5
    ax.yminorticksize = 5
    ax.xticklabelsize = 15
    ax.yticklabelsize = 15

    if legend
        axislegend(ax; position=legend_position, framevisible=true,
                   framewidth=0.8, framecolor=(:black, 0.35),
                   backgroundcolor=(:white, 0.88), padding=(6, 6, 6, 6))
    end
    return nothing
end


"""Format a single real value with exactly `ndigits` decimal places, including trailing zeros."""
function _format_fixed(v::Real, ndigits::Int)::String
    sign = v < 0 ? "-" : ""
    abs_v = abs(float(v))
    factor = 10.0^ndigits
    total = round(Int64, abs_v * factor)
    int_part = div(total, round(Int64, factor))
    dec_part = rem(total, round(Int64, factor))
    return "$(sign)$(int_part).$(lpad(dec_part, ndigits, '0'))"
end

"""Format tick labels uniformly across all values on an axis.
- If all values are integers: display without decimal point (e.g. 1700, 1750).
- Otherwise: use the minimum number of decimal places needed to make all
  rounded values unique, then pad every label with trailing zeros so that
  all labels share the same number of decimals (e.g. 0.68, 0.70 not 0.68, 0.7).
"""
function _compact_tick_labels(values::AbstractVector{<:Real})
    isempty(values) && return String[]

    # All values effectively integers -> no decimal point.
    if all(v -> isapprox(v, round(v); atol=1e-9 * max(abs(v), 1.0)), values)
        return [string(round(Int, v)) for v in values]
    end

    # Find the minimum number of decimal places that makes all rounded values unique.
    ndigits = 1
    while ndigits <= 8
        allunique(round.(values; digits=ndigits)) && break
        ndigits += 1
    end

    # Format every label with exactly ndigits decimal places (trailing zeros included).
    return [_format_fixed(v, ndigits) for v in values]
end


"""Return a readable major-step size using the 1-2-5 rule."""
function _nice_tick_step(span::Real, n_major::Int)::Float64
    target = max(float(span) / max(n_major - 1, 1), eps(Float64))
    exponent = floor(log10(target))
    fraction = target / 10.0^exponent
    nice_fraction = fraction < 1.5 ? 1.0 : (fraction < 3.0 ? 2.0 : (fraction < 7.0 ? 5.0 : 10.0))
    return nice_fraction * 10.0^exponent
end


"""Build readable colorbar ticks from automatic min/max values."""
function _colorbar_ticks_auto(vmin::Real, vmax::Real)::Vector{Float64}
    span = float(vmax - vmin)
    span <= 0 && return [float(vmin), float(vmax)]

    step = span <= 20 ? 1.0 : _nice_tick_step(span, 6)
    tick_start = ceil(float(vmin) / step) * step
    tick_end = ceil(float(vmax) / step) * step
    ticks = collect(tick_start:step:tick_end)
    return isempty(ticks) ? collect(range(float(vmin), float(vmax); length=5)) : ticks
end


"""Build rounded major ticks for time-like axes (e.g. 1700, 1750, ...)."""
function _rounded_major_ticks(xmin::Real, xmax::Real; n_major::Int=7)
    span = float(xmax - xmin)
    span <= 0 && return [float(xmin)]

    step = _nice_tick_step(span, n_major)
    start_tick = ceil(float(xmin) / step) * step
    end_tick = floor(float(xmax) / step) * step
    ticks = collect(start_tick:step:end_tick)

    # If the first pass is too sparse, refine once while keeping rounded values.
    if length(ticks) < 4
        step = _nice_tick_step(span, max(2 * n_major - 1, 3))
        start_tick = ceil(float(xmin) / step) * step
        end_tick = floor(float(xmax) / step) * step
        ticks = collect(start_tick:step:end_tick)
    end

    return isempty(ticks) ? collect(range(float(xmin), float(xmax); length=n_major)) : ticks
end


"""Force dense major ticks (at least 6) on both axes from plotted data ranges."""
function _set_dense_ticks!(ax,
                           x::AbstractVector{<:Real},
                           ys::AbstractVector{<:AbstractVector{<:Real}};
                           n_major::Int=7)
    xfinite = filter(isfinite, collect(x))
    yfinite = filter(isfinite, vcat((collect(y) for y in ys)...))
    if !isempty(xfinite)
        xmin, xmax = extrema(xfinite)
        if xmin == xmax
            delta = max(abs(xmin), 1.0) * 0.05
            xmin -= delta
            xmax += delta
        end
        ax.xticks = _rounded_major_ticks(xmin, xmax; n_major=n_major)
        ax.xtickformat = _compact_tick_labels
    end
    if !isempty(yfinite)
        ymin, ymax = extrema(yfinite)
        if ymin == ymax
            delta = max(abs(ymin), 1.0) * 0.05
            ymin -= delta
            ymax += delta
        end
        ax.yticks = _rounded_major_ticks(ymin, ymax; n_major=n_major)
        ax.ytickformat = _compact_tick_labels
    end
    return nothing
end


"""Create a visual subscript label, e.g. lsub("T", "agc") -> T_agc (rendered with subscript)."""
lsub(base::AbstractString, idx::AbstractString) = rich(base, subscript(idx))

end # module PlotHelpers


