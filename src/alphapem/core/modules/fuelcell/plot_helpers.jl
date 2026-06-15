# plot_helpers.jl
#
# Helper utilities shared by CairoMakie plotting functions.

module PlotHelpers
using CairoMakie

export _publication_colors,
       _macroscopic_color,
       _cell_current_color,
       _cell_voltage_color,
       _add_axis_toolbar!,
       _finalize_axis!,
       _label_with_rmse,
       _polarization_legend_base,
       _experimental_marker,
       _set_polarization_xlims!,
       _set_polarization_axis_limits!,
       _set_polarization_fixed_ticks!,
       _set_initial_temporal_limits!,
       _format_fixed,
       _compact_tick_labels,
       _nice_tick_step,
       _colorbar_ticks_auto,
       _rounded_major_ticks,
       _set_dense_ticks!,
       _plot_final_profile_along_gc!,
       _clear_dynamic_axes!,
       _clear_dynamic_legends!,
       saving_instructions!,
       gc_direction_labels,
       lsub

# Internal storage for axis metadata that cannot be attached directly to Makie objects.
const _axis_initial_limits = WeakKeyDict{Any, Observable{Any}}()

"""Helper to get or create the initial limits observable for an axis."""
function _get_initial_limits_obs(ax)
    return get!(_axis_initial_limits, ax) do
        Observable{Any}(nothing)
    end
end

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

"""Return a dedicated color for cell-wide (macroscopic) variables like current or voltage.
We use a vibrant but professional Deep Blue that is distinct from regional colors.
"""
function _macroscopic_color()
    return RGBf(0.1, 0.4, 0.7)
end

"""Return a dedicated color for the cell-wide current density."""
function _cell_current_color(nb_gc::Int)
    return _macroscopic_color()
end

"""Return a dedicated color for the cell-wide voltage, consistent across AlphaPEM."""
function _cell_voltage_color()
    return _macroscopic_color()
end

"""Add interactive 'Reset View' and 'Zoom' buttons to a specific axis."""
function _add_axis_toolbar!(ax)
    # Don't add toolbar for static backends (CairoMakie)
    # This prevents buttons from appearing in PDF/PNG exports.
    if string(Makie.current_backend()) == "CairoMakie"
        return nothing
    end

    # 1. Identify parent Figure
    fig = nothing
    curr = ax
    try
        while curr !== nothing && !(curr isa Figure)
            if hasproperty(curr, :parent)
                curr = curr.parent
            else
                break
            end
        end
        fig = curr
    catch
        return nothing
    end
    (fig === nothing || !(fig isa Figure)) && return nothing

    # 2. Identify position of the axis in the layout
    gc = nothing
    for content in fig.layout.content
        if content.content == ax
            gc = content
            break
        end
    end
    gc === nothing && return nothing

    # 3. Check for existing toolbar to avoid duplicates
    # We check if there's already a GridLayout in the Top() position of this cell
    for content in fig.layout.content
        if content.parent == gc.parent && 
           content.span == gc.span && 
           content.side == Top() && 
           content.content isa GridLayout
            return nothing
        end
    end

    # 4. Create the toolbar layout
    # Using Top() places buttons in the decoration area of the axis,
    # aligned with the title.
    gl = gc.parent[gc.span.rows, gc.span.cols, Top()] = GridLayout(
        halign = :right, 
        tellwidth = false, 
        padding = (0, 4, 1, 0)
    )
    
    # 5. Design: Styled background capsule
    Box(gl[1, 1:3], color = (:gray95, 0.7), strokewidth = 0.5, strokecolor = :gray80, cornerradius = 4)
    
    btn_attr = (
        width = 22, 
        height = 22,
        buttoncolor = :transparent,
        buttoncolor_hover = (:black, 0.08),
        buttoncolor_active = (:black, 0.15),
        labelcolor = :black,
        strokewidth = 0,
        cornerradius = 4,
        padding = (0, 0, 0, 0),
        fontsize = 14
    )
    
    btn_in  = Button(gl[1, 1]; label = "+", btn_attr...)
    btn_out = Button(gl[1, 2]; label = "-", btn_attr...)
    btn_res = Button(gl[1, 3]; label = "⌂", btn_attr...)
    
    colgap!(gl, 0)

    # 6. Interaction Logic
    obs = _get_initial_limits_obs(ax)
    on(btn_res.clicks) do _
        if obs[] !== nothing
            limits!(ax, obs[]...)
        else
            autolimits!(ax)
        end
    end

    on(btn_in.clicks) do _
        _zoom_axis!(ax, 0.8)
    end

    on(btn_out.clicks) do _
        _zoom_axis!(ax, 1.25)
    end

    return nothing
end

"""Helper to zoom into an axis (2D or 3D)."""
function _zoom_axis!(ax, factor)
    tname = string(typeof(ax))
    if ax isa Axis || occursin("Axis3", tname)
        lims = ax.finallimits[]
        center = lims.origin + lims.widths / 2
        new_widths = lims.widths * factor
        new_origin = center - new_widths / 2
        
        # Apply limits. Axis3 supports x/y/zlims!
        xlims!(ax, new_origin[1], new_origin[1] + new_widths[1])
        ylims!(ax, new_origin[2], new_origin[2] + new_widths[2])
        if occursin("Axis3", tname) && length(new_origin) >= 3
            zlims!(ax, new_origin[3], new_origin[3] + new_widths[3])
        end
    end
end

"""Apply common axis formatting for scientific figures."""
function _finalize_axis!(ax;
                         xlabel="",
                         ylabel="",
                         title="",
                         legend::Bool=false,
                         legend_position=:rb)
    # Apply labels: avoid wrapping already-rich text in another rich() call,
    # as deep nesting of RichText objects can cause SVG rendering failures.
    ax.xlabel = (xlabel isa Makie.RichText) ? xlabel : rich(xlabel; font=:bold)
    ax.ylabel = (ylabel isa Makie.RichText) ? ylabel : rich(ylabel; font=:bold)
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

    # Add interactive 'Reset View' (house) and zoom buttons to the axis.
    _get_initial_limits_obs(ax)
    _add_axis_toolbar!(ax)

    if legend
        axislegend(ax; position=legend_position, framevisible=true,
                   framewidth=0.8, framecolor=(:black, 0.35),
                   backgroundcolor=(:white, 0.88), padding=(6, 6, 6, 6))
    end
    return nothing
end

"""Append legacy RMSE text to a model legend label when available."""
function _label_with_rmse(base_label, sim_error)
    sim_error === nothing && return base_label
    return rich(base_label, " - ΔU", subscript("RMSE"), " = $(sim_error) %")
end

"""Return the legacy polarization legend base text used in historical Python plots."""
function _polarization_legend_base(type_fuel_cell::Symbol;
                                   simulation::Bool=true,
                                   calibration::Bool=false)
    prefix = simulation ? "Sim. - " : "Exp. - "

    if type_fuel_cell == :ZSW_GenStack
        suffix = calibration && simulation ? "nominal operating conditions" : "nominal"
        return prefix * suffix
    elseif type_fuel_cell == :ZSW_GenStack_Pa_1_61_Pc_1_41
        return rich(prefix, "P", subscript("a"), "/P", subscript("c"), " = 1.61/1.41 bar")
    elseif type_fuel_cell == :ZSW_GenStack_Pa_2_01_Pc_1_81
        return rich(prefix, "P", subscript("a"), "/P", subscript("c"), " = 2.01/1.81 bar")
    elseif type_fuel_cell == :ZSW_GenStack_Pa_2_4_Pc_2_2
        return rich(prefix, "P", subscript("a"), "/P", subscript("c"), " = 2.4/2.2 bar")
    elseif type_fuel_cell == :ZSW_GenStack_Pa_2_8_Pc_2_6
        return rich(prefix, "P", subscript("a"), "/P", subscript("c"), " = 2.8/2.6 bar")
    elseif type_fuel_cell == :ZSW_GenStack_T_62
        return prefix * "T = 62 °C"
    elseif type_fuel_cell == :ZSW_GenStack_T_76
        return prefix * "T = 76 °C"
    elseif type_fuel_cell == :ZSW_GenStack_T_84
        return prefix * "T = 84 °C"
    elseif type_fuel_cell == :EH_31_1_5
        return prefix * "P = 1.5 bar"
    elseif type_fuel_cell == :EH_31_2_0
        return prefix * "P = 2.0 bar"
    elseif type_fuel_cell == :EH_31_2_25
        return prefix * "P = 2.25 bar"
    elseif type_fuel_cell == :EH_31_2_5
        return prefix * "P = 2.5 bar"
    end

    return simulation ? "Simulation" : "Experiment"
end

"""Return an experimental marker style consistent with historical plotting conventions."""
function _experimental_marker(type_fuel_cell::Symbol)::Symbol
    if type_fuel_cell in (:ZSW_GenStack, :EH_31_1_5)
        return :rect
    elseif type_fuel_cell in (:ZSW_GenStack_Pa_1_61_Pc_1_41, :EH_31_2_0)
        return :utriangle
    elseif type_fuel_cell in (:ZSW_GenStack_Pa_2_01_Pc_1_81, :EH_31_2_25)
        return :dtriangle
    elseif type_fuel_cell in (:ZSW_GenStack_Pa_2_4_Pc_2_2, :EH_31_2_5)
        return :pentagon
    elseif type_fuel_cell == :ZSW_GenStack_Pa_2_8_Pc_2_6
        return :diamond
    elseif type_fuel_cell == :ZSW_GenStack_T_62
        return :cross
    elseif type_fuel_cell == :ZSW_GenStack_T_76
        return :xcross
    elseif type_fuel_cell == :ZSW_GenStack_T_84
        return :star5
    end
    return :rect
end

"""Ensure polarization plots keep the OCV point visible at x = 0."""
function _set_polarization_xlims!(ax, x_values::AbstractVector{<:Real})
    xfinite = filter(isfinite, collect(x_values))
    isempty(xfinite) && return nothing

    xmin, xmax = extrema(xfinite)
    xmin = min(0.0, xmin)
    if xmin == xmax
        delta = max(abs(xmin), 1.0) * 0.05
        xmin -= delta
        xmax += delta
    else
        margin = max((xmax - xmin) * 0.03, 1e-6)
        xmin -= margin
        xmax += margin
    end
    xlims!(ax, xmin, xmax)
    return nothing
end

"""Apply fixed axis limits for polarization charts."""
function _set_polarization_axis_limits!(ax)
    xlims!(ax, 0.0, 3.0)
    ylims!(ax, 0.4, 1.2)
    return nothing
end

"""Apply adaptive tick formatting for polarization charts.

Major ticks are placed at a fixed spacing (0.5 A·cm⁻² on x, 0.1 V on y) via
tick *functions* rather than fixed arrays.  Makie calls these functions with the
current view limits `(vmin, vmax)` at every pan/zoom event, so ticks extend
correctly outside the nominal polarization domain when the user moves the window."""
function _set_polarization_fixed_ticks!(ax)
    function _x_ticks(vmin, vmax)
        step  = 0.5
        start = ceil(vmin / step) * step
        stop  = floor(vmax / step) * step
        vals  = start <= stop ? collect(start:step:stop) : Float64[]
        return vals, _compact_tick_labels(vals)
    end
    function _y_ticks(vmin, vmax)
        step  = 0.1
        start = ceil(vmin / step) * step
        stop  = floor(vmax / step) * step
        vals  = start <= stop ? collect(start:step:stop) : Float64[]
        return vals, _compact_tick_labels(vals)
    end
    ax.xticks = _x_ticks
    ax.yticks = _y_ticks
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


"""Force adaptive dense major ticks (at least 6) on both axes.

Ticks are configured as functions that react to pan/zoom events, ensuring
they extend correctly even when the user moves outside the initial data domain.
"""
function _set_dense_ticks!(ax,
                           x::AbstractVector{<:Real},
                           ys::AbstractVector{<:AbstractVector{<:Real}};
                           n_major::Int=7)
    # Define adaptive tick functions using closures that capture n_major.
    function _adaptive_x_ticks(vmin, vmax)
        ticks = _rounded_major_ticks(vmin, vmax; n_major=n_major)
        return ticks, _compact_tick_labels(ticks)
    end
    function _adaptive_y_ticks(vmin, vmax)
        ticks = _rounded_major_ticks(vmin, vmax; n_major=n_major)
        return ticks, _compact_tick_labels(ticks)
    end

    ax.xticks = _adaptive_x_ticks
    ax.yticks = _adaptive_y_ticks
    
    # We also keep minor ticks visible (configured in _finalize_axis!).
    return nothing
end


"""Calculate and apply initial framing for temporal plots.

X limits are set to `t_range`, and Y limits are automatically computed from
data present within that range. This ensures the 'Home' button resets to a
clean, relevant view of the results after initialization.
"""
function _set_initial_temporal_limits!(ax,
                                       t::AbstractVector{<:Real},
                                       y_series::AbstractVector{<:AbstractVector{<:Real}},
                                       t_range::Tuple{Float64, Float64})
    t_start, t_end = t_range
    
    # Calculate Y limits based on data visible in the initial X range.
    mask = (t .>= t_start) .& (t .<= t_end)
    y_visible = Float64[]
    for y in y_series
        append!(y_visible, y[mask])
    end
    
    y_finite = filter(isfinite, y_visible)
    if !isempty(y_finite)
        ymin, ymax = extrema(y_finite)
        if ymin == ymax
            delta = max(abs(ymin), 1.0) * 0.05
            ymin -= delta
            ymax += delta
        else
            margin = (ymax - ymin) * 0.05
            ymin -= margin
            ymax += margin
        end
        
        # Apply limits and store them for the 'Home' button.
        obs = _get_initial_limits_obs(ax)
        obs[] = (t_start, t_end, ymin, ymax)
        limits!(ax, t_start, t_end, ymin, ymax)
    else
        xlims!(ax, t_start, t_end)
    end
    
    return nothing
end


"""Create a visual subscript label, e.g. lsub("T", "agc") -> T_agc (rendered with subscript)."""
function lsub(base::AbstractString, idx::AbstractString)
    # Avoid creating subscript("") which can cause SVG rendering issues.
    return isempty(idx) ? base : rich(base, subscript(idx))
end

"""Plot one final profile along the gas channel with the publication GC style."""
function _plot_final_profile_along_gc!(ax,
                                       x_gc::AbstractVector{<:Real},
                                       y_final::AbstractVector{<:Real},
                                       cfg;
                                       color,
                                       series_label,
                                       ylabel,
                                       reference_value=nothing,
                                       reference_label=nothing,
                                       legend_position=:rt)
    lines!(ax, x_gc, y_final; color=color, linewidth=2.8, label=series_label)
    scatter!(ax, x_gc, y_final; color=color, markersize=8)

    has_reference = reference_value !== nothing && isfinite(reference_value)
    if has_reference
        lines!(ax, [first(x_gc), last(x_gc)], [reference_value, reference_value];
               color=:black, linestyle=:dash, linewidth=2.0,
               label=reference_label)
    end

    y_finite = filter(isfinite, has_reference ? vcat(collect(y_final), [reference_value]) : collect(y_final))
    if !isempty(y_finite)
        ymin, ymax = extrema(y_finite)
        ymin == ymax && (delta = max(abs(ymin), 1.0) * 0.05; ymin -= delta; ymax += delta)
        ax.yticks = _rounded_major_ticks(ymin, ymax)
        ax.ytickformat = _compact_tick_labels
    end

    _finalize_axis!(ax;
                    xlabel=rich("Through gas-channel"; font=:bold),
                    ylabel=ylabel,
                    legend=true,
                    legend_position=legend_position)

    ax.xticks = (x_gc, gc_direction_labels(cfg, length(x_gc); triple_break_cathode_inlet=true))
    ax.xticklabelrotation = π / 10
    ax.xticklabelsize = 10
    ax.xminorticksvisible = false
    return nothing
end

"""Return GC-direction labels consistent with final postrun maps."""
function gc_direction_labels(cfg,
                             nb_gc::Int;
                             triple_break_cathode_inlet::Bool=false)
    labels = Any[string(i) for i in 1:nb_gc]
    if nb_gc == 1
        labels[1] = rich("1 layer")
    elseif cfg.type_flow == :counter_flow
        labels[1] = rich("cathode inlet\n(counter-flow)")
        labels[end] = rich("anode inlet\n(counter-flow)")
    else
        labels[1] = rich("inlet\n(co-flow)")
        labels[end] = rich("outlet\n(co-flow)")
    end
    return labels
end

"""Clear dynamic axes and attached legends for live plot refresh."""
function _clear_dynamic_axes!(items...)
    figures = Figure[]

    function _visit(item)
        item === nothing && return nothing
        if item isa Axis
            empty!(item)
            fig = item.parent
            fig in figures || push!(figures, fig)
        elseif item isa AbstractArray
            foreach(_visit, item)
        end
        return nothing
    end

    foreach(_visit, items)

    for fig in figures
        for block in reverse(copy(fig.content))
            block isa Legend && delete!(block)
        end
    end
    return nothing
end

"""Remove dynamic legends without clearing already plotted data."""
function _clear_dynamic_legends!(items...)
    figures = Figure[]

    function _visit(item)
        item === nothing && return nothing
        if item isa Axis
            fig = item.parent
            fig in figures || push!(figures, fig)
        elseif item isa AbstractArray
            foreach(_visit, item)
        end
        return nothing
    end

    foreach(_visit, items)

    for fig in figures
        for block in reverse(copy(fig.content))
            block isa Legend && delete!(block)
        end
    end
    return nothing
end

"""Return a coherent PDF export path with incremented index when needed."""
function _resolve_pdf_export_path(folder_path::String,
                                  filename::String)::String
    stem, _suffix = splitext(filename)
    candidate_stem = stem

    if isfile(joinpath(folder_path, candidate_stem * ".pdf"))
        stem_with_index = match(r"^(.*)_(\d+)$", stem)
        prefix = stem_with_index === nothing ? stem : stem_with_index.captures[1]
        counter = stem_with_index === nothing ? 1 : parse(Int, stem_with_index.captures[2]) + 1

        while true
            candidate_stem = "$(prefix)_$(counter)"
            !isfile(joinpath(folder_path, candidate_stem * ".pdf")) && break
            counter += 1
        end
    end

    return joinpath(folder_path, candidate_stem * ".pdf")
end

"""Save a figure to the project results directory."""
function saving_instructions!(_simu,
                              root_folder::String,
                              subfolder_name::String,
                              filename::String,
                              fig)
    fig === nothing && return nothing

    # Resolve current file and define repository markers to locate project root.
    cur = abspath(@__FILE__)
    markers = [".git", "Project.toml", "Manifest.toml", "README.md"]
    project_root = nothing
    parent = dirname(cur)
    while true
        any(ispath(joinpath(parent, m)) for m in markers) && (project_root = parent; break)
        new_parent = dirname(parent)
        new_parent == parent && break
        parent = new_parent
    end
    # Fallback to current working directory if no marker found.
    project_root === nothing && (project_root = pwd())

    # Build destination folder under project root and create it.
    folder_path = joinpath(project_root, root_folder, subfolder_name)
    mkpath(folder_path)

    pdf_path = _resolve_pdf_export_path(folder_path, filename)
    save(pdf_path, fig; backend=CairoMakie)
    return nothing
end

end # module PlotHelpers

