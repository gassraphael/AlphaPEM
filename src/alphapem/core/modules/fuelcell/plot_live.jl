# -*- coding: utf-8 -*-

# plot_live.jl
#
# In-place updates of the plots of a `:live` display. Included by `plot_helpers.jl`;
# nothing in here is used by the other display modes.
#
# Rebuilding every plot at each refresh is slow, and in a browser (WGLMakie) the axes are
# visibly empty between the removal of the old plots and the insertion of the new ones.
# Instead, the first refresh creates the plots and remembers them in call order; the next
# refreshes run the same plotting code, but `_lines!` only updates the data of the
# remembered plots, and the static decorations (labels, ticks, legend) are left alone.
#
# Usage, around the plotting calls of one refresh:
#     _live_begin!(axes...)    # starts a pass (the first one records, the others update)
#     ... plot functions calling `_lines!` instead of `lines!` ...
#     _live_end!(axes...)      # ends the pass

export _live_begin!,
       _live_end!,
       _live_updating,
       _lines!

"""
    _LiveAxisState

Plots remembered for one axis between the refreshes of a live display.

Fields
------
plots : Vector{Any}
    Plots created by the first (recording) pass, in call order.
cursor : Int
    Number of plotting calls made in the current pass.
active : Bool
    True between `_live_begin!` and `_live_end!`.
updating : Bool
    True when the current pass updates the remembered plots instead of creating them.
complete : Bool
    True once a recording pass went through to the end.
"""
mutable struct _LiveAxisState
    plots::Vector{Any}
    cursor::Int
    active::Bool
    updating::Bool
    complete::Bool
end

const _live_axes = WeakKeyDict{Any, _LiveAxisState}()

"""
    _live_updating(ax) -> Bool

Tell whether a live pass is currently updating the plots of `ax` in place.

Parameters
----------
ax : Axis
    Axis to test.

Returns
-------
Bool
    True while the plots of `ax` are updated rather than created, in which case the
    static decorations (labels, ticks, legend, annotations) must not be created again.
"""
function _live_updating(ax)
    st = get(_live_axes, ax, nothing)
    return st !== nothing && st.active && st.updating
end

"""
    _each_axis(f, items)

Call `f` on every axis found in `items`, ignoring `nothing`.

Parameters
----------
f : Function
    Function applied to each `Axis`.
items : Any
    Axes, or arrays (of any depth) of axes, as handed to `display!`.
"""
function _each_axis(f, items)
    for item in items
        item === nothing && continue
        if item isa Axis
            f(item)
        elseif item isa AbstractArray
            _each_axis(f, item)
        end
    end
    return nothing
end

"""
    _live_reset!(ax, st)

Remove the plots and the legend of `ax` so that its live state can be recorded again.

Parameters
----------
ax : Axis
    Axis to empty.
st : _LiveAxisState
    State of `ax`, whose remembered plots are forgotten.
"""
function _live_reset!(ax, st::_LiveAxisState)
    empty!(ax)
    if haskey(_axis_legends, ax)
        try
            delete!(_axis_legends[ax])
        catch
        end
        delete!(_axis_legends, ax)
    end
    empty!(st.plots)
    return nothing
end

"""
    _live_begin!(items...)

Start one live pass on the axes found in `items`.

The first pass records the plots that are created; the following ones update them in
place. A first pass that was interrupted (for instance by a closed browser tab) is
discarded and recorded again.

Parameters
----------
items : Any
    Axes, or arrays of axes, that the pass will draw on (`nothing` is ignored).
"""
function _live_begin!(items...)
    _each_axis(items) do ax
        st = get!(() -> _LiveAxisState(Any[], 0, false, false, false), _live_axes, ax)
        !st.complete && !isempty(st.plots) && _live_reset!(ax, st)
        st.cursor = 0
        st.active = true
        st.updating = st.complete
    end
    return nothing
end

"""
    _live_end!(items...)

End the live pass started by `_live_begin!`.

Plots that were recorded but not requested during an updating pass are removed, and the
recorded state is marked complete.

Parameters
----------
items : Any
    The same axes as in the matching `_live_begin!`.
"""
function _live_end!(items...)
    _each_axis(items) do ax
        st = get(_live_axes, ax, nothing)
        (st === nothing || !st.active) && return nothing
        if st.updating && st.cursor < length(st.plots)
            for plt in st.plots[st.cursor + 1:end]
                try
                    delete!(ax, plt)
                catch
                end
            end
            resize!(st.plots, st.cursor)
        end
        st.complete = true
        st.active = false
        st.updating = false
    end
    return nothing
end

# The simulation history is extended in place after each segment. A plot that kept a
# reference to one of its vectors would see it grow before its other vector is updated,
# and the browser renderer would then fail on arrays of different lengths.
_snapshot(x::AbstractVector) = copy(x)
_snapshot(x) = x

"""
    _live_plot!(create, PT, ax, args...; kw...)

Create a plot with `create`, or update the one remembered at the same position during a
live pass. Outside a live pass, it simply calls `create`.

Parameters
----------
create : Function
    Makie plotting function (`lines!`) used when a plot must be created.
PT : Type
    Plot type the remembered plot must have to be updated (`Makie.Lines`); a plot of
    another type is replaced.
ax : Axis
    Axis the plot belongs to.
args : Any
    Data arguments of the plot; they are copied (see `_snapshot`).
kw : Any
    Attributes, used only when a plot is created.

Returns
-------
Plot
    The created or updated plot.

References
----------
- Makie.jl documentation, https://docs.makie.org/stable/
"""
function _live_plot!(create::Function, ::Type{PT}, ax, args...; kw...) where {PT}
    st = get(_live_axes, ax, nothing)
    (st === nothing || !st.active) && return create(ax, args...; kw...)

    args = map(_snapshot, args)
    st.cursor += 1
    i = st.cursor
    if st.updating && i <= length(st.plots)
        if st.plots[i] isa PT
            Makie.update!(st.plots[i], args...)
            return st.plots[i]
        end
        # The plotting code changed since the recording: replace the stale plot.
        try
            delete!(ax, st.plots[i])
        catch
        end
        st.plots[i] = create(ax, args...; kw...)
        return st.plots[i]
    end

    plt = create(ax, args...; kw...)
    push!(st.plots, plt)
    return plt
end

"""`lines!` that updates the remembered plot in place during a live pass (see `_live_plot!`)."""
_lines!(ax, args...; kw...) = _live_plot!(lines!, Makie.Lines, ax, args...; kw...)
