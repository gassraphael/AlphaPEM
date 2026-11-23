# -*- coding: utf-8 -*-

"""
Interactive display of the figures in a web browser.

The figures are served by a Bonito server (WGLMakie), so no OpenGL context, X11 display or
GPU is required (Docker, HPC nodes, remote machines). One figure is served directly under
its name; several figures are listed in a menu. Every figure page offers PDF and PNG
downloads (exported with CairoMakie).

Sections:
1. Figure names and groups
2. Server and browser
3. Pages
4. Serving the figures
"""

import HTTP
using ..Core.Models.PlotHelpers: export_figure_bytes

const _PLOT_SERVER = Ref{Any}(nothing)
const _PLOT_SERVER_LOCK = ReentrantLock()
const _served_figure_routes = Dict{String, Any}()      # route => Bonito app or download handler.
const _prepared_figure_specs = Dict{String, String}[]  # name/title/group, same order as the prepared figures.
const _PLOT_LOGO_PATH = joinpath(@__DIR__, "..", "..", "..", "docs", "src", "images", "logo.png")

# Page stylesheet, read once at load time (and tracked so that edits trigger a recompilation).
const _PLOT_PAGE_CSS_PATH = joinpath(@__DIR__, "interactive_display.css")
include_dependency(_PLOT_PAGE_CSS_PATH)
const _PLOT_PAGE_CSS = read(_PLOT_PAGE_CSS_PATH, String)


# ═══════════════════════════════════════════════════════════════════════════════
#  1. Figure names and groups
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _interactive_plot_specs(cfg, nb_gc) -> Vector{Dict{String, String}}

Name the figures of a run, in display order.

Names are the file names (without extension) used when saving, and become the URL of each
figure. They match the figures returned by `figures_preparation` for the same configuration.

Parameters
----------
cfg : SimulationConfig
    Simulation configuration (current profile, display type and timing).
nb_gc : Integer
    Number of gas-channel nodes.

Returns
-------
Vector{Dict{String, String}}
    One dictionary per figure with the keys `"name"`, `"title"` and `"group"`.
"""
function _interactive_plot_specs(cfg::SimulationConfig, nb_gc::Integer)::Vector{Dict{String, String}}
    if cfg.type_display == :synthetic
        spec(name, title) = Dict("name" => name, "title" => title, "group" => "Synthetic view")
        if cfg.type_current isa StepParams
            specs = [spec("step_current_syn", "Step current: internal states")]
            cfg.display_timing == :postrun && nb_gc >= 3 &&
                push!(specs, spec("final_GC_profiles", "Final gas-channel profiles"))
            return specs
        elseif cfg.type_current isa PolarizationParams
            return [spec("pola_curve", "Polarization curve")]
        elseif cfg.type_current isa PolarizationCalibrationParams
            return [spec("pola_curve_cali", "Calibration polarization curve")]
        elseif cfg.type_current isa EISParams
            return [spec("Nyquist_plot_syn", "Impedance spectra")]
        end
    end
    return [Dict("name" => splitext(spec["filename"])[1][4:end], "title" => spec["title"], "group" => spec["group"])
            for spec in _web_plot_specs(cfg, nb_gc)]
end

"""
    _store_prepared_specs!(cfg, nb_gc)

Remember the names of the figures just prepared for `cfg`, for `_open_interactive_figures!`.

Parameters
----------
cfg : SimulationConfig
    Simulation configuration.
nb_gc : Integer
    Number of gas-channel nodes.
"""
function _store_prepared_specs!(cfg::SimulationConfig, nb_gc::Integer)
    empty!(_prepared_figure_specs)
    append!(_prepared_figure_specs, _interactive_plot_specs(cfg, nb_gc))
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  2. Server and browser
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _plot_server() -> Bonito.Server

Return the Bonito server hosting the interactive figures, starting it on first use.

The host and port can be set with the environment variables `ALPHAPEM_PLOT_HOST` (default
`127.0.0.1`; use `0.0.0.0` inside a container) and `ALPHAPEM_PLOT_PORT` (default `9385`).

Returns
-------
Bonito.Server
    The shared server.
"""
function _plot_server()
    lock(_PLOT_SERVER_LOCK) do
        if _PLOT_SERVER[] === nothing
            host = get(ENV, "ALPHAPEM_PLOT_HOST", "127.0.0.1")
            port = parse(Int, get(ENV, "ALPHAPEM_PLOT_PORT", "9385"))
            _PLOT_SERVER[] = WGLMakie.Bonito.Server(host, port)
        end
        return _PLOT_SERVER[]
    end
end

"""
    _plot_url(route) -> String

Build the URL under which `route` is reachable from the host machine.

Parameters
----------
route : AbstractString
    Route served by the plot server, such as `"/"` or `"/ucell_temporal"`.

Returns
-------
String
    Full URL (`0.0.0.0` is reported as `localhost`).
"""
function _plot_url(route::AbstractString)
    server = _plot_server()
    host = String(getfield(server, :url))
    host = host in ("0.0.0.0", "::") ? "localhost" : host
    return string("http://", host, ":", Int(getfield(server, :port)), route)
end

"""
    _open_browser(url)

Open `url` in the default browser. It silently does nothing when no browser is available
(for instance in Docker), and `ALPHAPEM_OPEN_BROWSER=false` disables it.

Parameters
----------
url : AbstractString
    Page to open.
"""
function _open_browser(url::AbstractString)
    lowercase(get(ENV, "ALPHAPEM_OPEN_BROWSER", "true")) == "false" && return nothing
    cmd = Sys.isapple()   ? `open $url` :
          Sys.iswindows() ? `cmd /c start "" $url` :
                            `xdg-open $url`
    try
        run(pipeline(cmd; stdout=devnull, stderr=devnull); wait=false)
    catch err
        @debug "Could not open a browser" exception=err
    end
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  3. Pages
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _plot_page(subtitle, content...; wide=false)

Wrap `content` in the AlphaPEM page layout (logo, subtitle, main area, footer).

Parameters
----------
subtitle : AbstractString
    Text shown under the logo.
content : Any
    Bonito DOM nodes (or figures) placed in the main area.
wide : Bool
    Use a wider container, for the figure pages.

Returns
-------
Bonito DOM node
    The complete page.

References
----------
- Bonito.jl, https://github.com/SimonDanisch/Bonito.jl
"""
function _plot_page(subtitle::AbstractString, content...; wide::Bool=false)
    Bonito = WGLMakie.Bonito
    DOM = Bonito.DOM
    logo = isfile(_PLOT_LOGO_PATH) ?
        [DOM.img(; src=Bonito.Asset(_PLOT_LOGO_PATH), alt="AlphaPEM logo")] : Any[DOM.h1("AlphaPEM")]
    return DOM.div(
        DOM.style(_PLOT_PAGE_CSS),
        DOM.div(
            DOM.header(logo..., DOM.p(subtitle)),
            DOM.main(content...),
            DOM.footer(DOM.p("AlphaPEM — interactive figures"));
            class=wide ? "ap-container ap-wide" : "ap-container");
        class="ap-page")
end

"""
    _plot_menu_page()

Build the menu page listing every served figure, grouped as in the AlphaPEM web interface.

Returns
-------
Bonito DOM node
    The menu page, with one card (title and name) per figure.
"""
function _plot_menu_page()
    DOM = WGLMakie.Bonito.DOM
    sections = Any[]
    groups = unique(spec["group"] for spec in _prepared_figure_specs)
    for group in groups
        cards = [DOM.a(DOM.strong(spec["title"]), DOM.span(spec["name"]); href="/" * spec["name"], class="ap-card")
                 for spec in _prepared_figure_specs if spec["group"] == group]
        push!(sections, DOM.h2(group), DOM.div(cards...; class="ap-grid"))
    end
    return _plot_page("Menu", sections...)
end


# ═══════════════════════════════════════════════════════════════════════════════
#  4. Serving the figures
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _open_interactive_figures!()

Serve every prepared figure (once) and open the browser on the entry page.

A single figure is served directly under its name; several figures are listed in a menu,
and only the menu URL is printed and opened. It does nothing when the figures are already
served, so it can be called after each refresh of a `:live` display.
"""
function _open_interactive_figures!()
    # Preliminary checks
    isempty(_served_figure_routes) || return nothing
    figs = _collect_web_figures(_prepared_interactive_figures...)
    isempty(figs) && return nothing
    specs = length(_prepared_figure_specs) == length(figs) ? _prepared_figure_specs :
            [Dict("name" => "figure_$i", "title" => "Figure $i", "group" => "Figures") for i in eachindex(figs)]
    Bonito = WGLMakie.Bonito
    server = _plot_server()
    multiple = length(figs) > 1

    # One page and two download endpoints per figure
    for (fig, spec) in zip(figs, specs)
        DataInspector(fig)
        name = spec["name"]
        route = "/" * name
        app = Bonito.App(title="AlphaPEM") do
            DOM = WGLMakie.Bonito.DOM
            back = multiple ? [DOM.a("← Menu"; href="/", class="ap-back")] : Any[]
            actions = DOM.div(
                DOM.a("Download PDF"; href=route * ".pdf", download=name * ".pdf", class="ap-btn"),
                DOM.a("Download PNG"; href=route * ".png", download=name * ".png", class="ap-btn");
                class="ap-actions")
            _plot_page(spec["title"], back..., DOM.div(fig; class="ap-figure"), actions; wide=true)
        end
        Bonito.route!(server, route => app)
        _served_figure_routes[route] = app

        # The figure is exported on request with CairoMakie.
        for (ext, mime) in (("pdf", "application/pdf"), ("png", "image/png"))
            handler = _ -> HTTP.Response(200,
                ["Content-Type" => mime, "Content-Disposition" => "attachment; filename=\"$(name).$(ext)\""],
                export_figure_bytes(fig, ext))
            Bonito.route!(server, "$(route).$(ext)" => handler)
            _served_figure_routes["$(route).$(ext)"] = handler
        end
    end

    # Entry page: the figure itself, or the menu
    entry = multiple ? "/" : "/" * first(specs)["name"]
    if multiple
        menu = Bonito.App(_plot_menu_page; title="AlphaPEM")
        Bonito.route!(server, "/" => menu)
        _served_figure_routes["/"] = menu
    end
    @info "Interactive figures available at $(_plot_url(entry))"
    _open_browser(_plot_url(entry))
    return nothing
end

"""
    _close_interactive_figures!()

Remove every served route so that the server can be reused by a later run.
"""
function _close_interactive_figures!()
    if _PLOT_SERVER[] !== nothing
        for route in keys(_served_figure_routes)
            try
                WGLMakie.Bonito.HTTPServer.delete_route!(_PLOT_SERVER[], route)
            catch err
                @debug "Could not delete Bonito route $route" exception=err
            end
        end
    end
    empty!(_served_figure_routes)
    return nothing
end

"""
    _plot_block_requested() -> Bool

Tell whether the run should keep serving the figures until Ctrl+C.

This is the default in a script (the process would otherwise exit and drop the server)
but not in a REPL; `ALPHAPEM_PLOT_BLOCK=true|false` overrides it.

Returns
-------
Bool
    True when the figures must stay served.
"""
_plot_block_requested() = lowercase(get(ENV, "ALPHAPEM_PLOT_BLOCK", string(!isinteractive()))) == "true"

"""
    _serve_until_interrupted!(message)

Keep serving the figures until the user presses Ctrl+C, then release the routes.

Parameters
----------
message : AbstractString
    Information shown to the user while waiting.
"""
function _serve_until_interrupted!(message::AbstractString)
    @info message
    try
        while true
            sleep(1)
        end
    catch err
        err isa InterruptException || rethrow()
    finally
        _close_interactive_figures!()
    end
    return nothing
end
