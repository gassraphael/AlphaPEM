# plot_eis.jl
#
# Frequency-domain EIS analysis and plots for cell voltage and internal states.
#
# The Fourier analysis is shared by voltage and states. Their physical conventions
# remain distinct and are implemented by `_eis_voltage_point` and `_eis_state_point`.


# ═══════════════════════════════════════════════════════════════════════════════
#  Shared EIS plot primitives
# ═══════════════════════════════════════════════════════════════════════════════

"""Draw one EIS series (Nyquist or Bode gain) with the shared marker-and-line style."""
function _plot_eis_series!(ax, x, y; color, linewidth=1.8, markersize=7,
                           strokecolor=:white, strokewidth=0.4,
                           linestyle=:solid, label=nothing)
    scatterlines!(ax, x, y; color, linewidth, markersize, strokecolor, strokewidth, linestyle, label)
    return nothing
end

"""Draw one EIS Bode-phase series, breaking connecting lines at ±180° wraps."""
function _plot_eis_phase_series!(ax, f, phase; color, linewidth=1.8, markersize=7,
                                 strokecolor=:white, strokewidth=0.4,
                                 linestyle=:solid, label=nothing)
    scatterlines!(ax, _break_phase_wraps(f, phase)...;
                  color, linewidth, markersize, strokecolor, strokewidth, linestyle, label)
    return nothing
end

"""Apply common Nyquist axis geometry and layout."""
function _finalize_eis_nyquist_axis!(ax; xlabel, ylabel, title, legend=false)
    ax.aspect = DataAspect()
    _finalize_axis!(ax; xlabel, ylabel, title, legend)
    return nothing
end

"""Apply common logarithmic-frequency Bode-gain axis layout."""
function _finalize_eis_gain_axis!(ax; xlabel, ylabel, title, legend=false, legend_position=:rt)
    ax.xscale = log10
    _finalize_axis!(ax; xlabel, ylabel, title, legend, legend_position)
    return nothing
end

"""Apply common logarithmic-frequency Bode-phase axis layout.

With `fixed_range=true` the phase axis is fixed to [-200°, 200°] with ticks every 90° (internal-state diagrams);
otherwise the limits and ticks are left automatic (impedance diagram).
"""
function _finalize_eis_phase_axis!(ax; xlabel, ylabel, title, legend=false, fixed_range=false)
    ax.xscale = log10
    fixed_range && (ax.yticks = -180:90:180)
    _finalize_axis!(ax; xlabel, ylabel, title, legend)
    fixed_range && ylims!(ax, -200, 200)
    return nothing
end

"""Shared frequency label used by the voltage and internal-state Bode diagrams."""
_eis_frequency_label() = rich("Frequency ", lsub("f", ""), " (Hz)")


# ═══════════════════════════════════════════════════════════════════════════════
#  Cell-voltage EIS curves (legacy public plotting API)
# ═══════════════════════════════════════════════════════════════════════════════

"""Plot the cell-voltage EIS Nyquist diagram using the legacy impedance conventions."""
function plot_EIS_curve_Nyquist(cd::AbstractCurrent,
                                Fourier_results::AbstractVector{FourierOutputs},
                                ax)
    Z_reals = Float64[]
    minus_Z_imags = Float64[]

    for res in Fourier_results
        Z_real, minus_Z_imag, _, _, _ = _eis_voltage_point(cd, res)
        if isfinite(Z_real) && isfinite(minus_Z_imag)
            push!(Z_reals, Z_real)
            push!(minus_Z_imags, minus_Z_imag)
        end
    end

    isempty(Z_reals) && return nothing

    _plot_eis_series!(ax, Z_reals, minus_Z_imags;
                              color=_publication_colors()[1], markersize=7)
    _finalize_eis_nyquist_axis!(ax;
                                xlabel=rich("Real impedance ", lsub("Z", "real"), " (mΩ·cm²)"),
                                ylabel=rich("Imaginary impedance ", rich("-", lsub("Z", "imag")), " (mΩ·cm²)"),
                                title="EIS - Nyquist")
    return nothing
end

"""Plot the cell-voltage EIS Bode amplitude using the legacy impedance conventions."""
function plot_EIS_curve_Bode_amplitude(cd::AbstractCurrent,
                                       Fourier_results::AbstractVector{FourierOutputs},
                                       ax)
    fs = Float64[]
    abs_Zs = Float64[]

    for res in Fourier_results
        _, _, f, abs_Z, _ = _eis_voltage_point(cd, res)
        if isfinite(f) && isfinite(abs_Z)
            push!(fs, f)
            push!(abs_Zs, abs_Z)
        end
    end

    isempty(fs) && return nothing

    _plot_eis_series!(ax, fs, abs_Zs;
                           color=_publication_colors()[2], markersize=7)
    _finalize_eis_gain_axis!(ax;
                             xlabel=_eis_frequency_label(),
                             ylabel=rich("Impedance amplitude ", lsub("|Z|", ""), " (mΩ·cm²)"),
                             title="EIS - Bode amplitude")
    return nothing
end

"""Plot the cell-voltage EIS Bode phase using the legacy impedance conventions."""
function plot_EIS_curve_Bode_angle(cd::AbstractCurrent,
                                   Fourier_results::AbstractVector{FourierOutputs},
                                   ax)
    fs = Float64[]
    phi_degs = Float64[]

    for res in Fourier_results
        _, _, f, _, phi_deg = _eis_voltage_point(cd, res)
        if isfinite(f) && isfinite(phi_deg)
            push!(fs, f)
            push!(phi_degs, phi_deg)
        end
    end

    isempty(fs) && return nothing

    _plot_eis_phase_series!(ax, fs, phi_degs;
                            color=_publication_colors()[3], markersize=7)
    _finalize_eis_phase_axis!(ax;
                              xlabel=_eis_frequency_label(),
                              ylabel=rich("Phase ", lsub("ϕ", ""), " (°)"),
                              title="EIS - Bode phase")
    return nothing
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Registry of the internal states displayed in the EIS frequency response
# ═══════════════════════════════════════════════════════════════════════════════

"""Internal states displayed as an EIS frequency response, in display order.

Each entry provides:
- `key`: dispatch key used by `_eis_state_signals`.
- `title`: human-readable quantity name.
- `unit`: unit of the amplitude.
- `symbol`: `(base, subscript)` pair used to typeset the quantity.

The cell voltage is intentionally absent: its amplitude and phase are already given by
the impedance Bode diagrams.
"""
const EIS_STATE_RESPONSES = (
    (key=:T,      title="Temperature",             unit="K",      symbol=("T", "")),
    (key=:Cv,     title="Vapor concentration",     unit="mol·m⁻³", symbol=("C", "v")),
    (key=:s,      title="Liquid saturation",       unit="–",      symbol=("s", "")),
    (key=:lambda, title="Water content",           unit="–",      symbol=("λ", "")),
    (key=:CH2,    title="Hydrogen concentration",  unit="mol·m⁻³", symbol=("C", "H₂")),
    (key=:CO2,    title="Oxygen concentration",    unit="mol·m⁻³", symbol=("C", "O₂")),
    (key=:P,      title="Pressure",                unit="bar",    symbol=("P", "")),
    (key=:CN2,    title="Nitrogen concentration",  unit="mol·m⁻³", symbol=("C", "N₂")),
    (key=:Phi_a,  title="Anode relative humidity", unit="–",      symbol=("Φ", "a")),
    (key=:Phi_c,  title="Cathode relative humidity", unit="–",    symbol=("Φ", "c")),
)

"""Return the registry entry associated with `key` (a `Symbol` or its string representation)."""
function eis_state_spec(key::Symbol)
    idx = findfirst(s -> s.key === key, EIS_STATE_RESPONSES)
    idx === nothing && throw(ArgumentError("Unknown EIS internal state `$key`."))
    return EIS_STATE_RESPONSES[idx]
end

function eis_state_spec(key::AbstractString)
    idx = findfirst(s -> string(s.key) == key, EIS_STATE_RESPONSES)
    idx === nothing && throw(ArgumentError("Unknown EIS internal state `$key`."))
    return EIS_STATE_RESPONSES[idx]
end


# ═══════════════════════════════════════════════════════════════════════════════
#  Signals: time series analysed for each internal state
# ═══════════════════════════════════════════════════════════════════════════════

"""One time series of an internal state, with its legend and line style."""
struct EISStateSignal
    base      :: String
    sub       :: String
    values    :: Vector{Float64}
    color     :: RGBf
    linestyle :: Symbol
end

EISStateSignal(base::AbstractString, sub::AbstractString, values::AbstractVector{<:Real},
               color; linestyle::Symbol=:solid) =
    EISStateSignal(String(base), String(sub), Float64.(values), RGBf(color), linestyle)

"""Return the typeset legend label of a signal."""
signal_label(s::EISStateSignal) = lsub(s.base, s.sub)

"""Return the plain-text name of a signal (used for data exports)."""
signal_name(s::EISStateSignal) = isempty(s.sub) ? s.base : string(s.base, "_", s.sub)

"""Return the MEA node of `mea` associated with `layer` (middle GDL/MPL node for porous layers)."""
function _layer_node(mea, layer::Symbol, gdl_mid::Int, mpl_mid::Int)
    layer === :agc  && return mea.agc
    layer === :agdl && return mea.agdl[gdl_mid]
    layer === :ampl && return mea.ampl[mpl_mid]
    layer === :acl  && return mea.acl
    layer === :mem  && return mea.mem
    layer === :ccl  && return mea.ccl
    layer === :cmpl && return mea.cmpl[mpl_mid]
    layer === :cgdl && return mea.cgdl[gdl_mid]
    layer === :cgc  && return mea.cgc
    throw(ArgumentError("Unknown MEA layer `$layer`."))
end

"""Build one signal per layer at the middle GC node, using the publication palette in layer order.

The subscript of each signal is `<prefix>,<layer>` (or `<layer>` when `prefix` is empty).
"""
function _layer_signals(outputs::SimulationOutputs,
                        cfg::SimulationConfig,
                        base::String,
                        prefix::String,
                        layers,
                        getter::Function)
    gdl_mid = middle_gdl_index(cfg)
    mpl_mid = middle_mpl_index(cfg)
    return [
        EISStateSignal(base, isempty(prefix) ? string(layer) : string(prefix, ",", layer),
                       extract_mid_mea_series(outputs, cfg, mea -> getter(_layer_node(mea, layer, gdl_mid, mpl_mid))),
                       _series_color(i))
        for (i, layer) in enumerate(layers)
    ]
end

"""Gas-channel nodes analysed along the channel: inlet, middle and outlet, as `(index, label)` pairs.

Duplicates are removed when the channel has fewer than three nodes.
"""
function _gc_probes(cfg::SimulationConfig)
    nb_gc = cfg.numerical_parameters.nb_gc
    mid_gc = middle_gas_channel_index(cfg)
    probes = [1 => "in"]
    mid_gc != 1 && mid_gc != nb_gc && push!(probes, mid_gc => "mid")
    nb_gc != 1 && push!(probes, nb_gc => "out")
    return probes
end

"""Build the inlet / middle / outlet signals of one quantity along the gas channel.

The middle position is only added when it differs from the inlet and the outlet, and the
outlet only when the channel has several nodes.
"""
function _gc_position_signals(outputs::SimulationOutputs,
                              cfg::SimulationConfig,
                              base::String,
                              sub::String,
                              color,
                              getter::Function)
    styles = Dict("in" => :dash, "mid" => :solid, "out" => :dot)
    return [EISStateSignal(base, string(sub, ",", label), extract_mea_series(outputs, i, getter), color;
                           linestyle=styles[label])
            for (i, label) in _gc_probes(cfg)]
end

"""Build the anode (solid) and cathode (dashed) signals of one per-GC-node quantity.

`anode` and `cathode` hold one series per GC node; only the inlet, middle and outlet nodes of
`probes` are kept. The subscript is `<prefix>,<label>`, with the anode and cathode prefixes given separately.
"""
function _anode_cathode_node_signals(base::String,
                                     anode_prefix::String,
                                     cathode_prefix::String,
                                     anode,
                                     cathode,
                                     probes)
    signals = EISStateSignal[]
    for (k, (i, label)) in enumerate(probes)
        push!(signals, EISStateSignal(base, "$(anode_prefix),$(label)", anode[i], _series_color(k); linestyle=:solid))
        push!(signals, EISStateSignal(base, "$(cathode_prefix),$(label)", cathode[i], _series_color(k); linestyle=:dash))
    end
    return signals
end

const _GAS_LAYERS_ANODE = (:agc, :agdl, :ampl, :acl)
const _GAS_LAYERS_CATHODE = (:ccl, :cmpl, :cgdl, :cgc)
const _GAS_LAYERS_ALL = (:agc, :agdl, :ampl, :acl, :ccl, :cmpl, :cgdl, :cgc)
const _THERMAL_LAYERS = (:agdl, :ampl, :acl, :mem, :ccl, :cmpl, :cgdl)

"""Return the time series analysed for the internal state `key`.

One method per entry of `EIS_STATE_RESPONSES`; the series follow the temporal displays
(middle GC node unless stated otherwise).
"""
_eis_state_signals(key::Symbol, args...) = _eis_state_signals(Val(key), args...)

_eis_state_signals(::Val{:T}, outputs, fc, cd, cfg) =
    _layer_signals(outputs, cfg, "T", "", _THERMAL_LAYERS, node -> node.T)

_eis_state_signals(::Val{:Cv}, outputs, fc, cd, cfg) =
    _layer_signals(outputs, cfg, "C", "v", _GAS_LAYERS_ALL, node -> node.C_v)

_eis_state_signals(::Val{:s}, outputs, fc, cd, cfg) =
    _layer_signals(outputs, cfg, "s", "", _GAS_LAYERS_ALL, node -> node.s)

function _eis_state_signals(::Val{:lambda}, outputs, fc, cd, cfg)
    signals = [
        EISStateSignal("λ", "acl", extract_mid_mea_series(outputs, cfg, mea -> mea.acl.lambda), _series_color(3)),
        EISStateSignal("λ", "mem", extract_mid_mea_series(outputs, cfg, mea -> mea.mem.lambda), _series_color(4)),
    ]
    append!(signals, _gc_position_signals(outputs, cfg, "λ", "ccl", _series_color(5), mea -> mea.ccl.lambda))
    return signals
end

_eis_state_signals(::Val{:CH2}, outputs, fc, cd, cfg) =
    _layer_signals(outputs, cfg, "C", "H₂", _GAS_LAYERS_ANODE, node -> node.C_H2)

function _eis_state_signals(::Val{:CO2}, outputs, fc, cd, cfg)
    signals = [
        EISStateSignal("C", "O₂,Pt", extract_mid_derived_gc_series(outputs, cfg, x -> x.C_O2_Pt), _series_color(10)),
        EISStateSignal("C", "O₂,ccl", extract_mid_mea_series(outputs, cfg, mea -> mea.ccl.C_O2), _series_color(5)),
    ]
    gdl_mid = middle_gdl_index(cfg)
    mpl_mid = middle_mpl_index(cfg)
    push!(signals,
          EISStateSignal("C", "O₂,cmpl", extract_mid_mea_series(outputs, cfg, mea -> mea.cmpl[mpl_mid].C_O2), _series_color(6)),
          EISStateSignal("C", "O₂,cgdl", extract_mid_mea_series(outputs, cfg, mea -> mea.cgdl[gdl_mid].C_O2), _series_color(7)))
    append!(signals, _gc_position_signals(outputs, cfg, "C", "O₂,cgc", _series_color(8), mea -> mea.cgc.C_O2))
    return signals
end

function _eis_state_signals(::Val{:P}, outputs, fc, cd, cfg)
    # Partial-pressure sum of the ideal-gas mixture at the middle GC node (bar).
    pressure(C_sum, T) = C_sum .* R .* T ./ 1e5
    agc(f) = extract_mid_mea_series(outputs, cfg, mea -> f(mea.agc))
    cgc(f) = extract_mid_mea_series(outputs, cfg, mea -> f(mea.cgc))

    return [
        EISStateSignal("P", "agc", pressure(agc(n -> n.C_v) .+ agc(n -> n.C_H2) .+ agc(n -> n.C_N2), agc(n -> n.T)), _series_color(1)),
        EISStateSignal("P", "cgc", pressure(cgc(n -> n.C_v) .+ cgc(n -> n.C_O2) .+ cgc(n -> n.C_N2), cgc(n -> n.T)), _series_color(2)),
        EISStateSignal("P", "a,in", extract_derived_series(outputs, x -> x.Pa_in) ./ 1e5, _series_color(3)),
        EISStateSignal("P", "c,in", extract_derived_series(outputs, x -> x.Pc_in) ./ 1e5, _series_color(5)),
    ]
end

function _eis_state_signals(::Val{:CN2}, outputs, fc, cd, cfg)
    nb_gc = cfg.numerical_parameters.nb_gc
    return _anode_cathode_node_signals(
        "C", "N₂,agc", "N₂,cgc",
        [extract_mea_series(outputs, i, mea -> mea.agc.C_N2) for i in 1:nb_gc],
        [extract_mea_series(outputs, i, mea -> mea.cgc.C_N2) for i in 1:nb_gc],
        _gc_probes(cfg))
end

_eis_state_signals(::Val{:Phi_a}, outputs, fc, cd, cfg) =
    _layer_signals(outputs, cfg, "Φ", "a", _GAS_LAYERS_ANODE, node -> node.C_v / C_v_sat(node.T))

_eis_state_signals(::Val{:Phi_c}, outputs, fc, cd, cfg) =
    _layer_signals(outputs, cfg, "Φ", "c", _GAS_LAYERS_CATHODE, node -> node.C_v / C_v_sat(node.T))


# ═══════════════════════════════════════════════════════════════════════════════
#  Frequency response computation and plotting
# ═══════════════════════════════════════════════════════════════════════════════

"""Safety factor κ: an oscillation is considered measurable only if its amplitude is at least `κ · atol`,
where `atol` is the absolute tolerance of the ODE solver (below, the response is buried in the solver error)."""
const EIS_MEASURABILITY_FACTOR = 10.0

"""Compute the signals of `key` and their Fourier response at every usable EIS frequency.

Points whose oscillation is below `EIS_MEASURABILITY_FACTOR · atol` keep their frequency but get a `NaN` gain and phase.

The returned gains (NamedTuples `f`, `gain`, `phase_deg`) are normalised by the imposed current amplitude (unit of the state per A·cm⁻²)
and the phases are the shift relative to the imposed current (positive when the state leads it),
without the π offset of the impedance Bode diagram.

Returns `(signals, responses)`; both are empty when no frequency segment is fully simulated.
"""
function _eis_state_responses(key::Symbol,
                              outputs::SimulationOutputs,
                              fc::AbstractFuelCell,
                              cd::AbstractCurrent,
                              cfg::SimulationConfig)
    t = time_history(outputs)
    signals = _eis_state_signals(key, outputs, fc, cd, cfg)
    A_min = EIS_MEASURABILITY_FACTOR * cfg.numerical_parameters.atol
    responses = map(signals) do s
        points = map(make_Fourier_transformation(t, s.values, cd, cfg)) do res
            p = _eis_state_point(cd, res)
            res.A >= A_min ? p : (p[1], NaN, NaN)
        end
        keep = filter(p -> isfinite(p[1]), points)
        (f = first.(keep), gain = getindex.(keep, 2), phase_deg = getindex.(keep, 3))
    end
    all(r -> isempty(r.f), responses) && return (EISStateSignal[], NamedTuple[])
    return (signals, responses)
end

"""Insert a `NaN` wherever the phase wraps around ±180° between two consecutive points.

This keeps the curve from drawing a near-vertical segment across the plot at each wrap.
"""
function _break_phase_wraps(f::AbstractVector, phase::AbstractVector)
    f_out = Float64[]
    phase_out = Float64[]
    for i in eachindex(f, phase)
        if i > firstindex(f) && abs(phase[i] - phase[i - 1]) > 180
            push!(f_out, NaN)
            push!(phase_out, NaN)
        end
        push!(f_out, f[i])
        push!(phase_out, phase[i])
    end
    return f_out, phase_out
end

"""Rich axis labels of an EIS state response, as a `NamedTuple`.

Fields: `frequency`, `gain`, `phase`, `real` and `imag` (Nyquist axes).
"""
function _eis_state_labels(spec)
    symbol = lsub(spec.symbol...)
    current = lsub("i", "fc")
    gain_unit = string(spec.unit, "/(A·cm⁻²)")
    transfer = rich("Δ", symbol, "/Δ", current)
    return (
        frequency = _eis_frequency_label(),
        gain      = rich(spec.title, " gain ", rich("|", transfer, "|"), " (", gain_unit, ")"),
        phase     = rich("Phase of ", symbol, " vs ", current, " (°)"),
        real      = rich("Real part of ", transfer, " (", gain_unit, ")"),
        imag      = rich("Imaginary part ", rich("-Im(", transfer, ")"), " (", gain_unit, ")"),
    )
end

"""
    plot_internal_state_EIS_response(key, outputs, fc, cd, cfg, ax_nyquist, ax_gain, ax_phase)

Plot the frequency response of one internal state during an EIS simulation.

The response of every series of the state is the complex gain `Δx/Δi_fc` (state oscillation over
imposed current oscillation), drawn at each excitation frequency as:
- a Nyquist diagram on `ax_nyquist` (minus the imaginary part versus the real part, as for the impedance);
- the gain modulus (amplitude divided by the imposed current amplitude) on `ax_gain`;
- the phase shift with respect to the imposed current (positive when the state leads it,
  no π offset) on `ax_phase`.

# Arguments
- `key::Symbol`: Internal state, one of the `key` entries of `EIS_STATE_RESPONSES`.
- `outputs::SimulationOutputs`: Simulation outputs.
- `fc::AbstractFuelCell`: Fuel-cell definition.
- `cd::AbstractCurrent`: EIS current profile.
- `cfg::SimulationConfig`: Simulation configuration.
- `ax_nyquist`, `ax_gain`, `ax_phase`: Makie axes receiving the three diagrams. The two Bode axes
  must use a logarithmic frequency scale.
"""
function plot_internal_state_EIS_response(key::Symbol,
                                          outputs::SimulationOutputs,
                                          fc::AbstractFuelCell,
                                          cd::AbstractCurrent,
                                          cfg::SimulationConfig,
                                          ax_nyquist,
                                          ax_gain,
                                          ax_phase)
    spec = eis_state_spec(key)
    signals, responses = _eis_state_responses(key, outputs, fc, cd, cfg)
    isempty(signals) && return nothing

    for (sig, resp) in zip(signals, responses)
        style = (color=sig.color, linestyle=sig.linestyle, linewidth=1.8, markersize=6,
                 strokecolor=:white, strokewidth=0.4, label=signal_label(sig))
        ok = findall(i -> isfinite(resp.gain[i]) && isfinite(resp.phase_deg[i]), eachindex(resp.f))
        isempty(ok) && continue
        f, gain, phase = resp.f[ok], resp.gain[ok], resp.phase_deg[ok]
        _plot_eis_series!(ax_nyquist, gain .* cosd.(phase), -gain .* sind.(phase); style...)
        _plot_eis_series!(ax_gain, f, gain; style...)
        _plot_eis_phase_series!(ax_phase, f, phase; style...)
    end

    labels = _eis_state_labels(spec)
    n_hidden = sum(count(g -> !isfinite(g), r.gain) for r in responses)

    # Same geometry convention as the impedance Nyquist diagram.
    _finalize_eis_nyquist_axis!(ax_nyquist;
                                xlabel=labels.real,
                                ylabel=labels.imag,
                                title="EIS - $(spec.title) Nyquist")

    _finalize_eis_gain_axis!(ax_gain;
                             xlabel=labels.frequency,
                             ylabel=labels.gain,
                             title="EIS - $(spec.title) gain",
                             legend=true,
                             legend_position=:rt)

    _finalize_eis_phase_axis!(ax_phase;
                              xlabel=labels.frequency,
                              ylabel=labels.phase,
                              title="EIS - $(spec.title) phase",
                              fixed_range=true)

    if n_hidden > 0
        for ax in (ax_gain, ax_phase)
            text!(ax, 0.02, 0.02; text="$n_hidden point(s) not shown:\nresponse below solver accuracy",
                  space=:relative, align=(:left, :bottom), fontsize=10, color=:gray40)
        end
    end
    return nothing
end

"""
    eis_state_response_table(key, outputs, fc, cd, cfg)

Return the `(headers, data)` table of the frequency response of one internal state.

The first column is the frequency (Hz), followed by the amplitude of every series and then
the phase of every series.
"""
function eis_state_response_table(key::Symbol,
                                  outputs::SimulationOutputs,
                                  fc::AbstractFuelCell,
                                  cd::AbstractCurrent,
                                  cfg::SimulationConfig)
    spec = eis_state_spec(key)
    signals, responses = _eis_state_responses(key, outputs, fc, cd, cfg)
    isempty(signals) && return (["Frequency (Hz)"], zeros(0, 1))

    headers = ["Frequency (Hz)";
               ["$(signal_name(s)) amplitude ($(spec.unit)/(A·cm⁻²))" for s in signals];
               ["$(signal_name(s)) phase (°)" for s in signals]]
    data = hcat(responses[1].f,
                (r.gain for r in responses)...,
                (r.phase_deg for r in responses)...)
    return (headers, data)
end
