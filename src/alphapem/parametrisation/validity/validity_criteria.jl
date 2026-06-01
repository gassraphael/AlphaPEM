# -*- coding: utf-8 -*-

"""
    ValidityCriteria

Module for validating AlphaPEM polarization curves.

A simulated polarization curve is classified as `:valid` if and only if all enabled
criteria are satisfied simultaneously:

- **start_in_range** — The first voltage value lies within a physically sensible range,
  typically (0 V, E₀) for a PEMFC (~1.23 V at standard conditions).
  Values outside this window indicate a non-physical equilibrium or a solver failure.
- **approx_monotonic** — Voltage decreases (approximately) with current density.
  Small upward bumps below a configurable tolerance are accepted; large rises are not.
- **positive_voltages** — The simulation reached the expected maximum current density
  (within a relative tolerance). When the solver triggers its safety stop before
  `i_max` is reached, it means cell voltage would have gone negative at that operating
  point — the criterion therefore fails. For the vector-based overload this criterion
  still checks the minimum voltage directly (legacy behaviour).

# Exports
- `ValidityCriteriaConfig`      — Configurable thresholds for each criterion
- `ValidationResult`            — Structured classification outcome with per-criterion flags
- `classify_polarization_curve` — Main entry point (`:valid` / `:invalid`)
- `check_start_voltage_range`   — Starting-voltage criterion
- `check_monotonicity`          — Approximate-monotonicity criterion
- `check_positive_voltages`     — No-negative-voltage criterion
"""
module ValidityCriteria

# Standard library only — no external packages needed for this module.

include(joinpath(@__DIR__, "../../utils/physics_constants.jl"))


# ─────────────────────────────────────────────────────────────────────────────
# INTERNAL HELPERS
# ─────────────────────────────────────────────────────────────────────────────

"""Return a consistently ordered `(Ucell, ifc)` pair.

The validity checks operate on the polarization curve as a function of current.
To make the functions robust to arbitrary input ordering, we sort the samples
by current density.
"""
function _sorted_curve(Ucell::AbstractVector{<:Real}, ifc::AbstractVector{<:Real})
    idx = sortperm(collect(ifc))
    return collect(Float64, Ucell[idx]), collect(Float64, ifc[idx])
end

export ValidityCriteriaConfig,
       ValidationResult,
       classify_polarization_curve,
       check_start_voltage_range,
       check_monotonicity,
       check_positive_voltages,
       _extract_polarization_sampling_indices

# ─────────────────────────────────────────────────────────────────────────────
# DATA STRUCTURES
# ─────────────────────────────────────────────────────────────────────────────

"""
    ValidityCriteriaConfig

Thresholds and activation flags controlling how a polarization curve is validated.

# Fields
- `voltage_range::Tuple{Float64,Float64}`: Acceptable window for the starting voltage (V).
  Default `(0.0, E0)` — the theoretical open-circuit voltage of a PEMFC.
- `monotonic_threshold::Float64`: Maximum tolerated upward bump between consecutive
  (current, voltage) points (V). Default `0.005` (5 mV).
- `min_voltage_threshold::Float64`: Minimum accepted voltage anywhere in the curve (V).
  Default `0.0`.
- `apply_start_range::Bool`: Enable the starting-voltage check. Default `true`.
- `apply_monotonicity::Bool`: Enable the monotonicity check. Default `true`.
- `apply_positive_voltages::Bool`: Enable the current-coverage (positive-voltage proxy) check.
  Default `true`.
- `current_reach_tol::Float64`: Relative tolerance for the current-coverage check.
  A simulation passes if `ifc_max_achieved ≥ ifc_max_expected × (1 − tol)`.
  Default `0.01` (1 %).  Only used by the simulator-based overload; the vector-based
  overload still uses `min_voltage_threshold` directly.

# Example
```julia
cfg = ValidityCriteriaConfig(
    voltage_range       = (0.0, E0),
    monotonic_threshold = 0.005,
)
```
"""
Base.@kwdef struct ValidityCriteriaConfig
    voltage_range::Tuple{Float64, Float64} = (0.0, E0)
    monotonic_threshold::Float64           = 0.005
    min_voltage_threshold::Float64         = 0.0
    current_reach_tol::Float64             = 0.01
    apply_start_range::Bool                = true
    apply_monotonicity::Bool               = true
    apply_positive_voltages::Bool          = true
end


"""
    ValidationResult

Structured outcome of a polarization-curve classification.

# Fields
- `classification::Symbol`: Overall verdict — `:valid` or `:invalid`.
- `start_in_range::Union{Bool,Nothing}`: Outcome of the starting-voltage check,
  or `nothing` when the criterion is disabled.
- `is_monotonic::Union{Bool,Nothing}`: Outcome of the monotonicity check,
  or `nothing` when disabled.
- `has_positive_voltages::Union{Bool,Nothing}`: Outcome of the no-negative-voltage check,
  or `nothing` when disabled.
- `details::String`: Human-readable summary (which criterion failed, or "all passed").
"""
struct ValidationResult
    classification::Symbol
    start_in_range::Union{Bool, Nothing}
    is_monotonic::Union{Bool, Nothing}
    has_positive_voltages::Union{Bool, Nothing}
    details::String
end


# ─────────────────────────────────────────────────────────────────────────────
# MAIN CLASSIFICATION
# ─────────────────────────────────────────────────────────────────────────────

"""
    classify_polarization_curve(Ucell, ifc,
                                cfg = ValidityCriteriaConfig()) -> ValidationResult

Classify a polarization curve as `:valid` or `:invalid`.

All enabled criteria are evaluated. The curve is valid only when every enabled
criterion passes. Individual results are stored in the returned `ValidationResult`
for diagnostic purposes.

# Arguments
- `Ucell::Vector{Float64}`: Cell voltages along the curve (V).
- `ifc::Vector{Float64}`: Corresponding current densities (A m⁻²).
- `cfg::ValidityCriteriaConfig`: Thresholds and activation flags (optional).

# Example
```julia
Ucell  = [1.05, 0.92, 0.80, 0.65, 0.50]
ifc    = [0.0,  5e3,  1e4,  1.5e4, 2e4]
result = classify_polarization_curve(Ucell, ifc)
println(result.classification)   # :valid or :invalid
```

"""
function classify_polarization_curve(Ucell::Vector{Float64},
                                     ifc::Vector{Float64},
                                     cfg::ValidityCriteriaConfig = ValidityCriteriaConfig())::ValidationResult
    if isempty(Ucell) || isempty(ifc)
        return ValidationResult(:invalid, false, false, false, "empty curve")
    end
    length(Ucell) == length(ifc) ||
        return ValidationResult(:invalid, false, false, false, "Ucell and ifc must have the same length")

    U_sorted, i_sorted = _sorted_curve(Ucell, ifc)

    start_ok = cfg.apply_start_range ? check_start_voltage_range(U_sorted, cfg.voltage_range) : nothing
    mono_ok = cfg.apply_monotonicity ? check_monotonicity(U_sorted, i_sorted, cfg.monotonic_threshold) : nothing
    pos_ok = cfg.apply_positive_voltages ? check_positive_voltages(U_sorted, cfg.min_voltage_threshold) : nothing

    all_ok = true
    (start_ok === false) && (all_ok = false)
    (mono_ok === false) && (all_ok = false)
    (pos_ok === false) && (all_ok = false)

    classification = all_ok ? :valid : :invalid

    details = if classification == :valid
        "all enabled criteria passed"
    else
        failed = String[]
        (start_ok === false) && push!(failed, "start_in_range")
        (mono_ok === false) && push!(failed, "approx_monotonic")
        (pos_ok === false) && push!(failed, "positive_voltages")
        isempty(failed) ? "invalid" : "failed: " * join(failed, ", ")
    end

    return ValidationResult(classification, start_ok, mono_ok, pos_ok, details)
end


"""
    classify_polarization_curve(simu, cfg = ValidityCriteriaConfig()) -> ValidationResult

Classify a polarization curve directly from an `AlphaPEM` simulator object.

Extracts discretized polarization points (sampled at stabilization times) from `simu.outputs`,
then delegates to the vector-based overload.

**Important**: This function only accepts polarization-type current profiles (standard or calibration).
If the current profile does not expose the required properties (`delta_t_ini`, `delta_t_break`,
and either `delta_i`/`i_max`/`v_load` or `i_exp`/`v_load`), the result will be `:invalid`.

For polarization-type current profiles, this function samples the curve at characteristic times:
- At `delta_t_ini` for the OCV point (i=0)
- At `delta_t_ini + k*(delta_t_load + delta_t_break)` for each current step

This ensures we only validate the stabilized points, not the transient behavior.
"""
function classify_polarization_curve(simu,   # ::AlphaPEM
                                     cfg::ValidityCriteriaConfig = ValidityCriteriaConfig())::ValidationResult
    # Check that simulator has outputs
    if !hasproperty(simu, :outputs) || getproperty(simu, :outputs) === nothing
        return ValidationResult(:invalid, false, false, false, "simulator has no outputs; run a simulation first")
    end

    outputs = getproperty(simu, :outputs)

    # Check that we have the derived outputs structure
    if !hasproperty(outputs, :derived) || !hasproperty(outputs.derived, :Ucell)
        return ValidationResult(:invalid, false, false, false, "simulator outputs do not contain derived.Ucell")
    end

    # Check that we have a current_density profile
    if !hasproperty(simu, :current_density)
        return ValidationResult(:invalid, false, false, false, "simulator does not have a current_density profile")
    end

    cd = getproperty(simu, :current_density)

    # Extract the full time history and voltage trajectory
    if !hasproperty(outputs, :solver) || !hasproperty(outputs.solver, :t)
        return ValidationResult(:invalid, false, false, false, "simulator outputs do not contain solver.t")
    end

    t_hist     = outputs.solver.t
    Ucell_full = outputs.derived.Ucell

    # Extract discretized polarization points (sampled at stabilization times).
    local sample_indices, ifc_sampled, sample_times
    try
        sample_indices, ifc_sampled, sample_times = _extract_polarization_sampling_indices(t_hist, cd)
    catch e
        # Not a polarization profile or unsupported structure.
        return ValidationResult(:invalid, nothing, nothing, nothing,
                                "not a supported polarization profile: $(sprint(showerror, e))")
    end

    # Extract stabilized cell voltages at the sampling times.
    Ucell_sampled = Float64[Ucell_full[i] for i in sample_indices]
    ifc_sampled_f = collect(Float64, ifc_sampled)

    U_sorted, i_sorted = _sorted_curve(Ucell_sampled, ifc_sampled_f)

    # ── Criterion 1: starting voltage ─────────────────────────────────────────
    start_ok = cfg.apply_start_range ?
               check_start_voltage_range(U_sorted, cfg.voltage_range) : nothing

    # ── Criterion 2: approximate monotonicity ─────────────────────────────────
    mono_ok = cfg.apply_monotonicity ?
              check_monotonicity(U_sorted, i_sorted, cfg.monotonic_threshold) : nothing

    # ── Criterion 3: current coverage (proxy for positive voltages) ───────────
    # The simulator stops before U_cell can go negative (safety stop).
    # So, it is verified if the simulation reached the expected maximum current
    # density: if it stopped early, the operating limit (U_cell → 0) was hit.
    pos_ok = if cfg.apply_positive_voltages
        t_end            = t_hist[end]
        # Last expected sample time that was actually reached by the solver.
        # sample_times are analytical; t_end is the exact solver stop time.
        k_last           = something(findlast(st -> st <= t_end, sample_times), 0)
        ifc_max_achieved = k_last > 0 ? ifc_sampled_f[k_last] : 0.0
        ifc_max_expected = ifc_sampled_f[end]
        check_positive_voltages(ifc_max_achieved, ifc_max_expected, cfg.current_reach_tol)
    else
        nothing
    end

    # ── Aggregate ─────────────────────────────────────────────────────────────
    all_ok = !any(x -> x === false, [start_ok, mono_ok, pos_ok])
    classification = all_ok ? :valid : :invalid

    details = if classification == :valid
        "all enabled criteria passed"
    else
        failed = String[]
        (start_ok === false) && push!(failed, "start_in_range")
        (mono_ok  === false) && push!(failed, "approx_monotonic")
        (pos_ok   === false) && push!(failed, "positive_voltages")
        isempty(failed) ? "invalid" : "failed: " * join(failed, ", ")
    end

    return ValidationResult(classification, start_ok, mono_ok, pos_ok, details)
end


"""
    _extract_polarization_sampling_indices(t_hist, cd)
        -> (indices::Vector{Int}, ifc_values::Vector{Float64}, sample_times::Vector{Float64})

Extract the indices in `t_hist` corresponding to characteristic polarization sampling
times, together with the analytically computed current density at each sample and the
sample times themselves.

Returns a tuple `(indices, ifc_values, sample_times)` so that callers never need to
call `current()` from the `Currents` module, keeping this module free of external
dependencies.  The `sample_times` vector is additionally used by
`classify_polarization_curve(simu, cfg)` to determine how far into the profile the
simulation actually ran (current-coverage check).
"""
function _extract_polarization_sampling_indices(t_hist::AbstractVector{<:Real},
                                                 cd)
    # Check required properties
    if !hasproperty(cd, :delta_t_ini)
        throw(ArgumentError("Current profile does not expose `delta_t_ini`"))
    end
    if !hasproperty(cd, :delta_t_break)
        throw(ArgumentError("Current profile does not expose `delta_t_break`"))
    end

    delta_t_ini   = Float64(getproperty(cd, :delta_t_ini))
    delta_t_break = Float64(getproperty(cd, :delta_t_break))

    # ── Standard polarization profile ─────────────────────────────────────────
    if hasproperty(cd, :delta_i) && hasproperty(cd, :i_max) && hasproperty(cd, :v_load)
        delta_i = Float64(getproperty(cd, :delta_i))
        v_load  = Float64(getproperty(cd, :v_load))
        i_max   = Float64(getproperty(cd, :i_max))

        delta_t_load = delta_i / v_load
        nb_loads     = floor(Int, i_max / delta_i)

        # k = 0 → OCV point (i=0);  k = 1..nb_loads → successive load steps
        sample_times = [delta_t_ini + k * (delta_t_load + delta_t_break)
                        for k in 0:nb_loads]
        ifc_values   = Float64[k * delta_i for k in 0:nb_loads]

        indices = [argmin(abs.(t_hist .- t)) for t in sample_times]
        return indices, ifc_values, sample_times

    # ── Calibration polarization profile ──────────────────────────────────────
    elseif hasproperty(cd, :i_exp) && hasproperty(cd, :v_load)
        i_exp  = getproperty(cd, :i_exp)
        v_load = Float64(getproperty(cd, :v_load))

        delta_t_load = abs(Float64(i_exp[1])) / v_load
        delta_t_cali = delta_t_load + delta_t_break

        sample_times = [delta_t_ini + k * delta_t_cali for k in 1:length(i_exp)]
        ifc_values   = collect(Float64, i_exp)

        indices = [argmin(abs.(t_hist .- t)) for t in sample_times]
        return indices, ifc_values, sample_times
    end

    throw(ArgumentError("Unsupported current profile for polarization sampling"))
end


# ─────────────────────────────────────────────────────────────────────────────
# INDIVIDUAL CRITERIA
# ─────────────────────────────────────────────────────────────────────────────

"""
    check_start_voltage_range(Ucell, voltage_range) -> Bool

Return `true` when `Ucell[begin]` lies within `voltage_range`.
"""
function check_start_voltage_range(Ucell::Vector{Float64},
                                   voltage_range::Tuple{Float64, Float64})::Bool
    isempty(Ucell) && return false
    u0 = Ucell[begin]
    isfinite(u0) || return false
    lo, hi = voltage_range
    return lo <= u0 <= hi
end


"""
    check_monotonicity(Ucell, ifc, threshold) -> Bool

Return `true` when every upward step in `Ucell` (as `ifc` increases) is smaller than
`threshold`.
"""
function check_monotonicity(Ucell::Vector{Float64},
                            ifc::Vector{Float64},
                            threshold::Float64)::Bool
    n = length(Ucell)
    (n == 0) && return false
    (n == 1) && return true
    length(ifc) == n || return false

    U_sorted, i_sorted = _sorted_curve(Ucell, ifc)

    prev_i = i_sorted[1]
    prev_u = U_sorted[1]
    (isfinite(prev_i) && isfinite(prev_u)) || return false

    for j in 2:n
        cur_i = i_sorted[j]
        cur_u = U_sorted[j]
        (isfinite(cur_i) && isfinite(cur_u)) || return false

        # Only enforce the constraint when current increases.
        if cur_i > prev_i
            (cur_u - prev_u > threshold) && return false
            prev_i = cur_i
            prev_u = cur_u
        else
            # Same (or decreasing) current: keep the most conservative voltage for the next step.
            prev_u = min(prev_u, cur_u)
        end
    end

    return true
end


"""
    check_positive_voltages(ifc_max_achieved, ifc_max_expected, tol = 0.01) -> Bool

Return `true` when the simulation reached the expected maximum current density.

`ifc_max_achieved` is the last current-density step that was fully simulated;
`ifc_max_expected` is the target value from the current profile.  The criterion
passes when `ifc_max_achieved ≥ ifc_max_expected × (1 − tol)`.

When `ifc_max_expected ≤ 0` the check is vacuously `true`.

This is the overload used by `classify_polarization_curve(simu, cfg)`: because the
simulator now triggers a safety stop before U_cell can go negative, the sign of the
sampled voltages is no longer a reliable indicator of the operating limit.  Comparing
achieved vs. expected current density is the correct proxy.
"""
function check_positive_voltages(ifc_max_achieved::Real,
                                  ifc_max_expected::Real,
                                  tol::Float64 = 0.01)::Bool
    ifc_max_expected <= 0.0 && return true
    return ifc_max_achieved >= ifc_max_expected * (1.0 - tol)
end


"""
    check_positive_voltages(Ucell, min_voltage = 0.0) -> Bool

Return `true` when every element of `Ucell` is ≥ `min_voltage`.

This legacy overload is used by the vector-based `classify_polarization_curve(Ucell, ifc, cfg)`
where no simulation metadata is available.
"""
function check_positive_voltages(Ucell::Vector{Float64},
                                  min_voltage::Float64 = 0.0)::Bool
    isempty(Ucell) && return false
    for u in Ucell
        (isfinite(u) && u >= min_voltage) || return false
    end
    return true
end

end # module ValidityCriteria
