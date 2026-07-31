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

"""
    _sorted_curve(Ucell, ifc) -> (U_sorted, i_sorted)

Sort polarization-curve points by increasing current density.

Validity checks operate on the E-I relationship, so curve points must be consistently ordered
by current. This internal helper ensures robustness to arbitrary input orderings by sorting
both vectors together by current density.

# Arguments
- `Ucell::AbstractVector{<:Real}`: Cell voltages (V). Converted to `Float64`.
- `ifc::AbstractVector{<:Real}`: Current densities (A m⁻²). Converted to `Float64`.

# Returns
- `(U_sorted, i_sorted)::Tuple{Vector{Float64}, Vector{Float64}}`: Voltages and currents
  sorted by increasing current density.

# Example
```julia
Ucell = [0.92, 1.05, 0.65]
ifc = [5e3, 0.0, 1.5e4]
U_s, i_s = _sorted_curve(Ucell, ifc)
# U_s ≈ [1.05, 0.92, 0.65], i_s ≈ [0.0, 5e3, 1.5e4]
```
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
    classify_polarization_curve(Ucell, ifc, cfg = ValidityCriteriaConfig()) -> ValidationResult

Classify a polarization curve sampled as discrete (voltage, current) points.

Evaluates all enabled criteria sequentially. A curve is `:valid` if and only if every
enabled criterion passes; otherwise it is `:invalid`. The result includes per-criterion
flags and a diagnostic message for troubleshooting.

# Arguments
- `Ucell::Vector{Float64}`: Cell voltages at each measurement point (V). Must be non-empty,
  same length as `ifc`, and contain only finite values.
- `ifc::Vector{Float64}`: Current densities at each measurement point (A m⁻²). Must be
  non-empty, same length as `Ucell`, and contain only finite values.
- `cfg::ValidityCriteriaConfig`: Configuration with thresholds and per-criterion activation
  flags. Default: `ValidityCriteriaConfig()` (all criteria enabled with standard thresholds).

# Returns
- `ValidationResult`: Structured outcome containing:
  - `classification::Symbol`: Overall verdict (`:valid` or `:invalid`).
  - `start_in_range::Union{Bool, Nothing}`: Starting-voltage check result (or `nothing`).
  - `is_monotonic::Union{Bool, Nothing}`: Monotonicity check result (or `nothing`).
  - `has_positive_voltages::Union{Bool, Nothing}`: Minimum-voltage check result (or `nothing`).
  - `details::String`: Human-readable summary.

# Notes
- Input vectors are automatically sorted by current density before validation.
- This is the vector-based overload for standalone curve data (no simulator metadata).
- For simulation outputs, use the simulator-based overload `classify_polarization_curve(simu, cfg)`.

# Example
```julia
Ucell  = [1.05, 0.92, 0.80, 0.65, 0.50]
ifc    = [0.0,  5e3,  1e4,  1.5e4, 2e4]
result = classify_polarization_curve(Ucell, ifc)
@assert result.classification == :valid
println(result.details)
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

Classify a polarization curve extracted from an AlphaPEM simulator object.

Automatically extracts discretized curve points (current and stabilized voltage) from the
simulator's solver history and outputs, then applies all enabled validation criteria.
Measurement points are sampled at the *end* of each loading + stabilisation period along
the descent (i_max → 0), ensuring only stabilized data is tested.

# Arguments
- `simu`: AlphaPEM simulator object with the following required properties:
  - `outputs`: Simulation results containing `solver.t` (time history, s) and
    `derived.Ucell` (voltage trajectory, V).
  - `current_density`: Current profile object exposing `t_starts`, `di_transitions`,
    `dt_loads`, `delta_t_breaks` (required for polarization profiles only).
- `cfg::ValidityCriteriaConfig`: Configuration with thresholds and per-criterion activation
  flags. Default: `ValidityCriteriaConfig()` (all criteria enabled with standard thresholds).

# Returns
- `ValidationResult`: Structured outcome (see `classify_polarization_curve(Ucell, ifc, cfg)`).

# Raises
- Returns `:invalid` with diagnostic message if:
  - Simulator has no outputs or simulation has not been run.
  - Current profile is not a polarization-type (missing transition properties).
  - Solver history is empty or malformed.

# Notes
- **Profile requirement**: This overload requires a polarization-type current profile
  (`PolarizationCurrent` or `PolarizationCalibrationCurrent`). Other profile types will
  result in `:invalid` with appropriate error message.
- **Safety mechanism**: The simulator triggers a safety stop if cell voltage would go negative,
  so the third criterion (`has_positive_voltages`) checks *coverage* (how far the descent
  progressed) rather than raw voltage signs.
- **Measurement schedule**: Points are sampled at planned stabilisation end times, but the
  exact indices are found by proximity to solver history, so results are robust to slight
  timing drift.

# Example
```julia
# Assume simu has been run with a polarization profile
result = classify_polarization_curve(simu)
if result.classification == :invalid
    println("Failed: \$(result.details)")
end
```
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
    # So, it is verified if the simulation completed the expected descent down to
    # i≈0: if it stopped early, the operating limit (U_cell → 0) was hit before the
    # full curve could be measured. `sample_times` are in temporal order regardless
    # of whether current is ascending or descending, so coverage is measured as the
    # fraction of descent steps reached rather than as a raw current value (the
    # last sampled current is ~0 by construction, not the peak).
    pos_ok = if cfg.apply_positive_voltages
        t_end        = t_hist[end]
        # Number of descent steps whose sample time was actually reached by the solver.
        # sample_times are analytical; t_end is the exact solver stop time.
        k_last       = something(findlast(st -> st <= t_end, sample_times), 0)
        n_expected   = length(sample_times)
        check_positive_voltages(Float64(k_last), Float64(n_expected), cfg.current_reach_tol)
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

Extract measurement times and current densities from a polarization-type current profile.

Polarization profiles (`PolarizationCurrent`, `PolarizationCalibrationCurrent`) arrange their
current schedule in four phases: (1) initial stabilisation (0 → i_ini), (2) fast ramp to peak
(i_ini → i_max), (3) measured descent (i_max → 0, sampled at each step end), and (4) hysteresis
check (0 → i_max, not used here). This function extracts indices and times for phase 3 only —
the actual polarization curve.

# Arguments
- `t_hist::AbstractVector{<:Real}`: Full time history from the solver (s). Used to locate
  indices nearest to planned sample times.
- `cd`: Current density profile object with properties `t_starts`, `di_transitions`, `dt_loads`,
  `delta_t_breaks`. These define the schedule directly (not reconstructed from i_max/di_step).

# Returns
- `indices::Vector{Int}`: Row indices into `t_hist` closest to each sample time.
- `ifc_values::Vector{Float64}`: Current density at each sampled point (A m⁻²), accumulated
  from `di_transitions` along the descent.
- `sample_times::Vector{Float64}`: Planned sampling times: `t_starts[k] + dt_loads[k] + delta_t_breaks[k]`
  for each descent step (s). Used by `classify_polarization_curve(simu, cfg)` to assess
  how far the solver progressed (current-coverage check).

# Raises
- `ArgumentError`: If the current profile is missing required properties, has fewer than 3
  transitions, or has no descent steps (first ascending transition encountered too early).

# Notes
- Returns only descent steps: stops at the first positive `di_transitions` (phase 4).
- Always indices ≤ t_hist end, even if solver stopped early before reaching planned times.
- This function avoids external module calls, keeping ValidityCriteria self-contained.

# Example
```julia
# Assumes cd has transitions: [i_ini, (i_max-i_ini), -(step1), -(step2), +(ascent)]
indices, ifc, times = _extract_polarization_sampling_indices(t_hist, cd)
# Returns descent times and currents; stops when first positive di_transitions appears
```
"""
function _extract_polarization_sampling_indices(t_hist::AbstractVector{<:Real},
                                                 cd)
    for prop in (:t_starts, :di_transitions, :dt_loads, :delta_t_breaks)
        hasproperty(cd, prop) ||
            throw(ArgumentError("Current profile does not expose `$prop`"))
    end

    t_starts       = getproperty(cd, :t_starts)
    di_transitions = getproperty(cd, :di_transitions)
    dt_loads       = getproperty(cd, :dt_loads)
    delta_t_breaks = getproperty(cd, :delta_t_breaks)

    length(di_transitions) >= 3 ||
        throw(ArgumentError("Polarization profile has no measurement step (expected: " *
                            "initial stabilisation, ramp to i_max, then measured steps)"))

    indices      = Int[]
    ifc_values   = Float64[]
    sample_times = Float64[]

    # i_curr tracks the current density reached after each transition; after steps 1-2
    # (0 -> i_ini -> i_max) it equals i_max, the starting point of the measured descent.
    i_curr = Float64(di_transitions[1]) + Float64(di_transitions[2])

    for k in 3:length(di_transitions)
        di_transitions[k] < 0 || break   # stop at the first ascending step (hysteresis segment)
        i_curr += Float64(di_transitions[k])
        t_sample = Float64(t_starts[k]) + Float64(dt_loads[k]) + Float64(delta_t_breaks[k])
        push!(ifc_values, i_curr)
        push!(sample_times, t_sample)
        push!(indices, argmin(abs.(t_hist .- t_sample)))
    end

    isempty(indices) &&
        throw(ArgumentError("Polarization profile has no descent (measurement) steps"))

    return indices, ifc_values, sample_times
end


# ─────────────────────────────────────────────────────────────────────────────
# INDIVIDUAL CRITERIA
# ─────────────────────────────────────────────────────────────────────────────

"""
    check_start_voltage_range(Ucell, voltage_range) -> Bool

Validate that the initial cell voltage lies within a physically sensible range.

The starting voltage (first point of the polarization curve) is a critical indicator of
numerical stability and physical correctness. Values outside the expected open-circuit range
suggest non-physical equilibrium or solver failure.

# Arguments
- `Ucell::Vector{Float64}`: Cell voltages along the polarization curve (V). Must be non-empty
  and contain finite values.
- `voltage_range::Tuple{Float64, Float64}`: Valid interval `(U_min, U_max)` for the starting
  voltage (V). Typically `(0.0, E₀)` where `E₀` ≈ 1.23 V is the theoretical open-circuit
  voltage of a PEMFC.

# Returns
- `Bool`: `true` if `Ucell[begin]` is finite and lies within `voltage_range` (inclusive);
  `false` otherwise. Returns `false` if `Ucell` is empty or its first value is non-finite.

# Example
```julia
Ucell = [1.05, 0.92, 0.80, 0.65]
result = check_start_voltage_range(Ucell, (0.0, 1.23))  # true
```
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

Validate that cell voltage decreases approximately monotonically with current density.

In a physical fuel cell, increasing current density causes voltage to decrease (E-I curve).
Small upward bumps due to numerical noise are tolerated; larger ones indicate hysteresis,
measurement error, or solver instability.

# Arguments
- `Ucell::Vector{Float64}`: Cell voltages along the polarization curve (V). All values must
  be finite.
- `ifc::Vector{Float64}`: Current densities (A m⁻²) corresponding to each voltage point.
  All values must be finite and same length as `Ucell`.
- `threshold::Float64`: Maximum tolerated upward step in voltage (V) as current increases.
  Typical value: 0.005 (5 mV).

# Returns
- `Bool`: `true` if all upward voltage steps (between consecutive increasing current points)
  are ≤ `threshold` and all values are finite; `false` if the constraint is violated,
  vectors have mismatched lengths, or contain non-finite values.

# Notes
- Vectors are automatically sorted by current density before checking.
- When multiple points have the same current, the minimum voltage is retained for comparison
  with the next distinct current value.

# Example
```julia
Ucell = [1.05, 0.92, 0.80, 0.65]
ifc = [0.0, 5e3, 1e4, 1.5e4]
result = check_monotonicity(Ucell, ifc, 0.005)  # true if no bump > 5 mV
```
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
    check_positive_voltages(achieved, expected, tol = 0.01) -> Bool

Validate that the simulation reached a sufficient fraction of the planned current descent.

In a real measurement, the polarization curve is traced by stepping current down from `i_max`
to 0. If the simulator stops early (safety limit triggered), it signals that the true
negative-voltage limit was reached before the full descent could be sampled. This overload
checks coverage as a fraction of descent *steps* (not raw current), since the descent always
ends near-zero by construction.

# Arguments
- `achieved::Real`: Number of descent steps (or time steps, or other extent measure) that
  the simulation actually completed.
- `expected::Real`: Total number of descent steps (or equivalent extent) that were planned.
- `tol::Float64`: Relative tolerance for step coverage. The check passes when
  `achieved ≥ expected × (1 − tol)`. Default: `0.01` (1 %).

# Returns
- `Bool`: `true` if `achieved ≥ expected × (1 − tol)` or `expected ≤ 0`; `false` otherwise.
  When `expected ≤ 0`, returns `true` (vacuously true condition).

# Notes
- This overload is used by `classify_polarization_curve(simu, cfg)` when full simulator
  metadata is available. The safety stop mechanism prevents voltages from going negative,
  so voltage range alone is insufficient; step coverage is the correct proxy.
- See also the vector-based overload for legacy behavior operating on raw voltages.

# Example
```julia
achieved = 18  # 18 out of 20 planned descent steps completed
expected = 20
result = check_positive_voltages(achieved, expected, 0.01)  # true if >= 19.8 steps
```
"""
function check_positive_voltages(achieved::Real,
                                  expected::Real,
                                  tol::Float64 = 0.01)::Bool
    expected <= 0.0 && return true
    return achieved >= expected * (1.0 - tol)
end


"""
    check_positive_voltages(Ucell, min_voltage = 0.0) -> Bool

Validate that all sampled cell voltages remain above a minimum threshold.

This is a simple bounds check: when operating a fuel cell, voltages should not go negative
under any circumstances. This overload is used when polarization-curve metadata is unavailable
(legacy vector-based path), so direct voltage inspection is the only available criterion.

# Arguments
- `Ucell::Vector{Float64}`: Cell voltages sampled along the polarization curve (V).
  All values must be finite.
- `min_voltage::Float64`: Minimum acceptable voltage (V). Default: `0.0`. Typically set to
  `0.0` (no negative voltages), but can be adjusted for application-specific thresholds.

# Returns
- `Bool`: `true` if all elements of `Ucell` are finite and ≥ `min_voltage`; `false` if
  `Ucell` is empty, any element is non-finite, or any element is < `min_voltage`.

# Notes
- This check is a legacy fallback for the vector-based `classify_polarization_curve(Ucell, ifc, cfg)`.
- When full simulator outputs are available, use the step-coverage overload
  `check_positive_voltages(achieved, expected, tol)` instead, which is more robust.

# Example
```julia
Ucell = [1.05, 0.92, 0.80, 0.65]
result = check_positive_voltages(Ucell, 0.0)  # true (all ≥ 0)

Ucell_bad = [1.05, 0.92, -0.05]
result_bad = check_positive_voltages(Ucell_bad, 0.0)  # false (one < 0)
```
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
