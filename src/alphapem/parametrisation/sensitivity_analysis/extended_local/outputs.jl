# -*- coding: utf-8 -*-

# ============================================================
# ASCENDING CURRENT PROFILE
# ============================================================

"""
    ascending_polarization_params(cfg::ExtendedLocalSensitivityConfig) -> AscendingPolarizationParams

Build the ascending current profile parameters from the analysis configuration.

The requested current densities come from `cfg.current_density_grid` (A/cm²) and are converted to
A/m², the unit used by AlphaPEM's current profiles. The durations come from `cfg`.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Configuration holding the current-density grid and profile durations.

Returns
-------
AscendingPolarizationParams
    Parameters of one ascending sweep, without return pass or hysteresis.
"""
function ascending_polarization_params(cfg::ExtendedLocalSensitivityConfig)::AscendingPolarizationParams
    return AscendingPolarizationParams(
        delta_t_ini = cfg.delta_t_ini,
        v_load = cfg.v_load,
        delta_t_break = cfg.delta_t_break,
        delta_t_break_lower_0_3 = cfg.delta_t_break_lower_0_3,
        delta_t_break_OCV = cfg.delta_t_break_OCV,
        delta_t_measurement = cfg.delta_t_measurement,
        i_grid = cfg.current_density_grid .* 1e4,
    )
end


# ============================================================
# VOLTAGE CURVE
# ============================================================

"""
    VoltageCurve

Cell voltage measured at the requested current densities, with its validity mask.

A point is valid only when the solver completed its whole stabilisation period and the measured
voltage is finite. The sweep is ascending and never skips a point: once a point is missing, all
higher current densities are missing too. Missing voltages are `NaN`; they are neither
interpolated nor imputed.

# Fields
- `current_density::Vector{Float64}`: Requested current densities (A/cm²).
- `voltage::Vector{Float64}`: Measured voltages (V), `NaN` where the point is missing.
- `valid::BitVector`: `true` where the voltage was measured.
- `last_reached_current_density::Float64`: Highest valid current density (A/cm²), `NaN` if none.
- `status::Symbol`: `:complete` (all points valid), `:partial` (a valid prefix) or `:invalid`
  (no valid point).
"""
struct VoltageCurve
    current_density::Vector{Float64}
    voltage::Vector{Float64}
    valid::BitVector
    last_reached_current_density::Float64
    status::Symbol

    function VoltageCurve(current_density::Vector{Float64}, voltage::Vector{Float64},
                          valid::BitVector)
        length(voltage) == length(current_density) == length(valid) ||
            throw(ArgumentError("current_density, voltage and valid must have the same length"))
        last_valid = findlast(valid)
        status = last_valid === nothing ? :invalid : (all(valid) ? :complete : :partial)
        last_reached = last_valid === nothing ? NaN : current_density[last_valid]
        return new(current_density, voltage, valid, last_reached, status)
    end
end


"""
    voltage_curve_from_history(t, U, profile::AscendingPolarizationCurrent) -> VoltageCurve

Average the cell voltage over the measurement window of each current density.

A window is used only if the solver reached its end: the stabilisation interrupted by a solver
stop, and all the following ones, are marked missing.

Parameters
----------
t : AbstractVector{<:Real}
    Increasing solver time history (s).
U : AbstractVector{<:Real}
    Cell voltage at each time of `t` (V).
profile : AscendingPolarizationCurrent
    Current profile that defines the requested current densities and measurement windows.

Returns
-------
VoltageCurve
    Voltages at the requested current densities, with their validity mask.

Raises
------
ArgumentError
    If `t` and `U` have different lengths.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function voltage_curve_from_history(t::AbstractVector{<:Real}, U::AbstractVector{<:Real},
                                    profile::AscendingPolarizationCurrent)::VoltageCurve
    length(t) == length(U) || throw(ArgumentError("t and U must have the same length"))

    n_points = length(profile.i_targets)
    voltage = fill(NaN, n_points)
    valid = falses(n_points)

    # Tolerance on window edges, since the solver stops exactly at the end of each profile step.
    tol = 1e-6
    t_reached = isempty(t) ? -Inf : t[end]

    for k in 1:n_points
        # The sweep never skips a missing point to resume at a higher current density.
        k > 1 && !valid[k - 1] && break

        # The solver stopped before (or during) this stabilisation: its window is incomplete.
        t_reached >= profile.t_meas_ends[k] - tol || break

        # Samples inside the measurement window (t is increasing).
        first_sample = searchsortedfirst(t, profile.t_meas_starts[k] - tol)
        last_sample = searchsortedlast(t, profile.t_meas_ends[k] + tol)
        first_sample <= last_sample || break
        window_voltage = @view U[first_sample:last_sample]
        all(isfinite, window_voltage) || break

        voltage[k] = mean(window_voltage)
        valid[k] = true
    end

    return VoltageCurve(profile.i_targets ./ 1e4, voltage, valid)
end


"""
    extract_voltage_curve(simu) -> VoltageCurve

Extract the cell voltage at the requested current densities from a finished simulation.

The voltages are averaged over each measurement window, without interpolation. A simulation
stopped by the solver (safety stop or timeout) keeps the points reached before the stop. A
simulation without outputs gives an `:invalid` curve.

Parameters
----------
simu : AlphaPEM
    Simulator run with an `AscendingPolarizationParams` current profile.

Returns
-------
VoltageCurve
    Voltages at the requested current densities, with their validity mask.

Raises
------
ArgumentError
    If the simulator does not use an ascending polarization current profile.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function extract_voltage_curve(simu)::VoltageCurve
    profile = simu.current_density
    profile isa AscendingPolarizationCurrent ||
        throw(ArgumentError("extract_voltage_curve requires an AscendingPolarizationParams simulation"))

    # A run without outputs has no measured point.
    if simu.outputs === nothing
        n_points = length(profile.i_targets)
        return VoltageCurve(profile.i_targets ./ 1e4, fill(NaN, n_points), falses(n_points))
    end

    return voltage_curve_from_history(
        time_history(simu.outputs), derived_outputs(simu.outputs).Ucell, profile)
end
