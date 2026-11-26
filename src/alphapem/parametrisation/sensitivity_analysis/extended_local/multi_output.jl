# -*- coding: utf-8 -*-

# ============================================================
# MULTI-OUTPUT MEASUREMENT SETS (GENERAL STRUCTURES ONLY)
# ============================================================
#
# Context. Only the cell voltage is analysed today. The reference studies how adding further
# measurements (high-frequency resistance, membrane water crossover) separates parameters that are
# collinear in voltage alone: the sensitivity vectors of the chosen outputs are concatenated and the
# collinearity index is recomputed on the longer vector. This file only fixes the vocabulary and the
# entry points of that extension; the computation is to be written once AlphaPEM exposes the extra
# outputs.

"""
    FUTURE_OUTPUTS

Outputs of the reference that AlphaPEM does not expose yet, so they cannot be simulated or analysed:
the high-frequency resistance (`:hfr`) and the membrane water crossover (`:crossover`). They may
already be named in a [`MeasurementSet`](@ref) but are not part of `SUPPORTED_OUTPUTS`.
"""
const FUTURE_OUTPUTS = (:hfr, :crossover)


"""
    MeasurementSet(name, outputs, weighting = :unit)

Group of outputs whose sensitivity vectors are concatenated before computing collinearity.

The weighting rescales each output's block of the concatenated vector, because outputs differ in
unit and in number of measurements:

- `:unit`: every output has weight 1, as in the reference (default).
- `:count_ratio`: each output is weighted by the ratio of the voltage measurement count to its own
  count, so that no output dominates through sheer number of points (an `N_V/N_HFR` weighting
  that is not part of the article).
- `Dict{Symbol, Float64}`: explicit positive weight per output.

Parameters
----------
name : Symbol
    Label of the set, used in figures and exports (for example `:V_HFR`).
outputs : Vector{Symbol}
    Distinct outputs, among `SUPPORTED_OUTPUTS` and `FUTURE_OUTPUTS`, always including `:voltage`.
weighting : Symbol or Dict{Symbol, Float64}, optional
    Weighting rule described above. Default is `:unit`.

Raises
------
ArgumentError
    If the outputs are empty, duplicated, unknown or lack `:voltage`, or the weighting is invalid.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
struct MeasurementSet
    name::Symbol
    outputs::Vector{Symbol}
    weighting::Union{Symbol, Dict{Symbol, Float64}}

    function MeasurementSet(name::Symbol, outputs::Vector{Symbol}, weighting = :unit)
        known = (SUPPORTED_OUTPUTS..., FUTURE_OUTPUTS...)
        !isempty(outputs) && length(unique(outputs)) == length(outputs) &&
            all(output -> output in known, outputs) && :voltage in outputs ||
            throw(ArgumentError("outputs must be distinct, among $(known), and include :voltage"))
        valid_weighting = weighting in (:unit, :count_ratio) ||
            (weighting isa AbstractDict && Set(keys(weighting)) == Set(outputs) &&
             all(w -> isfinite(w) && w > 0, values(weighting)))
        valid_weighting || throw(ArgumentError(
            "weighting must be :unit, :count_ratio, or a positive weight for every output"))
        return new(name, outputs,
                   weighting isa AbstractDict ? Dict{Symbol, Float64}(weighting) : weighting)
    end
end


"""
    DEFAULT_MEASUREMENT_SETS

Measurement sets compared by the reference: voltage alone, then voltage with HFR, with crossover,
and with both. Only the first is computable today.
"""
const DEFAULT_MEASUREMENT_SETS = (
    MeasurementSet(:V, [:voltage]),
    MeasurementSet(:V_HFR, [:voltage, :hfr]),
    MeasurementSet(:V_crossover, [:voltage, :crossover]),
    MeasurementSet(:V_HFR_crossover, [:voltage, :hfr, :crossover]),
)


"""
    combine_outputs(result, measurement_set) -> Array{Float64,3}

Concatenate and weight the sensitivity vectors of a measurement set (to be implemented).

Planned behaviour: for each candidate, stack the per-output sensitivity blocks of every parameter
along the measurement axis, apply the set's weighting, and respect each output's validity mask so
that missing values are never treated as zero. The result is shaped like the input of
`collinearity_indices` (`n_candidates × n_params × n_measurements`).

Parameters
----------
result : ExtendedLocalSimulationResult
    Completed sensitivity simulations holding every output of the set.
measurement_set : MeasurementSet
    Outputs to concatenate and their weighting.

Raises
------
ErrorException
    Always, until the extension is implemented.
"""
function combine_outputs(::ExtendedLocalSimulationResult, ::MeasurementSet)
    error("combine_outputs is not implemented: only the voltage output is available")
end


"""
    compare_measurement_sets(result, measurement_sets) -> Dict{Symbol, CollinearitySummary}

Compute the collinearity summary of each measurement set for comparison (to be implemented).

Planned behaviour: apply [`combine_outputs`](@ref) and the collinearity summary to every set so that
the reduction of a pair's collinearity when measurements are added can be plotted side by side,
keeping the same parameters and thresholds across sets.

Parameters
----------
result : ExtendedLocalSimulationResult
    Completed sensitivity simulations holding every output of the sets.
measurement_sets : AbstractVector{MeasurementSet}
    Sets to compare, in display order.

Raises
------
ErrorException
    Always, until the extension is implemented.
"""
function compare_measurement_sets(::ExtendedLocalSimulationResult, ::AbstractVector{MeasurementSet})
    error("compare_measurement_sets is not implemented: only the voltage output is available")
end
