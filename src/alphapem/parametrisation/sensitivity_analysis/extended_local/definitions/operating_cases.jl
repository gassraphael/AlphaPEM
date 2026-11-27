# -*- coding: utf-8 -*-

"""
    OperatingCase(label, operating_conditions)

One operating condition in the factorial library used by the extended local sensitivity analysis.

# Fields
- `label::String`: Stable identifier in factorial product order.
- `operating_conditions::OperatingConditions`: Operating conditions for this case.
"""
struct OperatingCase
    label::String
    operating_conditions::OperatingConditions
end


"""
    _effective_operating_condition_bounds(overrides) -> Dict

Merge partial operating-condition bound overrides with AlphaPEM's default bounds.

Parameters
----------
overrides : AbstractDict
    Partial mapping from operating-condition names to `(min, max, type)` tuples.

Returns
-------
Dict{Symbol,Tuple{Float64,Float64,Symbol}}
    A fresh dictionary containing the effective bounds.

Raises
------
ArgumentError
    If an override has an unknown name, invalid range, or incompatible type.
"""
function _effective_operating_condition_bounds(overrides::AbstractDict)
    # Keep the shared defaults intact: each call resolves its own effective bound set.
    bounds = copy(OPERATING_CONDITIONS_BOUNDS)

    # Overrides are partial, so every unspecified operating condition keeps its default range.
    for (name, override) in overrides
        haskey(bounds, name) ||
            throw(ArgumentError("Unknown operating condition in bounds override: $name"))
        override isa Tuple && length(override) == 3 ||
            throw(ArgumentError("Bounds override for $name must be a (min, max, type) tuple"))

        lo, hi, kind = override
        lo isa Real && hi isa Real && kind isa Symbol ||
            throw(ArgumentError("Bounds override for $name must contain two real bounds and a Symbol type"))
        lo_f = Float64(lo)
        hi_f = Float64(hi)
        isfinite(lo_f) && isfinite(hi_f) && lo_f < hi_f ||
            throw(ArgumentError("Bounds override for $name must be finite and have min < max"))
        kind == bounds[name][3] ||
            throw(ArgumentError("Bounds override for $name must retain type $(bounds[name][3])"))

        bounds[name] = (lo_f, hi_f, kind)
    end
    return bounds
end


"""
    _validate_operating_case_conditions(oc, bounds, pressure_constraint, label)

Validate one generated operating case against effective bounds and physical constraints.

Parameters
----------
oc : OperatingConditions
    Generated operating condition.
bounds : AbstractDict
    Effective bounds for bounded operating conditions.
pressure_constraint : OperatingConditionConstraint
    Exact anode-to-cathode pressure-gap constraint.
label : AbstractString
    Case label included in validation errors.

Returns
-------
Nothing
    Returns `nothing` when all checks pass.

Raises
------
ArgumentError
    If the case violates an effective bound, the pressure gap, or a basic physical constraint.
"""
function _validate_operating_case_conditions(
    oc::OperatingConditions,
    bounds::AbstractDict,
    pressure_constraint::OperatingConditionConstraint,
    label::AbstractString,
)
    # Reject invalid numeric data before applying range and physical checks.
    all(isfinite(getfield(oc, field)) for field in fieldnames(OperatingConditions)) ||
        throw(ArgumentError("Factorial case $label contains a non-finite operating condition"))

    # Check all effective bounds without silently clamping generated cases.
    for (name, (lo, hi, _)) in bounds
        value = getfield(oc, name)
        lo <= value <= hi ||
            throw(ArgumentError("Factorial case $label has $name=$value outside effective bounds [$lo, $hi]"))
    end

    # These physical invariants complement the configured ranges; bounds alone do not guarantee
    # that a generated operating point is meaningful for the fuel-cell model.
    oc.T_des > 273.15 ||
        throw(ArgumentError("Factorial case $label must have T_des > 273.15 K"))
    oc.Pa_des > 0.0 && oc.Pc_des > 0.0 ||
        throw(ArgumentError("Factorial case $label must have positive anode and cathode pressures"))
    oc.Sa >= 1.0 && oc.Sc >= 1.0 ||
        throw(ArgumentError("Factorial case $label must have stoichiometric ratios >= 1"))
    0.0 <= oc.Phi_a_des <= 1.0 && 0.0 <= oc.Phi_c_des <= 1.0 ||
        throw(ArgumentError("Factorial case $label must have relative humidities in [0, 1]"))
    0.0 <= oc.y_H2_in <= 1.0 ||
        throw(ArgumentError("Factorial case $label must have y_H2_in in [0, 1]"))
    oc.i_min_stoich >= 0.0 ||
        throw(ArgumentError("Factorial case $label must have i_min_stoich >= 0"))

    # Pc_des is dependent on Pa_des, so verify the exact relation rather than treating it as
    # an independent parameterized condition.
    expected_pc = pressure_constraint.fn(Dict(:Pa_des => oc.Pa_des))
    isapprox(oc.Pc_des, expected_pc; rtol = 0.0, atol = 1e-9) ||
        throw(ArgumentError("Factorial case $label does not satisfy the exact pressure-gap constraint"))
    return nothing
end


"""
    generate_operating_case_library(base::OperatingConditions; operating_bounds_overrides, pressure_gap)
        -> Vector{OperatingCase}

Generate all the three-level factorial operating cases from a base condition.

Each selected operating condition takes its effective lower bound, midpoint, and upper bound.
The cathode pressure is then derived from the exact anode-to-cathode pressure gap. Cases outside
the effective bounds are rejected rather than silently modified.

Parameters
----------
base : OperatingConditions
    Nominal operating conditions; fields not parameterized, such as `i_min_stoich`, are retained.
operating_bounds_overrides : AbstractDict
    Partial overrides to `OPERATING_CONDITIONS_BOUNDS`. Unspecified bounds retain their defaults.
pressure_gap : Real
    Exact difference `Pa_des - Pc_des`, in pascals. Default `0.2e5`.

Returns
-------
Vector{OperatingCase}
    The complete factorial library in deterministic product order.

Raises
------
ArgumentError
    If an override is invalid or any generated case violates its effective bounds or physical
    operating-condition constraints.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function generate_operating_case_library(
    base::OperatingConditions;
    operating_bounds_overrides::AbstractDict = Dict{Symbol,Tuple{Float64,Float64,Symbol}}(),
    pressure_gap::Real = 0.2e5,
)::Vector{OperatingCase}
    # === Resolve the bounds and pressure relationship for this library ===
    # Merge the caller's partial overrides with the default operating-condition bounds.
    bounds = _effective_operating_condition_bounds(operating_bounds_overrides)
    # Create the exact rule that derives cathode pressure from anode pressure.
    pressure_constraint = exact_pressure_gap_constraint(pressure_gap)

    # === Calculate the three design values for each parameterized condition ===
    # Prepare a mapping from each condition name to its three design values.
    operating_condition_levels = Dict{Symbol,NTuple{3,Float64}}()
    # Visit each condition that is varied in the factorial parameterization.
    for name in OPERATING_CONDITIONS_NAMES_FOR_PARAMETRISATION
        # Read this condition's effective lower and upper bounds; ignore the type field here.
        lo, hi, _ = bounds[name]
        # Store its lower bound, midpoint, and upper bound in that order.
        operating_condition_levels[name] = (lo, (lo + hi) / 2.0, hi)
    end

    # === Enumerate every combination of the three design values ===
    # Prepare the result vector that will hold cases in product order.
    cases = OperatingCase[]
    # Create one range of level positions (1, 2, 3) for every parameterized condition.
    # Position 1 selects the minimum, 2 the midpoint, and 3 the maximum.
    level_index_ranges = ntuple(
        _ -> 1:3,
        length(OPERATING_CONDITIONS_NAMES_FOR_PARAMETRISATION),
    )
    # `...` passes each range as a separate input to `Iterators.product`.
    # For two ranges `1:3`, it would yield `(1, 1)`, `(2, 1)`, `(3, 1)`, `(1, 2)`, and so on:
    # every possible choice of one position from each range.
    # With six ranges, it lazily yields 3^6 = 729 tuples, one for each combination of levels.
    level_indices = Iterators.product(level_index_ranges...)
    # Process each tuple and assign it a stable, one-based case number.
    for (case_index, indices) in enumerate(level_indices)
        # Start a mapping for this case's selected operating-condition values.
        operating_condition_values = Dict{Symbol,Float64}()
        # Visit conditions in the same order used to build the Cartesian product.
        for (condition_index, name) in enumerate(OPERATING_CONDITIONS_NAMES_FOR_PARAMETRISATION)
            # Use this condition's position to select its minimum, midpoint, or maximum value.
            operating_condition_values[name] =
                operating_condition_levels[name][indices[condition_index]]
        end

        # Derive Pc_des from the selected Pa_des; cathode pressure is not independently varied.
        operating_condition_values[:Pc_des] = pressure_constraint.fn(operating_condition_values)

        # Reconstruct a complete OperatingConditions value by visiting every struct field.
        conditions = OperatingConditions(;
            # Use generated values for selected or derived fields; otherwise retain the base value.
            (field => get(operating_condition_values, field, getfield(base, field))
             # Include all fields so the reconstructed object is complete.
             for field in fieldnames(OperatingConditions))...
        )
        # Format the sequential case number as a stable three-digit identifier.
        label = "factorial_" * lpad(string(case_index), 3, '0')

        # Reject the case if it violates effective bounds or physical constraints.
        _validate_operating_case_conditions(conditions, bounds, pressure_constraint, label)
        # Add the validated case to the library.
        push!(cases, OperatingCase(label, conditions))
    end

    # Return the complete library after all combinations have been generated and validated.
    return cases
end


"""
    generate_operating_case_library(cfg::ExtendedLocalSensitivityConfig) -> Vector{OperatingCase}

Build the complete factorial library using the selected fuel cell's nominal operating conditions
unless an explicit base override is supplied in `cfg`.

Parameters
----------
cfg : ExtendedLocalSensitivityConfig
    Configuration with the selected fuel cell, optional base and bound overrides, and exact pressure gap.

Returns
-------
Vector{OperatingCase}
    The 729-case factorial library.

Raises
------
ArgumentError
    If the configuration or any generated operating case is invalid.

References
----------
A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
fuel cell models—part i: Sensitivity analysis and parameter identifiability," Journal of The
Electrochemical Society, vol. 167, p. 044504, Mar. 2020. https://doi.org/10.1149/1945-7111/ab7091.
"""
function generate_operating_case_library(cfg::ExtendedLocalSensitivityConfig)::Vector{OperatingCase}
    # === Resolve base operating conditions from the configuration ===
    # Validate the configuration before using its fields to generate cases.
    validate_config(cfg)
    # Read the optional base conditions supplied by the caller.
    base = cfg.base_operating_conditions
    # Use the selected fuel cell's nominal conditions only when no base override was supplied.
    if base === nothing
        # Create the fuel cell selected by type, voltage zone, year, and gas-channel resolution.
        fuel_cell = create_fuelcell(
            cfg.fuel_cell_type, cfg.voltage_zone; year = cfg.year, nb_gc = cfg.nb_gc)
        # Extract its nominal operating conditions as the base for every generated case.
        base = fuel_cell.operating_conditions
    end

    # === Delegate generation to the method that accepts explicit base conditions ===
    # Call the base-conditions overload so there is only one implementation of case generation.
    return generate_operating_case_library(
        # Pass the resolved base operating conditions.
        base;
        # Forward the user's partial bound overrides.
        operating_bounds_overrides = cfg.operating_bounds_overrides,
        # Forward the exact anode-to-cathode pressure difference.
        pressure_gap = cfg.pressure_gap,
    )
end
