# -*- coding: utf-8 -*-

"""
Module SimulationConfigModule

This module defines the configuration structure for fuel cell simulations in AlphaPEM.

Structures:
- SimulationConfig: Structure containing the simulation configuration parameters.

Functions:
- validate_config(cfg::SimulationConfig): Validates the configuration parameter values.

Constants:
- ALLOWED_CURRENT_TYPES: Allowed current parameter types.
- ALLOWED_VOLTAGE_ZONE: Allowed voltage zones.
- ALLOWED_DISPLAY: Allowed display types.
- ALLOWED_DISPLAY_TIMING: Allowed display timing modes.
- ALLOWED_PURGE: Allowed purge types.
- ALLOWED_AUXILIARY: Allowed auxiliary system types.
"""

# _______________________________________________Simulation configuration_______________________________________________

module SimulationConfigModule

using ..Config: AbstractCurrentParams, StepParams, PolarizationParams,
                 PolarizationCalibrationParams, EISParams, NumericalParams
using ..StateScalingModule: StateScaling

export SimulationConfig, validate_config

Base.@kwdef mutable struct SimulationConfig{T<:AbstractCurrentParams}
    type_fuel_cell::Symbol = :ZSW_GenStack
    type_current::T = PolarizationParams()
    numerical_parameters::NumericalParams = NumericalParams()
    voltage_zone::Symbol = :full
    type_auxiliary::Symbol = :no_auxiliary
    type_flow::Symbol = :counter_flow
    type_purge::Symbol = :no_purge
    type_display::Symbol = :synthetic
    display_timing::Symbol=:postrun
    state_scaling::StateScaling = StateScaling()
end

# --- Allowed values (tu peux enrichir plus tard) ---

const ALLOWED_CURRENT_TYPES = (
    StepParams,
    PolarizationParams,
    PolarizationCalibrationParams,
    EISParams
)

const ALLOWED_VOLTAGE_ZONE = (
    :before_voltage_drop,
    :full
)

const ALLOWED_AUXILIARY = (
    :forced_convective_cathode_with_anodic_recirculation,
    :forced_convective_cathode_with_flow_through_anode,
    :no_auxiliary
)

const ALLOWED_FLOW = (
    :co_flow,
    :counter_flow,
)

const ALLOWED_PURGE = (
    :constant_purge,
    :periodic_purge,
    :no_purge
)

const ALLOWED_DISPLAY = (
    :multiple,
    :synthetic,
    :no_display
)

const ALLOWED_DISPLAY_TIMING = (
    :live,
    :postrun,
)

# --- Validation ---

function validate_config(cfg::SimulationConfig)

    any(T -> cfg.type_current isa T, ALLOWED_CURRENT_TYPES) ||
        error("Invalid type_current: $(typeof(cfg.type_current))")

    cfg.voltage_zone in ALLOWED_VOLTAGE_ZONE ||
        error("Invalid voltage_zone: $(cfg.voltage_zone)")

    cfg.type_auxiliary in ALLOWED_AUXILIARY ||
        error("Invalid type_auxiliary: $(cfg.type_auxiliary)")

    cfg.type_flow in ALLOWED_FLOW ||
        error("Invalid type_flow: $(cfg.type_flow)")

    cfg.type_purge in ALLOWED_PURGE ||
        error("Invalid type_purge: $(cfg.type_purge)")

    cfg.type_display in ALLOWED_DISPLAY ||
        error("Invalid type_display: $(cfg.type_display)")

    cfg.display_timing in ALLOWED_DISPLAY_TIMING ||
        error("Invalid display_timing: $(cfg.display_timing)")

    cell_scaling = cfg.state_scaling.cell
    manifold_scaling = cfg.state_scaling.manifold
    auxiliary_scaling = cfg.state_scaling.auxiliary

    for (name, value) in (
        ("numerical_parameters.nb_gc", cfg.numerical_parameters.nb_gc),
        ("numerical_parameters.nb_gdl", cfg.numerical_parameters.nb_gdl),
        ("numerical_parameters.nb_mpl", cfg.numerical_parameters.nb_mpl),
        ("numerical_parameters.nb_man", cfg.numerical_parameters.nb_man),
        ("numerical_parameters.purge_time", cfg.numerical_parameters.purge_time),
        ("numerical_parameters.delta_purge", cfg.numerical_parameters.delta_purge),
        ("numerical_parameters.delta_t_dyn_step", cfg.numerical_parameters.delta_t_dyn_step),
        ("numerical_parameters.rtol", cfg.numerical_parameters.rtol),
        ("numerical_parameters.atol", cfg.numerical_parameters.atol),
        ("state_scaling.cell.C_v", cell_scaling.C_v),
        ("state_scaling.cell.C_H2", cell_scaling.C_H2),
        ("state_scaling.cell.C_O2", cell_scaling.C_O2),
        ("state_scaling.cell.C_N2", cell_scaling.C_N2),
        ("state_scaling.cell.T", cell_scaling.T),
        ("state_scaling.cell.lambda", cell_scaling.lambda),
        ("state_scaling.cell.eta_c", cell_scaling.eta_c),
        ("state_scaling.cell.s", cell_scaling.s),
        ("state_scaling.manifold.P", manifold_scaling.P),
        ("state_scaling.manifold.Phi", manifold_scaling.Phi),
        ("state_scaling.auxiliary.Wcp", auxiliary_scaling.Wcp),
        ("state_scaling.auxiliary.Wa_inj", auxiliary_scaling.Wa_inj),
        ("state_scaling.auxiliary.Wc_inj", auxiliary_scaling.Wc_inj),
        ("state_scaling.auxiliary.Abp_a", auxiliary_scaling.Abp_a),
        ("state_scaling.auxiliary.Abp_c", auxiliary_scaling.Abp_c),
    )
        value > 0.0 || error("Invalid $name: $value (must be > 0)")
    end

    cfg.numerical_parameters.nb_gc >= 1 || error("Invalid numerical_parameters.nb_gc: must be >= 1")
    cfg.numerical_parameters.nb_gdl >= 1 || error("Invalid numerical_parameters.nb_gdl: must be >= 1")
    cfg.numerical_parameters.nb_mpl >= 1 || error("Invalid numerical_parameters.nb_mpl: must be >= 1")
    cfg.numerical_parameters.nb_man >= 1 || error("Invalid numerical_parameters.nb_man: must be >= 1")

    return cfg
end

end # module