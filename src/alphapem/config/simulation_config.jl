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
- ALLOWED_PLOT: Allowed plot types.
- ALLOWED_PURGE: Allowed purge types.
- ALLOWED_AUXILIARY: Allowed auxiliary system types.
"""

# _______________________________________________Simulation configuration_______________________________________________

module SimulationConfigModule

using ..Config: AbstractCurrentParams, StepParams, PolarizationParams,
                 PolarizationCalibrationParams, EISParams

export SimulationConfig, validate_config

Base.@kwdef mutable struct SimulationConfig{T<:AbstractCurrentParams}
    type_fuel_cell::Symbol = :ZSW_GenStack
    type_current::T = PolarizationParams()
    voltage_zone::Symbol = :full
    type_auxiliary::Symbol = :no_auxiliary
    type_purge::Symbol = :no_purge
    type_display::Symbol = :synthetic
    type_plot::Symbol = :fixed
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

const ALLOWED_PLOT = (
    :dynamic,
    :fixed
)

# --- Validation ---

function validate_config(cfg::SimulationConfig)

    any(T -> cfg.type_current isa T, ALLOWED_CURRENT_TYPES) ||
        error("Invalid type_current: $(typeof(cfg.type_current))")

    cfg.voltage_zone in ALLOWED_VOLTAGE_ZONE ||
        error("Invalid voltage_zone: $(cfg.voltage_zone)")

    cfg.type_auxiliary in ALLOWED_AUXILIARY ||
        error("Invalid type_auxiliary: $(cfg.type_auxiliary)")

    cfg.type_purge in ALLOWED_PURGE ||
        error("Invalid type_purge: $(cfg.type_purge)")

    cfg.type_display in ALLOWED_DISPLAY ||
        error("Invalid type_display: $(cfg.type_display)")

    cfg.type_plot in ALLOWED_PLOT ||
        error("Invalid type_plot: $(cfg.type_plot)")

    return cfg
end

end # module