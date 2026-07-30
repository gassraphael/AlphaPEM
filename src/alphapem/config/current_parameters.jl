# -*- coding: utf-8 -*-

"""This file contains the current density parameter structures."""

# ______________________________________________Current density parameters______________________________________________

abstract type AbstractCurrentParams end

"""
    StepParams

Parameters for the step current model.
"""
Base.@kwdef struct StepParams <: AbstractCurrentParams
    delta_t_ini::Float64   = 30.0 * 60.0 # (s). Initial time at zero current density for the stabilisation of the internal states (standard value).
    delta_t_load::Float64  = 30.0        # (s). Loading time for the step current density function, from 0 to i_step.
    delta_t_break::Float64 = 2.0 * 60.0  # (s). Time at i_step current density for the stabilisation of the internal states.
    i_ini::Float64         = 1.0e4       # (A.m-2). Initial current density for the step current density function.
    i_step::Float64        = 2.0e4       # (A.m-2). Current density for the step current density function.
end


"""
    PolarizationParams

Parameters for polarization curve simulation.
Priority for i_max (all values are rounded to nearest di_step multiple):
1. If user provides i_max explicitly (not NaN), use that value
2. Else if fuel_cell has experimental data, use max(i_exp)
3. Else default to 2.5 A/cm²
"""
Base.@kwdef struct PolarizationParams <: AbstractCurrentParams
    delta_t_ini::Float64             = 30 * 60.0 # (s). Initial time at zero current density for the stabilisation of the internal states.
    i_ini::Float64                   = 1.0e4     # (A.m-2). Initial current density for stabilization (before measurement steps begin).
    di_step::Float64                 = 0.05e4    # (A.m-2). Nominal current density step for the polarisation current density function.
    v_load::Float64                  = 0.01e4    # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break::Float64           = 300       # (s). Breaking time for one step current, for the stabilisation of the internal states and the measurement of the voltage.
    delta_t_break_lower_0_3::Float64 = 120.0     # (s). Breaking time for one step current when i_fc <= 0.3 A/cm².
    delta_t_break_OCV ::Float64      = 60.0      # (s). Breaking time for one step current when i_fc = 0 A/cm².
    delta_t_measurement::Float64     = 30.0      # (s). Time for the measurement of the voltage for each step current.
    i_max::Float64                   = NaN       # User-provided maximum current (A/m²). NaN means use experimental or default.
end


"""
    PolarizationCalibrationParams
"""
Base.@kwdef struct PolarizationCalibrationParams <: AbstractCurrentParams
    delta_t_ini::Float64             = 30 * 60.0  # (s). Initial time at zero current density for the stabilisation of the internal states.
    i_ini::Float64                   = 1.0e4      # (A.m-2). Initial current density for stabilization (before measurement steps begin).
    v_load::Float64                  = 0.01e4     # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break::Float64           = 300        # (s). Breaking time for one step current, for the stabilisation of the internal states and the measurement of the voltage.
    delta_t_break_lower_0_3::Float64 = 120.0      # (s). Breaking time for one step current when i_fc <= 0.3 A/cm².
    delta_t_break_OCV ::Float64      = 60.0       # (s). Breaking time for one step current when i_fc = 0 A/cm².
    delta_t_measurement::Float64     = 30.0       # (s). Time for the measurement of the voltage for each step current (except at i_fc < 0.3 V where a dedicated time is given for safety).
    i_exp::Vector{Float64}           = Float64[]  # Experimental current density values (A/m²) for calibration
end


"""
    EISParams
"""
Base.@kwdef struct EISParams <: AbstractCurrentParams
    i_EIS::Float64       = 1.5e4        # (A/m²). Parameters for the EIS curve.
    ratio::Float64       = 5.0 / 100.0  # (.). Parameters for the EIS curve.
    f_power_min::Float64 = -3.0         # (.). Power of the minimum frequency for the EIS current density function.
    f_power_max::Float64 = 5.0          # (.). Power of the maximum frequency for the EIS current density function.
    nb_f::Int            = 90           # (.). Number of frequencies tested for the EIS current density function.
    nb_points::Int       = 50           # (.). Number of points calculated per specific period for the EIS current density function.
end