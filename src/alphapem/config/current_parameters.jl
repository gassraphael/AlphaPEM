# -*- coding: utf-8 -*-

"""This file contains the current density parameter structures."""

# ______________________________________________Current density parameters______________________________________________

abstract type AbstractCurrentParams end

"""
    StepParams

Parameters for the step current model.
"""
Base.@kwdef struct StepParams <: AbstractCurrentParams
    delta_t_ini::Float64 = 30.0 * 60.0 # (s). Initial time at zero current density for the stabilisation of the internal states (standard value).
    delta_t_load::Float64 = 30.0 # (s). Loading time for the step current density function, from 0 to i_step.
    delta_t_break::Float64 = 2.0 * 60.0  # (s). Time at i_step current density for the stabilisation of the internal states.
    i_ini::Float64 = 1.0e4 # (A.m-2). Initial current density for the step current density function.
    i_step::Float64 = 2.0e4 # (A.m-2). Current density for the step current density function.
end


"""
    PolarizationParams
"""
Base.@kwdef struct PolarizationParams <: AbstractCurrentParams
    delta_t_ini::Float64 = 120.0 * 60.0 # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_i::Float64 = 0.05e4  # (A.m-2). Current density step for the polarisation current density function.
    v_load::Float64 = 0.01e4  # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break::Float64 = 15.0 * 60.0 # (s). Breaking time for one step current, for the stabilisation of the internal states.
    i_max::Float64 = 3.0e4 # Maximum current (default value, can be overridden by experimental current values if provided).
    type_fuel_cell::Union{Nothing, String} = nothing  # Optional
    voltage_zone::Union{Nothing, String} = nothing   # Optional
end


"""
    PolarizationCalibrationParams
"""
Base.@kwdef struct PolarizationCalibrationParams <: AbstractCurrentParams
    delta_t_ini::Float64 = 120.0 * 60.0 # (s). Initial time at zero current density for the stabilisation of the internal states.
    v_load::Float64 = 0.01e4  # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break::Float64 = 10.0 * 60.0  # (s). Breaking time for one step current, for the stabilisation of the internal states.
    type_fuel_cell::String
    voltage_zone::String
end


"""
    EISParams
"""
Base.@kwdef struct EISParams <: AbstractCurrentParams
    i_EIS::Float64 = 1.0e4  # (A/m²). Parameters for the EIS curve.
    ratio::Float64 = 5.0 / 100.0  # (.). Parameters for the EIS curve.
    f_power_min::Float64 = -3.0  # (.). Power of the minimum frequency for the EIS current density function.
    f_power_max::Float64 = 5.0 # (.). Power of the maximum frequency for the EIS current density function.
    nb_f::Int = 90 # (.). Number of frequencies tested for the EIS current density function.
    nb_points::Int = 50 # (.). Number of points calculated per specific period for the EIS current density function.
end