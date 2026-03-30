# -*- coding: utf-8 -*-

"""This file creates the current density objects."""

# _______________________________________________Current density factory_____________________________________________


function make_current(p::StepParams)
    return StepCurrent(
        p.delta_t_ini,
        p.delta_t_load,
        p.delta_t_break,
        p.i_ini,
        p.i_step
    )
end


function make_current(p::PolarizationParams)
    # If both type_fuel_cell and voltage_zone are provided, try to get experimental current values
    local i_max_val
    if p.type_fuel_cell !== nothing && p.voltage_zone !== nothing
        i_exp, _ = pola_exp_values_calibration(p.type_fuel_cell, p.voltage_zone)
        if length(i_exp) > 0
            i_max_val = maximum(i_exp)
        else
            i_max_val = p.i_max
        end
    else
        i_max_val = p.i_max
    end
    return PolarizationCurrent(
        p.delta_t_ini,
        p.v_load,
        p.delta_t_break,
        p.delta_i,
        i_max_val
    )
end


function make_current(p::PolarizationCalibrationParams)
    # type_fuel_cell and voltage_zone are mandatory
    i_exp, _ = pola_exp_values_calibration(
        p.type_fuel_cell,
        p.voltage_zone
    )
    return PolarizationCalibrationCurrent(
        p.delta_t_ini,
        p.v_load,
        p.delta_t_break,
        i_exp
    )
end


function make_current(p::EISParams)
    return EISCurrent(
        p.i_EIS,
        p.ratio,
        p.f_power_min,
        p.f_power_max,
        p.nb_f,
        p.nb_points
    )
end