# -*- coding: utf-8 -*-

"""This file creates the current density objects."""

# _______________________________________________Current density factory_____________________________________________

"""
    create_current(p::AbstractCurrentParams; fuel_cell=nothing)::AbstractCurrent

Create a current profile object according to the type of p.
If a `fuel_cell` object is provided, experimental values will be extracted from it for polarization and calibration profiles.
Otherwise, default parameter values are used.
"""

function create_current(p::StepParams, fuel_cell=nothing)::AbstractCurrent
    return StepCurrent(p)
end


function create_current(p::PolarizationParams, fuel_cell=nothing)::AbstractCurrent
    # If a FuelCell object is provided and contains experimental current data,
    # use the maximum value as i_max for the polarization profile.
    if fuel_cell !== nothing && hasproperty(fuel_cell, :pola_exp_data) &&
           hasproperty(fuel_cell.pola_exp_data, :i_exp) && !isempty(fuel_cell.pola_exp_data.i_exp)
            i_max_val = maximum(fuel_cell.pola_exp_data.i_exp)
            p = PolarizationParams(
                delta_t_ini = p.delta_t_ini,
                delta_i = p.delta_i,
                v_load = p.v_load,
                delta_t_break = p.delta_t_break,
                i_max = i_max_val,
            )
    end
    return PolarizationCurrent(p)
end


function create_current(p::PolarizationCalibrationParams, fuel_cell=nothing)::AbstractCurrent
    # If a FuelCell object is provided and contains calibration experimental data,
    # pass i_exp as a field in the parameter struct (requires constructor support).
    if fuel_cell !== nothing && hasproperty(fuel_cell, :pola_exp_data_cali) &&
       hasproperty(fuel_cell.pola_exp_data_cali, :i_exp) && !isempty(fuel_cell.pola_exp_data_cali.i_exp)
        p = PolarizationCalibrationParams(
            delta_t_ini = p.delta_t_ini,
            v_load = p.v_load,
            delta_t_break = p.delta_t_break,
            i_exp = Float64.(fuel_cell.pola_exp_data_cali.i_exp),
        )
    end
    return PolarizationCalibrationCurrent(p)
end



function create_current(p::EISParams, fuel_cell=nothing)::AbstractCurrent
    return EISCurrent(p)
end

