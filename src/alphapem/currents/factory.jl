# -*- coding: utf-8 -*-

"""This file creates the current density objects."""

# _______________________________________________Current density factory_____________________________________________

"""
    create_current(type_current::String; fuelcell=nothing)::AbstractCurrent

Create a current profile object according to the type string.
If a `fuelcell` object is provided, experimental values will be extracted from it for polarization and calibration profiles.
Otherwise, default parameter values are used.

Arguments:
    type_current::String: "step", "polarization", "polarization_for_cali", or "EIS"
    fuelcell: (optional) FuelCell object containing experimental current data

Returns:
    The corresponding current profile object.
"""
function create_current(type_current::String, fuelcell::AbstractFuelCell=nothing)::AbstractCurrent
    if type_current == "step"
        params = StepParams()
    elseif type_current == "polarization"
        params = PolarizationParams()
        # If a FuelCell object is provided and contains experimental current data,
        # use the maximum value as i_max for the polarization profile.
        if fuelcell !== nothing && hasproperty(fuelcell, :pola_exp_data) &&
           hasproperty(fuelcell.pola_exp_data, :i_exp) && !isempty(fuelcell.pola_exp_data.i_exp)
            i_max_val = maximum(fuelcell.pola_exp_data.i_exp)
            params = PolarizationParams(; NamedTuple(params)..., i_max=i_max_val)
        end
    elseif type_current == "polarization_for_cali"
        params = PolarizationCalibrationParams()
        # If a FuelCell object is provided and contains calibration experimental data,
        # pass i_exp as a field in the parameter struct (requires constructor support).
        if fuelcell !== nothing && hasproperty(fuelcell, :pola_exp_data_cali) &&
           hasproperty(fuelcell.pola_exp_data_cali, :i_exp) && !isempty(fuelcell.pola_exp_data_cali.i_exp)
            params = PolarizationCalibrationParams(; NamedTuple(params)..., i_exp=fuelcell.pola_exp_data_cali.i_exp)
        end
    elseif type_current == "EIS"
        params = EISParams()
    else
        throw(ArgumentError("Unknown type_current: $type_current"))
    end
    return create_current(params)
end


function create_current(p::StepParams)::AbstractCurrent
    return StepCurrent(
        p.delta_t_ini,
        p.delta_t_load,
        p.delta_t_break,
        p.i_ini,
        p.i_step
    )
end


function create_current(p::PolarizationParams)::AbstractCurrent
    return PolarizationCurrent(
        p.delta_t_ini,
        p.v_load,
        p.delta_t_break,
        p.delta_i,
        p.i_max_val
    )
end


function create_current(p::PolarizationCalibrationParams)::AbstractCurrent
    return PolarizationCalibrationCurrent(
        p.delta_t_ini,
        p.v_load,
        p.delta_t_break,
        p.i_exp
    )
end


function create_current(p::EISParams)::AbstractCurrent
    return EISCurrent(
        p.i_EIS,
        p.ratio,
        p.f_power_min,
        p.f_power_max,
        p.nb_f,
        p.nb_points
    )
end
