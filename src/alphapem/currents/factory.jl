# -*- coding: utf-8 -*-

"""This file creates the current density objects."""

# _______________________________________________Current density factory_____________________________________________

"""
    create_current(p::AbstractCurrentParams; fuelcell::AbstractFuelCell=nothing)::AbstractCurrent

Create a current profile object according to the type of p.
If a `fuelcell` object is provided, experimental values will be extracted from it for polarization and calibration profiles.
Otherwise, default parameter values are used.
"""

function create_current(p::StepParams)::AbstractCurrent
    return StepCurrent(p)
end


function create_current(p::PolarizationParams, fuelcell::AbstractFuelCell=nothing)::AbstractCurrent
    # If a FuelCell object is provided and contains experimental current data,
    # use the maximum value as i_max for the polarization profile.
    if fuelcell !== nothing && hasproperty(fuelcell, :pola_exp_data) &&
           hasproperty(fuelcell.pola_exp_data, :i_exp) && !isempty(fuelcell.pola_exp_data.i_exp)
            i_max_val = maximum(fuelcell.pola_exp_data.i_exp)
            p = PolarizationParams(; NamedTuple(p)..., i_max=i_max_val)
    end
    return PolarizationCurrent(p)
end


function create_current(p::PolarizationCalibrationParams, fuelcell::AbstractFuelCell=nothing)::AbstractCurrent
    # If a FuelCell object is provided and contains calibration experimental data,
    # pass i_exp as a field in the parameter struct (requires constructor support).
    if fuelcell !== nothing && hasproperty(fuelcell, :pola_exp_data_cali) &&
       hasproperty(fuelcell.pola_exp_data_cali, :i_exp) && !isempty(fuelcell.pola_exp_data_cali.i_exp)
        p = PolarizationCalibrationParams(; NamedTuple(p)..., i_exp=fuelcell.pola_exp_data_cali.i_exp)
    end
    return PolarizationCalibrationCurrent(p)
end


function create_current(p::EISParams)::AbstractCurrent
    return EISCurrent(p)
end
