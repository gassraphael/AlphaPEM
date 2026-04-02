# -*- coding: utf-8 -*-

"""
    DefaultFuelCell

A generic fallback FuelCell using default parameters, for unknown or manual types.
"""


struct DefaultFuelCell <: AbstractFuelCell
    physical_parameters::PhysicalParams
    operating_conditions::OperatingConditions
    pola_exp_data::PolaExperimentalData
    pola_exp_data_cali::PolaExperimentalData
    numerical_parameters::NumericalParams
end

function DefaultFuelCell()
    fc = DefaultFuelCell(
        PhysicalParams(),
        OperatingConditions(),
        PolaExperimentalData(),
        PolaExperimentalData(),
        NumericalParams()
    )
    return fc
end


