# -*- coding: utf-8 -*-

"""
    FuelCell factory for AlphaPEM

This module provides a factory function to instantiate the correct FuelCell type
(ZSWFuelCell, EH31FuelCell, or DefaultFuelCell) based on the type_fuel_cell string.
"""

include("abstract.jl")
include("zsw.jl")
include("eh31.jl")
include("default.jl")

const FUELCELL_TYPE_MAP = Dict(
    "ZSW-GenStack" => ZSWFuelCell,
    "ZSW-GenStack_Pa_1.61_Pc_1.41" => ZSWFuelCell,
    "ZSW-GenStack_Pa_2.01_Pc_1.81" => ZSWFuelCell,
    "ZSW-GenStack_Pa_2.4_Pc_2.2" => ZSWFuelCell,
    "ZSW-GenStack_Pa_2.8_Pc_2.6" => ZSWFuelCell,
    "ZSW-GenStack_T_62" => ZSWFuelCell,
    "ZSW-GenStack_T_76" => ZSWFuelCell,
    "ZSW-GenStack_T_84" => ZSWFuelCell,
    "EH-31_1.5" => EH31FuelCell,
    "EH-31_2.0" => EH31FuelCell,
    "EH-31_2.25" => EH31FuelCell,
    "EH-31_2.5" => EH31FuelCell
)

function create_fuelcell(type_fuel_cell::String, voltage_zone::String)::AbstractFuelCell
    simulator = get(FUELCELL_TYPE_MAP, type_fuel_cell, DefaultFuelCell)
    # DefaultFuelCell accepte seulement voltage_zone, les autres prennent (type_fuel_cell, voltage_zone)
    if simulator === DefaultFuelCell
        return DefaultFuelCell()
    else
        return simulator(type_fuel_cell, voltage_zone)
    end
end


