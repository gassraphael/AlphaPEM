# -*- coding: utf-8 -*-

"""
    FuelCell factory for AlphaPEM

This module provides a factory function to instantiate the correct FuelCell type
(ZSWFuelCell, EH31FuelCell, or DefaultFuelCell) based on the type_fuel_cell string.
"""


const FUELCELL_TYPE_MAP = Dict(
    :ZSW_GenStack => ZSWFuelCell,
    :ZSW_GenStack_Pa_1_61_Pc_1_41 => ZSWFuelCell,
    :ZSW_GenStack_Pa_2_01_Pc_1_81 => ZSWFuelCell,
    :ZSW_GenStack_Pa_2_4_Pc_2_2 => ZSWFuelCell,
    :ZSW_GenStack_Pa_2_8_Pc_2_6 => ZSWFuelCell,
    :ZSW_GenStack_T_62 => ZSWFuelCell,
    :ZSW_GenStack_T_76 => ZSWFuelCell,
    :ZSW_GenStack_T_84 => ZSWFuelCell,
    :EH_31_1_5 => EH31FuelCell,
    :EH_31_2_0 => EH31FuelCell,
    :EH_31_2_25 => EH31FuelCell,
    :EH_31_2_5 => EH31FuelCell
)

function create_fuelcell(type_fuel_cell::Symbol, voltage_zone::Symbol)::AbstractFuelCell
    simulator = get(FUELCELL_TYPE_MAP, type_fuel_cell, DefaultFuelCell)
    if simulator === DefaultFuelCell
        return DefaultFuelCell()
    else
        return simulator(type_fuel_cell, voltage_zone)
    end
end
