# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using PyCall

# Keep matplotlib usage via Python interop
const mpl = pyimport("matplotlib")
const plt = pyimport("matplotlib.pyplot")

# Importing constants' value and functions
include(joinpath(@__DIR__, "../fuelcell/Fuelcell.jl"))
include(joinpath(@__DIR__, "../currents/Currents.jl"))
include(joinpath(@__DIR__, "../core/models/AlphaPEM.jl"))
include(joinpath(@__DIR__, "run_simulation_modules.jl"))

# __________________________________________________AlphaPEM settings___________________________________________________
"""
- Users can select various preconfigured configurations for execution.
- Adjustments to these configurations can be made within settings.jl or current_densities.jl and their associated files.
- Selecting different type_fuel_cell during a single run results in simultaneous plots for various configurations.
"""

function run_simulation()::AlphaPEM
    # Fuel cell possibilities: "ZSW-GenStack"(2022), "ZSW-GenStack_Pa_1.61_Pc_1.41", "ZSW-GenStack_Pa_2.01_Pc_1.81",
    #                          "ZSW-GenStack_Pa_2.4_Pc_2.2", "ZSW-GenStack_Pa_2.8_Pc_2.6", "ZSW-GenStack_T_62",
    #                          "ZSW-GenStack_T_76", "ZSW-GenStack_T_84", "EH-31_1.5"(2021), "EH-31_2.0"(2021),
    #                          "EH-31_2.25"(2021), "EH-31_2.5"(2021), "manual_setup".
    # This parameter includes the fuel cell used in the model and the corresponding operating conditions.
    # - GenStack is a fuel cell developed in open source by ZSW (https://zenodo.org/records/14223364).
    # - EH-31 is a fuel cell developed by EH GROUP. 1.5, 2.0, 2.25 and 2.5 corresponds to the different pressure options.
    type_fuel_cells::Vector{String} = ["ZSW-GenStack"]

    # Current density possibilities: "step", "polarization", "polarization_for_cali", "EIS".
    type_current::String = "polarization"
    # Calibration zone: "before_voltage_drop", "full".
    # (only for "polarization" and "polarization_for_cali" current densities)
    voltage_zone::String = "full"
    # Auxiliary system possibilities: "forced-convective_cathode_with_anodic_recirculation",
    #                                 "forced-convective_cathode_with_flow-through_anode", "no_auxiliary".
    type_auxiliary::String = "no_auxiliary"
    # Purges possibilities: "constant_purge", "periodic_purge", "no_purge".
    type_purge::String = "no_purge"
    # Display possibilities: "multiple", "synthetic", "no_display".
    type_display::String = "synthetic"
    # Plot possibilities: "dynamic", "fixed". Using dynamic plot option enables real-time figure updates during
    # execution, albeit at the cost of decreased program speed.
    type_plot::String = "fixed"

    # __________________________________________________________Main____________________________________________________
    """This section is dedicated to ensuring the proper execution of the simulator, considering all the various
    possibilities including real-time figure updates and simultaneous plotting for different configurations.
    This should remain unaltered for regular program usage.
    """

    # Warning for potential issues with the selected configurations.
    if isempty(type_fuel_cells)
        throw(ArgumentError("type_fuel_cells must contain at least one fuel cell configuration."))
    end

    # Build Fuelcell objects for each configuration using the factory
    simulators = [create_fuelcell(type_fuel_cells[i], voltage_zone) for i in 1:length(type_fuel_cells)]

    # Build Current objects for each configuration using the factory
    current_densities = [create_current(type_current, simulators[i]) for i in 1:length(simulators)]

    # Create one simulator per selected fuel cell and current density configuration.
    simulators = [AlphaPEM(simulators[i], current_densities[i]) for i in 1:length(simulators)]

    # Check if type_current is valid and launch the simulation.
    Simulator::AlphaPEM = begin
        if type_current == "step"
            launch_AlphaPEM_for_step_current(simulators[1],
                                             select_nth_elements(operating_inputs, 1),
                                             select_nth_elements(current_parameters, 1),
                                             select_nth_elements(computing_parameters, 1))
        elseif type_current == "polarization"
            launch_AlphaPEM_for_polarization_current(simulators, operating_inputs, current_parameters,
                                                     computing_parameters)
        elseif type_current == "polarization_for_cali"
            launch_AlphaPEM_for_polarization_current_for_calibration(simulators, operating_inputs,
                                                                     current_parameters,
                                                                     computing_parameters)
        elseif type_current == "EIS"
            launch_AlphaPEM_for_EIS_current(simulators[1],
                                            select_nth_elements(operating_inputs, 1),
                                            select_nth_elements(current_parameters, 1),
                                            select_nth_elements(computing_parameters, 1))
        else
            throw(ArgumentError("You have to specify a type_current which is accepted."))
        end
    end

    # Disable interactive mode for non-blocking display.
    plt.ioff()
    # Ensure that the figures remain displayed after program execution.
    plt.show(; block=true)

    return Simulator
end

# Entry point: run_simulation() is called directly, as is standard practice for Julia scripts.
# This file is intended to be executed as a standalone script, not included as a module.
run_simulation()
