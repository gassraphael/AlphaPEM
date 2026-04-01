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

# __________________________________________________Run simulation___________________________________________________

function run_simulation(cfg::SimulationConfig)::AlphaPEM

    # Determine the number of simulators to create based on the length of type_fuel_cells.
    nb_fuel_cells = length(cfg.type_fuel_cells)

    # Build Fuelcell objects for each configuration using the factory
    fuel_cells = [create_fuelcell(cfg.type_fuel_cells[i], cfg.voltage_zone) for i in 1:nb_fuel_cells]

    # Build Current objects for each configuration using the factory
    current_densities = [create_current(cfg.type_current, fuel_cells[i]) for i in 1:nb_fuel_cells]

    # Create one simulator per selected fuel cell and current density configuration.
    AlphaPEM_simulators = [AlphaPEM(fuel_cells[i], current_densities[i], cfg) for i in 1:nb_fuel_cells]

    # Launch the simulation.
    if type_current == "step"
        launch_AlphaPEM_for_step_current(AlphaPEM_simulators[1])
    elseif type_current == "polarization"
        launch_AlphaPEM_for_polarization_current(AlphaPEM_simulators)
    elseif type_current == "polarization_for_cali"
        launch_AlphaPEM_for_polarization_current_for_calibration(AlphaPEM_simulators)
    elseif type_current == "EIS"
        launch_AlphaPEM_for_EIS_current(AlphaPEM_simulators[1])
    else
        throw(ArgumentError("You have to specify a type_current which is accepted."))
    end

    # Disable interactive mode for non-blocking display.
    plt.ioff()
    # Ensure that the figures remain displayed after program execution.
    plt.show(; block=true)

    return Simulator
end
