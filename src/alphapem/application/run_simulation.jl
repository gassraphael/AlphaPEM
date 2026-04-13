# -*- coding: utf-8 -*-

"""This file is designated for executing the AlphaPEM software package.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Shared package modules
using ..Config: SimulationConfig, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using ..Fuelcell: create_fuelcell
using ..Currents: create_current, step_duration
using ..Core.Models: AlphaPEM, simulate_model!, Display, Save_plot


include(joinpath(@__DIR__, "run_simulation_modules.jl"))

# __________________________________________________Run simulation___________________________________________________

function run_simulation(cfg::SimulationConfig)::AlphaPEM
    # Build a Fuelcell object with the configuration given using the factory
    fuel_cells = create_fuelcell(cfg.type_fuel_cell, cfg.voltage_zone)

    # Build a Current object with the configuration given using the factory
    current_density = create_current(cfg.type_current, fuel_cells)

    # Create a simulator
    AlphaPEM_simulator = AlphaPEM(fuel_cells, current_density, cfg)

    # Launch the simulation.
    if cfg.type_current isa StepParams
        launch_AlphaPEM_for_step_current(AlphaPEM_simulator)
    elseif cfg.type_current isa PolarizationParams
        launch_AlphaPEM_for_polarization_current(AlphaPEM_simulator)
    elseif cfg.type_current isa PolarizationCalibrationParams
        launch_AlphaPEM_for_polarization_current_for_calibration(AlphaPEM_simulator)
    elseif cfg.type_current isa EISParams
        launch_AlphaPEM_for_EIS_current(AlphaPEM_simulator)
    else
        throw(ArgumentError("You have to specify a type_current which is accepted."))
    end

    # Disable interactive mode for non-blocking display.
    if cfg.type_display != :no_display
        plt.ioff()
        # Ensure that the figures remain displayed after program execution.
        plt.show(; block=true)
    end

    return AlphaPEM_simulator
end

function run_simulation(cfgs::AbstractVector{<:SimulationConfig})::Vector{AlphaPEM}

    # Determine the number of simulators to create based on the length of type_fuel_cells.
    nb_fuel_cells = length(cfgs)

    # Build Fuelcell objects for each configuration using the factory
    fuel_cells = [create_fuelcell(cfgs[i].type_fuel_cell, cfgs[i].voltage_zone) for i in 1:nb_fuel_cells]

    # Build Current objects for each configuration using the factory
    current_densities = [create_current(cfgs[i].type_current, fuel_cells[i]) for i in 1:nb_fuel_cells]

    # Create one simulator per selected fuel cell and current density configuration.
    AlphaPEM_simulators = [AlphaPEM(fuel_cells[i], current_densities[i], cfgs[i]) for i in 1:nb_fuel_cells]

    # Launch the simulation.
    if all(cfg -> cfg.type_current isa PolarizationParams, cfgs)
        launch_AlphaPEM_for_polarization_current(AlphaPEM_simulators)
    elseif all(cfg -> cfg.type_current isa PolarizationCalibrationParams, cfgs)
        launch_AlphaPEM_for_polarization_current_for_calibration(AlphaPEM_simulators)
    else
        throw(ArgumentError("You have to specify homogeneous type_current parameters accepted for multi-simulator runs."))
    end

    # Disable interactive mode for non-blocking display.
    if any(cfg -> cfg.type_display != :no_display, cfgs)
        plt.ioff()
        # Ensure that the figures remain displayed after program execution.
        plt.show(; block=true)
    end

    return AlphaPEM_simulators
end

