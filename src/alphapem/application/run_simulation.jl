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
include(joinpath(@__DIR__, "../core/models/AlphaPEM.jl"))
include(joinpath(@__DIR__, "../config/parameters.jl"))
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

    # Retrieving parameters from settings.
    #   Imposed inputs
    (step_current_parameters, pola_current_parameters, pola_current_for_cali_parameters, i_EIS, ratio_EIS, f_EIS, t_EIS,
     current_density) = calculate_current_density_parameters(type_current)

    if isempty(type_fuel_cells)
        throw(ArgumentError("type_fuel_cells must contain at least one fuel cell configuration."))
    end

    #   Operating conditions (one set per selected fuel cell configuration)
    operating_conditions = [
        calculate_operating_inputs(deepcopy(pola_current_parameters), type_fuel_cell, voltage_zone)
        for type_fuel_cell in type_fuel_cells
    ]
    T_des = [oc[1] for oc in operating_conditions]
    Pa_des = [oc[2] for oc in operating_conditions]
    Pc_des = [oc[3] for oc in operating_conditions]
    Sa = [oc[4] for oc in operating_conditions]
    Sc = [oc[5] for oc in operating_conditions]
    Phi_a_des = [oc[6] for oc in operating_conditions]
    Phi_c_des = [oc[7] for oc in operating_conditions]
    y_H2_in = [oc[8] for oc in operating_conditions]
    pola_current_parameters_list = [oc[9] for oc in operating_conditions]

    #   Physical parameters
    (Hacl, Hccl, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc, Lgc, nb_channel_in_gc,
     Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem, Aact, nb_cell, e, K_l_ads, K_O2_ad_Pt, Re, i0_c_ref, kappa_co,
    kappa_c, C_scl) = calculate_physical_parameters(type_fuel_cells[1])

    #   Computing parameters
    nb_gc, nb_gdl, nb_mpl, purge_time, delta_purge, rtol, atol = calculate_computing_parameters(step_current_parameters)

    # Initialize the operating inputs and parameters dictionaries.
    # In Julia, vectors are naturally 1-based, so no [None] placeholder is needed.
    operating_inputs = Dict(
        "current_density" => current_density,
        "T_des" => T_des,
        "Pa_des" => Pa_des,
        "Pc_des" => Pc_des,
        "Sa" => Sa,
        "Sc" => Sc,
        "Phi_a_des" => Phi_a_des,
        "Phi_c_des" => Phi_c_des,
        "y_H2_in" => y_H2_in,
    )

    current_parameters = Dict(
        "step_current_parameters" => step_current_parameters,
        "pola_current_parameters" => pola_current_parameters_list,
        "pola_current_for_cali_parameters" => pola_current_for_cali_parameters,
        "i_EIS" => i_EIS,
        "ratio_EIS" => ratio_EIS,
        "t_EIS" => t_EIS,
        "f_EIS" => f_EIS,
    )

    accessible_physical_parameters = Dict(
        "Aact" => Aact, "nb_cell" => nb_cell, "Hagc" => Hagc, "Hcgc" => Hcgc, "Wagc" => Wagc,
        "Wcgc" => Wcgc, "Lgc" => Lgc, "nb_channel_in_gc" => nb_channel_in_gc, "Ldist" => Ldist,
        "Lm" => Lm, "A_T_a" => A_T_a, "A_T_c" => A_T_c, "Vasm" => Vasm, "Vcsm" => Vcsm,
        "Vaem" => Vaem, "Vcem" => Vcem,
    )

    undetermined_physical_parameters = Dict(
        "Hgdl" => Hgdl, "Hmpl" => Hmpl, "Hmem" => Hmem, "Hacl" => Hacl,
        "Hccl" => Hccl, "epsilon_gdl" => epsilon_gdl, "epsilon_mpl" => epsilon_mpl,
        "epsilon_c" => epsilon_c, "e" => e, "K_l_ads" => K_l_ads, "K_O2_ad_Pt" => K_O2_ad_Pt,
        "Re" => Re, "i0_c_ref" => i0_c_ref, "kappa_co" => kappa_co, "kappa_c" => kappa_c,
        "C_scl" => C_scl,
    )

    model_parameters = Dict(
        "nb_gc" => nb_gc, "nb_gdl" => nb_gdl, "nb_mpl" => nb_mpl, "purge_time" => purge_time,
        "delta_purge" => delta_purge, "rtol" => rtol, "atol" => atol,
    )

    computing_parameters = Dict(
        "type_fuel_cell" => type_fuel_cells,
        "type_current" => type_current,
        "voltage_zone" => voltage_zone,
        "type_auxiliary" => type_auxiliary,
        "type_purge" => type_purge,
        "type_display" => type_display,
        "type_plot" => type_plot,
    )

    # Create one simulator per selected fuel cell configuration.
    simulators = [AlphaPEM(accessible_physical_parameters, undetermined_physical_parameters, model_parameters)
                  for _ in 1:length(computing_parameters["type_fuel_cell"])]

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
