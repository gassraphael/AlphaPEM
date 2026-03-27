# -*- coding: utf-8 -*-

"""This module contains some of the functions required for the parameter calibration.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using PyCall
using Dates
using Statistics
using LinearAlgebra
using Interpolations
using Printf

# Keep colorama for colored output (optional, can be replaced with native Julia)
const Fore = pyimport("colorama").Fore
const Style = pyimport("colorama").Style

# Importing functions
include(joinpath(@__DIR__, "../config/current_densities.jl"))
include(joinpath(@__DIR__, "../config/pola_exp_values.jl"))

# Python interop for numpy operations
const np = pyimport("numpy")

# _________________________________________________Calibration modules__________________________________________________

"""Determines the parameter bounds of the fuel cell model for calibration.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `voltage_zone::String`: Zone of calibration: "before_voltage_drop", "full".
- `operating_inputs_1::Dict`: Operating inputs for the first fuel cell configuration.
- `operating_inputs_2::Dict`: Operating inputs for the second fuel cell configuration.

# Returns
- `varbound::Vector{Vector}`: List of the bounds on the parameters to calibrate. Each element is a vector
  containing the name, minimum value, maximum value, and type of the parameter.
- `gene_space::Vector{Dict}`: List of dictionaries used to define the bounds of the undetermined parameters for
  pygad. Each dictionary contains the 'low' and 'high' values for the parameter, and optionally a 'step' value
  for integer parameters.
"""
function parameter_bounds_for_calibration(type_fuel_cell::String,
                                          voltage_zone::String,
                                          operating_inputs_1::Dict,
                                          operating_inputs_2::Dict)

    Pc_des_1, Pc_des_2 = operating_inputs_1["Pc_des"], operating_inputs_2["Pc_des"]

    if type_fuel_cell in ("ZSW-GenStack", "ZSW-GenStack_Pa_1.61_Pc_1.41", "ZSW-GenStack_Pa_2.01_Pc_1.81",
                          "ZSW-GenStack_Pa_2.4_Pc_2.2", "ZSW-GenStack_Pa_2.8_Pc_2.6", "ZSW-GenStack_T_62",
                          "ZSW-GenStack_T_76", "ZSW-GenStack_T_84")
        # Fuel cell physical parameters
        Hacl_min, Hacl_max = 5e-6, 15e-6  # m. It is the thickness of the ACL.
        Hccl_min, Hccl_max = 5e-6, 20e-6  # m. It is the thickness of the CCL.
        Hmem_min, Hmem_max = 5e-6, 30e-6  # m. It is the thickness of the membrane.
        Hgdl_min, Hgdl_max = 100e-6, 150e-6  # m. It is the thickness of the gas diffusion layer.
        Hmpl_min, Hmpl_max = 40e-6, 100e-6  # m. It is the thickness of the microporous layer.
        epsilon_gdl_min, epsilon_gdl_max = 0.5, 0.9  # It is the anode/cathode GDL porosity, without units.
        # Constants based on the interaction between fluids and the structure
        e_min, e_max = 3, 5  # It is the capillary exponent, and should be an int number.
        K_l_ads_min, K_l_ads_max = 1, 100  # . It is the ratio between the liquid and vapor sorption rates of water in the membrane.
        K_O2_ad_Pt_min, K_O2_ad_Pt_max = 1, 10  # . It is the interfacial resistance coefficient of O2 adsorption on the Pt sites.
        # Voltage polarization
        Re_min, Re_max = 5e-8, 5e-6  # Ω.m². It is the electron conduction resistance of the circuit.
        i0_c_ref_min, i0_c_ref_max = 1e-1, 100  # A.m-2. It is the dry reference exchange current density at the cathode.
        kappa_co_min, kappa_co_max = 0.01, 40  # A.m-2. It is the crossover correction coefficient.
        kappa_c_min, kappa_c_max = 0.25, 4  # It is the overpotential correction exponent.
        # Undetermined parameter which is not considered yet (require the use of EIS curves to be calibrated)
        C_scl_min, C_sl_max = 2e7, 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.

        # Bounds gathering and type
        varbound = [
            ["Hacl", Hacl_min, Hacl_max, "real"],
            ["Hccl", Hccl_min, Hccl_max, "real"],
            ["Hmem", Hmem_min, Hmem_max, "real"],
            ["Hgdl", Hgdl_min, Hgdl_max, "real"],
            ["Hmpl", Hmpl_min, Hmpl_max, "real"],
            ["epsilon_gdl", epsilon_gdl_min, epsilon_gdl_max, "real"],
            ["e", e_min, e_max, "int"],
            ["Re", Re_min, Re_max, "real"],
            ["i0_d_c_ref", i0_c_ref_min, i0_c_ref_max, "real"],
            ["kappa_co", kappa_co_min, kappa_co_max, "real"],
            ["kappa_c", kappa_c_min, kappa_c_max, "real"],
        ]

        if voltage_zone == "full"
            push!(varbound, ["K_l_ads", K_l_ads_min, K_l_ads_max, "real"])
            push!(varbound, ["K_O2_ad_Pt", K_O2_ad_Pt_min, K_O2_ad_Pt_max, "real"])
        end

        gene_space = []  # List used to define the bounds of the undetermined parameters for pygad.
        for i in eachindex(varbound)
            name, min_val, max_val, type_val = varbound[i]
            if type_val == "int"
                push!(gene_space, Dict("low" => min_val, "high" => max_val, "step" => 1))
            else
                push!(gene_space, Dict("low" => min_val, "high" => max_val))
            end
        end

    elseif type_fuel_cell in ("EH-31_1.5", "EH-31_2.0", "EH-31_2.25", "EH-31_2.5")
        # Fuel cell physical parameters
        Hacl_min, Hacl_max = 8e-6, 20e-6  # m. It is the thickness of the ACL.
        Hmem_min, Hmem_max = 15e-6, 50e-6  # m. It is the thickness of the membrane.
        epsilon_gdl_min, epsilon_gdl_max = 0.40, 0.95  # It is the anode/cathode GDL porosity, without units.
        epsilon_c_min, epsilon_c_max = 0.15, 0.30  # It is the compression ratio of the GDL.
        # Constants based on the interaction between fluids and the structure
        e_min, e_max = 3, 5  # It is the capillary exponent, and should be an int number.
        K_O2_ad_Pt_min, K_O2_ad_Pt_max = 1, 10  # . It is the interfacial resistance coefficient of O2 adsorption on the Pt sites.
        # Voltage polarization
        Re_min, Re_max = 5e-7, 5e-6  # Ω.m². It is the electron conduction resistance of the circuit.
        i0_c_ref_min, i0_c_ref_max = 1e-1, 100  # A.m-2. It is the dry reference exchange current density at the cathode.
        kappa_co_min, kappa_co_max = 0.01, 40  # A.m-2. It is the crossover correction coefficient.
        kappa_c_min, kappa_c_max = 0.25, 4  # It is the overpotential correction exponent.
        # Undetermined parameter which is not considered yet (require the use of EIS curves to be calibrated)
        C_scl_min, C_sl_max = 2e7, 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.

        # Bounds gathering and type
        varbound = [
            ["Hacl", Hacl_min, Hacl_max, "real"],
            ["Hmem", Hmem_min, Hmem_max, "real"],
            ["epsilon_gdl", epsilon_gdl_min, epsilon_gdl_max, "real"],
            ["e", e_min, e_max, "int"],
            ["Re", Re_min, Re_max, "real"],
            ["i0_d_c_ref", i0_c_ref_min, i0_c_ref_max, "real"],
            ["kappa_co", kappa_co_min, kappa_co_max, "real"],
            ["kappa_c", kappa_c_min, kappa_c_max, "real"],
        ]

        if voltage_zone == "full"
            push!(varbound, ["K_O2_ad_Pt", K_O2_ad_Pt_min, K_O2_ad_Pt_max, "real"])
        end

        gene_space = []  # List used to define the bounds of the undetermined parameters for pygad.
        for i in eachindex(varbound)
            name, min_val, max_val, type_val = varbound[i]
            if type_val == "int"
                push!(gene_space, Dict("low" => min_val, "high" => max_val, "step" => 1))
            else
                push!(gene_space, Dict("low" => min_val, "high" => max_val))
            end
        end
    else
        throw(ArgumentError("A correct type_fuel_cell should be given."))
    end

    return varbound, gene_space
end


"""Determines the parameters of the fuel cell model for calibration.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `voltage_zone::String`: Zone of calibration: "before_voltage_drop", "full".

# Returns
- `operating_inputs::Dict`: Dictionary containing operating inputs (current_density, T_des, Pa_des, etc.).
- `current_parameters::Dict`: Dictionary containing current parameters.
- `accessible_physical_parameters::Dict`: Dictionary of accessible physical parameters.
- `undetermined_physical_parameters::Dict`: Dictionary of undetermined physical parameters.
- `model_parameters::Dict`: Dictionary of model parameters.
- `computing_parameters::Dict`: Dictionary of computing parameters.
- `i_exp::Vector{Float64}`: Experimental values of the current density.
- `U_exp::Vector{Float64}`: Experimental values of the voltage.
"""
function parameters_for_calibration(type_fuel_cell::String, voltage_zone::String)

    # Algorithm parameters for polarization curve generation
    type_auxiliary = "no_auxiliary"
    type_purge = "no_purge"
    type_display = "no_display"
    type_plot = "fixed"
    type_current = "polarization_for_cali"
    current_density = polarization_current_for_calibration
    delta_t_ini_step = 30 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states (standard value).
    delta_t_load_step = 1e-15  # (s). Loading time for the step current density function, from 0 to i_step.
    delta_t_break_step = 0  # (s). Time at i_step current density for the stabilisation of the internal states.
    i_step = 1.0e4  # (A.m-2). Current density for the step current density function.
    step_current_parameters = Dict("delta_t_ini_step" => delta_t_ini_step, "delta_t_load_step" => delta_t_load_step,
                                   "delta_t_break_step" => delta_t_break_step, "i_step" => i_step)
    delta_t_ini_pola = 30 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_pola = 30  # (s). Loading time for one step current of the polarisation current density function.
    delta_t_break_pola = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
    delta_i_pola = 0.05e4  # (A.m-2). Current density step for the polarisation current density function.
    pola_current_parameters = Dict("delta_t_ini_pola" => delta_t_ini_pola, "delta_t_load_pola" => delta_t_load_pola,
                                   "delta_t_break_pola" => delta_t_break_pola, "delta_i_pola" => delta_i_pola)
    delta_t_ini_pola_cali = 30 * 60  # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_t_load_pola_cali = 30  # (s). Loading time for one step current of the polarisation current density function.
    delta_t_break_pola_cali = 15 * 60  # (s). Breaking time for one step current, for the stabilisation of the internal states.
    pola_current_for_cali_parameters = Dict("delta_t_ini_pola_cali" => delta_t_ini_pola_cali,
                                            "delta_t_load_pola_cali" => delta_t_load_pola_cali,
                                            "delta_t_break_pola_cali" => delta_t_break_pola_cali)
    i_EIS, ratio_EIS = NaN, NaN  # (A/m², ). i_EIS is the current for which a ratio_EIS perturbation is added.
    f_EIS, t_EIS = NaN, NaN  # It is the EIS parameters.
    t_purge = (0.6, 15)  # s. It is the purge time and the distance between two purges.
    rtol = 1e-6  # Relative tolerance for the system of ODEs solver.
    atol = 1e-9  # Absolute tolerance for the system of ODEs solver.

    if type_fuel_cell in ("ZSW-GenStack", "ZSW-GenStack_Pa_1.61_Pc_1.41", "ZSW-GenStack_Pa_2.01_Pc_1.81",
                          "ZSW-GenStack_Pa_2.4_Pc_2.2", "ZSW-GenStack_Pa_2.8_Pc_2.6", "ZSW-GenStack_T_62",
                          "ZSW-GenStack_T_76", "ZSW-GenStack_T_84")
        # Given values by the author
        # Operating inputs
        if type_fuel_cell == "ZSW-GenStack_T_62"
            T_des = 62 + 273.15  # K. It is the temperature of the fuel cell.
        elseif type_fuel_cell == "ZSW-GenStack_T_76"
            T_des = 76 + 273.15  # K. It is the temperature of the fuel cell.
        elseif type_fuel_cell == "ZSW-GenStack_T_84"
            T_des = 84 + 273.15  # K. It is the temperature of the fuel cell.
        else
            T_des = 68 + 273.15  # K. It is the nominal temperature of the fuel cell.
        end
        Sa, Sc = 1.6, 1.6  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.398, 0.50  # It is the desired relative humidity.
        if type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
            Pa_des, Pc_des = 1.61e5, 1.41e5  # Pa. It is the desired pressure of the inlet fuel gas (at the anode/cathode).
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
            Pa_des, Pc_des = 2.01e5, 1.81e5  # Pa. It is the desired pressure of the inlet fuel gas (at the anode/cathode).
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
            Pa_des, Pc_des = 2.4e5, 2.2e5  # Pa. It is the desired pressure of the inlet fuel gas (at the anode/cathode).
        elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
            Pa_des, Pc_des = 2.8e5, 2.6e5  # Pa. It is the desired pressure of the inlet fuel gas (at the anode/cathode).
        else
            Pa_des, Pc_des = 2.2e5, 2.0e5  # Pa. It is the desired pressure of the inlet fuel gas (at the anode/cathode).
        end
        y_H2_in = 0.7  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola = 2.5e4  # (A.m-2). It is the maximum current density for the polarization curve.
        else  # voltage_zone == "before_voltage_drop"
            i_max_pola = 1.9e4
        end
        pola_current_parameters["i_max_pola"] = i_max_pola

        # Fuel cell physical parameters
        Aact = 283.87e-4  # m². It is the active area of the catalyst layer.
        nb_cell = 26  # . It is the number of cell in the stack.
        Hagc = 230e-6  # m. It is the thickness of the anode gas channel.
        Hcgc = 300e-6  # m. It is the thickness of the cathode gas channel.
        Wagc = 430e-6  # m. It is the width of the anode gas channel.
        Wcgc = 532e-6  # m. It is the width of the cathode gas channel.
        Lgc = 246.2e-3  # m. It is the length of the gas channel.
        nb_channel_in_gc = 105  # . It is the number of channels in the bipolar plate.
        Ldist = 71.1e-3  # m. It is the length of the distributor, which is the volume between the gas channel and the manifold.
        Lm = 25.8e-3  # m. It is the length of the manifold.
        A_T_a = 9.01e-4  # m². It is the inlet/exhaust anode manifold throttle area
        A_T_c = 22.61e-4  # m². It is the inlet/exhaust cathode manifold throttle area
        Vasm, Vcsm = Lm * A_T_a, Lm * A_T_c  # m3. It is the supply manifold volume.
        Vaem, Vcem = Vasm, Vcsm  # m-3. It is the exhaust manifold volume.

        # Fuel cell undetermined physical parameters.
        Hgdl = 127e-6  # m. It is the thickness of the gas diffusion layer.
        Hmpl = 70e-6  # m. It is the thickness of the microporous layer.
        epsilon_c = 0.2  # It is the compression ratio of the GDL.

        # Estimated undetermined parameters for the initialisation
        # Gas diffusion layer
        epsilon_gdl = 0.788  # It is the anode/cathode GDL porosity.
        epsilon_mpl = 0.425  # It is the porosity of the microporous layer.
        # Catalyst layer
        Hacl = 8e-6  # m. It is the thickness of the anode catalyst layer.
        Hccl = 17e-6  # m. It is the thickness of the cathode catalyst layer.
        # Membrane
        Hmem = 15e-6  # m. It is the thickness of the membrane.
        # Interaction parameters between water and PEMFC structure
        e = 4  # It is the capillary exponent
        K_l_ads = 1  # . It is an estimation of the ratio between the liquid and vapor sorption rates of water in the membrane. It should be in [10-1000] [shaoNewInsightsSteadystate2023].
        K_O2_ad_Pt = 5.4  # . It is the interfacial resistance coefficient of O2 adsorption on the Pt sites.
        # Voltage polarization
        Re = 1e-06  # ohm.m². It is the electron conduction resistance of the circuit.
        i0_c_ref = 14.43  # A.m-2. It is the dry reference exchange current density at the cathode.
        kappa_co = 5  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 1.026  # It is the overpotential correction exponent.
        C_scl = 2e7  # F.m-3. It is the volumetric space-charge layer capacitance.

        # Computing parameters
        nb_gc = 1  # It is the number of model points placed inside each gas channel.
        nb_gdl = 3  # It is the number of model points placed inside each GDL.
        nb_mpl = 2  # It is the number of model points placed inside each MPL.

    elseif type_fuel_cell in ("EH-31_1.5", "EH-31_2.0", "EH-31_2.25", "EH-31_2.5")
        # Given values by the author
        # Operating inputs
        T_des = 74 + 273.15  # K. It is the temperature of the fuel cell.
        Sa, Sc = 1.2, 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des, Phi_c_des = 0.4, 0.6  # It is the desired relative humidity.
        if type_fuel_cell == "EH-31_1.5"
            Pa_des, Pc_des = 1.5e5, 1.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        elseif type_fuel_cell == "EH-31_2.0"
            Pa_des, Pc_des = 2.0e5, 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        elseif type_fuel_cell == "EH-31_2.25"
            Pa_des, Pc_des = 2.25e5, 2.25e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        else  # type_fuel_cell == "EH-31_2.5"
            Pa_des, Pc_des = 2.5e5, 2.5e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        end
        y_H2_in = 1  # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola = 3.0e4  # (A.m-2). It is the maximum current density for the polarization curve.
        else  # voltage_zone == "before_voltage_drop"
            i_max_pola = 1.7e4
        end
        pola_current_parameters["i_max_pola"] = i_max_pola

        # Fuel cell physical parameters
        Aact = 85e-4  # m². It is the active area of the catalyst layer.
        Wagc = 450e-6  # m. It is the width of the anode gas channel.
        Wcgc = Wagc  # m. It is the width of the cathode gas channel.
        Lgc = 144e-3  # m. It is the length of one channel in the bipolar plate.
        nb_channel_in_gc = 67  # . It is the number of channels in the bipolar plate.

        # Extrapolated physical parameters
        nb_cell = 1  # . It is the number of cell in the stack.
        Hgdl = 200e-6  # m. It is the thickness of the gas diffusion layer.
        Hmpl = 30e-6  # m. It is the thickness of the microporous layer.
        epsilon_mpl = 0.4  # It is the porosity of the microporous layer.
        Hagc = 500e-6  # m. It is the thickness of the anode gas channel.
        Hcgc = Hagc  # m. It is the thickness of the cathode gas channel.
        Ldist = 5e-2  # m. It is the estimated length of the distributor, which is the volume between the gas channel and the manifold.
        Lm = 2.03e-3  # m. It is the length of the manifold.
        A_T_a = 11.8e-4  # m². It is the inlet/exhaust anode manifold throttle area
        A_T_c = 34.4e-4  # m². It is the inlet/exhaust cathode manifold throttle area
        Vasm, Vcsm = Lm * A_T_a, Lm * A_T_c  # m3. It is the supply manifold volume.
        Vaem, Vcem = Vasm, Vcsm  # m-3. It is the exhaust manifold volume.

        # Estimated undetermined parameters for the initialisation
        # Gas diffusion layer
        epsilon_gdl = 0.5002  # It is the anode/cathode GDL porosity.
        epsilon_c = 0.2  # It is the compression ratio of the GDL.
        # Catalyst layer
        Hacl = 8.593e-6  # m. It is the thickness of the anode catalyst layer.
        Hccl = Hacl  # m. It is the thickness of the cathode catalyst layer.
        # Membrane
        Hmem = 16.06e-6  # m. It is the thickness of the membrane.
        # Interaction parameters between water and PEMFC structure
        e = 4  # It is the capillary exponent
        K_l_ads = 1  # . It is an estimation of the ratio between the liquid and vapor sorption rates of water in the membrane. It should be in [10-1000] [shaoNewInsightsSteadystate2023].
        K_O2_ad_Pt = 5.4  # . It is the interfacial resistance coefficient of O2 adsorption on the Pt sites.
        # Voltage polarization
        Re = 1e-06  # ohm.m². It is the electron conduction resistance of the circuit.
        i0_c_ref = 14.43  # A.m-2. It is the reference exchange current density at the cathode.
        kappa_co = 30.42  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c = 0.4152  # It is the overpotential correction exponent.
        C_scl = 20e6  # F.m-3. It is the volumetric space-charge layer capacitance.

        # Computing parameters
        nb_gc = 1  # It is the number of model points placed inside each gas channel.
        nb_gdl = 3  # It is the number of model points placed inside each GDL.
        nb_mpl = 2  # It is the number of model points placed inside each MPL.

    else
        throw(ArgumentError("A correct type_fuel_cell should be given."))
    end

    # Initialize the operating inputs and parameters dictionaries.
    operating_inputs = Dict("current_density" => current_density, "T_des" => T_des, "Pa_des" => Pa_des,
                           "Pc_des" => Pc_des, "Sa" => Sa, "Sc" => Sc, "Phi_a_des" => Phi_a_des,
                           "Phi_c_des" => Phi_c_des, "y_H2_in" => y_H2_in)
    current_parameters = Dict("step_current_parameters" => step_current_parameters,
                             "pola_current_parameters" => pola_current_parameters,
                             "pola_current_for_cali_parameters" => pola_current_for_cali_parameters,
                             "i_EIS" => i_EIS, "ratio_EIS" => ratio_EIS, "t_EIS" => t_EIS, "f_EIS" => f_EIS)
    accessible_physical_parameters = Dict("Aact" => Aact, "nb_cell" => nb_cell, "Hagc" => Hagc, "Hcgc" => Hcgc,
                                         "Wagc" => Wagc, "Wcgc" => Wcgc, "Lgc" => Lgc,
                                         "nb_channel_in_gc" => nb_channel_in_gc, "Ldist" => Ldist, "Lm" => Lm,
                                         "A_T_a" => A_T_a, "A_T_c" => A_T_c, "Vasm" => Vasm, "Vcsm" => Vcsm,
                                         "Vaem" => Vaem, "Vcem" => Vcem)
    undetermined_physical_parameters = Dict("Hgdl" => Hgdl, "Hmpl" => Hmpl, "Hmem" => Hmem, "Hacl" => Hacl,
                                           "Hccl" => Hccl, "epsilon_gdl" => epsilon_gdl,
                                           "epsilon_mpl" => epsilon_mpl, "epsilon_c" => epsilon_c, "e" => e,
                                           "K_l_ads" => K_l_ads, "K_O2_ad_Pt" => K_O2_ad_Pt, "Re" => Re,
                                           "i0_c_ref" => i0_c_ref, "kappa_co" => kappa_co, "kappa_c" => kappa_c,
                                           "C_scl" => C_scl)
    model_parameters = Dict("nb_gc" => nb_gc, "nb_gdl" => nb_gdl, "nb_mpl" => nb_mpl, "t_purge" => t_purge,
                           "rtol" => rtol, "atol" => atol)
    computing_parameters = Dict("type_fuel_cell" => type_fuel_cell, "type_current" => type_current,
                               "voltage_zone" => voltage_zone, "type_auxiliary" => type_auxiliary,
                               "type_purge" => type_purge, "type_display" => type_display,
                               "type_plot" => type_plot)

    # Characteristic points of the experimental polarization curve
    i_exp, U_exp = pola_exp_values_calibration(type_fuel_cell, voltage_zone)

    return (operating_inputs, current_parameters, accessible_physical_parameters, undetermined_physical_parameters,
            model_parameters, computing_parameters, i_exp, U_exp)
end


"""Update the undetermined physical parameters dictionary with values from the solution.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `solution::Vector{Float64}`: List of parameter values obtained from the optimization algorithm.
- `varbound::Vector{Vector}`: List of parameter bounds and names. Each element contains the parameter name at index 1.
- `undetermined_physical_parameters::Dict`: Dictionary of undetermined physical parameters to be updated.

# Returns
- `Dict`: Updated dictionary of undetermined physical parameters.
"""
function update_undetermined_parameters(type_fuel_cell::String,
                                       solution::Vector{Float64},
                                       varbound::Vector{Vector},
                                       undetermined_physical_parameters::Dict)::Dict
    for i in eachindex(solution)
        param_name = varbound[i][1]
        if haskey(undetermined_physical_parameters, param_name)
            undetermined_physical_parameters[param_name] = solution[i]
        end
        if type_fuel_cell in ("EH-31_1.5", "EH-31_2.0", "EH-31_2.25", "EH-31_2.5")
            undetermined_physical_parameters["Hccl"] = undetermined_physical_parameters["Hacl"]
        end
    end
    return undetermined_physical_parameters
end


"""Calculate the simulation maximal error between the experimental and simulated polarization curves.

Two simulations on different operating conditions and on the same stack, and so two set of experimental data,
are considered as it is the minimum amount of data which is required for the calibration.

# Arguments
- `simulator_1::AlphaPEM`: PEM simulator which contains the simulation results for the first simulation.
- `U_exp_1::Vector{Float64}`: Experimental values of the voltage for the first simulation.
- `i_exp_1::Vector{Float64}`: Experimental values of the current density for the first simulation.
- `simulator_2::AlphaPEM`: PEM simulator which contains the simulation results for the second simulation.
- `U_exp_2::Vector{Float64}`: Experimental values of the voltage for the second simulation.
- `i_exp_2::Vector{Float64}`: Experimental values of the current density for the second simulation.

# Returns
- `sim_error::Float64`: Maximum error between the experimental and simulated polarization curves in percentage.
"""
function calculate_simulation_error(simulator_1::AlphaPEM,
                                   U_exp_1::Vector{Float64},
                                   i_exp_1::Vector{Float64},
                                   simulator_2::AlphaPEM,
                                   U_exp_2::Vector{Float64},
                                   i_exp_2::Vector{Float64})::Float64

    # Recovery of ifc_1
    t1 = collect(simulator_1.variables["t"])
    ifc_t_1 = [simulator_1.operating_inputs["current_density"](t, simulator_1.parameters) for t in t1]

    # Recovery of ifc_2
    t2 = collect(simulator_2.variables["t"])
    ifc_t_2 = [simulator_2.operating_inputs["current_density"](t, simulator_2.parameters) for t in t2]

    # Polarisation curve point recovery after stack stabilisation for Simulator1
    # Extraction of the parameters
    delta_t_ini_pola_cali_1 = simulator_1.parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]  # (s).
    delta_t_load_pola_cali_1 = simulator_1.parameters["pola_current_for_cali_parameters"]["delta_t_load_pola_cali"]  # (s).
    delta_t_break_pola_cali_1 = simulator_1.parameters["pola_current_for_cali_parameters"]["delta_t_break_pola_cali"]  # (s).

    # Calculation
    nb_loads1 = length(i_exp_1)  # Number of load which are made
    delta_t_cali_1 = delta_t_load_pola_cali_1 + delta_t_break_pola_cali_1  # s. It is the time of one load.
    ifc_discretized1 = [ifc_t_1[argmin(abs.(t1 .- (delta_t_ini_pola_cali_1 + i * delta_t_cali_1)))] for i in 1:nb_loads1]
    Ucell_discretized1 = [simulator_1.variables["Ucell"][argmin(abs.(t1 .- (delta_t_ini_pola_cali_1 + i * delta_t_cali_1)))] for i in 1:nb_loads1]

    # Polarisation curve point recovery after stack stabilisation for Simulator2
    # Extraction of the parameters
    delta_t_ini_pola_cali_2 = simulator_2.parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]  # (s).
    delta_t_load_pola_cali_2 = simulator_2.parameters["pola_current_for_cali_parameters"]["delta_t_load_pola_cali"]  # (s).
    delta_t_break_pola_cali_2 = simulator_2.parameters["pola_current_for_cali_parameters"]["delta_t_break_pola_cali"]  # (s).

    # Calculation
    nb_loads2 = length(i_exp_2)  # Number of load which are made
    delta_t_cali_2 = delta_t_load_pola_cali_2 + delta_t_break_pola_cali_2  # s. It is the time of one load.
    ifc_discretized2 = [ifc_t_2[argmin(abs.(t2 .- (delta_t_ini_pola_cali_2 + i * delta_t_cali_2)))] for i in 1:nb_loads2]
    Ucell_discretized2 = [simulator_2.variables["Ucell"][argmin(abs.(t2 .- (delta_t_ini_pola_cali_2 + i * delta_t_cali_2)))] for i in 1:nb_loads2]

    # Distance between the simulated and the experimental polarization curves (RMSE: root-mean-square error).
    rmse1 = sqrt(mean(((Ucell_discretized1 .- U_exp_1) ./ U_exp_1 .* 100) .^ 2))
    rmse2 = sqrt(mean(((Ucell_discretized2 .- U_exp_2) ./ U_exp_2 .* 100) .^ 2))
    sim_error = (rmse1 + rmse2) / 2  # in %.

    return sim_error
end


"""Print the calibration results by associating each optimized value with its parameter name.

# Arguments
- `convergence::Vector`: Convergence information of the genetic algorithm.
- `ga_instance`: Instance of PyGAD used for optimization.
- `solution::Vector{Float64}`: List of optimized parameter values.
- `varbound::Vector{Vector}`: List of parameter bounds and names.
- `sim_error::Float64`: Simulation error (RMSE) in percentage.

# Returns
- `Nothing`
"""
function print_calibration_results(convergence::Vector,
                                  ga_instance,
                                  solution::Vector{Float64},
                                  varbound::Vector{Vector},
                                  sim_error::Float64)::Nothing
    println("Convergence:\n", convergence)
    for idx in eachindex(solution)
        param_name = varbound[idx][1]
        println("Optimized parameter $param_name: $(solution[idx])")
    end
    println(Fore.RED, "\nSimulation error (RMSE): ", sim_error, " %")
    println(Style.RESET_ALL)
    best_solution_generation = hasattr(ga_instance, "best_solution_generation") ? py"int"(ga_instance.best_solution_generation) : -1
    if best_solution_generation != -1
        println("Best fitness value reached after $best_solution_generation generations.")
    end
    return nothing
end


"""Save the calibration results in a text file and a PyGAD file.

The optimized values are retrieved from the solution list and associated with their names via varbound.

# Arguments
- `convergence::Vector`: Convergence information from the genetic algorithm.
- `ga_instance`: Instance of PyGAD used for optimization.
- `solution::Vector{Float64}`: List of optimized parameter values.
- `varbound::Vector{Vector}`: List of parameter bounds and names.
- `sim_error::Float64`: Simulation error (RMSE) in percentage.
- `type_fuel_cell::String`: Type of fuel cell configuration.

# Returns
- `Nothing`
"""
function save_calibration_results(convergence::Vector,
                                 ga_instance,
                                 solution::Vector{Float64},
                                 varbound::Vector{Vector},
                                 sim_error::Float64,
                                 type_fuel_cell::String)::Nothing
    using Paths
    using Base.Filesystem

    # Resolve current file and locate project root using common repo markers
    cur = @__FILE__
    markers = [".git", "pyproject.toml", "setup.cfg", "requirements.txt", "Pipfile"]
    project_root = nothing
    parent = dirname(cur)
    while true
        if any(ispath(joinpath(parent, m)) for m in markers)
            project_root = parent
            break
        end
        new_parent = dirname(parent)
        new_parent == parent && break
        parent = new_parent
    end
    # Fallback to current working directory if no marker found.
    project_root === nothing && (project_root = pwd())

    # Prepare folder and base filename
    root_folder = "results"
    subfolder_name = if '_' in type_fuel_cell
        split(type_fuel_cell, '_')[1]
    else
        type_fuel_cell
    end
    folder_path = joinpath(project_root, root_folder, subfolder_name)
    mkpath(folder_path)

    # Find a non-colliding filename parameter_calibration_N.txt
    base_stem = "parameter_calibration"
    counter = 1
    while true
        txt_candidate = joinpath(folder_path, "$(base_stem)_$(counter).txt")
        if !isfile(txt_candidate)
            file_path = txt_candidate
            break
        end
        counter += 1
    end

    # Write the text file
    open(file_path, "w") do file
        write(file, "Convergence: " * string(convergence))
        for idx in eachindex(solution)
            param_name = varbound[idx][1]
            write(file, "\nOptimized parameter $param_name: $(solution[idx])")
        end
        write(file, "\nSimulation error (RMSE): " * string(sim_error) * "%")
        write(file, "\nAlgorithm works with " * type_fuel_cell * ".")
        best_solution_generation = hasattr(ga_instance, "best_solution_generation") ? py"int"(ga_instance.best_solution_generation) : -1
        if best_solution_generation != -1
            write(file, "\nBest fitness value reached after $best_solution_generation generations.")
        end
    end

    # Save the PyGAD instance in the same folder (without extension as PyGAD will add it)
    ga_instance.save(filename=joinpath(folder_path, "$(base_stem)_$(counter)"))

    # Remove ongoing pickle at project root if exists
    ongoing = joinpath(project_root, "parameter_calibration_ongoing.pkl")
    if isfile(ongoing)
        try
            rm(ongoing)
        catch
            # ignore deletion errors
        end
    end

    return nothing
end

