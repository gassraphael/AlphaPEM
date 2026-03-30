# -*- coding: utf-8 -*-

"""This file is used to set the parameters of the fuel cell system."""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries

# Importing functions
include(joinpath(@__DIR__, "parameters_specific.jl"))


# _______________________________________________________Settings_______________________________________________________

"""
    calculate_operating_inputs(pola_current_parameters::Dict{<:String, <:Real},
                               type_fuel_cell::Union{Nothing, String},
                               voltage_zone::String)

This function is used to set the operating inputs of the fuel cell system.

# Arguments
- `pola_current_parameters::Dict{String, Float64}`: Parameters for the polarization current density
  function.
- `type_fuel_cell::Union{Nothing, String}`: Type of fuel cell system.
- `voltage_zone::String`: Zone of the polarization curve which is considered. It can be `"full"` or
  `"before_voltage_drop"`.

# Returns
- `T_des`: Desired fuel cell temperature in Kelvin.
- `Pa_des: Desired anode pressure in Pascal.
- `Pc_des: Desired cathode pressure in Pascal.
- `Sa: Stoichiometric ratio of hydrogen.
- `Sc: Stoichiometric ratio of oxygen.
- `Phi_a_des: Desired anode relative humidity.
- `Phi_c_des: Desired cathode relative humidity.
- `y_H2_in: Molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the
  inlet.
- `pola_current_parameters::Dict`: Parameters for the polarization current density. It
  is a dictionary containing:
  - `"delta_t_ini_pola"`: the initial time (in seconds) at zero current density for the stabilisation of the internal
    states,
  - `"delta_t_load_pola"`: the loading time (in seconds) for one step current of the polarisation current density
    function,
  - `"delta_t_break_pola"`: the breaking time (in seconds) for one step current, for the stabilisation of the internal
    states,
  - `"delta_i_pola"`: the current density step (in A.m-2) for the polarisation current density function,
  - `"i_max_pola"`: the maximum current density (in A.m-2) for the polarization curve.
"""
function calculate_operating_inputs(pola_current_parameters::Dict{String, Float64},
                                    type_fuel_cell::String,
                                    voltage_zone::String)

    if type_fuel_cell == "manual_setup" # Setup which are not stored in "stored_operating_inputs".
        T_des::Float64 = 74.0 + 273.15  # K. It is the desired fuel cell temperature.
        Pa_des::Float64 = 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Pc_des::Float64 = 2.0e5  # Pa. It is the desired pressure of the fuel gas (at the anode/cathode).
        Sa::Float64 = 1.2  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Sc::Float64 = 2.0  # It is the stoichiometric ratio (of hydrogen and oxygen).
        Phi_a_des::Float64 = 0.4  # It is the desired relative humidity.
        Phi_c_des::Float64 = 0.6  # It is the desired relative humidity.
        y_H2_in::Float64 = 1.0 # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        i_max_pola::Float64 = 3.0e4  # A.m-2. It is the maximum current density for the polarization curve.
    else # Stored setup in "stored_operating_inputs".
        T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, i_max_pola = stored_operating_inputs(type_fuel_cell, voltage_zone)
    end

    pola_current_parameters["i_max_pola"] = i_max_pola  # Update the maximum current density for the polarization curve.
    return T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, pola_current_parameters
end


"""
    calculate_physical_parameters(type_fuel_cell::String)

This function is used to set the physical parameters of the fuel cell system.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell system.

# Returns
- `Hacl::Float64`: Thickness of the anode catalyst layer in meters.
- `Hccl::Float64`: Thickness of the cathode catalyst layer in meters.
- `Hmem::Float64`: Thickness of the membrane in meters.
- `Hgdl::Float64`: Thickness of the gas diffusion layer in meters.
- `epsilon_gdl::Float64`: Anode/cathode GDL porosity.
- `epsilon_c::Float64`: Compression ratio of the GDL.
- `Hmpl::Float64`: Thickness of the microporous layer in meters.
- `epsilon_mpl::Float64`: Porosity of the microporous layer.
- `Hagc::Float64`: Thickness of the anode gas channel in meters.
- `Hcgc::Float64`: Thickness of the cathode gas channel in meters.
- `Wagc::Float64`: Width of the anode gas channel in meters.
- `Wcgc::Float64`: Width of the cathode gas channel in meters.
- `Lgc::Float64`: Length of the gas channel in meters.
- `nb_channel_in_gc::Int64`: Number of channels in the bipolar plate.
- `Ldist::Float64`: Length of the distributor, which is the volume between the gas channel and the manifold, in m.
- `Lm::Float64`: Length of the manifold in m.
- `A_T_a::Float64`: Inlet/exhaust anode manifold throttle area in m².
- `A_T_c::Float64`: Inlet/exhaust cathode manifold throttle area in m².
- `Vasm::Float64`: Supply manifold volume at the anode in m³.
- `Vcsm::Float64`: Supply manifold volume at the cathode in m³.
- `Vaem::Float64`: Exhaust manifold volume at the anode in m³.
- `Vcem::Float64`: Exhaust manifold volume at the cathode in m³.
- `Aact::Float64`: Active area of the catalyst layer in m².
- `nb_cell::Int64`: Number of cell in the stack.
- `e::Int64`: Capillary exponent.
- `K_l_ads::Float64`: Ratio between the liquid and vapor sorption rates of water in the membrane.
- `K_O2_ad_Pt::Float64`: Interfacial resistance coefficient of O2 adsorption on the Pt sites.
- `Re::Float64`: Electron conduction resistance of the circuit in Ω.m².
- `i0_c_ref::Float64`: Reference exchange current density at the cathode in A.m-2.
- `kappa_co::Float64`: Crossover correction coefficient in mol.m-1.s-1.Pa-1.
- `kappa_c::Float64`: Overpotential correction exponent.
- `C_scl::Float64`: Volumetric space-charge layer capacitance in F.m-3.
"""
function calculate_physical_parameters(type_fuel_cell::String)

    if type_fuel_cell == "manual_setup" # Setup which are not stored in "stored_physical_parameters".
      # Fuel cell physical parameters: 𝜔 (which are not controllable by the system)
      # Global
      Aact::Float64 = 279.72e-4  # m². It is the MEA active area.
      nb_cell::Int64 = 1  # . It is the number of cell in the stack.
      #   Catalyst layer
      Hacl::Float64 = 8.089e-6  # m. It is the thickness of the anode catalyst layer.
      Hccl::Float64 = Hacl  # m. It is the thickness of the cathode catalyst layer.
      #   Membrane
      Hmem::Float64 = 2.0e-5  # m. It is the thickness of the membrane.
      #   Gas diffusion layer
      Hgdl::Float64 = 2.0e-4  # m. It is the thickness of the gas diffusion layer.
      epsilon_gdl::Float64 = 0.7011156494971454  # It is the anode/cathode GDL porosity.
      epsilon_c::Float64 = 0.27052745219052654  # It is the compression ratio of the GDL.
      #   Microporous layer
      Hmpl::Float64 = 3.0e-5  # m. It is the thickness of the microporous layer.
      epsilon_mpl::Float64 = 0.4  # It is the porosity of the microporous layer.
      #   Gas channel
      Hagc::Float64 = 5.0e-4  # m. It is the thickness of the anode gas channel.
      Hcgc::Float64 = Hagc  # m. It is the thickness of the cathode gas channel.
      Wagc::Float64 = 4.5e-4  # m. It is the width of the anode gas channel.
      Wcgc::Float64 = Wagc  # m. It is the width of the cathode gas channel.
      Lgc::Float64 = 144.0e-3  # m. It is the length of one channel in the bipolar plate.
      nb_channel_in_gc::Int64 = 67  # . It is the number of channels in the bipolar plate.
      Ldist::Float64 = 5.0e-2  # m. It is the estimated length of the distributor, which is the volume between the gas channel and the manifold.
      #   Auxiliaries
      Lm::Float64 = 25.8e-3  # m. It is the length of the manifold.
      L_endplate::Float64 = 46.8e-3  # m. It is the length of the endplate.
      A_T_a::Float64 = 11.8e-4  # m². It is the inlet/exhaust anode manifold throttle area
      A_T_c::Float64 = A_T_a  # m². It is the inlet/exhaust cathode manifold throttle area
      Vasm::Float64 = 7000.0e-6  # m3. It is the supply manifold volume.
      Vcsm::Float64 = 7000.0e-6  # m3. It is the supply manifold volume.
      Vaem::Float64 = 2400.0e-6  # m-3. It is the exhaust manifold volume.
      Vcem::Float64 = 2400.0e-6  # m-3. It is the exhaust manifold volume.
      V_endplate_a::Float64 = 33.6e-6  # m3. It is the anode endplate volume.
      V_endplate_c::Float64 = 86.6e-6  # m3. It is the cathode endplate volume.
      #   Interaction parameters between fluids and PEMFC structure
      e::Int64 = 5  # It is the capillary exponent
      K_l_ads::Float64 = 1.0  # . It is the ratio between the liquid and vapor sorption rates of water in the membrane. It should be in [10-1000] [shaoNewInsightsSteadystate2023].
      K_O2_ad_Pt::Float64 = 5.4  # . It is the interfacial resistance coefficient of O2 adsorption on the Pt sites.
      #   Voltage polarization
      Re::Float64 = 1.0e-06  # Ω.m². It is the electron conduction resistance of the circuit.
      i0_c_ref::Float64 = 14.43  # A.m-2. It is the dry reference exchange current density at the cathode.
      kappa_co::Float64 = 29.793535549174077  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
      kappa_c::Float64 = 1.6136446641573106  # It is the overpotential correction exponent.
      C_scl::Float64 = 2.0e7  # F.m-3. It is the volumetric space-charge layer capacitance.
    else # Stored setup in "stored_physical_parameters".
        (Hacl, Hccl, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc, Lgc,
         nb_channel_in_gc, Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem, Aact, nb_cell, e, K_l_ads, K_O2_ad_Pt, Re,
         i0_c_ref, kappa_co, kappa_c, C_scl) = stored_physical_parameters(type_fuel_cell)
    end

    return (Hacl, Hccl, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc, Lgc,
            nb_channel_in_gc, Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem, Aact, nb_cell, e, K_l_ads, K_O2_ad_Pt,
            Re, i0_c_ref, kappa_co, kappa_c, C_scl)
end


"""
    calculate_computing_parameters(step_current_parameters::Union{Nothing, Dict}=nothing)

This function is used to set the computing parameters of the fuel cell system.

# Arguments
- `step_current_parameters::Union{Nothing, Dict}`: Parameters for the step current
  density function.

# Returns
- `nb_gc::Int64`: Number of model nodes placed inside each gas channel.
- `nb_gdl::Int64`: Number of model nodes placed inside each GDL.
- `nb_mpl::Int64`: Number of model nodes placed inside each MPL.
- `purge_time::Float64`: The time for purging the system in seconds.
- `delta_purge::Float64`: The time between two purges in seconds.
- `rtol::Float64`: Relative tolerance for the system of ODEs solver.
- `atol::Float64`: Absolute tolerance for the system of ODEs solver.
"""
function calculate_computing_parameters(step_current_parameters::Union{Nothing, Dict}=nothing)

    # Setting the number of model points placed inside each layer:
    nb_gc::Int64 = 1  # It is the number of model points placed inside each gas channel.
    nb_gdl::Int64 = 3  # It is the number of model points placed inside each GDL.
    nb_mpl::Int64 = 2  # It is the number of model points placed inside each MPL.

    # Setting the purging parameters of the system and the dynamic display of the step current density function:
    purge_time::Float64 = 0.6  # (s). It is the purge time.
    delta_purge::Float64 = 15.0  # (s). It is the time between two purges.
    delta_t_dyn_step::Float64 = 0.1  # (s). Time for dynamic display of the step current density function.

    # Setting the tolerances for the system of ODEs solver:
    rtol::Float64 = 1.0e-6  # Relative tolerance for the system of ODEs solver.
    atol::Float64 = 1.0e-9  # Absolute tolerance for the system of ODEs solver.

    # Update the step current parameters.
    if !isnothing(step_current_parameters)
        step_current_parameters["delta_t_dyn_step"] = delta_t_dyn_step
    end
    return nb_gc, nb_gdl, nb_mpl, purge_time, delta_purge, rtol, atol
end


