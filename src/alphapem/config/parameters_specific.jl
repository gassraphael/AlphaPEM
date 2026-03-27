# -*- coding: utf-8 -*-

# This module contains some of the required functions for the settings.

# ___________________________________________________Settings modules___________________________________________________

"""
    stored_operating_inputs(type_fuel_cell::String, voltage_zone::String)
        -> Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}

This function gives the operating inputs which correspond to the given `type_fuel_cell`.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell configuration.
- `voltage_zone::String`:   Zone of the polarization curve which is considered.

# Returns
- `T_des::Float64`:       Desired fuel cell temperature in Kelvin.
- `Pa_des::Float64`:      Desired anode pressure in Pascal.
- `Pc_des::Float64`:      Desired cathode pressure in Pascal.
- `Sa::Float64`:          Stoichiometric ratio of hydrogen.
- `Sc::Float64`:          Stoichiometric ratio of oxygen.
- `Phi_a_des::Float64`:   Desired anode relative humidity.
- `Phi_c_des::Float64`:   Desired cathode relative humidity.
- `y_H2_in::Float64`:     Molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
- `i_max_pola::Float64`:  Maximum current density for the polarization curve in A⋅m⁻².
"""
function stored_operating_inputs(type_fuel_cell::String, voltage_zone::String)

    # For the ZSW Generic Stack fuel cell
    if type_fuel_cell == "ZSW-GenStack"
        T_des::Float64          = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des::Float64         = 2.2e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des::Float64         = 2.0e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa::Float64             = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc::Float64             = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des::Float64      = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des::Float64      = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in::Float64        = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola::Float64 = 2.500e4        # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.700e4        # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 1.61e5         # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 1.41e5         # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.150e4               # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 0.700e4               # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                = 2.01e5         # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 1.81e5         # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.450e4               # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.300e4               # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.4e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 2.2e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.500e4               # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.900e4               # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.8e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 2.6e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.500e4               # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.900e4               # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "ZSW-GenStack_T_62"
        T_des                   = 62.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.2e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.500e4   # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.500e4   # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "ZSW-GenStack_T_76"
        T_des                   = 76.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.2e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.500e4   # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.100e4   # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "ZSW-GenStack_T_84"
        T_des                   = 84.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.2e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.000e4   # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 0.700e4   # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end

    # For EH-31 fuel cell
    elseif type_fuel_cell == "EH-31_1.5"
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 1.5e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 1.5e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.250e4               # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.700e4               # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "EH-31_2.0"
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.500e4               # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.300e4               # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "EH-31_2.25"
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.25e5         # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.25e5         # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 2.800e4                # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.700e4                # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end
    elseif type_fuel_cell == "EH-31_2.5"
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.5e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.5e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
        if voltage_zone == "full"
            i_max_pola          = 3.000e4               # A.m-2. It is the maximum current density for the polarization curve.
        elseif voltage_zone == "before_voltage_drop"
            i_max_pola          = 1.600e4               # A.m-2. It is the maximum current density for the polarization curve.
        else
            error("the voltage_zone given is not valid.")
        end

    # For other fuel cells
    else
        error("the type_fuel_cell given is not valid.")
    end

    return T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in, i_max_pola
end


"""
    stored_physical_parameters(type_fuel_cell::String)
        -> Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64,
                 Float64, Float64, Float64, Float64, Float64, Int64,   Float64, Float64,
                 Float64, Float64, Float64, Float64, Float64, Float64, Float64, Int64,
                 Int64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}

This function gives the physical parameters which correspond to the given `type_fuel_cell`.

# Arguments
- `type_fuel_cell::String`: Type of fuel cell configuration.

# Returns
- `Hacl::Float64`:            Thickness of the anode catalyst layer in m.
- `Hccl::Float64`:            Thickness of the cathode catalyst layer in m.
- `Hmem::Float64`:            Thickness of the membrane in m.
- `Hgdl::Float64`:            Thickness of the gas diffusion layer in m.
- `epsilon_gdl::Float64`:     Anode/cathode GDL porosity.
- `epsilon_c::Float64`:       Compression ratio of the GDL.
- `Hmpl::Float64`:            Thickness of the microporous layer in m.
- `epsilon_mpl::Float64`:     Porosity of the microporous layer.
- `Hagc::Float64`:            Thickness of the anode gas channel in m.
- `Hcgc::Float64`:            Thickness of the cathode gas channel in m.
- `Wagc::Float64`:            Width of the anode gas channel in m.
- `Wcgc::Float64`:            Width of the cathode gas channel in m.
- `Lgc::Float64`:             Length of the gas channel in m.
- `nb_channel_in_gc::Int64`:  Number of channels in the bipolar plate.
- `Ldist::Float64`:           Length of the distributor (volume between the gas channel and the manifold) in m.
- `Lm::Float64`:              Length of the manifold in m.
- `A_T_a::Float64`:           Inlet/exhaust anode manifold throttle area in m².
- `A_T_c::Float64`:           Inlet/exhaust cathode manifold throttle area in m².
- `Vasm::Float64`:            Anode supply manifold volume in m³.
- `Vcsm::Float64`:            Cathode supply manifold volume in m³.
- `Vaem::Float64`:            Anode exhaust manifold volume in m³.
- `Vcem::Float64`:            Cathode exhaust manifold volume in m³.
- `Aact::Float64`:            Active area of the cell in m².
- `nb_cell::Int64`:           Number of cells in the stack.
- `e::Int64`:               Capillary exponent.
- `K_l_ads::Float64`:         Ratio between the liquid and vapor sorption rates of water in the membrane.
- `K_O2_ad_Pt::Float64`:      Interfacial resistance coefficient of O2 adsorption on the Pt sites.
- `Re::Float64`:              Electron conduction resistance of the circuit in Ω⋅m².
- `i0_c_ref::Float64`:        Reference exchange current density at the cathode in A⋅m⁻².
- `kappa_co::Float64`:        Crossover correction coefficient in mol⋅m⁻¹⋅s⁻¹⋅Pa⁻¹.
- `kappa_c::Float64`:         Overpotential correction exponent.
- `C_scl::Float64`:           Volumetric space-charge layer capacitance in F⋅m⁻³.
"""
function stored_physical_parameters(type_fuel_cell::String)

    # For the ZSW Generic Stack fuel cell
    if type_fuel_cell == "ZSW-GenStack" || type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41" ||
            type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81" || type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2" ||
            type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6" || type_fuel_cell == "ZSW-GenStack_T_62" ||
            type_fuel_cell == "ZSW-GenStack_T_76" || type_fuel_cell == "ZSW-GenStack_T_84"
        # Global
        Aact::Float64             = 283.87e-4              # m². It is the MEA active area.
        nb_cell::Int64            = 26                     # .   It is the number of cells in the stack.
        # Catalyst layer
        Hacl::Float64             = 8.112569325675836e-6   # m. It is the thickness of the anode catalyst layer.
        Hccl::Float64             = 7.605652607044295e-6   # m. It is the thickness of the cathode catalyst layer.
        # Membrane
        Hmem::Float64             = 13.646579963107156e-6  # m. It is the thickness of the membrane.
        # Gas diffusion layer
        Hgdl::Float64             = 121.28496643671034e-6  # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl::Float64      = 0.8436478459989776     # .  It is the anode/cathode GDL porosity.
        epsilon_c::Float64        = 0.2                    # .  It is the compression ratio of the GDL.
        #   Microporous layer
        Hmpl::Float64             = 43.98306893354156e-6   # m. It is the thickness of the microporous layer.
        epsilon_mpl::Float64      = 0.425                  # .  It is the porosity of the microporous layer.
        # Gas channel
        Hagc::Float64             = 230e-6                 # m. It is the thickness of the anode gas channel.
        Hcgc::Float64             = 300e-6                 # m. It is the thickness of the cathode gas channel.
        Wagc::Float64             = 430e-6                 # m. It is the width of the anode gas channel.
        Wcgc::Float64             = 532e-6                 # m. It is the width of the cathode gas channel.
        Lgc::Float64              = 246.2e-3               # m. It is the length of one channel in the bipolar plate.
        nb_channel_in_gc::Int64 = 105                      # .  It is the number of channels in the bipolar plate.
        Ldist::Float64            = 71.1e-3                # m. It is the length of the distributor, which is the volume between the gas channel and the manifold.
        #   Auxiliaries
        Lm::Float64               = 25.8e-3                # m.  It is the length of the manifold.
        A_T_a::Float64            = 9.01e-4                # m². It is the inlet/exhaust anode manifold throttle area.
        A_T_c::Float64            = 22.61e-4               # m². It is the inlet/exhaust cathode manifold throttle area.
        Vasm::Float64             = Lm * A_T_a             # m³. It is the supply manifold volume.
        Vcsm::Float64             = Lm * A_T_c             # m³. It is the supply manifold volume.
        Vaem::Float64             = Vasm                   # m³. It is the exhaust manifold volume.
        Vcem::Float64             = Vcsm                   # m³. It is the exhaust manifold volume.
        # Interaction parameters between fluids and PEMFC structure
        e::Int64                  = 3                      # .  It is the capillary exponent.
        K_l_ads::Float64          = 1.0                    # .  It is an estimation of the ratio between the liquid and vapor sorption rates of water in the membrane. It should be in [10-1000] [shaoNewInsightsSteadystate2023].
        K_O2_ad_Pt::Float64       = 7.346634385810734      # .  It is the interfacial resistance coefficient of O2 adsorption on the Pt sites.
        # Voltage polarization
        Re::Float64               = 1.545654084145453e-7   # Ω.m². It is the electron conduction resistance of the circuit.
        i0_c_ref::Float64         = 15.0                   # A.m-2. It is the dry reference exchange current density at the cathode.
        kappa_co::Float64         = 21.423681082096856     # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c::Float64          = 0.253020870903792      # .      It is the overpotential correction exponent.
        C_scl::Float64            = 2e7                    # F.m-3. It is the volumetric space-charge layer capacitance.

    # For EH-31 fuel cell
    elseif type_fuel_cell == "EH-31_1.5" || type_fuel_cell == "EH-31_2.0" ||
            type_fuel_cell == "EH-31_2.25" || type_fuel_cell == "EH-31_2.5"
        # Global
        Aact                      = 85e-4                  # m². It is the active area of the catalyst layer.
        nb_cell                   = 1                      # .   It is the number of cells in the stack.
        # Catalyst layer
        Hacl                      = 8.593e-6               # m. It is the thickness of the anode catalyst layer.
        Hccl                      = Hacl                   # m. It is the thickness of the cathode catalyst layer.
        # Membrane
        Hmem                      = 16.06e-6               # m. It is the thickness of the membrane.
        # Gas diffusion layer
        Hgdl                      = 200e-6                 # m. It is the thickness of the gas diffusion layer.
        epsilon_gdl               = 0.5002                 # .  It is the anode/cathode GDL porosity.
        epsilon_c                 = 0.2                    # .  It is the compression ratio of the GDL.
        #   Microporous layer
        Hmpl                      = 30e-6                  # m. It is the thickness of the microporous layer.
        epsilon_mpl               = 0.4                    # .  It is the porosity of the microporous layer.
        # Gas channel
        Hagc                      = 500e-6                 # m. It is the thickness of the anode gas channel.
        Hcgc                      = Hagc                   # m. It is the thickness of the cathode gas channel.
        Wagc                      = 450e-6                 # m. It is the width of the anode gas channel.
        Wcgc                      = Wagc                   # m. It is the width of the cathode gas channel.
        Lgc                       = 144e-3                 # m. It is the length of one channel in the bipolar plate.
        nb_channel_in_gc          = 67                     # .  It is the number of channels in the bipolar plate.
        Ldist                     = 5e-2                   # m. It is the estimated length of the distributor, which is the volume between the gas channel and the manifold.
        #   Auxiliaries
        Lm                        = 2.03                   # m.  It is the length of the manifold.
        A_T_a                     = 11.8e-4                # m². It is the inlet/exhaust anode manifold throttle area.
        A_T_c                     = 34.4e-4                # m². It is the inlet/exhaust cathode manifold throttle area.
        Vasm                      = Lm * A_T_a             # m³. It is the supply manifold volume.
        Vcsm                      = Lm * A_T_c             # m³. It is the supply manifold volume.
        Vaem                      = Vasm                   # m³. It is the exhaust manifold volume.
        Vcem                      = Vcsm                   # m³. It is the exhaust manifold volume.
        # Interaction parameters between fluids and PEMFC structure
        e                         = 4                      # .  It is the capillary exponent.
        K_l_ads                   = 1.0                    # .  It is an estimation of the ratio between the liquid and vapor sorption rates of water in the membrane. It should be in [10-1000] [shaoNewInsightsSteadystate2023].
        K_O2_ad_Pt                = 5.4                    # .  It is the interfacial resistance coefficient of O2 adsorption on the Pt sites.
        # Voltage polarization
        Re                        = 1e-6                   # Ω.m². It is the electron conduction resistance of the circuit.
        i0_c_ref                  = 14.43                  # A.m-2. It is the dry reference exchange current density at the cathode.
        kappa_co                  = 30.42                  # mol.m-1.s-1.Pa-1. It is the crossover correction coefficient.
        kappa_c                   = 0.4152                 # .      It is the overpotential correction exponent.
        C_scl                     = 20e6                   # F.m-3. It is the volumetric space-charge layer capacitance.

    # For other fuel cells
    else
        error("the type_fuel_cell given is not valid.")
    end

    return (Hacl, Hccl, Hmem, Hgdl, epsilon_gdl, epsilon_c, Hmpl, epsilon_mpl, Hagc, Hcgc, Wagc, Wcgc, Lgc,
            nb_channel_in_gc, Ldist, Lm, A_T_a, A_T_c, Vasm, Vcsm, Vaem, Vcem, Aact, nb_cell, e, K_l_ads,
            K_O2_ad_Pt, Re, i0_c_ref, kappa_co, kappa_c, C_scl)
end

