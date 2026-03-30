# src/alphapem/fuelcell/zsw.jl

"""
    ZSW

ZSW fuel cell model with parameters and experimental data.

Contains:
- physical, operating, numerical parameters (inherited from FuelCell)
- experimental polarization data (i_exp, U_exp)
"""
struct ZSW <: AbstractFuelCell
    physical_parameters::PhysicalParams
    operating_conditions::OperatingConditions
    experimental_values::ExperimentalValues
    numerical_parameters::NumericalParams
end

function physical_params(fc::ZSW)
    fc.physical_parameters = PhysicalParams(
        # Global
        Aact = 283.87e-4,                    # Active area of the catalyst layer in m²
        nb_cell = 26,                        # Number of cells in the stack
        # Catalyst layer
        Hacl = 8.112569325675836e-6,         # Thickness of the anode catalyst layer in meters
        Hccl = 7.605652607044295e-6,         # Thickness of the cathode catalyst layer in meters
        # Membrane
        Hmem = 13.646579963107156e-6,        # Thickness of the membrane in meters
        # Gas diffusion layer
        Hgdl = 121.28496643671034e-6,        # Thickness of the gas diffusion layer in meters
        epsilon_gdl = 0.8436478459989776,    # Anode/cathode GDL porosity
        epsilon_c = 0.2,                     # Compression ratio of the GDL
        #   Microporous layer
        Hmpl = 43.98306893354156e-6,         # Thickness of the microporous layer in meters
        epsilon_mpl = 0.425,                 # Porosity of the microporous layer
        # Gas channel
        Hagc = 230e-6,                       # Thickness of the anode gas channel in meters
        Hcgc = 300e-6,                       # Thickness of the cathode gas channel in meters
        Wagc = 430e-6,                       # Width of the anode gas channel in meters
        Wcgc = 532e-6,                       # Width of the cathode gas channel in meters
        Lgc = 246.2e-3,                      # Length of the gas channel in meters
        nb_channel_in_gc = 105,              # Number of channels in the bipolar plate
        Ldist = 71.1e-3,                     # Length of the distributor (between gas channel and manifold) in meters
        #   Auxiliaries
        Lm = 25.8e-3,                        # Length of the manifold in meters
        A_T_a = 9.01e-4,                     # Inlet/exhaust anode manifold throttle area in m²
        A_T_c = 22.61e-4,                    # Inlet/exhaust cathode manifold throttle area in m²
        Vasm = 25.8e-3 * 9.01e-4,            # Supply manifold volume at the anode in m³
        Vcsm = 25.8e-3 * 22.61e-4,           # Supply manifold volume at the cathode in m³
        Vaem = 25.8e-3 * 9.01e-4,            # Exhaust manifold volume at the anode in m³
        Vcem = 25.8e-3 * 22.61e-4,           # Exhaust manifold volume at the cathode in m³
        # Interaction parameters between fluids and PEMFC structure
        e = 3,                               # Capillary exponent
        K_l_ads = 1.0,                       # Ratio between the liquid and vapor sorption rates of water in the membrane
        K_O2_ad_Pt = 7.346634385810734,      # Interfacial resistance coefficient of O2 adsorption on the Pt sites
        # Voltage polarization
        Re = 1.545654084145453e-7,           # Electron conduction resistance of the circuit in Ω·m²
        i0_c_ref = 15.0,                     # Reference exchange current density at the cathode in A·m⁻²
        kappa_co = 21.423681082096856,       # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
        kappa_c = 0.253020870903792,         # Overpotential correction exponent
        C_scl = 2e7                          # Volumetric space-charge layer capacitance in F·m⁻³
    )
    return fc.physical_parameters
end


function operating_conditions(fc::ZSW, type_fuel_cell::String)
    if type_fuel_cell == "ZSW-GenStack"
        T_des::Float64          = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des::Float64         = 2.2e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des::Float64         = 2.0e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa::Float64             = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc::Float64             = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des::Float64      = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des::Float64      = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in::Float64        = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 1.61e5         # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 1.41e5         # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                = 2.01e5         # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 1.81e5         # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.4e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 2.2e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
        T_des                   = 68.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.8e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 2.6e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == "ZSW-GenStack_T_62"
        T_des                   = 62.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.2e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == "ZSW-GenStack_T_76"
        T_des                   = 76.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.2e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == "ZSW-GenStack_T_84"
        T_des                   = 84.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.2e5          # Pa. It is the desired pressure of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressure of the fuel gas at the cathode.
        Sa                      = 1.6            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 1.6            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.398          # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.50           # It is the desired relative humidity at the cathode.
        y_H2_in                 = 0.7            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    else
        error("Unknown type_fuel_cell: $type_fuel_cell")
    end

    fc.operating_conditions = OperatingConditions(T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in)
    return fc.operating_conditions
end


function experimental_values(fc::ZSW, type_fuel_cell::String, voltage_zone::String)
    if type_fuel_cell == "ZSW-GenStack"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.200, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp = [0.953, 0.864, 0.838, 0.819, 0.804, 0.778, 0.760, 0.743, 0.721, 0.703,
                       0.694, 0.685, 0.666, 0.644, 0.620, 0.593, 0.579, 0.546, 0.486]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.200, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700]
            U_exp = [0.953, 0.864, 0.838, 0.819, 0.804, 0.778, 0.760, 0.743, 0.721, 0.703,
                       0.694, 0.685, 0.666, 0.644, 0.620]
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_1.61_Pc_1.41"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700, 0.900,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.171]
            U_exp = [0.936, 0.835, 0.809, 0.795, 0.783, 0.759, 0.741, 0.725, 0.701, 0.670,
                       0.661, 0.633, 0.587, 0.541, 0.500, 0.457, 0.437, 0.402]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700]
            U_exp = [0.936, 0.835, 0.809, 0.795, 0.783, 0.759, 0.741, 0.725, 0.701]
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.01_Pc_1.81"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.415]
            U_exp = [0.946, 0.855, 0.830, 0.811, 0.795, 0.770, 0.752, 0.736, 0.717, 0.697,
                       0.685, 0.677, 0.655, 0.629, 0.599, 0.564, 0.545, 0.502, 0.450]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300]
            U_exp = [0.946, 0.855, 0.830, 0.811, 0.795, 0.770, 0.752, 0.736, 0.717, 0.697,
                       0.685, 0.677, 0.655]
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.4_Pc_2.2"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp = [0.949, 0.867, 0.841, 0.821, 0.807, 0.781, 0.763, 0.746, 0.725, 0.706,
                       0.697, 0.687, 0.670, 0.651, 0.630, 0.607, 0.595, 0.566, 0.514]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.200, 0.300, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900]
            U_exp = [0.949, 0.867, 0.841, 0.821, 0.807, 0.781, 0.763, 0.746, 0.725, 0.706,
                       0.697, 0.687, 0.670, 0.651, 0.630, 0.607]
        end
    elseif type_fuel_cell == "ZSW-GenStack_Pa_2.8_Pc_2.6"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp = [0.947, 0.872, 0.846, 0.827, 0.812, 0.787, 0.768, 0.752, 0.731, 0.711,
                       0.703, 0.694, 0.676, 0.659, 0.641, 0.622, 0.610, 0.588, 0.547]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900]
            U_exp = [0.947, 0.872, 0.846, 0.827, 0.812, 0.787, 0.768, 0.752, 0.731, 0.711,
                       0.703, 0.694, 0.676, 0.659, 0.641, 0.622]
        end
    elseif type_fuel_cell == "ZSW-GenStack_T_62"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp = [0.944, 0.855, 0.827, 0.808, 0.795, 0.771, 0.754, 0.739, 0.717, 0.696,
                       0.685, 0.675, 0.653, 0.631, 0.606, 0.581, 0.566, 0.532, 0.471]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500]
            U_exp = [0.944, 0.855, 0.827, 0.808, 0.795, 0.771, 0.754, 0.739, 0.717, 0.696,
                       0.685, 0.675, 0.653, 0.631]
        end
    elseif type_fuel_cell == "ZSW-GenStack_T_76"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.900,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000, 2.200, 2.500]
            U_exp = [0.946, 0.849, 0.825, 0.811, 0.799, 0.776, 0.759, 0.744, 0.724, 0.702,
                       0.691, 0.679, 0.652, 0.621, 0.587, 0.547, 0.527, 0.482, 0.406]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.299, 0.400, 0.498, 0.700, 0.900,
                       0.999, 1.099]
            U_exp = [0.946, 0.849, 0.825, 0.811, 0.799, 0.776, 0.759, 0.744, 0.724, 0.702,
                       0.691, 0.679]
        end
    elseif type_fuel_cell == "ZSW-GenStack_T_84"
        if voltage_zone == "full"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.300, 0.400, 0.498, 0.700, 0.901,
                       0.999, 1.099, 1.300, 1.500, 1.700, 1.900, 2.000]
            U_exp = [0.930, 0.847, 0.820, 0.805, 0.794, 0.772, 0.756, 0.741, 0.718, 0.686,
                       0.668, 0.650, 0.614, 0.575, 0.532, 0.486, 0.461]
        elseif voltage_zone == "before_voltage_drop"
            i_exp = [0.001, 0.050, 0.099, 0.150, 0.199, 0.300, 0.400, 0.498, 0.700]
            U_exp = [0.930, 0.847, 0.820, 0.805, 0.794, 0.772, 0.756, 0.741, 0.718]
        end
    else
        error("Unknown type_fuel_cell: $type_fuel_cell")
    end
    fc.experimental_values = ExperimentalValues(i_exp = i_exp, U_exp = U_exp)
    return fc.experimental_values
end


function numerical_params(fc::ZSW)
    fc.numerical_parameters = NumericalParams(
        # Setting the number of model points placed inside each layer:
        nb_gc = 1,                             # Number of model nodes placed inside each gas channel
        nb_gdl = 3,                            # Number of model nodes placed inside each GDL
        nb_mpl = 2,                            # Number of model nodes placed inside each MPL
        # Setting the purging parameters of the system and the dynamic display of the step current density function:
        purge_time = 0.6,                      # The time for purging the system in seconds
        delta_purge = 15.0,                    # The time between two purges in seconds
        # Setting the tolerances for the system of ODEs solver:
        rtol = 1e-6,                           # Relative tolerance for the system of ODEs solver
        atol = 1e-9                            # Absolute tolerance for the system of ODEs solver
    )
    return fc.numerical_parameters
end