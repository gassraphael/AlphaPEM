# src/alphapem/fuelcell/eh31.jl

"""
    EH31

EH31 fuel cell model with parameters and experimental data.

Contains:
- physical, operating, numerical parameters (inherited from FuelCell)
- experimental polarization data (i_exp, U_exp)
"""
mutable struct EH31FuelCell <: AbstractFuelCell
    physical_parameters::PhysicalParams
    operating_conditions::OperatingConditions
    pola_exp_data::PolaExperimentalData
    pola_exp_data_cali::PolaExperimentalData
    numerical_parameters::NumericalParams
end

# Simple constructor for EH31
function EH31FuelCell(type_fuel_cell::Symbol, voltage_zone::Symbol)
    # Create a temporary object with uninitialized fields
    fc = EH31FuelCell(
         PhysicalParams(),
         OperatingConditions(),
         PolaExperimentalData(),
         PolaExperimentalData(),
         NumericalParams()
    )
    fc.physical_parameters = eh31_physical_params()
    fc.operating_conditions = eh31_operating_conditions(type_fuel_cell)
    fc.pola_exp_data = eh31_pola_exp_data(type_fuel_cell, voltage_zone)
    fc.pola_exp_data_cali = eh31_pola_exp_data_calibration(type_fuel_cell, voltage_zone)
    fc.numerical_parameters = eh31_numerical_params()
    return fc
end

function eh31_physical_params()::PhysicalParams
    return PhysicalParams(
        # Global
        Aact = 85e-4,                        # Active area of the catalyst layer in m²
        nb_cell = 1,                         # Number of cells in the stack
        # Catalyst layer
        Hacl = 8.593e-6,                     # Thickness of the anode catalyst layer in meters
        Hccl = Hacl,                         # Thickness of the cathode catalyst layer in meters
        # Membrane
        Hmem = 16.06e-6,                     # Thickness of the membrane in meters
        # Gas diffusion layer
        Hgdl = 200e-6,                       # Thickness of the gas diffusion layer in meters
        epsilon_gdl = 0.5002,                # Anode/cathode GDL porosity
        epsilon_c = 0.2,                     # Compression ratio of the GDL
        #   Microporous layer
        Hmpl = 30e-6,                        # Thickness of the microporous layer in meters
        epsilon_mpl = 0.4,                   # Porosity of the microporous layer
        # Gas channel
        Hagc = 500e-6,                       # Thickness of the anode gas channel in meters
        Hcgc = Hagc,                         # Thickness of the cathode gas channel in meters
        Wagc = 450e-6,                       # Width of the anode gas channel in meters
        Wcgc = Wagc,                         # Width of the cathode gas channel in meters
        Lgc = 144e-3,                        # Length of the gas channel in meters
        nb_channel_in_gc = 67,               # Number of channels in the bipolar plate
        Ldist = 5e-2,                        # Length of the distributor (between gas channel and manifold) in meters
        #   Auxiliaries
        Lm = 2.03,                           # Length of the manifold in meters
        A_T_a = 11.8e-4,                     # Inlet/exhaust anode manifold throttle area in m²
        A_T_c = 34.4e-4,                     # Inlet/exhaust cathode manifold throttle area in m²
        Vasm = Lm * A_T_a,                   # Supply manifold volume at the anode in m³
        Vcsm = Lm * A_T_c,                   # Supply manifold volume at the cathode in m³
        Vaem = Vasm,                         # Exhaust manifold volume at the anode in m³
        Vcem = Vcsm,                         # Exhaust manifold volume at the cathode in m³
        # Interaction parameters between fluids and PEMFC structure
        e = 4,                               # Capillary exponent
        K_l_ads = 1.0,                       # Ratio between the liquid and vapor sorption rates of water in the membrane
        K_O2_ad_Pt = 5.4,                    # Interfacial resistance coefficient of O2 adsorption on the Pt sites
        # Voltage polarization
        Re = 1e-6,                           # Electron conduction resistance of the circuit in Ω·m²
        i0_c_ref = 14.43,                    # Reference exchange current density at the cathode in A·m⁻²
        kappa_co = 30.42,                    # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
        kappa_c = 0.4152,                    # Overpotential correction exponent
        C_scl = 2e7                          # Volumetric space-charge layer capacitance in F·m⁻³
    )
end


function eh31_operating_conditions(type_fuel_cell::Symbol)::OperatingConditions
    if type_fuel_cell == :EH_31_1_5
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 1.5e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 1.5e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == :EH_31_2_0
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.0e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 2.0e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == :EH_31_2_25
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.25e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 2.25e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    elseif type_fuel_cell == :EH_31_2_5
        T_des                   = 74.0 + 273.15  # K.  It is the desired fuel cell temperature.
        Pa_des                  = 2.5e5          # Pa. It is the desired pressures of the fuel gas at the anode.
        Pc_des                  = 2.5e5          # Pa. It is the desired pressures of the fuel gas at the cathode.
        Sa                      = 1.2            # It is the stoichiometric ratio of hydrogen at the anode.
        Sc                      = 2.0            # It is the stoichiometric ratio of oxygen at the cathode.
        Phi_a_des               = 0.4            # It is the desired relative humidity at the anode.
        Phi_c_des               = 0.6            # It is the desired relative humidity at the cathode.
        y_H2_in                 = 1.0            # It is the molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    else
        error("Unknown type_fuel_cell: $type_fuel_cell")
    end

    return OperatingConditions(T_des, Pa_des, Pc_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in)
end


function eh31_pola_exp_data(type_fuel_cell::Symbol, voltage_zone::Symbol)
    if type_fuel_cell == :EH_31_1_5  # at 1.5 bar
        if voltage_zone == :full
            i_exp_pola = [0.050, 0.068, 0.089, 0.110, 0.147, 0.185, 0.233, 0.293, 0.352, 0.395,
                          0.455, 0.510, 0.556, 0.620, 0.672, 0.738, 0.799, 0.850, 0.892, 0.942,
                          1.039, 1.139, 1.212, 1.269, 1.360, 1.432, 1.525, 1.604, 1.683, 1.765,
                          1.878, 1.966, 2.050, 2.109, 2.151, 2.188, 2.246]
            U_exp_pola = [0.900, 0.882, 0.865, 0.850, 0.834, 0.823, 0.811, 0.794, 0.781, 0.772,
                          0.761, 0.752, 0.745, 0.735, 0.728, 0.719, 0.712, 0.706, 0.700, 0.694,
                          0.681, 0.668, 0.660, 0.653, 0.641, 0.634, 0.622, 0.610, 0.599, 0.586,
                          0.570, 0.556, 0.540, 0.530, 0.521, 0.513, 0.500]
        elseif voltage_zone == :before_voltage_drop
            i_exp_pola = [0.050, 0.068, 0.089, 0.110, 0.147, 0.185, 0.233, 0.293, 0.352, 0.395,
                          0.455, 0.510, 0.556, 0.620, 0.672, 0.738, 0.799, 0.850, 0.892, 0.942,
                          1.039, 1.139, 1.212, 1.269, 1.360, 1.432, 1.525, 1.604, 1.683]
            U_exp_pola = [0.900, 0.882, 0.865, 0.850, 0.834, 0.823, 0.811, 0.794, 0.781, 0.772,
                          0.761, 0.752, 0.745, 0.735, 0.728, 0.719, 0.712, 0.706, 0.700, 0.694,
                          0.681, 0.668, 0.660, 0.653, 0.641, 0.634, 0.622, 0.610, 0.599]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end
    elseif type_fuel_cell == :EH_31_2_0  # at 2.0 bar
        if voltage_zone == :full
            i_exp_pola = [0.050, 0.057, 0.079, 0.106, 0.135, 0.171, 0.206, 0.242, 0.302, 0.346,
                          0.395, 0.434, 0.476, 0.531, 0.570, 0.623, 0.681, 0.731, 0.779, 0.822,
                          0.868, 0.930, 0.976, 1.031, 1.090, 1.134, 1.205, 1.242, 1.312, 1.358,
                          1.403, 1.453, 1.501, 1.569, 1.634, 1.725, 1.786, 1.857, 1.924, 1.979,
                          2.050, 2.125, 2.168, 2.214, 2.258, 2.308, 2.348, 2.413, 2.459]
            U_exp_pola = [0.900, 0.889, 0.874, 0.860, 0.853, 0.845, 0.837, 0.830, 0.817, 0.808,
                          0.800, 0.792, 0.786, 0.779, 0.772, 0.765, 0.759, 0.753, 0.747, 0.742,
                          0.737, 0.730, 0.726, 0.720, 0.714, 0.710, 0.702, 0.698, 0.690, 0.684,
                          0.679, 0.673, 0.668, 0.659, 0.651, 0.640, 0.631, 0.620, 0.608, 0.598,
                          0.586, 0.573, 0.565, 0.557, 0.548, 0.537, 0.528, 0.513, 0.502]
        elseif voltage_zone == :before_voltage_drop
            i_exp_pola = [0.050, 0.057, 0.079, 0.106, 0.135, 0.171, 0.206, 0.242, 0.302, 0.346,
                          0.395, 0.434, 0.476, 0.531, 0.570, 0.623, 0.681, 0.731, 0.779, 0.822,
                          0.868, 0.930, 0.976, 1.031, 1.090, 1.134, 1.205, 1.242]
            U_exp_pola = [0.900, 0.889, 0.874, 0.860, 0.853, 0.845, 0.837, 0.830, 0.817, 0.808,
                          0.800, 0.792, 0.786, 0.779, 0.772, 0.765, 0.759, 0.753, 0.747, 0.742,
                          0.737, 0.730, 0.726, 0.720, 0.714, 0.710, 0.702, 0.698]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end
    elseif type_fuel_cell == :EH_31_2_25  # at 2.25 bar
        if voltage_zone == :full
            i_exp_pola = [0.056, 0.095, 0.120, 0.138, 0.160, 0.183, 0.218, 0.248, 0.279, 0.315,
                          0.364, 0.409, 0.477, 0.536, 0.594, 0.641, 0.697, 0.748, 0.809, 0.866,
                          0.944, 1.011, 1.074, 1.142, 1.193, 1.252, 1.322, 1.381, 1.442, 1.496,
                          1.545, 1.599, 1.675, 1.746, 1.827, 1.868, 1.918, 2.004, 2.053, 2.114,
                          2.156, 2.209, 2.257, 2.310, 2.356, 2.403, 2.468, 2.513, 2.552, 2.600,
                          2.636, 2.679, 2.728, 2.794]
            U_exp_pola = [0.894, 0.882, 0.873, 0.867, 0.861, 0.854, 0.847, 0.840, 0.834, 0.827,
                          0.819, 0.812, 0.801, 0.793, 0.786, 0.781, 0.775, 0.771, 0.764, 0.759,
                          0.751, 0.746, 0.740, 0.734, 0.728, 0.723, 0.715, 0.709, 0.703, 0.698,
                          0.692, 0.686, 0.678, 0.670, 0.660, 0.654, 0.647, 0.635, 0.628, 0.618,
                          0.613, 0.604, 0.596, 0.587, 0.580, 0.570, 0.559, 0.551, 0.545, 0.536,
                          0.528, 0.520, 0.511, 0.497]
        elseif voltage_zone == :before_voltage_drop
            i_exp_pola = [0.056, 0.095, 0.120, 0.138, 0.160, 0.183, 0.218, 0.248, 0.279, 0.315,
                          0.364, 0.409, 0.477, 0.536, 0.594, 0.641, 0.697, 0.748, 0.809, 0.866,
                          0.944, 1.011, 1.074, 1.142, 1.193, 1.252, 1.322, 1.381, 1.442, 1.496,
                          1.545, 1.599, 1.675]
            U_exp_pola = [0.894, 0.882, 0.873, 0.867, 0.861, 0.854, 0.847, 0.840, 0.834, 0.827,
                          0.819, 0.812, 0.801, 0.793, 0.786, 0.781, 0.775, 0.771, 0.764, 0.759,
                          0.751, 0.746, 0.740, 0.734, 0.728, 0.723, 0.715, 0.709, 0.703, 0.698,
                          0.692, 0.686, 0.678]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end
    elseif type_fuel_cell == :EH_31_2_5  # at 2.5 bar
        if voltage_zone == :full
            i_exp_pola = [0.057, 0.070, 0.082, 0.101, 0.127, 0.145, 0.168, 0.200, 0.234, 0.267,
                          0.296, 0.331, 0.355, 0.388, 0.423, 0.467, 0.527, 0.577, 0.632, 0.685,
                          0.740, 0.789, 0.845, 0.898, 0.953, 1.030, 1.124, 1.192, 1.254, 1.314,
                          1.364, 1.434, 1.514, 1.587, 1.643, 1.707, 1.769, 1.826, 1.892, 1.972,
                          2.040, 2.124, 2.192, 2.265, 2.358, 2.429, 2.508, 2.572, 2.624, 2.691,
                          2.750, 2.822, 2.879, 2.918, 2.956, 2.988]
            U_exp_pola = [0.900, 0.892, 0.884, 0.875, 0.866, 0.861, 0.856, 0.850, 0.845, 0.840,
                          0.835, 0.829, 0.824, 0.820, 0.814, 0.807, 0.800, 0.793, 0.787, 0.783,
                          0.778, 0.775, 0.771, 0.767, 0.763, 0.758, 0.750, 0.744, 0.738, 0.732,
                          0.726, 0.719, 0.712, 0.703, 0.697, 0.691, 0.685, 0.679, 0.672, 0.663,
                          0.657, 0.648, 0.640, 0.632, 0.621, 0.610, 0.600, 0.591, 0.584, 0.575,
                          0.566, 0.555, 0.546, 0.537, 0.531, 0.524]
        elseif voltage_zone == :before_voltage_drop
            i_exp_pola = [0.057, 0.070, 0.082, 0.101, 0.127, 0.145, 0.168, 0.200, 0.234, 0.267,
                          0.296, 0.331, 0.355, 0.388, 0.423, 0.467, 0.527, 0.577, 0.632, 0.685,
                          0.740, 0.789, 0.845, 0.898, 0.953, 1.030, 1.124, 1.192, 1.254, 1.314,
                          1.364, 1.434, 1.514]
            U_exp_pola = [0.900, 0.892, 0.884, 0.875, 0.866, 0.861, 0.856, 0.850, 0.845, 0.840,
                          0.835, 0.829, 0.824, 0.820, 0.814, 0.807, 0.800, 0.793, 0.787, 0.783,
                          0.778, 0.775, 0.771, 0.767, 0.763, 0.758, 0.750, 0.744, 0.738, 0.732,
                          0.726, 0.719, 0.712]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end
    else
        error("Unknown type_fuel_cell: $type_fuel_cell")
    end

    return PolaExperimentalData(i_exp = i_exp_pola .* 1e4, U_exp = U_exp_pola)
end


"""
    eh31_pola_exp_data_calibration(type_fuel_cell, voltage_zone)

This function returns the experimental values of polarisation curves made on different fuel cells at different
operating conditions. The experimental values are used for calibrating the model and so are composed of a reduced
number of points compared to the polarization data function. These points are specifically chosen to be as few as
possible while still providing a good representation of the polarisation curve.
"""
function eh31_pola_exp_data_calibration(type_fuel_cell::Symbol, voltage_zone::Symbol)
    if type_fuel_cell == :EH_31_1_5  # at 1.5 bar
        if voltage_zone == :full
            i_exp_cali = [0.050, 0.110, 0.293, 1.039, 1.683, 1.966, 2.246]
            U_exp_cali = [0.900, 0.850, 0.794, 0.681, 0.599, 0.556, 0.500]
        elseif voltage_zone == :before_voltage_drop
            i_exp_cali = [0.050, 0.110, 0.293, 1.039, 1.683]
            U_exp_cali = [0.900, 0.850, 0.794, 0.681, 0.599]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end

    elseif type_fuel_cell == :EH_31_2_0  # at 2.0 bar
        if voltage_zone == :full
            i_exp_cali = [0.050, 0.106, 0.242, 0.681, 1.242, 1.501, 1.979, 2.459]
            U_exp_cali = [0.900, 0.860, 0.830, 0.759, 0.698, 0.668, 0.598, 0.502]
        elseif voltage_zone == :before_voltage_drop
            i_exp_cali = [0.050, 0.106, 0.242, 0.681, 1.242]
            U_exp_cali = [0.900, 0.860, 0.830, 0.759, 0.698]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end

    elseif type_fuel_cell == :EH_31_2_25  # at 2.25 bar
        if voltage_zone == :full
            i_exp_cali = [0.056, 0.183, 0.364, 1.011, 1.675, 1.918, 2.356, 2.794]
            U_exp_cali = [0.894, 0.854, 0.819, 0.746, 0.678, 0.647, 0.580, 0.497]
        elseif voltage_zone == :before_voltage_drop
            i_exp_cali = [0.056, 0.183, 0.364, 1.011, 1.675]
            U_exp_cali = [0.894, 0.854, 0.819, 0.746, 0.678]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end

    elseif type_fuel_cell == :EH_31_2_5  # at 2.5 bar
        if voltage_zone == :full
            i_exp_cali = [0.057, 0.127, 0.296, 0.527, 1.030, 1.514, 1.972, 2.358, 2.691, 2.988]
            U_exp_cali = [0.900, 0.866, 0.835, 0.800, 0.758, 0.712, 0.663, 0.621, 0.575, 0.524]
        elseif voltage_zone == :before_voltage_drop
            i_exp_cali = [0.057, 0.127, 0.296, 0.527, 1.030, 1.514]
            U_exp_cali = [0.900, 0.866, 0.835, 0.800, 0.758, 0.712]
        else
            throw(ArgumentError("The voltage_zone should be either :full or :before_voltage_drop."))
        end
    else
        throw(ArgumentError("Unknown type_fuel_cell: $type_fuel_cell"))
    end

    return PolaExperimentalData(i_exp = i_exp_cali .* 1e4, U_exp = U_exp_cali)
end


function eh31_numerical_params()
    return NumericalParams(
        # Setting the number of model points placed inside each layer:
        nb_gc = 1,                             # Number of model nodes placed inside each gas channel
        nb_gdl = 3,                            # Number of model nodes placed inside each GDL
        nb_mpl = 2,                            # Number of model nodes placed inside each MPL
        nb_man = 1,                            # Number of model nodes placed inside each manifold line
        # Setting the purging parameters of the system and the dynamic display of the step current density function:
        purge_time = 0.6,                      # The time for purging the system in seconds
        delta_purge = 15.0,                    # The time between two purges in seconds
        # Setting the tolerances for the system of ODEs solver:
        rtol = 1e-6,                           # Relative tolerance for the system of ODEs solver
        atol = 1e-9                            # Absolute tolerance for the system of ODEs solver
    )
end
