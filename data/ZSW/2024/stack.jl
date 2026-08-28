# ZSW GenStack (2024) master data
#
# This file contains the generalisable, fuel-cell-specific information:
#   - physical parameters,
#   - nominal operating conditions,
#   - undetermined parameter bounds used for calibration.
#
# Experimental polarization curves are stored separately under pola/.

using AlphaPEM.Config: PhysicalParams, OperatingConditions

const PHYSICAL_PARAMETERS = PhysicalParams(
    # Global
    Aact             = 283.87e-4,              # Active area of the catalyst layer in m²
    nb_cell          = 26,                     # Number of cells in the stack
    # Catalyst layer
    Hacl             = 3e-6,                   # Thickness of the anode catalyst layer in meters
    Hccl             = 10.179819442673712e-6,  # Thickness of the cathode catalyst layer in meters
    # Membrane
    Hmem             = 15e-6,                  # Thickness of the membrane in meters
    # Gas diffusion layer
    Hgdl             = 75.88436374521656e-6,   # Thickness of the gas diffusion layer in meters
    epsilon_gdl      = 0.8717693694526278,     # Anode/cathode GDL porosity
    # Microporous layer
    Hmpl             = 64.17561409938712e-6,   # Thickness of the microporous layer in meters
    epsilon_mpl      = 0.5163119718049662,     # Porosity of the microporous layer
    # Gas channel
    Hagc             = 230e-6,                 # Thickness of the anode gas channel in meters
    Hcgc             = 300e-6,                 # Thickness of the cathode gas channel in meters
    Wagc             = 430e-6,                 # Width of the anode gas channel in meters
    Wcgc             = 532e-6,                 # Width of the cathode gas channel in meters
    Lgc              = 246.2e-3,               # Length of the gas channel in meters
    nb_channel_in_gc = 105,                    # Number of channels in the bipolar plate
    Ldist            = 71.1e-3,                # Length of the distributor (between gas channel and manifold) in meters
    # Auxiliaries
    Lm               = 25.8e-3,                # Length of the manifold in meters
    A_T_a            = 9.01e-4,                # Inlet/exhaust anode manifold throttle area in m²
    A_T_c            = 22.61e-4,               # Inlet/exhaust cathode manifold throttle area in m²
    Vasm             = 25.8e-3 * 9.01e-4,      # Supply manifold volume at the anode in m³
    Vcsm             = 25.8e-3 * 22.61e-4,     # Supply manifold volume at the cathode in m³
    Vaem             = 25.8e-3 * 9.01e-4,      # Exhaust manifold volume at the anode in m³
    Vcem             = 25.8e-3 * 22.61e-4,     # Exhaust manifold volume at the cathode in m³
    # Interaction parameters between fluids and PEMFC structure
    e                = 3,                      # Capillary exponent
    # Volumic flow of O2 inside the CCL to the Pt sites
    IC_ccl           = 1.8766573409239173,     # Ionomer to carbon ratio in the cathode catalyst layer
    ECSA_0           = 90.18946251019264,      # Initial electrochemical surface area of the catalyst in cm2_Pt.cm-2_active_area
    r_carb           = 20.17782731935742e-9,   # Mean radius of the carbon particles in m
    # Voltage polarization
    Re               = 0,                      # Electron conduction resistance of the circuit in Ω·m²
    i0_c_ref         = 1.1319658074709191,     # Reference exchange current density at the cathode in A·m⁻²
    alpha_c          = 0.6552890241967357,     # Transfer coefficient of the cathode
    kappa_co         = 15.3655148859884,       # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
    kappa_c          = 0.25862840985329477,    # Overpotential correction exponent
    C_scl            = 2e7                     # Volumetric space-charge layer capacitance in F·m⁻³
)

const OPERATING_CONDITIONS = OperatingConditions(
    T_des        = 68.0 + 273.15,  # K. Desired cooling circuit temperature.
    Pa_des       = 2.2e5,          # Pa. Desired pressure of the fuel gas at the anode.
    Pc_des       = 2.0e5,          # Pa. Desired pressure of the fuel gas at the cathode.
    Sa           = 1.6,            # Stoichiometric ratio of hydrogen at the anode.
    Sc           = 1.6,            # Stoichiometric ratio of oxygen at the cathode.
    Phi_a_des    = 0.398,          # Desired relative humidity at the anode.
    Phi_c_des    = 0.50,           # Desired relative humidity at the cathode.
    y_H2_in      = 0.7,            # Molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    i_min_stoich = 0.5             # A/cm². Minimum current density used to compute the desired flows.
)

"""
    apply_operating_conditions_correction(oc::OperatingConditions) -> OperatingConditions

Apply stack-specific corrections to the operating conditions read from the data files.

For the ZSW GenStack, the cooling circuit is not modeled explicitly; the desired
stack temperature is therefore increased by 3 °C so that the simulated bipolar-plate
temperature matches the experimental cooling-circuit temperature.
"""
function apply_operating_conditions_correction(oc::OperatingConditions)::OperatingConditions
    return OperatingConditions(
        T_des       = oc.T_des + 3.0,
        Pa_des      = oc.Pa_des,
        Pc_des      = oc.Pc_des,
        Sa          = oc.Sa,
        Sc          = oc.Sc,
        Phi_a_des   = oc.Phi_a_des,
        Phi_c_des   = oc.Phi_c_des,
        y_H2_in     = oc.y_H2_in,
        i_min_stoich= oc.i_min_stoich,
    )
end

# Undetermined parameters for calibration: (name, min_bound, max_bound)
const UNDETERMINED_PARAMETERS_BEFORE_VOLTAGE_DROP = [
    (:Hccl,         10e-6, 15e-6),    # Cathode catalyst-layer thickness
    (:Hgdl,         70e-6, 88e-6),    # Gas-diffusion-layer thickness
    (:Hmpl,         60e-6, 79e-6),    # Microporous-layer thickness
    (:epsilon_gdl,  0.76, 0.88),      # GDL porosity
    (:epsilon_mpl,  0.41, 0.59),      # MPL porosity
    (:alpha_c,      0.62, 9.5),       # Cathode transfer coefficient
    (:e,            3, 5),            # Capillary exponent
    (:IC_ccl,       1.4, 1.9),        # Ionomer to carbon ratio in the cathode catalyst layer
    (:k_th_gdl,     0.15, 0.9),       # Thermal conductivity of the GDL
    (:wt_Pt_ccl,    0.5, 0.63),       # Weight fraction of platinum over carbon in the cathode catalyst layer
    (:L_Pt_ccl,     2.6e-3, 3.5e-3),  # Platinum loading in the cathode catalyst layer
    (:i0_c_ref,     1, 34.0),         # Reference cathode exchange current density
    (:kappa_co,     10, 20.0),        # Crossover correction coefficient
    (:kappa_c,      0.25, 1.3),       # Overpotential correction exponent
]

const UNDETERMINED_PARAMETERS_AFTER_VOLTAGE_DROP = [
    (:theta_c_cl,   90 * π / 180, 140 * π / 180), # CL contact angle
    (:r_carb,       10e-9, 21.55e-9),             # Mean radius of the carbon particles
    (:ECSA_0,       80.0, 190.0),                 # Initial electrochemical surface area of the catalyst
]

const UNDETERMINED_PARAMETERS_FULL =
    vcat(UNDETERMINED_PARAMETERS_BEFORE_VOLTAGE_DROP, UNDETERMINED_PARAMETERS_AFTER_VOLTAGE_DROP)
