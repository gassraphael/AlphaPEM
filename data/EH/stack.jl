# EH-31 (2022) master data
#
# This file contains the generalisable, fuel-cell-specific information:
#   - physical parameters,
#   - nominal operating conditions,
#   - undetermined parameter bounds used for calibration.
#
# Experimental polarization curves are stored separately under pola/.

using AlphaPEM.Config: PhysicalParams, OperatingConditions

# Define local variables for parameters used in multiple places or calculations
const Hcl = 8.593e-6                     # Thickness of the catalyst layers in meters
const Hgc = 500e-6                       # Thickness of the gas channels in meters
const Wgc = 450e-6                       # Width of the gas channels in meters

# Manifold parameters
const Lm = 2.03                          # Length of the manifold in meters
const A_T_a = 11.8e-4                    # Inlet/exhaust anode manifold throttle area in m²
const A_T_c = 34.4e-4                    # Inlet/exhaust cathode manifold throttle area in m²
const Vasm = Lm * A_T_a                  # Supply manifold volume at the anode in m³
const Vcsm = Lm * A_T_c                  # Supply manifold volume at the cathode in m³
const Vaem = Vasm                        # Exhaust manifold volume at the anode in m³
const Vcem = Vcsm                        # Exhaust manifold volume at the cathode in m³

const PHYSICAL_PARAMETERS = PhysicalParams(
    # Global
    Aact = 85e-4,                        # Active area of the catalyst layer in m²
    nb_cell = 1,                         # Number of cells in the stack
    # Catalyst layer
    Hacl = Hcl,                          # Thickness of the anode catalyst layer in meters
    Hccl = Hcl,                          # Thickness of the cathode catalyst layer in meters
    # Membrane
    Hmem = 16.06e-6,                     # Thickness of the membrane in meters
    # Gas diffusion layer
    Hgdl = 200e-6,                       # Thickness of the gas diffusion layer in meters
    epsilon_gdl = 0.5002,                # Anode/cathode GDL porosity
    # Microporous layer
    Hmpl = 30e-6,                        # Thickness of the microporous layer in meters
    epsilon_mpl = 0.4,                   # Porosity of the microporous layer
    # Gas channel
    Hagc = Hgc,                          # Thickness of the anode gas channel in meters
    Hcgc = Hgc,                          # Thickness of the cathode gas channel in meters
    Wagc = Wgc,                          # Width of the anode gas channel in meters
    Wcgc = Wgc,                          # Width of the cathode gas channel in meters
    Lgc = 144e-3,                        # Length of the gas channel in meters
    nb_channel_in_gc = 67,               # Number of channels in the bipolar plate
    Ldist = 5e-2,                        # Length of the distributor (between gas channel and manifold) in meters
    # Auxiliaries
    Lm = Lm,                             # Length of the manifold in meters
    A_T_a = A_T_a,                       # Inlet/exhaust anode manifold throttle area in m²
    A_T_c = A_T_c,                       # Inlet/exhaust cathode manifold throttle area in m²
    Vasm = Vasm,                         # Supply manifold volume at the anode in m³
    Vcsm = Vcsm,                         # Supply manifold volume at the cathode in m³
    Vaem = Vaem,                         # Exhaust manifold volume at the anode in m³
    Vcem = Vcem,                         # Exhaust manifold volume at the cathode in m³
    # Interaction parameters between fluids and PEMFC structure
    K_O2_ad_Pt = 5.4,                    # Interfacial resistance coefficient of O2 adsorption on the Pt sites
    # Voltage polarization
    Re = 1e-6,                           # Electron conduction resistance of the circuit in Ω·m²
    i0_c_ref = 14.43,                    # Reference exchange current density at the cathode in A·m⁻²
    kappa_co = 1,                        # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
    kappa_c = 0.4152,                    # Overpotential correction exponent
    C_scl = 2e7                          # Volumetric space-charge layer capacitance in F·m⁻³
)

const OPERATING_CONDITIONS = OperatingConditions(
    T_des        = 74.0 + 273.15,  # K. Desired fuel cell temperature.
    Pa_des       = 1.5e5,          # Pa. Desired pressure of the fuel gas at the anode.
    Pc_des       = 1.5e5,          # Pa. Desired pressure of the fuel gas at the cathode.
    Sa           = 1.2,            # Stoichiometric ratio of hydrogen at the anode.
    Sc           = 2.0,            # Stoichiometric ratio of oxygen at the cathode.
    Phi_a_des    = 0.4,            # Desired relative humidity at the anode.
    Phi_c_des    = 0.6,            # Desired relative humidity at the cathode.
    y_H2_in      = 1.0,            # Molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet.
    i_min_stoich = 0.5             # A/cm². Minimum current density used to compute the desired flows.
)

# Undetermined parameters for calibration: (name, min_bound, max_bound)
const UNDETERMINED_PARAMETERS_BEFORE_VOLTAGE_DROP = [
    (:Hacl,        1e-6, 15e-6),   # Anode/cathode catalyst-layer thickness (Hccl = Hacl)
    (:Hmem,        5e-6, 50e-6),   # Membrane thickness
    (:epsilon_gdl, 0.5, 0.9),      # GDL porosity
    (:Re,          5e-8, 5e-6),    # Electron-conduction resistance
    (:i0_c_ref,    0.1, 100.0),    # Reference cathode exchange current density
    (:kappa_co,    0.01, 40.0),    # Crossover correction coefficient
    (:kappa_c,     0.25, 4.0),     # Overpotential correction exponent
]

const UNDETERMINED_PARAMETERS_AFTER_VOLTAGE_DROP = [
    (:K_O2_ad_Pt,  0.1, 20.0),     # O₂ adsorption resistance coefficient
]

const UNDETERMINED_PARAMETERS_FULL =
    vcat(UNDETERMINED_PARAMETERS_BEFORE_VOLTAGE_DROP, UNDETERMINED_PARAMETERS_AFTER_VOLTAGE_DROP)
