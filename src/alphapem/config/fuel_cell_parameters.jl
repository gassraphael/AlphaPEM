# -*- coding: utf-8 -*-

"""This file contains the current density parameter structures.
Default values are provided for a typical PEMFC, but they can be modified by the user
when creating an instance of the structure.
The parameters are organized into three main categories: physical parameters, operating conditions,
and numerical parameters. Each category is defined as a separate structure that inherits from the
abstract type `AbstractFuelCellParams`. This allows for better organization and modularity in the code,
making it easier to manage and update the parameters as needed."""

# ______________________________________________Current density parameters______________________________________________

abstract type AbstractFuelCellParams end

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================

"""
    PhysicalParams

This structure contains the physical and geometrical parameters of the fuel cell.
"""
Base.@kwdef struct PhysicalParams <: AbstractFuelCellParams
    # Global
    Aact::Float64 = 300e-4          # MEA active area in meter squares.
    nb_cell::Int64 = 1              # Number of cells in the stack.
    #   Catalyst layer
    Hacl::Float64 = 10e-6           # Thickness of the anode catalyst layer in meters
    Hccl::Float64 = 10e-6           # Thickness of the cathode catalyst layer in meters
    #   Membrane
    Hmem::Float64 = 20e-6           # Thickness of the membrane in meters
    #   Gas diffusion layer
    Hgdl::Float64 = 200e-6          # Thickness of the gas diffusion layer in meters
    epsilon_gdl::Float64 = 0.7      # Anode/cathode GDL porosity
    epsilon_c::Float64 = 0.2        # Compression ratio of the GDL
    #   Microporous layer
    Hmpl::Float64 = 30e-6           # Thickness of the microporous layer in meters
    epsilon_mpl::Float64 = 0.4      # Porosity of the microporous layer
    #   Gas channel
    Hagc::Float64 = 500e-6          # Thickness of the anode gas channel in meters
    Hcgc::Float64 = 500e-6          # Thickness of the cathode gas channel in meters
    Wagc::Float64 = 450e-6          # Width of the anode gas channel in meters
    Wcgc::Float64 = 450e-6          # Width of the cathode gas channel in meters
    Lgc::Float64 = 0.144            # Length of the gas channel in meters
    nb_channel_in_gc::Int = 67      # Number of channels in the bipolar plate
    Ldist::Float64 = 0.05           # Length of the distributor (between gas channel and manifold) in meters
    #   Auxiliaries
    Lm::Float64 = 0.025             # Length of the manifold in meters
    L_endplate::Float64 = 45e-3     # Length of the endplate in meters.
    A_T_a::Float64 = 12e-4          # Inlet/exhaust anode manifold throttle area in m²
    A_T_c::Float64 = 12e-4          # Inlet/exhaust cathode manifold throttle area in m²
    Vasm::Float64 = 7e-3            # Supply manifold volume at the anode in m³
    Vcsm::Float64 = 7e-3            # Supply manifold volume at the cathode in m³
    Vaem::Float64 = 2.4e-3          # Exhaust manifold volume at the anode in m³
    Vcem::Float64 = 2.4e-3          # Exhaust manifold volume at the cathode in m³
    #   Interaction parameters between fluids and PEMFC structure
    e::Int64 = 5                    # Capillary exponent
    K_O2_ad_Pt::Float64 = 5.4       # Interfacial resistance coefficient of O2 adsorption on the Pt sites
    #   Voltage polarization
    Re::Float64 = 1e-6              # Electron conduction resistance of the circuit in Ω·m²
    i0_c_ref::Float64 = 14.43       # Reference exchange current density at the cathode in A·m⁻²
    kappa_co::Float64 = 30.0        # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
    kappa_c::Float64 = 1.0          # Overpotential correction exponent
    C_scl::Float64 = 2e7            # Volumetric space-charge layer capacitance in F·m⁻³
end


# ============================================================
# UNDETERMINED PARAMETER BOUNDS
# ============================================================

"""
    UNDETERMINED_PARAMETER_BOUNDS

Dictionary of default bounds for undetermined physical parameters.
These bounds serve as fallback values when no fuel-cell-specific bounds are available.

Each parameter is mapped to a tuple: (min::Float64, max::Float64, type::Symbol)
where type is either :real or :int.
"""
const UNDETERMINED_PARAMETER_BOUNDS = Dict{Symbol, Tuple{Float64, Float64, Symbol}}(
    :Hacl          => (5e-6, 20e-6, :real),           # Anode catalyst-layer thickness
    :Hccl          => (5e-6, 20e-6, :real),           # Cathode catalyst-layer thickness
    :Hmem          => (5e-6, 50e-6, :real),           # Membrane thickness
    :Hgdl          => (100e-6, 150e-6, :real),        # Gas-diffusion-layer thickness
    :Hmpl          => (40e-6, 100e-6, :real),         # Microporous-layer thickness
    :epsilon_gdl   => (0.5, 0.95, :real),             # GDL porosity
    :e             => (3.0, 5.0, :int),               # Capillary exponent
    :K_O2_ad_Pt    => (0.1, 10.0, :real),             # Resistance coefficient of O₂ adsorption on the Pt sites
    :Re            => (5e-8, 5e-6, :real),            # Electron-conduction resistance
    :i0_c_ref      => (0.1, 80.0, :real),             # Reference cathode exchange current density
    :kappa_co      => (0.01, 40.0, :real),            # Crossover correction coefficient
    :kappa_c       => (0.25, 4.0, :real),             # Overpotential correction exponent
)


"""
    PARAMETER_METADATA

Metadata for undetermined parameters: unit and description.
Used for displaying and exporting parameter bounds.
"""
const PARAMETER_METADATA = Dict{Symbol, Tuple{String, String}}(
    :Hacl          => ("m", "Anode catalyst-layer thickness"),
    :Hccl          => ("m", "Cathode catalyst-layer thickness"),
    :Hmem          => ("m", "Membrane thickness"),
    :Hgdl          => ("m", "Gas-diffusion-layer thickness"),
    :Hmpl          => ("m", "Microporous-layer thickness"),
    :epsilon_gdl   => ("—", "GDL porosity"),
    :e             => ("—", "Capillary exponent"),
    :K_O2_ad_Pt    => ("—", "O₂ adsorption resistance coefficient"),
    :Re            => ("Ω·m²", "Electron-conduction resistance"),
    :i0_c_ref      => ("A·m⁻²", "Reference cathode exchange current density"),
    :kappa_co      => ("—", "Crossover correction coefficient"),
    :kappa_c       => ("—", "Overpotential correction exponent"),
)



# ============================================================
# OPERATING PARAMETERS
# ============================================================

"""
    OperatingParams

This structure contains the operating parameters of the fuel cell system.
"""
Base.@kwdef struct OperatingConditions <: AbstractFuelCellParams
    T_des::Float64 = 74.0 + 273.15  # Desired fuel cell temperature in Kelvin
    Pa_des::Float64 = 2e5           # Desired anode pressure in Pascal
    Pc_des::Float64 = 2e5           # Desired cathode pressure in Pascal
    Sa::Float64 = 1.2               # Stoichiometric ratio of hydrogen
    Sc::Float64 = 2.0               # Stoichiometric ratio of oxygen
    Phi_a_des::Float64 = 0.4        # Desired anode relative humidity
    Phi_c_des::Float64 = 0.6        # Desired cathode relative humidity
    y_H2_in::Float64 = 1.0          # Molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet
end


# ============================================================
# EXPERIMENTAL DATA
# ============================================================

"""
    PolaExperimentalData

Structure for storing experimental polarization data of a fuel cell.
- `i_exp::Vector{Float64}`: Experimental current density values during polarization (A/m²).
- `U_exp::Vector{Float64}`: Experimental cell voltage values during polarization (V).
Both vectors are initialized as empty by default.
"""
Base.@kwdef struct PolaExperimentalData <: AbstractFuelCellParams
    i_exp::Vector{Float64} = Float64[] # Experimental current density values (A/m²)
    U_exp::Vector{Float64} = Float64[] # Experimental cell voltage values (V)
end


