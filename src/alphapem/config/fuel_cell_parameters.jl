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
    # ════════════════════════════════════════════════════════════════════════════
    # Geometrical parameters of the fuel cell
    # ════════════════════════════════════════════════════════════════════════════
    # Global
    Aact::Float64         = 300e-4           # MEA active area in meter squares.
    nb_cell::Int64        = 1                # Number of cells in the stack.
    #   Catalyst layer
    Hacl::Float64         = 10e-6            # Thickness of the anode catalyst layer in meters
    Hccl::Float64         = 10e-6            # Thickness of the cathode catalyst layer in meters
    #   Membrane
    Hmem::Float64         = 20e-6            # Thickness of the membrane in meters
    #   Gas diffusion layer
    Hgdl::Float64         = 200e-6           # Thickness of the gas diffusion layer in meters
    epsilon_gdl::Float64  = 0.7              # Anode/cathode GDL porosity
    #   Microporous layer
    Hmpl::Float64         = 30e-6            # Thickness of the microporous layer in meters
    epsilon_mpl::Float64  = 0.4              # Porosity of the microporous layer
    #   Gas channel
    Hagc::Float64         = 500e-6           # Thickness of the anode gas channel in meters
    Hcgc::Float64         = 500e-6           # Thickness of the cathode gas channel in meters
    Wagc::Float64         = 450e-6           # Width of the anode gas channel in meters
    Wcgc::Float64         = 450e-6           # Width of the cathode gas channel in meters
    Lgc::Float64          = 0.144            # Length of the gas channel in meters
    nb_channel_in_gc::Int = 67               # Number of channels in the bipolar plate
    Ldist::Float64        = 0.05             # Length of the distributor (between gas channel and manifold) in meters
    #   Auxiliaries
    Lm::Float64           = 0.025            # Length of the manifold in meters
    L_endplate::Float64   = 45e-3            # Length of the endplate in meters.
    A_T_a::Float64        = 12e-4            # Inlet/exhaust anode manifold throttle area in m²
    A_T_c::Float64        = 12e-4            # Inlet/exhaust cathode manifold throttle area in m²
    Vasm::Float64         = 7e-3             # Supply manifold volume at the anode in m³
    Vcsm::Float64         = 7e-3             # Supply manifold volume at the cathode in m³
    Vaem::Float64         = 2.4e-3           # Exhaust manifold volume at the anode in m³
    Vcem::Float64         = 2.4e-3           # Exhaust manifold volume at the cathode in m³

    # ════════════════════════════════════════════════════════════════════════════
    # Physical parameters of the model
    # ════════════════════════════════════════════════════════════════════════════
    # Model parameters for fluidic calculation
    rho_mem::Float64        = 1980           # kg.m-3. It is the density of the dry membrane.
    M_eq::Float64           = 1.1            # kg.mol-1. It is the equivalent molar mass of ionomer.
    tau_mpl::Float64        = 2              # It is the pore structure coefficient in the MPL, without units [Gen Inoue 2016 Journal Power Sources].
    tau_cl::Float64         = 4              # It is the pore structure coefficient in the CL, without units [Gen Inoue 2016 Journal Power Sources].
    r_s_gdl::Float64        = 2.0            # It is the exponent pore blockage in the GDL.
    r_s_mpl::Float64        = 2.5            # It is the exponent pore blockage in the MPL.
    r_s_cl::Float64         = 2.5            # It is the exponent pore blockage in the CL.
    Dp_gdl::Float64         = 33.2e-6        # m. It is the pore diameter of the GDL [ZSW GenStack].
    Dp_mpl::Float64         = 17.4e-6        # m. It is the pore diameter of the MPL [morganUnderstandingGasDiffusion2014].
    Dp_cl::Float64          = 0.15e-6        # m. It is the pore diameter of the CL [Ali Malekian 2019 International Journal of Hydrogen Energy].
    theta_c_gdl::Float64    = 120 * π / 180  # radian. It is the contact angle of GDL for liquid water.
    theta_c_mpl::Float64    = 135 * π / 180  # radian. It is the contact angle of MPL for liquid water.
    theta_c_cl::Float64     = 95 * π / 180   # radian. It is the contact angle of CL for liquid water.
    theta_l_rem::Float64    = 5e-5           # s/m. It is the coefficient of liquid water removal from the GDL to the GC [Ansys Fluent value from their User Guide].
    K_v_liq_gas::Float64    = 0.02           # . It is the liquid to gas velocity ratio in the GC [Ansys Fluent value from their User Guide].
    gamma_sorp_l::Float64   = 0.5            # s-1. Sorption rate of liquid water in the membrane [Ansys Fluent value from their User Guide].
    D_liq_dif::Float64      = 1e-5           # kg.m-1.s-1. It is the diffusion coefficient of liquid water in the GC [Ansys Fluent value from their User Guide].
    gamma_cond::Float64     = 1e8            # s-1. It is the overall condensation rate constant for water [Ansys Fluent value from their User Guide].
    gamma_evap::Float64     = 1e8            # s-1. It is the overall evaporation rate constant for water [Ansys Fluent value from their User Guide].
    epsilon_p::Float64      = 0.11           # . It is the percolation threshold porosity of the GDL.
    epsilon_c::Float64      = 0.2            # . It is the compression ratio of the GDL.
    alpha_p::Float64        = 0.785          # . It is a fitted value for the effective matter transfer in the GDL, for through plane direction.
    Eact_H2_cros_v::Float64 = 2.1e4          # J.mol-1. It is the activation energy of H2 for crossover in the under saturated membrane.
    Eact_H2_cros_l::Float64 = 1.8e4          # J.mol-1. It is the activation energy of H2 for crossover in the liquid-equilibrated membrane.
    Eact_O2_cros_v::Float64 = 2.2e4          # J.mol-1. It is the activation energy of oxygen for crossover in the under saturated membrane.
    Eact_O2_cros_l::Float64 = 2.0e4          # J.mol-1. It is the activation energy of oxygen for crossover in the liquid-equilibrated membrane.
    Kshape::Int64           = 2              # . Mathematical factor governing lambda_eq smoothing.
    # Interaction parameters between fluids and PEMFC structure
    e::Int64                = 3              # . Capillary exponent
    # Volumic flow of O2 inside the CCL to the Pt sites
    K_O2_ad_Pt::Float64    = 5.4             # Interfacial resistance coefficient of O2 adsorption on the Pt sites
    K_O2_dis_ion::Float64  = 8.5             # . It is the interfacial resistance coefficient of O2 dissolution inside the ionomer [haoModelingExperimentalValidation2015].
    K_O2_dis_l::Float64    = 1.0             # . It is the interfacial resistance coefficient of O2 dissolution inside the CL liquid water.
    IC::Float64            = 0.5             # . It is the ionomer to carbon ratio in the catalyst layer.
    ECSA_0::Float64        = 150             # cm2_Pt.cm-2_active_area. It is the initial electrochemical surface area of the catalyst.
    wt_Pt::Float64         = 0.5             # It is the weight fraction of platinum over carbon covered by platinum (Pt/C) in the cathode catalyst layer [haoModelingExperimentalValidation2015].
    L_Pt::Float64          = 0.3e-2          # kg.m-2. It is the platinum loading in the cathode catalyst layer.
    r_carb::Float64        = 40e-9           # m. It is the mean radius of the carbon particles.
    theta_Pt_0::Float64    = 0               # This is the initial platine-oxide coverage, assumed to be zero for simplification.
    rho_ion::Float64       = 1900            # kg.m-3. It is the density of the ionomer [haoModelingExperimentalValidation2015].
    #   Voltage calculation
    Re::Float64            = 1e-6            # Electron conduction resistance of the circuit in Ω·m²
    i0_c_ref::Float64      = 14.43           # Reference exchange current density at the cathode in A·m⁻²
    kappa_co::Float64      = 30.0            # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
    kappa_c::Float64       = 1.0             # Overpotential correction exponent
    C_scl::Float64         = 2e7             # Volumetric space-charge layer capacitance in F·m⁻³
    alpha_c::Float64       = 0.5             # It is the transfer coefficient of the cathode.

    # Model parameters for the heat transfer calculation
    #   Thermal conductivities
    k_th_gdl::Float64      = 0.3             # J.m-1.s-1.K-1. It is the thermal non-effective conductivity of the GDLs (non-effective ?) [ZSW].
    k_th_mpl::Float64      = 0.27            # J.m-1.s-1.K-1. It is the thermal conductivity of the MPLs [kotakaImpactInterfacialWater2014].
    k_th_cl::Float64       = 0.27            # J.m-1.s-1.K-1. It is the thermal conductivity of the CLs [vetterFreeOpenReference2019].
    k_th_mem::Float64      = 0.3             # J.m-1.s-1.K-1. It is the thermal conductivity of the membrane [vetterFreeOpenReference2019].
    #   Specific heat capacities
    Cp_gdl::Float64        = 568             # J.kg-1.K-1. It is the specific heat capacities of the GDLs [wangQuasi2DTransientModel2018].
    Cp_mpl::Float64        = 568             # J.kg-1.K-1. It is the specific heat capacities of the MPLs [yangEffectsOperatingConditions2019].
    Cp_cl::Float64         = 3300            # J.kg-1.K-1. It is the specific heat capacities the CLs [wangQuasi2DTransientModel2018].
    Cp_mem::Float64        = 833             # J.kg-1.K-1. It is the specific heat capacities of the membrane [wangQuasi2DTransientModel2018].
    #   Densities
    rho_gdl::Float64       = 1000            # kg.m-3. It is the density of the GDLs [wangQuasi2DTransientModel2018].
    rho_mpl::Float64       = 1000            # kg.m-3. It is the density of the MPLs [yangEffectsOperatingConditions2019].
    rho_cl::Float64        = 1000            # kg.m-3. It is the density of the CLs [wangQuasi2DTransientModel2018].
    #   Electrical conductivities
    sigma_e_gdl::Float64   = 1250            # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
    sigma_e_mpl::Float64   = 5000            # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [yangEffectsOperatingConditions2019].
    sigma_e_cl::Float64    = 350             # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
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
    :Hacl          => (1e-6, 15e-6, :real),                 # Anode catalyst-layer thickness
    :Hccl          => (5e-6, 20e-6, :real),                 # Cathode catalyst-layer thickness
    :Hmem          => (5e-6, 50e-6, :real),                 # Membrane thickness
    :Hgdl          => (70e-6, 150e-6, :real),               # Gas-diffusion-layer thickness
    :Hmpl          => (40e-6, 100e-6, :real),               # Microporous-layer thickness
    :epsilon_gdl   => (0.5, 0.9, :real),                    # GDL porosity
    :theta_c_gdl   => (90 * π / 180, 180 * π / 180, :real), # GDL contact angle
    :epsilon_mpl   => (0.3, 0.7, :real),                    # MPL porosity
    :IC            => (0.01, 2.0, :real),                   # Ionomer to carbon ratio in the catalyst layer
    :r_carb        => (10e-9, 100e-9, :real),               # Mean radius of the carbon particles
    :ECSA_0        => (10.0, 200.0, :real),                 # Initial electrochemical surface area of the catalyst
    :K_O2_dis_ion  => (0.1, 20.0, :real),                   # Interfacial resistance coefficient of O₂ dissolution inside the ionomer
    :K_O2_ad_Pt    => (0.1, 20.0, :real),                   # Interfacial resistance coefficient of O₂ adsorption on the Pt sites
    :alpha_c       => (0.01, 1.0, :real),                   # Cathode transfer coefficient
    :e             => (3, 5, :int),                         # Capillary exponent
    :Re            => (5e-8, 5e-6, :real),                  # Electron-conduction resistance
    :i0_c_ref      => (0.1, 1000.0, :real),                   # Reference cathode exchange current density
    :kappa_co      => (0.01, 40.0, :real),                  # Crossover correction coefficient
    :kappa_c       => (0.25, 4.0, :real),                   # Overpotential correction exponent
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
    :theta_c_gdl   => ("rad", "GDL contact angle"),
    :epsilon_mpl   => ("—", "MPL porosity"),
    :IC            => ("—", "Ionomer to carbon ratio in the catalyst layer"),
    :r_carb        => ("m", "Mean radius of the carbon particles"),
    :ECSA_0        => ("cm²_Pt.cm⁻²_active_area", "Initial electrochemical surface area of the catalyst"),
    :K_O2_dis_ion  => ("—", "O₂ dissolution interfacial resistance coefficient inside the ionomer"),
    :K_O2_ad_Pt    => ("—", "O₂ adsorption interfacial resistance coefficient on the Pt sites"),
    :alpha_c       => ("—", "Cathode transfer coefficient"),
    :e             => ("—", "Capillary exponent"),
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
    i_min_stoich::Float64 = 0.5     # Minimum current density used to compute the desired flows (A.cm-2)
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


