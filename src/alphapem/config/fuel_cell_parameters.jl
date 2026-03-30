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
    K_l_ads::Float64 = 1.0          # Ratio between the liquid and vapor sorption rates of water in the membrane
    K_O2_ad_Pt::Float64 = 5.4       # Interfacial resistance coefficient of O2 adsorption on the Pt sites
    #   Voltage polarization
    Re::Float64 = 1e-6              # Electron conduction resistance of the circuit in Ω·m²
    i0_c_ref::Float64 = 14.43       # Reference exchange current density at the cathode in A·m⁻²
    kappa_co::Float64 = 30.0        # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
    kappa_c::Float64 = 1.0          # Overpotential correction exponent
    C_scl::Float64 = 2e7            # Volumetric space-charge layer capacitance in F·m⁻³
end


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
# EXPERIMENTAL VALUES
# ============================================================

"""
    ExperimentalValues

Structure to store experimental current density and voltage values for a fuel cell.
- `i_exp::Vector{Float64}`: Experimental current density values (A/m²).
- `U_exp::Vector{Float64}`: Experimental cell voltage values (V).
By default, both vectors are initialized as empty.
"""
Base.@kwdef struct ExperimentalValues <: AbstractFuelCellParams
    i_exp::Vector{Float64} = Float64[] # Experimental current density values (A/m²)
    U_exp::Vector{Float64} = Float64[] # Experimental cell voltage values (V)
end


# ============================================================
# NUMERICAL PARAMETERS
# ============================================================

"""
    NumericalParams

This structure contains the numerical parameters for the fuel cell simulation.
"""
Base.@kwdef struct NumericalParams <: AbstractFuelCellParams
    nb_gc::Int = 1                  # Number of model nodes placed inside each gas channel
    nb_gdl::Int = 3                 # Number of model nodes placed inside each GDL
    nb_mpl::Int = 2                 # Number of model nodes placed inside each MPL
    purge_time::Float64 = 0.6       # The time for purging the system in seconds
    delta_purge::Float64 = 15.0     # The time between two purges in seconds
    delta_t_dyn_step::Float64 = 0.1 # Time for dynamic display of the step current density function in seconds
    rtol::Float64 = 1e-6            # Relative tolerance for the system of ODEs solver
    atol::Float64 = 1e-9            # Absolute tolerance for the system of ODEs solver
end
