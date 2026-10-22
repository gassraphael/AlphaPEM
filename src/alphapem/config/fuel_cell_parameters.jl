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
    # ── Global ──
    Aact::Float64         = 300e-4           # MEA active area in meter squares.
    nb_cell::Int64        = 1                # Number of cells in the stack.
    # ── Catalyst layer ──
    Hacl::Float64         = 5e-6             # Thickness of the anode catalyst layer in meters
    Hccl::Float64         = 15e-6            # Thickness of the cathode catalyst layer in meters
    # ── Membrane ──
    Hmem::Float64         = 20e-6            # Thickness of the membrane in meters
    # ── Gas diffusion layer ──
    Hgdl::Float64         = 125e-6           # Thickness of the gas diffusion layer in meters
    epsilon_gdl::Float64  = 0.7              # Anode/cathode GDL porosity
    # ── Microporous layer ──
    Hmpl::Float64         = 75e-6            # Thickness of the microporous layer in meters
    epsilon_mpl::Float64  = 0.4              # Porosity of the microporous layer
    # ── Gas channel ──
    Hagc::Float64         = 500e-6           # Thickness of the anode gas channel in meters
    Hcgc::Float64         = 500e-6           # Thickness of the cathode gas channel in meters
    Wagc::Float64         = 450e-6           # Width of the anode gas channel in meters
    Wcgc::Float64         = 450e-6           # Width of the cathode gas channel in meters
    Lgc::Float64          = 0.144            # Length of the gas channel in meters
    nb_channel_in_gc::Int = 67               # Number of channels in the bipolar plate
    Ldist::Float64        = 0.05             # Length of the distributor (between gas channel and manifold) in meters
    # ── Auxiliaries ──
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
    # ── Model parameters for fluidic calculation ──
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
    gamma_cond::Float64     = 1e8            # m-2. It is the effective gas-liquid phase-change geometry/closure factor for overall condensation rate of water [Ansys Fluent value from their User Guide https://de.scribd.com/document/344930548/ANSYS-Fluent-Advanced-Add-On-Modules (2016 version, page 238) / https://www.sciencedirect.com/science/article/pii/S0196890419303863 173 pp.].
    gamma_evap::Float64     = 1e8            # m-2. It is the effective gas-liquid phase-change geometry/closure factor for overall evaporation rate of water [Ansys Fluent value from their User Guide https://de.scribd.com/document/344930548/ANSYS-Fluent-Advanced-Add-On-Modules (2016 version, page 238) / https://www.sciencedirect.com/science/article/pii/S0196890419303863 173 pp.].
    epsilon_p::Float64      = 0.11           # . It is the percolation threshold porosity of the GDL.
    epsilon_c::Float64      = 0.2            # . It is the compression ratio of the GDL.
    alpha_p::Float64        = 0.785          # . It is a fitted value for the effective matter transfer in the GDL, for through plane direction.
    Eact_H2_cros_v::Float64 = 2.1e4          # J.mol-1. It is the activation energy of H2 for crossover in the under saturated membrane.
    Eact_H2_cros_l::Float64 = 1.8e4          # J.mol-1. It is the activation energy of H2 for crossover in the liquid-equilibrated membrane.
    Eact_O2_cros_v::Float64 = 2.2e4          # J.mol-1. It is the activation energy of oxygen for crossover in the under saturated membrane.
    Eact_O2_cros_l::Float64 = 2.0e4          # J.mol-1. It is the activation energy of oxygen for crossover in the liquid-equilibrated membrane.
    Kshape::Int64           = 2              # . Mathematical factor governing lambda_eq smoothing.
    # ── Interaction parameters between fluids and PEMFC structure ──
    e::Int64                = 3              # . Capillary exponent
    # ── Volumic flow of O2 inside the CCL to the Pt sites ──
    K_O2_ad_Pt::Float64    = 5.4             # Interfacial resistance coefficient of O2 adsorption on the Pt sites
    K_O2_dis_ion::Float64  = 8.5             # . It is the interfacial resistance coefficient of O2 dissolution inside the ionomer [haoModelingExperimentalValidation2015].
    K_O2_dis_l::Float64    = 1.0             # . It is the interfacial resistance coefficient of O2 dissolution inside the CL liquid water.
    IC_acl::Float64        = 0.6             # . It is the ionomer to carbon ratio in the anode catalyst layer [haoModelingExperimentalValidation2015].
    IC_ccl::Float64        = 0.95            # . It is the ionomer to carbon ratio in the cathode catalyst layer [haoModelingExperimentalValidation2015].
    ECSA_0::Float64        = 150             # cm2_Pt.cm-2_active_area. It is the initial electrochemical surface area of the catalyst.
    wt_Pt_acl::Float64     = 0.2             # It is the weight fraction of platinum over carbon covered by platinum (Pt/C) in the anode catalyst layer.
    wt_Pt_ccl::Float64     = 0.5             # It is the weight fraction of platinum over carbon covered by platinum (Pt/C) in the cathode catalyst layer [haoModelingExperimentalValidation2015].
    L_Pt_acl::Float64      = 0.5e-3          # kg.m-2. It is the platinum loading in the anode catalyst layer [haoModelingExperimentalValidation2015]..
    L_Pt_ccl::Float64      = 3e-3            # kg.m-2. It is the platinum loading in the cathode catalyst layer.
    r_carb::Float64        = 40e-9           # m. It is the mean radius of the carbon particles.
    theta_Pt_0::Float64    = 0               # This is the initial platine-oxide coverage, assumed to be zero for simplification.
    rho_ion::Float64       = 1900            # kg.m-3. It is the density of the ionomer [haoModelingExperimentalValidation2015].
    # ── Voltage calculation ──
    Re::Float64            = 1e-6            # Electron conduction resistance of the circuit in Ω·m²
    i0_c_ref::Float64      = 14.43           # Reference exchange current density at the cathode in A·m⁻²
    kappa_co::Float64      = 30.0            # Crossover correction coefficient in mol·m⁻¹·s⁻¹·Pa⁻¹
    kappa_c::Float64       = 1.0             # Overpotential correction exponent
    Eact_O2_red::Float64   = 27.7e3          # J.mol-1. Activation energy of oxygen reduction [futterPhysicalModelingPolymerelectrolyte2018].
    C_scl::Float64         = 2e7             # Volumetric space-charge layer capacitance in F·m⁻³
    alpha_c::Float64       = 0.5             # It is the transfer coefficient of the cathode.

    # ════════════════════════════════════════════════════════════════════════════
    # Model parameters for the heat transfer calculation
    # ════════════════════════════════════════════════════════════════════════════
    # ── Thermal conductivities ──
    k_th_gdl::Float64      = 0.3             # J.m-1.s-1.K-1. It is the thermal non-effective conductivity of the GDLs (non-effective ?) [ZSW].
    k_th_mpl::Float64      = 0.27            # J.m-1.s-1.K-1. It is the thermal conductivity of the MPLs [kotakaImpactInterfacialWater2014].
    k_th_cl::Float64       = 0.27            # J.m-1.s-1.K-1. It is the thermal conductivity of the CLs [vetterFreeOpenReference2019].
    k_th_mem::Float64      = 0.3             # J.m-1.s-1.K-1. It is the thermal conductivity of the membrane [vetterFreeOpenReference2019].
    # ── Specific heat capacities ──
    Cp_gdl::Float64        = 568             # J.kg-1.K-1. It is the specific heat capacities of the GDLs [wangQuasi2DTransientModel2018].
    Cp_mpl::Float64        = 568             # J.kg-1.K-1. It is the specific heat capacities of the MPLs [yangEffectsOperatingConditions2019].
    Cp_cl::Float64         = 3300            # J.kg-1.K-1. It is the specific heat capacities the CLs [wangQuasi2DTransientModel2018].
    Cp_mem::Float64        = 833             # J.kg-1.K-1. It is the specific heat capacities of the membrane [wangQuasi2DTransientModel2018].
    # ── Densities ──
    rho_gdl::Float64       = 1000            # kg.m-3. It is the density of the GDLs [wangQuasi2DTransientModel2018].
    rho_mpl::Float64       = 1000            # kg.m-3. It is the density of the MPLs [yangEffectsOperatingConditions2019].
    rho_cl::Float64        = 1000            # kg.m-3. It is the density of the CLs [wangQuasi2DTransientModel2018].
    # ── Electrical conductivities ──
    sigma_e_gdl::Float64   = 1250            # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
    sigma_e_mpl::Float64   = 5000            # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [yangEffectsOperatingConditions2019].
    sigma_e_cl::Float64    = 350             # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
end


# ============================================================
# PARAMETER BOUNDS
# ============================================================

"""
    PARAMETER_BOUNDS

Dictionary of the physical parameters whose value may fluctuate, with their fluctuation interval.
It is the reference list of bounds. The parameters that are actually left undetermined
(calibrated or sampled) are selected in the stack files of the `data` directory.

Each parameter is mapped to a tuple: (min::Float64, max::Float64, type::Symbol)
where type is either :real or :int.

Sources for some intervals:
- [Goshtasbi 2020] A. Goshtasbi et al., J. Electrochem. Soc. 167, 044504 (2020) (Part I and II, Table A·I).
- [Nobrega 2026] P. Affonso Nobrega et al., J. Electrochem. Soc. 173, 114503 (2026) (Table II), marapendi.
"""
const PARAMETER_BOUNDS = Dict{Symbol, Tuple{Float64, Float64, Symbol}}(
    # ════════════════════════════════════════════════════════════════════════════
    # Geometrical parameters of the fuel cell
    # ════════════════════════════════════════════════════════════════════════════
    # ── Catalyst layer ──
    :Hacl          => (1e-6, 15e-6, :real),                  # Anode catalyst-layer thickness
    :Hccl          => (5e-6, 20e-6, :real),                  # Cathode catalyst-layer thickness
    # ── Membrane ──
    :Hmem          => (5e-6, 50e-6, :real),                  # Membrane thickness
    # ── Gas diffusion layer ──
    :Hgdl          => (70e-6, 150e-6, :real),                # Gas-diffusion-layer thickness
    :epsilon_gdl   => (0.5, 0.9, :real),                     # GDL porosity
    # ── Microporous layer ──
    :Hmpl          => (40e-6, 100e-6, :real),                # Microporous-layer thickness
    :epsilon_mpl   => (0.3, 0.6, :real),                     # MPL porosity
    # ── Gas channel ──
    # :Hagc          => (1e-4, 1e-3, :real),                   # Anode gas-channel height. [Goshtasbi 2020] h_ch [0.01, 0.1] cm.
    # :Hcgc          => (1e-4, 1e-3, :real),                   # Cathode gas-channel height. [Goshtasbi 2020] h_ch [0.01, 0.1] cm.
    # :Wagc          => (1e-3, 4e-3, :real),                   # Anode gas-channel width. [Goshtasbi 2020] w_ch [0.1, 0.4] cm. NB: default value (450e-6) is below this range.
    # :Wcgc          => (1e-3, 4e-3, :real),                   # Cathode gas-channel width. [Goshtasbi 2020] w_ch [0.1, 0.4] cm. NB: default value (450e-6) is below this range.
    # :Lgc           => (0.15, 0.40, :real),                   # Gas-channel length. [Goshtasbi 2020] L_ch [15, 40] cm. NB: default value (0.144) is just below this range.

    # ════════════════════════════════════════════════════════════════════════════
    # Physical parameters of the model
    # ════════════════════════════════════════════════════════════════════════════
    # ── Model parameters for fluidic calculation ──
    # :rho_mem          => (?, ?, :real),   # ⚠️ A COMPLÉTER. Dry membrane density. No source.
    # :M_eq          => (0.6, 1.2, :real),                     # Equivalent molar mass of the ionomer. [Goshtasbi 2020] EW_mb [600, 1200] g/mol (wider than [Nobrega 2026] [700, 1100]).
    :tau_mpl       => (1.0, 4.0, :real),                     # Pore structure coefficient in the MPL
    # :tau_cl           => (?, ?, :real),   # ⚠️ A COMPLÉTER. Pore structure coefficient in the CL. Doubt: Goshtasbi n_Brugg (CL Bruggeman coefficient) [1, 3] may be related, but not the same definition.
    # :r_s_gdl          => (?, ?, :real),   # ⚠️ A COMPLÉTER. Pore-blockage exponent in the GDL. Doubt: Goshtasbi n_v [1, 3] / Nobrega n_s (1.5, 3) "saturation correction exponent" may correspond.
    # :r_s_mpl          => (?, ?, :real),   # ⚠️ A COMPLÉTER. Pore-blockage exponent in the MPL. Same doubt as r_s_gdl.
    # :r_s_cl           => (?, ?, :real),   # ⚠️ A COMPLÉTER. Pore-blockage exponent in the CL. Same doubt as r_s_gdl.
    # :Dp_gdl           => (?, ?, :real),   # ⚠️ A COMPLÉTER. GDL pore diameter. No source.
    # :Dp_mpl           => (?, ?, :real),   # ⚠️ A COMPLÉTER. MPL pore diameter. No source.
    # :Dp_cl         => (10e-9, 70e-9, :real),                 # CL pore diameter. [Nobrega 2026] d_p,CL (10, 70) nm. NB: default value (0.15e-6) is above this range.
    # :theta_c_gdl   => (90 * π / 180, 160 * π / 180, :real),  # GDL contact angle
    # :theta_c_mpl   => (110 * π / 180, 150 * π / 180, :real), # MPL contact angle. [Goshtasbi 2020] θ_MPL [110, 150]°.
    # :theta_c_cl    => (90 * π / 180, 140 * π / 180, :real),  # CL contact angle
    # :theta_l_rem      => (?, ?, :real),   # ⚠️ A COMPLÉTER. Liquid-water removal coefficient from the GDL to the GC. No source.
    # :K_v_liq_gas      => (?, ?, :real),   # ⚠️ A COMPLÉTER. Liquid-to-gas velocity ratio in the GC. No source.
    # :gamma_sorp_l     => (?, ?, :real),   # ⚠️ A COMPLÉTER. Sorption rate of liquid water in the membrane. No source (Nobrega k_abs is a different formulation).
    # :D_liq_dif        => (?, ?, :real),   # ⚠️ A COMPLÉTER. Diffusion coefficient of liquid water in the GC. No source.
    :gamma_cond    => (1e5, 1e9, :real),                     # Effective gas-liquid phase-change geometry/closure factor for overall condensation rate of water
    # :gamma_evap       => (?, ?, :real),   # ⚠️ A COMPLÉTER. Phase-change factor for evaporation (gamma_cond has bounds). No source.
    # :epsilon_p        => (?, ?, :real),   # ⚠️ A COMPLÉTER. Percolation threshold porosity of the GDL. No source.
    # :epsilon_c        => (?, ?, :real),   # ⚠️ A COMPLÉTER. Compression ratio of the GDL. No source.
    # :alpha_p          => (?, ?, :real),   # ⚠️ A COMPLÉTER. Fitted effective matter-transfer value in the GDL. No source.
    # :Eact_H2_cros_v   => (?, ?, :real),   # ⚠️ A COMPLÉTER. H2 crossover activation energy, vapor-equilibrated membrane. No source.
    # :Eact_H2_cros_l   => (?, ?, :real),   # ⚠️ A COMPLÉTER. H2 crossover activation energy, liquid-equilibrated membrane. No source.
    # :Eact_O2_cros_v   => (?, ?, :real),   # ⚠️ A COMPLÉTER. O2 crossover activation energy, vapor-equilibrated membrane. No source.
    # :Eact_O2_cros_l   => (?, ?, :real),   # ⚠️ A COMPLÉTER. O2 crossover activation energy, liquid-equilibrated membrane. No source.
    # ── Interaction parameters between fluids and PEMFC structure ──
    :e             => (3, 5, :int),                          # Capillary exponent
    # ── Volumic flow of O2 inside the CCL to the Pt sites ──
    # :K_O2_ad_Pt       => (?, ?, :real),   # ⚠️ A COMPLÉTER. O2 adsorption interfacial coefficient on Pt. Doubt: Nobrega k1, k2, k3 (2, 10) are interfacial resistance coefficients, but the mapping to K_O2_* is unclear.
    # :K_O2_dis_ion     => (?, ?, :real),   # ⚠️ A COMPLÉTER. O2 dissolution interfacial coefficient in the ionomer. Same doubt (Nobrega k1..k3 (2, 10)).
    # :K_O2_dis_l       => (?, ?, :real),   # ⚠️ A COMPLÉTER. O2 dissolution interfacial coefficient in liquid water. Same doubt (default 1.0 is outside (2, 10)).
    # :IC_acl        => (0.6, 1.3, :real),                     # Ionomer to carbon ratio in the anode catalyst layer. [Goshtasbi 2020] IC_Pt,an [0.6, 1.3].
    :IC_ccl        => (0.1, 2.0, :real),                     # Ionomer to carbon ratio in the cathode catalyst layer
    :ECSA_0        => (30.0, 200.0, :real),                  # Initial electrochemical surface area of the catalyst
    # :wt_Pt_acl     => (0.2, 0.7, :real),                     # Weight fraction of platinum over carbon in the anode catalyst layer. [Goshtasbi 2020] wt%Pt,an [0.2, 0.7].
    :wt_Pt_ccl     => (0.1, 0.7, :real),                     # Weight fraction of platinum over carbon in the cathode catalyst layer
    # :L_Pt_acl      => (5e-4, 3e-3, :real),                   # Platinum loading in the anode catalyst layer. [Goshtasbi 2020] L_Pt,an [0.05, 0.3] mg/cm² (1 mg/cm² = 1e-2 kg/m²).
    :L_Pt_ccl      => (1e-3, 5e-3, :real),                   # Platinum loading in the cathode catalyst layer
    :r_carb        => (10e-9, 100e-9, :real),                # Mean radius of the carbon particles
    # :rho_ion          => (?, ?, :real),   # ⚠️ A COMPLÉTER. Ionomer density. No source.
    # ── Voltage calculation ──
    :Re            => (5e-8, 5e-6, :real),                   # Electron-conduction resistance
    :i0_c_ref      => (0.1, 100.0, :real),                   # Reference cathode exchange current density
    :kappa_co      => (0.01, 40.0, :real),                   # Crossover correction coefficient
    :kappa_c       => (0.25, 4.0, :real),                    # Overpotential correction exponent
    # :Eact_O2_red   => (30e3, 80e3, :real),                   # Activation energy of oxygen reduction. [Goshtasbi 2020] and [Nobrega 2026] E_act,ca [30, 80] kJ/mol. NB: default value (27.7e3) is below this range.
    # :C_scl            => (?, ?, :real),   # ⚠️ A COMPLÉTER. Space-charge layer capacitance. No source.
    :alpha_c       => (0.5, 1.0, :real),                     # Cathode transfer coefficient

    # ════════════════════════════════════════════════════════════════════════════
    # Model parameters for the heat transfer calculation
    # ════════════════════════════════════════════════════════════════════════════
    # ── Thermal conductivities ──
    :k_th_gdl      => (0.1, 1.0, :real),                     # Thermal conductivity of the GDL
    # :k_th_mpl         => (?, ?, :real),   # ⚠️ A COMPLÉTER. MPL thermal conductivity. No source (k_th_gdl bounds exist; Goshtasbi k_T is a scaling factor).
    # :k_th_cl          => (?, ?, :real),   # ⚠️ A COMPLÉTER. CL thermal conductivity. No source.
    # :k_th_mem         => (?, ?, :real),   # ⚠️ A COMPLÉTER. Membrane thermal conductivity. No source.
    # ── Specific heat capacities ──
    # :Cp_gdl           => (?, ?, :real),   # ⚠️ A COMPLÉTER. GDL specific heat capacity. No source.
    # :Cp_mpl           => (?, ?, :real),   # ⚠️ A COMPLÉTER. MPL specific heat capacity. No source.
    # :Cp_cl            => (?, ?, :real),   # ⚠️ A COMPLÉTER. CL specific heat capacity. No source.
    # :Cp_mem           => (?, ?, :real),   # ⚠️ A COMPLÉTER. Membrane specific heat capacity. No source.
    # ── Densities ──
    # :rho_gdl          => (?, ?, :real),   # ⚠️ A COMPLÉTER. GDL density. No source.
    # :rho_mpl          => (?, ?, :real),   # ⚠️ A COMPLÉTER. MPL density. No source.
    # :rho_cl           => (?, ?, :real),   # ⚠️ A COMPLÉTER. CL density. No source.
    # ── Electrical conductivities ──
    # :sigma_e_gdl      => (?, ?, :real),   # ⚠️ A COMPLÉTER. GDL electrical conductivity. No source (Re already covers the total electronic resistance).
    # :sigma_e_mpl      => (?, ?, :real),   # ⚠️ A COMPLÉTER. MPL electrical conductivity. No source.
    # :sigma_e_cl       => (?, ?, :real),   # ⚠️ A COMPLÉTER. CL electrical conductivity. No source.
)


"""
    PARAMETER_METADATA

Metadata for the parameters of `PARAMETER_BOUNDS`: unit and description.
Used for displaying and exporting parameter bounds.
"""
const PARAMETER_METADATA = Dict{Symbol, Tuple{String, String}}(
    # ════════════════════════════════════════════════════════════════════════════
    # Geometrical parameters of the fuel cell
    # ════════════════════════════════════════════════════════════════════════════
    # ── Catalyst layer ──
    :Hacl          => ("m", "Anode catalyst-layer thickness"),
    :Hccl          => ("m", "Cathode catalyst-layer thickness"),
    # ── Membrane ──
    :Hmem          => ("m", "Membrane thickness"),
    # ── Gas diffusion layer ──
    :Hgdl          => ("m", "Gas-diffusion-layer thickness"),
    :epsilon_gdl   => ("—", "GDL porosity"),
    # ── Microporous layer ──
    :Hmpl          => ("m", "Microporous-layer thickness"),
    :epsilon_mpl   => ("—", "MPL porosity"),
    # ── Gas channel ──
    :Hagc          => ("m", "Anode gas-channel height"),
    :Hcgc          => ("m", "Cathode gas-channel height"),
    :Wagc          => ("m", "Anode gas-channel width"),
    :Wcgc          => ("m", "Cathode gas-channel width"),
    :Lgc           => ("m", "Gas-channel length"),

    # ════════════════════════════════════════════════════════════════════════════
    # Physical parameters of the model
    # ════════════════════════════════════════════════════════════════════════════
    # ── Model parameters for fluidic calculation ──
    :rho_mem       => ("kg·m⁻³", "Density of the dry membrane"),
    :M_eq          => ("kg·mol⁻¹", "Equivalent molar mass of the ionomer"),
    :tau_mpl       => ("—", "Pore structure coefficient in the MPL"),
    :tau_cl        => ("—", "Pore structure coefficient in the CL"),
    :r_s_gdl       => ("—", "Pore-blockage exponent in the GDL"),
    :r_s_mpl       => ("—", "Pore-blockage exponent in the MPL"),
    :r_s_cl        => ("—", "Pore-blockage exponent in the CL"),
    :Dp_gdl        => ("m", "GDL pore diameter"),
    :Dp_mpl        => ("m", "MPL pore diameter"),
    :Dp_cl         => ("m", "CL pore diameter"),
    :theta_c_gdl   => ("rad", "GDL contact angle"),
    :theta_c_mpl   => ("rad", "MPL contact angle"),
    :theta_c_cl    => ("rad", "CL contact angle"),
    :theta_l_rem   => ("s·m⁻¹", "Liquid-water removal coefficient from the GDL to the GC"),
    :K_v_liq_gas   => ("—", "Liquid-to-gas velocity ratio in the GC"),
    :gamma_sorp_l  => ("s⁻¹", "Sorption rate of liquid water in the membrane"),
    :D_liq_dif     => ("kg·m⁻¹·s⁻¹", "Diffusion coefficient of liquid water in the GC"),
    :gamma_cond    => ("m⁻²", "Effective gas-liquid phase-change geometry/closure factor for overall condensation rate of water"),
    :gamma_evap    => ("m⁻²", "Effective gas-liquid phase-change geometry/closure factor for overall evaporation rate of water"),
    :epsilon_p     => ("—", "Percolation threshold porosity of the GDL"),
    :epsilon_c     => ("—", "Compression ratio of the GDL"),
    :alpha_p       => ("—", "Fitted effective matter-transfer value in the GDL (through-plane)"),
    :Eact_H2_cros_v => ("J·mol⁻¹", "Activation energy of H2 crossover in the under-saturated membrane"),
    :Eact_H2_cros_l => ("J·mol⁻¹", "Activation energy of H2 crossover in the liquid-equilibrated membrane"),
    :Eact_O2_cros_v => ("J·mol⁻¹", "Activation energy of O2 crossover in the under-saturated membrane"),
    :Eact_O2_cros_l => ("J·mol⁻¹", "Activation energy of O2 crossover in the liquid-equilibrated membrane"),
    # ── Interaction parameters between fluids and PEMFC structure ──
    :e             => ("—", "Capillary exponent"),
    # ── Volumic flow of O2 inside the CCL to the Pt sites ──
    :K_O2_ad_Pt    => ("—", "Interfacial resistance coefficient of O2 adsorption on the Pt sites"),
    :K_O2_dis_ion  => ("—", "Interfacial resistance coefficient of O2 dissolution inside the ionomer"),
    :K_O2_dis_l    => ("—", "Interfacial resistance coefficient of O2 dissolution inside the CL liquid water"),
    :IC_acl        => ("—", "Ionomer to carbon ratio in the anode catalyst layer"),
    :IC_ccl        => ("—", "Ionomer to carbon ratio in the cathode catalyst layer"),
    :ECSA_0        => ("cm²_Pt.cm⁻²_active_area", "Initial electrochemical surface area of the catalyst"),
    :wt_Pt_acl     => ("—", "Weight fraction of platinum over carbon in the anode catalyst layer"),
    :wt_Pt_ccl     => ("—", "Weight fraction of platinum over carbon in the cathode catalyst layer"),
    :L_Pt_acl      => ("kg·m⁻²", "Platinum loading in the anode catalyst layer"),
    :L_Pt_ccl      => ("kg·m⁻²", "Platinum loading in the cathode catalyst layer"),
    :r_carb        => ("m", "Mean radius of the carbon particles"),
    :rho_ion       => ("kg·m⁻³", "Density of the ionomer"),
    # ── Voltage calculation ──
    :Re            => ("Ω·m²", "Electron-conduction resistance"),
    :i0_c_ref      => ("A·m⁻²", "Reference cathode exchange current density"),
    :kappa_co      => ("—", "Crossover correction coefficient"),
    :kappa_c       => ("—", "Overpotential correction exponent"),
    :Eact_O2_red   => ("J·mol⁻¹", "Activation energy of oxygen reduction"),
    :C_scl         => ("F·m⁻³", "Volumetric space-charge layer capacitance"),
    :alpha_c       => ("—", "Cathode transfer coefficient"),

    # ════════════════════════════════════════════════════════════════════════════
    # Model parameters for the heat transfer calculation
    # ════════════════════════════════════════════════════════════════════════════
    # ── Thermal conductivities ──
    :k_th_gdl      => ("W·m⁻¹·K⁻¹", "Thermal conductivity of the GDL"),
    :k_th_mpl      => ("W·m⁻¹·K⁻¹", "Thermal conductivity of the MPL"),
    :k_th_cl       => ("W·m⁻¹·K⁻¹", "Thermal conductivity of the CL"),
    :k_th_mem      => ("W·m⁻¹·K⁻¹", "Thermal conductivity of the membrane"),
    # ── Specific heat capacities ──
    :Cp_gdl        => ("J·kg⁻¹·K⁻¹", "Specific heat capacity of the GDL"),
    :Cp_mpl        => ("J·kg⁻¹·K⁻¹", "Specific heat capacity of the MPL"),
    :Cp_cl         => ("J·kg⁻¹·K⁻¹", "Specific heat capacity of the CL"),
    :Cp_mem        => ("J·kg⁻¹·K⁻¹", "Specific heat capacity of the membrane"),
    # ── Densities ──
    :rho_gdl       => ("kg·m⁻³", "Density of the GDL"),
    :rho_mpl       => ("kg·m⁻³", "Density of the MPL"),
    :rho_cl        => ("kg·m⁻³", "Density of the CL"),
    # ── Electrical conductivities ──
    :sigma_e_gdl   => ("Ω⁻¹·m⁻¹", "Electrical conductivity of the GDL"),
    :sigma_e_mpl   => ("Ω⁻¹·m⁻¹", "Electrical conductivity of the MPL"),
    :sigma_e_cl    => ("Ω⁻¹·m⁻¹", "Electrical conductivity of the CL"),
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
    y_H2_in::Float64 = 0.8          # Molar fraction of H2 in the dry anode gas mixture (H2/N2) injected at the inlet
    i_min_stoich::Float64 = 0.5     # Minimum current density used to compute the desired flows (A.cm-2)
end


"""
    OPERATING_CONDITIONS_BOUNDS

Dictionary of default bounds for operating conditions.
These bounds serve as fallback values when no fuel-cell-specific bounds are available.

Each parameter is mapped to a tuple: (min::Float64, max::Float64, type::Symbol)
where type is either :real or :int.
"""
const OPERATING_CONDITIONS_BOUNDS = Dict{Symbol, Tuple{Float64, Float64, Symbol}}(
    :T_des         => (40.0 + 273.15, 90.0 + 273.15, :real), # Desired fuel cell temperature in Kelvin
    :Pa_des        => (1.5e5, 3.5e5, :real),                 # Desired anode pressure in Pascal
    :Pc_des        => (1.5e5, 3.5e5, :real),                 # Desired cathode pressure in Pascal
    :Sa            => (1.5, 2.5, :real),                     # Stoichiometric ratio of hydrogen
    :Sc            => (1.5, 2.5, :real),                     # Stoichiometric ratio of oxygen
    :Phi_a_des     => (0.3, 0.9, :real),                     # Desired anode relative humidity
    :Phi_c_des     => (0.3, 0.9, :real),                     # Desired cathode relative humidity
)


"""
    OperatingConditionConstraint

A single constraint relating an operating condition `target` to one or more `sources`.

# Fields
- `target::Symbol`: The operating condition to enforce (e.g. `:Pc_des`).
- `sources::Vector{Symbol}`: The operating conditions the constraint depends on.
- `fn::Function`: Function `(sample_dict) -> value` that computes the bound or fixed value.
- `kind::Symbol`: Kind of constraint. One of:
  - `:(=)`: `target = fn(overrides)` (default).
  - `:(>=)`: `target >= fn(overrides)`; the sampled value is clamped from below.
  - `:(<=)`: `target <= fn(overrides)`; the sampled value is clamped from above.
- `active::Bool`: Whether the constraint is applied. Default `true`.

Example:
```julia
OperatingConditionConstraint(
    target = :Pc_des,
    sources = [:Pa_des],
    fn = d -> d[:Pa_des] - 0.5e5,
    kind = :(>=)
)
```
"""
struct OperatingConditionConstraint
    target::Symbol
    sources::Vector{Symbol}
    fn::Function
    kind::Symbol
    active::Bool
    name::String
end

function OperatingConditionConstraint(; target::Symbol,
                                        sources::Vector{Symbol},
                                        fn::Function,
                                        kind::Symbol = :(=),
                                        active::Bool = true,
                                        name::String = "")
    kind in (:(=), :(>=), :(<=)) ||
        throw(ArgumentError("OperatingConditionConstraint kind must be :(=), :(>=) or :(<=), got $kind"))
    return OperatingConditionConstraint(target, sources, fn, kind, active, name)
end


"""
    default_operating_condition_constraints()

Return the default set of operating-condition constraints used during Sobol sampling.

Default constraints:
- `Pc_des >= Pa_des - 0.5e5` (cathode pressure is at most 0.5e5 Pa below anode pressure)
- `Pc_des <= Pa_des + 0.5e5` (cathode pressure is at least 0.5e5 Pa above anode pressure)
"""
function default_operating_condition_constraints()::Vector{OperatingConditionConstraint}
    return [
        OperatingConditionConstraint(
            target = :Pc_des,
            sources = Symbol[:Pa_des],
            fn = d -> d[:Pa_des] - 0.5e5,
            kind = :(>=),
            name = "Pc_des_ge_Pa_des_minus_0.5e5"
        ),
        OperatingConditionConstraint(
            target = :Pc_des,
            sources = Symbol[:Pa_des],
            fn = d -> d[:Pa_des] + 0.5e5,
            kind = :(<=),
            name = "Pc_des_le_Pa_des_plus_0.5e5"
        )
    ]
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

