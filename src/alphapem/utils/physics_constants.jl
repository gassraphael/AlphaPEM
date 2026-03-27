# -*- coding: utf-8 -*-

"""This module contains physical constants which are used for modeling the PEM fuel cell."""

# __________________________________________________Physical constants__________________________________________________

# Physical constants
const F::Float64 = 96485.3321233  # C.mol-1. It is the Faraday constant.
const R::Float64 = 8.31446261815324  # J.mol-1.K-1. It is the universal gas constant.
const M_O2::Float64 = 3.2e-2  # kg.mol-1. It is the molar mass of O2.
const M_H2::Float64 = 2e-3  # kg.mol-1. It is the molar mass of H2.
const M_N2::Float64 = 2.8e-2  # kg.mol-1. It is the molar mass of N2.
const M_H2O::Float64 = M_H2 + M_O2 / 2  # kg.mol-1. It is the molar mass of H2O.
const gamma::Float64 = 1.401  # . It is the heat capacity ratio of dry air at 100°C.
const gamma_H2::Float64 = 1.404  # . It is the heat capacity ratio of H2 at 100°C.

# External environmental parameters
const Text::Float64 = 298  # K. It is the outside temperature.
const Pext::Float64 = 101325  # Pa. It is the outside pressure.
const Phi_ext::Float64 = 0.4  # It is the outside relative humidity.
const y_O2_ext::Float64 = 0.2095  # . It is the molar fraction of O2 in dry air.

# Model parameters for the cell
const rho_mem::Float64 = 1980  # kg.m-3. It is the density of the dry membrane.
const M_eq::Float64 = 1.1  # kg.mol-1. It is the equivalent molar mass of ionomer.
const tau_mpl::Float64 = 2  # It is the pore structure coefficient in the MPL, without units [Gen Inoue 2016 Journal Power Sources].
const tau_cl::Float64 = 4  # It is the pore structure coefficient in the CL, without units [Gen Inoue 2016 Journal Power Sources].
const r_s_gdl::Float64 = 2.0  # It is the exponent pore blockage in the GDL.
const r_s_mpl::Float64 = 2.5  # It is the exponent pore blockage in the MPL.
const r_s_cl::Float64 = 2.5  # It is the exponent pore blockage in the CL.
const Dp_gdl::Float64 = 33.2e-6  # m. It is the pore diameter of the GDL [ZSW GenStack].
const Dp_mpl::Float64 = 17.4e-6  # m. It is the pore diameter of the MPL [morganUnderstandingGasDiffusion2014].
const Dp_cl::Float64 = 0.15e-6  # m. It is the pore diameter of the CL [Ali Malekian 2019 International Journal of Hydrogen Energy].
const theta_c_gdl::Float64 = 120 * π / 180  # radian. It is the contact angle of GDL for liquid water.
const theta_c_mpl::Float64 = 135 * π / 180  # radian. It is the contact angle of MPL for liquid water.
const theta_c_cl::Float64 = 95 * π / 180   # radian. It is the contact angle of CL for liquid water.
const theta_l_rem::Float64 = 5e-5  # s/m. It is the coefficient of liquid water removal from the GDL to the GC [Ansys Fluent value from their User Guide].
const K_v_liq_gas::Float64 = 0.02  # . It is the liquid to gas velocity ratio in the GC [Ansys Fluent value from their User Guide].
const D_liq_dif::Float64 = 1e-5  # kg.m-1.s-1. It is the diffusion coefficient of liquid water in the GC [Ansys Fluent value from their User Guide].
const gamma_cond::Float64 = 1e8  # s-1. It is the overall condensation rate constant for water [Ansys Fluent value from their User Guide].
const gamma_evap::Float64 = 1e8  # s-1. It is the overall evaporation rate constant for water [Ansys Fluent value from their User Guide].
const epsilon_p::Float64 = 0.11  # . It is the percolation threshold porosity of the GDL.
const alpha_p::Float64 = 0.785  # . It is a fitted value for the effective matter transfer in the GDL, for through plane direction.
const Tref_cross::Float64 = 303.15  # K. It is the reference temperature for crossover.
const Eact_H2_cros_v::Float64 = 2.1e4  # J.mol-1. It is the activation energy of H2 for crossover in the under saturated membrane.
const Eact_H2_cros_l::Float64 = 1.8e4  # J.mol-1. It is the activation energy of H2 for crossover in the liquid-equilibrated membrane.
const Eact_O2_cros_v::Float64 = 2.2e4  # J.mol-1. It is the activation energy of oxygen for crossover in the under saturated membrane.
const Eact_O2_cros_l::Float64 = 2.0e4  # J.mol-1. It is the activation energy of oxygen for crossover in the liquid-equilibrated membrane.
const Kshape::Int64 = 2  # . Mathematical factor governing lambda_eq smoothing.
#   Volumic flow of O2 inside the CCL to the Pt sites
const K_O2_dis_ion::Float64 = 8.5  # . It is the interfacial resistance coefficient of O2 dissolution inside the ionomer [haoModelingExperimentalValidation2015].
const K_O2_dis_l::Float64 = 1.0  # . It is the interfacial resistance coefficient of O2 dissolution inside the CL liquid water.
const IC::Float64 = 0.5  # . It is the ionomer to carbon ratio in the catalyst layer.
const rho_ion::Float64 = 1900  # kg.m-3. It is the density of the ionomer [haoModelingExperimentalValidation2015].
const rho_carb::Float64 = 1950  # kg.m-3. It is the density of the carbon [haoModelingExperimentalValidation2015].
const rho_Pt::Float64 = 21450  # kg.m-3. It is the density of the platinum [haoModelingExperimentalValidation2015].
const r_carb::Float64 = 40e-9  # m. It is the radius of the carbon particles.
const theta_Pt_0::Float64 = 0  # This is the initial platine-oxide coverage, assumed to be zero for simplification.
const ECSA_0::Float64 = 150  # cm2_Pt.cm-2_active_area. It is the initial electrochemical surface area of the catalyst.
const wt_Pt::Float64 = 0.5  # It is the weight fraction of platinum over carbon covered by platinum (Pt/C) in the cathode catalyst layer [haoModelingExperimentalValidation2015].
const L_Pt::Float64 = 0.3e-2  # kg.m-2. It is the platinum loading in the cathode catalyst layer.
#   Voltage calculation
const C_O2ref_red::Float64 = 3.39  # mol.m-3. It is the reference concentration of oxygen for the reduction reaction.
const alpha_c::Float64 = 0.5  # It is the transfer coefficient of the cathode.
const E0::Float64 = 1.229  # V. It is the standard-state reversible voltage.
const Pref_eq::Float64 = 1e5  # Pa. It is the reference pressure for the equilibrium potential calculation.
const Eact_O2_red::Float64 = 27.7e3  # J.mol-1. It is the activation energy of oxygen reduction [futterPhysicalModelingPolymerelectrolyte2018].
const Tref_O2_red::Float64 = 323  # K. It is the reference temperature for the activation energy of oxygen reduction [futterPhysicalModelingPolymerelectrolyte2018].

# Model parameters for the heat transfer calculation
#   Thermal conductivities
const k_th_gdl::Float64 = 0.3   # J.m-1.s-1.K-1. It is the thermal non-effective conductivity of the GDLs (non-effective ?) [ZSW].
const k_th_mpl::Float64 = 0.27  # J.m-1.s-1.K-1. It is the thermal conductivity of the MPLs [kotakaImpactInterfacialWater2014].
const k_th_cl::Float64 = 0.27   # J.m-1.s-1.K-1. It is the thermal conductivity of the CLs [vetterFreeOpenReference2019].
const k_th_mem::Float64 = 0.3   # J.m-1.s-1.K-1. It is the thermal conductivity of the membrane [vetterFreeOpenReference2019].
#   Specific heat capacities
const Cp_gdl::Float64 = 568   # J.kg-1.K-1. It is the specific heat capacities of the GDLs [wangQuasi2DTransientModel2018].
const Cp_mpl::Float64 = 568   # J.kg-1.K-1. It is the specific heat capacities of the MPLs [yangEffectsOperatingConditions2019].
const Cp_cl::Float64 = 3300   # J.kg-1.K-1. It is the specific heat capacities the CLs [wangQuasi2DTransientModel2018].
const Cp_mem::Float64 = 833   # J.kg-1.K-1. It is the specific heat capacities of the membrane [wangQuasi2DTransientModel2018].
#   Densities
const rho_gdl::Float64 = 1000  # kg.m-3. It is the density of the GDLs [wangQuasi2DTransientModel2018].
const rho_mpl::Float64 = 1000  # kg.m-3. It is the density of the MPLs [yangEffectsOperatingConditions2019].
const rho_cl::Float64 = 1000   # kg.m-3. It is the density of the CLs [wangQuasi2DTransientModel2018].
#   Electrical conductivities
const sigma_e_gdl::Float64 = 1250  # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
const sigma_e_mpl::Float64 = 5000  # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [yangEffectsOperatingConditions2019].
const sigma_e_cl::Float64 = 350    # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
#   Molar entropy of reactions
const delta_s_HOR::Float64 = 0.104   # J.mol-1.K-1. It is the HOR molar reaction entropy [vetterFreeOpenReference2019].
const delta_s_ORR::Float64 = -163.3  # J.mol-1.K-1. It is the ORR molar reaction entropy [vetterFreeOpenReference2019].

# Model parameters for the balance of plant
const i_min_inlet_flows::Float64 = 0.3e4  # A.m-2. Minimum current density at which inlet mass flows are regulated to supply reactant flow.
const tau_cp::Float64 = 1   # s. It is the air compressor time constant.
const tau_hum::Float64 = 5  # s. It is the humidifier time constant.
const Kp_T::Float64 = 5e-8  # m².s-1.Pa-1. It is the proportional constant of the PD controller at the back pressure valve.
const Kd_T::Float64 = 1e-8  # m².Pa-1. It is the derivative constant of the PD controller at the back pressure valve.

