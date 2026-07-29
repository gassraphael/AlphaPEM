# -*- coding: utf-8 -*-

"""This module contains physical constants which are used for modeling the PEM fuel cell."""

# __________________________________________________Physical constants__________________________________________________

# Physical constants
const F::Float64           = 96485.3321233     # C.mol-1. It is the Faraday constant.
const R::Float64           = 8.31446261815324  # J.mol-1.K-1. It is the universal gas constant.
const M_O2::Float64        = 3.2e-2            # kg.mol-1. It is the molar mass of O2.
const M_H2::Float64        = 2e-3              # kg.mol-1. It is the molar mass of H2.
const M_N2::Float64        = 2.8e-2            # kg.mol-1. It is the molar mass of N2.
const M_H2O::Float64       = M_H2 + M_O2 / 2   # kg.mol-1. It is the molar mass of H2O.
const gamma::Float64       = 1.401             # . It is the heat capacity ratio of dry air at 100°C.
const gamma_H2::Float64    = 1.404             # . It is the heat capacity ratio of H2 at 100°C.

# External environmental parameters
const Text::Float64        = 298               # K. It is the outside temperature.
const Pext::Float64        = 101325            # Pa. It is the outside pressure.
const Phi_ext::Float64     = 0.4               # It is the outside relative humidity.
const y_O2_ext::Float64    = 0.2095            # . It is the molar fraction of O2 in dry air.

# Model parameters for the fluidic calculation
const Tref_cross::Float64  = 303.15            # K. It is the reference temperature for crossover.
#   Volumic flow of O2 inside the CCL to the Pt sites
const rho_carb::Float64    = 1950              # kg.m-3. It is the density of the carbon [haoModelingExperimentalValidation2015].
const rho_Pt::Float64      = 21450             # kg.m-3. It is the density of the platinum [haoModelingExperimentalValidation2015].
#   Voltage calculation
const C_O2ref_red::Float64 = 3.39              # mol.m-3. It is the reference concentration of oxygen for the reduction reaction.
const E0::Float64          = 1.229             # V. It is the standard-state reversible voltage.
const Pref_eq::Float64     = 1e5               # Pa. It is the reference pressure for the equilibrium potential calculation.
const Eact_O2_red::Float64 = 27.7e3            # J.mol-1. It is the activation energy of oxygen reduction [futterPhysicalModelingPolymerelectrolyte2018].
const Tref_O2_red::Float64 = 323               # K. It is the reference temperature for the activation energy of oxygen reduction [futterPhysicalModelingPolymerelectrolyte2018].

# Model parameters for the heat transfer calculation
#   Molar entropy of reactions
const delta_s_HOR::Float64 = 0.104             # J.mol-1.K-1. It is the HOR molar reaction entropy [vetterFreeOpenReference2019].
const delta_s_ORR::Float64 = -163.3            # J.mol-1.K-1. It is the ORR molar reaction entropy [vetterFreeOpenReference2019].

# Model parameters for the balance of plant
const i_min_inlet_flows::Float64 = 0.5e4       # A.m-2. Minimum current density at which inlet mass flows are regulated to supply reactant flow.
const tau_cp::Float64      = 1                 # s. It is the air compressor time constant.
const tau_hum::Float64     = 5                 # s. It is the humidifier time constant.
const Kp_T::Float64        = 5e-8              # m².s-1.Pa-1. It is the proportional constant of the PD controller at the back pressure valve.
const Kd_T::Float64        = 1e-8              # m².Pa-1. It is the derivative constant of the PD controller at the back pressure valve.

