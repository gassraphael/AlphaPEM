# -*- coding: utf-8 -*-

"""This module contains physical constants which are used for modeling the PEM fuel cell."""
import math

# __________________________________________________Physical constants__________________________________________________

# Physical constants
F = 96485.3321233  # C.mol-1. It is the Faraday constant.
R = 8.31446261815324  # J.mol-1.K-1. It is the universal gas constant.
M_O2 = 3.2e-2  # kg.mol-1. It is the molar mass of O2.
M_H2 = 2e-3  # kg.mol-1. It is the molar mass of H2.
M_N2 = 2.8e-2  # kg.mol-1. It is the molar mass of N2.
M_H2O = M_H2 + M_O2 / 2  # kg.mol-1. It is the molar mass of H2O.
gamma = 1.401  # . It is the heat capacity ratio of dry air at 100°C.
gamma_H2 = 1.404  # . It is the heat capacity ratio of H2 at 100°C.

# External environmental parameters
Text = 298  # K. It is the outside temperature.
Pext = 101325  # Pa. It is the outside pressure.
Phi_ext = 0.4  # It is the outside relative humidity.
y_O2_ext = 0.2095  # . It is the molar fraction of O2 in dry air.

# Model parameters for the cell
rho_mem = 1980  # kg.m-3. It is the density of the dry membrane.
M_eq = 1.1  # kg.mol-1. It is the equivalent molar mass of ionomer.
tau_mpl = 2  # It is the pore structure coefficient in the MPL, without units [Gen Inoue 2016 Journal Power Sources].
tau_cl = 4  # It is the pore structure coefficient in the CL, without units [Gen Inoue 2016 Journal Power Sources].
r_s_gdl = 2.0 # It is the exponent pore blockage in the GDL.
r_s_mpl = 2.5 # It is the exponent pore blockage in the MPL.
r_s_cl = 2.5 # It is the exponent pore blockage in the CL.
Dp_gdl = 33.2e-6  # m. It is the pore diameter of the GDL [ZSW GenStack].
Dp_mpl = 17.4e-6 # m. It is the pore diameter of the MPL [morganUnderstandingGasDiffusion2014].
Dp_cl = 0.15e-6  # m. It is the pore diameter of the CL [Ali Malekian 2019 International Journal of Hydrogen Energy].
theta_c_gdl = 120 * math.pi / 180  # radian. It is the contact angle of GDL for liquid water.
theta_c_mpl = 135 * math.pi / 180  # radian. It is the contact angle of MPL for liquid water.
theta_c_cl = 95 * math.pi / 180  # radian. It is the contact angle of CL for liquid water.
theta_l_rem = 5e-5 # s/m. It is the coefficient of liquid water removal from the GDL to the GC [Ansys Fluent value from their User Guide].
K_v_liq_gas = 0.02 # . It is the liquid to gas velocity ratio in the GC [Ansys Fluent value from their User Guide].
D_liq_dif = 1e-5  # kg.m-1.s-1. It is the diffusion coefficient of liquid water in the GC [Ansys Fluent value from their User Guide].
gamma_cond = 1e8  # s-1. It is the overall condensation rate constant for water [Ansys Fluent value from their User Guide].
gamma_evap = 1e8  # s-1. It is the overall evaporation rate constant for water [Ansys Fluent value from their User Guide].
epsilon_p = 0.11 #. It is the percolation threshold porosity of the GDL.
alpha_p = 0.785 #. It is a fitted value for the effective matter transfer in the GDL, for through plane direction.
Tref_cross = 303.15  # K. It is the reference temperature for crossover.
Eact_H2_cros_v = 2.1e4  # J.mol-1. It is the activation energy of H2 for crossover in the under saturated membrane.
Eact_H2_cros_l = 1.8e4  # J.mol-1. It is the activation energy of H2 for crossover in the liquid-equilibrated membrane.
Eact_O2_cros_v = 2.2e4  # J.mol-1. It is the activation energy of oxygen for crossover in the under saturated membrane.
Eact_O2_cros_l = 2.0e4  # J.mol-1. It is the activation energy of oxygen for crossover in the liquid-equilibrated membrane.
Kshape = 2  # . Mathematical factor governing lambda_eq smoothing.
#   Volumic flow of O2 inside the CCL to the Pt sites
K_O2_dis_ion = 8.5 # . It is the interfacial resistance coefficient of O2 dissolution inside the ionomer [haoModelingExperimentalValidation2015].
K_O2_dis_l = 1.0 # . It is the interfacial resistance coefficient of O2 dissolution inside the CL liquid water.
IC = 0.5 # . It is the ionomer to carbon ratio in the catalyst layer.
rho_ion = 1900 # kg.m-3. It is the density of the ionomer [haoModelingExperimentalValidation2015].
rho_carb = 1950 # kg.m-3. It is the density of the carbon [haoModelingExperimentalValidation2015].
rho_Pt = 21450 # kg.m-3. It is the density of the platinum [haoModelingExperimentalValidation2015].
r_carb = 40e-9 # m. It is the radius of the carbon particles.
theta_Pt_0 = 0 # This is the initial platine-oxide coverage, assumed to be zero for simplification.
ECSA_0 = 150 # cm2_Pt.cm-2_active_area. It is the initial electrochemical surface area of the catalyst.
wt_Pt = 0.5 # It is the weight fraction of platinum over carbon covered by platinum (Pt/C) in the cathode catalyst layer [haoModelingExperimentalValidation2015].
L_Pt = 0.3e-2 # kg.m-2. It is the platinum loading in the cathode catalyst layer.
#   Voltage calculation
C_O2ref_red = 3.39  # mol.m-3. It is the reference concentration of oxygen for the reduction reaction.
alpha_c = 0.5  # It is the transfer coefficient of the cathode.
E0 = 1.229  # V. It is the standard-state reversible voltage.
Pref_eq = 1e5  # Pa. It is the reference pressure for the equilibrium potential calculation.
Eact_O2_red = 27.7e3  # J.mol-1. It is the activation energy of oxygen reduction [futterPhysicalModelingPolymerelectrolyte2018].
Tref_O2_red = 323 # K. It is the reference temperature for the activation energy of oxygen reduction [futterPhysicalModelingPolymerelectrolyte2018].

# Model parameters for the heat transfer calculation
#   Thermal conductivities
k_th_gdl = 0.3 # J.m-1.s-1.K-1. It is the thermal non-effective conductivity of the GDLs (non-effective ?) [ZSW].
k_th_mpl = 0.27 # J.m-1.s-1.K-1. It is the thermal conductivity of the MPLs [kotakaImpactInterfacialWater2014].
k_th_cl = 0.27 # J.m-1.s-1.K-1. It is the thermal conductivity of the CLs [vetterFreeOpenReference2019].
k_th_mem = 0.3 # J.m-1.s-1.K-1. It is the thermal conductivity of the membrane [vetterFreeOpenReference2019].
#   Specific heat capacities
Cp_gdl = 568 # J.kg-1.K-1. It is the specific heat capacities of the GDLs [wangQuasi2DTransientModel2018].
Cp_mpl = 568 # J.kg-1.K-1. It is the specific heat capacities of the MPLs [yangEffectsOperatingConditions2019].
Cp_cl = 3300 # J.kg-1.K-1. It is the specific heat capacities the CLs [wangQuasi2DTransientModel2018].
Cp_mem = 833 # J.kg-1.K-1. It is the specific heat capacities of the membrane [wangQuasi2DTransientModel2018].
#   Densities
rho_gdl = 1000 # kg.m-3. It is the density of the GDLs [wangQuasi2DTransientModel2018].
rho_mpl = 1000 # kg.m-3. It is the density of the MPLs [yangEffectsOperatingConditions2019].
rho_cl = 1000 # kg.m-3. It is the density of the CLs [wangQuasi2DTransientModel2018].
#   Electrical conductivities
sigma_e_gdl = 1250 # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
sigma_e_mpl = 5000 # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [yangEffectsOperatingConditions2019].
sigma_e_cl = 350 # Ω-1.m-1. It is the electrical conductivity of the GDL (non-effective ?) [vetterFreeOpenReference2019].
#   Molar entropy of reactions
delta_s_HOR = 0.104  # J.mol-1.K-1. It is the HOR molar reaction entropy [vetterFreeOpenReference2019].
delta_s_ORR = -163.3  # J.mol-1.K-1. It is the ORR molar reaction entropy [vetterFreeOpenReference2019].

# Model parameters for the balance of plant
i_min_inlet_flows = 0.3e4  # A.m-2. Minimum current density at which inlet mass flows are regulated to supply reactant flow.
tau_cp = 1  # s. It is the air compressor time constant.
tau_hum = 5  # s. It is the humidifier time constant.
Kp_T = 5e-8  # m².s-1.Pa-1. It is the proportional constant of the PD controller at the back pressure valve.
Kd_T = 1e-8  # m².Pa-1. It is the derivative constant of the PD controller at the back pressure valve.
