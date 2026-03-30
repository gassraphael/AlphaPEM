# -*- coding: utf-8 -*-

"""This file contains the current density parameter structures."""

# ______________________________________________Current density parameters______________________________________________

abstract type AbstractFuelCellParams end

# ============================================================
# PHYSICAL PARAMETERS
# ============================================================

"""
    PhysicalParams

Physical and geometrical parameters of the fuel cell.
"""
Base.@kwdef struct PhysicalParams <: AbstractFuelCellParams

    # --- Catalyst layers ---
    Hacl::Float64 = 10e-6
    Hccl::Float64 = 10e-6

    # --- Membrane ---
    Hmem::Float64 = 20e-6

    # --- Gas diffusion layer ---
    Hgdl::Float64 = 200e-6
    epsilon_gdl::Float64 = 0.5
    epsilon_c::Float64 = 0.2

    # --- Microporous layer ---
    Hmpl::Float64 = 30e-6
    epsilon_mpl::Float64 = 0.4

    # --- Gas channels ---
    Hagc::Float64 = 500e-6
    Hcgc::Float64 = 500e-6
    Wagc::Float64 = 450e-6
    Wcgc::Float64 = 450e-6
    Lgc::Float64 = 0.144
    nb_channel_in_gc::Int = 67

    # --- Manifolds ---
    Ldist::Float64 = 0.05
    Lm::Float64 = 0.0258
    A_T_a::Float64 = 11.8e-4
    A_T_c::Float64 = 11.8e-4
    Vasm::Float64 = 7e-3
    Vcsm::Float64 = 7e-3
    Vaem::Float64 = 2.4e-3
    Vcem::Float64 = 2.4e-3

    # --- Global ---
    Aact::Float64 = 100e-4
    nb_cell::Int = 1

    # --- Electrochemistry ---
    e::Int = 5
    K_l_ads::Float64 = 1.0
    K_O2_ad_Pt::Float64 = 5.4
    Re::Float64 = 1e-6
    i0_c_ref::Float64 = 14.43
    kappa_co::Float64 = 30.0
    kappa_c::Float64 = 1.0
    C_scl::Float64 = 2e7
end


# ============================================================
# OPERATING PARAMETERS
# ============================================================

"""
    OperatingParams
"""
Base.@kwdef struct OperatingParams <: AbstractFuelCellParams
    T_des::Float64 = 353.15
    Pa_des::Float64 = 2e5
    Pc_des::Float64 = 2e5
    Sa::Float64 = 1.2
    Sc::Float64 = 2.0
    Phi_a_des::Float64 = 0.4
    Phi_c_des::Float64 = 0.6
    y_H2_in::Float64 = 1.0
end


# ============================================================
# NUMERICAL PARAMETERS
# ============================================================

"""
    NumericalParams
"""
Base.@kwdef struct NumericalParams <: AbstractFuelCellParams
    nb_gc::Int = 1
    nb_gdl::Int = 3
    nb_mpl::Int = 2

    purge_time::Float64 = 0.6
    delta_purge::Float64 = 15.0

    rtol::Float64 = 1e-6
    atol::Float64 = 1e-9
end