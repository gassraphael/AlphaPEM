# -*- coding: utf-8 -*-

"""This file represents all the matter flows inside the fuel cell system. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
include(joinpath(@__DIR__, "../../utils/physics_constants.jl"))
include(joinpath(@__DIR__, "../../utils/maths_functions.jl"))
include(joinpath(@__DIR__, "../modules/cell_voltage_modules.jl"))
include(joinpath(@__DIR__, "../modules/flows_1D_MEA_modules.jl"))


# ________________________________________________________Flows_________________________________________________________

"""Calculate the flows inside the fuel cell system.

Parameters
----------
sv_1D : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
i_fc :
    Fuel cell current density at time t (A.m-2).
v_a :
    Anode gas velocity at time t (m.s-1).
v_c :
    Cathode gas velocity at time t (m.s-1).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.

Returns
-------
Dict{String, Dict}
    Flows inside the fuel cell system. Julia vectors are 1-based and do not contain the
    dummy element formerly inserted in Python at index 0. Therefore, inter-node flow vectors
    have length `nb_gdl - 1` or `nb_mpl - 1`, while node-based source vectors have length
    `nb_gdl` or `nb_mpl`.
"""
function calculate_flows_1D_MEA(sv_1D::Dict,
                                i_fc,
                                v_a,
                                v_c,
                                fc::AbstractFuelCell)::Dict{String, Dict}

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of parameters
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des = fc.operating_conditions.T_des
    Hmem, Hacl, Hccl = pp.Hmem, pp.Hacl, pp.Hccl
    Wagc, Wcgc, Hagc, Hcgc = pp.Wagc, pp.Wcgc, pp.Hagc, pp.Hcgc
    epsilon_gdl, epsilon_mpl = pp.epsilon_gdl, pp.epsilon_mpl
    K_l_ads, kappa_co = pp.K_l_ads, pp.kappa_co
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv_1D["C_v_agc"], sv_1D["C_v_acl"], sv_1D["C_v_ccl"], sv_1D["C_v_cgc"]
    C_v_agdl, C_v_ampl = [sv_1D["C_v_agdl_$i"] for i in 1:nb_gdl], [sv_1D["C_v_ampl_$i"] for i in 1:nb_mpl]
    C_v_cmpl, C_v_cgdl = [sv_1D["C_v_cmpl_$i"] for i in 1:nb_mpl], [sv_1D["C_v_cgdl_$i"] for i in 1:nb_gdl]
    s_acl, s_ccl = sv_1D["s_acl"], sv_1D["s_ccl"]
    s_agdl, s_ampl = [sv_1D["s_agdl_$i"] for i in 1:nb_gdl], [sv_1D["s_ampl_$i"] for i in 1:nb_mpl]
    s_cmpl, s_cgdl = [sv_1D["s_cmpl_$i"] for i in 1:nb_mpl], [sv_1D["s_cgdl_$i"] for i in 1:nb_gdl]
    lambda_acl, lambda_mem, lambda_ccl = sv_1D["lambda_acl"], sv_1D["lambda_mem"], sv_1D["lambda_ccl"]
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv_1D["C_H2_agc"], sv_1D["C_H2_acl"], sv_1D["C_O2_ccl"], sv_1D["C_O2_cgc"]
    C_H2_agdl, C_H2_ampl = [sv_1D["C_H2_agdl_$i"] for i in 1:nb_gdl], [sv_1D["C_H2_ampl_$i"] for i in 1:nb_mpl]
    C_O2_cmpl, C_O2_cgdl = [sv_1D["C_O2_cmpl_$i"] for i in 1:nb_mpl], [sv_1D["C_O2_cgdl_$i"] for i in 1:nb_gdl]
    C_N2_agc, C_N2_cgc = sv_1D["C_N2_agc"], sv_1D["C_N2_cgc"]
    T_acl, T_mem, T_ccl = sv_1D["T_acl"], sv_1D["T_mem"], sv_1D["T_ccl"]
    T_agdl, T_ampl = [sv_1D["T_agdl_$i"] for i in 1:nb_gdl], [sv_1D["T_ampl_$i"] for i in 1:nb_mpl]
    T_cmpl, T_cgdl = [sv_1D["T_cmpl_$i"] for i in 1:nb_mpl], [sv_1D["T_cgdl_$i"] for i in 1:nb_gdl]

    # Intermediate values
    (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc, D_eff_EOD_acl_mem,
     D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem, D_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl,
     D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl,
     Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl,
     Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl) = flows_1D_MEA_int_values(sv_1D, i_fc, fc)

    # ________________________________________Dissolved water flows (mol.m-2.s-1)_______________________________________

    # Anode side
    J_lambda_acl_mem = D_eff_EOD_acl_mem * interpolate([lambda_acl, lambda_mem], [Hacl, Hmem]) -
                       rho_mem / M_eq * D_lambda_eff_acl_mem * d_dx(lambda_acl, lambda_mem, Hacl / 2, Hmem / 2)
    # Cathode side
    J_lambda_mem_ccl = D_eff_EOD_mem_ccl * interpolate([lambda_mem, lambda_ccl], [Hmem, Hccl]) -
                       rho_mem / M_eq * D_lambda_eff_mem_ccl * d_dx(lambda_mem, lambda_ccl, Hmem / 2, Hccl / 2)

    # _________________________________________Liquid water flows (kg.m-2.s-1)__________________________________________

    # Anode side
    Jl_agc_agdl = -theta_l_rem * epsilon_gdl * s_agdl[1] * max(Pcap_agdl + rho_agc * v_a^2 / 2, 0.0)
    Jl_agdl_agdl = [-D_cap_agdl_agdl[i] * d_dx(s_agdl[i], s_agdl[i + 1], H_gdl_node / 2)
                    for i in 1:(nb_gdl - 1)]
    Jl_agdl_ampl = -D_cap_agdl_ampl * d_dx(s_agdl[nb_gdl], s_ampl[1], H_gdl_node / 2, H_mpl_node / 2)
    Jl_ampl_ampl = [-D_cap_ampl_ampl[i] * d_dx(s_ampl[i], s_ampl[i + 1], H_mpl_node / 2)
                    for i in 1:(nb_mpl - 1)]
    Jl_ampl_acl = -D_cap_ampl_acl * d_dx(s_ampl[nb_mpl], s_acl, H_mpl_node / 2, Hacl / 2)

    # Cathode side
    Jl_ccl_cmpl = -D_cap_ccl_cmpl * d_dx(s_ccl, s_cmpl[1], Hccl / 2, H_mpl_node / 2)
    Jl_cmpl_cmpl = [-D_cap_cmpl_cmpl[i] * d_dx(s_cmpl[i], s_cmpl[i + 1], H_mpl_node / 2)
                    for i in 1:(nb_mpl - 1)]
    Jl_cmpl_cgdl = -D_cap_cmpl_cgdl * d_dx(s_cmpl[nb_mpl], s_cgdl[1], H_mpl_node / 2, H_gdl_node / 2)
    Jl_cgdl_cgdl = [-D_cap_cgdl_cgdl[i] * d_dx(s_cgdl[i], s_cgdl[i + 1], H_gdl_node / 2)
                    for i in 1:(nb_gdl - 1)]
    Jl_cgdl_cgc = theta_l_rem * epsilon_gdl * s_cgdl[nb_gdl] * max(Pcap_cgdl + rho_cgc * v_c^2 / 2, 0.0)

    # _____________________________________________Vapor flows (mol.m-2.s-1)____________________________________________

    # Conductive-convective vapor flows
    Jv_agc_agdl = h_a(Pagc, T_des, Wagc, Hagc) * (C_v_agc - C_v_agdl[1])  # Also calculated in velocity.jl
    Jv_cgdl_cgc = h_c(Pcgc, T_des, Wcgc, Hcgc) * (C_v_cgdl[nb_gdl] - C_v_cgc)  # Also calculated in velocity.jl

    # Conductive vapor flows
    #   Anode side
    Jv_agdl_agdl = [-Da_eff_agdl_agdl[i] * d_dx(C_v_agdl[i], C_v_agdl[i + 1], H_gdl_node / 2)
                    for i in 1:(nb_gdl - 1)]
    Jv_agdl_ampl = -Da_eff_agdl_ampl * d_dx(C_v_agdl[nb_gdl], C_v_ampl[1], H_gdl_node / 2, H_mpl_node / 2)
    Jv_ampl_ampl = [-Da_eff_ampl_ampl[i] * d_dx(C_v_ampl[i], C_v_ampl[i + 1], H_mpl_node / 2)
                    for i in 1:(nb_mpl - 1)]
    Jv_ampl_acl = -Da_eff_ampl_acl * d_dx(C_v_ampl[nb_mpl], C_v_acl, H_mpl_node / 2, Hacl / 2)

    #   Cathode side
    Jv_ccl_cmpl = -Dc_eff_ccl_cmpl * d_dx(C_v_ccl, C_v_cmpl[1], Hccl / 2, H_mpl_node / 2)
    Jv_cmpl_cmpl = [-Dc_eff_cmpl_cmpl[i] * d_dx(C_v_cmpl[i], C_v_cmpl[i + 1], H_mpl_node / 2)
                    for i in 1:(nb_mpl - 1)]
    Jv_cmpl_cgdl = -Dc_eff_cmpl_cgdl * d_dx(C_v_cmpl[nb_mpl], C_v_cgdl[1], H_mpl_node / 2, H_gdl_node / 2)
    Jv_cgdl_cgdl = [-Dc_eff_cgdl_cgdl[i] * d_dx(C_v_cgdl[i], C_v_cgdl[i + 1], H_gdl_node / 2)
                    for i in 1:(nb_gdl - 1)]

    # ______________________________H2 and O2 flows (mol.m-2.s-1 for J, mol.m-3.s-1 for S)______________________________

    # Hydrogen and oxygen consumption
    #   Anode side
    S_H2_reac = i_fc / (2 * F * Hacl)
    S_H2_cros = R * T_acl_mem_ccl / (Hmem * Hacl) * (k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl +
                                                     2 * k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl)
    #   Cathode side
    S_O2_reac = i_fc / (4 * F * Hccl)
    S_O2_cros = R * T_acl_mem_ccl / (Hmem * Hccl) * (k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl +
                                                     1 / 2 * k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl)

    # Conductive-convective H2 and O2 flows
    J_H2_agc_agdl = h_a(Pagc, T_des, Wagc, Hagc) * (C_H2_agc - C_H2_agdl[1])  # Also calculated in velocity.jl
    J_O2_cgdl_cgc = h_c(Pcgc, T_des, Wcgc, Hcgc) * (C_O2_cgdl[nb_gdl] - C_O2_cgc)  # Also calculated in velocity.jl

    # Conductive H2 and O2 flows
    #   Anode side
    J_H2_agdl_agdl = [-Da_eff_agdl_agdl[i] * d_dx(C_H2_agdl[i], C_H2_agdl[i + 1], H_gdl_node / 2)
                      for i in 1:(nb_gdl - 1)]
    J_H2_agdl_ampl = -Da_eff_agdl_ampl * d_dx(C_H2_agdl[nb_gdl], C_H2_ampl[1], H_gdl_node / 2, H_mpl_node / 2)
    J_H2_ampl_ampl = [-Da_eff_ampl_ampl[i] * d_dx(C_H2_ampl[i], C_H2_ampl[i + 1], H_mpl_node / 2)
                      for i in 1:(nb_mpl - 1)]
    J_H2_ampl_acl = -Da_eff_ampl_acl * d_dx(C_H2_ampl[nb_mpl], C_H2_acl, H_mpl_node / 2, Hacl / 2)

    #   Cathode side
    J_O2_ccl_cmpl = -Dc_eff_ccl_cmpl * d_dx(C_O2_ccl, C_O2_cmpl[1], Hccl / 2, H_mpl_node / 2)
    J_O2_cmpl_cmpl = [-Dc_eff_cmpl_cmpl[i] * d_dx(C_O2_cmpl[i], C_O2_cmpl[i + 1], H_mpl_node / 2)
                      for i in 1:(nb_mpl - 1)]
    J_O2_cmpl_cgdl = -Dc_eff_cmpl_cgdl * d_dx(C_O2_cmpl[nb_mpl], C_O2_cgdl[1], H_mpl_node / 2, H_gdl_node / 2)
    J_O2_cgdl_cgdl = [-Dc_eff_cgdl_cgdl[i] * d_dx(C_O2_cgdl[i], C_O2_cgdl[i + 1], H_gdl_node / 2)
                      for i in 1:(nb_gdl - 1)]

    # __________________________________________Water generated (mol.m-3.s-1)___________________________________________

    # Water produced in the membrane at the CL through the chemical reaction and crossover
    #   Anode side
    Sp_acl = 2 * k_O2(lambda_mem, T_mem, kappa_co) * R * T_acl_mem_ccl / (Hmem * Hacl) * C_O2_ccl
    #   Cathode side
    Sp_ccl = i_fc / (2 * F * Hccl) + k_H2(lambda_mem, T_mem, kappa_co) * R * T_acl_mem_ccl / (Hmem * Hccl) * C_H2_acl

    # Water absorption in the CL due to the contact between the ionomer and vapor or liquid water:
    #   Anode side
    Sv_abs_acl = (1 - s_acl) * gamma_sorp(C_v_acl, s_acl, lambda_acl, T_acl, Hacl) * rho_mem / M_eq *
                 (lambda_eq(C_v_acl, s_acl, T_acl) - lambda_acl)
    if s_acl > 0
        Sl_abs_acl = s_acl * K_l_ads * gamma_sorp(C_v_acl, s_acl, lambda_acl, T_acl, Hacl) * rho_mem / M_eq *
                     (lambda_eq(C_v_acl, s_acl, T_acl) - lambda_acl)
    else
        Sl_abs_acl = 0.0
    end

    #   Cathode side
    Sv_abs_ccl = (1 - s_ccl) * gamma_sorp(C_v_ccl, s_ccl, lambda_ccl, T_ccl, Hccl) * rho_mem / M_eq *
                 (lambda_eq(C_v_ccl, s_ccl, T_ccl) - lambda_ccl)
    if s_ccl > 0
        Sl_abs_ccl = s_ccl * K_l_ads * gamma_sorp(C_v_ccl, s_ccl, lambda_ccl, T_ccl, Hccl) * rho_mem / M_eq *
                     (lambda_eq(C_v_ccl, s_ccl, T_ccl) - lambda_ccl)
    else
        Sl_abs_ccl = 0.0
    end

    # Liquid water generated through vapor condensation or degenerated through evaporation
    #   Anode side
    Sl_agdl = [Svl("anode", s_agdl[i], C_v_agdl[i],
                   C_v_agdl[i] + C_H2_agdl[i] + C_N2_agc,
                   T_agdl[i], epsilon_gdl) for i in 1:nb_gdl]
    Sl_ampl = [Svl("anode", s_ampl[i], C_v_ampl[i],
                   C_v_ampl[i] + C_H2_ampl[i] + C_N2_agc,
                   T_ampl[i], epsilon_mpl) for i in 1:nb_mpl]
    Sl_acl = Svl("anode", s_acl, C_v_acl, C_v_acl + C_H2_acl + C_N2_agc, T_acl,
                 epsilon_cl(lambda_acl, T_acl, Hacl))

    #   Cathode side
    Sl_ccl = Svl("cathode", s_ccl, C_v_ccl, C_v_ccl + C_O2_ccl + C_N2_cgc, T_ccl,
                 epsilon_cl(lambda_ccl, T_ccl, Hccl))
    Sl_cmpl = [Svl("cathode", s_cmpl[i], C_v_cmpl[i],
                   C_v_cmpl[i] + C_O2_cmpl[i] + C_N2_cgc,
                   T_cmpl[i], epsilon_mpl) for i in 1:nb_mpl]
    Sl_cgdl = [Svl("cathode", s_cgdl[i], C_v_cgdl[i],
                   C_v_cgdl[i] + C_O2_cgdl[i] + C_N2_cgc,
                   T_cgdl[i], epsilon_gdl) for i in 1:nb_gdl]

    # Vapor generated through liquid water evaporation or degenerated through condensation
    #   Anode side
    Sv_agdl = [-x for x in Sl_agdl]
    Sv_ampl = [-x for x in Sl_ampl]
    Sv_acl = -Sl_acl

    #   Cathode side
    Sv_ccl = -Sl_ccl
    Sv_cmpl = [-x for x in Sl_cmpl]
    Sv_cgdl = [-x for x in Sl_cgdl]

    # _____________________________________Assemble and return the flow dictionary______________________________________

    return Dict{String, Dict}(
        "Jv" => Dict(
            "agc_agdl" => Jv_agc_agdl,
            "agdl_agdl" => Jv_agdl_agdl,
            "agdl_ampl" => Jv_agdl_ampl,
            "ampl_ampl" => Jv_ampl_ampl,
            "ampl_acl" => Jv_ampl_acl,
            "ccl_cmpl" => Jv_ccl_cmpl,
            "cmpl_cmpl" => Jv_cmpl_cmpl,
            "cmpl_cgdl" => Jv_cmpl_cgdl,
            "cgdl_cgdl" => Jv_cgdl_cgdl,
            "cgdl_cgc" => Jv_cgdl_cgc
        ),
        "Jl" => Dict(
            "agc_agdl" => Jl_agc_agdl,
            "agdl_agdl" => Jl_agdl_agdl,
            "agdl_ampl" => Jl_agdl_ampl,
            "ampl_ampl" => Jl_ampl_ampl,
            "ampl_acl" => Jl_ampl_acl,
            "ccl_cmpl" => Jl_ccl_cmpl,
            "cmpl_cmpl" => Jl_cmpl_cmpl,
            "cmpl_cgdl" => Jl_cmpl_cgdl,
            "cgdl_cgdl" => Jl_cgdl_cgdl,
            "cgdl_cgc" => Jl_cgdl_cgc
        ),
        "J_lambda" => Dict(
            "acl_mem" => J_lambda_acl_mem,
            "mem_ccl" => J_lambda_mem_ccl
        ),
        "J_H2" => Dict(
            "agc_agdl" => J_H2_agc_agdl,
            "agdl_agdl" => J_H2_agdl_agdl,
            "agdl_ampl" => J_H2_agdl_ampl,
            "ampl_ampl" => J_H2_ampl_ampl,
            "ampl_acl" => J_H2_ampl_acl
        ),
        "J_O2" => Dict(
            "ccl_cmpl" => J_O2_ccl_cmpl,
            "cmpl_cmpl" => J_O2_cmpl_cmpl,
            "cmpl_cgdl" => J_O2_cmpl_cgdl,
            "cgdl_cgdl" => J_O2_cgdl_cgdl,
            "cgdl_cgc" => J_O2_cgdl_cgc
        ),
        "S_abs" => Dict(
            "v_acl" => Sv_abs_acl,
            "l_acl" => Sl_abs_acl,
            "v_ccl" => Sv_abs_ccl,
            "l_ccl" => Sl_abs_ccl
        ),
        "Sp" => Dict(
            "acl" => Sp_acl,
            "ccl" => Sp_ccl
        ),
        "S_H2" => Dict(
            "reac" => S_H2_reac,
            "cros" => S_H2_cros
        ),
        "S_O2" => Dict(
            "reac" => S_O2_reac,
            "cros" => S_O2_cros
        ),
        "Sv" => Dict(
            "agdl" => Sv_agdl,
            "ampl" => Sv_ampl,
            "acl" => Sv_acl,
            "ccl" => Sv_ccl,
            "cmpl" => Sv_cmpl,
            "cgdl" => Sv_cgdl
        ),
        "Sl" => Dict(
            "agdl" => Sl_agdl,
            "ampl" => Sl_ampl,
            "acl" => Sl_acl,
            "ccl" => Sl_ccl,
            "cmpl" => Sl_cmpl,
            "cgdl" => Sl_cgdl
        )
    )
end

