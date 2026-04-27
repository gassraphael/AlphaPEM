# -*- coding: utf-8 -*-

"""This file represents all the matter flows inside the fuel cell system. It is a component of the fuel cell model.
"""

# ________________________________________________________Flows_________________________________________________________

"""Calculate the flows inside the fuel cell system.

Parameters
----------
sv_1D : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
i_fc : Float64
    Fuel cell current density at time t (A.m-2).
v_a : Float64
    Anode gas velocity at time t (m.s-1).
v_c : Float64
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
function calculate_flows_1D_MEA(sv_1D::CellState1D,
                                i_fc::Float64,
                                v_a::Float64,
                                v_c::Float64,
                                fc::AbstractFuelCell,
                                cfg::SimulationConfig)

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of parameters
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    T_des = fc.operating_conditions.T_des
    Hmem, Hacl, Hccl = pp.Hmem, pp.Hacl, pp.Hccl
    Wagc, Wcgc, Hagc, Hcgc = pp.Wagc, pp.Wcgc, pp.Hagc, pp.Hcgc
    epsilon_gdl, epsilon_mpl = pp.epsilon_gdl, pp.epsilon_mpl
    K_l_ads, kappa_co = pp.K_l_ads, pp.kappa_co
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv_1D.agc.C_v, sv_1D.acl.C_v, sv_1D.ccl.C_v, sv_1D.cgc.C_v
    C_v_agdl = [sv_1D.agdl[i].C_v for i in 1:nb_gdl]
    C_v_ampl = [sv_1D.ampl[i].C_v for i in 1:nb_mpl]
    C_v_cmpl = [sv_1D.cmpl[i].C_v for i in 1:nb_mpl]
    C_v_cgdl = [sv_1D.cgdl[i].C_v for i in 1:nb_gdl]
    s_acl, s_ccl = sv_1D.acl.s, sv_1D.ccl.s
    s_agdl = [sv_1D.agdl[i].s for i in 1:nb_gdl]
    s_ampl = [sv_1D.ampl[i].s for i in 1:nb_mpl]
    s_cmpl = [sv_1D.cmpl[i].s for i in 1:nb_mpl]
    s_cgdl = [sv_1D.cgdl[i].s for i in 1:nb_gdl]
    lambda_acl, lambda_mem, lambda_ccl = sv_1D.acl.lambda, sv_1D.mem.lambda, sv_1D.ccl.lambda
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv_1D.agc.C_H2, sv_1D.acl.C_H2, sv_1D.ccl.C_O2, sv_1D.cgc.C_O2
    C_H2_agdl = [sv_1D.agdl[i].C_H2 for i in 1:nb_gdl]
    C_H2_ampl = [sv_1D.ampl[i].C_H2 for i in 1:nb_mpl]
    C_O2_cmpl = [sv_1D.cmpl[i].C_O2 for i in 1:nb_mpl]
    C_O2_cgdl = [sv_1D.cgdl[i].C_O2 for i in 1:nb_gdl]
    C_N2_agc, C_N2_cgc = sv_1D.agc.C_N2, sv_1D.cgc.C_N2
    T_acl, T_mem, T_ccl = sv_1D.acl.T, sv_1D.mem.T, sv_1D.ccl.T
    T_agdl = [sv_1D.agdl[i].T for i in 1:nb_gdl]
    T_ampl = [sv_1D.ampl[i].T for i in 1:nb_mpl]
    T_cmpl = [sv_1D.cmpl[i].T for i in 1:nb_mpl]
    T_cgdl = [sv_1D.cgdl[i].T for i in 1:nb_gdl]

    # Intermediate values
    (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc, D_eff_EOD_acl_mem,
     D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem, D_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl,
     D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl,
     Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl,
     Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl) = flows_1D_MEA_int_values(sv_1D, i_fc, fc, cfg)

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

    # ____________________________________Assemble and return typed flow container____________________________________
    Jv = MEAVaporFluxes{nb_gdl, nb_mpl}(Jv_agc_agdl, Jv_agdl_agdl, Jv_agdl_ampl,
                                        Jv_ampl_ampl, Jv_ampl_acl,
                                        Jv_ccl_cmpl, Jv_cmpl_cmpl, Jv_cmpl_cgdl,
                                        Jv_cgdl_cgdl, Jv_cgdl_cgc)
    Jl = MEALiquidFluxes{nb_gdl, nb_mpl}(Jl_agc_agdl, Jl_agdl_agdl, Jl_agdl_ampl,
                                         Jl_ampl_ampl, Jl_ampl_acl,
                                         Jl_ccl_cmpl, Jl_cmpl_cmpl, Jl_cmpl_cgdl,
                                         Jl_cgdl_cgdl, Jl_cgdl_cgc)
    J_lambda = MEADissolvedWaterFlux(J_lambda_acl_mem, J_lambda_mem_ccl)
    J_H2 = MEAHydrogenFluxes{nb_gdl, nb_mpl}(J_H2_agc_agdl, J_H2_agdl_agdl, J_H2_agdl_ampl,
                                             J_H2_ampl_ampl, J_H2_ampl_acl)
    J_O2 = MEAOxygenFluxes{nb_gdl, nb_mpl}(J_O2_ccl_cmpl, J_O2_cmpl_cmpl, J_O2_cmpl_cgdl,
                                           J_O2_cgdl_cgdl, J_O2_cgdl_cgc)
    S_abs = MEASorptionSources(Sv_abs_acl, Sl_abs_acl, Sv_abs_ccl, Sl_abs_ccl)
    Sp = MEAWaterProductionSources(Sp_acl, Sp_ccl)
    S_H2 = MEAGasReactionSources(S_H2_reac, S_H2_cros)
    S_O2 = MEAGasReactionSources(S_O2_reac, S_O2_cros)
    Sv = MEAVaporSources{nb_gdl, nb_mpl}(Sv_agdl, Sv_ampl, Sv_acl, Sv_ccl, Sv_cmpl, Sv_cgdl)
    Sl = MEALiquidSources{nb_gdl, nb_mpl}(Sl_agdl, Sl_ampl, Sl_acl, Sl_ccl, Sl_cmpl, Sl_cgdl)

    return MEAFlows1D{nb_gdl, nb_mpl}(Jv, Jl, J_lambda, J_H2, J_O2, S_abs, Sp, S_H2, S_O2, Sv, Sl)
end

