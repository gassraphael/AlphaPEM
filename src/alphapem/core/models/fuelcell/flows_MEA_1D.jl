# -*- coding: utf-8 -*-

"""This file represents all the matter flows inside the fuel cell system. It is a component of the fuel cell model.
"""

# ________________________________________________________Flows_________________________________________________________

"""In-place MEA flow computation with reusable workspace.

Parameters
----------
flows_work : MEAFlowsWorkspace
    Workspace for the flows. It is used to avoid memory reallocation.
flows_int_work : MEAFlowsIntWorkspace
    Workspace for the intermediate flows. It is used to avoid memory reallocation.
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
MEAFlows1D{NB_GDL, NB_MPL}
    Typed flows inside the fuel cell system.
"""
function calculate_flows_1D_MEA!(flows_work::MEAFlowsWorkspace,
                                 flows_int_work::MEAFlowsIntWorkspace,
                                 sv_1D::CellState1D{NB_GDL, NB_MPL},
                                 i_fc::Float64,
                                 v_a::Float64,
                                 v_c::Float64,
                                 fc::AbstractFuelCell,
                                 cfg::SimulationConfig)::MEAFlows1D{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of parameters
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    T_des = fc.operating_conditions.T_des
    Hmem, Hacl, Hccl = pp.Hmem, pp.Hacl, pp.Hccl
    Wagc, Wcgc, Hagc, Hcgc = pp.Wagc, pp.Wcgc, pp.Hagc, pp.Hcgc
    epsilon_gdl, epsilon_mpl = pp.epsilon_gdl, pp.epsilon_mpl
    kappa_co = pp.kappa_co
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Extraction of the variables
    C_v_agc, C_v_agdl = sv_1D.agc.C_v, getproperty.(sv_1D.agdl, :C_v)
    C_v_ampl, C_v_acl = getproperty.(sv_1D.ampl, :C_v), sv_1D.acl.C_v
    C_v_ccl, C_v_cmpl,  = sv_1D.ccl.C_v, getproperty.(sv_1D.cmpl, :C_v)
    C_v_cgdl, C_v_cgc = getproperty.(sv_1D.cgdl, :C_v), sv_1D.cgc.C_v

    s_agc, s_agdl = sv_1D.agc.s, getproperty.(sv_1D.agdl, :s)
    s_ampl, s_acl = getproperty.(sv_1D.ampl, :s), sv_1D.acl.s
    s_ccl, s_cmpl = sv_1D.ccl.s, getproperty.(sv_1D.cmpl, :s)
    s_cgdl, s_cgc = getproperty.(sv_1D.cgdl, :s), sv_1D.cgc.s

    T_agc, T_agdl = sv_1D.agc.T, getproperty.(sv_1D.agdl, :T)
    T_ampl, T_acl = getproperty.(sv_1D.ampl, :T), sv_1D.acl.T
    T_mem = sv_1D.mem.T
    T_ccl, T_cmpl = sv_1D.ccl.T, getproperty.(sv_1D.cmpl, :T)
    T_cgdl, T_cgc = getproperty.(sv_1D.cgdl, :T), sv_1D.cgc.T

    C_H2_agc, C_H2_agdl = sv_1D.agc.C_H2, getproperty.(sv_1D.agdl, :C_H2)
    C_H2_ampl, C_H2_acl = getproperty.(sv_1D.ampl, :C_H2), sv_1D.acl.C_H2

    C_O2_ccl, C_O2_cmpl = sv_1D.ccl.C_O2, getproperty.(sv_1D.cmpl, :C_O2)
    C_O2_cgdl, C_O2_cgc = getproperty.(sv_1D.cgdl, :C_O2), sv_1D.cgc.C_O2

    C_N2_agc, C_N2_agdl = sv_1D.agc.C_N2, getproperty.(sv_1D.agdl, :C_N2)
    C_N2_ampl, C_N2_acl = getproperty.(sv_1D.ampl, :C_N2), sv_1D.acl.C_N2
    C_N2_ccl, C_N2_cmpl = sv_1D.ccl.C_N2, getproperty.(sv_1D.cmpl, :C_N2)
    C_N2_cgdl, C_N2_cgc = getproperty.(sv_1D.cgdl, :C_N2), sv_1D.cgc.C_N2

    lambda_acl, lambda_mem, lambda_ccl = sv_1D.acl.lambda, sv_1D.mem.lambda, sv_1D.ccl.lambda

    # Intermediate values
    (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc, D_eff_EOD_acl_mem,
     D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem, D_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl,
     D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl,
     Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl,
     Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl) = calculate_flows_1D_MEA_int_values!(flows_int_work, sv_1D, i_fc, fc, cfg)

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
    Jl_agdl_agdl = flows_work.Jl_agdl_agdl
    @inbounds for i in 1:(nb_gdl - 1)
        Jl_agdl_agdl[i] = -D_cap_agdl_agdl[i] * d_dx(s_agdl[i], s_agdl[i + 1], H_gdl_node / 2)
    end
    Jl_agdl_ampl = -D_cap_agdl_ampl * d_dx(s_agdl[nb_gdl], s_ampl[1], H_gdl_node / 2, H_mpl_node / 2)
    Jl_ampl_ampl = flows_work.Jl_ampl_ampl
    @inbounds for i in 1:(nb_mpl - 1)
        Jl_ampl_ampl[i] = -D_cap_ampl_ampl[i] * d_dx(s_ampl[i], s_ampl[i + 1], H_mpl_node / 2)
    end
    Jl_ampl_acl = -D_cap_ampl_acl * d_dx(s_ampl[nb_mpl], s_acl, H_mpl_node / 2, Hacl / 2)

    # Cathode side
    Jl_ccl_cmpl = -D_cap_ccl_cmpl * d_dx(s_ccl, s_cmpl[1], Hccl / 2, H_mpl_node / 2)
    Jl_cmpl_cmpl = flows_work.Jl_cmpl_cmpl
    @inbounds for i in 1:(nb_mpl - 1)
        Jl_cmpl_cmpl[i] = -D_cap_cmpl_cmpl[i] * d_dx(s_cmpl[i], s_cmpl[i + 1], H_mpl_node / 2)
    end
    Jl_cmpl_cgdl = -D_cap_cmpl_cgdl * d_dx(s_cmpl[nb_mpl], s_cgdl[1], H_mpl_node / 2, H_gdl_node / 2)
    Jl_cgdl_cgdl = flows_work.Jl_cgdl_cgdl
    @inbounds for i in 1:(nb_gdl - 1)
        Jl_cgdl_cgdl[i] = -D_cap_cgdl_cgdl[i] * d_dx(s_cgdl[i], s_cgdl[i + 1], H_gdl_node / 2)
    end
    Jl_cgdl_cgc = theta_l_rem * epsilon_gdl * s_cgdl[nb_gdl] * max(Pcap_cgdl + rho_cgc * v_c^2 / 2, 0.0)

    # _____________________________________________Vapor flows (mol.m-2.s-1)____________________________________________

    # Conductive-convective vapor flows
    Jv_agc_agdl = h_a(Pagc, T_des, Wagc, Hagc) * (C_v_agc - C_v_agdl[1])  # Also calculated in velocity.jl
    Jv_cgdl_cgc = h_c(Pcgc, T_des, Wcgc, Hcgc) * (C_v_cgdl[nb_gdl] - C_v_cgc)  # Also calculated in velocity.jl

    # Conductive vapor flows
    #   Anode side
    Jv_agdl_agdl = flows_work.Jv_agdl_agdl
    @inbounds for i in 1:(nb_gdl - 1)
        Jv_agdl_agdl[i] = -Da_eff_agdl_agdl[i] * d_dx(C_v_agdl[i], C_v_agdl[i + 1], H_gdl_node / 2)
    end
    Jv_agdl_ampl = -Da_eff_agdl_ampl * d_dx(C_v_agdl[nb_gdl], C_v_ampl[1], H_gdl_node / 2, H_mpl_node / 2)
    Jv_ampl_ampl = flows_work.Jv_ampl_ampl
    @inbounds for i in 1:(nb_mpl - 1)
        Jv_ampl_ampl[i] = -Da_eff_ampl_ampl[i] * d_dx(C_v_ampl[i], C_v_ampl[i + 1], H_mpl_node / 2)
    end
    Jv_ampl_acl = -Da_eff_ampl_acl * d_dx(C_v_ampl[nb_mpl], C_v_acl, H_mpl_node / 2, Hacl / 2)

    #   Cathode side
    Jv_ccl_cmpl = -Dc_eff_ccl_cmpl * d_dx(C_v_ccl, C_v_cmpl[1], Hccl / 2, H_mpl_node / 2)
    Jv_cmpl_cmpl = flows_work.Jv_cmpl_cmpl
    @inbounds for i in 1:(nb_mpl - 1)
        Jv_cmpl_cmpl[i] = -Dc_eff_cmpl_cmpl[i] * d_dx(C_v_cmpl[i], C_v_cmpl[i + 1], H_mpl_node / 2)
    end
    Jv_cmpl_cgdl = -Dc_eff_cmpl_cgdl * d_dx(C_v_cmpl[nb_mpl], C_v_cgdl[1], H_mpl_node / 2, H_gdl_node / 2)
    Jv_cgdl_cgdl = flows_work.Jv_cgdl_cgdl
    @inbounds for i in 1:(nb_gdl - 1)
        Jv_cgdl_cgdl[i] = -Dc_eff_cgdl_cgdl[i] * d_dx(C_v_cgdl[i], C_v_cgdl[i + 1], H_gdl_node / 2)
    end

    # __________________________________H2, O2 and N2 flows (mol.m-2.s-1)____________________________________

    # Hydrogen and oxygen consumption
    #   Anode side
    S_H2_reac = i_fc / (2 * F * Hacl)
    S_H2_cros = R * T_acl_mem_ccl / (Hmem * Hacl) * (k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl +
                                                     2 * k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl)
    #   Cathode side
    S_O2_reac = i_fc / (4 * F * Hccl)
    S_O2_cros = R * T_acl_mem_ccl / (Hmem * Hccl) * (k_O2(lambda_mem, T_mem, kappa_co) * C_O2_ccl +
                                                     1 / 2 * k_H2(lambda_mem, T_mem, kappa_co) * C_H2_acl)

    # Conductive-convective H2, O2 and N2 flows
    J_H2_agc_agdl = h_a(Pagc, T_des, Wagc, Hagc) * (C_H2_agc - C_H2_agdl[1])  # Also calculated in velocity.jl
    J_O2_cgdl_cgc = h_c(Pcgc, T_des, Wcgc, Hcgc) * (C_O2_cgdl[nb_gdl] - C_O2_cgc)  # Also calculated in velocity.jl
    J_N2_agc_agdl = h_a(Pagc, T_des, Wagc, Hagc) * (C_N2_agc - C_N2_agdl[1])
    J_N2_cgdl_cgc = h_c(Pcgc, T_des, Wcgc, Hcgc) * (C_N2_cgdl[nb_gdl] - C_N2_cgc)

    # Conductive H2, O2 and N2 flows
    #   Anode side
    J_H2_agdl_agdl = flows_work.J_H2_agdl_agdl
    J_N2_agdl_agdl = flows_work.J_N2_agdl_agdl
    @inbounds for i in 1:(nb_gdl - 1)
        J_H2_agdl_agdl[i] = -Da_eff_agdl_agdl[i] * d_dx(C_H2_agdl[i], C_H2_agdl[i + 1], H_gdl_node / 2)
        J_N2_agdl_agdl[i] = -Da_eff_agdl_agdl[i] * d_dx(C_N2_agdl[i], C_N2_agdl[i + 1], H_gdl_node / 2)
    end
    J_H2_agdl_ampl = -Da_eff_agdl_ampl * d_dx(C_H2_agdl[nb_gdl], C_H2_ampl[1], H_gdl_node / 2, H_mpl_node / 2)
    J_N2_agdl_ampl = -Da_eff_agdl_ampl * d_dx(C_N2_agdl[nb_gdl], C_N2_ampl[1], H_gdl_node / 2, H_mpl_node / 2)
    J_H2_ampl_ampl = flows_work.J_H2_ampl_ampl
    J_N2_ampl_ampl = flows_work.J_N2_ampl_ampl
    @inbounds for i in 1:(nb_mpl - 1)
        J_H2_ampl_ampl[i] = -Da_eff_ampl_ampl[i] * d_dx(C_H2_ampl[i], C_H2_ampl[i + 1], H_mpl_node / 2)
        J_N2_ampl_ampl[i] = -Da_eff_ampl_ampl[i] * d_dx(C_N2_ampl[i], C_N2_ampl[i + 1], H_mpl_node / 2)
    end
    J_H2_ampl_acl = -Da_eff_ampl_acl * d_dx(C_H2_ampl[nb_mpl], C_H2_acl, H_mpl_node / 2, Hacl / 2)
    J_N2_ampl_acl = -Da_eff_ampl_acl * d_dx(C_N2_ampl[nb_mpl], C_N2_acl, H_mpl_node / 2, Hacl / 2)

    #   Cathode side
    J_O2_ccl_cmpl = -Dc_eff_ccl_cmpl * d_dx(C_O2_ccl, C_O2_cmpl[1], Hccl / 2, H_mpl_node / 2)
    J_N2_ccl_cmpl = -Dc_eff_ccl_cmpl * d_dx(C_N2_ccl, C_N2_cmpl[1], Hccl / 2, H_mpl_node / 2)
    J_O2_cmpl_cmpl = flows_work.J_O2_cmpl_cmpl
    J_N2_cmpl_cmpl = flows_work.J_N2_cmpl_cmpl
    @inbounds for i in 1:(nb_mpl - 1)
        J_O2_cmpl_cmpl[i] = -Dc_eff_cmpl_cmpl[i] * d_dx(C_O2_cmpl[i], C_O2_cmpl[i + 1], H_mpl_node / 2)
        J_N2_cmpl_cmpl[i] = -Dc_eff_cmpl_cmpl[i] * d_dx(C_N2_cmpl[i], C_N2_cmpl[i + 1], H_mpl_node / 2)
    end
    J_O2_cmpl_cgdl = -Dc_eff_cmpl_cgdl * d_dx(C_O2_cmpl[nb_mpl], C_O2_cgdl[1], H_mpl_node / 2, H_gdl_node / 2)
    J_N2_cmpl_cgdl = -Dc_eff_cmpl_cgdl * d_dx(C_N2_cmpl[nb_mpl], C_N2_cgdl[1], H_mpl_node / 2, H_gdl_node / 2)
    J_O2_cgdl_cgdl = flows_work.J_O2_cgdl_cgdl
    J_N2_cgdl_cgdl = flows_work.J_N2_cgdl_cgdl
    @inbounds for i in 1:(nb_gdl - 1)
        J_O2_cgdl_cgdl[i] = -Dc_eff_cgdl_cgdl[i] * d_dx(C_O2_cgdl[i], C_O2_cgdl[i + 1], H_gdl_node / 2)
        J_N2_cgdl_cgdl[i] = -Dc_eff_cgdl_cgdl[i] * d_dx(C_N2_cgdl[i], C_N2_cgdl[i + 1], H_gdl_node / 2)
    end

    # __________________________________________Water generated (mol.m-3.s-1)___________________________________________

    # Water produced in the membrane at the CL through the chemical reaction and crossover
    #   Anode side
    Sp_acl = 2 * k_O2(lambda_mem, T_mem, kappa_co) * R * T_acl_mem_ccl / (Hmem * Hacl) * C_O2_ccl
    #   Cathode side
    Sp_ccl = i_fc / (2 * F * Hccl) + k_H2(lambda_mem, T_mem, kappa_co) * R * T_acl_mem_ccl / (Hmem * Hccl) * C_H2_acl

    # Water absorption in the CL due to the contact between the ionomer and vapor or liquid water:
    #   Anode side
    Sv_abs_acl = (1 - s_acl) * gamma_sorp_v(C_v_acl, s_acl, lambda_acl, T_acl, Hacl) * rho_mem / M_eq *
                 (lambda_eq(C_v_acl, s_acl, T_acl) - lambda_acl)
    if s_acl > 0
        Sl_abs_acl = s_acl * gamma_sorp_l * rho_mem / M_eq * (lambda_eq(C_v_acl, s_acl, T_acl) - lambda_acl)
    else
        Sl_abs_acl = 0.0
    end

    #   Cathode side
    Sv_abs_ccl = (1 - s_ccl) * gamma_sorp_v(C_v_ccl, s_ccl, lambda_ccl, T_ccl, Hccl) * rho_mem / M_eq *
                 (lambda_eq(C_v_ccl, s_ccl, T_ccl) - lambda_ccl)
    if s_ccl > 0
        Sl_abs_ccl = s_ccl * gamma_sorp_l * rho_mem / M_eq * (lambda_eq(C_v_ccl, s_ccl, T_ccl) - lambda_ccl)
    else
        Sl_abs_ccl = 0.0
    end

    # Liquid water generated through vapor condensation or degenerated through evaporation
    #   Anode side
    Sl_agdl = flows_work.Sl_agdl
    @inbounds for i in 1:nb_gdl
        Sl_agdl[i] = Svl(:anode, s_agdl[i], C_v_agdl[i],
                         C_v_agdl[i] + C_H2_agdl[i] + C_N2_agdl[i],
                         T_agdl[i], epsilon_gdl)
    end
    Sl_ampl = flows_work.Sl_ampl
    @inbounds for i in 1:nb_mpl
        Sl_ampl[i] = Svl(:anode, s_ampl[i], C_v_ampl[i],
                         C_v_ampl[i] + C_H2_ampl[i] + C_N2_ampl[i],
                         T_ampl[i], epsilon_mpl)
    end
    Sl_acl = Svl(:anode, s_acl, C_v_acl, C_v_acl + C_H2_acl + C_N2_acl, T_acl,
                 epsilon_cl(lambda_acl, T_acl, Hacl))

    #   Cathode side
    Sl_ccl = Svl(:cathode, s_ccl, C_v_ccl, C_v_ccl + C_O2_ccl + C_N2_ccl, T_ccl,
                 epsilon_cl(lambda_ccl, T_ccl, Hccl))
    Sl_cmpl = flows_work.Sl_cmpl
    @inbounds for i in 1:nb_mpl
        Sl_cmpl[i] = Svl(:cathode, s_cmpl[i], C_v_cmpl[i],
                         C_v_cmpl[i] + C_O2_cmpl[i] + C_N2_cmpl[i],
                         T_cmpl[i], epsilon_mpl)
    end
    Sl_cgdl = flows_work.Sl_cgdl
    @inbounds for i in 1:nb_gdl
        Sl_cgdl[i] = Svl(:cathode, s_cgdl[i], C_v_cgdl[i],
                         C_v_cgdl[i] + C_O2_cgdl[i] + C_N2_cgdl[i],
                         T_cgdl[i], epsilon_gdl)
    end

    # Vapor generated through liquid water evaporation or degenerated through condensation
    #   Anode side
    Sv_agdl = flows_work.Sv_agdl
    @inbounds for i in 1:nb_gdl
        Sv_agdl[i] = -Sl_agdl[i]
    end
    Sv_ampl = flows_work.Sv_ampl
    @inbounds for i in 1:nb_mpl
        Sv_ampl[i] = -Sl_ampl[i]
    end
    Sv_acl = -Sl_acl

    #   Cathode side
    Sv_ccl = -Sl_ccl
    Sv_cmpl = flows_work.Sv_cmpl
    @inbounds for i in 1:nb_mpl
        Sv_cmpl[i] = -Sl_cmpl[i]
    end
    Sv_cgdl = flows_work.Sv_cgdl
    @inbounds for i in 1:nb_gdl
        Sv_cgdl[i] = -Sl_cgdl[i]
    end

    # ____________________________________Assemble and return typed flow container____________________________________
    Jv = MEAVaporFluxes{NB_GDL, NB_MPL}(Jv_agc_agdl, Jv_agdl_agdl, Jv_agdl_ampl,
                                        Jv_ampl_ampl, Jv_ampl_acl,
                                        Jv_ccl_cmpl, Jv_cmpl_cmpl, Jv_cmpl_cgdl,
                                        Jv_cgdl_cgdl, Jv_cgdl_cgc)
    Jl = MEALiquidFluxes{NB_GDL, NB_MPL}(Jl_agc_agdl, Jl_agdl_agdl, Jl_agdl_ampl,
                                         Jl_ampl_ampl, Jl_ampl_acl,
                                         Jl_ccl_cmpl, Jl_cmpl_cmpl, Jl_cmpl_cgdl,
                                         Jl_cgdl_cgdl, Jl_cgdl_cgc)
    J_lambda = MEADissolvedWaterFlux(J_lambda_acl_mem, J_lambda_mem_ccl)
    J_H2 = MEAHydrogenFluxes{NB_GDL, NB_MPL}(J_H2_agc_agdl, J_H2_agdl_agdl, J_H2_agdl_ampl,
                                             J_H2_ampl_ampl, J_H2_ampl_acl)
    J_O2 = MEAOxygenFluxes{NB_GDL, NB_MPL}(J_O2_ccl_cmpl, J_O2_cmpl_cmpl, J_O2_cmpl_cgdl,
                                           J_O2_cgdl_cgdl, J_O2_cgdl_cgc)
    J_N2 = MEANitrogenFluxes{NB_GDL, NB_MPL}(J_N2_agc_agdl, J_N2_agdl_agdl, J_N2_agdl_ampl,
                                             J_N2_ampl_ampl, J_N2_ampl_acl,
                                             J_N2_ccl_cmpl, J_N2_cmpl_cmpl, J_N2_cmpl_cgdl,
                                             J_N2_cgdl_cgdl, J_N2_cgdl_cgc)
    S_abs = MEASorptionSources(Sv_abs_acl, Sl_abs_acl, Sv_abs_ccl, Sl_abs_ccl)
    Sp = MEAWaterProductionSources(Sp_acl, Sp_ccl)
    S_H2 = MEAGasReactionSources(S_H2_reac, S_H2_cros)
    S_O2 = MEAGasReactionSources(S_O2_reac, S_O2_cros)
    Sv = MEAVaporSources{NB_GDL, NB_MPL}(Sv_agdl, Sv_ampl, Sv_acl, Sv_ccl, Sv_cmpl, Sv_cgdl)
    Sl = MEALiquidSources{NB_GDL, NB_MPL}(Sl_agdl, Sl_ampl, Sl_acl, Sl_ccl, Sl_cmpl, Sl_cgdl)

    return MEAFlows1D{NB_GDL, NB_MPL}(Jv, Jl, J_lambda, J_H2, J_O2, J_N2, S_abs, Sp, S_H2, S_O2, Sv, Sl)
end


