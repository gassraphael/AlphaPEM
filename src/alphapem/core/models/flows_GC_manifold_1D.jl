# -*- coding: utf-8 -*-

"""This file represents all the flows passing through the auxiliaries. It is a component of the fuel cell model.
"""

# ______________________________________________________Auxiliaries_____________________________________________________

"""Calculate the flows passing through the auxiliaries.

Parameters
----------
sv_1D_cell : AbstractVector{<:CellState1D}
    Typed variables calculated by the solver (cell internal states).
sv_1D_manifold
    Typed variables calculated by the solver (manifold internal states).
sv_auxiliary
    Typed variables calculated by the solver (auxiliary internal states).
i_fc_cell : Float64
    Fuel cell current density at time t (A.m-2).
v_a : Vector{Float64}
    Velocity evolution at the anode side (m.s-1).
v_c : Vector{Float64}
    Velocity evolution at the cathode side (m.s-1).
Pa_in : Float64
    Inlet pressure at the anode side (Pa).
Pc_in : Float64
    Inlet pressure at the cathode side (Pa).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cfg : SimulationConfig
    Simulation configuration (provides `type_auxiliary`).

Returns
-------
GCManifoldFlows1D{NB_GC}
    Typed global and species-specific flows in the gas channels and auxiliaries.
"""
function calculate_flows_1D_GC_manifold(sv_1D_cell::AbstractVector{<:CellState1D},
                                        sv_1D_manifold,
                                        sv_auxiliary,
                                        i_fc_cell::Float64,
                                        v_a::Vector{Float64},
                                        v_c::Vector{Float64},
                                        Pa_in::Float64,
                                        Pc_in::Float64,
                                        fc::AbstractFuelCell,
                                        cfg::SimulationConfig)

    # __________________________________________________Preliminaries___________________________________________________

    # Extraction of the parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des, Phi_a_des, Phi_c_des, Sa, Sc, y_H2_in = oc.T_des, oc.Phi_a_des, oc.Phi_c_des, oc.Sa, oc.Sc, oc.y_H2_in
    Aact, nb_cell, Hagc, Hcgc = pp.Aact, pp.nb_cell, pp.Hagc, pp.Hcgc
    Wagc, Wcgc, Lgc, nb_channel_in_gc, A_T_a, A_T_c = pp.Wagc, pp.Wcgc, pp.Lgc, pp.nb_channel_in_gc, pp.A_T_a, pp.A_T_c
    nb_gc = np.nb_gc
    type_auxiliary = cfg.type_auxiliary
    agc_order = anode_gc_order(nb_gc, cfg.type_flow)

    # Extraction of the variables
    C_v_agc = [sv_1D_cell[i].agc.C_v for i in agc_order]
    C_v_cgc = [sv_1D_cell[i].cgc.C_v for i in 1:nb_gc]
    s_agc = [sv_1D_cell[i].agc.s for i in agc_order]
    s_cgc = [sv_1D_cell[i].cgc.s for i in 1:nb_gc]
    C_H2_agc = [sv_1D_cell[i].agc.C_H2 for i in agc_order]
    C_O2_cgc = [sv_1D_cell[i].cgc.C_O2 for i in 1:nb_gc]
    C_N2_cgc = [sv_1D_cell[i].cgc.C_N2 for i in 1:nb_gc]
    v_a_o = [v_a[i] for i in agc_order] # velocity oriented in the anode flow direction

    # Intermediate values
    (P_agc, P_cgc, Phi_agc, Phi_cgc, y_H2_agc, y_O2_cgc, M_agc, M_cgc, M_ext, M_H2_N2_in, rho_agc, rho_cgc, k_purge,
     Abp_a, Abp_c, mu_gaz_agc, mu_gaz_cgc) = flow_1D_GC_manifold_int_values(sv_1D_cell, sv_auxiliary, fc, cfg)
    W_des = desired_flows(sv_1D_cell, i_fc_cell, Pa_in, Pc_in, fc, cfg)

    # _________________________________________Inlet and outlet global flows____________________________________________
    """Global flows here refer to flows that integrate all the chemical species circulating together.
    Slight differences are to be noted in the expression of these flows depending on the type of auxiliary selected.
    """

    # Anode flow through the auxiliaries in mol.s-1
    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # pass
        # Wa_in = rho_asm_in_to_asm * v_a * A_T_a
        # Wasm_to_asm_out = rho_asm_to_asm_out * v_a * Hagc * Wagc
        # Wasm_out_to_agc = rho_asm_out_to_agc * v_a * Hagc * Wagc
        # Wagc_to_aem_in = rho_agc_to_aem_in * v_a * Hagc * Wagc
        # Waem_in_to_aem = rho_aem_in_to_aem * v_a * Hagc * Wagc
        # if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation # Attention: include a minimal flow for the pump, as for incoming flows.
        #     Ware = Maem_out_re * (Paem_out_re / (Paem_out_re - Phi_aem_out_re * Psat(T_des))) *
        #            (Sa - 1) * i_fc_cell / (2 * F) * (nb_cell * Aact)  # The pump exactly compensates the pressure drop.
        #     Wasm_in_re_to_asm = rho_asm_in_re_to_asm * v_a * A_T_a
        #     Waem_to_aem_out_re = rho_aem_to_aem_out_re * v_a * A_T_a
        #     Waem_to_aem_out = k_purge * rho_aem_to_aem_out * v_a * A_T_a
        #     Wa_out = k_purge * rho_aem_out_to_ext * v_a * A_T_a
        # else # type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        #     Waem_to_aem_out = rho_aem_to_aem_out * v_a * Abp_a
        #     Wa_out = rho_aem_out_to_ext * v_a * Abp_a
    else  # type_auxiliary == :no_auxiliary (only 1 cell)
        Wa_in = W_des.H2 + W_des.H2O_inj_a  # This expression is also present in calculate_velocity_evolution.
        Wa_out = P_agc[agc_order[end]] / (R * T_des) * v_a_o[end] * Hagc * Wagc * nb_cell * nb_channel_in_gc
    end

    # Anode flow entering/leaving the stack in mol.m-2.s-1
    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # pass
        # Ja_in = 0
        # Ja_out = 0
    else  # type_auxiliary == :no_auxiliary (only 1 cell)
        Ja_in = Wa_in / (Hagc * Wagc) / nb_cell / nb_channel_in_gc  # This expression is also present in calculate_velocity_evolution.
        Ja_out = Wa_out / (Hagc * Wagc) / nb_cell / nb_channel_in_gc
    end

    # Cathode flow through the auxiliaries in mol.s-1
    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # pass
        # Wc_in = rho_csm_in_to_csm * v_c * A_T_c
        # Wcsm_to_csm_out = rho_csm_to_csm_out * v_c * Hcgc * Wcgc
        # Wcsm_out_to_cgc = rho_csm_out_to_cgc * v_c * Hcgc * Wcgc
        # Wcgc_to_cem_in = rho_cgc_to_cem_in * v_c * Hcgc * Wcgc
        # Wcem_in_to_cem = rho_cem_in_to_cem * v_c * Hcgc * Wcgc
        # Wcem_to_cem_out = rho_cem_to_cem_out * v_c * Abp_c
        # Wc_out = rho_cem_out_to_ext * v_c * Abp_c
    else  # type_auxiliary == :no_auxiliary (only 1 cell)
        Wc_in = W_des.dry_air + W_des.H2O_inj_c  # This expression is also present in calculate_velocity_evolution.
        Wc_out = P_cgc[nb_gc] / (R * T_des) * v_c[nb_gc] * Hcgc * Wcgc * nb_cell * nb_channel_in_gc
    end

    # Cathode flow entering/leaving the stack in mol.m-2.s-1
    Jc_in = Wc_in / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc  # This expression is also present in calculate_velocity_evolution.
    Jc_out = Wc_out / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc

    # ________________________________________Inlet and outlet specific flows___________________________________________
    """Specific flows here refer to flows that integrate only a single chemical species within the ensemble of species
    circulating together. For example, only the water vapor flow within the ensemble of hydrogen and water vapor.
    """

    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # Jv_agc_in = Phi_asm * Psat(T_des) / Pasm * Ja_in
    else  # type_auxiliary == :no_auxiliary
        Jv_agc_in = Phi_a_des * Psat(T_des) / Pa_in * Ja_in
    end
    Jv_agc_agc = [C_v_agc[i] * v_a_o[i] for i in 1:(nb_gc - 1)]
    Jv_agc_out = C_v_agc[end] * R * T_des / P_agc[agc_order[end]] * Ja_out

    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # Jv_cgc_in = Phi_csm * Psat(T_des) / Pcsm * Jc_in
    else  # type_auxiliary == :no_auxiliary
        Jv_cgc_in = Phi_c_des * Psat(T_des) / Pc_in * Jc_in
    end
    Jv_cgc_cgc = [C_v_cgc[i] * v_c[i] for i in 1:(nb_gc - 1)]
    Jv_cgc_out = C_v_cgc[nb_gc] * R * T_des / P_cgc[nb_gc] * Jc_out

    # Liquid water flows at the GC (kg.m-2.s-1)
    #   At the anode side
    s_agc_outlet = 0.0  # Boundary condition at the outlet of the anode GC: no liquid water at the outlet.
    Jl_agc_agc_conv = [rho_H2O_l(T_des) * K_v_liq_gas * v_a_o[i] * s_agc[i] for i in 1:(nb_gc - 1)]
    Jl_agc_agc_dif = [-D_liq_dif * d_dx(s_agc[i], i == nb_gc ? s_agc_outlet : s_agc[i + 1], (Lgc / nb_gc) / 2)
                      for i in 1:(nb_gc - 1)]
    Jl_agc_agc = [Jl_agc_agc_conv[i] + Jl_agc_agc_dif[i] for i in 1:(nb_gc - 1)]
    Jl_agc_agc_conv_out = rho_H2O_l(T_des) * K_v_liq_gas * v_a_o[end] * s_agc[end]
    Jl_agc_agc_dif_out = -D_liq_dif * d_dx(s_agc[end], s_agc_outlet, (Lgc / nb_gc) / 2)
    Jl_agc_out = Jl_agc_agc_conv_out + Jl_agc_agc_dif_out

    #   At the cathode side
    s_cgc_outlet = 0.0  # Boundary condition at the outlet of the cathode GC: no liquid water at the outlet.
    Jl_cgc_cgc_conv = [rho_H2O_l(T_des) * K_v_liq_gas * v_c[i] * s_cgc[i] for i in 1:(nb_gc - 1)]
    Jl_cgc_cgc_dif = [-D_liq_dif * d_dx(s_cgc[i], i == nb_gc ? s_cgc_outlet : s_cgc[i + 1], (Lgc / nb_gc) / 2)
                      for i in 1:(nb_gc - 1)]
    Jl_cgc_cgc = [Jl_cgc_cgc_conv[i] + Jl_cgc_cgc_dif[i] for i in 1:(nb_gc - 1)]
    Jl_cgc_cgc_conv_out = rho_H2O_l(T_des) * K_v_liq_gas * v_c[nb_gc] * s_cgc[nb_gc]
    Jl_cgc_cgc_dif_out = -D_liq_dif * d_dx(s_cgc[nb_gc], s_cgc_outlet, (Lgc / nb_gc) / 2)
    Jl_cgc_out = Jl_cgc_cgc_conv_out + Jl_cgc_cgc_dif_out

    # H2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # pass
        # J_H2_agc_in = y_H2["asm_out"] * (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        # J_H2_agc_out = y_H2_agc * (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out
    else  # type_auxiliary == :forced_convective_cathode_with_anodic_recirculation or type_auxiliary == :no_auxiliary
        J_H2_agc_in = (1 - Phi_a_des * Psat(T_des) / Pa_in) * Ja_in
        J_H2_agc_agc = [C_H2_agc[i] * v_a_o[i] for i in 1:(nb_gc - 1)]
        J_H2_agc_out = C_H2_agc[end] * R * T_des / P_agc[agc_order[end]] * Ja_out
    end

    # O2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # pass
        # J_O2_cgc_in = y_O2_csm * (1 - Phi_csm * Psat(T_des) / Pcsm) * Jc_in
    else  # type_auxiliary == :no_auxiliary
        J_O2_cgc_in = y_O2_ext * (1 - Phi_c_des * Psat(T_des) / Pc_in) * Jc_in
    end
    J_O2_cgc_cgc = [C_O2_cgc[i] * v_c[i] for i in 1:(nb_gc - 1)]
    J_O2_cgc_out = C_O2_cgc[nb_gc] * R * T_des / P_cgc[nb_gc] * Jc_out

    # N2 flows at the GC (mol.m-2.s-1)
    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # pass
        # J_N2_agc_in = (1 - y_H2["asm_out"]) * (1 - Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc) * Ja_in
        # J_N2_agc_out = (1 - y_H2_agc) * (1 - Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in) * Ja_out
        # J_N2_cgc_in = (1 - y_O2_csm_out_to_cgc) * (1 - Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc) * Jc_in
        # J_N2_cgc_out = (1 - y_O2_cgc_to_cem_in) * (1 - Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in) * Jc_out
    else  # type_auxiliary == :no_auxiliary
        J_N2_agc_in = 0.0
        J_N2_agc_agc = fill(0.0, max(nb_gc - 1, 0))
        J_N2_agc_out = 0.0
        J_N2_cgc_in = (1 - y_O2_ext) * (1 - Phi_c_des * Psat(T_des) / Pc_in) * Jc_in
        J_N2_cgc_cgc = [C_N2_cgc[i] * v_c[i] for i in 1:(nb_gc - 1)]
        J_N2_cgc_out = C_N2_cgc[nb_gc] * R * T_des / P_cgc[nb_gc] * Jc_out
    end

    # Vapor flows at the manifold (mol.s-1)
    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # pass
        # Wv_asm_in_to_asm = Phi_asm_in_to_asm * Psat(T_des) / Pasm_in_to_asm * Wa_in
        # Wv_asm_to_asm_out = Phi_asm_to_asm_out * Psat(T_des) / Pasm_to_asm_out * Wasm_to_asm_out
        # Wv_asm_out_to_agc = Phi_asm_out_to_agc * Psat(T_des) / Pasm_out_to_agc * Wasm_out_to_agc
        # Wv_agc_to_aem_in = Phi_agc_to_aem_in * Psat(T_des) / Pagc_to_aem_in * Wagc_to_aem_in
        # Wv_aem_in_to_aem = Phi_aem_in_to_aem * Psat(T_des) / Paem_in_to_aem * Waem_in_to_aem
        # Wv_aem_to_aem_out = Phi_aem_to_aem_out * Psat(T_des) / Paem_to_aem_out * Waem_to_aem_out
        # Wv_a_out = Phi_aem_out * Psat(T_des) / Paem_out * Wa_out
        # if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
        #     # At the anode side
        #     Wv_asm_ext_to_in = 0
        #     Wv_asm_in_re_to_asm = Phi_asm_in_re_to_asm * Psat(T_des) / Pasm_in_re_to_asm * Wasm_in_re_to_asm
        #     Wv_aem_to_aem_out_re = Phi_aem_to_aem_out_re * Psat(T_des) / Paem_to_aem_out_re * Waem_to_aem_out_re
        #     Wv_are = Phi_aem_out_re * Psat(T_des) / Paem_out_re * (Ware / M["aem_out_re"])  # The pump exactly compensates the pressure drop.
        # else # type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        #     # At the anode side
        #     Wv_asm_ext_to_in = Wa_inj / M_H2O
        # # At the cathode side
        # Wv_csm_ext_to_in = Phi_ext * Psat(Text) / Pext * (Wcp / M["ext"]) + Wc_inj / M_H2O
        # Wv_csm_in_to_csm = Phi_csm_in_to_csm * Psat(T_des) / Pcsm_in_to_csm * Wc_in
        # Wv_csm_to_csm_out = Phi_csm_to_csm_out * Psat(T_des) / Pcsm_to_csm_out * Wcsm_to_csm_out
        # Wv_csm_out_to_cgc = Phi_csm_out_to_cgc * Psat(T_des) / Pcsm_out_to_cgc * Wcsm_out_to_cgc
        # Wv_cgc_to_cem_in = Phi_cgc_to_cem_in * Psat(T_des) / Pcgc_to_cem_in * Wcgc_to_cem_in
        # Wv_cem_in_to_cem = Phi_cem_in_to_cem * Psat(T_des) / Pcem_in_to_cem * Wcem_in_to_cem
        # Wv_cem_to_cem_out = Phi_cem_to_cem_out * Psat(T_des) / Pcem_to_cem_out * Wcem_to_cem_out
        # Wv_c_out = Phi_cem_out * Psat(T_des) / Pcem_out * Wc_out
    end

    if type_auxiliary == :no_auxiliary
        Jv = GCVaporFlows{nb_gc}(Jv_agc_in, Jv_agc_agc, Jv_agc_out,
                                 Jv_cgc_in, Jv_cgc_cgc, Jv_cgc_out)
        Jl = GCLiquidFlows{nb_gc}(Jl_agc_agc, Jl_agc_out, Jl_cgc_cgc, Jl_cgc_out)
        J_H2 = GCHydrogenFlows{nb_gc}(J_H2_agc_in, J_H2_agc_agc, J_H2_agc_out)
        J_O2 = GCOxygenFlows{nb_gc}(J_O2_cgc_in, J_O2_cgc_cgc, J_O2_cgc_out)
        J_N2 = GCNitrogenFlows{nb_gc}(J_N2_agc_in, J_N2_agc_agc, J_N2_agc_out,
                                      J_N2_cgc_in, J_N2_cgc_cgc, J_N2_cgc_out)
        W = GCMassFlows(Wa_in, Wa_out, Wc_in, Wc_out)
        return GCManifoldFlows1D{nb_gc}(Jv, Jl, J_H2, J_O2, J_N2, W)
    end
end