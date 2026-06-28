# -*- coding: utf-8 -*-

"""This file represents the calculation of the velocity over time. It is a component of the fuel cell model.
"""

# ________________________________________________________Velocity______________________________________________________

"""Build GC velocity/pressure profiles from inlet molar flows (in-place, zero allocation).

Parameters
----------
work : GCManifoldWorkspace
    Pre-allocated workspace.
sv : AbstractVector{<:CellState1D}
    Per-GC-node solver state.
J_a_in, J_c_in : Float64
    Inlet molar fluxes at the anode and cathode (mol·m⁻²·s⁻¹).
fc : AbstractFuelCell
    Fuel-cell instance providing geometric and operating parameters.
cfg : SimulationConfig
    Simulation configuration.
"""
function velocity_profiles_from_inlet_flows!(work::GCManifoldWorkspace,
                                             sv::AbstractVector{<:CellState1D},
                                             Ja_in::Float64,
                                             Jc_in::Float64,
                                             fc::AbstractFuelCell,
                                             cfg::SimulationConfig)

    # Extract parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    T_des, Pa_des, Pc_des = oc.T_des, oc.Pa_des, oc.Pc_des
    Hagc, Hcgc, Wagc, Wcgc = pp.Hagc, pp.Hcgc, pp.Wagc, pp.Wcgc
    Lgc, Ldist = pp.Lgc, pp.Ldist
    nb_gc, nb_gdl = np.nb_gc, np.nb_gdl
    L_node_gc = Lgc / nb_gc

    if cfg.type_auxiliary in (:forced_convective_cathode_with_anodic_recirculation,
                              :forced_convective_cathode_with_flow_through_anode)
        Pa_ext = Pext
        Pc_ext = Pext
    else
        Pa_ext = Pa_des
        Pc_ext = Pc_des
    end

    # ── Pass 1: Calculate GC thermodynamics, viscosities, and GC/GDL interface fluxes ────
    calculate_velocity_int_values!(work, sv, T_des, fc, cfg)

    # ── Pass 2: integrate molar flows from inlet to outlet (following flow direction) ──
    # `order[k]` is the physical GC index at flow step k; results stored by physical index.
    agc_order = anode_gc_order(nb_gc, cfg.type_flow)
    i_1   = agc_order[1]                               # physical index of the first GC node
    i_end = agc_order[end]                             # physical index of the last GC node
    @inbounds for k in 1:nb_gc
        i = agc_order[k]                               # physical index at numerical flow step k
        i_prev = k == 1 ? nothing : agc_order[k - 1]   # physical index at numerical flow step k-1
        Ja_prev = k == 1 ? Ja_in : work.J_a[i_prev]
        Jc_prev = k == 1 ? Jc_in : work.J_c[k - 1]     # cathode always co-flow (physical = flow order)
        work.J_a[i] = Ja_prev - work.J_tot_agc_agdl[i] * L_node_gc / Hagc
        work.J_c[k] = Jc_prev + work.J_tot_cgdl_cgc[k] * L_node_gc / Hcgc
    end

    # ── Pass 3: back-propagate pressure and velocity from outlet to inlet ────────
    # Outlet pressures and velocities.
    # Without GDL in the outlets, the stationary outlet flow is equal to the last GC node flow.
    Ja_out, Jc_out = work.J_a[i_end], work.J_c[end]
    Pa_out, Pc_out = Pa_ext, Pc_ext
    v_a_out, v_c_out = Ja_out / Pa_out * R * T_des, Jc_out / Pc_out * R * T_des

    # Last GC node pressure and velocity (with viscous drop along the last GC segment)
    work.P_a_chan[i_end] = Pa_out + 8 * π * work.mu_gaz_agc[i_end] * Ldist / (Hagc * Wagc) * v_a_out
    work.v_a[i_end]      = work.J_a[i_end] / work.P_a_chan[i_end] * R * T_des
    work.P_c_chan[end]   = Pc_out + 8 * π * work.mu_gaz_cgc[end] * Ldist / (Hcgc * Wcgc) * v_c_out
    work.v_c[end]        = work.J_c[end] / work.P_c_chan[end] * R * T_des

    # Back-propagate the pressure and velocity profiles from the last GC node to the first GC node.
    @inbounds for k in nb_gc:-1:2
        i      = agc_order[k]      # physical index at numerical flow step k
        i_prev = agc_order[k - 1]  # physical index at numerical flow step k-1
        work.P_a_chan[i_prev] = work.P_a_chan[i] + 8 * π * work.mu_gaz_agc[i] * L_node_gc / (Hagc * Wagc) *
                               (work.v_a[i] + work.J_tot_agc_agdl[i] / work.C_tot_agdl[i])
        work.v_a[i_prev]      = work.J_a[i_prev] / work.P_a_chan[i_prev] * R * T_des

        work.P_c_chan[k - 1] = work.P_c_chan[k] + 8 * π * work.mu_gaz_cgc[k] * L_node_gc / (Hcgc * Wcgc) *
                               (work.v_c[k] - work.J_tot_cgdl_cgc[k] / work.C_tot_cgdl[k])
        work.v_c[k - 1]      = work.J_c[k - 1] / work.P_c_chan[k - 1] * R * T_des
    end

    Pa_in = work.P_a_chan[i_1] + 8 * π * work.mu_gaz_agc[i_1] * L_node_gc / (Hagc * Wagc) *
             (work.v_a[i_1] + work.J_tot_agc_agdl[i_1] / work.C_tot_agdl[i_1])
    Pc_in = work.P_c_chan[1] + 8 * π * work.mu_gaz_cgc[1] * L_node_gc / (Hcgc * Wcgc) *
             (work.v_c[1] - work.J_tot_cgdl_cgc[1] / work.C_tot_cgdl[1])

    # work.v_a[i] is now indexed by physical position — no reordering needed.
    return work.v_a, work.v_c, Pa_in, Pc_in
end


"""Compute inlet-flow algebraic residuals for the DAE block in physical units."""
function velocity_inlet_flow_residuals!(gc_manifold_work::GCManifoldWorkspace,
                                        res::AbstractVector,
                                        J_a_in::Float64,
                                        J_c_in::Float64,
                                        sv::AbstractVector{<:CellState1D},
                                        i_fc_cell::Float64,
                                        fc::AbstractFuelCell,
                                        cfg::SimulationConfig)
    length(res) == 2 || throw(ArgumentError("res size mismatch in velocity_inlet_flow_residuals!."))
    cfg.type_auxiliary == :no_auxiliary ||
        throw(ArgumentError("velocity_inlet_flow_residuals! currently supports only :no_auxiliary."))

    _, _, P_a_in, P_c_in = velocity_profiles_from_inlet_flows!(gc_manifold_work, sv, J_a_in, J_c_in, fc, cfg)

    pp = fc.physical_parameters
    Hagc, Hcgc, Wagc, Wcgc = pp.Hagc, pp.Hcgc, pp.Wagc, pp.Wcgc
    nb_channel_in_gc, nb_cell = pp.nb_channel_in_gc, pp.nb_cell

    W_des = desired_flows(sv, i_fc_cell, P_a_in, P_c_in, fc, cfg)
    W_a_in = W_des.H2 + W_des.H2O_inj_a
    W_c_in = W_des.dry_air + W_des.H2O_inj_c
    J_a_in_calculated = W_a_in / (Hagc * Wagc) / nb_cell / nb_channel_in_gc
    J_c_in_calculated = W_c_in / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc

    res[1] = J_a_in_calculated - J_a_in
    res[2] = J_c_in_calculated - J_c_in
    return nothing
end


"""
    desired_flows(solver_variables, i_fc_cell, Pa_in, Pc_in, fc, cfg)

Calculate the desired flow for the air compressor and the humidifiers.

Parameters
----------
solver_variables : AbstractVector{<:CellState1D}
    Typed variables calculated by the solver at each gas-channel node.
i_fc_cell : Real
    Fuel cell current density (A.m-2).
Pa_in : Real
    Inlet pressure at the anode side (Pa).
Pc_in : Real
    Inlet pressure at the cathode side (Pa).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cfg : SimulationConfig
    Simulation configuration (provides `type_auxiliary`).

Returns
-------
DesiredInletFlows
    Desired hydrogen flow rate, dry-air flow rate, anode humidifier flow rate, and cathode humidifier flow rate.
"""
@inline function desired_flows(solver_variables::AbstractVector{<:CellState1D},
                       i_fc_cell::Float64,
                       Pa_in::Float64,
                       Pc_in::Float64,
                       fc::AbstractFuelCell,
                       cfg::SimulationConfig)::DesiredInletFlows

    # Extraction of the parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    T_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in = oc.T_des, oc.Sa, oc.Sc, oc.Phi_a_des, oc.Phi_c_des, oc.y_H2_in
    Hacl, Hmem, Hccl, Aact = pp.Hacl, pp.Hmem, pp.Hccl, pp.Aact
    kappa_co = pp.kappa_co
    nb_gc, nb_cell = np.nb_gc, pp.nb_cell

    # Physical quantities inside the stack
    #       The crossover current density i_n
    inv_H_total = 1.0 / (Hacl + Hmem + Hccl)
    w_acl = Hacl * inv_H_total
    w_mem = Hmem * inv_H_total
    w_ccl = Hccl * inv_H_total
    max_i_n = -Inf
    @inbounds for i in 1:nb_gc
        s_i = solver_variables[i]
        T_acl_i = s_i.acl.T
        T_mem_i = s_i.mem.T
        T_ccl_i = s_i.ccl.T
        lambda_mem_i = s_i.mem.lambda
        C_H2_acl_i = s_i.acl.C_H2
        C_O2_ccl_i = s_i.ccl.C_O2

        T_acl_mem_ccl = w_acl * T_acl_i + w_mem * T_mem_i + w_ccl * T_ccl_i
        i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl_i * k_H2(lambda_mem_i, T_mem_i, kappa_co)
        i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl_i * k_O2(lambda_mem_i, T_mem_i, kappa_co)
        i_n_i = i_H2 + i_O2
        max_i_n = max(max_i_n, i_n_i)
    end

    if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # The desired air compressor volume flow rate (mol.s-1). Warning: consider the minimum compressor flow!
        W_H2_des = 1 / y_H2_in * Sa * i_fc_cell / (2 * F) * (nb_cell * Aact)
        Wacp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc_cell, W_H2_des) # Adjust the desired compressor flow rate to ensure a minimum flow is maintained, based on the current density.
        W_air_ext_des = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * i_fc_cell / (4 * F) * (nb_cell * Aact)
        Wccp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc_cell, W_air_ext_des) # Adjust the desired compressor flow rate to ensure a minimum flow is maintained, based on the current density.

        # The desired humidifier volume flow rate at the anode side Wa_v_inj_des (mol.s-1). Warning: consider the minimum compressor flow!
        if cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
             Pasm = NaN
             Prd = Pasm
             W_H2_des        = 1 / y_H2_in * Sa * i_fc_cell / (2 * F) * (nb_cell * Aact)
             W_H2O_inj_a_des = Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des)) /
                               (1 - Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des))) * W_H2_des
        else  # type_auxiliary == "forced-convective_cathode_with_anodic_recirculation"
             W_H2O_inj_a_des = 0
        end

         # The desired humidifier volume flow rate at the cathode side Wc_v_inj_des (mol.s-1). Warning: consider the minimum compressor flow!
         Pcsm = NaN
         Wcp  = NaN
         Pcp  = Pcsm
         Wv_hum_in       = Phi_ext * Psat(Text) / Pext * Wcp  # Vapor flow rate from the outside.
         W_H20_c_des     = Phi_c_des * Psat(T_des) / Pcp * Wcp  # Desired vapor flow rate.
         W_H2O_inj_c_des = W_H20_c_des - Wv_hum_in  # Desired humidifier flow rate.

        _ = Wacp_des_adjusted
        _ = Wccp_des_adjusted
        # NOTE: When implementing desired_flows for auxiliary cases, remember that
        # :forced_convective_cathode_with_anodic_recirculation has no N2 at anode (pure H2 recirculation).
        throw(ErrorException("desired_flows is not yet implemented in Julia for the forced-convective auxiliary branches because the translated Python source still depends on unavailable variables such as Pasm, Pcsm, and Wcp."))
    else  # cfg.type_auxiliary == :no_auxiliary
        # At the anode side
        W_H2_des        = 1 / y_H2_in * Sa * (i_fc_cell + max_i_n) / (2 * F) * (nb_cell * Aact)
        W_H2O_inj_a_des = (Phi_a_des * Psat(T_des) / (Pa_in - Phi_a_des * Psat(T_des))) * W_H2_des

        # At the cathode side
        W_dry_air_des   = 1 / y_O2_ext * Sc * (i_fc_cell + max_i_n) / (4 * F) * (nb_cell * Aact)
        W_H2O_inj_c_des = (Phi_c_des * Psat(T_des) / (Pc_in - Phi_c_des * Psat(T_des))) * W_dry_air_des
    end

    return DesiredInletFlows(W_H2_des, W_dry_air_des, W_H2O_inj_a_des, W_H2O_inj_c_des)
end

