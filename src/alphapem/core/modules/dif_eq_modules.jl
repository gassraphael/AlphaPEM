# -*- coding: utf-8 -*-

"""This module is used to determine intermediate values for the calculation of the differential equations
and to implement integration events.
"""

# ____________________________________________Differential equations modules____________________________________________

"""Calculate typed intermediate values used by differential equations."""
function calculate_dif_eq_int_values(t::Float64,
                                     sv,
                                     fc::AbstractFuelCell,
                                     cfg::SimulationConfig,
                                     sv_manifold=nothing,
                                     sv_auxiliary=nothing,
                                     Ware=nothing)

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv.agc.C_v, sv.acl.C_v, sv.ccl.C_v, sv.cgc.C_v
    s_acl, s_ccl = sv.acl.s, sv.ccl.s
    lambda_acl, lambda_mem, lambda_ccl = sv.acl.lambda, sv.mem.lambda, sv.ccl.lambda
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv.agc.C_H2, sv.acl.C_H2, sv.ccl.C_O2, sv.cgc.C_O2
    C_N2_agc, C_N2_cgc = sv.agc.C_N2, sv.cgc.C_N2
    T_agc, T_acl, T_mem, T_ccl, T_cgc = sv.agc.T, sv.acl.T, sv.mem.T, sv.ccl.T, sv.cgc.T

    # Extraction of the parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des, y_H2_in = oc.T_des, oc.y_H2_in
    Hmem, Hacl, Hccl = pp.Hmem, pp.Hacl, pp.Hccl
    A_T_a, A_T_c = pp.A_T_a, pp.A_T_c
    epsilon_gdl, epsilon_mpl, kappa_co = pp.epsilon_gdl, pp.epsilon_mpl, pp.kappa_co
    purge_time, delta_purge = np.purge_time, np.delta_purge
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Physical quantities inside the stack
    #       Pressures
    P_agc = (C_v_agc + C_H2_agc + C_N2_agc) * R * T_agc
    P_cgc = (C_v_cgc + C_O2_cgc + C_N2_cgc) * R * T_cgc

    #       Total concentration
    C_tot_agc = C_v_agc + C_H2_agc + C_N2_agc
    C_tot_cgc = C_v_cgc + C_O2_cgc + C_N2_cgc

    #     H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2).
    y_O2_cgc = C_O2_cgc / (C_O2_cgc + C_N2_cgc)

    #       Molar masses
    M_agc = C_v_agc * R * T_des / P_agc * M_H2O +
            C_H2_agc * R * T_des / P_agc * M_H2 +
            C_N2_agc * R * T_des / P_agc * M_N2
    M_cgc = C_v_cgc * R * T_des / P_cgc * M_H2O +
            C_O2_cgc * R * T_des / P_cgc * M_O2 +
            C_N2_cgc * R * T_des / P_cgc * M_N2

    #       Density of the gas mixture.
    rho_agc = P_agc / (R * T_agc) * M_agc
    rho_cgc = P_cgc / (R * T_cgc) * M_cgc

    #       Vapor ratio over the gas mixture.
    x_H2O_v_agc = C_v_agc / C_tot_agc
    x_H2O_v_cgc = C_v_cgc / C_tot_cgc

    #       Dynamic viscosity of the gas mixture.
    mu_gaz_agc = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_agc, 1 - x_H2O_v_agc], T_agc)
    mu_gaz_cgc = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                  [x_H2O_v_cgc, y_O2_cgc * (1 - x_H2O_v_cgc), (1 - y_O2_cgc) * (1 - x_H2O_v_cgc)],
                                  T_cgc)

    #       Volumetric heat capacity (J.m-3.K-1)
    rho_Cp0_agdl = ntuple(i -> calculate_rho_Cp0("agdl", sv.agdl[i].T, sv.agdl[i].C_v,
                                                 sv.agdl[i].s, nothing, sv.agdl[i].C_H2, nothing, C_N2_agc,
                                                 epsilon_gdl), nb_gdl)
    rho_Cp0_ampl = ntuple(i -> calculate_rho_Cp0("ampl", sv.ampl[i].T, sv.ampl[i].C_v,
                                                 sv.ampl[i].s, nothing, sv.ampl[i].C_H2, nothing, C_N2_agc,
                                                 epsilon_mpl), nb_mpl)
    rho_Cp0_acl = calculate_rho_Cp0("acl", T_acl, C_v_acl, s_acl, lambda_acl, C_H2_acl, nothing, C_N2_agc,
                                    nothing, Hacl)
    rho_Cp0_mem = calculate_rho_Cp0("mem", T_mem, nothing, nothing, lambda_mem)
    rho_Cp0_ccl = calculate_rho_Cp0("ccl", T_ccl, C_v_ccl, s_ccl, lambda_ccl, nothing, C_O2_ccl,
                                    C_N2_cgc, nothing, Hccl)
    rho_Cp0_cmpl = ntuple(i -> calculate_rho_Cp0("cmpl", sv.cmpl[i].T, sv.cmpl[i].C_v,
                                                 sv.cmpl[i].s, nothing, nothing, sv.cmpl[i].C_O2, C_N2_cgc,
                                                 epsilon_mpl), nb_mpl)
    rho_Cp0_cgdl = ntuple(i -> calculate_rho_Cp0("cgdl", sv.cgdl[i].T, sv.cgdl[i].C_v,
                                                 sv.cgdl[i].s, nothing, nothing, sv.cgdl[i].C_O2, C_N2_cgc,
                                                 epsilon_gdl), nb_gdl)
    rho_Cp0 = MEAThermalIntermediates{nb_gdl, nb_mpl}(rho_Cp0_agdl, rho_Cp0_ampl, rho_Cp0_acl,
                                                       rho_Cp0_mem, rho_Cp0_ccl, rho_Cp0_cmpl, rho_Cp0_cgdl)

    #       Crossover current density
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                            [Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
    i_n = i_H2 + i_O2

    # Auxiliary/manifold branch now uses typed structures as well.
    v_re = nothing
    k_purge = nothing
    if cfg.type_auxiliary in (:forced_convective_cathode_with_anodic_recirculation,
                              :forced_convective_cathode_with_flow_through_anode) &&
       sv_manifold !== nothing && sv_auxiliary !== nothing

        # Extraction of the variables
        Pasm, Paem = sv_manifold.asm.nodes[1].P, sv_manifold.aem.nodes[1].P
        Pcsm, Pcem = sv_manifold.csm.nodes[1].P, sv_manifold.cem.nodes[1].P
        Phi_asm, Phi_aem = sv_manifold.asm.nodes[1].Phi, sv_manifold.aem.nodes[1].Phi
        Phi_csm, Phi_cem = sv_manifold.csm.nodes[1].Phi, sv_manifold.cem.nodes[1].Phi

        if cfg.type_purge == "no_purge"
            k_purge = 0.0
        elseif cfg.type_purge == "constant_purge"
            k_purge = 1.0
        elseif cfg.type_purge == "periodic_purge"
            period = purge_time + delta_purge
            phase = t - floor(Int, t / period) * period
            k_purge = phase <= purge_time ? 1.0 : 0.0
        else
            throw(ArgumentError("The type_purge variable should be correctly referenced."))
        end



        # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the exhaust manifolds
        y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_agc * R * T_des) /
                   (Paem - Phi_aem * Psat(T_des))
        y_O2_cem = (Pcem - Phi_cem * Psat(T_cgc) - C_N2_cgc * R * T_cgc) /
                   (Pcem - Phi_cem * Psat(T_cgc))

        # Molar masses at the anode side
        M_asm = if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            Phi_asm * Psat(T_des) / Pasm * M_H2O +
            (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
        else
            Phi_asm * Psat(T_des) / Pasm * M_H2O +
            y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 +
            (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
        end
        M_aem = if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            Phi_aem * Psat(T_des) / Paem * M_H2O +
            (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
        else
            Phi_aem * Psat(T_des) / Paem * M_H2O +
            y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 +
            (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
        end

        # Molar masses at the cathode side
        M_csm = Phi_csm * Psat(T_des) / Pcsm * M_H2O +
                y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 +
                (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
        M_cem = Phi_cem * Psat(T_des) / Pcem * M_H2O +
                y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 +
                (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2

        # Density of the gas mixture.
        rho_asm = Pasm / (R * T_des) * M_asm
        rho_aem = Paem / (R * T_des) * M_aem
        rho_csm = Pcsm / (R * T_des) * M_csm
        rho_cem = Pcem / (R * T_des) * M_cem

        # Boundary velocity at the anode recirculation inlet.
        if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            v_re = Ware !== nothing ? Ware / rho_aem_out_re / A_T_a : nothing
        else # cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
            v_re = nothing
        end

        # Vapor ratio over the gas mixture.
        x_H2O_v_asm = Phi_asm * Psat(T_des) / Pasm
        x_H2O_v_aem = Phi_aem * Psat(T_des) / Paem
        x_H2O_v_csm = Phi_csm * Psat(T_des) / Pcsm
        x_H2O_v_cem = Phi_cem * Psat(T_des) / Pcem

        # Dynamic viscosity of the gas mixture at the anode side.
        mu_asm = if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_asm, 1 - x_H2O_v_asm], T_des)
        else
            mu_mixture_gases(["H2O_v", "H2", "N2"],
                             [x_H2O_v_asm, y_H2_in * (1 - x_H2O_v_asm), (1 - y_H2_in) * (1 - x_H2O_v_asm)],
                             T_des)
        end
        mu_aem = if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_aem, 1 - x_H2O_v_aem], T_des)
        else
            mu_mixture_gases(["H2O_v", "H2", "N2"],
                             [x_H2O_v_aem, y_H2_aem * (1 - x_H2O_v_aem), (1 - y_H2_aem) * (1 - x_H2O_v_aem)],
                             T_des)
        end

        # Dynamic viscosity of the gas mixture at the cathode side.
        mu_csm = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                  [x_H2O_v_csm, y_O2_ext * (1 - x_H2O_v_csm), (1 - y_O2_ext) * (1 - x_H2O_v_csm)],
                                  T_des)
        mu_cem = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                  [x_H2O_v_cem, y_O2_cem * (1 - x_H2O_v_cem), (1 - y_O2_cem) * (1 - x_H2O_v_cem)],
                                  T_des)

        # Build typed manifold intermediates.
        manifold_intermediates = ManifoldIntermediates{1}(
            (ManifoldNodeIntermediates(Pasm, Phi_asm, M_asm, rho_asm, mu_asm),),
            (ManifoldNodeIntermediates(Paem, Phi_aem, M_aem, rho_aem, mu_aem),),
            (ManifoldNodeIntermediates(Pcsm, Phi_csm, M_csm, rho_csm, mu_csm),),
            (ManifoldNodeIntermediates(Pcem, Phi_cem, M_cem, rho_cem, mu_cem),)
        )

        Abp_a = clamp(sv_auxiliary.Abp_a, 0.0, A_T_a)
        Abp_c = clamp(sv_auxiliary.Abp_c, 0.0, A_T_c)
        auxiliary_intermediates = AuxiliaryIntermediates(v_re, k_purge, Abp_a, Abp_c)

        _ = manifold_intermediates
        _ = auxiliary_intermediates
    end

    agc_int = AnodeGCIntermediates(P_agc, C_tot_agc, rho_agc, mu_gaz_agc)
    cgc_int = CathodeGCIntermediates(P_cgc, C_tot_cgc, rho_cgc, mu_gaz_cgc)
    return MEAIntermediates1D{nb_gdl, nb_mpl}(agc_int, cgc_int, rho_Cp0, i_n, T_acl_mem_ccl, v_re, k_purge)
end


# Typed solver-vector helpers depend on model structs (`MEAState1D`, manifold and
# auxiliary containers). They are loaded only in the `Models` include context.
if @isdefined(MEAState1D)

"""Count the number of MEA/GC variables per gas-channel node in the solver vector."""
function _nb_solver_vars_per_gc(nb_gdl::Int, nb_mpl::Int)::Int
    return fieldcount(AnodeGCState) +
           nb_gdl * fieldcount(AnodeGDLState) +
           nb_mpl * fieldcount(AnodeMPLState) +
           fieldcount(AnodeCLState) +
           fieldcount(MembraneState) +
           fieldcount(CathodeCLState) +
           nb_mpl * fieldcount(CathodeMPLState) +
           nb_gdl * fieldcount(CathodeGDLState) +
           fieldcount(CathodeGCState)
end

"""Unpack one GC-node segment of the solver vector into a typed 1D MEA state."""
function _unpack_mea_state_1D(values::AbstractVector{<:Real}, nb_gdl::Int, nb_mpl::Int)
    idx = 1

    C_v_agc = Float64(values[idx]); idx += 1
    C_v_agdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    C_v_ampl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    C_v_acl = Float64(values[idx]); idx += 1
    C_v_ccl = Float64(values[idx]); idx += 1
    C_v_cmpl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    C_v_cgdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    C_v_cgc = Float64(values[idx]); idx += 1

    s_agc = Float64(values[idx]); idx += 1
    s_agdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    s_ampl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    s_acl = Float64(values[idx]); idx += 1
    s_ccl = Float64(values[idx]); idx += 1
    s_cmpl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    s_cgdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    s_cgc = Float64(values[idx]); idx += 1

    lambda_acl = Float64(values[idx]); idx += 1
    lambda_mem = Float64(values[idx]); idx += 1
    lambda_ccl = Float64(values[idx]); idx += 1

    C_H2_agc = Float64(values[idx]); idx += 1
    C_H2_agdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    C_H2_ampl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    C_H2_acl = Float64(values[idx]); idx += 1
    C_O2_ccl = Float64(values[idx]); idx += 1
    C_O2_cmpl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    C_O2_cgdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    C_O2_cgc = Float64(values[idx]); idx += 1
    C_N2_agc = Float64(values[idx]); idx += 1
    C_N2_cgc = Float64(values[idx]); idx += 1

    T_agc = Float64(values[idx]); idx += 1
    T_agdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    T_ampl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    T_acl = Float64(values[idx]); idx += 1
    T_mem = Float64(values[idx]); idx += 1
    T_ccl = Float64(values[idx]); idx += 1
    T_cmpl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    T_cgdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    T_cgc = Float64(values[idx]); idx += 1
    eta_c = Float64(values[idx]); idx += 1

    idx == length(values) + 1 ||
        throw(ArgumentError("Invalid 1D state segment length while unpacking solver vector."))

    agc = AnodeGCState(T_agc, C_v_agc, s_agc, C_H2_agc, C_N2_agc)
    agdl = ntuple(i -> AnodeGDLState(T_agdl[i], C_v_agdl[i], s_agdl[i], C_H2_agdl[i]), nb_gdl)
    ampl = ntuple(i -> AnodeMPLState(T_ampl[i], C_v_ampl[i], s_ampl[i], C_H2_ampl[i]), nb_mpl)
    acl = AnodeCLState(T_acl, C_v_acl, s_acl, lambda_acl, C_H2_acl)
    mem = MembraneState(T_mem, lambda_mem)
    ccl = CathodeCLState(T_ccl, C_v_ccl, s_ccl, lambda_ccl, C_O2_ccl, eta_c)
    cmpl = ntuple(i -> CathodeMPLState(T_cmpl[i], C_v_cmpl[i], s_cmpl[i], C_O2_cmpl[i]), nb_mpl)
    cgdl = ntuple(i -> CathodeGDLState(T_cgdl[i], C_v_cgdl[i], s_cgdl[i], C_O2_cgdl[i]), nb_gdl)
    cgc = CathodeGCState(T_cgc, C_v_cgc, s_cgc, C_O2_cgc, C_N2_cgc)

    return MEAState1D{nb_gdl, nb_mpl}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end

"""Create a derivative container initialized with NaN values."""
function _nan_mea_derivative_1D(nb_gdl::Int, nb_mpl::Int)
    z = NaN
    agc = AnodeGCDerivative(z, z, z, z, z)
    agdl = ntuple(_ -> AnodeGDLDerivative(z, z, z, z), nb_gdl)
    ampl = ntuple(_ -> AnodeMPLDerivative(z, z, z, z), nb_mpl)
    acl = AnodeCLDerivative(z, z, z, z, z)
    mem = MembraneDerivative(z, z)
    ccl = CathodeCLDerivative(z, z, z, z, z, z)
    cmpl = ntuple(_ -> CathodeMPLDerivative(z, z, z, z), nb_mpl)
    cgdl = ntuple(_ -> CathodeGDLDerivative(z, z, z, z), nb_gdl)
    cgc = CathodeGCDerivative(z, z, z, z, z)
    return MEACellDerivative1D{nb_gdl, nb_mpl}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end

"""Ensure all derivatives were assigned before returning to the solver."""
function _assert_derivative_complete(d::MEACellDerivative1D{nb_gdl, nb_mpl}) where {nb_gdl, nb_mpl}
    any(isnan, _pack_mea_derivative_1D(d)) &&
        throw(ArgumentError("At least one derivative entry is missing (NaN sentinel detected)."))
    return nothing
end

"""Repack one typed 1D derivative container into the solver ordering."""
function _pack_mea_derivative_1D(d::MEACellDerivative1D{nb_gdl, nb_mpl}) where {nb_gdl, nb_mpl}
    out = Float64[]
    sizehint!(out, _nb_solver_vars_per_gc(nb_gdl, nb_mpl))

    # C_v block
    push!(out, d.agc.C_v)
    append!(out, [d.agdl[i].C_v for i in 1:nb_gdl])
    append!(out, [d.ampl[i].C_v for i in 1:nb_mpl])
    push!(out, d.acl.C_v)
    push!(out, d.ccl.C_v)
    append!(out, [d.cmpl[i].C_v for i in 1:nb_mpl])
    append!(out, [d.cgdl[i].C_v for i in 1:nb_gdl])
    push!(out, d.cgc.C_v)

    # s block
    push!(out, d.agc.s)
    append!(out, [d.agdl[i].s for i in 1:nb_gdl])
    append!(out, [d.ampl[i].s for i in 1:nb_mpl])
    push!(out, d.acl.s)
    push!(out, d.ccl.s)
    append!(out, [d.cmpl[i].s for i in 1:nb_mpl])
    append!(out, [d.cgdl[i].s for i in 1:nb_gdl])
    push!(out, d.cgc.s)

    # lambda block
    push!(out, d.acl.lambda)
    push!(out, d.mem.lambda)
    push!(out, d.ccl.lambda)

    # H2/O2/N2 block
    push!(out, d.agc.C_H2)
    append!(out, [d.agdl[i].C_H2 for i in 1:nb_gdl])
    append!(out, [d.ampl[i].C_H2 for i in 1:nb_mpl])
    push!(out, d.acl.C_H2)
    push!(out, d.ccl.C_O2)
    append!(out, [d.cmpl[i].C_O2 for i in 1:nb_mpl])
    append!(out, [d.cgdl[i].C_O2 for i in 1:nb_gdl])
    push!(out, d.cgc.C_O2)
    push!(out, d.agc.C_N2)
    push!(out, d.cgc.C_N2)

    # T block + eta_c
    push!(out, d.agc.T)
    append!(out, [d.agdl[i].T for i in 1:nb_gdl])
    append!(out, [d.ampl[i].T for i in 1:nb_mpl])
    push!(out, d.acl.T)
    push!(out, d.mem.T)
    push!(out, d.ccl.T)
    append!(out, [d.cmpl[i].T for i in 1:nb_mpl])
    append!(out, [d.cgdl[i].T for i in 1:nb_gdl])
    push!(out, d.cgc.T)
    push!(out, d.ccl.eta_c)

    return out
end

"""Count manifold state variables in the solver vector."""
function _nb_solver_vars_manifolds(nb_man::Int)::Int
    manifold_lines = (:asm, :aem, :csm, :cem)
    return length(manifold_lines) * nb_man * fieldcount(ManifoldState)
end

"""Count auxiliary state variables currently present in the solver vector."""
function _nb_solver_vars_auxiliary(type_auxiliary::Symbol)::Int
    if type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                          :forced_convective_cathode_with_anodic_recirculation)
        return length(_packed_auxiliary_fields())
    end
    return 0
end

"""Fields currently packed in the solver vector for auxiliaries."""
function _packed_auxiliary_fields()
    return fieldnames(Auxiliary0DState)
end

"""Unpack manifold state from a solver-vector segment."""
function _unpack_manifold_state(values::AbstractVector{<:Real}, nb_manifold::Int)
    idx = 1

    P_asm = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold
    P_aem = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold
    P_csm = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold
    P_cem = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold

    Phi_asm = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold
    Phi_aem = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold
    Phi_csm = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold
    Phi_cem = ntuple(i -> Float64(values[idx + i - 1]), nb_manifold); idx += nb_manifold

    # Defensive check: ensure manifold segment length is exactly as expected.
    idx == length(values) + 1 ||
        throw(ArgumentError("Invalid manifold segment length while unpacking solver vector."))

    asm = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldState(P_asm[i], Phi_asm[i]), nb_manifold))
    aem = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldState(P_aem[i], Phi_aem[i]), nb_manifold))
    csm = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldState(P_csm[i], Phi_csm[i]), nb_manifold))
    cem = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldState(P_cem[i], Phi_cem[i]), nb_manifold))

    return _ManifoldStateBundle(asm, aem, csm, cem)
end

"""Create manifold derivative container initialised with NaN sentinels."""
function _nan_manifold_derivative_state(nb_manifold::Int)
    z = NaN
    mkline(n) = ManifoldLineDerivative{n}(ntuple(_ -> ManifoldDerivative(z, z), n))
    return _ManifoldDerivativeBundle(mkline(nb_manifold), mkline(nb_manifold),
                                     mkline(nb_manifold), mkline(nb_manifold))
end

"""Pack manifold derivatives into solver-vector ordering."""
function _pack_manifold_derivative_state(md)
    out = Float64[]

    append!(out, [md.asm.nodes[i].P for i in eachindex(md.asm.nodes)])
    append!(out, [md.aem.nodes[i].P for i in eachindex(md.aem.nodes)])
    append!(out, [md.csm.nodes[i].P for i in eachindex(md.csm.nodes)])
    append!(out, [md.cem.nodes[i].P for i in eachindex(md.cem.nodes)])

    append!(out, [md.asm.nodes[i].Phi for i in eachindex(md.asm.nodes)])
    append!(out, [md.aem.nodes[i].Phi for i in eachindex(md.aem.nodes)])
    append!(out, [md.csm.nodes[i].Phi for i in eachindex(md.csm.nodes)])
    append!(out, [md.cem.nodes[i].Phi for i in eachindex(md.cem.nodes)])

    return out
end

"""Ensure manifold derivatives are fully assigned."""
function _assert_manifold_derivative_complete(md)
    any(isnan, _pack_manifold_derivative_state(md)) &&
        throw(ArgumentError("At least one manifold derivative entry is missing (NaN sentinel detected)."))
    return nothing
end

"""Unpack auxiliary state from solver-vector segment."""
function _unpack_auxiliary_state(values::AbstractVector{<:Real})::Auxiliary0DState
    packed_fields = _packed_auxiliary_fields()
    # Keep strict consistency with the field-based packing schema.
    length(values) == length(packed_fields) ||
        throw(ArgumentError("Invalid auxiliary segment length while unpacking solver vector."))

    packed_values = NamedTuple{packed_fields}(Tuple(Float64(values[i]) for i in eachindex(values)))
    return Auxiliary0DState((getfield(packed_values, f) for f in packed_fields)...)
end

"""Create auxiliary derivative container initialised with NaN sentinels."""
function _nan_auxiliary_derivative()::Auxiliary0DDerivative
    z = NaN
    return Auxiliary0DDerivative(ntuple(_ -> z, fieldcount(Auxiliary0DDerivative))...)
end

"""Pack auxiliary derivatives into current solver-vector ordering."""
function _pack_auxiliary_derivative(d::Auxiliary0DDerivative)
    return Float64[getfield(d, f) for f in _packed_auxiliary_fields()]
end

"""Ensure auxiliary derivatives are fully assigned for active ordered fields."""
function _assert_auxiliary_derivative_complete(d::Auxiliary0DDerivative)
    any(isnan, _pack_auxiliary_derivative(d)) &&
        throw(ArgumentError("At least one auxiliary derivative entry is missing (NaN sentinel detected)."))
    return nothing
end

end


# ______________________________________Function which gives the integration event______________________________________

"""Create an integration event that stops the solver when a crucial variable (C_v, lambda, C_O2, C_H2)
becomes non-physical (below a threshold of 1e-5).

Parameters
----------
t :
    Time (s).
y : Vector
    Vector of the solver variables.
fc : AbstractFuelCell
    Fuel cell instance providing model parameters (used to retrieve `nb_gc` and
    the variable count per gas-channel node via `solver_variable_names`).
solver_variable_names : Vector{<:Vector{<:String}}
    Names of the solver variables.

Returns
-------
The difference between the minimum value of the crucial variables and 1e-5.
"""
function event_negative(t, y::Vector, fc::AbstractFuelCell,
                        solver_variable_names::Vector{<:Vector{<:String}})

    nb_gc = fc.numerical_parameters.nb_gc
    negative_solver_variables = Dict() # Dictionary to store the crucial variables.

    for (index, key) in enumerate(solver_variable_names[1])
        if startswith(key, "C_v_") || startswith(key, "lambda_") ||
           startswith(key, "C_O2_") || startswith(key, "C_H2_")
            for i in 1:nb_gc
                negative_solver_variables["$(key)_$i"] = y[index + (i - 1) * length(solver_variable_names[1])]
            end
        end
    end

    return minimum(values(negative_solver_variables)) - 1e-5 # 1e-5 is a control parameter to stop the program
                                                               # before having negative values.
end

