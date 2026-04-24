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
    return CellIntermediates1D{nb_gdl, nb_mpl}(agc_int, cgc_int, rho_Cp0, i_n, T_acl_mem_ccl, v_re, k_purge)
end


"""Return canonical Cell (MEA+GC) variable names for one GC node in solver ordering."""
function canonical_cell_solver_variable_names_1D(nb_gdl::Int, nb_mpl::Int)::Vector{String}
    return vcat(
        ["C_v_agc"], ["C_v_agdl_$(i)" for i in 1:nb_gdl], ["C_v_ampl_$(i)" for i in 1:nb_mpl],
        ["C_v_acl", "C_v_ccl"], ["C_v_cmpl_$(i)" for i in 1:nb_mpl], ["C_v_cgdl_$(i)" for i in 1:nb_gdl], ["C_v_cgc"],

        ["s_agc"], ["s_agdl_$(i)" for i in 1:nb_gdl], ["s_ampl_$(i)" for i in 1:nb_mpl],
        ["s_acl", "s_ccl"], ["s_cmpl_$(i)" for i in 1:nb_mpl], ["s_cgdl_$(i)" for i in 1:nb_gdl], ["s_cgc"],

        ["lambda_acl", "lambda_mem", "lambda_ccl"],

        ["C_H2_agc"], ["C_H2_agdl_$(i)" for i in 1:nb_gdl], ["C_H2_ampl_$(i)" for i in 1:nb_mpl], ["C_H2_acl", "C_O2_ccl"],
        ["C_O2_cmpl_$(i)" for i in 1:nb_mpl], ["C_O2_cgdl_$(i)" for i in 1:nb_gdl], ["C_O2_cgc", "C_N2_agc", "C_N2_cgc"],

        ["T_agc"], ["T_agdl_$(i)" for i in 1:nb_gdl], ["T_ampl_$(i)" for i in 1:nb_mpl],
        ["T_acl", "T_mem", "T_ccl"], ["T_cmpl_$(i)" for i in 1:nb_mpl], ["T_cgdl_$(i)" for i in 1:nb_gdl], ["T_cgc"],

        ["eta_c"],
    )
end

"""Backward-compatible alias for historical naming."""
canonical_mea_solver_variable_names(nb_gdl::Int, nb_mpl::Int) =
    canonical_cell_solver_variable_names_1D(nb_gdl, nb_mpl)

"""Return canonical algebraic variable names for the DAE.

The algebraic block is appended after the differential solver state and follows
the fixed ordering:
    [U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc], J_a_in, J_c_in]
"""
function canonical_algebraic_solver_variable_names(nb_gc::Int)::Vector{String}
    return vcat(
        ["U_cell"],
        ["i_fc_$(i)" for i in 1:nb_gc],
        ["C_O2_Pt_$(i)" for i in 1:nb_gc],
        ["J_a_in", "J_c_in"],
    )
end

"""Return the reference value for a canonical DAE algebraic variable name."""
function algebraic_state_scale(name::AbstractString, scaling::DAEAlgebraicScaling)::Float64
    name == "U_cell" && return scaling.U
    startswith(name, "i_fc_") && return scaling.i_fc
    startswith(name, "C_O2_Pt_") && return scaling.C_O2_Pt
    (name == "J_a_in" || name == "J_c_in") && return scaling.J_in
    throw(ArgumentError("Unknown canonical algebraic state variable name: \"$name\""))
end

"""Build scaling references for the DAE algebraic block."""
function build_algebraic_state_scaling(nb_gc::Int,
                                       scaling::DAEAlgebraicScaling)::Vector{Float64}
    names = canonical_algebraic_solver_variable_names(nb_gc)
    return [algebraic_state_scale(n, scaling) for n in names]
end

# ─────────────────────────────────────────────────────────────────────────────
# State vector scaling — build, scale, unscale
#
# Why this exists
# ----------------
# The solver state mixes variables with heterogeneous physical magnitudes:
#   - species concentrations are typically O(10)
#   - temperatures are typically O(10²)
#   - cathode overpotential is typically O(10⁻¹)
#   - liquid saturation can be much smaller than 1 in dry/transient regimes
#
# Such disparities do not change the physics, but they can make the numerical
# integration problem harder to condition and can make uniform tolerances less
# balanced across components. The solver therefore works on a dimensionless
# state vector whose entries are brought closer to O(1).
#
# Convention:  x̂ = x / x_ref   (dimensionless, order-1 solver variable)
#              dx̂/dt = (dx/dt) / x_ref
#
# Important: the physical equations are never rewritten in scaled variables.
# Scaling is applied only at the solver boundary, so all model kernels continue
# to read and produce physical quantities in their native units:
#   - y0_phys    → y0_scaled    before ODEProblem
#   - y_scaled   → y_phys       at the start of dydt!
#   - dy_phys    → dy_scaled    at the end of dydt!
#   - sol_scaled → sol_phys     in recovery!
# ─────────────────────────────────────────────────────────────────────────────

"""
    cell_state_scale(name, scaling) -> Float64

Return the reference value for a canonical cell state variable name.
The mapping is derived from the variable prefix, consistent with
`canonical_cell_solver_variable_names_1D`.

This keeps the canonical solver ordering as the single source of truth while
ensuring that each physical family (temperature, concentrations, saturation,
etc.) receives an appropriate reference magnitude.
"""
function cell_state_scale(name::AbstractString, scaling::CellStateScaling)::Float64
    startswith(name, "C_v_")     && return scaling.C_v
    startswith(name, "s_")       && return scaling.s
    startswith(name, "lambda_")  && return scaling.lambda
    startswith(name, "C_H2_")    && return scaling.C_H2
    startswith(name, "C_O2_")    && return scaling.C_O2
    startswith(name, "C_N2_")    && return scaling.C_N2
    startswith(name, "T_")       && return scaling.T
    name == "eta_c"              && return scaling.eta_c
    throw(ArgumentError("Unknown canonical cell state variable name: \"$name\""))
end

"""
    build_cell_state_scaling_1D(nb_gdl, nb_mpl, scaling) -> Vector{Float64}

Build the reference vector for one GC-node cell state segment, in canonical
solver ordering, derived from `canonical_cell_solver_variable_names_1D`.
No manual re-enumeration of the variable order.

This design avoids maintaining an additional hard-coded layout for the scaling,
which would be error-prone and could silently diverge from the true solver
ordering.
"""
function build_cell_state_scaling_1D(nb_gdl::Int, nb_mpl::Int,
                                     scaling::CellStateScaling)::Vector{Float64}
    names = canonical_cell_solver_variable_names_1D(nb_gdl, nb_mpl)
    return [cell_state_scale(n, scaling) for n in names]
end

"""
    build_solver_state_scaling(nb_gc, nb_gdl, nb_mpl, nb_man, type_auxiliary,
                               state_scaling) -> Vector{Float64}

Build the full solver state scaling vector, in the same three-segment ordering
as the ODE solver vector:
  1. Cell P2D  (nb_gc repetitions of the 1D cell block)
  2. Manifolds (identity placeholders — see warning below)
  3. Auxiliaries (identity placeholders — see warning below)

The purpose is to provide a single scaling vector aligned with the exact solver
layout so that scaling/unscaling can be done mechanically, without touching the
physical equations. For manifolds and auxiliaries, the current references are
temporary identity placeholders; warnings are emitted to make this explicit.
"""
function build_solver_state_scaling(nb_gc::Int, nb_gdl::Int, nb_mpl::Int, nb_man::Int,
                                    type_auxiliary::Symbol,
                                    state_scaling::StateScaling;
                                    include_algebraic::Bool=false)::Vector{Float64}
    scales = Float64[]

    # Segment 1 — Cell P2D: repeat the 1D block for each gas-channel node.
    cell_scales_1D = build_cell_state_scaling_1D(nb_gdl, nb_mpl, state_scaling.cell)
    append!(scales, repeat(cell_scales_1D, nb_gc))

    # Segment 2 — Manifolds (active only with auxiliary subsystems).
    has_auxiliary = type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                                       :forced_convective_cathode_with_anodic_recirculation)
    if has_auxiliary
        # TODO: dedicated manifold scaling references should be calibrated and
        #       set here once the manifold subsystem is fully reactivated.
        #       Current placeholder uses identity (P_ref=1, Phi_ref=1), which
        #       leaves manifold variables in physical units. This is acceptable
        #       as a temporary compatibility choice, but it does not provide the
        #       same order-1 normalization as the cell state vector.
        manifold_scaling_warning =
            "Manifold state scaling uses identity references (P_ref=1, Phi_ref=1). " *
            "This is a temporary placeholder. Dedicated manifold scaling should be " *
            "adjusted when the manifold subsystem is reactivated."
        @warn manifold_scaling_warning maxlog=1
        m = state_scaling.manifold
        # Order mirrors _pack_manifold_derivative_state! / _unpack_manifold_state:
        # P_asm, P_aem, P_csm, P_cem, Phi_asm, Phi_aem, Phi_csm, Phi_cem
        append!(scales, fill(m.P,   4 * nb_man))
        append!(scales, fill(m.Phi, 4 * nb_man))

        # TODO: dedicated auxiliary scaling references should be calibrated and
        #       set here once the BoP auxiliary subsystem is fully reactivated.
        #       Current placeholder uses identity, leaving auxiliary variables
        #       (Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c) in physical units. This keeps
        #       the implementation simple today, but the same conditioning logic
        #       should eventually be extended to these states as well.
        auxiliary_scaling_warning =
            "Auxiliary state scaling uses identity references (all refs = 1). " *
            "This is a temporary placeholder. Dedicated auxiliary scaling should be " *
            "adjusted when the auxiliary subsystem is reactivated."
        @warn auxiliary_scaling_warning maxlog=1
        a = state_scaling.auxiliary
        # Order mirrors _packed_auxiliary_fields() / Auxiliary0DState field order.
        append!(scales, [a.Wcp, a.Wa_inj, a.Wc_inj, a.Abp_a, a.Abp_c])
    end

    # Segment 4 — DAE algebraic variables appended after the differential state.
    include_algebraic && append!(scales, build_algebraic_state_scaling(nb_gc, state_scaling.dae_algebraic))

    return scales
end

"""Build solver scaling from any simulation-like container exposing `fuel_cell` and `cfg`.

This helper intentionally avoids a concrete `AlphaPEM` type annotation to prevent
cross-module load-order coupling (`core/modules` is loaded before `core/models`).
"""
function build_solver_state_scaling(simu;
                                    include_algebraic::Bool=true)::Vector{Float64}
    return build_internal_solver_state_scaling(simu.fuel_cell, simu.cfg;
                                               include_algebraic=include_algebraic)
end

"""
    build_internal_solver_state_scaling(fc, cfg) -> Vector{Float64}

Build the fixed internal solver scaling vector used during integration.

This helper lives in `core/modules` rather than `core/models` so that the
solver-scaling infrastructure remains grouped with the packing/unpacking helpers
and can be reused by both the RHS and the simulation orchestration code.
"""
function build_internal_solver_state_scaling(fc::AbstractFuelCell,
                                             cfg::SimulationConfig;
                                             include_algebraic::Bool=false)::Vector{Float64}
    np = fc.numerical_parameters
    return build_solver_state_scaling(np.nb_gc, np.nb_gdl, np.nb_mpl, np.nb_man,
                                      cfg.type_auxiliary, StateScaling();
                                      include_algebraic=include_algebraic)
end

"""
    scale_values(values, scales) -> Vector{Float64}

Convert a physical state vector to a dimensionless (scaled) solver vector.
    x̂[i] = x[i] / scales[i]

This transformation is used only for the solver-facing representation, not for
the model equations themselves.
"""
@inline function scale_values(values::AbstractVector{Float64},
                              scales::AbstractVector{Float64})::Vector{Float64}
    return values ./ scales
end

"""
    unscale_values(values_scaled, scales) -> Vector{Float64}

Convert a dimensionless (scaled) solver vector back to physical units.
    x[i] = x̂[i] * scales[i]

This inverse transformation is what allows the RHS and post-processing code to
continue working entirely with physical variables.
"""
@inline function unscale_values(values_scaled::AbstractVector{Float64},
                                scales::AbstractVector{Float64})::Vector{Float64}
    return values_scaled .* scales
end

"""Count the number of CellState1D variables (MEA+GC) per gas-channel node in the solver vector."""
function _nb_solver_vars_cell_1D(nb_gdl::Int, nb_mpl::Int)::Int
    return length(canonical_cell_solver_variable_names_1D(nb_gdl, nb_mpl))
end

"""Count algebraic variables appended in the DAE solver vector."""
function _nb_solver_vars_algebraic(nb_gc::Int)::Int
    return length(canonical_algebraic_solver_variable_names(nb_gc))
end

"""Build the solver differential/algebraic mask for DAEProblem construction.

Why this mask matters
---------------------
IDA solves a mixed DAE system `F(t, y, ydot) = 0` where only part of `y`
corresponds to true differential states. The remaining states are algebraic
constraints and must be flagged explicitly.

This `BitVector` is used in two places with the exact same meaning:
1) `DAEProblem(...; differential_vars=mask)` so Sundials knows which rows are
   differential vs algebraic;
2) `_build_consistent_initial_solver_derivatives(...)` so only differential
   entries are assigned in `ydot0` while algebraic entries are left untouched.

Keeping one canonical constructor for this mask avoids subtle inconsistencies
between problem definition and initialisation.
"""
function build_solver_differential_vars(nb_gc::Int, nb_gdl::Int, nb_mpl::Int, nb_man::Int,
                                        type_auxiliary::Symbol;
                                        include_algebraic::Bool=false)::BitVector
    n_vars_cell_1D = _nb_solver_vars_cell_1D(nb_gdl, nb_mpl)
    has_auxiliary = type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                                       :forced_convective_cathode_with_anodic_recirculation)
    n_vars_manifold = has_auxiliary ? _nb_solver_vars_manifolds(nb_man) : 0
    n_vars_auxiliary = _nb_solver_vars_auxiliary(type_auxiliary)
    n_diff = nb_gc * n_vars_cell_1D + n_vars_manifold + n_vars_auxiliary

    if include_algebraic
        n_alg = _nb_solver_vars_algebraic(nb_gc)
        # Differential block first, algebraic block appended at the end.
        # This mirrors the canonical DAE state layout used everywhere else.
        return vcat(trues(n_diff), falses(n_alg))
    end
    return trues(n_diff)
end

"""Build a consistent initial `dy/dt` guess in scaled solver coordinates for IDA.

The DAE residual uses the convention:
    F_diff = ydot_IDA - f(y)

At `t0`, if we started with `ydot0 = 0`, the differential residual would be
`F_diff(t0) = -f(y0)`, which can be large and lead to poor Newton starts.
This helper computes one residual evaluation at `(t0, y0, ydot=0)` and sets:
    ydot0[i] = -F_i(t0, y0, 0)
for differential rows only, i.e. `ydot0 ≈ f(y0)`.

Algebraic rows are intentionally left at zero because those equations do not
define a time derivative.
"""
function _build_consistent_initial_solver_derivatives(residual!,
                                                      packed,
                                                      initial_solver_values::Vector{Float64},
                                                      t0::Float64,
                                                      differential_vars::BitVector)::Vector{Float64}
    n = length(initial_solver_values)
    n == length(differential_vars) ||
        throw(ArgumentError("differential_vars size mismatch in _build_consistent_initial_solver_derivatives."))

    # Start from a neutral guess: no derivative anywhere.
    dydt0 = zeros(Float64, n)

    # One residual call at (t0, y0, ydot=0) gives the mismatch F(t0, y0, 0).
    # For differential rows, we then negate this mismatch to enforce
    # F_diff(t0, y0, ydot0) ≈ 0 before IDA starts iterating.
    res0 = zeros(Float64, n)
    residual!(res0, dydt0, initial_solver_values, packed, t0)

    @inbounds for i in eachindex(differential_vars)
        differential_vars[i] || continue
        dydt0[i] = -res0[i]
    end

    return dydt0
end

# Typed solver-vector helpers depend on model structs (`CellState1D`, manifold and
# auxiliary containers). They are loaded only in the `Models` include context.
if @isdefined(CellState1D)

"""Unpack one GC-node segment of the solver vector into a typed `CellState1D`."""
@inline function _unpack_cell_state_1D(values::AbstractVector{<:Real}, nb_gdl::Int, nb_mpl::Int)
    return _unpack_cell_state_1D(values, Val(nb_gdl), Val(nb_mpl))
end

@inline function _unpack_cell_state_1D(values::AbstractVector{<:Real}, ::Val{N_GDL}, ::Val{N_MPL}) where {N_GDL, N_MPL}
    idx = Ref(1)

    @inline read_scalar!() = begin
        @inbounds v = Float64(values[idx[]])
        idx[] += 1
        return v
    end
    @inline read_block!(::Val{N}) where {N} = ntuple(_ -> read_scalar!(), Val(N))

    C_v_agc = read_scalar!()
    C_v_agdl = read_block!(Val(N_GDL))
    C_v_ampl = read_block!(Val(N_MPL))
    C_v_acl = read_scalar!()
    C_v_ccl = read_scalar!()
    C_v_cmpl = read_block!(Val(N_MPL))
    C_v_cgdl = read_block!(Val(N_GDL))
    C_v_cgc = read_scalar!()

    s_agc = read_scalar!()
    s_agdl = read_block!(Val(N_GDL))
    s_ampl = read_block!(Val(N_MPL))
    s_acl = read_scalar!()
    s_ccl = read_scalar!()
    s_cmpl = read_block!(Val(N_MPL))
    s_cgdl = read_block!(Val(N_GDL))
    s_cgc = read_scalar!()

    lambda_acl = read_scalar!()
    lambda_mem = read_scalar!()
    lambda_ccl = read_scalar!()

    C_H2_agc = read_scalar!()
    C_H2_agdl = read_block!(Val(N_GDL))
    C_H2_ampl = read_block!(Val(N_MPL))
    C_H2_acl = read_scalar!()
    C_O2_ccl = read_scalar!()
    C_O2_cmpl = read_block!(Val(N_MPL))
    C_O2_cgdl = read_block!(Val(N_GDL))
    C_O2_cgc = read_scalar!()
    C_N2_agc = read_scalar!()
    C_N2_cgc = read_scalar!()

    T_agc = read_scalar!()
    T_agdl = read_block!(Val(N_GDL))
    T_ampl = read_block!(Val(N_MPL))
    T_acl = read_scalar!()
    T_mem = read_scalar!()
    T_ccl = read_scalar!()
    T_cmpl = read_block!(Val(N_MPL))
    T_cgdl = read_block!(Val(N_GDL))
    T_cgc = read_scalar!()
    eta_c = read_scalar!()

    idx[] == length(values) + 1 ||
        throw(ArgumentError("Invalid 1D state segment length while unpacking solver vector."))

    agc = AnodeGCState(T_agc, C_v_agc, s_agc, C_H2_agc, C_N2_agc)
    agdl = ntuple(i -> AnodeGDLState(T_agdl[i], C_v_agdl[i], s_agdl[i], C_H2_agdl[i]), Val(N_GDL))
    ampl = ntuple(i -> AnodeMPLState(T_ampl[i], C_v_ampl[i], s_ampl[i], C_H2_ampl[i]), Val(N_MPL))
    acl = AnodeCLState(T_acl, C_v_acl, s_acl, lambda_acl, C_H2_acl)
    mem = MembraneState(T_mem, lambda_mem)
    ccl = CathodeCLState(T_ccl, C_v_ccl, s_ccl, lambda_ccl, C_O2_ccl, eta_c)
    cmpl = ntuple(i -> CathodeMPLState(T_cmpl[i], C_v_cmpl[i], s_cmpl[i], C_O2_cmpl[i]), Val(N_MPL))
    cgdl = ntuple(i -> CathodeGDLState(T_cgdl[i], C_v_cgdl[i], s_cgdl[i], C_O2_cgdl[i]), Val(N_GDL))
    cgc = CathodeGCState(T_cgc, C_v_cgc, s_cgc, C_O2_cgc, C_N2_cgc)

    return CellState1D{N_GDL, N_MPL}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end

"""Create a derivative container initialized with NaN values."""
function _nan_cell_derivative_1D(nb_gdl::Int, nb_mpl::Int)
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
    return CellDerivative1D{nb_gdl, nb_mpl}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end

"""Ensure all scalar values in one typed 1D derivative are assigned (no NaN sentinel left)."""
function _assert_cell_derivative_complete!(d::CellDerivative1D{nb_gdl, nb_mpl}) where {nb_gdl, nb_mpl}
    fail() = throw(ArgumentError("At least one derivative entry is missing (NaN sentinel detected)."))

    isnan(d.agc.C_v) && fail()
    isnan(d.agc.s) && fail()
    isnan(d.agc.C_H2) && fail()
    isnan(d.agc.C_N2) && fail()
    isnan(d.agc.T) && fail()

    for i in 1:nb_gdl
        isnan(d.agdl[i].C_v) && fail()
        isnan(d.agdl[i].s) && fail()
        isnan(d.agdl[i].C_H2) && fail()
        isnan(d.agdl[i].T) && fail()
        isnan(d.cgdl[i].C_v) && fail()
        isnan(d.cgdl[i].s) && fail()
        isnan(d.cgdl[i].C_O2) && fail()
        isnan(d.cgdl[i].T) && fail()
    end

    for i in 1:nb_mpl
        isnan(d.ampl[i].C_v) && fail()
        isnan(d.ampl[i].s) && fail()
        isnan(d.ampl[i].C_H2) && fail()
        isnan(d.ampl[i].T) && fail()
        isnan(d.cmpl[i].C_v) && fail()
        isnan(d.cmpl[i].s) && fail()
        isnan(d.cmpl[i].C_O2) && fail()
        isnan(d.cmpl[i].T) && fail()
    end

    isnan(d.acl.C_v) && fail()
    isnan(d.acl.s) && fail()
    isnan(d.acl.lambda) && fail()
    isnan(d.acl.C_H2) && fail()
    isnan(d.acl.T) && fail()

    isnan(d.mem.lambda) && fail()
    isnan(d.mem.T) && fail()

    isnan(d.ccl.C_v) && fail()
    isnan(d.ccl.s) && fail()
    isnan(d.ccl.lambda) && fail()
    isnan(d.ccl.C_O2) && fail()
    isnan(d.ccl.T) && fail()
    isnan(d.ccl.eta_c) && fail()

    isnan(d.cgc.C_v) && fail()
    isnan(d.cgc.s) && fail()
    isnan(d.cgc.C_O2) && fail()
    isnan(d.cgc.C_N2) && fail()
    isnan(d.cgc.T) && fail()
    return nothing
end

"""Ensure all MEA node derivatives in the P2D container are fully assigned."""
function _assert_fuelcell_derivative_complete!(d::FuelCellDerivativeP2D{nb_gdl, nb_mpl, nb_gc}) where {nb_gdl, nb_mpl, nb_gc}
    for i in 1:nb_gc
        _assert_cell_derivative_complete!(d.nodes[i])
    end
    return nothing
end

"""Write one typed 1D derivative container directly into `dy` starting at `offset` (1-based, no allocation)."""
function _pack_cell_derivative_1D!(dy::AbstractVector{Float64}, offset::Int,
                                   d::CellDerivative1D{nb_gdl, nb_mpl}) where {nb_gdl, nb_mpl}
    idx = offset
    # C_v block
    dy[idx] = d.agc.C_v; idx += 1
    for i in 1:nb_gdl; dy[idx] = d.agdl[i].C_v; idx += 1; end
    for i in 1:nb_mpl; dy[idx] = d.ampl[i].C_v; idx += 1; end
    dy[idx] = d.acl.C_v; idx += 1
    dy[idx] = d.ccl.C_v; idx += 1
    for i in 1:nb_mpl; dy[idx] = d.cmpl[i].C_v; idx += 1; end
    for i in 1:nb_gdl; dy[idx] = d.cgdl[i].C_v; idx += 1; end
    dy[idx] = d.cgc.C_v; idx += 1
    # s block
    dy[idx] = d.agc.s; idx += 1
    for i in 1:nb_gdl; dy[idx] = d.agdl[i].s; idx += 1; end
    for i in 1:nb_mpl; dy[idx] = d.ampl[i].s; idx += 1; end
    dy[idx] = d.acl.s; idx += 1
    dy[idx] = d.ccl.s; idx += 1
    for i in 1:nb_mpl; dy[idx] = d.cmpl[i].s; idx += 1; end
    for i in 1:nb_gdl; dy[idx] = d.cgdl[i].s; idx += 1; end
    dy[idx] = d.cgc.s; idx += 1
    # lambda block
    dy[idx] = d.acl.lambda; idx += 1
    dy[idx] = d.mem.lambda; idx += 1
    dy[idx] = d.ccl.lambda; idx += 1
    # H2/O2/N2 block
    dy[idx] = d.agc.C_H2; idx += 1
    for i in 1:nb_gdl; dy[idx] = d.agdl[i].C_H2; idx += 1; end
    for i in 1:nb_mpl; dy[idx] = d.ampl[i].C_H2; idx += 1; end
    dy[idx] = d.acl.C_H2; idx += 1
    dy[idx] = d.ccl.C_O2; idx += 1
    for i in 1:nb_mpl; dy[idx] = d.cmpl[i].C_O2; idx += 1; end
    for i in 1:nb_gdl; dy[idx] = d.cgdl[i].C_O2; idx += 1; end
    dy[idx] = d.cgc.C_O2; idx += 1
    dy[idx] = d.agc.C_N2; idx += 1
    dy[idx] = d.cgc.C_N2; idx += 1
    # T block + eta_c
    dy[idx] = d.agc.T; idx += 1
    for i in 1:nb_gdl; dy[idx] = d.agdl[i].T; idx += 1; end
    for i in 1:nb_mpl; dy[idx] = d.ampl[i].T; idx += 1; end
    dy[idx] = d.acl.T; idx += 1
    dy[idx] = d.mem.T; idx += 1
    dy[idx] = d.ccl.T; idx += 1
    for i in 1:nb_mpl; dy[idx] = d.cmpl[i].T; idx += 1; end
    for i in 1:nb_gdl; dy[idx] = d.cgdl[i].T; idx += 1; end
    dy[idx] = d.cgc.T; idx += 1
    dy[idx] = d.ccl.eta_c
    return nothing
end

"""Write all GC-node typed derivatives directly into `dy` with zero allocations (SciML iip convention)."""
function _pack_fuelcell_derivative_p2d!(dydt::AbstractVector{Float64},
                                         d::FuelCellDerivativeP2D{nb_gdl, nb_mpl, nb_gc},
                                         n_vars_cell_1D::Int) where {nb_gdl, nb_mpl, nb_gc}
    for i in 1:nb_gc
        _pack_cell_derivative_1D!(dydt, (i - 1) * n_vars_cell_1D + 1, d.nodes[i])
    end
    return nothing
end

"""Write manifold derivatives directly into `dy` starting at `offset` (no allocation)."""
function _pack_manifold_derivative_state!(dydt::AbstractVector{Float64}, offset::Int, md)
    idx = offset
    for i in eachindex(md.asm.nodes); dydt[idx] = md.asm.nodes[i].P;   idx += 1; end
    for i in eachindex(md.aem.nodes); dydt[idx] = md.aem.nodes[i].P;   idx += 1; end
    for i in eachindex(md.csm.nodes); dydt[idx] = md.csm.nodes[i].P;   idx += 1; end
    for i in eachindex(md.cem.nodes); dydt[idx] = md.cem.nodes[i].P;   idx += 1; end
    for i in eachindex(md.asm.nodes); dydt[idx] = md.asm.nodes[i].Phi; idx += 1; end
    for i in eachindex(md.aem.nodes); dydt[idx] = md.aem.nodes[i].Phi; idx += 1; end
    for i in eachindex(md.csm.nodes); dydt[idx] = md.csm.nodes[i].Phi; idx += 1; end
    for i in eachindex(md.cem.nodes); dydt[idx] = md.cem.nodes[i].Phi; idx += 1; end
    return nothing
end

"""Write auxiliary derivatives directly into `dy` starting at `offset` (no allocation)."""
function _pack_auxiliary_derivative!(dydt::AbstractVector{Float64}, offset::Int, d::Auxiliary0DDerivative)
    idx = offset
    for f in _packed_auxiliary_fields()
        dydt[idx] = getfield(d, f); idx += 1
    end
    return nothing
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