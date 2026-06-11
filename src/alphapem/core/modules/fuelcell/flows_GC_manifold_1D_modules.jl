# -*- coding: utf-8 -*-

"""
This module is used to calculate intermediate values for auxiliary flow calculations.
"""

# _____________________________________________________Preliminaries____________________________________________________

"""
    anode_gc_order(nb_gc, type_flow)

Return the anode GC traversal order (from physical inlet to physical outlet).

- The returned sequence gives the **physical GC index** at each **numerical flow position**.
  In other words, `order[p] = i` means:
  - `p` is the position in flow-oriented numerical arrays (inlet -> outlet),
  - `i` is the corresponding physical GC index in the model state arrays (`sv[i]`).

- `:co_flow`      -> `1:nb_gc`
- `:counter_flow` -> `nb_gc:-1:1`
"""
@inline function anode_gc_order(nb_gc::Int, type_flow::Symbol)
    if type_flow == :co_flow
        return 1:nb_gc
    elseif type_flow == :counter_flow
        return nb_gc:-1:1
    else
        throw(ArgumentError("Unsupported type_flow for GC orientation: $(type_flow)"))
    end
end

"""
    anode_gc_pos_map(nb_gc, type_flow)

Return a lookup vector `pos` such that `pos[i]` gives the position of
physical GC index `i` along the anode flow direction (inlet -> outlet).

- `i` is the **physical GC index** used by model state arrays (`sv[i]`).
- `pos[i]` is the **numerical position** in flow-oriented arrays (inlet -> outlet).

This is the inverse mapping of `anode_gc_order`:
- `order[p] = i`  <=>  `pos[i] = p`
"""
@inline function anode_gc_pos_map(nb_gc::Int, type_flow::Symbol)::Vector{Int}
    order = anode_gc_order(nb_gc, type_flow)
    pos = zeros(Int, nb_gc)
    @inbounds for p in 1:nb_gc
        pos[order[p]] = p
    end
    return pos
end

"""
    flow_1D_GC_manifold_int_values(sv_1D_cell, sv_auxiliary, fc, cfg)

Calculate intermediate values used for GC/manifold flow computations.

Arguments
---------
sv_1D_cell : AbstractVector{<:CellState1D}
    Typed variables computed by the solver (internal stack states).
sv_auxiliary
    Typed variables computed by the auxiliary system (internal auxiliary states).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cfg : SimulationConfig
    Simulation configuration (provides `type_auxiliary`).

Returns
-------
Tuple
    `(P_agc, P_cgc, Phi_agc, Phi_cgc, y_H2_agc, y_O2_cgc, M_agc, M_cgc,
      M_ext, M_H2_N2_in, rho_agc, rho_cgc, k_purge, Abp_a, Abp_c,
      mu_gaz_agc, mu_gaz_cgc)`
"""
function flow_1D_GC_manifold_int_values(sv_1D_cell::AbstractVector{<:CellState1D},
                                        sv_auxiliary,
                                        fc::AbstractFuelCell,
                                        cfg::SimulationConfig)
    # Extract parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    T_des, y_H2_in = oc.T_des, oc.y_H2_in
    Hmem, Hacl, Hccl, kappa_co = pp.Hmem, pp.Hacl, pp.Hccl, pp.kappa_co
    nb_gc, purge_time, delta_purge = np.nb_gc, np.purge_time, np.delta_purge
    nb_gc = np.nb_gc

    # Initialise variables
    k_purge = nothing
    Abp_a = sv_auxiliary === nothing ? nothing : getproperty(sv_auxiliary, :Abp_a)
    Abp_c = sv_auxiliary === nothing ? nothing : getproperty(sv_auxiliary, :Abp_c)
    P_agc = Vector{Float64}(undef, nb_gc)
    P_cgc = Vector{Float64}(undef, nb_gc)
    Phi_agc = Vector{Float64}(undef, nb_gc)
    Phi_cgc = Vector{Float64}(undef, nb_gc)
    y_H2_agc = Vector{Float64}(undef, nb_gc)
    y_O2_cgc = Vector{Float64}(undef, nb_gc)
    M_agc = Vector{Float64}(undef, nb_gc)
    M_cgc = Vector{Float64}(undef, nb_gc)
    rho_agc = Vector{Float64}(undef, nb_gc)
    rho_cgc = Vector{Float64}(undef, nb_gc)
    mu_gaz_agc = Vector{Float64}(undef, nb_gc)
    mu_gaz_cgc = Vector{Float64}(undef, nb_gc)

    # Physical quantities outside the stack
    #       Molar masses
    Psat_Text = Psat(Text)
    Psat_Tdes = Psat(T_des)
    M_ext = Phi_ext * Psat_Text / Pext * M_H2O +
            y_O2_ext * (1 - Phi_ext * Psat_Text / Pext) * M_O2 +
            (1 - y_O2_ext) * (1 - Phi_ext * Psat_Text / Pext) * M_N2
    M_H2_N2_in = y_H2_in * M_H2 + (1 - y_H2_in) * M_N2

    # Physical quantities inside the stack
    #       Pressures, humidity ratios, dry-gas composition ratios,
    #       molar masses, densities, and dynamic viscosities
    for i in 1:nb_gc
        sv_i = sv_1D_cell[i]
        C_v_agc_i = sv_i.agc.C_v
        C_v_cgc_i = sv_i.cgc.C_v
        C_H2_agc_i = sv_i.agc.C_H2
        C_O2_cgc_i = sv_i.cgc.C_O2
        C_N2_agc_i = sv_i.agc.C_N2
        C_N2_cgc_i = sv_i.cgc.C_N2
        T_agc_i = sv_i.agc.T
        T_cgc_i = sv_i.cgc.T

        P_agc_i = (C_v_agc_i + C_H2_agc_i + C_N2_agc_i) * R * T_agc_i
        P_cgc_i = (C_v_cgc_i + C_O2_cgc_i + C_N2_cgc_i) * R * T_cgc_i
        P_agc[i] = P_agc_i
        P_cgc[i] = P_cgc_i

        Phi_agc_i = C_v_agc_i / C_v_sat(T_agc_i)
        Phi_cgc_i = C_v_cgc_i / C_v_sat(T_cgc_i)
        Phi_agc[i] = Phi_agc_i
        Phi_cgc[i] = Phi_cgc_i

        y_H2_agc_i = C_H2_agc_i / (C_H2_agc_i + C_N2_agc_i)
        y_O2_cgc_i = C_O2_cgc_i / (C_O2_cgc_i + C_N2_cgc_i)
        y_H2_agc[i] = y_H2_agc_i
        y_O2_cgc[i] = y_O2_cgc_i

        M_agc_i = C_v_agc_i * R * T_des / P_agc_i * M_H2O +
                  C_H2_agc_i * R * T_des / P_agc_i * M_H2 +
                  C_N2_agc_i * R * T_des / P_agc_i * M_N2
        vapor_share_cgc_i = Phi_cgc_i * Psat_Tdes / P_cgc_i
        M_cgc_i = vapor_share_cgc_i * M_H2O +
                  y_O2_cgc_i * (1 - vapor_share_cgc_i) * M_O2 +
                  (1 - y_O2_cgc_i) * (1 - vapor_share_cgc_i) * M_N2
        M_agc[i] = M_agc_i
        M_cgc[i] = M_cgc_i

        rho_agc[i] = P_agc_i / (R * T_agc_i) * M_agc_i
        rho_cgc[i] = P_cgc_i / (R * T_cgc_i) * M_cgc_i

        x_H2O_v_agc_i = C_v_agc_i / (C_v_agc_i + C_H2_agc_i + C_N2_agc_i)
        x_H2O_v_cgc_i = C_v_cgc_i / (C_v_cgc_i + C_O2_cgc_i + C_N2_cgc_i)
        mu_gaz_agc[i] = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_agc_i, 1 - x_H2O_v_agc_i], T_agc_i)
        mu_gaz_cgc[i] = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                         [x_H2O_v_cgc_i, y_O2_cgc_i * (1 - x_H2O_v_cgc_i),
                                          (1 - y_O2_cgc_i) * (1 - x_H2O_v_cgc_i)], T_cgc_i)
    end

    # Physical quantities in the auxiliary system
    if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # The typed auxiliary/manifold refactor is not finished yet in this helper.
        # Keep the branch syntactically valid for compilation, preserve the legacy
        # lines below for later recovery, and only retain safe placeholder values.
        k_purge = nothing
        Abp_a = Abp_a === nothing ? nothing : clamp(Abp_a, 0.0, pp.A_T_a)
        Abp_c = Abp_c === nothing ? nothing : clamp(Abp_c, 0.0, pp.A_T_c)

#= Legacy auxiliary block kept for later refactor.
         # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at EM
         y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_a * R * T_des) / (Paem - Phi_aem * Psat(T_des))
         y_O2_cem = (Pcem - Phi_cem * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem - Phi_cem * Psat(T_cgc))

         # Molar masses
         if parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation"
             Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O +
                    (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
             Maem = Phi_aem * Psat(T_des) / Paem * M_H2O +
                    (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
         else  # parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode"
             Masm = Phi_asm * Psat(T_des) / Pasm * M_H2O +
                    y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 +
                    (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
             Maem = Phi_aem * Psat(T_des) / Paem * M_H2O +
                    y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 +
                    (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
         # Cathode-side molar masses
         Mcsm = Phi_csm * Psat(T_des) / Pcsm * M_H2O +
                y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 +
                (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
         Mcem = Phi_cem * Psat(T_des) / Pcem * M_H2O +
                y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 +
                (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2

         # Gas-mixture density
         rho_asm = Pasm / (R * T_des) * Masm
         rho_aem = Paem / (R * T_des) * Maem
         rho_csm = Pcsm / (R * T_des) * Mcsm
         rho_cgc = Pcgc / (R * T_cgc) * Mcgc
         rho_cem = Pcem / (R * T_cgc) * Mcem

         # Purge
         if type_purge == "no_purge"
             k_purge = 0
         elseif type_purge == "constant_purge"
             k_purge = 1
         elseif type_purge == "periodic_purge"
             if (t - Int(floor(t / (purge_time + delta_purge))) * (purge_time + delta_purge)) <= purge_time
                 k_purge = 1
             else
                 k_purge = 0
             end
         else
             error("The type_purge variable must be correctly referenced.")
         end
         # Back-pressure valve area
         if Abp_a > A_T_a
             Abp_a = A_T_a
         elseif Abp_a < 0
             Abp_a = 0
         end
         if Abp_c > A_T_c
             Abp_c = A_T_c
         elseif Abp_c < 0
             Abp_c = 0
         end
=#
    else  # cfg.type_auxiliary == :no_auxiliary
        k_purge = nothing
        Abp_a = nothing
        Abp_c = nothing
    end

    return P_agc, P_cgc, Phi_agc, Phi_cgc, y_H2_agc, y_O2_cgc, M_agc, M_cgc, M_ext, M_H2_N2_in, rho_agc, rho_cgc,
           k_purge, Abp_a, Abp_c, mu_gaz_agc, mu_gaz_cgc
end

