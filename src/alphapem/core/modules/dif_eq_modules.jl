# -*- coding: utf-8 -*-

"""This module is used to determine intermediate values for the calculation of the differential equations
and to implement integration events.
"""

# ____________________________________________Differential equations modules____________________________________________

"""This function calculates intermediate values for the calculation of the differential equations.

Parameters
----------
t : Float64
    Time (s).
sv : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cfg : SimulationConfig
    Simulation configuration.

Returns
-------
Dict
    Dictionary containing all intermediate values used by the differential equations.
"""
function calculate_dif_eq_int_values(t::Float64, sv::Dict, fc::AbstractFuelCell, cfg::SimulationConfig)::Dict

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv["C_v_agc"], sv["C_v_acl"], sv["C_v_ccl"], sv["C_v_cgc"]
    s_acl, s_ccl = sv["s_acl"], sv["s_ccl"]
    lambda_acl, lambda_mem, lambda_ccl = sv["lambda_acl"], sv["lambda_mem"], sv["lambda_ccl"]
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv["C_H2_agc"], sv["C_H2_acl"], sv["C_O2_ccl"], sv["C_O2_cgc"]
    C_N2_agc, C_N2_cgc = sv["C_N2_agc"], sv["C_N2_cgc"]
    T_agc, T_acl, T_mem, T_ccl, T_cgc = sv["T_agc"], sv["T_acl"], sv["T_mem"], sv["T_ccl"], sv["T_cgc"]
    Pasm, Paem, Pcsm, Pcem = get(sv, "Pasm", nothing), get(sv, "Paem", nothing), get(sv, "Pcsm", nothing), get(sv, "Pcem", nothing)
    Phi_asm, Phi_aem = get(sv, "Phi_asm", nothing), get(sv, "Phi_aem", nothing)
    Phi_csm, Phi_cem = get(sv, "Phi_csm", nothing), get(sv, "Phi_cem", nothing)

    # Extraction of the parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des, y_H2_in = oc.T_des, oc.y_H2_in
    Hmem, Hacl, Hccl = pp.Hmem, pp.Hacl, pp.Hccl
    Lgc, nb_channel_in_gc, Lm = pp.Lgc, pp.nb_channel_in_gc, pp.Lm
    epsilon_gdl, epsilon_mpl, kappa_co = pp.epsilon_gdl, pp.epsilon_mpl, pp.kappa_co
    nb_gdl, nb_mpl, purge_time, delta_purge = np.nb_gdl, np.nb_mpl, np.purge_time, np.delta_purge
    type_purge = cfg.type_purge

    # Physical quantities outside the stack
    #       Molar masses
    M_ext = Phi_ext * Psat(Text) / Pext * M_H2O +
            y_O2_ext * (1 - Phi_ext * Psat(Text) / Pext) * M_O2 +
            (1 - y_O2_ext) * (1 - Phi_ext * Psat(Text) / Pext) * M_N2
    M_H2_N2_in = y_H2_in * M_H2 + (1 - y_H2_in) * M_N2

    # Physical quantities inside the stack
    #       Pressures
    P_agc = (C_v_agc + C_H2_agc + C_N2_agc) * R * T_agc
    P_cgc = (C_v_cgc + C_O2_cgc + C_N2_cgc) * R * T_cgc

    #       Total concentration
    C_tot_agc = C_v_agc + C_H2_agc + C_N2_agc
    C_tot_cgc = C_v_cgc + C_O2_cgc + C_N2_cgc

    #       Humidities
    Phi_cgc = C_v_cgc / C_v_sat(T_cgc)

    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC
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
    x_H2O_v_agc = C_v_agc / (C_v_agc + C_H2_agc + C_N2_agc)
    x_H2O_v_cgc = C_v_cgc / (C_v_cgc + C_O2_cgc + C_N2_cgc)

    #       Dynamic viscosity of the gas mixture.
    mu_gaz_agc = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_agc, 1 - x_H2O_v_agc], T_agc)
    mu_gaz_cgc = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                  [x_H2O_v_cgc, y_O2_cgc * (1 - x_H2O_v_cgc), (1 - y_O2_cgc) * (1 - x_H2O_v_cgc)],
                                  T_cgc)

    #       Volumetric heat capacity (J.m-3.K-1)
    rho_Cp0 = Dict()

    for i in 1:nb_gdl
        rho_Cp0["agdl_$i"] = calculate_rho_Cp0("agdl", sv["T_agdl_$i"], sv["C_v_agdl_$i"],
                               sv["s_agdl_$i"], nothing, sv["C_H2_agdl_$i"], nothing, C_N2_agc,
                               epsilon_gdl)
    end

    for i in 1:nb_mpl
        rho_Cp0["ampl_$i"] = calculate_rho_Cp0("ampl", sv["T_ampl_$i"], sv["C_v_ampl_$i"],
                               sv["s_ampl_$i"], nothing, sv["C_H2_ampl_$i"], nothing, C_N2_agc,
                               epsilon_mpl)
    end

    rho_Cp0["acl"] = calculate_rho_Cp0("acl", T_acl, C_v_acl, s_acl, lambda_acl, C_H2_acl, nothing, C_N2_agc,
                                       nothing, Hacl)
    rho_Cp0["mem"] = calculate_rho_Cp0("mem", T_mem, nothing, nothing, lambda_mem)
    rho_Cp0["ccl"] = calculate_rho_Cp0("ccl", T_ccl, C_v_ccl, s_ccl, lambda_ccl, nothing, C_O2_ccl,
                                       C_N2_cgc, nothing, Hccl)

    for i in 1:nb_mpl
        rho_Cp0["cmpl_$i"] = calculate_rho_Cp0("cmpl", sv["T_cmpl_$i"], sv["C_v_cmpl_$i"],
                               sv["s_cmpl_$i"], nothing, nothing, sv["C_O2_cmpl_$i"], C_N2_cgc,
                               epsilon_mpl)
    end

    for i in 1:nb_gdl
        rho_Cp0["cgdl_$i"] = calculate_rho_Cp0("cgdl", sv["T_cgdl_$i"], sv["C_v_cgdl_$i"],
                               sv["s_cgdl_$i"], nothing, nothing, sv["C_O2_cgdl_$i"], C_N2_cgc,
                               epsilon_gdl)
    end

    #       The crossover current density i_n
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                            [Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
    i_n = i_H2 + i_O2

    if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode

        #=
        # Purge
        if type_purge == "no_purge"
            k_purge = 0
        elseif type_purge == "constant_purge"
            k_purge = 1
        elseif type_purge == "periodic_purge"
            if (t - floor(Int, t / (purge_time + delta_purge)) * (purge_time + delta_purge)) <= purge_time
                k_purge = 1
            else
                k_purge = 0
            end
        else
            throw(ArgumentError("The type_purge variable should be correctly referenced."))
        end

        # H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the EM
        y_H2_aem = (Paem - Phi_aem * Psat(T_des) - C_N2_a * R * T_des) / (Paem - Phi_aem * Psat(T_des))
        y_O2_cem = (Pcem - Phi_cem * Psat(T_cgc) - C_N2_c * R * T_cgc) / (Pcem - Phi_cem * Psat(T_cgc))

        # Molar masses at the anode side
        if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            M["asm"] = Phi_asm * Psat(T_des) / Pasm * M_H2O +
                       (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2
            M["aem"] = Phi_aem * Psat(T_des) / Paem * M_H2O +
                       (1 - Phi_aem * Psat(T_des) / Paem) * M_H2
        else # cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
            M["asm"] = Phi_asm * Psat(T_des) / Pasm * M_H2O +
                       y_H2_in * (1 - Phi_asm * Psat(T_des) / Pasm) * M_H2 +
                       (1 - y_H2_in) * (1 - Phi_asm * Psat(T_des) / Pasm) * M_N2
            M["aem"] = Phi_aem * Psat(T_des) / Paem * M_H2O +
                       y_H2_aem * (1 - Phi_aem * Psat(T_des) / Paem) * M_H2 +
                       (1 - y_H2_aem) * (1 - Phi_aem * Psat(T_des) / Paem) * M_N2
        end

        # Molar masses at the cathode side
        M["csm"] = Phi_csm * Psat(T_des) / Pcsm * M_H2O +
                   y_O2_ext * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_O2 +
                   (1 - y_O2_ext) * (1 - Phi_csm * Psat(T_des) / Pcsm) * M_N2
        M["cem"] = Phi_cem * Psat(T_des) / Pcem * M_H2O +
                   y_O2_cem * (1 - Phi_cem * Psat(T_des) / Pcem) * M_O2 +
                   (1 - y_O2_cem) * (1 - Phi_cem * Psat(T_des) / Pcem) * M_N2

        # Density/concentration of the gas mixture.
        C_tot_a_in = Pasm_in / (R * T_des)
        rho_asm = Pasm / (R * T_des) * Masm
        rho_agc = P["agc_$i"] / (R * sv["T_agc_$i"]) * Magc
        rho_aem = Paem / (R * T_des) * Maem
        if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            rho_asm_in_re = Pasm_in_re / (R * T_des) * Masm_in_re
            rho_aem_out_re = Paem_out_re / (R * T_des) * Maem_out_re
        else
            rho_asm_in_re, rho_aem_out_re = nothing, nothing
        end
        rho_a_ext = Pext / (R * T_des) * Maem_out
        C_tot_a_ext = Pext / (R * T_des) # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed.
        C_tot_c_in = Pcsm_in / (R * T_des)
        rho_csm = Pcsm / (R * T_des) * Mcsm
        rho_cgc = P["cgc_$i"] / (R * sv["T_cgc_$i"]) * Mcgc
        rho_cem = Pcem / (R * T_cgc) * Mcem
        rho_c_ext = Pext / (R * T_des) * Mcem_out
        C_tot_c_ext = Pext * Mcem_out / (R * T_des) # Boundary condition: at the exit, pressure and temperature are fixed. So, the total concentration is fixed.

        # Vapor ratio over the gas mixture.
        x_H2O_v_asm = Phi_asm * Psat(T_des) / Pasm
        x_H2O_v_agc = C_v_agc / (C_v_agc + C_H2_agc + C_N2_a)
        x_H2O_v_aem = Phi_aem * Psat(T_des) / Paem
        x_H2O_v_a_ext = Phi_a_ext * Psat(T_des) / Pext
        x_H2O_v_csm = Phi_csm * Psat(T_des) / Pcsm
        x_H2O_v_cgc = C_v_cgc / (C_v_cgc + C_O2_cgc + C_N2_c)
        x_H2O_v_cem = Phi_cem * Psat(T_des) / Pcem
        x_H2O_v_c_ext = Phi_c_ext * Psat(T_des) / Pext

        # Molar fraction of H2 in the dry gas mixture (H2/N2)
        y_H2_agc = C_H2_agc / (C_H2_agc + C_N2_a)
        y_O2_cgc = C_O2_cgc / (C_O2_cgc + C_N2_c)

        # Dynamic viscosity of the gas mixture at the anode side.
        if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            mu_gaz_asm = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_asm, 1 - x_H2O_v_asm], T_des)
            mu_gaz_agc = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_agc, 1 - x_H2O_v_agc], T_agc)
            mu_gaz_aem = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_aem, 1 - x_H2O_v_aem], T_des)
            mu_gaz_a_ext = mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_a_ext, 1 - x_H2O_v_a_ext], T_des)
        else # cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
            mu_gaz_asm = mu_mixture_gases(["H2O_v", "H2", "N2"],
                                          [x_H2O_v_asm, y_H2_in * (1 - x_H2O_v_asm), (1 - y_H2_in) * (1 - x_H2O_v_asm)],
                                          T_des)
            mu_gaz_agc = mu_mixture_gases(["H2O_v", "H2", "N2"],
                                          [x_H2O_v_agc, y_H2_agc * (1 - x_H2O_v_agc), (1 - y_H2_agc) * (1 - x_H2O_v_agc)],
                                          T_agc)
            mu_gaz_aem = mu_mixture_gases(["H2O_v", "H2", "N2"],
                                          [x_H2O_v_aem, y_H2_aem * (1 - x_H2O_v_aem), (1 - y_H2_aem) * (1 - x_H2O_v_aem)],
                                          T_des)
            mu_gaz_a_ext = mu_mixture_gases(["H2O_v", "H2", "N2"],
                                            [x_H2O_v_a_ext, y_H2_aem_out * (1 - x_H2O_v_a_ext), (1 - y_H2_aem_out) * (1 - x_H2O_v_a_ext)],
                                            T_des)
        end

        # Dynamic viscosity of the gas mixture at the cathode side.
        mu_gaz_csm = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                      [x_H2O_v_csm, y_O2_ext * (1 - x_H2O_v_csm), (1 - y_O2_ext) * (1 - x_H2O_v_csm)],
                                      T_des)
        mu_gaz_cgc = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                      [x_H2O_v_cgc, y_O2_cgc * (1 - x_H2O_v_cgc), (1 - y_O2_cgc) * (1 - x_H2O_v_cgc)],
                                      T_cgc)
        mu_gaz_cem = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                      [x_H2O_v_cem, y_O2_cem * (1 - x_H2O_v_cem), (1 - y_O2_cem) * (1 - x_H2O_v_cem)],
                                      T_des)
        mu_gas_c_ext = mu_mixture_gases(["H2O_v", "O2", "N2"],
                                        [x_H2O_v_c_ext, y_O2_cem_out * (1 - x_H2O_v_c_ext), (1 - y_O2_cem_out) * (1 - x_H2O_v_c_ext)],
                                        T_des)

        # Boundary velocities
        if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
            v_re = Ware / rho_aem_out_re / A_T_a
        else # cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
            v_re = nothing
        end
        =#

    else # cfg.type_auxiliary == :no_auxiliary
        # Set to `nothing` the variables not used when "no_auxiliary" system is considered.
        v_re, Lman_to_endplate, Lman_to_man_gc, k_purge = (nothing, nothing, nothing, nothing)
    end

    return Dict(
        "rho_Cp0" => rho_Cp0,
        "v_re" => v_re,
        "k_purge" => k_purge,
        "rho_agc" => rho_agc,
        "rho_cgc" => rho_cgc,
        "C_tot_agc" => C_tot_agc,
        "C_tot_cgc" => C_tot_cgc,
        "mu_gaz_agc" => mu_gaz_agc,
        "mu_gaz_cgc" => mu_gaz_cgc,
        "P_agc" => P_agc,
        "P_cgc" => P_cgc,
        "i_n" => i_n,
    )
end


# Typed solver-vector helpers depend on model structs (`MEAState1D`, manifold and
# auxiliary containers). They are loaded only in the `Models` include context.
if @isdefined(MEAState1D)

"""Count the number of MEA/GC variables per gas-channel node in the solver vector.

The ordering matches `create_initial_variable_values` in `AlphaPEM.jl`.
"""
function _nb_solver_vars_per_gc(nb_gdl::Int, nb_mpl::Int)::Int
    return fieldcount(AnodeGCNode) +
           nb_gdl * fieldcount(AnodeGDLNode) +
           nb_mpl * fieldcount(AnodeMPLNode) +
           fieldcount(AnodeCLNode) +
           fieldcount(MembraneNode) +
           fieldcount(CathodeCLNode) +
           nb_mpl * fieldcount(CathodeMPLNode) +
           nb_gdl * fieldcount(CathodeGDLNode) +
           fieldcount(CathodeGCNode)
end


"""Unpack one GC-node segment of the solver vector into a typed 1D MEA state."""
function _unpack_mea_state_1D(values::AbstractVector{<:Real}, nb_gdl::Int, nb_mpl::Int)
    idx = 1

    # C_v block
    C_v_agc = Float64(values[idx]); idx += 1
    C_v_agdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    C_v_ampl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    C_v_acl = Float64(values[idx]); idx += 1
    C_v_ccl = Float64(values[idx]); idx += 1
    C_v_cmpl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    C_v_cgdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    C_v_cgc = Float64(values[idx]); idx += 1

    # s block
    s_agc = Float64(values[idx]); idx += 1
    s_agdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    s_ampl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    s_acl = Float64(values[idx]); idx += 1
    s_ccl = Float64(values[idx]); idx += 1
    s_cmpl = [Float64(values[idx + i - 1]) for i in 1:nb_mpl]; idx += nb_mpl
    s_cgdl = [Float64(values[idx + i - 1]) for i in 1:nb_gdl]; idx += nb_gdl
    s_cgc = Float64(values[idx]); idx += 1

    # lambda block
    lambda_acl = Float64(values[idx]); idx += 1
    lambda_mem = Float64(values[idx]); idx += 1
    lambda_ccl = Float64(values[idx]); idx += 1

    # H2/O2/N2 block
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

    # T block + eta_c
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

    # Defensive check: ensure the full segment was consumed exactly once.
    idx == length(values) + 1 ||
        throw(ArgumentError("Invalid 1D state segment length while unpacking solver vector."))

    agc = AnodeGCNode(T_agc, C_v_agc, s_agc, C_H2_agc, C_N2_agc)
    agdl = ntuple(i -> AnodeGDLNode(T_agdl[i], C_v_agdl[i], s_agdl[i], C_H2_agdl[i]), nb_gdl)
    ampl = ntuple(i -> AnodeMPLNode(T_ampl[i], C_v_ampl[i], s_ampl[i], C_H2_ampl[i]), nb_mpl)
    acl = AnodeCLNode(T_acl, C_v_acl, s_acl, lambda_acl, C_H2_acl)
    mem = MembraneNode(T_mem, lambda_mem)
    ccl = CathodeCLNode(T_ccl, C_v_ccl, s_ccl, lambda_ccl, C_O2_ccl, eta_c)
    cmpl = ntuple(i -> CathodeMPLNode(T_cmpl[i], C_v_cmpl[i], s_cmpl[i], C_O2_cmpl[i]), nb_mpl)
    cgdl = ntuple(i -> CathodeGDLNode(T_cgdl[i], C_v_cgdl[i], s_cgdl[i], C_O2_cgdl[i]), nb_gdl)
    cgc = CathodeGCNode(T_cgc, C_v_cgc, s_cgc, C_O2_cgc, C_N2_cgc)

    return MEAState1D{nb_gdl, nb_mpl}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end


"""Create a derivative container initialized with NaN values.

NaN sentinels make missing derivative assignments fail fast.
"""
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
    return MEADerivative1D{nb_gdl, nb_mpl}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end


"""Ensure all derivatives were assigned before returning to the solver."""
function _assert_derivative_complete(d::MEADerivative1D{nb_gdl, nb_mpl}) where {nb_gdl, nb_mpl}
    any(isnan, _pack_mea_derivative_1D(d)) &&
        throw(ArgumentError("At least one derivative entry is missing (NaN sentinel detected)."))
    return nothing
end


"""Repack one typed 1D derivative container into the solver ordering."""
function _pack_mea_derivative_1D(d::MEADerivative1D{nb_gdl, nb_mpl}) where {nb_gdl, nb_mpl}
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
    return length(manifold_lines) * nb_man * fieldcount(ManifoldNode)
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


"""Unpack manifold state from a solver-vector segment.

Ordering is:
P(asm, aem, csm, cem) then Phi(asm, aem, csm, cem), each by local manifold node index.
"""
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

    asm = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldNode(P_asm[i], Phi_asm[i]), nb_manifold))
    aem = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldNode(P_aem[i], Phi_aem[i]), nb_manifold))
    csm = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldNode(P_csm[i], Phi_csm[i]), nb_manifold))
    cem = ManifoldLine{nb_manifold}(ntuple(i -> ManifoldNode(P_cem[i], Phi_cem[i]), nb_manifold))

    return _ManifoldStateBundle{nb_manifold}(asm, aem, csm, cem)
end


"""Create manifold derivative container initialised with NaN sentinels."""
function _nan_manifold_derivative_state(nb_manifold::Int)
    z = NaN
    mkline(n) = ManifoldLineDerivative{n}(ntuple(_ -> ManifoldDerivative(z, z), n))
    return _ManifoldDerivativeBundle{nb_manifold}(mkline(nb_manifold), mkline(nb_manifold),
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

