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

