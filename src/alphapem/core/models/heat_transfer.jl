# -*- coding: utf-8 -*-

"""This file represents all the heat transfers occurring inside the fuel cell system.
It is a component of the fuel cell model.
"""

# ____________________________________________________Heat transfers____________________________________________________

"""Calculate the heat transfers occurring inside the fuel cell system.

Parameters
----------
sv_1D : Dict
    Variables calculated by the solver (fuel cell internal states).
    `sv` is a contraction of solver_variables for enhanced readability.
i_fc :
    Fuel cell current density at time t (A.m-2).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
S_abs : Dict
    Water absorption rates from the CL to the membrane (mol.m-3.s-1).
Sl : Dict
    Liquid water absorption rates (mol.m-3.s-1).

Returns
-------
Dict{String, Dict} where
    Heat transfers occurring inside the fuel cell system.
"""
function calculate_heat_transfers(sv_1D::Dict, i_fc, fc::AbstractFuelCell, S_abs::Dict, Sl::Dict)::Dict{String, Dict}

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of the variables
    T_acl, T_mem, T_ccl = sv_1D["T_acl"], sv_1D["T_mem"], sv_1D["T_ccl"]
    lambda_acl, lambda_mem, lambda_ccl = sv_1D["lambda_acl"], sv_1D["lambda_mem"], sv_1D["lambda_ccl"]
    s_acl, s_ccl, eta_c = sv_1D["s_acl"], sv_1D["s_ccl"], sv_1D["eta_c"]

    # Extraction of the parameters
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des = fc.operating_conditions.T_des
    Hmem, Hgdl, Hmpl, Hacl, Hccl = pp.Hmem, pp.Hgdl, pp.Hmpl, pp.Hacl, pp.Hccl
    epsilon_gdl, epsilon_mpl, epsilon_c = pp.epsilon_gdl, pp.epsilon_mpl, pp.epsilon_c
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Intermediate values
    (Hgdl_node, Hmpl_node, k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_ampl, k_th_eff_ampl_ampl,
     k_th_eff_ampl_acl, k_th_eff_acl_mem, k_th_eff_mem_ccl, k_th_eff_ccl_cmpl, k_th_eff_cmpl_cmpl,
     k_th_eff_cmpl_cgdl, k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc) = heat_transfer_int_values(sv_1D, fc)

    # ______________________________________________Heat flows (J.m-2.s-1)______________________________________________

    # Anode side
    T_agc_mean = T_des
    T_cgc_mean = T_des
    Jt_agc_agdl = -k_th_eff_agc_agdl * d_dx(T_agc_mean, sv_1D["T_agdl_1"], Hgdl_node / 2)
    Jt_agdl_agdl = Dict(
        "agdl_agdl_$i" => -k_th_eff_agdl_agdl[i] * d_dx(sv_1D["T_agdl_$i"], sv_1D["T_agdl_$(i + 1)"], Hgdl_node / 2)
        for i in 1:(nb_gdl - 1)
    )

    Jt_agdl_ampl = -k_th_eff_agdl_ampl * d_dx(sv_1D["T_agdl_$nb_gdl"], sv_1D["T_ampl_1"], Hgdl_node / 2, Hmpl_node / 2)
    Jt_ampl_ampl = Dict(
        "ampl_ampl_$i" => -k_th_eff_ampl_ampl[i] * d_dx(sv_1D["T_ampl_$i"], sv_1D["T_ampl_$(i + 1)"], Hmpl_node / 2)
        for i in 1:(nb_mpl - 1)
    )
    Jt_ampl_acl = -k_th_eff_ampl_acl * d_dx(sv_1D["T_ampl_$nb_mpl"], T_acl, Hmpl_node / 2, Hacl / 2)

    # Membrane side
    Jt_acl_mem = -k_th_eff_acl_mem * d_dx(T_acl, T_mem, Hacl / 2, Hmem / 2)
    Jt_mem_ccl = -k_th_eff_mem_ccl * d_dx(T_mem, T_ccl, Hmem / 2, Hccl / 2)

    # Cathode side
    Jt_ccl_cmpl = -k_th_eff_ccl_cmpl * d_dx(T_ccl, sv_1D["T_cmpl_1"], Hccl / 2, Hmpl_node / 2)
    Jt_cmpl_cmpl = Dict(
        "cmpl_cmpl_$i" => -k_th_eff_cmpl_cmpl[i] * d_dx(sv_1D["T_cmpl_$i"], sv_1D["T_cmpl_$(i + 1)"], Hmpl_node / 2)
        for i in 1:(nb_mpl - 1)
    )
    Jt_cmpl_cgdl = -k_th_eff_cmpl_cgdl * d_dx(sv_1D["T_cmpl_$nb_mpl"], sv_1D["T_cgdl_1"], Hmpl_node, Hgdl_node)
    Jt_cgdl_cgdl = Dict(
        "cgdl_cgdl_$i" => -k_th_eff_cgdl_cgdl[i] * d_dx(sv_1D["T_cgdl_$i"], sv_1D["T_cgdl_$(i + 1)"], Hgdl_node / 2)
        for i in 1:(nb_gdl - 1)
    )

    Jt_cgdl_cgc = -k_th_eff_cgdl_cgc * d_dx(sv_1D["T_cgdl_$nb_gdl"], T_cgc_mean, Hgdl_node / 2)

    Jt = merge(
        Dict("agc_agdl" => Jt_agc_agdl),
        Jt_agdl_agdl,
        Dict("agdl_ampl" => Jt_agdl_ampl),
        Jt_ampl_ampl,
        Dict("ampl_acl" => Jt_ampl_acl,
                       "acl_mem" => Jt_acl_mem,
                       "mem_ccl" => Jt_mem_ccl,
                       "ccl_cmpl" => Jt_ccl_cmpl),
        Jt_cmpl_cmpl,
        Dict("cmpl_cgdl" => Jt_cmpl_cgdl),
        Jt_cgdl_cgdl,
        Dict("cgdl_cgc" => Jt_cgdl_cgc)
    )

    # ____________________________________________Heat generated (J.m-3.s-1)____________________________________________

    # The heat dissipated by the electrochemical reaction 2*H2 + O2 -> 2*H2O, in J.m-3.s-1.
    #    It is given by the sum of Peltier and activation heats [vetterFreeOpenReference2019].
    S_r_acl = i_fc / (2 * F * Hacl)  # mol.m-3.s-1. It is the amount of hydrogen consumed at the ACL.
    S_r_ccl = i_fc / (4 * F * Hccl)  # mol.m-3.s-1. It is the amount of oxygen consumed at the CCL.
    Q_r = Dict(
        "acl" => S_r_acl * T_acl * (-delta_s_HOR),
        "ccl" => S_r_ccl * T_ccl * (-delta_s_ORR) + i_fc * eta_c / Hccl
    )

    # The heat dissipated by the absorption of water from the CL to the membrane, in J.m-3.s-1.
    Q_sorp = Dict(
        "v_acl" => S_abs["v_acl"] * (-delta_h_abs(T_acl)),
        "l_acl" => S_abs["l_acl"] * (-delta_h_abs(T_acl)),
        "v_ccl" => S_abs["v_ccl"] * (-delta_h_abs(T_ccl)),
        "l_ccl" => S_abs["l_ccl"] * (-delta_h_abs(T_ccl))
    )

    # The heat dissipated by the liquefaction of vapor water, in J.m-3.s-1.
    Q_liq = merge(
        Dict(
            "agdl_$i" => Sl["agdl"][i] * (-delta_h_liq(sv_1D["T_agdl_$i"]))
            for i in 1:nb_gdl
        ),
        Dict(
            "cgdl_$i" => Sl["cgdl"][i] * (-delta_h_liq(sv_1D["T_cgdl_$i"]))
            for i in 1:nb_gdl
        ),
        Dict(
            "ampl_$i" => Sl["ampl"][i] * (-delta_h_liq(sv_1D["T_ampl_$i"]))
            for i in 1:nb_mpl
        ),
        Dict(
            "cmpl_$i" => Sl["cmpl"][i] * (-delta_h_liq(sv_1D["T_cmpl_$i"]))
            for i in 1:nb_mpl
        ),
        Dict(
            "acl" => Sl["acl"] * (-delta_h_liq(T_acl)),
            "ccl" => Sl["ccl"] * (-delta_h_liq(T_ccl))
        )
    )

    # The heat dissipated by the ionic currents (Joule heating + Ohm's law), in J.m-3.s-1.
    Q_p = Dict(
        "mem" => i_fc^2 / sigma_p_eff("mem", lambda_mem, T_mem),
        "ccl" => i_fc^2 / (3 * sigma_p_eff("ccl", lambda_ccl, T_ccl, Hccl))
    )

    # The heat dissipated by the electric currents (Joule heating + Ohm's law), in J.m-3.s-1.
    sigma_e_eff_gdl = sigma_e_eff("gdl", epsilon_gdl, epsilon_c)
    sigma_e_eff_mpl = sigma_e_eff("mpl", epsilon_mpl)
    Q_e = merge(
        Dict(
            "agdl_$i" => i_fc^2 / sigma_e_eff_gdl
            for i in 1:nb_gdl
        ),
        Dict(
            "ampl_$i" => i_fc^2 / sigma_e_eff_mpl
            for i in 1:nb_mpl
        ),
        Dict(
            "acl" => i_fc^2 / sigma_e_eff("cl", nothing, nothing, lambda_acl, T_acl, Hacl),
            "ccl" => i_fc^2 / (3 * sigma_e_eff("cl", nothing, nothing, lambda_ccl, T_ccl, Hccl))
        ),
        Dict(
            "cmpl_$i" => i_fc^2 / sigma_e_eff_mpl
            for i in 1:nb_mpl
        ),
        Dict(
            "cgdl_$i" => i_fc^2 / sigma_e_eff_gdl
            for i in 1:nb_gdl
        )
    )

    return Dict(
        "Jt" => Jt,
        "Q_r" => Q_r,
        "Q_sorp" => Q_sorp,
        "Q_liq" => Q_liq,
        "Q_p" => Q_p,
        "Q_e" => Q_e
    )
end
