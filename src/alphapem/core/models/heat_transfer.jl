# -*- coding: utf-8 -*-

"""This file represents all the heat transfers occurring inside the fuel cell system.
It is a component of the fuel cell model.
"""

# ____________________________________________________Heat transfers____________________________________________________

"""Calculate the heat transfers occurring inside the fuel cell system.

Parameters
----------
sv_1D : MEAState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
i_fc : Float64
    Fuel cell current density at time t (A.m-2).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
S_abs : MEASorptionSources
    Typed water absorption rates from the CL to the membrane (mol.m-3.s-1).
Sl : MEALiquidSources{NB_GDL, NB_MPL}
    Typed liquid water phase-change source terms at each porous-layer node (mol.m-3.s-1).

Returns
-------
MEAHeatFlows1D{NB_GDL, NB_MPL}
    Typed heat transfer outputs (conductive fluxes and volumetric heat sources)
    for one gas-channel column.
"""
function calculate_heat_transfers(sv_1D::MEAState1D{NB_GDL, NB_MPL},
                                  i_fc::Float64,
                                  fc::AbstractFuelCell,
                                  S_abs::MEASorptionSources,
                                  Sl::MEALiquidSources{NB_GDL, NB_MPL}
                                  )::MEAHeatFlows1D{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of the variables (typed access via MEAState1D struct fields)
    T_acl, T_mem, T_ccl = sv_1D.acl.T, sv_1D.mem.T, sv_1D.ccl.T
    lambda_acl, lambda_mem, lambda_ccl = sv_1D.acl.lambda, sv_1D.mem.lambda, sv_1D.ccl.lambda
    s_acl, s_ccl, eta_c = sv_1D.acl.s, sv_1D.ccl.s, sv_1D.ccl.eta_c

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
    Jt_agc_agdl = -k_th_eff_agc_agdl * d_dx(T_agc_mean, sv_1D.agdl[1].T, Hgdl_node / 2)
    Jt_agdl_agdl = [-k_th_eff_agdl_agdl[i] * d_dx(sv_1D.agdl[i].T, sv_1D.agdl[i + 1].T, Hgdl_node / 2)
                    for i in 1:(nb_gdl - 1)]
    Jt_agdl_ampl = -k_th_eff_agdl_ampl * d_dx(sv_1D.agdl[nb_gdl].T, sv_1D.ampl[1].T, Hgdl_node / 2, Hmpl_node / 2)
    Jt_ampl_ampl = [-k_th_eff_ampl_ampl[i] * d_dx(sv_1D.ampl[i].T, sv_1D.ampl[i + 1].T, Hmpl_node / 2)
                    for i in 1:(nb_mpl - 1)]
    Jt_ampl_acl = -k_th_eff_ampl_acl * d_dx(sv_1D.ampl[nb_mpl].T, T_acl, Hmpl_node / 2, Hacl / 2)

    # Membrane side
    Jt_acl_mem = -k_th_eff_acl_mem * d_dx(T_acl, T_mem, Hacl / 2, Hmem / 2)
    Jt_mem_ccl = -k_th_eff_mem_ccl * d_dx(T_mem, T_ccl, Hmem / 2, Hccl / 2)

    # Cathode side
    Jt_ccl_cmpl = -k_th_eff_ccl_cmpl * d_dx(T_ccl, sv_1D.cmpl[1].T, Hccl / 2, Hmpl_node / 2)
    Jt_cmpl_cmpl = [-k_th_eff_cmpl_cmpl[i] * d_dx(sv_1D.cmpl[i].T, sv_1D.cmpl[i + 1].T, Hmpl_node / 2)
                    for i in 1:(nb_mpl - 1)]
    Jt_cmpl_cgdl = -k_th_eff_cmpl_cgdl * d_dx(sv_1D.cmpl[nb_mpl].T, sv_1D.cgdl[1].T, Hmpl_node, Hgdl_node)
    Jt_cgdl_cgdl = [-k_th_eff_cgdl_cgdl[i] * d_dx(sv_1D.cgdl[i].T, sv_1D.cgdl[i + 1].T, Hgdl_node / 2)
                    for i in 1:(nb_gdl - 1)]
    Jt_cgdl_cgc = -k_th_eff_cgdl_cgc * d_dx(sv_1D.cgdl[nb_gdl].T, T_cgc_mean, Hgdl_node / 2)

    Jt = MEAThermalFluxes{NB_GDL, NB_MPL}(
        Jt_agc_agdl, Jt_agdl_agdl, Jt_agdl_ampl, Jt_ampl_ampl, Jt_ampl_acl,
        Jt_acl_mem, Jt_mem_ccl, Jt_ccl_cmpl, Jt_cmpl_cmpl, Jt_cmpl_cgdl,
        Jt_cgdl_cgdl, Jt_cgdl_cgc
    )

    # ____________________________________________Heat generated (J.m-3.s-1)____________________________________________

    # The heat dissipated by the electrochemical reaction 2*H2 + O2 -> 2*H2O, in J.m-3.s-1.
    #    It is given by the sum of Peltier and activation heats [vetterFreeOpenReference2019].
    S_r_acl = i_fc / (2 * F * Hacl)  # mol.m-3.s-1. It is the amount of hydrogen consumed at the ACL.
    S_r_ccl = i_fc / (4 * F * Hccl)  # mol.m-3.s-1. It is the amount of oxygen consumed at the CCL.
    Q_r = MEAReactionHeat(
        S_r_acl * T_acl * (-delta_s_HOR), # Q_r_acl
        S_r_ccl * T_ccl * (-delta_s_ORR) + i_fc * eta_c / Hccl # Q_r_ccl (Peltier + activation heat)
    )

    # The heat dissipated by the absorption of water from the CL to the membrane, in J.m-3.s-1.
    Q_sorp = MEASorptionHeat(
        S_abs.v_acl * (-delta_h_abs(T_acl)), # Q_sorp_v_acl
        S_abs.l_acl * (-delta_h_abs(T_acl)), # Q_sorp_l_acl
        S_abs.v_ccl * (-delta_h_abs(T_ccl)), # Q_sorp_v_ccl
        S_abs.l_ccl * (-delta_h_abs(T_ccl)) # Q_sorp_l_ccl
    )

    # The heat dissipated by the liquefaction of vapor water, in J.m-3.s-1.
    Q_liq = MEALiquidHeat{NB_GDL, NB_MPL}(
        ntuple(i -> Sl.agdl[i] * (-delta_h_liq(sv_1D.agdl[i].T)), NB_GDL), # Q_liq_agdl
        ntuple(i -> Sl.ampl[i] * (-delta_h_liq(sv_1D.ampl[i].T)), NB_MPL), # Q_liq_ampl
        Sl.acl * (-delta_h_liq(T_acl)), # Q_liq_acl
        Sl.ccl * (-delta_h_liq(T_ccl)), # Q_liq_ccl
        ntuple(i -> Sl.cmpl[i] * (-delta_h_liq(sv_1D.cmpl[i].T)), NB_MPL), # Q_liq_cmpl
        ntuple(i -> Sl.cgdl[i] * (-delta_h_liq(sv_1D.cgdl[i].T)), NB_GDL) # Q_liq_cgdl
    )

    # The heat dissipated by the ionic currents (Joule heating + Ohm's law), in J.m-3.s-1.
    Q_p = MEAProtonHeat(
        i_fc^2 / sigma_p_eff("mem", lambda_mem, T_mem), # Q_p_mem
        i_fc^2 / (3 * sigma_p_eff("ccl", lambda_ccl, T_ccl, Hccl)) # Q_p_ccl
    )

    # The heat dissipated by the electric currents (Joule heating + Ohm's law), in J.m-3.s-1.
    sigma_e_eff_gdl = sigma_e_eff("gdl", epsilon_gdl, epsilon_c)
    sigma_e_eff_mpl = sigma_e_eff("mpl", epsilon_mpl)
    Q_e = MEAElectricHeat{NB_GDL, NB_MPL}(
        ntuple(_ -> i_fc^2 / sigma_e_eff_gdl, NB_GDL), # Q_e_agdl
        ntuple(_ -> i_fc^2 / sigma_e_eff_mpl, NB_MPL), # Q_e_ampl
        i_fc^2 / sigma_e_eff("cl", nothing, nothing, lambda_acl, T_acl, Hacl), # Q_e_acl
        i_fc^2 / (3 * sigma_e_eff("cl", nothing, nothing, lambda_ccl, T_ccl, Hccl)), # Q_e_ccl
        ntuple(_ -> i_fc^2 / sigma_e_eff_mpl, NB_MPL), # Q_e_cmpl
        ntuple(_ -> i_fc^2 / sigma_e_eff_gdl, NB_GDL) # Q_e_cgdl
    )

    return MEAHeatFlows1D{NB_GDL, NB_MPL}(Jt, Q_r, Q_sorp, Q_liq, Q_p, Q_e)
end
