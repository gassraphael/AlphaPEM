# -*- coding: utf-8 -*-

"""This file represents all the heat transfers occurring inside the fuel cell system.
It is a component of the fuel cell model.
"""

# ____________________________________________________Heat transfers____________________________________________________

"""In-place heat-transfer computation with reusable workspaces.

Parameters
----------
heat_work : MEAHeatWorkspace
    Pre-allocated workspace for thermal flux vectors.
heat_int_work : MEAHeatIntWorkspace
    Pre-allocated workspace for intermediate thermal conductivity values.
sv_1D : CellState1D{NB_GDL, NB_MPL}
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
function calculate_heat_transfers!(heat_work::MEAHeatWorkspace,
                                   heat_int_work::MEAHeatIntWorkspace,
                                   sv_1D::CellState1D{NB_GDL, NB_MPL},
                                   i_fc::Float64,
                                   fc::AbstractFuelCell,
                                   cfg::SimulationConfig,
                                   S_abs::MEASorptionSources,
                                   Sl::MEALiquidSources{NB_GDL, NB_MPL}
                                   )::MEAHeatFlows1D{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # ___________________________________________________Preliminaries__________________________________________________

    # Extraction of the variables (typed access via CellState1D struct fields)
    T_agc, T_agdl = sv_1D.agc.T, getproperty.(sv_1D.agdl, :T)
    T_ampl, T_acl = getproperty.(sv_1D.ampl, :T), sv_1D.acl.T
    T_mem = sv_1D.mem.T
    T_ccl, T_cmpl = sv_1D.ccl.T, getproperty.(sv_1D.cmpl, :T)
    T_cgdl, T_cgc = getproperty.(sv_1D.cgdl, :T), sv_1D.cgc.T

    s_agc, s_agdl = sv_1D.agc.s, getproperty.(sv_1D.agdl, :s)
    s_ampl, s_acl = getproperty.(sv_1D.ampl, :s), sv_1D.acl.s
    s_ccl, s_cmpl = sv_1D.ccl.s, getproperty.(sv_1D.cmpl, :s)
    s_cgdl, s_cgc = getproperty.(sv_1D.cgdl, :s), sv_1D.cgc.s

    lambda_acl, lambda_mem, lambda_ccl = sv_1D.acl.lambda, sv_1D.mem.lambda, sv_1D.ccl.lambda
    eta_c = sv_1D.ccl.eta_c

    # Extraction of the parameters
    pp = fc.physical_parameters
    T_des = fc.operating_conditions.T_des
    Hmem, Hgdl, Hmpl, Hacl, Hccl = pp.Hmem, pp.Hgdl, pp.Hmpl, pp.Hacl, pp.Hccl
    epsilon_gdl, epsilon_mpl = pp.epsilon_gdl, pp.epsilon_mpl

    # Intermediate values
    (Hgdl_node, Hmpl_node, k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_ampl, k_th_eff_ampl_ampl,
    k_th_eff_ampl_acl, k_th_eff_acl_mem, k_th_eff_mem_ccl, k_th_eff_ccl_cmpl, k_th_eff_cmpl_cmpl,
    k_th_eff_cmpl_cgdl, k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc) = calculate_heat_int_values!(heat_int_work, sv_1D, fc, cfg)

    # ______________________________________________Heat flows (J.m-2.s-1)______________________________________________

    # Anode side
    T_agc_mean = T_des
    T_cgc_mean = T_des
    Jt_agc_agdl = -k_th_eff_agc_agdl * d_dx(T_agc_mean, T_agdl[1], Hgdl_node / 2)
    Jt_agdl_agdl = heat_work.Jt_agdl_agdl
    @inbounds for i in 1:(NB_GDL - 1)
        Jt_agdl_agdl[i] = -k_th_eff_agdl_agdl[i] * d_dx(T_agdl[i], T_agdl[i + 1], Hgdl_node / 2)
    end
    Jt_agdl_ampl = -k_th_eff_agdl_ampl * d_dx(T_agdl[NB_GDL], T_ampl[1], Hgdl_node / 2, Hmpl_node / 2)
    Jt_ampl_ampl = heat_work.Jt_ampl_ampl
    @inbounds for i in 1:(NB_MPL - 1)
        Jt_ampl_ampl[i] = -k_th_eff_ampl_ampl[i] * d_dx(T_ampl[i], T_ampl[i + 1], Hmpl_node / 2)
    end
    Jt_ampl_acl = -k_th_eff_ampl_acl * d_dx(T_ampl[NB_MPL], T_acl, Hmpl_node / 2, Hacl / 2)

    # Membrane side
    Jt_acl_mem = -k_th_eff_acl_mem * d_dx(T_acl, T_mem, Hacl / 2, Hmem / 2)
    Jt_mem_ccl = -k_th_eff_mem_ccl * d_dx(T_mem, T_ccl, Hmem / 2, Hccl / 2)

    # Cathode side
    Jt_ccl_cmpl = -k_th_eff_ccl_cmpl * d_dx(T_ccl, T_cmpl[1], Hccl / 2, Hmpl_node / 2)
    Jt_cmpl_cmpl = heat_work.Jt_cmpl_cmpl
    @inbounds for i in 1:(NB_MPL - 1)
        Jt_cmpl_cmpl[i] = -k_th_eff_cmpl_cmpl[i] * d_dx(T_cmpl[i], T_cmpl[i + 1], Hmpl_node / 2)
    end
    Jt_cmpl_cgdl = -k_th_eff_cmpl_cgdl * d_dx(T_cmpl[NB_MPL], T_cgdl[1], Hmpl_node, Hgdl_node)
    Jt_cgdl_cgdl = heat_work.Jt_cgdl_cgdl
    @inbounds for i in 1:(NB_GDL - 1)
        Jt_cgdl_cgdl[i] = -k_th_eff_cgdl_cgdl[i] * d_dx(T_cgdl[i], T_cgdl[i + 1], Hgdl_node / 2)
    end
    Jt_cgdl_cgc = -k_th_eff_cgdl_cgc * d_dx(T_cgdl[NB_GDL], T_cgc_mean, Hgdl_node / 2)

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
        ntuple(i -> Sl.agdl[i] * (-delta_h_liq(T_agdl[i])), NB_GDL), # Q_liq_agdl
        ntuple(i -> Sl.ampl[i] * (-delta_h_liq(T_ampl[i])), NB_MPL), # Q_liq_ampl
        Sl.acl * (-delta_h_liq(T_acl)), # Q_liq_acl
        Sl.ccl * (-delta_h_liq(T_ccl)), # Q_liq_ccl
        ntuple(i -> Sl.cmpl[i] * (-delta_h_liq(T_cmpl[i])), NB_MPL), # Q_liq_cmpl
        ntuple(i -> Sl.cgdl[i] * (-delta_h_liq(T_cgdl[i])), NB_GDL) # Q_liq_cgdl
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

