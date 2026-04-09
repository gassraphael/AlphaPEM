# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# ____________________________________________________Assembly__________________________________________________________

"""Assemble one complete MEA derivative from per-physics derivative contributions."""
function assemble_mea_derivative_1D(dw::MEADissolvedWaterDerivative,
                                    lw::MEALiquidWaterDerivative{NB_GDL, NB_MPL},
                                    vw::MEAVaporDerivative{NB_GDL, NB_MPL},
                                    sd::MEAH2O2SpeciesDerivative{NB_GDL, NB_MPL},
                                    vd::MEAVoltageDerivative,
                                    td::MEATemperatureDerivative{NB_GDL, NB_MPL}) where {NB_GDL, NB_MPL}

    z = NaN

    # Gas-channel derivatives are still filled by gas-channel equations later.
    agc = AnodeGCDerivative(z, z, z, z, z)
    cgc = CathodeGCDerivative(z, z, z, z, z)

    agdl = ntuple(NB_GDL) do j
        AnodeGDLDerivative(vw.agdl_C_v[j], lw.agdl_s[j], sd.agdl_C_H2[j], td.agdl_T[j])
    end
    ampl = ntuple(NB_MPL) do j
        AnodeMPLDerivative(vw.ampl_C_v[j], lw.ampl_s[j], sd.ampl_C_H2[j], td.ampl_T[j])
    end
    acl = AnodeCLDerivative(vw.acl_C_v, lw.acl_s, sd.acl_C_H2, dw.acl_lambda, td.acl_T)
    mem = MembraneDerivative(dw.mem_lambda, td.mem_T)
    ccl = CathodeCLDerivative(vw.ccl_C_v, lw.ccl_s, sd.ccl_C_O2, dw.ccl_lambda, td.ccl_T, vd.ccl_eta_c)
    cmpl = ntuple(NB_MPL) do j
        CathodeMPLDerivative(vw.cmpl_C_v[j], lw.cmpl_s[j], sd.cmpl_C_O2[j], td.cmpl_T[j])
    end
    cgdl = ntuple(NB_GDL) do j
        CathodeGDLDerivative(vw.cgdl_C_v[j], lw.cgdl_s[j], sd.cgdl_C_O2[j], td.cgdl_T[j])
    end

    return MEACellDerivative1D{NB_GDL, NB_MPL}(agc, agdl, ampl, acl, mem, ccl, cmpl, cgdl, cgc)
end

# ____________________________________________________Main functions____________________________________________________

"""Calculate dissolved-water (lambda) dynamics contribution.

Parameters
----------
sv : MEAState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
S_abs : MEASorptionSources
    Water absorption rates at the CLs (mol·m⁻³·s⁻¹).
J_lambda : MEADissolvedWaterFlux
    Dissolved-water inter-layer fluxes (kg·m⁻²·s⁻¹).
Sp : MEAWaterProductionSources
    Water production rates at the CLs (mol·m⁻³·s⁻

Returns
-------
MEADissolvedWaterDerivative
    Container with the lambda derivatives for the ACL, membrane and CCL.
"""
function calculate_dyn_dissoved_water_evolution_inside_MEA(
        sv::MEAState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        S_abs::MEASorptionSources,
        J_lambda::MEADissolvedWaterFlux,
        Sp::MEAWaterProductionSources
)::MEADissolvedWaterDerivative where {NB_GDL, NB_MPL}

    dlambda_acl = M_eq / (rho_mem * epsilon_mc(sv.acl.lambda, sv.acl.T, pp.Hacl)) *
                  (-J_lambda.acl_mem / pp.Hacl + S_abs.v_acl + S_abs.l_acl + Sp.acl)
    dlambda_mem = M_eq / rho_mem * (J_lambda.acl_mem - J_lambda.mem_ccl) / pp.Hmem
    dlambda_ccl = M_eq / (rho_mem * epsilon_mc(sv.ccl.lambda, sv.ccl.T, pp.Hccl)) *
                  (J_lambda.mem_ccl / pp.Hccl + S_abs.v_ccl + S_abs.l_ccl + Sp.ccl)

    return MEADissolvedWaterDerivative(dlambda_acl, dlambda_mem, dlambda_ccl)
end


"""Calculate the dynamic evolution of liquid water in the porous layers.

Parameters
----------
sv : MEAState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
Jl : MEALiquidFluxes{NB_GDL, NB_MPL}
    Liquid-water inter-layer fluxes (kg·m⁻²·s⁻¹).
S_abs : MEASorptionSources
    Water absorption rates at the CLs (mol·m⁻³·s⁻¹).
Sl : MEALiquidSources{NB_GDL, NB_MPL}
    Liquid-water phase-change source terms (mol·m⁻³·s⁻¹).

Returns
-------
MEACellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with s (liquid saturation) derivatives filled in.
"""
function calculate_dyn_liquid_water_evolution_inside_MEA(
        sv::MEAState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jl::MEALiquidFluxes{NB_GDL, NB_MPL},
        S_abs::MEASorptionSources,
        Sl::MEALiquidSources{NB_GDL, NB_MPL}
)::MEALiquidWaterDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # Surface-area reduction factors due to the ribs between the GC and GDL.
    Jl_agc_agdl_red = Jl.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    Jl_cgdl_cgc_red = Jl.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    H_gdl_node = pp.Hgdl / NB_GDL   # thickness of one GDL node
    H_mpl_node = pp.Hmpl / NB_MPL   # thickness of one MPL node

    # Anode GDL
    agdl_s = ntuple(NB_GDL) do j
        Jl_in = j == 1 ? Jl_agc_agdl_red : Jl.agdl_agdl[j - 1]
        Jl_out = j == NB_GDL ? Jl.agdl_ampl : Jl.agdl_agdl[j]
        1 / (rho_H2O_l(sv.agdl[j].T) * pp.epsilon_gdl) * ((Jl_in - Jl_out) / H_gdl_node + M_H2O * Sl.agdl[j])
    end

    # Anode MPL
    ampl_s = ntuple(NB_MPL) do j
        Jl_in = j == 1 ? Jl.agdl_ampl : Jl.ampl_ampl[j - 1]
        Jl_out = j == NB_MPL ? Jl.ampl_acl : Jl.ampl_ampl[j]
        1 / (rho_H2O_l(sv.ampl[j].T) * pp.epsilon_mpl) * ((Jl_in - Jl_out) / H_mpl_node + M_H2O * Sl.ampl[j])
    end

    # Anode and cathode CLs
    acl_s = 1 / (rho_H2O_l(sv.acl.T) * epsilon_cl(sv.acl.lambda, sv.acl.T, pp.Hacl)) *
            (Jl.ampl_acl / pp.Hacl - M_H2O * S_abs.l_acl + M_H2O * Sl.acl)
    ccl_s = 1 / (rho_H2O_l(sv.ccl.T) * epsilon_cl(sv.ccl.lambda, sv.ccl.T, pp.Hccl)) *
            (-Jl.ccl_cmpl / pp.Hccl - M_H2O * S_abs.l_ccl + M_H2O * Sl.ccl)

    # Cathode MPL
    cmpl_s = ntuple(NB_MPL) do j
        Jl_in = j == 1 ? Jl.ccl_cmpl : Jl.cmpl_cmpl[j - 1]
        Jl_out = j == NB_MPL ? Jl.cmpl_cgdl : Jl.cmpl_cmpl[j]
        1 / (rho_H2O_l(sv.cmpl[j].T) * pp.epsilon_mpl) * ((Jl_in - Jl_out) / H_mpl_node + M_H2O * Sl.cmpl[j])
    end

    # Cathode GDL
    cgdl_s = ntuple(NB_GDL) do j
        Jl_in = j == 1 ? Jl.cmpl_cgdl : Jl.cgdl_cgdl[j - 1]
        Jl_out = j == NB_GDL ? Jl_cgdl_cgc_red : Jl.cgdl_cgdl[j]
        1 / (rho_H2O_l(sv.cgdl[j].T) * pp.epsilon_gdl) * ((Jl_in - Jl_out) / H_gdl_node + M_H2O * Sl.cgdl[j])
    end

    return MEALiquidWaterDerivative{NB_GDL, NB_MPL}(agdl_s, ampl_s, acl_s, ccl_s, cmpl_s, cgdl_s)
end


"""Calculate the dynamic evolution of water vapour in the porous layers and CLs.

Parameters
----------
sv : MEAState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
Jv : MEAVaporFluxes{NB_GDL, NB_MPL}
    Water-vapour inter-layer fluxes (mol·m⁻²·s⁻¹).
Sv : MEAVaporSources{NB_GDL, NB_MPL}
    Vapour phase-change source terms (mol·m⁻³·s⁻¹).
S_abs : MEASorptionSources
    Water absorption rates at the CLs (mol·m⁻³·s⁻¹).

Returns
-------
MEACellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with C_v derivatives filled in.
"""
function calculate_dyn_vapor_evolution_inside_MEA(
        sv::MEAState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jv::MEAVaporFluxes{NB_GDL, NB_MPL},
        Sv::MEAVaporSources{NB_GDL, NB_MPL},
        S_abs::MEASorptionSources
)::MEAVaporDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    Jv_agc_agdl_red = Jv.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    Jv_cgdl_cgc_red = Jv.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    agdl_C_v = ntuple(NB_GDL) do j
        Jv_in = j == 1 ? Jv_agc_agdl_red : Jv.agdl_agdl[j - 1]
        Jv_out = j == NB_GDL ? Jv.agdl_ampl : Jv.agdl_agdl[j]
        1 / (pp.epsilon_gdl * (1 - sv.agdl[j].s)) * ((Jv_in - Jv_out) / H_gdl_node + Sv.agdl[j])
    end

    ampl_C_v = ntuple(NB_MPL) do j
        Jv_in = j == 1 ? Jv.agdl_ampl : Jv.ampl_ampl[j - 1]
        Jv_out = j == NB_MPL ? Jv.ampl_acl : Jv.ampl_ampl[j]
        1 / (pp.epsilon_mpl * (1 - sv.ampl[j].s)) * ((Jv_in - Jv_out) / H_mpl_node + Sv.ampl[j])
    end

    acl_C_v = 1 / (epsilon_cl(sv.acl.lambda, sv.acl.T, pp.Hacl) * (1 - sv.acl.s)) *
              (Jv.ampl_acl / pp.Hacl - S_abs.v_acl + Sv.acl)
    ccl_C_v = 1 / (epsilon_cl(sv.ccl.lambda, sv.ccl.T, pp.Hccl) * (1 - sv.ccl.s)) *
              (-Jv.ccl_cmpl / pp.Hccl - S_abs.v_ccl + Sv.ccl)

    cmpl_C_v = ntuple(NB_MPL) do j
        Jv_in = j == 1 ? Jv.ccl_cmpl : Jv.cmpl_cmpl[j - 1]
        Jv_out = j == NB_MPL ? Jv.cmpl_cgdl : Jv.cmpl_cmpl[j]
        1 / (pp.epsilon_mpl * (1 - sv.cmpl[j].s)) * ((Jv_in - Jv_out) / H_mpl_node + Sv.cmpl[j])
    end

    cgdl_C_v = ntuple(NB_GDL) do j
        Jv_in = j == 1 ? Jv.cmpl_cgdl : Jv.cgdl_cgdl[j - 1]
        Jv_out = j == NB_GDL ? Jv_cgdl_cgc_red : Jv.cgdl_cgdl[j]
        1 / (pp.epsilon_gdl * (1 - sv.cgdl[j].s)) * ((Jv_in - Jv_out) / H_gdl_node + Sv.cgdl[j])
    end

    return MEAVaporDerivative{NB_GDL, NB_MPL}(agdl_C_v, ampl_C_v, acl_C_v, ccl_C_v, cmpl_C_v, cgdl_C_v)
end


"""Calculate the dynamic evolution of H₂ (anode) and O₂ (cathode) in the porous layers.

Parameters
----------
sv : MEAState1D{NB_GDL, NB_MPL}
    Typed 1D internal state for one gas-channel column.
pp : PhysicalParams
    Fuel-cell physical parameters container (geometry, thicknesses and porous properties).
J_H2 : MEAHydrogenFluxes{NB_GDL, NB_MPL}
    Hydrogen inter-layer fluxes (mol·m⁻²·s⁻¹).
J_O2 : MEAOxygenFluxes{NB_GDL, NB_MPL}
    Oxygen inter-layer fluxes (mol·m⁻²·s⁻¹).
S_H2 : MEAGasReactionSources
    Hydrogen reaction / crossover source terms (mol·m⁻³·s⁻¹).
S_O2 : MEAGasReactionSources
    Oxygen reaction / crossover source terms (mol·m⁻³·s⁻¹).

Returns
-------
MEACellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with C_H2 and C_O2 derivatives filled in.
"""
function calculate_dyn_H2_O2_N2_evolution_inside_MEA(
        sv::MEAState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        J_H2::MEAHydrogenFluxes{NB_GDL, NB_MPL},
        J_O2::MEAOxygenFluxes{NB_GDL, NB_MPL},
        S_H2::MEAGasReactionSources,
        S_O2::MEAGasReactionSources
)::MEAH2O2SpeciesDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    J_H2_agc_agdl_red = J_H2.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_O2_cgdl_cgc_red = J_O2.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    agdl_C_H2 = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_H2_agc_agdl_red : J_H2.agdl_agdl[j - 1]
        J_out = j == NB_GDL ? J_H2.agdl_ampl : J_H2.agdl_agdl[j]
        1 / (pp.epsilon_gdl * (1 - sv.agdl[j].s)) * (J_in - J_out) / H_gdl_node
    end

    ampl_C_H2 = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_H2.agdl_ampl : J_H2.ampl_ampl[j - 1]
        J_out = j == NB_MPL ? J_H2.ampl_acl : J_H2.ampl_ampl[j]
        1 / (pp.epsilon_mpl * (1 - sv.ampl[j].s)) * (J_in - J_out) / H_mpl_node
    end

    acl_C_H2 = 1 / (epsilon_cl(sv.acl.lambda, sv.acl.T, pp.Hacl) * (1 - sv.acl.s)) *
               (J_H2.ampl_acl / pp.Hacl - S_H2.reac - S_H2.cros)

    ccl_C_O2 = 1 / (epsilon_cl(sv.ccl.lambda, sv.ccl.T, pp.Hccl) * (1 - sv.ccl.s)) *
               (-J_O2.ccl_cmpl / pp.Hccl - S_O2.reac - S_O2.cros)

    cmpl_C_O2 = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_O2.ccl_cmpl : J_O2.cmpl_cmpl[j - 1]
        J_out = j == NB_MPL ? J_O2.cmpl_cgdl : J_O2.cmpl_cmpl[j]
        1 / (pp.epsilon_mpl * (1 - sv.cmpl[j].s)) * (J_in - J_out) / H_mpl_node
    end

    cgdl_C_O2 = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_O2.cmpl_cgdl : J_O2.cgdl_cgdl[j - 1]
        J_out = j == NB_GDL ? J_O2_cgdl_cgc_red : J_O2.cgdl_cgdl[j]
        1 / (pp.epsilon_gdl * (1 - sv.cgdl[j].s)) * (J_in - J_out) / H_gdl_node
    end

    return MEAH2O2SpeciesDerivative{NB_GDL, NB_MPL}(agdl_C_H2, ampl_C_H2, acl_C_H2, ccl_C_O2, cmpl_C_O2, cgdl_C_O2)
end


"""Calculate the dynamic evolution of the cathode overpotential eta_c.

Parameters
----------
i_fc
    Fuel cell current density (A·m⁻²).
C_O2_Pt
    Oxygen concentration at the platinum surface (mol·m⁻³).
T_ccl : Float64
    Temperature in the cathode catalyst layer (K).
eta_c : Float64
    Cathode overpotential (V).
pp : PhysicalParams
    Fuel-cell physical parameters container.
i_n
    Crossover current density (A·m⁻²).

Returns
-------
MEACellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with eta_c derivative filled in (ccl).
"""
function calculate_dyn_voltage_evolution(
        i_fc,
        C_O2_Pt,
        T_ccl::Float64,
        eta_c::Float64,
        pp::PhysicalParams,
        i_n
)::MEAVoltageDerivative

    deta_c = 1 / (pp.C_scl * pp.Hccl) * ((i_fc + i_n) -
             pp.i0_c_ref * (C_O2_Pt / C_O2ref_red)^pp.kappa_c *
             exp(-Eact_O2_red / (R * T_ccl) * (1 / T_ccl - 1 / Tref_O2_red)) *
             exp(alpha_c * F / (R * T_ccl) * eta_c))

    return MEAVoltageDerivative(deta_c)
end


"""Calculate the dynamic evolution of temperature throughout the MEA.

Parameters
----------
pp : PhysicalParams
    Fuel-cell physical parameters container (layer thicknesses come from `pp`).
nb_gdl, nb_mpl : Int64
    Node counts (match NB_GDL and NB_MPL).
rho_Cp0 : MEAThermalIntermediates{NB_GDL, NB_MPL}
    Volumetric heat capacities at each node (J·m⁻³·K⁻¹).
Jt : MEAThermalFluxes{NB_GDL, NB_MPL}
    Conductive thermal fluxes through the MEA (W·m⁻²).
Q_r : MEAReactionHeat
    Electrochemical reaction heat sources (W·m⁻³).
Q_sorp : MEASorptionHeat
    Sorption heat sources (W·m⁻³).
Q_liq : MEALiquidHeat{NB_GDL, NB_MPL}
    Liquefaction / evaporation heat sources (W·m⁻³).
Q_p : MEAProtonHeat
    Ionic (protonic) Joule-heating sources (W·m⁻³).
Q_e : MEAElectricHeat{NB_GDL, NB_MPL}
    Electronic Joule-heating sources (W·m⁻³).

Returns
-------
MEACellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with T derivatives filled in for all MEA layers.
"""
function calculate_dyn_temperature_evolution_inside_MEA(
        rho_Cp0::MEAThermalIntermediates{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jt::MEAThermalFluxes{NB_GDL, NB_MPL},
        Q_r::MEAReactionHeat,
        Q_sorp::MEASorptionHeat,
        Q_liq::MEALiquidHeat{NB_GDL, NB_MPL},
        Q_p::MEAProtonHeat,
        Q_e::MEAElectricHeat{NB_GDL, NB_MPL}
)::MEATemperatureDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    agdl_T = ntuple(NB_GDL) do j
        Jt_in = j == 1 ? Jt.agc_agdl : Jt.agdl_agdl[j - 1]
        Jt_out = j == NB_GDL ? Jt.agdl_ampl : Jt.agdl_agdl[j]
        (1 / rho_Cp0.agdl[j]) * ((Jt_in - Jt_out) / H_gdl_node + Q_liq.agdl[j] + Q_e.agdl[j])
    end

    ampl_T = ntuple(NB_MPL) do j
        Jt_in = j == 1 ? Jt.agdl_ampl : Jt.ampl_ampl[j - 1]
        Jt_out = j == NB_MPL ? Jt.ampl_acl : Jt.ampl_ampl[j]
        (1 / rho_Cp0.ampl[j]) * ((Jt_in - Jt_out) / H_mpl_node + Q_liq.ampl[j] + Q_e.ampl[j])
    end

    acl_T = (1 / rho_Cp0.acl) * ((Jt.ampl_acl - Jt.acl_mem) / pp.Hacl +
            Q_r.acl + Q_sorp.v_acl + Q_sorp.l_acl + Q_liq.acl + Q_e.acl)
    mem_T = (1 / rho_Cp0.mem) * ((Jt.acl_mem - Jt.mem_ccl) / pp.Hmem + Q_p.mem)
    ccl_T = (1 / rho_Cp0.ccl) * ((Jt.mem_ccl - Jt.ccl_cmpl) / pp.Hccl +
            Q_r.ccl + Q_sorp.v_ccl + Q_sorp.l_ccl + Q_liq.ccl + Q_p.ccl + Q_e.ccl)

    cmpl_T = ntuple(NB_MPL) do j
        Jt_in = j == 1 ? Jt.ccl_cmpl : Jt.cmpl_cmpl[j - 1]
        Jt_out = j == NB_MPL ? Jt.cmpl_cgdl : Jt.cmpl_cmpl[j]
        (1 / rho_Cp0.cmpl[j]) * ((Jt_in - Jt_out) / H_mpl_node + Q_liq.cmpl[j] + Q_e.cmpl[j])
    end

    cgdl_T = ntuple(NB_GDL) do j
        Jt_in = j == 1 ? Jt.cmpl_cgdl : Jt.cgdl_cgdl[j - 1]
        Jt_out = j == NB_GDL ? Jt.cgdl_cgc : Jt.cgdl_cgdl[j]
        (1 / rho_Cp0.cgdl[j]) * ((Jt_in - Jt_out) / H_gdl_node + Q_liq.cgdl[j] + Q_e.cgdl[j])
    end

    return MEATemperatureDerivative{NB_GDL, NB_MPL}(agdl_T, ampl_T, acl_T, mem_T, ccl_T, cmpl_T, cgdl_T)
end
