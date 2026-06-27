# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""


# ____________________________________________________Main functions____________________________________________________

"""Calculate dissolved-water (lambda) dynamics contribution.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
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
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        S_abs::MEASorptionSources,
        J_lambda::MEADissolvedWaterFlux,
        Sp::MEAWaterProductionSources
)::MEADissolvedWaterDerivative where {NB_GDL, NB_MPL}

    # Extraction of the variables
    T_acl, T_ccl = sv.acl.T, sv.ccl.T
    lambda_acl, lambda_mem, lambda_ccl = sv.acl.lambda, sv.mem.lambda, sv.ccl.lambda

    # Differential equations
    d_lambda_acl_dt = M_eq / (rho_mem * epsilon_mc(lambda_acl, T_acl, pp.Hacl)) *
                  (-J_lambda.acl_mem / pp.Hacl + S_abs.v_acl + S_abs.l_acl + Sp.acl)
    d_lambda_mem_dt = M_eq / rho_mem * (J_lambda.acl_mem - J_lambda.mem_ccl) / pp.Hmem
    d_lambda_ccl_dt = M_eq / (rho_mem * epsilon_mc(lambda_ccl, T_ccl, pp.Hccl)) *
                  (J_lambda.mem_ccl / pp.Hccl + S_abs.v_ccl + S_abs.l_ccl + Sp.ccl)

    return MEADissolvedWaterDerivative(d_lambda_acl_dt, d_lambda_mem_dt, d_lambda_ccl_dt)
end


"""Calculate the dynamic evolution of liquid water in the porous layers.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
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
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with s (liquid saturation) derivatives filled in.
"""
function calculate_dyn_liquid_water_evolution_inside_MEA(
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jl::MEALiquidFluxes{NB_GDL, NB_MPL},
        S_abs::MEASorptionSources,
        Sl::MEALiquidSources{NB_GDL, NB_MPL}
)::MEALiquidWaterDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # Extraction of the variables
    T_agdl, T_ampl = getproperty.(sv.agdl, :T), getproperty.(sv.ampl, :T)
    T_cmpl, T_cgdl = getproperty.(sv.cmpl, :T), getproperty.(sv.cgdl, :T)
    T_acl, T_ccl = sv.acl.T, sv.ccl.T
    lambda_acl, lambda_ccl = sv.acl.lambda, sv.ccl.lambda

    # Surface-area reduction factors due to the ribs between the GC and GDL.
    Jl_agc_agdl_red = Jl.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    Jl_cgdl_cgc_red = Jl.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    H_gdl_node = pp.Hgdl / NB_GDL   # thickness of one GDL node
    H_mpl_node = pp.Hmpl / NB_MPL   # thickness of one MPL node

    # Differential equations
    #   Anode GDL
    d_s_agdl_dt = ntuple(NB_GDL) do j
        Jl_in = j == 1 ? Jl_agc_agdl_red : Jl.agdl_agdl[j - 1]
        Jl_out = j == NB_GDL ? Jl.agdl_ampl : Jl.agdl_agdl[j]
        1 / (rho_H2O_l(T_agdl[j]) * pp.epsilon_gdl) * ((Jl_in - Jl_out) / H_gdl_node + M_H2O * Sl.agdl[j])
    end

    #   Anode MPL
    d_s_ampl_dt = ntuple(NB_MPL) do j
        Jl_in = j == 1 ? Jl.agdl_ampl : Jl.ampl_ampl[j - 1]
        Jl_out = j == NB_MPL ? Jl.ampl_acl : Jl.ampl_ampl[j]
        1 / (rho_H2O_l(T_ampl[j]) * pp.epsilon_mpl) * ((Jl_in - Jl_out) / H_mpl_node + M_H2O * Sl.ampl[j])
    end

    #   Anode and cathode CLs
    d_s_acl_dt = 1 / (rho_H2O_l(T_acl) * epsilon_cl(lambda_acl, T_acl, pp.Hacl)) *
            (Jl.ampl_acl / pp.Hacl - M_H2O * S_abs.l_acl + M_H2O * Sl.acl)
    d_s_ccl_dt = 1 / (rho_H2O_l(T_ccl) * epsilon_cl(lambda_ccl, T_ccl, pp.Hccl)) *
            (-Jl.ccl_cmpl / pp.Hccl - M_H2O * S_abs.l_ccl + M_H2O * Sl.ccl)

    #   Cathode MPL
    d_s_cmpl_dt = ntuple(NB_MPL) do j
        Jl_in = j == 1 ? Jl.ccl_cmpl : Jl.cmpl_cmpl[j - 1]
        Jl_out = j == NB_MPL ? Jl.cmpl_cgdl : Jl.cmpl_cmpl[j]
        1 / (rho_H2O_l(T_cmpl[j]) * pp.epsilon_mpl) * ((Jl_in - Jl_out) / H_mpl_node + M_H2O * Sl.cmpl[j])
    end

    #   Cathode GDL
    d_s_cgdl_dt = ntuple(NB_GDL) do j
        Jl_in = j == 1 ? Jl.cmpl_cgdl : Jl.cgdl_cgdl[j - 1]
        Jl_out = j == NB_GDL ? Jl_cgdl_cgc_red : Jl.cgdl_cgdl[j]
        1 / (rho_H2O_l(T_cgdl[j]) * pp.epsilon_gdl) * ((Jl_in - Jl_out) / H_gdl_node + M_H2O * Sl.cgdl[j])
    end

    return MEALiquidWaterDerivative{NB_GDL, NB_MPL}(d_s_agdl_dt, d_s_ampl_dt, d_s_acl_dt, d_s_ccl_dt, d_s_cmpl_dt, d_s_cgdl_dt)
end


"""Calculate the dynamic evolution of water vapour in the porous layers and CLs.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
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
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with C_v derivatives filled in.
"""
function calculate_dyn_vapor_evolution_inside_MEA(
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        Jv::MEAVaporFluxes{NB_GDL, NB_MPL},
        Sv::MEAVaporSources{NB_GDL, NB_MPL},
        S_abs::MEASorptionSources
)::MEAVaporDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # Extraction of the variables
    s_agc, s_agdl = sv.agc.s, getproperty.(sv.agdl, :s)
    s_ampl, s_acl = getproperty.(sv.ampl, :s), sv.acl.s
    s_ccl, s_cmpl = sv.ccl.s, getproperty.(sv.cmpl, :s)
    s_cgdl, s_cgc = getproperty.(sv.cgdl, :s), sv.cgc.s

    T_agc, T_agdl = sv.agc.T, getproperty.(sv.agdl, :T)
    T_ampl, T_acl = getproperty.(sv.ampl, :T), sv.acl.T
    T_ccl, T_cmpl = sv.ccl.T, getproperty.(sv.cmpl, :T)
    T_cgdl, T_cgc = getproperty.(sv.cgdl, :T), sv.cgc.T

    lambda_acl, lambda_mem, lambda_ccl = sv.acl.lambda, sv.mem.lambda, sv.ccl.lambda

    # Surface-area reduction factors due to the ribs between the GC and GDL.
    Jv_agc_agdl_red = Jv.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    Jv_cgdl_cgc_red = Jv.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    # Differential equations
    #   Anode GDL
    d_C_v_agdl_dt = ntuple(NB_GDL) do j
        Jv_in = j == 1 ? Jv_agc_agdl_red : Jv.agdl_agdl[j - 1]
        Jv_out = j == NB_GDL ? Jv.agdl_ampl : Jv.agdl_agdl[j]
        1 / (pp.epsilon_gdl * (1 - s_agdl[j])) * ((Jv_in - Jv_out) / H_gdl_node + Sv.agdl[j])
    end

    #   Anode MPL
    d_C_v_ampl_dt = ntuple(NB_MPL) do j
        Jv_in = j == 1 ? Jv.agdl_ampl : Jv.ampl_ampl[j - 1]
        Jv_out = j == NB_MPL ? Jv.ampl_acl : Jv.ampl_ampl[j]
        1 / (pp.epsilon_mpl * (1 - s_ampl[j])) * ((Jv_in - Jv_out) / H_mpl_node + Sv.ampl[j])
    end

    #    Anode and cathode CLs
    d_C_v_acl_dt = 1 / (epsilon_cl(lambda_acl, T_acl, pp.Hacl) * (1 - s_acl)) *
              (Jv.ampl_acl / pp.Hacl - S_abs.v_acl + Sv.acl)
    d_C_v_ccl_dt = 1 / (epsilon_cl(lambda_ccl, T_ccl, pp.Hccl) * (1 - s_ccl)) *
              (-Jv.ccl_cmpl / pp.Hccl - S_abs.v_ccl + Sv.ccl)

    #   Cathode MPL
    d_C_v_cmpl_dt = ntuple(NB_MPL) do j
        Jv_in = j == 1 ? Jv.ccl_cmpl : Jv.cmpl_cmpl[j - 1]
        Jv_out = j == NB_MPL ? Jv.cmpl_cgdl : Jv.cmpl_cmpl[j]
        1 / (pp.epsilon_mpl * (1 - s_cmpl[j])) * ((Jv_in - Jv_out) / H_mpl_node + Sv.cmpl[j])
    end

    #   Cathode GDL
    d_C_v_cgdl_dt = ntuple(NB_GDL) do j
        Jv_in = j == 1 ? Jv.cmpl_cgdl : Jv.cgdl_cgdl[j - 1]
        Jv_out = j == NB_GDL ? Jv_cgdl_cgc_red : Jv.cgdl_cgdl[j]
        1 / (pp.epsilon_gdl * (1 - s_cgdl[j])) * ((Jv_in - Jv_out) / H_gdl_node + Sv.cgdl[j])
    end

    return MEAVaporDerivative{NB_GDL, NB_MPL}(d_C_v_agdl_dt, d_C_v_ampl_dt, d_C_v_acl_dt, d_C_v_ccl_dt, d_C_v_cmpl_dt, d_C_v_cgdl_dt)
end


"""Calculate the dynamic evolution of H₂ (anode) and O₂ (cathode) in the porous layers.

Parameters
----------
sv : CellState1D{NB_GDL, NB_MPL}
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
CellDerivative1D{NB_GDL, NB_MPL}
    Updated derivative container with C_H2 and C_O2 derivatives filled in.
"""
function calculate_dyn_H2_O2_N2_evolution_inside_MEA(
        sv::CellState1D{NB_GDL, NB_MPL},
        pp::PhysicalParams,
        J_H2::MEAHydrogenFluxes{NB_GDL, NB_MPL},
        J_O2::MEAOxygenFluxes{NB_GDL, NB_MPL},
        J_N2::MEANitrogenFluxes{NB_GDL, NB_MPL},
        S_H2::MEAGasReactionSources,
        S_O2::MEAGasReactionSources
)::MEAGasSpeciesDerivative{NB_GDL, NB_MPL} where {NB_GDL, NB_MPL}

    # Extraction of the variables
    s_agdl, s_ampl, s_acl = getproperty.(sv.agdl, :s), getproperty.(sv.ampl, :s), sv.acl.s
    s_ccl, s_cmpl, s_cgdl = sv.ccl.s, getproperty.(sv.cmpl, :s), getproperty.(sv.cgdl, :s)

    T_acl, T_ccl = sv.acl.T, sv.ccl.T
    lambda_acl, lambda_ccl = sv.acl.lambda, sv.ccl.lambda

    # Surface-area reduction factors due to the ribs between the GC and GDL.
    J_H2_agc_agdl_red = J_H2.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_O2_cgdl_cgc_red = J_O2.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_N2_agc_agdl_red = J_N2.agc_agdl * (pp.Wagc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)
    J_N2_cgdl_cgc_red = J_N2.cgdl_cgc * (pp.Wcgc * pp.Lgc) / (pp.Aact / pp.nb_channel_in_gc)

    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    # Differential equations
    #   Anode GDL
    d_C_H2_agdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_H2_agc_agdl_red : J_H2.agdl_agdl[j - 1]
        J_out = j == NB_GDL ? J_H2.agdl_ampl : J_H2.agdl_agdl[j]
        1 / (pp.epsilon_gdl * (1 - s_agdl[j])) * (J_in - J_out) / H_gdl_node
    end
    d_C_N2_agdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_N2_agc_agdl_red : J_N2.agdl_agdl[j - 1]
        J_out = j == NB_GDL ? J_N2.agdl_ampl : J_N2.agdl_agdl[j]
        1 / (pp.epsilon_gdl * (1 - s_agdl[j])) * (J_in - J_out) / H_gdl_node
    end

    #   Anode MPL
    d_C_H2_ampl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_H2.agdl_ampl : J_H2.ampl_ampl[j - 1]
        J_out = j == NB_MPL ? J_H2.ampl_acl : J_H2.ampl_ampl[j]
        1 / (pp.epsilon_mpl * (1 - s_ampl[j])) * (J_in - J_out) / H_mpl_node
    end
    d_C_N2_ampl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_N2.agdl_ampl : J_N2.ampl_ampl[j - 1]
        J_out = j == NB_MPL ? J_N2.ampl_acl : J_N2.ampl_ampl[j]
        1 / (pp.epsilon_mpl * (1 - s_ampl[j])) * (J_in - J_out) / H_mpl_node
    end

    #   Anode CL
    d_C_H2_acl_dt = 1 / (epsilon_cl(lambda_acl, T_acl, pp.Hacl) * (1 - s_acl)) *
               (J_H2.ampl_acl / pp.Hacl - S_H2.reac - S_H2.cros)
    d_C_N2_acl_dt = 1 / (epsilon_cl(lambda_acl, T_acl, pp.Hacl) * (1 - s_acl)) *
               (J_N2.ampl_acl / pp.Hacl)

    #   Cathode CL
    d_C_O2_ccl_dt = 1 / (epsilon_cl(lambda_ccl, T_ccl, pp.Hccl) * (1 - s_ccl)) *
               (-J_O2.ccl_cmpl / pp.Hccl - S_O2.reac - S_O2.cros)
    d_C_N2_ccl_dt = 1 / (epsilon_cl(lambda_ccl, T_ccl, pp.Hccl) * (1 - s_ccl)) *
               (-J_N2.ccl_cmpl / pp.Hccl)

    #   Cathode MPL
    d_C_O2_cmpl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_O2.ccl_cmpl : J_O2.cmpl_cmpl[j - 1]
        J_out = j == NB_MPL ? J_O2.cmpl_cgdl : J_O2.cmpl_cmpl[j]
        1 / (pp.epsilon_mpl * (1 - s_cmpl[j])) * (J_in - J_out) / H_mpl_node
    end
    d_C_N2_cmpl_dt = ntuple(NB_MPL) do j
        J_in = j == 1 ? J_N2.ccl_cmpl : J_N2.cmpl_cmpl[j - 1]
        J_out = j == NB_MPL ? J_N2.cmpl_cgdl : J_N2.cmpl_cmpl[j]
        1 / (pp.epsilon_mpl * (1 - s_cmpl[j])) * (J_in - J_out) / H_mpl_node
    end

    #   Cathode GDL
    d_C_O2_cgdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_O2.cmpl_cgdl : J_O2.cgdl_cgdl[j - 1]
        J_out = j == NB_GDL ? J_O2_cgdl_cgc_red : J_O2.cgdl_cgdl[j]
        1 / (pp.epsilon_gdl * (1 - s_cgdl[j])) * (J_in - J_out) / H_gdl_node
    end
    d_C_N2_cgdl_dt = ntuple(NB_GDL) do j
        J_in = j == 1 ? J_N2.cmpl_cgdl : J_N2.cgdl_cgdl[j - 1]
        J_out = j == NB_GDL ? J_N2_cgdl_cgc_red : J_N2.cgdl_cgdl[j]
        1 / (pp.epsilon_gdl * (1 - s_cgdl[j])) * (J_in - J_out) / H_gdl_node
    end

    return MEAGasSpeciesDerivative{NB_GDL, NB_MPL}(
        d_C_H2_agdl_dt, d_C_H2_ampl_dt, d_C_H2_acl_dt,
        d_C_N2_agdl_dt, d_C_N2_ampl_dt, d_C_N2_acl_dt,
        d_C_O2_ccl_dt, d_C_O2_cmpl_dt, d_C_O2_cgdl_dt,
        d_C_N2_ccl_dt, d_C_N2_cmpl_dt, d_C_N2_cgdl_dt
    )
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
CellDerivative1D{NB_GDL, NB_MPL}
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

    # During nonlinear/DAE iterations the algebraic unknown `C_O2_Pt` can briefly
    # step outside its physical domain (e.g. slightly negative). Protect the
    # Butler–Volmer concentration term which uses fractional powers.
    C_O2_Pt_safe = _positive_concentration_value(C_O2_Pt)

    # Differential equation
    d_eta_c_ccl_dt = 1 / (pp.C_scl * pp.Hccl) * ((i_fc + i_n) -
             pp.i0_c_ref * (C_O2_Pt_safe / C_O2ref_red)^pp.kappa_c *
             exp(-Eact_O2_red / (R * T_ccl) * (1 / T_ccl - 1 / Tref_O2_red)) *
             exp(alpha_c * F / (R * T_ccl) * eta_c))

    return MEAVoltageDerivative(d_eta_c_ccl_dt)
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
CellDerivative1D{NB_GDL, NB_MPL}
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

    # Extraction of the parameters
    H_gdl_node = pp.Hgdl / NB_GDL
    H_mpl_node = pp.Hmpl / NB_MPL

    # Differential equations
    #   Anode GDL
    d_T_agdl_dt = ntuple(NB_GDL) do j
        Jt_in = j == 1 ? Jt.agc_agdl : Jt.agdl_agdl[j - 1]
        Jt_out = j == NB_GDL ? Jt.agdl_ampl : Jt.agdl_agdl[j]
        (1 / rho_Cp0.agdl[j]) * ((Jt_in - Jt_out) / H_gdl_node + Q_liq.agdl[j] + Q_e.agdl[j])
    end

    #   Anode MPL
    d_T_ampl_dt = ntuple(NB_MPL) do j
        Jt_in = j == 1 ? Jt.agdl_ampl : Jt.ampl_ampl[j - 1]
        Jt_out = j == NB_MPL ? Jt.ampl_acl : Jt.ampl_ampl[j]
        (1 / rho_Cp0.ampl[j]) * ((Jt_in - Jt_out) / H_mpl_node + Q_liq.ampl[j] + Q_e.ampl[j])
    end

    #  Anode and cathode CLs + membrane
    d_T_acl_dt = (1 / rho_Cp0.acl) * ((Jt.ampl_acl - Jt.acl_mem) / pp.Hacl +
            Q_r.acl + Q_sorp.v_acl + Q_sorp.l_acl + Q_liq.acl + Q_e.acl)
    d_T_mem_dt = (1 / rho_Cp0.mem) * ((Jt.acl_mem - Jt.mem_ccl) / pp.Hmem + Q_p.mem)
    d_T_ccl_dt = (1 / rho_Cp0.ccl) * ((Jt.mem_ccl - Jt.ccl_cmpl) / pp.Hccl +
            Q_r.ccl + Q_sorp.v_ccl + Q_sorp.l_ccl + Q_liq.ccl + Q_p.ccl + Q_e.ccl)

    #   Cathode MPL
    d_T_cmpl_dt = ntuple(NB_MPL) do j
        Jt_in = j == 1 ? Jt.ccl_cmpl : Jt.cmpl_cmpl[j - 1]
        Jt_out = j == NB_MPL ? Jt.cmpl_cgdl : Jt.cmpl_cmpl[j]
        (1 / rho_Cp0.cmpl[j]) * ((Jt_in - Jt_out) / H_mpl_node + Q_liq.cmpl[j] + Q_e.cmpl[j])
    end

    #  Cathode GDL
    d_T_cgdl_dt = ntuple(NB_GDL) do j
        Jt_in = j == 1 ? Jt.cmpl_cgdl : Jt.cgdl_cgdl[j - 1]
        Jt_out = j == NB_GDL ? Jt.cgdl_cgc : Jt.cgdl_cgdl[j]
        (1 / rho_Cp0.cgdl[j]) * ((Jt_in - Jt_out) / H_gdl_node + Q_liq.cgdl[j] + Q_e.cgdl[j])
    end

    return MEATemperatureDerivative{NB_GDL, NB_MPL}(d_T_agdl_dt, d_T_ampl_dt, d_T_acl_dt, d_T_mem_dt, d_T_ccl_dt, d_T_cmpl_dt, d_T_cgdl_dt)
end
