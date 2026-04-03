# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model.
"""

# ____________________________________________________Main functions____________________________________________________

"""This function calculates the dynamic evolution of the dissolved water in the membrane and the catalyst layers.

Parameters
----------
dif_eq_1D : Dict
    Dictionary used for saving the differential equations.
sv_1D : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
Hmem : Float64
    Thickness of the membrane (m).
Hacl : Float64
    Thickness of the anode catalyst layer (m).
Hccl : Float64
    Thickness of the cathode catalyst layer (m).
S_abs : Dict
    Source terms of absorbed water in the anode and cathode catalyst layers (mol.m-3.s-1).
J_lambda : Dict
    Water fluxes between the different layers of the membrane (mol.m-2.s-1).
Sp : Dict
    Source terms of produced water in the anode and cathode catalyst layers (mol.m-3.s-1).
"""
function calculate_dyn_dissoved_water_evolution_inside_MEA(dif_eq_1D::Dict, sv_1D::Dict, Hmem::Float64,
                                                           Hacl::Float64, Hccl::Float64, S_abs::Dict,
                                                           J_lambda::Dict, Sp::Dict)
    dif_eq_1D["dlambda_acl / dt"] = M_eq / (rho_mem * epsilon_mc(sv_1D["lambda_acl"], sv_1D["T_acl"], Hacl)) *
                                     (-J_lambda["acl_mem"] / Hacl + S_abs["v_acl"] + S_abs["l_acl"] + Sp["acl"])
    dif_eq_1D["dlambda_mem / dt"] = M_eq / rho_mem * (J_lambda["acl_mem"] - J_lambda["mem_ccl"]) / Hmem
    dif_eq_1D["dlambda_ccl / dt"] = M_eq / (rho_mem * epsilon_mc(sv_1D["lambda_ccl"], sv_1D["T_ccl"], Hccl)) *
                                     (J_lambda["mem_ccl"] / Hccl + S_abs["v_ccl"] + S_abs["l_ccl"] + Sp["ccl"])

    return nothing
end


"""This function calculates the dynamic evolution of the liquid water in the gas diffusion and catalyst layers.

Parameters
----------
dif_eq_1D : Dict
    Dictionary used for saving the differential equations.
sv_1D : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
Aact : Float64
    Active area of one cell (m²).
Wagc : Float64
    Width of the anode gas channel (m).
Wcgc : Float64
    Width of the cathode gas channel (m).
Lgc : Float64
    Length of one channel of the gas channel (m).
nb_channel_in_gc : Int64
    Number of channels in the gas channel.
Hgdl : Float64
    Thickness of the gas diffusion layer (m).
Hmpl : Float64
    Thickness of the microporous layer (m).
Hacl : Float64
    Thickness of the anode catalyst layer (m).
Hccl : Float64
    Thickness of the cathode catalyst layer (m).
epsilon_gdl : Float64
    Anode/cathode GDL porosity.
epsilon_mpl : Float64
    Anode/cathode MPL porosity.
nb_gc : Int64
    Number of model nodes placed inside each GC.
nb_gdl : Int64
    Number of model nodes placed inside each GDL.
nb_mpl : Int64
    Number of model nodes placed inside each MPL.
Jl : Dict
    Liquid water flow between the different layers (mol.m-2.s-1).
S_abs : Dict
    Water absorption in the CLs (mol.m-3.s-1).
Sl : Dict
    Liquid water source terms inside the different layers (mol.m-3.s-1).
"""
function calculate_dyn_liquid_water_evolution_inside_MEA(dif_eq_1D::Dict, sv_1D::Dict, Aact::Float64, Wagc::Float64,
                                                         Wcgc::Float64, Lgc::Float64, nb_channel_in_gc::Int64,
                                                         Hgdl::Float64, Hmpl::Float64, Hacl::Float64, Hccl::Float64,
                                                         epsilon_gdl::Float64, epsilon_mpl::Float64, nb_gc::Int64,
    nb_gdl::Int64,
    nb_mpl::Int64,
    Jl::Dict,
    S_abs::Dict,
    Sl::Dict
)

    # At the anode side
    #       Inside the AGDL
    Jl_agc_agdl_red = Jl["agc_agdl"] * (Wagc * Lgc) / (Aact / nb_channel_in_gc)  # There is a surface reduction due to the presence of the ribs, which is taken into account here.
    if nb_gdl == 1
        dif_eq_1D["ds_agdl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_agdl_1"]) * epsilon_gdl) *
                                       ((Jl_agc_agdl_red - Jl["agdl_ampl"]) / Hgdl + M_H2O * Sl["agdl"][1])
    elseif nb_gdl == 2
        dif_eq_1D["ds_agdl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_agdl_1"]) * epsilon_gdl) *
                                       ((Jl_agc_agdl_red - Jl["agdl_agdl"][1]) / (Hgdl / nb_gdl) + M_H2O * Sl["agdl"][1])
        dif_eq_1D["ds_agdl_2 / dt"] = 1 / (rho_H2O_l(sv_1D["T_agdl_2"]) * epsilon_gdl) *
                                       ((Jl["agdl_agdl"][1] - Jl["agdl_ampl"]) / (Hgdl / nb_gdl) + M_H2O * Sl["agdl"][2])
    else  # n_gdl > 2
        dif_eq_1D["ds_agdl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_agdl_1"]) * epsilon_gdl) *
                                       ((Jl_agc_agdl_red - Jl["agdl_agdl"][1]) / (Hgdl / nb_gdl) + M_H2O * Sl["agdl"][1])
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["ds_agdl_$(i) / dt"] = 1 / (rho_H2O_l(sv_1D["T_agdl_$(i)"]) * epsilon_gdl) *
                                              ((Jl["agdl_agdl"][i - 1] - Jl["agdl_agdl"][i]) / (Hgdl / nb_gdl) +
                                               M_H2O * Sl["agdl"][i])
        end
        dif_eq_1D["ds_agdl_$(nb_gdl) / dt"] = 1 / (rho_H2O_l(sv_1D["T_agdl_$(nb_gdl)"]) * epsilon_gdl) *
                                               ((Jl["agdl_agdl"][nb_gdl - 1] - Jl["agdl_ampl"]) / (Hgdl / nb_gdl) +
                                                M_H2O * Sl["agdl"][nb_gdl])
    end

    #      Inside the AMPL
    if nb_mpl == 1
        dif_eq_1D["ds_ampl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_ampl_1"]) * epsilon_mpl) *
                                       ((Jl["agdl_ampl"] - Jl["ampl_acl"]) / Hmpl + M_H2O * Sl["ampl"][1])
    elseif nb_mpl == 2
        dif_eq_1D["ds_ampl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_ampl_1"]) * epsilon_mpl) *
                                       ((Jl["agdl_ampl"] - Jl["ampl_ampl"][1]) / (Hmpl / nb_mpl) + M_H2O * Sl["ampl"][1])
        dif_eq_1D["ds_ampl_2 / dt"] = 1 / (rho_H2O_l(sv_1D["T_ampl_2"]) * epsilon_mpl) *
                                       ((Jl["ampl_ampl"][1] - Jl["ampl_acl"]) / (Hmpl / nb_mpl) + M_H2O * Sl["ampl"][2])
    else  # n_mpl > 2
        dif_eq_1D["ds_ampl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_ampl_1"]) * epsilon_mpl) *
                                       ((Jl["agdl_ampl"] - Jl["ampl_ampl"][1]) / (Hmpl / nb_mpl) + M_H2O * Sl["ampl"][1])
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["ds_ampl_$(i) / dt"] = 1 / (rho_H2O_l(sv_1D["T_ampl_$(i)"]) * epsilon_mpl) *
                                              ((Jl["ampl_ampl"][i - 1] - Jl["ampl_ampl"][i]) / (Hmpl / nb_mpl) +
                                               M_H2O * Sl["ampl"][i])
        end
        dif_eq_1D["ds_ampl_$(nb_mpl) / dt"] = 1 / (rho_H2O_l(sv_1D["T_ampl_$(nb_mpl)"]) * epsilon_mpl) *
                                               ((Jl["ampl_ampl"][nb_mpl - 1] - Jl["ampl_acl"]) / (Hmpl / nb_mpl) +
                                                M_H2O * Sl["ampl"][nb_mpl])
    end

    #      Inside the ACL
    dif_eq_1D["ds_acl / dt"] = 1 / (rho_H2O_l(sv_1D["T_acl"]) * epsilon_cl(sv_1D["lambda_acl"], sv_1D["T_acl"], Hacl)) *
                                (Jl["ampl_acl"] / Hacl - M_H2O * S_abs["l_acl"] + M_H2O * Sl["acl"])

    # At the cathode side
    #       Inside the CCL
    dif_eq_1D["ds_ccl / dt"] = 1 / (rho_H2O_l(sv_1D["T_ccl"]) * epsilon_cl(sv_1D["lambda_ccl"], sv_1D["T_ccl"], Hccl)) *
                                (-Jl["ccl_cmpl"] / Hccl - M_H2O * S_abs["l_ccl"] + M_H2O * Sl["ccl"])

    #       Inside the CMPL
    if nb_mpl == 1
        dif_eq_1D["ds_cmpl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cmpl_1"]) * epsilon_mpl) *
                                       ((Jl["ccl_cmpl"] - Jl["cmpl_cgdl"]) / Hmpl + M_H2O * Sl["cmpl"][1])
    elseif nb_mpl == 2
        dif_eq_1D["ds_cmpl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cmpl_1"]) * epsilon_mpl) *
                                       ((Jl["ccl_cmpl"] - Jl["cmpl_cmpl"][1]) / (Hmpl / nb_mpl) + M_H2O * Sl["cmpl"][1])
        dif_eq_1D["ds_cmpl_2 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cmpl_2"]) * epsilon_mpl) *
                                       ((Jl["cmpl_cmpl"][1] - Jl["cmpl_cgdl"]) / (Hmpl / nb_mpl) + M_H2O * Sl["cmpl"][2])
    else  # n_mpl > 2
        dif_eq_1D["ds_cmpl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cmpl_1"]) * epsilon_mpl) *
                                       ((Jl["ccl_cmpl"] - Jl["cmpl_cmpl"][1]) / (Hmpl / nb_mpl) + M_H2O * Sl["cmpl"][1])
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["ds_cmpl_$(i) / dt"] = 1 / (rho_H2O_l(sv_1D["T_cmpl_$(i)"]) * epsilon_mpl) *
                                              ((Jl["cmpl_cmpl"][i - 1] - Jl["cmpl_cmpl"][i]) / (Hmpl / nb_mpl) +
                                               M_H2O * Sl["cmpl"][i])
        end
        dif_eq_1D["ds_cmpl_$(nb_mpl) / dt"] = 1 / (rho_H2O_l(sv_1D["T_cmpl_$(nb_mpl)"]) * epsilon_mpl) *
                                               ((Jl["cmpl_cmpl"][nb_mpl - 1] - Jl["cmpl_cgdl"]) / (Hmpl / nb_mpl) +
                                                M_H2O * Sl["cmpl"][nb_mpl])
    end

    #       Inside the CGDL
    Jl_cgdl_cgc_red = Jl["cgdl_cgc"] * (Wcgc * Lgc) / (Aact / nb_channel_in_gc)  # There is a surface reduction due to the presence of the ribs, which is taken into account here.
    if nb_gdl == 1
        dif_eq_1D["ds_cgdl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cgdl_1"]) * epsilon_gdl) *
                                       ((Jl["cmpl_cgdl"] - Jl_cgdl_cgc_red) / Hgdl + M_H2O * Sl["cgdl"][1])
    elseif nb_gdl == 2
        dif_eq_1D["ds_cgdl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cgdl_1"]) * epsilon_gdl) *
                                       ((Jl["cmpl_cgdl"] - Jl["cgdl_cgdl"][1]) / (Hgdl / nb_gdl) + M_H2O * Sl["cgdl"][1])
        dif_eq_1D["ds_cgdl_2 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cgdl_2"]) * epsilon_gdl) *
                                       ((Jl["cgdl_cgdl"][1] - Jl_cgdl_cgc_red) / (Hgdl / nb_gdl) + M_H2O * Sl["cgdl"][2])
    else
        dif_eq_1D["ds_cgdl_1 / dt"] = 1 / (rho_H2O_l(sv_1D["T_cgdl_1"]) * epsilon_gdl) *
                                       ((Jl["cmpl_cgdl"] - Jl["cgdl_cgdl"][1]) / (Hgdl / nb_gdl) + M_H2O * Sl["cgdl"][1])
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["ds_cgdl_$(i) / dt"] = 1 / (rho_H2O_l(sv_1D["T_cgdl_$(i)"]) * epsilon_gdl) *
                                              ((Jl["cgdl_cgdl"][i - 1] - Jl["cgdl_cgdl"][i]) / (Hgdl / nb_gdl) + M_H2O * Sl["cgdl"][i])
        end
        dif_eq_1D["ds_cgdl_$(nb_gdl) / dt"] = 1 / (rho_H2O_l(sv_1D["T_cgdl_$(nb_gdl)"]) * epsilon_gdl) *
                                               ((Jl["cgdl_cgdl"][nb_gdl - 1] - Jl_cgdl_cgc_red) / (Hgdl / nb_gdl) + M_H2O * Sl["cgdl"][nb_gdl])
    end

    return nothing
end


"""This function calculates the dynamic evolution of the vapor in the gas diffusion layers, the microporous layers,
and the catalyst layers.

Parameters
----------
dif_eq_1D : Dict
    Dictionary used for saving the differential equations.
sv_1D : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
Aact : Float64
    Active area of one cell (m²).
Wagc : Float64
    Width of the anode gas channel (m).
Wcgc : Float64
    Width of the cathode gas channel (m).
Lgc : Float64
    Length of one channel of the gas channel (m).
nb_channel_in_gc : Int64
    Number of channels in the gas channel.
Hgdl : Float64
    Thickness of the gas diffusion layer (m).
Hmpl : Float64
    Thickness of the microporous layer (m).
Hacl : Float64
    Thickness of the anode catalyst layer (m).
Hccl : Float64
    Thickness of the cathode catalyst layer (m).
epsilon_gdl : Float64
    Anode/cathode GDL porosity.
epsilon_mpl : Float64
    Anode/cathode MPL porosity.
nb_gc : Int64
    Number of model nodes placed inside each GC.
nb_gdl : Int64
    Number of model nodes placed inside each GDL.
nb_mpl : Int64
    Number of model nodes placed inside each MPL.
Jv : Dict
    Vapor flow between the different layers (mol.m-2.s-1).
Sv : Dict
    Vapor source terms inside the different layers (mol.m-3.s-1).
S_abs : Dict
    Water absorption in the CLs (mol.m-3.s-1).
"""
function calculate_dyn_vapor_evolution_inside_MEA(
    dif_eq_1D::Dict,
    sv_1D::Dict,
    Aact::Float64,
    Wagc::Float64,
    Wcgc::Float64,
    Lgc::Float64,
    nb_channel_in_gc::Int64,
    Hgdl::Float64,
    Hmpl::Float64,
    Hacl::Float64,
    Hccl::Float64,
    epsilon_gdl::Float64,
    epsilon_mpl::Float64,
    nb_gc::Int64,
    nb_gdl::Int64,
    nb_mpl::Int64,
    Jv::Dict,
    Sv::Dict,
    S_abs::Dict
)

    # At the anode side
    #       Inside the AGDL
    Jv_agc_agdl_red = Jv["agc_agdl"] * (Wagc * Lgc) / (Aact / nb_channel_in_gc)  # There is a surface reduction due to the presence of the ribs, which is taken into account here.
    if nb_gdl == 1
        dif_eq_1D["dC_v_agdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_1"])) *
                                         ((Jv_agc_agdl_red - Jv["agdl_ampl"]) / Hgdl + Sv["agdl"][1])
    elseif nb_gdl == 2
        dif_eq_1D["dC_v_agdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_1"])) *
                                         ((Jv_agc_agdl_red - Jv["agdl_agdl"][1]) / (Hgdl / nb_gdl) + Sv["agdl"][1])
        dif_eq_1D["dC_v_agdl_2 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_2"])) *
                                         ((Jv["agdl_agdl"][1] - Jv["agdl_ampl"]) / (Hgdl / nb_gdl) + Sv["agdl"][2])
    else  # n_gdl > 2
        dif_eq_1D["dC_v_agdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_1"])) *
                                         ((Jv_agc_agdl_red - Jv["agdl_agdl"][1]) / (Hgdl / nb_gdl) + Sv["agdl"][1])
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["dC_v_agdl_$(i) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_$(i)"])) *
                                                ((Jv["agdl_agdl"][i - 1] - Jv["agdl_agdl"][i]) / (Hgdl / nb_gdl) + Sv["agdl"][i])
        end
        dif_eq_1D["dC_v_agdl_$(nb_gdl) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_$(nb_gdl)"])) *
                                                 ((Jv["agdl_agdl"][nb_gdl - 1] - Jv["agdl_ampl"]) / (Hgdl / nb_gdl) + Sv["agdl"][nb_gdl])
    end

    #       Inside the AMPL
    if nb_mpl == 1
        dif_eq_1D["dC_v_ampl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_1"])) *
                                         ((Jv["agdl_ampl"] - Jv["ampl_acl"]) / Hmpl + Sv["ampl"][1])
    elseif nb_mpl == 2
        dif_eq_1D["dC_v_ampl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_1"])) *
                                         ((Jv["agdl_ampl"] - Jv["ampl_ampl"][1]) / (Hmpl / nb_mpl) + Sv["ampl"][1])
        dif_eq_1D["dC_v_ampl_2 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_2"])) *
                                         ((Jv["ampl_ampl"][1] - Jv["ampl_acl"]) / (Hmpl / nb_mpl) + Sv["ampl"][2])
    else  # n_mpl > 2
        dif_eq_1D["dC_v_ampl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_1"])) *
                                         ((Jv["agdl_ampl"] - Jv["ampl_ampl"][1]) / (Hmpl / nb_mpl) + Sv["ampl"][1])
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["dC_v_ampl_$(i) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_$(i)"])) *
                                                ((Jv["ampl_ampl"][i - 1] - Jv["ampl_ampl"][i]) / (Hmpl / nb_mpl) + Sv["ampl"][i])
        end
        dif_eq_1D["dC_v_ampl_$(nb_mpl) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_$(nb_mpl)"])) *
                                                 ((Jv["ampl_ampl"][nb_mpl - 1] - Jv["ampl_acl"]) / (Hmpl / nb_mpl) + Sv["ampl"][nb_mpl])
    end

    #       Inside the ACL
    dif_eq_1D["dC_v_acl / dt"] = 1 / (epsilon_cl(sv_1D["lambda_acl"], sv_1D["T_acl"], Hacl) * (1 - sv_1D["s_acl"])) *
                                  (Jv["ampl_acl"] / Hacl - S_abs["v_acl"] + Sv["acl"])

    # At the cathode side
    #       Inside the CCL
    dif_eq_1D["dC_v_ccl / dt"] = 1 / (epsilon_cl(sv_1D["lambda_ccl"], sv_1D["T_ccl"], Hccl) * (1 - sv_1D["s_ccl"])) *
                                  (-Jv["ccl_cmpl"] / Hccl - S_abs["v_ccl"] + Sv["ccl"])

    #       Inside the CMPL
    if nb_mpl == 1
        dif_eq_1D["dC_v_cmpl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_1"])) *
                                         ((Jv["ccl_cmpl"] - Jv["cmpl_cgdl"]) / Hmpl + Sv["cmpl"][1])
    elseif nb_mpl == 2
        dif_eq_1D["dC_v_cmpl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_1"])) *
                                         ((Jv["ccl_cmpl"] - Jv["cmpl_cmpl"][1]) / (Hmpl / nb_mpl) + Sv["cmpl"][1])
        dif_eq_1D["dC_v_cmpl_2 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_2"])) *
                                         ((Jv["cmpl_cmpl"][1] - Jv["cmpl_cgdl"]) / (Hmpl / nb_mpl) + Sv["cmpl"][2])
    else  # n_mpl > 2
        dif_eq_1D["dC_v_cmpl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_1"])) *
                                         ((Jv["ccl_cmpl"] - Jv["cmpl_cmpl"][1]) / (Hmpl / nb_mpl) + Sv["cmpl"][1])
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["dC_v_cmpl_$(i) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_$(i)"])) *
                                                ((Jv["cmpl_cmpl"][i - 1] - Jv["cmpl_cmpl"][i]) / (Hmpl / nb_mpl) + Sv["cmpl"][i])
        end
        dif_eq_1D["dC_v_cmpl_$(nb_mpl) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_$(nb_mpl)"])) *
                                                 ((Jv["cmpl_cmpl"][nb_mpl - 1] - Jv["cmpl_cgdl"]) / (Hmpl / nb_mpl) + Sv["cmpl"][nb_mpl])
    end

    #       Inside the CGDL
    Jv_cgdl_cgc_red = Jv["cgdl_cgc"] * (Wcgc * Lgc) / (Aact / nb_channel_in_gc)  # There is a surface reduction due to the presence of the ribs, which is taken into account here.
    if nb_gdl == 1
        dif_eq_1D["dC_v_cgdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_1"])) *
                                         ((Jv["cmpl_cgdl"] - Jv_cgdl_cgc_red) / Hgdl + Sv["cgdl"][1])
    elseif nb_gdl == 2
        dif_eq_1D["dC_v_cgdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_1"])) *
                                         ((Jv["cmpl_cgdl"] - Jv["cgdl_cgdl"][1]) / (Hgdl / nb_gdl) + Sv["cgdl"][1])
        dif_eq_1D["dC_v_cgdl_2 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_2"])) *
                                         ((Jv["cgdl_cgdl"][1] - Jv_cgdl_cgc_red) / (Hgdl / nb_gdl) + Sv["cgdl"][2])
    else  # n_gdl > 2
        dif_eq_1D["dC_v_cgdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_1"])) *
                                         ((Jv["cmpl_cgdl"] - Jv["cgdl_cgdl"][1]) / (Hgdl / nb_gdl) + Sv["cgdl"][1])
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["dC_v_cgdl_$(i) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_$(i)"])) *
                                                ((Jv["cgdl_cgdl"][i - 1] - Jv["cgdl_cgdl"][i]) / (Hgdl / nb_gdl) + Sv["cgdl"][i])
        end
        dif_eq_1D["dC_v_cgdl_$(nb_gdl) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_$(nb_gdl)"])) *
                                                 ((Jv["cgdl_cgdl"][nb_gdl - 1] - Jv_cgdl_cgc_red) / (Hgdl / nb_gdl) + Sv["cgdl"][nb_gdl])
    end

    return nothing
end


"""This function calculates the dynamic evolution of the hydrogen, oxygen and nitrogen in the gas diffusion layers,
the microporous layers, and the catalyst layers.

Parameters
----------
dif_eq_1D : Dict
    Dictionary used for saving the differential equations.
sv_1D : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
Aact : Float64
    Active area of one cell (m²).
Wagc : Float64
    Width of the anode gas channel (m).
Wcgc : Float64
    Width of the cathode gas channel (m).
Lgc : Float64
    Length of one channel of the gas channel (m).
nb_channel_in_gc : Int64
    Number of channels in the gas channel.
Hgdl : Float64
    Thickness of the gas diffusion layer (m).
Hmpl : Float64
    Thickness of the microporous layer (m).
Hacl : Float64
    Thickness of the anode catalyst layer (m).
Hccl : Float64
    Thickness of the cathode catalyst layer (m).
epsilon_gdl : Float64
    Anode/cathode GDL porosity.
epsilon_mpl : Float64
    Anode/cathode MPL porosity.
nb_gdl : Int64
    Number of model nodes placed inside each GDL.
nb_mpl : Int64
    Number of model nodes placed inside each MPL.
nb_gc : Int64
    Number of model nodes placed inside each GC.
J_H2 : Dict
    Hydrogen flow between the different layers (mol.m-2.s-1).
J_O2 : Dict
    Oxygen flow between the different layers (mol.m-2.s-1).
S_H2 : Dict
    Hydrogen source terms inside the different layers (mol.m-3.s-1).
S_O2 : Dict
    Oxygen source terms inside the different layers (mol.m-3.s-1).
"""
function calculate_dyn_H2_O2_N2_evolution_inside_MEA(
    dif_eq_1D::Dict,
    sv_1D::Dict,
    Aact::Float64,
    Wagc::Float64,
    Wcgc::Float64,
    Lgc::Float64,
    nb_channel_in_gc::Int64,
    Hgdl::Float64,
    Hmpl::Float64,
    Hacl::Float64,
    Hccl::Float64,
    epsilon_gdl::Float64,
    epsilon_mpl::Float64,
    nb_gdl::Int64,
    nb_mpl::Int64,
    nb_gc::Int64,
    J_H2::Dict,
    J_O2::Dict,
    S_H2::Dict,
    S_O2::Dict
)

    # At the anode side
    #      Inside the AGDL
    J_H2_agc_agdl_red = J_H2["agc_agdl"] * (Wagc * Lgc) / (Aact / nb_channel_in_gc)  # There is a surface reduction due to the presence of the ribs, which is taken into account here.
    if nb_gdl == 1
        dif_eq_1D["dC_H2_agdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_1"])) *
                                          (J_H2_agc_agdl_red - J_H2["agdl_ampl"]) / Hgdl
    elseif nb_gdl == 2
        dif_eq_1D["dC_H2_agdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_1"])) *
                                          (J_H2_agc_agdl_red - J_H2["agdl_agdl"][1]) / (Hgdl / nb_gdl)
        dif_eq_1D["dC_H2_agdl_2 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_2"])) *
                                          (J_H2["agdl_agdl"][1] - J_H2["agdl_ampl"]) / (Hgdl / nb_gdl)
    else  # n_gdl > 2
        dif_eq_1D["dC_H2_agdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_1"])) *
                                          (J_H2_agc_agdl_red - J_H2["agdl_agdl"][1]) / (Hgdl / nb_gdl)
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["dC_H2_agdl_$(i) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_$(i)"])) *
                                                 (J_H2["agdl_agdl"][i - 1] - J_H2["agdl_agdl"][i]) / (Hgdl / nb_gdl)
        end
        dif_eq_1D["dC_H2_agdl_$(nb_gdl) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_agdl_$(nb_gdl)"])) *
                                                  (J_H2["agdl_agdl"][nb_gdl - 1] - J_H2["agdl_ampl"]) / (Hgdl / nb_gdl)
    end

    #      Inside the AMPL
    if nb_mpl == 1
        dif_eq_1D["dC_H2_ampl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_1"])) * (J_H2["agdl_ampl"] - J_H2["ampl_acl"]) / Hmpl
    elseif nb_mpl == 2
        dif_eq_1D["dC_H2_ampl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_1"])) *
                                          (J_H2["agdl_ampl"] - J_H2["ampl_ampl"][1]) / (Hmpl / nb_mpl)
        dif_eq_1D["dC_H2_ampl_2 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_2"])) *
                                          (J_H2["ampl_ampl"][1] - J_H2["ampl_acl"]) / (Hmpl / nb_mpl)
    else  # n_mpl > 2
        dif_eq_1D["dC_H2_ampl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_1"])) *
                                          (J_H2["agdl_ampl"] - J_H2["ampl_ampl"][1]) / (Hmpl / nb_mpl)
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["dC_H2_ampl_$(i) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_$(i)"])) *
                                                 (J_H2["ampl_ampl"][i - 1] - J_H2["ampl_ampl"][i]) / (Hmpl / nb_mpl)
        end
        dif_eq_1D["dC_H2_ampl_$(nb_mpl) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_ampl_$(nb_mpl)"])) *
                                                  (J_H2["ampl_ampl"][nb_mpl - 1] - J_H2["ampl_acl"]) / (Hmpl / nb_mpl)
    end

    #      Inside the ACL
    dif_eq_1D["dC_H2_acl / dt"] = 1 / (epsilon_cl(sv_1D["lambda_acl"], sv_1D["T_acl"], Hacl) * (1 - sv_1D["s_acl"])) *
                                   (J_H2["ampl_acl"] / Hacl - S_H2["reac"] - S_H2["cros"])

    # At the cathode side
    #      Inside the CCL
    dif_eq_1D["dC_O2_ccl / dt"] = 1 / (epsilon_cl(sv_1D["lambda_ccl"], sv_1D["T_ccl"], Hccl) * (1 - sv_1D["s_ccl"])) *
                                   (-J_O2["ccl_cmpl"] / Hccl - S_O2["reac"] - S_O2["cros"])

    #      Inside the CMPL
    if nb_mpl == 1
        dif_eq_1D["dC_O2_cmpl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_1"])) * (J_O2["ccl_cmpl"] - J_O2["cmpl_cgdl"]) / Hmpl
    elseif nb_mpl == 2
        dif_eq_1D["dC_O2_cmpl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_1"])) *
                                          (J_O2["ccl_cmpl"] - J_O2["cmpl_cmpl"][1]) / (Hmpl / nb_mpl)
        dif_eq_1D["dC_O2_cmpl_2 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_2"])) *
                                          (J_O2["cmpl_cmpl"][1] - J_O2["cmpl_cgdl"]) / (Hmpl / nb_mpl)
    else  # n_mpl > 2
        dif_eq_1D["dC_O2_cmpl_1 / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_1"])) *
                                          (J_O2["ccl_cmpl"] - J_O2["cmpl_cmpl"][1]) / (Hmpl / nb_mpl)
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["dC_O2_cmpl_$(i) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_$(i)"])) *
                                                 (J_O2["cmpl_cmpl"][i - 1] - J_O2["cmpl_cmpl"][i]) / (Hmpl / nb_mpl)
        end
        dif_eq_1D["dC_O2_cmpl_$(nb_mpl) / dt"] = 1 / (epsilon_mpl * (1 - sv_1D["s_cmpl_$(nb_mpl)"])) *
                                                  (J_O2["cmpl_cmpl"][nb_mpl - 1] - J_O2["cmpl_cgdl"]) / (Hmpl / nb_mpl)
    end

    #      Inside the CGDL
    J_O2_cgdl_cgc_red = J_O2["cgdl_cgc"] * (Wcgc * Lgc) / (Aact / nb_channel_in_gc)  # There is a surface reduction due to the presence of the ribs, which is taken into account here.
    if nb_gdl == 1
        dif_eq_1D["dC_O2_cgdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_1"])) * (J_O2["cmpl_cgdl"] - J_O2_cgdl_cgc_red) / Hgdl
    elseif nb_gdl == 2
        dif_eq_1D["dC_O2_cgdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_1"])) *
                                          (J_O2["cmpl_cgdl"] - J_O2["cgdl_cgdl"][1]) / (Hgdl / nb_gdl)
        dif_eq_1D["dC_O2_cgdl_2 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_2"])) *
                                          (J_O2["cgdl_cgdl"][1] - J_O2_cgdl_cgc_red) / (Hgdl / nb_gdl)
    else
        dif_eq_1D["dC_O2_cgdl_1 / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_1"])) *
                                          (J_O2["cmpl_cgdl"] - J_O2["cgdl_cgdl"][1]) / (Hgdl / nb_gdl)
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["dC_O2_cgdl_$(i) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_$(i)"])) *
                                                 (J_O2["cgdl_cgdl"][i - 1] - J_O2["cgdl_cgdl"][i]) / (Hgdl / nb_gdl)
        end
        dif_eq_1D["dC_O2_cgdl_$(nb_gdl) / dt"] = 1 / (epsilon_gdl * (1 - sv_1D["s_cgdl_$(nb_gdl)"])) *
                                                  (J_O2["cgdl_cgdl"][nb_gdl - 1] - J_O2_cgdl_cgc_red) / (Hgdl / nb_gdl)
    end

    return nothing
end


"""This function calculates the dynamic evolution of the temperature in the fuel cell.

Parameters
----------
dif_eq_1D : Dict
    Dictionary used for saving the differential equations.
Hgdl : Float64
    Thickness of the gas diffusion layer, in m.
Hmpl : Float64
    Thickness of the microporous layer, in m.
Hacl : Float64
    Thickness of the anode catalyst layer, in m.
Hccl : Float64
    Thickness of the cathode catalyst layer, in m.
Hmem : Float64
    Thickness of the membrane, in m.
nb_gdl : Int64
    Number of model nodes placed inside each GDL.
nb_mpl : Int64
    Number of model nodes placed inside each MPL.
rho_Cp0 : Dict
    Volumetric heat capacity of the different components of the fuel cell system, in J.m-3.K-1.
Jt : Dict
    Heat flows occuring inside the fuel cell system, J.m-2.s-1.
Q_r : Dict
    Heat dissipated by the electrochemical reaction 2*H2 + O2 -> 2*H2O, in J.m-3.s-1.
Q_sorp : Dict
    Heat dissipated by the absorption of water from the CL to the membrane, in J.m-3.s-1.
Q_liq : Dict
    Heat dissipated by the evaporation of liquid water, in J.m-3.s-1.
Q_p : Dict
    Heat dissipated by the ionic currents (Joule heating + Ohm's law), in J.m-3.s-1.
Q_e : Dict
    Heat dissipated by the electric currents (Joule heating + Ohm's law), in J.m-3.s-1.
"""
function calculate_dyn_temperature_evolution_inside_MEA(
    dif_eq_1D::Dict,
    Hgdl::Float64,
    Hmpl::Float64,
    Hacl::Float64,
    Hccl::Float64,
    Hmem::Float64,
    nb_gdl::Int64,
    nb_mpl::Int64,
    rho_Cp0::Dict,
    Jt::Dict,
    Q_r::Dict,
    Q_sorp::Dict,
    Q_liq::Dict,
    Q_p::Dict,
    Q_e::Dict
)

    # At the anode side
    #       Inside the AGDL
    if nb_gdl == 1
        dif_eq_1D["dT_agdl_1 / dt"] = (1 / rho_Cp0["agdl_1"]) * ((Jt["agc_agdl"] - Jt["agdl_ampl"]) / Hgdl +
                                                                  Q_liq["agdl_1"] + Q_e["agdl_1"])
    elseif nb_gdl == 2
        dif_eq_1D["dT_agdl_1 / dt"] = (1 / rho_Cp0["agdl_1"]) * ((Jt["agc_agdl"] - Jt["agdl_agdl_1"]) / (Hgdl / nb_gdl) +
                                                                  Q_liq["agdl_1"] + Q_e["agdl_1"])
        dif_eq_1D["dT_agdl_2 / dt"] = (1 / rho_Cp0["agdl_2"]) * ((Jt["agdl_agdl_1"] - Jt["agdl_ampl"]) / (Hgdl / nb_gdl) +
                                                                  Q_liq["agdl_2"] + Q_e["agdl_2"])
    else  # n_gdl > 2
        dif_eq_1D["dT_agdl_1 / dt"] = (1 / rho_Cp0["agdl_1"]) * ((Jt["agc_agdl"] - Jt["agdl_agdl_1"]) / (Hgdl / nb_gdl) +
                                                                  Q_liq["agdl_1"] + Q_e["agdl_1"])
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["dT_agdl_$(i) / dt"] = (1 / rho_Cp0["agdl_$(i)"]) *
                                              ((Jt["agdl_agdl_$(i - 1)"] - Jt["agdl_agdl_$(i)"]) / (Hgdl / nb_gdl) +
                                               Q_liq["agdl_$(i)"] + Q_e["agdl_$(i)"])
        end
        dif_eq_1D["dT_agdl_$(nb_gdl) / dt"] = (1 / rho_Cp0["agdl_$(nb_gdl)"]) *
                                               ((Jt["agdl_agdl_$(nb_gdl - 1)"] - Jt["agdl_ampl"]) / (Hgdl / nb_gdl) +
                                                Q_liq["agdl_$(nb_gdl)"] + Q_e["agdl_$(nb_gdl)"])
    end

    #      Inside the AMPL
    if nb_mpl == 1
        dif_eq_1D["dT_ampl_1 / dt"] = (1 / rho_Cp0["ampl_1"]) * ((Jt["agdl_ampl"] - Jt["ampl_acl"]) / Hmpl +
                                                                  Q_liq["ampl_1"] + Q_e["ampl_1"])
    elseif nb_mpl == 2
        dif_eq_1D["dT_ampl_1 / dt"] = (1 / rho_Cp0["ampl_1"]) * ((Jt["agdl_ampl"] - Jt["ampl_ampl_1"]) / (Hmpl / nb_mpl) +
                                                                  Q_liq["ampl_1"] + Q_e["ampl_1"])
        dif_eq_1D["dT_ampl_2 / dt"] = (1 / rho_Cp0["ampl_2"]) * ((Jt["ampl_ampl_1"] - Jt["ampl_acl"]) / (Hmpl / nb_mpl) +
                                                                  Q_liq["ampl_2"] + Q_e["ampl_2"])
    else  # n_mpl > 2
        dif_eq_1D["dT_ampl_1 / dt"] = (1 / rho_Cp0["ampl_1"]) * ((Jt["agdl_ampl"] - Jt["ampl_ampl_1"]) / (Hmpl / nb_mpl) +
                                                                  Q_liq["ampl_1"] + Q_e["ampl_1"])
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["dT_ampl_$(i) / dt"] = (1 / rho_Cp0["ampl_$(i)"]) *
                                              ((Jt["ampl_ampl_$(i - 1)"] - Jt["ampl_ampl_$(i)"]) / (Hmpl / nb_mpl) +
                                               Q_liq["ampl_$(i)"] + Q_e["ampl_$(i)"])
        end
        dif_eq_1D["dT_ampl_$(nb_mpl) / dt"] = (1 / rho_Cp0["ampl_$(nb_mpl)"]) *
                                               ((Jt["ampl_ampl_$(nb_mpl - 1)"] - Jt["ampl_acl"]) / (Hmpl / nb_mpl) +
                                                Q_liq["ampl_$(nb_mpl)"] + Q_e["ampl_$(nb_mpl)"])
    end

    #      Inside the ACL
    dif_eq_1D["dT_acl / dt"] = (1 / rho_Cp0["acl"]) *
                                ((Jt["ampl_acl"] - Jt["acl_mem"]) / Hacl +
                                 Q_r["acl"] + Q_sorp["v_acl"] + Q_sorp["l_acl"] + Q_liq["acl"] + Q_e["acl"])

    # Inside the membrane
    dif_eq_1D["dT_mem / dt"] = (1 / rho_Cp0["mem"]) * ((Jt["acl_mem"] - Jt["mem_ccl"]) / Hmem + Q_p["mem"])

    # At the cathode side
    #       Inside the CCL
    dif_eq_1D["dT_ccl / dt"] = (1 / rho_Cp0["ccl"]) *
                                ((Jt["mem_ccl"] - Jt["ccl_cmpl"]) / Hccl +
                                 Q_r["ccl"] + Q_sorp["v_ccl"] + Q_sorp["l_ccl"] + Q_liq["ccl"] + Q_p["ccl"] + Q_e["ccl"])

    #      Inside the CMPL
    if nb_mpl == 1
        dif_eq_1D["dT_cmpl_1 / dt"] = (1 / rho_Cp0["cmpl_1"]) * ((Jt["ccl_cmpl"] - Jt["cmpl_cgdl"]) / Hmpl +
                                                                  Q_liq["cmpl_1"] + Q_e["cmpl_1"])
    elseif nb_mpl == 2
        dif_eq_1D["dT_cmpl_1 / dt"] = (1 / rho_Cp0["cmpl_1"]) *
                                       ((Jt["ccl_cmpl"] - Jt["cmpl_cmpl_1"]) / (Hmpl / nb_mpl) +
                                        Q_liq["cmpl_1"] + Q_e["cmpl_1"])
        dif_eq_1D["dT_cmpl_2 / dt"] = (1 / rho_Cp0["cmpl_2"]) *
                                       ((Jt["cmpl_cmpl_1"] - Jt["cmpl_cgdl"]) / (Hmpl / nb_mpl) +
                                        Q_liq["cmpl_2"] + Q_e["cmpl_2"])
    else  # n_mpl > 2
        dif_eq_1D["dT_cmpl_1 / dt"] = (1 / rho_Cp0["cmpl_1"]) *
                                       ((Jt["ccl_cmpl"] - Jt["cmpl_cmpl_1"]) / (Hmpl / nb_mpl) +
                                        Q_liq["cmpl_1"] + Q_e["cmpl_1"])
        for i in 2:(nb_mpl - 1)
            dif_eq_1D["dT_cmpl_$(i) / dt"] = (1 / rho_Cp0["cmpl_$(i)"]) *
                                              ((Jt["cmpl_cmpl_$(i - 1)"] - Jt["cmpl_cmpl_$(i)"]) / (Hmpl / nb_mpl) +
                                               Q_liq["cmpl_$(i)"] + Q_e["cmpl_$(i)"])
        end
        dif_eq_1D["dT_cmpl_$(nb_mpl) / dt"] = (1 / rho_Cp0["cmpl_$(nb_mpl)"]) *
                                               ((Jt["cmpl_cmpl_$(nb_mpl - 1)"] - Jt["cmpl_cgdl"]) / (Hmpl / nb_mpl) +
                                                Q_liq["cmpl_$(nb_mpl)"] + Q_e["cmpl_$(nb_mpl)"])
    end

    #       Inside the CGDL
    if nb_gdl == 1
        dif_eq_1D["dT_cgdl_1 / dt"] = (1 / rho_Cp0["cgdl_1"]) * ((Jt["cmpl_cgdl"] - Jt["cgdl_cgc"]) / Hgdl +
                                                                  Q_liq["cgdl_1"] + Q_e["cgdl_1"])
    elseif nb_gdl == 2
        dif_eq_1D["dT_cgdl_1 / dt"] = (1 / rho_Cp0["cgdl_1"]) *
                                       ((Jt["cmpl_cgdl"] - Jt["cgdl_cgdl_1"]) / (Hgdl / nb_gdl) +
                                        Q_liq["cgdl_1"] + Q_e["cgdl_1"])
        dif_eq_1D["dT_cgdl_2 / dt"] = (1 / rho_Cp0["cgdl_2"]) *
                                       ((Jt["cgdl_cgdl_1"] - Jt["cgdl_cgc"]) / (Hgdl / nb_gdl) +
                                        Q_liq["cgdl_2"] + Q_e["cgdl_2"])
    else  # n_gdl > 2
        dif_eq_1D["dT_cgdl_1 / dt"] = (1 / rho_Cp0["cgdl_1"]) *
                                       ((Jt["cmpl_cgdl"] - Jt["cgdl_cgdl_1"]) / (Hgdl / nb_gdl) +
                                        Q_liq["cgdl_1"] + Q_e["cgdl_1"])
        for i in 2:(nb_gdl - 1)
            dif_eq_1D["dT_cgdl_$(i) / dt"] = (1 / rho_Cp0["cgdl_$(i)"]) *
                                              ((Jt["cgdl_cgdl_$(i - 1)"] - Jt["cgdl_cgdl_$(i)"]) / (Hgdl / nb_gdl) +
                                               Q_liq["cgdl_$(i)"] + Q_e["cgdl_$(i)"])
        end
        dif_eq_1D["dT_cgdl_$(nb_gdl) / dt"] = (1 / rho_Cp0["cgdl_$(nb_gdl)"]) *
                                               ((Jt["cgdl_cgdl_$(nb_gdl - 1)"] - Jt["cgdl_cgc"]) / (Hgdl / nb_gdl) +
                                                Q_liq["cgdl_$(nb_gdl)"] + Q_e["cgdl_$(nb_gdl)"])
    end

    return nothing
end


"""This function calculates the dynamic evolution of the cell overpotential eta_c.

Parameters
----------
dif_eq_1D : Dict
    Dictionary used for saving the differential equations.
i_fc
    Fuel cell current density (A.m-2).
C_O2_Pt
    Oxygen concentration at the platinum surface (mol.m-3).
T_ccl : Float64
    Fuel cell temperature in the cathode catalyst layer (K).
eta_c : Float64
    Cell overpotential (V).
Hccl : Float64
    Thickness of the cathode catalyst layer (m).
i0_c_ref : Float64
    Reference exchange current density at the cathode (A.m-2).
kappa_c : Float64
    Overpotential correction exponent.
C_scl : Float64
    Volumetric space-charge layer capacitance (F.m-3).
i_n
    Crossover current density (A.m-2).
"""
function calculate_dyn_voltage_evolution(
    dif_eq_1D::Dict,
    i_fc,
    C_O2_Pt,
    T_ccl::Float64,
    eta_c::Float64,
    Hccl::Float64,
    i0_c_ref::Float64,
    kappa_c::Float64,
    C_scl::Float64,
    i_n
)
    dif_eq_1D["deta_c / dt"] = 1 / (C_scl * Hccl) * ((i_fc + i_n) - i0_c_ref * (C_O2_Pt / C_O2ref_red)^kappa_c *
                                                      exp(-Eact_O2_red / (R * T_ccl) * (1 / T_ccl - 1 / Tref_O2_red)) *
                                                      exp(alpha_c * F / (R * T_ccl) * eta_c))

    # Expression used when theta_Pt is considered as a dynamic variable
    # theta_Pt_0 = 0  # This is the initial platine-oxide coverage, assumed to be zero for simplification.
    # omega = 3.0e3 # J.mol-1, It is the energy parameter for the Temkin isotherm [Hao 2015]
    #
    # dif_eq["deta_c / dt"] = 1 / (C_scl * Hccl) * ((i_fc + i_n) - i0_d_c_ref * ECSA_0 * (1 - theta_Pt_0) * (C_O2_Pt / C_O2ref_red)^kappa_c *
    #                                                exp(-Eact_O2_red / (R * T_ccl) * (1 / T_ccl - 1 / Tref_O2_red)) *
    #                                                exp(alpha_c * F / (R * T_ccl) * eta_c + omega * theta_Pt_0 / (R * T_ccl)))  # sign to be checked for omega.

    return nothing
end

