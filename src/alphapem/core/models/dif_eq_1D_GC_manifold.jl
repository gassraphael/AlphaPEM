# This file represents all the differential equations used for the fuel cell model.

# _____________________________________________________Preliminaries____________________________________________________
# Importing constants' value and functions
include(joinpath(@__DIR__, "../../utils/physics_functions.jl"))


# ____________________________________________________Main functions____________________________________________________

"""Calculate dynamic gas evolution in the gas channels.

Parameters
----------
dif_eq : AbstractVector{Dict}
    Dictionaries used for saving the differential equations, one per GC node.
sv : AbstractVector{Dict}
    Solver variables, one dictionary per GC node.
Hagc : Float64
    Thickness of the anode gas channel (m).
Hcgc : Float64
    Thickness of the cathode gas channel (m).
Lgc : Float64
    Length of the gas channel (m).
nb_gc : Int64
    Number of gas channels.
type_auxiliary : String
    Type of auxiliary components used in the fuel cell system.
Jv : Dict
    Vapor flow between the different layers (mol.m-2.s-1).
J_H2 : Dict
    Hydrogen flow between the different layers (mol.m-2.s-1).
J_O2 : Dict
    Oxygen flow between the different layers (mol.m-2.s-1).
J_N2 : Dict
    Nitrogen flow between the different layers (mol.m-2.s-1).
"""
function calculate_dyn_gas_evolution_inside_gas_channel(
    dif_eq::AbstractVector{Dict},
    sv::AbstractVector{Dict},
    Hagc::Float64,
    Hcgc::Float64,
    Lgc::Float64,
    nb_gc::Int64,
    type_auxiliary::Symbol,
    Jv::Dict,
    J_H2::Dict,
    J_O2::Dict,
    J_N2::Dict
)

    # At the anode side, inside the AGC
    if nb_gc == 1
        dif_eq[1]["dC_v_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                     (Jv["agc_in"] - Jv["agc_out"]) / Lgc - Jv["agc_agdl"][1] / Hagc
    elseif nb_gc == 2
        dif_eq[1]["dC_v_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                     (Jv["agc_in"] - Jv["agc_agc"][1]) / (Lgc / nb_gc) - Jv["agc_agdl"][1] / Hagc
        dif_eq[2]["dC_v_agc / dt"] = 1 / (1 - sv[2]["s_agc"]) *
                                     (Jv["agc_agc"][1] - Jv["agc_out"]) / (Lgc / nb_gc) - Jv["agc_agdl"][2] / Hagc
    else # nb_gc > 2
        dif_eq[1]["dC_v_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                     (Jv["agc_in"] - Jv["agc_agc"][1]) / (Lgc / nb_gc) - Jv["agc_agdl"][1] / Hagc
        for i in 2:(nb_gc - 1)
            dif_eq[i]["dC_v_agc / dt"] = 1 / (1 - sv[i]["s_agc"]) *
                                         (Jv["agc_agc"][i - 1] - Jv["agc_agc"][i]) / (Lgc / nb_gc) - Jv["agc_agdl"][i] / Hagc
        end
        dif_eq[nb_gc]["dC_v_agc / dt"] = 1 / (1 - sv[nb_gc]["s_agc"]) *
                                         (Jv["agc_agc"][nb_gc - 1] - Jv["agc_out"]) / (Lgc / nb_gc) - Jv["agc_agdl"][nb_gc] / Hagc
    end

    if nb_gc == 1
        dif_eq[1]["dC_H2_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                      (J_H2["agc_in"] - J_H2["agc_out"]) / Lgc - J_H2["agc_agdl"][1] / Hagc
    elseif nb_gc == 2
        dif_eq[1]["dC_H2_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                      (J_H2["agc_in"] - J_H2["agc_agc"][1]) / (Lgc / nb_gc) - J_H2["agc_agdl"][1] / Hagc
        dif_eq[2]["dC_H2_agc / dt"] = 1 / (1 - sv[2]["s_agc"]) *
                                      (J_H2["agc_agc"][1] - J_H2["agc_out"]) / (Lgc / nb_gc) - J_H2["agc_agdl"][2] / Hagc
    else # nb_gc > 2
        dif_eq[1]["dC_H2_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                      (J_H2["agc_in"] - J_H2["agc_agc"][1]) / (Lgc / nb_gc) - J_H2["agc_agdl"][1] / Hagc
        for i in 2:(nb_gc - 1)
            dif_eq[i]["dC_H2_agc / dt"] = 1 / (1 - sv[i]["s_agc"]) *
                                          (J_H2["agc_agc"][i - 1] - J_H2["agc_agc"][i]) / (Lgc / nb_gc) - J_H2["agc_agdl"][i] / Hagc
        end
        dif_eq[nb_gc]["dC_H2_agc / dt"] = 1 / (1 - sv[nb_gc]["s_agc"]) *
                                          (J_H2["agc_agc"][nb_gc - 1] - J_H2["agc_out"]) / (Lgc / nb_gc) - J_H2["agc_agdl"][nb_gc] / Hagc
    end

    if type_auxiliary == :forced_convective_cathode_with_flow_through_anode  # Test bench: simulated H2 recirculation which leads to N2 in the anode.
        if nb_gc == 1
            dif_eq[1]["dC_N2_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                          (J_N2["agc_in"] - J_N2["agc_out"]) / Lgc
        elseif nb_gc == 2
            dif_eq[1]["dC_N2_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                          (J_N2["agc_in"] - J_N2["agc_agc"][1]) / (Lgc / nb_gc)
            dif_eq[2]["dC_N2_agc / dt"] = 1 / (1 - sv[2]["s_agc"]) *
                                          (J_N2["agc_agc"][1] - J_N2["agc_out"]) / (Lgc / nb_gc)
        else # nb_gc > 2
            dif_eq[1]["dC_N2_agc / dt"] = 1 / (1 - sv[1]["s_agc"]) *
                                          (J_N2["agc_in"] - J_N2["agc_agc"][1]) / (Lgc / nb_gc)
            for i in 2:(nb_gc - 1)
                dif_eq[i]["dC_N2_agc / dt"] = 1 / (1 - sv[i]["s_agc"]) *
                                              (J_N2["agc_agc"][i - 1] - J_N2["agc_agc"][i]) / (Lgc / nb_gc)
            end
            dif_eq[nb_gc]["dC_N2_agc / dt"] = 1 / (1 - sv[nb_gc]["s_agc"]) *
                                              (J_N2["agc_agc"][nb_gc - 1] - J_N2["agc_out"]) / (Lgc / nb_gc)
        end
    else
        for i in 1:nb_gc
            dif_eq[i]["dC_N2_agc / dt"] = 0
        end
    end

    # At the cathode side, inside the CGC
    if nb_gc == 1
        dif_eq[1]["dC_v_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                     (Jv["cgc_in"] - Jv["cgc_out"]) / Lgc + Jv["cgdl_cgc"][1] / Hcgc
    elseif nb_gc == 2
        dif_eq[1]["dC_v_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                     (Jv["cgc_in"] - Jv["cgc_cgc"][1]) / (Lgc / nb_gc) + Jv["cgdl_cgc"][1] / Hcgc
        dif_eq[2]["dC_v_cgc / dt"] = 1 / (1 - sv[2]["s_cgc"]) *
                                     (Jv["cgc_cgc"][1] - Jv["cgc_out"]) / (Lgc / nb_gc) + Jv["cgdl_cgc"][2] / Hcgc
    else # nb_gc > 2
        dif_eq[1]["dC_v_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                     (Jv["cgc_in"] - Jv["cgc_cgc"][1]) / (Lgc / nb_gc) + Jv["cgdl_cgc"][1] / Hcgc
        for i in 2:(nb_gc - 1)
            dif_eq[i]["dC_v_cgc / dt"] = 1 / (1 - sv[i]["s_cgc"]) *
                                         (Jv["cgc_cgc"][i - 1] - Jv["cgc_cgc"][i]) / (Lgc / nb_gc) + Jv["cgdl_cgc"][i] / Hcgc
        end
        dif_eq[nb_gc]["dC_v_cgc / dt"] = 1 / (1 - sv[nb_gc]["s_cgc"]) *
                                         (Jv["cgc_cgc"][nb_gc - 1] - Jv["cgc_out"]) / (Lgc / nb_gc) + Jv["cgdl_cgc"][nb_gc] / Hcgc
    end

    if nb_gc == 1
        dif_eq[1]["dC_O2_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                      (J_O2["cgc_in"] - J_O2["cgc_out"]) / Lgc + J_O2["cgdl_cgc"][1] / Hcgc
    elseif nb_gc == 2
        dif_eq[1]["dC_O2_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                      (J_O2["cgc_in"] - J_O2["cgc_cgc"][1]) / (Lgc / nb_gc) + J_O2["cgdl_cgc"][1] / Hcgc
        dif_eq[2]["dC_O2_cgc / dt"] = 1 / (1 - sv[2]["s_cgc"]) *
                                      (J_O2["cgc_cgc"][1] - J_O2["cgc_out"]) / (Lgc / nb_gc) + J_O2["cgdl_cgc"][2] / Hcgc
    else # nb_gc > 2
        dif_eq[1]["dC_O2_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                      (J_O2["cgc_in"] - J_O2["cgc_cgc"][1]) / (Lgc / nb_gc) + J_O2["cgdl_cgc"][1] / Hcgc
        for i in 2:(nb_gc - 1)
            dif_eq[i]["dC_O2_cgc / dt"] = 1 / (1 - sv[i]["s_cgc"]) *
                                          (J_O2["cgc_cgc"][i - 1] - J_O2["cgc_cgc"][i]) / (Lgc / nb_gc) + J_O2["cgdl_cgc"][i] / Hcgc
        end
        dif_eq[nb_gc]["dC_O2_cgc / dt"] = 1 / (1 - sv[nb_gc]["s_cgc"]) *
                                          (J_O2["cgc_cgc"][nb_gc - 1] - J_O2["cgc_out"]) / (Lgc / nb_gc) + J_O2["cgdl_cgc"][nb_gc] / Hcgc
    end

    if nb_gc == 1
        dif_eq[1]["dC_N2_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                      (J_N2["cgc_in"] - J_N2["cgc_out"]) / Lgc
    elseif nb_gc == 2
        dif_eq[1]["dC_N2_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                      (J_N2["cgc_in"] - J_N2["cgc_cgc"][1]) / (Lgc / nb_gc)
        dif_eq[2]["dC_N2_cgc / dt"] = 1 / (1 - sv[2]["s_cgc"]) *
                                      (J_N2["cgc_cgc"][1] - J_N2["cgc_out"]) / (Lgc / nb_gc)
    else # nb_gc > 2
        dif_eq[1]["dC_N2_cgc / dt"] = 1 / (1 - sv[1]["s_cgc"]) *
                                      (J_N2["cgc_in"] - J_N2["cgc_cgc"][1]) / (Lgc / nb_gc)
        for i in 2:(nb_gc - 1)
            dif_eq[i]["dC_N2_cgc / dt"] = 1 / (1 - sv[i]["s_cgc"]) *
                                          (J_N2["cgc_cgc"][i - 1] - J_N2["cgc_cgc"][i]) / (Lgc / nb_gc)
        end
        dif_eq[nb_gc]["dC_N2_cgc / dt"] = 1 / (1 - sv[nb_gc]["s_cgc"]) *
                                          (J_N2["cgc_cgc"][nb_gc - 1] - J_N2["cgc_out"]) / (Lgc / nb_gc)
    end
    return nothing
end


"""Calculate dynamic liquid-water evolution in the gas channels.

Parameters
----------
dif_eq : AbstractVector{Dict}
    Dictionaries used for saving the differential equations.
T_des : Float64
    Desired fuel cell temperature (K).
Hagc : Float64
    Thickness of the anode gas channel (m).
Hcgc : Float64
    Thickness of the cathode gas channel (m).
Lgc : Float64
    Length of the gas channel (m).
nb_gc : Int64
    Number of gas channels.
Jl : Dict
    Liquid water flow between the different layers (mol.m-2.s-1).
"""
function calculate_dyn_liq_evolution_inside_gas_channel(
    dif_eq::AbstractVector{Dict},
    T_des::Float64,
    Hagc::Float64,
    Hcgc::Float64,
    Lgc::Float64,
    nb_gc::Int64,
    Jl::Dict
)

    # At the anode side, inside the AGC
    if nb_gc == 1
        dif_eq[1]["ds_agc / dt"] = 1 / rho_H2O_l(T_des) * (- Jl["agc_out"] / Lgc - Jl["agc_agdl"][1] / Hagc)
    elseif nb_gc == 2
        dif_eq[1]["ds_agc / dt"] = 1 / rho_H2O_l(T_des) * (- Jl["agc_agc"][1] / (Lgc / nb_gc) - Jl["agc_agdl"][1] / Hagc)
        dif_eq[2]["ds_agc / dt"] = 1 / rho_H2O_l(T_des) * ((Jl["agc_agc"][1] - Jl["agc_out"]) / (Lgc / nb_gc) - Jl["agc_agdl"][2] / Hagc)
    else # nb_gc > 2
        dif_eq[1]["ds_agc / dt"] = 1 / rho_H2O_l(T_des) * (- Jl["agc_agc"][1] / (Lgc / nb_gc) - Jl["agc_agdl"][1] / Hagc)
        for i in 2:(nb_gc - 1)
            dif_eq[i]["ds_agc / dt"] = 1 / rho_H2O_l(T_des) * ((Jl["agc_agc"][i - 1] - Jl["agc_agc"][i]) / (Lgc / nb_gc) - Jl["agc_agdl"][i] / Hagc)
        end
        dif_eq[nb_gc]["ds_agc / dt"] = 1 / rho_H2O_l(T_des) * ((Jl["agc_agc"][nb_gc - 1] - Jl["agc_out"]) / (Lgc / nb_gc) - Jl["agc_agdl"][nb_gc] / Hagc)
    end

    # At the cathode side, inside the CGC
    if nb_gc == 1
        dif_eq[1]["ds_cgc / dt"] = 1 / rho_H2O_l(T_des) * (- Jl["cgc_out"] / Lgc + Jl["cgdl_cgc"][1] / Hcgc)
    elseif nb_gc == 2
        dif_eq[1]["ds_cgc / dt"] = 1 / rho_H2O_l(T_des) * (- Jl["cgc_cgc"][1] / (Lgc / nb_gc) + Jl["cgdl_cgc"][1] / Hcgc)
        dif_eq[2]["ds_cgc / dt"] = 1 / rho_H2O_l(T_des) * ((Jl["cgc_cgc"][1] - Jl["cgc_out"]) / (Lgc / nb_gc) + Jl["cgdl_cgc"][2] / Hcgc)
    else # nb_gc > 2
        dif_eq[1]["ds_cgc / dt"] = 1 / rho_H2O_l(T_des) * (- Jl["cgc_cgc"][1] / (Lgc / nb_gc) + Jl["cgdl_cgc"][1] / Hcgc)
        for i in 2:(nb_gc - 1)
            dif_eq[i]["ds_cgc / dt"] = 1 / rho_H2O_l(T_des) * ((Jl["cgc_cgc"][i - 1] - Jl["cgc_cgc"][i]) / (Lgc / nb_gc) + Jl["cgdl_cgc"][i] / Hcgc)
        end
        dif_eq[nb_gc]["ds_cgc / dt"] = 1 / rho_H2O_l(T_des) * ((Jl["cgc_cgc"][nb_gc - 1] - Jl["cgc_out"]) / (Lgc / nb_gc) + Jl["cgdl_cgc"][nb_gc] / Hcgc)
    end
    return nothing
end


"""Calculate dynamic temperature evolution in the gas channels.

Parameters
----------
dif_eq : AbstractVector{Dict}
    Dictionaries used for saving the differential equations.
nb_gc : Int64
    Number of gas channels.
"""
function calculate_dyn_temperature_evolution_inside_gas_channel(
    dif_eq::AbstractVector{Dict},
    nb_gc::Int64
)

    # At the anode side, inside the AGC
    for i in 1:nb_gc
        dif_eq[i]["dT_agc / dt"] = 0  # Dirichlet boundary condition. T_agc is initialized to T_fc and remains constant.
    end

    # At the cathode side, inside the CGC
    for i in 1:nb_gc
        dif_eq[i]["dT_cgc / dt"] = 0  # Dirichlet boundary condition. T_cgc is initialized to T_fc and remains constant.
    end
    return nothing
end


"""Calculate pressure and humidity dynamics in the manifolds.

Parameters
----------
dif_eq : Dict
    Dictionary used for saving the differential equations.
T_des : Float64
    Fuel cell temperature (K).
nb_cell : Int64
    Number of cells in the fuel cell stack.
Vasm : Float64
    Volume of the anode supply manifold (m³).
Vcsm : Float64
    Volume of the cathode supply manifold (m³).
Vaem : Float64
    Volume of the anode exhaust manifold (m³).
Vcem : Float64
    Volume of the cathode exhaust manifold (m³).
type_auxiliary : String
    Type of auxiliary components used in the fuel cell model.
W : Dict
    Matter flows between the different layers (mol.s-1).
Wv : Dict
    Vapor flows between the different layers (mol.s-1).
"""
function calculate_dyn_manifold_pressure_and_humidity_evolution(
    dif_eq::Dict,
    T_des::Float64,
    nb_cell::Int64,
    Vasm::Float64,
    Vcsm::Float64,
    Vaem::Float64,
    Vcem::Float64,
    type_auxiliary::Symbol,
    W::Dict,
    Wv::Dict
)

    # # Pressure evolution inside the manifolds
    # if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
    #    type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #     # At the anode side
    #     if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
    #         dif_eq["dPasm / dt"] = (W["a_in"] + W["asm_in_re_to_asm"] - nb_cell * Wasm_to_asm_out) / Vasm * R * T_des
    #         dif_eq["dPaem / dt"] = (nb_cell * Waem_in_to_aem - Waem_to_aem_out - Waem_to_aem_out_re) / Vaem * R * T_des
    #     else  # type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #         dif_eq["dPasm / dt"] = (W["a_in"] - nb_cell * Wasm_to_asm_out) / Vasm * R * T_des
    #         dif_eq["dPaem / dt"] = (nb_cell * Waem_in_to_aem - Waem_to_aem_out) / Vaem * R * T_des
    #     end
    #     # At the cathode side
    #     dif_eq["dPcsm / dt"] = (Wc_in - nb_cell * Wcsm_to_csm_out) / Vcsm * R * T_des
    #     dif_eq["dPcem / dt"] = (nb_cell * Wcem_in_to_cem - Wcem_to_cem_out) / Vcem * R * T_des
    # end
    #
    # # Humidity evolution inside the manifolds
    # if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
    #    type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #     # At the anode side
    #     if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
    #         dif_eq["dPhi_asm / dt"] = (Wv_asm_in_to_asm + Wv_asm_in_re_to_asm - nb_cell * Wv_asm_to_asm_out) / Vasm * R * T_des / Psat(T_des)
    #         dif_eq["dPhi_aem / dt"] = (nb_cell * Wv_aem_in_to_aem - Wv_aem_to_aem_out_re - Wv_aem_to_aem_out) / Vaem * R * T_des / Psat(T_des)
    #     else  # type_auxiliary == :forced_convective_cathode_with_flow_through_anode
    #         dif_eq["dPhi_asm / dt"] = (Wv_asm_in_to_asm - nb_cell * Wv_asm_to_asm_out) / Vasm * R * T_des / Psat(T_des)
    #         dif_eq["dPhi_aem / dt"] = (nb_cell * Wv_aem_in_to_aem - Wv_aem_to_aem_out) / Vaem * R * T_des / Psat(T_des)
    #     end
    #     # At the cathode side
    #     dif_eq["dPhi_csm / dt"] = (Wv_csm_in_to_csm - nb_cell * Wv_csm_to_csm_out) / Vcsm * R * T_des / Psat(T_des)
    #     dif_eq["dPhi_cem / dt"] = 
    # end
    return nothing
end

