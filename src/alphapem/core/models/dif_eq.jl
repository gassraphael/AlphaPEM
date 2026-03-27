# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model. It is a master file which converts
the 1D+1D+1D model to several 1D models in order to ease the coding.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
include(joinpath(@__DIR__, "../../utils/physics_constants.jl"))
include(joinpath(@__DIR__, "velocity.jl"))
include(joinpath(@__DIR__, "flows_1D_MEA.jl"))
include(joinpath(@__DIR__, "flows_1D_GC_manifold.jl"))
include(joinpath(@__DIR__, "current_distribution_1D_GC.jl"))
include(joinpath(@__DIR__, "heat_transfer.jl"))
include(joinpath(@__DIR__, "dif_eq_1D_MEA.jl"))
include(joinpath(@__DIR__, "dif_eq_1D_GC_manifold.jl"))
include(joinpath(@__DIR__, "dif_eq_auxiliaries.jl"))
include(joinpath(@__DIR__, "../modules/cell_voltage_modules.jl"))
include(joinpath(@__DIR__, "../modules/dif_eq_modules.jl"))


# ______________________Objective function to solve. It gives the system of differential equations______________________

"""This function gives the system of differential equations to solve.

Parameters
----------
t :
    Time (s).
y : Vector
    Vector of the solver variables.
operating_inputs : Dict
    Operating inputs of the fuel cell.
parameters : Dict
    Parameters of the fuel cell model.
solver_variable_names : Vector{Vector{String}}
    Names of the solver variables. The first element contains the MEA and GC variable names,
    the second element contains the manifold variable names, and the third element contains
    the auxiliary variable names.

Returns
-------
Vector
    Vector containing the derivative of the solver variables.
"""
function dydt(t,
              y::Vector{<:Number},
              operating_inputs::Dict,
              parameters::Dict,
              solver_variable_names::Vector{Vector{String}})::Vector

    # Extraction of frequently used parameters
    T_des, Pa_des, Pc_des = operating_inputs["T_des"], operating_inputs["Pa_des"], operating_inputs["Pc_des"]
    Hacl, Hccl, Hmem = parameters["Hacl"], parameters["Hccl"], parameters["Hmem"]
    Hgdl, Hmpl = parameters["Hgdl"], parameters["Hmpl"]
    Aact, Wagc, Wcgc, Lgc = parameters["Aact"], parameters["Wagc"], parameters["Wcgc"], parameters["Lgc"]
    Hagc, Hcgc = parameters["Hagc"], parameters["Hcgc"]
    Vasm, Vcsm, Vaem, Vcem = parameters["Vasm"], parameters["Vcsm"], parameters["Vaem"], parameters["Vcem"]
    A_T_a, A_T_c = parameters["A_T_a"], parameters["A_T_c"]
    nb_cell, nb_channel_in_gc = parameters["nb_cell"], parameters["nb_channel_in_gc"]
    epsilon_gdl, epsilon_mpl = parameters["epsilon_gdl"], parameters["epsilon_mpl"]
    i0_c_ref, kappa_c, C_scl = parameters["i0_c_ref"], parameters["kappa_c"], parameters["C_scl"]
    nb_gc, nb_gdl, nb_mpl = parameters["nb_gc"], parameters["nb_gdl"], parameters["nb_mpl"]
    type_auxiliary = parameters["type_auxiliary"]

    # Creation of the sv (solver variables) and dif_eq dictionaries. They are intermediates to simplify the code.
    sv_1D_cell = [Dict() for _ in 1:nb_gc] # List of dicts: index i -> solver variables at GC node i
    sv_1D_manifold = Dict()                # Solver variables inside the gas channel and manifold 1D line
    sv_auxiliary = Dict()                  # Auxiliary solver variables
    dif_eq_1D_cell = [Dict() for _ in 1:nb_gc] # List of dicts: index i -> differential equations at GC node i
    dif_eq_1D_manifold = Dict()            # Differential equations inside the gas channel 1D line
    dif_eq_auxiliary = Dict()              # Differential equations of the auxiliary variables
    for i in 1:nb_gc  # Each dictionary in this loop corresponds to one gas channel node
        for (index, variable) in enumerate(solver_variable_names[1])  # Concern only the MEA and GC variables
            sv_1D_cell[i][variable] = y[index + (i - 1) * length(solver_variable_names[1])]
        end
    end
    if parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode" ||
       parameters["type_auxiliary"] == "forced-convective_cathode_with_anodic_recirculation"
        for (index, variable) in enumerate(solver_variable_names[2])  # Concern only the manifold variables
            sv_1D_manifold[variable] = y[index + nb_gc * length(solver_variable_names[1])]
        end
        for (index, variable) in enumerate(solver_variable_names[3])  # Concern only the auxiliary variables
            sv_auxiliary[variable] = y[index + nb_gc * length(solver_variable_names[1]) + length(solver_variable_names[2])]
        end
    end

    # Conditions to pursue the calculations
    for i in 1:nb_gc
        if sv_1D_cell[i]["eta_c"] > E0
            throw(ArgumentError("The cathode overpotential is higher than the open circuit voltage at time t = " *
                                string(t) * " s. It means that the voltage is negative, which is not possible."))
        end
    end

    # Intermediate values (one dict per GC node, naturally 1-based)
    dif_eq_int_values = [calculate_dif_eq_int_values(t, sv_1D_cell[i], operating_inputs, parameters)
                         for i in 1:nb_gc]

    # Calculate the local current density at each node of the GC.
    i_fc_cell = operating_inputs["current_density"](t, parameters)
    i_fc = calculate_1D_GC_current_density(i_fc_cell, sv_1D_cell, parameters)

    # Calculation of the oxygen concentration at the platinum surface in the cathode catalyst layer
    C_O2_Pt = [calculate_C_O2_Pt(i_fc[i], sv_1D_cell[i], parameters) for i in 1:nb_gc]

    # Calculation of the velocities inside the GC and the manifolds
    v_a, v_c, Pa_in, Pc_in = calculate_velocity_evolution(sv_1D_cell, i_fc_cell, operating_inputs, parameters)

    # Calculation of the flows. These elements are vectors of dictionaries. Each index corresponds to one GC node.
    # Each element of the dictionary corresponds to one flow inside this GC node.
    flows_1D_MEA = [calculate_flows_1D_MEA(sv_1D_cell[i], i_fc[i], v_a[i], v_c[i], operating_inputs, parameters)
                    for i in 1:nb_gc]
    flows_1D_GC_manifold = calculate_flows_1D_GC_manifold(sv_1D_cell, sv_1D_manifold, sv_auxiliary, i_fc_cell,
                                                           v_a, v_c, Pa_in, Pc_in, operating_inputs, parameters)
    for (dif_eq_key, flow_names) in [("agc_agdl", ("Jl", "Jv", "J_H2")), ("cgdl_cgc", ("Jl", "Jv", "J_O2"))]
        for flow_name in flow_names
            flows_1D_GC_manifold[flow_name][dif_eq_key] = [flows_1D_MEA[i][flow_name][dif_eq_key]
                                                            for i in 1:nb_gc]
        end
    end
    heat_flows_global = [calculate_heat_transfers(sv_1D_cell[i], i_fc[i], operating_inputs, parameters,
                                                   flows_1D_MEA[i]["S_abs"], flows_1D_MEA[i]["Sl"])
                         for i in 1:nb_gc]

    # Calculation of the dynamic evolutions
    for i in 1:nb_gc
        sv_i = sv_1D_cell[i]
        dif_eq_i = dif_eq_1D_cell[i]
        flows_i = flows_1D_MEA[i]
        heat_i = heat_flows_global[i]
        dif_eq_int_values_i = dif_eq_int_values[i]

    #       Inside the MEA
        calculate_dyn_dissoved_water_evolution_inside_MEA(dif_eq_i, sv_i, Hmem, Hacl, Hccl,
                                                          flows_i["S_abs"], flows_i["J_lambda"], flows_i["Sp"])
        calculate_dyn_liquid_water_evolution_inside_MEA(dif_eq_i, sv_i,
                                                         Aact, Wagc, Wcgc, Lgc, nb_channel_in_gc,
                                                         Hgdl, Hmpl, Hacl, Hccl,
                                                         epsilon_gdl, epsilon_mpl,
                                                         nb_gc, nb_gdl, nb_mpl,
                                                         flows_i["Jl"], flows_i["S_abs"], flows_i["Sl"])
        calculate_dyn_vapor_evolution_inside_MEA(dif_eq_i, sv_i,
                                                  Aact, Wagc, Wcgc, Lgc, nb_channel_in_gc,
                                                  Hgdl, Hmpl, Hacl, Hccl,
                                                  epsilon_gdl, epsilon_mpl,
                                                  nb_gc, nb_gdl, nb_mpl,
                                                  flows_i["Jv"], flows_i["Sv"], flows_i["S_abs"])
        calculate_dyn_H2_O2_N2_evolution_inside_MEA(dif_eq_i, sv_i,
                                                      Aact, Wagc, Wcgc, Lgc, nb_channel_in_gc,
                                                      Hgdl, Hmpl, Hacl, Hccl,
                                                      epsilon_gdl, epsilon_mpl,
                                                      nb_gdl, nb_mpl, nb_gc,
                                                      flows_i["J_H2"], flows_i["J_O2"],
                                                      flows_i["S_H2"], flows_i["S_O2"])
        calculate_dyn_voltage_evolution(dif_eq_i, i_fc[i], C_O2_Pt[i],             # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                         sv_i["T_ccl"], sv_i["eta_c"],
                                         Hccl, i0_c_ref, kappa_c, C_scl,
                                         dif_eq_int_values_i["i_n"])
        calculate_dyn_temperature_evolution_inside_MEA(dif_eq_i,
                                                        Hgdl, Hmpl, Hacl, Hccl, Hmem,
                                                        nb_gdl, nb_mpl,
                                                        dif_eq_int_values_i["rho_Cp0"],
                                                        heat_i["Jt"], heat_i["Q_r"],
                                                        heat_i["Q_sorp"], heat_i["Q_liq"],
                                                        heat_i["Q_p"], heat_i["Q_e"])
    end
    #       Inside the gas channels and the manifolds
    calculate_dyn_gas_evolution_inside_gas_channel(dif_eq_1D_cell, sv_1D_cell,
                                                    Hagc, Hcgc, Lgc, nb_gc, type_auxiliary,
                                                    flows_1D_GC_manifold["Jv"], flows_1D_GC_manifold["J_H2"],
                                                    flows_1D_GC_manifold["J_O2"], flows_1D_GC_manifold["J_N2"])
    calculate_dyn_liq_evolution_inside_gas_channel(dif_eq_1D_cell, T_des, Hagc, Hcgc, Lgc, nb_gc, flows_1D_GC_manifold["Jl"])
    calculate_dyn_temperature_evolution_inside_gas_channel(dif_eq_1D_cell, nb_gc)
    if type_auxiliary != "no_auxiliary"
        calculate_dyn_manifold_pressure_and_humidity_evolution(dif_eq_1D_manifold,
                                                                T_des, nb_cell, Vasm, Vcsm, Vaem, Vcem, type_auxiliary,
                                                                flows_1D_GC_manifold["W"],
                                                                flows_1D_GC_manifold["Wv"])
    end
    #       Inside the auxiliaries
    if type_auxiliary != "no_auxiliary"
        calculate_dyn_air_compressor_evolution(dif_eq_auxiliary,
                                                parameters["Pacp_des"],
                                                sv_auxiliary["Pasm_out"],
                                                parameters["Pccp_des"],
                                                sv_auxiliary["Pcsm_out"],
                                                type_auxiliary)
        # Note: in the original Python code, `calculate_dyn_humidifier_evolution` is called with
        # `**dif_eq_int_values` where `dif_eq_int_values` is a list, not a dict. This appears to
        # be an inconsistency in the Python source. Here, `parameters` is used as the source of
        # desired injection rates (Wa_inj_des, Wc_inj_des), which is the most plausible intent.
        calculate_dyn_humidifier_evolution(dif_eq_auxiliary,
                                            sv_auxiliary["Wa_inj"],
                                            sv_auxiliary["Wc_inj"],
                                            type_auxiliary,
                                            parameters["Wa_inj_des"],
                                            parameters["Wc_inj_des"])
        calculate_dyn_throttle_area_controler(dif_eq_auxiliary, sv_auxiliary,
                                               Pa_des, Pc_des, A_T_a, A_T_c, type_auxiliary,
                                               flows_1D_GC_manifold["Pagc"], flows_1D_GC_manifold["Pcgc"])
    end

    # All the dif_eq dictionaries are converted because the solver requires an ordered list to work
    dif_eq_global = [dif_eq_1D_cell[i]["d" * key * " / dt"]
                     for i in 1:nb_gc for key in solver_variable_names[1]]
    if type_auxiliary == "forced-convective_cathode_with_flow-through_anode" ||
       type_auxiliary == "forced-convective_cathode_with_anodic_recirculation"
        append!(dif_eq_global,
                [dif_eq_1D_manifold["d" * key * " / dt"] for key in solver_variable_names[2]])
        append!(dif_eq_global,
                [dif_eq_auxiliary["d" * key * " / dt"] for key in solver_variable_names[3]])
    end
    return dif_eq_global
end

