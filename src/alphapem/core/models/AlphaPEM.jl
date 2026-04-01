# -*- coding: utf-8 -*-

"""
- Objectives: Create an open-source software package to simulate the PEM fuel cell for control system applications.
- Authors until V1.0: Raphael GASS, Zhongliang LI, Rachid OUTBIB, Samir JEMEI and Daniel HISSEL.
- Authors from V1.0: Raphael GASS, Dominique GRONDIN, Cedric DAMOUR, Samir JEMEI and Michel BENNE.
---
This file describes the AlphaPEM class, which is a PEM fuel cell system simulator.
The model is one-dimensional, dynamic, biphasic, and isothermal. It has been published in the following articles:
- Gass et al 2024 J. Electrochem. Soc. https://doi.org/10.1149/1945-7111/ad305a
- Gass et al 2024 SSRN http://dx.doi.org/10.2139/ssrn.4812343
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries.
using DifferentialEquations

# Importing constants' value and functions.
include(joinpath(@__DIR__, "dif_eq.jl"))
include(joinpath(@__DIR__, "../modules/cell_voltage_modules.jl"))
include(joinpath(@__DIR__, "../modules/dif_eq_modules.jl"))
include(joinpath(@__DIR__, "../modules/flows_1D_MEA_modules.jl"))
include(joinpath(@__DIR__, "../modules/display_modules.jl"))


# _______________________________________________________AlphaPEM_______________________________________________________


mutable struct AlphaPEM
    fuel_cell::AbstractFuelCell
    current_density::AbstractCurrent
    cfg::SimulationConfig
    solver_variable_names::Vector{Vector{String}}
    all_variable_names::Vector{String}
    variables::Dict
    time_interval::Tuple{Float64, Float64}
    initial_variable_values::Vector
    sol
end


"""Initialise all parameters defining a fuel cell stack operation.

Parameters
----------


Returns
-------
AlphaPEM
    Initialised simulator instance.
"""
function AlphaPEM(fuel_cell::AbstractFuelCell, current_density::AbstractCurrent, cfg::SimulationConfig)::AlphaPEM

    # Initialize the variables' dictionary.
    solver_variable_names = [[
        "C_v_agc", "C_v_agdl", "C_v_ampl", "C_v_acl", "C_v_ccl", "C_v_cmpl", "C_v_cgdl", "C_v_cgc",
        "s_agc", "s_agdl", "s_ampl", "s_acl", "s_ccl", "s_cmpl", "s_cgdl", "s_cgc",
        "lambda_acl", "lambda_mem", "lambda_ccl",
        "C_H2_agc", "C_H2_agdl", "C_H2_ampl", "C_H2_acl",
        "C_O2_ccl", "C_O2_cmpl", "C_O2_cgdl", "C_O2_cgc",
        "C_N2_agc", "C_N2_cgc",
        "T_agc", "T_agdl", "T_ampl", "T_acl", "T_mem", "T_ccl", "T_cmpl", "T_cgdl", "T_cgc",
        "eta_c"
    ]]

    simu = AlphaPEM(
        fuel_cell, #
        current_density, #
        cfg, #
        solver_variable_names, # solver_variable_names::Vector{Vector{String}}
        String[], # all_variable_names::Vector{String}
        Dict{String, Vector{Number}}(), # variables::Dict{String, Vector{Number}}
        (0.0, 0.0), # time_interval::Tuple{Float64, Float64}
        [], # initial_variable_values::Vector{Number}
        nothing, # sol
    )
    # Several points are considered in each GC, GDL and MPL. This must be inserted into the solver_variable_names.
    solver_variable_names_extension!(simu)
    return simu
end


"""Extend `solver_variable_names` for GDL/MPL spatial discretisation.

Parameters
----------
simu : AlphaPEM
    Fuel cell simulator instance.

Returns
-------
Nothing
    The function updates `simu.solver_variable_names` in place.
"""
function solver_variable_names_extension!(simu::AlphaPEM)
    # Several points are considered in each GDL, MPL and GC. They must be inserted into the solver_variable_names.
    new_points_location = ["C_v_agdl", "C_v_ampl", "C_v_cmpl", "C_v_cgdl",
                           "s_agdl", "s_ampl", "s_cmpl", "s_cgdl",
                           "C_H2_agdl", "C_H2_ampl", "C_O2_cmpl", "C_O2_cgdl",
                           "T_agdl", "T_ampl", "T_cmpl", "T_cgdl"]

    for variable in new_points_location
        index = findfirst(==(variable), simu.solver_variable_names[1])
        index === nothing && continue
        # Delete the previous points.
        deleteat!(simu.solver_variable_names[1], index)
        # Increase the number of points.
        if endswith(variable, "gdl")
            splice!(simu.solver_variable_names[1], index:index-1,
                ["$(variable)_$(i)" for i in 1:simu.fuel_cell.physical_parameters.nb_gdl])
        elseif endswith(variable, "mpl")
            splice!(simu.solver_variable_names[1], index:index-1,
                ["$(variable)_$(i)" for i in 1:simu.fuel_cell.physical_parameters.nb_mpl])
        end
    end
    return nothing
end


"""Simulate the model with the given operating inputs.

Parameters
----------
simu : AlphaPEM
    Fuel cell simulator instance.

initial_variable_values : Union{Nothing, Vector}, optional
    Initial values of the solver variables. If `nothing`, values are generated
    from a no-current equilibrium.
time_interval : Union{Nothing, Tuple{Float64, Float64}}, optional
    Time interval for numerical resolution. If `nothing`, it is generated
    according to the chosen current profile.

Returns
-------
Nothing
    The function updates `simu` in place.
"""
function simulate_model!(simu::AlphaPEM,
                         cfg::SimulationConfig,
                         initial_variable_values:: Union{Nothing, Vector}=nothing,
                         time_interval:: Union{Nothing, Tuple{Float64, Float64}}=nothing)

    # General warnings.
    if cfg.type_fuel_cell in (:EH_31_1_5, :EH_31_2_0, :E_H31_2_25, :EH_31_2_5)
        println("Warning: EH-Group fuel cell examples may be outdated. Using ZSW-GenStack is recommended.\n")
    end

    if cfg.voltage_zone == :EIS
        throw(ArgumentError("The EIS generation is currently undergoing maintenance."))
    end

    if cfg.type_auxiliary in (:forced_convective_cathode_with_anodic_recirculation,
                              :forced_convective_cathode_with_flow_through_anode)
        cfg.type_auxiliary = :no_auxiliary
        println("Warning: auxiliaries were temporarily removed; \"no_auxiliary\" is automatically used.\n")
    end

    if simu.fuel_cell.operating_conditions.Pa_des < Pext || simu.fuel_cell.operating_conditions.Pc_des < Pext
        throw(ArgumentError("The desired pressure is too low. It cannot be lower than the pressure outside the stack."))
    end

    # Initialize the variables' dictionaries.
    if cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                              :forced_convective_cathode_with_anodic_recirculation)
        push!(simu.solver_variable_names,
              ["Pasm", "Paem", "Pcsm", "Pcem", "Phi_asm", "Phi_aem", "Phi_csm", "Phi_cem"])
        push!(simu.solver_variable_names, ["Wcp", "Wa_inj", "Wc_inj", "Abp_a", "Abp_c"])
    end
    simu.all_variable_names = vcat(
        reduce(vcat, simu.solver_variable_names),
        ["t", "i_fc", "C_O2_Pt", "Ucell", "v_a", "v_c", "Pa_in", "Pc_in"],
        ["Phi_a_des", "Phi_c_des"],
    )
    simu.variables = Dict{String, Vector{Number}}(k => Number[] for k in simu.all_variable_names)

    # Create the dynamic evolution.
    #       Create time intervals
    simu.time_interval = time_interval === nothing ? simu.current_density.time_interval : time_interval
    #       Create the initial variable values
    simu.initial_variable_values = initial_variable_values === nothing ?
                                    create_initial_variable_values(simu) : initial_variable_values

    #       Solve the differential equation system.
    #           Pack external data passed to the ODE right\-hand side.
    packed = (fuel_cell=simu.fuel_cell, solver_variable_names=simu.solver_variable_names)
    #           Define RHS in SciML signature: f(y, p, t) -> dy/dt.
    rhs = (y, p, t) -> dydt(t, y, p.fuel_cell, p.solver_variable_names)
    #           Stop integration as soon as a critical variable becomes negative.
    #           The callback monitors the event condition continuously in time.
    condition = (y, t, integ) -> event_negative(t, y, integ.p.fuel_cell, integ.p.solver_variable_names)
    cb_negative = ContinuousCallback(condition, integ -> terminate!(integ))
    #           Build and solve the ODE problem with FBDF for stiff dynamics.
    prob = ODEProblem(rhs, simu.initial_variable_values, simu.time_interval, packed)
    simu.sol = solve(prob, FBDF(); reltol=simu.parameters["rtol"], abstol=simu.parameters["atol"],
                     callback=cb_negative)

    #       Recover the variable values calculated by the solver into the dictionary.
    recovery!(simu)
    return nothing
end


"""Calculate the time intervals for numerical resolution, according to the current chosen,
if it is not provided.

Parameters
----------


Returns
-------
Tuple{Float64, Float64}
    Time interval for numerical resolution. It is used when `time_interval == nothing`.
"""
function create_time_interval(current_density::AbstractCurrent)::Tuple{Float64, Float64}

    # Recovery of the good time interval.
    if type_current == "step"
        t0_interval = 0.0 # s.
        tf_interval = step_current_parameters["delta_t_ini_step"] +
                      step_current_parameters["delta_t_load_step"] +
                      step_current_parameters["delta_t_break_step"] # s.
    elseif type_current == "polarization"
        # Extraction of the parameters.
        delta_t_ini_pola = pola_current_parameters["delta_t_ini_pola"]
        delta_t_load_pola = pola_current_parameters["delta_t_load_pola"]
        delta_t_break_pola = pola_current_parameters["delta_t_break_pola"]
        delta_i_pola = pola_current_parameters["delta_i_pola"]
        i_max_pola = pola_current_parameters["i_max_pola"]
        # Calculation.
        t0_interval = 0.0 # s.
        tf_interval = delta_t_ini_pola + Int(i_max_pola / delta_i_pola) * (delta_t_load_pola + delta_t_break_pola)
    elseif type_current == "polarization_for_cali"
        # Extraction of the parameters.
        delta_t_ini_pola_cali = pola_current_for_cali_parameters["delta_t_ini_pola_cali"]
        delta_t_load_pola_cali = pola_current_for_cali_parameters["delta_t_load_pola_cali"]
        delta_t_break_pola_cali = pola_current_for_cali_parameters["delta_t_break_pola_cali"]
        i_exp_cali_t, _ = simu.fuel_cell.pola_exp_data_cali
        # Calculation.
        delta_t_pola_cali = delta_t_load_pola_cali + delta_t_break_pola_cali
        t0_interval = 0.0
        tf_interval = delta_t_ini_pola_cali + length(i_exp_cali_t) * delta_t_pola_cali # s.
    else
        throw(ArgumentError("Please enter a recognized type_current option for calculating the time interval."))
    end
    return (t0_interval, tf_interval)
end


"""Create the initial values of the solver variables if they are not provided.
It is generated considering an equilibrium inside the fuel cell with H2, O2 and N2,
at the external pressure, humidity and temperature, without flow or current.

Parameters
----------
simu : AlphaPEM
    Fuel cell simulator instance.

Returns
-------
Vector
    Initial values of the solver variables. It is used when
    `initial_variable_values == nothing`.
"""
function create_initial_variable_values(simu::AlphaPEM)::Vector
    # Extraction of the operating inputs and parameters.
    current_density, T_des = simu.operating_inputs["current_density"], simu.operating_inputs["T_des"]
    Pa_des, Pc_des = simu.fuel_cell.operating_conditions.Pa_des, simu.fuel_cell.operating_conditions.Pc_des
    Phi_a_des, Phi_c_des = simu.operating_inputs["Phi_a_des"], simu.operating_inputs["Phi_c_des"]
    y_H2_in = simu.operating_inputs["y_H2_in"]
    Hmem, kappa_co, kappa_c = simu.parameters["Hmem"], simu.parameters["kappa_co"], simu.parameters["kappa_c"]
    i0_c_ref = simu.parameters["i0_c_ref"]
    nb_gc, nb_gdl, nb_mpl = simu.parameters["nb_gc"], simu.parameters["nb_gdl"], simu.parameters["nb_mpl"]

    # Initial fuel cell states.
    #   Intermediate values.
    T_ini = T_des
    if cfg.type_auxiliary in (:forced_convective_cathode_with_anodic_recirculation,
                              :forced_convective_cathode_with_flow_through_anode)
        Pa_ini, Pc_ini = Pext, Pext
        Phi_a_ini, Phi_c_ini = Phi_ext, Phi_ext
    else # cfg.type_auxiliary == "no_auxiliaries".
        Pa_ini, Pc_ini = Pa_des, Pc_des
        Phi_a_ini, Phi_c_ini = Phi_a_des, Phi_c_des
    end

    Psat_ini = 101325 * 10 ^ (-2.1794 + 0.02953 * (T_ini - 273.15) - 9.1837e-5 * (T_ini - 273.15)^2 +
                              1.4454e-7 * (T_ini - 273.15)^3)
    #   Initial fuel cell states.
    C_v_a_ini = Phi_a_ini * Psat_ini / (R * T_ini)
    C_v_c_ini = Phi_c_ini * Psat_ini / (R * T_ini)
    C_O2_ini = y_O2_ext * (Pc_ini - Phi_c_ini * Psat_ini) / (R * T_ini)
    C_N2_cgc_ini = (1 - y_O2_ext) * (Pc_ini - Phi_c_ini * Psat_ini) / (R * T_ini)

    if cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        C_H2_ini = y_H2_in * (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)
        C_N2_agc_ini = (1 - y_H2_in) * (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)
    else
        C_H2_ini = (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)
        C_N2_agc_ini = 0.0
    end

    s_ini = 0.0
    lambda_mem_ini = lambda_eq(C_v_c_ini, s_ini, T_ini)
    i_fc_ini = current_density(simu.time_interval[1], simu.parameters)
    i_n_ini = 2 * F * R * T_ini / Hmem * C_H2_ini * k_H2(lambda_mem_ini, T_ini, kappa_co) +
              4 * F * R * T_ini / Hmem * C_O2_ini * k_O2(lambda_mem_ini, T_ini, kappa_co)
    eta_c_ini = R * T_ini / (alpha_c * F) * log((i_fc_ini + i_n_ini) / i0_c_ref *
                                                 1 / exp(-Eact_O2_red / (R * T_ini) * (1 / T_ini - 1 / Tref_O2_red)) *
                                                 (C_O2ref_red / C_O2_ini)^kappa_c)

    # Initial auxiliary system state.
    Wcp_ini = 0.0
    Wa_inj_ini = 0.0
    Wc_inj_ini = 0.0
    Abp_a_ini = 0.0
    Abp_c_ini = 0.0

    # Main variable initialization.
    C_v_agc, C_v_agdl, C_v_ampl, C_v_acl = (C_v_a_ini, C_v_a_ini, C_v_a_ini, C_v_a_ini)
    C_v_ccl, C_v_cmpl, C_v_cgdl, C_v_cgc = (C_v_c_ini, C_v_c_ini, C_v_c_ini, C_v_c_ini)
    s_agc, s_agdl, s_ampl, s_acl, s_ccl, s_cmpl, s_cgdl, s_cgc = (s_ini, s_ini, s_ini, s_ini, s_ini, s_ini, s_ini, s_ini)
    lambda_acl, lambda_mem, lambda_ccl = (lambda_mem_ini, lambda_mem_ini, lambda_mem_ini)
    C_H2_agc, C_H2_agdl, C_H2_ampl, C_H2_acl = (C_H2_ini, C_H2_ini, C_H2_ini, C_H2_ini)
    C_O2_ccl, C_O2_cmpl, C_O2_cgdl, C_O2_cgc = (C_O2_ini, C_O2_ini, C_O2_ini, C_O2_ini)
    if cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        C_N2_agc, C_N2_cgc = C_N2_agc_ini, C_N2_cgc_ini
    else
        C_N2_agc, C_N2_cgc = 0.0, C_N2_cgc_ini
    end
    T_agc, T_agdl, T_ampl, T_acl, T_mem, T_ccl, T_cmpl, T_cgdl, T_cgc = (T_ini, T_ini, T_ini, T_ini, T_ini, T_ini, T_ini, T_ini, T_ini)
    eta_c = eta_c_ini
    Pasm, Paem = Pa_ini, Pa_ini
    Pcsm, Pcem = Pc_ini, Pc_ini
    Phi_asm, Phi_aem = Phi_a_ini, Phi_a_ini
    Phi_csm, Phi_cem = Phi_c_ini, Phi_c_ini
    Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c = Wcp_ini, Wa_inj_ini, Wc_inj_ini, Abp_a_ini, Abp_c_ini

    # Gathering of the variables initial value into one list, only for one gas channel node.
    initial_variable_values_1D = vcat(
        [C_v_agc], fill(C_v_agdl, nb_gdl), fill(C_v_ampl, nb_mpl), [C_v_acl, C_v_ccl],
        fill(C_v_cmpl, nb_mpl), fill(C_v_cgdl, nb_gdl), [C_v_cgc],
        [s_agc], fill(s_agdl, nb_gdl), fill(s_ampl, nb_mpl), [s_acl, s_ccl],
        fill(s_cmpl, nb_mpl), fill(s_cgdl, nb_gdl), [s_cgc],
        [lambda_acl, lambda_mem, lambda_ccl],
        [C_H2_agc], fill(C_H2_agdl, nb_gdl), fill(C_H2_ampl, nb_mpl), [C_H2_acl, C_O2_ccl],
        fill(C_O2_cmpl, nb_mpl), fill(C_O2_cgdl, nb_gdl), [C_O2_cgc], [C_N2_agc], [C_N2_cgc],
        [T_agc], fill(T_agdl, nb_gdl), fill(T_ampl, nb_mpl), [T_acl], [T_mem], [T_ccl],
        fill(T_cmpl, nb_mpl), fill(T_cgdl, nb_gdl), [T_cgc], [eta_c],
    )
    # Replication for each gas channel node.
    initial_variable_values = repeat(initial_variable_values_1D, nb_gc)

    # Addition of the auxiliary system initial states.
    if cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                              :forced_convective_cathode_with_anodic_recirculation)
        append!(initial_variable_values, [Pasm, Paem, Pcsm, Pcem, Phi_asm, Phi_aem, Phi_csm, Phi_cem,
                                          Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c])
    end
    return initial_variable_values
end


"""Recover the values calculated by the solver and add them into the variables dictionary.
Some internal states are rebuilt manually because they are not all directly exported
by the numerical solver.

Parameters
----------
simu : AlphaPEM
    Fuel cell simulator instance.

Returns
-------
Nothing
    The function updates `simu.variables` in place.
"""
function recovery!(simu::AlphaPEM)
    # Recovery of the time span.
    simu.variables["t"] = collect(simu.sol.t)

    # Recovery of the main variables dynamic evolution.
    nb_gc = simu.parameters["nb_gc"]
    for (index, key) in enumerate(simu.solver_variable_names[1]) # recovery of the MEA and GC variables
        simu.variables[key] = [[simu.sol.u[j][index + (i - 1) * length(simu.solver_variable_names[1])] for j in eachindex(simu.sol.u)]
                                for i in 1:nb_gc]
    end
    if cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                              :forced_convective_cathode_with_anodic_recirculation)
        for (index, key) in enumerate(simu.solver_variable_names[2]) # recovery of the manifold variables
            simu.variables[key] = [simu.sol.u[j][index + nb_gc * length(simu.solver_variable_names[1])]
                                    for j in eachindex(simu.sol.u)]
        end
        for (index, key) in enumerate(simu.solver_variable_names[3]) # recovery of the auxiliary variables
            simu.variables[key] = [simu.sol.u[j][index + nb_gc * length(simu.solver_variable_names[1]) +
                                                  length(simu.solver_variable_names[2])]
                                    for j in eachindex(simu.sol.u)]
        end
    end

    # Recovery of more variables.
    simu.variables["v_a"] = [[] for _ in 1:nb_gc]
    simu.variables["v_c"] = [[] for _ in 1:nb_gc]
    simu.variables["C_O2_Pt"] = [[] for _ in 1:nb_gc]
    simu.variables["i_fc"] = [[] for _ in 1:nb_gc]

    for (j, t_j) in enumerate(simu.variables["t"])
        # ... recovery of the variables inside the MEA 1D line.
        solver_variables_1D_MEA = [Dict() for _ in 1:nb_gc]
        for k in 1:nb_gc
            for (index, variable) in enumerate(simu.solver_variable_names[1])
                solver_variables_1D_MEA[k][variable] = simu.sol.u[j][index + (k - 1) * length(simu.solver_variable_names[1])]
            end
        end

        # ... recovery of i_fc and C_O2_Pt.
        i_fc_cell = simu.operating_inputs["current_density"](t_j, simu.parameters)
        i_fc = calculate_1D_GC_current_density(i_fc_cell, solver_variables_1D_MEA, simu.parameters)
        for k in 1:nb_gc
            push!(simu.variables["i_fc"][k], i_fc[k])
            push!(simu.variables["C_O2_Pt"][k], calculate_C_O2_Pt(i_fc[k], solver_variables_1D_MEA[k], simu.parameters))
        end

        # ... recovery of Ucell, v_a, v_c, Pa_in and Pc_in.
        push!(simu.variables["Ucell"], calculate_cell_voltage(i_fc[1], simu.variables["C_O2_Pt"][1][end],
                                                                solver_variables_1D_MEA[1], simu.parameters))
        v_a, v_c, Pa_in, Pc_in = calculate_velocity_evolution(solver_variables_1D_MEA, i_fc_cell,
                                                              simu.operating_inputs, simu.parameters)
        for k in 1:nb_gc
            push!(simu.variables["v_a"][k], v_a[k])
            push!(simu.variables["v_c"][k], v_c[k])
        end
        push!(simu.variables["Pa_in"], Pa_in)
        push!(simu.variables["Pc_in"], Pc_in)
    end
    return nothing
end


"""Display the plots of the program.

Parameters
----------
simu : AlphaPEM
    Fuel cell simulator instance.
ax1 : optional
    Axes for the first set of plots. Default is `nothing`.
ax2 : optional
    Axes for the second set of plots. Default is `nothing`.
ax3 : optional
    Axes for the third set of plots. Default is `nothing`.

Returns
-------
Nothing
"""
function Display(simu::AlphaPEM, ax1=nothing, ax2=nothing, ax3=nothing)
    # Folder name.
    subfolder_name = split(cfg.type_fuel_cell, '_')[1]

    # Display.
    if type_current == "step"
        if type_display == "multiple"
            figs_axes = [plt.subplots(figsize=(8, 8)) for _ in 1:11]
            figs = [fa[1] for fa in figs_axes]
            axes = [fa[2] for fa in figs_axes]
            plot_ifc_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, axes[1])
            plot_C_v_1D_temporal(simu.variables, simu.parameters, axes[2])
            plot_lambda_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, axes[3])
            plot_s_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, axes[4])
            plot_C_O2_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, axes[5])
            plot_C_H2_1D_temporal(simu.variables, simu.parameters, axes[6])
            plot_C_N2_1D_temporal(simu.variables, simu.parameters, axes[7])
            plot_T_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, axes[8])
            plot_Ucell(simu.variables, simu.parameters, axes[9])
            plot_P_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, axes[10])
            plot_v_1D_temporal(simu.variables, simu.parameters, axes[11])

            # Keep the same special handling as Python for step+multiple outputs.
            Saving_instructions(simu, "results", subfolder_name, "step_current_ifc_1.pdf", figs[1])
            Saving_instructions(simu, "results", subfolder_name, "step_current_Cv_1.pdf", figs[2])
            Saving_instructions(simu, "results", subfolder_name, "step_current_lambda_1.pdf", figs[3])
            Saving_instructions(simu, "results", subfolder_name, "step_current_s_1.pdf", figs[4])
            Saving_instructions(simu, "results", subfolder_name, "step_current_C_O2_1.pdf", figs[5])
            Saving_instructions(simu, "results", subfolder_name, "step_current_C_H2_1.pdf", figs[6])
            Saving_instructions(simu, "results", subfolder_name, "step_current_C_N2_1.pdf", figs[7])
            Saving_instructions(simu, "results", subfolder_name, "step_current_T_1.pdf", figs[8])
            Saving_instructions(simu, "results", subfolder_name, "step_current_Ucell_1.pdf", figs[9])
            Saving_instructions(simu, "results", subfolder_name, "step_current_P_1.pdf", figs[10])
            Saving_instructions(simu, "results", subfolder_name, "step_current_v_1.pdf", figs[11])

            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif type_display == "synthetic"
            plot_ifc_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[1, 1])
            plot_Ucell(simu.variables, simu.parameters, ax1[1, 2])
            plot_T_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[1, 3])

            plot_C_v_1D_temporal(simu.variables, simu.parameters, ax1[2, 1])
            plot_s_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[2, 2])
            plot_lambda_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[2, 3])

            plot_C_H2_1D_temporal(simu.variables, simu.parameters, ax1[3, 1])
            plot_C_O2_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[3, 2])
            plot_P_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[3, 3])

            if simu.parameters["type_plot"] == "fixed"
                plot_T_pseudo_2D_final(simu.variables, simu.operating_inputs, simu.parameters, ax2)
            end

            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(1.0)
        end
    elseif type_current == "polarization"
        if type_display == "multiple"
            plot_polarisation_curve(simu.variables, simu.operating_inputs, simu.parameters, ax1[1])
            plot_power_density_curve(simu.variables, simu.operating_inputs, simu.parameters, length(simu.variables["t"]), ax1[2])
            plot_cell_efficiency(simu.variables, simu.operating_inputs, simu.parameters, length(simu.variables["t"]), ax1[3])

            plot_lambda_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax2[2])
            plot_s_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax2[3])
            plot_T_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax2[4])
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif type_display == "synthetic"
            plot_polarisation_curve(simu.variables, simu.operating_inputs, simu.parameters, ax1)
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif type_display == "no_display"
            plot_polarisation_curve(simu.variables, simu.operating_inputs, simu.parameters, ax1, false)
        end
    elseif type_current == "polarization_for_cali"
        if type_display == "multiple"
            plot_polarisation_curve_for_cali(simu.variables, simu.operating_inputs, simu.parameters, ax1[1])
            plot_lambda_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[2])
            plot_s_1D_temporal(simu.variables, simu.operating_inputs, simu.parameters, ax1[3])
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif type_display == "synthetic"
            plot_polarisation_curve_for_cali(simu.variables, simu.operating_inputs, simu.parameters, ax1)
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        end
    elseif type_current == "EIS"
        Fourier_results = make_Fourier_transformation(simu.variables, simu.operating_inputs, simu.parameters)
        if type_display == "multiple"
            plot_EIS_curve_Nyquist(simu.parameters, Fourier_results, ax1)
            plot_EIS_curve_Bode_amplitude(simu.parameters, Fourier_results, ax2)
            plot_EIS_curve_Bode_angle(simu.parameters, Fourier_results, ax3)
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif type_display == "synthetic"
            plot_EIS_curve_Nyquist(simu.parameters, Fourier_results, ax1[1])
            plot_EIS_curve_Bode_amplitude(simu.parameters, Fourier_results, ax1[2])
            plot_EIS_curve_Bode_angle(simu.parameters, Fourier_results, ax1[3])
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        end
    end
    return nothing
end


"""Save the plots. The filenames are generated according to `type_current`
and `type_display`.

Parameters
----------
simu : AlphaPEM
    Fuel cell simulator instance.
fig1 : optional
    Figure for the first plot. Default is `nothing`.
fig2 : optional
    Figure for the second plot. Default is `nothing`.
fig3 : optional
    Figure for the third plot. Default is `nothing`.

Returns
-------
Nothing
"""
function Save_plot(simu::AlphaPEM, fig1=nothing, fig2=nothing, fig3=nothing)
    # Folder name.
    subfolder_name = split(cfg.type_fuel_cell, '_')[1]

    # For the step current.
    if type_current == "step"
        if type_display == "synthetic"
            Saving_instructions(simu, "results", subfolder_name, "step_current_syn_1.pdf", fig1)
            type_plot == "fixed" && Saving_instructions(simu, "results", subfolder_name, "final_temperature_dist_1.pdf", fig2)
        end
    # For the polarization curve.
    elseif type_current == "polarization"
        if type_display == "multiple"
            Saving_instructions(simu, "results", subfolder_name, "global_indicators_1.pdf", fig1)
            Saving_instructions(simu, "results", subfolder_name, "pola_curve_syn_1.pdf", fig2)
        elseif type_display == "synthetic"
            Saving_instructions(simu, "results", subfolder_name, "pola_curve_1.pdf", fig1)
        end
    # For the EIS curve.
    elseif type_current == "EIS"
        if type_display == "multiple"
            Saving_instructions(simu, "results", subfolder_name, "Nyquist_plot_1.pdf", fig1)
            Saving_instructions(simu, "results", subfolder_name, "Bode_amplitude_curve_1.pdf", fig2)
            Saving_instructions(simu, "results", subfolder_name, "Bode_angle_curve_1.pdf", fig3)
        elseif type_display == "synthetic"
            Saving_instructions(simu, "results", subfolder_name, "Nyquist_plot_syn_1.pdf", fig1)
        end
    # For the polarization curve for calibration.
    elseif type_current == "polarization_for_cali"
        if type_display == "multiple"
            Saving_instructions(simu, "results", subfolder_name, "impact_cali_on_internal_state_1.pdf", fig1)
        elseif type_display == "synthetic"
            Saving_instructions(simu, "results", subfolder_name, "pola_curve_cali_1.pdf", fig1)
        end
    end
    return nothing
end


"""Give the saving instructions for figures.

Parameters
----------
simu : AlphaPEM
    Fuel cell simulator instance.
root_folder : String
    Root folder for saving.
subfolder_name : String
    Subfolder name for saving.
filename : String
    Filename for saving.
fig :
    Figure object to save.

Returns
-------
Nothing
"""
function Saving_instructions(simu::AlphaPEM,
                             root_folder::String,
                             subfolder_name::String,
                             filename::String,
                             fig)
    # Resolve current file and define repository markers to locate project root.
    cur = abspath(@__FILE__)
    markers = [".git", "pyproject.toml", "setup.cfg", "requirements.txt", "Pipfile"]
    project_root = nothing
    parent = dirname(cur)
    while true
        any(ispath(joinpath(parent, m)) for m in markers) && (project_root = parent; break)
        new_parent = dirname(parent)
        new_parent == parent && break
        parent = new_parent
    end
    # Fallback to current working directory if no marker found.
    project_root === nothing && (project_root = pwd())

    # Build destination folder under project root and create it.
    folder_path = joinpath(project_root, root_folder, subfolder_name)
    mkpath(folder_path)

    # Prepare a safe filename: keep extension and append _N before extension if collision.
    file_path = joinpath(folder_path, filename)
    if isfile(file_path)
        stem, suffix = splitext(filename)
        suffix = isempty(suffix) ? ".pdf" : suffix
        counter = 1
        while true
            candidate = joinpath(folder_path, "$(stem)_$(counter)$(suffix)")
            !isfile(candidate) && (file_path = candidate; break)
            counter += 1
        end
    end

    # Save the figure with the same parameters as before.
    fig.savefig(file_path, dpi=900, transparent=false, bbox_inches="tight")
    return nothing
end

