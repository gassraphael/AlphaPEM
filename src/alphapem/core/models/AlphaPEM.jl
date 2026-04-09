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

# Creating the variable names
const MANIFOLD_SOLVER_VARIABLE_NAMES  = ["Pasm", "Paem", "Pcsm", "Pcem", "Phi_asm", "Phi_aem", "Phi_csm", "Phi_cem"]
const AUXILIARY_SOLVER_VARIABLE_NAMES = ["Wcp", "Wa_inj", "Wc_inj", "Abp_a", "Abp_c"]
const DERIVED_VARIABLE_NAMES          = ["t", "i_fc", "C_O2_Pt", "Ucell", "v_a", "v_c", "Pa_in", "Pc_in"]

# _______________________________________________________AlphaPEM_______________________________________________________


mutable struct AlphaPEM
    fuel_cell::AbstractFuelCell
    current_density::AbstractCurrent
    cfg::SimulationConfig
    variables::Dict
    time_interval::Tuple{Float64, Float64}
    initial_variable_values::Vector
    sol
end


"""Create an `AlphaPEM` simulator from typed fuel-cell, current-profile, and configuration objects.

# Arguments
- `fuel_cell::AbstractFuelCell`: Fuel-cell definition with physical, operating, and numerical parameters.
- `current_density::AbstractCurrent`: Current-density profile used during the simulation.
- `cfg::SimulationConfig`: Global simulation configuration.

# Returns
- `AlphaPEM`: Initialised simulator instance.
"""
function AlphaPEM(fuel_cell::AbstractFuelCell, current_density::AbstractCurrent, cfg::SimulationConfig)::AlphaPEM

    simu = AlphaPEM(
        fuel_cell, #
        current_density, #
        cfg, #
        Dict{String, Any}(), # variables::Dict
        (0.0, 0.0), # time_interval::Tuple{Float64, Float64}
        [], # initial_variable_values::Vector{Number}
        nothing, # sol
    )
    return simu
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
                         initial_variable_values:: Union{Nothing, Vector}=nothing,
                         time_interval:: Union{Nothing, Tuple{Float64, Float64}}=nothing)

    # General warnings.
    if simu.cfg.type_fuel_cell in (:EH_31_1_5, :EH_31_2_0, :E_H31_2_25, :EH_31_2_5)
        println("Warning: EH-Group fuel cell examples may be outdated. Using ZSW-GenStack is recommended.\n")
    end

    if simu.cfg.voltage_zone == :EIS
        throw(ArgumentError("The EIS generation is currently undergoing maintenance."))
    end

    if simu.cfg.type_auxiliary in (:forced_convective_cathode_with_anodic_recirculation,
                              :forced_convective_cathode_with_flow_through_anode)
        simu.cfg.type_auxiliary = :no_auxiliary
        println("Warning: auxiliaries were temporarily removed; \"no_auxiliary\" is automatically used.\n")
    end

    if simu.fuel_cell.operating_conditions.Pa_des < Pext || simu.fuel_cell.operating_conditions.Pc_des < Pext
        throw(ArgumentError("The desired pressure is too low. It cannot be lower than the pressure outside the stack."))
    end

    # Initialize the variables' dictionaries.
    has_auxiliary = simu.cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                                                :forced_convective_cathode_with_anodic_recirculation)
    canonical_names = canonical_mea_solver_variable_names(simu.fuel_cell.numerical_parameters.nb_gdl,
                                                          simu.fuel_cell.numerical_parameters.nb_mpl)
    manifold_names = has_auxiliary ? MANIFOLD_SOLVER_VARIABLE_NAMES : String[]
    auxiliary_names = has_auxiliary ? AUXILIARY_SOLVER_VARIABLE_NAMES : String[]
    all_variable_names = vcat(
        canonical_names,
        manifold_names,
        auxiliary_names,
        DERIVED_VARIABLE_NAMES,
    )
    simu.variables = Dict{String, Any}(k => Number[] for k in all_variable_names)

    # Create the dynamic evolution.
    #       Create time intervals
    simu.time_interval = time_interval === nothing ? simu.current_density.time_interval : time_interval
    #       Create the initial variable values
    simu.initial_variable_values = initial_variable_values === nothing ?
                                    create_initial_variable_values(simu) : initial_variable_values

    #       Solve the differential equation system.
    #           Pack external data passed to the ODE right\-hand side.
    packed = (fuel_cell=simu.fuel_cell, current_density=simu.current_density, cfg=simu.cfg)
    #           Define RHS in SciML signature: f(y, p, t) -> dy/dt.
    rhs = (y, p, t) -> dydt(t, y, p.fuel_cell, p.current_density, p.cfg)
    #           Build and solve the ODE problem with FBDF for stiff dynamics.
    prob = ODEProblem(rhs, simu.initial_variable_values, simu.time_interval, packed)
    simu.sol = solve(prob, FBDF(autodiff=false); reltol=simu.fuel_cell.numerical_parameters.rtol,
                     abstol=simu.fuel_cell.numerical_parameters.atol)

    #       Recover the variable values calculated by the solver into the dictionary.
    recovery!(simu)
    return nothing
end


"""Build the default initial state vector for the ODE solver.

The initial state corresponds to an equilibrium condition inside the fuel cell with
hydrogen, oxygen, and nitrogen at the external pressure, humidity, and temperature,
without flow or load current.

# Arguments
- `simu::AlphaPEM`: Fuel-cell simulator instance.

# Returns
- `Vector`: Initial values of the solver variables.
"""
function create_initial_variable_values(simu::AlphaPEM)::Vector
    # Extraction of the parameter classes for better readability.
    oc = simu.fuel_cell.operating_conditions
    pp = simu.fuel_cell.physical_parameters
    np = simu.fuel_cell.numerical_parameters
    # Extraction of frequently used parameters
    T_des, Pa_des, Pc_des = oc.T_des, oc.Pa_des, oc.Pc_des
    Phi_a_des, Phi_c_des, y_H2_in = oc.Phi_a_des, oc.Phi_c_des, oc.y_H2_in
    Hmem, kappa_co, kappa_c, i0_c_ref = pp.Hmem, pp.kappa_co, pp.kappa_c, pp.i0_c_ref
    nb_gc, nb_gdl, nb_mpl = np.nb_gc, np.nb_gdl, np.nb_mpl

    # Initial fuel cell states.
    #   Intermediate values.
    T_ini = T_des
    if simu.cfg.type_auxiliary in (:forced_convective_cathode_with_anodic_recirculation,
                              :forced_convective_cathode_with_flow_through_anode)
        Pa_ini, Pc_ini = Pext, Pext
        Phi_a_ini, Phi_c_ini = Phi_ext, Phi_ext
    else
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

    if simu.cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        C_H2_ini = y_H2_in * (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)
        C_N2_agc_ini = (1 - y_H2_in) * (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)
    else
        C_H2_ini = (Pa_ini - Phi_a_ini * Psat_ini) / (R * T_ini)
        C_N2_agc_ini = 0.0
    end

    s_ini = 0.0
    lambda_mem_ini = lambda_eq(C_v_c_ini, s_ini, T_ini)
    i_fc_ini = current(simu.current_density, simu.time_interval[1])
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
    if simu.cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
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

    # Gather initial values in the canonical typed solver ordering.
    names_1D = canonical_mea_solver_variable_names(nb_gdl, nb_mpl)
    values_1D = Dict{String, Float64}(
        "C_v_agc" => C_v_agc, "C_v_acl" => C_v_acl, "C_v_ccl" => C_v_ccl, "C_v_cgc" => C_v_cgc,
        "s_agc" => s_agc, "s_acl" => s_acl, "s_ccl" => s_ccl, "s_cgc" => s_cgc,
        "lambda_acl" => lambda_acl, "lambda_mem" => lambda_mem, "lambda_ccl" => lambda_ccl,
        "C_H2_agc" => C_H2_agc, "C_H2_acl" => C_H2_acl,
        "C_O2_ccl" => C_O2_ccl, "C_O2_cgc" => C_O2_cgc,
        "C_N2_agc" => C_N2_agc, "C_N2_cgc" => C_N2_cgc,
        "T_agc" => T_agc, "T_acl" => T_acl, "T_mem" => T_mem, "T_ccl" => T_ccl, "T_cgc" => T_cgc,
        "eta_c" => eta_c,
    )
    for i in 1:nb_gdl
        values_1D["C_v_agdl_$(i)"] = C_v_agdl
        values_1D["C_v_cgdl_$(i)"] = C_v_cgdl
        values_1D["s_agdl_$(i)"] = s_agdl
        values_1D["s_cgdl_$(i)"] = s_cgdl
        values_1D["C_H2_agdl_$(i)"] = C_H2_agdl
        values_1D["C_O2_cgdl_$(i)"] = C_O2_cgdl
        values_1D["T_agdl_$(i)"] = T_agdl
        values_1D["T_cgdl_$(i)"] = T_cgdl
    end
    for i in 1:nb_mpl
        values_1D["C_v_ampl_$(i)"] = C_v_ampl
        values_1D["C_v_cmpl_$(i)"] = C_v_cmpl
        values_1D["s_ampl_$(i)"] = s_ampl
        values_1D["s_cmpl_$(i)"] = s_cmpl
        values_1D["C_H2_ampl_$(i)"] = C_H2_ampl
        values_1D["C_O2_cmpl_$(i)"] = C_O2_cmpl
        values_1D["T_ampl_$(i)"] = T_ampl
        values_1D["T_cmpl_$(i)"] = T_cmpl
    end
    initial_variable_values_1D = [values_1D[name] for name in names_1D]
    # Replication for each gas channel node.
    initial_variable_values = repeat(initial_variable_values_1D, nb_gc)

    # Addition of the auxiliary system initial states.
    if simu.cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                              :forced_convective_cathode_with_anodic_recirculation)
        append!(initial_variable_values, [Pasm, Paem, Pcsm, Pcem, Phi_asm, Phi_aem, Phi_csm, Phi_cem,
                                          Wcp, Wa_inj, Wc_inj, Abp_a, Abp_c])
    end
    return initial_variable_values
end


"""Populate `simu.variables` from the solver output.

Some derived internal quantities are rebuilt manually because they are not directly
stored in the numerical solution object.

# Arguments
- `simu::AlphaPEM`: Fuel-cell simulator instance.

# Returns
- `Nothing`: The function updates `simu.variables` in place.
"""
function recovery!(simu::AlphaPEM)
    # Recovery of the time span.
    simu.variables["t"] = collect(simu.sol.t)

    # Recovery of the main variables dynamic evolution.
    np = simu.fuel_cell.numerical_parameters
    nb_gc, nb_gdl, nb_mpl = np.nb_gc, np.nb_gdl, np.nb_mpl
    n_vars_mea_1D = _nb_solver_vars_per_gc(nb_gdl, nb_mpl)
    canonical_names = canonical_mea_solver_variable_names(nb_gdl, nb_mpl)
    length(canonical_names) == n_vars_mea_1D ||
        throw(ArgumentError("Canonical MEA layout size mismatch in recovery!."))

    for (index, key) in enumerate(canonical_names) # recovery of MEA and GC variables in canonical order
        simu.variables[key] = [[simu.sol.u[j][index + (i - 1) * n_vars_mea_1D] for j in eachindex(simu.sol.u)]
                               for i in 1:nb_gc]
    end
    if simu.cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                              :forced_convective_cathode_with_anodic_recirculation)
        for (index, key) in enumerate(MANIFOLD_SOLVER_VARIABLE_NAMES) # recovery of the manifold variables
            simu.variables[key] = [simu.sol.u[j][index + nb_gc * n_vars_mea_1D] for j in eachindex(simu.sol.u)]
        end
        n_vars2 = length(MANIFOLD_SOLVER_VARIABLE_NAMES)
        for (index, key) in enumerate(AUXILIARY_SOLVER_VARIABLE_NAMES) # recovery of the auxiliary variables
            simu.variables[key] = [simu.sol.u[j][index + nb_gc * n_vars_mea_1D + n_vars2] for j in eachindex(simu.sol.u)]
        end
    end

    # Recovery of more variables.
    simu.variables["v_a"] = [[] for _ in 1:nb_gc]
    simu.variables["v_c"] = [[] for _ in 1:nb_gc]
    simu.variables["C_O2_Pt"] = [[] for _ in 1:nb_gc]
    simu.variables["i_fc"] = [[] for _ in 1:nb_gc]


    for (j, t_j) in enumerate(simu.variables["t"])
        # ... recovery of the variables inside the MEA 1D line.
        solver_variables_1D_MEA = [_unpack_mea_state_1D(@view(simu.sol.u[j][(k - 1) * n_vars_mea_1D + 1:k * n_vars_mea_1D]),
                                                        nb_gdl, nb_mpl)
                                   for k in 1:nb_gc]

        # ... recovery of i_fc and C_O2_Pt.
        i_fc_cell = current(simu.current_density, t_j)
        i_fc = calculate_1D_GC_current_density(i_fc_cell, solver_variables_1D_MEA, simu.fuel_cell)
        for k in 1:nb_gc
            push!(simu.variables["i_fc"][k], i_fc[k])
            push!(simu.variables["C_O2_Pt"][k], calculate_C_O2_Pt(i_fc[k], solver_variables_1D_MEA[k], simu.fuel_cell))
        end

        # ... recovery of Ucell, v_a, v_c, Pa_in and Pc_in.
        push!(simu.variables["Ucell"], calculate_cell_voltage(i_fc[1], simu.variables["C_O2_Pt"][1][end],
                                                                solver_variables_1D_MEA[1], simu.fuel_cell))
        v_a, v_c, Pa_in, Pc_in = calculate_velocity_evolution(solver_variables_1D_MEA, i_fc_cell, simu.fuel_cell, simu.cfg)
        for k in 1:nb_gc
            push!(simu.variables["v_a"][k], v_a[k])
            push!(simu.variables["v_c"][k], v_c[k])
        end
        push!(simu.variables["Pa_in"], Pa_in)
        push!(simu.variables["Pc_in"], Pc_in)
    end
    return nothing
end


"""Display the plots associated with the current simulation.

# Arguments
- `simu::AlphaPEM`: Fuel-cell simulator instance.
- `ax1`: Axes for the first set of plots. Defaults to `nothing`.
- `ax2`: Axes for the second set of plots. Defaults to `nothing`.
- `ax3`: Axes for the third set of plots. Defaults to `nothing`.

# Returns
- `Nothing`
"""
function Display(simu::AlphaPEM, ax1=nothing, ax2=nothing, ax3=nothing)
    # Folder name.
    subfolder_name = String(split(String(simu.cfg.type_fuel_cell), '_')[1])

    # Display.
    if simu.cfg.type_current isa StepParams
        if simu.cfg.type_display == :multiple
            figs_axes = [plt.subplots(figsize=(8, 8)) for _ in 1:11]
            figs = [fa[1] for fa in figs_axes]
            axes = [fa[2] for fa in figs_axes]
            plot_ifc_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[1])
            plot_C_v_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[2])
            plot_lambda_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[3])
            plot_s_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[4])
            plot_C_O2_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[5])
            plot_C_H2_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[6])
            plot_C_N2_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[7])
            plot_T_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[8])
            plot_Ucell(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[9])
            plot_P_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[10])
            plot_v_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, axes[11])

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
        elseif simu.cfg.type_display == :synthetic
            plot_ifc_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[1, 1])
            plot_Ucell(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[1, 2])
            plot_T_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[1, 3])

            plot_C_v_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[2, 1])
            plot_s_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[2, 2])
            plot_lambda_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[2, 3])

            plot_C_H2_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[3, 1])
            plot_C_O2_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[3, 2])
            plot_P_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[3, 3])

            if simu.cfg.type_plot == :fixed
                plot_T_pseudo_2D_final(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax2)
            end

            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(1.0)
        end
    elseif simu.cfg.type_current isa PolarizationParams
        if simu.cfg.type_display == :multiple
            plot_polarisation_curve(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[1])
            plot_power_density_curve(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, length(simu.variables["t"]), ax1[2])
            plot_cell_efficiency(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, length(simu.variables["t"]), ax1[3])

            plot_lambda_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax2[2])
            plot_s_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax2[3])
            plot_T_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax2[4])
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif simu.cfg.type_display == :synthetic
            plot_polarisation_curve(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1)
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif simu.cfg.type_display == :no_display
            plot_polarisation_curve(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1, false)
        end
    elseif simu.cfg.type_current isa PolarizationCalibrationParams
        if simu.cfg.type_display == :multiple
            plot_polarisation_curve_for_cali(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[1])
            plot_lambda_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[2])
            plot_s_1D_temporal(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1[3])
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif simu.cfg.type_display == :synthetic
            plot_polarisation_curve_for_cali(simu.variables, simu.fuel_cell, simu.current_density, simu.cfg, ax1)
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        end
    elseif simu.cfg.type_current isa EISParams
        Fourier_results = make_Fourier_transformation(simu.variables, simu.current_density, simu.cfg)
        if simu.cfg.type_display == :multiple
            plot_EIS_curve_Nyquist(simu.fuel_cell, simu.current_density, simu.cfg, Fourier_results, ax1)
            plot_EIS_curve_Bode_amplitude(simu.fuel_cell, simu.current_density, simu.cfg, Fourier_results, ax2)
            plot_EIS_curve_Bode_angle(simu.fuel_cell, simu.current_density, simu.cfg, Fourier_results, ax3)
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        elseif simu.cfg.type_display == :synthetic
            plot_EIS_curve_Nyquist(simu.fuel_cell, simu.current_density, simu.cfg, Fourier_results, ax1[1])
            plot_EIS_curve_Bode_amplitude(simu.fuel_cell, simu.current_density, simu.cfg, Fourier_results, ax1[2])
            plot_EIS_curve_Bode_angle(simu.fuel_cell, simu.current_density, simu.cfg, Fourier_results, ax1[3])
            # A break is necessary to plot the new points in dynamic mode.
            plt.pause(0.1)
        end
    end
    return nothing
end


"""Save the plots generated for the current simulation.

The output filenames depend on the current profile and on the selected display mode.

# Arguments
- `simu::AlphaPEM`: Fuel-cell simulator instance.
- `fig1`: Figure for the first plot. Defaults to `nothing`.
- `fig2`: Figure for the second plot. Defaults to `nothing`.
- `fig3`: Figure for the third plot. Defaults to `nothing`.

# Returns
- `Nothing`
"""
function Save_plot(simu::AlphaPEM, fig1=nothing, fig2=nothing, fig3=nothing)
    # Folder name.
    subfolder_name = String(split(String(simu.cfg.type_fuel_cell), '_')[1])

    # For the step current.
    if simu.cfg.type_current isa StepParams
        if simu.cfg.type_display == :synthetic
            Saving_instructions(simu, "results", subfolder_name, "step_current_syn_1.pdf", fig1)
            simu.cfg.type_plot == :fixed && Saving_instructions(simu, "results", subfolder_name, "final_temperature_dist_1.pdf", fig2)
        end
    # For the polarization curve.
    elseif simu.cfg.type_current isa PolarizationParams
        if simu.cfg.type_display == :multiple
            Saving_instructions(simu, "results", subfolder_name, "global_indicators_1.pdf", fig1)
            Saving_instructions(simu, "results", subfolder_name, "pola_curve_syn_1.pdf", fig2)
        elseif simu.cfg.type_display == :synthetic
            Saving_instructions(simu, "results", subfolder_name, "pola_curve_1.pdf", fig1)
        end
    # For the EIS curve.
    elseif simu.cfg.type_current isa EISParams
        if simu.cfg.type_display == :multiple
            Saving_instructions(simu, "results", subfolder_name, "Nyquist_plot_1.pdf", fig1)
            Saving_instructions(simu, "results", subfolder_name, "Bode_amplitude_curve_1.pdf", fig2)
            Saving_instructions(simu, "results", subfolder_name, "Bode_angle_curve_1.pdf", fig3)
        elseif simu.cfg.type_display == :synthetic
            Saving_instructions(simu, "results", subfolder_name, "Nyquist_plot_syn_1.pdf", fig1)
        end
    # For the polarization curve for calibration.
    elseif simu.cfg.type_current isa PolarizationCalibrationParams
        if simu.cfg.type_display == :multiple
            Saving_instructions(simu, "results", subfolder_name, "impact_cali_on_internal_state_1.pdf", fig1)
        elseif simu.cfg.type_display == :synthetic
            Saving_instructions(simu, "results", subfolder_name, "pola_curve_cali_1.pdf", fig1)
        end
    end
    return nothing
end


"""Save a figure to the project results directory.

# Arguments
- `simu::AlphaPEM`: Fuel-cell simulator instance.
- `root_folder::String`: Root folder for saving.
- `subfolder_name::String`: Subfolder name for saving.
- `filename::String`: Target filename.
- `fig`: Figure object to save.

# Returns
- `Nothing`
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

