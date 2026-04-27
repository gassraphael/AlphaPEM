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
using SparseArrays: sparse, SparseMatrixCSC, rowvals, nzrange, nonzeros
using SciMLBase: DAEFunction, successful_retcode
using Sundials: IDA
using .PlotHelpers: _clear_dynamic_axes!, saving_instructions!


# _______________________________________________________AlphaPEM_______________________________________________________


mutable struct AlphaPEM
    fuel_cell::AbstractFuelCell
    current_density::AbstractCurrent
    cfg::SimulationConfig
    time_interval::Tuple{Float64, Float64}
    initial_variable_values::Vector{Float64}
    outputs::Union{Nothing, SimulationOutputs}
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
        (0.0, 0.0), # time_interval::Tuple{Float64, Float64}
        Float64[], # initial_variable_values::Vector{Float64}
        nothing, # outputs::Union{Nothing, SimulationOutputs}
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
     Initial values of the physical solver variables. If `nothing`, values are
     generated from a no-current equilibrium. These values remain stored in
     physical units on `simu.initial_variable_values`, then are internally
     scaled before the DAE problem is solved.
time_interval : Union{Nothing, Tuple{Float64, Float64}}, optional
    Time interval for numerical resolution. If `nothing`, it is generated
    according to the chosen current profile.

Returns
-------
Nothing
    The function updates `simu` in place.
"""
function simulate_model!(simu::AlphaPEM,
                         initial_variable_values::Union{Nothing, AbstractVector{<:Real}}=nothing,
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

    # Flow configuration warning for ZSW fuel cells with co_flow.
    if contains(string(typeof(simu.fuel_cell)), "ZSWFuelCell") && simu.cfg.type_flow == :co_flow
        @warn "ZSWFuelCell with standard operating conditions typically requires counter-flow configuration for optimal performance. " *
              "Consider setting type_flow = :counter_flow in SimulationConfig."
    end

    if simu.fuel_cell.operating_conditions.Pa_des < Pext || simu.fuel_cell.operating_conditions.Pc_des < Pext
        throw(ArgumentError("The desired pressure is too low. It cannot be lower than the pressure outside the stack."))
    end

    # Initialize the outputs.
    simu.outputs = nothing

    # Create the dynamic evolution.
    #       Create time intervals
    simu.time_interval = time_interval === nothing ? simu.current_density.time_interval : time_interval
    #       Create the initial variable values
    simu.initial_variable_values = initial_variable_values === nothing ?
                                    create_initial_variable_values(simu) : initial_variable_values

     #       Solve the differential-algebraic equation system.
     #           Pre-calculate constant solver vector dimensions to avoid recomputation in dae_residual!.
     np = simu.cfg.numerical_parameters
     nb_gdl, nb_mpl, nb_gc, nb_man = np.nb_gdl, np.nb_mpl, np.nb_gc, np.nb_man
     n_vars_cell_1D = _nb_solver_vars_cell_1D(nb_gdl, nb_mpl)
     has_auxiliary = simu.cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                                                  :forced_convective_cathode_with_anodic_recirculation)
     n_vars_manifold = has_auxiliary ? _nb_solver_vars_manifolds(nb_man) : 0
     n_vars_auxiliary = _nb_solver_vars_auxiliary(simu.cfg.type_auxiliary)
     differential_vars = build_solver_differential_vars(nb_gc, nb_gdl, nb_mpl, nb_man,
                                                        simu.cfg.type_auxiliary;
                                                        include_algebraic=true)

     simu.initial_variable_values = _ensure_dae_initial_values!(simu, simu.initial_variable_values,
                                                                n_vars_cell_1D, n_vars_manifold,
                                                                n_vars_auxiliary)

     #           Build the fixed internal scaling vector and scale the initial state.
     #           IDA only sees the dimensionless scaled state vector,
     #           while the rest of the code keeps physical units.
     solver_state_scaling = build_solver_state_scaling(simu; include_algebraic=true)
     length(solver_state_scaling) == length(simu.initial_variable_values) ||
         throw(ArgumentError("Internal solver scaling size mismatch in simulate_model!."))
     initial_solver_values = scale_values(simu.initial_variable_values, solver_state_scaling)
     atol_scaled = np.atol ./ solver_state_scaling # preserves the same physical absolute precision after scaling.
     #           Pack external data passed to the DAE residual.
     packed = (fuel_cell=simu.fuel_cell, current_density=simu.current_density, cfg=simu.cfg,
               n_vars_cell_1D=n_vars_cell_1D, n_vars_manifold=n_vars_manifold,
               n_vars_auxiliary=n_vars_auxiliary, solver_state_scaling=solver_state_scaling,
               differential_vars=differential_vars)
     #           Define DAE residual in SciML iip=true signature: F!(res, dydt, y, p, t) -> nothing.
     #           The pre-allocated buffers are managed by the solver — zero output allocation per call.
     residual! = (res, dydt_IDA, y, p, t) -> dae_residual!(res, dydt_IDA, y, t, p.fuel_cell, p.current_density,
                                                            p.cfg, p.n_vars_cell_1D, p.n_vars_manifold,
                                                            p.n_vars_auxiliary, p.solver_state_scaling)

     #           Build a consistent initial derivative vector for differential rows.
     initial_solver_derivatives = _build_consistent_initial_solver_derivatives(residual!, packed,
                                                                                initial_solver_values,
                                                                                simu.time_interval[1],
                                                                                differential_vars)

     #           Build a sparse Jacobian structure for IDA/KLU.
     #           SciML then uses this pattern to compute sparse Jacobians robustly.
     jac_prototype = _build_dae_jacobian_prototype(residual!, packed,
                                                    initial_solver_derivatives,
                                                    initial_solver_values,
                                                    simu.time_interval[1],
                                                    differential_vars)
     jacobian! = (J, dydt_IDA, y, p, gamma, t) -> _dae_jacobian_fd!(J, dydt_IDA, y, p, gamma, t,
                                                                    residual!, differential_vars)
      #           Build and solve the DAE problem with IDA (Sundials).
     dae_fun = DAEFunction(residual!; jac=jacobian!, jac_prototype=jac_prototype)
     prob = DAEProblem(dae_fun, initial_solver_derivatives, initial_solver_values, simu.time_interval, packed;
                       differential_vars=differential_vars)
     simu.sol = solve(prob, IDA(linear_solver=:KLU);
                      reltol=np.rtol, abstol=atol_scaled,
                      initializealg=NoInit())
     #           Check that the solver converged successfully; raise an informative error otherwise.
     successful_retcode(simu.sol.retcode) ||
     throw(ErrorException("IDA solve failed in simulate_model!: retcode = $(simu.sol.retcode)"))

    #       Recover the variable values calculated by the solver into outputs.
    recovery!(simu)
    return nothing
end


"""Build the default initial state vector for the DAE solver.

The initial state represents a no-flow thermodynamic equilibrium: the cell has been
maintained at the desired operating temperature, pressure, and relative humidity long
enough for full relaxation, with no liquid water (s = 0) and no load current.

# Arguments
- `simu::AlphaPEM`: Fuel-cell simulator instance.

# Returns
- `Vector{Float64}`: Initial values of the differential solver variables in canonical ordering.
"""
function create_initial_variable_values(simu::AlphaPEM)::Vector{Float64}
    # Extraction of the parameter classes for better readability.
    oc = simu.fuel_cell.operating_conditions
    pp = simu.fuel_cell.physical_parameters
    np = simu.cfg.numerical_parameters
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

    Psat_ini = Psat(T_ini)
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
    # Ionomer water content: each CL equilibrates with its own local humidity at rest.
    lambda_acl_ini = lambda_eq(C_v_a_ini, s_ini, T_ini)   # Anode CL ↔ anode humidity
    lambda_ccl_ini = lambda_eq(C_v_c_ini, s_ini, T_ini)   # Cathode CL ↔ cathode humidity
    lambda_mem_ini = (lambda_acl_ini + lambda_ccl_ini) / 2 # Membrane: linear-profile average
    # Cathode overpotential: quasi-static value consistent with the initial current.
    i_fc_ini = current(simu.current_density, simu.time_interval[1])
    #       Crossover current density (H2 and O2 membrane permeation).
    i_n_ini = 2 * F * R * T_ini / Hmem * C_H2_ini * k_H2(lambda_mem_ini, T_ini, kappa_co) +
              4 * F * R * T_ini / Hmem * C_O2_ini * k_O2(lambda_mem_ini, T_ini, kappa_co)
    #       Compute the true oxygen concentration at Pt directly from local CCL variables.
    C_O2_Pt_ini = calculate_C_O2_Pt(i_fc_ini, s_ini, lambda_ccl_ini, C_O2_ini, T_ini, simu.fuel_cell)
    #       Butler-Volmer quasi-static inversion: solve dη_c/dt = 0.
    eta_c_ini = R * T_ini / (alpha_c * F) * log((i_fc_ini + i_n_ini) / i0_c_ref *
                                                 1 / exp(-Eact_O2_red / (R * T_ini) * (1 / T_ini - 1 / Tref_O2_red)) *
                                                 (C_O2ref_red / C_O2_Pt_ini)^kappa_c)

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
    lambda_acl, lambda_mem, lambda_ccl = lambda_acl_ini, lambda_mem_ini, lambda_ccl_ini
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
    names_1D = canonical_cell_solver_variable_names_1D(nb_gdl, nb_mpl)
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


"""Append missing DAE algebraic initial values in physical units when needed."""
function _ensure_dae_initial_values!(simu::AlphaPEM,
                                     initial_values::Vector{Float64},
                                     n_vars_cell_1D::Int,
                                     n_vars_manifold::Int,
                                     n_vars_auxiliary::Int)::Vector{Float64}
    np = simu.cfg.numerical_parameters
    pp = simu.fuel_cell.physical_parameters
    nb_gc = np.nb_gc
    n_diff = nb_gc * n_vars_cell_1D + n_vars_manifold + n_vars_auxiliary
    n_alg = _nb_solver_vars_algebraic(nb_gc)

    length(initial_values) == n_diff + n_alg && return initial_values
    length(initial_values) == n_diff ||
        throw(ArgumentError("Initial state size mismatch in _ensure_dae_initial_values!."))

    # Build typed differential states from the physical initialization.
    sv_cell_1D = [_unpack_cell_state_1D(@view(initial_values[(k - 1) * n_vars_cell_1D + 1:k * n_vars_cell_1D]),
                                        np.nb_gdl, np.nb_mpl)
                  for k in 1:nb_gc]
    t0 = simu.time_interval[1]
    i_fc_cell_0 = current(simu.current_density, t0)

    # Initialize current-distribution algebraic states with existing robust kernels.
    i_fc_0 = calculate_1D_GC_current_density(i_fc_cell_0, sv_cell_1D, simu.cfg, simu.fuel_cell)
    C_O2_Pt_0 = [calculate_C_O2_Pt(i_fc_0[k], sv_cell_1D[k], simu.fuel_cell) for k in 1:nb_gc]
    U_cell_0 = calculate_cell_voltage(i_fc_0[1], C_O2_Pt_0[1], sv_cell_1D[1], simu.fuel_cell)

    # Initialize inlet-flow algebraic states via pressure-consistent desired flows.
    _, _, Pa_in_0, Pc_in_0 = calculate_velocity_evolution(sv_cell_1D, i_fc_cell_0, simu.fuel_cell, simu.cfg)
    W_des_0 = desired_flows(sv_cell_1D, i_fc_cell_0, Pa_in_0, Pc_in_0, simu.fuel_cell, simu.cfg)
    J_a_in_0 = (W_des_0.H2 + W_des_0.H2O_inj_a) / (pp.Hagc * pp.Wagc) / pp.nb_cell / pp.nb_channel_in_gc
    J_c_in_0 = (W_des_0.dry_air + W_des_0.H2O_inj_c) / (pp.Hcgc * pp.Wcgc) / pp.nb_cell / pp.nb_channel_in_gc

    append!(initial_values, [U_cell_0])
    append!(initial_values, i_fc_0)
    append!(initial_values, C_O2_Pt_0)
    append!(initial_values, [J_a_in_0, J_c_in_0])
    return initial_values
end


"""Populate `simu.outputs` from the solver output.

# Arguments
- `simu::AlphaPEM`: Fuel-cell simulator instance.

# Returns
- `Nothing`: The function updates `simu.outputs` in place.
"""
function recovery!(simu::AlphaPEM)
    # Recovery of the time span.
    t_hist = collect(simu.sol.t)

    # Recovery of the main variables dynamic evolution.
    np = simu.cfg.numerical_parameters
    nb_gc, nb_gdl, nb_mpl, nb_man = np.nb_gc, np.nb_gdl, np.nb_mpl, np.nb_man
    n_vars_cell_1D = _nb_solver_vars_cell_1D(nb_gdl, nb_mpl)
    has_auxiliary = simu.cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                                                 :forced_convective_cathode_with_anodic_recirculation)
    n_vars_manifold = has_auxiliary ? _nb_solver_vars_manifolds(nb_man) : 0
    n_vars_auxiliary = _nb_solver_vars_auxiliary(simu.cfg.type_auxiliary)
    n_diff = nb_gc * n_vars_cell_1D + n_vars_manifold + n_vars_auxiliary

    solver_state_scaling = build_solver_state_scaling(simu; include_algebraic=true)
    canonical_names = canonical_cell_solver_variable_names_1D(nb_gdl, nb_mpl)
    length(canonical_names) == n_vars_cell_1D ||
        throw(ArgumentError("Canonical MEA layout size mismatch in recovery!."))
    length(solver_state_scaling) == length(simu.sol.u[1]) ||
        throw(ArgumentError("Internal solver scaling size mismatch in recovery!."))

    # Typed output buffers.
    n_t = length(t_hist)
    solver_states = Vector{FuelCellStateP2D{nb_gdl, nb_mpl, nb_gc}}(undef, n_t)
    i_fc_hist = [Float64[] for _ in 1:nb_gc]
    C_O2_Pt_hist = [Float64[] for _ in 1:nb_gc]
    v_a_hist = [Float64[] for _ in 1:nb_gc]
    v_c_hist = [Float64[] for _ in 1:nb_gc]
    Ucell_hist = Float64[]
    Pa_in_hist = Float64[]
    Pc_in_hist = Float64[]

    for j in eachindex(t_hist)
        y_phys_j = unscale_values(simu.sol.u[j], solver_state_scaling)

        # Recovery of the variables inside the MEA 1D line.
        sv_cell_1D = [_unpack_cell_state_1D(@view(y_phys_j[(k - 1) * n_vars_cell_1D + 1:k * n_vars_cell_1D]),
                                            nb_gdl, nb_mpl)
                      for k in 1:nb_gc]
        solver_states[j] = FuelCellStateP2D{nb_gdl, nb_mpl, nb_gc}(Tuple(sv_cell_1D))

        # Recovery of the DAE algebraic states in canonical ordering:
        # [U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc], J_a_in, J_c_in]
        alg_offset = n_diff
        U_cell = y_phys_j[alg_offset + 1]
        i_fc = @view(y_phys_j[alg_offset + 2:alg_offset + 1 + nb_gc])
        C_O2_Pt = @view(y_phys_j[alg_offset + 2 + nb_gc:alg_offset + 1 + 2 * nb_gc])
        J_a_in = y_phys_j[alg_offset + 2 * nb_gc + 2]
        J_c_in = y_phys_j[alg_offset + 2 * nb_gc + 3]
        # Store these algebraic states for each time step.
        push!(Ucell_hist, U_cell)
        for k in 1:nb_gc
            push!(i_fc_hist[k], i_fc[k])
            push!(C_O2_Pt_hist[k], C_O2_Pt[k])
        end

        # Recover velocity profiles and inlet pressures from the algebraic inlet flow densities.
        v_a, v_c, Pa_in, Pc_in = velocity_profiles_from_inlet_flows(sv_cell_1D, J_a_in, J_c_in,
                                                                     simu.fuel_cell, simu.cfg)
        # Store anode and cathode velocity profiles for each gas channel node.
        for k in 1:nb_gc
            push!(v_a_hist[k], v_a[k])
            push!(v_c_hist[k], v_c[k])
        end
        # Store inlet pressures for anode and cathode sides.
        push!(Pa_in_hist, Pa_in)
        push!(Pc_in_hist, Pc_in)
    end

    simu.outputs = SimulationOutputs{nb_gdl, nb_mpl, nb_gc}(
        SolverTrajectory{nb_gdl, nb_mpl, nb_gc}(t_hist, solver_states),
        DerivedOutputs{nb_gc}(Ucell_hist, i_fc_hist, C_O2_Pt_hist, v_a_hist, v_c_hist, Pa_in_hist, Pc_in_hist),
    )
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

function display!(simu::AlphaPEM, _ax1=nothing, _ax2=nothing, _ax3=nothing)
    outputs = simu.outputs
    outputs === nothing && throw(ArgumentError("display! requires available simulation outputs. Run simulate_model! first."))
    simu.cfg.type_display == :no_display && return nothing
    nb_gc = simu.cfg.numerical_parameters.nb_gc
    has_extended_gc_profiles = nb_gc >= 3
    is_postrun_display = simu.cfg.display_timing == :postrun
    can_plot_extended_gc_postrun = has_extended_gc_profiles && is_postrun_display

    simu.cfg.display_timing == :live && !(simu.cfg.type_current isa EISParams) &&
        _clear_dynamic_axes!(_ax1, _ax2, _ax3)

    if simu.cfg.type_current isa StepParams
        if simu.cfg.type_display == :synthetic && _ax1 !== nothing
            plot_ifc_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[1, 1])
            plot_Ucell(outputs, simu.current_density, simu.cfg, _ax1[1, 2])
            plot_T_1D_temporal(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1[1, 3])

            plot_C_v_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[2, 1])
            plot_s_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[2, 2])
            plot_lambda_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[2, 3])

            plot_C_H2_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[3, 1])
            plot_C_O2_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[3, 2])
            plot_P_1D_temporal(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1[3, 3])

            if (can_plot_extended_gc_postrun &&
                _ax2 isa AbstractMatrix && size(_ax2, 1) >= 2 && size(_ax2, 2) >= 2)
                plot_T_pseudo_2D_final(outputs, simu.fuel_cell, _ax2[1, 1].parent, _ax2[1, 1], simu.cfg)
                plot_ifc_GC_final(outputs, simu.current_density, simu.cfg, _ax2[1, 2])
                plot_C_O2_Pt_GC_final(outputs, simu.cfg, _ax2[2, 1])
                plot_lambda_mem_GC_final(outputs, simu.cfg, _ax2[2, 2])
            end
        elseif simu.cfg.type_display == :multiple && _ax1 isa AbstractVector && length(_ax1) >= 14
            # Multiple mode: one figure per internal-state plot.
            plot_ifc_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[1])
            plot_Ucell(outputs, simu.current_density, simu.cfg, _ax1[2])
            plot_T_1D_temporal(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1[3])
            plot_C_v_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[4])
            plot_s_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[5])
            plot_lambda_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[6])
            plot_C_H2_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[7])
            plot_C_O2_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[8])
            plot_P_1D_temporal(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1[9])
            plot_C_N2_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[10])
            plot_Phi_a_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[11])
            plot_Phi_c_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[12])
            plot_v_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[13])
            plot_Re_nb_1D_temporal(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1[14])

            if can_plot_extended_gc_postrun && _ax2 !== nothing
                plot_T_pseudo_2D_final(outputs, simu.fuel_cell, _ax2.parent, _ax2, simu.cfg)
                if _ax3 isa AbstractVector && length(_ax3) >= 3
                    plot_ifc_GC_final(outputs, simu.current_density, simu.cfg, _ax3[1])
                    plot_C_O2_Pt_GC_final(outputs, simu.cfg, _ax3[2])
                    plot_lambda_mem_GC_final(outputs, simu.cfg, _ax3[3])
                elseif _ax3 !== nothing
                    plot_ifc_GC_final(outputs, simu.current_density, simu.cfg, _ax3)
                end
            end
        end
    elseif simu.cfg.type_current isa PolarizationParams
        if simu.cfg.type_display == :synthetic && _ax1 !== nothing
            plot_polarization_curve(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1)
        elseif simu.cfg.type_display == :multiple && _ax1 isa AbstractVector && length(_ax1) >= 5 && _ax2 !== nothing
            # Multiple mode: internal states and derived polarization curves in individual figures.
            plot_ifc_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[1])
            plot_Ucell(outputs, simu.current_density, simu.cfg, _ax1[2])
            plot_T_1D_temporal(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1[3])
            plot_power_density_curve(outputs, simu.current_density, simu.cfg, _ax1[4])
            plot_cell_efficiency(outputs, simu.current_density, simu.cfg, _ax1[5])
            if !(_ax2 isa AbstractVector)
                plot_polarization_curve(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax2)
            end
        end
    elseif simu.cfg.type_current isa PolarizationCalibrationParams
        if simu.cfg.type_display == :synthetic && _ax1 !== nothing
            plot_polarization_curve_for_cali(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1)
        elseif simu.cfg.type_display == :multiple && _ax1 isa AbstractVector && length(_ax1) >= 5 && _ax2 !== nothing
            plot_ifc_1D_temporal(outputs, simu.current_density, simu.cfg, _ax1[1])
            plot_Ucell(outputs, simu.current_density, simu.cfg, _ax1[2])
            plot_T_1D_temporal(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax1[3])
            plot_power_density_curve(outputs, simu.current_density, simu.cfg, _ax1[4])
            plot_cell_efficiency(outputs, simu.current_density, simu.cfg, _ax1[5])
            if !(_ax2 isa AbstractVector)
                plot_polarization_curve_for_cali(outputs, simu.fuel_cell, simu.current_density, simu.cfg, _ax2)
            end
        end
    elseif simu.cfg.type_current isa EISParams
        # EIS display uses one Fourier point per dynamic segment.
        Fourier_results = make_Fourier_transformation(outputs, simu.current_density, simu.cfg)
        if simu.cfg.type_display == :synthetic && _ax1 !== nothing
            plot_EIS_curve_Nyquist(simu.current_density, Fourier_results, _ax1[1])
            plot_EIS_curve_Bode_amplitude(simu.current_density, Fourier_results, _ax1[2])
            plot_EIS_curve_Bode_angle(simu.current_density, Fourier_results, _ax1[3])
        elseif simu.cfg.type_display == :multiple && _ax1 !== nothing && _ax2 !== nothing && _ax3 !== nothing
            plot_EIS_curve_Nyquist(simu.current_density, Fourier_results, _ax1)
            plot_EIS_curve_Bode_amplitude(simu.current_density, Fourier_results, _ax2)
            plot_EIS_curve_Bode_angle(simu.current_density, Fourier_results, _ax3)
        end
    else
        @warn "display!: plotting migration to CairoMakie is in progress; this simulation/display mode is not ported yet." maxlog=1
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
function save_plot!(simu::AlphaPEM, _fig1=nothing, _fig2=nothing, _fig3=nothing)
    simu.outputs === nothing && throw(ArgumentError("save_plot! requires available simulation outputs. Run simulate_model! first."))
    nb_gc = simu.cfg.numerical_parameters.nb_gc
    has_extended_gc_profiles = nb_gc >= 3
    is_postrun_display = simu.cfg.display_timing == :postrun
    can_save_extended_gc_postrun = has_extended_gc_profiles && is_postrun_display

    # Folder name.
    subfolder_name = String(split(String(simu.cfg.type_fuel_cell), '_')[1])

    # For the step current.
    if simu.cfg.type_current isa StepParams
        if simu.cfg.type_display == :synthetic
            saving_instructions!(simu, "results", subfolder_name, "step_current_syn_1.pdf", _fig1)
            if can_save_extended_gc_postrun
                saving_instructions!(simu, "results", subfolder_name, "final_GC_profiles_1.pdf", _fig2)
            end
        elseif simu.cfg.type_display == :multiple && _fig1 isa AbstractVector
            # Multiple mode: one file per internal-state plot.
            step_files = [
                "ifc_1D_temporal_1.pdf",
                "Ucell_temporal_1.pdf",
                "T_1D_temporal_1.pdf",
                "Cv_1D_temporal_1.pdf",
                "s_1D_temporal_1.pdf",
                "lambda_1D_temporal_1.pdf",
                "CH2_1D_temporal_1.pdf",
                "CO2_1D_temporal_1.pdf",
                "P_1D_temporal_1.pdf",
                "CN2_1D_temporal_1.pdf",
                "Phi_a_1D_temporal_1.pdf",
                "Phi_c_1D_temporal_1.pdf",
                "v_1D_temporal_1.pdf",
                "Re_1D_temporal_1.pdf",
            ]
            for (fig_i, filename) in zip(_fig1, step_files)
                saving_instructions!(simu, "results", subfolder_name, filename, fig_i)
            end
            if can_save_extended_gc_postrun
                saving_instructions!(simu, "results", subfolder_name, "final_temperature_map_and_current_GC_1.pdf", _fig2)
            end
            if can_save_extended_gc_postrun && _fig3 isa AbstractVector && length(_fig3) >= 3
                saving_instructions!(simu, "results", subfolder_name, "ifc_GC_final_1.pdf", _fig3[1])
                saving_instructions!(simu, "results", subfolder_name, "CO2_Pt_GC_final_1.pdf", _fig3[2])
                saving_instructions!(simu, "results", subfolder_name, "lambda_mem_GC_final_1.pdf", _fig3[3])
            elseif can_save_extended_gc_postrun
                saving_instructions!(simu, "results", subfolder_name, "ifc_GC_final_1.pdf", _fig3)
            end
        end
    # For the polarization curve.
    elseif simu.cfg.type_current isa PolarizationParams
        if simu.cfg.type_display == :multiple
            if _fig1 isa AbstractVector && length(_fig1) >= 5
                saving_instructions!(simu, "results", subfolder_name, "ifc_1D_temporal_1.pdf",     _fig1[1])
                saving_instructions!(simu, "results", subfolder_name, "Ucell_temporal_1.pdf",      _fig1[2])
                saving_instructions!(simu, "results", subfolder_name, "T_1D_temporal_1.pdf",       _fig1[3])
                saving_instructions!(simu, "results", subfolder_name, "power_density_curve_1.pdf", _fig1[4])
                saving_instructions!(simu, "results", subfolder_name, "efficiency_curve_1.pdf",    _fig1[5])
            end
            saving_instructions!(simu, "results", subfolder_name, "pola_curve_1.pdf", _fig2)
        elseif simu.cfg.type_display == :synthetic
            saving_instructions!(simu, "results", subfolder_name, "pola_curve_1.pdf", _fig1)
        end
    # For the EIS curve.
    elseif simu.cfg.type_current isa EISParams
        if simu.cfg.type_display == :multiple
            saving_instructions!(simu, "results", subfolder_name, "Nyquist_plot_1.pdf", _fig1)
            saving_instructions!(simu, "results", subfolder_name, "Bode_amplitude_curve_1.pdf", _fig2)
            saving_instructions!(simu, "results", subfolder_name, "Bode_angle_curve_1.pdf", _fig3)
        elseif simu.cfg.type_display == :synthetic
            saving_instructions!(simu, "results", subfolder_name, "Nyquist_plot_syn_1.pdf", _fig1)
        end
    # For the polarization curve for calibration.
    elseif simu.cfg.type_current isa PolarizationCalibrationParams
        if simu.cfg.type_display == :multiple
            if _fig1 isa AbstractVector && length(_fig1) >= 5
                saving_instructions!(simu, "results", subfolder_name, "ifc_1D_temporal_cali_1.pdf",     _fig1[1])
                saving_instructions!(simu, "results", subfolder_name, "Ucell_temporal_cali_1.pdf",      _fig1[2])
                saving_instructions!(simu, "results", subfolder_name, "T_1D_temporal_cali_1.pdf",       _fig1[3])
                saving_instructions!(simu, "results", subfolder_name, "power_density_curve_cali_1.pdf", _fig1[4])
                saving_instructions!(simu, "results", subfolder_name, "efficiency_curve_cali_1.pdf",    _fig1[5])
            end
            saving_instructions!(simu, "results", subfolder_name, "pola_curve_cali_1.pdf", _fig2)
        elseif simu.cfg.type_display == :synthetic
            saving_instructions!(simu, "results", subfolder_name, "pola_curve_cali_1.pdf", _fig1)
        end
    end

    return nothing
end

