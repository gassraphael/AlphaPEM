# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model. It is a master file which converts
the 1D+1D+1D model to several 1D models in order to ease the coding.
"""

# ______________________Objective function to solve. It gives the system of differential equations______________________

# Typed packing/unpacking and derivative consistency helpers are defined in
# `src/alphapem/core/modules/dif_eq_modules.jl`.


"""In-place DAE residual for IDA (SciML iip=true convention: F!(res, dy, y, p, t)).

The solver vectors `y` and `dydt` are dimensionless/scaled. The residual is built as:
- differential block: `res_diff = dydt_IDA - dydt_model`
- algebraic block: residuals of current-distribution and inlet-flow constraints.

Parameters
----------
res : Vector{Float64}
    Pre-allocated residual vector written in place by the DAE residual function.
    IDA drives this vector to zero (`F(t, y, dy) = 0`).
dydt_IDA : Vector{Float64}
    Time-derivative vector `dy/dt` proposed by IDA at the current Newton
    iteration (scaled solver coordinates).
y : Vector{Float64}
    Current state vector provided by the solver (scaled solver coordinates),
    containing both differential and algebraic variables.
n_vars_cell_1D : Int
    Number of state variables for a single gas-channel (GC) column (MEA + GC).
n_vars_manifold : Int
    Number of manifold state variables (0 when no auxiliary system is present).
n_vars_auxiliary : Int
    Number of auxiliary-system state variables (0 when no auxiliary system is present).
solver_state_scaling : Vector{Float64}
    Reference magnitudes used to non-dimensionalise `y` and `dydt_IDA`.
    Built once by `build_solver_state_scaling` and reused at every residual call.
y_phys_work : Vector{Float64}
    Pre-allocated work buffer that receives the un-scaled physical state
    `y_phys = y .* solver_state_scaling`.  Avoids a heap allocation per call.
flows_work : Vector{MEAFlowsWorkspace}
    One pre-allocated workspace per GC node.  Passed to `calculate_flows_1D_MEA!`
    so that the in-place variant can reuse its internal arrays without allocating
    new memory at each residual evaluation.
heat_work : MEAHeatWorkspace
    Single pre-allocated workspace shared across all GC nodes for
    `calculate_heat_transfers!`.  Heat flows for each node are computed
    sequentially, so one workspace is sufficient.
current_res_scales : Vector{Float64}
    Reference magnitudes used to non-dimensionalise the current-distribution
    residuals.  Built once by `build_current_residual_scaling` and reused at every residual call.
j_in_scale : Float64
    Reference magnitude used to non-dimensionalise the inlet-flow residuals.
    Built once by `build_j_in_residual_scaling` and reused at every residual call.
"""
function dae_residual!(res::Vector{Float64}, dydt_IDA::Vector{Float64}, y::Vector{Float64}, t::Float64,
                       fc::AbstractFuelCell, cd::AbstractCurrent, cfg::SimulationConfig,
                       n_vars_cell_1D::Int, n_vars_manifold::Int,
                       n_vars_auxiliary::Int, solver_state_scaling::Vector{Float64},
                       y_phys_work::Vector{Float64},
                       flows_work::Vector{MEAFlowsWorkspace},
                       heat_work::MEAHeatWorkspace,
                       current_res_scales::Vector{Float64},
                       j_in_scale::Float64)

    # Extraction of frequently used parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    T_des = oc.T_des
    nb_gc, nb_gdl, nb_mpl = np.nb_gc, np.nb_gdl, np.nb_mpl
    type_auxiliary = cfg.type_auxiliary
    has_auxiliary = type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                                       :forced_convective_cathode_with_anodic_recirculation)

    # Build typed local state and derivative containers for each GC node, manifold, and auxiliary system.
    n_diff = nb_gc * n_vars_cell_1D + n_vars_manifold + n_vars_auxiliary
    n_alg = _nb_solver_vars_algebraic(nb_gc)
    expected_len = n_diff + n_alg
    length(y) == expected_len ||
        throw(ArgumentError("Unexpected solver vector size in dae_residual!."))
    length(dydt_IDA) == expected_len ||
        throw(ArgumentError("Unexpected derivative vector size in dae_residual!."))
    length(res) == expected_len ||
        throw(ArgumentError("Unexpected residual vector size in dae_residual!."))
    length(solver_state_scaling) == expected_len ||
        throw(ArgumentError("Unexpected solver scaling vector size in dae_residual!."))
    length(y_phys_work) == expected_len ||
        throw(ArgumentError("Unexpected physical-state work buffer size in dae_residual!."))
    length(flows_work) == nb_gc ||
        throw(ArgumentError("Unexpected MEA flow workspace count in dae_residual!."))

    # Convert the solver state back to physical units before evaluating the
    # transport, electrochemical and thermal equations.
    unscale_values!(y_phys_work, y, solver_state_scaling)
    y_phys = y_phys_work
    cell_state_type = typeof(_unpack_cell_state_1D(@view(y_phys[1:n_vars_cell_1D]), nb_gdl, nb_mpl))
    sv_cell_1D = Vector{cell_state_type}(undef, nb_gc)
    @inbounds for i in 1:nb_gc
        i_start = (i - 1) * n_vars_cell_1D + 1
        i_stop = i * n_vars_cell_1D
        sv_cell_1D[i] = _unpack_cell_state_1D(@view(y_phys[i_start:i_stop]), nb_gdl, nb_mpl)
    end
    dif_eq_cell_1D = Vector{typeof(_nan_cell_derivative_1D(nb_gdl, nb_mpl))}(undef, nb_gc)

    if has_auxiliary
        manifold_offset = nb_gc * n_vars_cell_1D
        aux_offset = manifold_offset + n_vars_manifold
        sv_manifold_1D = _unpack_manifold_state(@view(y_phys[manifold_offset + 1:aux_offset]), np.nb_man)
        sv_auxiliary = _unpack_auxiliary_state(@view(y_phys[aux_offset + 1:n_diff]))
        dif_eq_manifold_1D = _nan_manifold_derivative_state(np.nb_man)
        dif_eq_auxiliary = _nan_auxiliary_derivative()
    else
        sv_manifold_1D = nothing
        sv_auxiliary = nothing
        dif_eq_manifold_1D = nothing
        dif_eq_auxiliary = nothing
    end

    # Algebraic DAE states (canonical ordering):
    # [U_cell, i_fc[1:nb_gc], C_O2_Pt[1:nb_gc], J_a_in, J_c_in]
    alg_offset = n_diff
    U_cell = y_phys[alg_offset + 1]
    i_fc = @view(y_phys[alg_offset + 2:alg_offset + 1 + nb_gc])
    C_O2_Pt = @view(y_phys[alg_offset + 2 + nb_gc:alg_offset + 1 + 2 * nb_gc])
    J_a_in = y_phys[alg_offset + 2 * nb_gc + 2]
    J_c_in = y_phys[alg_offset + 2 * nb_gc + 3]

    # Algebraic block A: current-distribution residuals in physical units.
    i_fc_cell = current(cd, t)
    # Length is 2*nb_gc + 1 for [U_cell + i_fc + C_O2_Pt].
    # We use start + length - 1 explicitly to avoid off-by-one mistakes.
    alg_current_start = alg_offset + 1
    alg_current_len = 2 * nb_gc + 1
    alg_current_end = alg_current_start + alg_current_len - 1
    res_current = @view(res[alg_current_start:alg_current_end])
    gc_current_distribution_residuals!(res_current, U_cell, i_fc, C_O2_Pt, i_fc_cell, sv_cell_1D, fc, cfg)
    res_current ./= current_res_scales

    # Algebraic block B: inlet-flow residuals in physical units.
    alg_flow_start = alg_current_end + 1
    alg_flow_end = alg_flow_start + 1
    res_flow = @view(res[alg_flow_start:alg_flow_end])
    velocity_inlet_flow_residuals!(res_flow, J_a_in, J_c_in, sv_cell_1D, i_fc_cell, fc, cfg)
    res_flow ./= j_in_scale

    # Rebuild flow-dependent fields from algebraic states.
    v_a, v_c, Pa_in, Pc_in = velocity_profiles_from_inlet_flows(sv_cell_1D, J_a_in, J_c_in, fc, cfg)

    # Conditions to pursue the calculations
    for i in 1:nb_gc
        if sv_cell_1D[i].ccl.eta_c > E0
            throw(ArgumentError("The cathode overpotential is higher than the open circuit voltage at time t = " *
                                string(t) * " s. It means that the voltage is negative, which is not possible."))
        end
    end

    # Calculation of the flows for each GC node.
    # The in-place variant reuses the pre-allocated workspace to avoid allocations inside the solver loop.
    first_flow = calculate_flows_1D_MEA!(flows_work[1], sv_cell_1D[1], i_fc[1], v_a[1], v_c[1], fc, cfg)
    flows_1D_MEA = Vector{typeof(first_flow)}(undef, nb_gc)
    flows_1D_MEA[1] = first_flow
    @inbounds for i in 2:nb_gc
        flows_1D_MEA[i] = calculate_flows_1D_MEA!(flows_work[i], sv_cell_1D[i], i_fc[i], v_a[i], v_c[i], fc, cfg)
    end
    flows_1D_GC_manifold = calculate_flows_1D_GC_manifold(sv_cell_1D, sv_auxiliary, i_fc_cell,
                                                           v_a, v_c, Pa_in, Pc_in, fc, cfg)
    # Calculation of the dynamic evolutions inside the MEA.
    @inbounds for i in 1:nb_gc
        dif_eq_int_values_i = calculate_dif_eq_int_values(t, sv_cell_1D[i], fc, cfg, sv_manifold_1D, sv_auxiliary)
        # heat_work is shared across iterations (sequential loop) — one workspace is sufficient.
        heat_flows_i = calculate_heat_transfers!(heat_work, sv_cell_1D[i], i_fc[i], fc, cfg, flows_1D_MEA[i].S_abs,
                                                 flows_1D_MEA[i].Sl)

        dif_eq_mea_diss_water_i = calculate_dyn_dissoved_water_evolution_inside_MEA(sv_cell_1D[i], pp,
                                                                                      flows_1D_MEA[i].S_abs,
                                                                                      flows_1D_MEA[i].J_lambda,
                                                                                      flows_1D_MEA[i].Sp)
        dif_eq_mea_liq_water_i = calculate_dyn_liquid_water_evolution_inside_MEA(sv_cell_1D[i],
                                                                                   pp,
                                                                                   flows_1D_MEA[i].Jl,
                                                                                   flows_1D_MEA[i].S_abs,
                                                                                   flows_1D_MEA[i].Sl)
        dif_eq_mea_vapor_water_i = calculate_dyn_vapor_evolution_inside_MEA(sv_cell_1D[i],
                                                                             pp,
                                                                             flows_1D_MEA[i].Jv,
                                                                             flows_1D_MEA[i].Sv,
                                                                             flows_1D_MEA[i].S_abs)
        dif_eq_mea_species_i = calculate_dyn_H2_O2_N2_evolution_inside_MEA(sv_cell_1D[i],
                                                                            pp,
                                                                            flows_1D_MEA[i].J_H2,
                                                                            flows_1D_MEA[i].J_O2,
                                                                            flows_1D_MEA[i].S_H2,
                                                                            flows_1D_MEA[i].S_O2)
        dif_eq_voltage_i = calculate_dyn_voltage_evolution(i_fc[i], C_O2_Pt[i],
                                                            sv_cell_1D[i].ccl.T,
                                                            sv_cell_1D[i].ccl.eta_c,
                                                            pp,
                                                            dif_eq_int_values_i.i_n)
        dif_eq_mea_temperature_i = calculate_dyn_temperature_evolution_inside_MEA(dif_eq_int_values_i.rho_Cp0,
                                                                                   pp,
                                                                                   heat_flows_i.Jt,
                                                                                   heat_flows_i.Q_r,
                                                                                   heat_flows_i.Q_sorp,
                                                                                   heat_flows_i.Q_liq,
                                                                                   heat_flows_i.Q_p,
                                                                                   heat_flows_i.Q_e)

        dif_eq_cell_1D[i] = assemble_mea_derivative_1D(dif_eq_mea_diss_water_i,
                                                       dif_eq_mea_liq_water_i,
                                                       dif_eq_mea_vapor_water_i,
                                                       dif_eq_mea_species_i,
                                                       dif_eq_voltage_i,
                                                       dif_eq_mea_temperature_i)
    end

    #       Inside the gas channels: compute independent GC contributions, then assemble once.
    dif_eq_gc_gas = calculate_dyn_gas_evolution_inside_gas_channel(sv_cell_1D,
                                                                    pp,
                                                                    cfg,
                                                                    flows_1D_GC_manifold,
                                                                    flows_1D_MEA)
    dif_eq_gc_liq = calculate_dyn_liq_evolution_inside_gas_channel(T_des,
                                                                    pp,
                                                                    cfg,
                                                                    flows_1D_GC_manifold,
                                                                    flows_1D_MEA)
    dif_eq_gc_temperature = calculate_dyn_temperature_evolution_inside_gas_channel(nb_gc)
    dif_eq_cell_1D = assemble_gc_derivative_1D(dif_eq_cell_1D,
                                               dif_eq_gc_gas,
                                               dif_eq_gc_liq,
                                               dif_eq_gc_temperature)

    #       Inside the manifolds and auxiliaries
    if has_auxiliary
        calculate_dyn_manifold_pressure_and_humidity_evolution(dif_eq_manifold_1D,
                                                                T_des, type_auxiliary,
                                                                flows_1D_GC_manifold.W,
                                                                flows_1D_GC_manifold.Wv)
        dif_eq_auxiliary = calculate_dyn_air_compressor_evolution(dif_eq_auxiliary, sv_auxiliary, cfg)
        dif_eq_auxiliary = calculate_dyn_humidifier_evolution(dif_eq_auxiliary, sv_auxiliary, cfg)
        dif_eq_auxiliary = calculate_dyn_throttle_area_controler(dif_eq_auxiliary, sv_auxiliary, cfg)
    end

    # Differential block residual: F_diff = dydt_IDA - dydt_model
    # IDA expects F(t, y, dy) = 0.
    # For differential states, this is exactly dy - f(y) = 0.
    # Pack model derivative into a work view (no allocation).
    dydt_model = @view(res[1:n_diff])
    fuelcell_derivative = FuelCellDerivativeP2D{nb_gdl, nb_mpl, nb_gc}(Tuple(dif_eq_cell_1D))
    _assert_fuelcell_derivative_complete!(fuelcell_derivative)
    _pack_fuelcell_derivative_p2d!(dydt_model, fuelcell_derivative, n_vars_cell_1D)

    if has_auxiliary
        manifold_offset = nb_gc * n_vars_cell_1D
        _assert_manifold_derivative_complete(dif_eq_manifold_1D)
        _pack_manifold_derivative_state!(dydt_model, manifold_offset + 1, dif_eq_manifold_1D)
        aux_offset = manifold_offset + n_vars_manifold
        _assert_auxiliary_derivative_complete(dif_eq_auxiliary)
        _pack_auxiliary_derivative!(dydt_model, aux_offset + 1, dif_eq_auxiliary)
    end

    # Convert physical f(y) to scaled solver-space derivative.
    dydt_model ./= @view(solver_state_scaling[1:n_diff])

    # Final differential residual.
    res_diff = @view(res[1:n_diff])
    @views res_diff .= dydt_IDA[1:n_diff] .- dydt_model

    return nothing
end

