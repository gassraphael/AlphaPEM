# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model. It is a master file which converts
the 1D+1D+1D model to several 1D models in order to ease the coding.
"""

# ______________________Objective function to solve. It gives the system of differential equations______________________

# Typed packing/unpacking and derivative consistency helpers are defined in
# `src/alphapem/core/modules/dif_eq_modules.jl`.

"""In-place RHS for the ODE solver (SciML iip=true convention: f!(dy, y, p, t)).

Writes the derivative directly into the pre-allocated vector `dy` managed by the solver,
eliminating all output-vector allocations.

The solver vector `y` is dimensionless and scaled so that most state variables are
of order 1. The physical model itself is still evaluated in physical units:
`y` is first unscaled into a physical state vector, then `dy` is rescaled before
being returned to the ODE solver.

Parameters
----------
dy : Vector{Float64}
    Pre-allocated output vector (managed by the ODE solver). Written in place.
y : Vector{Float64}
    Current solver state vector.
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cd : AbstractCurrent
    Current profile instance (prescribes current as a function of time).
cfg : SimulationConfig
    Simulation configuration (type_auxiliary, etc.).
n_vars_cell_1D : Int
    Pre-calculated number of solver variables per gas-channel node.
n_vars_manifold : Int
    Pre-calculated number of solver variables in manifolds.
n_vars_auxiliary : Int
    Pre-calculated number of solver variables in auxiliary systems.
 * solver_state_scaling : Vector{Float64}
     Scaling factors aligned with the full solver state vector ordering.
"""
function dydt!(dy::Vector{Float64}, t::Float64, y::Vector{Float64}, fc::AbstractFuelCell, cd::AbstractCurrent,
               cfg::SimulationConfig, n_vars_cell_1D::Int, n_vars_manifold::Int,
               n_vars_auxiliary::Int, solver_state_scaling::Vector{Float64})

    # Extraction of frequently used parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des = oc.T_des
    Lgc, Hagc, Hcgc = pp.Lgc, pp.Hagc, pp.Hcgc
    nb_gc, nb_gdl, nb_mpl, nb_man = np.nb_gc, np.nb_gdl, np.nb_mpl, np.nb_man

    # Global model composition
    type_auxiliary = cfg.type_auxiliary
    has_auxiliary = type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                                       :forced_convective_cathode_with_anodic_recirculation)

    # Build typed local state and derivative containers for each GC node, manifold, and auxiliary system.
    n_vars_cell_P2D = nb_gc * n_vars_cell_1D
    expected_len = n_vars_cell_P2D + n_vars_manifold + n_vars_auxiliary
    length(y) == expected_len ||
        throw(ArgumentError("Unexpected solver vector size for typed path."))
    length(solver_state_scaling) == expected_len ||
        throw(ArgumentError("Unexpected solver scaling vector size for typed path."))

    # Convert the solver state back to physical units before evaluating the
    # transport, electrochemical and thermal equations.
    y_phys = unscale_values(y, solver_state_scaling)

    sv_cell_1D = [_unpack_cell_state_1D(@view(y_phys[(i - 1) * n_vars_cell_1D + 1:i * n_vars_cell_1D]),
                                      nb_gdl, nb_mpl)
                 for i in 1:nb_gc]
    dif_eq_cell_1D = Vector{typeof(_nan_cell_derivative_1D(nb_gdl, nb_mpl))}(undef, nb_gc)

    if has_auxiliary
        manifold_offset = n_vars_cell_P2D
        aux_offset = n_vars_cell_P2D + n_vars_manifold
        sv_manifold_1D = _unpack_manifold_state(@view(y_phys[manifold_offset + 1:aux_offset]), nb_man)
        sv_auxiliary = _unpack_auxiliary_state(@view(y_phys[aux_offset + 1:end]))
        dif_eq_manifold_1D = _nan_manifold_derivative_state(nb_man)
        dif_eq_auxiliary = _nan_auxiliary_derivative()
    else
        sv_manifold_1D = nothing
        sv_auxiliary = nothing
        dif_eq_manifold_1D = nothing
        dif_eq_auxiliary = nothing
    end

    # Conditions to pursue the calculations
    for i in 1:nb_gc
        if sv_cell_1D[i].ccl.eta_c > E0
            throw(ArgumentError("The cathode overpotential is higher than the open circuit voltage at time t = " *
                                string(t) * " s. It means that the voltage is negative, which is not possible."))
        end
    end

    # Intermediate values (one container per GC node)
    dif_eq_int_values = [calculate_dif_eq_int_values(t, sv_cell_1D[i], fc, cfg, sv_manifold_1D, sv_auxiliary)
                         for i in 1:nb_gc]

    # Calculate the local current density at each node of the GC.
    i_fc_cell = current(cd, t)
    i_fc = calculate_1D_GC_current_density(i_fc_cell, sv_cell_1D, fc)

    # Calculation of the oxygen concentration at the platinum surface in the cathode catalyst layer
    C_O2_Pt = [calculate_C_O2_Pt(i_fc[i], sv_cell_1D[i], fc) for i in 1:nb_gc]

    # Calculation of the velocities inside the GC and the manifolds
    v_a, v_c, Pa_in, Pc_in = calculate_velocity_evolution(sv_cell_1D, i_fc_cell, fc, cfg)

    # Calculation of the flows for each GC node.
    flows_1D_MEA = [calculate_flows_1D_MEA(sv_cell_1D[i], i_fc[i], v_a[i], v_c[i], fc)
                    for i in 1:nb_gc]
    flows_1D_GC_manifold = calculate_flows_1D_GC_manifold(sv_cell_1D, sv_manifold_1D, sv_auxiliary, i_fc_cell,
                                                           v_a, v_c, Pa_in, Pc_in, fc, cfg)
    heat_flows_global = [calculate_heat_transfers(sv_cell_1D[i], i_fc[i], fc, flows_1D_MEA[i].S_abs,
                                                  flows_1D_MEA[i].Sl)
                         for i in 1:nb_gc]

    # Calculation of the dynamic evolutions inside the MEA.
    dif_eq_mea_diss_water = [calculate_dyn_dissoved_water_evolution_inside_MEA(sv_cell_1D[i], pp,
                                                                               flows_1D_MEA[i].S_abs,
                                                                               flows_1D_MEA[i].J_lambda,
                                                                               flows_1D_MEA[i].Sp)
                             for i in 1:nb_gc]
    dif_eq_mea_liq_water = [calculate_dyn_liquid_water_evolution_inside_MEA(sv_cell_1D[i],
                                                                            pp,
                                                                            flows_1D_MEA[i].Jl,
                                                                            flows_1D_MEA[i].S_abs,
                                                                            flows_1D_MEA[i].Sl)
                            for i in 1:nb_gc]
    dif_eq_mea_vapor_water = [calculate_dyn_vapor_evolution_inside_MEA(sv_cell_1D[i],
                                                                       pp,
                                                                       flows_1D_MEA[i].Jv,
                                                                       flows_1D_MEA[i].Sv,
                                                                       flows_1D_MEA[i].S_abs)
                              for i in 1:nb_gc]
    dif_eq_mea_species = [calculate_dyn_H2_O2_N2_evolution_inside_MEA(sv_cell_1D[i],
                                                                      pp,
                                                                      flows_1D_MEA[i].J_H2,
                                                                      flows_1D_MEA[i].J_O2,
                                                                      flows_1D_MEA[i].S_H2,
                                                                      flows_1D_MEA[i].S_O2)
                          for i in 1:nb_gc]
    dif_eq_voltage = [calculate_dyn_voltage_evolution(i_fc[i], C_O2_Pt[i],
                                                      sv_cell_1D[i].ccl.T,
                                                      sv_cell_1D[i].ccl.eta_c,
                                                      pp,
                                                      dif_eq_int_values[i].i_n)
                      for i in 1:nb_gc]
    dif_eq_mea_temperature = [calculate_dyn_temperature_evolution_inside_MEA(dif_eq_int_values[i].rho_Cp0,
                                                                             pp,
                                                                             heat_flows_global[i].Jt,
                                                                             heat_flows_global[i].Q_r,
                                                                             heat_flows_global[i].Q_sorp,
                                                                             heat_flows_global[i].Q_liq,
                                                                             heat_flows_global[i].Q_p,
                                                                             heat_flows_global[i].Q_e)
                               for i in 1:nb_gc]

    for i in 1:nb_gc
        dif_eq_cell_1D[i] = assemble_mea_derivative_1D(dif_eq_mea_diss_water[i],
                                                      dif_eq_mea_liq_water[i],
                                                      dif_eq_mea_vapor_water[i],
                                                      dif_eq_mea_species[i],
                                                      dif_eq_voltage[i],
                                                      dif_eq_mea_temperature[i])
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

    # Pack typed derivatives directly into the pre-allocated solver buffer dy (no allocation).
    fuelcell_derivative = FuelCellDerivativeP2D{nb_gdl, nb_mpl, nb_gc}(Tuple(dif_eq_cell_1D))
    _assert_fuelcell_derivative_complete!(fuelcell_derivative)
    _pack_fuelcell_derivative_p2d!(dy, fuelcell_derivative, n_vars_cell_1D)

    if has_auxiliary
        manifold_offset = n_vars_cell_P2D
        _assert_manifold_derivative_complete(dif_eq_manifold_1D)
        _pack_manifold_derivative_state!(dy, manifold_offset + 1, dif_eq_manifold_1D)
        aux_offset = n_vars_cell_P2D + n_vars_manifold
        _assert_auxiliary_derivative_complete(dif_eq_auxiliary)
        _pack_auxiliary_derivative!(dy, aux_offset + 1, dif_eq_auxiliary)
    end

    # Convert physical derivatives back to solver-space derivatives.
    dy ./= solver_state_scaling

    return nothing
end


