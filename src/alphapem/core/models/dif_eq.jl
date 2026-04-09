# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model. It is a master file which converts
the 1D+1D+1D model to several 1D models in order to ease the coding.
"""

# ______________________Objective function to solve. It gives the system of differential equations______________________

# Typed packing/unpacking and derivative consistency helpers are defined in
# `src/alphapem/core/modules/dif_eq_modules.jl`.

"""This function gives the system of differential equations to solve.

Parameters
----------
t : Float64
    Time (s).
y : Vector{Float64}
    Vector of the solver variables.
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cd : AbstractCurrent
    Current profile instance (prescribes current as a function of time).
cfg : SimulationConfig
    Simulation configuration (type_auxiliary, etc.).

Returns
-------
Vector{Float64}
    Vector containing the derivative of the solver variables.
"""
function dydt(t::Float64, y::Vector{Float64}, fc::AbstractFuelCell, cd::AbstractCurrent,
              cfg::SimulationConfig)::Vector{Float64}

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

    # Build typed local state and derivative containers for each GC node.
    n_vars_per_gc = _nb_solver_vars_per_gc(nb_gdl, nb_mpl)
    n_vars_mea = nb_gc * n_vars_per_gc
    n_vars_manifold = has_auxiliary ? _nb_solver_vars_manifolds(nb_man) : 0
    n_vars_auxiliary = _nb_solver_vars_auxiliary(type_auxiliary)
    expected_len = n_vars_mea + n_vars_manifold + n_vars_auxiliary
    length(y) == expected_len ||
        throw(ArgumentError("Unexpected solver vector size for typed path."))

    sv_mea_1D = [_unpack_mea_state_1D(@view(y[(i - 1) * n_vars_per_gc + 1:i * n_vars_per_gc]),
                                      nb_gdl, nb_mpl)
                 for i in 1:nb_gc]
    dif_eq_mea_1D = Vector{typeof(_nan_mea_derivative_1D(nb_gdl, nb_mpl))}(undef, nb_gc)

    if has_auxiliary
        manifold_offset = n_vars_mea
        aux_offset = n_vars_mea + n_vars_manifold
        sv_manifold_1D = _unpack_manifold_state(@view(y[manifold_offset + 1:aux_offset]), nb_man)
        sv_auxiliary = _unpack_auxiliary_state(@view(y[aux_offset + 1:end]))
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
        if sv_mea_1D[i].ccl.eta_c > E0
            throw(ArgumentError("The cathode overpotential is higher than the open circuit voltage at time t = " *
                                string(t) * " s. It means that the voltage is negative, which is not possible."))
        end
    end

    # Intermediate values (one container per GC node)
    dif_eq_int_values = [calculate_dif_eq_int_values(t, sv_mea_1D[i], fc, cfg, sv_manifold_1D, sv_auxiliary)
                         for i in 1:nb_gc]

    # Calculate the local current density at each node of the GC.
    i_fc_cell = current(cd, t)
    i_fc = calculate_1D_GC_current_density(i_fc_cell, sv_mea_1D, fc)

    # Calculation of the oxygen concentration at the platinum surface in the cathode catalyst layer
    C_O2_Pt = [calculate_C_O2_Pt(i_fc[i], sv_mea_1D[i], fc) for i in 1:nb_gc]

    # Calculation of the velocities inside the GC and the manifolds
    v_a, v_c, Pa_in, Pc_in = calculate_velocity_evolution(sv_mea_1D, i_fc_cell, fc, cfg)

    # Calculation of the flows for each GC node.
    flows_1D_MEA = [calculate_flows_1D_MEA(sv_mea_1D[i], i_fc[i], v_a[i], v_c[i], fc)
                    for i in 1:nb_gc]
    flows_1D_GC_manifold = calculate_flows_1D_GC_manifold(sv_mea_1D, sv_manifold_1D, sv_auxiliary, i_fc_cell,
                                                           v_a, v_c, Pa_in, Pc_in, fc, cfg)
    heat_flows_global = [calculate_heat_transfers(sv_mea_1D[i], i_fc[i], fc, flows_1D_MEA[i].S_abs,
                                                  flows_1D_MEA[i].Sl)
                         for i in 1:nb_gc]

    # Calculation of the dynamic evolutions inside the MEA.
    # Each calculate_dyn_* builds one physics-specific derivative contribution over all GC nodes.
    dif_eq_mea_diss_water = [calculate_dyn_dissoved_water_evolution_inside_MEA(sv_mea_1D[i], pp,
                                                                               flows_1D_MEA[i].S_abs,
                                                                               flows_1D_MEA[i].J_lambda,
                                                                               flows_1D_MEA[i].Sp)
                             for i in 1:nb_gc]
    dif_eq_mea_liq_water = [calculate_dyn_liquid_water_evolution_inside_MEA(sv_mea_1D[i],
                                                                            pp,
                                                                            flows_1D_MEA[i].Jl,
                                                                            flows_1D_MEA[i].S_abs,
                                                                            flows_1D_MEA[i].Sl)
                            for i in 1:nb_gc]
    dif_eq_mea_vapor_water = [calculate_dyn_vapor_evolution_inside_MEA(sv_mea_1D[i],
                                                                       pp,
                                                                       flows_1D_MEA[i].Jv,
                                                                       flows_1D_MEA[i].Sv,
                                                                       flows_1D_MEA[i].S_abs)
                              for i in 1:nb_gc]
    dif_eq_mea_species = [calculate_dyn_H2_O2_N2_evolution_inside_MEA(sv_mea_1D[i],
                                                                      pp,
                                                                      flows_1D_MEA[i].J_H2,
                                                                      flows_1D_MEA[i].J_O2,
                                                                      flows_1D_MEA[i].S_H2,
                                                                      flows_1D_MEA[i].S_O2)
                          for i in 1:nb_gc]
    dif_eq_voltage = [calculate_dyn_voltage_evolution(i_fc[i], C_O2_Pt[i],
                                                      sv_mea_1D[i].ccl.T,
                                                      sv_mea_1D[i].ccl.eta_c,
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
        dif_eq_mea_1D[i] = assemble_mea_derivative_1D(dif_eq_mea_diss_water[i],
                                                      dif_eq_mea_liq_water[i],
                                                      dif_eq_mea_vapor_water[i],
                                                      dif_eq_mea_species[i],
                                                      dif_eq_voltage[i],
                                                      dif_eq_mea_temperature[i])
    end

    #       Inside the gas channels: compute independent GC contributions, then assemble once.
    dif_eq_gc_gas = calculate_dyn_gas_evolution_inside_gas_channel(sv_mea_1D,
                                                                   pp,
                                                                   cfg,
                                                                   flows_1D_GC_manifold,
                                                                   flows_1D_MEA)
    dif_eq_gc_liq = calculate_dyn_liq_evolution_inside_gas_channel(T_des,
                                                                   pp,
                                                                   flows_1D_GC_manifold,
                                                                   flows_1D_MEA)
    dif_eq_gc_temperature = calculate_dyn_temperature_evolution_inside_gas_channel(nb_gc)
    dif_eq_mea_1D = assemble_gc_derivative_1D(dif_eq_mea_1D,
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

    # Repack typed derivatives into solver ordering.
    fuelcell_derivative = FuelCellDerivativeP2D{nb_gdl, nb_mpl, nb_gc}(Tuple(dif_eq_mea_1D))
    _assert_fuelcell_derivative_complete(fuelcell_derivative)
    dif_eq_global = _pack_fuelcell_derivative_p2d(fuelcell_derivative)

    if has_auxiliary
        _assert_manifold_derivative_complete(dif_eq_manifold_1D)
        append!(dif_eq_global, _pack_manifold_derivative_state(dif_eq_manifold_1D))
        _assert_auxiliary_derivative_complete(dif_eq_auxiliary)
        append!(dif_eq_global, _pack_auxiliary_derivative(dif_eq_auxiliary))
    end

    length(dif_eq_global) == expected_len ||
        throw(ArgumentError("Unexpected derivative vector size for typed path."))

    return dif_eq_global
end
