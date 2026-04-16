# -*- coding: utf-8 -*-

"""This file represents the calculation of the velocity over time. It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using NonlinearSolve
using LineSearches
using SciMLBase


# ________________________________________________________Velocity______________________________________________________


"""
    calculate_velocity_evolution(sv, i_fc_cell, fc, cfg)

Calculate the gas velocities at the anode and cathode.
This function finds the velocities `v_a` and `v_c` that make the pressures computed from the Hagen-Poiseuille
relation equal to the pressures inferred from the ideal gas law applied to the desired molar flows.
A nonlinear root solver from `NonlinearSolve.jl` is used.

Parameters
----------
sv : AbstractVector{<:CellState1D}
    Typed solver state variables. `sv[i]` provides the state of gas-channel node `i`.
i_fc_cell : Float64
    Fuel cell current density at time t (A.m-2).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cfg : SimulationConfig
    Simulation configuration (provides `type_auxiliary`).

Returns
-------
Tuple(Vector{Float64}, Vector{Float64}, Float64, Float64)
    Tuple of (v_a, v_c, Pa_in, Pc_in):
    - Anode gas velocity profile (m.s-1)
    - Cathode gas velocity profile (m.s-1)
    - Anode inlet pressure (Pa)
    - Cathode inlet pressure (Pa)

Raises
------
ErrorException
    If the nonlinear solver does not converge.
"""
function calculate_velocity_evolution(sv::AbstractVector{<:CellState1D}, i_fc_cell::Float64, fc::AbstractFuelCell, cfg::SimulationConfig)::Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64}

    # Extraction of the parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des, Pa_des, Pc_des = oc.T_des, oc.Pa_des, oc.Pc_des
    Hagc, Hcgc, Wagc, Wcgc = pp.Hagc, pp.Hcgc, pp.Wagc, pp.Wcgc
    Lgc, Ldist, nb_channel_in_gc, nb_cell = pp.Lgc, pp.Ldist, pp.nb_channel_in_gc, pp.nb_cell
    nb_gc, nb_gdl = np.nb_gc, np.nb_gdl

    # Extraction of the variables
    C_v_agc = [sv[i].agc.C_v for i in 1:nb_gc]
    C_v_agdl_1 = [sv[i].agdl[1].C_v for i in 1:nb_gc]
    C_v_cgdl_nb_gdl = [sv[i].cgdl[nb_gdl].C_v for i in 1:nb_gc]
    C_v_cgc = [sv[i].cgc.C_v for i in 1:nb_gc]
    C_H2_agc = [sv[i].agc.C_H2 for i in 1:nb_gc]
    C_H2_agdl_1 = [sv[i].agdl[1].C_H2 for i in 1:nb_gc]
    C_O2_cgdl_nb_gdl = [sv[i].cgdl[nb_gdl].C_O2 for i in 1:nb_gc]
    C_O2_cgc = [sv[i].cgc.C_O2 for i in 1:nb_gc]
    C_N2_agc = [sv[i].agc.C_N2 for i in 1:nb_gc]
    C_N2_cgc = [sv[i].cgc.C_N2 for i in 1:nb_gc]
    T_agc = [sv[i].agc.T for i in 1:nb_gc]
    T_cgc = [sv[i].cgc.T for i in 1:nb_gc]

    # Intermediate calculation
    #       Length of one gas channel node (m).
    L_node_gc = Lgc / nb_gc

    #       Pressures (Pa).
    if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
        cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        Pa_ext = Pext
        Pc_ext = Pext
    else  # cfg.type_auxiliary == :no_auxiliary
        Pa_ext = Pa_des
        Pc_ext = Pc_des
    end
    Pagc = [(C_v_agc[i] + C_H2_agc[i] + C_N2_agc[i]) * R * T_agc[i] for i in 1:nb_gc]
    Pcgc = [(C_v_cgc[i] + C_O2_cgc[i] + C_N2_cgc[i]) * R * T_cgc[i] for i in 1:nb_gc]

    #       H2/O2 ratio in the dry anode/cathode gas mixture (H2/N2 or O2/N2) at the GC.
    y_H2_agc = [C_H2_agc[i] / (C_H2_agc[i] + C_N2_agc[i]) for i in 1:nb_gc]
    y_O2_cgc = [C_O2_cgc[i] / (C_O2_cgc[i] + C_N2_cgc[i]) for i in 1:nb_gc]

    #       Vapor ratio over the gas mixture.
    x_H2O_v_agc = [C_v_agc[i] / (C_v_agc[i] + C_H2_agc[i] + C_N2_agc[i]) for i in 1:nb_gc]
    x_H2O_v_cgc = [C_v_cgc[i] / (C_v_cgc[i] + C_O2_cgc[i] + C_N2_cgc[i]) for i in 1:nb_gc]

    #       Dynamic viscosity of the gas mixture.
    mu_gaz_agc = [mu_mixture_gases(["H2O_v", "H2"], [x_H2O_v_agc[i], 1 - x_H2O_v_agc[i]], T_agc[i])
                  for i in 1:nb_gc]
    mu_gaz_cgc = [mu_mixture_gases(["H2O_v", "O2", "N2"],
                                   [x_H2O_v_cgc[i],
                                    y_O2_cgc[i] * (1 - x_H2O_v_cgc[i]),
                                    (1 - y_O2_cgc[i]) * (1 - x_H2O_v_cgc[i])],
                                   T_cgc[i]) for i in 1:nb_gc]

    # Calculation of the boundary conditions at the GC/GDL interface
    C_tot_agdl = [C_v_agdl_1[i] + C_H2_agdl_1[i] + C_N2_agc[i] for i in 1:nb_gc]
    C_tot_cgdl = [C_v_cgdl_nb_gdl[i] + C_O2_cgdl_nb_gdl[i] + C_N2_cgc[i] for i in 1:nb_gc]
    J_tot_agc_agdl = zeros(Float64, nb_gc)
    J_tot_cgdl_cgc = zeros(Float64, nb_gc)
    @inbounds for i in 1:nb_gc
        Jv_agc_agdl = h_a(Pagc[i], T_des, Wagc, Hagc) * (C_v_agc[i] - C_v_agdl_1[i])  # Also calculated in flows_1D_MEA.jl
        J_H2_agc_agdl = h_a(Pagc[i], T_des, Wagc, Hagc) * (C_H2_agc[i] - C_H2_agdl_1[i])  # Also calculated in flows_1D_MEA.jl
        Jv_cgdl_cgc = h_c(Pcgc[i], T_des, Wcgc, Hcgc) * (C_v_cgdl_nb_gdl[i] - C_v_cgc[i])  # Also calculated in flows_1D_MEA.jl
        J_O2_cgdl_cgc = h_c(Pcgc[i], T_des, Wcgc, Hcgc) * (C_O2_cgdl_nb_gdl[i] - C_O2_cgc[i])  # Also calculated in flows_1D_MEA.jl
        Jl_agc_agdl = 0.0 / M_H2O  # Should be added later, knowing that it requires the knowledge of v_agc...
        Jl_cgdl_cgc = 0.0 / M_H2O  # Should be added later, knowing that it requires the knowledge of v_cgc...
        J_tot_agc_agdl[i] = Jv_agc_agdl + J_H2_agc_agdl + Jl_agc_agdl  # Total molar flow from the AGC to the AGDL (mol.m-2.s-1).
        J_tot_cgdl_cgc[i] = Jv_cgdl_cgc + J_O2_cgdl_cgc + Jl_cgdl_cgc  # Total molar flow from the CGDL to the CGC (mol.m-2.s-1).
    end

    # Pre-allocated Float64 buffers for the residual closure.
    # Re-used across Newton iterations and line-search evaluations (Float64 path only).
    # The ForwardDiff Jacobian path (Dual element type) still allocates, as it needs typed arrays.
    _res_J_a = Vector{Float64}(undef, nb_gc)
    _res_J_c = Vector{Float64}(undef, nb_gc)
    _res_P_a = Vector{Float64}(undef, nb_gc)
    _res_P_c = Vector{Float64}(undef, nb_gc)
    _res_v_a = Vector{Float64}(undef, nb_gc)
    _res_v_c = Vector{Float64}(undef, nb_gc)

    # Shared profile kernel used both by Newton residuals and final output extraction.
    function compute_profiles!(J_a, J_c, P_a, P_c, v_a, v_c, J_a_in_guessed, J_c_in_guessed)
        # Continuity equations are used for calculating the molar flows at all points along the GC at stationary state.
        @inbounds for i in 1:nb_gc
            J_a_previous = i == 1 ? J_a_in_guessed : J_a[i - 1]
            J_c_previous = i == 1 ? J_c_in_guessed : J_c[i - 1]
            J_a[i] = J_a_previous - J_tot_agc_agdl[i] * L_node_gc / Hagc
            J_c[i] = J_c_previous + J_tot_cgdl_cgc[i] * L_node_gc / Hcgc
        end
        J_a_out = J_a[end]  # Inside the distributor, there are no mass transfer with the GDL.
        J_c_out = J_c[end]  # Inside the distributor, there are no mass transfer with the GDL.

        # Velocities at the outlets of the GCs (m/s).
        P_a_out = Pa_ext
        P_c_out = Pc_ext
        v_a_out = J_a_out / P_a_out * R * T_des
        v_c_out = J_c_out / P_c_out * R * T_des

        # Backward calculation of pressures and velocities along the GC using modified Hagen-Poiseuille equation.
        P_a[end] = P_a_out + 8 * π * mu_gaz_agc[end] * Ldist / (Hagc * Wagc) * v_a_out  # Hagen-Poiseuille equation inside the distributor.
        v_a[end] = J_a[end] / P_a[end] * R * T_des  # It is considered that the outlet gas viscosity at static state is close to mu_gaz_agc[end] in order to simplify the calculations.
        P_c[end] = P_c_out + 8 * π * mu_gaz_cgc[end] * Ldist / (Hcgc * Wcgc) * v_c_out
        v_c[end] = J_c[end] / P_c[end] * R * T_des

        @inbounds for i in nb_gc:-1:2
            # At the node side
            P_a[i - 1] = P_a[i] + 8 * π * mu_gaz_agc[i] * L_node_gc / (Hagc * Wagc) *
                         (v_a[i] + J_tot_agc_agdl[i] / C_tot_agdl[i])  # Modified Hagen-Poiseuille equation, considering the convective vapor flow from the GC to the GDL.
            v_a[i - 1] = J_a[i - 1] / P_a[i - 1] * R * T_des  # Velocity calculated using the known molar flow and pressure at node i-1.

            # At the cathode side
            P_c[i - 1] = P_c[i] + 8 * π * mu_gaz_cgc[i] * L_node_gc / (Hcgc * Wcgc) *
                         (v_c[i] - J_tot_cgdl_cgc[i] / C_tot_cgdl[i])  # It is considered that the gas viscosity at static state is close to mu_gaz_cgc[i] in order to simplify the calculations.
            v_c[i - 1] = J_c[i - 1] / P_c[i - 1] * R * T_des
        end

        P_a_in = P_a[1] + 8 * π * mu_gaz_agc[1] * L_node_gc / (Hagc * Wagc) *
                 (v_a[1] + J_tot_agc_agdl[1] / C_tot_agdl[1])
        P_c_in = P_c[1] + 8 * π * mu_gaz_cgc[1] * L_node_gc / (Hcgc * Wcgc) *
                 (v_c[1] - J_tot_cgdl_cgc[1] / C_tot_cgdl[1])

        return P_a_in, P_c_in
    end

    # In-place residual function for the nonlinear solver applied to inlet molar flows.
    function residuals!(res, J_in_guessed, _)
        # Intermediate values
        FT = eltype(J_in_guessed)  # Numeric type of the guessed inlet molar flows (Float64 or Dual for autodiff).
        J_a_in_guessed, J_c_in_guessed = J_in_guessed[1], J_in_guessed[2]
        if FT === Float64
            # Reuse pre-allocated buffers — no heap allocation on the hot path.
            J_a, J_c = _res_J_a, _res_J_c
            P_a, P_c = _res_P_a, _res_P_c
            v_a, v_c = _res_v_a, _res_v_c
        else
            # ForwardDiff Jacobian evaluation: element type is Dual — must allocate.
            J_a, J_c = Vector{FT}(undef, nb_gc), Vector{FT}(undef, nb_gc)
            P_a, P_c = Vector{FT}(undef, nb_gc), Vector{FT}(undef, nb_gc)
            v_a, v_c = Vector{FT}(undef, nb_gc), Vector{FT}(undef, nb_gc)
        end
        P_a_in, P_c_in = compute_profiles!(J_a, J_c, P_a, P_c, v_a, v_c, J_a_in_guessed, J_c_in_guessed)

        # Desired molar flows at anode and cathode
        if cfg.type_auxiliary == :no_auxiliary
            W_des_calculated = desired_flows(sv, i_fc_cell, P_a_in, P_c_in, fc, cfg)
            Wa_in_calculated = W_des_calculated.H2 + W_des_calculated.H2O_inj_a  # Expression also present in flow calculations
            J_a_in_calculated = Wa_in_calculated / (Hagc * Wagc) / nb_cell / nb_channel_in_gc
            Wc_in_calculated = W_des_calculated.dry_air + W_des_calculated.H2O_inj_c  # This expression is also present in calculate_velocity_evolution.
            J_c_in_calculated = Wc_in_calculated / (Hcgc * Wcgc) / nb_cell / nb_channel_in_gc
        end

        # Residuals: difference between the calculated and guessed inlet molar flows.
        res[1] = J_a_in_calculated - J_a_in_guessed
        res[2] = J_c_in_calculated - J_c_in_guessed
        return nothing
    end

    # Calculation of the initial molar flow using NonlinearSolve and the pressure-drop relation
    #       Initial guesses
    v_medium = 10.0  # Initial guess for the velocity (m/s).
    x0 = [v_medium * Pa_des / (R * T_des), v_medium * Pc_des / (R * T_des)]

    #       Solver call
    prob = NonlinearProblem(residuals!, x0)
    solver = NewtonRaphson(linesearch = LineSearchesJL(; method=LineSearches.BackTracking()))
    sol = solve(prob, solver; abstol=1e-10, reltol=1e-10, maxiters=200)

    #       Check for convergence
    if !SciMLBase.successful_retcode(sol)
        throw(ErrorException("Convergence failed in calculate_velocity_evolution: $(sol.retcode)"))
    end

    #       Extract initial flow rates
    J_a_in, J_c_in = sol.u[1], sol.u[2]

    # Reuse pre-allocated Float64 buffers to build final profiles (no extra vector allocation).
    P_a_in, P_c_in = compute_profiles!(_res_J_a, _res_J_c, _res_P_a, _res_P_c, _res_v_a, _res_v_c, J_a_in, J_c_in)

    return _res_v_a, _res_v_c, P_a_in, P_c_in
end


"""
    desired_flows(solver_variables, i_fc_cell, Pa_in, Pc_in, fc, cfg)

Calculate the desired flow for the air compressor and the humidifiers.

Parameters
----------
solver_variables : AbstractVector{<:CellState1D}
    Typed variables calculated by the solver at each gas-channel node.
i_fc_cell : Real
    Fuel cell current density (A.m-2).
Pa_in : Real
    Inlet pressure at the anode side (Pa).
Pc_in : Real
    Inlet pressure at the cathode side (Pa).
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.
cfg : SimulationConfig
    Simulation configuration (provides `type_auxiliary`).

Returns
-------
DesiredInletFlows
    Desired hydrogen flow rate, dry-air flow rate, anode humidifier flow rate, and cathode humidifier flow rate.
"""
@inline function desired_flows(solver_variables::AbstractVector{<:CellState1D},
                       i_fc_cell::Real,
                       Pa_in::Real,
                       Pc_in::Real,
                       fc::AbstractFuelCell,
                       cfg::SimulationConfig)::DesiredInletFlows

    # Extraction of the parameters
    oc = fc.operating_conditions
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    T_des, Sa, Sc, Phi_a_des, Phi_c_des, y_H2_in = oc.T_des, oc.Sa, oc.Sc, oc.Phi_a_des, oc.Phi_c_des, oc.y_H2_in
    Hacl, Hmem, Hccl, Aact = pp.Hacl, pp.Hmem, pp.Hccl, pp.Aact
    kappa_co = pp.kappa_co
    nb_gc, nb_cell = np.nb_gc, pp.nb_cell

    # Physical quantities inside the stack
    #       The crossover current density i_n
    inv_H_total = 1.0 / (Hacl + Hmem + Hccl)
    w_acl = Hacl * inv_H_total
    w_mem = Hmem * inv_H_total
    w_ccl = Hccl * inv_H_total
    max_i_n = -Inf
    @inbounds for i in 1:nb_gc
        s_i = solver_variables[i]
        T_acl_i = s_i.acl.T
        T_mem_i = s_i.mem.T
        T_ccl_i = s_i.ccl.T
        lambda_mem_i = s_i.mem.lambda
        C_H2_acl_i = s_i.acl.C_H2
        C_O2_ccl_i = s_i.ccl.C_O2

        T_acl_mem_ccl = w_acl * T_acl_i + w_mem * T_mem_i + w_ccl * T_ccl_i
        i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl_i * k_H2(lambda_mem_i, T_mem_i, kappa_co)
        i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl_i * k_O2(lambda_mem_i, T_mem_i, kappa_co)
        i_n_i = i_H2 + i_O2
        max_i_n = max(max_i_n, i_n_i)
    end

    if cfg.type_auxiliary == :forced_convective_cathode_with_anodic_recirculation ||
       cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        # The desired air compressor volume flow rate (mol.s-1). Warning: consider the minimum compressor flow!
        W_H2_des = 1 / y_H2_in * Sa * i_fc_cell / (2 * F) * (nb_cell * Aact)
        Wacp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc_cell, W_H2_des) # Adjust the desired compressor flow rate to ensure a minimum flow is maintained, based on the current density.
        W_air_ext_des = Pext / (Pext - Phi_ext * Psat(Text)) * 1 / y_O2_ext * Sc * i_fc_cell / (4 * F) * (nb_cell * Aact)
        Wccp_des_adjusted = adjust_compressor_flow_with_minimum(i_fc_cell, W_air_ext_des) # Adjust the desired compressor flow rate to ensure a minimum flow is maintained, based on the current density.

        # The desired humidifier volume flow rate at the anode side Wa_v_inj_des (mol.s-1). Warning: consider the minimum compressor flow!
        if cfg.type_auxiliary == :forced_convective_cathode_with_flow_through_anode
             Pasm = NaN
             Prd = Pasm
             W_H2_des = 1 / y_H2_in * Sa * i_fc_cell / (2 * F) * (nb_cell * Aact)
             W_H2O_inj_a_des = Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des)) /
                               (1 - Phi_a_des * Psat(T_des) / (Prd + Phi_a_des * Psat(T_des))) * W_H2_des
        else  # type_auxiliary == "forced-convective_cathode_with_anodic_recirculation"
             W_H2O_inj_a_des = 0
        end

         # The desired humidifier volume flow rate at the cathode side Wc_v_inj_des (mol.s-1). Warning: consider the minimum compressor flow!
         Pcsm = NaN
         Wcp = NaN
         Pcp = Pcsm
         Wv_hum_in = Phi_ext * Psat(Text) / Pext * Wcp  # Vapor flow rate from the outside.
         W_H20_c_des = Phi_c_des * Psat(T_des) / Pcp * Wcp  # Desired vapor flow rate.
         W_H2O_inj_c_des = W_H20_c_des - Wv_hum_in  # Desired humidifier flow rate.

        _ = Wacp_des_adjusted
        _ = Wccp_des_adjusted
        throw(ErrorException("desired_flows is not yet implemented in Julia for the forced-convective auxiliary branches because the translated Python source still depends on unavailable variables such as Pasm, Pcsm, and Wcp."))
    else  # cfg.type_auxiliary == :no_auxiliary
        # At the anode side
        W_H2_des = Sa * (i_fc_cell + max_i_n) / (2 * F) * (nb_cell * Aact)
        W_H2O_inj_a_des = (Phi_a_des * Psat(T_des) / (Pa_in - Phi_a_des * Psat(T_des))) * W_H2_des

        # At the cathode side
        W_dry_air_des = 1 / y_O2_ext * Sc * (i_fc_cell + max_i_n) / (4 * F) * (nb_cell * Aact)
        W_H2O_inj_c_des = (Phi_c_des * Psat(T_des) / (Pc_in - Phi_c_des * Psat(T_des))) * W_dry_air_des
    end

    return DesiredInletFlows(W_H2_des, W_dry_air_des, W_H2O_inj_a_des, W_H2O_inj_c_des)
end


"""
    adjust_compressor_flow_with_minimum(i_fc_cell, Wcp_des)

Adjust the desired compressor flow rate to ensure a minimum flow is maintained, based on the current density.

Parameters
----------
i_fc_cell : Real
    Actual fuel cell current density (A.m-2).
Wcp_des : Real
    Desired compressor flow rate (mol.s-1).

Returns
-------
Real
    Adjusted compressor flow rate (mol.s-1) ensuring the minimum flow.
"""
@inline function adjust_compressor_flow_with_minimum(i_fc_cell::Real, Wcp_des::Real)::Real

    # Parameters for minimum current density adjustment
    i_cp_min = 0.3e4  # (A/m²) Minimum current density for compressor flow.
    delta_i_load_step = 0.01e4  # (A/m²) Minimum current density step for reaching the minimum compressor flow.

    if i_fc_cell <= i_cp_min + 3 * delta_i_load_step
        Wcp_des_adjusted = (
            Wcp_des * i_cp_min / i_fc_cell * (1.0 + tanh(4 * (i_fc_cell - (delta_i_load_step / 2)) / (delta_i_load_step / 2))) / 2 +
            Wcp_des * (1 - i_cp_min / i_fc_cell) * (1.0 + tanh(4 * (i_fc_cell - i_cp_min - (delta_i_load_step / 2)) / (delta_i_load_step / 2))) / 2
        )
        return Wcp_des_adjusted
    else  # For higher current densities, the compressor flow is not adjusted, and so it is faster to return the original value.
        return Wcp_des
    end
end

