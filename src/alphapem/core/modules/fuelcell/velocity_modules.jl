# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the velocity calculation.
"""

# _____________________________________________________Velocity modules____________________________________________________

"""Calculate intermediate values for velocity computation: GC thermodynamics, viscosities, and GC/GDL interface fluxes (in-place, zero allocation).

Parameters
----------
work : GCManifoldWorkspace
    Pre-allocated workspace. Results are written into `work` fields.
sv : AbstractVector{<:CellState1D}
    Per-GC-node solver state.
T_des : Float64
    Design temperature (K).
fc : AbstractFuelCell
    Fuel-cell instance providing geometric and operating parameters.
cfg : SimulationConfig
    Simulation configuration.

Returns
-------
Nothing
    Results are stored in-place in `work` fields: P_agc, x_H2O_v_agc, x_H2_agc, x_N2_agc, mu_gaz_agc,
    C_tot_agdl, J_tot_agc_agdl, P_cgc, x_H2O_v_cgc, x_O2_cgc, x_N2_cgc, mu_gaz_cgc,
    C_tot_cgdl, J_tot_cgdl_cgc.
"""
function calculate_velocity_int_values!(work::GCManifoldWorkspace,
                                        sv::AbstractVector{<:CellState1D},
                                        T_des::Float64,
                                        fc::AbstractFuelCell,
                                        cfg::SimulationConfig)

    # Extract parameters
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    Hagc, Hcgc, Wagc, Wcgc = pp.Hagc, pp.Hcgc, pp.Wagc, pp.Wcgc
    nb_gc, nb_gdl = np.nb_gc, np.nb_gdl

    @inbounds for (k, i) in enumerate(anode_gc_order(nb_gc, cfg.type_flow)) # In the anode, the flow is reversed for counter-flow.
        # Extract the GC state variables in the good order (counter-flow or co-flow)
        sv_i = sv[i]
        C_v_agc, C_H2_agc, C_N2_agc, T_agc = sv_i.agc.C_v, sv_i.agc.C_H2, sv_i.agc.C_N2, sv_i.agc.T

        # Calculate the AGC pressure and store it in natural order (not counter-flow order)
        C_tot_agc     = C_v_agc + C_H2_agc + C_N2_agc
        P_agc         = C_tot_agc * R * T_agc
        work.P_agc[k] = P_agc

        # Calculate the mole fractions and viscosities, then store them in natural order (not counter-flow order)
        x_H2O_v_agc, x_H2_agc, x_N2_agc = C_v_agc  / C_tot_agc, C_H2_agc / C_tot_agc, C_N2_agc / C_tot_agc
        work.x_H2O_v_agc[k], work.x_H2_agc[k], work.x_N2_agc[k] = x_H2O_v_agc, x_H2_agc, x_N2_agc
        work.mu_gaz_agc[k]  = mu_mixture_gases(:H2O_v, x_H2O_v_agc, :H2, x_H2_agc, :N2, x_N2_agc, T_agc)

        # Calculate the total molar concentration at the AGC/AGDL interface and the net molar flux from AGC to AGDL,
        # then store them in natural order (not counter-flow order)
        C_v_agdl_1, C_H2_agdl_1, C_N2_agdl_1 = sv_i.agdl[1].C_v, sv_i.agdl[1].C_H2, sv_i.agdl[1].C_N2
        h_agc                    = h_a(P_agc, T_des, Wagc, Hagc)
        work.C_tot_agdl[k]       = C_v_agdl_1 + C_H2_agdl_1 + C_N2_agdl_1
        work.J_tot_agc_agdl[k]   = h_agc * (C_v_agc - C_v_agdl_1) +
                                   h_agc * (C_H2_agc - C_H2_agdl_1) +
                                   h_agc * (C_N2_agc - C_N2_agdl_1)
    end

    @inbounds for i in 1:nb_gc # In the cathode
        # Extract the GC state variables
        sv_i   = sv[i]
        C_v_cgc, C_O2_cgc, C_N2_cgc, T_cgc = sv_i.cgc.C_v, sv_i.cgc.C_O2, sv_i.cgc.C_N2, sv_i.cgc.T

        # Calculate the CGC pressure
        C_tot_cgc     = C_v_cgc + C_O2_cgc + C_N2_cgc
        P_cgc         = C_tot_cgc * R * T_cgc
        work.P_cgc[i] = P_cgc

        # Calculate the mole fractions and viscosities
        x_H2O_v_cgc, x_O2_cgc, x_N2_cgc = C_v_cgc / C_tot_cgc, C_O2_cgc / C_tot_cgc, C_N2_cgc / C_tot_cgc
        work.x_H2O_v_cgc[i], work.x_O2_cgc[i], work.x_N2_cgc[i] = x_H2O_v_cgc, x_O2_cgc, x_N2_cgc
        work.mu_gaz_cgc[i]  = mu_mixture_gases(:H2O_v, x_H2O_v_cgc, :O2, x_O2_cgc, :N2, x_N2_cgc, T_cgc)

        # Calculate the total molar concentration at the CGDL/CGC interface and the net molar flux from CGDL to CGC
        C_v_cgdl_last, C_O2_cgdl_last, C_N2_cgdl_last = sv_i.cgdl[nb_gdl].C_v, sv_i.cgdl[nb_gdl].C_O2, sv_i.cgdl[nb_gdl].C_N2
        h_cgc                          = h_c(P_cgc, T_des, Wcgc, Hcgc)
        work.C_tot_cgdl[i]             = C_v_cgdl_last + C_O2_cgdl_last + C_N2_cgdl_last
        work.J_tot_cgdl_cgc[i]         = h_cgc * (C_v_cgdl_last - C_v_cgc) +
                                         h_cgc * (C_O2_cgdl_last - C_O2_cgc) +
                                         h_cgc * (C_N2_cgdl_last - C_N2_cgc)
    end

    return nothing
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
@inline function adjust_compressor_flow_with_minimum(i_fc_cell::Float64, Wcp_des::Float64)::Float64

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
