# -*- coding: utf-8 -*-

"""This file represents all the differential equations used for the fuel cell model."""

# ____________________________________________________Main functions____________________________________________________

"""This function calculates the dynamic evolution of the air compressor.

# Arguments
- `dif_eq::Dict`: Dictionary used for saving the differential equations.
- `Pacp_des`: Desired pressure at the anode compressor outlet (Pa).
- `Pasm_out`: Pressure at the anode supply manifold outlet (Pa).
- `Pccp_des`: Desired pressure at the cathode compressor outlet (Pa).
- `Pcsm_out`: Pressure at the cathode supply manifold outlet (Pa).
- `type_auxiliary::Symbol`: Type of auxiliary components used in the fuel cell model.

# Returns
- `Nothing`: `dif_eq` is modified in place.
"""
function calculate_dyn_air_compressor_evolution()
    # This section is intentionally left unchanged for now.
    # dif_eq["dPasm_out / dt"] = (Pacp_des - Pasm_out) / tau_cp # Estimation at the first order.
    return nothing
end


"""This function calculates the dynamic evolution of the humidifiers.

# Arguments
- `dif_eq::Dict`: Dictionary used for saving the differential equations.
- `Wa_inj`: Injected water mass flow at the anode side (kg.s-1).
- `Wc_inj`: Injected water mass flow at the cathode side (kg.s-1).
- `type_auxiliary::Symbol`: Type of auxiliary components used in the fuel cell model.
- `Wa_inj_des`: Desired injected water mass flow at the anode side (kg.s-1).
- `Wc_inj_des`: Desired injected water mass flow at the cathode side (kg.s-1).

# Returns
- `Nothing`: `dif_eq` is modified in place.
"""
function calculate_dyn_humidifier_evolution(dif_eq::Dict, Wa_inj, Wc_inj, type_auxiliary::Symbol,
                                            Wa_inj_des, Wc_inj_des)::Nothing
    # Anode and cathode humidifiers evolution
    if type_auxiliary == :forced_convective_cathode_with_anodic_recirculation
        dif_eq["dWa_inj / dt"] = 0.0
        dif_eq["dWc_inj / dt"] = (Wc_inj_des - Wc_inj) / tau_hum  # Estimation at the first order.
    elseif type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        dif_eq["dWa_inj / dt"] = (Wa_inj_des - Wa_inj) / tau_hum  # Estimation at the first order.
        dif_eq["dWc_inj / dt"] = (Wc_inj_des - Wc_inj) / tau_hum  # Estimation at the first order.
    end

    return nothing
end


"""This function calculates the dynamic evolution of the throttle area inside the anode and cathode auxiliaries.
This function has to be executed after `calculate_dyn_vapor_evolution` and `calculate_dyn_H2_O2_N2_evolution`.

# Arguments
- `dif_eq::Dict`: Dictionary used for saving the differential equations.
- `sv::Dict`: Dictionary containing the solver variables.
- `Pa_des`: Desired pressure inside the anode gas channel (Pa).
- `Pc_des`: Desired pressure inside the cathode gas channel (Pa).
- `A_T_a`: Maximum throttle area at the anode side (m²).
- `A_T_c`: Maximum throttle area at the cathode side (m²).
- `type_auxiliary::Symbol`: Type of auxiliary components used in the fuel cell model.
- `Pagc:`: Pressure inside the anode gas channel (Pa).
- `Pcgc`: Pressure inside the cathode gas channel (Pa).

# Returns
- `Nothing`: `dif_eq` is modified in place.
"""
function calculate_dyn_throttle_area_controler(dif_eq::Dict, sv::Dict, Pa_des, Pc_des, A_T_a, A_T_c,
                                               type_auxiliary::Symbol, Pagc, Pcgc)::Nothing

    # Extraction of the variables
    T_agc, T_cgc = sv["T_agc"], sv["T_cgc"]
    Abp_a, Abp_c = sv["Abp_a"], sv["Abp_c"]

    # Calculation of the pressure derivative inside the gas channels
    dPagcdt = (dif_eq["dC_v_agc / dt"] + dif_eq["dC_H2_agc / dt"] + dif_eq["dC_N2_agc / dt"]) * R * T_agc
    dPcgcdt = (dif_eq["dC_v_cgc / dt"] + dif_eq["dC_O2_cgc / dt"] + dif_eq["dC_N2_cgc / dt"]) * R * T_cgc

    # Throttle area evolution inside the anode auxiliaries
    if type_auxiliary == :forced_convective_cathode_with_flow_through_anode
        dif_eq["dAbp_a / dt"] = -Kp_T * (Pa_des - Pagc) + Kd_T * dPagcdt  # PD controller
        if Abp_a > A_T_a && dif_eq["dAbp_a / dt"] > 0.0  # The throttle area cannot be higher than the maximum value
            dif_eq["dAbp_a / dt"] = 0.0
        elseif Abp_a < 0.0 && dif_eq["dAbp_a / dt"] < 0.0  # The throttle area cannot be lower than 0
            dif_eq["dAbp_a / dt"] = 0.0
        end
    end

    # Throttle area evolution inside the cathode auxiliaries
    dif_eq["dAbp_c / dt"] = -Kp_T * (Pc_des - Pcgc) + Kd_T * dPcgcdt  # PD controller
    if Abp_c > A_T_c && dif_eq["dAbp_c / dt"] > 0.0  # The throttle area cannot be higher than the maximum value
        dif_eq["dAbp_c / dt"] = 0.0
    elseif Abp_c < 0.0 && dif_eq["dAbp_c / dt"] < 0.0  # The throttle area cannot be lower than 0
        dif_eq["dAbp_c / dt"] = 0.0
    end

    return nothing
end

