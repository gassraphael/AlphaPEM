# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the flows calculation.
"""

# _____________________________________________________Flow modules_____________________________________________________

"""Calculate intermediate values for the flows calculation.

Parameters
----------
sv : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
i_fc : Float64
    Current density of the fuel cell (A/m²).
fc : AbstractFuelCell
    The fuel cell instance providing model parameters.
cfg : SimulationConfig
    Simulation configuration (provides numerical parameters).

Returns
-------
Tuple(29 elements)
    Tuple containing all intermediate values used by the flows calculation.
    Elements: (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc,
    D_eff_EOD_acl_mem, D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem, D_lambda_eff_mem_ccl,
    D_cap_agdl_agdl, D_cap_agdl_ampl, D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl,
    D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl, Da_eff_agdl_agdl, Da_eff_agdl_ampl,
    Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl, Dc_eff_cmpl_cgdl,
    Dc_eff_cgdl_cgdl, T_acl_mem_ccl)
"""
function flows_1D_MEA_int_values(sv::CellState1D, i_fc::Float64, fc::AbstractFuelCell, cfg::SimulationConfig)::Tuple

    # Extraction of the parameters
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    Hacl, Hccl, Hmem, Hgdl, Hmpl = pp.Hacl, pp.Hccl, pp.Hmem, pp.Hgdl, pp.Hmpl
    Wagc, Wcgc = pp.Wagc, pp.Wcgc
    epsilon_gdl, epsilon_mpl, epsilon_c, e = pp.epsilon_gdl, pp.epsilon_mpl, pp.epsilon_c, pp.e
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv.agc.C_v, sv.acl.C_v, sv.ccl.C_v, sv.cgc.C_v
    C_v_agdl = [sv.agdl[i].C_v for i in 1:nb_gdl]
    C_v_ampl = [sv.ampl[i].C_v for i in 1:nb_mpl]
    C_v_cmpl = [sv.cmpl[i].C_v for i in 1:nb_mpl]
    C_v_cgdl = [sv.cgdl[i].C_v for i in 1:nb_gdl]
    s_acl, s_ccl = sv.acl.s, sv.ccl.s
    s_agdl = [sv.agdl[i].s for i in 1:nb_gdl]
    s_ampl = [sv.ampl[i].s for i in 1:nb_mpl]
    s_cmpl = [sv.cmpl[i].s for i in 1:nb_mpl]
    s_cgdl = [sv.cgdl[i].s for i in 1:nb_gdl]
    lambda_acl, lambda_mem, lambda_ccl = sv.acl.lambda, sv.mem.lambda, sv.ccl.lambda
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv.agc.C_H2, sv.acl.C_H2, sv.ccl.C_O2, sv.cgc.C_O2
    C_H2_agdl = [sv.agdl[i].C_H2 for i in 1:nb_gdl]
    C_H2_ampl = [sv.ampl[i].C_H2 for i in 1:nb_mpl]
    C_O2_cmpl = [sv.cmpl[i].C_O2 for i in 1:nb_mpl]
    C_O2_cgdl = [sv.cgdl[i].C_O2 for i in 1:nb_gdl]
    C_N2_agc, C_N2_cgc = sv.agc.C_N2, sv.cgc.C_N2
    T_agc, T_acl, T_mem, T_ccl, T_cgc = sv.agc.T, sv.acl.T, sv.mem.T, sv.ccl.T, sv.cgc.T
    T_agdl = [sv.agdl[i].T for i in 1:nb_gdl]
    T_ampl = [sv.ampl[i].T for i in 1:nb_mpl]
    T_cmpl = [sv.cmpl[i].T for i in 1:nb_mpl]
    T_cgdl = [sv.cgdl[i].T for i in 1:nb_gdl]

    # Transitory parameter
    H_gdl_node = Hgdl / nb_gdl
    H_mpl_node = Hmpl / nb_mpl

    # Pressures in the stack
    Pagc  = (C_v_agc + C_H2_agc + C_N2_agc) * R * T_agc
    Pagdl = [(C_v_agdl[i] + C_H2_agdl[i] + C_N2_agc) * R * T_agdl[i] for i in 1:nb_gdl]
    Pampl = [(C_v_ampl[i] + C_H2_ampl[i] + C_N2_agc) * R * T_ampl[i] for i in 1:nb_mpl]
    Pacl  = (C_v_acl + C_H2_acl + C_N2_agc) * R * T_acl
    Pccl  = (C_v_ccl + C_O2_ccl + C_N2_cgc) * R * T_ccl
    Pcmpl = [(C_v_cmpl[i] + C_O2_cmpl[i] + C_N2_cgc) * R * T_cmpl[i] for i in 1:nb_mpl]
    Pcgdl = [(C_v_cgdl[i] + C_O2_cgdl[i] + C_N2_cgc) * R * T_cgdl[i] for i in 1:nb_gdl]
    Pcgc  = (C_v_cgc + C_O2_cgc + C_N2_cgc) * R * T_cgc

    # Capillary pressures in the stack
    Pcap_agdl = Pcap("gdl", s_agdl[1],      T_agdl[1],      epsilon_gdl, epsilon_c)
    Pcap_cgdl = Pcap("gdl", s_cgdl[nb_gdl], T_cgdl[nb_gdl], epsilon_gdl, epsilon_c)

    # Densities in the GC
    rho_agc = C_H2_agc * M_H2 + C_v_agc * M_H2O + C_N2_agc * M_N2
    rho_cgc = C_O2_cgc * M_O2 + C_v_cgc * M_H2O + C_N2_cgc * M_N2

    # Weighted mean values ...
    #       ... of the EOD flow of water in the membrane
    D_eff_EOD_acl_mem = hmean([D_EOD_eff(i_fc, lambda_acl, T_acl, Hacl), D_EOD(lambda_mem)],
                              [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_eff_EOD_mem_ccl = hmean([D_EOD(lambda_mem), D_EOD_eff(i_fc, lambda_ccl, T_ccl, Hccl)],
                              [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])

    #       ... of the diffusion coefficient of water in the membrane
    D_lambda_eff_acl_mem = hmean([D_lambda_eff(lambda_acl, T_acl, Hacl), D_lambda(lambda_mem)],
                                 [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_lambda_eff_mem_ccl = hmean([D_lambda(lambda_mem), D_lambda_eff(lambda_ccl, T_ccl, Hccl)],
                                 [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])

    # Pre-computed inter-layer CL porosities and weight factors (avoid repeated calls and divisions)
    epsilon_acl = epsilon_cl(lambda_acl, T_acl, Hacl)  # CL porosity at the anode side.
    epsilon_ccl = epsilon_cl(lambda_ccl, T_ccl, Hccl)  # CL porosity at the cathode side.
    H_gdl_mpl  = H_gdl_node + H_mpl_node                  # Sum of GDL and MPL node thicknesses.
    H_mpl_acl  = H_mpl_node + Hacl                        # Sum of MPL and ACL thicknesses.
    H_ccl_mpl  = Hccl + H_mpl_node                        # Sum of CCL and MPL thicknesses.
    w_gdl_mpl  = H_gdl_node / H_gdl_mpl                   # GDL-side weight at GDL/MPL interface.
    w_mpl_gdl  = H_mpl_node / H_gdl_mpl                   # MPL-side weight at GDL/MPL interface.
    w_mpl_acl  = H_mpl_node / H_mpl_acl                   # MPL-side weight at MPL/ACL interface.
    w_acl_mpl  = Hacl       / H_mpl_acl                   # ACL-side weight at MPL/ACL interface.
    w_ccl_mpl  = Hccl       / H_ccl_mpl                   # CCL-side weight at CCL/MPL interface.
    w_mpl_ccl  = H_mpl_node / H_ccl_mpl                   # MPL-side weight at CCL/MPL interface.

    #       ... of the capillary coefficient
    D_cap_agdl_agdl = Vector(undef, nb_gdl - 1)
    @inbounds for i in 1:(nb_gdl - 1)
        D_cap_agdl_agdl[i] = hmean([Dcap("gdl", s_agdl[i],     T_agdl[i],     epsilon_gdl, e, epsilon_c),
                         Dcap("gdl", s_agdl[i + 1], T_agdl[i + 1], epsilon_gdl, e, epsilon_c)])
    end

    D_cap_agdl_ampl = hmean([Dcap("gdl", s_agdl[nb_gdl], T_agdl[nb_gdl], epsilon_gdl, e, epsilon_c),
                             Dcap("mpl", s_ampl[1],      T_ampl[1],      epsilon_mpl, e)],
                            [w_gdl_mpl, w_mpl_gdl])

    D_cap_ampl_ampl = Vector(undef, nb_mpl - 1)
    @inbounds for i in 1:(nb_mpl - 1)
        D_cap_ampl_ampl[i] = hmean([Dcap("mpl", s_ampl[i],     T_ampl[i],     epsilon_mpl, e),
                                     Dcap("mpl", s_ampl[i + 1], T_ampl[i + 1], epsilon_mpl, e)])
    end

    D_cap_ampl_acl = hmean([Dcap("mpl", s_ampl[nb_mpl], T_ampl[nb_mpl], epsilon_mpl, e),
                            Dcap("cl",  s_acl,          T_acl,          epsilon_acl, e)],
                           [w_mpl_acl, w_acl_mpl])

    D_cap_ccl_cmpl = hmean([Dcap("cl",  s_ccl,    T_ccl,    epsilon_ccl, e),
                            Dcap("mpl", s_cmpl[1], T_cmpl[1], epsilon_mpl, e)],
                           [w_ccl_mpl, w_mpl_ccl])

    D_cap_cmpl_cmpl = Vector(undef, nb_mpl - 1)
    @inbounds for i in 1:(nb_mpl - 1)
        D_cap_cmpl_cmpl[i] = hmean([Dcap("mpl", s_cmpl[i],     T_cmpl[i],     epsilon_mpl, e),
                                     Dcap("mpl", s_cmpl[i + 1], T_cmpl[i + 1], epsilon_mpl, e)])
    end

    D_cap_cmpl_cgdl = hmean([Dcap("mpl", s_cmpl[nb_mpl], T_cmpl[nb_mpl], epsilon_mpl, e),
                             Dcap("gdl", s_cgdl[1],      T_cgdl[1],      epsilon_gdl, e, epsilon_c)],
                            [w_mpl_gdl, w_gdl_mpl])

    D_cap_cgdl_cgdl = Vector(undef, nb_gdl - 1)
    @inbounds for i in 1:(nb_gdl - 1)
        D_cap_cgdl_cgdl[i] = hmean([Dcap("gdl", s_cgdl[i],     T_cgdl[i],     epsilon_gdl, e, epsilon_c),
                         Dcap("gdl", s_cgdl[i + 1], T_cgdl[i + 1], epsilon_gdl, e, epsilon_c)])
    end

    #       ... of the effective diffusion coefficient
    Da_eff_agdl_agdl = Vector(undef, nb_gdl - 1)
    @inbounds for i in 1:(nb_gdl - 1)
        Da_eff_agdl_agdl[i] = hmean([Da_eff("gdl", s_agdl[i],     T_agdl[i],     Pagdl[i],     epsilon_gdl, epsilon_c),
                          Da_eff("gdl", s_agdl[i + 1], T_agdl[i + 1], Pagdl[i + 1], epsilon_gdl, epsilon_c)])
    end

    Da_eff_agdl_ampl = hmean([Da_eff("gdl", s_agdl[nb_gdl], T_agdl[nb_gdl], Pagdl[nb_gdl], epsilon_gdl, epsilon_c),
                              Da_eff("mpl", s_ampl[1],      T_ampl[1],      Pampl[1],      epsilon_mpl)],
                             [w_gdl_mpl, w_mpl_gdl])

    Da_eff_ampl_ampl = Vector(undef, nb_mpl - 1)
    @inbounds for i in 1:(nb_mpl - 1)
        Da_eff_ampl_ampl[i] = hmean([Da_eff("mpl", s_ampl[i],     T_ampl[i],     Pampl[i],     epsilon_mpl),
                                      Da_eff("mpl", s_ampl[i + 1], T_ampl[i + 1], Pampl[i + 1], epsilon_mpl)])
    end

    Da_eff_ampl_acl = hmean([Da_eff("mpl", s_ampl[nb_mpl], T_ampl[nb_mpl], Pampl[nb_mpl], epsilon_mpl),
                             Da_eff("cl",  s_acl,          T_acl,          Pacl,          epsilon_acl)],
                            [w_mpl_acl, w_acl_mpl])

    Dc_eff_ccl_cmpl = hmean([Dc_eff("cl",  s_ccl,    T_ccl,    Pccl,    epsilon_ccl),
                             Dc_eff("mpl", s_cmpl[1], T_cmpl[1], Pcmpl[1], epsilon_mpl)],
                            [w_ccl_mpl, w_mpl_ccl])

    Dc_eff_cmpl_cmpl = Vector(undef, nb_mpl - 1)
    @inbounds for i in 1:(nb_mpl - 1)
        Dc_eff_cmpl_cmpl[i] = hmean([Dc_eff("mpl", s_cmpl[i],     T_cmpl[i],     Pcmpl[i],     epsilon_mpl),
                                      Dc_eff("mpl", s_cmpl[i + 1], T_cmpl[i + 1], Pcmpl[i + 1], epsilon_mpl)])
    end

    Dc_eff_cmpl_cgdl = hmean([Dc_eff("mpl", s_cmpl[nb_mpl], T_cmpl[nb_mpl], Pcmpl[nb_mpl], epsilon_mpl),
                              Dc_eff("gdl", s_cgdl[1],      T_cgdl[1],      Pcgdl[1],      epsilon_gdl, epsilon_c)],
                             [w_mpl_gdl, w_gdl_mpl])

    Dc_eff_cgdl_cgdl = Vector(undef, nb_gdl - 1)
    @inbounds for i in 1:(nb_gdl - 1)
        Dc_eff_cgdl_cgdl[i] = hmean([Dc_eff("gdl", s_cgdl[i],     T_cgdl[i],     Pcgdl[i],     epsilon_gdl, epsilon_c),
                          Dc_eff("gdl", s_cgdl[i + 1], T_cgdl[i + 1], Pcgdl[i + 1], epsilon_gdl, epsilon_c)])
    end

    #       ... of the temperature
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                            [Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])

    return (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc, D_eff_EOD_acl_mem,
            D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem, D_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl,
            D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl,
            Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl,
            Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl)
end


""" This function calculates the capillary coefficient at the GDL, the MPL or the CL, in kg.m.s-1, considering
GDL compression.

Parameters
----------
element : String
    Specifies the element for which the capillary coefficient is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
s :
    Liquid water saturation variable.
T :
    Temperature in K.
epsilon : Float64
    Porosity.
e : Int64
    Capillary exponent.
epsilon_c : Union{Float64, Nothing}, optional
    Compression ratio of the GDL.

Returns
-------
Real
    Capillary coefficient at the GDL, MPL or CL in kg.m.s-1.
"""
function Dcap(element::String,
              s,
              T,
              epsilon::Float64,
              e::Int64,
              epsilon_c::Union{Float64, Nothing}=nothing)

    K0_value = K0(element, epsilon, epsilon_c)
    if element == "gdl"
        theta_c_value = theta_c_gdl
    elseif element == "mpl"
        theta_c_value = theta_c_mpl
    elseif element == "cl"
        theta_c_value = theta_c_cl
    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end

    return sigma(T) * K0_value / nu_l(T) * abs(cos(theta_c_value)) *
           (epsilon / K0_value)^0.5 * (s^e + 1e-7) * (1.417 - 4.24 * s + 3.789 * s^2)
end


""" This function calculates the capillary pressure at the GDL, the MPL or the CL, in kg.m.s-1.

Parameters
----------
element : String
    Specifies the element for which the capillary pressure is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
s :
    Liquid water saturation variable.
T :
    Temperature in K.
epsilon : Float64
    Porosity.
epsilon_c : Union{Float64, Nothing}, optional
    Compression ratio of the GDL.

Returns
-------
Real
    Capillary pressure at the selected element.
"""
function Pcap(element::String,
              s,
              T,
              epsilon::Float64,
              epsilon_c::Union{Float64, Nothing}=nothing)

    K0_value = K0(element, epsilon, epsilon_c)
    if element == "gdl"
        theta_c_value = theta_c_gdl
    elseif element == "mpl"
        theta_c_value = theta_c_mpl
    elseif element == "cl"
        theta_c_value = theta_c_cl
    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end

    s_num = s + 1e-7 # To avoid numerical issues when s is slightly negative due to numerical noise.

    return sigma(T) * abs(cos(theta_c_value)) * (epsilon / K0_value)^0.5 *
           (1.417 * s_num - 2.12 * s_num^2 + 1.263 * s_num^3)
end


"""This function calculates the diffusion coefficient at the anode, in m².s-1.

Parameters
----------
P :
    Pressure in Pa.
T :
    Temperature in K.

Returns
-------
Da :
    Diffusion coefficient at the anode in m².s-1.
"""
function Da(P, T)
    return 1.644e-4 * (T / 333)^2.334 * (101325 / P)
end


"""This function calculates the diffusion coefficient at the cathode, in m².s-1.

Parameters
----------
P :
    Pressure in Pa.
T :
    Temperature in K.

Returns
-------
Dc
    Diffusion coefficient at the cathode in m².s-1.
"""
function Dc(P, T)
    return 3.242e-5 * (T / 333)^2.334 * (101325 / P)
end


"""This function calculates the effective diffusion coefficient at the GDL, MPL or CL and at the anode,
in m².s-1, considering GDL compression.

Parameters
----------
element : String
    Specifies the element for which the effective diffusion coefficient is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
s :
    Liquid water saturation variable.
T :
    Temperature in K.
P :
    Pressure in Pa.
epsilon : Float64
    Porosity.
epsilon_c : Union{Float64, Nothing}, optional
    Compression ratio of the GDL.

Returns
-------
Da_eff
    Effective diffusion coefficient at the anode in m².s-1.
"""
function Da_eff(element::String,
                s,
                T,
                P,
                epsilon::Float64,
                epsilon_c::Union{Float64, Nothing}=nothing)

    if element == "gdl" # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67
            beta2 = -1.59
        else
            beta2 = -0.90
        end
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p))^alpha_p)
        return epsilon / tau_gdl * exp(beta2 * epsilon_c) * (1 - s)^r_s_gdl * Da(P, T)

    elseif element == "mpl" # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s)^r_s_mpl * Da(P, T)

    elseif element == "cl" # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_cl * (1 - s)^r_s_cl * Da(P, T)

    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end
end


"""This function calculates the effective diffusion coefficient at the GDL, MPL or CL and at the cathode,
in m².s-1, considering GDL compression.

Parameters
----------
element : String
    Specifies the element for which the effective diffusion coefficient is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
s :
    Liquid water saturation variable.
T :
    Temperature in K.
P :
    Pressure in Pa.
epsilon : Float64
    Porosity.
epsilon_c : Union{Float64, Nothing}, optional
    Compression ratio of the GDL.

Returns
-------
Dc_eff
    Effective diffusion coefficient at the cathode in m².s-1.
"""
function Dc_eff(element::String,
                s,
                T,
                P,
                epsilon::Float64,
                epsilon_c::Union{Float64, Nothing}=nothing)

    if element == "gdl" # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67
            beta2 = -1.59
        else
            beta2 = -0.90
        end
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p))^alpha_p)
        return epsilon / tau_gdl * exp(beta2 * epsilon_c) * (1 - s)^r_s_gdl * Dc(P, T)

    elseif element == "mpl" # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s)^r_s_mpl * Dc(P, T)

    elseif element == "cl" # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_cl * (1 - s)^r_s_cl * Dc(P, T)

    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end
end


"""This function calculates the effective convective-conductive mass transfer coefficient at the anode, in m.s-1.

Parameters
----------
P :
    Pressure in Pa.
T :
    Temperature in K.
Wgc : Float64
    Width of the gas channel in m.
Hgc : Float64
    Thickness of the gas channel in m.

Returns
-------
h_a
    Effective convective-conductive mass transfer coefficient at the anode in m.s-1.
"""
function h_a(P, T, Wgc::Float64, Hgc::Float64)
    Sh = 0.9247 * NaNMath.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.
    return Sh * Da(P, T) / Hgc
end


"""This function calculates the effective convective-conductive mass transfer coefficient at the cathode, in m.s-1.

Parameters
----------
P :
    Pressure in Pa.
T :
    Temperature in K.
Wgc : Float64
    Width of the gas channel in m.
Hgc : Float64
    Thickness of the gas channel in m.

Returns
-------
h_c
    Effective convective-conductive mass transfer coefficient at the cathode in m.s-1.
"""
function h_c(P, T, Wgc::Float64, Hgc::Float64)
    Sh = 0.9247 * NaNMath.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.y
    return Sh * Dc(P, T) / Hgc
end


"""This function calculates the equilibrium water content in the membrane from the vapor phase. Hinatsu's expression
has been selected.

Parameters
----------
a_w :
    Water activity.

Returns
-------
lambda_v_eq
    Equilibrium water content in the membrane from the vapor phase.
"""
function lambda_v_eq(a_w)
    return 0.300 + 10.8 * a_w - 16.0 * a_w^2 + 14.1 * a_w^3
end


"""This function calculates the equilibrium water content in the membrane from the liquid phase.
Hinatsu's expression has been selected. It is valid for N-form membranes for 25 to 100 degC.

Parameters
----------
T :
    Temperature in K.

Returns
-------
lambda_l_eq
    Equilibrium water content in the membrane from the liquid phase.
"""
function lambda_l_eq(T)
    return 10.0 * 1.84e-2 * (T - 273.15) + 9.90e-4 * (T - 273.15)^2
end


"""This function calculates the equilibrium water content in the membrane. Hinatsu's expression modified with
Bao's formulation has been selected.

Parameters
----------
C_v :
    Water concentration variable in mol.m-3.
s :
    Liquid water saturation variable.
T :
    Temperature in K.

Returns
-------
lambda_eq
    Equilibrium water content in the membrane.
"""
function lambda_eq(C_v, s, T)
    a_w = C_v / C_v_sat(T) + 2 * s  # Water activity.
    return 0.5 * lambda_v_eq(a_w)                                          * (1 - tanh(100 * (a_w - 1))) +
           0.5 * (lambda_v_eq(1) + ((lambda_l_eq(T) - lambda_v_eq(1)) / 2) * (1 - exp(-Kshape * (a_w - 1)))) *
                                                                             (1 + tanh(100 * (a_w - 1)))
end


"""This function calculates the diffusion coefficient of water in the bulk membrane, in m².s-1.

Parameters
----------
lambdaa :
    Water content in the membrane.

Returns
-------
D_lambda
    Diffusion coefficient of water in the membrane in m².s-1.
"""
function D_lambda(lambdaa)
    return 4.1e-10 * (lambdaa / 25.0)^0.15 * (1.0 + tanh((lambdaa - 2.5) / 1.4))
end


"""This function calculates the effective diffusion coefficient of water in the ionomer of the catalyst layers,
in m².s-1.

Parameters
----------
lambdaa :
    Water content in the catalyst layer.
T :
    Temperature in K.
Hcl : Float64
    Thickness of the CL layer.

Returns
-------
D_lambda_eff
    Effective diffusion coefficient of water in the catalyst layer in m².s-1.
"""
function D_lambda_eff(lambdaa, T, Hcl::Float64)
    return epsilon_mc(lambdaa, T, Hcl) / tau_cl * D_lambda(lambdaa)
end


"""This function calculates the electro-osmotic drag diffusion coefficient of water in the membrane, in mol.m-2.s-1.

Parameters
----------
i_fc :
    Fuel cell current density in A.m-2.

Returns
-------
D_EOD
    Electro-osmotic drag diffusion coefficient of water in the membrane in mol.m-2.s-1.
"""
function D_EOD(i_fc)
    return 2.5 / 22 * i_fc / F
end


"""This function calculates the effective electro-osmotic drag diffusion coefficient of water in the ionomer of the
catalyst layers, in mol.m-2.s-1.

Parameters
----------
i_fc :
    Fuel cell current density in A.m-2.
lambdaa :
    Water content in the catalyst layer.
T :
    Temperature in K.
Hcl : Float64
    Thickness of the CL layer.

Returns
-------

    Effective electro-osmotic drag diffusion coefficient of water in the catalyst layer in mol.m-2.s-1.
"""
function D_EOD_eff(i_fc, lambdaa, T, Hcl::Float64)
    return epsilon_mc(lambdaa, T, Hcl) / tau_cl * D_EOD(i_fc)
end


"""This function calculates the water volume fraction of the membrane.

Parameters
----------
lambdaa :
    Water content in the membrane.
T :
    Temperature in K.

Returns
-------
fv
    Water volume fraction of the membrane.
"""
function fv(lambdaa, T)
    return (lambdaa * M_H2O / rho_H2O_l(T)) / (M_eq / rho_mem + lambdaa * M_H2O / rho_H2O_l(T))
end


"""This function calculates the sorption rate of water in the membrane, in s-1.

Parameters
----------
C_v :
    Water concentration variable in mol.m-3.
s :
    Liquid water saturation variable.
lambdaa :
    Water content in the membrane.
T :
    Temperature in K.
Hcl : Float64
    Thickness of the CL layer.

Returns
-------
gamma_sorp
    Sorption rate of water in the membrane in s-1.
"""
function gamma_sorp(C_v, s, lambdaa, T, Hcl::Float64)

    fv_value = fv(lambdaa, T)
    gamma_abs = (1.14e-5 * fv_value) / Hcl * exp(2416 * (1 / 303 - 1 / T))
    gamma_des = (4.59e-5 * fv_value) / Hcl * exp(2416 * (1 / 303 - 1 / T))

    # Transition function between absorption and desorption
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (lambda_eq(C_v, s, T) - lambdaa))) # Transition function.

    return w * gamma_abs + (1 - w) * gamma_des # Interpolation between absorption and desorption.
end


"""This function calculates the phase transfer rate of water condensation or evaporation, in mol.m-3.s-1.
It is positive for condensation and negative for evaporation.

Parameters
----------
element : String
    Specifies the element for which the phase transfer rate is calculated.
s :
    Liquid water saturation variable.
C_v :
    Water concentration variable in mol.m-3.
Ctot :
    Total gas concentration in mol.m-3.
T :
    Temperature in K.
epsilon : Float64
    Porosity.

Returns
-------
Svl :
    Phase transfer rate of water condensation or evaporation in mol.m-3.s-1.
"""
function Svl(element::String,
             s,
             C_v,
             Ctot,
             T,
             epsilon::Float64)

    # Calculation of the total and partial pressures
    Ptot = Ctot * R * T # Total pressure.
    P_v = C_v * R * T # Partial pressure of vapor.

    # Determination of the diffusion coefficient at the anode or the cathode
    if element == "anode"
        D_value = Da(Ptot, T)  # Diffusion coefficient at the anode.
    else  # element == "cathode"
        D_value = Dc(Ptot, T)  # Diffusion coefficient at the cathode.
    end

    Svl_cond = gamma_cond * M_H2O / (R * T) * epsilon * (1 - s) * D_value * Ptot * NaNMath.log((Ptot - Psat(T)) / (Ptot - P_v))
    Svl_evap = gamma_evap * M_H2O / (R * T) * epsilon * s * D_value * Ptot * NaNMath.log((Ptot - Psat(T)) / (Ptot - P_v))

    # Transition function between condensation and evaporation
    K_transition = 3e-3 # This is a constant that defines the sharpness of the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (Psat(T) - P_v))) # Transition function.

    return w * Svl_evap + (1 - w) * Svl_cond # Interpolation between condensation and evaporation.
end


"""This function calculates the water surface tension, in N.m-1, as a function of the temperature.

Parameters
----------
T :
    Temperature in K.

Returns
-------
sigma :
    Water surface tension in N.m-1.
"""
function sigma(T)
    return 235.8e-3 * ((647.15 - T) / 647.15)^1.256 * (1 - 0.625 * (647.15 - T) / 647.15)
end


"""This function calculates the intrinsic permeability, in m², considering GDL compression.

Parameters
----------
element : String
    Specifies the element for which the intrinsic permeability is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
epsilon : Float64
    Porosity.
epsilon_c : Union{Float64, Nothing}, optional
    Compression ratio of the GDL.

Returns
-------
K0 : Float64
    Intrinsic permeability in m².

Sources
-------
1. Qin Chen 2020 - Two-dimensional multi-physics modeling of porous transport layer in polymer electrolyte membrane
   electrolyzer for water splitting - for the Blake-Kozeny equation.
2. M.L. Stewart 2005 - A study of pore geometry effects on anisotropy in hydraulic permeability using the
   lattice-Boltzmann method - for the Blake-Kozeny equation.
"""
function K0(element::String,
            epsilon::Float64,
            epsilon_c::Union{Float64, Nothing}=nothing)::Float64

    if element == "gdl"
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67
            beta1 = -3.60
        else
            beta1 = -2.60
        end
        return epsilon / (8 * NaNMath.log(epsilon)^2) * (epsilon - epsilon_p)^(alpha_p + 2) *
               4.6e-6^2 / ((1 - epsilon_p)^alpha_p * ((alpha_p + 1) * epsilon - epsilon_p)^2) * exp(beta1 * epsilon_c)

    elseif element == "mpl"
        return (Dp_mpl^2 / 150) * (epsilon^3 / ((1 - epsilon)^2)) # Using the Blake-Kozeny equation.

    elseif element == "cl"
        return (Dp_cl^2 / 150) * (epsilon^3 / ((1 - epsilon)^2)) # Using the Blake-Kozeny equation.

    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end
end


"""This function calculates the permeability coefficient of the membrane for hydrogen, in mol.m-1.s-1.Pa-1.

Parameters
----------
lambdaa :
    Water content in the membrane.
T :
    Temperature in K.
kappa_co : Float64
    Crossover correction coefficient in mol.m-1.s-1.Pa-1.

Returns
-------
k_H2
    Permeability coefficient of the membrane for hydrogen in mol.m-1.s-1.Pa-1.
"""
function k_H2(lambdaa, T, kappa_co::Float64)

    # Calculation of the permeability coefficient of the membrane for hydrogen
    k_H2_d = kappa_co * (0.29 + 2.2 * fv(lambdaa, T)) * 1e-14 * exp(Eact_H2_cros_v / R * (1 / Tref_cross - 1 / T))
    k_H2_l = kappa_co * 1.8 * 1e-14 * exp(Eact_H2_cros_l / R * (1 / Tref_cross - 1 / T))

    # Transition function between under-saturated and liquid-saturated states
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (lambda_l_eq(T) - lambdaa)))  # Transition function.

    return w * k_H2_d + (1 - w) * k_H2_l  # Interpolation between under-saturated and liquid-equilibrated H2 crossover.
end


"""This function calculates the permeability coefficient of the membrane for oxygen, in mol.m-1.s-1.Pa-1.

Parameters
----------
lambdaa :
    Water content in the membrane.
T :
    Temperature in K.
kappa_co : Float64
    Crossover correction coefficient in mol.m-1.s-1.Pa-1.

Returns
-------
k_O2
    Permeability coefficient of the membrane for oxygen in mol.m-1.s-1.Pa-1.
"""
function k_O2(lambdaa, T, kappa_co::Float64)

    # Calculation of the permeability coefficient of the membrane for oxygen
    k_O2_v = kappa_co * (0.11 + 1.9 * fv(lambdaa, T)) * 1e-14 * exp(Eact_O2_cros_v / R * (1 / Tref_cross - 1 / T))
    k_O2_l = kappa_co * 1.2 * 1e-14 * exp(Eact_O2_cros_l / R * (1 / Tref_cross - 1 / T))

    # Transition function between under-saturated and liquid-saturated states
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (lambda_l_eq(T) - lambdaa)))  # Transition function.

    return w * k_O2_v + (1 - w) * k_O2_l  # Interpolation between under-saturated and liquid-equilibrated O2 crossover.
end

