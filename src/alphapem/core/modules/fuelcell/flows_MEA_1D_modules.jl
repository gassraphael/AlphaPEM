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
    D_EOD_acl_mem, D_EOD_mem_ccl, K_lambda_eff_acl_mem, K_lambda_eff_mem_ccl,
    D_cap_agdl_agdl, D_cap_agdl_ampl, D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl,
    D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl, Da_eff_agdl_agdl, Da_eff_agdl_ampl,
    Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl, Dc_eff_cmpl_cgdl,
    Dc_eff_cgdl_cgdl, T_acl_mem_ccl)
"""
function calculate_flows_1D_MEA_int_values!(flows_int_work::MEAFlowsIntWorkspace, sv::CellState1D, i_fc::Float64,
                                            fc::AbstractFuelCell, cfg::SimulationConfig)::Tuple

    # Extraction of the parameters
    pp = fc.physical_parameters
    np = cfg.numerical_parameters
    Hacl, Hccl, Hmem, Hgdl, Hmpl = pp.Hacl, pp.Hccl, pp.Hmem, pp.Hgdl, pp.Hmpl
    Wagc, Wcgc = pp.Wagc, pp.Wcgc
    epsilon_gdl, epsilon_mpl, epsilon_c = pp.epsilon_gdl, pp.epsilon_mpl, pp.epsilon_c
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Extraction of the variables
    C_v_agc, C_v_agdl = sv.agc.C_v, getproperty.(sv.agdl, :C_v)
    C_v_ampl, C_v_acl = getproperty.(sv.ampl, :C_v), sv.acl.C_v
    C_v_ccl, C_v_cmpl,  = sv.ccl.C_v, getproperty.(sv.cmpl, :C_v)
    C_v_cgdl, C_v_cgc = getproperty.(sv.cgdl, :C_v), sv.cgc.C_v

    s_agc, s_agdl = sv.agc.s, getproperty.(sv.agdl, :s)
    s_ampl, s_acl = getproperty.(sv.ampl, :s), sv.acl.s
    s_ccl, s_cmpl = sv.ccl.s, getproperty.(sv.cmpl, :s)
    s_cgdl, s_cgc = getproperty.(sv.cgdl, :s), sv.cgc.s

    T_agc, T_agdl = sv.agc.T, getproperty.(sv.agdl, :T)
    T_ampl, T_acl = getproperty.(sv.ampl, :T), sv.acl.T
    T_mem = sv.mem.T
    T_ccl, T_cmpl = sv.ccl.T, getproperty.(sv.cmpl, :T)
    T_cgdl, T_cgc = getproperty.(sv.cgdl, :T), sv.cgc.T

    C_H2_agc, C_H2_agdl = sv.agc.C_H2, getproperty.(sv.agdl, :C_H2)
    C_H2_ampl, C_H2_acl = getproperty.(sv.ampl, :C_H2), sv.acl.C_H2

    C_O2_ccl, C_O2_cmpl = sv.ccl.C_O2, getproperty.(sv.cmpl, :C_O2)
    C_O2_cgdl, C_O2_cgc = getproperty.(sv.cgdl, :C_O2), sv.cgc.C_O2

    C_N2_agc, C_N2_agdl = sv.agc.C_N2, getproperty.(sv.agdl, :C_N2)
    C_N2_ampl, C_N2_acl = getproperty.(sv.ampl, :C_N2), sv.acl.C_N2
    C_N2_ccl, C_N2_cmpl = sv.ccl.C_N2, getproperty.(sv.cmpl, :C_N2)
    C_N2_cgdl, C_N2_cgc = getproperty.(sv.cgdl, :C_N2), sv.cgc.C_N2

    lambda_acl, lambda_mem, lambda_ccl = sv.acl.lambda, sv.mem.lambda, sv.ccl.lambda

    # Transitory parameter
    H_gdl_node = Hgdl / nb_gdl
    H_mpl_node = Hmpl / nb_mpl

    # Pressures in the stack
    Pagc  = (C_v_agc + C_H2_agc + C_N2_agc) * R * T_agc
    Pagdl = [(C_v_agdl[i] + C_H2_agdl[i] + C_N2_agdl[i]) * R * T_agdl[i] for i in 1:nb_gdl]
    Pampl = [(C_v_ampl[i] + C_H2_ampl[i] + C_N2_ampl[i]) * R * T_ampl[i] for i in 1:nb_mpl]
    Pacl  = (C_v_acl + C_H2_acl + C_N2_acl) * R * T_acl
    Pccl  = (C_v_ccl + C_O2_ccl + C_N2_ccl) * R * T_ccl
    Pcmpl = [(C_v_cmpl[i] + C_O2_cmpl[i] + C_N2_cmpl[i]) * R * T_cmpl[i] for i in 1:nb_mpl]
    Pcgdl = [(C_v_cgdl[i] + C_O2_cgdl[i] + C_N2_cgdl[i]) * R * T_cgdl[i] for i in 1:nb_gdl]
    Pcgc  = (C_v_cgc + C_O2_cgc + C_N2_cgc) * R * T_cgc

    # Capillary pressures in the stack
    Pcap_agdl = Pcap(:gdl, s_agdl[1],      T_agdl[1],      epsilon_gdl, epsilon_c, pp)
    Pcap_cgdl = Pcap(:gdl, s_cgdl[nb_gdl], T_cgdl[nb_gdl], epsilon_gdl, epsilon_c, pp)

    # Densities in the GC
    rho_agc = C_H2_agc * M_H2 + C_v_agc * M_H2O + C_N2_agc * M_N2
    rho_cgc = C_O2_cgc * M_O2 + C_v_cgc * M_H2O + C_N2_cgc * M_N2

    # EOD interface flux coefficients per lambda unit.
    # The local segment current density i_fc, defined per geometric active area, is assigned to the proton current
    # density at both CL/membrane interfaces in this formulation.
    D_EOD_acl_mem = D_EOD(i_fc)
    D_EOD_mem_ccl = D_EOD(i_fc)

    # Weighted harmonic means of the dissolved-water back-diffusion conductance per lambda gradient.

    K_lambda_acl = cl_dry_ionomer_storage_capacity(:acl, Hacl, pp) * D_lambda_eff(:acl, lambda_acl, T_acl, Hacl, pp)
    K_lambda_mem = pp.rho_mem / pp.M_eq * D_lambda(lambda_mem)
    K_lambda_ccl = cl_dry_ionomer_storage_capacity(:ccl, Hccl, pp) * D_lambda_eff(:ccl, lambda_ccl, T_ccl, Hccl, pp)
    K_lambda_eff_acl_mem = hmean(K_lambda_acl, K_lambda_mem,
                                 Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem))
    K_lambda_eff_mem_ccl = hmean(K_lambda_mem, K_lambda_ccl,
                                 Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl))

    # Pre-computed inter-layer CL porosities and weight factors (avoid repeated calls and divisions)
    epsilon_acl = epsilon_cl(:acl, lambda_acl, T_acl, Hacl, pp)  # CL porosity at the anode side.
    epsilon_ccl = epsilon_cl(:ccl, lambda_ccl, T_ccl, Hccl, pp)  # CL porosity at the cathode side.
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
    D_cap_agdl_agdl = flows_int_work.D_cap_agdl_agdl
    @inbounds for i in 1:(nb_gdl - 1)
        D_cap_agdl_agdl[i] = hmean(Dcap(:gdl, s_agdl[i],     T_agdl[i],     epsilon_gdl, epsilon_c, pp),
                                    Dcap(:gdl, s_agdl[i + 1], T_agdl[i + 1], epsilon_gdl, epsilon_c, pp))
    end

    D_cap_agdl_ampl = hmean(Dcap(:gdl, s_agdl[nb_gdl], T_agdl[nb_gdl], epsilon_gdl, epsilon_c, pp),
                             Dcap(:mpl, s_ampl[1],      T_ampl[1],      epsilon_mpl, nothing, pp),
                             w_gdl_mpl, w_mpl_gdl)

    D_cap_ampl_ampl = flows_int_work.D_cap_ampl_ampl
    @inbounds for i in 1:(nb_mpl - 1)
        D_cap_ampl_ampl[i] = hmean(Dcap(:mpl, s_ampl[i],     T_ampl[i],     epsilon_mpl, nothing, pp),
                                    Dcap(:mpl, s_ampl[i + 1], T_ampl[i + 1], epsilon_mpl, nothing, pp))
    end

    D_cap_ampl_acl = hmean(Dcap(:mpl, s_ampl[nb_mpl], T_ampl[nb_mpl], epsilon_mpl, nothing, pp),
                           Dcap(:cl,  s_acl,          T_acl,          epsilon_acl, nothing, pp),
                           w_mpl_acl, w_acl_mpl)

    D_cap_ccl_cmpl = hmean(Dcap(:cl,  s_ccl,    T_ccl,    epsilon_ccl, nothing, pp),
                           Dcap(:mpl, s_cmpl[1], T_cmpl[1], epsilon_mpl, nothing, pp),
                           w_ccl_mpl, w_mpl_ccl)

    D_cap_cmpl_cmpl = flows_int_work.D_cap_cmpl_cmpl
    @inbounds for i in 1:(nb_mpl - 1)
        D_cap_cmpl_cmpl[i] = hmean(Dcap(:mpl, s_cmpl[i],     T_cmpl[i],     epsilon_mpl, nothing, pp),
                                    Dcap(:mpl, s_cmpl[i + 1], T_cmpl[i + 1], epsilon_mpl, nothing, pp))
    end

    D_cap_cmpl_cgdl = hmean(Dcap(:mpl, s_cmpl[nb_mpl], T_cmpl[nb_mpl], epsilon_mpl, nothing, pp),
                             Dcap(:gdl, s_cgdl[1],      T_cgdl[1],      epsilon_gdl, epsilon_c, pp),
                             w_mpl_gdl, w_gdl_mpl)

    D_cap_cgdl_cgdl = flows_int_work.D_cap_cgdl_cgdl
    @inbounds for i in 1:(nb_gdl - 1)
        D_cap_cgdl_cgdl[i] = hmean(Dcap(:gdl, s_cgdl[i],     T_cgdl[i],     epsilon_gdl, epsilon_c, pp),
                                    Dcap(:gdl, s_cgdl[i + 1], T_cgdl[i + 1], epsilon_gdl, epsilon_c, pp))
    end

    #       ... of the effective diffusion coefficient
    Da_eff_agdl_agdl = flows_int_work.Da_eff_agdl_agdl
    @inbounds for i in 1:(nb_gdl - 1)
        Da_eff_agdl_agdl[i] = hmean(Da_eff(:gdl, s_agdl[i],     T_agdl[i],     Pagdl[i],     epsilon_gdl, epsilon_c, pp),
                                     Da_eff(:gdl, s_agdl[i + 1], T_agdl[i + 1], Pagdl[i + 1], epsilon_gdl, epsilon_c, pp))
    end

    Da_eff_agdl_ampl = hmean(Da_eff(:gdl, s_agdl[nb_gdl], T_agdl[nb_gdl], Pagdl[nb_gdl], epsilon_gdl, epsilon_c, pp),
                              Da_eff(:mpl, s_ampl[1],      T_ampl[1],      Pampl[1],      epsilon_mpl, nothing, pp),
                              w_gdl_mpl, w_mpl_gdl)

    Da_eff_ampl_ampl = flows_int_work.Da_eff_ampl_ampl
    @inbounds for i in 1:(nb_mpl - 1)
        Da_eff_ampl_ampl[i] = hmean(Da_eff(:mpl, s_ampl[i],     T_ampl[i],     Pampl[i],     epsilon_mpl, nothing, pp),
                                     Da_eff(:mpl, s_ampl[i + 1], T_ampl[i + 1], Pampl[i + 1], epsilon_mpl, nothing, pp))
    end

    Da_eff_ampl_acl = hmean(Da_eff(:mpl, s_ampl[nb_mpl], T_ampl[nb_mpl], Pampl[nb_mpl], epsilon_mpl, nothing, pp),
                             Da_eff(:cl,  s_acl,          T_acl,          Pacl,          epsilon_acl, nothing, pp),
                             w_mpl_acl, w_acl_mpl)

    Dc_eff_ccl_cmpl = hmean(Dc_eff(:cl,  s_ccl,    T_ccl,    Pccl,    epsilon_ccl, nothing, pp),
                             Dc_eff(:mpl, s_cmpl[1], T_cmpl[1], Pcmpl[1], epsilon_mpl, nothing, pp),
                             w_ccl_mpl, w_mpl_ccl)

    Dc_eff_cmpl_cmpl = flows_int_work.Dc_eff_cmpl_cmpl
    @inbounds for i in 1:(nb_mpl - 1)
        Dc_eff_cmpl_cmpl[i] = hmean(Dc_eff(:mpl, s_cmpl[i],     T_cmpl[i],     Pcmpl[i],     epsilon_mpl, nothing, pp),
                                     Dc_eff(:mpl, s_cmpl[i + 1], T_cmpl[i + 1], Pcmpl[i + 1], epsilon_mpl, nothing, pp))
    end

    Dc_eff_cmpl_cgdl = hmean(Dc_eff(:mpl, s_cmpl[nb_mpl], T_cmpl[nb_mpl], Pcmpl[nb_mpl], epsilon_mpl, nothing, pp),
                              Dc_eff(:gdl, s_cgdl[1],      T_cgdl[1],      Pcgdl[1],      epsilon_gdl, epsilon_c, pp),
                              w_mpl_gdl, w_gdl_mpl)

    Dc_eff_cgdl_cgdl = flows_int_work.Dc_eff_cgdl_cgdl
    @inbounds for i in 1:(nb_gdl - 1)
        Dc_eff_cgdl_cgdl[i] = hmean(Dc_eff(:gdl, s_cgdl[i],     T_cgdl[i],     Pcgdl[i],     epsilon_gdl, epsilon_c, pp),
                                     Dc_eff(:gdl, s_cgdl[i + 1], T_cgdl[i + 1], Pcgdl[i + 1], epsilon_gdl, epsilon_c, pp))
    end

    #       ... of the temperature
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                            [Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])

    return (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc, D_EOD_acl_mem,
            D_EOD_mem_ccl, K_lambda_eff_acl_mem, K_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl,
            D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl,
            Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl,
            Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl)
end


""" This function calculates the capillary coefficient at the GDL, the MPL or the CL, in kg.m.s-1, considering
GDL compression.

Parameters
----------
element : Symbol
    Specifies the element for which the capillary coefficient is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
s :
    Liquid water saturation variable.
T :
    Temperature in K.
epsilon : Float64
    Porosity.
epsilon_c : Union{Float64, Nothing}
    Compression ratio of the GDL.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
Real
    Capillary coefficient at the GDL, MPL or CL in kg.m.s-1.
"""
function Dcap(element::Symbol,
              s,
              T,
              epsilon::Float64,
              epsilon_c::Union{Float64, Nothing},
              pp::PhysicalParams)

    # Extraction of the parameters
    e = pp.e  # Capillary exponent.
    theta_c_gdl, theta_c_mpl, theta_c_cl = pp.theta_c_gdl, pp.theta_c_mpl, pp.theta_c_cl

    K0_value = K0(element, epsilon, epsilon_c, pp)
    s_eff = _clamped_fraction_value(s)
    if element == :gdl
        theta_c_value = theta_c_gdl
    elseif element == :mpl
        theta_c_value = theta_c_mpl
    elseif element == :cl
        theta_c_value = theta_c_cl
    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end

    return sigma(T) * K0_value / nu_l(T) * abs(cos(theta_c_value)) *
           (epsilon / K0_value)^0.5 * (s_eff^e + 1e-7) *
           (1.417 - 4.24 * s_eff + 3.789 * s_eff^2)
end


""" This function calculates the capillary pressure at the GDL, the MPL or the CL, in kg.m.s-1.

Parameters
----------
element : Symbol
    Specifies the element for which the capillary pressure is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
s :
    Liquid water saturation variable.
T :
    Temperature in K.
epsilon : Float64
    Porosity.
epsilon_c : Union{Float64, Nothing}
    Compression ratio of the GDL.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
Real
    Capillary pressure at the selected element.
"""
function Pcap(element::Symbol,
              s,
              T,
              epsilon::Float64,
              epsilon_c::Union{Float64, Nothing},
              pp::PhysicalParams)

    # Extraction of the parameters
    theta_c_gdl, theta_c_mpl, theta_c_cl = pp.theta_c_gdl, pp.theta_c_mpl, pp.theta_c_cl

    K0_value = K0(element, epsilon, epsilon_c, pp)
    s_eff = _clamped_fraction_value(s)
    if element == :gdl
        theta_c_value = theta_c_gdl
    elseif element == :mpl
        theta_c_value = theta_c_mpl
    elseif element == :cl
        theta_c_value = theta_c_cl
    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end

    return sigma(T) * abs(cos(theta_c_value)) * (epsilon / K0_value)^0.5 *
           (1.417 * s_eff - 2.12 * s_eff^2 + 1.263 * s_eff^3)
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
    T_eff = _positive_temperature_value(T)
    P_eff = _positive_pressure_value(P)
    return 1.644e-4 * (T_eff / 333)^2.334 * (101325 / P_eff)
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
    T_eff = _positive_temperature_value(T)
    P_eff = _positive_pressure_value(P)
    return 3.242e-5 * (T_eff / 333)^2.334 * (101325 / P_eff)
end


"""This function calculates the effective diffusion coefficient at the GDL, MPL or CL and at the anode,
in m².s-1, considering GDL compression.

Parameters
----------
element : Symbol
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
epsilon_c : Union{Float64, Nothing}
    Compression ratio of the GDL.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
Da_eff
    Effective diffusion coefficient at the anode in m².s-1.
"""
function Da_eff(element::Symbol,
                s,
                T,
                P,
                epsilon::Float64,
                epsilon_c::Union{Float64, Nothing},
                pp::PhysicalParams)

    # Extraction of the parameters
    epsilon_p, alpha_p = pp.epsilon_p, pp.alpha_p
    r_s_gdl, r_s_mpl, r_s_cl = pp.r_s_gdl, pp.r_s_mpl, pp.r_s_cl
    tau_mpl, tau_void_cl = pp.tau_mpl, pp.tau_void_cl

    s_eff = _clamped_fraction_value(s)
    if element == :gdl # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67
            beta2 = -1.59
        else
            beta2 = -0.90
        end
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p))^alpha_p)
        return epsilon / tau_gdl * exp(beta2 * epsilon_c) * (1 - s_eff)^r_s_gdl * Da(P, T)

    elseif element == :mpl # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s_eff)^r_s_mpl * Da(P, T)

    elseif element == :cl # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_void_cl * (1 - s_eff)^r_s_cl * Da(P, T)

    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end
end


"""This function calculates the effective diffusion coefficient at the GDL, MPL or CL and at the cathode,
in m².s-1, considering GDL compression.

Parameters
----------
element : Symbol
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
epsilon_c : Union{Float64, Nothing}
    Compression ratio of the GDL.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
Dc_eff
    Effective diffusion coefficient at the cathode in m².s-1.
"""
function Dc_eff(element::Symbol,
                s,
                T,
                P,
                epsilon::Float64,
                epsilon_c::Union{Float64, Nothing},
                pp::PhysicalParams)

    # Extraction of the parameters
    epsilon_p, alpha_p = pp.epsilon_p, pp.alpha_p
    r_s_gdl, r_s_mpl, r_s_cl = pp.r_s_gdl, pp.r_s_mpl, pp.r_s_cl
    tau_mpl, tau_void_cl = pp.tau_mpl, pp.tau_void_cl

    s_eff = _clamped_fraction_value(s)
    if element == :gdl # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67
            beta2 = -1.59
        else
            beta2 = -0.90
        end
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p))^alpha_p)
        return epsilon / tau_gdl * exp(beta2 * epsilon_c) * (1 - s_eff)^r_s_gdl * Dc(P, T)

    elseif element == :mpl # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s_eff)^r_s_mpl * Dc(P, T)

    elseif element == :cl # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_void_cl * (1 - s_eff)^r_s_cl * Dc(P, T)

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
Hinatsu's expression has been selected. It is valid for Nafion®117 S-form membranes for 25 to 130 degC.

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
    return 10.0 + 1.84e-2 * (T - 273.15) + 9.90e-4 * (T - 273.15)^2
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
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
lambda_eq
    Equilibrium water content in the membrane.
"""
function lambda_eq(C_v, s, T, pp::PhysicalParams)
    Kshape = pp.Kshape  # Mathematical factor governing lambda_eq smoothing.
    # Sanitise inputs: during nonlinear iterations the solver may probe unphysical
    # states (C_v < 0, s > 1, T very low).  There, C_v / C_v_sat(T) can diverge to
    # -Inf, making exp(-Kshape*(a_w-1)) overflow to +Inf and the product with
    # (1 + tanh(...)) = 0 evaluate to NaN.  Bounding a_w's components keeps the
    # residual finite so the solver can reject the step instead of crashing.
    C_v_eff = _nonnegative_value(C_v)
    s_eff = _clamped_fraction_value(s)
    a_w = C_v_eff / C_v_sat(T) + 2 * s_eff  # Water activity.
    return 0.5 * lambda_v_eq(a_w)                                          * (1 - tanh(100 * (a_w - 1))) +
           0.5 * (lambda_v_eq(1) + (lambda_l_eq(T) - lambda_v_eq(1)) * (1 - exp(-Kshape * (a_w - 1)))) *
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
    lambda_eff = _nonnegative_value(lambdaa)
    return 4.1e-10 * (lambda_eff / 25.0)^0.15 * (1.0 + tanh((lambda_eff - 2.5) / 1.4))
end


"""This function calculates the effective diffusion coefficient of dissolved water in the catalyst layer ionomer phase,
in m².s-1.

Parameters
----------
element : Symbol
    Either `:acl` (anode) or `:ccl` (cathode).
lambdaa :
    Water content in the catalyst layer ionomer, defined as the number of water molecules per fixed sulfonic-acid site.
T :
    Temperature in K.
Hcl : Float64
    Thickness of the CL layer in m.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
D_lambda_eff
    Effective diffusion coefficient of dissolved water in the catalyst layer ionomer phase in m².s-1.

Notes
-----
The fixed-site storage capacity of the CL ionomer is handled separately through K_lambda = C_fix * D_lambda_eff on a bulk CL-volume basis. Therefore this coefficient applies only the ionomer-phase tortuosity tau_ion to the material diffusion coefficient and does not multiply by the wet ionomer volume fraction epsilon_mc.

tau_void_cl is not used here because it is the pore-structure coefficient for gas transport through the CL pore space, whereas tau_ion describes the effective CL ionomer-network factor. Applying this proton-conduction-derived factor to water diffusion is a modeling assumption documented in tau_ion.
"""
function D_lambda_eff(element::Symbol, lambdaa, T, Hcl::Float64, pp::PhysicalParams)
    return D_lambda(lambdaa) / tau_ion(element, lambdaa, T, Hcl, pp)
end


"""Calculate the signed electro-osmotic drag flux coefficient per lambda unit, in mol.m-2.s-1.

Parameters
----------
i_fc :
    Local segment current density per geometric active area in A.m-2.
    In this formulation, it is assigned to the proton current density at both CL/membrane interfaces.

Returns
-------
D_EOD
    Flux coefficient in mol.m-2.s-1.
    Multiplication by the reconstructed interface lambda gives the unprotected EOD flux.

The caller assigns this coefficient to D_EOD_acl_mem and D_EOD_mem_ccl. These are coefficients per lambda unit,
not water fluxes. Multiplication by lambda_acl_mem_eod and lambda_mem_ccl_eod gives the respective protected EOD
contributions; equal coefficients therefore do not imply equal EOD fluxes. The final net fluxes additionally include
their respective back-diffusion contributions and net-flux protections.

Notes
-----
The local transport law is J_EOD = n_d(lambda) * i_p / F, with n_d(lambda) = 2.5 * lambda / 22.
This is a current-driven flux law, not a diffusion coefficient. The harmonic averaging of back-diffusion conductances
does not by itself justify harmonic averaging of this EOD coefficient.

The dissolved-water inventory balance uses fluxes at the CL/membrane interfaces.
With i_p(interface) = i_fc, the unprotected EOD contributions follow Gass et al., Table 1:
ACL/membrane: (2.5 / 22) * i_fc / F * lambda_acl_mem.
Membrane/CCL: (2.5 / 22) * i_fc / F * lambda_mem_ccl.

With this current basis and the adopted drag law, no additional CL ionomer volume-fraction or tortuosity factor
multiplies the EOD coefficient. Such factors affect transport properties used to calculate currents from potential
gradients; they do not independently reduce the prescribed proton current in the water-per-proton flux relation.

This interface-current assignment is consistent with the layer-integrated protonic
Joule heating for the same linear CL current profile: the interface current is
i_fc, whereas the layer average of i_p^2 is i_fc^2 / 3.

lambda_acl_mem and lambda_mem_ccl denote the interface water contents corresponding to the paper's lambda_acl,mem
and lambda_mem,ccl. They are reconstructed by distance-weighted linear interpolation of the adjacent layer states.
The paper describes arithmetic averaging between nodes but does not explicitly prescribe these distance weights
for interface hydration; this reconstruction follows the later repository implementation. It is an approximation,
not an exact solution of the coupled through-plane EOD/back-diffusion problem. Hydration is a state variable,
not a transport conductance, so the harmonic resistance average used for back diffusion is not applied to hydration.

The implementation additionally applies donor-hydration protections, yielding lambda_acl_mem_eod and
lambda_mem_ccl_eod, and then limits the combined EOD/back-diffusion flux according to its donor direction.
These protections are additional numerical closures, not part of the cited published EOD expressions.
Consequently, the protected flux is not identical to the unprotected expression published by Gass et al.

Sources
-------
1. Gass et al. (2024), An advanced 1D physics-based model for PEM hydrogen fuel cells with enhanced overvoltage
   prediction, Table 1 (ACL/membrane and membrane/CCL dissolved-water fluxes).
   Version-specific preprint: https://arxiv.org/pdf/2404.07508v1.
2. Gass et al. (2024), A critical review of proton exchange membrane fuel cells matter transports and voltage
   polarisation for modelling, Section 2.3 and Eq. 3 (drag law and current density per active area).
   Version-specific preprint: https://arxiv.org/pdf/2410.13323v1.
3. Vetter and Schumacher (2019; preprint 2018), Free open reference implementation of a two-phase PEM fuel cell model,
   Computer Physics Communications, Table 1, Table 3 and Eq. 23. DOI: 10.1016/j.cpc.2018.07.023.
4. Kulikovsky and McIntyre (2011), Heat flux from the catalyst layer of a fuel cell,
   Electrochimica Acta 56, 9172-9179 (CL heat-transport reference). DOI: 10.1016/j.electacta.2011.07.113.
"""
function D_EOD(i_fc)
    return 2.5 / 22 * i_fc / F
end


const DEFAULT_LAMBDA_CONSTITUTIVE_EPS = 1.0e-8
const DEFAULT_NEGATIVE_EVENT_FLOOR = 1.0e-5
const DEFAULT_LAMBDA_INVENTORY_EPS = 1.0e-4
const DEFAULT_EOD_DONOR_LAMBDA_SCALE = 0.25


"""This function returns a smooth non-negative continuation of max(value, 0).

This helper is used only for dissolved-water constitutive factors. It does not clamp the ODE state itself.
The default smoothing width follows a previous implementation based on the AlphaPEM V1.3 version in Python, where it has been robust for the previous model formulation. It still needs to be re-verified in the integrated (V2.0) Julia-model.
"""
function _lambda_smooth_positive_part(value, eps_value=DEFAULT_LAMBDA_CONSTITUTIVE_EPS)
    return 0.5 * (value + sqrt(value * value + eps_value * eps_value))
end


"""This function returns a smooth donor-inventory availability factor. The factor is 0 when the donor-side dissolved-water inventory is depleted, 1 when sufficient donor inventory is available, and changes smoothly between both limits to avoid a kink in the model equations.

The limiter is applied only to water-removing dissolved-water fluxes. It preserves the raw lambda states and smoothly
reduces additional removal when the donor inventory approaches the lower model floor.

The default floor and ramp width follow the implementation based on the AlphaPEM V1.3 version in Python (`1e-5` and `1e-4`).
They are not newly validated (V2.0) Julia-model parameters; they have been robust in the previous model and must be checked again in the new integrated model.
"""
function _dissolved_inventory_limiter(lambdaa,
                                      inventory_floor=DEFAULT_NEGATIVE_EVENT_FLOOR,
                                      inventory_eps=DEFAULT_LAMBDA_INVENTORY_EPS)
    if inventory_eps <= 0.0
        throw(ArgumentError("The dissolved-water inventory epsilon must be positive."))
    end

    lambda_available = lambdaa - inventory_floor
    if lambda_available <= 0.0
        return 0.0
    elseif lambda_available >= inventory_eps
        return 1.0
    end

    normalized_availability = lambda_available / inventory_eps
    return normalized_availability^2 * (3.0 - 2.0 * normalized_availability)
end


"""This function limits a dissolved-water flux with the inventory of the donor side.

Positive flux follows the local minus-to-plus convention. Negative flux is a return flux and is therefore limited by the plus-side dissolved-water inventory.
The default limiter parameters follow the implementation based on the AlphaPEM V1.3 version in Python and must be re-verified in the new integrated (V2.0) Julia-model.
"""
function _limit_directed_dissolved_flux(flux,
                                        lambda_minus,
                                        lambda_plus,
                                        inventory_floor=DEFAULT_NEGATIVE_EVENT_FLOOR,
                                        inventory_eps=DEFAULT_LAMBDA_INVENTORY_EPS)
    if flux > 0.0
        return flux * _dissolved_inventory_limiter(lambda_minus, inventory_floor, inventory_eps)
    elseif flux < 0.0
        return flux * _dissolved_inventory_limiter(lambda_plus, inventory_floor, inventory_eps)
    else
        return flux
    end
end


"""This function returns a smooth EOD lambda factor limited by the donor-side hydration.

The interface lambda is first reconstructed geometrically from the adjacent finite-volume states. The caller selects
the donor according to the EOD current direction. This helper applies a smooth minimum with the donor hydration to
reduce the gross EOD term when the donor is less hydrated than the reconstructed interface.
It is a state-dependent hydration factor, not an inventory-per-time-step bound; by itself it does not guarantee
nonnegative inventories after a numerical time step. The final net flux is also subject to the directed inventory ramp.

The default smoothing scale (`0.25` in lambda units) follows the implementation based on the AlphaPEM V1.3 version in Python, where it has been robust for the previous model.
This parameterization is a numerical closure for the inventory limiter and still needs to be verified for the new integrated (V2.0) Julia-model.
"""
function _donor_limited_eod_lambda(lambda_interface,
                                   lambda_donor,
                                   smoothing_scale=DEFAULT_EOD_DONOR_LAMBDA_SCALE)
    if smoothing_scale <= 0.0
        throw(ArgumentError("The EOD donor lambda smoothing scale must be positive."))
    end

    lambda_interface_pos = _lambda_smooth_positive_part(lambda_interface)
    lambda_donor_pos = _lambda_smooth_positive_part(lambda_donor)
    lambda_eff_raw = 0.5 * (
        lambda_interface_pos + lambda_donor_pos -
        sqrt((lambda_interface_pos - lambda_donor_pos)^2 + smoothing_scale^2)
    )
    return _lambda_smooth_positive_part(lambda_eff_raw)
end


"""This function calculates the water volume fraction of the membrane.

Parameters
----------
lambdaa :
    Water content in the membrane.
T :
    Temperature in K.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
fv
    Water volume fraction of the membrane.
"""
function fv(lambdaa, T, pp::PhysicalParams)
    M_eq, rho_mem = pp.M_eq, pp.rho_mem  # Equivalent molar mass and density of the dry membrane.
    lambda_eff = _nonnegative_value(lambdaa)
    T_eff = _positive_temperature_value(T)
    return _clamped_fraction_value( (lambda_eff * M_H2O / rho_H2O_l(T_eff)) /
                                    (M_eq / rho_mem + lambda_eff * M_H2O / rho_H2O_l(T_eff)) )
end


"""This function calculates the absorption/desorption rate of vapor in the ionomer, in s-1.

Parameters
----------
C_v :
    Water concentration variable in mol.m-3.
s :
    Liquid water saturation variable.
lambdaa :
    Water content in the ionomer.
T :
    Temperature in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
gamma_sorp
    Absorption/desorption rate of vapor in the ionomer in s-1.
"""
function gamma_sorp_v(C_v, s, lambdaa, T, Hcl::Float64, pp::PhysicalParams)

    T_eff = _positive_temperature_value(T)
    fv_value = fv(lambdaa, T_eff, pp)
    gamma_abs = (1.14e-5 * fv_value) / Hcl * exp(2416 * (1 / 303 - 1 / T_eff))
    gamma_des = (4.59e-5 * fv_value) / Hcl * exp(2416 * (1 / 303 - 1 / T_eff))

    # Transition function between absorption and desorption
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (lambda_eq(C_v, s, T, pp) - lambdaa))) # Transition function.

    return w * gamma_abs + (1 - w) * gamma_des # Interpolation between absorption and desorption.
end


"""This function calculates the phase transfer rate of water condensation or evaporation, in mol.m-3.s-1.
It is positive for condensation and negative for evaporation.

Parameters
----------
element : Symbol
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
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
Svl :
    Phase transfer rate of water condensation or evaporation in mol.m-3.s-1.
"""
function Svl(element::Symbol,
             s,
             C_v,
             Ctot,
             T,
             epsilon::Float64,
             pp::PhysicalParams)

    # Extraction of the parameters
    gamma_cond, gamma_evap = pp.gamma_cond, pp.gamma_evap

    s_eff = _clamped_fraction_value(s)
    C_v_eff = _nonnegative_value(C_v)
    T_eff = _positive_temperature_value(T)
    # Calculation of the total and partial pressures
    Ptot = _positive_pressure_value(Ctot * R * T_eff) # Total pressure.
    P_v = _bounded_vapor_pressure_value(C_v_eff * R * T_eff, Ptot)
    Psat_eff = min(Psat(T_eff), prevfloat(Ptot))

    # Determination of the diffusion coefficient at the anode or the cathode
    if element == :anode
        D_value = Da(Ptot, T_eff)  # Diffusion coefficient at the anode.
    else  # element == :cathode
        D_value = Dc(Ptot, T_eff)  # Diffusion coefficient at the cathode.
    end

    Svl_cond = gamma_cond / (R * T_eff) * epsilon * (1 - s_eff) * D_value * Ptot * log((Ptot - Psat_eff) / (Ptot - P_v))
    Svl_evap = gamma_evap / (R * T_eff) * epsilon * s_eff * D_value * Ptot * log((Ptot - Psat_eff) / (Ptot - P_v))

    # Transition function between condensation and evaporation
    K_transition = 3e-3 # This is a constant that defines the sharpness of the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (Psat_eff - P_v))) # Transition function.

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
    T_eff = _liquid_water_temperature_value(T)
    return 235.8e-3 * ((647.15 - T_eff) / 647.15)^1.256 * (1 - 0.625 * (647.15 - T_eff) / 647.15)
end


"""This function calculates the intrinsic permeability, in m², considering GDL compression.

Parameters
----------
element : Symbol
    Specifies the element for which the intrinsic permeability is calculated.
    Must be either "gdl" (gas diffusion layer), "mpl" (micro-porous layer) or "cl" (catalyst layer).
epsilon : Float64
    Porosity.
epsilon_c : Union{Float64, Nothing}
    Compression ratio of the GDL.
pp : PhysicalParams
    Physical parameters of the fuel cell.

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
function K0(element::Symbol,
            epsilon::Float64,
            epsilon_c::Union{Float64, Nothing},
            pp::PhysicalParams)::Float64

    # Extraction of the parameters
    epsilon_p, alpha_p = pp.epsilon_p, pp.alpha_p
    Dp_mpl, Dp_cl = pp.Dp_mpl, pp.Dp_cl

    if element == :gdl
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67
            beta1 = -3.60
        else
            beta1 = -2.60
        end
        return epsilon / (8 * log(epsilon)^2) * (epsilon - epsilon_p)^(alpha_p + 2) *
               4.6e-6^2 / ((1 - epsilon_p)^alpha_p * ((alpha_p + 1) * epsilon - epsilon_p)^2) * exp(beta1 * epsilon_c)

    elseif element == :mpl
        return (Dp_mpl^2 / 150) * (epsilon^3 / ((1 - epsilon)^2)) # Using the Blake-Kozeny equation.

    elseif element == :cl
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
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
k_H2
    Permeability coefficient of the membrane for hydrogen in mol.m-1.s-1.Pa-1.
"""
function k_H2(lambdaa, T, kappa_co::Float64, pp::PhysicalParams)

    # Extraction of the parameters
    Eact_H2_cros_v, Eact_H2_cros_l = pp.Eact_H2_cros_v, pp.Eact_H2_cros_l

    T_eff = _positive_temperature_value(T)
    # Calculation of the permeability coefficient of the membrane for hydrogen
    k_H2_d = kappa_co * (0.29 + 2.2 * fv(lambdaa, T, pp)) * 1e-14 * exp(Eact_H2_cros_v / R * (1 / Tref_cross - 1 / T_eff))
    k_H2_l = kappa_co * 1.8 * 1e-14 * exp(Eact_H2_cros_l / R * (1 / Tref_cross - 1 / T_eff))

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
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
k_O2
    Permeability coefficient of the membrane for oxygen in mol.m-1.s-1.Pa-1.
"""
function k_O2(lambdaa, T, kappa_co::Float64, pp::PhysicalParams)

    # Extraction of the parameters
    Eact_O2_cros_v, Eact_O2_cros_l = pp.Eact_O2_cros_v, pp.Eact_O2_cros_l

    T_eff = _positive_temperature_value(T)
    # Calculation of the permeability coefficient of the membrane for oxygen
    k_O2_v = kappa_co * (0.11 + 1.9 * fv(lambdaa, T, pp)) * 1e-14 * exp(Eact_O2_cros_v / R * (1 / Tref_cross - 1 / T_eff))
    k_O2_l = kappa_co * 1.2 * 1e-14 * exp(Eact_O2_cros_l / R * (1 / Tref_cross - 1 / T_eff))

    # Transition function between under-saturated and liquid-saturated states
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (lambda_l_eq(T) - lambdaa)))  # Transition function.

    return w * k_O2_v + (1 - w) * k_O2_l  # Interpolation between under-saturated and liquid-equilibrated O2 crossover.
end
