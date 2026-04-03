# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the heat transfer calculation.
"""

# _________________________________________________Heat transfer modules________________________________________________

"""Calculate intermediate values for the heat transfer calculation.

Parameters
----------
sv : Dict
    Variables calculated by the solver. They correspond to the fuel cell internal states.
    `sv` is a contraction of solver_variables for enhanced readability.
fc : AbstractFuelCell
    The fuel cell instance providing model parameters.

Returns
-------
Tuple
    Tuple containing the intermediate values used by the heat calculation.
"""
function heat_transfer_int_values(sv::Dict,
                                  fc::AbstractFuelCell)::Tuple

    # Extraction of the parameters
    pp = fc.physical_parameters
    np = fc.numerical_parameters
    Hgdl, Hmpl, Hacl, Hccl = pp.Hgdl, pp.Hmpl, pp.Hacl, pp.Hccl
    Hmem, epsilon_gdl, epsilon_mpl, epsilon_c = pp.Hmem, pp.epsilon_gdl, pp.epsilon_mpl, pp.epsilon_c
    nb_gdl, nb_mpl = np.nb_gdl, np.nb_mpl

    # Extraction of the variables
    C_v_acl, C_v_ccl = sv["C_v_acl"], sv["C_v_ccl"]
    C_v_agdl, C_v_ampl = [sv["C_v_agdl_$i"] for i in 1:nb_gdl], [sv["C_v_ampl_$i"] for i in 1:nb_mpl]
    C_v_cmpl, C_v_cgdl = [sv["C_v_cmpl_$i"] for i in 1:nb_mpl], [sv["C_v_cgdl_$i"] for i in 1:nb_gdl]
    s_acl, s_ccl = sv["s_acl"], sv["s_ccl"]
    s_agdl, s_ampl = [sv["s_agdl_$i"] for i in 1:nb_gdl], [sv["s_ampl_$i"] for i in 1:nb_mpl]
    s_cmpl, s_cgdl = [sv["s_cmpl_$i"] for i in 1:nb_mpl], [sv["s_cgdl_$i"] for i in 1:nb_gdl]
    lambda_acl, lambda_mem, lambda_ccl = sv["lambda_acl"], sv["lambda_mem"], sv["lambda_ccl"]
    C_H2_acl, C_O2_ccl = sv["C_H2_acl"], sv["C_O2_ccl"]
    C_H2_agdl, C_H2_ampl = [sv["C_H2_agdl_$i"] for i in 1:nb_gdl], [sv["C_H2_ampl_$i"] for i in 1:nb_mpl]
    C_O2_cmpl, C_O2_cgdl = [sv["C_O2_cmpl_$i"] for i in 1:nb_mpl], [sv["C_O2_cgdl_$i"] for i in 1:nb_gdl]
    C_N2_agc, C_N2_cgc = sv["C_N2_agc"], sv["C_N2_cgc"]
    T_acl, T_mem, T_ccl = sv["T_acl"], sv["T_mem"], sv["T_ccl"]
    T_agdl, T_ampl = [sv["T_agdl_$i"] for i in 1:nb_gdl], [sv["T_ampl_$i"] for i in 1:nb_mpl]
    T_cmpl, T_cgdl = [sv["T_cmpl_$i"] for i in 1:nb_mpl], [sv["T_cgdl_$i"] for i in 1:nb_gdl]

    # Calculation of intermediate values
    Hgdl_node = Hgdl / nb_gdl
    Hmpl_node = Hmpl / nb_mpl

    # Weighted harmonic means of the effective thermal conductivity
    k_th_eff_agc_agdl = k_th_eff("agdl", T_agdl[1], C_v_agdl[1], s_agdl[1], nothing,
                                 C_H2_agdl[1], nothing, C_N2_agc, epsilon_gdl, nothing, epsilon_c)

    k_th_eff_agdl_agdl = Vector{Number}(undef, max(nb_gdl - 1, 0))
    @inbounds for i in 1:(nb_gdl - 1)
        k_th_eff_agdl_agdl[i] = hmean([
            k_th_eff("agdl", T_agdl[i], C_v_agdl[i], s_agdl[i], nothing,
                     C_H2_agdl[i], nothing, C_N2_agc, epsilon_gdl, nothing, epsilon_c),
            k_th_eff("agdl", T_agdl[i + 1], C_v_agdl[i + 1], s_agdl[i + 1], nothing,
                     C_H2_agdl[i + 1], nothing, C_N2_agc, epsilon_gdl, nothing, epsilon_c)
        ])
    end

    k_th_eff_agdl_ampl = hmean([
        k_th_eff("agdl", T_agdl[nb_gdl], C_v_agdl[nb_gdl], s_agdl[nb_gdl], nothing,
             C_H2_agdl[nb_gdl], nothing, C_N2_agc, epsilon_gdl, nothing, epsilon_c),
        k_th_eff("ampl", T_ampl[1], C_v_ampl[1], s_ampl[1], nothing, C_H2_ampl[1],
             nothing, C_N2_agc, epsilon_mpl)
    ], [Hgdl_node / 2, Hmpl_node / 2])

    k_th_eff_ampl_ampl = Vector{Number}(undef, max(nb_mpl - 1, 0))
    @inbounds for i in 1:(nb_mpl - 1)
        k_th_eff_ampl_ampl[i] = hmean([
            k_th_eff("ampl", T_ampl[i], C_v_ampl[i], s_ampl[i], nothing,
                     C_H2_ampl[i], nothing, C_N2_agc, epsilon_mpl),
            k_th_eff("ampl", T_ampl[i + 1], C_v_ampl[i + 1], s_ampl[i + 1], nothing,
                     C_H2_ampl[i + 1], nothing, C_N2_agc, epsilon_mpl)
        ])
    end

    k_th_eff_ampl_acl = hmean([
        k_th_eff("ampl", T_ampl[nb_mpl], C_v_ampl[nb_mpl], s_ampl[nb_mpl], nothing,
             C_H2_ampl[nb_mpl], nothing, C_N2_agc, epsilon_mpl),
        k_th_eff("acl", T_acl, C_v_acl, s_acl, lambda_acl, C_H2_acl, nothing, C_N2_agc, nothing, Hacl)
    ], [Hmpl_node / 2, Hacl / 2])

    k_th_eff_acl_mem = hmean([
        k_th_eff("acl", T_acl, C_v_acl, s_acl, lambda_acl, C_H2_acl, nothing, C_N2_agc, nothing, Hacl),
        k_th_eff("mem", T_mem, nothing, nothing, lambda_mem)
    ], [Hacl / 2, Hmem / 2])

    k_th_eff_mem_ccl = hmean([
        k_th_eff("mem", T_mem, nothing, nothing, lambda_mem),
        k_th_eff("ccl", T_ccl, C_v_ccl, s_ccl, lambda_ccl, nothing, C_O2_ccl, C_N2_cgc, nothing, Hccl)
    ], [Hmem / 2, Hccl / 2])

    k_th_eff_ccl_cmpl = hmean([
        k_th_eff("ccl", T_ccl, C_v_ccl, s_ccl, lambda_ccl, nothing, C_O2_ccl, C_N2_cgc, nothing, Hccl),
        k_th_eff("cmpl", T_cmpl[1], C_v_cmpl[1], s_cmpl[1], nothing, nothing, C_O2_cmpl[1],
             C_N2_cgc, epsilon_mpl)
    ], [Hccl / 2, Hmpl_node / 2])

    k_th_eff_cmpl_cmpl = Vector{Number}(undef, max(nb_mpl - 1, 0))
    @inbounds for i in 1:(nb_mpl - 1)
        k_th_eff_cmpl_cmpl[i] = hmean([
            k_th_eff("cmpl", T_cmpl[i], C_v_cmpl[i], s_cmpl[i], nothing,
                     nothing, C_O2_cmpl[i], C_N2_cgc, epsilon_mpl),
            k_th_eff("cmpl", T_cmpl[i + 1], C_v_cmpl[i + 1], s_cmpl[i + 1], nothing,
                     nothing, C_O2_cmpl[i + 1], C_N2_cgc, epsilon_mpl)
        ])
    end

    k_th_eff_cmpl_cgdl = hmean([
        k_th_eff("cmpl", T_cmpl[nb_mpl], C_v_cmpl[nb_mpl], s_cmpl[nb_mpl], nothing,
             nothing, C_O2_cmpl[nb_mpl], C_N2_cgc, epsilon_mpl),
        k_th_eff("cgdl", T_cgdl[1], C_v_cgdl[1], s_cgdl[1], nothing, nothing, C_O2_cgdl[1],
             C_N2_cgc, epsilon_gdl, nothing, epsilon_c)
    ], [Hmpl_node / 2, Hgdl_node / 2])

    k_th_eff_cgdl_cgdl = Vector{Number}(undef, max(nb_gdl - 1, 0))
    @inbounds for i in 1:(nb_gdl - 1)
        k_th_eff_cgdl_cgdl[i] = hmean([
            k_th_eff("cgdl", T_cgdl[i], C_v_cgdl[i], s_cgdl[i], nothing,
                     nothing, C_O2_cgdl[i], C_N2_cgc, epsilon_gdl, nothing, epsilon_c),
            k_th_eff("cgdl", T_cgdl[i + 1], C_v_cgdl[i + 1], s_cgdl[i + 1], nothing,
                     nothing, C_O2_cgdl[i + 1], C_N2_cgc, epsilon_gdl, nothing, epsilon_c)
        ])
    end

    k_th_eff_cgdl_cgc = k_th_eff("cgdl", T_cgdl[nb_gdl], C_v_cgdl[nb_gdl],
                                 s_cgdl[nb_gdl], nothing, nothing, C_O2_cgdl[nb_gdl], C_N2_cgc,
                                 epsilon_gdl, nothing, epsilon_c)

    return (Hgdl_node, Hmpl_node, k_th_eff_agc_agdl, k_th_eff_agdl_agdl, k_th_eff_agdl_ampl, k_th_eff_ampl_ampl,
            k_th_eff_ampl_acl, k_th_eff_acl_mem, k_th_eff_mem_ccl, k_th_eff_ccl_cmpl, k_th_eff_cmpl_cmpl,
            k_th_eff_cmpl_cgdl, k_th_eff_cgdl_cgdl, k_th_eff_cgdl_cgc)
end


"""This function calculates the effective proton conductivity, in Ω-1.m-1, in either the membrane or the CCL.

Parameters
----------
element : String
    Specifies the element for which the proton conductivity is calculated.
    Must be either `"mem"` (membrane) or `"ccl"` (cathode catalyst layer).
lambdaa :
    Water content in the membrane.
T :
    Temperature in K.
Hcl : Union{Float64, Nothing}
    Thickness of the CL layer.

Returns
-------
sigma_p_eff :
    Proton conductivity in Ω-1.m-1.
"""
function sigma_p_eff(element::String,
                     lambdaa,
                     T,
                     Hcl::Union{Float64, Nothing}=nothing)

    lambda_transition = 1.0

    if element == "mem"  # The proton conductivity at the membrane.
        sigma_p_eff_low = 0.1879 * exp(1268 * (1 / 303.15 - 1 / T))
        sigma_p_eff_high = (0.5139 * lambdaa - 0.326) * exp(1268 * (1 / 303.15 - 1 / T))
    elseif element == "ccl"  # The effective proton conductivity at the cathode catalyst layer.
        if Hcl === nothing
            throw(ArgumentError("Hcl must be provided when element == 'ccl'."))
        end
        sigma_p_eff_low = epsilon_mc(lambdaa, T, Hcl) * 0.1879 * exp(1268 * (1 / 303.15 - 1 / T))
        sigma_p_eff_high = epsilon_mc(lambdaa, T, Hcl) * (0.5139 * lambdaa - 0.326) * exp(1268 * (1 / 303.15 - 1 / T))
    else
        throw(ArgumentError("The element should be either 'mem' or 'ccl'."))
    end

    # Transition function between low and high lambda.
    K_transition = 10.0  # The higher it is, the sharper the transition between two states.
    w = 0.5 * (1 + tanh(K_transition * (lambda_transition - lambdaa)))

    return w * sigma_p_eff_low + (1 - w) * sigma_p_eff_high
end


"""This function calculates the effective electrical conductivity, in Ω-1.m-1, in either the GDL, the MPL or the CL,
considering GDL compression.

Parameters
----------
element : String
    Specifies the element for which the electrical conductivity is calculated.
    Must be either `"gdl"`, `"mpl"`, or `"cl"`.
epsilon : Union{Float64, Nothing}
    Porosity.
epsilon_c : Union{Float64, Nothing}
    Compression ratio of the GDL.
lambda_cl :
    Water content in the CL.
T_cl :
    Temperature inside the CL in K.
Hcl : Union{Real, Nothing}
    Thickness of the CL layer.

Returns
-------
sigma_e_eff
    Effective electrical conductivity in Ω-1.m-1.
"""
function sigma_e_eff(element::String,
                     epsilon::Union{Float64, Nothing}=nothing,
                     epsilon_c::Union{Float64, Nothing}=nothing,
                     lambda_cl=nothing,
                     T_cl=nothing,
                     Hcl::Union{Float64, Nothing}=nothing)

    if element == "gdl"  # The effective electrical conductivity at the GDL.
        if epsilon === nothing || epsilon_c === nothing
            throw(ArgumentError("epsilon and epsilon_c must be provided when element == 'gdl'."))
        end
        beta3 = epsilon < 0.67 ? 4.04 : 4.40
        return (1 - epsilon) * sigma_e_gdl * exp(beta3 * epsilon_c)

    elseif element == "mpl"  # The effective electrical conductivity at the MPL.
        if epsilon === nothing
            throw(ArgumentError("epsilon must be provided when element == 'mpl'."))
        end
        return (1 - epsilon) * sigma_e_mpl

    elseif element == "cl"  # The effective electrical conductivity at the CL.
        if lambda_cl === nothing || T_cl === nothing || Hcl === nothing
            throw(ArgumentError("lambda_cl, T_cl and Hcl must be provided when element == 'cl'."))
        end
        return (1 - epsilon_cl(lambda_cl, T_cl, Hcl) - epsilon_mc(lambda_cl, T_cl, Hcl)) * sigma_e_cl

    else
        throw(ArgumentError("The element should be either 'gdl', 'mpl' or 'cl'."))
    end
end


"""This function calculates the thermal conductivity of fluids, in J.m-1.s-1.K-1, as a function of temperature.

Parameters
----------
component : String
    Specifies the gas for which the thermal conductivity is calculated.
    Must be either `"H2O_l"` (liquid water), `"H2O_v"` (vapor), `"H2"` (hydrogen), `"O2"` (oxygen), or `"N2"` (nitrogen).
T :
    Temperature in K.

Returns
-------
k_th
    Thermal conductivity of the selected fluid in J.m-1.s-1.K-1.
"""
function k_th(component::String, T)

    if component == "H2O_l"  # For T >= 273.16 and T <= 633.15 K.
        return -0.2987 + 4.7054e-3 * T - 5.6209e-6 * T^2
    elseif component == "H2O_v"  # For T >= 150 K and T <= 1500 K.
        return 5.6199e-3 + 1.5699e-5 * T + 1.0106e-7 * T^2 - 2.4282e-11 * T^3
    elseif component == "H2"  # For T >= 14 K and T <= 1500 K.
        return 1.0979e-2 + 6.6411e-4 * T - 3.4378e-7 * T^2 + 9.7283e-11 * T^3
    elseif component == "O2"  # For T >= 80 K and T <= 2000 K.
        return 1.5475e-4 + 9.4153e-5 * T - 2.7529e-8 * T^2 + 5.2069e-12 * T^3
    elseif component == "N2"  # For T >= 63 K and T <= 1500 K.
        return -2.2678e-4 + 1.0275e-4 * T - 6.0151e-8 * T^2 + 2.2332e-11 * T^3
    else
        throw(ArgumentError("The element should be either 'H2O_l', 'H2O_v', 'H2', 'O2' or 'N2'."))
    end
end


"""This function calculates the thermal conductivity of a gas mixture, in J.m-1.s-1.K-1.
The Lindsay-Bromley (Wassiljewa) method is used.

Parameters
----------
k_th_g : Vector
    Thermal conductivities of each pure gas component, in J.m-1.s-1.K-1, at the same temperature.
mu_g : Vector
    Viscosity of each pure gas component, in Pa.s, at the same temperature.
x : Vector
    Mole fractions of each gas component in the mixture (must sum to 1).
M : Vector
    Molar masses of each gas component (in kg.mol-1).

Returns
-------
k_th_gaz_mixture
    Thermal conductivity of the gas mixture, in J.m-1.s-1.K-1.
"""
function k_th_gaz_mixture(k_th_g::Vector,
                          mu_g::Vector,
                          x::Vector,
                          M::Vector)

    n = length(k_th_g)
    if length(mu_g) != n || length(x) != n || length(M) != n
        throw(ArgumentError("All input vectors must have the same length."))
    end

    total_x = 0.0
    @inbounds for xi in x
        total_x += xi
    end
    if abs(total_x - 1.0) > 1e-6
        throw(ArgumentError("The sum of the molar fractions should be 1."))
    end

    epsilon_TS = 0.85  # Value suggested by Tandon and Saxena in 1965.

    # Calculation of A_W using the Maxon and Saxena suggestion.
    A_W = Matrix{Float64}(undef, n, n)
    @inbounds for i in 1:n
        for j in 1:n
            if i == j
                A_W[i, j] = 1.0
            else
                A_W[i, j] = (epsilon_TS * (1 + sqrt(mu_g[i] / mu_g[j]) * (M[j] / M[i])^0.25)^2) /
                            sqrt(8 * (1 + M[i] / M[j]))
            end
        end
    end

    # Calculation of the thermal conductivity of the gas mixture.
    k_th_gaz_mixture = 0.0
    @inbounds for i in 1:n
        prod_x_A_w = 0.0
        for j in 1:n
            prod_x_A_w += x[j] * A_W[i, j]
        end
        k_th_gaz_mixture += x[i] * k_th_g[i] / prod_x_A_w
    end

    return k_th_gaz_mixture
end


"""This function calculates the effective thermal conductivity, in J.m-1.s-1.K-1, in either the GDL, the MPL,
the CL, or the membrane.

Parameters
----------
element : String
    Must be either `"agdl"`, `"cgdl"`, `"ampl"`, `"cmpl"`, `"acl"`, `"ccl"`, or `"mem"`.
T :
    Temperature in K.
C_v :
    Water concentration variable in mol.m-3.
s :
    Liquid water saturation variable.
lambdaa :
    Water content in the membrane.
C_H2 :
    Concentration of hydrogen in the AGDL or ACL.
C_O2 :
    Concentration of oxygen in the CGDL or CCL.
C_N2 :
    Concentration of nitrogen in the AGDL, ACL, CGDL or CCL.
epsilon : Union{Float64, Nothing}
    Porosity.
Hcl : Union{Float64, Nothing}
    Thickness of the CL layer.
epsilon_c : Union{Float64, Nothing}
    Compression ratio of the GDL.

Returns
-------
k_th_eff
    Effective thermal conductivity in J.m-1.s-1.K-1.
"""
function k_th_eff(element::String,
                  T,
                  C_v=nothing,
                  s=nothing,
                  lambdaa=nothing,
                  C_H2=nothing,
                  C_O2=nothing,
                  C_N2=nothing,
                  epsilon::Union{Float64, Nothing}=nothing,
                  Hcl::Union{Float64, Nothing}=nothing,
                  epsilon_c::Union{Float64, Nothing}=nothing)::Number

    if element == "agdl" || element == "cgdl"  # The effective thermal conductivity at the GDL.
        if C_v === nothing || s === nothing || C_N2 === nothing || epsilon === nothing || epsilon_c === nothing
            throw(ArgumentError("C_v, s, C_N2, epsilon and epsilon_c must be provided for GDL elements."))
        end
        beta3 = epsilon < 0.67 ? 4.04 : 4.40
        if element == "agdl"  # The thermal conductivity of the gas mixture in the AGDL.
            if C_H2 === nothing
                throw(ArgumentError("C_H2 must be provided for 'agdl'."))
            end
            sum_C = C_v + C_H2 + C_N2
            x_v, x_h2, x_n2 = C_v / sum_C, C_H2 / sum_C, C_N2 / sum_C
            k_th_gaz = k_th_gaz_mixture([k_th("H2O_v", T), k_th("H2", T), k_th("N2", T)],
                                        [mu_gaz("H2O_v", T), mu_gaz("H2", T), mu_gaz("N2", T)],
                                        [x_v, x_h2, x_n2],
                                        [M_H2O, M_H2, M_N2])
        else  # The thermal conductivity of the gas mixture in the CGDL.
            if C_O2 === nothing
                throw(ArgumentError("C_O2 must be provided for 'cgdl'."))
            end
            sum_C = C_v + C_O2 + C_N2
            x_v, x_o2, x_n2 = C_v / sum_C, C_O2 / sum_C, C_N2 / sum_C
            k_th_gaz = k_th_gaz_mixture([k_th("H2O_v", T), k_th("O2", T), k_th("N2", T)],
                                        [mu_gaz("H2O_v", T), mu_gaz("O2", T), mu_gaz("N2", T)],
                                        [x_v, x_o2, x_n2],
                                        [M_H2O, M_O2, M_N2])
        end
        return hmean([k_th_gdl * exp(beta3 * epsilon_c), k_th("H2O_l", T), k_th_gaz],
                     [1 - epsilon, epsilon * s, epsilon * (1 - s)])

    elseif element == "ampl" || element == "cmpl"  # The effective thermal conductivity at the MPL.
        if C_v === nothing || s === nothing || C_N2 === nothing || epsilon === nothing
            throw(ArgumentError("C_v, s, C_N2 and epsilon must be provided for MPL elements."))
        end
        if element == "ampl"  # The thermal conductivity of the gas mixture in the AMPL.
            if C_H2 === nothing
                throw(ArgumentError("C_H2 must be provided for 'ampl'."))
            end
            sum_C = C_v + C_H2 + C_N2
            x_v, x_h2, x_n2 = C_v / sum_C, C_H2 / sum_C, C_N2 / sum_C
            k_th_gaz = k_th_gaz_mixture([k_th("H2O_v", T), k_th("H2", T), k_th("N2", T)],
                                        [mu_gaz("H2O_v", T), mu_gaz("H2", T), mu_gaz("N2", T)],
                                        [x_v, x_h2, x_n2],
                                        [M_H2O, M_H2, M_N2])
        else  # The thermal conductivity of the gas mixture in the CMPL.
            if C_O2 === nothing
                throw(ArgumentError("C_O2 must be provided for 'cmpl'."))
            end
            sum_C = C_v + C_O2 + C_N2
            x_v, x_o2, x_n2 = C_v / sum_C, C_O2 / sum_C, C_N2 / sum_C
            k_th_gaz = k_th_gaz_mixture([k_th("H2O_v", T), k_th("O2", T), k_th("N2", T)],
                                        [mu_gaz("H2O_v", T), mu_gaz("O2", T), mu_gaz("N2", T)],
                                        [x_v, x_o2, x_n2],
                                        [M_H2O, M_O2, M_N2])
        end
        return hmean([k_th_mpl, k_th("H2O_l", T), k_th_gaz],
                     [1 - epsilon, epsilon * s, epsilon * (1 - s)])

    elseif element == "acl" || element == "ccl"  # The effective thermal conductivity at the CL.
        if C_v === nothing || s === nothing || lambdaa === nothing || C_N2 === nothing || Hcl === nothing
            throw(ArgumentError("C_v, s, lambdaa, C_N2 and Hcl must be provided for CL elements."))
        end

        fv_val = fv(lambdaa, T)
        epsilon_mc_val = epsilon_mc(lambdaa, T, Hcl)
        epsilon_cl_val = epsilon_cl(lambdaa, T, Hcl)
        k_th_eff_mem = hmean([k_th_mem, k_th("H2O_l", T)], [1 - fv_val, fv_val])

        if element == "acl"  # The thermal conductivity of the gas mixture in the ACL.
            if C_H2 === nothing
                throw(ArgumentError("C_H2 must be provided for 'acl'."))
            end
            sum_C = C_v + C_H2 + C_N2
            x_v, x_h2, x_n2 = C_v / sum_C, C_H2 / sum_C, C_N2 / sum_C
            k_th_gaz = k_th_gaz_mixture([k_th("H2O_v", T), k_th("H2", T), k_th("N2", T)],
                                        [mu_gaz("H2O_v", T), mu_gaz("H2", T), mu_gaz("N2", T)],
                                        [x_v, x_h2, x_n2],
                                        [M_H2O, M_H2, M_N2])
        else  # The thermal conductivity of the gas mixture in the CCL.
            if C_O2 === nothing
                throw(ArgumentError("C_O2 must be provided for 'ccl'."))
            end
            sum_C = C_v + C_O2 + C_N2
            x_v, x_o2, x_n2 = C_v / sum_C, C_O2 / sum_C, C_N2 / sum_C
            k_th_gaz = k_th_gaz_mixture([k_th("H2O_v", T), k_th("O2", T), k_th("N2", T)],
                                        [mu_gaz("H2O_v", T), mu_gaz("O2", T), mu_gaz("N2", T)],
                                        [x_v, x_o2, x_n2],
                                        [M_H2O, M_O2, M_N2])
        end

        return hmean([k_th_cl, k_th_eff_mem, k_th("H2O_l", T), k_th_gaz],
                     [1 - epsilon_cl_val - epsilon_mc_val, epsilon_mc_val, epsilon_cl_val * s, epsilon_cl_val * (1 - s)])

    elseif element == "mem"  # The effective thermal conductivity at the membrane.
        if lambdaa === nothing
            throw(ArgumentError("lambdaa must be provided for 'mem'."))
        end
        fv_val = fv(lambdaa, T)
        return hmean([k_th_mem, k_th("H2O_l", T)], [1 - fv_val, fv_val])

    else
        throw(ArgumentError("The element should be either 'agdl', 'cgdl', 'ampl', 'cmpl', 'acl', 'ccl' or 'mem'."))
    end
end


"""This function calculates the specific heat capacity of fluids, in J.kg-1.K-1, as a function of temperature.

Parameters
----------
component : String
    Must be either `"H2O_l"` (liquid water), `"H2O_v"` (vapor), `"H2"` (hydrogen), `"O2"` (oxygen), or `"N2"` (nitrogen).
T :
    Temperature in K.

Returns
-------
Cp0
    Specific heat capacity of the selected fluid in J.kg-1.K-1.
"""
function Cp0(component::String, T)

    if component == "H2O_l"  # For T >= 298 and T <= 500 K.
        return 1 / M_H2O * (-203.6060 + 1523.290 * (T / 1000) - 3196.413 * (T / 1000)^2 + 2474.455 * (T / 1000)^3 +
                            3.855326 / (T / 1000)^2)
    elseif component == "H2O_v"  # For T = 350 K.
        return 1880
    elseif component == "H2"  # For T >= 298 K and T <= 1000 K.
        return 1 / M_H2 * (33.066178 - 11.363417 * (T / 1000) + 11.432816 * (T / 1000)^2 - 2.772874 * (T / 1000)^3 -
                           0.158558 / (T / 1000)^2)
    elseif component == "O2"  # For T >= 100 K and T <= 700 K.
        return 1 / M_O2 * (31.32234 - 20.23531 * (T / 1000) + 57.86644 * (T / 1000)^2 - 36.50624 * (T / 1000)^3 -
                           0.007374 / (T / 1000)^2)
    elseif component == "N2"  # For T >= 100 K and T <= 500 K.
        return 1 / M_N2 * (28.98641 + 1.853978 * (T / 1000) - 9.647459 * (T / 1000)^2 + 16.63537 * (T / 1000)^3 +
                           0.000117 / (T / 1000)^2)
    else
        throw(ArgumentError("The element should be either 'H2O_l', 'H2O_v', 'H2', 'O2' or 'N2'."))
    end
end


"""This function calculates the standard enthalpy of fluids, in J.mol-1, as a function of temperature.

Parameters
----------
component : String
    Must be either `"H2O_l"` (liquid water) or `"H2O_v"` (vapor).
T :
    Temperature in K.

Returns
-------
h0
    Standard enthalpy of the selected fluid in J.mol-1.
"""
function h0(component::String, T)

    if component == "H2O_l"  # For T >= 298 and T <= 500 K.
        return (-285.83 - 203.6060 * (T / 1000) + 1523.290 * (T / 1000)^2 / 2 - 3196.413 * (T / 1000)^3 / 3 +
                2474.455 * (T / 1000)^4 / 4 - 3.855326 / (T / 1000) - 256.5478 + 285.8304) * 1e3
    elseif component == "H2O_v"  # For T = 298.15 K.
        return -241.83 * 1e3 + Cp0("H2O_v", T) * M_H2O * (T - 298.15)
    else
        throw(ArgumentError("The element should be either 'H2O_l' or 'H2O_v'."))
    end
end


"""This function calculates the volumetric heat capacity, in J.m-3.K-1, in either the GDL, the MPL, the CL or
the membrane.

Parameters
----------
element : String
    Must be either `"agdl"`, `"cgdl"`, `"ampl"`, `"cmpl"`, `"acl"`, `"ccl"`, or `"mem"`.
T :
    Temperature in K.
C_v :
    Water concentration variable in mol.m-3.
s :
    Liquid water saturation variable.
lambdaa :
    Water content in the membrane.
C_H2 :
    Concentration of hydrogen in the AGDL or ACL.
C_O2 :
    Concentration of oxygen in the CGDL or CCL.
C_N2 :
    Concentration of nitrogen in the AGDL, ACL, CGDL or CCL.
epsilon : Union{Float64, Nothing}
    Porosity.
Hcl : Union{Float64, Nothing}
    Thickness of the CL layer.

Returns
-------
rho_Cp0
    Volumetric heat capacity in J.m-3.K-1.
"""
function calculate_rho_Cp0(element::String,
                           T,
                           C_v=nothing,
                           s=nothing,
                           lambdaa=nothing,
                           C_H2=nothing,
                           C_O2=nothing,
                           C_N2=nothing,
                           epsilon::Union{Float64, Nothing}=nothing,
                           Hcl::Union{Float64, Nothing}=nothing)

    if element == "agdl" || element == "cgdl" || element == "ampl" || element == "cmpl"
        if C_v === nothing || s === nothing || C_N2 === nothing || epsilon === nothing
            throw(ArgumentError("C_v, s, C_N2 and epsilon must be provided for GDL/MPL elements."))
        end

        if element == "agdl" || element == "ampl"  # In the anode.
            if C_H2 === nothing
                throw(ArgumentError("C_H2 must be provided for anode elements."))
            end
            sum_C = C_v + C_H2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0("H2O_v", T), M_H2 * C_H2 * Cp0("H2", T), M_N2 * C_N2 * Cp0("N2", T)],
                                  [C_v / sum_C, C_H2 / sum_C, C_N2 / sum_C])
        else  # In the cathode.
            if C_O2 === nothing
                throw(ArgumentError("C_O2 must be provided for cathode elements."))
            end
            sum_C = C_v + C_O2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0("H2O_v", T), M_O2 * C_O2 * Cp0("O2", T), M_N2 * C_N2 * Cp0("N2", T)],
                                  [C_v / sum_C, C_O2 / sum_C, C_N2 / sum_C])
        end

        if element == "agdl" || element == "cgdl"  # In the GDLs.
            return average([rho_gdl * Cp_gdl, rho_H2O_l(T) * Cp0("H2O_l", T), rho_Cp0_gaz],
                           [1 - epsilon, epsilon * s, epsilon * (1 - s)])
        else  # In the MPLs.
            return average([rho_mpl * Cp_mpl, rho_H2O_l(T) * Cp0("H2O_l", T), rho_Cp0_gaz],
                           [1 - epsilon, epsilon * s, epsilon * (1 - s)])
        end

    elseif element == "acl" || element == "ccl"  # The volumetric heat capacity at the CL.
        if C_v === nothing || s === nothing || lambdaa === nothing || C_N2 === nothing || Hcl === nothing
            throw(ArgumentError("C_v, s, lambdaa, C_N2 and Hcl must be provided for CL elements."))
        end

        epsilon_mc_val = epsilon_mc(lambdaa, T, Hcl)
        epsilon_cl_val = epsilon_cl(lambdaa, T, Hcl)

        if element == "acl"  # The heat capacity of the gas mixture in the ACL.
            if C_H2 === nothing
                throw(ArgumentError("C_H2 must be provided for 'acl'."))
            end
            sum_C = C_v + C_H2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0("H2O_v", T), M_H2 * C_H2 * Cp0("H2", T), M_N2 * C_N2 * Cp0("N2", T)],
                                  [C_v / sum_C, C_H2 / sum_C, C_N2 / sum_C])
        else  # The heat capacity of the gas mixture in the CCL.
            if C_O2 === nothing
                throw(ArgumentError("C_O2 must be provided for 'ccl'."))
            end
            sum_C = C_v + C_O2 + C_N2
            rho_Cp0_gaz = average([M_H2O * C_v * Cp0("H2O_v", T), M_O2 * C_O2 * Cp0("O2", T), M_N2 * C_N2 * Cp0("N2", T)],
                                  [C_v / sum_C, C_O2 / sum_C, C_N2 / sum_C])
        end

        return average([rho_cl * Cp_cl, rho_mem * Cp_mem, rho_H2O_l(T) * Cp0("H2O_l", T), rho_Cp0_gaz],
                       [1 - epsilon_cl_val - epsilon_mc_val, epsilon_mc_val, epsilon_cl_val * s, epsilon_cl_val * (1 - s)])

    elseif element == "mem"  # The volumetric heat capacity at the membrane.
        if lambdaa === nothing
            throw(ArgumentError("lambdaa must be provided for 'mem'."))
        end
        fv_val = fv(lambdaa, T)
        return average([rho_mem * Cp_mem, rho_H2O_l(T) * Cp0("H2O_l", T)], [1 - fv_val, fv_val])

    else
        throw(ArgumentError("The element should be either 'agdl', 'cgdl', 'ampl', 'cmpl', 'acl', 'ccl' or 'mem'."))
    end
end


"""This function computes the molar enthalpy of liquefaction of water at a given temperature, in J.mol-1.
It is calculated as the difference in molar enthalpy between liquid water (H2O_l) and water vapor (H2O_v).

Parameters
----------
T :
    Temperature in K.

Returns
-------
delta_h_liq
    Molar enthalpy of liquefaction in J.mol-1.
"""
function delta_h_liq(T)
    return h0("H2O_l", T) - h0("H2O_v", T)
end


"""This function computes the molar enthalpy of absorption of water at a given temperature, in J.mol-1.
This reaction is exothermic.

Parameters
----------
T :
    Temperature in K.

Returns
-------
delta_h_abs
    Molar enthalpy of absorption in the CL in J.mol-1.
"""
function delta_h_abs(T)
    return delta_h_liq(T)
end

