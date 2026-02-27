# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the flows calculation.
"""
import math
from functools import lru_cache

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from alphapem.utils.physics_constants import (R, M_H2, M_O2, M_N2, M_H2O, theta_c_gdl, theta_c_mpl, theta_c_cl,
                                              epsilon_p, alpha_p, r_s_gdl, tau_mpl, r_s_mpl, tau_cl, r_s_cl, Kshape, F,
                                              M_eq, rho_mem, gamma_cond, gamma_evap, Dp_mpl, Dp_cl, Eact_H2_cros_v,
                                              Tref_cross, Eact_H2_cros_l, Eact_O2_cros_v, Eact_O2_cros_l)
from alphapem.utils.maths_functions import (hmean, average)
from alphapem.core.modules.cell_voltage_modules import epsilon_mc, epsilon_cl
from alphapem.utils.physics_functions import nu_l, C_v_sat, rho_H2O_l, Psat


# _____________________________________________________Flow modules_____________________________________________________

def flows_1D_MEA_int_values(sv, i_fc, parameters):
    """This functions calculates intermediate values for the flows calculation.

    Parameters
    ----------
    sv : dict
        Variables calculated by the solver. They correspond to the fuel cell internal states.
        sv is a contraction of solver_variables for enhanced readability.
    i_fc : float
        Current density of the fuel cell (A/m²).
    parameters : dict
        Parameters of the fuel cell model.

    Returns
    -------
    Pagc : float
        Global pressure in the anode gas channel (Pa).
    Pcgc : float
        Global pressure in the cathode gas channel (Pa).
    lambda_acl_mem : float
        Water content in the ACL and the membrane (kg/kg).
    """

    # Extraction of the variables
    C_v_agc, C_v_acl, C_v_ccl, C_v_cgc = sv['C_v_agc'], sv['C_v_acl'], sv['C_v_ccl'], sv['C_v_cgc']
    s_acl, s_ccl = sv['s_acl'], sv['s_ccl']
    lambda_acl, lambda_mem, lambda_ccl = sv['lambda_acl'], sv['lambda_mem'], sv['lambda_ccl']
    C_H2_agc, C_H2_acl, C_O2_ccl, C_O2_cgc = sv['C_H2_agc'], sv['C_H2_acl'], sv['C_O2_ccl'], sv['C_O2_cgc']
    C_N2_agc, C_N2_cgc = sv['C_N2_agc'], sv['C_N2_cgc']
    T_agc, T_acl, T_mem, T_ccl, T_cgc = sv['T_agc'], sv['T_acl'], sv['T_mem'], sv['T_ccl'], sv['T_cgc']
    # Extraction of the operating inputs and the parameters
    epsilon_gdl, epsilon_mpl, epsilon_c = parameters['epsilon_gdl'], parameters['epsilon_mpl'], parameters['epsilon_c']
    e, Hacl, Hccl, Hmem = parameters['e'], parameters['Hacl'], parameters['Hccl'], parameters['Hmem']
    Hgdl, Hmpl, Wagc, Wcgc = parameters['Hgdl'], parameters['Hmpl'], parameters['Wagc'], parameters['Wcgc']
    nb_gdl, nb_mpl = parameters['nb_gdl'], parameters['nb_mpl']

    # Transitory parameter
    H_gdl_node = Hgdl / nb_gdl
    H_mpl_node = Hmpl / nb_mpl

    # Pressures in the stack
    Pagc = (C_v_agc + C_H2_agc + C_N2_agc) * R * T_agc
    Pagdl = [None] + [(sv[f'C_v_agdl_{i}'] + sv[f'C_H2_agdl_{i}'] + C_N2_agc) * R * sv[f'T_agdl_{i}'] for i in range(1, nb_gdl + 1)]
    Pampl = [None] + [(sv[f'C_v_ampl_{i}'] + sv[f'C_H2_ampl_{i}'] + C_N2_agc) * R * sv[f'T_ampl_{i}'] for i in range(1, nb_mpl + 1)]
    Pacl = (C_v_acl + C_H2_acl + C_N2_agc) * R * T_acl
    Pccl = (C_v_ccl + C_O2_ccl + C_N2_cgc) * R * T_ccl
    Pcmpl = [None] + [(sv[f'C_v_cmpl_{i}'] + sv[f'C_O2_cmpl_{i}'] + C_N2_cgc) * R * sv[f'T_cmpl_{i}'] for i in range(1, nb_mpl + 1)]
    Pcgdl = [None] + [(sv[f'C_v_cgdl_{i}'] + sv[f'C_O2_cgdl_{i}'] + C_N2_cgc) * R * sv[f'T_cgdl_{i}'] for i in range(1, nb_gdl + 1)]
    Pcgc = (C_v_cgc + C_O2_cgc + C_N2_cgc) * R * T_cgc

    # Capillary pressures in the stack
    Pcap_agdl = Pcap('gdl', sv['s_agdl_1'], sv['T_agdl_1'], epsilon_gdl, epsilon_c=epsilon_c)
    Pcap_cgdl = Pcap('gdl', sv[f's_cgdl_{nb_gdl}'], sv[f'T_cgdl_{nb_gdl}'], epsilon_gdl, epsilon_c=epsilon_c)

    # Densities in the GC
    rho_agc = C_H2_agc * M_H2 + C_v_agc * M_H2O + C_N2_agc * M_N2
    rho_cgc = C_O2_cgc * M_O2 + C_v_cgc * M_H2O + C_N2_cgc * M_N2

    # Weighted mean values ...
    #       ... of the EOD flow of water in the membrane
    D_eff_EOD_acl_mem = hmean([D_EOD_eff(i_fc, lambda_acl, T_acl, Hacl), D_EOD(lambda_mem)],
                          weights = [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_eff_EOD_mem_ccl = hmean([D_EOD(lambda_mem), D_EOD_eff(i_fc, lambda_ccl, T_ccl, Hccl)],
                          weights = [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])
    #       ... of the diffusion coefficient of water in the membrane
    D_lambda_eff_acl_mem = hmean([D_lambda_eff(lambda_acl, T_acl, Hacl), D_lambda(lambda_mem)],
                          weights = [Hacl / (Hacl + Hmem), Hmem / (Hacl + Hmem)])
    D_lambda_eff_mem_ccl = hmean([D_lambda(lambda_mem), D_lambda_eff(lambda_ccl, T_ccl, Hccl)],
                          weights = [Hmem / (Hmem + Hccl), Hccl / (Hmem + Hccl)])
    #       ... of the capillary coefficient
    D_cap_agdl_agdl = [None] + [hmean([Dcap('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                         Dcap('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, nb_gdl)]
    D_cap_agdl_ampl = hmean([Dcap('gdl', sv[f's_agdl_{nb_gdl}'], sv[f'T_agdl_{nb_gdl}'], epsilon_gdl, e,
                                   epsilon_c=epsilon_c),
                                  Dcap('mpl', sv['s_ampl_1'], sv['T_ampl_1'], epsilon_mpl, e)],
                             weights=[H_gdl_node / (H_gdl_node + H_mpl_node), H_mpl_node / (H_gdl_node + H_mpl_node)])
    D_cap_ampl_ampl = [None] + [hmean([Dcap('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], epsilon_mpl, e),
                                         Dcap('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], epsilon_mpl, e)])
                                for i in range(1, nb_mpl)]
    D_cap_ampl_acl = hmean([Dcap('mpl', sv[f's_ampl_{nb_mpl}'], sv[f'T_ampl_{nb_mpl}'], epsilon_mpl, e),
                                    Dcap('cl', s_acl, T_acl, epsilon_cl(lambda_acl, T_acl, Hacl), e)],
                             weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    D_cap_ccl_cmpl = hmean([Dcap('cl', s_ccl, T_ccl, epsilon_cl(lambda_ccl, T_ccl, Hccl), e),
                            Dcap('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], epsilon_mpl, e)],
                           weights=[Hccl / (Hccl + H_mpl_node), H_mpl_node / (Hccl + H_mpl_node)])
    D_cap_cmpl_cmpl = [None] + [hmean([Dcap('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], epsilon_mpl, e),
                                       Dcap('mpl', sv[f's_cmpl_{i + 1}'], sv[f'T_cmpl_{i + 1}'], epsilon_mpl, e)])
                                for i in range(1, nb_mpl)]
    D_cap_cmpl_cgdl = hmean([Dcap('mpl', sv[f's_cmpl_{nb_mpl}'], sv[f'T_cmpl_{nb_mpl}'], epsilon_mpl, e),
                             Dcap('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], epsilon_gdl, e,
                                  epsilon_c=epsilon_c)],
                            weights=[H_mpl_node / (H_mpl_node + H_gdl_node), H_gdl_node / (H_mpl_node + H_gdl_node)])
    D_cap_cgdl_cgdl = [None] + [hmean([Dcap('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c),
                                        Dcap('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], epsilon_gdl, e,
                                            epsilon_c=epsilon_c)]) for i in range(1, nb_gdl)]
    #       ... of the effective diffusion coefficient
    Da_eff_agdl_agdl = [None] + [hmean([Da_eff('gdl', sv[f's_agdl_{i}'], sv[f'T_agdl_{i}'], Pagdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                          Da_eff('gdl', sv[f's_agdl_{i+1}'], sv[f'T_agdl_{i+1}'], Pagdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, nb_gdl)]
    Da_eff_agdl_ampl = hmean([Da_eff('gdl', sv[f's_agdl_{nb_gdl}'], sv[f'T_agdl_{nb_gdl}'], Pagdl[nb_gdl],
                                          epsilon_gdl, epsilon_c = epsilon_c),
                                      Da_eff('mpl', sv['s_ampl_1'], sv['T_ampl_1'], Pampl[1], epsilon_mpl)],
                               weights = [H_gdl_node / (H_gdl_node + H_mpl_node), H_mpl_node / (H_gdl_node + H_mpl_node)])
    Da_eff_ampl_ampl = [None] + [hmean([Da_eff('mpl', sv[f's_ampl_{i}'], sv[f'T_ampl_{i}'], Pampl[i],
                                              epsilon_mpl),
                                        Da_eff('mpl', sv[f's_ampl_{i+1}'], sv[f'T_ampl_{i+1}'], Pampl[i+1],
                                              epsilon_mpl)]) for i in range(1, nb_mpl)]
    Da_eff_ampl_acl = hmean([Da_eff('mpl', sv[f's_ampl_{nb_mpl}'], sv[f'T_ampl_{nb_mpl}'], Pampl[nb_mpl],
                                            epsilon_mpl),
                                     Da_eff('cl', s_acl, T_acl, Pacl, epsilon_cl(lambda_acl, T_acl, Hacl))],
                              weights=[H_mpl_node / (H_mpl_node + Hacl), Hacl / (H_mpl_node + Hacl)])
    Dc_eff_ccl_cmpl = hmean([Dc_eff('cl', s_ccl, T_ccl, Pccl, epsilon_cl(lambda_ccl, T_ccl, Hccl)),
                               Dc_eff('mpl', sv['s_cmpl_1'], sv['T_cmpl_1'], Pcmpl[1], epsilon_mpl)],
                              weights=[Hccl / (H_mpl_node + Hccl), H_mpl_node / (H_mpl_node + Hccl)])
    Dc_eff_cmpl_cmpl = [None] + [hmean([Dc_eff('mpl', sv[f's_cmpl_{i}'], sv[f'T_cmpl_{i}'], Pcmpl[i],
                                                 epsilon_mpl),
                                          Dc_eff('mpl', sv[f's_cmpl_{i + 1}'], sv[f'T_cmpl_{i + 1}'], Pcmpl[i + 1],
                                                 epsilon_mpl)]) for i in range(1, nb_mpl)]
    Dc_eff_cmpl_cgdl = hmean([Dc_eff('mpl', sv[f's_cmpl_{nb_mpl}'], sv[f'T_cmpl_{nb_mpl}'], Pcmpl[nb_mpl],
                                       epsilon_mpl),
                                    Dc_eff('gdl', sv['s_cgdl_1'], sv['T_cgdl_1'], Pcgdl[1], epsilon_gdl,
                                           epsilon_c=epsilon_c)],
                               weights=[H_mpl_node / (H_mpl_node + H_gdl_node), H_gdl_node / (H_mpl_node + H_gdl_node)])
    Dc_eff_cgdl_cgdl = [None] + [hmean([Dc_eff('gdl', sv[f's_cgdl_{i}'], sv[f'T_cgdl_{i}'], Pcgdl[i],
                                              epsilon_gdl, epsilon_c = epsilon_c),
                                        Dc_eff('gdl', sv[f's_cgdl_{i+1}'], sv[f'T_cgdl_{i+1}'], Pcgdl[i+1],
                                              epsilon_gdl, epsilon_c = epsilon_c)]) for i in range(1, nb_gdl)]
    #       ... of the temperature
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                        weights=[Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])

    return (H_gdl_node, H_mpl_node, Pagc, Pcgc, Pcap_agdl, Pcap_cgdl, rho_agc, rho_cgc, D_eff_EOD_acl_mem,
            D_eff_EOD_mem_ccl, D_lambda_eff_acl_mem, D_lambda_eff_mem_ccl, D_cap_agdl_agdl, D_cap_agdl_ampl,
            D_cap_ampl_ampl, D_cap_ampl_acl, D_cap_ccl_cmpl, D_cap_cmpl_cmpl, D_cap_cmpl_cgdl, D_cap_cgdl_cgdl,
            Da_eff_agdl_agdl, Da_eff_agdl_ampl, Da_eff_ampl_ampl, Da_eff_ampl_acl, Dc_eff_ccl_cmpl, Dc_eff_cmpl_cmpl,
            Dc_eff_cmpl_cgdl, Dc_eff_cgdl_cgdl, T_acl_mem_ccl)


def Dcap(element, s, T, epsilon, e, epsilon_c=None):
    """ This function calculates the capillary coefficient at the GDL, the MPL or the CL, in kg.m.s-1, considering
     GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the capillary coefficient is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.
    epsilon : float
        Porosity.
    e : float
        Capillary exponent.
    epsilon_c : float, optional
        Compression ratio of the GDL.

    Returns
    -------
    float
        Capillary coefficient at the GDL, MPL or CL in kg.m.s-1.
    """

    K0_value = K0(element, epsilon, epsilon_c)
    if element == 'gdl':
        theta_c_value = theta_c_gdl
    elif element == 'mpl':
        theta_c_value = theta_c_mpl
    elif element == 'cl':
        theta_c_value = theta_c_cl
    else:
        raise ValueError("The element should be either 'gdl', 'mpl' or 'cl'.")

    return sigma(T) * K0_value / nu_l(T) * abs(math.cos(theta_c_value)) * \
           (epsilon / K0_value) ** 0.5 * (s ** e + 1e-7) * (1.417 - 4.24 * s + 3.789 * s ** 2)


def Pcap(element, s, T, epsilon, epsilon_c=None):
    """ This function calculates the capillary pressure at the GDL, the MPL or the CL, in kg.m.s-1.

    Parameters
    ----------
    element : str
        Specifies the element for which the capillary coefficient is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.

    Returns
    -------
    float
        Capillary pressure at the selected element in kg.m.s-1.
    """

    K0_value = K0(element, epsilon, epsilon_c)
    if element == 'gdl':
        theta_c_value = theta_c_gdl
    elif element == 'mpl':
        theta_c_value = theta_c_mpl
    elif element == 'cl':
        theta_c_value = theta_c_cl
    else:
        raise ValueError("The element should be either 'gdl', 'mpl' or 'cl'.")

    s_num = s + 1e-7 # To avoid numerical issues when s = 0.

    return sigma(T) * abs(math.cos(theta_c_value)) * (epsilon / K0_value) ** 0.5 * \
           (1.417 * s_num - 2.12 * s_num ** 2 + 1.263 * s_num ** 3)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def Da(P, T):
    """This function calculates the diffusion coefficient at the anode, in m².s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Diffusion coefficient at the anode in m².s-1.
    """
    return 1.644e-4 * (T / 333) ** 2.334 * (101325 / P)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def Dc(P, T):
    """This function calculates the diffusion coefficient at the cathode, in m².s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Diffusion coefficient at the cathode in m².s-1.
    """
    return 3.242e-5 * (T / 333) ** 2.334 * (101325 / P)


def Da_eff(element, s, T, P, epsilon, epsilon_c=None):
    """This function calculates the effective diffusion coefficient at the GDL, TL, MPL or the CL and at the anode,
     in m².s-1, considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the effective diffusion coefficient is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.

    Returns
    -------
    float
        Effective diffusion coefficient at the anode in m².s-1.
    """

    if element == 'gdl': # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta2 = -1.59
        else:
            beta2 = -0.90
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p)) ** alpha_p)
        return epsilon / tau_gdl * math.exp(beta2 * epsilon_c) * (1 - s) ** r_s_gdl * Da(P, T)

    elif element == 'mpl': # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s) ** r_s_mpl * Da(P, T)

    elif element == 'cl': # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_cl * (1 - s) ** r_s_cl * Da(P, T)

    else:
        raise ValueError("The element should be either 'gdl', 'tl', 'mpl' or 'cl'.")


def Dc_eff(element, s, T, P, epsilon, epsilon_c=None):
    """This function calculates the effective diffusion coefficient at the GDL, MPL, TL or the CL and at the cathode,
     in m².s-1, considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the effective diffusion coefficient is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.
    P : float
        Pressure in Pa.
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.

    Returns
    -------
    float
        Effective diffusion coefficient at the cathode in m².s-1.
    """

    if element == 'gdl': # The effective diffusion coefficient at the GDL using Tomadakis and Sotirchos model.
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta2 = -1.59
        else:
            beta2 = -0.90
        tau_gdl = 1 / (((epsilon - epsilon_p) / (1 - epsilon_p)) ** alpha_p)
        return epsilon / tau_gdl * math.exp(beta2 * epsilon_c) * (1 - s) ** r_s_gdl * Dc(P, T)

    elif element == 'mpl': # The effective diffusion coefficient at the MPL using Bruggeman model.
        return epsilon / tau_mpl * (1 - s) ** r_s_mpl * Dc(P, T)

    elif element == 'cl': # The effective diffusion coefficient at the CL using Bruggeman model.
        return epsilon / tau_cl * (1 - s) ** r_s_cl * Dc(P, T)

    else:
        raise ValueError("The element should be either 'gdl', 'tl', 'mpl' or 'cl'.")


def h_a(P, T, Wgc, Hgc):
    """This function calculates the effective convective-conductive mass transfer coefficient at the anode, in m.s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    Wgc : float
        Width of the gas channel in m.
    Hgc : float
        Thickness of the gas channel in m.

    Returns
    -------
    float
        Effective convective-conductive mass transfer coefficient at the anode in m.s-1.
    """
    Sh = 0.9247 * math.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.
    return Sh * Da(P, T) / Hgc


def h_c(P, T, Wgc, Hgc):
    """This function calculates the effective convective-conductive mass transfer coefficient at the cathode, in m.s-1.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T : float
        Temperature in K.
    Wgc : float
        Width of the gas channel in m.
    Hgc : float
        Thickness of the gas channel in m.

    Returns
    -------
    float
        Effective convective-conductive mass transfer coefficient at the cathode in m.s-1.
    """
    Sh = 0.9247 * math.log(Wgc / Hgc) + 2.3787  # Sherwood coefficient.
    return Sh * Dc(P, T) / Hgc


def lambda_v_eq(a_w):
    """This function calculates the equilibrium water content in the membrane from the vapor phase. Hinatsu's expression
     has been selected.

    Parameters
    ----------
    a_w: float
        Water activity.

    Returns
    -------
    float
        Equilibrium water content in the membrane from the vapor phase.
    """
    return 0.300 + 10.8 * a_w - 16.0 * a_w ** 2 + 14.1 * a_w ** 3


def lambda_l_eq(T):
    """This function calculates the equilibrium water content in the membrane from the liquid phase.
    Hinatsu's expression has been selected. It is valid for N-form membranes for 25 to 100 °C.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Equilibrium water content in the membrane from the liquid phase.
    """
    return 10.0 * 1.84e-2 * (T - 273.15) + 9.90e-4 * (T - 273.15)**2


def lambda_eq(C_v, s, T):
    """This function calculates the equilibrium water content in the membrane. Hinatsu's expression modified with
    Bao's formulation has been selected.

    Parameters
    ----------
    C_v : float
        Water concentration variable in mol.m-3.
    s : float
        Liquid water saturation variable.
    T : float
        Temperature in K.

    Returns
    -------
    float
        Equilibrium water content in the membrane.
    """
    a_w = C_v / C_v_sat(T) + 2 * s  # water activity
    return 0.5 * lambda_v_eq(a_w)                                          * (1 - math.tanh(100 * (a_w - 1))) + \
           0.5 * (lambda_v_eq(1) + ((lambda_l_eq(T) - lambda_v_eq(1)) / 2) * (1 - math.exp(-Kshape * (a_w - 1)))) * \
                                                                             (1 + math.tanh(100 * (a_w - 1)))


def D_lambda(lambdaa):
    """This function calculates the diffusion coefficient of water in the bulk membrane, in m².s-1.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.

    Returns
    -------
    float
        Diffusion coefficient of water in the membrane in m².s-1.
    """
    return 4.1e-10 * (lambdaa / 25.0) ** 0.15 * (1.0 + math.tanh((lambdaa - 2.5) / 1.4))


def D_lambda_eff(lambdaa, T, Hcl):
    """This function calculates the effective diffusion coefficient of water in the ionomer of the catalyst layers,
    in m².s-1.

    Parameters
    ----------
    lambdaa : float
        Water content in the catalyst layer.
    T : float
        Temperature in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Effective diffusion coefficient of water in the catalyst layer in m².s-1.
    """
    return epsilon_mc(lambdaa, T, Hcl) / tau_cl * D_lambda(lambdaa)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def D_EOD(i_fc):
    """This function calculates the electro-osmotic drag diffusion coefficient of water in the membrane, in mol.m-2.s-1.

    Parameters
    ----------
    i_fc : float
        Fuel cell current density in A.m-2.

    Returns
    -------
    float
        Electro-osmotic drag diffusion coefficient of water in the membrane in mol.m-2.s-1.
    """
    return 2.5 / 22 * i_fc / F


def D_EOD_eff(i_fc, lambdaa, T, Hcl):
    """This function calculates the effective electro-osmotic drag diffusion coefficient of water in the ionomer of the
    catalyst layers, in mol.m-2.s-1.

    Parameters
    ----------
    i_fc : float
        Fuel cell current density in A.m-2.
    lambdaa : float
        Water content in the catalyst layer.
    T : float
        Temperature in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Effective electro-osmotic drag diffusion coefficient of water in the catalyst layer in mol.m-2.s-1.
    """
    return epsilon_mc(lambdaa, T, Hcl) / tau_cl * D_EOD(i_fc)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def fv(lambdaa, T):
    """This function calculates the water volume fraction of the membrane.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.

    Returns
    -------
    float
        Water volume fraction of the membrane.
    """

    return (lambdaa * M_H2O / rho_H2O_l(T)) / (M_eq / rho_mem + lambdaa * M_H2O / rho_H2O_l(T))


def gamma_sorp(C_v, s, lambdaa, T, Hcl):
    """This function calculates the sorption rate of water in the membrane, in s-1.

    Parameters
    ----------
    C_v : float
        Water concentration variable in mol.m-3.
    s : float
        Liquid water saturation variable.
    lambdaa : float
        Water content in the membrane.
    T : float
        Temperature in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Sorption rate of water in the membrane in s-1.
    """

    fv_value = fv(lambdaa, T)
    gamma_abs = (1.14e-5 * fv_value) / Hcl * math.exp(2416 * (1 / 303 - 1 / T))
    gamma_des = (4.59e-5 * fv_value) / Hcl * math.exp(2416 * (1 / 303 - 1 / T))

    # Transition function between absorption and desorption
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_eq(C_v, s, T) - lambdaa))) # transition function

    return w * gamma_abs + (1 - w) * gamma_des # interpolation between absorption and desorption


@lru_cache(maxsize=None) # Cache the results to optimize performance
def Svl(element, s, C_v, Ctot, T, epsilon):
    """This function calculates the phase transfer rate of water condensation or evaporation, in mol.m-3.s-1.
    It is positive for condensation and negative for evaporation.

    Parameters
    ----------
    element : str
        Specifies the element for which the phase transfer rate is calculated.
    s : float
        Liquid water saturation variable.
    C_v : float
        Water concentration variable in mol.m-3.
    Ctot : float
        Total gas concentration in mol.m-3.
    T : float
        Temperature in K.
    epsilon : float
        Porosity.

    Returns
    -------
    float
        Phase transfer rate of water condensation or evaporation in mol.m-3.s-1.
    """

    # Calculation of the total and partial pressures
    Ptot = Ctot * R * T # Total pressure
    P_v = C_v * R * T # Partial pressure of vapor

    # Determination of the diffusion coefficient at the anode or the cathode
    if element == 'anode':
        D_value = Da(Ptot, T)  # Diffusion coefficient at the anode
    else:  # element == 'cathode'
        D_value = Dc(Ptot, T)  # Diffusion coefficient at the cathode

    Svl_cond = gamma_cond * M_H2O / (R * T) * epsilon * (1 - s) * D_value * Ptot * math.log((Ptot - Psat(T)) / (Ptot - P_v))
    Svl_evap = gamma_evap * M_H2O / (R * T) * epsilon * s * D_value * Ptot * math.log((Ptot - Psat(T)) / (Ptot - P_v))

    # Transition function between condensation and evaporation
    K_transition = 3e-3 # This is a constant that defines the sharpness of the transition between two states.  The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (Psat(T) - P_v))) # transition function

    return w * Svl_evap + (1 - w) * Svl_cond # interpolation between condensation and evaporation


def sigma(T):
    """This function calculates the water surface tension, in N.m-1, as a function of the temperature.

    Parameters
    ----------
    T : float
        Temperature in K.

    Returns
    -------
    float
        Water surface tension in N.m-1.
    """
    return 235.8e-3 * ((647.15 - T) / 647.15) ** 1.256 * (1 - 0.625 * (647.15 - T) / 647.15)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def K0(element, epsilon, epsilon_c=None):
    """This function calculates the intrinsic permeability, in m², considering GDL compression.

    Parameters
    ----------
    element : str
        Specifies the element for which the intrinsic permeability is calculated.
        Must be either 'gdl' (gas diffusion layer) or 'cl' (catalyst layer).
    epsilon : float
        Porosity.
    epsilon_c : float, optional
        Compression ratio of the GDL.

    Returns
    -------
    float
        Intrinsic permeability in m².

    Sources
    -------
    1. Qin Chen 2020 - Two-dimensional multi-physics modeling of porous transport layer in polymer electrolyte membrane
    electrolyzer for water splitting - for the Blake-Kozeny equation.
    2. M.L. Stewart 2005 - A study of pore geometry effects on anisotropy in hydraulic permeability using the
    lattice-Boltzmann method - for the Blake-Kozeny equation.
    """

    if element == 'gdl':
        # According to the GDL porosity, the GDL compression effect is different.
        if epsilon < 0.67:
            beta1 = -3.60
        else:
            beta1 = -2.60
        return epsilon / (8 * math.log(epsilon) ** 2) * (epsilon - epsilon_p) ** (alpha_p + 2) * \
            4.6e-6 ** 2 / ((1 - epsilon_p) ** alpha_p * ((alpha_p + 1) * epsilon - epsilon_p) ** 2) * math.exp(beta1 * epsilon_c)

    elif element == 'mpl':
        return (Dp_mpl**2 / 150) * (epsilon**3 / ((1-epsilon)**2)) # Using the Blake-Kozeny equation

    elif element == 'cl':
        return (Dp_cl**2 / 150) * (epsilon**3 / ((1-epsilon)**2)) # Using the Blake-Kozeny equation

    else:
        raise ValueError("The element should be either 'gdl', 'mpl' or 'cl'.")


def k_H2(lambdaa, T, kappa_co):
    """This function calculates the permeability coefficient of the membrane for hydrogen, in mol.m−1.s−1.Pa−1.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.
    T : float
        Temperature in K.
    kappa_co : float
        Crossover correction coefficient in mol.m-1.s-1.Pa-1.

    Returns
    -------
    float
        Permeability coefficient of the membrane for hydrogen in mol.m−1.s−1.Pa−1.
    """

    # Calculation of the permeability coefficient of the membrane for hydrogen
    k_H2_d = kappa_co * (0.29 + 2.2 * fv(lambdaa, T)) * 1e-14 * math.exp(Eact_H2_cros_v / R * (1 / Tref_cross - 1 / T))
    k_H2_l = kappa_co * 1.8 * 1e-14 * math.exp(Eact_H2_cros_l / R * (1 / Tref_cross - 1 / T))

    # Transition function between under-saturated and liquid-saturated states
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_l_eq(T) - lambdaa)))  # transition function

    return w * k_H2_d + (1 - w) * k_H2_l  # interpolation between under-saturated and liquid-equilibrated H2 crossover


def k_O2(lambdaa, T, kappa_co):
    """This function calculates the permeability coefficient of the membrane for oxygen, in mol.m−1.s−1.Pa−1.

    Parameters
    ----------
    lambdaa : float
        Water content in the membrane.
    T : float
        Temperature in K.
    kappa_co : float
        Crossover correction coefficient in mol.m-1.s-1.Pa-1.

    Returns
    -------
    float
        Permeability coefficient of the membrane for oxygen in mol.m−1.s−1.Pa−1.
    """

    # Calculation of the permeability coefficient of the membrane for oxygen
    k_O2_v = kappa_co * (0.11 + 1.9 * fv(lambdaa, T)) * 1e-14 * math.exp(Eact_O2_cros_v / R * (1 / Tref_cross - 1 / T))
    k_O2_l = kappa_co * 1.2 * 1e-14 * math.exp(Eact_O2_cros_l / R * (1 / Tref_cross - 1 / T))

    # Transition function between under-saturated and liquid-saturated states
    K_transition = 10  # It is a constant that defines the sharpness of the transition between two states. The higher it is, the sharper the transition is.
    w = 0.5 * (1 + math.tanh(K_transition * (lambda_l_eq(T) - lambdaa)))  # transition function

    return w * k_O2_v + (1 - w) * k_O2_l  # interpolation between under-saturated and liquid-equilibrated O2 crossover
