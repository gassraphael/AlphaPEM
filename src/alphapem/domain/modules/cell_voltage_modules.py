# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the voltage calculation.
"""
import math
from functools import lru_cache

# _____________________________________________________Preliminaries____________________________________________________

# Importing constants' value and functions
from alphapem.utils.physics_constants import (F, K_O2_dis_l, r_carb, K_O2_dis_ion, theta_Pt_0, rho_Pt, rho_carb, wt_Pt,
                                              ECSA_0, L_Pt, IC, rho_ion, M_H2O, M_eq)
from alphapem.utils.physics_functions import rho_H2O_l


# _________________________________________________Cell voltage modules_________________________________________________

def calculate_C_O2_Pt(i_fc, sv_1D, parameters):
    """This function calculates the oxygen concentration at the platinum surface in the cathode catalyst layer.
    Parameters
    ----------
    i_fc : float
        The current density (A/m²).
    sv_1D : dict
        The dictionary containing the variables calculated by the solver.
    parameters : dict
        The dictionary containing the parameters.

    Returns
    -------
    C_O2_Pt : float
        The oxygen concentration at the platinum surface in the cathode catalyst layer (mol/m³).

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    # Extraction of the variables
    s_ccl, lambda_ccl = sv_1D['s_ccl'], sv_1D['lambda_ccl']
    C_O2_ccl = sv_1D['C_O2_ccl']
    T_ccl = sv_1D['T_ccl']
    # Extraction of the operating inputs and the parameters
    Hccl, K_O2_ad_Pt = parameters['Hccl'], parameters['K_O2_ad_Pt']

    C_O2_Pt = C_O2_ccl - i_fc / (4 * F * Hccl) * R_T_O2_Pt(s_ccl, lambda_ccl, T_ccl, Hccl, K_O2_ad_Pt) / a_c(lambda_ccl, T_ccl, Hccl)

    if C_O2_Pt <= 0:
        raise ValueError("Calculated C_O2_Pt is non-physical (negative or zero). Check input parameters and variables.")

    return C_O2_Pt


@lru_cache(maxsize=None) # Cache the results to optimize performance
def R_T_O2_Pt(s, lambdaa, T, Hcl, K_O2_ad_Pt):
    """This function calculates the total resistance of oxygen to the platinium particules inside the CCL, defined as the
     sum of the different dissolution, diffusion and adsorption resistances.

    Parameters
    ----------
    s : float
        Liquid water saturation in the CL.
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.
    K_O2_ad_Pt : float
        Interfacial resistance coefficient of O2 adsorption on the Pt sites, without units.

    Returns
    -------
    float
        Total resistance of O2 inside the CCL to the Pt particules in s.m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    2. Georg A. Futter - Article 2018 - Physical modeling of polymer-electrolyte membrane fuel cells - Understanding
    water management and impedance spectra.
    3. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
    """

    return R_O2_dis_l(s, lambdaa, T, Hcl) + R_O2_dif_l(s, lambdaa, T, Hcl) + \
           R_O2_dis_ion(lambdaa, T, Hcl) + R_O2_dif_ion_eff(lambdaa, T, Hcl) + \
           R_O2_ad_Pt_eff(lambdaa, T, Hcl, K_O2_ad_Pt)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def R_O2_dis_l(s, lambdaa, T, Hcl):
    """This function calculates the dissolution resistance of oxygen in the CCL liquid water film, in s.m-1.
    The assumption to make R_02_dis_l proportional to R_O2_dif_l is strong.

    Parameters
    ----------
    s : float
        Liquid water saturation in the CL.
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Dissolution resistance of O2 in the liquid water film, in s.m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return K_O2_dis_l * R_O2_dif_l(s, lambdaa, T, Hcl)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def R_O2_dif_l(s, lambdaa, T, Hcl):
    """This function calculates the diffusion resistance of oxygen inside the CCL liquid water film, in s.m-1.

    Parameters
    ----------
    s : float
        Liquid water saturation in the CL.
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Diffusion resistance of O2 inside the CCL liquid water film, in s.m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    2. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells
    3. Ping Han - Article 1996 - Temperature dependence of oxygen diffusion in H20 and D20
    """

    delta_ion_val = delta_ion(lambdaa, T, Hcl)
    delta_H2O_l = (s * epsilon_cl(lambdaa, T, Hcl) * r_carb**3 / epsilon_carb(Hcl) + (r_carb + delta_ion_val)**3)**(1/3) - \
                  (r_carb + delta_ion_val) # The liquid water film thickness in the CL, in m.

    D_O2_dif_l =  10 ** (-8.410 + 773.8 / T - (506.4 / T)**2) # The effective diffusion coefficient of O2 in the liquid water film, in m².s-1.

    return delta_H2O_l / D_O2_dif_l


def R_O2_dis_ion(lambdaa, T, Hcl):
    """This function calculates the dissolution resistance of oxygen in the CCL ionomer film, in s.m-1.
    The assumption to make R_02_dis_ion proportional to R_02_dif_ion is strong.

    Parameters
    ----------
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Dissolution resistance of O2 in the CCL ionomer film, in s.m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return K_O2_dis_ion * R_O2_dif_ion(lambdaa, T, Hcl)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def R_O2_dif_ion(lambdaa, T, Hcl):
    """This function calculates the diffusion resistance of oxygen inside the CCL ionomer film, in s.m-1.

    Parameters
    ----------
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Diffusion resistance of O2 inside the CCL ionomer film, in s.m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    2. Georg A. Futter - Article 2018 - Physical modeling of polymer-electrolyte membrane fuel cells - Understanding
    water management and impedance spectra.
    """

    D_O2_dif_ion = 17.45e-10 * math.exp(-1514 / T) # This is the effective diffusion coefficient of O2 in the ionomer film, in m².s-1.

    return delta_ion(lambdaa, T, Hcl) / D_O2_dif_ion


@lru_cache(maxsize=None) # Cache the results to optimize performance
def R_O2_dif_ion_eff(lambdaa, T, Hcl):
    """This function calculates the effective diffusion resistance of oxygen inside the CCL ionomer film, in s.m-1.

    Parameters
    ----------
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Effective diffusion resistance of O2 inside the CCL ionomer film, in s.m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return (r_carb + delta_ion(lambdaa, T, Hcl))**2 / (r_Pt()**2 * (1 - theta_Pt_0)) * \
           rho_Pt / rho_carb * (r_Pt() / r_carb)**3 * (1 - wt_Pt) / wt_Pt * \
           R_O2_dif_ion(lambdaa, T, Hcl)


def R_O2_ad_Pt(lambdaa, T, Hcl, K_O2_ad_Pt):
    """This function calculates the adsorption resistance of oxygen on the Pt particules inside the CCL, in s.m-1.
    The assumption to make R_O2_ad_Pt proportional to R_O2_dif_ion is strong.

    Parameters
    ----------
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.
    K_O2_ad_Pt : float
        Interfacial resistance coefficient of O2 adsorption on the Pt sites, without units.

    Returns
    -------
    float
        Adsorption resistance of O2 on the Pt particules inside the CCL, in s

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return K_O2_ad_Pt * R_O2_dif_ion(lambdaa, T, Hcl)


def R_O2_ad_Pt_eff(lambdaa, T, Hcl, K_O2_ad_Pt):
    """This function calculates the effective adsorption resistance of oxygen on the Pt particules inside the CCL, in s.m-1.
    Parameters
    ----------
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.
    K_O2_ad_Pt : float
        Interfacial resistance coefficient of O2 adsorption on the Pt sites, without units.

    Returns
    -------
    float
        Effective adsorption resistance of O2 on the Pt particules inside the CCL, in s.m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return (r_carb + delta_ion(lambdaa, T, Hcl))**2 / (r_Pt()**2 * (1 - theta_Pt_0)) * \
           rho_Pt / rho_carb * (r_Pt() / r_carb)**3 * (1 - wt_Pt) / wt_Pt * \
           R_O2_ad_Pt(lambdaa, T, Hcl, K_O2_ad_Pt)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def r_Pt():
    """This function calculates the platine particle radius, in m.

    Returns
    -------
    float
        Platine particle radius in m.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return 3 / (rho_Pt * ECSA_0 / L_Pt) # This is the platine particle radius, in m.


def delta_ion(lambdaa, T, Hcl):
    """This function calculates the ionomer film thickness in the CL, in m. It should be in [7-9] nm.

    Parameters
    ----------
    lambdaa : float
        Water content in the CL.
    T : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Ionomer film thickness in the CL in m.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    2. Georg A. Futter - Article 2018 - Physical modeling of polymer-electrolyte membrane fuel cells - Understanding
    water management and impedance spectra.
    """

    return r_carb * ((epsilon_mc(lambdaa, T, Hcl) / epsilon_carb(Hcl) + 1)**(1/3) - 1)


@lru_cache(maxsize=None) # Cache the results to optimize performance
def epsilon_carb(Hccl):
    """This function calculates the carbon volume fraction in the CCL.

    Parameters
    ----------
    Hccl : float
        Thickness of the CCL layer.

    Returns
    -------
    float
        Carbon volume fraction in the CCL.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    L_carb = L_Pt * (1 - wt_Pt) / wt_Pt  # This is the carbon loading in the CCL, in kg.m-2.
    epsilon_carb = L_carb / (rho_carb * Hccl) # This is the volume fraction of carbon in the CCL.

    if epsilon_carb >= 1:
        print("epsilon_carb: ", epsilon_carb, " Hccl: ", Hccl, " wt_Pt: ", wt_Pt)
        raise ValueError("The calculated carbon volume fraction in the CCL is greater than or equal to 1. "
                         "Please check the inputs Hccl and wt_Pt.")
    return epsilon_carb


@lru_cache(maxsize=None) # Cache the results to optimize performance
def epsilon_Pt(Hccl):
    """This function calculates the Pt volume fraction in the CCL.

    Parameters
    ----------
    Hccl : float
        Thickness of the CCL layer.

    Returns
    -------
    float
        Carbon volume fraction in the CCL.

    Sources
    -------
    1. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
    """

    epsilon_Pt = L_Pt / (rho_Pt * Hccl)  # This is the volume fraction of carbon in the CCL.

    if epsilon_Pt >= 1:
        print("epsilon_Pt: ", epsilon_Pt, " Hccl: ", Hccl, " wt_Pt: ", wt_Pt)
        raise ValueError("The calculated Pt volume fraction in the CCL is greater than or equal to 1. "
                         "Please check the inputs Hccl and wt_Pt.")
    return epsilon_Pt


@lru_cache(maxsize=None) # Cache the results to optimize performance
def a_c(lambdaa, T_cl, Hccl):
    """This function calculates the volumetric surface area of the ionomer in the CL, in m-1.
    Parameters
    ----------
    lambdaa : float
        Water content in the CL.
    T_cl : float
        Temperature inside the CL in K.
    Hccl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Specific surface area of the ionomer in the CL in m-1.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    return 3 * epsilon_carb(Hccl) / r_carb**3 * (r_carb + delta_ion(lambdaa, T_cl, Hccl))**2


def epsilon_mc(lambda_cl, T_cl, Hcl):
    """This function calculates the ionomer volume fraction in the CL.

    Parameters
    ----------
    lambda_cl : float
        Water content in the CL.
    T_cl : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        Ionomer volume fraction in the CL.

    Sources
    -------
    1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
    in PEM Fuel Cells.
    """

    epsilon_mc = IC * epsilon_carb(Hcl) * rho_carb / rho_ion * (1 + (M_H2O * rho_ion) / (
                rho_H2O_l(T_cl) * M_eq) * lambda_cl)

    if epsilon_mc >= 1:
        print("epsilon_mc: ", epsilon_mc, " Hcl: ", Hcl, " IC: ", IC, " wt_Pt: ", wt_Pt)
        raise ValueError("The calculated ionomer volume fraction in the CCL is greater than or equal to 1. "
                         "Please check the inputs Hcl, IC and wt_Pt.")
    return epsilon_mc


def epsilon_cl(lambda_cl, T_cl, Hcl):
    """This function calculates the CL porosity.

    Parameters
    ----------
    lambda_cl : float
        Water content in the CL.
    T_cl : float
        Temperature inside the CL in K.
    Hcl : float
        Thickness of the CL layer.

    Returns
    -------
    float
        CL porosity.

    Sources
    -------
    1. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
    """

    epsilon_cl = 1 - epsilon_carb(Hcl) - epsilon_Pt(Hcl) - epsilon_mc(lambda_cl, T_cl, Hcl)

    if epsilon_cl <= 0:
        print("epsilon_cl: ", epsilon_cl, " Hcl: ", Hcl, " wt_Pt: ", wt_Pt)
        raise ValueError("The calculated porosity in the CCL is less than or equal to 0. "
                         "Please check the inputs Hcl and wt_Pt.")
    return epsilon_cl
