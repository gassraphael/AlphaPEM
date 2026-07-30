# -*- coding: utf-8 -*-

"""This module is used to calculate intermediate values for the voltage calculation.
"""

# _________________________________________________Cell voltage modules_________________________________________________


"""Calculate the oxygen concentration at the platinum surface in the cathode catalyst layer.

Parameters
----------
i_fc : Float64
    The current density (A/m²).
sv : CellState1D
    The typed 1D cell-column state (MEA+GC) for one gas-channel position.
fc : AbstractFuelCell
    The fuel cell instance providing model parameters.

Returns
-------
Float64
    The oxygen concentration at the platinum surface in the cathode catalyst layer (mol/m³).

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function calculate_C_O2_Pt(i_fc::Real,
                           sv::CellState1D,
                           fc::AbstractFuelCell)

    # Extraction of the variables (typed access via CellState1D struct fields)
    C_O2_ccl, C_O2_cmpl = sv.ccl.C_O2, getproperty.(sv.cmpl, :C_O2)
    s_ccl, s_cmpl = sv.ccl.s, getproperty.(sv.cmpl, :s)
    T_ccl, T_cmpl = sv.ccl.T, getproperty.(sv.cmpl, :T)
    lambda_ccl = sv.ccl.lambda

    return calculate_C_O2_Pt(i_fc, s_ccl, lambda_ccl, C_O2_ccl, T_ccl, fc)
end


"""Calculate oxygen concentration at Pt surface from explicit cathode CL state variables.

This overload avoids constructing a full `CellState1D` when only local CCL
variables are needed (e.g. during initialization).
"""
function calculate_C_O2_Pt(i_fc::Real,
                           s_ccl::Real,
                           lambda_ccl::Real,
                           C_O2_ccl::Real,
                           T_ccl::Real,
                           fc::AbstractFuelCell)

    # Extraction of the parameters
    pp = fc.physical_parameters
    Hccl, K_O2_ad_Pt = pp.Hccl, pp.K_O2_ad_Pt
    C_O2_Pt_raw = C_O2_ccl - i_fc / (4 * F * Hccl) *
                  R_T_O2_Pt(s_ccl, lambda_ccl, T_ccl, Hccl, K_O2_ad_Pt, pp) /
                  a_c(:ccl, lambda_ccl, T_ccl, Hccl, pp)
    return _positive_concentration_value(C_O2_Pt_raw) # This avoids an unphysical negative C_O2_Pt and keeps Newton stable.
end


"""This function calculates the total resistance of oxygen to the platinium particules inside the CCL, defined as the
 sum of the different dissolution, diffusion and adsorption resistances.

Parameters
----------
s :
    Liquid water saturation in the CL.
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
K_O2_ad_Pt : Float64
    Interfacial resistance coefficient of O2 adsorption on the Pt sites, without units.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
Real
    Total resistance of O2 inside the CCL to the Pt particules in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
2. Georg A. Futter - Article 2018 - Physical modeling of polymer-electrolyte membrane fuel cells - Understanding
water management and impedance spectra.
3. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
"""
function R_T_O2_Pt(s, lambdaa, T, Hcl::Float64, K_O2_ad_Pt::Float64, pp::PhysicalParams)
    return R_O2_dis_l(s, lambdaa, T, Hcl, pp) + R_O2_dif_l(s, lambdaa, T, Hcl, pp) +
           R_O2_dis_ion(lambdaa, T, Hcl, pp)  + R_O2_dif_ion_eff(lambdaa, T, Hcl, pp) +
           R_O2_ad_Pt_eff(lambdaa, T, Hcl, K_O2_ad_Pt, pp)
end


"""This function calculates the dissolution resistance of oxygen in the CCL liquid water film, in s.m-1.
The assumption to make R_02_dis_l proportional to R_O2_dif_l is strong.

Parameters
----------
s :
    Liquid water saturation in the CL.
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
R_O2_dis_l
    Dissolution resistance of O2 in the liquid water film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_dis_l(s, lambdaa, T, Hcl::Float64, pp::PhysicalParams)
    K_O2_dis_l = pp.K_O2_dis_l  # Interfacial resistance coefficient of O2 dissolution inside the CL liquid water.
    return K_O2_dis_l * R_O2_dif_l(s, lambdaa, T, Hcl, pp)
end


"""This function calculates the diffusion resistance of oxygen inside the CCL liquid water film, in s.m-1.

Parameters
----------
s :
    Liquid water saturation in the CL.
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
R_O2_dif_l
    Diffusion resistance of O2 inside the CCL liquid water film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
2. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells
3. Ping Han - Article 1996 - Temperature dependence of oxygen diffusion in H20 and D20
"""
function R_O2_dif_l(s, lambdaa, T, Hcl::Float64, pp::PhysicalParams)
    r_carb = pp.r_carb  # Mean radius of the carbon particles.

    T_eff = _positive_temperature_value(T)
    delta_ion_val = delta_ion(:ccl, lambdaa, T_eff, Hcl, pp)
    s_num = _nonnegative_value(s)
    delta_H2O_l = (s_num * epsilon_cl(:ccl, lambdaa, T_eff, Hcl, pp) * r_carb^3 / epsilon_carb(:ccl, Hcl, pp) +
                  (r_carb + delta_ion_val)^3)^(1 / 3) -
                  (r_carb + delta_ion_val) # The liquid water film thickness in the CL, in m.

    D_O2_dif_l = 10^(-8.410 + 773.8 / T_eff - (506.4 / T_eff)^2) # The effective diffusion coefficient of O2 in the liquid water film, in m².s-1.

    return delta_H2O_l / D_O2_dif_l
end


"""This function calculates the dissolution resistance of oxygen in the CCL ionomer film, in s.m-1.
The assumption to make R_02_dis_ion proportional to R_02_dif_ion is strong.

Parameters
----------
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
R_O2_dis_ion
    Dissolution resistance of O2 in the CCL ionomer film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_dis_ion(lambdaa, T, Hcl::Float64, pp::PhysicalParams)
    K_O2_dis_ion = pp.K_O2_dis_ion  # Interfacial resistance coefficient of O2 dissolution inside the ionomer.
    return K_O2_dis_ion * R_O2_dif_ion(lambdaa, T, Hcl, pp)
end


"""This function calculates the diffusion resistance of oxygen inside the CCL ionomer film, in s.m-1.

Parameters
----------
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
R_O2_dif_ion :
    Diffusion resistance of O2 inside the CCL ionomer film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
2. Georg A. Futter - Article 2018 - Physical modeling of polymer-electrolyte membrane fuel cells - Understanding
water management and impedance spectra.
"""
function R_O2_dif_ion(lambdaa, T, Hcl::Float64, pp::PhysicalParams)

    T_eff = _positive_temperature_value(T)
    D_O2_dif_ion = 17.45e-10 * exp(-1514 / T_eff) # This is the effective diffusion coefficient of O2 in the ionomer film, in m².s-1.

    return delta_ion(:ccl, lambdaa, T_eff, Hcl, pp) / D_O2_dif_ion
end


"""This function calculates the effective diffusion resistance of oxygen inside the CCL ionomer film, in s.m-1.

Parameters
----------
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
R_O2_dif_ion_eff
    Effective diffusion resistance of O2 inside the CCL ionomer film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_dif_ion_eff(lambdaa, T, Hcl::Float64, pp::PhysicalParams)
    r_carb, theta_Pt_0, wt_Pt_ccl = pp.r_carb, pp.theta_Pt_0, pp.wt_Pt_ccl

    delta_ion_val = delta_ion(:ccl, lambdaa, T, Hcl, pp)
    r_Pt_val = r_Pt(pp)
    geom_factor = (r_carb + delta_ion_val)^2 / (r_Pt_val^2 * (1 - theta_Pt_0)) *
                  rho_Pt / rho_carb * (r_Pt_val / r_carb)^3 * (1 - wt_Pt_ccl) / wt_Pt_ccl

    return geom_factor * R_O2_dif_ion(lambdaa, T, Hcl, pp)
end


"""This function calculates the adsorption resistance of oxygen on the Pt particules inside the CCL, in s.m-1.
The assumption to make R_O2_ad_Pt proportional to R_O2_dif_ion is strong.

Parameters
----------
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
K_O2_ad_Pt : Float64
    Interfacial resistance coefficient of O2 adsorption on the Pt sites, without units.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
R_O2_ad_Pt :
    Adsorption resistance of O2 on the Pt particules inside the CCL, in s

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_ad_Pt(lambdaa, T, Hcl::Float64, K_O2_ad_Pt::Float64, pp::PhysicalParams)
    return K_O2_ad_Pt * R_O2_dif_ion(lambdaa, T, Hcl, pp)
end


"""This function calculates the effective adsorption resistance of oxygen on the Pt particules inside the CCL, in s.m-1.
Parameters
----------
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
K_O2_ad_Pt : Float64
    Interfacial resistance coefficient of O2 adsorption on the Pt sites, without units.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
R_O2_ad_Pt_eff
    Effective adsorption resistance of O2 on the Pt particules inside the CCL, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_ad_Pt_eff(lambdaa, T, Hcl::Float64, K_O2_ad_Pt::Float64, pp::PhysicalParams)
    r_carb, theta_Pt_0, wt_Pt_ccl = pp.r_carb, pp.theta_Pt_0, pp.wt_Pt_ccl
    delta_ion_val = delta_ion(:ccl, lambdaa, T, Hcl, pp)
    r_Pt_val = r_Pt(pp)
    geom_factor = (r_carb + delta_ion_val)^2 / (r_Pt_val^2 * (1 - theta_Pt_0)) *
                  rho_Pt / rho_carb * (r_Pt_val / r_carb)^3 * (1 - wt_Pt_ccl) / wt_Pt_ccl

    return geom_factor * R_O2_ad_Pt(lambdaa, T, Hcl, K_O2_ad_Pt, pp)
end


"""This function calculates the platine particle radius, in m.

Parameters
----------
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
r_Pt
    Platine particle radius in m.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function r_Pt(pp::PhysicalParams)
    ECSA_0, L_Pt_ccl = pp.ECSA_0, pp.L_Pt_ccl  # Initial electrochemical surface area and Pt loading of the catalyst.
    return 3 / (rho_Pt * ECSA_0 / L_Pt_ccl) # This is the platine particle radius, in m.
end


"""This function calculates the ionomer film thickness in the CL, in m. It should be in [7-9] nm.

Parameters
----------
element : Symbol
    Either `:acl` (anode) or `:ccl` (cathode) — selects the electrode-specific Pt loading and Pt/C weight fraction.
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
delta_ion
    Ionomer film thickness in the CL in m.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
2. Georg A. Futter - Article 2018 - Physical modeling of polymer-electrolyte membrane fuel cells - Understanding
water management and impedance spectra.
"""
function delta_ion(element::Symbol, lambdaa, T, Hcl::Float64, pp::PhysicalParams)
    r_carb = pp.r_carb  # Mean radius of the carbon particles.
    return r_carb * ((epsilon_mc(element, lambdaa, T, Hcl, pp) / epsilon_carb(element, Hcl, pp) + 1)^(1 / 3) - 1)
end


"""This function calculates the carbon volume fraction in the catalyst layer (ACL or CCL).

Parameters
----------
element : Symbol
    Either `:acl` (anode) or `:ccl` (cathode) — selects the electrode-specific Pt loading and Pt/C weight fraction.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
epsilon_carb :
    Carbon volume fraction in the CL.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function epsilon_carb(element::Symbol, Hcl::Float64, pp::PhysicalParams)
    L_Pt, wt_Pt = element == :acl ? (pp.L_Pt_acl, pp.wt_Pt_acl) :
                  element == :ccl ? (pp.L_Pt_ccl, pp.wt_Pt_ccl) :
                  throw(ArgumentError("The element should be either 'acl' or 'ccl'."))
    L_carb = L_Pt * (1 - wt_Pt) / wt_Pt  # This is the carbon loading in the CL, in kg.m-2.
    epsilon_carb_val = L_carb / (rho_carb * Hcl) # This is the volume fraction of carbon in the CL.
    if epsilon_carb_val >= 1
        println("epsilon_carb: ", epsilon_carb_val, " element: ", element, " Hcl: ", Hcl, " wt_Pt: ", wt_Pt)
        throw(ArgumentError("The calculated carbon volume fraction in the $(element) is greater than or equal to 1. Please check the inputs Hcl and wt_Pt."))
    end
    return epsilon_carb_val
end


"""This function calculates the Pt volume fraction in the catalyst layer (ACL or CCL).

Parameters
----------
element : Symbol
    Either `:acl` (anode) or `:ccl` (cathode) — selects the electrode-specific Pt loading.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
epsilon_Pt :
    Pt volume fraction in the CL.

Sources
-------
1. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
"""
function epsilon_Pt(element::Symbol, Hcl::Float64, pp::PhysicalParams)
    L_Pt, wt_Pt = element == :acl ? (pp.L_Pt_acl, pp.wt_Pt_acl) :
                  element == :ccl ? (pp.L_Pt_ccl, pp.wt_Pt_ccl) :
                  throw(ArgumentError("The element should be either 'acl' or 'ccl'."))
    epsilon_Pt_val = L_Pt / (rho_Pt * Hcl)  # This is the volume fraction of Pt in the CL.
    if epsilon_Pt_val >= 1
        println("epsilon_Pt: ", epsilon_Pt_val, " element: ", element, " Hcl: ", Hcl, " wt_Pt: ", wt_Pt)
        throw(ArgumentError("The calculated Pt volume fraction in the $(element) is greater than or equal to 1. Please check the inputs Hcl and wt_Pt."))
    end
    return epsilon_Pt_val
end


"""This function calculates the volumetric surface area of the ionomer in the CL, in m-1.
Parameters
----------
element : Symbol
    Either `:acl` (anode) or `:ccl` (cathode).
lambdaa :
    Water content in the CL.
T_cl :
    Temperature inside the CL in K.
Hccl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
a_c :
    Specific surface area of the ionomer in the CL in m⁻¹.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function a_c(element::Symbol, lambdaa, T_cl, Hccl::Float64, pp::PhysicalParams)
    r_carb = pp.r_carb  # Mean radius of the carbon particles.
    return 3 * epsilon_carb(element, Hccl, pp) / r_carb^3 * (r_carb + delta_ion(element, lambdaa, T_cl, Hccl, pp))^2
end


"""This function calculates the ionomer volume fraction in the CL.

Parameters
----------
element : Symbol
    Either `:acl` (anode) or `:ccl` (cathode) — selects the electrode-specific Pt loading and Pt/C weight fraction.
lambda_cl :
    Water content in the CL.
T_cl :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
epsilon_mc :
    Ionomer volume fraction in the CL.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function epsilon_mc(element::Symbol, lambda_cl, T_cl, Hcl::Float64, pp::PhysicalParams)

    IC = element == :acl ? pp.IC_acl :
         element == :ccl ? pp.IC_ccl :
         throw(ArgumentError("The element should be either 'acl' or 'ccl'."))
    rho_ion, M_eq = pp.rho_ion, pp.M_eq

    lambda_eff = _nonnegative_value(lambda_cl)
    epsilon_mc_val = IC * epsilon_carb(element, Hcl, pp) * rho_carb / rho_ion *
                     (1 + (M_H2O * rho_ion) / (rho_H2O_l(T_cl) * M_eq) * lambda_eff)

    if epsilon_mc_val >= 1
        println("epsilon_mc: ", epsilon_mc_val, " element: ", element, " Hcl: ", Hcl, " IC: ", IC)
        throw(ArgumentError("The calculated ionomer volume fraction in the $(element) is greater than or equal to 1. Please check the inputs Hcl and IC."))
    end
    return epsilon_mc_val
end


"""This function calculates the CL porosity.

Parameters
----------
element : Symbol
    Either `:acl` (anode) or `:ccl` (cathode) — selects the electrode-specific Pt loading and Pt/C weight fraction.
lambda_cl :
    Water content in the CL.
T_cl :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.
pp : PhysicalParams
    Physical parameters of the fuel cell.

Returns
-------
epsilon_cl
    CL porosity.

Sources
-------
1. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
"""
function epsilon_cl(element::Symbol, lambda_cl, T_cl, Hcl::Float64, pp::PhysicalParams)

    epsilon_cl_val = 1 - epsilon_carb(element, Hcl, pp) - epsilon_Pt(element, Hcl, pp) -
                      epsilon_mc(element, lambda_cl, T_cl, Hcl, pp)

    if epsilon_cl_val <= 0
        println("epsilon_cl: ", epsilon_cl_val, " element: ", element, " Hcl: ", Hcl)
        throw(ArgumentError("The calculated porosity in the $(element) is less than or equal to 0. Please check the inputs Hcl and, for the $(element), wt_Pt."))
    end
    return epsilon_cl_val
end
