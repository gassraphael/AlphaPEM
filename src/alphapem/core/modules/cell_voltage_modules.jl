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
function calculate_C_O2_Pt(i_fc::Float64,
                           sv::CellState1D,
                           fc::AbstractFuelCell)::Float64

    # Extraction of the variables (typed access via CellState1D struct fields)
    s_ccl      = sv.ccl.s
    lambda_ccl = sv.ccl.lambda
    C_O2_ccl   = sv.ccl.C_O2
    T_ccl      = sv.ccl.T
    # Extraction of the parameters
    Hccl, K_O2_ad_Pt = fc.physical_parameters.Hccl, fc.physical_parameters.K_O2_ad_Pt

    C_O2_Pt = C_O2_ccl - i_fc / (4 * F * Hccl) *
              R_T_O2_Pt(s_ccl, lambda_ccl, T_ccl, Hccl, K_O2_ad_Pt) /
              a_c(lambda_ccl, T_ccl, Hccl)

    return max(C_O2_Pt, 0.0)
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
function R_T_O2_Pt(s, lambdaa, T, Hcl::Float64, K_O2_ad_Pt::Float64)
    return R_O2_dis_l(s, lambdaa, T, Hcl) + R_O2_dif_l(s, lambdaa, T, Hcl) +
           R_O2_dis_ion(lambdaa, T, Hcl)  + R_O2_dif_ion_eff(lambdaa, T, Hcl) +
           R_O2_ad_Pt_eff(lambdaa, T, Hcl, K_O2_ad_Pt)
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

Returns
-------
R_O2_dis_l
    Dissolution resistance of O2 in the liquid water film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_dis_l(s, lambdaa, T, Hcl::Float64)
    return K_O2_dis_l * R_O2_dif_l(s, lambdaa, T, Hcl)
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
function R_O2_dif_l(s, lambdaa, T, Hcl::Float64)

    delta_ion_val = delta_ion(lambdaa, T, Hcl)
    delta_H2O_l = (s * epsilon_cl(lambdaa, T, Hcl) * r_carb^3 / epsilon_carb(Hcl) +
                  (r_carb + delta_ion_val)^3)^(1 / 3) -
                  (r_carb + delta_ion_val) # The liquid water film thickness in the CL, in m.

    D_O2_dif_l = 10^(-8.410 + 773.8 / T - (506.4 / T)^2) # The effective diffusion coefficient of O2 in the liquid water film, in m².s-1.

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

Returns
-------
R_O2_dis_ion
    Dissolution resistance of O2 in the CCL ionomer film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_dis_ion(lambdaa, T, Hcl::Float64)
    return K_O2_dis_ion * R_O2_dif_ion(lambdaa, T, Hcl)
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
function R_O2_dif_ion(lambdaa, T, Hcl::Float64)

    D_O2_dif_ion = 17.45e-10 * exp(-1514 / T) # This is the effective diffusion coefficient of O2 in the ionomer film, in m².s-1.

    return delta_ion(lambdaa, T, Hcl) / D_O2_dif_ion
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

Returns
-------
R_O2_dif_ion_eff
    Effective diffusion resistance of O2 inside the CCL ionomer film, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_dif_ion_eff(lambdaa, T, Hcl::Float64)
    delta_ion_val = delta_ion(lambdaa, T, Hcl)
    r_Pt_val = r_Pt()
    geom_factor = (r_carb + delta_ion_val)^2 / (r_Pt_val^2 * (1 - theta_Pt_0)) *
                  rho_Pt / rho_carb * (r_Pt_val / r_carb)^3 * (1 - wt_Pt) / wt_Pt

    return geom_factor * R_O2_dif_ion(lambdaa, T, Hcl)
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

Returns
-------
R_O2_ad_Pt :
    Adsorption resistance of O2 on the Pt particules inside the CCL, in s

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_ad_Pt(lambdaa, T, Hcl::Float64, K_O2_ad_Pt::Float64)
    return K_O2_ad_Pt * R_O2_dif_ion(lambdaa, T, Hcl)
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

Returns
-------
R_O2_ad_Pt_eff
    Effective adsorption resistance of O2 on the Pt particules inside the CCL, in s.m-1.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function R_O2_ad_Pt_eff(lambdaa, T, Hcl::Float64, K_O2_ad_Pt::Float64)
    delta_ion_val = delta_ion(lambdaa, T, Hcl)
    r_Pt_val = r_Pt()
    geom_factor = (r_carb + delta_ion_val)^2 / (r_Pt_val^2 * (1 - theta_Pt_0)) *
                  rho_Pt / rho_carb * (r_Pt_val / r_carb)^3 * (1 - wt_Pt) / wt_Pt

    return geom_factor * R_O2_ad_Pt(lambdaa, T, Hcl, K_O2_ad_Pt)
end


"""This function calculates the platine particle radius, in m.

Returns
-------
r_Pt
    Platine particle radius in m.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function r_Pt()
    return 3 / (rho_Pt * ECSA_0 / L_Pt) # This is the platine particle radius, in m.
end


"""This function calculates the ionomer film thickness in the CL, in m. It should be in [7-9] nm.

Parameters
----------
lambdaa :
    Water content in the CL.
T :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.

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
function delta_ion(lambdaa, T, Hcl::Float64)
    return r_carb * ((epsilon_mc(lambdaa, T, Hcl) / epsilon_carb(Hcl) + 1)^(1 / 3) - 1)
end


"""This function calculates the carbon volume fraction in the CCL.

Parameters
----------
Hccl : Float64
    Thickness of the CCL layer.

Returns
-------
epsilon_carb :
    Carbon volume fraction in the CCL.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function epsilon_carb(Hccl::Float64)
    L_carb = L_Pt * (1 - wt_Pt) / wt_Pt  # This is the carbon loading in the CCL, in kg.m-2.
    epsilon_carb_val = L_carb / (rho_carb * Hccl) # This is the volume fraction of carbon in the CCL.
    if epsilon_carb_val >= 1
        println("epsilon_carb: ", epsilon_carb_val, " Hccl: ", Hccl, " wt_Pt: ", wt_Pt)
        throw(ArgumentError("The calculated carbon volume fraction in the CCL is greater than or equal to 1. Please check the inputs Hccl and wt_Pt."))
    end
    return epsilon_carb_val
end


"""This function calculates the Pt volume fraction in the CCL.

Parameters
----------
Hccl : Float64
    Thickness of the CCL layer.

Returns
-------
epsilon_Pt :
    Carbon volume fraction in the CCL.

Sources
-------
1. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
"""
function epsilon_Pt(Hccl::Float64)
    epsilon_Pt_val = L_Pt / (rho_Pt * Hccl)  # This is the volume fraction of Pt in the CCL.
    if epsilon_Pt_val >= 1
        println("epsilon_Pt: ", epsilon_Pt_val, " Hccl: ", Hccl, " wt_Pt: ", wt_Pt)
        throw(ArgumentError("The calculated Pt volume fraction in the CCL is greater than or equal to 1. Please check the inputs Hccl and wt_Pt."))
    end
    return epsilon_Pt_val
end


"""This function calculates the volumetric surface area of the ionomer in the CL, in m-1.
Parameters
----------
lambdaa :
    Water content in the CL.
T_cl :
    Temperature inside the CL in K.
Hccl : Float64
    Thickness of the CL layer.

Returns
-------
a_c :
    Specific surface area of the ionomer in the CL in m⁻¹.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function a_c(lambdaa, T_cl, Hccl::Float64)
    return 3 * epsilon_carb(Hccl) / r_carb^3 * (r_carb + delta_ion(lambdaa, T_cl, Hccl))^2
end


"""This function calculates the ionomer volume fraction in the CL.

Parameters
----------
lambda_cl :
    Water content in the CL.
T_cl :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.

Returns
-------
epsilon_mc :
    Ionomer volume fraction in the CL.

Sources
-------
1. Liang Hao - Article 2015 - Modeling and Experimental Validation of Pt Loading and Electrode Composition Effects
in PEM Fuel Cells.
"""
function epsilon_mc(lambda_cl, T_cl, Hcl::Float64)

    epsilon_mc_val = IC * epsilon_carb(Hcl) * rho_carb / rho_ion *
                     (1 + (M_H2O * rho_ion) / (rho_H2O_l(T_cl) * M_eq) * lambda_cl)

    if epsilon_mc_val >= 1
        println("epsilon_mc: ", epsilon_mc_val, " Hcl: ", Hcl, " IC: ", IC, " wt_Pt: ", wt_Pt)
        throw(ArgumentError("The calculated ionomer volume fraction in the CCL is greater than or equal to 1. Please check the inputs Hcl, IC and wt_Pt."))
    end
    return epsilon_mc_val
end


"""This function calculates the CL porosity.

Parameters
----------
lambda_cl :
    Water content in the CL.
T_cl :
    Temperature inside the CL in K.
Hcl : Float64
    Thickness of the CL layer.

Returns
-------
epsilon_cl
    CL porosity.

Sources
-------
1. Alireza Goshtasbi - Article 2020 - A Mathematical Model toward Real-Time Monitoring of Automotive PEM Fuel Cells.
"""
function epsilon_cl(lambda_cl, T_cl, Hcl::Float64)

    epsilon_cl_val = 1 - epsilon_carb(Hcl) - epsilon_Pt(Hcl) - epsilon_mc(lambda_cl, T_cl, Hcl)

    if epsilon_cl_val <= 0
        println("epsilon_cl: ", epsilon_cl_val, " Hcl: ", Hcl, " wt_Pt: ", wt_Pt)
        throw(ArgumentError("The calculated porosity in the CCL is less than or equal to 0. Please check the inputs Hcl and wt_Pt."))
    end
    return epsilon_cl_val
end

