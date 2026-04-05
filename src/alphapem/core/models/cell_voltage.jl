# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the cell voltage. It is a component of the fuel cell model.
"""

# _____________________________________________________Cell voltage_____________________________________________________

"""Calculate the cell voltage in volt.

Parameters
----------
i_fc : Float64
    The current density (A/m²).
C_O2_Pt : Float64
    The oxygen concentration at the platinum surface in the cathode catalyst layer (mol/m³).
sv : Dict
    The dictionary containing the variables calculated by the solver.
fc : AbstractFuelCell
    The fuel cell instance providing model parameters.

Returns
-------
Float64
    The cell voltage in volt.
"""
function calculate_cell_voltage(i_fc::Float64, C_O2_Pt::Float64, sv::Dict, fc::AbstractFuelCell)::Float64

    # Extraction of the variables
    s_ccl, lambda_mem, lambda_ccl = sv["s_ccl"], sv["lambda_mem"], sv["lambda_ccl"]
    C_H2_acl, C_O2_ccl = sv["C_H2_acl"], sv["C_O2_ccl"]
    eta_c = sv["eta_c"]
    T_acl, T_mem, T_ccl = sv["T_acl"], sv["T_mem"], sv["T_ccl"]
    # Extraction of the parameters
    Hmem, Hacl, Hccl = fc.physical_parameters.Hmem, fc.physical_parameters.Hacl, fc.physical_parameters.Hccl
    Re, kappa_co = fc.physical_parameters.Re, fc.physical_parameters.kappa_co

    # The equilibrium potential
    Ueq = E0 - 8.5e-4 * (T_ccl - 298.15) + R * T_ccl / (2 * F) * (log(R * T_acl * C_H2_acl / Pref_eq) +
                                                                  0.5 * log(R * T_ccl * C_O2_Pt / Pref_eq))

    # The crossover current density
    T_acl_mem_ccl = average([T_acl, T_mem, T_ccl],
                            [Hacl / (Hacl + Hmem + Hccl), Hmem / (Hacl + Hmem + Hccl), Hccl / (Hacl + Hmem + Hccl)])
    i_H2 = 2 * F * R * T_acl_mem_ccl / Hmem * C_H2_acl * k_H2(lambda_mem, T_mem, kappa_co)
    i_O2 = 4 * F * R * T_acl_mem_ccl / Hmem * C_O2_ccl * k_O2(lambda_mem, T_mem, kappa_co)
    i_n = i_H2 + i_O2

    # The proton resistance
    #       The proton resistance at the membrane : Rmem
    Rmem = Hmem / sigma_p_eff("mem", lambda_mem, T_mem)
    #       The proton resistance at the cathode catalyst layer : Rccl
    Rccl = Hccl / sigma_p_eff("ccl", lambda_ccl, T_ccl, Hccl)
    #       The total proton resistance
    Rp = Rmem + Rccl  # Its value is around [4-7]e-6 ohm.m².

    # The cell voltage
    Ucell = Ueq - eta_c - (i_fc + i_n) * (Rp + Re)

    return Ucell
end
