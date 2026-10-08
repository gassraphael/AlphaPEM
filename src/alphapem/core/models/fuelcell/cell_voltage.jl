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
sv : CellState1D
    The typed 1D cell-column state (MEA+GC) for one gas-channel position.
fc : AbstractFuelCell
    The fuel cell instance providing model parameters.

Returns
-------
Float64
    The cell voltage in volt.
"""
function calculate_cell_voltage(i_fc::Real, C_O2_Pt::Real, sv::CellState1D, fc::AbstractFuelCell)

    # Extraction of the variables
    T_acl = _positive_temperature_value(sv.acl.T)
    T_mem = _positive_temperature_value(sv.mem.T)
    T_ccl = _positive_temperature_value(sv.ccl.T)

    lambda_mem, lambda_ccl = sv.mem.lambda, sv.ccl.lambda

    C_H2_acl = _nonnegative_value(sv.acl.C_H2)

    eta_c = sv.ccl.eta_c
    C_O2_Pt_safe = _positive_concentration_value(C_O2_Pt)

    # Extraction of the parameters
    pp = fc.physical_parameters
    Hmem, Hccl = pp.Hmem, pp.Hccl
    Re = pp.Re

    # The equilibrium potential
    Ueq = E0 - 8.5e-4 * (T_ccl - 298.15) + R * T_ccl / (2 * F) *
          (log(R * T_acl * C_H2_acl / Pref_eq) +
           0.5 * log(R * T_ccl * C_O2_Pt_safe / Pref_eq))

    # The proton resistance
    #       The proton resistance at the membrane : Rmem
    Rmem = Hmem / sigma_p_eff(:mem, lambda_mem, T_mem, nothing, pp)
    #       The proton resistance at the cathode catalyst layer : Rccl
    Rccl = Hccl / sigma_p_eff(:ccl, lambda_ccl, T_ccl, Hccl, pp)
    #       The total proton resistance
    Rp = Rmem + Rccl  # Its value is around [4-7]e-6 ohm.m².

    # Ohmic voltage losses use the external current density only.
    # Gass et al., IJHE (2025), Eq. (49), doi:10.1016/j.ijhydene.2024.11.374.
    # The critical review, Eq. (64) and Sec. 8.3.3, arXiv:2410.13323v1
    # (doi:10.1149/1945-7111/ad305a), explicitly follows O'Hayre and argues
    # against adding the crossover-equivalent current to the ohmic term
    # as in the combined voltage formulation of Dicks and Rand.
    # Primary sources:
    # - O'Hayre, Cha, Colella and Prinz, Fuel Cell Fundamentals, 3rd ed.,
    #   Wiley, 2016, Sec. 6.1, p. 206, Eqs. (6.3)-(6.4):
    #   activation/concentration losses use j + j_leak; ohmic losses use j.
    # - Dicks and Rand, Fuel Cell Systems Explained, 3rd ed., Wiley, 2018,
    #   Sec. 3.8, p. 57, Eq. (3.22): the combined equation uses (i + i_n) * r.
    # Gas crossover remains in the overpotential dynamics and species balances.
    Ucell = Ueq - eta_c - i_fc * (Rp + Re)

    return Ucell
end
