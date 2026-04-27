# -*- coding: utf-8 -*-

"""This file contains the numerical parameter structure for the fuel cell simulation.
Default values are provided for a typical PEMFC, but they can be modified by the user
when creating an instance of the structure."""

# ============================================================
# NUMERICAL PARAMETERS
# ============================================================

"""
    NumericalParams

This structure contains the numerical parameters for the fuel cell simulation.
"""
Base.@kwdef struct NumericalParams <: AbstractFuelCellParams
    nb_gc::Int = 5                  # Number of model nodes placed inside each gas channel
    nb_gdl::Int = 3                 # Number of model nodes placed inside each GDL
    nb_mpl::Int = 2                 # Number of model nodes placed inside each MPL
    nb_man::Int = 1                 # Number of model nodes placed inside each manifold line
    purge_time::Float64 = 0.6       # The time for purging the system in seconds
    delta_purge::Float64 = 15.0     # The time between two purges in seconds
    delta_t_dyn_step::Float64 = 0.1 # Time for dynamic display of the step current density function in seconds
    rtol::Float64 = 1e-4            # Relative tolerance for the system of ODEs solver
    atol::Float64 = 5e-7            # Absolute tolerance for the system of ODEs solver
end

