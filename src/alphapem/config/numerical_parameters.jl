# -*- coding: utf-8 -*-

"""This file contains the numerical parameter structure for the fuel cell simulation.
Default values are provided for a typical PEMFC, but they can be modified by the user
when creating an instance of the structure."""

# ============================================================
# NUMERICAL PARAMETERS
# ============================================================

# Pa. Development default for a numerical regularization, not a material constant.
const DEFAULT_PHASE_CHANGE_PRESSURE_WIDTH = 10.0

# mol.m-3. Provisional donor regularization; independent of atol and time step.
const DEFAULT_SORPTION_VAPOR_CONCENTRATION_WIDTH = 1.0e-6

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
    rtol::Float64 = 1e-3            # Relative tolerance for the system of ODEs solver
    atol::Float64 = 1e-6            # Absolute tolerance for the system of ODEs solver
    phase_change_pressure_width::Float64 = DEFAULT_PHASE_CHANGE_PRESSURE_WIDTH # Pa; 0 selects piecewise phase transfer.
    sorption_vapor_concentration_width::Float64 = DEFAULT_SORPTION_VAPOR_CONCENTRATION_WIDTH # mol.m-3; positive uptake ramp near an empty vapor donor.
    maxiters::Int = 5e5             # Maximum number of internal time steps that the IDA solver is allowed to take.
    max_run_time_s::Float64 = 5*60  # Maximum simulation runtime in seconds (0 = no limit)
    save_freq::Float64 = 10.0       # Maximum output save frequency in Hz (caps save density, does not force a sampling grid)
end

