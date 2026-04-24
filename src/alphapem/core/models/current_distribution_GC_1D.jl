# -*- coding: utf-8 -*-

"""This file represents the equations for calculating the current distribution through the GC.
It is a component of the fuel cell model.
"""

# _____________________________________________________Preliminaries____________________________________________________

# Importing the necessary libraries
using NonlinearSolve
using SciMLBase: successful_retcode


# _________________________________________________Current distribution_________________________________________________

"""Calculate the local current density distribution in the 1D direction of the GC.

Parameters
----------
i_fc_cell : Float64
    Fuel cell current density at time t (A.m-2).
sv : AbstractVector{<:CellState1D}
    Typed internal states calculated by the solver.
    sv is a contraction of solver_variables for enhanced readability.
    sv[i] is the typed 1D cell-column state (MEA+GC) associated with gas channel i.
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.

Returns
-------
Vector{Float64}
    Local current density distribution in the 1D direction of the GC (A.m-2).
    Julia vectors are naturally 1-based, so no dummy element is stored at index 0.
    i_fc[i] corresponds to gas channel i, for i in 1:nb_gc.
"""
function calculate_1D_GC_current_density(i_fc_cell::Float64, sv::AbstractVector{<:CellState1D}, fc::AbstractFuelCell)::Vector{Float64}

    # Extraction of the parameters
    nb_gc = fc.numerical_parameters.nb_gc
        # Fast path: if there is only one gas channel, the local current density equals the cell current density.
    nb_gc == 1 && return [Float64(i_fc_cell)]
    # Extraction of the variables
    C_O2_ccl = [sv[i].ccl.C_O2 for i in 1:nb_gc]

    # Internal scaling improves conditioning of this nonlinear solve (mixed
    # voltage/current/concentration magnitudes) while keeping a physical API.
    x_scales, res_scales = _build_gc_current_density_scaling(nb_gc)

    # Residual function for NonlinearSolve solver applied on the local current density
    function residuals!(res, x, _)

        # Convert solver unknowns from scaled to physical values.
        x_phys = x .* x_scales

        # Recovery of the guessed variable values
        U_cell_guessed  = x_phys[1]
        @views i_fc_guessed    = x_phys[2:nb_gc+1]        # view: no copy allocated
        @views C_O2_Pt_guessed = x_phys[nb_gc+2:2*nb_gc+1]  # view: no copy allocated

        # Residuals: difference between guessed and calculated algebraic states.
        gc_current_distribution_residuals!(res, U_cell_guessed, i_fc_guessed, C_O2_Pt_guessed, i_fc_cell, sv, fc)

        # Return dimensionless residuals to balance equation blocks.
        res ./= res_scales
    end

    # Calculation of the 1D GC current density by solving the system of equations defined by the residuals function
    # using NonlinearSolve with Levenberg-Marquardt
    #       Initial guesses
    x0 = Vector{Float64}(undef, 2 * nb_gc + 1)
    x0[1]            = calculate_cell_voltage(i_fc_cell, C_O2_ccl[1], sv[1], fc)
    x0[2:nb_gc+1]   .= i_fc_cell
    x0[nb_gc+2:2*nb_gc+1] = C_O2_ccl
    x0_scaled = scale_values(x0, x_scales)
    #       Solver call
    prob = NonlinearProblem(residuals!, x0_scaled, nothing)
    sol  = solve(prob, LevenbergMarquardt(); abstol=1e-6, reltol=1e-6, maxiters=400)

    #       Check for convergence
    if !successful_retcode(sol.retcode)
        error("Convergence failed in calculate_1D_GC_current_density: retcode = $(sol.retcode)")
    end
    #       Extract the results in physical units
    sol_phys = unscale_values(sol.u, x_scales)
    i_fc = Vector(sol_phys[2:nb_gc+1])

    return i_fc
end

