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
sv : AbstractVector{<:MEAState1D}
    Typed internal states calculated by the solver.
    sv is a contraction of solver_variables for enhanced readability.
    sv[i] is the typed 1D MEA state associated with gas channel i.
fc : AbstractFuelCell
    Fuel cell instance providing model parameters.

Returns
-------
Vector{Float64}
    Local current density distribution in the 1D direction of the GC (A.m-2).
    Julia vectors are naturally 1-based, so no dummy element is stored at index 0.
    i_fc[i] corresponds to gas channel i, for i in 1:nb_gc.
"""
function calculate_1D_GC_current_density(i_fc_cell::Float64, sv::AbstractVector{<:MEAState1D}, fc::AbstractFuelCell)::Vector{Float64}

    # Extraction of the parameters
    nb_gc = fc.numerical_parameters.nb_gc
        # Fast path: if there is only one gas channel, the local current density equals the cell current density.
    nb_gc == 1 && return [Float64(i_fc_cell)]
    # Extraction of the variables
    C_O2_ccl = [sv[i].ccl.C_O2 for i in 1:nb_gc]

    # Residual function for NonlinearSolve solver applied on the local current density
    function residuals!(res, x, _)

        # Recovery of the guessed variable values
        U_cell_guessed  = x[1]
        @views i_fc_guessed    = x[2:nb_gc+1]        # view: no copy allocated
        @views C_O2_Pt_guessed = x[nb_gc+2:2*nb_gc+1]  # view: no copy allocated

        # Residuals: difference between guessed and calculated values
        #   Equation set 1 – cell voltage consistency across all GC positions (nb_gc equations)
        @inbounds for i in 1:nb_gc
            res[i] = calculate_cell_voltage(i_fc_guessed[i], C_O2_Pt_guessed[i], sv[i], fc) - U_cell_guessed
        end
        #   Equation set 2 – average current density conservation (1 equation)
        res[nb_gc+1] = i_fc_cell - average(i_fc_guessed)
        #   Equation set 3 – oxygen concentration at the Pt surface (nb_gc equations)
        @inbounds for i in 1:nb_gc
            res[nb_gc+1+i] = calculate_C_O2_Pt(i_fc_guessed[i], sv[i], fc) - C_O2_Pt_guessed[i]
        end
    end

    # Calculation of the 1D GC current density by solving the system of equations defined by the residuals function
    # using NonlinearSolve with Levenberg-Marquardt
    #       Initial guesses
    x0 = Vector(undef, 2 * nb_gc + 1)
    x0[1]            = calculate_cell_voltage(i_fc_cell, C_O2_ccl[1], sv[1], fc)
    x0[2:nb_gc+1]   .= i_fc_cell
    x0[nb_gc+2:2*nb_gc+1] = C_O2_ccl
    #       Solver call
    prob = NonlinearProblem(residuals!, x0, nothing)
    sol  = solve(prob, LevenbergMarquardt())

    #       Check for convergence
    if !successful_retcode(sol.retcode)
        error("Convergence failed in calculate_1D_GC_current_density: retcode = $(sol.retcode)")
    end
    #       Extract the results (indices 2 to nb_gc+1 of the solution vector correspond to i_fc[1:nb_gc])
    i_fc = Vector(sol.u[2:nb_gc+1])

    return i_fc
end

