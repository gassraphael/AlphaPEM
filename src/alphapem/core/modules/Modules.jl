"""
    AlphaPEM.Core.Modules
Computational modules for the PEMFC simulator (finite-volume physics kernels).
Includes fluid flow dynamics, heat transfer, electrochemical voltage calculations,
ODE right-hand-side helpers, and display/post-processing utilities.
"""
module Modules
using ...Utils   # AlphaPEM.Utils: constants + maths/physics functions
using ...Config  # AlphaPEM.Config: experimental data (pola_exp_values, ...)
using PyCall
using Interpolations
using FFTW
using Statistics

# Include in dependency order
include("cell_voltage_modules.jl")
include("flows_1D_MEA_modules.jl")
include("heat_modules.jl")
include("flows_1D_GC_manifold_modules.jl")
include("dif_eq_modules.jl")
include("display_calc_modules.jl")
include("display_modules.jl")

end  # module Modules

