"""
    AlphaPEM.Core.Modules
Computational modules for the PEMFC simulator (finite-volume physics kernels).
Includes fluid flow dynamics, heat transfer, electrochemical voltage calculations,
ODE right-hand-side helpers, and display/post-processing utilities.
"""
module Modules
using ...Utils   # AlphaPEM.Utils: constants + maths/physics functions
using ...Config  # AlphaPEM.Config: experimental data (pola_exp_values, ...)
using ...Config: SimulationConfig, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using ...Config: StateScaling, CellStateScaling, ManifoldStateScaling, AuxiliaryStateScaling, DAEAlgebraicScaling
using ...Fuelcell: AbstractFuelCell
using ...Currents: AbstractCurrent, current
using ..Types    # AlphaPEM.Core.Types: all domain structs (MEAThermalIntermediates, …)
using Interpolations
using FFTW
using Statistics

# Include in dependency order
# Fuel Cell specific modules
include("fuelcell/cell_voltage_modules.jl")
include("fuelcell/flows_MEA_1D_modules.jl")
include("fuelcell/heat_modules.jl")
include("fuelcell/flows_GC_manifold_1D_modules.jl")
include("fuelcell/dif_eq_modules.jl")
include("fuelcell/outputs_accessors.jl")
include("fuelcell/plot_postprocess.jl")

# Electrolyzer specific modules (to be added)
# include("electrolyzer/...")

end  # module Modules

