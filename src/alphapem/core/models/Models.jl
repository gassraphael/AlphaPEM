"""
    AlphaPEM.Core.Models
Main fuel cell model class and provides access to all model-related functionality
for PEMFC simulation.

Types:
    AlphaPEM: Main fuel cell simulator class implementing physics-based modeling
"""
module Models

using DifferentialEquations
using CairoMakie
using NaNMath    # NaNMath avoids hard DomainError on rare non-physical transients (log argument <= 0).
using ...Utils   # AlphaPEM.Utils: constants + maths/physics functions
using ...Config  # AlphaPEM.Config: experimental data (pola_exp_values, ...)
using ...Config: SimulationConfig, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using ...Fuelcell: AbstractFuelCell
using ...Currents: AbstractCurrent, current, delta_t_load, solver_tstops, solver_dtmax
using ..Types    # AlphaPEM.Core.Types: all domain structs (cell_state, cell_derivative, …)

# Explicit aliases for utility constants that may clash with Base exports.
const Text = Utils.Text
const Pext = Utils.Pext
const Phi_ext = Utils.Phi_ext
const y_O2_ext = Utils.y_O2_ext
const interpolate = Utils.interpolate

# ── Include shared helpers exactly once ─────
# Fuel Cell helpers
include("../modules/fuelcell/cell_voltage_modules.jl")
include("../modules/fuelcell/current_distribution_modules.jl")
include("../modules/fuelcell/flows_MEA_1D_modules.jl")
include("../modules/fuelcell/heat_modules.jl")
include("../modules/fuelcell/flows_GC_manifold_1D_modules.jl")
include("../modules/fuelcell/outputs_accessors.jl")
include("../modules/fuelcell/plot_postprocess.jl")
include("../modules/fuelcell/plot_helpers.jl")

# Electrolyzer helpers (to be added)
# include("../modules/electrolyzer/...")

# ── Include model files in dependency order ─────────────────────────────────
# Fuel Cell specific models
include("../modules/fuelcell/dif_eq_modules.jl")
include("fuelcell/plot.jl")
include("fuelcell/velocity.jl")
include("fuelcell/cell_voltage.jl")
include("fuelcell/current_distribution_GC_1D.jl")
include("fuelcell/flows_MEA_1D.jl")
include("fuelcell/flows_GC_manifold_1D.jl")
include("fuelcell/heat_transfer.jl")
include("fuelcell/dif_eq/dif_eq_MEA_1D.jl")
include("fuelcell/dif_eq/dif_eq_GC_manifold_1D.jl")
include("fuelcell/dif_eq/dif_eq_auxiliaries.jl")
include("fuelcell/dif_eq.jl")

# Electrolyzer specific models (to be added)
# include("electrolyzer/...")

# Common components
include("dae_jacobian.jl")
include("AlphaPEM.jl")

# Public API
export AlphaPEM

end  # module Models

