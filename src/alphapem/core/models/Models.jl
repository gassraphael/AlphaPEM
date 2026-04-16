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
using ...Utils   # AlphaPEM.Utils: constants + maths/physics functions
using ...Config  # AlphaPEM.Config: experimental data (pola_exp_values, ...)
using ...Config: SimulationConfig, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using ...Fuelcell: AbstractFuelCell
using ...Currents: AbstractCurrent, current, delta_t_load
using ..Types    # AlphaPEM.Core.Types: all domain structs (cell_state, cell_derivative, …)

# Explicit aliases for utility constants that may clash with Base exports.
const Text = Utils.Text
const Pext = Utils.Pext
const Phi_ext = Utils.Phi_ext
const y_O2_ext = Utils.y_O2_ext
const interpolate = Utils.interpolate

# ── Include shared helpers exactly once ─────
include("../modules/cell_voltage_modules.jl")
include("../modules/flows_MEA_1D_modules.jl")
include("../modules/heat_modules.jl")
include("../modules/flows_GC_manifold_1D_modules.jl")
include("../modules/outputs_accessors.jl")
include("../modules/plot_postprocess.jl")
include("../modules/plot_helpers.jl")

# ── Include model files in dependency order ─────────────────────────────────
include("../modules/dif_eq_modules.jl")
include("plot.jl")
include("velocity.jl")
include("cell_voltage.jl")
include("current_distribution_GC_1D.jl")
include("flows_MEA_1D.jl")
include("flows_GC_manifold_1D.jl")
include("heat_transfer.jl")
include("dif_eq_MEA_1D.jl")
include("dif_eq_GC_manifold_1D.jl")
include("dif_eq_auxiliaries.jl")
include("dif_eq.jl")
include("AlphaPEM.jl")

# Public API
export AlphaPEM

end  # module Models

