"""
    AlphaPEM.Core.Models
Main fuel cell model class and provides access to all model-related functionality
for PEMFC simulation.

Types:
    AlphaPEM: Main fuel cell simulator class implementing physics-based modeling
"""
module Models

using DifferentialEquations
using PyCall: PyNULL, pyimport, copy!
using ...Utils   # AlphaPEM.Utils: constants + maths/physics functions
using ...Config  # AlphaPEM.Config: experimental data (pola_exp_values, ...)
using ...Config: SimulationConfig, StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using ...Fuelcell: AbstractFuelCell
using ...Currents: AbstractCurrent, current, delta_t_load
using ..Types    # AlphaPEM.Core.Types: all domain structs (cell_state, cell_derivative, …)

const np = PyNULL()
const mpl = PyNULL()
const plt = PyNULL()
const LogLocator = PyNULL()
const colors = PyNULL()

# Explicit aliases for utility constants that may clash with Base exports.
const Text = Utils.Text
const Pext = Utils.Pext
const Phi_ext = Utils.Phi_ext
const y_O2_ext = Utils.y_O2_ext
const interpolate = Utils.interpolate

# ── Include shared helpers exactly once ─────
include("../modules/cell_voltage_modules.jl")
include("../modules/flows_1D_MEA_modules.jl")
include("../modules/heat_modules.jl")
include("../modules/flows_1D_GC_manifold_modules.jl")
include("../modules/display_calc_modules.jl")
include("../modules/display_modules.jl")

# ── Include model files in dependency order ─────────────────────────────────
include("../modules/dif_eq_modules.jl")
include("velocity.jl")
include("cell_voltage.jl")
include("current_distribution_1D_GC.jl")
include("flows_1D_MEA.jl")
include("flows_1D_GC_manifold.jl")
include("heat_transfer.jl")
include("dif_eq_1D_MEA.jl")
include("dif_eq_1D_GC_manifold.jl")
include("dif_eq_auxiliaries.jl")
include("dif_eq.jl")
include("AlphaPEM.jl")

# Public API
export AlphaPEM

function __init__()
    # Rebind matplotlib/numpy Python handles after precompilation.
    copy!(np, pyimport("numpy"))
    copy!(mpl, pyimport("matplotlib"))
    copy!(plt, pyimport("matplotlib.pyplot"))
    copy!(LogLocator, pyimport("matplotlib.ticker").LogLocator)
    copy!(colors, plt.get_cmap("tab20"))
end

end  # module Models

