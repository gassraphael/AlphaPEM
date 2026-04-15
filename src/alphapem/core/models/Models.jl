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
include("../modules/flows_MEA_1D_modules.jl")
include("../modules/heat_modules.jl")
include("../modules/flows_GC_manifold_1D_modules.jl")
include("../modules/outputs_accessors.jl")
include("../modules/plot_postprocessing_modules.jl")
include("../modules/plotting_modules.jl")

# ── Include model files in dependency order ─────────────────────────────────
include("../modules/dif_eq_modules.jl")
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

function __init__()
    # Rebind matplotlib/numpy Python handles after precompilation.
    copy!(np, pyimport("numpy"))
    copy!(mpl, pyimport("matplotlib"))
    copy!(plt, pyimport("matplotlib.pyplot"))
    copy!(LogLocator, pyimport("matplotlib.ticker").LogLocator)
    publication_palette = [
        "#0072B2", # blue
        "#D55E00", # vermillion
        "#009E73", # bluish green
        "#CC79A7", # reddish purple
        "#56B4E9", # sky blue
        "#E69F00", # orange
        "#882255", # wine
        "#44AA99", # teal
        "#999933", # olive
        "#332288", # indigo
        "#DDCC77", # sand
    ]
    copy!(colors, mpl.colors.ListedColormap(publication_palette; name="alphapem_publication"))
end

end  # module Models

