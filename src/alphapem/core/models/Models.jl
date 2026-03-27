"""
    AlphaPEM.Core.Models
Main fuel cell model class and provides access to all model-related functionality
for PEMFC simulation.

Types:
    AlphaPEM: Main fuel cell simulator class implementing physics-based modeling
"""
module Models

using DifferentialEquations
using ...Utils   # AlphaPEM.Utils: constants + maths/physics functions
using ...Config  # AlphaPEM.Config: experimental data (pola_exp_values, ...)

# ── Include in dependency order ──────────────────────────────────────────────
include("dif_eq.jl")
include("velocity.jl")
include("cell_voltage.jl")
include("current_distribution_1D_GC.jl")
include("AlphaPEM.jl")

# Public API
export AlphaPEM

end  # module Models

