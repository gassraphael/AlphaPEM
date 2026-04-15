# -*- coding: utf-8 -*-

"""
    AlphaPEM.Application

Application entry point for simulation execution and runtime plotting setup.

File responsibilities:
- `run_simulation.jl`: public API + orchestration (launch, dispatch, final display)
- `run_simulation_modules.jl`: internal helpers (figure preparation, internal-state utilities)
"""
module Application

using CairoMakie
import GLMakie

include("run_simulation.jl")

export run_simulation

end  # module Application
