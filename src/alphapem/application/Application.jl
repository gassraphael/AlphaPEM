# -*- coding: utf-8 -*-

"""
    AlphaPEM.Application

Application entry point for simulation execution and runtime plotting setup.

File responsibilities:
- `run_simulation.jl`: public API + orchestration (launch, dispatch, final display)
- `run_simulation_modules.jl`: internal helpers (figure preparation, internal-state utilities)
- `interactive_display.jl`: serving the figures in a web browser (menu, downloads)
- `run_simulation_live.jl`: segmented `:live` runs, a debugging aid kept apart
"""
module Application

using CairoMakie
import WGLMakie
import HTTP

include("run_simulation.jl")

export run_simulation, prepare_web_figures,
       serve_interactive_figure!, wait_for_interactive_figures!

end  # module Application
