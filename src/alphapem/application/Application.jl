# -*- coding: utf-8 -*-

"""
    AlphaPEM.Application

Application entry point for simulation execution and runtime plotting setup.

File responsibilities:
- `run_simulation.jl`: public API + orchestration (launch, dispatch, final display)
- `run_simulation_modules.jl`: internal helpers (figure preparation, internal-state utilities)
"""
module Application

using PyCall: PyNULL, pyimport, copy!

const mpl = PyNULL()
const plt = PyNULL()

function __init__()
    # Rebind Python objects at runtime to avoid NULL PyObject after precompilation.
    copy!(mpl, pyimport("matplotlib"))
    copy!(plt, pyimport("matplotlib.pyplot"))
end

include("run_simulation.jl")

export run_simulation

end  # module Application
