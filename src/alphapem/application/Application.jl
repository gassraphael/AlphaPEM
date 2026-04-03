# -*- coding: utf-8 -*-

"""
    AlphaPEM.Application

This module contains the application entry points and runtime management for AlphaPEM,
including the main simulation execution interface and orchestration logic.

Modules:
    - run_simulation: Main execution module
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

