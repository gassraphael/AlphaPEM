# -*- coding: utf-8 -*-

"""
    AlphaPEM.Application

This module contains the application entry points and runtime management for AlphaPEM,
including the main simulation execution interface and orchestration logic.

Modules:
    - run_simulation: Main execution module
"""
module Application

include("run_simulation.jl")

export run_simulation

end  # module Application

