# -*- coding: utf-8 -*-

"""
Helper functions for simple sensitivity analysis.

This file acts as a lightweight loader: it includes the focused sub-modules stored in
the `simple_sensitivity/` directory.
"""

include(joinpath(@__DIR__, "simple_sensitivity", "regions.jl"))
include(joinpath(@__DIR__, "simple_sensitivity", "polarization.jl"))
include(joinpath(@__DIR__, "simple_sensitivity", "impacts.jl"))
include(joinpath(@__DIR__, "simple_sensitivity", "output.jl"))
