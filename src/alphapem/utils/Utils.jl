# -*- coding: utf-8 -*-

"""
    AlphaPEM.Utils

This module provides mathematical and physical utility functions used throughout
AlphaPEM, including constants and transport-property helpers.

Modules:
    - maths_functions: Mathematical helpers for model discretization
    - physics_functions: Physical property calculations and closures
"""
module Utils

include("physics_constants.jl")
include("maths_functions.jl")
include("physics_functions.jl")

# Re-export all symbols defined by the included utility files.
for _name in names(@__MODULE__; all=false, imported=false)
    if _name in (:Utils, :eval, :include) || startswith(String(_name), "#")
        continue
    end
    @eval export $_name
end

end  # module Utils

