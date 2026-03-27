# -*- coding: utf-8 -*-

"""
    AlphaPEM.Core

This module contains the core domain models and computational modules that define
PEMFC physics and simulation behavior.

Modules:
    - Models: Core domain models (main `AlphaPEM` fuel cell model)
    - Modules: Physics kernels (flows, heat transfer, voltage, ODE helpers, display)
"""
module Core

include("models/Models.jl")
include("modules/Modules.jl")

using .Models: AlphaPEM

export Models, Modules, AlphaPEM

end  # module Core

