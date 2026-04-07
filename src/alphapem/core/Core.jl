# -*- coding: utf-8 -*-

"""
    AlphaPEM.Core

This module contains the core domain models and computational modules that define
PEMFC physics and simulation behavior.

Modules:
    - Types:   Pure data types (structs) — loaded first, no dependencies
    - Models:  Core domain models (main `AlphaPEM` fuel cell model)
    - Modules: Physics kernels (flows, heat transfer, voltage, ODE helpers, display)
"""
module Core

include("types/Types.jl")    # ① types first — no deps
include("models/Models.jl")  # ② models depend on Types
include("modules/Modules.jl")

using .Models: AlphaPEM

export Types, Models, Modules, AlphaPEM

end  # module Core

