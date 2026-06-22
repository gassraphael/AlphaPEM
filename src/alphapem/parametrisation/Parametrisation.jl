# -*- coding: utf-8 -*-

"""
    AlphaPEM.Parametrisation

Parameter identification tools for AlphaPEM.

## Recommended two-stage workflow

Calibrating AlphaPEM's undetermined parameters (`Hacl`, `Re`, `epsilon_gdl`, …)
directly over the full prior ranges is expensive and brittle, because a large
fraction of that space produces non-physical polarization curves.
The recommended approach is therefore:

**Stage 1 — Restrict the space (ValidParameterRegion)**

 - in construction.

**Stage 2 — Calibrate within the restricted space** *(calibration module, work in progress)*

Load the restricted bounds produced in Stage 1 and pass them to the genetic-algorithm
calibration, focusing the search on the high-validity region:

```julia
using AlphaPEM.Parametrisation.ValidParameterRegion: load_restricted_bounds
    bounds = load_restricted_bounds("results/model_validity/restricted_bounds_PRIM.yaml")
# → Dict(:Hacl => (7.0e-6, 1.2e-5), :Re => (2.0e-7, 1.8e-6), ...)
# Pass `bounds` to the calibration configuration.
```

## Sub-systems

### Valid Parameter Region (`ValidParameterRegion`)
 - in construction

### Calibration (`Calibration`)
GA-based parameter calibration system in pure Julia.

## Exports

- `ValidParameterRegion` — validity analysis module
- `Calibration` — genetic algorithm calibration module
"""
module Parametrisation

include("common/ParametrisationCommon.jl")
include("calibration.jl")

export ParametrisationCommon, Calibration

end  # module Parametrisation
