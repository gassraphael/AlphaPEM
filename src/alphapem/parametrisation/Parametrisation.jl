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

Run `run_validity_analysis` to identify a compact hyperbox where AlphaPEM
consistently produces physically meaningful curves.  Results are saved as a
YAML file that can be reused in Stage 2:

```julia
using AlphaPEM.Parametrisation.ValidParameterRegion

cfg = ValidityAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    n_samples      = 2000,
    output_dir     = "results/model_validity",
)
result = run_validity_analysis(cfg)   # generates original_bounds.yaml + classified_configurations.csv

# With PRIM (requires R + IRD package):
prim_cfg = PRIMConfig(
    ird_package_dir       = "external/IRD_method_2023/irdpackage",
    reference_config_path = "results/model_validity/reference_config.yaml",
    output_dir            = "results/model_validity",
)
result = run_validity_analysis(cfg, prim_cfg)  # also generates restricted_bounds_PRIM.yaml
```

**Stage 2 — Calibrate within the restricted space** *(calibration module, work in progress)*

Load the restricted bounds produced in Stage 1 and pass them to the genetic-algorithm
calibration, focusing the search on the high-validity region:

```julia
using AlphaPEM.Parametrisation.ValidParameterRegion: load_restricted_bounds
bounds = load_restricted_bounds("results/model_validity/restricted_bounds_PRIM.yaml")
# → Dict(:Hacl => (7.0e-6, 1.2e-5), :Re => (2.0e-7, 1.8e-6), ...)
# Pass `bounds` to the calibration configuration once that module is re-integrated.
```

## Sub-systems

### Valid Parameter Region (`ValidParameterRegion`)
Identifies the region of the undetermined-parameter space where AlphaPEM produces
physically meaningful polarization curves, using Latin Hypercube Sampling, batch
simulation, and the PRIM/IRD method.  See `examples/run_parameter_validity.jl`.

### Calibration *(work in progress)*
⚠ `calibration_modules.jl` and `calibration.jl` are currently being refactored
and are not yet included here.  They will be re-integrated once the PyCall /
PyGAD dependency is resolved.

## Exports

- `ValidParameterRegion` — main module (contains `ValidityCriteria`,
  `ConfigurationSampling`, `PRIMInterface`, `ResultsExport`)
"""
module Parametrisation

# NOTE: calibration_modules.jl is currently broken (depends on PyCall / PyGAD
# which are not yet ported to pure Julia) and is intentionally excluded until
# the calibration sub-system is refactored.
# include("calibration_modules.jl")

include("valid_parameter_region.jl")

export ValidParameterRegion

end  # module Parametrisation
