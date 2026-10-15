# -*- coding: utf-8 -*-

"""
    ExtendedLocalSensitivityAnalysis

Extended local sensitivity and parameter identifiability analysis of AlphaPEM.

Local sensitivities (normalised finite differences) are evaluated at many points
sampled globally in the parameter space with a Sobol sequence, then analysed through
their magnitude, collinearity and Hessian-based identifiability.

# References
- A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano and T. Ersal, "Effective Parameterization of
  PEM Fuel Cell Models – Part I: Sensitivity Analysis and Parameter Identifiability",
  J. Electrochem. Soc. 167(4), 044504, 2020. DOI: 10.1149/1945-7111/ab7091.
- P. Affonso Nobrega, *marapendi* (PERSEE – Mines Paris – PSL University), MIT licence,
  <https://github.com/affonsop/marapendi>.

The module also hosts the simulation-infrastructure helpers shared with
`SobolSensitivityAnalysis` (warning logger, hashing, checkpoint I/O).
"""
module ExtendedLocalSensitivityAnalysis

using DataFrames
using JLD2
using Logging
using SHA

include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "warning_logger.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "hashing.jl"))
include(joinpath(@__DIR__, "sensitivity_analysis", "shared", "checkpoint_io.jl"))

end  # module ExtendedLocalSensitivityAnalysis
