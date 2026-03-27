# -*- coding: utf-8 -*-

__precompile__(false)

"""
	AlphaPEM

1D+1D dynamic simulator of PEM fuel cells for embedded applications.

This module provides package-level access to configuration, core models,
application entry points, calibration helpers, and interfaces.
"""
module AlphaPEM

const VERSION = "1.4.0"
const AUTHOR = "Raphael Gass"
const EMAIL = "gassraphael@proton.me"
const LICENSE = "GPLv3"

include("alphapem/utils/Utils.jl")
include("alphapem/config/Config.jl")
include("alphapem/core/Core.jl")
include("alphapem/application/Application.jl")
include("alphapem/parametrisation/Parametrisation.jl")
include("alphapem/interfaces/Interfaces.jl")

export VERSION,
	   AUTHOR,
	   EMAIL,
	   LICENSE,
	   Utils,
	   Config,
	   Core,
	   Application,
	   Parametrisation,
	   Interfaces

end  # module AlphaPEM

