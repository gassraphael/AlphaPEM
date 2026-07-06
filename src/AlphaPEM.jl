# -*- coding: utf-8 -*-

"""
	AlphaPEM

1D+1D dynamic simulator of PEM fuel cells for embedded applications.

This module provides package-level access to configuration, core models,
application entry points, calibration helpers, and interfaces.
"""
module AlphaPEM

const VERSION = "2.0.0"
const AUTHOR = "Raphaël Gass"
const EMAIL = "gassraphael@proton.me"
const LICENSE = "BSD-3-Clause"

# Suppress verbose output from CondaPkg and dependencies
ENV["JULIA_CONDAPKG_VERBOSITY"] = "-1"
ENV["PYTHONCALL_EXE_PATH"] = "" # Suppress PythonCall startup messages

import CondaPkg

include("alphapem/utils/Utils.jl")
include("alphapem/config/Config.jl")
include("alphapem/fuelcell/Fuelcell.jl")
include("alphapem/currents/Currents.jl")
include("alphapem/core/Core.jl")
include("alphapem/application/Application.jl")
include("alphapem/interfaces/Interfaces.jl")
include("alphapem/parametrisation/Parametrisation.jl")

export VERSION,
	   AUTHOR,
	   EMAIL,
	   LICENSE,
	   Utils,
	   Config,
	   Fuelcell,
	   Currents,
	   Core,
	   Application,
	   Interfaces,
	   Parametrisation

end  # module AlphaPEM

