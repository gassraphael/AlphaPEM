# -*- coding: utf-8 -*-

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

# Suppress CondaPkg verbose initialization messages
ENV["PIXI_DISABLE_PROGRESS_BAR"] = "1"
let
    original_stdout = stdout
    redirect_stdout(devnull)
    try
        import CondaPkg
    finally
        redirect_stdout(original_stdout)
    end
end

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

