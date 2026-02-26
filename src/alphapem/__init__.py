# -*- coding: utf-8 -*-

"""AlphaPEM: 1D+1D Dynamic Simulator of PEM Fuel Cells for Embedded Applications

This package provides a comprehensive simulation framework for modeling proton exchange membrane
fuel cell (PEMFC) systems. It includes physics-based, finite-volume, pseudo-two-dimensional models
with dynamic two-phase and non-isothermal capabilities.

Main Features:
    - Dynamic simulation of PEMFC internal states and voltage
    - Support for various current density profiles (step, polarization, EIS)
    - Automatic parameter calibration
    - Graphical user interface (GUI)
    - Finite volume method for accurate spatial discretization

Modules:
    - domain: Core data models and parameters
    - application: Simulation execution interface
    - interfaces: GUI and user interfaces
    - utils: Utility functions for mathematics and physics

See Also:
    - Documentation: https://gassraphael.github.io/AlphaPEM/
    - Repository: https://github.com/gassraphael/AlphaPEM
"""

__version__ = "1.3.0"
__author__ = "Raphaël Gass"
__email__ = "gassraphael@proton.me"
__license__ = "GPLv3"

from alphapem.domain.models import AlphaPEM

__all__ = [
    "AlphaPEM",
    "__version__",
    "__author__",
    "__email__",
    "__license__",
]
