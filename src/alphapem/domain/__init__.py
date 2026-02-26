# -*- coding: utf-8 -*-

"""AlphaPEM Domain

This module contains the core domain models, parameters, and computational modules that define
the physics and behavior of the PEMFC simulator.

Submodules:
    - models: Core domain models (AlphaPEM fuel cell model)
    - modules: Physical simulation modules (flows, heat transfer, voltage, etc.)
    - parameters: Configuration and parameter definitions
"""

from alphapem.domain.models import AlphaPEM

__all__ = [
    "AlphaPEM",
]

