# -*- coding: utf-8 -*-

"""AlphaPEM Application Layer

This module contains the application entry points and runtime management for AlphaPEM,
including the main simulation execution interface and orchestration logic.

Modules:
    - run_simulation: Main execution module with CLI entry point
"""

from alphapem.application.run_simulation import run_simulation

__all__ = [
    "run_simulation",
]

