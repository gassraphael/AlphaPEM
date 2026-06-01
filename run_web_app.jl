# -*- coding: utf-8 -*-

"""
AlphaPEM Web Simulator - Launch Script

This is the entry point for users to launch the AlphaPEM web interface.

Usage:
    julia run_web_app.jl

After starting, your browser will automatically open to:
    http://localhost:8000/

To stop the server, press Ctrl+C in the terminal.
"""

include(joinpath(@__DIR__, "src", "alphapem", "interfaces", "web", "WebApp.jl"))

