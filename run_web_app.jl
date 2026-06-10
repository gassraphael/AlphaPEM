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

# Suppress GTK warnings and errors
ENV["GTK_MODULES"] = ""
ENV["GTK_DEBUG"] = ""

# Handle optional port argument
if length(ARGS) > 0
    try
        # Validate that the argument is a number
        port = parse(Int, ARGS[1])
        ENV["ALPHAPEM_PORT"] = string(port)
    catch
        println("Warning: Invalid port argument '$(ARGS[1])'. Using default port 8000.")
    end
end

include(joinpath(@__DIR__, "src", "alphapem", "interfaces", "web", "WebApp.jl"))

