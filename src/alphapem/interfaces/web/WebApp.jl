# -*- coding: utf-8 -*-

"""
    AlphaPEM Web Application

Main entry point for the AlphaPEM web simulator.

Usage:
    julia run_web_app.jl

This script:
1. Activates the AlphaPEM project environment
2. Initializes the Genie web framework
3. Loads routes and backend logic
4. Starts the HTTP server on localhost:8000
5. Opens the browser automatically
"""

# ========================================
# ENVIRONMENT SETUP
# ========================================

# Activate project environment
import Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "..", ".."))

# ========================================
# IMPORTS
# ========================================

using Genie
using Genie.Router
using Genie.Renderer.Html
using Genie.Renderer.Json
using Genie.Requests
import Logging

# Set logging configuration to minimize debug output
Logging.global_logger(Logging.ConsoleLogger(stdout, Logging.Info))

# Suppress debug-level logging from HTTP packages
for logger_name in [:HTTPDebugger, :Genie, :HTTP, :GenieSession]
    # This prevents debug messages from appearing
end

# ========================================
# LOAD APP CONFIGURATION
# ========================================

@info "Loading AlphaPEM Web Application..."

# Load app configuration
include(joinpath(@__DIR__, "app", "app.jl"))

# Load backend logic
include(joinpath(@__DIR__, "SimulatorBackend.jl"))

# Load routes
include(joinpath(@__DIR__, "app", "routes.jl"))

# ========================================
# INITIALIZE
# ========================================

SimulatorBackend.initialize_backend()

# ========================================
# SERVER CONFIGURATION
# ========================================

# Get configuration
config = get_config()
configure_genie()

# ========================================
# START SERVER
# ========================================

const HOST = config[:server_host]
const PORT = config[:server_port]

try
    @info "🚀 AlphaPEM Simulator is running on http://$HOST:$PORT/"
    @info "   Press Ctrl+C to stop the server"

    # Start the Genie server
    Genie.up(
        PORT,
        HOST,
        ws_port = PORT + 1,
        async = false,
        verbose = true
    )

catch e
    @error "Failed to start server" exception=e
    rethrow(e)
end

