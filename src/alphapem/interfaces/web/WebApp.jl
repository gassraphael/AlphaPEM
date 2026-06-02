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

# Custom logger filter to suppress Info messages from HTTP/WebSocket packages
struct FilterLogger <: Logging.AbstractLogger
    logger::Logging.AbstractLogger
end

function Logging.handle_message(logger::FilterLogger, level, message, _module, group, id, file, line; kwargs...)
    # Suppress Info messages from specific packages
    if level == Logging.Info
        module_name = string(_module)
        # Filter out messages from HTTP, WebSocket, and Genie packages during startup
        if any(occursin(pkg, module_name) for pkg in ["HTTP", "OpenSSL", "Genie"])
            return
        end
    end
    Logging.handle_message(logger.logger, level, message, _module, group, id, file, line; kwargs...)
end

Logging.shouldlog(logger::FilterLogger, level, _module, group, id) =
    Logging.shouldlog(logger.logger, level, _module, group, id)

Logging.min_enabled_level(logger::FilterLogger) = Logging.min_enabled_level(logger.logger)

# Apply the custom logger
Logging.global_logger(FilterLogger(Logging.ConsoleLogger(stdout, Logging.Info)))

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
    @info "   AlphaPEM Simulator is running on http://$HOST:$PORT/"
    @info "   Press Ctrl+C to stop the server"

    # Start the Genie server
    Genie.up(
        PORT,
        HOST,
        ws_port = PORT + 1,
        async = false,
        verbose = false
    )

catch e
    @error "Failed to start server" exception=e
    rethrow(e)
end

