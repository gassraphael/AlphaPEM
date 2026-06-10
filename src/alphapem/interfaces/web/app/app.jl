# -*- coding: utf-8 -*-

"""
    AlphaPEM Web Application Configuration

This module configures the Genie web framework for the AlphaPEM simulator.
All settings for the web server, routes, and middleware are defined here.
"""

using Genie, Stipple, StippleUI
import Logging

# ------ GENIE SERVER CONFIGURATION ------

"""
Configure the Genie server with appropriate settings for AlphaPEM simulator.

# Settings:
- server_host: Listen address (0.0.0.0 = all interfaces, localhost = local only)
- server_port: Port number (8000 = standard)
- server_workers: Number of worker threads
- log_level: Logging verbosity
- websockets_enforced: Security setting for WebSocket connections
"""

# Configure Genie settings
const GENIE_CONFIG = Dict(
    :server_host => "127.0.0.1",           # localhost only (secure)
    :server_port => parse(Int, get(ENV, "ALPHAPEM_PORT", "8000")), # Port from env or default
    :server_workers => 4,                  # Worker threads
    :log_level => Logging.Warn,            # Verbosity level
    :env => "dev",                         # Development mode
    :websockets_enforced => false,         # No WebSocket enforcement
)

"""
Initialize Genie configuration based on environment.
"""
function configure_genie()
    for (key, value) in GENIE_CONFIG
        if key == :server_host
            Genie.config.server_host = value
        elseif key == :server_port
            Genie.config.server_port = value
        elseif key == :log_level
            Genie.config.log_level = value
        end
    end

    return nothing
end

"""
Get the Genie configuration dictionary.
"""
function get_config()::Dict
    return GENIE_CONFIG
end

export configure_genie, get_config

