"""
    AlphaPEM.Interfaces

This module provides user-interface entry points for AlphaPEM.

Available interfaces:
- Web interface: Modern web-based GUI accessible via browser (localhost:8000)
"""
module Interfaces

"""
    launch_web_simulator()

Launch the AlphaPEM web simulator.

This function starts a local Genie web server and opens the simulator in your default browser.

# Usage

```julia
using AlphaPEM
AlphaPEM.Interfaces.launch_web_simulator()
```

Or from the command line:
```bash
julia run_web_app.jl
```

The web interface will be available at: http://localhost:8000/

# Notes

- The server runs on localhost:8000 to ensure local execution and security
- All simulations run on your machine (not on a remote server)
- To stop the server, press Ctrl+C in the terminal
- Data and results remain on your computer
"""
function launch_web_simulator()
    @info "Launching AlphaPEM Web Simulator..."

    # Include the web app entry point
    include(joinpath(@__DIR__, "web", "WebApp.jl"))
end

export launch_web_simulator

end  # module Interfaces

