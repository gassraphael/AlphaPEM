# -*- coding: utf-8 -*-

"""
    AlphaPEM Web Application Routes

This module defines all HTTP endpoints for the AlphaPEM simulator web interface.
Routes handle:
- Static page serving (HTML/CSS/JS)
- API endpoints for simulations
- Parameter validation and configuration
"""

using Genie, Genie.Router, Genie.Renderer, Genie.Requests, HTTP
using JSON

include("../SimulatorBackend.jl")

# ------ PAGES DIR ------
const PAGES_DIR = joinpath(@__DIR__, "..", "pages")
const PUBLIC_DIR = joinpath(@__DIR__, "..", "public")

# ------ HELPER FUNCTIONS ------

"""
Return JSON response with proper HTTP headers.
"""
function json(data::Dict; status=200)
    json_string = JSON.json(data)
    HTTP.Response(status, ["Content-Type" => "application/json; charset=utf-8"], json_string)
end

function json(data; status=200)
    json_string = JSON.json(data)
    HTTP.Response(status, ["Content-Type" => "application/json; charset=utf-8"], json_string)
end

# ------ STATIC PAGES ------

# Serve the main index page (landing page of the simulator).
route("/"; method="GET") do
    page_path = joinpath(PAGES_DIR, "index.html")
    page_content = read(page_path, String)
    HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"], page_content)
end

# Serve the simulator main page with tabs and controls.
route("/simulator"; method="GET") do
    try
        page_content = read(joinpath(PAGES_DIR, "simulator.html"), String)
        HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"], page_content)
    catch e
        @error "Error serving simulator page" error=e
        json(Dict(:error => "Failed to load simulator page: $(string(e))", :status => 500))
    end
end

# Serve results page after simulation completion.
route("/results/:result_id"; method="GET") do
    result_id::String = Router.params()[:result_id]
    page_content = read(joinpath(PAGES_DIR, "results.html"), String)
    HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"], page_content)
end

# ------ CSS STYLESHEET ------

# Serve CSS stylesheet.
route("/css/style.css"; method="GET") do
    css_path = joinpath(PUBLIC_DIR, "css", "style.css")
    if isfile(css_path)
        css_content = read(css_path, String)
        HTTP.Response(200, ["Content-Type" => "text/css; charset=utf-8"], css_content)
    else
        json(Dict(:error => "CSS stylesheet not found"))
    end
end

# ------ API ENDPOINTS FOR SIMULATION ------

# API endpoint: Get available fuel cell presets.
# Returns a JSON dictionary with fuel cell options and their default parameters.
route("/api/fuel_cells"; method="GET") do
    fuel_cells = SimulatorBackend.get_available_fuel_cells()
    json(fuel_cells)
end

# API endpoint: Get default parameters for a selected fuel cell.
#
# Parameters:
# - fuel_cell_type: Symbol or String identifying the fuel cell model
#
# Response:
# JSON object containing all default parameters for the selected fuel cell.
route("/api/fuel_cell_defaults/:fuel_cell_type"; method="GET") do
    fuel_cell_type::String = Router.params()[:fuel_cell_type]

    try
        defaults = SimulatorBackend.get_fuel_cell_defaults(fuel_cell_type)
        json(defaults)
    catch e
        json(Dict(:error => "Invalid fuel cell type: $(string(e))"))
    end
end

# API endpoint: Validate user input parameters.
#
# Before launching a simulation, parameters must be validated to ensure:
# - All values are within physically acceptable ranges
# - Required parameters are provided
# - Unit conversions are correct
#
# Request body: JSON object with all parameters
#
# Response:
# - If valid: {:valid => true}
# - If invalid: {:valid => false, :errors => [...]}
route("/api/validate_parameters"; method="POST") do
    try
        params = JSON.parse(String(Requests.rawpayload()))
        validation_result = SimulatorBackend.validate_parameters(params)

        json(validation_result)
    catch e
        json(Dict(:error => "Parameter validation failed: $(string(e))", :valid => false))
    end
end

# API endpoint: Launch a STEP CURRENT simulation.
#
# Simulates a step current density change and monitors cell response.
#
# Request body:
# - operating_conditions: Dict with T, Pa, Pc, Sa, Sc, Phi_a, Phi_c, y_H2_in
# - accessible_parameters: Dict with geometric/volume parameters
# - undetermined_parameters: Dict with physical model parameters
# - computing_parameters: Dict with numerical solver settings
# - current_params: Dict with step current profile (Δt_ini, Δt_load, Δt_break, i_step)
#
# Response:
# - If successful: {:success => true, :result_id => "...", :message => "..."}
# - If error: {:success => false, :error => "..."}
route("/api/simulate/step"; method="POST") do
    try
        params = JSON.parse(String(Requests.rawpayload()))

        # Validate before running
        validation = SimulatorBackend.validate_parameters(params)
        if !validation[:valid]
            return json(Dict(:success => false, :error => validation[:errors]))
        end

        # Build simulation configuration
        config = SimulatorBackend.build_simulation_config(params, :step)

        # Run simulation
        result = SimulatorBackend.run_step_simulation(config)

        json(Dict(:success => true, :result_id => result[:id], :message => "Step simulation completed"))
    catch e
        json(Dict(:success => false, :error => "Step simulation failed: $(string(e))"))
    end
end

# API endpoint: Launch a POLARIZATION CURVE simulation.
#
# Generates a polarization curve by sweeping current density.
#
# Request body: Same structure as step simulation
#
# Additional parameters:
# - delta_i_pola: Current density step size
# - i_max_pola: Maximum current density to reach
#
# Response:
# - If successful: {:success => true, :result_id => "...", :message => "..."}
# - If error: {:success => false, :error => "..."}
route("/api/simulate/polarization"; method="POST") do
    try
        params = JSON.parse(String(Requests.rawpayload()))

        validation = SimulatorBackend.validate_parameters(params)
        if !validation[:valid]
            return json(Dict(:success => false, :error => validation[:errors]))
        end

        config = SimulatorBackend.build_simulation_config(params, :polarization)
        result = SimulatorBackend.run_polarization_simulation(config)

        json(Dict(:success => true, :result_id => result[:id], :message => "Polarization curve completed"))
    catch e
        json(Dict(:success => false, :error => "Polarization simulation failed: $(string(e))"))
    end
end

# API endpoint: Launch an EIS (Electrochemical Impedance Spectroscopy) simulation.
#
# Performs frequency sweep impedance measurement.
#
# Request body: Same structure as other simulations
#
# Additional parameters:
# - i_EIS: DC current density for EIS
# - ratio_EIS: AC current amplitude as ratio of i_EIS
# - f_power_min_EIS: Minimum frequency (as power of 10)
# - f_power_max_EIS: Maximum frequency (as power of 10)
# - nb_f_EIS: Number of frequency points
# - nb_points_EIS: Number of data points per frequency
#
# Response:
# - If successful: {:success => true, :result_id => "...", :message => "..."}
# - If error: {:success => false, :error => "..."}
route("/api/simulate/eis"; method="POST") do
    try
        params = JSON.parse(String(Requests.rawpayload()))

        validation = SimulatorBackend.validate_parameters(params)
        if !validation[:valid]
            return json(Dict(:success => false, :error => validation[:errors]))
        end

        config = SimulatorBackend.build_simulation_config(params, :eis)
        result = SimulatorBackend.run_eis_simulation(config)

        json(Dict(:success => true, :result_id => result[:id], :message => "EIS simulation completed"))
    catch e
        json(Dict(:success => false, :error => "EIS simulation failed: $(string(e))"))
    end
end

# API endpoint: Get simulation status/progress.
#
# Allows frontend to check if a simulation is still running.
#
# Parameters:
# - result_id: The simulation result ID to check
#
# Response:
# {
#     :status => "running" | "completed" | "failed",
#     :progress => 0-100 (percentage),
#     :message => "Current status message"
# }
route("/api/status/:result_id"; method="GET") do
    result_id::String = Router.params()[:result_id]

    try
        status = SimulatorBackend.get_simulation_status(result_id)
        json(status)
    catch e
        json(Dict(:status => "error", :message => "Status check failed: $(string(e))"))
    end
end

# API endpoint: Get results data for a completed simulation.
#
# Parameters:
# - result_id: The simulation result ID
#
# Response:
# JSON object containing:
# - simulation_type: "step" | "polarization" | "eis"
# - status: Current status
# - start_time: When simulation started
# - elapsed_time: How long it took
# - data: Time series and other data
# - kpis: Key performance indicators
route("/api/results/:result_id"; method="GET") do
    result_id::String = Router.params()[:result_id]

    try
        results = SimulatorBackend.get_detailed_results(result_id)
        json(results)
    catch e
        json(Dict(:error => "Result retrieval failed: $(string(e))", :status => 500))
    end
end

# API endpoint: Download result in specified format.
#
# Parameters:
# - result_id: The simulation result ID
# - format: Export format (csv, json, xlsx) - optional, defaults to json
#
# Response:
# File with appropriate content type
route("/api/download/:result_id"; method="GET") do
    result_id::String = Router.params()[:result_id]
    format = get(Router.params(), :format, "json")  # Default to JSON

    try
        # Get export data
        export_data = SimulatorBackend.export_results(result_id, format)

        # Determine content type and extension
        content_type, ext = if format == "csv"
            ("text/csv; charset=utf-8", "csv")
        elseif format in ["xlsx", "excel", "xls"]
            ("application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "xlsx")
        else  # json by default
            ("application/json; charset=utf-8", "json")
        end

        # Return file with appropriate headers
        HTTP.Response(
            200,
            ["Content-Type" => content_type,
             "Content-Disposition" => "attachment;filename=alphapem_results_$(result_id).$(ext)"],
            export_data
        )
    catch e
        json(Dict(:error => "Download failed: $(string(e))", :status => 500))
    end
end

# ------ ERROR HANDLING ------

# Note: Catch-all 404 route commented out for now - may conflict with specific routes
# route("/*"; method="GET") do
#     json(:error => "Endpoint not found")
# end

