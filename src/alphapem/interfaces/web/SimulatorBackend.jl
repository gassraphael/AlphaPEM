# -*- coding: utf-8 -*-

"""
    SimulatorBackend

Backend logic for the AlphaPEM web simulator.

This module manages:
- Parameter management and validation
- Simulation configuration building
- Interface between web routes and AlphaPEM core
- Result storage and retrieval
"""

module SimulatorBackend

export PLOTS_DIR, initialize_backend, run_step_simulation, run_polarization_simulation, run_eis_simulation, get_detailed_results

using AlphaPEM
using AlphaPEM.Config: SimulationConfig, StepParams, PolarizationParams, EISParams, NumericalParams
using AlphaPEM.Config: PhysicalParams, OperatingConditions
using AlphaPEM.Application: run_simulation, prepare_web_figures, figures_preparation, display!
using JSON
using Dates
using Logging
using WGLMakie
import CairoMakie
using Statistics: mean
using XLSX

# Import AlphaPEM internals for data extraction
import AlphaPEM.Core.Models: middle_gdl_index, middle_mpl_index,
                             time_history, derived_outputs,
                             make_Fourier_transformation, C_v_sat, extract_mea_series, 
                             extract_mid_mea_series, extract_derived_series,
                             extract_derived_gc_series, extract_mid_derived_gc_series,
                             polarisation_sampling_indices
import AlphaPEM.Core.Modules: _eis_point, calculate_reynolds_numbers, final_temperature_matrix_celsius
import AlphaPEM.Application: _web_plot_specs
import AlphaPEM.Currents: current, EISCurrent, PolarizationCurrent, PolarizationCalibrationCurrent
import AlphaPEM.Utils: R

const RESULTS_DIR = joinpath(pwd(), "results")
const PLOTS_DIR = joinpath(RESULTS_DIR, "web_plots")
const SIMULATION_RESULTS = Dict{String, Dict}()  # In-memory result cache

# ------ INTERACTIVE PLOT SERVER (Bonito) ------
#
# Option A architecture: we keep a single live Bonito server in the Julia
# process. Each simulation registers one route per figure under
# `/plots/<result_id>/<plot_key>` returning a fresh `Bonito.App(() -> fig)`.
# The frontend iframe points to those URLs, so every plot switch opens a
# brand new Bonito session, which is what guarantees true interactivity
# (zoom, pan, hover, ...) without the spinner / "Plot initialized multiple
# times" issue that any pre-export to HTML would cause.

const BONITO_HOST = "127.0.0.1"
const BONITO_PORT = 9384
const BONITO_SERVER = Ref{Any}(nothing)
const BONITO_LOCK   = ReentrantLock()
# Keep figures alive so the closures registered as Bonito apps don't get GC'd.
const BONITO_FIGURES = Dict{String, Vector{Any}}()
const BONITO_ROUTES  = Dict{String, Vector{String}}()

"""
Lazy-start (and return) the global Bonito server used to host interactive plots.
"""
function _bonito_server()
    if BONITO_SERVER[] === nothing
        Bonito = WGLMakie.Bonito
        lock(BONITO_LOCK) do
            if BONITO_SERVER[] === nothing
                @info "    Starting interactive plot server on http://$(BONITO_HOST):$(BONITO_PORT)/"
                BONITO_SERVER[] = Bonito.Server(BONITO_HOST, BONITO_PORT)
            end
        end
    end
    return BONITO_SERVER[]
end

"""
Build the public base URL of the Bonito server (as seen by the browser).

Reads the host/port from the live server because Bonito may pick a different
port than requested if the preferred one is already in use.
"""
function _bonito_base_url()::String
    server = _bonito_server()
    proto = try
        String(getfield(server, :protocol))
    catch
        "http://"
    end
    host = try
        String(getfield(server, :url))
    catch
        BONITO_HOST
    end
    port = try
        Int(getfield(server, :port))
    catch
        BONITO_PORT
    end
    # `proto` from Bonito already ends with "://"
    return string(proto, host, ":", port)
end

"""
Register one Bonito route per figure for the given simulation `result_id`.

`specs` is the list of `_web_plot_specs(...)` dicts; `figures` is the matching
list of WGLMakie `Figure` objects. Routes are added to the global Bonito
server and the figures are stored in `BONITO_FIGURES` so they stay alive.

Returns an updated copy of `specs` where each entry has an absolute `"url"`
field pointing to the Bonito server. The original `"filename"` / `"title"` /
`"group"` / `"key"` fields are preserved.
"""
function register_bonito_plots!(result_id::String, specs, figures)
    isempty(specs) && return specs

    Bonito = WGLMakie.Bonito
    server = _bonito_server()

    # Drop any previously registered routes/figures for this result_id so a
    # re-run cleanly overwrites the old ones.
    _unregister_bonito_plots!(result_id)

    BONITO_FIGURES[result_id] = collect(figures)
    BONITO_ROUTES[result_id]  = String[]

    enriched = Dict{String, String}[]
    for (i, spec) in enumerate(specs)
        fig = figures[i]
        key = haskey(spec, "key") ? String(spec["key"]) : "plot_$(i)"
        title = haskey(spec, "title") ? String(spec["title"]) : "Plot $(i)"
        route_path = "/plots/$(result_id)/$(key)"

        # Capture `fig` by value so each closure returns its own figure.
        app = Bonito.App(() -> fig; title=title)
        try
            Bonito.route!(server, route_path => app)
        catch e
            @error "    Failed to register Bonito route $route_path" exception=e
            rethrow(e)
        end
        push!(BONITO_ROUTES[result_id], route_path)

        new_spec = Dict{String, String}()
        for (k, v) in spec
            new_spec[String(k)] = String(v)
        end
        new_spec["url"] = string(_bonito_base_url(), route_path)
        new_spec["key"] = key
        push!(enriched, new_spec)
    end

    return enriched
end

"""
Tear down previously registered Bonito routes/figures for `result_id`.
"""
function _unregister_bonito_plots!(result_id::String)
    haskey(BONITO_ROUTES, result_id) || return nothing
    Bonito = WGLMakie.Bonito
    server = BONITO_SERVER[]
    if server !== nothing
        for route_path in BONITO_ROUTES[result_id]
            try
                Bonito.HTTPServer.delete_route!(server, route_path)
            catch e
                @debug "    Could not delete Bonito route $route_path: $e"
            end
        end
    end
    delete!(BONITO_ROUTES,  result_id)
    delete!(BONITO_FIGURES, result_id)
    return nothing
end

"""
Initialize backend directories and caches.
"""
function initialize_backend()
    # Create results directory if it doesn't exist
    if !isdir(RESULTS_DIR)
        mkpath(RESULTS_DIR)
    end
    if !isdir(PLOTS_DIR)
        mkpath(PLOTS_DIR)
    end

    # Ensure WGLMakie is active for the web interface
    WGLMakie.activate!()

    # Eagerly start the Bonito plot server so the first plot switch is fast.
    try
        _bonito_server()
    catch e
        @warn "    Could not start Bonito plot server eagerly: $e"
    end
end

# ------ FUEL CELL PRESETS ------

"""
Metadata for parameters: factor to convert from SI (m, m², etc.) to UI units,
and the target unit string.
"""
const PARAM_UI_CONVERSION = Dict(
    :Aact => (factor=1e4, unit="cm²", precision=1),
    :Hacl => (factor=1e6, unit="μm", precision=1),
    :Hccl => (factor=1e6, unit="μm", precision=1),
    :Hmem => (factor=1e6, unit="μm", precision=1),
    :Hgdl => (factor=1e6, unit="μm", precision=1),
    :Hmpl => (factor=1e6, unit="μm", precision=1),
    :Hagc => (factor=1e6, unit="μm", precision=1),
    :Hcgc => (factor=1e6, unit="μm", precision=1),
    :Wagc => (factor=1e6, unit="μm", precision=1),
    :Wcgc => (factor=1e6, unit="μm", precision=1),
    :Lgc  => (factor=1e3, unit="mm", precision=1),
    :Ldist => (factor=1e3, unit="mm", precision=1),
    :Lm   => (factor=1e3, unit="mm", precision=1),
    :A_T_a => (factor=1e4, unit="cm²", precision=1),
    :A_T_c => (factor=1e4, unit="cm²", precision=1),
    :Vasm => (factor=1e6, unit="cm³", precision=1),
    :Vcsm => (factor=1e6, unit="cm³", precision=1),
    :Vaem => (factor=1e6, unit="cm³", precision=1),
    :Vcem => (factor=1e6, unit="cm³", precision=1),
    :Re   => (factor=1e6, unit="μΩ·m²", precision=3),
    :epsilon_gdl => (factor=1.0, unit="", precision=3),
    :K_O2_ad_Pt => (factor=1.0, unit="", precision=3),
    :kappa_c => (factor=1.0, unit="", precision=3),
    :i0_c_ref => (factor=1.0, unit="A/m²", precision=2),
    :kappa_co => (factor=1.0, unit="mol·m⁻¹·s⁻¹·Pa⁻¹", precision=2),
    :C_scl => (factor=1e-6, unit="MF/m³", precision=0),
    # Run Parameters
    :delta_t_ini => (factor=1/60, unit="min", precision=1),
    :delta_t_load => (factor=1.0, unit="s", precision=1),
    :delta_t_break => (factor=1/60, unit="min", precision=1),
    :i_ini => (factor=1e-4, unit="A/cm²", precision=2),
    :i_step => (factor=1e-4, unit="A/cm²", precision=2),
    :di_step => (factor=1e-4, unit="A/cm²", precision=2),
    :i_max => (factor=1e-4, unit="A/cm²", precision=2),
    :i_EIS => (factor=1e-4, unit="A/cm²", precision=2),
    :i_static => (factor=1e-4, unit="A/cm²", precision=2),
    :v_load => (factor=1e-4, unit="A/cm²·s", precision=3),
)

"""
Get list of available fuel cell models.

Returns a dictionary with fuel cell types and their descriptions.
"""
function get_available_fuel_cells()::Dict
    fuel_cells = Dict(
        :ZSW_GenStack => Dict(
            :name => "ZSW",
            :description => "ZSW GenStack",
        ),
        :EH31_2022 => Dict(
            :name => "EH-31 (2022)",
            :description => "EH-31 fuel cell",
        ),
    )

    return fuel_cells
end

"""
Get default parameters for a specific fuel cell type.

# Arguments:
- fuel_cell_type: Symbol or String identifying the fuel cell

# Returns:
Dictionary with all default operating conditions and parameters.
"""
function get_fuel_cell_defaults(fuel_cell_type::String)::Dict
    # Handle mapping to internal symbols
    mapped_type = if fuel_cell_type == "ZSW_GenStack"
        :ZSW_GenStack
    elseif fuel_cell_type == "EH31_2022"
        :EH31_2022
    else
        Symbol(fuel_cell_type)
    end

    # Load parameters using AlphaPEM functions
    try
        # Use the factory to create a fuel cell instance
        # This ensures we get exactly the parameters defined in fuelcell/*.jl
        fc = AlphaPEM.Fuelcell.create_fuelcell(mapped_type, :before_voltage_drop)
        
        # Get undetermined parameters list for this fuel cell
        und_params_list = AlphaPEM.Fuelcell.undetermined_parameters(fc, :before_voltage_drop)
        und_param_keys = [p[1] for p in und_params_list]

        ap = fc.physical_parameters
        oc = fc.operating_conditions
        
        # Helper to convert struct to Dict
        function struct_to_dict(obj)
            d = Dict{Symbol, Any}()
            for field in fieldnames(typeof(obj))
                d[field] = getfield(obj, field)
            end
            return d
        end
        
        ap_dict = struct_to_dict(ap)
        
        accessible = Dict()
        undetermined = Dict()
        
        for (k, v) in ap_dict
            # Use the keys from undetermined_parameters() if available, 
            # otherwise fallback to generic detection
            if k in und_param_keys
                undetermined[k] = v
            else
                accessible[k] = v
            end
        end

        return Dict(
            :operating_conditions => Dict(
                :T_fc => oc.T_des - 273.15,
                :Pa => oc.Pa_des / 1e5,
                :Pc => oc.Pc_des / 1e5,
                :Sa => oc.Sa,
                :Sc => oc.Sc,
                :Phi_a => oc.Phi_a_des,
                :Phi_c => oc.Phi_c_des,
                :y_H2_in => oc.y_H2_in,
            ),
            :accessible_parameters => accessible,
            :undetermined_parameters => undetermined,
            :undetermined_list => und_param_keys, # Explicitly send the list for UI highlighting
            :param_metadata => PARAM_UI_CONVERSION,
            :computing_parameters => Dict(
                :nb_gc => 1,
                :nb_gdl => 3,
                :nb_mpl => 2,
                :rtol => 1e-3,
                :atol => 1e-6,
                :maxiters => 100000,
            ),
            :model_config => Dict(
                :voltage_zone => "before_voltage_drop",
                :type_auxiliary => "no_auxiliary",
                :type_flow => "counter_flow",
                :type_purge => "no_purge",
            ),
            :step_parameters => struct_to_dict(StepParams()),
            :polarization_parameters => struct_to_dict(PolarizationParams()),
            :eis_parameters => struct_to_dict(EISParams()),
            :param_metadata => PARAM_UI_CONVERSION,
        )
    catch e
        @warn "Error loading fuel cell defaults for $mapped_type" exception=e
    end

    # Fallback to generic defaults (units matched to UI expectations)
    defaults = Dict(
        :operating_conditions => Dict(
            :T_fc => 60.0,              # °C - Cell temperature
            :Pa => 1.5,                 # bar - Anode pressure
            :Pc => 1.5,                 # bar - Cathode pressure
            :Sa => 2.0,                 # Anode stoichiometry
            :Sc => 2.0,                 # Cathode stoichiometry
            :Phi_a => 0.6,              # Anode humidity
            :Phi_c => 0.6,              # Cathode humidity
            :y_H2_in => 0.95,           # Anode inlet H2 ratio
        ),
        :step_parameters => Dict(
            :delta_t_ini => 1800.0,
            :delta_t_load => 30.0,
            :delta_t_break => 120.0,
            :i_ini => 10000.0,
            :i_step => 20000.0,
        ),
        :polarization_parameters => Dict(
            :delta_t_ini => 7200.0,
            :di_step => 500.0,
            :v_load => 100.0,
            :delta_t_break => 900.0,
            :i_max => 25000.0,
        ),
        :eis_parameters => Dict(
            :i_EIS => 10000.0,
            :ratio => 0.05,
            :f_power_min => -3.0,
            :f_power_max => 5.0,
            :nb_f => 90,
            :nb_points => 50,
        ),
        :model_config => Dict(
            :voltage_zone => "before_voltage_drop",
            :type_auxiliary => "no_auxiliary",
            :type_flow => "counter_flow",
            :type_purge => "no_purge",
        ),
        :accessible_parameters => Dict(
            :Aact => 0.03,              # m² - Active area
            :nb_cell => 1,              # Number of cells
            :Hagc => 500e-6,            # m - Anode GC thickness
            :Hcgc => 500e-6,            # m - Cathode GC thickness
            :Wagc => 450e-6,            # m - Anode GC width
            :Wcgc => 450e-6,            # m - Cathode GC width
            :Lgc => 0.144,              # m - GC length
            :nb_channel_in_gc => 50,    # Number of channels
            :Ldist => 0.05,             # m - Distributor length
            :Lm => 0.025,               # m - Manifold length
            :A_T_a => 1e-3,             # m² - Anode throttle area
            :A_T_c => 1e-3,             # m² - Cathode throttle area
            :Vasm => 5e-5,              # m³ - Supply anode manifold
            :Vcsm => 5e-5,              # m³ - Supply cathode manifold
            :Vaem => 5e-5,              # m³ - Exhaust anode manifold
            :Vcem => 5e-5,              # m³ - Exhaust cathode manifold
        ),
        :undetermined_parameters => Dict(
            :Hgdl => 200e-6,            # m - GDL thickness
            :Hmpl => 50e-6,             # m - MPL thickness
            :Hacl => 10e-6,             # m - Anode CL thickness
            :Hccl => 10e-6,             # m - Cathode CL thickness
            :Hmem => 15e-6,             # m - Membrane thickness
            :epsilon_gdl => 0.7,        # GDL porosity
            :epsilon_mpl => 0.5,        # MPL porosity
            :e => 4,                    # Capillary exponent
            :gamma_sorp_l => 0.5,       # Sorption rate
            :K_O2_ad_Pt => 5.4,         # O2 adsorption coefficient
            :Re => 1.5e-7,              # Ω.m² - Electron resistance
            :i0_c_ref => 15.0,          # A/m² - Exchange current
            :kappa_co => 20.0,          # mol/(m.s.Pa) - Crossover
            :kappa_c => 1.0,            # Overpotential exponent
            :C_scl => 2e7,              # F/m³ - Double layer capacitance
        ),
        :undetermined_list => [:Hacl, :Hccl, :Hmem, :Hgdl, :Hmpl, :epsilon_gdl, :e, :Re, :i0_c_ref, :kappa_co, :kappa_c, :K_O2_ad_Pt, :C_scl],
        :computing_parameters => Dict(
            :nb_gc => 1,                # GC nodes
            :nb_gdl => 5,               # GDL nodes
            :nb_mpl => 5,               # MPL nodes
            :t_purge => 0.6,            # s - Purge time
            :delta_t_purge => 15.0,     # s - Time between purges
            :rtol => 1e-3,              # Solver relative tolerance
            :atol => 1e-6,              # Solver absolute tolerance
        ),
        :step_parameters => Dict(
            :delta_t_ini => 30.0 * 60.0, # s - Initial stabilization
            :delta_t_load => 30.0,       # s - Loading time
            :delta_t_break => 2.0 * 60.0,# s - Break time
            :i_ini => 1.0e4,             # A/m² - Initial current
            :i_step => 1.5e4,            # A/m² - Step current
        ),
        :polarization_parameters => Dict(
            :delta_t_ini => 30.0 * 60.0, # s
            :delta_t_load => 30.0,       # s
            :delta_t_break => 2.0 * 60.0,# s
            :di_step => 0.5e4,           # A/m² - Step size
            :i_max => 2.0e4,             # A/m² - Max current
        ),
        :eis_parameters => Dict(
            :i_static => 1.0e4,          # A/m² - DC current
            :ratio => 0.05,              # AC amplitude ratio
            :f_power_min => 0,           # Min frequency power
            :f_power_max => 4,           # Max frequency power
            :nb_frequencies => 10,       # Number of frequencies
            :nb_points => 20,            # Points per frequency
        ),
        :param_metadata => PARAM_UI_CONVERSION,
    )

    return defaults
end

# ------ PARAMETER VALIDATION ------

"""
Validate all user-provided parameters.

Checks:
- Physical feasibility (ranges, signs, consistency)
- Required fields presence
- Unit consistency

# Returns:
Dict with :valid (bool) and :errors (array) fields
"""
function validate_parameters(params::Dict)::Dict
    errors = String[]

    # Auto-fill missing sections from fuel cell defaults if available
    if haskey(params, :fuel_cell_type)
        try
            defaults = get_fuel_cell_defaults(string(params[:fuel_cell_type]))
            for section in [:operating_conditions, :accessible_parameters,
                           :undetermined_parameters, :computing_parameters]
                if !haskey(params, section)
                    params[section] = get(defaults, section, Dict())
                    @info "Auto-filled missing section: $section from defaults"
                end
            end
        catch e
            @warn "Could not auto-fill parameters from defaults: $e"
        end
    end

    # Check structure
    required_sections = [:operating_conditions, :accessible_parameters,
                        :undetermined_parameters, :computing_parameters]

    for section in required_sections
        if !haskey(params, section)
            push!(errors, "Missing section: $section")
        end
    end

    if !isempty(errors)
        @error "Parameter validation failed" errors=errors
        return Dict(:valid => false, :errors => errors)
    end

    # Validate operating conditions
    oc = params[:operating_conditions]

    if haskey(oc, :T_fc) && (oc[:T_fc] < -273.15 || oc[:T_fc] > 200)
        push!(errors, "Temperature out of range: $oc[:T_fc] °C")
    end

    if haskey(oc, :Pa) && (oc[:Pa] < 0 || oc[:Pa] > 10)
        push!(errors, "Anode pressure out of range: $(oc[:Pa]) bar")
    end

    if haskey(oc, :Pc) && (oc[:Pc] < 0 || oc[:Pc] > 10)
        push!(errors, "Cathode pressure out of range: $(oc[:Pc]) bar")
    end

    if haskey(oc, :Sa) && (oc[:Sa] < 1 || oc[:Sa] > 10)
        push!(errors, "Anode stoichiometry out of range: $(oc[:Sa])")
    end

    if haskey(oc, :Sc) && (oc[:Sc] < 1 || oc[:Sc] > 10)
        push!(errors, "Cathode stoichiometry out of range: $(oc[:Sc])")
    end

    if haskey(oc, :Phi_a) && (oc[:Phi_a] < 0 || oc[:Phi_a] > 1)
        push!(errors, "Anode humidity out of range: $(oc[:Phi_a])")
    end

    if haskey(oc, :Phi_c) && (oc[:Phi_c] < 0 || oc[:Phi_c] > 1)
        push!(errors, "Cathode humidity out of range: $(oc[:Phi_c])")
    end

    # Validate accessible parameters
    ap = params[:accessible_parameters]

    if haskey(ap, :Aact) && ap[:Aact] <= 0
        push!(errors, "Active area must be positive: $(ap[:Aact]) cm²")
    end

    if haskey(ap, :nb_cell) && ap[:nb_cell] < 1
        push!(errors, "Number of cells must be at least 1: $(ap[:nb_cell])")
    end

    # Validate computing parameters
    cp = params[:computing_parameters]

    if haskey(cp, :nb_gc) && cp[:nb_gc] < 1
        push!(errors, "Number of GC nodes must be at least 1: $(cp[:nb_gc])")
    end

    if haskey(cp, :rtol) && (cp[:rtol] <= 0 || cp[:rtol] > 1e-3)
        push!(errors, "Solver rtol out of range: $(cp[:rtol])")
    end

    if haskey(cp, :atol) && (cp[:atol] <= 0 || cp[:atol] > 1e-3)
        push!(errors, "Solver atol out of range: $(cp[:atol])")
    end

    if !isempty(errors)
        return Dict(:valid => false, :errors => errors)
    end

    return Dict(:valid => true, :errors => String[])
end

# ------ SIMULATION CONFIGURATION ------

"""
Build SimulationConfig from web parameters.

# Arguments:
- params: Dictionary with all parameters from web form
- sim_type: Symbol (:step, :polarization, :eis)

# Returns:
AlphaPEM.Config.SimulationConfig object
"""
function build_simulation_config(params::Dict, sim_type::Symbol)::SimulationConfig

    # Build current profile parameters based on simulation type
    current_params = if sim_type == :step
        sp = get(params, :step_parameters, Dict())
        # Merge with defaults if available
        if haskey(params, :fuel_cell_type)
            try
                defaults = get_fuel_cell_defaults(string(params[:fuel_cell_type]))
                for (key, val) in get(defaults, :step_parameters, Dict())
                    if !haskey(sp, key)
                        sp[key] = val
                    end
                end
            catch e
                @warn "Could not merge step parameter defaults: $e"
            end
        end
        StepParams(
            delta_t_ini = sp[:delta_t_ini],
            delta_t_load = sp[:delta_t_load],
            delta_t_break = sp[:delta_t_break],
            i_ini = sp[:i_ini],
            i_step = sp[:i_step],
        )
    elseif sim_type == :polarization
        pp = get(params, :polarization_parameters, Dict())
        # Merge with defaults if available
        if haskey(params, :fuel_cell_type)
            try
                defaults = get_fuel_cell_defaults(string(params[:fuel_cell_type]))
                for (key, val) in get(defaults, :polarization_parameters, Dict())
                    if !haskey(pp, key)
                        pp[key] = val
                    end
                end
            catch e
                @warn "Could not merge polarization parameter defaults: $e"
            end
        end
        PolarizationParams(
            delta_t_ini = pp[:delta_t_ini],
            di_step = pp[:di_step],
            v_load = pp[:v_load],
            delta_t_break = pp[:delta_t_break],
            i_max = pp[:i_max],
        )
    elseif sim_type == :eis
        ep = get(params, :eis_parameters, Dict())
        # Merge with defaults if available
        if haskey(params, :fuel_cell_type)
            try
                defaults = get_fuel_cell_defaults(string(params[:fuel_cell_type]))
                for (key, val) in get(defaults, :eis_parameters, Dict())
                    if !haskey(ep, key)
                        ep[key] = val
                    end
                end
            catch e
                @warn "Could not merge EIS parameter defaults: $e"
            end
        end
        EISParams(
            i_EIS = ep[:i_static],
            ratio = ep[:ratio],
            f_power_min = ep[:f_power_min],
            f_power_max = ep[:f_power_max],
            nb_f = ep[:nb_frequencies],
            nb_points = ep[:nb_points],
        )
    else
        error("Unknown simulation type: $sim_type")
    end

    # Build numerical parameters (from computing_parameters)
    cp = params[:computing_parameters]
    num_params = NumericalParams(
        nb_gc = cp[:nb_gc],
        nb_gdl = cp[:nb_gdl],
        nb_mpl = cp[:nb_mpl],
        rtol = cp[:rtol],
        atol = cp[:atol],
        maxiters = get(cp, :maxiters, 100000), # Increased default to match NumericalParams
    )

    # Build custom physical parameters
    ap_data = get(params, :accessible_parameters, Dict())
    up_data = get(params, :undetermined_parameters, Dict())
    
    # 1. Determine base fuel cell type for defaults
    fc_type_str = string(get(params, :fuel_cell_type, "ZSW_GenStack"))
    v_zone_str = string(get(params, :voltage_zone, "before_voltage_drop"))
    
    # Strip "custom_" prefix if present to get the base model defaults
    base_type_str = replace(fc_type_str, "custom_" => "")
    base_type = Symbol(base_type_str)
    v_zone = Symbol(v_zone_str)

    # 2. Get base physical parameters for this fuel cell type
    base_fc = try
        AlphaPEM.Fuelcell.create_fuelcell(base_type, v_zone)
    catch
        # Fallback to ZSW
        AlphaPEM.Fuelcell.create_fuelcell(:ZSW_GenStack, :before_voltage_drop)
    end
    base_ap = base_fc.physical_parameters
    
    # 3. Start with base parameters as a Dict
    filtered_phys_dict = Dict{Symbol, Any}()
    for field in fieldnames(PhysicalParams)
        filtered_phys_dict[field] = getfield(base_ap, field)
    end
    
    # 4. Override with web data (which is already in SI units)
    for (k, v) in ap_data
        key = Symbol(k)
        if haskey(filtered_phys_dict, key)
            filtered_phys_dict[key] = Float64(v)
        end
    end
    for (k, v) in up_data
        key = Symbol(k)
        if haskey(filtered_phys_dict, key)
            filtered_phys_dict[key] = Float64(v)
        end
    end
    
    custom_ap = PhysicalParams(; filtered_phys_dict...)
    
    # Build custom operating conditions
    oc_data = get(params, :operating_conditions, Dict())
    custom_oc = OperatingConditions(
        T_des     = Float64(get(oc_data, "T_fc", get(oc_data, :T_fc, 60.0))) + 273.15,
        Pa_des    = Float64(get(oc_data, "Pa", get(oc_data, :Pa, 1.5))) * 1e5,
        Pc_des    = Float64(get(oc_data, "Pc", get(oc_data, :Pc, 1.5))) * 1e5,
        Sa        = Float64(get(oc_data, "Sa", get(oc_data, :Sa, 2.0))),
        Sc        = Float64(get(oc_data, "Sc", get(oc_data, :Sc, 2.0))),
        Phi_a_des = Float64(get(oc_data, "Phi_a", get(oc_data, :Phi_a, 0.6))),
        Phi_c_des = Float64(get(oc_data, "Phi_c", get(oc_data, :Phi_c, 0.6))),
        y_H2_in   = Float64(get(oc_data, "y_H2_in", get(oc_data, :y_H2_in, 0.95)))
    )

    # Build model configuration (accept both flattened payload and nested model_config payload)
    model_cfg = get(params, :model_config, Dict())
    voltage_zone = Symbol(get(params, :voltage_zone, get(model_cfg, :voltage_zone, "before_voltage_drop")))
    type_auxiliary = Symbol(get(params, :type_auxiliary, get(model_cfg, :type_auxiliary, "no_auxiliary")))
    type_flow = Symbol(get(params, :type_flow, get(model_cfg, :type_flow, "counter_flow")))
    type_purge = Symbol(get(params, :type_purge, get(model_cfg, :type_purge, "no_purge")))

    # Resolve base fuel-cell symbol to ensure factory compatibility. Synthetic IDs
    # like "custom_*" or "default" are mapped back to known types (ZSW/EH31)
    # to avoid falling back to DefaultFuelCell (which causes NaN derivatives).
    resolved_fuel_cell_type = if base_fc isa AlphaPEM.Fuelcell.ZSWFuelCell ||
                                 base_fc isa AlphaPEM.Fuelcell.EH31FuelCell
        base_type
    else
        :ZSW_GenStack
    end

    config = SimulationConfig(
        type_fuel_cell = resolved_fuel_cell_type,
        type_current = current_params,
        numerical_parameters = num_params,
        voltage_zone = voltage_zone,
        type_auxiliary = type_auxiliary,
        type_flow = type_flow,
        type_purge = type_purge,
        type_display = :no_display,     # Server-side, don't display plots
        display_timing = :postrun,      # Display after simulation
        physical_parameters = custom_ap,
        operating_conditions = custom_oc,
    )

    return config
end

# ------ SIMULATION EXECUTION ------

"""
Run a step current simulation.

# Arguments:
- config: SimulationConfig object

# Returns:
Dictionary with :id and :status
"""
function run_step_simulation(config::SimulationConfig; params=nothing)::Dict
    result_id = generate_result_id("step")
    start_time = time()
    if !isnothing(params) && haskey(params, :req_start_time)
        start_time = params[:req_start_time]
    end

    try
        # Capture the actual output from the simulation
        output = run_simulation(config)

         # Prepare WGLMakie figures and host them via the live Bonito server.
         plot_files = try
             specs, figs = prepare_web_figures(output)
             enriched = register_bonito_plots!(result_id, specs, figs)
             @debug "    Registered $(length(enriched)) Bonito plot routes for $result_id"
             enriched
         catch e
             @error "    Failed to prepare/register web plots for step simulation" exception=e
              Dict{String, String}[]
         end

        elapsed = time() - start_time

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "step",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => config,
            :input_params => params,  # Store original input parameters
            :output => output,  # Store the actual simulation results
            :plot_files => plot_files,
        )

        @info "    Step simulation completed: $result_id in $(round(elapsed, digits=2))s"

        return Dict(:id => result_id, :status => "completed")

    catch e
        elapsed = time() - start_time
        @error "    Step simulation failed: $result_id" exception=e

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "step",
            :status => "failed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :error => string(e),
        )

        return Dict(:id => result_id, :status => "failed", :error => string(e))
    end
end

"""
Run a polarization curve simulation.
"""
function run_polarization_simulation(config::SimulationConfig; params=nothing)::Dict
    result_id = generate_result_id("polarization")
    start_time = time()
    if !isnothing(params) && haskey(params, :req_start_time)
        start_time = params[:req_start_time]
    end

    try
        # Capture the actual output from the simulation
        output = run_simulation(config)

         # Prepare WGLMakie figures and host them via the live Bonito server.
         plot_files = try
             specs, figs = prepare_web_figures(output)
             enriched = register_bonito_plots!(result_id, specs, figs)
             @debug "    Registered $(length(enriched)) Bonito plot routes for $result_id"
             enriched
         catch e
             @error "    Failed to prepare/register web plots for polarization simulation" exception=e
              Dict{String, String}[]
         end

        elapsed = time() - start_time

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "polarization",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => config,
            :input_params => params,  # Store original input parameters
            :output => output,  # Store the actual simulation results
            :plot_files => plot_files,
        )

        @info "    Polarization simulation completed: $result_id in $(round(elapsed, digits=2))s"

        return Dict(:id => result_id, :status => "completed")

    catch e
        elapsed = time() - start_time
        @error "    Polarization simulation failed: $result_id" exception=e

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "polarization",
            :status => "failed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :error => string(e),
        )

        return Dict(:id => result_id, :status => "failed", :error => string(e))
    end
end

"""
Run an EIS simulation.
"""
function run_eis_simulation(config::SimulationConfig; params=nothing)::Dict
    result_id = generate_result_id("eis")
    start_time = time()
    if !isnothing(params) && haskey(params, :req_start_time)
        start_time = params[:req_start_time]
    end

    try
        eis_config = SimulationConfig(
            type_fuel_cell     = config.type_fuel_cell,
            type_current       = config.type_current,
            numerical_parameters = config.numerical_parameters,
            voltage_zone       = config.voltage_zone,
            type_auxiliary     = config.type_auxiliary,
            type_flow          = config.type_flow,
            type_purge         = config.type_purge,
            type_display       = :no_display,
            display_timing     = :postrun,
        )
        # Capture the actual output from the simulation
        output = run_simulation(eis_config)

         # Prepare WGLMakie figures and host them via the live Bonito server.
         plot_files = try
             specs, figs = prepare_web_figures(output)
             enriched = register_bonito_plots!(result_id, specs, figs)
             @debug "    Registered $(length(enriched)) Bonito plot routes for $result_id"
             enriched
         catch e
             @error "    Failed to prepare/register web plots for EIS simulation" exception=e
             Dict{String, String}[]
         end

        elapsed = time() - start_time

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "eis",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => eis_config,
            :input_params => params,  # Store original input parameters
            :output => output,  # Store the actual simulation results
            :plot_files => plot_files,
        )

        @info "    EIS simulation completed: $result_id in $(round(elapsed, digits=2))s"

        return Dict(:id => result_id, :status => "completed")

    catch e
        elapsed = time() - start_time
        @error "    EIS simulation failed: $result_id" exception=e

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "eis",
            :status => "failed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :error => string(e),
        )

        return Dict(:id => result_id, :status => "failed", :error => string(e))
    end
end

# ------ RESULT RETRIEVAL ------

"""
Get the status of a simulation.
"""
function get_simulation_status(result_id::String)::Dict
    if haskey(SIMULATION_RESULTS, result_id)
        result = SIMULATION_RESULTS[result_id]
        return Dict(
            :status => result[:status],
            :type => result[:type],
            :message => "$(result[:type]) simulation: $(result[:status])",
            :elapsed_time => get(result, :elapsed_time, -1.0),
        )
    else
        return Dict(
            :status => "unknown",
            :message => "Simulation result not found",
        )
    end
end


function _legacy_plot_title(filename::AbstractString, index::Integer)::String
    if occursin("plot_main", filename)
        return "Main Results"
    elseif occursin("plot_gc", filename)
        return "Gas Channel Profiles"
    end
    return "Plot $(index)"
end


function _get_plot_field(entry, key::AbstractString, fallback::String)::String
    if entry isa AbstractDict
        if haskey(entry, key)
            return String(entry[key])
        elseif haskey(entry, Symbol(key))
            return String(entry[Symbol(key)])
        end
    end
    return fallback
end


"""Normalize stored plot metadata so the API always returns `{url, title, group, key}` objects."""
function _normalize_plot_entries(result_id::String, entries)::Vector{Dict{String, String}}
    plots = Dict{String, String}[]

    for (index, entry) in enumerate(entries)
        if entry isa AbstractString
            filename = String(entry)
            push!(plots, Dict(
                "url"   => "/api/plots/$(result_id)/$(filename)",
                "title" => _legacy_plot_title(filename, index),
                "group" => "Plots",
                "key"   => "plot_$(index)",
            ))
            continue
        end

        # When the spec already carries a fully-qualified URL (e.g. the live
        # Bonito plot server), use it directly. Otherwise fall back to the
        # legacy filename-based `/api/plots/...` route.
        explicit_url = _get_plot_field(entry, "url", "")
        url = if !isempty(explicit_url)
            explicit_url
        else
            filename = _get_plot_field(entry, "filename", "")
            isempty(filename) ? "" : "/api/plots/$(result_id)/$(filename)"
        end
        isempty(url) && continue

        push!(plots, Dict(
            "url"   => url,
            "title" => _get_plot_field(entry, "title", "Plot $(index)"),
            "group" => _get_plot_field(entry, "group", "Plots"),
            "key"   => _get_plot_field(entry, "key", "plot_$(index)"),
        ))
    end

    return plots
end

"""
Get detailed results data with actual simulation data.

Returns full dataset including time series, KPIs, etc.
"""
function get_detailed_results(result_id::String)::Dict
    if !haskey(SIMULATION_RESULTS, result_id)
        error("Result not found: $result_id")
    end

    result = SIMULATION_RESULTS[result_id]
    sim_type = result[:type]

    # Use actual simulation output if available, otherwise fall back to synthetic data
    if haskey(result, :output) && result[:output] !== nothing
        try
            data = extract_simulation_data(result[:output], sim_type)
        catch e
            @warn "Could not extract simulation data, using synthetic data: $e"
            data = generate_simulation_data(sim_type)
        end
    else
        data = generate_simulation_data(sim_type)
    end

    # Build normalized metadata for the generated interactive plot files.
    plot_entries = get(result, :plot_files, Any[])
    plot_urls = _normalize_plot_entries(result_id, plot_entries)

    return Dict(
        :simulation_type => sim_type,
        :status => result[:status],
        :start_time => haskey(result, :start_time) ? Dates.format(result[:start_time], "yyyy-mm-dd HH:MM:SS") : "N/A",
        :elapsed_time => get(result, :elapsed_time, 0.0),
        :error => get(result, :error, ""),
        :input_params => get(result, :input_params, nothing),
        :data => data,
        :plots => plot_urls,
    )
end

"""
Extract a scalar time series from solver states using the provided accessor function.
The accessor receives one CellState1D and returns a Float64.
Values are averaged across all GC nodes at each timestep.
"""
function _extract_state_series(solver::Any, accessor::Function)::Vector{Float64}
    result = Float64[]
    for state in solver.states
        vals = Float64[]
        for node in state.nodes
            try
                push!(vals, Float64(accessor(node)))
            catch
            end
        end
        push!(result, isempty(vals) ? NaN : mean(vals))
    end
    return result
end


"""
Extract actual simulation data from SimulationOutputs object.
"""
function extract_simulation_data(output::Any, sim_type::String)::Dict
    # output is an AlphaPEM object, we need output.outputs to get SimulationOutputs
    # Then access: output.outputs.solver and output.outputs.derived

    try
        # Check if output has .outputs field
        if !hasfield(typeof(output), :outputs) || output.outputs === nothing
            error("No simulation outputs available in AlphaPEM object")
        end

        sim_outputs = output.outputs
        solver = sim_outputs.solver
        derived = sim_outputs.derived

        time_vec = collect(solver.t)  # Convert to Vector{Float64}
        voltage_vec = collect(derived.Ucell)  # Voltage in V

        # For data stored per GC node, calculate mean across all GC nodes at each timestep
        current_vec = mean_across(derived.i_fc)  # Current in A/m²
        O2_conc = mean_across(derived.C_O2_Pt)  # O2 concentration at catalyst
        v_anode = mean_across(derived.v_a)  # Anode velocity
        v_cathode = mean_across(derived.v_c)  # Cathode velocity

        # Handle empty vectors
        if isempty(current_vec)
            current_vec = zeros(length(time_vec))
        end

        if sim_type == "step"
            # For step simulation, extract complete time series data
            power_vec = voltage_vec .* current_vec

            # Extract state-based variables (mean across GC nodes) to match :synthetic display
            T_vec     = _extract_state_series(solver, node -> node.ccl.T)
            C_v_vec   = _extract_state_series(solver, node -> node.agc.C_v)
            lambda_vec = _extract_state_series(solver, node -> node.mem.lambda)
            C_H2_vec  = _extract_state_series(solver, node -> node.agc.C_H2)
            s_vec     = _extract_state_series(solver, node -> begin
                agdl = node.agdl
                isempty(agdl) ? node.agc.s : agdl[1].s
            end)
            C_O2_vec  = _extract_state_series(solver, node -> begin
                cgdl = node.cgdl
                isempty(cgdl) ? node.cgc.C_O2 : cgdl[1].C_O2
            end)

            result = Dict(
                :time => time_vec,
                :voltage => voltage_vec,
                :current => current_vec,
                :power => power_vec,
                :Pa_in => derived.Pa_in,
                :Pc_in => derived.Pc_in,
                :O2_concentration => O2_conc,
                :anode_velocity => v_anode,
                :cathode_velocity => v_cathode,
                # State-based variables matching :synthetic panels
                :T => T_vec,
                :C_v => C_v_vec,
                :s => s_vec,
                :lambda => lambda_vec,
                :C_H2 => C_H2_vec,
                :C_O2 => C_O2_vec,
                :data_points => length(time_vec),
            )

            return result

        elseif sim_type == "polarization"
            # For polarization, voltage and current form the V-I curve
            power_vec = voltage_vec .* current_vec

            result = Dict(
                :current => current_vec,
                :voltage => voltage_vec,
                :power => power_vec,
                :Pa_in => derived.Pa_in,
                :Pc_in => derived.Pc_in,
                :O2_concentration => O2_conc,
                :data_points => length(current_vec),
            )

            return result

        elseif sim_type == "eis"
            # Extract frequency/impedance points via the unified post-processing
            f_results = make_Fourier_transformation(sim_outputs, output.current_density, output.cfg)

            # Map points for structured JSON response
            freqs = Float64[]
            Z_real = Float64[]
            Z_imag = Float64[]
            abs_Z = Float64[]
            phi_deg = Float64[]

            for res in f_results
                p = _eis_point(output.current_density, res)
                if isfinite(p[3])
                    # p = (Z_real, minus_Z_imag, f, abs_Z, phi)
                    push!(freqs,   p[3])
                    push!(Z_real,  p[1])
                    push!(Z_imag, -p[2])
                    push!(abs_Z,   p[4])
                    push!(phi_deg, p[5])
                end
            end

            result = Dict(
                :voltage             => voltage_vec,
                :current             => current_vec,
                :Pa_in               => derived.Pa_in,
                :Pc_in               => derived.Pc_in,
                :frequency           => freqs,
                :impedance_real      => Z_real,
                :impedance_imag      => Z_imag,
                :impedance_amplitude => abs_Z,
                :phase               => phi_deg,
                :data_points         => length(freqs)
            )

            return result
        else
            return Dict(:data_points => 0)
        end
    catch e
        @error "Error extracting simulation data: $e"
        rethrow(e)
    end
end

"""
Generate synthetic simulation data for visualization.
"""
function generate_simulation_data(sim_type::String)::Dict
    if sim_type == "step"
        # Generate step current response
        time = collect(0:0.5:50)
        voltage = 0.85 .- 0.001 .* time
        current = [5000; fill(15000, length(time)-1)]
        power = voltage .* current
        efficiency = 80 .- 0.1 .* time

        return Dict(
            :time => time,
            :voltage => voltage,
            :current => current,
            :power => power,
            :efficiency => efficiency,
            :data_points => length(time),
        )

    elseif sim_type == "polarization"
        # Generate polarization curve
        current = collect(0:1000:20000)
        voltage = 0.85 .- 0.015 .* (current ./ 1000)
        power = voltage .* current
        efficiency = 80 .- 2 .* (current ./ 1000)

        return Dict(
            :current => current,
            :voltage => voltage,
            :power => power,
            :efficiency => efficiency,
            :data_points => length(current),
        )

    elseif sim_type == "eis"
        # Generate EIS impedance data
        frequency = [10^(i*0.2) for i in 0:19]
        impedance_real = 0.1 .+ 0.005 .* collect(1:20)
        impedance_imag = 0.05 .* sin.(collect(1:20) .* 0.3)

        return Dict(
            :frequency => frequency,
            :impedance_real => impedance_real,
            :impedance_imag => impedance_imag,
            :data_points => length(frequency),
        )
    else
        return Dict(:data_points => 0)
    end
end



"""
Calculate mean of a vector of vectors (across nodes for each timestep).
"""
function mean_across(vec_of_vecs::Vector{Vector{T}}) where T
    if isempty(vec_of_vecs) || isempty(vec_of_vecs[1])
        return Float64[]
    end
    n_nodes = length(vec_of_vecs)
    n_timesteps = length(vec_of_vecs[1])
    return [mean([vec_of_vecs[node][t] for node in 1:n_nodes]) for t in 1:n_timesteps]
end

"""
Export simulation results in various formats.
Returns a tuple (data, suggested_filename).
"""
function export_results(result_id::String, format::String="json"; index=nothing)::Any
    if !haskey(SIMULATION_RESULTS, result_id)
        error("Result not found: $result_id")
    end

    format_lower = lowercase(format)
    
    data = if format_lower == "json"
        JSON.json(SIMULATION_RESULTS[result_id])
    elseif format_lower == "csv" || format_lower == "xlsx" || format_lower == "excel" || format_lower == "xls"
        export_results_xlsx(result_id)
    elseif format_lower == "pdf"
        export_results_pdf(result_id, index)
    else
        error("Unsupported export format: $format")
    end
    
    filename = get_export_filename(result_id, format_lower, index=index)
    
    return (data=data, filename=filename)
end

"""
Get plot data for a specific key.
Returns a tuple (headers, data_matrix).
"""
function _get_plot_data(simu, key, config)
    outputs = simu.outputs
    cd = simu.current_density
    nb_gc = config.numerical_parameters.nb_gc
    
    # ── Temporal plots ────────────────────────────────────────────────────────
    if key in ["ifc_1d_temporal", "ifc_1d_temporal_cali"]
        t = time_history(outputs)
        headers = ["Time (s)", "i_fc_cell (A/cm²)"]
        for i in 1:nb_gc
            push!(headers, "i_fc_$(i) (A/cm²)")
        end
        data = hcat(t, current(cd, t) ./ 1e4)
        for i in 1:nb_gc
            data = hcat(data, extract_derived_gc_series(outputs, i, x -> x.i_fc) ./ 1e4)
        end
        return headers, data
        
    elseif key in ["Ucell_temporal", "Ucell_temporal_cali"]
        t = time_history(outputs)
        headers = ["Time (s)", "Cell Voltage (V)"]
        data = hcat(t, extract_derived_series(outputs, x -> x.Ucell))
        return headers, data
        
    elseif key in ["T_1D_temporal", "T_1D_temporal_cali"]
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        # Average T for each layer
        layers = [
            (node -> node.agc.T, "T_agc (K)"),
            (node -> node.agdl[middle_gdl_index(config)].T, "T_agdl (K)"),
            (node -> node.ampl[middle_mpl_index(config)].T, "T_ampl (K)"),
            (node -> node.acl.T, "T_acl (K)"),
            (node -> node.mem.T, "T_mem (K)"),
            (node -> node.ccl.T, "T_ccl (K)"),
            (node -> node.cmpl[middle_mpl_index(config)].T, "T_cmpl (K)"),
            (node -> node.cgdl[middle_gdl_index(config)].T, "T_cgdl (K)"),
            (node -> node.cgc.T, "T_cgc (K)")
        ]
        for (acc, name) in layers
            push!(headers, name)
            data = hcat(data, extract_mid_mea_series(outputs, config, acc))
        end
        return headers, data

    elseif key == "Cv_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        layers = [
            (node -> node.agc.C_v, "Cv_agc (mol/m³)"),
            (node -> node.agdl[middle_gdl_index(config)].C_v, "Cv_agdl (mol/m³)"),
            (node -> node.ampl[middle_mpl_index(config)].C_v, "Cv_ampl (mol/m³)"),
            (node -> node.acl.C_v, "Cv_acl (mol/m³)"),
            (node -> node.ccl.C_v, "Cv_ccl (mol/m³)"),
            (node -> node.cmpl[middle_mpl_index(config)].C_v, "Cv_cmpl (mol/m³)"),
            (node -> node.cgdl[middle_gdl_index(config)].C_v, "Cv_cgdl (mol/m³)"),
            (node -> node.cgc.C_v, "Cv_cgc (mol/m³)")
        ]
        for (acc, name) in layers
            push!(headers, name)
            data = hcat(data, extract_mid_mea_series(outputs, config, acc))
        end
        push!(headers, "Cv_sat_ccl (mol/m³)")
        T_ccl = extract_mid_mea_series(outputs, config, mea -> mea.ccl.T)
        data = hcat(data, [C_v_sat(T) for T in T_ccl])
        return headers, data

    elseif key == "s_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        layers = [
            (node -> node.agc.s, "s_agc (-)"),
            (node -> node.agdl[middle_gdl_index(config)].s, "s_agdl (-)"),
            (node -> node.ampl[middle_mpl_index(config)].s, "s_ampl (-)"),
            (node -> node.acl.s, "s_acl (-)"),
            (node -> node.ccl.s, "s_ccl (-)"),
            (node -> node.cmpl[middle_mpl_index(config)].s, "s_cmpl (-)"),
            (node -> node.cgdl[middle_gdl_index(config)].s, "s_cgdl (-)"),
            (node -> node.cgc.s, "s_cgc (-)")
        ]
        for (acc, name) in layers
            push!(headers, name)
            data = hcat(data, extract_mid_mea_series(outputs, config, acc))
        end
        return headers, data

    elseif key == "lambda_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)", "lambda_acl (-)", "lambda_mem (-)", "lambda_ccl (-)"]
        data = hcat(t, 
            extract_mid_mea_series(outputs, config, node -> node.acl.lambda),
            extract_mid_mea_series(outputs, config, node -> node.mem.lambda),
            extract_mid_mea_series(outputs, config, node -> node.ccl.lambda)
        )
        return headers, data

    elseif key == "CH2_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        for i in 1:nb_gc
            push!(headers, "CH2_$(i) (mol/m³)")
            data = hcat(data, extract_mea_series(outputs, i, mea -> mea.agc.C_H2))
        end
        return headers, data

    elseif key == "CO2_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        for i in 1:nb_gc
            push!(headers, "CO2_$(i) (mol/m³)")
            data = hcat(data, extract_mea_series(outputs, i, mea -> mea.cgc.C_O2))
        end
        return headers, data

    elseif key == "P_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        for i in 1:nb_gc
            push!(headers, "Pa_$(i) (Pa)")
            data = hcat(data, extract_mea_series(outputs, i, mea -> (mea.agc.C_H2 + mea.agc.C_v + mea.agc.C_N2) * R * mea.agc.T))
            push!(headers, "Pc_$(i) (Pa)")
            data = hcat(data, extract_mea_series(outputs, i, mea -> (mea.cgc.C_O2 + mea.cgc.C_v + mea.cgc.C_N2) * R * mea.cgc.T))
        end
        return headers, data

    elseif key == "CN2_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        for i in 1:nb_gc
            push!(headers, "CN2_$(i) (mol/m³)")
            data = hcat(data, extract_mea_series(outputs, i, mea -> mea.agc.C_N2))
            push!(headers, "CN2_cgc_$(i) (mol/m³)")
            data = hcat(data, extract_mea_series(outputs, i, mea -> mea.cgc.C_N2))
        end
        return headers, data

    elseif key == "Phi_a_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        layers = [
            (node -> node.agc.C_v / C_v_sat(node.agc.T), "Phi_a,agc (-)"),
            (node -> node.agdl[middle_gdl_index(config)].C_v / C_v_sat(node.agdl[middle_gdl_index(config)].T), "Phi_a,agdl (-)"),
            (node -> node.ampl[middle_mpl_index(config)].C_v / C_v_sat(node.ampl[middle_mpl_index(config)].T), "Phi_a,ampl (-)"),
            (node -> node.acl.C_v / C_v_sat(node.acl.T), "Phi_a,acl (-)")
        ]
        for (acc, name) in layers
            push!(headers, name)
            data = hcat(data, extract_mid_mea_series(outputs, config, acc))
        end
        return headers, data

    elseif key == "Phi_c_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        layers = [
            (node -> node.ccl.C_v / C_v_sat(node.ccl.T), "Phi_c,ccl (-)"),
            (node -> node.cmpl[middle_mpl_index(config)].C_v / C_v_sat(node.cmpl[middle_mpl_index(config)].T), "Phi_c,cmpl (-)"),
            (node -> node.cgdl[middle_gdl_index(config)].C_v / C_v_sat(node.cgdl[middle_gdl_index(config)].T), "Phi_c,cgdl (-)"),
            (node -> node.cgc.C_v / C_v_sat(node.cgc.T), "Phi_c,cgc (-)")
        ]
        for (acc, name) in layers
            push!(headers, name)
            data = hcat(data, extract_mid_mea_series(outputs, config, acc))
        end
        return headers, data

    elseif key == "v_1D_temporal"
        t = time_history(outputs)
        headers = ["Time (s)"]
        data = reshape(t, :, 1)
        for i in 1:nb_gc
            push!(headers, "v_a_$(i) (m/s)")
            data = hcat(data, extract_derived_gc_series(outputs, i, x -> x.v_a))
            push!(headers, "v_c_$(i) (m/s)")
            data = hcat(data, extract_derived_gc_series(outputs, i, x -> x.v_c))
        end
        return headers, data

    elseif key == "Re_1D_temporal"
        Re_a_full, Re_c_full = calculate_reynolds_numbers(outputs, simu.fuel_cell)
        t = time_history(outputs)

        headers = ["Time (s)"]
        data = reshape(t, :, 1)

        for i in 1:nb_gc
            push!(headers, "Re_a_$(i) (-)")
            data = hcat(data, Re_a_full[i])
            push!(headers, "Re_c_$(i) (-)")
            data = hcat(data, Re_c_full[i])
        end
        return headers, data

    # ── Performance curves ───────────────────────────────────────────────────
    elseif key == "hysteresis_curve"
        # Hysteresis curve shows all points (not just stationary points)
        t_all = time_history(outputs)
        Ucell = extract_derived_series(outputs, x -> x.Ucell)
        i_fc = [current(cd, tt) for tt in t_all] ./ 1e4
        headers = ["Current Density (A/cm²)", "Cell Voltage (V)"]
        return headers, hcat(i_fc, Ucell)

    elseif key in ["polarization_curve", "polarization_curve_cali",
                   "power_density_curve", "power_density_curve_cali",
                   "efficiency_curve", "efficiency_curve_cali"]

        t_all = time_history(outputs)
        if cd isa Union{PolarizationCurrent, PolarizationCalibrationCurrent}
            # Only export stationary points for polarization curves as requested
            indices = polarisation_sampling_indices(outputs, cd)
            t = t_all[indices]
            Ucell = extract_derived_series(outputs, x -> x.Ucell)[indices]
        else
            t = t_all
            Ucell = extract_derived_series(outputs, x -> x.Ucell)
        end
        i_fc = [current(cd, tt) for tt in t] ./ 1e4

        if key in ["polarization_curve", "polarization_curve_cali"]
            headers = ["Current Density (A/cm²)", "Cell Voltage (V)"]
            return headers, hcat(i_fc, Ucell)
        elseif key in ["power_density_curve", "power_density_curve_cali"]
            headers = ["Current Density (A/cm²)", "Power Density (W/cm²)"]
            return headers, hcat(i_fc, Ucell .* i_fc)
        else # efficiency_curve
            headers = ["Current Density (A/cm²)", "Efficiency (%)"]
            return headers, hcat(i_fc, (Ucell ./ 1.25) .* 100)
        end

    # ── EIS ──────────────────────────────────────────────────────────────────
    elseif key in ["Nyquist_plot", "Bode_amplitude_curve", "Bode_angle_curve"]
        if cd isa EISCurrent
            # Select headers and data extraction based on the plot type
            # p from _eis_point is (Z_real, minus_Z_imag, f, abs_Z, phi)
            
            # Use our unified post-processing to get all analyzed points.
            # (Note: make_Fourier_transformation is called for each sheet, but it's 
            # consistent with how other plots are handled here).
            fourier_results = make_Fourier_transformation(outputs, cd, config)
            rows = []
            
            for res in fourier_results
                p = _eis_point(cd, res)
                if isfinite(p[3])
                    if key == "Nyquist_plot"
                        # Only Z_real and -Z_imag as displayed on screen
                        push!(rows, [p[1], p[2]])
                    elseif key == "Bode_amplitude_curve"
                        # Frequency and Amplitude
                        push!(rows, [p[3], p[4]])
                    else # Bode_angle_curve
                        # Frequency and Phase
                        push!(rows, [p[3], p[5]])
                    end
                end
            end

            if key == "Nyquist_plot"
                headers = ["Z_real (mΩ·cm²)", "-Z_imag (mΩ·cm²)"]
            elseif key == "Bode_amplitude_curve"
                headers = ["Frequency (Hz)", "|Z| (mΩ·cm²)"]
            else
                headers = ["Frequency (Hz)", "Phase (°)"]
            end
            
            if isempty(rows)
                return headers, zeros(0, 2)
            end
            return headers, vcat([r' for r in rows]...)
        end

    # ── Final profiles ───────────────────────────────────────────────────────
    elseif key == "ifc_GC_final"
        headers = ["Node Index", "Final Current Density (A/cm²)"]
        x_gc = collect(1:nb_gc)
        i_fc_hist = derived_outputs(outputs).i_fc
        i_fc_final = [i_fc_hist[k][end] for k in 1:nb_gc] ./ 1e4
        return headers, hcat(x_gc, i_fc_final)

    elseif key == "CO2_Pt_GC_final"
        headers = ["Node Index", "Final Pt O2 Concentration (mol/m³)"]
        x_gc = collect(1:nb_gc)
        C_O2_Pt_hist = derived_outputs(outputs).C_O2_Pt
        C_O2_Pt_final = [C_O2_Pt_hist[k][end] for k in 1:nb_gc]
        return headers, hcat(x_gc, C_O2_Pt_final)

    elseif key == "lambda_mem_GC_final"
        headers = ["Node Index", "Final Membrane Water Content (-)"]
        x_gc = collect(1:nb_gc)
        lambda_hist = [extract_mea_series(outputs, k, node -> node.mem.lambda) for k in 1:nb_gc]
        lambda_final = [h[end] for h in lambda_hist]
        return headers, hcat(x_gc, lambda_final)
        
    elseif key == "T_pseudo_2D_final"
        headers = ["Node Index", "Final Temperature (K)"]
        x_gc = collect(1:nb_gc)
        T_hist = [extract_mea_series(outputs, k, node -> node.mem.T) for k in 1:nb_gc]
        T_final = [h[end] for h in T_hist]
        return headers, hcat(x_gc, T_final)
    end
    
    return [], zeros(0,0)
end


"""
Export a specific simulation plot as PDF.
"""
function export_results_pdf(result_id::String, index::Any)::Vector{UInt8}
    if !haskey(SIMULATION_RESULTS, result_id)
        error("Result not found: $result_id")
    end
    result = SIMULATION_RESULTS[result_id]
    output = get(result, :output, nothing)
    if output === nothing
        error("Simulation output not available for $result_id")
    end
    
    config = result[:config]
    nb_gc = config.numerical_parameters.nb_gc
    
    # Temporarily switch to :multiple / :postrun for rendering.
    orig_display = config.type_display
    orig_timing  = config.display_timing
    config.type_display   = :multiple
    config.display_timing = :postrun
    
    try
        # Activate CairoMakie for PDF generation
        CairoMakie.activate!()
        
        # Prepare figures (identical to how they are prepared for the web)
        fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(config, nb_gc; backend=:cairo)
        
        # Populate with data
        display!(output, ax1, ax2, ax3)
        
        # Collect all figures
        figures = Any[]
        for fig in (fig1, fig2, fig3)
            fig === nothing && continue
            if fig isa AbstractVector
                append!(figures, fig)
            else
                push!(figures, fig)
            end
        end
        
        if isempty(figures)
            error("No figures generated for export (check if display is enabled in config)")
        end
        
        # Determine which figure to export
        idx = 1
        if index !== nothing
            try
                # index is passed as a string from the URL
                idx = parse(Int, string(index)) + 1 # JavaScript index is 0-based
            catch
                idx = 1
            end
        end
        
        if idx < 1 || idx > length(figures)
            idx = 1
        end
        
        target_fig = figures[idx]
        
        # Export to PDF via temporary file
        data = mktempdir() do tmpdir
            path = joinpath(tmpdir, "plot.pdf")
            save(path, target_fig)
            return read(path)
        end
        
        return data
    finally
        # Restore original config
        config.type_display   = orig_display
        config.display_timing = orig_timing
        
        # Restore WGLMakie for the web interface
        WGLMakie.activate!()
    end
end

"""
Export simulation results as XLSX with multiple sheets.
"""
function export_results_xlsx(result_id::String)::Vector{UInt8}
    if !haskey(SIMULATION_RESULTS, result_id)
        error("Result not found: $result_id")
    end

    result = SIMULATION_RESULTS[result_id]
    output = result[:output]
    config = result[:config]
    nb_gc = config.numerical_parameters.nb_gc
    
    # Save original timing mode
    orig_timing = config.display_timing
    
    # Use full history for XLSX export
    config.display_timing = :full_history
    
    try
        # Get the same plot list as shown in the UI
        plot_specs = _web_plot_specs(config, nb_gc)
        
        # Create a temporary file to save the XLSX
        xl_data = mktempdir() do tmpdir
            xl_path = joinpath(tmpdir, "results.xlsx")
            
            XLSX.openxlsx(xl_path, mode="w") do xf
                # The 'w' mode creates one sheet named "Sheet1" by default
                sheet = xf[1]
                first_sheet = true
                
                for spec in plot_specs
                    key = spec["key"]
                    title = spec["title"]
                    
                    headers, data = _get_plot_data(output, key, config)
                    
                    if !isempty(headers)
                        # Sanitize title for sheet name (max 31 chars, no special chars)
                        sheet_name = replace(title, r"[\\/*?:\[\]]" => "_")
                        if length(sheet_name) > 31
                            sheet_name = sheet_name[1:31]
                        end
                        
                        if first_sheet
                            XLSX.rename!(sheet, sheet_name)
                            first_sheet = false
                        else
                            sheet = XLSX.addsheet!(xf, sheet_name)
                        end
                        
                        # Write headers in the first row
                        for (c, h) in enumerate(headers)
                            sheet[1, c] = h
                        end
                        
                        # Write data starting from second row
                        if !isempty(data)
                            # data is a Matrix{Float64}
                            # XLSX can write the whole matrix at once
                            sheet[2, 1] = data
                        end
                    end
                end
                
                # If no sheets were written (unlikely), keep at least one
                if first_sheet
                    XLSX.rename!(sheet, "No Data")
                end
            end
            return read(xl_path)
        end
        
        return xl_data
    finally
        # Restore original timing mode
        config.display_timing = orig_timing
    end
end

"""
Get suggested filename for export.
"""
function get_export_filename(result_id::String, format::String; index=nothing)::String
    if !haskey(SIMULATION_RESULTS, result_id)
        return "alphapem_export.$(format)"
    end
    
    result = SIMULATION_RESULTS[result_id]
    config = result[:config]
    nb_gc = config.numerical_parameters.nb_gc
    
    # result_id already follows format AlphaPEM_YYYY_MM_DD_type_vX
    base_name = result_id
    
    if format == "pdf"
        plot_specs = _web_plot_specs(config, nb_gc)
        idx = 1
        if index !== nothing
            try
                idx = parse(Int, string(index)) + 1
            catch
                idx = 1
            end
        end
        
        if idx >= 1 && idx <= length(plot_specs)
            plot_title = plot_specs[idx]["title"]
            # Sanitize plot title for filename (lowercase, underscores)
            plot_slug = replace(lowercase(plot_title), r"[^a-z0-9]+" => "_")
            plot_slug = strip(plot_slug, '_')
            
            # Insert plot_slug before _vX
            m = match(r"^(.*)(_v\d+)$", base_name)
            if m !== nothing
                return "$(m.captures[1])_$(plot_slug)$(m.captures[2]).pdf"
            else
                return "$(base_name)_$(plot_slug).pdf"
            end
        end
        return "$(base_name).pdf"
    elseif format in ["csv", "xlsx", "excel", "xls"]
        # Both "csv" button and XLSX button should now provide XLSX with sheets
        return "$(base_name).xlsx"
    else
        return "$(base_name).$(format)"
    end
end

# ------ UTILITY FUNCTIONS ------

"""
Generate unique result ID following the format: AlphaPEM_YYYY_MM_DD_type_vX
"""
function generate_result_id(sim_type::String)::String
    date_str = Dates.format(now(), "yyyy_mm_dd")
    
    # Simple versioning: find existing results for the same day and type
    base_id = "AlphaPEM_$(date_str)_$(sim_type)"
    
    version = 1
    id = "$(base_id)_v$(version)"
    
    while haskey(SIMULATION_RESULTS, id)
        version += 1
        id = "$(base_id)_v$(version)"
    end
    
    return id
end

end  # end module SimulatorBackend

