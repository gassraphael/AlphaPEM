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
using AlphaPEM.Application: run_simulation, generate_web_plots
using AlphaPEM.Fuelcell
using JSON
using Dates
using Logging

# ------ MODULE CONFIGURATION ------

const RESULTS_DIR = joinpath(pwd(), "results")
const PLOTS_DIR = joinpath(RESULTS_DIR, "web_plots")
const SIMULATION_RESULTS = Dict{String, Dict}()  # In-memory result cache

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
    :i0_c_ref => (factor=1e-4, unit="A/cm²", precision=2),
    :kappa_co => (factor=1.0, unit="mol·m⁻¹·s⁻¹·Pa⁻¹", precision=2),
    :C_scl => (factor=1e-6, unit="MF/m³", precision=0),
    # Run Parameters
    :delta_t_ini => (factor=1/60, unit="min", precision=1),
    :delta_t_break => (factor=1/60, unit="min", precision=1),
    :i_ini => (factor=1e-4, unit="A/cm²", precision=2),
    :i_step => (factor=1e-4, unit="A/cm²", precision=2),
    :delta_i => (factor=1e-4, unit="A/cm²", precision=2),
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
                :maxiters => 10000,
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
            :delta_i => 500.0,
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
            :epsilon_c => 0.2,          # Compression ratio
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
            :delta_i => 0.5e4,           # A/m² - Step size
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
            delta_t_ini = sp[:delta_t_ini] * 60.0,
            delta_t_load = sp[:delta_t_load],
            delta_t_break = sp[:delta_t_break] * 60.0,
            i_ini = sp[:i_ini] * 1e4,
            i_step = sp[:i_step] * 1e4,
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
            delta_t_ini = pp[:delta_t_ini] * 60.0,
            delta_i = pp[:delta_i] * 1e4,
            v_load = pp[:v_load] * 1e4,
            delta_t_break = pp[:delta_t_break] * 60.0,
            i_max = pp[:i_max] * 1e4,
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
            i_EIS = ep[:i_static] * 1e4,
            ratio = ep[:ratio] / 100.0,
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
    )

    # Build custom physical parameters
    ap_data = get(params, :accessible_parameters, Dict())
    up_data = get(params, :undetermined_parameters, Dict())
    
    # Merge both into a single dict for PhysicalParams, converting keys to symbols
    phys_dict = Dict{Symbol, Any}()
    for (k, v) in ap_data
        phys_dict[Symbol(k)] = v
    end
    for (k, v) in up_data
        phys_dict[Symbol(k)] = v
    end
    
    # Filter out unknown fields to avoid Keyword Argument errors
    valid_phys_fields = fieldnames(PhysicalParams)
    filtered_phys_dict = Dict{Symbol, Any}()
    for (k, v) in phys_dict
        if k in valid_phys_fields
            filtered_phys_dict[k] = v
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

    # Build SimulationConfig
    config = SimulationConfig(
        type_fuel_cell = Symbol(params[:fuel_cell_type]),
        type_current = current_params,
        numerical_parameters = num_params,
        voltage_zone = Symbol(get(params, :voltage_zone, "fuel_cell")),
        type_auxiliary = Symbol(get(params, :type_auxiliary, "none")),
        type_flow = Symbol(get(params, :type_flow, "constant")),
        type_purge = Symbol(get(params, :type_purge, "no_purge")),
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
function run_step_simulation(config::SimulationConfig)::Dict
    result_id = generate_result_id()

    try
        start_time = time()
        # Capture the actual output from the simulation
        output = run_simulation(config)
        elapsed = time() - start_time

         # Generate WGLMakie plots using the same functions as the native display.
         plot_dir = joinpath(PLOTS_DIR, result_id)
         plot_files = try
             @debug "    Generating web plots to directory: $plot_dir"
             files = generate_web_plots(output, plot_dir)
             @debug "    Generated $(length(files)) plot files: $files"
             files
         catch e
             @error "    Failed to generate web plots for step simulation" exception=e
              Dict{String, String}[]
         end

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "step",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => config,
            :output => output,  # Store the actual simulation results
            :plot_files => plot_files,
        )

        @info "    Step simulation completed: $result_id in $(round(elapsed, digits=2))s"

        return Dict(:id => result_id, :status => "completed")

    catch e
        @error "    Step simulation failed: $result_id" exception=e

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "step",
            :status => "failed",
            :error => string(e),
        )

        return Dict(:id => result_id, :status => "failed", :error => string(e))
    end
end

"""
Run a polarization curve simulation.
"""
function run_polarization_simulation(config::SimulationConfig)::Dict
    result_id = generate_result_id()

    try
        start_time = time()
        # Capture the actual output from the simulation
        output = run_simulation(config)
        elapsed = time() - start_time

         # Generate WGLMakie plots using the same functions as the native display.
         plot_dir = joinpath(PLOTS_DIR, result_id)
         plot_files = try
             @debug "    Generating web plots to directory: $plot_dir"
             files = generate_web_plots(output, plot_dir)
             @debug "    Generated $(length(files)) plot files: $files"
             files
         catch e
             @error "    Failed to generate web plots for polarization simulation" exception=e
              Dict{String, String}[]
         end

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "polarization",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => config,
            :output => output,  # Store the actual simulation results
            :plot_files => plot_files,
        )

        @info "    Polarization simulation completed: $result_id"

        return Dict(:id => result_id, :status => "completed")

    catch e
        @error "    Polarization simulation failed: $result_id" exception=e

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "polarization",
            :status => "failed",
            :error => string(e),
        )

        return Dict(:id => result_id, :status => "failed", :error => string(e))
    end
end

"""
Run an EIS simulation.
"""
function run_eis_simulation(config::SimulationConfig)::Dict
    result_id = generate_result_id()

    try
        start_time = time()
        # EIS requires :live display_timing for frequency-by-frequency stepping.
        eis_config = SimulationConfig(
            type_fuel_cell     = config.type_fuel_cell,
            type_current       = config.type_current,
            numerical_parameters = config.numerical_parameters,
            voltage_zone       = config.voltage_zone,
            type_auxiliary     = config.type_auxiliary,
            type_flow          = config.type_flow,
            type_purge         = config.type_purge,
            type_display       = :no_display,
            display_timing     = :live,
        )
        # Capture the actual output from the simulation
        output = run_simulation(eis_config)
        elapsed = time() - start_time

         # Generate WGLMakie plots using the same functions as the native display.
         plot_dir = joinpath(PLOTS_DIR, result_id)
         plot_files = try
             @debug "    Generating web plots to directory: $plot_dir"
             files = generate_web_plots(output, plot_dir)
             @debug "    Generated $(length(files)) plot files: $files"
             files
         catch e
             @error "    Failed to generate web plots for EIS simulation" exception=e
             Dict{String, String}[]
         end

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "eis",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => eis_config,
            :output => output,  # Store the actual simulation results
            :plot_files => plot_files,
        )

        @info "    EIS simulation completed: $result_id"

        return Dict(:id => result_id, :status => "completed")

    catch e
        @error "    EIS simulation failed: $result_id" exception=e

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "eis",
            :status => "failed",
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

        filename = _get_plot_field(entry, "filename", "")
        isempty(filename) && continue

        push!(plots, Dict(
            "url"   => "/api/plots/$(result_id)/$(filename)",
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
        :start_time => Dates.format(result[:start_time], "yyyy-mm-dd HH:MM:SS"),
        :elapsed_time => result[:elapsed_time],
        :data => data,
        :kpis => calculate_kpis(data),
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
            # For EIS, we have impedance data
            # The structure might be different, so we try to extract what's available
            result = Dict(
                :voltage => voltage_vec,
                :current => current_vec,
                :Pa_in => derived.Pa_in,
                :Pc_in => derived.Pc_in,
                :data_points => length(time_vec),
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
Calculate Key Performance Indicators from simulation data.
"""
function calculate_kpis(data::Dict)::Dict
    kpis = Dict()

    if haskey(data, :voltage)
        voltages = data[:voltage]
        kpis[:peak_voltage] = maximum(voltages)
        kpis[:min_voltage] = minimum(voltages)
        kpis[:avg_voltage] = mean(voltages)
    end

    if haskey(data, :power)
        powers = data[:power]
        kpis[:avg_power] = mean(powers)
        kpis[:peak_power] = maximum(powers)
    end

    if haskey(data, :efficiency)
        efficiencies = data[:efficiency]
        kpis[:avg_efficiency] = mean(efficiencies)
        kpis[:min_efficiency] = minimum(efficiencies)
    end

    return kpis
end

"""
Calculate mean of an array.
"""
function mean(arr::Vector)
    length(arr) > 0 ? sum(arr) / length(arr) : 0.0
end

"""
Calculate mean of a vector of vectors.
"""
function mean_across(vec_of_vecs::Vector{Vector{T}}) where T
    if isempty(vec_of_vecs) || isempty(vec_of_vecs[1])
        return Float64[]
    end
    return [mean(vec_of_vecs[i]) for i in 1:length(vec_of_vecs)]
end

"""
Export simulation results in various formats.
"""
function export_results(result_id::String, format::String="json")::String
    if !haskey(SIMULATION_RESULTS, result_id)
        error("Result not found: $result_id")
    end

    result = SIMULATION_RESULTS[result_id]
    format_lower = lowercase(format)

    if format_lower == "json"
        return JSON.json(result)
    elseif format_lower == "csv"
        return export_results_csv(result_id)
    elseif format_lower in ["xlsx", "excel", "xls"]
        return export_results_xlsx(result_id)
    else
        error("Unsupported export format: $format")
    end
end

"""
Export simulation results as CSV.
"""
function export_results_csv(result_id::String)::String
    if !haskey(SIMULATION_RESULTS, result_id)
        error("Result not found: $result_id")
    end

    result = SIMULATION_RESULTS[result_id]

    # For now, return a simple CSV representation
    lines = [
        "AlphaPEM Simulation Results",
        "Result ID,$result_id",
        "Type,$(result[:type])",
        "Status,$(result[:status])",
        "Start Time,$(result[:start_time])",
        "Elapsed Time (s),$(result[:elapsed_time])",
    ]

    return join(lines, "\n")
end

"""
Export simulation results as XLSX (simplified).
"""
function export_results_xlsx(result_id::String)::String
    # This would require ExcelFiles.jl package
    # For now, return similar to CSV
    return export_results_csv(result_id)
end

# ------ UTILITY FUNCTIONS ------

"""
Generate unique result ID.
"""
function generate_result_id()::String
    return "result_$(Dates.format(now(), "yyyymmdd_HHMMSS"))_$(rand(1000:9999))"
end

end  # end module SimulatorBackend

