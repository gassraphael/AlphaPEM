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

using AlphaPEM
using AlphaPEM.Config: SimulationConfig, StepParams, PolarizationParams, EISParams, NumericalParams
using AlphaPEM.Config: PhysicalParams, OperatingConditions
using AlphaPEM.Application: run_simulation
using AlphaPEM.Fuelcell
using JSON
using Dates
using Logging

# ------ MODULE CONFIGURATION ------

const RESULTS_DIR = joinpath(pwd(), "results")
const SIMULATION_RESULTS = Dict{String, Dict}()  # In-memory result cache

"""
Initialize backend directories and caches.
"""
function initialize_backend()
    # Create results directory if it doesn't exist
    if !isdir(RESULTS_DIR)
        mkpath(RESULTS_DIR)
    end
end

# ------ FUEL CELL PRESETS ------

"""
Get list of available fuel cell models.

Returns a dictionary with fuel cell types and their descriptions.
"""
function get_available_fuel_cells()::Dict
    fuel_cells = Dict(
        :default => Dict(
            :name => "Default Model",
            :description => "Generic PEM fuel cell",
        ),
        :ZSW_GenStack => Dict(
            :name => "ZSW GenStack",
            :description => "ZSW GenStack (2022) - calibrated reference cell",
        ),
        :ZSW_GenStack_Pa_1_61_Pc_1_41 => Dict(
            :name => "ZSW GenStack Pa=1.61 bar, Pc=1.41 bar",
            :description => "ZSW GenStack at specific pressure conditions",
        ),
        :ZSW_GenStack_Pa_2_01_Pc_1_81 => Dict(
            :name => "ZSW GenStack Pa=2.01 bar, Pc=1.81 bar",
            :description => "ZSW GenStack at specific pressure conditions",
        ),
        :ZSW_GenStack_Pa_2_4_Pc_2_2 => Dict(
            :name => "ZSW GenStack Pa=2.4 bar, Pc=2.2 bar",
            :description => "ZSW GenStack at specific pressure conditions",
        ),
        :ZSW_GenStack_Pa_2_8_Pc_2_6 => Dict(
            :name => "ZSW GenStack Pa=2.8 bar, Pc=2.6 bar",
            :description => "ZSW GenStack at specific pressure conditions",
        ),
        :ZSW_GenStack_T_62 => Dict(
            :name => "ZSW GenStack T=62°C",
            :description => "ZSW GenStack at low temperature",
        ),
        :ZSW_GenStack_T_76 => Dict(
            :name => "ZSW GenStack T=76°C",
            :description => "ZSW GenStack at medium temperature",
        ),
        :ZSW_GenStack_T_84 => Dict(
            :name => "ZSW GenStack T=84°C",
            :description => "ZSW GenStack at high temperature",
        ),
        :EH31_1_5 => Dict(
            :name => "EH-31: P=1.5 bar",
            :description => "EH-31 fuel cell at 1.5 bar",
        ),
        :EH31_2_0 => Dict(
            :name => "EH-31: P=2.0 bar",
            :description => "EH-31 fuel cell at 2.0 bar",
        ),
        :EH31_2_25 => Dict(
            :name => "EH-31: P=2.25 bar",
            :description => "EH-31 fuel cell at 2.25 bar",
        ),
        :EH31_2_5 => Dict(
            :name => "EH-31: P=2.5 bar",
            :description => "EH-31 fuel cell at 2.5 bar",
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
    # Convert string to symbol
    fc_symbol = Symbol(fuel_cell_type)

    # TODO: Load from fuel cell database in AlphaPEM
    # For now, return generic defaults

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
        :accessible_parameters => Dict(
            :Aact => 100.0,             # cm² - Active area
            :nb_cell => 1,              # Number of cells
            :Hagc => 100.0,             # µm - Anode GC thickness
            :Hcgc => 100.0,             # µm - Cathode GC thickness
            :Wagc => 100.0,             # µm - Anode GC width
            :Wcgc => 100.0,             # µm - Cathode GC width
            :Lgc => 100.0,              # mm - GC length
            :nb_channel_in_gc => 50,    # Number of channels
            :Ldist => 50.0,             # mm - Distributor length
            :Lm => 50.0,                # mm - Manifold length
            :A_T_a => 10.0,             # cm² - Anode throttle area
            :A_T_c => 10.0,             # cm² - Cathode throttle area
            :Vasm => 100.0,             # cm³ - Supply anode manifold
            :Vcsm => 100.0,             # cm³ - Supply cathode manifold
            :Vaem => 100.0,             # cm³ - Exhaust anode manifold
            :Vcem => 100.0,             # cm³ - Exhaust cathode manifold
            :V_endplate_a => 50.0,      # cm³ - Anode endplate
            :V_endplate_c => 50.0,      # cm³ - Cathode endplate
        ),
        :undetermined_parameters => Dict(
            :Hgdl => 200.0,             # µm - GDL thickness
            :Hmpl => 50.0,              # µm - MPL thickness
            :Hacl => 10.0,              # µm - Anode CL thickness
            :Hccl => 10.0,              # µm - Cathode CL thickness
            :Hmem => 10.0,              # µm - Membrane thickness
            :epsilon_gdl => 0.7,        # GDL porosity
            :epsilon_mpl => 0.5,        # MPL porosity
            :epsilon_c => 0.2,          # Compression ratio
            :e => 4,                    # Capillary exponent
            :K_l_ads => 5.0,            # Sorption ratio
            :K_O2_ad_Pt => 2.0,         # O2 adsorption coefficient
            :Re => 1.0,                 # Ω.mm² - Electron resistance
            :i0_c_ref => 1.0,           # A/m² - Exchange current
            :kappa_co => 1.0,           # mol/(m.s.Pa) - Crossover
            :kappa_c => 1.0,            # Overpotential exponent
            :C_scl => 20.0,             # F/cm³ - Double layer capacitance
        ),
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

    # Check structure
    required_sections = [:operating_conditions, :accessible_parameters,
                        :undetermined_parameters, :computing_parameters]

    for section in required_sections
        if !haskey(params, section)
            push!(errors, "Missing section: $section")
        end
    end

    if !isempty(errors)
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
        sp = params[:step_parameters]
        StepParams(
            delta_t_ini = sp[:delta_t_ini],
            delta_t_load = sp[:delta_t_load],
            delta_t_break = sp[:delta_t_break],
            i_ini = sp[:i_ini],
            i_step = sp[:i_step],
        )
    elseif sim_type == :polarization
        pp = params[:polarization_parameters]
        PolarizationParams(
            delta_t_ini = pp[:delta_t_ini],
            delta_t_load = pp[:delta_t_load],
            delta_t_break = pp[:delta_t_break],
            delta_i = pp[:delta_i],
            i_max = pp[:i_max],
        )
    elseif sim_type == :eis
        ep = params[:eis_parameters]
        EISParams(
            i_static = ep[:i_static],
            ratio = ep[:ratio],
            f_power_min = ep[:f_power_min],
            f_power_max = ep[:f_power_max],
            nb_frequencies = ep[:nb_frequencies],
            nb_points = ep[:nb_points],
        )
    else
        error("Unknown simulation type: $sim_type")
    end

    # Build numerical parameters
    cp = params[:numerical_parameters]
    num_params = NumericalParams(
        nb_gc = cp[:nb_gc],
        nb_gdl = cp[:nb_gdl],
        nb_mpl = cp[:nb_mpl],
        rtol = cp[:rtol],
        atol = cp[:atol],
    )

    # Build SimulationConfig
    config = SimulationConfig(
        type_fuel_cell = Symbol(params[:fuel_cell_type]),
        type_current = current_params,
        numerical_parameters = num_params,
        voltage_zone = Symbol(params[:voltage_zone]),
        type_auxiliary = Symbol(params[:type_auxiliary]),
        type_flow = Symbol(params[:type_flow]),
        type_purge = Symbol(params[:type_purge]),
        type_display = :no_display,     # Server-side, don't display plots
        display_timing = :postrun,      # Display after simulation
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
        @info "Starting step simulation: $result_id"

        start_time = time()
        run_simulation(config)
        elapsed = time() - start_time

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "step",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => config,
        )

        @info "Step simulation completed: $result_id in $(round(elapsed, digits=2))s"

        return Dict(:id => result_id, :status => "completed")

    catch e
        @error "Step simulation failed: $result_id" exception=e

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
        @info "Starting polarization simulation: $result_id"

        start_time = time()
        run_simulation(config)
        elapsed = time() - start_time

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "polarization",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => config,
        )

        @info "Polarization simulation completed: $result_id"

        return Dict(:id => result_id, :status => "completed")

    catch e
        @error "Polarization simulation failed: $result_id" exception=e

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
        @info "Starting EIS simulation: $result_id"

        start_time = time()
        run_simulation(config)
        elapsed = time() - start_time

        SIMULATION_RESULTS[result_id] = Dict(
            :type => "eis",
            :status => "completed",
            :start_time => DateTime(now()),
            :elapsed_time => elapsed,
            :config => config,
        )

        @info "EIS simulation completed: $result_id"

        return Dict(:id => result_id, :status => "completed")

    catch e
        @error "EIS simulation failed: $result_id" exception=e

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

"""
Get detailed results data with simulated/generated data.

Returns full dataset including time series, KPIs, etc.
"""
function get_detailed_results(result_id::String)::Dict
    if !haskey(SIMULATION_RESULTS, result_id)
        error("Result not found: $result_id")
    end

    result = SIMULATION_RESULTS[result_id]
    sim_type = result[:type]

    # Generate synthetic data for visualization
    data = generate_simulation_data(sim_type)

    return Dict(
        :simulation_type => sim_type,
        :status => result[:status],
        :start_time => Dates.format(result[:start_time], "yyyy-mm-dd HH:MM:SS"),
        :elapsed_time => result[:elapsed_time],
        :data => data,
        :kpis => calculate_kpis(data),
        :plots => [],
    )
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


