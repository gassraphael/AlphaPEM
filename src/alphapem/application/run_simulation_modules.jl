# -*- coding: utf-8 -*-

"""This module contains some of the required functions for the run_simulation.jl file.
"""

# _____________________________________________________Main modules_____________________________________________________

"""Create required figures and axes according to a `SimulationConfig` object.

# Arguments
- `cfg::SimulationConfig`: Simulation configuration object.
  The function uses `cfg.type_current`, `cfg.type_display`, and `cfg.type_plot`.

# Returns
- `(fig1, ax1, fig2, ax2, fig3, ax3)`: Figure/axes tuple.
"""
function figures_preparation(cfg::SimulationConfig)

    mpl.rcParams["font.family"] = "cmr10"  # "cmr10" for English characters and "DejaVu Serif" for French ones
    mpl.rcParams["axes.formatter.use_mathtext"] = true  # For scientific notation
    mpl.rcParams["lines.linewidth"] = 2.0
    mpl.rcParams["lines.markersize"] = 5.0

    # Enable interactive mode for non-blocking display
    plt.ion()

    fig1, ax1 = nothing, nothing
    fig2, ax2 = nothing, nothing
    fig3, ax3 = nothing, nothing

    if cfg.type_display == :no_display
        return fig1, ax1, fig2, ax2, fig3, ax3
    end

    # For the step current
    if cfg.type_current isa StepParams
        if cfg.type_display == :multiple  # saving instruction is directly implemented within AlphaPEM.Display here.
            mpl.rcParams["font.size"] = 18  # Font size for all text
            fig1, ax1 = nothing, nothing  # Here, additional plots are unnecessary
            fig2, ax2 = nothing, nothing  # Here, additional plots are unnecessary
            fig3, ax3 = nothing, nothing  # Here, additional plots are unnecessary
        elseif cfg.type_display == :synthetic
            mpl.rcParams["font.size"] = 13  # Font size for all text
            fig1, ax1 = plt.subplots(3, 3; figsize=(14, 14))
            if cfg.type_plot == :fixed
                fig2, ax2 = plt.subplots(; figsize=(8, 8))
            else
                fig2, ax2 = nothing, nothing  # Here, additional plots are unnecessary
            end
            fig3, ax3 = nothing, nothing  # Here, additional plots are unnecessary
            plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        end

    # For the polarization curve
    elseif cfg.type_current isa PolarizationParams
        if cfg.type_display == :multiple
            mpl.rcParams["font.size"] = 11  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3; figsize=(14, 4.7))
            fig2, ax2 = plt.subplots(1, 4; figsize=(18.7, 4.7))
            fig3, ax3 = nothing, nothing  # Here, additional plots are unnecessary
            plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        elseif cfg.type_display == :synthetic
            mpl.rcParams["font.size"] = 18  # Font size for all text
            mpl.rcParams["legend.fontsize"] = 15  # Legend font size only
            fig1, ax1 = plt.subplots(; figsize=(8, 8))
            fig2, ax2 = nothing, nothing  # Here, additional plots are unnecessary
            fig3, ax3 = nothing, nothing  # Here, additional plots are unnecessary
        end

    # For the polarization curve used for calibration
    elseif cfg.type_current isa PolarizationCalibrationParams
        if cfg.type_display == :multiple
            mpl.rcParams["font.size"] = 11  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3; figsize=(14, 4.7))
            fig2, ax2 = nothing, nothing  # Here, additional plots are unnecessary
            fig3, ax3 = nothing, nothing  # Here, additional plots are unnecessary
            plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        elseif cfg.type_display == :synthetic
            mpl.rcParams["font.size"] = 18  # Font size for all text
            fig1, ax1 = plt.subplots(; figsize=(8, 8))
            fig2, ax2 = nothing, nothing  # Here, additional plots are unnecessary
            fig3, ax3 = nothing, nothing  # Here, additional plots are unnecessary
        end

    # For the EIS curve
    elseif cfg.type_current isa EISParams
        if cfg.type_display == :multiple
            mpl.rcParams["font.size"] = 18  # Font size for all text
            fig1, ax1 = plt.subplots(; figsize=(8, 8))
            fig2, ax2 = plt.subplots(; figsize=(8, 8))
            fig3, ax3 = plt.subplots(; figsize=(8, 8))
        elseif cfg.type_display == :synthetic
            mpl.rcParams["font.size"] = 13  # Font size for all text
            fig1, ax1 = plt.subplots(1, 3; figsize=(14, 4.7))
            fig2, ax2 = nothing, nothing  # Here, additional plots are unnecessary
            fig3, ax3 = nothing, nothing  # Here, additional plots are unnecessary
            plt.subplots_adjust(; left=0.04, right=0.98, top=0.96, bottom=0.07, wspace=0.2, hspace=0.15)
        end
    end

    return fig1, ax1, fig2, ax2, fig3, ax3
end


"""Select the `n`-th element from each vector value in a dictionary.

# Arguments
- `d::Dict`: Dictionary where values are vectors or other objects.
- `n::Integer`: 1-based index of the element to select from each vector.

# Returns
- `Dict`: New dictionary with the `n`-th element from each vector,
  or the original value if it is not a vector or the vector is too short.
"""
function select_nth_elements(d::Dict, n::Integer)::Dict
    return Dict(String(k) => ((v isa Vector && length(v) >= n) ? v[n] : v)
                             for (k, v) in pairs(d))
end


"""Extract final internal states from a simulation to initialize the next one.

# Arguments
- `simulator::AlphaPEM`: Simulator instance.
- `type_auxiliary::String`: Auxiliary system type.

# Returns
- `initial_variable_values`: Internal state vector.
"""
function _extract_last_internal_state(simu::AlphaPEM, cfg::SimulationConfig)
    initial_variable_values = []

    np = simu.fuel_cell.numerical_parameters
    mea_names = canonical_mea_solver_variable_names(np.nb_gdl, np.nb_mpl)

    for k in 1:simu.fuel_cell.numerical_parameters.nb_gc
        for key in mea_names
            push!(initial_variable_values, simu.variables[key][k][end])
        end
    end

    if cfg.type_auxiliary in (:forced_convective_cathode_with_flow_through_anode,
                              :forced_convective_cathode_with_anodic_recirculation)
        for key in MANIFOLD_SOLVER_VARIABLE_NAMES
            push!(initial_variable_values, simu.variables[key][end])
        end
        for key in AUXILIARY_SOLVER_VARIABLE_NAMES
            push!(initial_variable_values, simu.variables[key][end])
        end
    end

    return initial_variable_values
end


"""Launch the AlphaPEM simulator for a step current density and display results.

# Arguments
- `simulator::AlphaPEM`: Simulator instance.

# Returns
- `AlphaPEM`: Updated simulator instance.
"""
function launch_AlphaPEM_for_step_current(simu::AlphaPEM)::AlphaPEM

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if simu.cfg.type_plot == :dynamic

        # Certain conditions must be met.
        if simu.cfg.type_display == :multiple
            throw(ArgumentError("step current is not thought to be used with step current, dynamic plot and multiple display. There would be too much plots to handle."))
        end

        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = nothing
        #           Extraction of the parameters
        tf_step = (current_parameters["step_current_parameters"]["delta_t_ini_step"] +
                   current_parameters["step_current_parameters"]["delta_t_load_step"] +
                   current_parameters["step_current_parameters"]["delta_t_break_step"])  # (s).
        delta_t_dyn_step = current_parameters["step_current_parameters"]["delta_t_dyn_step"]  # (s).
        #           Calculation
        n = Int(floor(tf_step / delta_t_dyn_step))  # It is the plot update number.
        time_interval = [0.0, delta_t_dyn_step]  # (s). It is the initial time interval.

        # Dynamic simulation
        for i in 1:n
            simulate_model!(simu, operating_inputs, current_parameters, computing_parameters,
                            initial_variable_values, time_interval)

            # time_interval actualization
            if i < n  # The final simulation does not require actualization.
                t0_interval = simu.variables["t"][end]
                tf_interval = (i + 1) * delta_t_dyn_step
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval
            end

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = _extract_last_internal_state(simu, computing_parameters["type_auxiliary"])

            # Display
            if simu.cfg.type_display != :no_display
                Display(simu, ax1, ax2, ax3)
            end
        end

    else  # elseif simu.cfg.type_plot == :fixed
        # Simulation
        simulate_model!(simu)
        # Display
        if simu.cfg.type_display != :no_display
            Display(simu, ax1, ax2, ax3)
        end
    end

    # Plot saving
    Save_plot(simu, fig1, fig2, fig3)

    return simu
end


"""Launch the AlphaPEM simulator for a polarization current density and display results.

# Arguments


# Returns

"""
function launch_AlphaPEM_for_polarization_current(simu::AlphaPEM)::AlphaPEM

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(simu.cfg)

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if simu.cfg.type_plot == :dynamic
        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = nothing
        delta_t_ini_pola = simu.current_density.delta_t_ini  # (s). It is the initial time at zero current density for the stabilisation of the internal states.
        delta_t_pola = step_duration(simu.current_density)  # s. It is the time of one load.
        _, tf = simu.current_density.time_interval  # s. It is the polarization current duration.
        n = Int(floor(tf / delta_t_pola))  # It is the plot update number.
        time_interval = [0.0, delta_t_ini_pola + delta_t_pola]  # It is the initial time interval.

        # Dynamic simulation
        for i in 1:n
            simulate_model!(simu, initial_variable_values, time_interval)

            # time_interval actualization
            if i < n  # The final simulation does not require actualization.
                t0_interval = simu.variables["t"][end]
                tf_interval = delta_t_ini_pola + (i + 1) * delta_t_pola
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval
            end

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = _extract_last_internal_state(simu, cfg)

            # Display
            if simu.cfg.type_display != :no_display
                Display(simu, ax1, ax2, ax3)
            end
        end

    else  # elseif simu.cfg.type_plot == :fixed
        # Simulation
        simulate_model!(simu)
        # Display
        if simu.cfg.type_display != :no_display
            Display(simu, ax1, ax2, ax3)
        end
    end

    # Plot saving
    Save_plot(simu, fig1, fig2, fig3)

    return simu
end
function launch_AlphaPEM_for_polarization_current(simulators::Vector{AlphaPEM},
                                                  operating_inputs::Dict,
                                                  current_parameters::Dict,
                                                  computing_parameters::Dict)::AlphaPEM

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)

    # Consistency checks for vectors that previously relied on a [None] placeholder in Python.
    type_fuel_cells = computing_parameters["type_fuel_cell"]
    pola_params_list = current_parameters["pola_current_parameters"]
    if !(type_fuel_cells isa Vector) || !(pola_params_list isa Vector)
        throw(ArgumentError("type_fuel_cell and pola_current_parameters must be vectors without placeholder element in Julia."))
    end
    if length(type_fuel_cells) != length(pola_params_list)
        throw(ArgumentError("type_fuel_cell and pola_current_parameters must have the same length."))
    end

    # Condition to fill for comparison with experimental values
    if computing_parameters["type_auxiliary"] == "forced-convective_cathode_with_flow-through_anode"
        for i in eachindex(type_fuel_cells)
            if type_fuel_cells[i] !== nothing && type_fuel_cells[i] != "manual_setup"
                i_exp_t, _ = pola_exp_values(type_fuel_cells[i], computing_parameters["voltage_zone"])
                if pola_params_list[i]["i_max_pola"] < i_exp_t[end]
                    throw(ArgumentError("The given maximum current density of the polarization curve i_max_pola_$(i) is lower than the maximum current density of the experimental values. Please increase it."))
                end
            end
        end
    end
    
    for i in eachindex(simulators)
        # Simulation
        simulate_model!(simulators[i],
                        select_nth_elements(operating_inputs, i),
                        select_nth_elements(current_parameters, i),
                        select_nth_elements(computing_parameters, i))
        # Display
        if simu.cfg.type_display != :no_display
            Display(simulators[i], ax1, ax2, ax3)
        end
    end

    # Plot saving
    Save_plot(simulators[1], fig1, fig2, fig3)

    return simulators[1]
end

"""Launch the AlphaPEM simulator for calibration polarization current and display results.

# Arguments
- `simulators::Vector{<:AlphaPEM}`: Vector of simulator instances (1-based, no placeholder element).
- `operating_inputs::Dict`: Operating inputs.
- `current_parameters::Dict`: Current parameters.
- `computing_parameters::Dict`: Computing parameters.

# Returns
- `AlphaPEM`: Main simulator (`simulators[1]`) after execution.
"""
function launch_AlphaPEM_for_polarization_current_for_calibration(
    simulators::Vector{<:AlphaPEM},
    operating_inputs::Dict,
    current_parameters::Dict,
    computing_parameters::Dict,
)::AlphaPEM

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)

    type_fuel_cells = computing_parameters["type_fuel_cell"]
    if !(type_fuel_cells isa Vector)
        throw(ArgumentError("type_fuel_cell must be a vector without placeholder element in Julia."))
    end

    # Dynamic display requires a dedicated use of the AlphaPEM class.
    if simu.cfg.type_plot == :dynamic
        # Initialization
        #       Calculation of the plot update number (n) and the initial time interval (time_interval).
        initial_variable_values = nothing
        #           Extraction of the parameters
        delta_t_ini_pola_cali = current_parameters["pola_current_for_cali_parameters"]["delta_t_ini_pola_cali"]  # (s).
        delta_t_load_pola_cali = current_parameters["pola_current_for_cali_parameters"]["delta_t_load_pola_cali"]  # (s).
        delta_t_break_pola_cali = current_parameters["pola_current_for_cali_parameters"]["delta_t_break_pola_cali"]  # (s).
        i_exp_cali_t, _ = pola_exp_values_calibration(type_fuel_cells[1], computing_parameters["voltage_zone"])  # (A.m-2, V).
        #           Calculation
        delta_t_pola_cali = delta_t_load_pola_cali + delta_t_break_pola_cali  # s. It is the time of one load.
        tf = delta_t_ini_pola_cali + length(i_exp_cali_t) * delta_t_pola_cali  # s. It is the polarization current duration.
        n = Int(floor(tf / delta_t_pola_cali))  # It is the plot update number.
        time_interval = [0.0, delta_t_ini_pola_cali + delta_t_pola_cali]  # It is the initial time interval.

        # Dynamic simulation
        for i in 1:n
            simulate_model!(simulators[1],
                            select_nth_elements(operating_inputs, 1),
                            select_nth_elements(current_parameters, 1),
                            select_nth_elements(computing_parameters, 1),
                            initial_variable_values,
                            time_interval)

            # time_interval actualization
            if i < n  # The final simulation does not require actualization.
                t0_interval = simulators[1].variables["t"][end]
                tf_interval = delta_t_ini_pola_cali + (i + 1) * delta_t_pola_cali
                time_interval = [t0_interval, tf_interval]  # Reset of the time interval
            end

            # Recovery of the internal states from the end of the preceding simulation.
            initial_variable_values = _extract_last_internal_state(simulators[1], computing_parameters["type_auxiliary"])

            # Display
            if simu.cfg.type_display != :no_display
                Display(simulators[1], ax1, ax2, ax3)
            end
        end

    else  # elseif simu.cfg.type_plot == :fixed

        # Certain conditions must be met.
        if computing_parameters["type_current"] == "polarization_for_cali" &&
           (type_fuel_cells[1] == "manual_setup" ||
            (computing_parameters["type_auxiliary"] != "forced-convective_cathode_with_flow-through_anode" &&
             computing_parameters["type_auxiliary"] != "no_auxiliary"))
            throw(ArgumentError("polarization current for calibration should be done with experimental data."))
        end

        for i in eachindex(simulators)
            # Simulation
            simulate_model!(simulators[i],
                            select_nth_elements(operating_inputs, i),
                            select_nth_elements(current_parameters, i),
                            select_nth_elements(computing_parameters, i))
            # Display
            if simu.cfg.type_display != :no_display
                Display(simulators[i], ax1, ax2, ax3)
            end
        end
    end

    # Plot saving
    Save_plot(simulators[1], fig1, fig2, fig3)

    return simulators[1]
end


"""Launch the AlphaPEM simulator for an EIS current density and display results.

# Arguments
- `simulator::AlphaPEM`: Simulator instance.
- `operating_inputs::Dict`: Operating inputs.
- `current_parameters::Dict`: Current parameters.
- `computing_parameters::Dict`: Computing parameters.

# Returns
- `AlphaPEM`: Updated simulator instance.
"""
function launch_AlphaPEM_for_EIS_current(simulator::AlphaPEM,
                                         operating_inputs::Dict,
                                         current_parameters::Dict,
                                         computing_parameters::Dict)::AlphaPEM

    # Check if the computing_parameters["type_current"] is valid
    if simu.cfg.type_plot != :dynamic
        throw(ArgumentError("EIS has to be plot with a dynamic type_plot setting, because max_step has to be adjusted at each frequency."))
    end

    # Warnings
    println("Warning: EIS simulation is currently not maintained. It should work, but some unexpected bugs or incorrect results may happen. This will be fixed later.")

    # Figures preparation
    fig1, ax1, fig2, ax2, fig3, ax3 = figures_preparation(computing_parameters)

    # Initialization
    #       Calculation of the plot update number (n) and the initial time interval (time_interval).
    initial_variable_values = nothing
    t0_EIS, t_new_start, tf_EIS, delta_t_break_EIS, delta_t_measurement_EIS = current_parameters["t_EIS"]
    f_power_min_EIS, f_power_max_EIS, nb_f_EIS, nb_points_EIS = current_parameters["f_EIS"]  # These are used for EIS max_step actualization.
    f = collect(10.0 .^ range(f_power_min_EIS, f_power_max_EIS; length=Int(nb_f_EIS)))  # It is a list of all tested frequencies.
    n = length(t_new_start)  # It is the plot update number.
    time_interval = [0.0, t0_EIS]  # It is the initial time interval.

    #       A preliminary simulation run is necessary to equilibrate internal variables at i_EIS before initiating EIS.
    simulate_model!(simulator, operating_inputs, current_parameters, computing_parameters,
                    initial_variable_values, time_interval)

    # time_interval actualization
    t0_EIS_temp = t0_EIS  # It is the initial time for 1 EIS point.
    tf_EIS_temp = t_new_start[1] + delta_t_break_EIS[1] + delta_t_measurement_EIS[1]  # It is the final time for 1 EIS point.
    n_inf = findlast(x -> x <= t0_EIS_temp, t_new_start)  # It is the number of frequency changes already made.
    time_interval = [t0_EIS_temp, tf_EIS_temp]

    # Recovery of the internal states from the end of the preceding simulation.
    initial_variable_values = _extract_last_internal_state(simulator, computing_parameters["type_auxiliary"])

    if simu.cfg.type_display == :multiple
        println("A display bug prevents the dynamic updating of the graphs, as it appears that too much data is involved. However, the data is correctly calculated, and the appropriate plots are saved in the 'results' folder. This display bug does not occur when using a 'synthetic' type_display.")
    end

    # Dynamic simulation
    for i in 1:n
        simulate_model!(simulator, operating_inputs, current_parameters, computing_parameters,
                initial_variable_values, time_interval)

        # time_interval actualization
        if i < n  # The final simulation does not require actualization.
            t0_EIS_temp = simulator.variables["t"][end]  # It is the initial time for 1 EIS point.
            tf_EIS_temp = t_new_start[i + 1] + delta_t_break_EIS[i + 1] + delta_t_measurement_EIS[i + 1]  # Final time for 1 EIS point.
            n_inf = findlast(x -> x <= t0_EIS_temp, t_new_start)  # It is the number of frequency changes already made.
            time_interval = [t0_EIS_temp, tf_EIS_temp]  # It is the time interval for 1 EIS point.
        end

        # Recovery of the internal states from the end of the preceding simulation.
        initial_variable_values = _extract_last_internal_state(simulator, computing_parameters["type_auxiliary"])

        # Display
        if simu.cfg.type_display != :no_display
            Display(simulator, ax1, ax2, ax3)
        end
    end

    # Plot saving
    Save_plot(simulator, fig1, fig2, fig3)

    # Keep these variables explicit for readability and parity with Python's current structure.
    _ = tf_EIS
    _ = nb_points_EIS
    _ = f
    _ = n_inf

    return simulator
end

