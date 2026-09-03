"""
    _plot_calibration_results(result, output_dir)

Generate and save a figure showing the experimental vs. simulated polarization curves
for all calibration conditions using the best identified parameters.
"""
function _plot_calibration_results(result, output_dir::String)
    configs = result.config.simulation_configs
    n_configs = length(configs)
    n_configs == 0 && return nothing

    @info "Generating calibration results figure for $n_configs condition(s)..."

    # 1. Setup Figure
    # Adjust figure size based on the number of subplots
    fig = Figure(size = (800, 400 * n_configs))

    try
        for (i, sc) in enumerate(configs)
            # 2. Re-simulate each condition with the best identified parameters
            # Use the specified voltage zone for the final plot
            full_sc = SimulationConfig(
                type_fuel_cell = sc.type_fuel_cell,
                year           = sc.year,
                voltage_zone   = sc.voltage_zone,
                display_timing = :postrun # Ensure we get discretized points
            )

            fc = create_fuelcell(full_sc.type_fuel_cell, full_sc.voltage_zone; year=full_sc.year, nb_gc=sc.numerical_parameters.nb_gc)
            fc.physical_parameters = result.best_params # Use calibrated parameters

            # Use the complete polarization current profile
            cd = create_current(PolarizationParams(), fc)

            # Execute simulation
            sim = AlphaPEM(fc, cd, full_sc)
            simulate_model!(sim)

            # 3. Plot on a new axis
            ax = Axis(fig[i, 1])
            plot_polarization_curve(sim.outputs, fc, cd, full_sc, ax)

            # Add a subtitle or prefix to the axis title if multiple conditions exist
            if n_configs > 1
                ax.title = "Condition $i: " * ax.title[]
            end
        end

        # 4. Save the figure
        out_path = joinpath(output_dir, "calibration_polarization_curves.png")
        save(out_path, fig)
    catch e
        @warn "Failed to generate calibration plots: $e"
        # Don't rethrow, as this is a post-processing step that shouldn't crash the calibration
    end

    return nothing
end
