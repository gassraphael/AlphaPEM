# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :polarization
You can modify any SimulationConfig parameter below.
"""

include("../src/alphapem/config/simulation_config.jl")
include("../src/alphapem/application/run_simulation.jl")
include("../src/alphapem/config/current_parameters.jl")

using .SimulationConfigModule
using .Config: PolarizationParams

current_params = PolarizationParams(
    delta_t_ini = 120.0 * 60.0,  # (s). Initial time at zero current density for the stabilisation of the internal states.
    delta_i = 0.05e4,            # (A.m-2). Current density step for the polarisation current density function.
    v_load = 0.01e4,             # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break = 15.0 * 60.0, # (s). Breaking time for one step current, for the stabilisation of the internal states.
    i_max = 3.0e4                # Maximum current (default value, can be overridden by experimental current values if provided).
)

# List of fuel cell types to simulate. You can add more fuel cell types to the list,
# but make sure they are implemented in the create_fuelcell factory function.
type_fuel_cell_list = [:ZSW_GenStack]

# If only one fuel cell type is selected, run a single simulation. Otherwise, run multiple simulations in parallel.
if length(type_fuel_cell_list) == 1
    cfg = SimulationConfig(
        type_fuel_cell = type_fuel_cell_list[1],
        type_current = current_params,
        voltage_zone = :full,
        type_auxiliary = :no_auxiliary,
        type_purge = :no_purge,
        type_display = :synthetic,
        type_plot = :fixed
    )
    start_time = time() # Starting time
    run_simulation(cfg)
    algo_time = time() - start_time # Ending time
    println("Time of the algorithm in second : ", algo_time)
else
    cfgs = Vector{SimulationConfig}(undef, length(type_fuel_cell_list))
    for i in 1:length(type_fuel_cell_list)
        cfgs[i] = SimulationConfig(
            type_fuel_cell = type_fuel_cell_list[i],
            type_current = current_params,
            voltage_zone = :full,
            type_auxiliary = :no_auxiliary,
            type_purge = :no_purge,
            type_display = :synthetic,
            type_plot = :fixed
        )
    end
    start_time = time() # Starting time
    run_simulation(cfgs)
    algo_time = time() - start_time # Ending time
    println("Time of the algorithm in second : ", algo_time)
end
