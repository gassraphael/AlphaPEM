# -*- coding: utf-8 -*-

"""
Example: Run AlphaPEM with type_current = :polarization
You can modify any SimulationConfig parameter below.
"""

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# Must be loaded before AlphaPEM for interactive display (skipped silently on headless servers).
try; import GLMakie; catch; end

using AlphaPEM.Config: SimulationConfig, PolarizationParams, NumericalParams
using AlphaPEM.Application: run_simulation

current_params = PolarizationParams(
    delta_t_ini = 30 * 60.0,  # (s). Initial time at zero current density for the stabilisation of the internal states.
    di_step = 0.05e4,            # (A.m-2). Current density step for the polarisation current density function.
    v_load = 0.01e4,             # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break = 5 * 60.0, # (s). Breaking time for one step current, for the stabilisation of the internal states.
    i_max = 3.0e4                # Maximum current (default value, can be overridden by experimental current values if provided).
)

# List of fuel cell types to simulate. You can add more fuel cell types to the list,
# but make sure they are implemented in the create_fuelcell factory function.
# Examples:
# :ZSW_GenStack, :ZSW_GenStack_Pa_1_61_Pc_1_41, :ZSW_GenStack_Pa_2_01_Pc_1_81,
# :ZSW_GenStack_Pa_2_4_Pc_2_2, :ZSW_GenStack_Pa_2_8_Pc_2_6, :ZSW_GenStack_T_62,
# :ZSW_GenStack_T_76, :ZSW_GenStack_T_84
type_fuel_cell_list = [:ZSW_GenStack]
nb_gc_pola = 1

# If only one fuel cell type is selected, run a single simulation. Otherwise, run multiple simulations in parallel.
if length(type_fuel_cell_list) == 1
    cfg = SimulationConfig(
        type_fuel_cell = type_fuel_cell_list[1],
        type_current = current_params,
        numerical_parameters = NumericalParams(nb_gc = nb_gc_pola),
        voltage_zone = :full, # :before_voltage_drop, :full.
        type_auxiliary = :no_auxiliary, # :forced_convective_cathode_with_anodic_recirculation, :forced_convective_cathode_with_flow_through_anode, :no_auxiliary.
        type_purge = :no_purge, # :constant_purge, :periodic_purge, :no_purge.
        type_display = :synthetic, # :multiple, :synthetic, :no_display.
        display_timing = :postrun # :live, :postrun.
    )
    start_time = time() # Starting time
    run_simulation(cfg)
    algo_time = time() - start_time # Ending time
    println("Time of the algorithm in second : ", algo_time)
else
    cfgs = [
        SimulationConfig(
            type_fuel_cell = type_fuel_cell_list[i],
            type_current = current_params,
            numerical_parameters = NumericalParams(nb_gc = nb_gc_pola),
            voltage_zone = :full, # :before_voltage_drop, :full.
            type_auxiliary = :no_auxiliary, # :forced_convective_cathode_with_anodic_recirculation, :forced_convective_cathode_with_flow_through_anode, :no_auxiliary.
            type_purge = :no_purge, # :constant_purge, :periodic_purge, :no_purge.
            type_display = :synthetic, # :multiple, :synthetic, :no_display.
            display_timing = :postrun # :live, :postrun.
        )
        for i in 1:length(type_fuel_cell_list)
    ]
    start_time = time() # Starting time
    run_simulation(cfgs)
    algo_time = time() - start_time # Ending time
    println("Time of the algorithm in second : ", algo_time)
end
