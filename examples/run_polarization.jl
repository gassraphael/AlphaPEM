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
    delta_t_ini   = 30 * 60.0,   # (s). Initial time at zero current density for the stabilisation of the internal states.
    di_step       = 0.05e4,      # (A.m-2). Current density step for the polarisation current density function.
    v_load        = 0.01e4,      # (A.m-2.s-1). Loading rate for one step current of the polarisation current density function.
    delta_t_break = 5 * 60.0,    # (s). Breaking time for one step current, for the stabilisation of the internal states.
    i_max         = NaN          # Maximum current density.
)

# Fuel-cell types and line styles to simulate. Add types implemented by create_fuelcell
# and choose either :solid or :dot for each curve.
# Examples:
# :ZSW_nominal, :ZSW_Pa_1_61_Pc_1_41, :ZSW_Pa_2_01_Pc_1_81,
# :ZSW_Pa_2_4_Pc_2_2, :ZSW_Pa_2_8_Pc_2_6, :ZSW_T_62,
# :ZSW_T_76, :ZSW_T_84
fuel_cell_cases = [
    (type = :ZSW_nominal,         linestyle = :dot),
    (type = :ZSW_Pa_1_61_Pc_1_41, linestyle = :dot),
    (type = :ZSW_Pa_2_01_Pc_1_81, linestyle = :solid),
    (type = :ZSW_Pa_2_4_Pc_2_2,   linestyle = :solid),
    (type = :ZSW_Pa_2_8_Pc_2_6,   linestyle = :solid),
]
year_pola = 2024
nb_gc_pola = 5

# If only one fuel cell type is selected, run a single simulation. Otherwise, run multiple simulations in parallel.
if length(fuel_cell_cases) == 1
    cfg = SimulationConfig(
        type_fuel_cell       = fuel_cell_cases[1].type,
        year                 = year_pola,
        type_current         = current_params,
        numerical_parameters = NumericalParams(nb_gc = nb_gc_pola),
        voltage_zone         = :full,         # :before_voltage_drop, :full.
        type_auxiliary       = :no_auxiliary, # :forced_convective_cathode_with_anodic_recirculation, :forced_convective_cathode_with_flow_through_anode, :no_auxiliary.
        type_flow            = :counter_flow, # :co_flow, :counter_flow.
        type_purge           = :no_purge,     # :constant_purge, :periodic_purge, :no_purge.
        type_display         = :synthetic,    # :multiple, :synthetic, :no_display.
        display_timing       = :postrun       # :live, :postrun.
    )
    start_time = time() # Starting time
    run_simulation(cfg; linestyle=fuel_cell_cases[1].linestyle)
    algo_time = time() - start_time # Ending time
    println("Time of the algorithm in second : ", algo_time)
else
    cfgs = [
        SimulationConfig(
            type_fuel_cell       = case.type,
            type_current         = current_params,
            numerical_parameters = NumericalParams(nb_gc = nb_gc_pola),
            voltage_zone         = :full,                # :before_voltage_drop, :full.
            type_auxiliary       = :no_auxiliary,        # :forced_convective_cathode_with_anodic_recirculation, :forced_convective_cathode_with_flow_through_anode, :no_auxiliary.
            type_flow            = :counter_flow,        # :co_flow, :counter_flow.
            type_purge           = :no_purge,            # :constant_purge, :periodic_purge, :no_purge.
            type_display         = :synthetic,           # :multiple, :synthetic, :no_display.
            display_timing       = :postrun              # :live, :postrun.
        )
        for case in fuel_cell_cases
    ]
    start_time = time() # Starting time
    run_simulation(cfgs; linestyles=[case.linestyle for case in fuel_cell_cases])
    algo_time = time() - start_time # Ending time
    println("Time of the algorithm in second : ", algo_time)
end
