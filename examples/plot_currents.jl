# examples/plot_currents.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# --- Load AlphaPEM from the project environment ---
using AlphaPEM.Config: StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using AlphaPEM.Currents: StepCurrent, PolarizationCurrent, PolarizationCalibrationCurrent, EISCurrent, current
using CairoMakie

# --- Time resolution ---
n_points = 1000  # points for plotting

# === Step Current ===
step_p = StepParams(
    delta_t_ini = 30*60.0,  # delta_t_ini = 30 min
    delta_t_load = 30.0,    # delta_t_load = 30 s
    delta_t_break = 2*60.0, # delta_t_break = 2 min
    i_ini = 1.0e4,          # i_ini = 1.0 A/cm²
    i_step = 2.0e4          # i_step = 2.0 A/cm²
)
step_c = StepCurrent(step_p)

t0_step, tf_step = step_c.time_interval
t_step = range(t0_step, tf_step, length=n_points)
i_step_vals = [current(step_c, t) for t in t_step]

# === Polarization Current ===
pola_p = PolarizationParams(
    delta_t_ini = 120*60.0,  # delta_t_ini = 120 min
    v_load = 0.01e4,         # v_load = 0.01 A.cm-2.s-1
    delta_t_break = 15*60.0, # delta_t_break = 15 min
    delta_i = 0.05e4,        # delta_i = 0.05 A/cm²
    i_max = 3.0e4            # i_max = 3.0 A/cm²
)
pola_c = PolarizationCurrent(pola_p)

t0_pola, tf_pola = pola_c.time_interval
t_pola = range(t0_pola, tf_pola, length=n_points)
i_pola = [current(pola_c, t) for t in t_pola]

# === Polarization Calibration Current (from experimental data) ===
pola_cali_p = PolarizationCalibrationParams(
    delta_t_ini = 120*60.0,                                    # delta_t_ini = 120 min
    v_load = 0.01e4,                                           # v_load = 0.01 A/cm²/s
    delta_t_break = 15*60.0,                                   # delta_t_break = 15 min
    i_exp = [0.0, 0.5e4, 1.0e4, 2.0e4, 2.5e4, 3.0e4]          # i_exp = experimental current densities
)
pola_cali = PolarizationCalibrationCurrent(pola_cali_p)

t0_pola_cali, tf_pola_cali = pola_cali.time_interval
t_pola_cali = range(t0_pola_cali, tf_pola_cali, length=n_points)
i_pola_cali = [current(pola_cali, t) for t in t_pola_cali]

# === EIS Current ===
eis_p = EISParams(
    i_EIS = 1.0e4,      # i_EIS = 1.0 A/cm²
    ratio = 0.05,       # ratio = 0.05
    f_power_min = -3.0, # f_power_min = 10^-3 Hz
    f_power_max = 5.0,  # f_power_max = 10^5 Hz
    nb_f = 90,          # nb_f = 90 frequencies
    nb_points = 50      # nb_points = 50 points per frequency
)
eis_c = EISCurrent(eis_p)

t0_eis, tf_eis = eis_c.time_interval
t_eis = range(t0_eis, tf_eis, length=n_points)
i_eis = [current(eis_c, t) for t in t_eis]

# === Plot: all current profiles on a 2×2 figure ===
fig = Figure(size = (1200, 800))

ax1 = Axis(fig[1, 1],
    title  = "Step Current Profile",
    xlabel = "Time (min)",
    ylabel = "Current density (A m⁻²)")
lines!(ax1, collect(t_step) ./ 60, i_step_vals, linewidth = 2)

ax2 = Axis(fig[1, 2],
    title  = "Polarization Current Profile",
    xlabel = "Time (min)",
    ylabel = "Current density (A m⁻²)")
lines!(ax2, collect(t_pola) ./ 60, i_pola, linewidth = 2)

ax3 = Axis(fig[2, 1],
    title  = "Polarization Calibration Current Profile",
    xlabel = "Time (min)",
    ylabel = "Current density (A m⁻²)")
lines!(ax3, collect(t_pola_cali) ./ 60, i_pola_cali, linewidth = 2)

ax4 = Axis(fig[2, 2],
    title  = "EIS Current Profile",
    xlabel = "Time (min)",
    ylabel = "Current density (A m⁻²)")
lines!(ax4, collect(t_eis) ./ 60, i_eis, linewidth = 2)

display(fig)
