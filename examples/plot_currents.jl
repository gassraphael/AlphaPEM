# examples/plot_currents.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# --- Load AlphaPEM from the project environment ---
using AlphaPEM.Config: StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using AlphaPEM.Currents: StepCurrent, PolarizationCurrent, PolarizationCalibrationCurrent, EISCurrent, current
using AlphaPEM.Core.Models.PlotHelpers: _publication_colors, _finalize_axis!, lsub
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
i_step_vals = [current(step_c, t) for t in t_step] ./ 1e4

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
i_pola = [current(pola_c, t) for t in t_pola] ./ 1e4

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
i_pola_cali = [current(pola_cali, t) for t in t_pola_cali] ./ 1e4

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

# Uniform sampling on [0, tf] causes strong aliasing for high frequencies.
# Build an adaptive grid: `nb_points` samples per local period on each EIS segment.
t_eis = Float64[]

# Include the initial stabilization/ramp with a moderate resolution.
append!(t_eis, collect(range(0.0, stop=t0_eis, length=2000)))

for k in eachindex(eis_c.f)
    t_start = eis_c.t_new_start[k]
    t_stop = t_start + eis_c.delta_t_break[k] + eis_c.delta_t_measurement[k]
    dt = 1.0 / (eis_c.f[k] * eis_p.nb_points)

    t_segment = collect(t_start:dt:t_stop)
    if !isempty(t_eis) && !isempty(t_segment) && t_segment[1] == t_eis[end]
        t_segment = t_segment[2:end]
    end
    append!(t_eis, t_segment)
end

i_eis = current(eis_c, t_eis) ./ 1e4

# === Plot: all current profiles on a 2×2 figure ===
fig = Figure(size = (1200, 800))
palette = _publication_colors()

ax1 = Axis(fig[1, 1], title = "Step current")
lines!(ax1, collect(t_step), i_step_vals, color=:black, linewidth=2.8)
_finalize_axis!(ax1;
    xlabel = rich("Time ", lsub("t", ""), " (s)"),
    ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))

ax2 = Axis(fig[1, 2], title = "Polarization current")
lines!(ax2, collect(t_pola), i_pola, color=palette[1], linewidth=2.8)
_finalize_axis!(ax2;
    xlabel = rich("Time ", lsub("t", ""), " (s)"),
    ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))

ax3 = Axis(fig[2, 1], title = "Polarization calibration current")
lines!(ax3, collect(t_pola_cali), i_pola_cali, color=palette[2], linewidth=2.8)
_finalize_axis!(ax3;
    xlabel = rich("Time ", lsub("t", ""), " (s)"),
    ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))

ax4 = Axis(fig[2, 2], title = "EIS current")
lines!(ax4, t_eis, i_eis, color=palette[3], linewidth=2.8)
_finalize_axis!(ax4;
    xlabel = rich("Time ", lsub("t", ""), " (s)"),
    ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))

display(fig)
