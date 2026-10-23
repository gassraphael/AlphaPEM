# examples/plot_currents.jl

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# --- Load AlphaPEM from the project environment ---
using AlphaPEM.Config: StepParams, PolarizationParams, PolarizationCalibrationParams, EISParams
using AlphaPEM.Currents: StepCurrent, PolarizationCurrent, PolarizationCalibrationCurrent, EISCurrent, current
using AlphaPEM.Core.Models.PlotHelpers: _publication_colors, _cell_current_color, _finalize_axis!, lsub, saving_instructions!
using GLMakie
GLMakie.activate!()  # Enable interactive window

# --- Time resolution ---
n_points = 10000  # points for plotting
show_polarization_calibration = false

# === Step Current ===
step_p = StepParams(
    delta_t_ini   = 30*60.0, # delta_t_ini = 30 min
    delta_t_load  = 30.0,    # delta_t_load = 30 s
    delta_t_break = 2*60.0,  # delta_t_break = 2 min
    i_ini         = 1.0e4,   # i_ini = 1.0 A/cm²
    i_step        = 1.5e4    # i_step = 2.0 A/cm²
)
step_c = StepCurrent(step_p)

t0_step, tf_step = step_c.time_interval
t_step = range(t0_step, tf_step, length=n_points)
t_step_min = t_step ./ 60
i_step_vals = [current(step_c, t) for t in t_step] ./ 1e4

# === Polarization Current ===
pola_p = PolarizationParams(
    delta_t_ini   = 30*60.0,  # delta_t_ini = 120 min
    v_load        = 0.01e4,   # v_load = 0.01 A.cm-2.s-1
    delta_t_break = 300,      # delta_t_break = 15 min
    di_step       = 0.05e4,   # di_step = 0.05 A/cm²
    i_max         = 2.5e4     # i_max = 2.5 A/cm²
)
pola_c = PolarizationCurrent(pola_p)

t0_pola, tf_pola = pola_c.time_interval
t_pola = range(t0_pola, tf_pola, length=n_points)
t_pola_h = t_pola ./ 3600
i_pola = [current(pola_c, t) for t in t_pola] ./ 1e4

# === Polarization Calibration Current (from experimental data) ===
pola_cali_p = PolarizationCalibrationParams(
    delta_t_ini   = 30*60.0,                                  # delta_t_ini = 120 min
    v_load        = 0.01e4,                                   # v_load = 0.01 A/cm²/s
    delta_t_break = 300,                                      # delta_t_break = 15 min
    i_exp         = [0.0, 0.25e4, 0.5e4, 1.0e4, 2.0e4, 2.5e4] # i_exp = experimental current densities
)
pola_cali = PolarizationCalibrationCurrent(pola_cali_p)

t0_pola_cali, tf_pola_cali = pola_cali.time_interval
t_pola_cali = range(t0_pola_cali, tf_pola_cali, length=n_points)
t_pola_cali_min = t_pola_cali ./ 60
i_pola_cali = [current(pola_cali, t) for t in t_pola_cali] ./ 1e4

# === EIS Current ===
eis_p = EISParams(
    i_EIS       = 1.0e4, # i_EIS = 1.0 A/cm²
    ratio       = 0.05,  # ratio = 0.05
    f_power_min = -3.0,  # f_power_min = 10^-3 Hz
    f_power_max = 5.0,   # f_power_max = 10^5 Hz
    nb_f        = 90,    # nb_f = 90 frequencies
    nb_points   = 50     # nb_points = 50 points per frequency
)
eis_c = EISCurrent(eis_p)

t0_eis, tf_eis = eis_c.time_interval

# Uniform sampling on [0, tf] causes strong aliasing for high frequencies.
# Build an adaptive grid: `nb_points` samples per local period on each EIS segment.
t_eis = Float64[]

# Include the initial stabilization/ramp with a moderate resolution.
append!(t_eis, collect(range(0.0, stop=eis_c.t0, length=2000)))

for k in eachindex(eis_c.f)
    t_start = eis_c.t_new_start[k]
    t_stop = t_start + eis_c.delta_t_break[k] + eis_c.delta_t_measurement[k]
    dt = 1.0 / (eis_c.f[k] * eis_p.nb_points)

    t_segment = collect(t_start:dt:t_stop)
    append!(t_eis, t_segment)
end

t_eis_h = t_eis ./ 3600
i_eis = [current(eis_c, t) for t in t_eis] ./ 1e4

if show_polarization_calibration || length(eis_c.f) < 2
    n_zoom_periods = 3
    zoom_start = eis_c.t_new_start[1]
    zoom_stop = zoom_start + n_zoom_periods / eis_c.f[1]
    zoom_title = "First $n_zoom_periods oscillations"
else
    n_zoom_periods = 5
    frequency_change = eis_c.t_new_start[2]
    zoom_start = frequency_change - n_zoom_periods / eis_c.f[1]
    zoom_stop = frequency_change + n_zoom_periods / eis_c.f[2]
    zoom_title = "$n_zoom_periods oscillations before and after one frequency change"
end
zoom_indices = findall(t -> zoom_start <= t <= zoom_stop, t_eis)
t_eis_zoom_h = t_eis[zoom_indices] ./ 3600
i_eis_zoom = i_eis[zoom_indices]

_one_decimal_tick_labels(values) = [string(round(value; digits=1)) for value in values]

# === Plot: all current profiles on a 2×2 figure (EIS with a zoom inset) ===
fig = Figure(size = (1200, 800))
palette = _publication_colors()

ax1 = Axis(fig[1, 1], title = "Step current")
lines!(ax1, t_step_min, i_step_vals, color=_cell_current_color(1), linewidth=2.8)
_finalize_axis!(ax1;
    xlabel = rich("Time ", lsub("t", ""), " (min)"),
    ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))
ax1.xticks = range(first(t_step_min), last(t_step_min); length=6)
ax1.xtickformat = _one_decimal_tick_labels

ax2 = Axis(fig[1, 2], title = "Polarization current")
lines!(ax2, t_pola_h, i_pola, color=palette[1], linewidth=2.8)
_finalize_axis!(ax2;
    xlabel = rich("Time ", lsub("t", ""), " (h)"),
    ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))
ax2.xticks = range(first(t_pola_h), last(t_pola_h); length=6)
ax2.xtickformat = _one_decimal_tick_labels

if show_polarization_calibration
    ax3 = Axis(fig[2, 1], title = "Polarization calibration current")
    lines!(ax3, t_pola_cali_min, i_pola_cali, color=palette[2], linewidth=2.8)
    _finalize_axis!(ax3;
        xlabel = rich("Time ", lsub("t", ""), " (min)"),
        ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))
    ax3.xticks = range(first(t_pola_cali_min), last(t_pola_cali_min); length=6)
    ax3.xtickformat = _one_decimal_tick_labels
end

eis_position = show_polarization_calibration ? fig[2, 2] : fig[2, 1:2]
ax4 = Axis(eis_position, title = "EIS current")
lines!(ax4, t_eis_h, i_eis, color=palette[3], linewidth=2.8)
_finalize_axis!(ax4;
    xlabel = rich("Time ", lsub("t", ""), " (h)"),
    ylabel = rich("Current density ", lsub("i", "fc"), " (A·cm⁻²)"))
ax4.xticks = range(first(t_eis_h), last(t_eis_h); length=6)
ax4.xtickformat = _one_decimal_tick_labels

# Frame around the zoomed region on the main EIS curve
zoom_x = [zoom_start, zoom_stop, zoom_stop, zoom_start, zoom_start] ./ 3600
zoom_y = [minimum(i_eis_zoom), minimum(i_eis_zoom), maximum(i_eis_zoom), maximum(i_eis_zoom), minimum(i_eis_zoom)]
lines!(ax4, zoom_x, zoom_y; color=:orange, linewidth=1.5)

# Inset axis showing the first oscillations or the first frequency change
ax_zoom = Axis(eis_position,
    width=Relative(show_polarization_calibration ? 0.5 : 0.42),
    height=Relative(show_polarization_calibration ? 0.4 : 0.5),
    halign=0.95, valign=0.4, tellwidth=false, tellheight=false,
    backgroundcolor=:white, title = zoom_title, titlesize=13)
translate!(ax_zoom.blockscene, 0, 0, 100)
lines!(ax_zoom, t_eis_zoom_h, i_eis_zoom, color=palette[3], linewidth=2.0)
ax_zoom.xticks = range(first(t_eis_zoom_h), last(t_eis_zoom_h); length=4)
ax_zoom.xtickformat = _one_decimal_tick_labels
ax_zoom.xminorticks = IntervalsBetween(5)
ax_zoom.xminorticksvisible = true
ax_zoom.yminorticks = IntervalsBetween(5)
ax_zoom.yminorticksvisible = true
ax_zoom.xgridvisible = false
ax_zoom.ygridvisible = false

# Connectors from the zoom frame to the inset corners (pixel space, kept in sync with the layout)
connector_points = lift(ax4.scene.viewport, ax_zoom.scene.viewport, ax4.finallimits) do vp4, vpz, lims
    to_px(x, y) = Point2f((x - lims.origin[1]) / lims.widths[1] * vp4.widths[1],
                          (y - lims.origin[2]) / lims.widths[2] * vp4.widths[2])
    inset_left   = vpz.origin[1] - vp4.origin[1]
    inset_bottom = vpz.origin[2] - vp4.origin[2]
    inset_top    = inset_bottom + vpz.widths[2]
    return [to_px(zoom_x[2], zoom_y[3]), Point2f(inset_left, inset_top),
            Point2f(NaN, NaN),
            to_px(zoom_x[2], zoom_y[1]), Point2f(inset_left, inset_bottom)]
end
lines!(ax4, connector_points; space=:pixel, color=:orange, linestyle=:dash, linewidth=1.2)

for spine in (:leftspinecolor, :rightspinecolor, :topspinecolor, :bottomspinecolor)
    setproperty!(ax_zoom, spine, :orange)
end

# === Display and Save ===
# Open an interactive window (GLMakie)
DataInspector(fig) # Add hover inspection
screen = GLMakie.Screen()
display(screen, fig)

# Save to the results directory (CairoMakie PDF)
saving_instructions!(nothing, "results", "currents", "plot_currents.pdf", fig)
