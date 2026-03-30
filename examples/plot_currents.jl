# examples/plot_currents.jl

using Plots

# --- Include your current models ---
include("../src/alphapem/currents/step.jl")
include("../src/alphapem/currents/polarization.jl")
include("../src/alphapem/currents/eis.jl")
include("../src/alphapem/currents/polarization_for_cali.jl")

# --- Time resolution ---
n_points = 1000  # points for plotting

# === Step Current ===
step_c = StepCurrent(
    30*60.0,     # delta_t_ini = 30 min
    30.0,        # delta_t_load = 30 s
    2*60.0,      # delta_t_break = 2 min
    1.0e4,       # i_ini = 1.0 A/cm²
    2.0e4        # i_step = 2.0 A/cm²
)

t0_step, tf_step = time_interval(step_c)
t_step = range(t0_step, tf_step, length=n_points)
i_step = [current(step_c, t) for t in t_step]

# === Polarization Current ===
pola_c = PolarizationCurrent(
    120*60.0,    # delta_t_ini = 120 min
    0.01e4,      # v_load = 0.01 A.cm-2.s-1
    15*60.0,     # delta_t_break = 15 min
    0.05e4,      # delta_i = 0.05 A/cm²
    3.0e4        # i_max = 3.0 A/cm²
)

t0_pola, tf_pola = time_interval(pola_c)
t_pola = range(t0_pola, tf_pola, length=n_points)
i_pola = [current(pola_c, t) for t in t_pola]

# === Polarization Calibration Current (from experimental data) ===
pola_cali = PolarizationCalibrationCurrent(
    120*60.0,                                # delta_t_ini = 120 min
    0.01e4,                                  # v_load = 0.01 A/cm²/s
    15*60.0,                                 # delta_t_break = 15 min
    [0.0, 0.5e4, 1.0e4, 2.0e4, 2.5e4, 3.0e4] # i_exp = experimental current densities
)

t0_pola_cali, tf_pola_cali = time_interval(pola_cali)
t_pola_cali = range(t0_pola_cali, tf_pola_cali, length=n_points)
i_pola_cali = [current(pola_cali, t) for t in t_pola_cali]

# === EIS Current ===
eis_c = EISCurrent(
    1.0e4,        # i_EIS = 1.0 A/cm²
    0.05,         # ratio = 0.05
    -3.0,         # f_power_min = 10^-3 Hz
    5.0,          # f_power_max = 10^5 Hz
    90,           # nb_f = 90 frequencies
    50            # nb_points = 50 points per frequency
)

t0_eis, tf_eis = time_interval(eis_c)
t_eis = range(t0_eis, tf_eis, length=n_points)
i_eis = [current(eis_c, t) for t in t_eis]

# === Plot Step Current ===
plt_step = plot(
    t_step ./ 60, i_step,
    xlabel = "Time (min)",
    ylabel = "Current density (A/m²)",
    label = "Step Current",
    lw = 2,
    grid = true,
    title = "Step Current Profile"
)
display(plt_step)

# === Plot Polarization Current ===
plt_pola = plot(
    t_pola ./ 60, i_pola,
    xlabel = "Time (min)",
    ylabel = "Current density (A/m²)",
    label = "Polarization Current",
    lw = 2,
    grid = true,
    title = "Polarization Current Profile"
)
display(plt_pola)

# === Plot Polarization Calibration Current ===
plt_pola_cali = plot(
    t_pola_cali ./ 60, i_pola_cali,
    xlabel = "Time (min)",
    ylabel = "Current density (A/m²)",
    label = "Polarization Calibration Current",
    lw = 2,
    grid = true,
    title = "Polarization Calibration Current Profile"
)
display(plt_pola_cali)

 # === Plot EIS Current (commented for now) ===
 plt_eis = plot(
     t_eis ./ 60, i_eis,
     xlabel = "Time (min)",
     ylabel = "Current density (A/m²)",
     label = "EIS Current",
     lw = 2,
     grid = true,
     title = "EIS Current Profile"
 )
 display(plt_eis)
