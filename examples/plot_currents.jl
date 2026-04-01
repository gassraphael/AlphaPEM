# examples/plot_currents.jl

using Plots

# --- Include your current models ---
include("../src/alphapem/config/current_parameters.jl")
include("../src/alphapem/currents/step.jl")
include("../src/alphapem/currents/polarization.jl")
include("../src/alphapem/currents/eis.jl")
include("../src/alphapem/currents/polarization_for_cali.jl")

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
i_step = [current(step_c, t) for t in t_step]

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
