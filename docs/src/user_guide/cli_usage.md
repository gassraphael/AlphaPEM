# Command line interface

For advanced users and batch processing, AlphaPEM provides full programmatic access via Julia interface. 
This guide covers command-line simulations and integration into custom workflows.

## Running example scripts

The `examples/` directory contains ready-to-run Julia scripts demonstrating all simulation types and features.

### Step response simulation

Observe the dynamic response to a constant current step:

```bash
julia --project=. examples/run_step.jl
```

**Configuration in `run_step.jl`:**

```julia
config = SimulationConfig(
    type_fuel_cell=:ZSW_GenStack,
    type_current=StepParams(i_ini = 1.0e4, i_step = 1.5e4),
    voltage_zone=:before_voltage_drop,
    # ... additional parameters
)
```

**Output:** Time-series data (voltage, temperature, pressures, water content, ...)

### Polarization curve

Generate steady-state voltage vs. current characteristic:

```bash
julia --project=. examples/run_polarization.jl
```

**Configuration in `run_polarization.jl`:**

```julia
config = SimulationConfig(
    type_fuel_cell=:ZSW_GenStack,
    type_current=PolarizationParams(
        di_step = 0.05e4,
        v_load = 0.01e4,
        i_max = 3.0e4
    ),
    # ... additional parameters
)
```

**Output:** V-I curve (voltage and power vs. current density)

### Electrochemical impedance spectroscopy

Characterize internal impedance across frequency range:

```bash
julia --project=. examples/run_EIS.jl
```

**Configuration in `run_EIS.jl`:**

```julia
config = SimulationConfig(
    type_fuel_cell=:ZSW_GenStack,
    type_current=EISParams(
        i_EIS = 1.0e4,
        ratio = 5.0 / 100.0,
        f_power_min = -3.0,     
        f_power_max = 5.0 
    ),
    # ... additional parameters
)
```

**Output:** Nyquist plot data, Bode magnitude and phase

### Plotting current profiles

Visualize current density profiles before simulation:

```bash
julia --project=. examples/plot_currents.jl
```

Useful for verifying profile shapes and timing parameters.

## Valid parameter region analysis

Identify physically valid undetermined parameter bounds (requires optional R installation):

```bash
julia --project=. examples/run_parameter_validity.jl
```

**Workflow:**
1. Latin hypercube sampling of parameter space
2. Batch simulation and validation of results
3. Identify compact valid regions using PRIM/MaxBox
4. Export bounds for calibration

**Output:**
- Valid parameter region bounds (YAML file)
- Report (TXT file)



## Parameter calibration (Genetic Algorithm)

Automatically adjust undetermined parameters to match experimental data using a genetic algorithm:

```bash
julia --project=. examples/run_calibration.jl
```

**Configuration in `run_calibration.jl`:**

```julia
calibration_config = CalibrationConfig(
    calibration_conditions=[SimulationConfig(
                                type_fuel_cell = :ZSW_GenStack,
                                voltage_zone   = :before_voltage_drop,)],
    ga_config=GAConfig(
        num_generations = 500,
        population_size = 128,
        target_error    = 1/100,),
    # ... additional parameters
)
```

**Output:**
- Calibrated parameter set (YAML file)
- Report with history (YAML file)
- Plot comparing simulation vs. experimental data (PNG files)

See [Calibration Guide](../advanced/calibration.md) for detailed workflows.

## Programmatic access (Julia)

Use AlphaPEM as a library in your own Julia code:

```julia
using AlphaPEM
using AlphaPEM: SimulationConfig, StepParams, run_simulation
using AlphaPEM.Core.Models: extract_mid_mea_series

# Define configuration
config = SimulationConfig(
    type_fuel_cell=:ZSW_GenStack,
    type_current=StepParams(i_ini = 1.0e4, i_step = 1.5e4),
)

# Run simulation
simu = run_simulation(config)

# Access results programmatically
U = simu.outputs.derived.Ucell      # Voltage time series (V)
i = simu.outputs.derived.i_fc       # Current density (A/m²)
T = extract_mid_mea_series(simu.outputs, config, mea -> mea.acl.T) # Temperature (K)
# ... additional fields
```

## Integration with Python

Call AlphaPEM from Python using `juliacall`:

```bash
pip install juliacall
```

```python
from juliacall import Main as jl

# Load AlphaPEM
jl.seval("using AlphaPEM")

# Configure simulation (via Julia)
jl.seval("""
config = SimulationConfig(
    type_fuel_cell=:ZSW_GenStack,
    type_current=StepParams(i_ini = 1.0e4, i_step = 1.5e4),
)
""")

# Run simulation
simu = jl.run_simulation(jl.config)

# Access results in Python
voltage = simu.outputs.derived.Ucell
current = simu.outputs.derived.i_fc
print(f"Peak voltage: {max(voltage)} V")
```

This approach enables seamless polyglot workflows combining Julia's performance with Python's ecosystem.

## Output data organization

Simulation results are saved to the `results/` directory, organized by operation type and fuel cell:

```
results/
├── [FUEL_CELL_TYPE]/              # Step, polarization, EIS simulations
│   ├── step_current_syn_1.pdf
│   ├── pola_curve_1.pdf
│   ├── Nyquist_plot_1.pdf
│   └── ... (various PDF plots)
├── calibration/
│   └── [FUEL_CELL_TYPE]/           # Calibration results
│       ├── calibration_report.yaml  # Calibration summary
│       ├── final_population.yaml    # Best parameters found
│       ├── calibrated_bounds.yaml   # Restricted parameter bounds
│       ├── calibration_checkpoint.yaml
│       ├── calibration_checkpoint_population.yaml
│       └── calibration_polarization_curves.png
├── model_validity/
│   └── [EXPERIMENT_NAME]/           # Parameter validity analysis
│       ├── bounds_initial.yaml
│       ├── bounds_restricted.yaml
│       ├── parameter_classification.csv
│       └── generated_curves.csv
├── benchmark/                        # Benchmark simulations
├── currents/                         # Current profile plots
└── web_plots/
    └── [RESULT_ID]/                  # Web app downloaded plots
        ├── *.pdf                     # Downloaded simulation plots
        └── *.xlsx                    # Downloaded data tables
```


**Plot formats:** Simulation results are saved as PDF files.
**YAML data files:** Contain parameter sets, bounds, and metadata (used for calibration initialization and 
result analysis).

## Troubleshooting

| Issue | Solution                                              |
|-------|-------------------------------------------------------|
| Slow execution on first run | Julia JIT compilation is normal; first run is slower  |
| Out-of-memory errors | Reduce mesh nodes or simulation resolution            |
| Parameter calibration not converging | Widen parameter number/bounds or increase generations |

## Next steps

- **Advanced calibration**: Read [Calibration guide](../advanced/calibration.md)

---

**Questions?** Contact [raphael.gass@univ-reunion.fr](mailto:raphael.gass@univ-reunion.fr) or visit [GitHub issues](https://github.com/gassraphael/AlphaPEM/issues).
