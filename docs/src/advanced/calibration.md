# Parameter Calibration

AlphaPEM includes a genetic algorithm-based calibration system that automatically adjusts undetermined 
physical parameters to match experimental fuel cell measurements.

## Overview

Calibration solves the inverse problem: given experimental data (typically a polarization curve), 
find parameter values that make the model predict that data.

**When to calibrate:**
- You have experimental V-I curves for a specific fuel cell
- Default model parameters don't match your cell's behavior
- You need a cell-specific predictive model for control or design

**Process:**
1. Define experimental data (one or more operating conditions)
2. Set parameter bounds and GA hyperparameters
3. Run calibration (genetic algorithm optimizes parameters)
4. Export calibrated parameters to YAML
5. Use calibrated model in future simulations

## Quick Start

Create `calibration_config.jl`:

```julia
using AlphaPEM.Config: SimulationConfig, PolarizationCalibrationParams, CalibrationConfig, GAConfig
using AlphaPEM.Application: run_calibration

# Experimental current densities (A/m²)
i_exp = [0.1e4, 0.5e4, 1.0e4, 1.5e4, 2.0e4, 2.5e4]

calibration_config = CalibrationConfig(
    simulation_configs=[
        SimulationConfig(
            type_fuel_cell=:ZSW_GenStack,
            type_current=PolarizationCalibrationParams(i_exp=i_exp),
            voltage_zone=:before_voltage_drop,
        )
    ],
    ga_config=GAConfig(
        num_generations=200,
        pop_size=100,
        target_error=0.01,  # Stop if RMSE < 1%
    ),
    output_dir="results/calibration"
)

run_calibration(calibration_config)
```

## Configuration Structure

### CalibrationConfig

Top-level configuration combining simulation and GA settings.

| Field                     | Type                       | Default                 | Purpose                                                                                                                |
|---------------------------|----------------------------|-------------------------|------------------------------------------------------------------------------------------------------------------------|
| `simulation_configs`      | `Vector{SimulationConfig}` | One default config      | Experimental conditions (fuel cell, operating conditions) and experimental current data |
| `ga_config`               | `GAConfig`                 | Default GA settings     | Genetic algorithm hyperparameters                                                                                      |
| `parallel`                | `Bool`                     | `true`                  | Use multi-threading for population evaluation                                                                          |
| `initial_population_file` | `String or Nothing`        | `nothing`               | Warm-start from previous calibration checkpoint                                                                        |
| `output_dir`              | `String`                   | `"results/calibration"` | Where to save results and logs                                                                                         |
| `save_frequency`          | `Int`                      | `1`                     | Save checkpoint every N generations                                                                                    |

### PolarizationCalibrationParams

Current profile for calibration (experimental V-I curve).

| Field           | Type              | Default       | Purpose                               |
|-----------------|-------------------|---------------|---------------------------------------|
| `delta_t_ini`   | `Float64`         | 1800.0 s      | Initial stabilization time (zero current) |
| `v_load`        | `Float64`         | 0.01e4 A/m²/s | Loading rate (ramp current smoothly)  |
| `delta_t_break` | `Float64`         | 300.0 s       | Hold time at each current step        |
| `i_exp`         | `Vector{Float64}` | `[]`          | Experimental current densities (A/m²) |

Example with real experimental points:

```julia
PolarizationCalibrationParams(
    i_exp=[0.1e4, 0.3e4, 0.5e4, 1.0e4, 1.5e4, 2.0e4, 2.5e4],
    delta_t_ini=30*60,  # 30 min stabilization
)
```

### GAConfig

Genetic algorithm hyperparameters.

| Field                | Type             | Default            | Purpose                                       |
|----------------------|------------------|--------------------|-----------------------------------------------|
| `num_generations`    | `Int`            | 1000               | Total generations to evolve                   |
| `pop_size`           | `Int`            | 128                | Population size per generation                |
| `num_parents_mating` | `Int`            | `⌊0.2 × pop_size⌋` | Parents selected for reproduction             |
| `mutation_num_genes` | `Int`            | 1                  | Genes mutated per individual (1 = best found) |
| `elitism`            | `Int`            | 1                  | Best individuals carried forward              |
| `target_error`       | `Float64`        | 0.005 (0.5%)       | Stop early if RMSE < threshold                |
| `seed`               | `Int or Nothing` | `nothing`          | Random seed (`nothing` = stochastic runs)     |

### Multi-Condition Calibration

Calibrate against multiple operating conditions simultaneously (stronger constraint):

```julia
calibration_config = CalibrationConfig(
    simulation_configs=[
        # Condition 1: 70°C, low pressure
        SimulationConfig(
            type_fuel_cell=:ZSW_GenStack,
            temperature=70.0,
            pressure_anode=1.0,
            pressure_cathode=1.0,
            type_current=PolarizationCalibrationParams(
                i_exp=[0.5e4, 1.0e4, 1.5e4, 2.0e4]
            ),
        ),
        # Condition 2: 70°C, high pressure (same fuel cell type)
        SimulationConfig(
            type_fuel_cell=:ZSW_GenStack,
            temperature=70.0,
            pressure_anode=3.0,
            pressure_cathode=3.0,
            type_current=PolarizationCalibrationParams(
                i_exp=[0.5e4, 1.0e4, 2.0e4, 2.5e4]
            ),
        ),
    ],
    ga_config=GAConfig(num_generations=300),
)
```

**Important:** All conditions must use the same `type_fuel_cell` and `voltage_zone` 
(they calibrate the same parameter set).

## Warm-Start from Checkpoint

Resume a previous calibration run:

```julia
calibration_config = CalibrationConfig(
    simulation_configs=[...],
    ga_config=GAConfig(num_generations=500),  # Run 500 more generations
    initial_population_file="results/calibration/checkpoint_gen_200.yaml"
)

run_calibration(calibration_config)
```

Checkpoint files are auto-saved based on `save_frequency` and contain the population state.

## Output and Results

After calibration, inspect results in the output directory:

```
results/calibration/
├── calibration_report.yaml          # Summary (best params, fitness, time)
├── final_population.yaml            # All individuals from last generation
├── calibration_checkpoint.yaml      # Latest checkpoint for resume
├── calibrated_bounds.yaml           # Parameter bounds tightened by GA
└── calibration_polarization_curves.png  # Simulated vs. experimental comparison
```

### Calibration Report

Extract and use calibrated parameters:

```julia
using YAML

report = YAML.load(open("results/calibration/calibration_report.yaml"))
bounds  = YAML.load(open("results/calibration/calibrated_bounds.yaml"))

# Results are stored under the "results" and "parameters" keys
best_rmse = report["results"]["min_rmse_achieved"]
best_params = bounds["parameters"]

# Metadata contains execution details
execution_time = report["metadata"]["execution_time_seconds"]

println("Best RMSE: $(best_rmse)%")
println("Execution time: $(execution_time) s")
```

## Practical Guidance

### Preparing Experimental Data

**Format:** Vector of experimental current densities in A/m².

```julia
i_exp = [0.1e4, 0.5e4, 1.0e4, 1.5e4, 2.0e4, 2.5e4]
```

**Density:** Use 5–10 points across the operating range. 
More points provide better constraint but increase computation time.

**Noise:** Experimental data should be smoothed (simple moving average acceptable).

### Choosing Parameter Bounds

Before calibration, the code selects bounds from `fuel_cell_parameters.jl` based on `type_fuel_cell`. 
Bounds define the search space for the GA.

**Strategy:** First run the [parameter validity analysis](validity_analysis.md) to identify feasible ranges, 
then set bounds accordingly.

### Tuning GA Hyperparameters

| Goal                          | Adjustment                                              | Rationale                                                                           |
|-------------------------------|---------------------------------------------------------|-------------------------------------------------------------------------------------|
| **Faster convergence (time)** | Decrease `pop_size`                                     | Fewer simulations per generation (useful if bounds are already tight).              |
| **Better exploration**        | Increase `pop_size` or increase `mutation_num_genes`    | Explores a wider parameter space to avoid local minima.                             |
| **Strict early stopping**     | Increase `target_error` (e.g., from `0.005` to `0.015`) | Stops the algorithm earlier when a satisfactory (but less precise) error is reached. |
| **Prevent overfitting**       | Use multi-condition calibration                         | Calibrates against multiple operating conditions simultaneously.                    |

### Expected Performance

Runtime depends on the number of nodes used for each simulation, amount of simulation and population considered. 
Typical runs are 2-3 days on clusters.

## Troubleshooting

| Issue                             | Diagnosis                                                        | Fix                                                                                  |
|-----------------------------------|------------------------------------------------------------------|--------------------------------------------------------------------------------------|
| **No convergence**                | GA stuck at high error (>4% RMSE)                                | Widen parameter bounds; increase `pop_size` or `num_generations`                     |
| **Oscillating fitness**           | Best fitness fluctuates or degrades between generations          | Increase `elitism` (at least 1) or decrease mutation rate                            |
| **Memory errors / Crashes**       | Memory limit exceeded due to concurrent evaluations              | Reduce the number of parallel threads or reduce mesh nodes in `numerical_parameters` |

## References

- [Valid Parameter Region Analysis](validity_analysis.md) — Identify valid parameter space before calibration
- [Command Line Usage](../user_guide/cli_usage.md) — Running `run_calibration.jl`

---

**Questions?** Contact [raphael.gass@univ-reunion.fr](mailto:raphael.gass@univ-reunion.fr) or see [GitHub Issues](https://github.com/gassraphael/AlphaPEM/issues).
