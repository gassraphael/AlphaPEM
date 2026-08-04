# Sobol Global Sensitivity Analysis

Variance-based global sensitivity analysis of AlphaPEM using the Sobol method.
This module quantifies how undetermined physical parameters and operating conditions
influence the polarization curve in three characteristic regions: activation, ohmic,
and mass transport.

## Overview

The analysis runs AlphaPEM directly for each Sobol sample.

**What it computes:**
- First-order indices `S1`: direct influence of each parameter.
- Total-order indices `ST`: direct influence plus all interactions with other parameters.
- Optional second-order indices `S2`: pairwise interaction effects (very expensive).

**Output aggregation:**
Indices are computed separately for the AUC (area under the curve) of each polarization
region, so you can see which parameters matter where.

## When to Use

Use when:
- You want to rank parameters by influence on cell voltage.
- You need to know whether activation, ohmic, or mass-transport losses are most sensitive
to a given parameter.
- You want to reduce the parameter space before calibration.
- You want to study the effect of operating conditions (temperature, pressure, humidity,
stoichiometry) alongside physical parameters.

Skip if:
- You already have tight domain bounds and a clear idea of the dominant parameters.
- Compute budget is too small for the required number of simulations.

## Output Metric: Area Under the Curve (AUC)

Each simulation produces a full polarization curve `U = f(i)`. To apply variance-based
sensitivity analysis, this curve must be reduced to one or more scalar outputs.

The chosen scalar is the **area under the polarization curve** (AUC) inside each
characteristic current-density region:

```
AUC_region = ∫ U(i) di   for i in [i_min, i_max]
```

AUC is computed by trapezoidal integration after clipping negative voltages to zero if needed.

**Why AUC?**
- It condenses the voltage behaviour over an entire region into a single number.
- It is smooth with respect to parameter changes, which improves Sobol convergence.
- It avoids relying on a single current-density point, which can be noisy or unrepresentative.

If you need point-wise sensitivities, you can adapt the output matrix in the module,
but AUC is the recommended default.

## Workflow

1. **Define inputs:** physical undetermined parameters + optional operating conditions.
2. **Generate Sobol matrices:** `A` and `B` using `QuasiMonteCarlo.SobolSample` with a
   randomised shift.
3. **Fuse designs:** `GlobalSensitivity.fuse_designs(A, B)` produces the full set of
   parameter vectors to evaluate.
4. **Validation:** discard physically inconsistent operating conditions.
5. **Batch simulation:** run AlphaPEM for every valid column of the fused design.
6. **Imputation:** replace failed curves by the average of the `k` nearest valid neighbours
   in normalized parameter space.
7. **Compute indices:** call `GlobalSensitivity.gsa` with the pre-computed AUC matrix.
8. **Export and plot:** CSV tables, YAML summary, bar plots, and ranking plots.

## Running the Analysis

```bash
julia --project=. examples/run_sobol_sensitivity_analysis.jl
```

### Minimal configuration

```julia
using AlphaPEM.Parametrisation.SobolSensitivityAnalysis

cfg = SobolAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    voltage_zone   = :full,
    N              = 1024,
)

result = run_sobol_analysis(cfg)
```

### Fix operating conditions to a single value

Use `excluded_operating_conditions` to keep parameters at their nominal value and remove
them from the Sobol sampling. This is useful when you want to study physical parameters
under a fixed experimental condition, e.g. fixed stoichiometries:

```julia
cfg = SobolAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    voltage_zone   = :full,
    N              = 1024,
    excluded_operating_conditions = [:Sa, :Sc],
)
```

### Enforce relationships between operating conditions

Use `OperatingConditionConstraint` to express post-sampling relationships. The default constraint enforces
for example `Pc_des = Pa_des - 0.5e5`. You can add your own:

```julia
using AlphaPEM.Config: OperatingConditionConstraint

my_constraints = [
    OperatingConditionConstraint(
        target  = :Pc_des,
        sources = [:Pa_des],
        fn      = d -> d[:Pa_des] - 0.5e5
    ),
]

cfg = SobolAnalysisConfig(
    fuel_cell_type = :ZSW_GenStack,
    N              = 1024,
    operating_condition_constraints = my_constraints,
)
```

### Adjust region thresholds

The thresholds separating activation, ohmic, and mass transport are configurable:

```julia
cfg = SobolAnalysisConfig(
    fuel_cell_type    = :ZSW_GenStack,
    region_thresholds = (0.4, 1.6),  # A/cm²
)
```

These values should be chosen based on the experimental polarization curve. 


## Output Interpretation

Results are saved to `results/sobol_sensitivity/`:

```
results/sobol_sensitivity/
├── sobol_indices_activation.csv
├── sobol_indices_ohmic.csv
├── sobol_indices_mass_transport.csv
├── sobol_summary.yaml
├── sobol_indices_ST.png
├── sobol_ranking_ST.png
└── raw_curves.csv          # if save_curves = true
```

### `sobol_indices_<region>.csv`

| Column | Meaning |
|--------|---------|
| `parameter` | Input parameter name |
| `S1` | First-order Sobol index |
| `S1_conf` | 95 % confidence interval for `S1` |
| `ST` | Total-order Sobol index |
| `ST_conf` | 95 % confidence interval for `ST` |

Interpretation guidelines:
- `ST ≈ 0`: the parameter has negligible influence.
- `S1 ≈ ST`: the parameter acts mostly independently.
- `ST ≫ S1`: the parameter has strong interactions with others.
- `ST > 0.1`: typically considered an influential parameter.

### `sobol_summary.yaml`

Contains metadata (fuel cell, thresholds, execution time), the list of input parameters,
and the top five parameters per region ranked by `ST`.

## Execution Details

AlphaPEM uses Julia's multi-threading for parallel batch simulation. The number of threads
is controlled when launching Julia:

```bash
julia --threads=auto --project=. examples/run_sobol_sensitivity_analysis.jl
```

**Simulation count:**
- For `S1` and `ST` only: `N × (2 + n_params)` simulations.
- With `second_order = true`: much larger, roughly `N × (2 + 2 × n_params + n_params × (n_params - 1) / 2)`.

With the default `N = 10_000` and about 20 input parameters, the S1/ST analysis requires
roughly 220 000 simulations. Start with `N = 512` or `N = 1024` to validate the workflow,
then increase `N` for publication-quality indices.

**Performance guidelines:**
- Use `nb_gc = 1` in `numerical_params` for fastest batch runs.
- Set `max_run_time_s` low enough to avoid hanging on divergent samples, but high enough
  for normal convergence.
- If many simulations fail, tighten parameter bounds or run a [valid parameter region analysis](../advanced/validity_analysis.md) first.

## Troubleshooting

| Issue | Cause | Fix |
|-------|-------|-----|
| Many failed simulations | Parameter combinations are non-physical or too extreme | Tighten bounds, run validity analysis first, or increase `max_run_time_s` |
| KNN imputation fails | Too few valid neighbours | Reduce `knn_k` or increase `N` |
| Sobol indices do not converge | `N` too small | Increase `N`; check that confidence intervals shrink |
| Out of memory | Large `N` with `save_curves = true` | Set `save_curves = false` or reduce `N` |
| Second-order analysis too slow | Combinatorial cost | Use `second_order = false` for screening, then enable S2 only for a reduced parameter set |

## Integration with Calibration

The recommended parametrisation workflow is:

1. **[Valid Parameter Region Analysis](../advanced/validity_analysis.md)** — restrict the parameter space.
2. **Sobol Sensitivity Analysis** — identify the most influential parameters.
3. **[Calibration](../advanced/calibration.md)** — focus the genetic algorithm on the restricted, high-influence parameters.

---

**Questions?** Contact [raphael.gass@univ-reunion.fr](mailto:raphael.gass@univ-reunion.fr).
