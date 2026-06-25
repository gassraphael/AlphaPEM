# Valid Parameter Region Analysis

Identifies physically valid regions of the parameter space before calibration. 
This prevents the genetic algorithm from wasting computation on non-physical parameter combinations.

## When to Use

**Optional feature** (requires R + IRD package; see [Installation](../getting_started/installation.md)).

Use when:
- You need to understand which parameter combinations produce valid fuel cell behavior
- You plan multi-condition or complex calibrations and want to restrict the search space
- You want to quantify parameter sensitivity and robustness

Skip if:
- You have trusted bounds from domain experts
- Your parameter space is well-understood

## Workflow

1. **Latin Hypercube Sampling (LHS):** Draw N configurations uniformly from parameter space
2. **Batch Simulation:** Run all configurations in parallel
3. **Validation:** Check each result against physical criteria:
   - Voltage monotonically decreases with current
   - Voltage stays within realistic range 
4. **PRIM/MaxBox:** Identify compact hyperboxes with ≥80% valid configurations
5. **Export:** Save restricted bounds for calibration

## Running Analysis

Example configuration in `run_parameter_validity.jl`:

```julia
using AlphaPEM.Parametrisation.ValidParameterRegion
using AlphaPEM.Parametrisation.ValidParameterRegion: IRDConfig

analysis_cfg = ValidityAnalysisConfig(
    fuel_cell_type=:ZSW_GenStack,
    voltage_zone=:before_voltage_drop,
    n_samples=10_000,
    hyperbox_finder_method = [:PRIM, :MaxBox],
)
ird_cfg = IRDConfig(
    ird_package_dir   = "external/IRD_method_2023/irdpackage",
    probability_range = (0.8, 1.0), # valid region must contain ≥80% valid samples
)
results = run_validity_analysis(analysis_cfg, ird_cfg)
```

## Output Interpretation

Results are saved to `results/model_validity/[FUEL_CELL]/`:

```
results/model_validity/ZSW/
├── bounds_initial.yaml          # Original parameter ranges
├── bounds_restricted.yaml       # Tightened ranges (valid region)
├── parameter_classification.csv # Each sample: valid/invalid flag
├── generated_curves.csv         # Simulated polarization curves
└── final_report.txt             # Summary 
```

### Example: bounds_restricted.yaml

```yaml
Hgdl:
    min: 0.0001
    max: 0.0001496
epsilon_gdl:
    min: 0.7
    max: 0.9
i0_c_ref:
    min: 5.63438
    max: 100.0
```


## Understanding Validity Criteria

**Voltage Monotonicity:**
Cell voltage must decrease (or stay flat) as current increases. 

**Voltage Range:**
Valid range 0.0-1.23 V. Outside this range indicates parameter combinations incompatible with PEMFC physics.

**Convergence:**
Solver must converge (not diverge or hit time limits). Non-convergent solutions are marked invalid.

## Execution Details

AlphaPEM uses Julia's multi-threading for parallel batch simulation. 
Julia automatically uses all available CPU cores by default.

**Performance guidelines:**
- Optimal on machines with 16+ CPU cores or clusters.
- Typical runtime: 8-10 hours for 10_000 samples on modern hardware.
- If memory-constrained, reduce `n_samples` (start with 500–1000, scale up gradually).


## Troubleshooting

| Issue                          | Cause                             | Fix                                                                            |
|--------------------------------|-----------------------------------|--------------------------------------------------------------------------------|
| Very small valid region (<10%) | Parameter ranges too broad        | Use domain knowledge to tighten bounds before analysis                         |
| Analysis hangs                 | Solver divergence on some samples | Increase solver tolerances or reduce `numerical_parameters.nb_gc` (mesh nodes) |
| Memory errors                  | Large `n_samples` × mesh × cores  | Reduce `n_samples` or number of Julia threads                                  |
| No valid configurations found  | Fundamental incompatibility       | Reconsider parameter ranges or model assumptions                               |

## References

- Credits: [README Parametrisation Section](../../README.md)
- IRD Package: [supplementary_2023_ird](https://github.com/slds-lmu/supplementary_2023_ird)
- [Calibration Guide](calibration.md) — Next step after identifying valid region
- [Installation](../getting_started/installation.md) — R + IRD setup

---

**Questions?** Contact [raphael.gass@univ-reunion.fr](mailto:raphael.gass@univ-reunion.fr).
