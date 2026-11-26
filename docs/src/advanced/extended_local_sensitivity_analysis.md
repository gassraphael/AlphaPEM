# Extended Local Sensitivity Analysis

Local sensitivity analysis repeated at many points of the parameter space, followed by the selection of
the parameters that can be identified from the data. The method follows Goshtasbi et al. (2020).

## Overview

A **local sensitivity** measures how much the simulated voltage changes when one parameter is slightly
modified (normalised finite differences). Because it depends on where the other parameters are, it is
evaluated at many candidate parameter sets drawn with a Sobol sequence, and for a factorial library of
729 operating conditions.

**What it computes:**
- Sensitivities of the voltage to every parameter, for every candidate and operating condition.
- Their magnitude (how much a parameter matters) and their collinearity (whether two parameters have
  the same effect and cannot be told apart).
- A reduced set of parameters to calibrate.

The analysis is **practical**, not structural: its conclusions depend on the model, the parameter bounds
and the measurements, and do not by themselves establish the physical importance of a parameter.

## Workflow

`run_extended_local_sensitivity_analysis(cfg)` coordinates the following steps.

1. **Define the parameter space:** choose the parameters to analyse and fix integer-valued parameters.
2. **Generate and screen candidates:** draw Sobol parameter sets and reject those that do not reproduce the experimental polarization curves.
3. **Build the operating-condition set:** create the factorial library of operating cases.
4. **Prepare live figures:** create the convergence, voltage-sensitivity and collinearity figures for updates during the run.
5. **Compute local sensitivities:** evaluate every accepted candidate across all operating cases, saving checkpoints as work units finish.
6. **Summarize and select parameters:** assess collinearity, apply the parameter-selection thresholds, and export the figures and summary tables.

### Parameter selection

The selection is made in two stages among the analysed parameters:

1. **Collinearity filter** (`psi_th`): when the median collinearity of a pair of parameters exceeds
   `psi_th`, the less sensitive of the two is removed.
2. **Sensitivity threshold** (`sensitivity_threshold`): among the remaining parameters, those whose
   median sensitivity magnitude reaches the threshold are selected.

!!! note
    The threshold selection is a deliberate simplification. It is **not** the method of Goshtasbi et al.,
    who choose the subset size from singular values and maximise a D-optimality criterion with a genetic
    search. The threshold selection judges each parameter on its own and cannot reveal redundancy among
    three or more parameters that no single pair shows.

## Running the Analysis

```bash
julia --project=. examples/run_extended_local_sensitivity_analysis.jl
```

The example script only sets the analysis parameters and displays the outcome:

```julia
using AlphaPEM.Parametrisation.ExtendedLocalSensitivityAnalysis

cfg = ExtendedLocalSensitivityConfig(
    n_candidates          = 2^9,
    sensitivity_threshold = 0.1,   # minimum median sensitivity to select a parameter
    collinearity_threshold = 0.707, # pairs highlighted in the collinearity figure
    psi_th                = 0.9,   # median collinearity above which a parameter is removed
)

outcome = run_extended_local_sensitivity_analysis(cfg)
outcome.analysis.selection.selected   # selected parameters
```

## Long runs: checkpoint and result archive

The simulations are by far the most expensive part, so two files protect and reuse them. Both are written
in `output_dir`.

**Checkpoint** (`extended_local_sensitivity_checkpoint.jld2`). Written every `checkpoint_frequency`
finished work units. If the run is interrupted (or reaches the walltime of a cluster), launch it again with
`resume = true`: finished work and the screening are restored and only the missing simulations are run.
The checkpoint is deleted once the final archive has been written. It is ignored, with a warning, if the
configuration of the simulations changed.

**Result archive** (`extended_local_sensitivity_result.jld2`, if `save_result = true`). Contains every
computed quantity: sensitivity blocks, validity masks, screening records, candidates, operating cases,
configuration and timings. It lets the analysis be repeated without simulating again:

```julia
cfg = ExtendedLocalSensitivityConfig(psi_th = 0.8, sensitivity_threshold = 0.05,
                                     output_dir = "results/extended_local_sensitivity")
outcome = reanalyze_extended_local_sensitivity(cfg)
```

Only the analysis settings (`sensitivity_threshold`, `collinearity_threshold`, `psi_th`, `output_dir`) may
differ from the original run. In the example script, set `LOAD_SAVED_RESULT = true`.

## Outputs

| File | Content |
|---|---|
| `convergence.pdf`, `voltage_sensitivity.pdf`, `collinearity.pdf` | Figures, overwritten at every update. |
| `collinear_pairs.csv` | Collinear pairs whose two parameters are above `sensitivity_threshold`. |
| `sensitivity_summary.csv` | One row per parameter: median sensitivity, removal by the collinearity filter, final selection. |
| `extended_local_sensitivity_summary.yaml` | Julia/AlphaPEM versions, seed, configuration hash, design sizes, simulation outcomes, thresholds, and the parameter selection (removed as collinear, below threshold, selected). |
| `extended_local_sensitivity_result.jld2` | Archive of the complete result. |
| `extended_local_sensitivity_checkpoint.jld2` | Resume state, present only while a run is unfinished. |

## Running on a cluster

```bash
qsub scripts/run_extended_local_sensitivity_analysis_cluster.sh
```

The job synchronises its results, including the checkpoint, back to the submission directory
periodically and on SIGTERM. After a walltime stop, set `RESUME = true` in the example script and
submit again.

## References

- A. Goshtasbi, J. Chen, J. R. Waldecker, S. Hirano, and T. Ersal, "Effective parameterization of pem
  fuel cell models—part i: Sensitivity analysis and parameter identifiability," *J. Electrochem. Soc.*,
  vol. 167, p. 044504, 2020. <https://doi.org/10.1149/1945-7111/ab7091>
- A. Goshtasbi et al., "Effective parameterization of pem fuel cell models—part ii: Robust parameter
  subset selection, robust optimal experimental design, and multi-step parameter identification
  algorithm," *J. Electrochem. Soc.*, vol. 167, p. 044505, 2020.
