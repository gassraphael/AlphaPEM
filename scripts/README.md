# AlphaPEM Scripts
This directory contains utility scripts for running AlphaPEM, particularly on HPC (High Performance Computing) clusters.

## `run_calibration_cluster.sh`
PBS (Portable Batch System) script that submits an AlphaPEM calibration job on a cluster. It prepares a temporary workspace, configures Python dependencies, runs the calibration module, and copies results back to the original project directory.

## `run_extended_local_sensitivity_analysis_cluster.sh`
PBS script that submits the extended local sensitivity analysis (`examples/run_extended_local_sensitivity_analysis.jl`). It prepares a temporary workspace, runs the analysis, and synchronises the results, including the resume checkpoint, back to the project directory periodically and on SIGTERM. After a walltime stop, set `RESUME = true` in the example script and submit again.

## License
This script is part of AlphaPEM (BSD-3-Clause). The calibration module uses PyGAD (BSD-3-Clause); see `src/alphapem/parametrisation/LICENSE-BSD-3-CLAUSE` for details.
