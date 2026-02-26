# AlphaPEM Scripts
This directory contains utility scripts for running AlphaPEM, particularly on HPC (High Performance Computing) clusters.

## `run_calibration_cluster.sh`
PBS (Portable Batch System) script that submits an AlphaPEM calibration job on a cluster. It prepares a temporary workspace, configures Python dependencies, runs the calibration module, and copies results back to the original project directory.

## License
This script is part of AlphaPEM (GPLv3). The calibration module uses PyGAD (BSD-3-Clause); see `src/alphapem/parametrisation/LICENSE-BSD-3-CLAUSE` for details.
