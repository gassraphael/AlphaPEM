#!/bin/bash

#################################################################################
# AlphaPEM - HPC Cluster Parameter-Validity Script
#################################################################################
#
# Description:
#   This PBS script launches the AlphaPEM parameter-validity-region analysis
#   (examples/run_parameter_validity.jl) on an HPC (High Performance Computing)
#   cluster. It configures the environment, installs dependencies, and executes
#   the full pipeline: LHS sampling, batch simulation & classification, and
#   IRD-based restricted-region analysis (PRIM/MaxBox, via R).
#
# Usage:
#   qsub scripts/run_parameter_validity_cluster.sh
#
# Note:
#   - Adapt PBS parameters according to your cluster (ncpus, mem, walltime)
#   - Script copies entire project to temporary workspace for execution
#   - Results are automatically copied back to original directory
#   - The IRD package must already be cloned once on the login node:
#       git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/IRD_method_2023
#     See README.md § Installation for full instructions (R + required packages).
#
# Author: Raphaël Gass
# License: GPLv3
#################################################################################


# **PBS (Portable Batch System) Configuration:**
# select: number of nodes to allocate
# ncpus: number of CPUs per node (adapt according to your cluster; the batch
#        simulation loop scales with Threads.nthreads(), set below from ncpus)
#PBS -l select=1:ncpus=32:mem=8gb

# Maximum execution time on cluster (format HH:MM:SS)
#PBS -l walltime=90:00:00

# Standard output and error files
#PBS -o out_parameter_validity
#PBS -e err_parameter_validity

# Job name displayed in queue
#PBS -N AlphaPEM_parameter_validity

# Export environment variables to job
#PBS -V

# Queue name
#PBS -q longq

# Email address for notifications
#PBS -M raphael.gass@univ-reunion.fr

# Notification types: (b)eginning, (e)nd, (a)bortion
# Email sent for each of these states
#PBS -m bea


# **Job Information Display (for debugging):**
echo "================================================================================"
echo "         AlphaPEM Parameter Validity Region Analysis - Job Information"
echo "================================================================================"
echo -n 'Job is running on node: '; cat $PBS_NODEFILE
echo "--------------------------------------------------------------------------------"
echo "PBS: qsub launched from:       $PBS_O_HOST"
echo "PBS: originating queue:        $PBS_O_QUEUE"
echo "PBS: executing queue:          $PBS_QUEUE"
echo "PBS: working directory:        $PBS_O_WORKDIR"
echo "PBS: execution mode:           $PBS_ENVIRONMENT"
echo "PBS: job identifier:           $PBS_JOBID"
echo "PBS: job name:                 $PBS_JOBNAME"
echo "PBS: node file:                $PBS_NODEFILE"
echo "PBS: home directory:           $PBS_O_HOME"
echo "PBS: PATH:                     $PBS_O_PATH"
echo "================================================================================"
echo ""


# **Environment Preparation:**
# Change to submission directory
cd "$PBS_O_WORKDIR"

# Identify project structure
PROJECT_ROOT=$(pwd)                           # AlphaPEM project root
PROJECT_NAME=$(basename "$PROJECT_ROOT")      # Project name (AlphaPEM)

echo "[INFO] Project: $PROJECT_NAME"
echo "[INFO] Root: $PROJECT_ROOT"
echo ""

# Create temporary working directory
# (local home directory is often size-limited on clusters)
export PBS_TMPDIR=/gpfs/scratch/$USER/$PBS_JOBID
mkdir -p "$PBS_TMPDIR"
echo "[INFO] Temporary workspace created: $PBS_TMPDIR"

# Copy entire project to temporary directory, excluding .git (read-only objects cause cp errors)
echo "[INFO] Copying project to temporary workspace..."
rsync -a --exclude='.git' "$PROJECT_ROOT/" "$PBS_TMPDIR/$PROJECT_NAME/"
echo "[INFO] Copy completed."
echo ""

# Move to project copy
cd "$PBS_TMPDIR/$PROJECT_NAME"


# **Conda/glibc Compatibility Workaround:**
# This cluster runs glibc 2.17 (CentOS 7), but conda-forge packages now require glibc >= 2.28.
# Setting CONDA_OVERRIDE_GLIBC=2.28 bypasses the virtual package version check so pixi can
# install modern Python (CondaPkg.toml dependency, unrelated to the R/IRD step below).
export CONDA_OVERRIDE_GLIBC=2.28


# **Julia Environment Configuration:**
echo "================================================================================"
echo "              Julia Environment Configuration"
echo "================================================================================"

# Use Julia directly from PATH
julia_version=$(julia --version)
echo "[INFO] Using: $julia_version"

# Setup Julia environment
echo "[INFO] Instantiating Julia project..."
julia --project -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
echo "[INFO] Project instantiated successfully"

# Precompile packages, ignoring display-dependent packages (GLMakie/WGLMakie)
# that fail on headless servers without an X11 display.
echo "[INFO] Precompiling packages (non-strict, headless-safe)..."
julia --project -e 'using Pkg; Pkg.precompile(strict=false)'
echo "[INFO] Precompilation completed"
echo ""

# **R / IRD Environment Check:**
# STEP 3 of run_parameter_validity.jl (PRIM/MaxBox restricted-region analysis)
# shells out to Rscript. Check upfront rather than failing after paying for the
# full batch simulation (STEP 1-2, by far the most expensive part of the job).
echo "================================================================================"
echo "              R / IRD Environment Check"
echo "================================================================================"

IRD_PKG_DIR="$PBS_TMPDIR/$PROJECT_NAME/external/IRD_method_2023/irdpackage"

if ! command -v Rscript &> /dev/null; then
    echo "[ERROR] Rscript not found in PATH. Install R before submitting this job."
    exit 1
fi
echo "[INFO] Using: $(Rscript --version 2>&1)"

if [ ! -d "$IRD_PKG_DIR" ]; then
    echo "[ERROR] IRD package directory not found: $IRD_PKG_DIR"
    echo "[ERROR] Run once on the login node (inside the repository, not the scratch copy):"
    echo "[ERROR]   git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/IRD_method_2023"
    exit 1
fi
echo "[INFO] IRD package directory found: $IRD_PKG_DIR"

echo "[INFO] Checking required R packages (devtools, mlr3, mlr3learners, mlr3pipelines, iml, ranger, yaml, jsonlite, data.table, optparse)..."
Rscript -e '
required <- c("devtools", "mlr3", "mlr3learners", "mlr3pipelines", "iml", "ranger",
              "yaml", "jsonlite", "data.table", "optparse")
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    cat("[ERROR] Missing R packages:", paste(missing, collapse = ", "), "\n")
    quit(status = 1)
} else {
    cat("[INFO] All required R packages are available.\n")
}
'
R_STATUS=$?
if [ $R_STATUS -ne 0 ]; then
    echo "[ERROR] R/IRD check failed. Aborting job."
    exit $R_STATUS
fi
echo ""

# **Parameter-Validity Analysis Execution:**
echo "================================================================================"
echo "         Launching AlphaPEM Parameter Validity Region Analysis"
echo "================================================================================"
echo "[INFO] Script: examples/run_parameter_validity.jl"
echo "[INFO] Execution start: $(date)"
echo "--------------------------------------------------------------------------------"
echo ""

# Launch parameter-validity script
# Count CPUs from PBS node file; the script itself also auto-restarts with the
# right thread count if launched with too few threads (see PARALLEL/N_THREADS
# constants at the top of run_parameter_validity.jl), so this is a safety net.
NCPUS=$(wc -l < $PBS_NODEFILE)
julia --threads=$NCPUS --project examples/run_parameter_validity.jl

echo ""
echo "--------------------------------------------------------------------------------"
echo "[INFO] Execution end: $(date)"
echo "================================================================================"
echo ""


# **Results Backup:**
echo "================================================================================"
echo "              Results Backup"
echo "================================================================================"

# Copy results back to original directory (exclude .git to avoid permission errors)
echo "[INFO] Copying results to original directory..."
cd "$PBS_TMPDIR/$PROJECT_NAME"
rsync -a --exclude='.git' . "$PROJECT_ROOT/"
echo "[INFO] Results copied successfully"

# Return to original directory
cd "$PBS_O_WORKDIR"

# Clean up temporary directory
echo "[INFO] Cleaning temporary workspace..."
rm -rf "$PBS_TMPDIR"
echo "[INFO] Cleanup completed"

echo "================================================================================"
echo "              Job completed successfully"
echo "================================================================================"
