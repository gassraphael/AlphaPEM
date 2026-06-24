#!/bin/bash

#################################################################################
# AlphaPEM - HPC Cluster Calibration Script
#################################################################################
#
# Description:
#   This PBS script launches automated parameter calibration of AlphaPEM on an
#   HPC (High Performance Computing) cluster. It configures the environment,
#   installs dependencies, and executes the genetic algorithm calibration module.
#
# Usage:
#   qsub scripts/run_calibration_cluster.sh
#
# Note:
#   - Adapt PBS parameters according to your cluster (ncpus, mem, walltime)
#   - Script copies entire project to temporary workspace for execution
#   - Results are automatically copied back to original directory
#
# Author: Raphaël Gass
# License: GPLv3
#################################################################################


# **PBS (Portable Batch System) Configuration:**
# select: number of nodes to allocate
# ncpus: number of CPUs per node (adapt according to your cluster)
# mem: RAM memory allocated per node
#PBS -l select=1:ncpus=32:mem=1gb

# Maximum execution time on cluster (format HH:MM:SS)
#PBS -l walltime=90:00:00

# Standard output and error files
#PBS -o out_calibration
#PBS -e err_calibration

# Job name displayed in queue
#PBS -N AlphaPEM_calibration

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
echo "              AlphaPEM Parameter Calibration - Job Information"
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
# install modern Python and PyGAD (both are pure-Python at runtime and do not call glibc 2.28 APIs).
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

# Check that pygad is available (PythonCall/CondaPkg manage the Python environment automatically;
# pygad is declared in CondaPkg.toml and installed on first use — no manual setup needed).
echo "[INFO] Checking Python/pygad availability..."
julia --project << 'JULIA_EOF'
    using PythonCall
    try
        pygad = pyimport("pygad")
        pyver = pyconvert(String, pyimport("sys").version.split(" ")[0])
        println("[INFO] pygad OK (Python: $pyver)")
    catch e
        println(stderr, "[ERROR] pygad not found. Run once on the login node:")
        println(stderr, "  julia --project=. -e 'using CondaPkg; CondaPkg.resolve()'")
        exit(1)
    end
JULIA_EOF
PYGAD_STATUS=$?
if [ $PYGAD_STATUS -ne 0 ]; then
    echo "[ERROR] Python/pygad check failed. Aborting job."
    exit $PYGAD_STATUS
fi
echo ""

# **Calibration Execution:**
echo "================================================================================"
echo "              Launching AlphaPEM Calibration"
echo "================================================================================"
echo "[INFO] Script: examples/run_calibration.jl"
echo "[INFO] Execution start: $(date)"
echo "--------------------------------------------------------------------------------"
echo ""

# Launch calibration script
# Count CPUs from PBS node file and explicitly set threads
NCPUS=$(wc -l < $PBS_NODEFILE)
julia --threads=$NCPUS --project examples/run_calibration.jl

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

