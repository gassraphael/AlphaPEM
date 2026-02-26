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

# Copy entire project to temporary directory
echo "[INFO] Copying project to temporary workspace..."
cp -r "$PROJECT_ROOT" "$PBS_TMPDIR"
echo "[INFO] Copy completed."
echo ""

# Move to project copy
cd "$PBS_TMPDIR/$PROJECT_NAME"


# **Python Environment Configuration:**
echo "================================================================================"
echo "              Python Environment Configuration"
echo "================================================================================"

# Load cluster modules
# Note: Adapt these versions to match your cluster's available modules
module purge
module load mpi/intel/2019_update2
module load tools/python/3.7.2
echo "[INFO] Modules loaded: MPI Intel 2019, Python 3.7.2"

# Clean PYTHONPATH to avoid conflicts
unset PYTHONPATH
echo "[INFO] PYTHONPATH cleared"

# Create Python virtual environment
echo "[INFO] Creating virtual environment..."
python3 -m venv env
source env/bin/activate
echo "[INFO] Virtual environment activated"

# Install dependencies
echo "[INFO] Installing Python dependencies..."
python3 -m pip install --upgrade --force-reinstall --no-cache-dir pip setuptools wheel
python3 -m pip install --upgrade --force-reinstall --no-cache-dir -e .
echo "[INFO] Dependencies installed successfully"
echo ""

# Add project root to PYTHONPATH for imports
export PYTHONPATH="$PBS_TMPDIR/$PROJECT_NAME:$PYTHONPATH"
echo "[INFO] PYTHONPATH configured: $PYTHONPATH"
echo ""


# **Calibration Execution:**
echo "================================================================================"
echo "              Launching AlphaPEM Calibration"
echo "================================================================================"
echo "[INFO] Module: alphapem.parametrisation.calibration"
echo "[INFO] Execution start: $(date)"
echo "--------------------------------------------------------------------------------"
echo ""

# Launch calibration module with MPI
# Note: Use -np 1 for serial execution, or adjust according to your needs
mpiexec -np 1 python3 -m alphapem.parametrisation.calibration

echo ""
echo "--------------------------------------------------------------------------------"
echo "[INFO] Execution end: $(date)"
echo "================================================================================"
echo ""


# **Results Backup:**
echo "================================================================================"
echo "              Results Backup"
echo "================================================================================"

# Copy results back to original directory
echo "[INFO] Copying results to original directory..."
cd "$PBS_TMPDIR/$PROJECT_NAME"
cp -ru . "$PROJECT_ROOT"
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

