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
#   - R + a C++ compiler must already be installed once on the login node in a
#     "r-env" conda environment (no sudo available on the cluster):
#       curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o miniforge.sh
#       bash miniforge.sh -b -p "$HOME/miniforge3"
#       export PATH="$HOME/miniforge3/bin:$PATH"
#       conda create -y -n r-env -c conda-forge \
#           r-base=4.3.2 cmake compilers \
#           libcurl openssl libxml2 libwebp libpng libtiff libjpeg-turbo \
#           freetype fontconfig harfbuzz fribidi
#     See README.md § Installation, step 4, for full instructions.
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


# **Headless-Node Workarounds:**
# `#PBS -V` copies the submitting shell's whole environment to the compute node, DISPLAY
# included (typically DISPLAY=localhost:10.0, left over from an X11-forwarded SSH session).
# That display does not exist on the node, so GLFW aborts with "X11: Failed to open
# display". This job never draws on screen, so drop the variable.
unset DISPLAY

# Julia runs an automatic, whole-environment precompilation pass every time a package is
# loaded from a cold cache — including on `using AlphaPEM` inside the example script. That
# implicit pass has no error tolerance and dies on GLMakie (OpenGL/X11 unavailable here),
# which is what aborts the job before any AlphaPEM code runs. Disable it; precompilation
# is driven explicitly below instead.
export JULIA_PKG_PRECOMPILE_AUTO=0


# **R / IRD Environment Setup:**
# STEP 3 of run_parameter_validity.jl (PRIM/MaxBox restricted-region analysis)
# shells out to Rscript. Set this up upfront rather than failing after paying
# for the full batch simulation (STEP 1-2, by far the most expensive part of the job).
# This must also happen BEFORE the Julia setup below: AlphaPEM loads RCall.jl, which
# postpones its own precompilation when no R installation is on PATH.
echo "================================================================================"
echo "              R / IRD Environment Setup"
echo "================================================================================"

# Activate the "r-env" conda environment (R + C++ compiler), created once on the
# login node as described in README.md § Installation, step 4.
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate r-env

if ! command -v Rscript &> /dev/null; then
    echo "[ERROR] Rscript not found in PATH after 'conda activate r-env'."
    echo "[ERROR] Create the environment once on the login node (README.md § Installation, step 4):"
    echo "[ERROR]   conda create -y -n r-env -c conda-forge r-base=4.3.2 cmake compilers \\"
    echo "[ERROR]       libcurl openssl libxml2 libwebp libpng libtiff libjpeg-turbo \\"
    echo "[ERROR]       freetype fontconfig harfbuzz fribidi"
    exit 1
fi
echo "[INFO] Using: $(Rscript --version 2>&1)"

# Clone the IRD package if not already present (README.md § Installation, step 4.b).
IRD_PKG_DIR="$PBS_TMPDIR/$PROJECT_NAME/external/IRD_method_2023/irdpackage"
if [ ! -d "$IRD_PKG_DIR" ]; then
    echo "[INFO] IRD package not found. Cloning it..."
    git clone https://github.com/slds-lmu/supplementary_2023_ird.git external/IRD_method_2023
fi
echo "[INFO] IRD package directory found: $IRD_PKG_DIR"

# Install the required R packages if not already present (README.md § Installation, step 4.c).
echo "[INFO] Checking required R packages (devtools, mlr3, mlr3learners, mlr3pipelines, iml, ranger, yaml, jsonlite, data.table, optparse)..."
Rscript -e '
required <- c("devtools", "mlr3", "mlr3learners", "mlr3pipelines", "iml", "ranger",
              "yaml", "jsonlite", "data.table", "optparse")
missing  <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
    cat("[INFO] Missing R packages:", paste(missing, collapse = ", "), "-- installing them.\n")
    quit(status = 1)
} else {
    cat("[INFO] All required R packages are available.\n")
}
'
if [ $? -ne 0 ]; then
    Rscript src/alphapem/parametrisation/validity/R/install_r_packages.R
    R_INSTALL_STATUS=$?
    if [ $R_INSTALL_STATUS -ne 0 ]; then
        echo "[ERROR] R package installation failed. Aborting job."
        exit $R_INSTALL_STATUS
    fi
fi
echo ""

# **Julia Environment Configuration:**
# Runs after the R setup above so that RCall.jl finds R and precompiles for real.
echo "================================================================================"
echo "              Julia Environment Configuration"
echo "================================================================================"

# Use Julia directly from PATH
julia_version=$(julia --version)
echo "[INFO] Using: $julia_version"

# Refresh the package registries before instantiating. The committed Manifest.toml pins
# versions (e.g. PureKLU 1.1.0) that a stale General clone on the cluster does not know
# about, which fails with "Unsatisfiable requirements detected".
echo "[INFO] Updating Julia registries..."
julia --project -e 'using Pkg; Pkg.Registry.update()'

# Instantiate strictly from the committed Manifest.toml.
# Do NOT call Pkg.resolve() here: it re-resolves the whole dependency graph against the
# registry and discards the known-good pinned versions shipped with the project.
echo "[INFO] Instantiating Julia project..."
julia --project -e 'using Pkg; Pkg.instantiate()'
if [ $? -ne 0 ]; then
    echo "[ERROR] Instantiation failed. Aborting job."
    exit 1
fi
echo "[INFO] Project instantiated successfully"

# Precompile every direct dependency EXCEPT GLMakie, which requires OpenGL/X11 and can
# never precompile on a compute node.
# Note: `Pkg.precompile(strict=false)` does NOT help here — `strict` only downgrades
# failures of *indirect* dependencies to warnings, whereas GLMakie is a direct dependency
# and still raises. Excluding it by name keeps every genuine failure fatal.
echo "[INFO] Precompiling packages (headless-safe, GLMakie excluded)..."
julia --project -e '
    using Pkg
    headless_unsafe = ["GLMakie"]
    pkgs = sort([n for n in keys(Pkg.project().dependencies) if !(n in headless_unsafe)])
    Pkg.precompile(pkgs)
'
if [ $? -ne 0 ]; then
    echo "[ERROR] Precompilation failed. Aborting job."
    exit 1
fi

# Warm AlphaPEM's own cache and check the environment actually loads (RCall included),
# before paying for the full LHS batch simulation.
echo "[INFO] Warming AlphaPEM cache..."
julia --project -e 'using AlphaPEM'
if [ $? -ne 0 ]; then
    echo "[ERROR] AlphaPEM failed to load. Aborting job."
    exit 1
fi
echo "[INFO] Precompilation completed"
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
