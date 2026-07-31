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
# Pin HOME/USER to the values PBS captured at submission time.
#
# WORKAROUND for a cluster misconfiguration: two accounts share UID/GID 1149
# ("getent passwd 1149" returns both rgass and bthomont). PBS's MoM rebuilds HOME and
# USER on the compute node via getpwuid(), and which of the two entries wins depends on
# the node's resolution order — so HOME randomly lands on the *other* account's home
# directory (and "id" warns "cannot find name for GID 1149" because that GID has no
# group entry at all). PBS_O_HOME/PBS_O_LOGNAME were resolved on the login node, where
# the ordering is correct, so they are the trustworthy values.
#
# This matters far beyond cosmetics: the Julia depot (~/.julia), the pixi/CondaPkg cache,
# $HOME/miniforge3 (the "r-env" conda environment used below) and /gpfs/scratch/$USER all
# derive from these variables.
if [ -n "$PBS_O_HOME" ]; then
    export HOME="$PBS_O_HOME"
fi
if [ -n "$PBS_O_LOGNAME" ]; then
    export USER="$PBS_O_LOGNAME"
    export LOGNAME="$PBS_O_LOGNAME"
fi
echo "[INFO] HOME: $HOME"
echo "[INFO] USER: $USER"

# Fail fast rather than 150 lines further down if HOME is still unusable.
if [ ! -d "$HOME" ] || [ ! -w "$HOME" ]; then
    echo "[ERROR] HOME ('$HOME') is not a writable directory. Aborting job."
    exit 1
fi

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

# Suppress the implicit precompilation pass that Pkg runs after each of its operations, so
# the Julia setup below stays explicit and ordered. Note this does NOT cover the pass Base
# runs at load time — see the GLMakie discussion in "Julia Environment Configuration".
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

# Locate the Miniforge installation holding the "r-env" environment (R + C++ compiler),
# created once on the login node as described in README.md § Installation, step 4.
# Several candidates are probed instead of hard-coding "$HOME/miniforge3": see the UID
# collision documented in "Environment Preparation" above. Set ALPHAPEM_CONDA_ROOT to
# override if Miniforge lives elsewhere.
CONDA_ROOT=""
for candidate in "$ALPHAPEM_CONDA_ROOT" "$HOME/miniforge3" "$PBS_O_HOME/miniforge3" \
                 "$HOME/miniconda3" "$PBS_O_HOME/miniconda3"; do
    if [ -n "$candidate" ] && [ -f "$candidate/etc/profile.d/conda.sh" ]; then
        CONDA_ROOT="$candidate"
        break
    fi
done
if [ -z "$CONDA_ROOT" ] && command -v conda &> /dev/null; then
    CONDA_ROOT=$(conda info --base 2>/dev/null)
fi

if [ -z "$CONDA_ROOT" ] || [ ! -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]; then
    echo "[ERROR] No conda installation found (looked for etc/profile.d/conda.sh under"
    echo "[ERROR]   \$ALPHAPEM_CONDA_ROOT, \$HOME/miniforge3, \$PBS_O_HOME/miniforge3, ...)."
    echo "[ERROR] Install Miniforge once on the login node (README.md § Installation, step 4),"
    echo "[ERROR] or point ALPHAPEM_CONDA_ROOT at an existing installation."
    exit 1
fi
echo "[INFO] Conda base: $CONDA_ROOT"

source "$CONDA_ROOT/etc/profile.d/conda.sh"
if ! conda activate r-env; then
    echo "[ERROR] 'conda activate r-env' failed under $CONDA_ROOT."
    echo "[ERROR] Available environments:"
    conda env list 2>&1 | sed 's/^/[ERROR]   /'
    exit 1
fi

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

# Drop GLMakie from THIS COPY of the project (the scratch copy is disposable, and
# Project.toml/Manifest.toml are excluded from the results rsync at the end of the job).
#
# This is the only thing that actually works. GLMakie needs OpenGL/X11 and can never
# precompile on a compute node, and in Julia >= 1.12 load-time precompilation is done by
# Base, not Pkg (loading.jl -> Precompilation.precompilepkgs(...; _from_loading=true)),
# which:
#   - ignores JULIA_PKG_PRECOMPILE_AUTO (only Pkg.should_autoprecompile() reads it, for
#     Pkg operations such as the Pkg.rm below — hence the export is still worth keeping);
#   - walks the WHOLE environment graph whatever package was requested;
#   - raises on any failure of a package in `project_deps`, i.e. anything listed in
#     Project.toml [deps] (precompilation.jl: `if strict || (dep in project_deps)`).
# So `Pkg.precompile(strict=false)` and excluding GLMakie from an explicit precompile
# call both fail to protect the eventual `using AlphaPEM`: as long as GLMakie sits in
# [deps], every load from this environment tries it and dies.
#
# AlphaPEM never imports GLMakie (src/alphapem/application/run_simulation_modules.jl
# resolves it lazily through Base.loaded_modules and falls back to CairoMakie), so
# removing it changes nothing for a headless batch run.
echo "[INFO] Removing GLMakie from the scratch project copy (headless node, no OpenGL/X11)..."
julia --project -e 'using Pkg; haskey(Pkg.project().dependencies, "GLMakie") && Pkg.rm("GLMakie")'
if [ $? -ne 0 ]; then
    echo "[ERROR] Could not remove GLMakie from the project copy. Aborting job."
    exit 1
fi

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

# Rebuild RCall against the R of "r-env".
# RCall 0.14 resolves R at *build* time and bakes the result into deps/deps.jl; it does
# NOT read ENV["R_HOME"] when loaded. Having R on PATH is therefore not enough: if RCall
# was ever built without R, it keeps `Rhome = ""` and declares `__precompile__(false)`
# ("No R installation found by RCall.jl ... Importing RCall will fail"), which propagates
# to AlphaPEM since src/.../validity/ird_interface.jl does `using RCall`.
export R_HOME=$(R RHOME)
echo "[INFO] R_HOME: $R_HOME"
echo "[INFO] Building RCall against r-env..."
julia --project -e 'using Pkg; Pkg.build("RCall")'
if [ $? -ne 0 ]; then
    echo "[ERROR] Pkg.build(\"RCall\") failed. Aborting job."
    exit 1
fi

# Precompile the whole environment. GLMakie is gone, so this can now be strict: any
# failure left is a real one and must stop the job.
echo "[INFO] Precompiling packages..."
julia --project -e 'using Pkg; Pkg.precompile()'
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

# Copy results back to original directory.
# Excluded: .git (read-only objects cause permission errors) and the root Project.toml /
# Manifest.toml, which this job mutated on purpose to drop GLMakie — that edit is local to
# the scratch copy and must not leak back into the repository.
echo "[INFO] Copying results to original directory..."
cd "$PBS_TMPDIR/$PROJECT_NAME"
rsync -a --exclude='.git' --exclude='/Project.toml' --exclude='/Manifest.toml' . "$PROJECT_ROOT/"
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
