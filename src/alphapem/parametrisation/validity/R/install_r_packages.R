#!/usr/bin/env Rscript
# install_r_packages.R
#
# One-time setup script: installs all R packages required by AlphaPEM's
# PRIM-based valid parameter region analysis (run_parameter_validity.jl).
#
# This file is the R equivalent of `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
# All required packages are declared here so that the installation is reproducible
# and independent of the calling process's permissions.
#
# PREREQUISITES — Before running this script, install the required system libraries:
#
#   Linux (Debian/Ubuntu):
#     sudo apt update && sudo apt install -y \
#         r-base build-essential cmake \
#         libcurl4-openssl-dev libssl-dev libxml2-dev \
#         libwebp-dev libpng-dev libtiff5-dev libjpeg-dev \
#         libfreetype6-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
#
#   macOS: brew install r   (system headers included via Xcode Command Line Tools)
#   Windows: install Rtools from https://cran.r-project.org/bin/windows/Rtools/
#
# Usage (run once, from the AlphaPEM root directory):
#
#   Linux / macOS:
#     sudo Rscript src/alphapem/parametrisation/validity/R/install_r_packages.R
#
#   Windows (run as Administrator):
#     Rscript src\alphapem\parametrisation\validity\R\install_r_packages.R

required <- c(
  "R6",             # R6 classes used by irdpackage (Prim, MaxBox, PostProcessing)
  "checkmate",      # fast argument checking (irdpackage import)
  "paradox",        # parameter sets used by irdpackage (ParamDbl, ParamInt, …)
  "optparse",       # CLI argument parsing
  "data.table",     # fast CSV reading
  "mlr3",           # machine-learning framework
  "mlr3learners",   # learner implementations (Random Forest, …)
  "mlr3pipelines",  # pre-processing pipelines
  "iml",            # interpretable ML (Predictor, PRIM, MaxBox wrappers)
  "ranger",         # fast Random Forest implementation (requires C++ compiler)
  "yaml",           # YAML read/write
  "jsonlite"        # JSON read/write
)

already    <- required[vapply(required, requireNamespace, logical(1), quietly = TRUE)]
to_install <- setdiff(required, already)

if (length(already)) {
  message("Already installed: ", paste(already, collapse = ", "))
}

if (!length(to_install)) {
  message("All required R packages are already installed. Nothing to do.")
  quit(status = 0)
}

message("Installing: ", paste(to_install, collapse = ", "))
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(to_install)

# Verify that every package loaded successfully
still_missing <- to_install[!vapply(to_install, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("Installation failed for: ", paste(still_missing, collapse = ", "),
       "\nCheck the error messages above and retry with sudo / as Administrator.")
}

message("\nAll R packages installed successfully.")
message("You can now run:  RUN_PRIM=true julia --project=. --threads=auto examples/run_parameter_validity.jl")

