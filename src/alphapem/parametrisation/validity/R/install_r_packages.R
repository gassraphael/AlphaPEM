#!/usr/bin/env Rscript
# install_r_packages.R
#
# One-time setup script: installs all R packages required by AlphaPEM's
# PRIM-based valid parameter region analysis (run_parameter_validity.jl).
#
# This file is the R equivalent of `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
# Specific versions are pinned to ensure compatibility with the irdpackage
# (external/IRD_method_2023/), which uses the paradox 0.x API
# (ParamDbl$new, ParamInt$new, ParamFct$new, ParamSet$new, SamplerUnif).
# paradox 1.0+ dropped these classes; mlr3 0.19+ requires paradox >= 1.0.1.
# Therefore mlr3 must be pinned to 0.18.0 (last version accepting paradox 0.x).
#
# PREREQUISITES — Before running this script, install the required system libraries:
#
#   Linux (Debian/Ubuntu):
#     sudo apt update && sudo apt install -y \
#         r-base build-essential cmake \
#         libcurl4-openssl-dev libssl-dev libxml2-dev \
#         libpng-dev libtiff5-dev libjpeg-dev \
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

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- Helper: install a specific version from CRAN archive -------------------
install_version_archive <- function(pkg, version) {
  installed <- tryCatch(packageVersion(pkg), error = function(e) NULL)
  if (!is.null(installed) && installed == version) {
    message(sprintf("  %s %s already installed — skipping.", pkg, version))
    return(invisible(NULL))
  }
  if (!is.null(installed)) {
    message(sprintf("  Downgrading %s %s → %s …", pkg, installed, version))
  } else {
    message(sprintf("  Installing %s %s …", pkg, version))
  }
  url <- sprintf(
    "https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz",
    pkg, pkg, version
  )
  install.packages(url, repos = NULL, type = "source")
}

# ---- Pinned versions (paradox 0.x ecosystem) --------------------------------
# These versions form a consistent, tested set compatible with the irdpackage.
pinned <- list(
  # mlr3misc first (required by paradox and mlr3)
  # 0.14.0 is the minimum actually enforced at load-time by mlr3 0.18.0
  mlr3misc     = "0.14.0",
  # mlr3measures (required by mlr3)
  mlr3measures = "0.5.0",
  # paradox 0.11.1: last release with ParamDbl$new / ParamSet$new / SamplerUnif
  paradox      = "0.11.1",
  # mlr3 0.18.0: last release that accepts paradox >= 0.10.0 (not 1.0+)
  mlr3         = "0.18.0",
  # mlr3learners 0.6.0: last release requiring only mlr3 >= 0.17.1
  mlr3learners = "0.6.0"
)

# ---- Latest-version packages (no paradox constraints) -----------------------
latest <- c(
  "R6",            # R6 classes used by irdpackage
  "checkmate",     # fast argument checking
  "optparse",      # CLI argument parsing
  "data.table",    # fast CSV reading
  "iml",           # interpretable ML (Predictor, used as bridge to irdpackage)
  "ranger",        # fast Random Forest (requires C++ compiler)
  "yaml",          # YAML read/write
  "jsonlite",      # JSON read/write
  "RhpcBLASctl"   # BLAS thread control (required by mlr3 0.18.0)
)

# ---- Install pinned versions first ------------------------------------------
message("\n=== Installing pinned versions (paradox 0.x ecosystem) ===")
for (pkg in names(pinned)) {
  install_version_archive(pkg, pinned[[pkg]])
}

# ---- Install / update latest-version packages -------------------------------
message("\n=== Installing/updating remaining packages ===")
already  <- latest[vapply(latest, requireNamespace, logical(1), quietly = TRUE)]
to_inst  <- setdiff(latest, already)
if (length(already)) message("Already installed: ", paste(already, collapse = ", "))
if (length(to_inst)) {
  message("Installing: ", paste(to_inst, collapse = ", "))
  install.packages(to_inst)
}

# ---- Verify all packages load -----------------------------------------------
message("\n=== Verifying installation ===")
all_pkgs <- c(names(pinned), latest)
missing  <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Installation failed for: ", paste(missing, collapse = ", "),
    "\nCheck the error messages above and retry with sudo / as Administrator."
  )
}

# Final version report
for (pkg in names(pinned)) {
  v <- as.character(packageVersion(pkg))
  ok <- v == pinned[[pkg]]
  message(sprintf("  %-14s %s  %s", pkg, v, if (ok) "[OK]" else paste0("[MISMATCH — expected ", pinned[[pkg]], "]")))
}
for (pkg in latest) {
  message(sprintf("  %-14s %s", pkg, as.character(packageVersion(pkg))))
}

message("\nAll R packages installed successfully.")
message("You can now run:")
message("  RUN_PRIM=true julia --project=. --threads=auto examples/run_parameter_validity.jl")

