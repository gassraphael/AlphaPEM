#!/usr/bin/env Rscript
# install_r_packages.R
#
# One-time setup script: installs all R packages required by AlphaPEM's
# PRIM-based valid parameter region analysis (run_parameter_validity.jl).
#
# This file is the R equivalent of `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.
#
# WHY A FROZEN CRAN SNAPSHOT INSTEAD OF "LATEST":
# irdpackage (external/IRD_method_2023/) is unmaintained and hard-pinned to the
# paradox 0.x API (ParamDbl$new, ParamInt$new, ParamFct$new, ParamSet$new,
# SamplerUnif). paradox 1.0+ (June 2024) dropped those classes, and mlr3 0.19+
# requires paradox >= 1.0.1 — so the whole mlr3/paradox ecosystem must resolve
# to versions that predate that break. Pinning a handful of package *versions*
# by hand (as this script used to do) doesn't stay correct: every other
# package installed as "latest" keeps drifting forward, and a future
# transitive dependency bump can silently break the pinned packages again.
#
# Instead, every package below is installed from a single frozen CRAN mirror
# snapshot (Posit Package Manager, https://packagemanager.posit.co), dated
# 2024-03-15 — three days after mlr3 0.18.0 (the last release accepting
# paradox 0.x) and three months before paradox 1.0.0. Because the snapshot is
# frozen in time, install.packages() always resolves the exact same package
# versions (mlr3 0.18.0, paradox 0.11.1, and every dependency as it existed
# that day), no matter what is released on CRAN afterwards. There is nothing
# left to re-break as time passes: this script's behavior is fixed forever
# unless SNAPSHOT_DATE below is deliberately changed.
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

# Frozen CRAN snapshot date — see the header comment for why this must not
# simply track "latest". Change this only after re-validating irdpackage
# against the new date's paradox/mlr3 versions.
SNAPSHOT_DATE <- "2024-03-15"
options(repos = c(CRAN = sprintf("https://packagemanager.posit.co/cran/%s", SNAPSHOT_DATE)))

# Force serial installation (Ncpus = 1): parallel install.packages() runs
# concurrent `R CMD INSTALL` processes against the same library folder,
# which can deadlock on file locks. Safe to parallelize later, at runtime
# (see run_parameter_validity.jl), just not during package installation.
options(Ncpus = 1)

# ---- Packages required by irdpackage and run_parameter_validity.jl ----------
# Expected versions resolved from the SNAPSHOT_DATE mirror above (recorded
# here only so the final report can flag an unexpected mismatch — they are
# not fetched individually, install.packages() resolves them all at once).
expected <- list(
  R6            = "2.5.1",
  checkmate     = "2.3.1",
  optparse      = "1.7.4",
  data.table    = "1.15.2",
  iml           = "0.11.1",
  ranger        = "0.16.0",
  yaml          = "2.3.8",
  jsonlite      = "1.8.8",
  RhpcBLASctl   = "0.23.42", # packageVersion() normalizes "0.23-42" to "0.23.42"
  mlr3misc      = "0.14.0",
  mlr3measures  = "0.5.0",
  paradox       = "0.11.1", # last release with ParamDbl$new / ParamSet$new / SamplerUnif
  mlr3          = "0.18.0", # last release that accepts paradox 0.x (0.19+ requires paradox >= 1.0.1)
  mlr3learners  = "0.6.0"
)
pkgs <- names(expected)

# ---- Install everything from the frozen snapshot in one resolution pass -----
# Letting install.packages() resolve the whole dependency graph at once (rather
# than installing package-by-package) is what guarantees a consistent,
# mutually-compatible set — the snapshot pins every transitive dependency too.
message("\n=== Installing R packages from CRAN snapshot ", SNAPSHOT_DATE, " ===")
install.packages(pkgs)

# ---- Verify all packages load -----------------------------------------------
message("\n=== Verifying installation ===")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Installation failed for: ", paste(missing, collapse = ", "),
    "\nCheck the error messages above and retry with sudo / as Administrator."
  )
}

# Final version report
for (pkg in pkgs) {
  v  <- as.character(packageVersion(pkg))
  ok <- v == expected[[pkg]]
  message(sprintf(
    "  %-14s %s  %s", pkg, v,
    if (ok) "[OK]" else paste0("[UNEXPECTED — snapshot resolved ", v, ", script recorded ", expected[[pkg]], "]")
  ))
}

message("\nAll R packages installed successfully.")

# ---- Warm Julia's precompilation cache (headless-safe) ----------------------
# GLMakie is a direct Julia dependency of AlphaPEM that requires OpenGL/X11 and
# therefore can never precompile on a headless server (e.g. an HPC cluster). On
# a fresh or changed Manifest, Julia's automatic precompile pass dies on it and
# aborts run_parameter_validity.jl before any AlphaPEM code runs.
#
# Two things are needed to avoid this, and BOTH matter:
#
#  1. JULIA_PKG_PRECOMPILE_AUTO=0 — `Pkg.instantiate()` and every `using` from a
#     cold cache trigger their OWN implicit whole-environment precompile pass,
#     which has no error tolerance at all. This disables it so the explicit call
#     below is the only precompilation that happens.
#
#  2. Excluding GLMakie by name — `Pkg.precompile(strict=false)` does NOT work:
#     `strict` only downgrades failures of *indirect* dependencies to warnings,
#     while GLMakie is a *direct* dependency and still raises
#     "ERROR: The following 1 direct dependency failed to precompile: GLMakie".
#     Naming the exclusion also keeps any genuine failure fatal.
#
# Running this once, right here, warms the on-disk compile cache so the failure
# never happens again afterwards — whether Julia is later launched through the
# cluster scripts OR directly, e.g.:
#   julia --project=. --threads=auto examples/run_parameter_validity.jl
julia_warm_expr <- paste0(
  "using Pkg; Pkg.instantiate(); ",
  "skip = [\"GLMakie\"]; ",
  "Pkg.precompile(sort([n for n in keys(Pkg.project().dependencies) if !(n in skip)]))"
)
message("\n=== Warming Julia precompilation cache (headless-safe) ===")
if (nzchar(Sys.which("julia"))) {
  julia_status <- system2(
    "julia",
    c("--project=.", "-e", shQuote(julia_warm_expr)),
    env = "JULIA_PKG_PRECOMPILE_AUTO=0"
  )
  if (julia_status == 0) {
    message("  Julia packages precompiled (GLMakie is deliberately skipped: it needs OpenGL/X11).")
  } else {
    message(sprintf("  [WARNING] Julia precompilation exited with status %d. Rerun manually if needed:", julia_status))
    message(sprintf("    JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. -e '%s'", julia_warm_expr))
  }
} else {
  message("  [INFO] `julia` not found on PATH — skipping. Before running AlphaPEM, precompile manually:")
  message(sprintf("    JULIA_PKG_PRECOMPILE_AUTO=0 julia --project=. -e '%s'", julia_warm_expr))
}

message("\nYou can now run:")
message("  RUN_PRIM=true julia --project=. --threads=auto examples/run_parameter_validity.jl")

