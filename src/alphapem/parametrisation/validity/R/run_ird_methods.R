#!/usr/bin/env Rscript
# run_prim.R
#
# Interpretable Regional Descriptors (IRD) runner for AlphaPEM.
#
# What this script does:
#   1. Load a CSV with parameter columns and a classification column
#   2. Train a Random Forest to predict validity
#   3. Run PRIM and/or MaxBox (and optionally Maire) via the irdpackage
#   4. Post-process the resulting box
#   5. Save outputs to --outdir:
#        IRD_bounds_<method>_<run_name>.yaml
#        IRD_report_<method>_<run_name>.txt
#        RF_metrics_<run_name>.txt
#
# CLI arguments:
#   --data                 Path to CSV (features + target)
#   --target               Target column name           [default: classification]
#   --positive             Positive class label         [default: valid]
#   --xinterest            Path to YAML/JSON reference config
#   --range                Probability range 'low,high' [default: 0.8,1.0]
#   --outdir               Output directory             [default: results/model_validity]
#   --methods              Comma-sep: PRIM,MaxBox,Maire [default: PRIM,MaxBox]
#   --run_name             Tag appended to output filenames
#   --categorical_overrides Comma-sep columns to treat as categorical [default: e]
#   --seed                 Random seed                  [default: 42]
#   --ird_pkg_dir          Path to the local irdpackage folder
#   --helpers_path         Absolute path to ird_helpers.R
#                          (defaults to <script_dir>/ird_helpers.R)

# ---- Check required packages -----------------------------------------------
# Packages must be pre-installed via:
#   sudo Rscript src/alphapem/parametrisation/validity/R/install_r_packages.R
# (see README.md § Installation — step 8c)
required <- c("optparse", "data.table", "mlr3", "mlr3learners",
              "iml", "yaml", "jsonlite", "tools", "R6", "checkmate", "paradox")
miss <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop(
    "Missing R packages: ", paste(miss, collapse = ", "), "\n",
    "Please run the following command from the AlphaPEM root directory:\n",
    "  sudo Rscript src/alphapem/parametrisation/validity/R/install_r_packages.R\n"
  )
}

library(optparse)
library(data.table)
library(R6)
library(checkmate)
library(paradox)
library(mlr3)
library(mlr3learners)
library(iml)
library(jsonlite)

# ---- Helper: load irdpackage ------------------------------------------------
use_ird_dev <- function(pkg_dir) {
  message("[irdpackage] Loading from: ", pkg_dir)

  r_dir <- file.path(pkg_dir, "R")
  if (!dir.exists(r_dir)) {
    stop("irdpackage R/ directory not found: ", r_dir,
         "\nMake sure the IRD package is cloned at: ", pkg_dir)
  }

  # Collect and sort R files (sorted = deterministic load order)
  r_files <- sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE,
                              recursive = FALSE))
  if (!length(r_files)) {
    stop("No .R files found in: ", r_dir)
  }

  for (f in r_files) {
    # Patch Maire.R on-the-fly to make tensorflow optional (PRIM/MaxBox don't
    # need it, but an eager tf call at parse time would crash the session).
    if (basename(f) == "Maire.R") {
      txt <- readLines(f, warn = FALSE)
      if (length(txt) >= 1 && grepl("^tf\\s*=\\s*tf\\$compat\\$v1", txt[[1]])) {
        tmp <- tempfile(fileext = ".R")
        txt[[1]] <- "tf <- tryCatch(tensorflow::tf$compat$v1, error = function(e) NULL)"
        writeLines(txt, tmp)
        source(tmp, local = FALSE)
        message("[irdpackage] Patched Maire.R for optional tensorflow loading")
        next
      }
    }
    source(f, local = FALSE)
  }
  invisible(TRUE)
}

# ---- Parse CLI --------------------------------------------------------------
option_list <- list(
  make_option(c("-d", "--data"),
    type = "character",
    help = "Path to CSV with features + target.",
    metavar = "FILE"),
  make_option(c("-t", "--target"),
    type = "character", default = "classification",
    help = "Target column [default: %default]."),
  make_option(c("-p", "--positive"),
    type = "character", default = "valid",
    help = "Positive class label [default: %default]."),
  make_option(c("-x", "--xinterest"),
    type = "character",
    default = file.path("configs", "IRD_reference_config.yaml"),
    help = "Path to YAML or JSON with x_interest mapping."),
  make_option(c("-r", "--range"),
    type = "character", default = "0.8,1.0",
    help = "Desired probability range 'low,high' [default: %default]."),
  make_option(c("-o", "--outdir"),
    type = "character", default = file.path("results", "model_validity"),
    help = "Output directory [default: %default]."),
  make_option(c("-m", "--methods"),
    type = "character", default = "PRIM,MaxBox",
    help = "Methods to run: PRIM, MaxBox, Maire [default: %default]."),
  make_option(c("-n", "--run_name"),
    type = "character", default = NULL,
    help = "Run name tag. Default is a timestamp."),
  make_option(c("-c", "--categorical_overrides"),
    type = "character", default = "e",
    help = "Comma-sep feature names to force as categorical [default: %default]."),
  make_option(c("-s", "--seed"),
    type = "integer", default = 42L,
    help = "Random seed [default: %default]."),
  make_option("--ird_pkg_dir",
    type = "character",
    default = file.path("external", "IRD_method_2023", "irdpackage"),
    help = "Path to local irdpackage directory."),
  make_option("--helpers_path",
    type = "character", default = NULL,
    help = "Absolute path to ird_helpers.R. Defaults to same dir as this script.")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Source ird_helpers.R ---------------------------------------------------
if (is.null(opt$helpers_path) || !nchar(opt$helpers_path)) {
  # Try to detect the script's own directory
  script_dir <- tryCatch({
    args     <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg)) dirname(normalizePath(sub("--file=", "", file_arg)))
    else getwd()
  }, error = function(e) getwd())
  opt$helpers_path <- file.path(script_dir, "ird_helpers.R")
}
if (!file.exists(opt$helpers_path)) {
  stop("ird_helpers.R not found at: ", opt$helpers_path,
       "\nPass --helpers_path <absolute_path> to fix this.")
}
source(opt$helpers_path)

# ---- Validate inputs --------------------------------------------------------
if (is.null(opt$data) || !file.exists(opt$data)) {
  stop("--data file not found: ", opt$data)
}
if (is.null(opt$xinterest) || !file.exists(opt$xinterest)) {
  stop("--xinterest file not found: ", opt$xinterest)
}

desired_range <- as.numeric(strsplit(trimws(opt$range), "[, ]+")[[1]])
if (length(desired_range) != 2 || any(is.na(desired_range))) {
  stop("Invalid --range. Expected format: '0.8,1.0' (got: ", opt$range, ")")
}

methods <- trimws(strsplit(opt$methods, ",")[[1]])
methods <- methods[methods %in% c("PRIM", "MaxBox", "Maire")]
if (!length(methods)) stop("No valid methods selected. Choose: PRIM, MaxBox, Maire")

categorical_overrides <- character(0)
if (!is.null(opt$categorical_overrides) && nchar(opt$categorical_overrides)) {
  categorical_overrides <- unique(trimws(strsplit(opt$categorical_overrides, ",")[[1]]))
}

if (is.null(opt$run_name) || !nchar(opt$run_name)) {
  opt$run_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Load irdpackage --------------------------------------------------------
use_ird_dev(opt$ird_pkg_dir)

# ---- Read x_interest --------------------------------------------------------
read_mapping <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "json")        return(jsonlite::fromJSON(path))
  if (ext %in% c("yml", "yaml")) return(yaml::yaml.load_file(path))
  stop("x_interest file must be JSON or YAML. Got: ", path)
}
x_interest_list <- read_mapping(opt$xinterest)

# ---- Load data --------------------------------------------------------------
dt <- fread(opt$data, data.table = FALSE)
df <- as.data.frame(dt)

if (!opt$target %in% names(df)) {
  stop("Target column not found in data: '", opt$target, "'")
}
names(df)[names(df) == opt$target] <- "validity"

# Coerce target to binary factor (positive class vs. all others → "invalid")
# Multiple negative classes (e.g. "failed", "invalid") are merged into one.
if (!is.factor(df$validity)) df$validity <- factor(df$validity)
labs <- levels(df$validity)
if (!(opt$positive %in% labs)) {
  stop("Positive class '", opt$positive, "' not in target levels: ",
       paste(labs, collapse = ", "))
}
neg_labels <- setdiff(labs, opt$positive)
neg        <- "invalid"  # canonical label for the merged negative class
if (length(neg_labels) > 1) {
  message("Merging ", length(neg_labels), " negative classes (",
          paste(neg_labels, collapse = ", "), ") into '", neg, "'")
} else if (length(neg_labels) == 1) {
  neg <- neg_labels
}
# Recode: positive stays, everything else → neg
df$validity <- ifelse(df$validity == opt$positive, opt$positive, neg)
df$validity <- factor(df$validity, levels = c(neg, opt$positive))

# ---- Select features (intersection of x_interest keys and data columns) -----
all_features <- setdiff(names(df), "validity")
xi_keys      <- names(x_interest_list)

if (!length(xi_keys)) stop("x_interest has no keys.")

features <- intersect(all_features, xi_keys)
if (!length(features)) stop("None of the x_interest keys match columns in the data.")

dropped <- setdiff(all_features, features)
if (length(dropped)) {
  message("Dropping ", length(dropped), " columns not in x_interest: ",
          paste(dropped, collapse = ", "))
}
missing_keys <- setdiff(xi_keys, all_features)
if (length(missing_keys)) {
  message("Ignoring x_interest keys not in data: ",
          paste(missing_keys, collapse = ", "))
}

x_interest           <- as.data.frame(as.list(x_interest_list[features]),
                                        stringsAsFactors = FALSE)
x_interest$validity  <- factor(opt$positive, levels = levels(df$validity))

# ---- Type harmonization ------------------------------------------------------
harmonize_types <- function(df, x_interest, features, categorical_overrides) {
  for (f in features) {
    if (f %in% categorical_overrides) {
      if (!is.factor(df[[f]])) df[[f]] <- factor(df[[f]])
      xi_val <- as.character(x_interest[[f]][1])
      levs   <- levels(df[[f]])
      if (!xi_val %in% levs) levels(df[[f]]) <- union(levs, xi_val)
      x_interest[[f]] <- factor(x_interest[[f]], levels = levels(df[[f]]))
    } else {
      suppressWarnings({
        df[[f]]         <- as.numeric(df[[f]])
        x_interest[[f]] <- as.numeric(x_interest[[f]])
      })
      if (!is.numeric(df[[f]]) || !is.numeric(x_interest[[f]])) {
        stop("Feature '", f, "' must be numeric but could not be coerced. ",
             "Use --categorical_overrides if it is categorical.")
      }
    }
  }
  # Sanity check: types must match
  df_cls <- vapply(features, function(f) class(df[[f]])[1], "")
  xi_cls <- vapply(features, function(f) class(x_interest[[f]])[1], "")
  mism   <- features[df_cls != xi_cls]
  if (length(mism)) {
    stop("Type mismatch after harmonization for: ", paste(mism, collapse = ", "))
  }
  list(df = df, x_interest = x_interest)
}

hz          <- harmonize_types(df, x_interest, features, categorical_overrides)
df          <- hz$df
x_interest  <- hz$x_interest

# Keep only selected features + target
df <- df[, c(features, "validity"), drop = FALSE]

# ---- Train Random Forest ----------------------------------------------------
set.seed(opt$seed)

task  <- TaskClassif$new(id = "pem", backend = df, target = "validity")

if (requireNamespace("ranger", quietly = TRUE)) {
  learner_id <- "classif.ranger"
} else if (requireNamespace("randomForest", quietly = TRUE)) {
  learner_id <- "classif.randomForest"
  message("ranger not available; falling back to classif.randomForest")
} else {
  stop("No supported RF learner available. Install 'ranger' or 'randomForest'.")
}

mod <- lrn(learner_id, predict_type = "prob")

split      <- partition(task, ratio = 0.8)
task_train <- task$clone()$filter(split$train)
task_test  <- task$clone()$filter(split$test)

mod$train(task_train)

pred_test  <- mod$predict(task_test)
met_test   <- pred_test$score(msrs(c("classif.acc", "classif.precision",
                                      "classif.recall", "classif.fbeta",
                                      "classif.auc")))
cat("Test set performance:\n")
print(met_test)

mod$train(task)
pred_full  <- mod$predict(task)
met_full   <- pred_full$score(msrs(c("classif.acc", "classif.precision",
                                      "classif.recall", "classif.fbeta",
                                      "classif.auc")))
cat("Full dataset performance:\n")
print(met_full)

# Save RF metrics
capture_tbl   <- function(x) paste(capture.output(print(x)), collapse = "\n")
metrics_text  <- paste0(
  "Random Forest metrics\n",
  "Seed: ",          opt$seed,     "\n",
  "Data: ",          opt$data,     "\n",
  "Target: ",        opt$target,   "    Positive class: ", opt$positive, "\n\n",
  "[Test set]\n",    capture_tbl(met_test),  "\n\n",
  "[Full dataset]\n", capture_tbl(met_full), "\n"
)
metrics_path <- file.path(opt$outdir, paste0("RF_metrics_", opt$run_name, ".txt"))
dir.create(dirname(metrics_path), showWarnings = FALSE, recursive = TRUE)
writeLines(metrics_text, metrics_path)
message("Saved RF metrics: ", metrics_path)

# ---- Build iml Predictor ---------------------------------------------------
pred <- Predictor$new(
  model = mod,
  data  = df,
  y     = "validity",
  type  = "classification",
  class = opt$positive
)

# ---- IRD runner functions ---------------------------------------------------
run_prim <- function() {
  prim  <- Prim$new(predictor = pred)
  bx    <- prim$find_box(x_interest = x_interest, desired_range = desired_range)
  post  <- PostProcessing$new(predictor = pred)$find_box(
    x_interest    = x_interest,
    desired_range = desired_range,
    box_init      = bx$box
  )
  post
}

run_maxbox <- function() {
  mb   <- MaxBox$new(predictor = pred, quiet = FALSE, strategy = "traindata")
  bx   <- mb$find_box(x_interest = x_interest, desired_range = desired_range)
  post <- PostProcessing$new(predictor = pred)$find_box(
    x_interest    = x_interest,
    desired_range = desired_range,
    box_init      = bx$box
  )
  post
}

run_maire <- function() {
  if (!requireNamespace("tensorflow", quietly = TRUE)) {
    message("tensorflow not available — skipping MAIRE.")
    return(NULL)
  }
  tensorflow::tf$compat$v1$disable_eager_execution()
  mair <- Maire$new(
    predictor         = pred,
    num_of_iterations = 100L,
    convergence       = TRUE,
    quiet             = FALSE,
    strategy          = "traindata"
  )
  bx   <- mair$find_box(x_interest = x_interest, desired_range = desired_range)
  post <- PostProcessing$new(predictor = pred, subbox_relsize = 0.3)$find_box(
    x_interest    = x_interest,
    desired_range = desired_range,
    box_init      = bx$box
  )
  post
}

# ---- Execute selected methods and export ------------------------------------
for (m in methods) {
  cat("\n=== Running", m, "===\n")

  post_box <- NULL
  if      (m == "PRIM")   post_box <- tryCatch(run_prim(),   error = function(e) { message("PRIM failed: ",   e$message); NULL })
  else if (m == "MaxBox") post_box <- tryCatch(run_maxbox(), error = function(e) { message("MaxBox failed: ", e$message); NULL })
  else if (m == "Maire")  post_box <- tryCatch(run_maire(),  error = function(e) { message("Maire failed: ",  e$message); NULL })

  if (is.null(post_box)) next

  base      <- paste0(m, "_", opt$run_name)
  yaml_path <- file.path(opt$outdir, paste0("IRD_bounds_",  base, ".yaml"))
  txt_path  <- file.path(opt$outdir, paste0("IRD_report_",  base, ".txt"))

  write_ird_yaml(
    post_box_obj         = post_box,
    data_pc              = df,
    outfile              = yaml_path,
    categorical_overrides = categorical_overrides,
    round_digits         = 6
  )

  write_ird_text_report(
    post_box_obj  = post_box,
    data_pc       = df,
    outfile       = txt_path,
    method        = m,
    desired_range = desired_range,
    desired_class = opt$positive,
    postprocessed = TRUE
  )

  cat("Saved:\n  ", yaml_path, "\n  ", txt_path, "\n")
}

message("run_prim.R — done.")



