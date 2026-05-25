# ---------------------------------------------------------------------------
#  ird_helpers.R
#
#  Provides:
#    - write_ird_yaml():        save a PostProcessing box to a YAML file
#    - write_ird_text_report(): save a human-readable report to a text file
# ---------------------------------------------------------------------------

suppressWarnings({
  if (!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
})

library(yaml)

# ----------------------------------------------------------
# Boolean handler so YAML shows True/False instead of yes/no
# ----------------------------------------------------------

yaml_bool_handler <- function(x) {
  structure(ifelse(x, "True", "False"), class = "verbatim")
}

# ----------------------------------------------------------
# write_ird_yaml()
# Write an IRD PostProcessing result box to a YAML file.
#
# The output format is:
#   parameters:
#     - name: <feature>
#       type: continuous | categorical
#       low: <number>      # continuous only
#       high: <number>     # continuous only
#       values: [...]      # categorical only
#       fixed: True|False
# ----------------------------------------------------------
write_ird_yaml <- function(
  post_box_obj,
  data_pc,
  outfile,
  categorical_overrides = character(0),
  round_digits = 6
) {
  feature_names <- names(post_box_obj$x_interest)

  lower_rd  <- post_box_obj$box$lower
  upper_rd  <- post_box_obj$box$upper
  levels_rd <- post_box_obj$box$levels
  if (is.null(levels_rd)) levels_rd <- vector("list", length(feature_names))
  if (is.null(names(levels_rd))) names(levels_rd) <- feature_names

  is_categorical <- function(f) {
    f %in% categorical_overrides ||
      is.factor(data_pc[[f]]) || is.character(data_pc[[f]])
  }

  safe_round <- function(x) {
    if (is.numeric(x)) return(signif(x, round_digits))
    suppressWarnings({
      xx <- as.numeric(x)
      if (!any(is.na(xx))) return(signif(xx, round_digits))
    })
    x
  }

  normalize_levels <- function(vals) {
    if (is.null(vals)) return(NULL)
    v <- unlist(vals, use.names = FALSE)
    suppressWarnings({
      vv <- as.numeric(v)
      if (!any(is.na(vv))) return(signif(vv, round_digits))
    })
    as.character(v)
  }

  params <- lapply(feature_names, function(f) {
    entry <- list(name = f)
    if (is_categorical(f)) {
      entry$type <- "categorical"
      vals <- normalize_levels(levels_rd[[f]])
      if (is.null(vals)) {
        if (is.factor(data_pc[[f]])) {
          vals <- levels(data_pc[[f]])
        } else {
          vals <- unique(as.character(data_pc[[f]]))
        }
      }
      entry$values <- vals
      entry$fixed  <- length(vals) == 1
    } else {
      entry$type <- "continuous"
      low  <- as.numeric(lower_rd[[f]])
      high <- as.numeric(upper_rd[[f]])
      if (is.na(low) || is.na(high)) {
        low  <- suppressWarnings(min(as.numeric(data_pc[[f]]), na.rm = TRUE))
        high <- suppressWarnings(max(as.numeric(data_pc[[f]]), na.rm = TRUE))
      }
      entry$low   <- safe_round(low)
      entry$high  <- safe_round(high)
      entry$fixed <- isTRUE(all.equal(as.numeric(low), as.numeric(high)))
    }
    entry
  })

  out_list  <- list(parameters = params)
  yaml_text <- as.yaml(
    out_list,
    indent   = 2,
    line.sep = "\n",
    handlers = list(logical = yaml_bool_handler)
  )
  writeLines(yaml_text, con = outfile)
  invisible(outfile)
}


# ----------------------------------------------------------
# write_ird_text_report()
# Save a console-style IRD report to a plain-text file.
# ----------------------------------------------------------
write_ird_text_report <- function(
  post_box_obj,
  data_pc,
  outfile,
  method,
  desired_range,
  desired_class  = "valid",
  postprocessed  = TRUE,
  digits_main    = 2,
  digits_scientific = 3,
  equal_tol      = 1e-12
) {
  stopifnot(length(desired_range) == 2)

  fmt_num <- function(x) {
    if (!is.numeric(x) || is.na(x)) return(as.character(x))
    ax <- abs(x)
    if ((ax >= 1e5) || (ax > 0 && ax < 1e-3)) {
      return(format(x, scientific = TRUE,  digits = digits_scientific, trim = TRUE))
    } else {
      return(format(round(x, digits_main), nsmall = digits_main, trim = TRUE, scientific = FALSE))
    }
  }

  fmt_set <- function(vals) {
    if (is.null(vals)) return("{}")
    v <- unlist(vals, use.names = FALSE)
    suppressWarnings({
      vv <- as.numeric(v)
      if (!any(is.na(vv)))
        return(paste0("{", paste(vapply(vv, fmt_num, character(1)), collapse = ", "), "}"))
    })
    paste0("{", paste(as.character(v), collapse = ", "), "}")
  }

  fmt_interval <- function(lo, hi) {
    if (is.na(lo) || is.na(hi)) return("[NA, NA]")
    if (isTRUE(abs(hi - lo) <= equal_tol)) return(paste0("{", fmt_num(lo), "}"))
    paste0("[", fmt_num(lo), ", ", fmt_num(hi), "]")
  }

  pad_left <- function(s, w) sprintf("%-*s", w, s)

  feature_names <- names(post_box_obj$x_interest)
  lower_rd  <- post_box_obj$box$lower
  upper_rd  <- post_box_obj$box$upper
  lower_1d  <- post_box_obj$box_single$lower
  upper_1d  <- post_box_obj$box_single$upper
  levels_rd <- post_box_obj$box$levels
  levels_1d <- post_box_obj$box_single$levels
  metrics   <- post_box_obj$evaluate()

  col_xint <- vapply(feature_names, function(f) {
    xi <- post_box_obj$x_interest[[f]][1]
    if (is.numeric(xi)) fmt_num(xi) else as.character(xi)
  }, character(1))

  is_cat <- vapply(feature_names, function(f) {
    has_lev <- !is.null(levels_rd) && !is.null(levels_rd[[f]]) && length(levels_rd[[f]]) > 0
    has_lev || is.factor(data_pc[[f]]) || is.character(data_pc[[f]])
  }, logical(1))

  col_rd <- vapply(seq_along(feature_names), function(i) {
    f <- feature_names[i]
    if (is_cat[i]) fmt_set(if (!is.null(levels_rd[[f]])) levels_rd[[f]] else unique(data_pc[[f]]))
    else fmt_interval(suppressWarnings(as.numeric(lower_rd[[f]])),
                      suppressWarnings(as.numeric(upper_rd[[f]])))
  }, character(1))

  col_1d <- vapply(seq_along(feature_names), function(i) {
    f <- feature_names[i]
    if (is_cat[i]) fmt_set(if (!is.null(levels_1d[[f]])) levels_1d[[f]] else unique(data_pc[[f]]))
    else fmt_interval(suppressWarnings(as.numeric(lower_1d[[f]])),
                      suppressWarnings(as.numeric(upper_1d[[f]])))
  }, character(1))

  col_range <- vapply(seq_along(feature_names), function(i) {
    f <- feature_names[i]
    if (is_cat[i]) {
      vals <- if (is.factor(data_pc[[f]])) levels(data_pc[[f]]) else unique(as.character(data_pc[[f]]))
      fmt_set(vals)
    } else {
      suppressWarnings({
        lo <- min(as.numeric(data_pc[[f]]), na.rm = TRUE)
        hi <- max(as.numeric(data_pc[[f]]), na.rm = TRUE)
      })
      fmt_interval(lo, hi)
    }
  }, character(1))

  header <- c("feature", "x_interest", "regional descriptor", "1-dim descriptor", "range")
  M      <- cbind(feature_names, col_xint, col_rd, col_1d, col_range)
  widths <- pmax(nchar(header),
                 apply(M, 2, function(col) max(nchar(col), na.rm = TRUE)))

  fmt_row <- function(v) paste(mapply(pad_left, v, widths), collapse = "  ")

  table_lines <- c(
    fmt_row(header),
    vapply(seq_len(nrow(M)), function(i) fmt_row(M[i, ]), character(1))
  )

  if (is.null(names(metrics))) names(metrics) <- paste0("metric_", seq_along(metrics))
  met_names  <- names(metrics)
  met_vals   <- vapply(metrics, function(v)
    if (is.numeric(v)) format(v, digits = 8, scientific = FALSE, trim = TRUE)
    else as.character(v), character(1))
  met_widths <- pmax(nchar(met_names), nchar(met_vals)) + 2
  met_line1  <- paste(mapply(function(n, w) sprintf("%-*s", w, n), met_names, met_widths), collapse = "")
  met_line2  <- paste(mapply(function(v, w) sprintf("%-*s", w, v), met_vals,  met_widths), collapse = "")

  lines <- c(
    "Regional Descriptors",
    "",
    paste0("Method: ",          method),
    paste0("Post-processed: ",  if (postprocessed) "True" else "False"),
    paste0("Desired class: ",   desired_class),
    paste0("Desired range: [",  fmt_num(desired_range[1]), ", ", fmt_num(desired_range[2]), "]"),
    "",
    "Descriptor:",
    table_lines,
    "",
    met_line1,
    met_line2,
    ""
  )
  writeLines(lines, con = outfile)
  invisible(outfile)
}

