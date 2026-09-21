# ==============================================================================
# generate_validation_summary.R
#
# Consolidates every per-script `*_summary.txt` narrative under
# `<output>/results/` into one file, `results/VALIDATION_SUMMARY.txt`, listing
# every PASS/FAIL line found and an overall verdict. This is the single
# validation log meant to travel with the bundle handed to collaborators.
#
# It only reads text these scripts already wrote; it does not re-run or
# re-check anything itself. Most individual checks would already have
# stopped the pipeline via stop() if they failed (see docs/PIPELINE_MAP.md ->
# "Known operational caveat" for the one documented exception, Script 01's
# legacy QC, which is reported as a warning and so can appear here as FAIL
# even though the pipeline continued).
#
# Usage (from an R console or RStudio, after running some or all of the
# pipeline):
#
#   source("generate_validation_summary.R")
#   generate_validation_summary("/full/path/to/your_data_folder/output")
#
# Or, with no argument, it looks for ./output relative to the current
# working directory.
# ==============================================================================

generate_validation_summary <- function(output_dir = "output") {

  results_dir <- file.path(output_dir, "results")

  if (!dir.exists(results_dir)) {
    stop(
      "No results/ directory found at: ", results_dir, "\n",
      "Run the pipeline (or at least some scripts) first, or pass the ",
      "correct output_dir."
    )
  }

  summary_files <- sort(
    list.files(
      results_dir,
      pattern = "_summary\\.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
  )

  if (length(summary_files) == 0L) {
    stop(
      "No *_summary.txt files found under: ", results_dir, "\n",
      "Run at least one script before generating the validation summary."
    )
  }

  report <- c(
    "==============================================================",
    "CONSOLIDATED VALIDATION SUMMARY",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(
      "Scanned ", length(summary_files),
      " summary file(s) under ", results_dir
    ),
    "=============================================================="
  )

  n_pass <- 0L
  n_fail <- 0L
  failing_files <- character()

  for (f in summary_files) {
    lines <- readLines(f, warn = FALSE)
    pass_lines <- grep("PASS", lines, value = TRUE)
    fail_lines <- grep("FAIL(ED)?", lines, value = TRUE)

    rel <- sub(results_dir, "", f, fixed = TRUE)
    rel <- sub("^[/\\\\]", "", rel)

    report <- c(report, "", paste0("-- ", rel, " --"))

    if (length(fail_lines) > 0L) {
      failing_files <- c(failing_files, rel)
      n_fail <- n_fail + length(fail_lines)
      report <- c(report, paste0("  [FAIL] ", fail_lines))
    }
    if (length(pass_lines) > 0L) {
      n_pass <- n_pass + length(pass_lines)
      report <- c(report, paste0("  [PASS] ", pass_lines))
    }
    if (length(pass_lines) == 0L && length(fail_lines) == 0L) {
      report <- c(
        report,
        "  (no PASS/FAIL line found in this file -- inspect manually)"
      )
    }
  }

  report <- c(
    report,
    "",
    "==============================================================",
    sprintf(
      "TOTAL: %d PASS line(s), %d FAIL line(s) across %d file(s).",
      n_pass, n_fail, length(summary_files)
    )
  )

  if (n_fail > 0L) {
    report <- c(
      report,
      "",
      "OVERALL: FAIL -- inspect the file(s) listed above before trusting results.",
      paste0("Files with a FAIL line: ", paste(failing_files, collapse = ", "))
    )
  } else {
    report <- c(
      report,
      "",
      "OVERALL: PASS -- every scanned summary reported PASS, no FAIL/FAILED lines found."
    )
  }

  out_path <- file.path(results_dir, "VALIDATION_SUMMARY.txt")
  writeLines(report, out_path)

  cat(paste(report, collapse = "\n"), "\n")
  cat("\nWritten to:", out_path, "\n")

  invisible(list(path = out_path, n_pass = n_pass, n_fail = n_fail))
}
