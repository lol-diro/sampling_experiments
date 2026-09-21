# ==============================================================================
# verify_inputs.R
#
# Pure-R replacement for a shell checksum check. Verifies the MD5 fingerprints
# of the restricted/controlled inputs BEFORE running the pipeline. This does
# not modify or duplicate any check inside 00_config.R -- it is an independent,
# additive convenience you can run from a plain R console or RStudio.
#
# Usage (from an R console or RStudio, working directory does not matter):
#
#   source("verify_inputs.R")
#   verify_inputs("/full/path/to/your_data_folder")
#
# nsclc_tracerx_2017.tar.gz is expected in that same folder, alongside the
# three .RData files, and is checked automatically. Pass
# tracerx_archive_path = NULL to skip that check, or an explicit path to
# check a copy kept somewhere else.
#
# Every line printed should read OK. If any line reads MISMATCH or MISSING,
# stop -- you do not have the exact dataset this pipeline was validated
# against.
# ==============================================================================

# Same frozen fingerprints declared in 00_config.R (EXPECTED_INPUT_MD5 /
# EXPECTED_TRACERX_ARCHIVE_MD5), reproduced here for an independent check that
# does not require 00_config.R to already be configured and sourced.
.expected_input_md5 <- c(
  wgs_samples_matching.RData = "b48b33d718bde6cc0d9dae2c6bbac033",
  somatic_mutations.RData    = "30b6aa77fcf4885741cb703ea2a034cc",
  clinical_data.RData        = "f26c9f0f5f43560d525f56ac03e10526"
)

.expected_tracerx_archive_md5 <- "c4006c662ece2e125b2f05babcac8345"

verify_inputs <- function(
    data_dir,
    tracerx_archive_path = file.path(data_dir, "nsclc_tracerx_2017.tar.gz")) {

  all_ok <- TRUE

  cat("== Verifying HCC input files in:", data_dir, "==\n")

  for (fname in names(.expected_input_md5)) {
    fpath <- file.path(data_dir, fname)

    if (!file.exists(fpath)) {
      cat(sprintf("  MISSING   %s\n", fname))
      all_ok <- FALSE
      next
    }

    observed <- unname(tools::md5sum(fpath))
    expected <- unname(.expected_input_md5[[fname]])

    if (identical(observed, expected)) {
      cat(sprintf("  OK        %s\n", fname))
    } else {
      cat(sprintf(
        "  MISMATCH  %s (expected %s, observed %s)\n",
        fname, expected, observed
      ))
      all_ok <- FALSE
    }
  }

  if (!is.null(tracerx_archive_path)) {
    cat("\n== Verifying TRACERx100 archive:", tracerx_archive_path, "==\n")

    if (!file.exists(tracerx_archive_path)) {
      cat(sprintf("  MISSING   %s\n", basename(tracerx_archive_path)))
      all_ok <- FALSE
    } else {
      observed <- unname(tools::md5sum(tracerx_archive_path))
      if (identical(observed, .expected_tracerx_archive_md5)) {
        cat(sprintf("  OK        %s\n", basename(tracerx_archive_path)))
      } else {
        cat(sprintf(
          "  MISMATCH  %s (expected %s, observed %s)\n",
          basename(tracerx_archive_path),
          .expected_tracerx_archive_md5, observed
        ))
        all_ok <- FALSE
      }
    }
  } else {
    cat("\n(No tracerx_archive_path supplied: skipping the TRACERx100 archive check.)\n")
  }

  cat("\n== Verifying legacy regression file (required by scripts 03, 04) ==\n")
  legacy_path <- file.path(data_dir, "patient_level_summary.tsv")
  if (file.exists(legacy_path)) {
    cat("  OK        patient_level_summary.tsv (present; no reference checksum available for this file)\n")
  } else {
    cat("  MISSING   patient_level_summary.tsv\n")
    all_ok <- FALSE
  }

  cat("\n")
  if (all_ok) {
    cat("All checked inputs match the frozen fingerprints recorded in 00_config.R.\n")
  } else {
    cat("At least one input is missing or does not match. Do not proceed until resolved.\n")
  }

  invisible(all_ok)
}
