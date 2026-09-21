# ==============================================================================
# 04_hcc_exhaustive_sampling.R
#
# Purpose
#   Apply the validated sampling engine to the complete HCC primary event set.
#
# This script:
#   1. Enumerates EVERY non-empty sector subset for each patient with n >= 2.
#   2. Computes the validated subset-level endpoints.
#   3. Builds patient-level summaries for uniform unconstrained sampling:
#        - MEAN across all k-sector subsets = expected performance under a
#          uniformly random k-sector selection (primary estimand).
#        - MEDIAN across all subsets = legacy/typical-subset sensitivity.
#   4. Validates the cohort-wide results independently using exact
#      hypergeometric identities.
#   5. Regresses ALL legacy-equivalent metrics against the original pipeline
#      for EVERY patient x k combination.
#
# No cohort-level biological inference is performed here.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ------------------------------ locate files ---------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) > 0L) {
  script_dir <- dirname(
    normalizePath(
      sub("^--file=", "", file_arg[[1]]),
      winslash = "/",
      mustWork = TRUE
    )
  )
} else {
  script_dir <- getwd()
}

source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "functions_sampling.R"))

cat("Script 04: exhaustive primary sampling\n")
cat("=====================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
depth_dir <- file.path(PATHS$output_dir, "results", "hcc_depth")
dir.create(depth_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- constants -----------------------------------
# These are deterministic QC expectations derived from the already validated
# sector-count distribution after excluding the single n=1 patient:
#
# total non-empty subsets =
#   sum_p (2^n_p - 1) = 4590
#
# patient x k rows =
#   sum_p n_p = 489
EXPECTED_PATIENTS_N_GE_2 <- 122L
EXPECTED_TOTAL_SUBSETS <- 4590L
EXPECTED_PATIENT_K_ROWS <- 489L

TOLERANCE <- 1e-9

# ------------------------------- helpers -------------------------------------
safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

max_abs_diff_na_safe <- function(a, b) {
  if (length(a) != length(b)) return(Inf)

  a_na <- is.na(a)
  b_na <- is.na(b)

  if (any(xor(a_na, b_na))) return(Inf)

  keep <- !a_na & !b_na
  if (!any(keep)) return(0)

  max(abs(as.numeric(a[keep]) - as.numeric(b[keep])))
}

all_defined_close <- function(x, target, tolerance = TOLERANCE) {
  keep <- !is.na(x)
  if (!any(keep)) return(TRUE)
  all(abs(x[keep] - target) <= tolerance)
}

# ---------------------------- require inputs ---------------------------------
sample_map_path <- file.path(PATHS$intermediate_dir, "hcc_sample_map.rds")
event_presence_path <- file.path(
  PATHS$intermediate_dir, "hcc_event_presence.rds"
)
reference_events_path <- file.path(
  PATHS$intermediate_dir, "hcc_reference_events.rds"
)

required_paths <- c(
  sample_map_path,
  event_presence_path,
  reference_events_path
)

missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing validated Script 01/02 intermediate file(s):\n",
    paste0("  ", missing_paths, collapse = "\n")
  )
}

sample_map <- as.data.table(readRDS(sample_map_path))
event_presence <- as.data.table(readRDS(event_presence_path))
reference_events <- as.data.table(readRDS(reference_events_path))

required_sample_cols <- c("patient_id", "sample_id", "n_sectors")
required_presence_cols <- c("patient_id", "sample_id", "event_key")
required_ref_cols <- c(
  "patient_id", "event_key", "full_count",
  "is_ubiquitous", "is_nonubiquitous", "is_private"
)

if (!all(required_sample_cols %in% names(sample_map))) {
  stop("hcc_sample_map.rds is missing required columns.")
}
if (!all(required_presence_cols %in% names(event_presence))) {
  stop("hcc_event_presence.rds is missing required columns.")
}
if (!all(required_ref_cols %in% names(reference_events))) {
  stop("hcc_reference_events.rds is missing required columns.")
}

sample_map[, patient_id := as.character(patient_id)]
sample_map[, sample_id := as.character(sample_id)]
event_presence[, patient_id := as.character(patient_id)]
event_presence[, sample_id := as.character(sample_id)]
reference_events[, patient_id := as.character(patient_id)]

# Key by patient for efficient extraction from the multi-million-row objects.
setkey(sample_map, patient_id)
setkey(event_presence, patient_id)
setkey(reference_events, patient_id)

patient_table <- unique(
  sample_map[, .(patient_id, n = as.integer(n_sectors))]
)
setorder(patient_table, patient_id)

analysis_patients <- patient_table[n >= 2L]

if (nrow(analysis_patients) != EXPECTED_PATIENTS_N_GE_2) {
  stop(
    "Expected ", EXPECTED_PATIENTS_N_GE_2,
    " patients with n>=2, observed ", nrow(analysis_patients), "."
  )
}

# Independent combinatorial QC before any genomic computation.
expected_subsets_from_metadata <- sum(2^analysis_patients$n - 1)
expected_patient_k_from_metadata <- sum(analysis_patients$n)

if (expected_subsets_from_metadata != EXPECTED_TOTAL_SUBSETS) {
  stop(
    "Validated metadata no longer imply ", EXPECTED_TOTAL_SUBSETS,
    " total non-empty subsets; observed expectation = ",
    expected_subsets_from_metadata, "."
  )
}
if (expected_patient_k_from_metadata != EXPECTED_PATIENT_K_ROWS) {
  stop(
    "Validated metadata no longer imply ", EXPECTED_PATIENT_K_ROWS,
    " patient x k rows; observed expectation = ",
    expected_patient_k_from_metadata, "."
  )
}

cat(
  "Analysis population: ", nrow(analysis_patients),
  " patients with n>=2\n", sep = ""
)
cat(
  "Expected exhaustive subset rows: ",
  format(EXPECTED_TOTAL_SUBSETS, big.mark = ","),
  "\n\n", sep = ""
)

# ------------------------------------------------------------------------------
# 1. Cohort-wide exhaustive subset enumeration
# ------------------------------------------------------------------------------

subset_results_by_patient <- vector("list", nrow(analysis_patients))
analytic_results_by_patient <- vector("list", nrow(analysis_patients))

for (ii in seq_len(nrow(analysis_patients))) {
  pid <- analysis_patients$patient_id[[ii]]
  n <- analysis_patients$n[[ii]]

  cat(
    sprintf(
      "[%3d/%3d] %s: n=%d",
      ii, nrow(analysis_patients), pid, n
    )
  )

  sm <- sample_map[J(pid), .(sample_id), nomatch = 0L]
  if (nrow(sm) != n) {
    stop(
      "Patient ", pid, ": sample-map n mismatch (metadata n=", n,
      ", sample rows=", nrow(sm), ")."
    )
  }

  # Sample order is deliberately deterministic. Exhaustive unconstrained
  # results are invariant to this ordering.
  setorder(sm, sample_id)
  sample_ids <- sm$sample_id

  pres <- event_presence[
    J(pid),
    .(sample_id, event_key),
    nomatch = 0L
  ]

  ref <- reference_events[
    J(pid),
    .(
      event_key,
      full_count,
      is_ubiquitous,
      is_nonubiquitous,
      is_private
    ),
    nomatch = 0L
  ]

  if (nrow(ref) == 0L) {
    stop("Patient ", pid, " has no reference events.")
  }

  pdata <- build_patient_incidence(
    patient_presence = pres,
    patient_reference = ref,
    sample_ids = sample_ids
  )

  if (pdata$n_samples != n) {
    stop("Patient ", pid, ": incidence-matrix sector count mismatch.")
  }

  patient_subset_list <- vector("list", n)

  # Independent analytic expectations for the NEW and legacy-compatible
  # event-recovery endpoints under uniform exhaustive sampling.
  nonubiq_m <- ref$full_count[ref$is_nonubiquitous]
  n_nonubiq <- length(nonubiq_m)
  n_private <- sum(ref$is_private)

  analytic_rows <- vector("list", n)

  for (k in seq_len(n)) {
    zz <- as.data.table(compute_all_subset_metrics(pdata, k))

    expected_n_subsets <- choose(n, k)
    if (nrow(zz) != expected_n_subsets) {
      stop(
        "Patient ", pid, ", k=", k,
        ": expected ", expected_n_subsets,
        " subsets but computed ", nrow(zz), "."
      )
    }

    zz[, patient_id := pid]
    zz[, n := n]
    zz[, sampling_fraction := k / n]

    setcolorder(
      zz,
      c(
        "patient_id", "n", "k", "sampling_fraction",
        "subset_id", "subset_indices", "subset_samples",
        setdiff(
          names(zz),
          c(
            "patient_id", "n", "k", "sampling_fraction",
            "subset_id", "subset_indices", "subset_samples"
          )
        )
      )
    )

    patient_subset_list[[k]] <- zz

    # Exact event-level hypergeometric expectations.
    if (n_nonubiq > 0L) {
      denom <- choose(n, k)

      p_not_detected <- choose(n - nonubiq_m, k) / denom
      p_apparent_ubiquity <- choose(nonubiq_m, k) / denom

      expected_detection <- mean(1 - p_not_detected)
      expected_apparent_ubiquity <- mean(p_apparent_ubiquity)
      expected_heterogeneity_classification <- mean(
        1 - p_not_detected - p_apparent_ubiquity
      )
    } else {
      expected_detection <- NA_real_
      expected_apparent_ubiquity <- NA_real_
      expected_heterogeneity_classification <- NA_real_
    }

    expected_private <- if (n_private > 0L) k / n else NA_real_

    analytic_rows[[k]] <- data.table(
      patient_id = pid,
      n = n,
      k = k,
      analytic_recall_nonubiquitous_detection = expected_detection,
      analytic_recall_heterogeneity_classification =
        expected_heterogeneity_classification,
      analytic_apparent_ubiquity_error = expected_apparent_ubiquity,
      analytic_recall_private = expected_private
    )
  }

  subset_results_by_patient[[ii]] <- rbindlist(
    patient_subset_list,
    use.names = TRUE,
    fill = TRUE
  )
  analytic_results_by_patient[[ii]] <- rbindlist(analytic_rows)

  cat(
    sprintf(
      " -> %d reference events, %d subset rows\n",
      nrow(ref),
      nrow(subset_results_by_patient[[ii]])
    )
  )

  rm(
    sm, sample_ids, pres, ref, pdata,
    patient_subset_list, analytic_rows,
    nonubiq_m
  )

  # Important for the very high-mutation-burden patient(s).
  invisible(gc(FALSE))
}

subset_results <- rbindlist(
  subset_results_by_patient,
  use.names = TRUE,
  fill = TRUE
)
analytic_expectations <- rbindlist(
  analytic_results_by_patient,
  use.names = TRUE,
  fill = TRUE
)

rm(subset_results_by_patient, analytic_results_by_patient)
invisible(gc())

setorder(subset_results, patient_id, k, subset_id)
setorder(analytic_expectations, patient_id, k)

if (nrow(subset_results) != EXPECTED_TOTAL_SUBSETS) {
  stop(
    "Cohort-wide subset row count mismatch: expected ",
    EXPECTED_TOTAL_SUBSETS, ", observed ", nrow(subset_results), "."
  )
}

# ------------------------------------------------------------------------------
# 2. Patient-level unconstrained summaries
# ------------------------------------------------------------------------------

# IMPORTANT:
#   MEAN across all k-sector subsets = expected performance if each subset is
#   sampled with equal probability.
#
#   MEDIAN is retained as a legacy/typical-subset sensitivity, not as the
#   mathematical expectation of uniform random selection.

patient_summary <- subset_results[, .(
  n_subsets = .N,

  n_events_reference = n_events_reference[1L],
  n_nonubiquitous_reference = n_nonubiquitous_reference[1L],
  n_private_reference = n_private_reference[1L],

  recall_nonubiquitous_detection_mean =
    safe_mean(recall_nonubiquitous_detection),
  recall_nonubiquitous_detection_median =
    safe_median(recall_nonubiquitous_detection),

  recall_heterogeneity_classification_mean =
    safe_mean(recall_heterogeneity_classification),
  recall_heterogeneity_classification_median =
    safe_median(recall_heterogeneity_classification),

  apparent_ubiquity_error_mean =
    safe_mean(apparent_ubiquity_error),
  apparent_ubiquity_error_median =
    safe_median(apparent_ubiquity_error),

  conditional_apparent_ubiquity_rate_mean =
    safe_mean(conditional_apparent_ubiquity_rate),
  conditional_apparent_ubiquity_rate_median =
    safe_median(conditional_apparent_ubiquity_rate),

  recall_private_mean =
    safe_mean(recall_private),
  recall_private_median =
    safe_median(recall_private),

  nonubiquitous_fraction_abs_error_mean =
    safe_mean(nonubiquitous_fraction_abs_error),
  nonubiquitous_fraction_abs_error_median =
    safe_median(nonubiquitous_fraction_abs_error),

  jaccard_ith_abs_error_mean =
    safe_mean(jaccard_ith_abs_error),
  jaccard_ith_abs_error_median =
    safe_median(jaccard_ith_abs_error),

  full_nonubiquitous_fraction =
    full_nonubiquitous_fraction[1L],
  full_jaccard_ith =
    full_jaccard_ith[1L]
), by = .(patient_id, n, k, sampling_fraction)]

setorder(patient_summary, patient_id, k)

if (nrow(patient_summary) != EXPECTED_PATIENT_K_ROWS) {
  stop(
    "Patient-level summary row count mismatch: expected ",
    EXPECTED_PATIENT_K_ROWS, ", observed ", nrow(patient_summary), "."
  )
}

# Exact number of possible subsets for every patient x k.
subset_count_qc <- patient_summary[
  ,
  .(
    patient_id,
    n,
    k,
    observed_n_subsets = n_subsets,
    expected_n_subsets = choose(n, k)
  )
]
subset_count_qc[
  ,
  pass := observed_n_subsets == expected_n_subsets
]

if (!all(subset_count_qc$pass)) {
  if (DEBUG) {
    fwrite(
      subset_count_qc[pass == FALSE],
      file.path(PATHS$qc_dir, "04_subset_count_failures.tsv"),
      sep = "\t"
    )
  }
  stop("At least one patient x k combination has an incorrect subset count.")
}

# ------------------------------------------------------------------------------
# 3. Global mathematical validation of patient-level means
# ------------------------------------------------------------------------------

analytic_validation <- merge(
  patient_summary,
  analytic_expectations,
  by = c("patient_id", "n", "k"),
  all = TRUE,
  sort = TRUE
)

if (nrow(analytic_validation) != EXPECTED_PATIENT_K_ROWS) {
  stop("Analytic-validation merge produced an unexpected number of rows.")
}

analytic_validation[
  ,
  detection_abs_error :=
    abs(
      recall_nonubiquitous_detection_mean -
        analytic_recall_nonubiquitous_detection
    )
]
analytic_validation[
  ,
  heterogeneity_classification_abs_error :=
    abs(
      recall_heterogeneity_classification_mean -
        analytic_recall_heterogeneity_classification
    )
]
analytic_validation[
  ,
  apparent_ubiquity_abs_error :=
    abs(
      apparent_ubiquity_error_mean -
        analytic_apparent_ubiquity_error
    )
]
analytic_validation[
  ,
  private_abs_error :=
    abs(
      recall_private_mean -
        analytic_recall_private
    )
]

max_detection_error <- max(
  analytic_validation$detection_abs_error,
  na.rm = TRUE
)
max_classification_error <- max(
  analytic_validation$heterogeneity_classification_abs_error,
  na.rm = TRUE
)
max_apparent_ubiquity_error <- max(
  analytic_validation$apparent_ubiquity_abs_error,
  na.rm = TRUE
)
max_private_error <- max(
  analytic_validation$private_abs_error,
  na.rm = TRUE
)

analytic_pass <-
  max_detection_error <= TOLERANCE &&
  max_classification_error <= TOLERANCE &&
  max_apparent_ubiquity_error <= TOLERANCE &&
  max_private_error <= TOLERANCE

if (DEBUG) {
  fwrite(
    analytic_validation[
      ,
      .(
        patient_id, n, k,
        recall_nonubiquitous_detection_mean,
        analytic_recall_nonubiquitous_detection,
        detection_abs_error,
        recall_heterogeneity_classification_mean,
        analytic_recall_heterogeneity_classification,
        heterogeneity_classification_abs_error,
        apparent_ubiquity_error_mean,
        analytic_apparent_ubiquity_error,
        apparent_ubiquity_abs_error,
        recall_private_mean,
        analytic_recall_private,
        private_abs_error
      )
    ],
    file.path(PATHS$qc_dir, "04_analytic_expectation_validation.tsv"),
    sep = "\t"
  )
}

if (!analytic_pass) {
  stop(
    "Cohort-wide analytic expectation validation FAILED. ",
    "See 04_analytic_expectation_validation.tsv."
  )
}

# ------------------------------------------------------------------------------
# 4. Full-sampling (k=n) identities across the complete cohort
# ------------------------------------------------------------------------------

full_sampling <- patient_summary[k == n]

full_sampling_pass <-
  nrow(full_sampling) == EXPECTED_PATIENTS_N_GE_2 &&
  all_defined_close(
    full_sampling$recall_nonubiquitous_detection_mean, 1
  ) &&
  all_defined_close(
    full_sampling$recall_heterogeneity_classification_mean, 1
  ) &&
  all_defined_close(
    full_sampling$apparent_ubiquity_error_mean, 0
  ) &&
  all_defined_close(
    full_sampling$recall_private_mean, 1
  ) &&
  all_defined_close(
    full_sampling$nonubiquitous_fraction_abs_error_mean, 0
  ) &&
  all_defined_close(
    full_sampling$jaccard_ith_abs_error_mean, 0
  )

if (!full_sampling_pass) {
  if (DEBUG) {
    fwrite(
      full_sampling,
      file.path(PATHS$qc_dir, "04_full_sampling_identity_failure.tsv"),
      sep = "\t"
    )
  }
  stop("At least one cohort-wide k=n identity failed.")
}

# ------------------------------------------------------------------------------
# 5. ALL-patient / ALL-k regression against original legacy pipeline
# ------------------------------------------------------------------------------

legacy_summary_path <- file.path(
  PROJECT_ROOT,
  "patient_level_summary.tsv"
)

if (!file.exists(legacy_summary_path)) {
  stop(
    "Legacy patient_level_summary.tsv not found at:\n",
    legacy_summary_path
  )
}

legacy <- fread(legacy_summary_path)

required_legacy_cols <- c(
  "patient_id", "k", "n_subsets",
  "n_events_full", "n_nontruncal_full", "n_private_full",
  "recall_nontruncal_mean", "recall_nontruncal_median",
  "recall_private_mean", "recall_private_median",
  "branch_fraction_abs_error_mean",
  "branch_fraction_abs_error_median",
  "jaccard_ith_abs_error_mean",
  "jaccard_ith_abs_error_median",
  "full_branch_fraction",
  "full_jaccard_ith"
)

if (!all(required_legacy_cols %in% names(legacy))) {
  stop("Legacy patient_level_summary.tsv schema is not as expected.")
}

legacy[, patient_id := as.character(patient_id)]
legacy[, k := as.integer(k)]

# Prefix every non-key legacy column BEFORE merging. Several legacy columns have
# the same names as their current-pipeline counterparts (e.g. recall_private_mean); explicit
# names prevent merge suffixes from making the regression ambiguous.
legacy_cmp <- copy(legacy[, ..required_legacy_cols])
legacy_nonkeys <- setdiff(names(legacy_cmp), c("patient_id", "k"))
setnames(
  legacy_cmp,
  legacy_nonkeys,
  paste0("legacy__", legacy_nonkeys)
)

legacy_regression <- merge(
  patient_summary,
  legacy_cmp,
  by = c("patient_id", "k"),
  all = TRUE,
  sort = TRUE
)

if (nrow(legacy_regression) != EXPECTED_PATIENT_K_ROWS ||
    anyNA(legacy_regression$patient_id) ||
    anyNA(legacy_regression$k) ||
    anyNA(legacy_regression$legacy__n_subsets) ||
    anyNA(legacy_regression$n_subsets)) {
  stop(
    "Legacy regression merge failed to produce exactly ",
    EXPECTED_PATIENT_K_ROWS, " fully matched patient x k rows."
  )
}

# Exact integer/count comparisons.
count_checks <- data.table(
  metric = c(
    "n_subsets",
    "n_events_reference",
    "n_nonubiquitous_reference",
    "n_private_reference"
  ),
  legacy_metric = c(
    "legacy__n_subsets",
    "legacy__n_events_full",
    "legacy__n_nontruncal_full",
    "legacy__n_private_full"
  )
)

count_checks[
  ,
  n_mismatches := vapply(
    seq_len(.N),
    function(i) {
      current_col <- metric[[i]]
      old_col <- legacy_metric[[i]]
      sum(
        as.numeric(legacy_regression[[current_col]]) !=
          as.numeric(legacy_regression[[old_col]])
      )
    },
    integer(1)
  )
]
count_checks[, pass := n_mismatches == 0L]

# Numeric comparisons. The new names deliberately replace legacy terminology
# ("nontruncal", "branch fraction") without changing the underlying legacy
# metric definitions.
numeric_map <- data.table(
  metric = c(
    "recall_nonubiquitous_detection_mean",
    "recall_nonubiquitous_detection_median",
    "recall_private_mean",
    "recall_private_median",
    "nonubiquitous_fraction_abs_error_mean",
    "nonubiquitous_fraction_abs_error_median",
    "jaccard_ith_abs_error_mean",
    "jaccard_ith_abs_error_median",
    "full_nonubiquitous_fraction",
    "full_jaccard_ith"
  ),
  legacy_metric = c(
    "legacy__recall_nontruncal_mean",
    "legacy__recall_nontruncal_median",
    "legacy__recall_private_mean",
    "legacy__recall_private_median",
    "legacy__branch_fraction_abs_error_mean",
    "legacy__branch_fraction_abs_error_median",
    "legacy__jaccard_ith_abs_error_mean",
    "legacy__jaccard_ith_abs_error_median",
    "legacy__full_branch_fraction",
    "legacy__full_jaccard_ith"
  )
)

numeric_map[
  ,
  max_abs_error := vapply(
    seq_len(.N),
    function(i) {
      max_abs_diff_na_safe(
        legacy_regression[[metric[[i]]]],
        legacy_regression[[legacy_metric[[i]]]]
      )
    },
    numeric(1)
  )
]
numeric_map[, pass := max_abs_error <= TOLERANCE]

legacy_regression_summary <- rbindlist(
  list(
    count_checks[
      ,
      .(
        metric,
        legacy_metric,
        comparison_type = "exact_count",
        error = as.numeric(n_mismatches),
        pass
      )
    ],
    numeric_map[
      ,
      .(
        metric,
        legacy_metric,
        comparison_type = "numeric",
        error = max_abs_error,
        pass
      )
    ]
  ),
  use.names = TRUE
)

if (DEBUG) {
  fwrite(
    legacy_regression_summary,
    file.path(PATHS$qc_dir, "04_legacy_regression_summary.tsv"),
    sep = "\t"
  )
}

legacy_pass <- all(legacy_regression_summary$pass)

if (!legacy_pass) {
  if (DEBUG) {
    fwrite(
      legacy_regression,
      file.path(PATHS$qc_dir, "04_legacy_regression_full_failure.tsv"),
      sep = "\t"
    )
  }
  stop(
    "At least one cohort-wide legacy regression check FAILED. ",
    "See 04_legacy_regression_summary.tsv."
  )
}

max_legacy_numeric_error <- max(numeric_map$max_abs_error)

# ------------------------------------------------------------------------------
# 6. Save canonical primary sampling outputs
# ------------------------------------------------------------------------------

subset_rds_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_exhaustive_subsets_all_filtered.rds"
)
patient_summary_rds_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_unconstrained_patient_summary_all_filtered.rds"
)

saveRDS(subset_results, subset_rds_path)
saveRDS(patient_summary, patient_summary_rds_path)

fwrite(
  subset_results,
  file.path(
    depth_dir,
    "04_hcc_exhaustive_subsets_all_filtered.tsv"
  ),
  sep = "\t"
)

fwrite(
  patient_summary,
  file.path(
    depth_dir,
    "04_hcc_unconstrained_patient_summary_all_filtered.tsv"
  ),
  sep = "\t"
)

# Compact metadata/QC table by original n.
by_n_qc <- analysis_patients[
  ,
  .(
    n_patients = .N,
    total_exhaustive_subsets = sum(2^n - 1),
    total_patient_k_rows = sum(n)
  ),
  by = n
]
setorder(by_n_qc, n)

if (DEBUG) {
  fwrite(
    by_n_qc,
    file.path(PATHS$qc_dir, "04_exhaustive_sampling_by_n.tsv"),
    sep = "\t"
  )
}

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  "Script 04: exhaustive primary sampling",
  "=====================================================",
  "",
  sprintf(
    "Patients analysed (n>=2): %d",
    nrow(analysis_patients)
  ),
  sprintf(
    "Patient x k rows: %d",
    nrow(patient_summary)
  ),
  sprintf(
    "Exhaustive subset rows: %d",
    nrow(subset_results)
  ),
  "",
  "Primary unconstrained estimand:",
  "  patient-level MEAN across all equally likely k-sector subsets",
  "Legacy sensitivity:",
  "  patient-level MEDIAN across all k-sector subsets",
  "",
  sprintf(
    "Exact subset-count QC: PASS (%d total subsets)",
    EXPECTED_TOTAL_SUBSETS
  ),
  sprintf(
    "k=n full-reference identities: %s",
    if (full_sampling_pass) "PASS" else "FAIL"
  ),
  sprintf(
    "Analytic hypergeometric validation: %s",
    if (analytic_pass) "PASS" else "FAIL"
  ),
  sprintf(
    "  max detection-recall error: %.3g",
    max_detection_error
  ),
  sprintf(
    "  max heterogeneity-classification error: %.3g",
    max_classification_error
  ),
  sprintf(
    "  max apparent-ubiquity error: %.3g",
    max_apparent_ubiquity_error
  ),
  sprintf(
    "  max private-recall error: %.3g",
    max_private_error
  ),
  sprintf(
    "All-patient/all-k legacy regression: %s",
    if (legacy_pass) "PASS" else "FAIL"
  ),
  sprintf(
    "  max legacy numeric error: %.3g",
    max_legacy_numeric_error
  ),
  "",
  sprintf(
    "Elapsed time: %.1f seconds",
    elapsed
  ),
  "",
  "No cohort-level biological inference was performed in Script 04.",
  "Proceed to depth/cohort aggregation only if every QC above is PASS."
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "04_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
