# ==============================================================================
# 11_tracerx100_exhaustive_sampling.R
#
# Purpose
#   Apply the already validated HCC exhaustive sampling engine to the
#   frozen TRACERx100 protein-altering primary-tumour event reference.
#
# Canonical TRACERx100 input (Scripts 09-10):
#   - 100 patients
#   - 323 primary tumour regions (-R#)
#   - region-count distribution:
#       n=2:40, n=3:23, n=4:22, n=5:9, n=6:1, n=7:5
#   - 31,919 patient-specific protein-altering exact variants
#   - 8,946 full-reference non-ubiquitous events
#   - 5,484 private events
#   - four n=2 patients have zero non-ubiquitous/private events:
#       CRUK0040, CRUK0059, CRUK0061, CRUK0090
#
# Scientific estimand
#   For each patient and each k, enumerate EVERY choose(n,k) primary-region
#   subset. Thus the patient-level mean across subsets is the exact expectation
#   under uniformly random/unconstrained selection of k primary tumour regions.
#
# Metrics are identical to the validated HCC engine:
#   - non-ubiquitous detection recall
#   - heterogeneity-classification recall
#   - apparent-ubiquity error
#   - conditional apparent-ubiquity rate
#   - private recall
#   - non-ubiquitous-fraction absolute error
#   - Jaccard ITH absolute error
#
# Zero-denominator rule
#   If a patient's full-reference non-ubiquitous count is zero, detection,
#   classification and apparent-ubiquity metrics remain NA (never zero).
#   If private count is zero, private recall remains NA.
#
# This script performs exhaustive sampling + mathematical/QC validation only.
# Cohort-level learning-curve inference and HCC comparison are deferred to
# Script 12.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ------------------------------ locate script --------------------------------
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

cat("Script 11: TRACERx100 exhaustive sampling\n")
cat("======================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
tracerx_result_dir <- file.path(PATHS$output_dir, "results", "tracerx100")
dir.create(tracerx_result_dir, recursive = TRUE, showWarnings = FALSE)

sample_map_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_primary_sample_map.rds"
)
presence_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_event_presence_protein_altering.rds"
)
reference_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_reference_events_protein_altering.rds"
)
reference_summary_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_patient_reference_summary_protein_altering.rds"
)

for (p in c(
  sample_map_path,
  presence_path,
  reference_path,
  reference_summary_path
)) {
  if (!file.exists(p)) {
    stop("Missing validated Script 09/10 input:\n", p)
  }
}

# ------------------------------- constants -----------------------------------
EXPECTED_N_DISTRIBUTION <- c(
  `2` = 40L,
  `3` = 23L,
  `4` = 22L,
  `5` = 9L,
  `6` = 1L,
  `7` = 5L
)

EXPECTED_PATIENTS <- 100L
EXPECTED_REGIONS <- 323L
EXPECTED_TOTAL_SUBSETS <- 1588L
EXPECTED_PATIENT_K_ROWS <- 323L

EXPECTED_REFERENCE_EVENTS <- 31919L
EXPECTED_NONUBIQUITOUS_EVENTS <- 8946L
EXPECTED_PRIVATE_EVENTS <- 5484L

EXPECTED_ZERO_NONUBIQ <- sort(c(
  "CRUK0040", "CRUK0059", "CRUK0061", "CRUK0090"
))

TOLERANCE <- 1e-10

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

# ------------------------------- 1. load --------------------------------------
sample_map <- as.data.table(readRDS(sample_map_path))
presence <- as.data.table(readRDS(presence_path))
reference <- as.data.table(readRDS(reference_path))
reference_summary <- as.data.table(readRDS(reference_summary_path))

required_sample_cols <- c(
  "patient_id", "sample_id", "region_number"
)
required_presence_cols <- c(
  "patient_id", "sample_id", "event_key"
)
required_reference_cols <- c(
  "patient_id", "event_key", "full_count", "n_regions",
  "is_ubiquitous", "is_nonubiquitous", "is_private"
)
required_refsum_cols <- c(
  "patient_id", "n_regions",
  "n_events_reference",
  "n_nonubiquitous_reference",
  "n_private_reference"
)

if (!all(required_sample_cols %in% names(sample_map))) {
  stop("TRACERx100 sample map schema is not as expected.")
}
if (!all(required_presence_cols %in% names(presence))) {
  stop("TRACERx100 event-presence schema is not as expected.")
}
if (!all(required_reference_cols %in% names(reference))) {
  stop("TRACERx100 event-reference schema is not as expected.")
}
if (!all(required_refsum_cols %in% names(reference_summary))) {
  stop("TRACERx100 patient reference-summary schema is not as expected.")
}

sample_map[
  ,
  `:=`(
    patient_id = as.character(patient_id),
    sample_id = as.character(sample_id),
    region_number = as.integer(region_number)
  )
]
presence[
  ,
  `:=`(
    patient_id = as.character(patient_id),
    sample_id = as.character(sample_id),
    event_key = as.character(event_key)
  )
]
reference[
  ,
  `:=`(
    patient_id = as.character(patient_id),
    event_key = as.character(event_key)
  )
]
reference_summary[
  ,
  patient_id := as.character(patient_id)
]

# --------------------------- 2. input fingerprints ----------------------------
patient_table <- sample_map[
  ,
  .(n = .N),
  by = patient_id
]
setorder(patient_table, patient_id)

n_distribution <- patient_table[
  ,
  .N,
  by = n
]
setorder(n_distribution, n)

expected_n_dt <- data.table(
  n = as.integer(names(EXPECTED_N_DISTRIBUTION)),
  expected_N = as.integer(EXPECTED_N_DISTRIBUTION)
)
n_distribution_qc <- merge(
  expected_n_dt,
  n_distribution[, .(n, observed_N = N)],
  by = "n",
  all = TRUE,
  sort = TRUE
)
n_distribution_qc[
  ,
  pass := expected_N == observed_N
]

if (nrow(patient_table) != EXPECTED_PATIENTS ||
    nrow(sample_map) != EXPECTED_REGIONS ||
    !all(n_distribution_qc$pass)) {
  if (DEBUG) {
    fwrite(
      n_distribution_qc,
      file.path(
        PATHS$qc_dir,
        "11_tracerx100_region_count_distribution_failure.tsv"
      ),
      sep = "\t"
    )
  }
  stop("TRACERx100 cohort fingerprint QC FAILED.")
}

if (nrow(reference) != EXPECTED_REFERENCE_EVENTS ||
    sum(reference$is_nonubiquitous) != EXPECTED_NONUBIQUITOUS_EVENTS ||
    sum(reference$is_private) != EXPECTED_PRIVATE_EVENTS) {
  stop("TRACERx100 event-reference fingerprint QC FAILED.")
}

zero_nonubiq_ids <- sort(
  reference_summary[
    n_nonubiquitous_reference == 0L,
    patient_id
  ]
)

if (!identical(zero_nonubiq_ids, EXPECTED_ZERO_NONUBIQ)) {
  stop("TRACERx100 zero-denominator patient fingerprint QC FAILED.")
}

# Pure combinatorial expectations from the frozen region-count distribution.
expected_total_subsets <- sum(2^patient_table$n - 1L)
expected_patient_k_rows <- sum(patient_table$n)

if (expected_total_subsets != EXPECTED_TOTAL_SUBSETS) {
  stop(
    "Expected ", EXPECTED_TOTAL_SUBSETS,
    " total non-empty subsets; metadata imply ",
    expected_total_subsets, "."
  )
}
if (expected_patient_k_rows != EXPECTED_PATIENT_K_ROWS) {
  stop(
    "Expected ", EXPECTED_PATIENT_K_ROWS,
    " patient x k rows; metadata imply ",
    expected_patient_k_rows, "."
  )
}

cat(
  "Analysis population: ", EXPECTED_PATIENTS,
  " patients / ", EXPECTED_REGIONS, " primary regions\n", sep = ""
)
cat(
  "Expected exhaustive subset rows: ",
  EXPECTED_TOTAL_SUBSETS, "\n\n", sep = ""
)

setkey(sample_map, patient_id)
setkey(presence, patient_id)
setkey(reference, patient_id)

# ------------------------------------------------------------------------------
# 3. Exhaustive sampling for every patient / every k
# ------------------------------------------------------------------------------

subset_results_by_patient <- vector("list", nrow(patient_table))
analytic_results_by_patient <- vector("list", nrow(patient_table))

for (ii in seq_len(nrow(patient_table))) {
  pid <- patient_table$patient_id[[ii]]
  n <- patient_table$n[[ii]]

  cat(
    sprintf(
      "[%3d/%3d] %s: n=%d",
      ii, nrow(patient_table), pid, n
    )
  )

  sm <- sample_map[
    J(pid),
    .(sample_id, region_number),
    nomatch = 0L
  ]

  if (nrow(sm) != n) {
    stop(pid, ": sample-map n mismatch.")
  }

  # Deterministic order by primary region number. Exhaustive unconstrained
  # results are invariant to ordering, but this makes subset IDs reproducible.
  setorder(sm, region_number, sample_id)
  sample_ids <- sm$sample_id

  pres <- presence[
    J(pid),
    .(sample_id, event_key),
    nomatch = 0L
  ]

  ref <- reference[
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
    stop(pid, ": no protein-altering reference events.")
  }

  pdata <- build_patient_incidence(
    patient_presence = pres,
    patient_reference = ref,
    sample_ids = sample_ids
  )

  if (pdata$n_samples != n) {
    stop(pid, ": incidence-matrix region count mismatch.")
  }

  patient_subset_list <- vector("list", n)
  analytic_rows <- vector("list", n)

  nonubiq_m <- ref$full_count[ref$is_nonubiquitous]
  n_nonubiq <- length(nonubiq_m)
  n_private <- sum(ref$is_private)

  for (kval in seq_len(n)) {
    zz <- as.data.table(
      compute_all_subset_metrics(pdata, kval)
    )

    if (nrow(zz) != choose(n, kval)) {
      stop(
        pid, " / k=", kval,
        ": subset-count mismatch."
      )
    }

    zz[
      ,
      `:=`(
        patient_id = pid,
        n = n,
        k = kval,
        sampling_fraction = kval / n
      )
    ]

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

    patient_subset_list[[kval]] <- zz

    # Independent exact hypergeometric expectations.
    if (n_nonubiq > 0L) {
      denom <- choose(n, kval)
      p_not_detected <- choose(n - nonubiq_m, kval) / denom
      p_apparent_ubiquity <- choose(nonubiq_m, kval) / denom

      expected_detection <- mean(1 - p_not_detected)
      expected_apparent <- mean(p_apparent_ubiquity)
      expected_classification <- mean(
        1 - p_not_detected - p_apparent_ubiquity
      )
    } else {
      expected_detection <- NA_real_
      expected_apparent <- NA_real_
      expected_classification <- NA_real_
    }

    expected_private <- if (n_private > 0L) {
      kval / n
    } else {
      NA_real_
    }

    analytic_rows[[kval]] <- data.table(
      patient_id = pid,
      n = n,
      k = kval,
      analytic_recall_nonubiquitous_detection =
        expected_detection,
      analytic_recall_heterogeneity_classification =
        expected_classification,
      analytic_apparent_ubiquity_error =
        expected_apparent,
      analytic_recall_private =
        expected_private
    )
  }

  subset_results_by_patient[[ii]] <- rbindlist(
    patient_subset_list,
    use.names = TRUE,
    fill = TRUE
  )

  analytic_results_by_patient[[ii]] <- rbindlist(
    analytic_rows,
    use.names = TRUE,
    fill = TRUE
  )

  cat(
    sprintf(
      " -> %d reference events, %d subset rows\n",
      nrow(ref),
      nrow(subset_results_by_patient[[ii]])
    )
  )

  rm(
    sm, sample_ids, pres, ref, pdata,
    patient_subset_list, analytic_rows, nonubiq_m
  )
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

setorder(subset_results, patient_id, k, subset_id)
setorder(analytic_expectations, patient_id, k)

if (nrow(subset_results) != EXPECTED_TOTAL_SUBSETS) {
  stop(
    "Exhaustive row-count QC FAILED: expected ",
    EXPECTED_TOTAL_SUBSETS, ", observed ",
    nrow(subset_results), "."
  )
}

# ------------------------------------------------------------------------------
# 4. Patient-level unconstrained summaries
# ------------------------------------------------------------------------------

patient_summary <- subset_results[
  ,
  .(
    n_subsets = .N,

    n_events_reference = n_events_reference[1L],
    n_nonubiquitous_reference =
      n_nonubiquitous_reference[1L],
    n_private_reference =
      n_private_reference[1L],

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
  ),
  by = .(
    patient_id,
    n,
    k,
    sampling_fraction
  )
]

setorder(patient_summary, patient_id, k)

if (nrow(patient_summary) != EXPECTED_PATIENT_K_ROWS) {
  stop(
    "Patient-summary row-count QC FAILED: expected ",
    EXPECTED_PATIENT_K_ROWS, ", observed ",
    nrow(patient_summary), "."
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
    expected_n_subsets = choose(n, k),
    pass = n_subsets == choose(n, k)
  )
]

if (!all(subset_count_qc$pass)) {
  if (DEBUG) {
    fwrite(
      subset_count_qc[pass == FALSE],
      file.path(
        PATHS$qc_dir,
        "11_tracerx100_subset_count_failures.tsv"
      ),
      sep = "\t"
    )
  }
  stop("At least one patient x k has an incorrect subset count.")
}

# ------------------------------------------------------------------------------
# 5. Independent analytic validation of patient-level means
# ------------------------------------------------------------------------------

analytic_validation <- merge(
  patient_summary,
  analytic_expectations,
  by = c("patient_id", "n", "k"),
  all = TRUE,
  sort = TRUE
)

if (nrow(analytic_validation) != EXPECTED_PATIENT_K_ROWS) {
  stop("Analytic-validation merge row count is incorrect.")
}

analytic_validation[
  ,
  detection_abs_error := abs(
    recall_nonubiquitous_detection_mean -
      analytic_recall_nonubiquitous_detection
  )
]
analytic_validation[
  ,
  classification_abs_error := abs(
    recall_heterogeneity_classification_mean -
      analytic_recall_heterogeneity_classification
  )
]
analytic_validation[
  ,
  apparent_ubiquity_abs_error := abs(
    apparent_ubiquity_error_mean -
      analytic_apparent_ubiquity_error
  )
]
analytic_validation[
  ,
  private_abs_error := abs(
    recall_private_mean -
      analytic_recall_private
  )
]

max_detection_error <- max(
  analytic_validation$detection_abs_error,
  na.rm = TRUE
)
max_classification_error <- max(
  analytic_validation$classification_abs_error,
  na.rm = TRUE
)
max_apparent_error <- max(
  analytic_validation$apparent_ubiquity_abs_error,
  na.rm = TRUE
)
max_private_error <- max(
  analytic_validation$private_abs_error,
  na.rm = TRUE
)

analytic_pass <- all(c(
  max_detection_error,
  max_classification_error,
  max_apparent_error,
  max_private_error
) <= TOLERANCE)

if (!analytic_pass) {
  if (DEBUG) {
    fwrite(
      analytic_validation,
      file.path(
        PATHS$qc_dir,
        "11_tracerx100_analytic_validation_failure.tsv"
      ),
      sep = "\t"
    )
  }
  stop("TRACERx100 analytic hypergeometric validation FAILED.")
}

# Explicitly verify that the four zero-denominator patients remain NA for the
# relevant metrics at EVERY k.
zero_metric_rows <- patient_summary[
  patient_id %in% EXPECTED_ZERO_NONUBIQ
]

zero_denominator_pass <-
  nrow(zero_metric_rows) == 8L &&
  all(is.na(
    zero_metric_rows$recall_nonubiquitous_detection_mean
  )) &&
  all(is.na(
    zero_metric_rows$recall_heterogeneity_classification_mean
  )) &&
  all(is.na(
    zero_metric_rows$apparent_ubiquity_error_mean
  )) &&
  all(is.na(
    zero_metric_rows$recall_private_mean
  ))

if (!zero_denominator_pass) {
  if (DEBUG) {
    fwrite(
      zero_metric_rows,
      file.path(
        PATHS$qc_dir,
        "11_tracerx100_zero_denominator_failure.tsv"
      ),
      sep = "\t"
    )
  }
  stop("Zero-denominator NA propagation QC FAILED.")
}

# ------------------------------------------------------------------------------
# 6. Full-reference k=n identities
# ------------------------------------------------------------------------------

full_sampling <- patient_summary[k == n]

if (nrow(full_sampling) != EXPECTED_PATIENTS) {
  stop("Expected one k=n row per patient.")
}

defined_nonubiq <- full_sampling$n_nonubiquitous_reference > 0L
defined_private <- full_sampling$n_private_reference > 0L

full_identity_pass <-
  all(
    abs(
      full_sampling[
        defined_nonubiq,
        recall_nonubiquitous_detection_mean
      ] - 1
    ) <= TOLERANCE
  ) &&
  all(
    abs(
      full_sampling[
        defined_nonubiq,
        recall_heterogeneity_classification_mean
      ] - 1
    ) <= TOLERANCE
  ) &&
  all(
    abs(
      full_sampling[
        defined_nonubiq,
        apparent_ubiquity_error_mean
      ]
    ) <= TOLERANCE
  ) &&
  all(
    abs(
      full_sampling[
        defined_private,
        recall_private_mean
      ] - 1
    ) <= TOLERANCE
  ) &&
  all(
    abs(
      full_sampling$nonubiquitous_fraction_abs_error_mean
    ) <= TOLERANCE
  ) &&
  all(
    abs(
      full_sampling$jaccard_ith_abs_error_mean
    ) <= TOLERANCE
  )

if (!full_identity_pass) {
  if (DEBUG) {
    fwrite(
      full_sampling,
      file.path(
        PATHS$qc_dir,
        "11_tracerx100_full_sampling_identity_failure.tsv"
      ),
      sep = "\t"
    )
  }
  stop("TRACERx100 k=n full-reference identity QC FAILED.")
}

# ------------------------------------------------------------------------------
# 7. Regression against Script 10 reference counts
# ------------------------------------------------------------------------------

regression <- merge(
  patient_summary[
    k == n,
    .(
      patient_id,
      n,
      n_events_reference,
      n_nonubiquitous_reference,
      n_private_reference
    )
  ],
  reference_summary[
    ,
    .(
      patient_id,
      n_script10 = n_regions,
      n_events_script10 = n_events_reference,
      n_nonubiquitous_script10 =
        n_nonubiquitous_reference,
      n_private_script10 =
        n_private_reference
    )
  ],
  by = "patient_id",
  all = TRUE,
  sort = TRUE
)

regression[
  ,
  pass := (
    n == n_script10 &
    n_events_reference == n_events_script10 &
    n_nonubiquitous_reference ==
      n_nonubiquitous_script10 &
    n_private_reference ==
      n_private_script10
  )
]

if (nrow(regression) != EXPECTED_PATIENTS ||
    anyNA(regression$pass) ||
    !all(regression$pass)) {
  if (DEBUG) {
    fwrite(
      regression,
      file.path(
        PATHS$qc_dir,
        "11_tracerx100_script10_reference_regression_failure.tsv"
      ),
      sep = "\t"
    )
  }
  stop("Regression against Script 10 reference counts FAILED.")
}

# ------------------------------------------------------------------------------
# 8. Additional monotonicity QC
# ------------------------------------------------------------------------------

monotonic_cols <- c(
  "recall_nonubiquitous_detection_mean",
  "recall_heterogeneity_classification_mean",
  "recall_private_mean"
)

monotonic_failures <- list()
mfi <- 0L

for (pid in unique(patient_summary$patient_id)) {
  d <- patient_summary[patient_id == pid][order(k)]

  for (cc in monotonic_cols) {
    vv <- d[[cc]]

    # Zero-denominator NA series are valid and skipped.
    if (all(is.na(vv))) next

    dd <- diff(vv)
    bad <- which(!is.na(dd) & dd < -TOLERANCE)

    if (length(bad) > 0L) {
      mfi <- mfi + 1L
      monotonic_failures[[mfi]] <- data.table(
        patient_id = pid,
        metric = cc,
        k_from = d$k[bad],
        k_to = d$k[bad + 1L],
        value_from = vv[bad],
        value_to = vv[bad + 1L]
      )
    }
  }
}

if (length(monotonic_failures) > 0L) {
  if (DEBUG) {
    fwrite(
      rbindlist(monotonic_failures),
      file.path(
        PATHS$qc_dir,
        "11_tracerx100_monotonicity_failures.tsv"
      ),
      sep = "\t"
    )
  }
  stop("TRACERx100 patient-level monotonicity QC FAILED.")
}

# ------------------------------------------------------------------------------
# 9. Save canonical Script 11 outputs
# ------------------------------------------------------------------------------


saveRDS(
  patient_summary,
  file.path(
    PATHS$intermediate_dir,
    "tracerx100_unconstrained_patient_summary_protein_altering.rds"
  )
)

fwrite(
  subset_results,
  file.path(
    tracerx_result_dir,
    "11_tracerx100_exhaustive_subsets_protein_altering.tsv"
  ),
  sep = "\t"
)

fwrite(
  patient_summary,
  file.path(
    tracerx_result_dir,
    "11_tracerx100_unconstrained_patient_summary_protein_altering.tsv"
  ),
  sep = "\t"
)

if (DEBUG) {
  fwrite(
    n_distribution_qc,
    file.path(
      PATHS$qc_dir,
      "11_tracerx100_region_count_distribution_qc.tsv"
    ),
    sep = "\t"
  )
}

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
        classification_abs_error,
        apparent_ubiquity_error_mean,
        analytic_apparent_ubiquity_error,
        apparent_ubiquity_abs_error,
        recall_private_mean,
        analytic_recall_private,
        private_abs_error
      )
    ],
    file.path(
      PATHS$qc_dir,
      "11_tracerx100_analytic_expectation_validation.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    regression,
    file.path(
      PATHS$qc_dir,
      "11_tracerx100_script10_reference_regression.tsv"
    ),
    sep = "\t"
  )
}

# Descriptive N(k) table only; no biological inference.
available_n_by_k <- patient_summary[
  ,
  .(
    n_patients = uniqueN(patient_id),
    n_metric_detection =
      sum(!is.na(recall_nonubiquitous_detection_mean))
  ),
  by = k
]
setorder(available_n_by_k, k)

strict_n_by_k <- rbindlist(
  lapply(
    sort(unique(patient_summary$k)),
    function(kval) {
      d <- patient_summary[k == kval & n > kval]
      data.table(
        k = kval,
        n_patients = uniqueN(d$patient_id),
        n_metric_detection =
          sum(!is.na(d$recall_nonubiquitous_detection_mean))
      )
    }
  )
)

if (DEBUG) {
  fwrite(
    available_n_by_k,
    file.path(
      PATHS$qc_dir,
      "11_tracerx100_available_N_by_k.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    strict_n_by_k,
    file.path(
      PATHS$qc_dir,
      "11_tracerx100_strict_N_by_k.tsv"
    ),
    sep = "\t"
  )
}

# ------------------------------------------------------------------------------
# 10. Human-readable QC summary
# ------------------------------------------------------------------------------

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  "Script 11: TRACERx100 exhaustive sampling",
  "======================================================",
  "",
  sprintf(
    "Patients analysed: %d",
    uniqueN(patient_summary$patient_id)
  ),
  sprintf(
    "Primary tumour regions: %d",
    sum(unique(patient_table[, .(patient_id, n)])$n)
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
  "  patient-level MEAN across all equally likely k-region subsets",
  "Legacy/typical-subset sensitivity:",
  "  patient-level MEDIAN across all k-region subsets",
  "",
  "QC:",
  sprintf(
    "  exact subset-count QC: PASS (%d total subsets)",
    EXPECTED_TOTAL_SUBSETS
  ),
  "  Script 10 full-reference count regression: PASS",
  "  zero-denominator NA propagation: PASS",
  "  patient-level monotonicity: PASS",
  "  k=n full-reference identities: PASS",
  "  analytic hypergeometric validation: PASS",
  sprintf(
    "    max detection-recall error: %.3g",
    max_detection_error
  ),
  sprintf(
    "    max heterogeneity-classification error: %.3g",
    max_classification_error
  ),
  sprintf(
    "    max apparent-ubiquity error: %.3g",
    max_apparent_error
  ),
  sprintf(
    "    max private-recall error: %.3g",
    max_private_error
  ),
  "",
  "Zero-denominator patients retained:",
  paste0(
    "  ",
    paste(EXPECTED_ZERO_NONUBIQ, collapse = ", ")
  ),
  "  detection/classification/apparent-ubiquity/private recall remain NA.",
  "",
  "Available-cohort N by k:",
  paste0(
    "  k=", available_n_by_k$k,
    ": N=", available_n_by_k$n_patients,
    ", defined detection N=",
    available_n_by_k$n_metric_detection,
    collapse = "\n"
  ),
  "",
  "Strict true-downsampling N by k (n>k):",
  paste0(
    "  k=", strict_n_by_k$k,
    ": N=", strict_n_by_k$n_patients,
    ", defined detection N=",
    strict_n_by_k$n_metric_detection,
    collapse = "\n"
  ),
  "",
  "No cohort-level learning-curve inference was performed in Script 11.",
  "Proceed to TRACERx100/HCC external-validation analysis only if all QC PASS.",
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed)
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "11_tracerx100_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
