# ==============================================================================
# 12_hcc_tracerx100_external_validation_exact_bootstrap.R
#
# Purpose
#   Perform the prespecified cohort-level external validation of the HCC
#   protein-altering sampling-depth learning curve in TRACERx100.
#
# IMPORTANT SCIENTIFIC FRAMING
#   This is NOT a search for a new "optimal k" in TRACERx.
#   It asks whether the qualitative sampling-depth behavior observed in HCC
#   persists in an independent multiregion NSCLC cohort generated with WES.
#
# Primary cross-cohort comparison
#   SAME TRUE-DOWNSAMPLING DESIGN:
#     HCC       : fixed n >= 5, N = 55, k = 1..4
#     TRACERx100: fixed n >= 5, N = 15, k = 1..4
#
#   SAME event concept:
#     protein-altering exact variants
#
#   SAME sampling estimand:
#     within-patient MEAN across all exhaustive k-region subsets
#
#   SAME main endpoints:
#     1. non-ubiquitous detection recall
#     2. heterogeneity-classification recall
#
# Cohort-level estimator
#   Median across patients, with exact patient-level nonparametric percentile-bootstrap CI.
#
# Cross-cohort contrasts
#   Descriptive difference in medians:
#       TRACERx100 - HCC
#   using the exact independent empirical bootstrap distributions of the two cohort medians.
#
# Marginal gains
#   Within each fixed cohort:
#       R(k) - R(k-1), paired within patient
#   Cross-cohort difference compares the medians of these patient-level gains.
#
# Operational 80% threshold
#   Retained ONLY for continuity with the legacy HCC analysis.
#   It is not used as a biological optimum or as a formal validation criterion.
#
# Reference-architecture context
#   Because the cohorts differ biologically and technologically (WGS HCC versus
#   WES TRACERx100), we also compare full-reference protein-altering event
#   architecture in the fixed n>=5 cohorts:
#     - total protein-altering event burden
#     - non-ubiquitous fraction
#     - private-event fraction
#
# No p-values, equivalence tests, or post-hoc optimization are performed.
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
source(file.path(script_dir, "functions_bootstrap_exact.R"))

cat("Script 12: HCC / TRACERx100 external validation\n")
cat("==========================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
cross_dir <- file.path(PATHS$output_dir, "results", "cross_cohort")
dir.create(cross_dir, recursive = TRUE, showWarnings = FALSE)

# Canonical patient-level sampling outputs.
hcc_patient_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_unconstrained_patient_summary_protein_altering.rds"
)
trx_patient_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_unconstrained_patient_summary_protein_altering.rds"
)

# Full-reference architecture sources.
hcc_reference_summary_path <- file.path(
  PATHS$qc_dir,
  "02_reference_summary_by_patient.tsv"
)
trx_reference_summary_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_patient_reference_summary_protein_altering.rds"
)

# Regression target for the already validated HCC protein-altering summary.
hcc_depth_summary_path <- file.path(
  PATHS$output_dir,
  "results",
  "hcc_sensitivity",
  "08_hcc_protein_altering_depth_summary.tsv"
)

required_paths <- c(
  hcc_patient_path,
  trx_patient_path,
  hcc_reference_summary_path,
  trx_reference_summary_path,
  hcc_depth_summary_path
)

missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing validated prerequisite file(s):\n",
    paste0("  ", missing_paths, collapse = "\n")
  )
}

# ------------------------------- constants -----------------------------------
B <- NA_integer_  # exact bootstrap: no Monte-Carlo replicates
BASE_SEED <- as.integer(PARAMS$seed)
CI_LEVEL <- 0.95
THRESHOLD <- as.numeric(PARAMS$adequate_recovery_threshold)

FIXED_K <- 1:4
GAIN_TO_K <- 2:4

EXPECTED_HCC_FIXED_N <- 55L
EXPECTED_TRX_FIXED_N <- 15L
EXPECTED_HCC_N_GE6 <- 7L
EXPECTED_TRX_N_GE6 <- 6L

EXPECTED_TRX_AVAILABLE_N <- c(
  `1` = 100L,
  `2` = 100L,
  `3` = 60L,
  `4` = 37L,
  `5` = 15L,
  `6` = 6L,
  `7` = 5L
)

EXPECTED_TRX_STRICT_N <- c(
  `1` = 100L,
  `2` = 60L,
  `3` = 37L,
  `4` = 15L,
  `5` = 6L,
  `6` = 5L
)

MAIN_METRICS <- data.table(
  metric = c(
    "nonubiquitous_detection_recall",
    "heterogeneity_classification_recall"
  ),
  column = c(
    "recall_nonubiquitous_detection_mean",
    "recall_heterogeneity_classification_mean"
  ),
  role = c(
    "primary",
    "key_secondary"
  )
)

SUPPORTIVE_METRICS <- data.table(
  metric = c(
    "apparent_ubiquity_error",
    "private_recall",
    "nonubiquitous_fraction_abs_error",
    "jaccard_ith_abs_error"
  ),
  column = c(
    "apparent_ubiquity_error_mean",
    "recall_private_mean",
    "nonubiquitous_fraction_abs_error_mean",
    "jaccard_ith_abs_error_mean"
  )
)

ALL_METRICS <- rbindlist(
  list(
    MAIN_METRICS[, .(metric, column)],
    SUPPORTIVE_METRICS[, .(metric, column)]
  )
)

# ------------------------------- helpers -------------------------------------
safe_mean <- function(v) {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  mean(v)
}

safe_median <- function(v) {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  stats::median(v)
}

safe_quantile <- function(v, p) {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  as.numeric(
    stats::quantile(
      v,
      probs = p,
      names = FALSE,
      type = 7
    )
  )
}

bootstrap_median_ci <- function(v, B, seed, level = 0.95) {
  exact_bootstrap_median_ci(
    v = v,
    level = level,
    min_n = 1L
  )
}

bootstrap_independent_median_difference <- function(
    x_trx,
    x_hcc,
    B,
    seed,
    level = 0.95) {

  exact_bootstrap_independent_median_difference_ci(
    x = x_trx,
    y = x_hcc,
    level = level
  )
}

summarise_curve <- function(
    d,
    cohort,
    analysis,
    cohort_definition,
    k_values,
    selector,
    seed_offset) {

  rows <- list()
  ri <- 0L

  for (kval in k_values) {
    keep <- selector(d, kval)
    dk <- d[which(keep)]

    if (nrow(dk) == 0L) {
      stop(
        "No rows for ", cohort,
        " / ", analysis,
        " / k=", kval, "."
      )
    }

    if (any(dk$k != kval) ||
        anyDuplicated(dk$patient_id) > 0L) {
      stop(
        "Invalid patient x k structure for ",
        cohort, " / ", analysis,
        " / k=", kval, "."
      )
    }

    for (mi in seq_len(nrow(ALL_METRICS))) {
      metric_id <- ALL_METRICS$metric[[mi]]
      col_name <- ALL_METRICS$column[[mi]]

      if (!col_name %in% names(dk)) {
        stop(
          "Required metric column absent: ",
          col_name
        )
      }

      v <- dk[[col_name]]
      ci <- bootstrap_median_ci(
        v,
        B = B,
        seed = BASE_SEED +
          seed_offset +
          kval * 100L +
          mi,
        level = CI_LEVEL
      )
      ci_defined <- all(is.finite(ci))

      ri <- ri + 1L

      rows[[ri]] <- data.table(
        cohort = cohort,
        analysis = analysis,
        cohort_definition = cohort_definition,
        k = as.integer(kval),
        n_patients = nrow(dk),
        n_metric = sum(!is.na(v)),
        n_equal_k = sum(dk$n == kval),
        fraction_equal_k = mean(dk$n == kval),
        metric = metric_id,
        patient_subset_estimand = "mean",
        cohort_estimand = "median",
        median = safe_median(v),
        q1 = safe_quantile(v, 0.25),
        q3 = safe_quantile(v, 0.75),
        bootstrap_ci_lower = ci[[1L]],
        bootstrap_ci_upper = ci[[2L]],
        bootstrap_replicates = NA_integer_,
        bootstrap_method = if (ci_defined) {
          "exact_nonparametric_percentile"
        } else {
          NA_character_
        },
        bootstrap_level = if (ci_defined) CI_LEVEL else NA_real_
      )
    }
  }

  rbindlist(rows)
}

make_fixed_gain_table <- function(
    d,
    cohort,
    n_min,
    seed_offset) {

  fixed <- d[n >= n_min & k %in% FIXED_K]

  patient_sets <- lapply(
    FIXED_K,
    function(kval) {
      sort(fixed[k == kval, patient_id])
    }
  )

  if (!all(vapply(
    patient_sets[-1L],
    function(z) identical(z, patient_sets[[1L]]),
    logical(1)
  ))) {
    stop(
      cohort,
      ": fixed n>=", n_min,
      " patient set changes across k."
    )
  }

  rows <- list()
  ri <- 0L

  for (mi in seq_len(nrow(MAIN_METRICS))) {
    metric_id <- MAIN_METRICS$metric[[mi]]
    col_name <- MAIN_METRICS$column[[mi]]

    wide <- dcast(
      fixed[
        ,
        .(
          patient_id,
          k,
          value = get(col_name)
        )
      ],
      patient_id ~ k,
      value.var = "value"
    )

    if (!all(as.character(FIXED_K) %in% names(wide))) {
      stop(
        cohort,
        ": incomplete fixed-cohort gain matrix."
      )
    }

    for (k_to in GAIN_TO_K) {
      k_from <- k_to - 1L

      delta <- (
        wide[[as.character(k_to)]] -
          wide[[as.character(k_from)]]
      )

      ci <- bootstrap_median_ci(
        delta,
        B = B,
        seed = BASE_SEED +
          seed_offset +
          mi * 1000L +
          k_to,
        level = CI_LEVEL
      )

      ri <- ri + 1L

      rows[[ri]] <- data.table(
        cohort = cohort,
        analysis = paste0("fixed_n_ge", n_min),
        metric = metric_id,
        k_from = k_from,
        k_to = k_to,
        n_patients = length(delta),
        median_patient_gain = safe_median(delta),
        q1_patient_gain = safe_quantile(delta, 0.25),
        q3_patient_gain = safe_quantile(delta, 0.75),
        mean_patient_gain = safe_mean(delta),
        bootstrap_ci_lower = ci[[1L]],
        bootstrap_ci_upper = ci[[2L]],
        bootstrap_replicates = NA_integer_,
        bootstrap_method = "exact_nonparametric_percentile",
        bootstrap_level = CI_LEVEL
      )
    }
  }

  rbindlist(rows)
}

# ------------------------------- 1. load --------------------------------------
hcc <- as.data.table(readRDS(hcc_patient_path))
trx <- as.data.table(readRDS(trx_patient_path))

hcc_ref <- fread(hcc_reference_summary_path)
trx_ref <- as.data.table(readRDS(trx_reference_summary_path))
hcc_depth_validated <- fread(hcc_depth_summary_path)

required_patient_cols <- c(
  "patient_id", "n", "k",
  ALL_METRICS$column
)

if (!all(required_patient_cols %in% names(hcc))) {
  stop(
    "HCC protein patient summary schema mismatch."
  )
}
if (!all(required_patient_cols %in% names(trx))) {
  stop(
    "TRACERx100 patient summary schema mismatch."
  )
}

hcc[
  ,
  `:=`(
    patient_id = as.character(patient_id),
    n = as.integer(n),
    k = as.integer(k)
  )
]

trx[
  ,
  `:=`(
    patient_id = as.character(patient_id),
    n = as.integer(n),
    k = as.integer(k)
  )
]

# ----------------------- 2. fixed-cohort structural QC ------------------------
hcc_fixed_ids <- sort(
  hcc[n >= 5L & k == 1L, patient_id]
)
trx_fixed_ids <- sort(
  trx[n >= 5L & k == 1L, patient_id]
)

if (length(hcc_fixed_ids) != EXPECTED_HCC_FIXED_N) {
  stop("HCC fixed n>=5 N QC FAILED.")
}
if (length(trx_fixed_ids) != EXPECTED_TRX_FIXED_N) {
  stop("TRACERx100 fixed n>=5 N QC FAILED.")
}

for (kval in FIXED_K) {
  if (!identical(
    sort(hcc[n >= 5L & k == kval, patient_id]),
    hcc_fixed_ids
  )) {
    stop(
      "HCC fixed n>=5 patient set changes at k=",
      kval, "."
    )
  }

  if (!identical(
    sort(trx[n >= 5L & k == kval, patient_id]),
    trx_fixed_ids
  )) {
    stop(
      "TRACERx100 fixed n>=5 patient set changes at k=",
      kval, "."
    )
  }
}

# Every k<=4 fixed-cohort observation must be true downsampling.
if (hcc[n >= 5L & k %in% FIXED_K, any(k >= n)]) {
  stop("HCC fixed n>=5 contains k=n at k<=4.")
}
if (trx[n >= 5L & k %in% FIXED_K, any(k >= n)]) {
  stop("TRACERx100 fixed n>=5 contains k=n at k<=4.")
}

# Main metrics must be defined for every fixed-cohort patient.
for (cc in MAIN_METRICS$column) {
  if (hcc[n >= 5L & k %in% FIXED_K, any(is.na(get(cc)))]) {
    stop("HCC fixed-cohort main endpoint contains NA.")
  }

  if (trx[n >= 5L & k %in% FIXED_K, any(is.na(get(cc)))]) {
    stop("TRACERx100 fixed-cohort main endpoint contains NA.")
  }
}

# -------------------- 3. TRACERx cohort-view summaries ------------------------
trx_available <- summarise_curve(
  d = trx,
  cohort = "TRACERx100",
  analysis = "available_cohort",
  cohort_definition =
    "all patients with n>=k; descriptive and ceiling-sensitive",
  k_values = 1:7,
  selector = function(d, kval) {
    d$k == kval
  },
  seed_offset = 10000L
)

trx_strict <- summarise_curve(
  d = trx,
  cohort = "TRACERx100",
  analysis = "strict_downsampling",
  cohort_definition =
    "patients with n>k; true downsampling, cohort changes with k",
  k_values = 1:6,
  selector = function(d, kval) {
    d$k == kval & d$n > kval
  },
  seed_offset = 20000L
)

trx_fixed5 <- summarise_curve(
  d = trx,
  cohort = "TRACERx100",
  analysis = "fixed_n_ge5",
  cohort_definition =
    "same 15 patients with n>=5 at k=1..4; primary external-validation view",
  k_values = FIXED_K,
  selector = function(d, kval) {
    d$k == kval & d$n >= 5L
  },
  seed_offset = 30000L
)

trx_nge6 <- summarise_curve(
  d = trx,
  cohort = "TRACERx100",
  analysis = "n_ge6_exploratory",
  cohort_definition =
    "same 6 patients with n>=6 at k=1..5; exploratory deeper-sampling view",
  k_values = 1:5,
  selector = function(d, kval) {
    d$k == kval & d$n >= 6L
  },
  seed_offset = 40000L
)

trx_curve_summary <- rbindlist(
  list(
    trx_available,
    trx_strict,
    trx_fixed5,
    trx_nge6
  ),
  use.names = TRUE
)

# N(k) regression against Script 11.
trx_available_N <- trx_curve_summary[
  analysis == "available_cohort" &
    metric == "nonubiquitous_detection_recall",
  .(
    k,
    observed_N = n_patients,
    defined_N = n_metric
  )
]
setorder(trx_available_N, k)

if (!identical(
  as.integer(trx_available_N$observed_N),
  as.integer(EXPECTED_TRX_AVAILABLE_N)
)) {
  stop("TRACERx100 available N(k) regression FAILED.")
}

trx_strict_N <- trx_curve_summary[
  analysis == "strict_downsampling" &
    metric == "nonubiquitous_detection_recall",
  .(
    k,
    observed_N = n_patients,
    defined_N = n_metric
  )
]
setorder(trx_strict_N, k)

if (!identical(
  as.integer(trx_strict_N$observed_N),
  as.integer(EXPECTED_TRX_STRICT_N)
)) {
  stop("TRACERx100 strict N(k) regression FAILED.")
}

# --------------------- 4. HCC fixed-cohort reconstruction ---------------------
hcc_fixed5 <- summarise_curve(
  d = hcc,
  cohort = "HCC",
  analysis = "fixed_n_ge5",
  cohort_definition =
    "same 55 patients with n>=5 at k=1..4; HCC protein-altering comparator",
  k_values = FIXED_K,
  selector = function(d, kval) {
    d$k == kval & d$n >= 5L
  },
  seed_offset = 50000L
)

# Numerical regression of reconstructed HCC fixed-n>=5 medians against the
# already validated Script 08 summary.
hcc_regression_target <- hcc_depth_validated[
  analysis == "fixed_n_ge5" &
    metric %in% MAIN_METRICS$metric,
  .(
    k,
    metric,
    validated_median = median
  )
]

hcc_regression_observed <- hcc_fixed5[
  metric %in% MAIN_METRICS$metric,
  .(
    k,
    metric,
    reconstructed_median = median
  )
]

hcc_regression <- merge(
  hcc_regression_target,
  hcc_regression_observed,
  by = c("k", "metric"),
  all = TRUE,
  sort = TRUE
)

hcc_regression[
  ,
  abs_error := abs(
    reconstructed_median -
      validated_median
  )
]

if (nrow(hcc_regression) != 8L ||
    anyNA(hcc_regression$abs_error) ||
    any(hcc_regression$abs_error > 1e-12)) {
  if (DEBUG) {
    fwrite(
      hcc_regression,
      file.path(
        PATHS$qc_dir,
        "12_hcc_fixed_curve_regression_failure.tsv"
      ),
      sep = "\t"
    )
  }
  stop("HCC Script 08 fixed-curve regression FAILED.")
}

# ------------------ 5. fixed n>=5 cross-cohort comparison --------------------
fixed_curve_both <- rbindlist(
  list(
    hcc_fixed5,
    trx_fixed5
  ),
  use.names = TRUE
)

curve_comparison_rows <- list()
cc_i <- 0L

for (mi in seq_len(nrow(MAIN_METRICS))) {
  metric_id <- MAIN_METRICS$metric[[mi]]
  col_name <- MAIN_METRICS$column[[mi]]

  for (kval in FIXED_K) {
    hv <- hcc[
      n >= 5L & k == kval,
      get(col_name)
    ]
    tv <- trx[
      n >= 5L & k == kval,
      get(col_name)
    ]

    diff_ci <- bootstrap_independent_median_difference(
      x_trx = tv,
      x_hcc = hv,
      B = B,
      seed = BASE_SEED +
        60000L +
        mi * 1000L +
        kval,
      level = CI_LEVEL
    )

    cc_i <- cc_i + 1L

    curve_comparison_rows[[cc_i]] <- data.table(
      metric = metric_id,
      endpoint_role = MAIN_METRICS$role[[mi]],
      k = as.integer(kval),

      hcc_N = length(hv),
      hcc_median = safe_median(hv),
      hcc_q1 = safe_quantile(hv, 0.25),
      hcc_q3 = safe_quantile(hv, 0.75),

      tracerx_N = length(tv),
      tracerx_median = safe_median(tv),
      tracerx_q1 = safe_quantile(tv, 0.25),
      tracerx_q3 = safe_quantile(tv, 0.75),

      tracerx_minus_hcc_median =
        safe_median(tv) - safe_median(hv),

      difference_bootstrap_ci_lower =
        diff_ci[[1L]],
      difference_bootstrap_ci_upper =
        diff_ci[[2L]],

      bootstrap_replicates = NA_integer_,
      bootstrap_method = "exact_nonparametric_percentile",
      bootstrap_level = CI_LEVEL,

      interpretation =
        "descriptive external-validation contrast; not equivalence testing"
    )
  }
}

fixed_curve_comparison <- rbindlist(curve_comparison_rows)

# -------------------------- 6. paired marginal gains --------------------------
hcc_gains <- make_fixed_gain_table(
  hcc,
  cohort = "HCC",
  n_min = 5L,
  seed_offset = 70000L
)

trx_gains <- make_fixed_gain_table(
  trx,
  cohort = "TRACERx100",
  n_min = 5L,
  seed_offset = 80000L
)

gain_comparison_rows <- list()
gc_i <- 0L

for (mi in seq_len(nrow(MAIN_METRICS))) {
  metric_id <- MAIN_METRICS$metric[[mi]]
  col_name <- MAIN_METRICS$column[[mi]]

  hwide <- dcast(
    hcc[
      n >= 5L & k %in% FIXED_K,
      .(
        patient_id,
        k,
        value = get(col_name)
      )
    ],
    patient_id ~ k,
    value.var = "value"
  )

  twide <- dcast(
    trx[
      n >= 5L & k %in% FIXED_K,
      .(
        patient_id,
        k,
        value = get(col_name)
      )
    ],
    patient_id ~ k,
    value.var = "value"
  )

  for (k_to in GAIN_TO_K) {
    k_from <- k_to - 1L

    hg <- (
      hwide[[as.character(k_to)]] -
        hwide[[as.character(k_from)]]
    )

    tg <- (
      twide[[as.character(k_to)]] -
        twide[[as.character(k_from)]]
    )

    diff_ci <- bootstrap_independent_median_difference(
      x_trx = tg,
      x_hcc = hg,
      B = B,
      seed = BASE_SEED +
        90000L +
        mi * 1000L +
        k_to,
      level = CI_LEVEL
    )

    gc_i <- gc_i + 1L

    gain_comparison_rows[[gc_i]] <- data.table(
      metric = metric_id,
      endpoint_role = MAIN_METRICS$role[[mi]],
      k_from = k_from,
      k_to = k_to,

      hcc_N = length(hg),
      hcc_median_gain = safe_median(hg),
      hcc_q1_gain = safe_quantile(hg, 0.25),
      hcc_q3_gain = safe_quantile(hg, 0.75),

      tracerx_N = length(tg),
      tracerx_median_gain = safe_median(tg),
      tracerx_q1_gain = safe_quantile(tg, 0.25),
      tracerx_q3_gain = safe_quantile(tg, 0.75),

      tracerx_minus_hcc_median_gain =
        safe_median(tg) - safe_median(hg),

      difference_bootstrap_ci_lower =
        diff_ci[[1L]],
      difference_bootstrap_ci_upper =
        diff_ci[[2L]],

      bootstrap_replicates = NA_integer_,
      bootstrap_method = "exact_nonparametric_percentile",
      bootstrap_level = CI_LEVEL
    )
  }
}

gain_comparison <- rbindlist(gain_comparison_rows)

# ---------------------- 7. operational threshold comparison ------------------
threshold_rows <- list()
th_i <- 0L

for (cohort_id in c("HCC", "TRACERx100")) {
  d <- if (cohort_id == "HCC") hcc else trx

  for (kval in FIXED_K) {
    v <- d[
      n >= 5L & k == kval,
      recall_nonubiquitous_detection_mean
    ]

    th_i <- th_i + 1L

    threshold_rows[[th_i]] <- data.table(
      cohort = cohort_id,
      analysis = "fixed_n_ge5",
      k = kval,
      n_patients = length(v),
      median_detection_recall = safe_median(v),
      threshold = THRESHOLD,
      median_reaches_threshold =
        safe_median(v) >= THRESHOLD,
      n_patients_ge_threshold =
        sum(v >= THRESHOLD, na.rm = TRUE),
      proportion_patients_ge_threshold =
        mean(v >= THRESHOLD, na.rm = TRUE)
    )
  }
}

threshold_comparison <- rbindlist(threshold_rows)

first_threshold <- threshold_comparison[
  ,
  {
    kk <- k[median_reaches_threshold %in% TRUE]

    list(
      first_k_with_median_ge_threshold =
        if (length(kk) == 0L) {
          NA_integer_
        } else {
          min(kk)
        }
    )
  },
  by = cohort
]

threshold_comparison <- merge(
  threshold_comparison,
  first_threshold,
  by = "cohort",
  all.x = TRUE,
  sort = FALSE
)

# -------------------- 8. reference-architecture comparison -------------------
if (!all(c(
  "event_set",
  "patient_id",
  "n_sectors",
  "n_events_reference",
  "n_nonubiquitous_reference",
  "n_private_reference"
) %in% names(hcc_ref))) {
  stop("HCC reference-summary schema mismatch.")
}

if (!all(c(
  "patient_id",
  "n_regions",
  "n_events_reference",
  "n_nonubiquitous_reference",
  "n_private_reference"
) %in% names(trx_ref))) {
  stop("TRACERx100 reference-summary schema mismatch.")
}

hcc_ref_fixed <- hcc_ref[
  event_set == "protein_altering" &
    n_sectors >= 5L
]

trx_ref_fixed <- trx_ref[
  n_regions >= 5L
]

if (nrow(hcc_ref_fixed) != EXPECTED_HCC_FIXED_N) {
  stop("HCC fixed reference-architecture N QC FAILED.")
}
if (nrow(trx_ref_fixed) != EXPECTED_TRX_FIXED_N) {
  stop("TRACERx fixed reference-architecture N QC FAILED.")
}

hcc_ref_fixed[
  ,
  `:=`(
    nonubiquitous_fraction =
      n_nonubiquitous_reference /
        n_events_reference,
    private_fraction =
      n_private_reference /
        n_events_reference
  )
]

trx_ref_fixed[
  ,
  `:=`(
    nonubiquitous_fraction =
      n_nonubiquitous_reference /
        n_events_reference,
    private_fraction =
      n_private_reference /
        n_events_reference
  )
]

reference_variables <- c(
  "n_events_reference",
  "n_nonubiquitous_reference",
  "nonubiquitous_fraction",
  "private_fraction"
)

reference_architecture_rows <- list()
ra_i <- 0L

for (variable in reference_variables) {
  hv <- hcc_ref_fixed[[variable]]
  tv <- trx_ref_fixed[[variable]]

  diff_ci <- bootstrap_independent_median_difference(
    x_trx = tv,
    x_hcc = hv,
    B = B,
    seed = BASE_SEED +
      100000L +
      match(variable, reference_variables),
    level = CI_LEVEL
  )

  ra_i <- ra_i + 1L

  reference_architecture_rows[[ra_i]] <- data.table(
    variable = variable,

    hcc_N = length(hv),
    hcc_median = safe_median(hv),
    hcc_q1 = safe_quantile(hv, 0.25),
    hcc_q3 = safe_quantile(hv, 0.75),

    tracerx_N = length(tv),
    tracerx_median = safe_median(tv),
    tracerx_q1 = safe_quantile(tv, 0.25),
    tracerx_q3 = safe_quantile(tv, 0.75),

    tracerx_minus_hcc_median =
      safe_median(tv) - safe_median(hv),

    difference_bootstrap_ci_lower =
      diff_ci[[1L]],
    difference_bootstrap_ci_upper =
      diff_ci[[2L]],

    interpretation =
      "contextual cohort architecture; biology and assay differences are not separable"
  )
}

reference_architecture <- rbindlist(
  reference_architecture_rows
)

# ----------------------- 9. n>=6 k=5 exploratory check -----------------------
hcc_nge6_ids <- hcc[n >= 6L & k == 1L, patient_id]
trx_nge6_ids <- trx[n >= 6L & k == 1L, patient_id]

if (length(hcc_nge6_ids) != EXPECTED_HCC_N_GE6) {
  stop("HCC n>=6 exploratory N QC FAILED.")
}
if (length(trx_nge6_ids) != EXPECTED_TRX_N_GE6) {
  stop("TRACERx n>=6 exploratory N QC FAILED.")
}

k5_exploratory_rows <- list()
ke_i <- 0L

for (mi in seq_len(nrow(MAIN_METRICS))) {
  metric_id <- MAIN_METRICS$metric[[mi]]
  col_name <- MAIN_METRICS$column[[mi]]

  hv <- hcc[
    n >= 6L & k == 5L,
    get(col_name)
  ]

  tv <- trx[
    n >= 6L & k == 5L,
    get(col_name)
  ]

  ke_i <- ke_i + 1L

  k5_exploratory_rows[[ke_i]] <- data.table(
    metric = metric_id,
    hcc_N = length(hv),
    hcc_median = safe_median(hv),
    hcc_q1 = safe_quantile(hv, 0.25),
    hcc_q3 = safe_quantile(hv, 0.75),
    tracerx_N = length(tv),
    tracerx_median = safe_median(tv),
    tracerx_q1 = safe_quantile(tv, 0.25),
    tracerx_q3 = safe_quantile(tv, 0.75),
    tracerx_minus_hcc_median =
      safe_median(tv) - safe_median(hv),
    inferential_status =
      "exploratory_descriptive_only_small_N"
  )
}

k5_exploratory <- rbindlist(k5_exploratory_rows)

# ------------------------------ 10. outputs -----------------------------------
if (DEBUG) {
  fwrite(
    trx_curve_summary,
    file.path(
      cross_dir,
      "12_tracerx100_depth_summary_protein_altering.tsv"
    ),
    sep = "\t"
  )
}

fwrite(
  fixed_curve_both,
  file.path(
    cross_dir,
    "12_hcc_tracerx100_fixed_n_ge5_curve_summary.tsv"
  ),
  sep = "\t"
)

if (DEBUG) {
  fwrite(
    fixed_curve_comparison,
    file.path(
      cross_dir,
      "12_hcc_tracerx100_fixed_n_ge5_curve_comparison.tsv"
    ),
    sep = "\t"
  )
}

fwrite(
  rbindlist(
    list(hcc_gains, trx_gains),
    use.names = TRUE
  ),
  file.path(
    cross_dir,
    "12_hcc_tracerx100_fixed_n_ge5_paired_gains.tsv"
  ),
  sep = "\t"
)

if (DEBUG) {
  fwrite(
    gain_comparison,
    file.path(
      cross_dir,
      "12_hcc_tracerx100_fixed_n_ge5_gain_comparison.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    threshold_comparison,
    file.path(
      cross_dir,
      "12_hcc_tracerx100_operational_threshold_comparison.tsv"
    ),
    sep = "\t"
  )
}

fwrite(
  reference_architecture,
  file.path(
    cross_dir,
    "12_hcc_tracerx100_reference_architecture_comparison.tsv"
  ),
  sep = "\t"
)

if (DEBUG) {
  fwrite(
    k5_exploratory,
    file.path(
      cross_dir,
      "12_hcc_tracerx100_n_ge6_k5_exploratory.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    hcc_regression,
    file.path(
      PATHS$qc_dir,
      "12_hcc_script08_fixed_curve_regression.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    trx_available_N,
    file.path(
      PATHS$qc_dir,
      "12_tracerx100_available_N_regression.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    trx_strict_N,
    file.path(
      PATHS$qc_dir,
      "12_tracerx100_strict_N_regression.tsv"
    ),
    sep = "\t"
  )
}


# ------------------------ 11. human-readable summary ---------------------------
get_curve_value <- function(cohort_id, metric_id, kval) {
  z <- fixed_curve_both[
    cohort == cohort_id &
      metric == metric_id &
      k == kval
  ]

  if (nrow(z) != 1L) {
    stop(
      "Could not retrieve unique curve row for ",
      cohort_id, " / ", metric_id,
      " / k=", kval, "."
    )
  }

  z
}

get_gain_value <- function(cohort_id, metric_id, k_to_value) {
  z <- rbindlist(
    list(hcc_gains, trx_gains),
    use.names = TRUE
  )[
    cohort == cohort_id &
      metric == metric_id &
      k_to == k_to_value
  ]

  if (nrow(z) != 1L) {
    stop(
      "Could not retrieve unique gain row for ",
      cohort_id, " / ", metric_id,
      " / k_to=", k_to_value, "."
    )
  }

  z
}

h4d <- get_curve_value(
  "HCC",
  "nonubiquitous_detection_recall",
  4L
)
t4d <- get_curve_value(
  "TRACERx100",
  "nonubiquitous_detection_recall",
  4L
)
h4c <- get_curve_value(
  "HCC",
  "heterogeneity_classification_recall",
  4L
)
t4c <- get_curve_value(
  "TRACERx100",
  "heterogeneity_classification_recall",
  4L
)

hcc_first_threshold <- unique(
  threshold_comparison[
    cohort == "HCC",
    first_k_with_median_ge_threshold
  ]
)
trx_first_threshold <- unique(
  threshold_comparison[
    cohort == "TRACERx100",
    first_k_with_median_ge_threshold
  ]
)

hcc_ref_nonfrac <- reference_architecture[
  variable == "nonubiquitous_fraction"
]
trx_ref_priv <- reference_architecture[
  variable == "private_fraction"
]

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  "Script 12: HCC / TRACERx100 external validation",
  "==========================================================",
  "",
  "Primary fixed true-downsampling comparison:",
  sprintf(
    "  HCC: N=%d patients with n>=5, k=1..4",
    EXPECTED_HCC_FIXED_N
  ),
  sprintf(
    "  TRACERx100: N=%d patients with n>=5, k=1..4",
    EXPECTED_TRX_FIXED_N
  ),
  "",
  "Protein-altering non-ubiquitous detection recall:",
  paste0(
    "  HCC medians k=1..4: ",
    paste(
      sprintf(
        "%.4f",
        fixed_curve_both[
          cohort == "HCC" &
            metric == "nonubiquitous_detection_recall",
          median
        ]
      ),
      collapse = ", "
    )
  ),
  paste0(
    "  TRACERx100 medians k=1..4: ",
    paste(
      sprintf(
        "%.4f",
        fixed_curve_both[
          cohort == "TRACERx100" &
            metric == "nonubiquitous_detection_recall",
          median
        ]
      ),
      collapse = ", "
    )
  ),
  "",
  "Protein-altering heterogeneity-classification recall:",
  paste0(
    "  HCC medians k=1..4: ",
    paste(
      sprintf(
        "%.4f",
        fixed_curve_both[
          cohort == "HCC" &
            metric == "heterogeneity_classification_recall",
          median
        ]
      ),
      collapse = ", "
    )
  ),
  paste0(
    "  TRACERx100 medians k=1..4: ",
    paste(
      sprintf(
        "%.4f",
        fixed_curve_both[
          cohort == "TRACERx100" &
            metric == "heterogeneity_classification_recall",
          median
        ]
      ),
      collapse = ", "
    )
  ),
  "",
  "At k=4:",
  sprintf(
    "  detection: HCC %.4f vs TRACERx100 %.4f (TRACERx-HCC = %.4f)",
    h4d$median,
    t4d$median,
    t4d$median - h4d$median
  ),
  sprintf(
    "  classification: HCC %.4f vs TRACERx100 %.4f (TRACERx-HCC = %.4f)",
    h4c$median,
    t4c$median,
    t4c$median - h4c$median
  ),
  "",
  "Paired marginal detection gains in fixed n>=5 cohorts:",
  paste0(
    "  HCC: ",
    paste(
      sprintf(
        "%.4f",
        hcc_gains[
          metric == "nonubiquitous_detection_recall",
          median_patient_gain
        ]
      ),
      collapse = ", "
    ),
    " for 1->2, 2->3, 3->4"
  ),
  paste0(
    "  TRACERx100: ",
    paste(
      sprintf(
        "%.4f",
        trx_gains[
          metric == "nonubiquitous_detection_recall",
          median_patient_gain
        ]
      ),
      collapse = ", "
    ),
    " for 1->2, 2->3, 3->4"
  ),
  "",
  sprintf(
    "Operational %.0f%% threshold first reached by cohort median:",
    100 * THRESHOLD
  ),
  sprintf(
    "  HCC: k=%s",
    hcc_first_threshold
  ),
  sprintf(
    "  TRACERx100: k=%s",
    trx_first_threshold
  ),
  "  This difference is descriptive; the threshold is not a validation criterion.",
  "",
  "Reference-architecture context in fixed n>=5 cohorts:",
  sprintf(
    "  median non-ubiquitous fraction: HCC %.4f vs TRACERx100 %.4f",
    hcc_ref_nonfrac$hcc_median,
    hcc_ref_nonfrac$tracerx_median
  ),
  sprintf(
    "  median private fraction: HCC %.4f vs TRACERx100 %.4f",
    trx_ref_priv$hcc_median,
    trx_ref_priv$tracerx_median
  ),
  "  TRACERx100 has lower observed non-ubiquitous and private-event fractions",
  "  in this dataset; biology and WES/WGS calling differences cannot be disentangled.",
  "",
  "QC:",
  "  HCC Script 08 fixed-curve regression: PASS",
  "  TRACERx Script 11 available/strict N(k) regression: PASS",
  "  fixed patient sets constant across k=1..4: PASS",
  "  every fixed k=1..4 point is true downsampling: PASS",
  "  no zero-denominator main-endpoint patient in either fixed n>=5 cohort: PASS",
  "",
  "Interpretation guardrails:",
  "  - external validation is based on concordant learning-curve behavior, not equality;",
  "  - no cross-cohort p-value or equivalence claim is made;",
  "  - early-depth differences are interpreted in light of different reference ITH architecture;",
  "  - k=5 comparison is exploratory only because n>=6 sample sizes are HCC N=7 and TRACERx100 N=6;",
  "  - no spatial analysis is attempted in TRACERx100 because HCC-like spatial coordinates are unavailable.",
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed)
)

writeLines(
  summary_lines,
  file.path(
    PATHS$qc_dir,
    "12_hcc_tracerx100_qc_summary.txt"
  )
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
