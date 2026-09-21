# ==============================================================================
# 05_hcc_depth_analysis_exact_bootstrap.R
#
# Purpose
#   Cohort-level analysis of exhaustive UNCONSTRAINED HCC downsampling using the
#   validated all-filtered SNV/indel event set.
#
# Scientific estimand
#   Within patient:
#     mean across all k-sector subsets = exact expected performance under
#     uniformly random/unconstrained selection of k sectors.
#
#   Across patients:
#     median of patient-level expected performance, with IQR and patient-level
#     exact nonparametric percentile-bootstrap 95% CI.
#
# Prespecified analysis views
#   A. available_cohort:
#        all patients with n >= k; descriptive continuity with legacy analysis.
#   B. strict_downsampling:
#        only patients with n > k; removes k=n ceiling points at each k but the
#        patient population changes with k.
#   C. fixed_n_ge5:
#        the same 55 patients with n >= 5, analysed at k=1..4; primary robustness
#        analysis because every point is true downsampling (k<n).
#   D. n_ge6_exploratory:
#        the same 7 patients with n >= 6, analysed at k=1..5; exploratory check
#        of five-sector performance when information exists beyond sector 5.
#   E. original-n strata:
#        n=2,3,4,5 and n>=6, descriptive only.
#
# Additional diagnostics
#   - legacy subset-MEDIAN versus primary subset-MEAN sensitivity
#   - operational 80% adequate-recovery threshold
#   - paired marginal gains in fixed_n_ge5
#   - ceiling diagnostics (fraction with n=k)
#   - leave-one-patient-out median influence, including D010 if eligible
#
# This script does NOT analyse spatial design.
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
source(file.path(script_dir, "functions_bootstrap_exact.R"))

cat("Script 05: HCC depth analysis\n")
cat("============================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
depth_dir <- file.path(PATHS$output_dir, "results", "hcc_depth")
dir.create(depth_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_unconstrained_patient_summary_all_filtered.rds"
)

if (!file.exists(input_path)) {
  stop("Missing validated Script 04 input:\n", input_path)
}

x <- as.data.table(readRDS(input_path))

required_cols <- c(
  "patient_id", "n", "k", "sampling_fraction", "n_subsets",
  "recall_nonubiquitous_detection_mean",
  "recall_nonubiquitous_detection_median",
  "recall_heterogeneity_classification_mean",
  "recall_heterogeneity_classification_median",
  "apparent_ubiquity_error_mean",
  "apparent_ubiquity_error_median",
  "recall_private_mean",
  "recall_private_median",
  "nonubiquitous_fraction_abs_error_mean",
  "nonubiquitous_fraction_abs_error_median",
  "jaccard_ith_abs_error_mean",
  "jaccard_ith_abs_error_median"
)

if (!all(required_cols %in% names(x))) {
  stop(
    "Script 04 patient summary lacks required columns:\n",
    paste(setdiff(required_cols, names(x)), collapse = ", ")
  )
}

x[, patient_id := as.character(patient_id)]
x[, n := as.integer(n)]
x[, k := as.integer(k)]
setorder(x, patient_id, k)

# ------------------------------- constants -----------------------------------
B <- NA_integer_  # exact bootstrap: no Monte-Carlo replicates
BASE_SEED <- as.integer(PARAMS$seed)
THRESHOLD <- as.numeric(PARAMS$adequate_recovery_threshold)
CI_LEVEL <- 0.95
ALPHA <- 1 - CI_LEVEL

EXPECTED_AVAILABLE_N <- c(`1` = 122L, `2` = 122L, `3` = 100L, `4` = 76L, `5` = 55L)
EXPECTED_STRICT_N <- c(`1` = 122L, `2` = 100L, `3` = 76L, `4` = 55L, `5` = 7L)
EXPECTED_FIXED_N_GE5 <- 55L
EXPECTED_N_GE6 <- 7L

metric_map <- data.table(
  metric = c(
    "nonubiquitous_detection_recall",
    "heterogeneity_classification_recall",
    "apparent_ubiquity_error",
    "private_recall",
    "nonubiquitous_fraction_abs_error",
    "jaccard_ith_abs_error"
  ),
  subset_mean_column = c(
    "recall_nonubiquitous_detection_mean",
    "recall_heterogeneity_classification_mean",
    "apparent_ubiquity_error_mean",
    "recall_private_mean",
    "nonubiquitous_fraction_abs_error_mean",
    "jaccard_ith_abs_error_mean"
  ),
  subset_median_column = c(
    "recall_nonubiquitous_detection_median",
    "recall_heterogeneity_classification_median",
    "apparent_ubiquity_error_median",
    "recall_private_median",
    "nonubiquitous_fraction_abs_error_median",
    "jaccard_ith_abs_error_median"
  ),
  direction = c(
    "higher_better",
    "higher_better",
    "lower_better",
    "higher_better",
    "lower_better",
    "lower_better"
  )
)

# ------------------------------- helpers -------------------------------------
safe_median <- function(v) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  stats::median(v)
}

safe_quantile <- function(v, p) {
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  as.numeric(stats::quantile(v, probs = p, names = FALSE, type = 7))
}

bootstrap_median_ci <- function(v, B, seed, level = 0.95) {
  exact_bootstrap_median_ci(
    v = v,
    level = level,
    min_n = 1L
  )
}

summarise_one_k <- function(d, analysis_id, cohort_definition, k_value,
                            seed_base) {
  if (nrow(d) == 0L) {
    stop("No rows supplied to summarise_one_k for ", analysis_id,
         ", k=", k_value, ".")
  }

  if (anyDuplicated(d$patient_id) > 0L) {
    stop("Duplicate patient rows in ", analysis_id, ", k=", k_value, ".")
  }

  n_patients <- nrow(d)
  n_equal_k <- sum(d$n == k_value)

  rows <- vector("list", nrow(metric_map))

  for (mi in seq_len(nrow(metric_map))) {
    metric_id <- metric_map$metric[[mi]]
    col_name <- metric_map$subset_mean_column[[mi]]
    v <- d[[col_name]]
    n_metric <- sum(!is.na(v))

    ci <- bootstrap_median_ci(
      v,
      B = B,
      seed = seed_base + mi,
      level = CI_LEVEL
    )
    ci_defined <- all(is.finite(ci))

    rows[[mi]] <- data.table(
      analysis = analysis_id,
      cohort_definition = cohort_definition,
      k = as.integer(k_value),
      n_patients = as.integer(n_patients),
      n_metric = as.integer(n_metric),
      n_equal_k = as.integer(n_equal_k),
      fraction_equal_k = n_equal_k / n_patients,
      metric = metric_id,
      patient_level_subset_estimator = "mean",
      cohort_estimator = "median",
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
      bootstrap_level = if (ci_defined) CI_LEVEL else NA_real_,
      direction = metric_map$direction[[mi]]
    )
  }

  rbindlist(rows)
}

summarise_analysis <- function(analysis_id, cohort_definition, selector_by_k,
                               k_values, seed_offset) {
  ans <- vector("list", length(k_values))

  for (jj in seq_along(k_values)) {
    kval <- as.integer(k_values[[jj]])

    # selector_by_k returns a logical vector of length nrow(x).
    keep <- selector_by_k(kval)
    if (length(keep) != nrow(x) || anyNA(keep)) {
      stop("Invalid selector for ", analysis_id, ", k=", kval, ".")
    }

    d <- x[which(keep)]

    # Every patient must appear exactly once at this k.
    if (any(d$k != kval)) {
      stop("Selector returned a row with incorrect k in ", analysis_id, ".")
    }

    ans[[jj]] <- summarise_one_k(
      d = d,
      analysis_id = analysis_id,
      cohort_definition = cohort_definition,
      k_value = kval,
      seed_base = BASE_SEED + seed_offset + jj * 100L
    )
  }

  rbindlist(ans)
}

curve_detection_table <- function(summary_long) {
  summary_long[
    metric == "nonubiquitous_detection_recall",
    .(
      analysis,
      cohort_definition,
      k,
      n_patients,
      n_equal_k,
      fraction_equal_k,
      median_detection_recall = median,
      q1_detection_recall = q1,
      q3_detection_recall = q3,
      bootstrap_ci_lower,
      bootstrap_ci_upper
    )
  ]
}

# ------------------------------------------------------------------------------
# 0. Patient-level mathematical/QC properties before cohort aggregation
# ------------------------------------------------------------------------------

# For the exact expectation under exhaustive sampling, detection, correct
# heterogeneity classification, and private recall must be non-decreasing as k
# increases within every patient.
monotonic_cols <- c(
  "recall_nonubiquitous_detection_mean",
  "recall_heterogeneity_classification_mean",
  "recall_private_mean"
)

monotonic_failures <- list()
mf_i <- 0L

patient_ids <- unique(x$patient_id)
for (pid in patient_ids) {
  idx <- which(x$patient_id == pid)
  d <- x[idx][order(k)]

  for (cc in monotonic_cols) {
    vv <- d[[cc]]
    dd <- diff(vv)
    bad <- which(!is.na(dd) & dd < -1e-10)
    if (length(bad) > 0L) {
      mf_i <- mf_i + 1L
      monotonic_failures[[mf_i]] <- data.table(
        patient_id = pid,
        metric_column = cc,
        k_from = d$k[bad],
        k_to = d$k[bad + 1L],
        value_from = vv[bad],
        value_to = vv[bad + 1L]
      )
    }
  }
}

if (length(monotonic_failures) > 0L) {
  monotonic_failure_table <- rbindlist(monotonic_failures)
  if (DEBUG) {
    fwrite(
      monotonic_failure_table,
      file.path(PATHS$qc_dir, "05_patient_monotonicity_failures.tsv"),
      sep = "\t"
    )
  }
  stop("Patient-level monotonicity QC FAILED.")
}

# ------------------------------------------------------------------------------
# 1. Four prespecified depth-analysis views
# ------------------------------------------------------------------------------

available_summary <- summarise_analysis(
  analysis_id = "available_cohort",
  cohort_definition = "all patients with n >= k; descriptive/legacy-continuity view",
  selector_by_k = function(kval) {
    x$k == kval
  },
  k_values = 1:5,
  seed_offset = 10000L
)

strict_summary <- summarise_analysis(
  analysis_id = "strict_downsampling",
  cohort_definition = "patients with n > k; removes k=n ceiling points, changing cohort by k",
  selector_by_k = function(kval) {
    x$k == kval & x$n > kval
  },
  k_values = 1:5,
  seed_offset = 20000L
)

fixed5_summary <- summarise_analysis(
  analysis_id = "fixed_n_ge5",
  cohort_definition = "same patients with n >= 5, k=1..4; primary depth robustness",
  selector_by_k = function(kval) {
    x$k == kval & x$n >= 5L
  },
  k_values = 1:4,
  seed_offset = 30000L
)

nge6_summary <- summarise_analysis(
  analysis_id = "n_ge6_exploratory",
  cohort_definition = "same patients with n >= 6, k=1..5; exploratory five-sector check",
  selector_by_k = function(kval) {
    x$k == kval & x$n >= 6L
  },
  k_values = 1:5,
  seed_offset = 40000L
)

core_summary_long <- rbindlist(
  list(
    available_summary,
    strict_summary,
    fixed5_summary,
    nge6_summary
  ),
  use.names = TRUE
)

# Deterministic cohort-size QC.
observed_available_N <- curve_detection_table(available_summary)$n_patients
names(observed_available_N) <- as.character(1:5)

observed_strict_N <- curve_detection_table(strict_summary)$n_patients
names(observed_strict_N) <- as.character(1:5)

if (!identical(as.integer(observed_available_N),
               as.integer(EXPECTED_AVAILABLE_N))) {
  stop(
    "Available-cohort N(k) QC FAILED. Observed: ",
    paste(observed_available_N, collapse = ", ")
  )
}

if (!identical(as.integer(observed_strict_N),
               as.integer(EXPECTED_STRICT_N))) {
  stop(
    "Strict-downsampling N(k) QC FAILED. Observed: ",
    paste(observed_strict_N, collapse = ", ")
  )
}

fixed5_N <- unique(curve_detection_table(fixed5_summary)$n_patients)
if (length(fixed5_N) != 1L || fixed5_N != EXPECTED_FIXED_N_GE5) {
  stop("fixed_n_ge5 cohort-size QC FAILED.")
}

nge6_N <- unique(curve_detection_table(nge6_summary)$n_patients)
if (length(nge6_N) != 1L || nge6_N != EXPECTED_N_GE6) {
  stop("n_ge6 exploratory cohort-size QC FAILED.")
}

# Every fixed_n_ge5 and n_ge6 point must be genuine downsampling.
if (any(fixed5_summary$fraction_equal_k != 0)) {
  stop("fixed_n_ge5 unexpectedly contains k=n points.")
}
if (any(nge6_summary$fraction_equal_k != 0)) {
  stop("n_ge6_exploratory unexpectedly contains k=n points.")
}
if (any(strict_summary$fraction_equal_k != 0)) {
  stop("strict_downsampling unexpectedly contains k=n points.")
}

# ------------------------------------------------------------------------------
# 2. Ceiling diagnostic for available-cohort curve
# ------------------------------------------------------------------------------

ceiling_diagnostic <- unique(
  available_summary[
    metric == "nonubiquitous_detection_recall",
    .(
      k,
      n_patients,
      n_equal_k,
      fraction_equal_k
    )
  ]
)
setorder(ceiling_diagnostic, k)

# ------------------------------------------------------------------------------
# 3. Operational 80% adequate-recovery threshold
# ------------------------------------------------------------------------------

core_detection <- core_summary_long[
  metric == "nonubiquitous_detection_recall"
]

threshold_by_k <- core_detection[
  ,
  .(
    analysis,
    cohort_definition,
    k,
    n_patients,
    n_equal_k,
    fraction_equal_k,
    median_detection_recall = median,
    threshold = THRESHOLD,
    median_reaches_threshold = median >= THRESHOLD
  )
]

threshold_first <- threshold_by_k[
  ,
  {
    kk <- k[median_reaches_threshold %in% TRUE]
    list(
      threshold = THRESHOLD,
      first_k_with_median_ge_threshold =
        if (length(kk) == 0L) NA_integer_ else min(kk)
    )
  },
  by = .(analysis, cohort_definition)
]

# Add patient-level proportion reaching threshold at each k.
threshold_prop_rows <- list()
tp_i <- 0L

for (analysis_id in unique(threshold_by_k$analysis)) {
  ks <- threshold_by_k[analysis == analysis_id, k]

  for (kval in ks) {
    if (analysis_id == "available_cohort") {
      keep <- x$k == kval
    } else if (analysis_id == "strict_downsampling") {
      keep <- x$k == kval & x$n > kval
    } else if (analysis_id == "fixed_n_ge5") {
      keep <- x$k == kval & x$n >= 5L
    } else if (analysis_id == "n_ge6_exploratory") {
      keep <- x$k == kval & x$n >= 6L
    } else {
      stop("Unknown threshold analysis: ", analysis_id)
    }

    d <- x[which(keep)]
    vv <- d$recall_nonubiquitous_detection_mean

    tp_i <- tp_i + 1L
    threshold_prop_rows[[tp_i]] <- data.table(
      analysis = analysis_id,
      k = as.integer(kval),
      n_patients = nrow(d),
      n_patients_ge_threshold = sum(vv >= THRESHOLD, na.rm = TRUE),
      proportion_patients_ge_threshold =
        mean(vv >= THRESHOLD, na.rm = TRUE)
    )
  }
}

threshold_proportions <- rbindlist(threshold_prop_rows)

threshold_by_k <- merge(
  threshold_by_k,
  threshold_proportions,
  by = c("analysis", "k", "n_patients"),
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------------------------
# 4. Paired marginal information gains in the fixed n>=5 cohort
# ------------------------------------------------------------------------------

fixed_rows <- x[which(x$n >= 5L & x$k <= 4L)]
fixed_patient_ids <- sort(unique(fixed_rows$patient_id))

if (length(fixed_patient_ids) != EXPECTED_FIXED_N_GE5) {
  stop("Unexpected number of fixed_n_ge5 patient IDs.")
}

gain_metric_map <- data.table(
  metric = c(
    "nonubiquitous_detection_recall",
    "heterogeneity_classification_recall"
  ),
  column = c(
    "recall_nonubiquitous_detection_mean",
    "recall_heterogeneity_classification_mean"
  )
)

gain_rows <- list()
gi <- 0L

for (mi in seq_len(nrow(gain_metric_map))) {
  metric_id <- gain_metric_map$metric[[mi]]
  col_name <- gain_metric_map$column[[mi]]

  d <- fixed_rows[, .(patient_id, k, value = get(col_name))]
  wide <- dcast(d, patient_id ~ k, value.var = "value")

  if (!all(as.character(1:4) %in% names(wide))) {
    stop("Fixed-cohort gain table lacks one or more k columns.")
  }
  if (nrow(wide) != EXPECTED_FIXED_N_GE5) {
    stop("Fixed-cohort gain table has unexpected patient count.")
  }

  for (k_to in 2:4) {
    k_from <- k_to - 1L
    delta <- wide[[as.character(k_to)]] - wide[[as.character(k_from)]]

    ci <- bootstrap_median_ci(
      delta,
      B = B,
      seed = BASE_SEED + 50000L + mi * 1000L + k_to,
      level = CI_LEVEL
    )

    gi <- gi + 1L
    gain_rows[[gi]] <- data.table(
      analysis = "fixed_n_ge5_paired_gain",
      metric = metric_id,
      k_from = k_from,
      k_to = k_to,
      n_patients = length(delta),
      median_patient_gain = safe_median(delta),
      q1_patient_gain = safe_quantile(delta, 0.25),
      q3_patient_gain = safe_quantile(delta, 0.75),
      bootstrap_ci_lower = ci[[1L]],
      bootstrap_ci_upper = ci[[2L]],
      bootstrap_replicates = NA_integer_,
      bootstrap_method = "exact_nonparametric_percentile",
      bootstrap_level = CI_LEVEL
    )
  }
}

fixed5_paired_gains <- rbindlist(gain_rows)

# ------------------------------------------------------------------------------
# 5. Original sampling-depth strata
# ------------------------------------------------------------------------------

strata <- list(
  n_exact_2 = function() x$n == 2L,
  n_exact_3 = function() x$n == 3L,
  n_exact_4 = function() x$n == 4L,
  n_exact_5 = function() x$n == 5L,
  n_ge_6 = function() x$n >= 6L
)

strata_k <- list(
  n_exact_2 = 1:2,
  n_exact_3 = 1:3,
  n_exact_4 = 1:4,
  n_exact_5 = 1:5,
  n_ge_6 = 1:5
)

stratified_rows <- list()
si <- 0L

for (stratum_name in names(strata)) {
  stratum_keep_all <- strata[[stratum_name]]()

  for (kval in strata_k[[stratum_name]]) {
    keep <- stratum_keep_all & x$k == kval
    d <- x[which(keep)]

    if (nrow(d) == 0L) next

    si <- si + 1L
    z <- summarise_one_k(
      d = d,
      analysis_id = paste0("stratum_", stratum_name),
      cohort_definition = paste0("original sampling-depth stratum: ", stratum_name),
      k_value = kval,
      seed_base = BASE_SEED + 60000L + si * 100L
    )
    z[, stratum := stratum_name]
    z[, is_full_reference_point := all(d$n == kval)]
    stratified_rows[[si]] <- z
  }
}

stratified_summary <- rbindlist(stratified_rows, use.names = TRUE, fill = TRUE)

# ------------------------------------------------------------------------------
# 6. Primary subset-MEAN vs legacy subset-MEDIAN sensitivity
# ------------------------------------------------------------------------------

sensitivity_analyses <- c("available_cohort", "fixed_n_ge5")
sensitivity_rows <- list()
ssi <- 0L

for (analysis_id in sensitivity_analyses) {
  if (analysis_id == "available_cohort") {
    k_values <- 1:5
  } else {
    k_values <- 1:4
  }

  for (kval in k_values) {
    if (analysis_id == "available_cohort") {
      keep <- x$k == kval
    } else {
      keep <- x$k == kval & x$n >= 5L
    }

    d <- x[which(keep)]

    for (mi in seq_len(nrow(metric_map))) {
      metric_id <- metric_map$metric[[mi]]
      mean_col <- metric_map$subset_mean_column[[mi]]
      median_col <- metric_map$subset_median_column[[mi]]

      primary_cohort_median <- safe_median(d[[mean_col]])
      legacy_cohort_median <- safe_median(d[[median_col]])

      ssi <- ssi + 1L
      sensitivity_rows[[ssi]] <- data.table(
        analysis = analysis_id,
        k = kval,
        n_patients = nrow(d),
        metric = metric_id,
        cohort_median_of_subset_mean = primary_cohort_median,
        cohort_median_of_subset_median = legacy_cohort_median,
        legacy_minus_primary =
          legacy_cohort_median - primary_cohort_median
      )
    }
  }
}

mean_vs_median_sensitivity <- rbindlist(sensitivity_rows)

# ------------------------------------------------------------------------------
# 7. Leave-one-patient-out median influence
# ------------------------------------------------------------------------------

influence_analyses <- c("available_cohort", "fixed_n_ge5")
influence_metrics <- gain_metric_map
influence_rows <- list()
ii <- 0L

for (analysis_id in influence_analyses) {
  k_values <- if (analysis_id == "available_cohort") 1:5 else 1:4

  for (kval in k_values) {
    if (analysis_id == "available_cohort") {
      keep <- x$k == kval
    } else {
      keep <- x$k == kval & x$n >= 5L
    }

    d <- x[which(keep)]

    for (mi in seq_len(nrow(influence_metrics))) {
      metric_id <- influence_metrics$metric[[mi]]
      col_name <- influence_metrics$column[[mi]]
      vv <- d[[col_name]]

      baseline <- safe_median(vv)
      loo_values <- numeric(nrow(d))

      for (jj in seq_len(nrow(d))) {
        loo_values[[jj]] <- safe_median(vv[-jj])
      }

      abs_change <- abs(loo_values - baseline)
      max_idx <- which.max(abs_change)

      d010_idx <- which(d$patient_id == "D010")
      d010_change <- if (length(d010_idx) == 1L) {
        loo_values[[d010_idx]] - baseline
      } else {
        NA_real_
      }

      ii <- ii + 1L
      influence_rows[[ii]] <- data.table(
        analysis = analysis_id,
        k = kval,
        metric = metric_id,
        n_patients = nrow(d),
        baseline_median = baseline,
        max_abs_leave_one_out_change = abs_change[[max_idx]],
        patient_causing_max_abs_change = d$patient_id[[max_idx]],
        d010_is_in_cohort = length(d010_idx) == 1L,
        d010_leave_one_out_change = d010_change
      )
    }
  }
}

leave_one_out_influence <- rbindlist(influence_rows)

# ------------------------------------------------------------------------------
# 8. Save canonical outputs
# ------------------------------------------------------------------------------

fwrite(
  core_summary_long,
  file.path(depth_dir, "05_hcc_depth_core_summary_all_filtered.tsv"),
  sep = "\t"
)

if (DEBUG) {
  fwrite(
    ceiling_diagnostic,
    file.path(depth_dir, "05_hcc_available_cohort_ceiling_diagnostic.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    threshold_by_k,
    file.path(depth_dir, "05_hcc_adequate_recovery_threshold_by_k.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    threshold_first,
    file.path(depth_dir, "05_hcc_adequate_recovery_threshold_first_k.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    fixed5_paired_gains,
    file.path(depth_dir, "05_hcc_fixed_n_ge5_paired_gains.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    stratified_summary,
    file.path(depth_dir, "05_hcc_depth_stratified_by_original_n.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    mean_vs_median_sensitivity,
    file.path(depth_dir, "05_hcc_subset_mean_vs_median_sensitivity.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    leave_one_out_influence,
    file.path(PATHS$qc_dir, "05_leave_one_out_median_influence.tsv"),
    sep = "\t"
  )
}


# ------------------------------------------------------------------------------
# 9. Compact human-readable QC/result summary
# ------------------------------------------------------------------------------

get_detection_row <- function(summary_table, analysis_id, kval) {
  idx <- which(
    summary_table$analysis == analysis_id &
      summary_table$metric == "nonubiquitous_detection_recall" &
      summary_table$k == kval
  )
  if (length(idx) != 1L) {
    stop("Could not identify unique detection row for ", analysis_id,
         ", k=", kval, ".")
  }
  summary_table[idx]
}

get_classification_row <- function(summary_table, analysis_id, kval) {
  idx <- which(
    summary_table$analysis == analysis_id &
      summary_table$metric == "heterogeneity_classification_recall" &
      summary_table$k == kval
  )
  if (length(idx) != 1L) {
    stop("Could not identify unique classification row for ", analysis_id,
         ", k=", kval, ".")
  }
  summary_table[idx]
}

avail4 <- get_detection_row(core_summary_long, "available_cohort", 4L)
avail5 <- get_detection_row(core_summary_long, "available_cohort", 5L)
fixed4 <- get_detection_row(core_summary_long, "fixed_n_ge5", 4L)
fixed4_class <- get_classification_row(core_summary_long, "fixed_n_ge5", 4L)
nge6_5 <- get_detection_row(core_summary_long, "n_ge6_exploratory", 5L)
nge6_5_class <- get_classification_row(
  core_summary_long, "n_ge6_exploratory", 5L
)

threshold_fixed <- threshold_first[
  analysis == "fixed_n_ge5",
  first_k_with_median_ge_threshold
]
threshold_available <- threshold_first[
  analysis == "available_cohort",
  first_k_with_median_ge_threshold
]

max_d010_detection_effect <- max(
  abs(
    leave_one_out_influence[
      metric == "nonubiquitous_detection_recall",
      d010_leave_one_out_change
    ]
  ),
  na.rm = TRUE
)

max_loo_detection_effect <- max(
  leave_one_out_influence[
    metric == "nonubiquitous_detection_recall",
    max_abs_leave_one_out_change
  ],
  na.rm = TRUE
)

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  "Script 05: HCC depth analysis",
  "============================================",
  "",
  "Bootstrap method: exact empirical nonparametric distribution of the patient-level median",
  sprintf("Bootstrap CI: %.0f%% percentile; no Monte-Carlo simulation error",
          100 * CI_LEVEL),
  sprintf("Operational adequate-recovery threshold: %.2f", THRESHOLD),
  "",
  "Cohort-size QC: PASS",
  sprintf(
    "  available N(k=1..5): %s",
    paste(EXPECTED_AVAILABLE_N, collapse = ", ")
  ),
  sprintf(
    "  strict n>k N(k=1..5): %s",
    paste(EXPECTED_STRICT_N, collapse = ", ")
  ),
  sprintf("  fixed n>=5: N=%d at k=1..4", EXPECTED_FIXED_N_GE5),
  sprintf("  n>=6 exploratory: N=%d at k=1..5", EXPECTED_N_GE6),
  "Patient-level monotonicity QC: PASS",
  "",
  "Key all-filtered SNV/indel depth results:",
  sprintf(
    "  available cohort k=4: median expected non-ubiquitous detection recall = %.4f (N=%d)",
    avail4$median, avail4$n_patients
  ),
  sprintf(
    "  available cohort k=5: median expected non-ubiquitous detection recall = %.4f (N=%d; %.1f%% have n=k)",
    avail5$median,
    avail5$n_patients,
    100 * avail5$fraction_equal_k
  ),
  sprintf(
    "  fixed n>=5 k=4: median expected non-ubiquitous detection recall = %.4f (N=%d)",
    fixed4$median, fixed4$n_patients
  ),
  sprintf(
    "  fixed n>=5 k=4: median expected heterogeneity-classification recall = %.4f",
    fixed4_class$median
  ),
  sprintf(
    "  n>=6 exploratory k=5: median expected non-ubiquitous detection recall = %.4f (N=%d)",
    nge6_5$median, nge6_5$n_patients
  ),
  sprintf(
    "  n>=6 exploratory k=5: median expected heterogeneity-classification recall = %.4f",
    nge6_5_class$median
  ),
  "",
  sprintf(
    "First k with cohort median detection recall >= %.2f:",
    THRESHOLD
  ),
  sprintf("  available cohort: k=%s", threshold_available),
  sprintf("  fixed n>=5 (k<=4): k=%s", threshold_fixed),
  "",
  "Influence QC (patient-weighted medians):",
  sprintf(
    "  maximum absolute D010 deletion effect on detection-recall median = %.6f",
    max_d010_detection_effect
  ),
  sprintf(
    "  maximum leave-one-patient-out effect on detection-recall median across core analyses = %.6f",
    max_loo_detection_effect
  ),
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed),
  "",
  "Interpretation guardrails:",
  "  - available k=5 contains substantial k=n ceiling contribution and is descriptive;",
  "  - fixed n>=5, k=1..4 is the primary robustness view for depth;",
  "  - n>=6, k=5 is exploratory because N=7;",
  "  - the 0.80 threshold is an operational legacy threshold, not a biological optimum;",
  "  - no spatial sampling comparison is performed in Script 05."
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "05_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
