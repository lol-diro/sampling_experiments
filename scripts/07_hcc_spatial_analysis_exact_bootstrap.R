# ==============================================================================
# 07_hcc_spatial_analysis_exact_bootstrap.R
#
# Purpose
#   Cohort-level paired analysis of HCC spatial sampling designs constructed and
#   validated in Script 06.
#
# Primary scientific contrast
#   grid-based dispersed sampling versus clustered/adjacent sampling
#
# Sign convention
#   ALL "advantage" quantities are coded so that:
#       delta > 0  => the left-hand / more dispersed design performs better
#
# Primary spatial inference
#   strict spatial cohort, all available patients with n > k, k = 2,3,4
#
# Prespecified robustness / sensitivity
#   - strict spatial fixed cohort n >= 5, same 44 patients at k = 2,3,4
#   - order-only spatial sensitivity, k = 2,3,4
#   - order-only fixed n >= 5 sensitivity, same 45 patients at k = 2,3,4
#
# Exploratory only
#   k = 5 (N=4 in the strict spatial cohort)
#
# Primary endpoint
#   non-ubiquitous detection recall
#
# Key secondary endpoint
#   heterogeneity-classification recall
#
# Supportive endpoints
#   private recall
#   apparent-ubiquity error
#   non-ubiquitous-fraction absolute error
#   Jaccard ITH absolute error
#
# Statistical unit
#   patient
#
# Main effect estimator
#   median of paired patient-level grid-minus-cluster advantages
#   with exact patient-level nonparametric percentile-bootstrap 95% CI.
#
# No subset is treated as an independent biological replicate.
# No p-value is required for the primary analysis; effect size, uncertainty,
# direction consistency, and robustness are reported explicitly.
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

cat("Script 07: spatial analysis\n")
cat("==========================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
spatial_dir <- file.path(PATHS$output_dir, "results", "hcc_spatial")
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_spatial_designs_all_filtered.rds"
)

if (!file.exists(input_path)) {
  stop("Missing validated Script 06 input:\n", input_path)
}

obj <- readRDS(input_path)

required_top <- c("primary", "order_only_sensitivity")
if (!all(required_top %in% names(obj))) {
  stop("Script 06 RDS does not have the expected top-level structure.")
}

primary_design <- as.data.table(obj$primary$design_summary)
primary_delta <- as.data.table(obj$primary$patient_deltas)

order_design <- as.data.table(obj$order_only_sensitivity$design_summary)
order_delta <- as.data.table(obj$order_only_sensitivity$patient_deltas)

# ------------------------------- constants -----------------------------------
B <- NA_integer_  # exact bootstrap: no Monte-Carlo replicates
BASE_SEED <- as.integer(PARAMS$seed)
CI_LEVEL <- 0.95
MIN_BOOTSTRAP_N <- 20L
TOLERANCE <- 1e-12

PRIMARY_K <- 2:4
EXPLORATORY_K <- 5L

EXPECTED_N <- list(
  primary_strict_all = c(`2` = 80L, `3` = 60L, `4` = 44L),
  primary_strict_fixed_n_ge5 = c(`2` = 44L, `3` = 44L, `4` = 44L),
  order_only_all = c(`2` = 86L, `3` = 65L, `4` = 45L),
  order_only_fixed_n_ge5 = c(`2` = 45L, `3` = 45L, `4` = 45L),
  primary_k5 = c(`5` = 4L),
  order_k5 = c(`5` = 4L)
)

metric_table <- data.table(
  metric = c(
    "recall_nonubiquitous_detection",
    "recall_heterogeneity_classification",
    "recall_private",
    "apparent_ubiquity_error",
    "nonubiquitous_fraction_abs_error",
    "jaccard_ith_abs_error"
  ),
  endpoint_role = c(
    "primary",
    "key_secondary",
    "supportive",
    "supportive",
    "supportive",
    "supportive"
  ),
  raw_direction = c(
    "higher_better",
    "higher_better",
    "higher_better",
    "lower_better",
    "lower_better",
    "lower_better"
  )
)

contrast_table <- data.table(
  contrast = c(
    "grid_vs_cluster",
    "grid_vs_unconstrained",
    "unconstrained_vs_cluster"
  ),
  column_prefix = c(
    "grid_advantage__",
    "grid_minus_unconstrained__",
    "unconstrained_minus_cluster__"
  ),
  role = c(
    "primary_spatial_contrast",
    "supportive_ordering",
    "supportive_ordering"
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
  as.numeric(stats::quantile(v, probs = p, names = FALSE, type = 7))
}

safe_min <- function(v) {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  min(v)
}

safe_max <- function(v) {
  v <- as.numeric(v)
  v <- v[!is.na(v)]
  if (length(v) == 0L) return(NA_real_)
  max(v)
}

bootstrap_median_ci <- function(v, B, seed, level = 0.95,
                                min_n = MIN_BOOTSTRAP_N) {
  exact_bootstrap_median_ci(
    v = v,
    level = level,
    min_n = min_n
  )
}

validate_N_by_k <- function(d, expected, label) {
  observed <- d[
    ,
    .(n_patients = uniqueN(patient_id)),
    by = k
  ]
  observed <- observed[k %in% as.integer(names(expected))]
  setorder(observed, k)

  exp_dt <- data.table(
    k = as.integer(names(expected)),
    expected_n_patients = as.integer(expected)
  )

  z <- merge(exp_dt, observed, by = "k", all = TRUE, sort = TRUE)
  z[, analysis := label]
  z[, pass := expected_n_patients == n_patients]

  if (!all(z$pass)) {
    if (DEBUG) {
      fwrite(
        z,
        file.path(
          PATHS$qc_dir,
          paste0("07_", label, "_N_failure.tsv")
        ),
        sep = "\t"
      )
    }
    stop(label, ": cohort-size QC FAILED.")
  }

  z
}

summarise_paired_contrasts <- function(d, analysis_id,
                                       cohort_definition,
                                       inferential_status,
                                       seed_offset) {
  if (nrow(d) == 0L) {
    stop("No data in analysis ", analysis_id, ".")
  }

  if (anyDuplicated(d[, .(patient_id, k)]) > 0L) {
    stop("Duplicate patient x k rows in ", analysis_id, ".")
  }

  k_values <- sort(unique(d$k))
  rows <- list()
  ri <- 0L

  for (kval in k_values) {
    dk <- d[k == kval]

    for (ci in seq_len(nrow(contrast_table))) {
      contrast_id <- contrast_table$contrast[[ci]]
      prefix <- contrast_table$column_prefix[[ci]]
      contrast_role <- contrast_table$role[[ci]]

      for (mi in seq_len(nrow(metric_table))) {
        metric_id <- metric_table$metric[[mi]]
        endpoint_role <- metric_table$endpoint_role[[mi]]
        col_name <- paste0(prefix, metric_id)

        if (!col_name %in% names(dk)) {
          stop(
            "Required paired-delta column missing: ", col_name,
            " in ", analysis_id, "."
          )
        }

        v <- as.numeric(dk[[col_name]])
        keep <- !is.na(v)
        vv <- v[keep]
        ids <- dk$patient_id[keep]

        if (length(vv) == 0L) {
          next
        }

        ci_boot <- bootstrap_median_ci(
          vv,
          B = B,
          seed = BASE_SEED + seed_offset +
            kval * 1000L + ci * 100L + mi,
          level = CI_LEVEL
        )

        n_pos <- sum(vv > TOLERANCE)
        n_neg <- sum(vv < -TOLERANCE)
        n_zero <- length(vv) - n_pos - n_neg

        ri <- ri + 1L
        rows[[ri]] <- data.table(
          analysis = analysis_id,
          cohort_definition = cohort_definition,
          inferential_status = inferential_status,
          contrast = contrast_id,
          contrast_role = contrast_role,
          metric = metric_id,
          endpoint_role = endpoint_role,
          k = as.integer(kval),

          n_patients = nrow(dk),
          n_metric = length(vv),

          median_paired_advantage = safe_median(vv),
          q1_paired_advantage = safe_quantile(vv, 0.25),
          q3_paired_advantage = safe_quantile(vv, 0.75),
          mean_paired_advantage = safe_mean(vv),
          min_paired_advantage = safe_min(vv),
          max_paired_advantage = safe_max(vv),

          n_positive = n_pos,
          n_negative = n_neg,
          n_zero = n_zero,
          proportion_positive = n_pos / length(vv),
          proportion_negative = n_neg / length(vv),
          proportion_zero = n_zero / length(vv),

          bootstrap_ci_lower = ci_boot[[1L]],
          bootstrap_ci_upper = ci_boot[[2L]],
          bootstrap_replicates = NA_integer_,
          bootstrap_method =
            if (length(vv) >= MIN_BOOTSTRAP_N) "exact_nonparametric_percentile" else NA_character_,
          bootstrap_level =
            if (length(vv) >= MIN_BOOTSTRAP_N) CI_LEVEL else NA_real_,

          sign_convention =
            "positive = left/more-dispersed design performs better"
        )
      }
    }
  }

  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

summarise_absolute_design_performance <- function(d, analysis_id,
                                                  cohort_definition,
                                                  inferential_status) {
  if (nrow(d) == 0L) {
    stop("No design-summary data in ", analysis_id, ".")
  }

  required_designs <- c(
    "grid_dispersed",
    "unconstrained",
    "clustered_adjacent"
  )

  if (!all(required_designs %in% unique(d$design))) {
    stop("Not all three designs are present in ", analysis_id, ".")
  }

  rows <- list()
  ri <- 0L

  for (kval in sort(unique(d$k))) {
    dk <- d[k == kval]

    for (design_id in required_designs) {
      dd <- dk[design == design_id]

      for (mi in seq_len(nrow(metric_table))) {
        metric_id <- metric_table$metric[[mi]]
        endpoint_role <- metric_table$endpoint_role[[mi]]
        col_name <- paste0(metric_id, "_mean")

        if (!col_name %in% names(dd)) {
          stop(
            "Required design metric missing: ", col_name,
            " in ", analysis_id, "."
          )
        }

        v <- dd[[col_name]]

        ri <- ri + 1L
        rows[[ri]] <- data.table(
          analysis = analysis_id,
          cohort_definition = cohort_definition,
          inferential_status = inferential_status,
          k = as.integer(kval),
          design = design_id,
          metric = metric_id,
          endpoint_role = endpoint_role,
          n_patients = nrow(dd),
          n_metric = sum(!is.na(v)),
          cohort_median = safe_median(v),
          cohort_q1 = safe_quantile(v, 0.25),
          cohort_q3 = safe_quantile(v, 0.75),
          cohort_mean = safe_mean(v)
        )
      }
    }
  }

  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

leave_one_out_median_influence <- function(d, analysis_id) {
  main_metrics <- c(
    "recall_nonubiquitous_detection",
    "recall_heterogeneity_classification"
  )

  rows <- list()
  ri <- 0L

  for (kval in sort(unique(d$k))) {
    dk <- d[k == kval]

    for (metric_id in main_metrics) {
      col_name <- paste0("grid_advantage__", metric_id)
      v <- as.numeric(dk[[col_name]])
      keep <- !is.na(v)
      vv <- v[keep]
      ids <- dk$patient_id[keep]

      if (length(vv) < 2L) next

      baseline <- stats::median(vv)
      loo <- vapply(
        seq_along(vv),
        function(jj) stats::median(vv[-jj]),
        numeric(1)
      )
      changes <- loo - baseline
      abs_changes <- abs(changes)
      max_abs_change <- max(abs_changes)
      tied <- which(
        abs(abs_changes - max_abs_change) <=
          sqrt(.Machine$double.eps)
      )

      ri <- ri + 1L
      rows[[ri]] <- data.table(
        analysis = analysis_id,
        k = as.integer(kval),
        metric = metric_id,
        n_patients = length(vv),
        baseline_median_paired_advantage = baseline,
        max_abs_leave_one_out_change = max_abs_change,
        n_patients_tied_for_max_abs_change = length(tied),
        patients_tied_for_max_abs_change =
          paste(ids[tied], collapse = ";"),
        signed_leave_one_out_changes_at_max =
          paste(
            format(
              changes[tied],
              digits = 17,
              scientific = FALSE,
              trim = TRUE
            ),
            collapse = ";"
          )
      )
    }
  }

  rbindlist(rows)
}

rank_extreme_patient_effects <- function(d, analysis_id, n_each_side = 5L) {
  main_metrics <- c(
    "recall_nonubiquitous_detection",
    "recall_heterogeneity_classification"
  )

  rows <- list()
  ri <- 0L

  for (kval in sort(unique(d$k))) {
    dk <- d[k == kval]

    for (metric_id in main_metrics) {
      col_name <- paste0("grid_advantage__", metric_id)
      z <- dk[
        !is.na(get(col_name)),
        .(
          patient_id,
          n,
          k,
          paired_advantage = get(col_name)
        )
      ]
      setorder(z, paired_advantage, patient_id)

      n_side <- min(as.integer(n_each_side), nrow(z))

      low <- head(z, n_side)
      low[, extreme_direction := "most_negative"]
      low[, rank_within_direction := seq_len(.N)]

      high <- tail(z, n_side)
      setorder(high, -paired_advantage, patient_id)
      high[, extreme_direction := "most_positive"]
      high[, rank_within_direction := seq_len(.N)]

      ri <- ri + 1L
      zz <- rbindlist(list(low, high))
      zz[, analysis := analysis_id]
      zz[, metric := metric_id]
      rows[[ri]] <- zz
    }
  }

  rbindlist(rows, use.names = TRUE)
}

# ------------------------------------------------------------------------------
# 1. Build prespecified analysis datasets
# ------------------------------------------------------------------------------

# Primary strict: Script 06 already contains only k<n rows.
primary_all_delta <- primary_delta[k %in% PRIMARY_K]
primary_all_design <- primary_design[k %in% PRIMARY_K]

primary_fixed_delta <- primary_delta[k %in% PRIMARY_K & n >= 5L]
primary_fixed_design <- primary_design[k %in% PRIMARY_K & n >= 5L]

order_all_delta <- order_delta[k %in% PRIMARY_K]
order_all_design <- order_design[k %in% PRIMARY_K]

order_fixed_delta <- order_delta[k %in% PRIMARY_K & n >= 5L]
order_fixed_design <- order_design[k %in% PRIMARY_K & n >= 5L]

primary_k5_delta <- primary_delta[k == EXPLORATORY_K]
primary_k5_design <- primary_design[k == EXPLORATORY_K]

order_k5_delta <- order_delta[k == EXPLORATORY_K]
order_k5_design <- order_design[k == EXPLORATORY_K]

# ------------------------------------------------------------------------------
# 2. Deterministic cohort-size and fixed-cohort QC
# ------------------------------------------------------------------------------

N_qc <- rbindlist(
  list(
    validate_N_by_k(
      primary_all_delta,
      EXPECTED_N$primary_strict_all,
      "primary_strict_all"
    ),
    validate_N_by_k(
      primary_fixed_delta,
      EXPECTED_N$primary_strict_fixed_n_ge5,
      "primary_strict_fixed_n_ge5"
    ),
    validate_N_by_k(
      order_all_delta,
      EXPECTED_N$order_only_all,
      "order_only_all"
    ),
    validate_N_by_k(
      order_fixed_delta,
      EXPECTED_N$order_only_fixed_n_ge5,
      "order_only_fixed_n_ge5"
    ),
    validate_N_by_k(
      primary_k5_delta,
      EXPECTED_N$primary_k5,
      "primary_k5"
    ),
    validate_N_by_k(
      order_k5_delta,
      EXPECTED_N$order_k5,
      "order_k5"
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

# Fixed cohort must contain exactly the same patient set at k=2,3,4.
validate_fixed_patient_set <- function(d, expected_n, label) {
  sets <- lapply(
    PRIMARY_K,
    function(kval) sort(d[k == kval, unique(patient_id)])
  )
  identical_sets <- all(vapply(
    sets[-1L],
    function(z) identical(z, sets[[1L]]),
    logical(1)
  ))

  if (!identical_sets || length(sets[[1L]]) != expected_n) {
    stop(label, ": fixed patient-set QC FAILED.")
  }
  TRUE
}

fixed_primary_pass <- validate_fixed_patient_set(
  primary_fixed_delta, 44L, "primary_strict_fixed_n_ge5"
)
fixed_order_pass <- validate_fixed_patient_set(
  order_fixed_delta, 45L, "order_only_fixed_n_ge5"
)

# ------------------------------------------------------------------------------
# 3. Paired spatial contrasts
# ------------------------------------------------------------------------------

contrast_summaries <- rbindlist(
  list(
    summarise_paired_contrasts(
      primary_all_delta,
      analysis_id = "primary_strict_all_n_gt_k",
      cohort_definition =
        "strict spatial cohort; all available patients with n>k; k=2..4",
      inferential_status = "primary",
      seed_offset = 10000L
    ),
    summarise_paired_contrasts(
      primary_fixed_delta,
      analysis_id = "primary_strict_fixed_n_ge5",
      cohort_definition =
        "strict spatial cohort; same patients with n>=5 at k=2..4",
      inferential_status = "robustness",
      seed_offset = 20000L
    ),
    summarise_paired_contrasts(
      order_all_delta,
      analysis_id = "order_only_all_n_gt_k",
      cohort_definition =
        "order-only sensitivity cohort; all available patients with n>k; k=2..4",
      inferential_status = "sensitivity",
      seed_offset = 30000L
    ),
    summarise_paired_contrasts(
      order_fixed_delta,
      analysis_id = "order_only_fixed_n_ge5",
      cohort_definition =
        "order-only sensitivity cohort; same patients with n>=5 at k=2..4",
      inferential_status = "sensitivity",
      seed_offset = 40000L
    ),
    summarise_paired_contrasts(
      primary_k5_delta,
      analysis_id = "primary_strict_k5",
      cohort_definition =
        "strict spatial cohort; k=5 and n>5",
      inferential_status = "exploratory_N4_no_CI",
      seed_offset = 50000L
    ),
    summarise_paired_contrasts(
      order_k5_delta,
      analysis_id = "order_only_k5",
      cohort_definition =
        "order-only sensitivity cohort; k=5 and n>5",
      inferential_status = "exploratory_N4_no_CI",
      seed_offset = 60000L
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------------------------
# 4. Absolute design performance
# ------------------------------------------------------------------------------

absolute_summaries <- rbindlist(
  list(
    summarise_absolute_design_performance(
      primary_all_design,
      "primary_strict_all_n_gt_k",
      "strict spatial cohort; all available patients with n>k; k=2..4",
      "primary"
    ),
    summarise_absolute_design_performance(
      primary_fixed_design,
      "primary_strict_fixed_n_ge5",
      "strict spatial cohort; same patients with n>=5 at k=2..4",
      "robustness"
    ),
    summarise_absolute_design_performance(
      order_all_design,
      "order_only_all_n_gt_k",
      "order-only sensitivity cohort; all available patients with n>k; k=2..4",
      "sensitivity"
    ),
    summarise_absolute_design_performance(
      order_fixed_design,
      "order_only_fixed_n_ge5",
      "order-only sensitivity cohort; same patients with n>=5 at k=2..4",
      "sensitivity"
    ),
    summarise_absolute_design_performance(
      primary_k5_design,
      "primary_strict_k5",
      "strict spatial cohort; k=5 and n>5",
      "exploratory_N4"
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

# ------------------------------------------------------------------------------
# 5. Influence and patient-level heterogeneity diagnostics
# ------------------------------------------------------------------------------

influence <- rbindlist(
  list(
    leave_one_out_median_influence(
      primary_all_delta,
      "primary_strict_all_n_gt_k"
    ),
    leave_one_out_median_influence(
      primary_fixed_delta,
      "primary_strict_fixed_n_ge5"
    ),
    leave_one_out_median_influence(
      order_all_delta,
      "order_only_all_n_gt_k"
    ),
    leave_one_out_median_influence(
      order_fixed_delta,
      "order_only_fixed_n_ge5"
    )
  ),
  use.names = TRUE
)

extremes <- rbindlist(
  list(
    rank_extreme_patient_effects(
      primary_all_delta,
      "primary_strict_all_n_gt_k"
    ),
    rank_extreme_patient_effects(
      primary_fixed_delta,
      "primary_strict_fixed_n_ge5"
    )
  ),
  use.names = TRUE
)

# ------------------------------------------------------------------------------
# 6. Cross-sensitivity comparison of the MAIN paired effect
# ------------------------------------------------------------------------------

main_effect <- contrast_summaries[
  contrast == "grid_vs_cluster" &
    metric %in% c(
      "recall_nonubiquitous_detection",
      "recall_heterogeneity_classification"
    )
]

strict_main <- main_effect[
  analysis == "primary_strict_all_n_gt_k",
  .(
    k,
    metric,
    strict_N = n_patients,
    strict_median = median_paired_advantage,
    strict_ci_lower = bootstrap_ci_lower,
    strict_ci_upper = bootstrap_ci_upper,
    strict_prop_positive = proportion_positive
  )
]

order_main <- main_effect[
  analysis == "order_only_all_n_gt_k",
  .(
    k,
    metric,
    order_N = n_patients,
    order_median = median_paired_advantage,
    order_ci_lower = bootstrap_ci_lower,
    order_ci_upper = bootstrap_ci_upper,
    order_prop_positive = proportion_positive
  )
]

strict_vs_order <- merge(
  strict_main,
  order_main,
  by = c("k", "metric"),
  all = TRUE,
  sort = TRUE
)
strict_vs_order[
  ,
  order_minus_strict_median :=
    order_median - strict_median
]

fixed_main <- main_effect[
  analysis == "primary_strict_fixed_n_ge5",
  .(
    k,
    metric,
    fixed_N = n_patients,
    fixed_median = median_paired_advantage,
    fixed_ci_lower = bootstrap_ci_lower,
    fixed_ci_upper = bootstrap_ci_upper,
    fixed_prop_positive = proportion_positive
  )
]

strict_vs_fixed <- merge(
  strict_main,
  fixed_main,
  by = c("k", "metric"),
  all = TRUE,
  sort = TRUE
)
strict_vs_fixed[
  ,
  fixed_minus_varying_median :=
    fixed_median - strict_median
]

# ------------------------------------------------------------------------------
# 7. Explicit QC: CIs for primary/main endpoint must be generated at k=2..4
# ------------------------------------------------------------------------------

primary_main_rows <- contrast_summaries[
  analysis == "primary_strict_all_n_gt_k" &
    contrast == "grid_vs_cluster" &
    metric %in% c(
      "recall_nonubiquitous_detection",
      "recall_heterogeneity_classification"
    )
]

if (nrow(primary_main_rows) != 6L) {
  stop("Expected six primary grid-vs-cluster main-endpoint rows.")
}

if (anyNA(primary_main_rows$bootstrap_ci_lower) ||
    anyNA(primary_main_rows$bootstrap_ci_upper)) {
  stop("Primary k=2..4 bootstrap CI generation FAILED.")
}

# k=5 is intentionally descriptive without a bootstrap CI.
k5_main_rows <- contrast_summaries[
  analysis == "primary_strict_k5" &
    contrast == "grid_vs_cluster"
]

if (any(!is.na(k5_main_rows$bootstrap_ci_lower)) ||
    any(!is.na(k5_main_rows$bootstrap_ci_upper))) {
  stop("k=5 exploratory rows unexpectedly received inferential CIs.")
}

# ------------------------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------------------------

if (DEBUG) {
  fwrite(
    contrast_summaries,
    file.path(
      spatial_dir,
      "07_hcc_spatial_paired_contrast_summary_all_filtered.tsv"
    ),
    sep = "\t"
  )
}

fwrite(
  absolute_summaries,
  file.path(
    spatial_dir,
    "07_hcc_spatial_absolute_design_summary_all_filtered.tsv"
  ),
  sep = "\t"
)

fwrite(
  main_effect,
  file.path(
    spatial_dir,
    "07_hcc_spatial_main_effect_summary.tsv"
  ),
  sep = "\t"
)

if (DEBUG) {
  fwrite(
    strict_vs_order,
    file.path(
      spatial_dir,
      "07_hcc_spatial_strict_vs_order_sensitivity.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    strict_vs_fixed,
    file.path(
      spatial_dir,
      "07_hcc_spatial_varying_vs_fixed_n_ge5.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    extremes,
    file.path(
      spatial_dir,
      "07_hcc_spatial_extreme_patient_effects.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    N_qc,
    file.path(PATHS$qc_dir, "07_spatial_N_by_k_qc.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    influence,
    file.path(PATHS$qc_dir, "07_spatial_leave_one_out_influence.tsv"),
    sep = "\t"
  )
}


# ------------------------------------------------------------------------------
# 9. Human-readable summary
# ------------------------------------------------------------------------------

get_main <- function(analysis_id, metric_id, kval) {
  z <- contrast_summaries[
    analysis == analysis_id &
      contrast == "grid_vs_cluster" &
      metric == metric_id &
      k == kval
  ]
  if (nrow(z) != 1L) {
    stop(
      "Could not retrieve unique summary row for ",
      analysis_id, " / ", metric_id, " / k=", kval
    )
  }
  z
}

summary_lines <- c(
  "Script 07: spatial analysis",
  "==========================================",
  "",
  "Bootstrap method for inferential k=2..4 analyses: exact empirical nonparametric median distribution",
  sprintf(
    "Bootstrap CI: %.0f%% percentile; no Monte-Carlo simulation error",
    100 * CI_LEVEL
  ),
  "",
  "Primary spatial contrast:",
  "  grid-based dispersed minus clustered/adjacent",
  "  positive advantage = grid/dispersed performs better",
  "",
  "Primary strict spatial cohort:"
)

for (kval in PRIMARY_K) {
  a <- get_main(
    "primary_strict_all_n_gt_k",
    "recall_nonubiquitous_detection",
    kval
  )
  b <- get_main(
    "primary_strict_all_n_gt_k",
    "recall_heterogeneity_classification",
    kval
  )

  summary_lines <- c(
    summary_lines,
    sprintf(
      paste0(
        "  k=%d, N=%d: detection median paired advantage = %.4f ",
        "[95%% CI %.4f, %.4f]; %.1f%% patients positive"
      ),
      kval,
      a$n_patients,
      a$median_paired_advantage,
      a$bootstrap_ci_lower,
      a$bootstrap_ci_upper,
      100 * a$proportion_positive
    ),
    sprintf(
      paste0(
        "             heterogeneity-classification median paired advantage = %.4f ",
        "[95%% CI %.4f, %.4f]; %.1f%% patients positive"
      ),
      b$median_paired_advantage,
      b$bootstrap_ci_lower,
      b$bootstrap_ci_upper,
      100 * b$proportion_positive
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "Strict fixed n>=5 robustness:"
)

for (kval in PRIMARY_K) {
  a <- get_main(
    "primary_strict_fixed_n_ge5",
    "recall_nonubiquitous_detection",
    kval
  )
  b <- get_main(
    "primary_strict_fixed_n_ge5",
    "recall_heterogeneity_classification",
    kval
  )

  summary_lines <- c(
    summary_lines,
    sprintf(
      "  k=%d, N=%d: detection median advantage = %.4f; classification median advantage = %.4f",
      kval,
      a$n_patients,
      a$median_paired_advantage,
      b$median_paired_advantage
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "Order-only sensitivity:"
)

for (kval in PRIMARY_K) {
  a <- get_main(
    "order_only_all_n_gt_k",
    "recall_nonubiquitous_detection",
    kval
  )
  b <- get_main(
    "order_only_all_n_gt_k",
    "recall_heterogeneity_classification",
    kval
  )

  summary_lines <- c(
    summary_lines,
    sprintf(
      "  k=%d, N=%d: detection median advantage = %.4f; classification median advantage = %.4f",
      kval,
      a$n_patients,
      a$median_paired_advantage,
      b$median_paired_advantage
    )
  )
}

a5 <- get_main(
  "primary_strict_k5",
  "recall_nonubiquitous_detection",
  5L
)
b5 <- get_main(
  "primary_strict_k5",
  "recall_heterogeneity_classification",
  5L
)

max_loo_detection <- max(
  influence[
    analysis == "primary_strict_all_n_gt_k" &
      metric == "recall_nonubiquitous_detection",
    max_abs_leave_one_out_change
  ],
  na.rm = TRUE
)

max_loo_classification <- max(
  influence[
    analysis == "primary_strict_all_n_gt_k" &
      metric == "recall_heterogeneity_classification",
    max_abs_leave_one_out_change
  ],
  na.rm = TRUE
)

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  summary_lines,
  "",
  "k=5 exploratory only:",
  sprintf(
    "  N=%d: detection median advantage = %.4f; classification median advantage = %.4f",
    a5$n_patients,
    a5$median_paired_advantage,
    b5$median_paired_advantage
  ),
  "  No bootstrap CI is reported for k=5 because N=4.",
  "",
  "Influence QC in primary strict analysis:",
  sprintf(
    "  maximum leave-one-patient-out change in detection median advantage = %.6f",
    max_loo_detection
  ),
  sprintf(
    "  maximum leave-one-patient-out change in classification median advantage = %.6f",
    max_loo_classification
  ),
  "",
  "Cohort-size QC: PASS",
  "Fixed-patient-set QC: PASS",
  "",
  "Interpretation guardrails:",
  "  - patient is the statistical unit; subset windows are not biological replicates;",
  "  - the primary effect is the median WITHIN-PATIENT paired grid-vs-cluster advantage;",
  "  - mean paired advantage is reported as a heterogeneity diagnostic, not substituted for the primary median;",
  "  - k=5 is descriptive/exploratory only (N=4);",
  "  - physical distances are not inferred from T-label differences;",
  "  - no continuous-distance model is run in Script 07.",
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed)
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "07_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
