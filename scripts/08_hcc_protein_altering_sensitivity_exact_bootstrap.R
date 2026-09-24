# ==============================================================================
# 08_hcc_protein_altering_sensitivity_exact_bootstrap.R
#
# Purpose
#   Prespecified protein-altering sensitivity analysis using the event flag
#   constructed and QC-validated in Script 02:
#
#       protein_altering == TRUE  <=>  IMPACT in {HIGH, MODERATE}
#
# This is NOT a new primary analysis. It asks whether the main conclusions from
# the all-filtered SNV/indel analysis remain qualitatively consistent when the
# event universe is restricted to protein-altering variants.
#
# The script:
#   1. Reuses the validated event-presence/reference objects from Script 02.
#   2. Enumerates all 4,590 HCC subsets exactly for the restricted event set.
#   3. Repeats the key depth summaries:
#        - available cohort, k=1..5
#        - strict n>k, k=1..5
#        - fixed n>=5, k=1..4
#        - n>=6 exploratory, k=1..5
#   4. Reuses the EXACT spatial selections frozen in Script 06 by joining their
#      canonical subset signatures to the protein-altering exhaustive universe.
#   5. Repeats the key grid-vs-cluster spatial contrasts.
#   6. Directly compares the sensitivity results with the all-filtered results
#      from Scripts 05 and 07.
#
# No variant re-calling, driver-only filtering, CNA/SV analysis, or additional
# annotation model is introduced here.
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
source(file.path(script_dir, "functions_bootstrap_exact.R"))

cat("Script 08: protein-altering sensitivity\n")
cat("======================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
sensitivity_dir <- file.path(PATHS$output_dir, "results", "hcc_sensitivity")
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- constants -----------------------------------
B <- NA_integer_  # exact bootstrap: no Monte-Carlo replicates
BASE_SEED <- as.integer(PARAMS$seed)
CI_LEVEL <- 0.95
TOLERANCE <- 1e-9

EXPECTED_PATIENTS_N_GE2 <- 122L
EXPECTED_TOTAL_SUBSETS <- 4590L
EXPECTED_PATIENT_K_ROWS <- 489L
EXPECTED_PROTEIN_EVENTS <- 36264L

EXPECTED_AVAILABLE_N <- c(`1` = 122L, `2` = 122L, `3` = 100L, `4` = 76L, `5` = 55L)
EXPECTED_STRICT_N <- c(`1` = 122L, `2` = 100L, `3` = 76L, `4` = 55L, `5` = 7L)
EXPECTED_PRIMARY_SPATIAL_N <- c(`2` = 80L, `3` = 60L, `4` = 44L, `5` = 4L)
EXPECTED_ORDER_SPATIAL_N <- c(`2` = 86L, `3` = 65L, `4` = 45L, `5` = 4L)

METRICS <- c(
  "recall_nonubiquitous_detection",
  "recall_heterogeneity_classification",
  "apparent_ubiquity_error",
  "recall_private",
  "nonubiquitous_fraction_abs_error",
  "jaccard_ith_abs_error"
)

# Script 05 uses manuscript-facing metric labels that differ slightly from the
# internal function column names. Keep the protein-altering DEPTH output on the
# exact same naming convention so it can be joined to Script 05 without any
# implicit renaming or post-hoc guesswork.
DEPTH_METRIC_LABELS <- c(
  recall_nonubiquitous_detection = "nonubiquitous_detection_recall",
  recall_heterogeneity_classification = "heterogeneity_classification_recall",
  apparent_ubiquity_error = "apparent_ubiquity_error",
  recall_private = "private_recall",
  nonubiquitous_fraction_abs_error = "nonubiquitous_fraction_abs_error",
  jaccard_ith_abs_error = "jaccard_ith_abs_error"
)

HIGHER_BETTER <- c(
  "recall_nonubiquitous_detection",
  "recall_heterogeneity_classification",
  "recall_private"
)

LOWER_BETTER <- c(
  "apparent_ubiquity_error",
  "nonubiquitous_fraction_abs_error",
  "jaccard_ith_abs_error"
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
  as.numeric(stats::quantile(v, p, names = FALSE, type = 7))
}

bootstrap_median_ci <- function(v, B, seed, level = 0.95, min_n = 20L) {
  exact_bootstrap_median_ci(
    v = v,
    level = level,
    min_n = min_n
  )
}

canonical_sample_signature <- function(ids) {
  ids <- as.character(ids)
  if (length(ids) == 0L || anyNA(ids) || any(!nzchar(ids))) {
    stop("Invalid sample IDs in subset signature.")
  }
  if (anyDuplicated(ids) > 0L) {
    stop("Duplicate sample IDs within subset signature.")
  }
  paste(sort(ids), collapse = "|")
}

subset_string_to_signature <- function(z) {
  ids <- strsplit(as.character(z), ",", fixed = TRUE)[[1L]]
  canonical_sample_signature(ids)
}

max_abs_diff_na_safe <- function(a, b) {
  if (length(a) != length(b)) return(Inf)
  aa <- is.na(a)
  bb <- is.na(b)
  if (any(xor(aa, bb))) return(Inf)
  keep <- !aa & !bb
  if (!any(keep)) return(0)
  max(abs(as.numeric(a[keep]) - as.numeric(b[keep])))
}

# Cohort-level median of patient-level expected performance.
summarise_depth_view <- function(patient_summary, analysis_id,
                                 cohort_definition, k_values, selector,
                                 seed_offset) {
  rows <- list()
  ri <- 0L

  for (kval in k_values) {
    keep <- selector(patient_summary, kval)
    d <- patient_summary[which(keep)]

    if (nrow(d) == 0L || any(d$k != kval)) {
      stop("Invalid depth selector for ", analysis_id, ", k=", kval, ".")
    }
    if (anyDuplicated(d$patient_id) > 0L) {
      stop("Duplicate patient in depth view ", analysis_id, ", k=", kval, ".")
    }

    for (mi in seq_along(METRICS)) {
      metric <- METRICS[[mi]]
      metric_label <- unname(DEPTH_METRIC_LABELS[[metric]])
      if (is.null(metric_label) || is.na(metric_label)) {
        stop("No canonical depth metric label defined for: ", metric)
      }

      col_name <- paste0(metric, "_mean")

      if (!col_name %in% names(d)) {
        stop("Missing patient-summary metric: ", col_name)
      }

      v <- d[[col_name]]
      ci <- bootstrap_median_ci(
        v,
        B = B,
        seed = BASE_SEED + seed_offset + kval * 100L + mi,
        level = CI_LEVEL
      )

      ri <- ri + 1L
      rows[[ri]] <- data.table(
        event_set = "protein_altering",
        analysis = analysis_id,
        cohort_definition = cohort_definition,
        k = as.integer(kval),
        n_patients = nrow(d),
        n_metric = sum(!is.na(v)),
        n_equal_k = sum(d$n == kval),
        fraction_equal_k = mean(d$n == kval),
        metric = metric_label,
        median = safe_median(v),
        q1 = safe_quantile(v, 0.25),
        q3 = safe_quantile(v, 0.75),
        bootstrap_ci_lower = ci[[1L]],
        bootstrap_ci_upper = ci[[2L]],
        bootstrap_replicates = NA_integer_,
        bootstrap_method =
          if (sum(!is.na(v)) >= 20L) "exact_nonparametric_percentile" else NA_character_,
        bootstrap_level =
          if (sum(!is.na(v)) >= 20L) CI_LEVEL else NA_real_
      )
    }
  }

  rbindlist(rows)
}

summarise_selected_design <- function(z, cohort_name) {
  if (nrow(z) == 0L) stop("Empty selected spatial table for ", cohort_name)

  z[
    ,
    c(
      list(
        n_design_subsets = .N
      ),
      lapply(.SD, safe_mean)
    ),
    by = .(patient_id, n, k, design),
    .SDcols = METRICS
  ][
    ,
    spatial_cohort := cohort_name
  ][]
}

make_spatial_deltas <- function(design_summary) {
  grid <- design_summary[design == "grid_dispersed"]
  cluster <- design_summary[design == "clustered_adjacent"]

  key <- c("spatial_cohort", "patient_id", "n", "k")

  g <- grid[, c(key, METRICS), with = FALSE]
  c <- cluster[, c(key, METRICS), with = FALSE]

  setnames(g, METRICS, paste0(METRICS, "__grid"))
  setnames(c, METRICS, paste0(METRICS, "__cluster"))

  z <- merge(g, c, by = key, all = TRUE, sort = TRUE)

  if (anyNA(z$patient_id) || nrow(z) != nrow(grid) || nrow(z) != nrow(cluster)) {
    stop("Protein spatial grid/cluster merge failed.")
  }

  for (metric in HIGHER_BETTER) {
    z[[paste0("grid_advantage__", metric)]] <-
      z[[paste0(metric, "__grid")]] -
      z[[paste0(metric, "__cluster")]]
  }

  for (metric in LOWER_BETTER) {
    z[[paste0("grid_advantage__", metric)]] <-
      z[[paste0(metric, "__cluster")]] -
      z[[paste0(metric, "__grid")]]
  }

  z
}

summarise_spatial_effect <- function(d, analysis_id, cohort_definition,
                                     k_values, selector, seed_offset) {
  rows <- list()
  ri <- 0L

  for (kval in k_values) {
    keep <- selector(d, kval)
    dk <- d[which(keep)]

    if (nrow(dk) == 0L || any(dk$k != kval)) {
      stop("Invalid spatial selector for ", analysis_id, ", k=", kval, ".")
    }
    if (anyDuplicated(dk$patient_id) > 0L) {
      stop("Duplicate patient in spatial analysis ", analysis_id, ".")
    }

    for (mi in seq_along(METRICS)) {
      metric <- METRICS[[mi]]
      col_name <- paste0("grid_advantage__", metric)
      v <- dk[[col_name]]
      vv <- v[!is.na(v)]

      ci <- bootstrap_median_ci(
        vv,
        B = B,
        seed = BASE_SEED + seed_offset + kval * 100L + mi,
        level = CI_LEVEL
      )

      ri <- ri + 1L
      rows[[ri]] <- data.table(
        event_set = "protein_altering",
        analysis = analysis_id,
        cohort_definition = cohort_definition,
        contrast = "grid_vs_cluster",
        k = as.integer(kval),
        n_patients = nrow(dk),
        n_metric = length(vv),
        metric = metric,
        median_paired_advantage = safe_median(vv),
        q1_paired_advantage = safe_quantile(vv, 0.25),
        q3_paired_advantage = safe_quantile(vv, 0.75),
        mean_paired_advantage = safe_mean(vv),
        n_positive = sum(vv > TOLERANCE),
        n_negative = sum(vv < -TOLERANCE),
        n_zero = sum(abs(vv) <= TOLERANCE),
        proportion_positive = mean(vv > TOLERANCE),
        bootstrap_ci_lower = ci[[1L]],
        bootstrap_ci_upper = ci[[2L]],
        bootstrap_replicates = NA_integer_,
        bootstrap_method =
          if (length(vv) >= 20L) "exact_nonparametric_percentile" else NA_character_,
        sign_convention = "positive = grid/dispersed performs better"
      )
    }
  }

  rbindlist(rows)
}

# ---------------------------- require inputs ---------------------------------
sample_map_path <- file.path(
  PATHS$intermediate_dir, "hcc_sample_map.rds"
)
presence_path <- file.path(
  PATHS$intermediate_dir, "hcc_event_presence.rds"
)
reference_path <- file.path(
  PATHS$intermediate_dir, "hcc_reference_events.rds"
)

primary_selection_path <- file.path(
  PATHS$output_dir, "results", "hcc_spatial",
  "06_hcc_spatial_selected_subsets_primary_strict.tsv"
)
order_selection_path <- file.path(
  PATHS$output_dir, "results", "hcc_spatial",
  "06_hcc_spatial_selected_subsets_order_sensitivity.tsv"
)

all_filtered_depth_path <- file.path(
  PATHS$output_dir, "results", "hcc_depth",
  "05_hcc_depth_core_summary_all_filtered.tsv"
)
all_filtered_spatial_path <- file.path(
  PATHS$output_dir, "results", "hcc_spatial",
  "07_hcc_spatial_main_effect_summary.tsv"
)

required_paths <- c(
  sample_map_path,
  presence_path,
  reference_path,
  primary_selection_path,
  order_selection_path,
  all_filtered_depth_path,
  all_filtered_spatial_path
)

missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing validated prerequisite file(s):\n",
    paste0("  ", missing_paths, collapse = "\n")
  )
}

sample_map <- as.data.table(readRDS(sample_map_path))
event_presence_all <- as.data.table(readRDS(presence_path))
reference_all <- as.data.table(readRDS(reference_path))

# ---------------------- protein-altering reference/presence -------------------
if (!"protein_altering" %in% names(reference_all)) {
  stop("Reference object lacks protein_altering flag from Script 02.")
}

reference_protein <- reference_all[protein_altering == TRUE]

if (nrow(reference_protein) != EXPECTED_PROTEIN_EVENTS) {
  stop(
    "Protein-altering event-count QC failed: expected ",
    EXPECTED_PROTEIN_EVENTS, ", observed ", nrow(reference_protein), "."
  )
}

protein_keys <- unique(reference_protein$event_key)
event_presence_protein <- event_presence_all[event_key %in% protein_keys]

if (uniqueN(reference_protein$patient_id) != 123L) {
  stop("Protein-altering reference no longer includes all 123 source patients.")
}

# All 122 patients with n>=2 must have at least one non-ubiquitous event, based
# on the already audited Script 02 reference. Verify rather than assume.
protein_ref_by_patient <- reference_protein[
  ,
  .(
    n_events_reference = .N,
    n_nonubiquitous_reference = sum(is_nonubiquitous),
    n_private_reference = sum(is_private)
  ),
  by = patient_id
]

sample_patient_n <- unique(
  sample_map[, .(patient_id, n = as.integer(n_sectors))]
)

protein_ref_by_patient <- merge(
  protein_ref_by_patient,
  sample_patient_n,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

if (protein_ref_by_patient[n >= 2L & n_nonubiquitous_reference == 0L, .N] > 0L) {
  stop("At least one n>=2 patient has zero protein-altering non-ubiquitous events.")
}
if (protein_ref_by_patient[n >= 2L & n_private_reference == 0L, .N] > 0L) {
  stop("At least one n>=2 patient has zero protein-altering private events.")
}

setkey(sample_map, patient_id)
setkey(event_presence_protein, patient_id)
setkey(reference_protein, patient_id)

analysis_patients <- unique(
  sample_map[n_sectors >= 2L, .(patient_id, n = as.integer(n_sectors))]
)
setorder(analysis_patients, patient_id)

if (nrow(analysis_patients) != EXPECTED_PATIENTS_N_GE2) {
  stop("Expected 122 n>=2 patients in protein sensitivity.")
}

# ------------------------------------------------------------------------------
# 1. Exact exhaustive protein-altering sampling
# ------------------------------------------------------------------------------

cat(
  "Enumerating protein-altering subsets for ",
  nrow(analysis_patients), " patients...\n", sep = ""
)

subset_list <- vector("list", nrow(analysis_patients))
patient_summary_list <- vector("list", nrow(analysis_patients))
analytic_max_errors <- c(
  detection = 0,
  classification = 0,
  apparent_ubiquity = 0,
  private = 0
)

for (ii in seq_len(nrow(analysis_patients))) {
  pid <- analysis_patients$patient_id[[ii]]
  n <- analysis_patients$n[[ii]]

  sm <- sample_map[J(pid), .(sample_id), nomatch = 0L]
  setorder(sm, sample_id)
  sample_ids <- sm$sample_id

  pres <- event_presence_protein[
    J(pid),
    .(sample_id, event_key),
    nomatch = 0L
  ]

  ref <- reference_protein[
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

  pdata <- build_patient_incidence(
    patient_presence = pres,
    patient_reference = ref,
    sample_ids = sample_ids
  )

  patient_subset <- vector("list", n)
  patient_summary <- vector("list", n)

  m_nonubiq <- ref$full_count[ref$is_nonubiquitous]
  n_nonubiq <- length(m_nonubiq)
  n_private <- sum(ref$is_private)

  for (kval in seq_len(n)) {
    zz <- as.data.table(compute_all_subset_metrics(pdata, kval))
    zz[, patient_id := pid]
    zz[, n := n]
    zz[
      ,
      subset_signature := vapply(
        subset_samples,
        subset_string_to_signature,
        character(1)
      )
    ]

    patient_subset[[kval]] <- zz

    # Exact analytic expectations for the new restricted event set.
    denom <- choose(n, kval)
    p_not_detected <- choose(n - m_nonubiq, kval) / denom
    p_apparent <- choose(m_nonubiq, kval) / denom

    analytic_detection <- mean(1 - p_not_detected)
    analytic_apparent <- mean(p_apparent)
    analytic_classification <- mean(1 - p_not_detected - p_apparent)
    analytic_private <- if (n_private > 0L) kval / n else NA_real_

    observed_detection <- safe_mean(zz$recall_nonubiquitous_detection)
    observed_apparent <- safe_mean(zz$apparent_ubiquity_error)
    observed_classification <- safe_mean(
      zz$recall_heterogeneity_classification
    )
    observed_private <- safe_mean(zz$recall_private)

    analytic_max_errors[["detection"]] <- max(
      analytic_max_errors[["detection"]],
      abs(observed_detection - analytic_detection)
    )
    analytic_max_errors[["classification"]] <- max(
      analytic_max_errors[["classification"]],
      abs(observed_classification - analytic_classification)
    )
    analytic_max_errors[["apparent_ubiquity"]] <- max(
      analytic_max_errors[["apparent_ubiquity"]],
      abs(observed_apparent - analytic_apparent)
    )
    if (!is.na(observed_private) && !is.na(analytic_private)) {
      analytic_max_errors[["private"]] <- max(
        analytic_max_errors[["private"]],
        abs(observed_private - analytic_private)
      )
    }

    ps <- data.table(
      patient_id = pid,
      n = n,
      k = kval,
      n_subsets = nrow(zz)
    )

    for (metric in METRICS) {
      ps[[paste0(metric, "_mean")]] <- safe_mean(zz[[metric]])
      ps[[paste0(metric, "_median")]] <- safe_median(zz[[metric]])
    }

    patient_summary[[kval]] <- ps
  }

  subset_list[[ii]] <- rbindlist(patient_subset)
  patient_summary_list[[ii]] <- rbindlist(patient_summary)

  rm(
    sm, sample_ids, pres, ref, pdata,
    patient_subset, patient_summary, m_nonubiq
  )
  invisible(gc(FALSE))
}

protein_subsets <- rbindlist(subset_list, use.names = TRUE, fill = TRUE)
protein_patient_summary <- rbindlist(
  patient_summary_list,
  use.names = TRUE,
  fill = TRUE
)

rm(subset_list, patient_summary_list)
invisible(gc())

setorder(protein_subsets, patient_id, k, subset_id)
setorder(protein_patient_summary, patient_id, k)

if (nrow(protein_subsets) != EXPECTED_TOTAL_SUBSETS) {
  stop(
    "Protein exhaustive subset count QC FAILED: expected ",
    EXPECTED_TOTAL_SUBSETS, ", observed ", nrow(protein_subsets), "."
  )
}
if (nrow(protein_patient_summary) != EXPECTED_PATIENT_K_ROWS) {
  stop("Protein patient x k row count QC FAILED.")
}
if (max(analytic_max_errors) > TOLERANCE) {
  stop(
    "Protein analytic hypergeometric validation FAILED; max error = ",
    format(max(analytic_max_errors), scientific = TRUE)
  )
}

# Patient-level monotonicity for exact expectations.
for (metric_col in c(
  "recall_nonubiquitous_detection_mean",
  "recall_heterogeneity_classification_mean",
  "recall_private_mean"
)) {
  bad <- protein_patient_summary[
    order(k),
    any(diff(get(metric_col)) < -TOLERANCE),
    by = patient_id
  ][V1 == TRUE]

  if (nrow(bad) > 0L) {
    stop("Protein patient-level monotonicity failed for ", metric_col, ".")
  }
}

# ------------------------------------------------------------------------------
# 2. Protein-altering depth summaries
# ------------------------------------------------------------------------------

depth_summary <- rbindlist(
  list(
    summarise_depth_view(
      protein_patient_summary,
      "available_cohort",
      "all patients with n>=k; descriptive continuity view",
      1:5,
      function(d, kval) d$k == kval,
      10000L
    ),
    summarise_depth_view(
      protein_patient_summary,
      "strict_downsampling",
      "patients with n>k; no k=n ceiling points",
      1:5,
      function(d, kval) d$k == kval & d$n > kval,
      20000L
    ),
    summarise_depth_view(
      protein_patient_summary,
      "fixed_n_ge5",
      "same 55 patients with n>=5 at k=1..4",
      1:4,
      function(d, kval) d$k == kval & d$n >= 5L,
      30000L
    ),
    summarise_depth_view(
      protein_patient_summary,
      "n_ge6_exploratory",
      "same 7 patients with n>=6 at k=1..5",
      1:5,
      function(d, kval) d$k == kval & d$n >= 6L,
      40000L
    )
  )
)

# N(k) QC for key depth views.
get_N <- function(analysis_id) {
  z <- depth_summary[
    analysis == analysis_id &
      metric == "nonubiquitous_detection_recall",
    .(k, n_patients)
  ]
  setorder(z, k)
  z
}

if (!identical(
  as.integer(get_N("available_cohort")$n_patients),
  as.integer(EXPECTED_AVAILABLE_N)
)) {
  stop("Protein available-cohort N(k) QC FAILED.")
}
if (!identical(
  as.integer(get_N("strict_downsampling")$n_patients),
  as.integer(EXPECTED_STRICT_N)
)) {
  stop("Protein strict-downsampling N(k) QC FAILED.")
}

# ------------------------------------------------------------------------------
# 3. Reuse EXACT frozen Script 06 spatial selections
# ------------------------------------------------------------------------------

primary_sel <- fread(primary_selection_path)
order_sel <- fread(order_selection_path)

selection_required <- c(
  "patient_id", "n", "k", "design", "subset_signature"
)

if (!all(selection_required %in% names(primary_sel)) ||
    !all(selection_required %in% names(order_sel))) {
  stop("Script 06 selected-subset schema is not as expected.")
}

# Only grid and clustered selections are stored in Script 06.
if (!all(unique(primary_sel$design) %in%
         c("grid_dispersed", "clustered_adjacent"))) {
  stop("Unexpected design in primary Script 06 selection table.")
}
if (!all(unique(order_sel$design) %in%
         c("grid_dispersed", "clustered_adjacent"))) {
  stop("Unexpected design in order-only Script 06 selection table.")
}

# Select columns explicitly by NAME. The previous c(list(patient_id, ...),
# lapply(.SD, identity)) expression retained the values but did not guarantee
# the key-column names, which caused patient_id to be unavailable downstream.
protein_subset_for_join <- protein_subsets[
  ,
  c(
    "patient_id", "n", "k", "subset_signature",
    METRICS
  ),
  with = FALSE
]

if (anyDuplicated(
  protein_subset_for_join[
    ,
    .(patient_id, k, subset_signature)
  ]
) > 0L) {
  stop("Duplicate protein subset signatures before spatial join.")
}

join_spatial_selection <- function(sel, cohort_name) {
  sel2 <- copy(sel)
  sel2[, patient_id := as.character(patient_id)]

  z <- merge(
    sel2[
      ,
      .(
        patient_id,
        n,
        k,
        design,
        design_selection_id,
        subset_signature
      )
    ],
    protein_subset_for_join,
    by = c("patient_id", "n", "k", "subset_signature"),
    all.x = TRUE,
    sort = TRUE
  )

  if (nrow(z) != nrow(sel2)) {
    stop(cohort_name, ": protein spatial join changed row count.")
  }

  # All primary endpoint values should exist; all n>=2 patients have at least
  # one protein-altering non-ubiquitous event.
  if (anyNA(z$recall_nonubiquitous_detection)) {
    stop(cohort_name, ": unmatched/missing protein spatial subset detected.")
  }

  z[, spatial_cohort := cohort_name]
  z
}

primary_selected_protein <- join_spatial_selection(
  primary_sel, "primary_strict"
)
order_selected_protein <- join_spatial_selection(
  order_sel, "order_only_sensitivity"
)

primary_design_protein <- summarise_selected_design(
  primary_selected_protein, "primary_strict"
)
order_design_protein <- summarise_selected_design(
  order_selected_protein, "order_only_sensitivity"
)

primary_delta_protein <- make_spatial_deltas(primary_design_protein)
order_delta_protein <- make_spatial_deltas(order_design_protein)

# Spatial N QC.
get_spatial_N <- function(d) {
  z <- d[, .(n_patients = uniqueN(patient_id)), by = k]
  setorder(z, k)
  z
}

pN <- get_spatial_N(primary_delta_protein)
oN <- get_spatial_N(order_delta_protein)

if (!identical(as.integer(pN$n_patients),
               as.integer(EXPECTED_PRIMARY_SPATIAL_N))) {
  stop("Protein primary spatial N(k) QC FAILED.")
}
if (!identical(as.integer(oN$n_patients),
               as.integer(EXPECTED_ORDER_SPATIAL_N))) {
  stop("Protein order-only spatial N(k) QC FAILED.")
}

# ------------------------------------------------------------------------------
# 4. Protein-altering spatial effect summaries
# ------------------------------------------------------------------------------

spatial_summary <- rbindlist(
  list(
    summarise_spatial_effect(
      primary_delta_protein,
      "primary_strict_all_n_gt_k",
      "strict spatial cohort; all available n>k patients",
      2:4,
      function(d, kval) d$k == kval,
      50000L
    ),
    summarise_spatial_effect(
      primary_delta_protein,
      "primary_strict_fixed_n_ge5",
      "strict spatial cohort; same 44 n>=5 patients",
      2:4,
      function(d, kval) d$k == kval & d$n >= 5L,
      60000L
    ),
    summarise_spatial_effect(
      order_delta_protein,
      "order_only_all_n_gt_k",
      "order-only sensitivity cohort; all available n>k patients",
      2:4,
      function(d, kval) d$k == kval,
      70000L
    ),
    summarise_spatial_effect(
      order_delta_protein,
      "order_only_fixed_n_ge5",
      "order-only sensitivity cohort; same 45 n>=5 patients",
      2:4,
      function(d, kval) d$k == kval & d$n >= 5L,
      80000L
    ),
    summarise_spatial_effect(
      primary_delta_protein,
      "primary_strict_k5",
      "strict spatial cohort; k=5 and n>5, exploratory N=4",
      5L,
      function(d, kval) d$k == kval,
      90000L
    )
  )
)

# ------------------------------------------------------------------------------
# 4b. Core sensitivity checkpoint
# ------------------------------------------------------------------------------

# Save the validated restricted exhaustive objects BEFORE cross-event-set
# comparison/reporting. If a later presentation-level join fails, the expensive
# core sensitivity calculation remains recoverable and inspectable.

# ------------------------------------------------------------------------------
# 5. Direct comparison with all-filtered results
# ------------------------------------------------------------------------------

all_depth <- fread(all_filtered_depth_path)
all_spatial <- fread(all_filtered_spatial_path)

# Fail early if depth naming conventions are no longer aligned.
expected_depth_labels <- unname(DEPTH_METRIC_LABELS)
if (!all(expected_depth_labels %in% unique(depth_summary$metric))) {
  stop(
    "Protein depth-summary metric naming QC FAILED. Missing: ",
    paste(
      setdiff(expected_depth_labels, unique(depth_summary$metric)),
      collapse = ", "
    )
  )
}
if (!all(
  c(
    "nonubiquitous_detection_recall",
    "heterogeneity_classification_recall"
  ) %in% unique(all_depth$metric)
)) {
  stop("All-filtered Script 05 depth metric naming QC FAILED.")
}

depth_compare_metrics <- c(
  "nonubiquitous_detection_recall",
  "heterogeneity_classification_recall"
)

protein_depth_main <- depth_summary[
  analysis %in% c("fixed_n_ge5", "n_ge6_exploratory") &
    metric %in% depth_compare_metrics,
  .(
    analysis,
    k,
    metric,
    protein_N = n_patients,
    protein_median = median,
    protein_ci_lower = bootstrap_ci_lower,
    protein_ci_upper = bootstrap_ci_upper
  )
]

all_depth_main <- all_depth[
  analysis %in% c("fixed_n_ge5", "n_ge6_exploratory") &
    metric %in% depth_compare_metrics,
  .(
    analysis,
    k,
    metric,
    all_filtered_N = n_patients,
    all_filtered_median = median,
    all_filtered_ci_lower = bootstrap_ci_lower,
    all_filtered_ci_upper = bootstrap_ci_upper
  )
]

depth_comparison <- merge(
  protein_depth_main,
  all_depth_main,
  by = c("analysis", "k", "metric"),
  all = TRUE,
  sort = TRUE
)
depth_comparison[
  ,
  protein_minus_all_filtered :=
    protein_median - all_filtered_median
]

if (nrow(depth_comparison) != nrow(protein_depth_main) ||
    nrow(depth_comparison) != nrow(all_depth_main) ||
    anyNA(depth_comparison$protein_N) ||
    anyNA(depth_comparison$all_filtered_N) ||
    any(depth_comparison$protein_N != depth_comparison$all_filtered_N)) {

  if (DEBUG) {
    fwrite(
      depth_comparison,
      file.path(
        PATHS$qc_dir,
        "08_depth_comparison_alignment_failure.tsv"
      ),
      sep = "\t"
    )
  }

  stop(
    "Protein-vs-all-filtered depth comparison alignment QC FAILED. ",
    "Diagnostic table written to 08_depth_comparison_alignment_failure.tsv."
  )
}

protein_spatial_main <- spatial_summary[
  analysis %in% c(
    "primary_strict_all_n_gt_k",
    "primary_strict_fixed_n_ge5",
    "order_only_all_n_gt_k"
  ) &
    metric %in% c(
      "recall_nonubiquitous_detection",
      "recall_heterogeneity_classification"
    ),
  .(
    analysis,
    k,
    metric,
    protein_N = n_patients,
    protein_median_advantage = median_paired_advantage,
    protein_ci_lower = bootstrap_ci_lower,
    protein_ci_upper = bootstrap_ci_upper,
    protein_prop_positive = proportion_positive
  )
]

all_spatial_main2 <- all_spatial[
  analysis %in% c(
    "primary_strict_all_n_gt_k",
    "primary_strict_fixed_n_ge5",
    "order_only_all_n_gt_k"
  ) &
    metric %in% c(
      "recall_nonubiquitous_detection",
      "recall_heterogeneity_classification"
    ),
  .(
    analysis,
    k,
    metric,
    all_filtered_N = n_patients,
    all_filtered_median_advantage = median_paired_advantage,
    all_filtered_ci_lower = bootstrap_ci_lower,
    all_filtered_ci_upper = bootstrap_ci_upper,
    all_filtered_prop_positive = proportion_positive
  )
]

spatial_comparison <- merge(
  protein_spatial_main,
  all_spatial_main2,
  by = c("analysis", "k", "metric"),
  all = TRUE,
  sort = TRUE
)
spatial_comparison[
  ,
  protein_minus_all_filtered_advantage :=
    protein_median_advantage - all_filtered_median_advantage
]

if (nrow(spatial_comparison) != nrow(protein_spatial_main) ||
    anyNA(spatial_comparison$protein_N) ||
    anyNA(spatial_comparison$all_filtered_N) ||
    any(spatial_comparison$protein_N != spatial_comparison$all_filtered_N)) {
  stop("Protein-vs-all-filtered spatial comparison alignment QC FAILED.")
}

# ------------------------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------------------------


saveRDS(
  protein_patient_summary,
  file.path(
    PATHS$intermediate_dir,
    "hcc_unconstrained_patient_summary_protein_altering.rds"
  )
)


if (DEBUG) {
  fwrite(
    protein_patient_summary,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_altering_patient_summary.tsv"
    ),
    sep = "\t"
  )
}
fwrite(
  depth_summary,
  file.path(
    sensitivity_dir,
    "08_hcc_protein_altering_depth_summary.tsv"
  ),
  sep = "\t"
)
if (DEBUG) {
  fwrite(
    spatial_summary,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_altering_spatial_summary.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    primary_design_protein,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_altering_spatial_patient_design_summary_primary_strict.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    primary_delta_protein,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_altering_spatial_patient_deltas_primary_strict.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    order_design_protein,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_altering_spatial_patient_design_summary_order_sensitivity.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    order_delta_protein,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_altering_spatial_patient_deltas_order_sensitivity.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    depth_comparison,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_vs_all_filtered_depth_comparison.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    spatial_comparison,
    file.path(
      sensitivity_dir,
      "08_hcc_protein_vs_all_filtered_spatial_comparison.tsv"
    ),
    sep = "\t"
  )
}

# QC audit.
analytic_qc <- data.table(
  endpoint = names(analytic_max_errors),
  max_absolute_error = as.numeric(analytic_max_errors),
  pass = as.numeric(analytic_max_errors) <= TOLERANCE
)
if (DEBUG) {
  fwrite(
    analytic_qc,
    file.path(PATHS$qc_dir, "08_protein_analytic_validation.tsv"),
    sep = "\t"
  )
}

# ------------------------------------------------------------------------------
# 7. Human-readable summary
# ------------------------------------------------------------------------------

get_depth <- function(analysis_id, metric_id, kval) {
  z <- depth_summary[
    analysis == analysis_id &
      metric == metric_id &
      k == kval
  ]
  if (nrow(z) != 1L) {
    stop("Could not retrieve unique protein depth summary row.")
  }
  z
}

get_spatial <- function(analysis_id, metric_id, kval) {
  z <- spatial_summary[
    analysis == analysis_id &
      metric == metric_id &
      k == kval
  ]
  if (nrow(z) != 1L) {
    stop("Could not retrieve unique protein spatial summary row.")
  }
  z
}

d4 <- get_depth(
  "fixed_n_ge5",
  "nonubiquitous_detection_recall",
  4L
)
h4 <- get_depth(
  "fixed_n_ge5",
  "heterogeneity_classification_recall",
  4L
)
d5 <- get_depth(
  "n_ge6_exploratory",
  "nonubiquitous_detection_recall",
  5L
)

summary_lines <- c(
  "Script 08: protein-altering sensitivity",
  "======================================================",
  "",
  sprintf(
    "Protein-altering patient-specific events: %d",
    nrow(reference_protein)
  ),
  "Definition: Script 02 IMPACT in {HIGH, MODERATE}",
  "",
  sprintf(
    "Exhaustive subset rows: %d (expected %d)",
    nrow(protein_subsets), EXPECTED_TOTAL_SUBSETS
  ),
  sprintf(
    "Patient x k rows: %d (expected %d)",
    nrow(protein_patient_summary), EXPECTED_PATIENT_K_ROWS
  ),
  "Analytic hypergeometric validation: PASS",
  sprintf(
    "  maximum endpoint error: %.3g",
    max(analytic_max_errors)
  ),
  "Patient-level monotonicity QC: PASS",
  "",
  "Protein-altering depth sensitivity:",
  sprintf(
    paste0(
      "  fixed n>=5 k=4: detection median = %.4f ",
      "[95%% CI %.4f, %.4f], N=%d"
    ),
    d4$median, d4$bootstrap_ci_lower, d4$bootstrap_ci_upper,
    d4$n_patients
  ),
  sprintf(
    paste0(
      "  fixed n>=5 k=4: heterogeneity-classification median = %.4f ",
      "[95%% CI %.4f, %.4f]"
    ),
    h4$median, h4$bootstrap_ci_lower, h4$bootstrap_ci_upper
  ),
  sprintf(
    "  n>=6 exploratory k=5: detection median = %.4f, N=%d",
    d5$median, d5$n_patients
  ),
  "",
  "Protein-altering primary strict spatial sensitivity:"
)

for (kval in 2:4) {
  a <- get_spatial(
    "primary_strict_all_n_gt_k",
    "recall_nonubiquitous_detection",
    kval
  )
  b <- get_spatial(
    "primary_strict_all_n_gt_k",
    "recall_heterogeneity_classification",
    kval
  )

  summary_lines <- c(
    summary_lines,
    sprintf(
      paste0(
        "  k=%d, N=%d: grid-vs-cluster detection median advantage = %.4f ",
        "[95%% CI %.4f, %.4f]; %.1f%% positive"
      ),
      kval, a$n_patients, a$median_paired_advantage,
      a$bootstrap_ci_lower, a$bootstrap_ci_upper,
      100 * a$proportion_positive
    ),
    sprintf(
      paste0(
        "             classification median advantage = %.4f ",
        "[95%% CI %.4f, %.4f]; %.1f%% positive"
      ),
      b$median_paired_advantage,
      b$bootstrap_ci_lower, b$bootstrap_ci_upper,
      100 * b$proportion_positive
    )
  )
}

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  summary_lines,
  "",
  "Comparison tables against the all-filtered primary analysis were written.",
  "Interpret qualitative concordance only after reviewing those tables;",
  "this script does not impose a post-hoc concordance threshold.",
  "",
  "If the main depth and spatial conclusions remain qualitatively concordant,",
  "no driver-only or additional functional filtering is required.",
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed)
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "08_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
