# ==============================================================================
# 06_hcc_spatial_designs.R
#
# Purpose
#   Construct and validate the HCC spatial sampling designs using the already
#   validated exhaustive all-filtered subset table from Script 04.
#
# IMPORTANT: this script does NOT recompute genomic metrics. It selects exact
# subsets from the canonical exhaustive universe, which guarantees that every
# grid/clustered result is evaluated with the same event definitions and metric
# functions as the unconstrained analysis.
#
# Spatial definitions frozen before outcome analysis
#
#   PRIMARY spatial cohort:
#     - patient-level spatial_primary_eligible == TRUE from Script 01
#     - no MF label
#     - exactly one T coordinate per sector
#     - no duplicate/composite/missing T coordinate
#     - consecutive T coordinates
#
#   ORDER-ONLY sensitivity cohort:
#     - spatial_order_eligible == TRUE
#     - same exclusions, but clean gaps in T numbering are allowed
#     - selection still uses ORDER/RANK only; no metric physical distance is
#       inferred from differences in T labels.
#
#   GRID-BASED DISPERSED:
#     - ordered sector ranks 1..n
#     - k >= 2
#     - span the two endpoints
#     - among endpoint-spanning subsets, minimize squared deviation from k
#       equally spaced ideal rank positions
#     - retain/average all mathematically tied grid solutions
#
#   CLUSTERED/ADJACENT:
#     - every contiguous k-sector window in ordered sector rank
#     - mean across all n-k+1 possible windows = location-agnostic expected
#       performance of an adjacent/clustered design
#
#   UNCONSTRAINED:
#     - all choose(n,k) subsets
#     - mean across all subsets, already validated in Script 04
#
# Spatial contrast is evaluated only when k < n.
# Planned k:
#   k=2,3,4 = primary spatial range
#   k=5     = exploratory only (very small N)
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

cat("Script 06: spatial design construction\n")
cat("=====================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
spatial_dir <- file.path(PATHS$output_dir, "results", "hcc_spatial")
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- constants -----------------------------------
K_VALUES <- 2:5
TOLERANCE <- 1e-10

# Fixed QC expectations derived from the already validated Script 01 metadata.
EXPECTED_PRIMARY_N_BY_K <- c(`2` = 80L, `3` = 60L, `4` = 44L, `5` = 4L)
EXPECTED_ORDER_N_BY_K <- c(`2` = 86L, `3` = 65L, `4` = 45L, `5` = 4L)

# Metrics inherited verbatim from the validated exhaustive subset table.
METRIC_COLS <- c(
  "recall_nonubiquitous_detection",
  "recall_heterogeneity_classification",
  "apparent_ubiquity_error",
  "conditional_apparent_ubiquity_rate",
  "recall_private",
  "nonubiquitous_fraction_abs_error",
  "jaccard_ith_abs_error"
)

# Positive paired deltas will always mean "grid performs better".
HIGHER_BETTER_METRICS <- c(
  "recall_nonubiquitous_detection",
  "recall_heterogeneity_classification",
  "recall_private"
)

LOWER_BETTER_METRICS <- c(
  "apparent_ubiquity_error",
  "nonubiquitous_fraction_abs_error",
  "jaccard_ith_abs_error"
)

# conditional_apparent_ubiquity_rate is retained in the design summaries but is
# not promoted to the paired-delta table because it is a conditional supportive
# quantity and can be undefined when no non-ubiquitous event is detected.

# ------------------------------- helpers -------------------------------------
safe_mean <- function(v) {
  if (length(v) == 0L || all(is.na(v))) return(NA_real_)
  mean(v, na.rm = TRUE)
}

safe_median <- function(v) {
  if (length(v) == 0L || all(is.na(v))) return(NA_real_)
  stats::median(v, na.rm = TRUE)
}

safe_min <- function(v) {
  if (length(v) == 0L || all(is.na(v))) return(NA_real_)
  min(v, na.rm = TRUE)
}

safe_max <- function(v) {
  if (length(v) == 0L || all(is.na(v))) return(NA_real_)
  max(v, na.rm = TRUE)
}

canonical_sample_signature <- function(ids) {
  ids <- as.character(ids)
  if (length(ids) == 0L || anyNA(ids) || any(!nzchar(ids))) {
    stop("Invalid sample ID while constructing a subset signature.")
  }
  if (anyDuplicated(ids) > 0L) {
    stop("Duplicate sample ID within one subset signature.")
  }
  paste(sort(ids), collapse = "|")
}

subset_string_to_signature <- function(z) {
  ids <- strsplit(as.character(z), ",", fixed = TRUE)[[1L]]
  canonical_sample_signature(ids)
}

rank_matrix_to_selection_table <- function(rank_matrix, ordered_samples,
                                           ordered_t, design) {
  if (is.null(dim(rank_matrix))) {
    rank_matrix <- matrix(rank_matrix, ncol = 1L)
  }

  ans <- vector("list", ncol(rank_matrix))

  for (jj in seq_len(ncol(rank_matrix))) {
    rr <- as.integer(rank_matrix[, jj])
    samples <- ordered_samples[rr]
    tvals <- ordered_t[rr]

    ans[[jj]] <- data.table(
      design = design,
      design_selection_id = jj,
      rank_indices = paste(rr, collapse = ","),
      selected_samples_spatial_order = paste(samples, collapse = ","),
      subset_signature = canonical_sample_signature(samples),
      rank_span = max(rr) - min(rr),
      selected_t_values = paste(tvals, collapse = ","),
      t_label_span = max(tvals) - min(tvals)
    )
  }

  rbindlist(ans)
}

summarise_design_rows <- function(selected_rows, patient_id, n, k, design,
                                  spatial_cohort) {
  if (nrow(selected_rows) == 0L) {
    stop(
      "No selected rows supplied for ", patient_id,
      ", k=", k, ", design=", design, "."
    )
  }

  ans <- data.table(
    spatial_cohort = spatial_cohort,
    patient_id = patient_id,
    n = as.integer(n),
    k = as.integer(k),
    design = design,
    n_design_subsets = nrow(selected_rows)
  )

  for (cc in METRIC_COLS) {
    ans[[paste0(cc, "_mean")]] <- safe_mean(selected_rows[[cc]])
    ans[[paste0(cc, "_median")]] <- safe_median(selected_rows[[cc]])
    ans[[paste0(cc, "_min")]] <- safe_min(selected_rows[[cc]])
    ans[[paste0(cc, "_max")]] <- safe_max(selected_rows[[cc]])
  }

  ans
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

# ---------------------------- require inputs ---------------------------------
sample_map_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_sample_map.rds"
)
subset_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_exhaustive_subsets_all_filtered.rds"
)
unconstrained_summary_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_unconstrained_patient_summary_all_filtered.rds"
)

required_paths <- c(
  sample_map_path,
  subset_path,
  unconstrained_summary_path
)

missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing validated input file(s):\n",
    paste0("  ", missing_paths, collapse = "\n")
  )
}

sample_map <- as.data.table(readRDS(sample_map_path))
subset_results <- as.data.table(readRDS(subset_path))
unconstrained_summary <- as.data.table(readRDS(unconstrained_summary_path))

required_sample_cols <- c(
  "patient_id", "sample_id", "sector_num", "n_sectors",
  "spatial_primary_eligible", "spatial_order_eligible"
)

required_subset_cols <- c(
  "patient_id", "n", "k", "subset_id", "subset_samples",
  METRIC_COLS
)

required_unconstrained_cols <- c(
  "patient_id", "n", "k",
  "recall_nonubiquitous_detection_mean",
  "recall_heterogeneity_classification_mean",
  "apparent_ubiquity_error_mean",
  "conditional_apparent_ubiquity_rate_mean",
  "recall_private_mean",
  "nonubiquitous_fraction_abs_error_mean",
  "jaccard_ith_abs_error_mean"
)

if (!all(required_sample_cols %in% names(sample_map))) {
  stop(
    "Sample map lacks required columns: ",
    paste(setdiff(required_sample_cols, names(sample_map)), collapse = ", ")
  )
}
if (!all(required_subset_cols %in% names(subset_results))) {
  stop(
    "Exhaustive subset table lacks required columns: ",
    paste(setdiff(required_subset_cols, names(subset_results)), collapse = ", ")
  )
}
if (!all(required_unconstrained_cols %in% names(unconstrained_summary))) {
  stop(
    "Unconstrained patient summary lacks required columns: ",
    paste(
      setdiff(required_unconstrained_cols, names(unconstrained_summary)),
      collapse = ", "
    )
  )
}

sample_map[, patient_id := as.character(patient_id)]
sample_map[, sample_id := as.character(sample_id)]
subset_results[, patient_id := as.character(patient_id)]
unconstrained_summary[, patient_id := as.character(patient_id)]

# A comma is used only as a serialized delimiter in Script 04.
if (any(grepl(",", sample_map$sample_id, fixed = TRUE))) {
  stop("At least one sample ID contains a comma; subset parsing would be unsafe.")
}

# Add canonical sample-set signature to the exhaustive universe.
subset_results[
  ,
  subset_signature := vapply(
    subset_samples,
    subset_string_to_signature,
    character(1)
  )
]

if (anyDuplicated(
  subset_results[, .(patient_id, k, subset_signature)]
) > 0L) {
  stop(
    "Duplicate patient x k x subset-signature rows detected in exhaustive table."
  )
}

setkey(subset_results, patient_id, k)
setkey(unconstrained_summary, patient_id, k)

# ------------------------------------------------------------------------------
# 1. Construct one spatial cohort
# ------------------------------------------------------------------------------

construct_spatial_cohort <- function(
    cohort_name,
    eligibility_column,
    expected_n_by_k,
    require_consecutive) {

  cat("Constructing spatial cohort: ", cohort_name, "\n", sep = "")

  patient_flags <- unique(
    sample_map[
      ,
      .(
        patient_id,
        n = as.integer(n_sectors),
        eligible = as.logical(get(eligibility_column))
      )
    ]
  )

  eligible_patients <- patient_flags[eligible == TRUE & n >= 2L]
  setorder(eligible_patients, patient_id)

  selected_subset_rows <- list()
  patient_design_rows <- list()
  selection_qc_rows <- list()

  sr_i <- 0L
  pd_i <- 0L
  qc_i <- 0L

  for (pp in seq_len(nrow(eligible_patients))) {
    pid <- eligible_patients$patient_id[[pp]]
    n <- eligible_patients$n[[pp]]

    sm <- sample_map[
      patient_id == pid,
      .(
        sample_id,
        sector_num = as.integer(sector_num)
      )
    ]

    if (nrow(sm) != n) {
      stop(
        cohort_name, " / ", pid,
        ": sample count does not equal validated n."
      )
    }
    if (anyNA(sm$sector_num) || anyDuplicated(sm$sector_num) > 0L) {
      stop(
        cohort_name, " / ", pid,
        ": eligible patient has missing/duplicate sector coordinates."
      )
    }

    setorder(sm, sector_num, sample_id)

    if (require_consecutive) {
      expected_positions <- seq.int(min(sm$sector_num), length.out = n)
      if (!identical(as.integer(sm$sector_num),
                     as.integer(expected_positions))) {
        stop(
          cohort_name, " / ", pid,
          ": primary-eligible patient does not have consecutive T coordinates."
        )
      }
    }

    ordered_samples <- sm$sample_id
    ordered_t <- sm$sector_num

    if (anyDuplicated(ordered_samples) > 0L) {
      stop(cohort_name, " / ", pid, ": duplicate sample ID.")
    }

    max_k <- min(max(K_VALUES), n - 1L)
    patient_k_values <- K_VALUES[K_VALUES <= max_k]

    if (length(patient_k_values) == 0L) next

    for (kval in patient_k_values) {
      all_rows <- subset_results[
        .(pid, as.integer(kval)),
        nomatch = 0L
      ]

      if (nrow(all_rows) != choose(n, kval)) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": exhaustive subset count mismatch."
        )
      }

      # -------- unconstrained --------
      pd_i <- pd_i + 1L
      patient_design_rows[[pd_i]] <- summarise_design_rows(
        selected_rows = all_rows,
        patient_id = pid,
        n = n,
        k = kval,
        design = "unconstrained",
        spatial_cohort = cohort_name
      )

      # -------- grid-based dispersed --------
      grid_ranks <- grid_subsets_rank(n, kval)
      grid_selection <- rank_matrix_to_selection_table(
        rank_matrix = grid_ranks,
        ordered_samples = ordered_samples,
        ordered_t = ordered_t,
        design = "grid_dispersed"
      )

      if (anyDuplicated(grid_selection$subset_signature) > 0L) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": duplicated tied grid selection."
        )
      }

      grid_idx <- match(
        grid_selection$subset_signature,
        all_rows$subset_signature
      )
      if (anyNA(grid_idx)) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": at least one grid subset was not found in exhaustive universe."
        )
      }

      grid_rows <- all_rows[grid_idx]

      # Grid k>=2 must span first and last ordered sector.
      grid_rank_lists <- strsplit(
        grid_selection$rank_indices,
        ",",
        fixed = TRUE
      )
      grid_endpoint_ok <- all(vapply(
        grid_rank_lists,
        function(z) {
          rr <- as.integer(z)
          min(rr) == 1L && max(rr) == n
        },
        logical(1)
      ))

      if (!grid_endpoint_ok) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": grid selection does not span both endpoints."
        )
      }

      pd_i <- pd_i + 1L
      patient_design_rows[[pd_i]] <- summarise_design_rows(
        selected_rows = grid_rows,
        patient_id = pid,
        n = n,
        k = kval,
        design = "grid_dispersed",
        spatial_cohort = cohort_name
      )

      # Store exact selected grid subsets + genomic outcomes.
      sr_i <- sr_i + 1L
      gz <- cbind(
        grid_selection,
        grid_rows[, ..METRIC_COLS]
      )
      gz[, spatial_cohort := cohort_name]
      gz[, patient_id := pid]
      gz[, n := n]
      gz[, k := kval]
      setcolorder(
        gz,
        c(
          "spatial_cohort", "patient_id", "n", "k",
          "design", "design_selection_id",
          "rank_indices", "selected_t_values",
          "selected_samples_spatial_order",
          "subset_signature",
          "rank_span", "t_label_span",
          METRIC_COLS
        )
      )
      selected_subset_rows[[sr_i]] <- gz

      # -------- clustered / adjacent --------
      cluster_ranks <- clustered_subsets_rank(n, kval)
      cluster_selection <- rank_matrix_to_selection_table(
        rank_matrix = cluster_ranks,
        ordered_samples = ordered_samples,
        ordered_t = ordered_t,
        design = "clustered_adjacent"
      )

      expected_windows <- n - kval + 1L
      if (nrow(cluster_selection) != expected_windows) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": clustered window count mismatch."
        )
      }

      if (anyDuplicated(cluster_selection$subset_signature) > 0L) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": duplicate clustered selection."
        )
      }

      # Consecutive in ORDER/RANK, by definition.
      cluster_rank_lists <- strsplit(
        cluster_selection$rank_indices,
        ",",
        fixed = TRUE
      )
      cluster_contiguous_ok <- all(vapply(
        cluster_rank_lists,
        function(z) {
          rr <- as.integer(z)
          identical(rr, seq.int(min(rr), length.out = length(rr)))
        },
        logical(1)
      ))

      if (!cluster_contiguous_ok) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": clustered selector produced a non-contiguous rank window."
        )
      }

      cluster_idx <- match(
        cluster_selection$subset_signature,
        all_rows$subset_signature
      )
      if (anyNA(cluster_idx)) {
        stop(
          cohort_name, " / ", pid, " / k=", kval,
          ": at least one clustered subset was not found in exhaustive universe."
        )
      }

      cluster_rows <- all_rows[cluster_idx]

      pd_i <- pd_i + 1L
      patient_design_rows[[pd_i]] <- summarise_design_rows(
        selected_rows = cluster_rows,
        patient_id = pid,
        n = n,
        k = kval,
        design = "clustered_adjacent",
        spatial_cohort = cohort_name
      )

      sr_i <- sr_i + 1L
      cz <- cbind(
        cluster_selection,
        cluster_rows[, ..METRIC_COLS]
      )
      cz[, spatial_cohort := cohort_name]
      cz[, patient_id := pid]
      cz[, n := n]
      cz[, k := kval]
      setcolorder(
        cz,
        c(
          "spatial_cohort", "patient_id", "n", "k",
          "design", "design_selection_id",
          "rank_indices", "selected_t_values",
          "selected_samples_spatial_order",
          "subset_signature",
          "rank_span", "t_label_span",
          METRIC_COLS
        )
      )
      selected_subset_rows[[sr_i]] <- cz

      # -------- selection-level QC row --------
      qc_i <- qc_i + 1L
      selection_qc_rows[[qc_i]] <- data.table(
        spatial_cohort = cohort_name,
        patient_id = pid,
        n = n,
        k = kval,
        n_all_subsets = nrow(all_rows),
        expected_all_subsets = choose(n, kval),
        n_grid_solutions = nrow(grid_selection),
        n_cluster_windows = nrow(cluster_selection),
        expected_cluster_windows = expected_windows,
        grid_endpoint_span_ok = grid_endpoint_ok,
        cluster_rank_contiguous_ok = cluster_contiguous_ok,
        primary_positions_consecutive =
          if (require_consecutive) TRUE else NA
      )
    }
  }

  selected_table <- rbindlist(
    selected_subset_rows,
    use.names = TRUE,
    fill = TRUE
  )
  design_summary <- rbindlist(
    patient_design_rows,
    use.names = TRUE,
    fill = TRUE
  )
  selection_qc <- rbindlist(
    selection_qc_rows,
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    selected_table,
    patient_id, k, design, design_selection_id
  )
  setorder(
    design_summary,
    patient_id, k, design
  )
  setorder(
    selection_qc,
    patient_id, k
  )

  # Exactly 3 design summaries per patient x k.
  patient_k_counts <- unique(
    design_summary[, .(patient_id, n, k)]
  )
  design_counts <- design_summary[, .N, by = .(patient_id, k)]
  if (any(design_counts$N != 3L)) {
    stop(
      cohort_name,
      ": not every patient x k has exactly three design summaries."
    )
  }

  observed_n_by_k <- patient_k_counts[
    ,
    .(n_patients = uniqueN(patient_id)),
    by = k
  ]
  setorder(observed_n_by_k, k)

  expected_table <- data.table(
    k = as.integer(names(expected_n_by_k)),
    expected_n_patients = as.integer(expected_n_by_k)
  )
  observed_n_by_k <- merge(
    expected_table,
    observed_n_by_k,
    by = "k",
    all = TRUE,
    sort = TRUE
  )
  observed_n_by_k[
    ,
    pass := expected_n_patients == n_patients
  ]

  if (!all(observed_n_by_k$pass)) {
    if (DEBUG) {
      fwrite(
        observed_n_by_k,
        file.path(
          PATHS$qc_dir,
          paste0("06_", cohort_name, "_N_by_k_failure.tsv")
        ),
        sep = "\t"
      )
    }
    stop(cohort_name, ": N-by-k QC FAILED.")
  }

  if (!all(selection_qc$n_all_subsets ==
           selection_qc$expected_all_subsets)) {
    stop(cohort_name, ": exhaustive subset-count QC FAILED.")
  }
  if (!all(selection_qc$n_cluster_windows ==
           selection_qc$expected_cluster_windows)) {
    stop(cohort_name, ": clustered-window-count QC FAILED.")
  }
  if (!all(selection_qc$grid_endpoint_span_ok)) {
    stop(cohort_name, ": grid endpoint-span QC FAILED.")
  }
  if (!all(selection_qc$cluster_rank_contiguous_ok)) {
    stop(cohort_name, ": cluster-contiguity QC FAILED.")
  }

  # In the strict primary cohort with consecutive coordinates:
  # grid rank span = n-1 and clustered rank span = k-1.
  if (require_consecutive) {
    strict_selected <- selected_table
    grid_bad <- strict_selected[
      design == "grid_dispersed" &
        rank_span != (n - 1L)
    ]
    cluster_bad <- strict_selected[
      design == "clustered_adjacent" &
        rank_span != (k - 1L)
    ]

    if (nrow(grid_bad) > 0L || nrow(cluster_bad) > 0L) {
      stop(
        cohort_name,
        ": strict-cohort spatial-span identity FAILED."
      )
    }
  }

  list(
    selected_table = selected_table,
    design_summary = design_summary,
    selection_qc = selection_qc,
    observed_n_by_k = observed_n_by_k
  )
}

# ------------------------------------------------------------------------------
# 2. Build PRIMARY and ORDER-ONLY sensitivity spatial datasets
# ------------------------------------------------------------------------------

primary <- construct_spatial_cohort(
  cohort_name = "primary_strict",
  eligibility_column = "spatial_primary_eligible",
  expected_n_by_k = EXPECTED_PRIMARY_N_BY_K,
  require_consecutive = TRUE
)

order_sensitivity <- construct_spatial_cohort(
  cohort_name = "order_only_sensitivity",
  eligibility_column = "spatial_order_eligible",
  expected_n_by_k = EXPECTED_ORDER_N_BY_K,
  require_consecutive = FALSE
)

# ------------------------------------------------------------------------------
# 3. Regression of spatial "unconstrained" rows against Script 04 summaries
# ------------------------------------------------------------------------------

all_design_summary <- rbindlist(
  list(
    primary$design_summary,
    order_sensitivity$design_summary
  ),
  use.names = TRUE
)

u <- all_design_summary[design == "unconstrained"]

# One row per spatial cohort x patient x k; compare every unconstrained mean to
# the canonical Script 04 patient-level expected performance.
regression <- merge(
  u,
  unconstrained_summary,
  by = c("patient_id", "n", "k"),
  all.x = TRUE,
  sort = TRUE,
  suffixes = c("__spatial", "__script04")
)

regression_metric_map <- data.table(
  spatial_col = paste0(METRIC_COLS, "_mean"),
  script04_col = c(
    "recall_nonubiquitous_detection_mean__script04",
    "recall_heterogeneity_classification_mean__script04",
    "apparent_ubiquity_error_mean__script04",
    "conditional_apparent_ubiquity_rate_mean__script04",
    "recall_private_mean__script04",
    "nonubiquitous_fraction_abs_error_mean__script04",
    "jaccard_ith_abs_error_mean__script04"
  )
)

# After merge, only overlapping column names receive suffixes. Every spatial
# metric mean overlaps its Script 04 counterpart, so normalize the spatial names.
for (cc in METRIC_COLS) {
  old_name <- paste0(cc, "_mean__spatial")
  new_name <- paste0(cc, "_mean")
  if (old_name %in% names(regression)) {
    setnames(regression, old_name, new_name)
  }
}

regression_errors <- regression_metric_map[
  ,
  .(
    metric = gsub("_mean$", "", spatial_col),
    max_abs_error = vapply(
      seq_len(.N),
      function(i) {
        max_abs_diff_na_safe(
          regression[[spatial_col[[i]]]],
          regression[[script04_col[[i]]]]
        )
      },
      numeric(1)
    )
  )
]
regression_errors[, pass := max_abs_error <= TOLERANCE]

if (!all(regression_errors$pass)) {
  if (DEBUG) {
    fwrite(
      regression,
      file.path(PATHS$qc_dir, "06_unconstrained_regression_failure.tsv"),
      sep = "\t"
    )
  }
  stop("Spatial unconstrained regression against Script 04 FAILED.")
}

# ------------------------------------------------------------------------------
# 4. Patient-level grid-vs-cluster deltas (no cohort inference yet)
# ------------------------------------------------------------------------------

make_patient_deltas <- function(design_summary, cohort_name) {
  g <- design_summary[design == "grid_dispersed"]
  c <- design_summary[design == "clustered_adjacent"]
  u <- design_summary[design == "unconstrained"]

  key_cols <- c("spatial_cohort", "patient_id", "n", "k")

  g_keep <- copy(g)
  c_keep <- copy(c)
  u_keep <- copy(u)

  # Retain means only for paired outcome contrasts.
  mean_cols <- paste0(METRIC_COLS, "_mean")

  g_keep <- g_keep[, c(key_cols, "n_design_subsets", mean_cols), with = FALSE]
  c_keep <- c_keep[, c(key_cols, "n_design_subsets", mean_cols), with = FALSE]
  u_keep <- u_keep[, c(key_cols, "n_design_subsets", mean_cols), with = FALSE]

  setnames(
    g_keep,
    c("n_design_subsets", mean_cols),
    paste0(c("n_design_subsets", mean_cols), "__grid")
  )
  setnames(
    c_keep,
    c("n_design_subsets", mean_cols),
    paste0(c("n_design_subsets", mean_cols), "__cluster")
  )
  setnames(
    u_keep,
    c("n_design_subsets", mean_cols),
    paste0(c("n_design_subsets", mean_cols), "__unconstrained")
  )

  z <- Reduce(
    function(a, b) merge(a, b, by = key_cols, all = TRUE, sort = TRUE),
    list(g_keep, c_keep, u_keep)
  )

  if (anyNA(z$patient_id) || nrow(z) != nrow(g)) {
    stop(cohort_name, ": paired-design merge FAILED.")
  }

  for (metric in HIGHER_BETTER_METRICS) {
    gcol <- paste0(metric, "_mean__grid")
    ccol <- paste0(metric, "_mean__cluster")
    ucol <- paste0(metric, "_mean__unconstrained")

    z[[paste0("grid_advantage__", metric)]] <-
      z[[gcol]] - z[[ccol]]
    z[[paste0("grid_minus_unconstrained__", metric)]] <-
      z[[gcol]] - z[[ucol]]
    z[[paste0("unconstrained_minus_cluster__", metric)]] <-
      z[[ucol]] - z[[ccol]]
  }

  for (metric in LOWER_BETTER_METRICS) {
    gcol <- paste0(metric, "_mean__grid")
    ccol <- paste0(metric, "_mean__cluster")
    ucol <- paste0(metric, "_mean__unconstrained")

    # Reverse error contrasts so positive = grid better.
    z[[paste0("grid_advantage__", metric)]] <-
      z[[ccol]] - z[[gcol]]
    z[[paste0("grid_minus_unconstrained__", metric)]] <-
      z[[ucol]] - z[[gcol]]
    z[[paste0("unconstrained_minus_cluster__", metric)]] <-
      z[[ccol]] - z[[ucol]]
  }

  z
}

primary_deltas <- make_patient_deltas(
  primary$design_summary,
  "primary_strict"
)

order_deltas <- make_patient_deltas(
  order_sensitivity$design_summary,
  "order_only_sensitivity"
)

# ------------------------------------------------------------------------------
# 5. Save canonical Script 06 outputs
# ------------------------------------------------------------------------------

fwrite(
  primary$selected_table,
  file.path(
    spatial_dir,
    "06_hcc_spatial_selected_subsets_primary_strict.tsv"
  ),
  sep = "\t"
)
if (DEBUG) {
  fwrite(
    primary$design_summary,
    file.path(
      spatial_dir,
      "06_hcc_spatial_patient_design_summary_primary_strict.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    primary_deltas,
    file.path(
      spatial_dir,
      "06_hcc_spatial_patient_deltas_primary_strict.tsv"
    ),
    sep = "\t"
  )
}

fwrite(
  order_sensitivity$selected_table,
  file.path(
    spatial_dir,
    "06_hcc_spatial_selected_subsets_order_sensitivity.tsv"
  ),
  sep = "\t"
)
if (DEBUG) {
  fwrite(
    order_sensitivity$design_summary,
    file.path(
      spatial_dir,
      "06_hcc_spatial_patient_design_summary_order_sensitivity.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    order_deltas,
    file.path(
      spatial_dir,
      "06_hcc_spatial_patient_deltas_order_sensitivity.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    rbindlist(
      list(primary$selection_qc, order_sensitivity$selection_qc),
      use.names = TRUE
    ),
    file.path(PATHS$qc_dir, "06_spatial_selection_qc.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    rbindlist(
      list(primary$observed_n_by_k, order_sensitivity$observed_n_by_k),
      idcol = "cohort_index"
    ),
    file.path(PATHS$qc_dir, "06_spatial_N_by_k_qc.tsv"),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    regression_errors,
    file.path(PATHS$qc_dir, "06_unconstrained_regression_summary.tsv"),
    sep = "\t"
  )
}

saveRDS(
  list(
    primary = list(
      design_summary = primary$design_summary,
      patient_deltas = primary_deltas
    ),
    order_only_sensitivity = list(
      design_summary = order_sensitivity$design_summary,
      patient_deltas = order_deltas
    )
  ),
  file.path(
    PATHS$intermediate_dir,
    "hcc_spatial_designs_all_filtered.rds"
  )
)

# ------------------------------------------------------------------------------
# 6. Human-readable QC summary
# ------------------------------------------------------------------------------

primary_n <- primary$observed_n_by_k
order_n <- order_sensitivity$observed_n_by_k

n_tied_grid_primary <- primary$selection_qc[n_grid_solutions > 1L, .N]
max_grid_ties_primary <- max(primary$selection_qc$n_grid_solutions)
max_regression_error <- max(regression_errors$max_abs_error)

fixed_primary_n_ge5 <- unique(
  primary_deltas[n >= 5L & k %in% 2:4, patient_id]
)

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  "Script 06: spatial design construction",
  "=====================================================",
  "",
  "Spatial design definitions:",
  "  grid_dispersed     = endpoint-spanning, approximately equally spaced ordered ranks",
  "  clustered_adjacent = all contiguous k-sector rank windows, averaged over window location",
  "  unconstrained      = all k-sector subsets, uniformly weighted",
  "",
  "No metric physical distance is inferred from T-label differences.",
  "",
  "Primary strict spatial cohort N by k (requires n>k):",
  paste0(
    "  k=", primary_n$k, ": N=", primary_n$n_patients,
    collapse = "\n"
  ),
  "",
  "Order-only spatial sensitivity N by k (requires n>k):",
  paste0(
    "  k=", order_n$k, ": N=", order_n$n_patients,
    collapse = "\n"
  ),
  "",
  sprintf(
    "Primary strict fixed n>=5 patient count for k=2..4 robustness: %d",
    length(fixed_primary_n_ge5)
  ),
  "",
  "Selection QC: PASS",
  sprintf(
    "  primary patient x k combinations with >1 tied grid solution: %d",
    n_tied_grid_primary
  ),
  sprintf(
    "  maximum number of tied grid solutions in primary cohort: %d",
    max_grid_ties_primary
  ),
  "  all grid selections span ordered endpoints: PASS",
  "  all clustered selections are contiguous ordered-rank windows: PASS",
  "  all selected subsets matched the validated exhaustive universe: PASS",
  "",
  "Unconstrained regression against Script 04: PASS",
  sprintf(
    "  maximum absolute metric error: %.3g",
    max_regression_error
  ),
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed),
  "",
  "No cohort-level grid-vs-cluster inference was performed in Script 06.",
  "Proceed to Script 07 only if every QC above is PASS."
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "06_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
