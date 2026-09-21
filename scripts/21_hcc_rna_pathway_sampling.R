# ==============================================================================
# 21_hcc_rna_pathway_sampling.R
#
# Purpose
#   Minimal, rigorous analysis of sampling depth for transcriptomic
#   heterogeneity using the 50 Hallmark GSVA scores already published for each
#   tumour RNA-seq sector in Supplementary Table mmc4.xlsx.
#
# Why use the published GSVA scores rather than recompute from raw counts?
#   mmc4.xlsx already contains the exact
#   50 Hallmark GSVA scores used in the published study. Recomputing them would
#   add version-, gene-set-, normalization-, and implementation-dependence
#   without improving the answer to that question.
#
# Primary fixed cohort
#   Patients with >=5 tumour RNA-seq sectors in the published supplement.
#   Evaluate k=1..4 only, so every point is true downsampling (k<n).
#
# Exhaustive subset analysis
#   Every possible subset of k RNA sectors is enumerated; no Monte-Carlo sector
#   sampling is used. For each patient and k, the estimand is the exact mean
#   across all possible subsets.
#
# Endpoints
#   1) landscape_spearman
#      Spearman correlation across the 50 Hallmark pathways between the mean
#      pathway-score vector of a subset and the mean vector from all available
#      RNA sectors of that patient. Higher = better preservation.
#
#   2) landscape_rmse_z
#      RMSE between subset and full-patient pathway centroids after each pathway
#      is z-standardized across all 406 published tumour RNA sectors.
#      Lower = better preservation. The z-standardization prevents pathways
#      with larger numerical GSVA ranges from dominating the distance.
#
#   3) pathway_ith_absolute_error
#      Absolute error in pathway-level ITH. Pathway ITH is the mean pairwise
#      RMS Euclidean distance between sector vectors in the globally
#      z-standardized 50-pathway space. This endpoint is undefined for k=1.
#      Lower = better recovery of transcriptomic heterogeneity magnitude.
#
# Statistical unit
#   Patient. Cohort summaries are medians of patient-level exact expectations.
#   95% CIs use the validated exact-bootstrap empirical nonparametric
#   bootstrap distribution of the median.
#
# Important guardrails
#   - Full available RNA-sector set is a reference, not whole-tumour truth.
#   - No arbitrary "optimal k" threshold is imposed for RNA.
#   - RNA and DNA have different available-sector cohorts; this script does not
#     perform an invalid unpaired cross-modality hypothesis test.
#
# Required
#   - 00_config.R
#   - validated the intermediate/ output directory/hcc_sample_map.rds
#   - mmc4.xlsx
#
# Optional environment variable
#   MMC4_XLSX_PATH=/full/path/mmc4.xlsx
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "Package 'readxl' is required.\n",
    "Install once with: install.packages('readxl')"
  )
}

# ------------------------------ locate script --------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)

if (length(file_arg) > 0L) {
  script_dir <- dirname(
    normalizePath(
      sub("^--file=", "", file_arg[[1L]]),
      winslash = "/",
      mustWork = TRUE
    )
  )
} else {
  script_dir <- getwd()
}

source(file.path(script_dir, "00_config.R"))
source(file.path(script_dir, "functions_bootstrap_exact.R"))

cat("Script 21: HCC RNA pathway sampling\n")
cat("=======================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------- constants -----------------------------------
CI_LEVEL <- 0.95
PRIMARY_K <- 1:4
PRIMARY_MIN_N <- 5L

EXPECTED_MMC4_MD5 <- "33c14fc9e52e32cca6ea8c9c771d8985"
DNA_SHEET <- "merged_DNA_sample_table"
RNA_SHEET <- "merged_RNA_sample_table"

EXPECTED_DNA_PATIENTS <- 123L
EXPECTED_DNA_SECTORS <- 490L

EXPECTED_RNA_PATIENTS <- 111L
EXPECTED_RNA_SECTORS <- 406L
EXPECTED_GSVA_PATHWAYS <- 50L
EXPECTED_FIXED_RNA_N_GE5 <- 39L

EXPECTED_RNA_SECTOR_DISTRIBUTION <- c(
  `1` = 2L,
  `2` = 24L,
  `3` = 25L,
  `4` = 21L,
  `5` = 38L,
  `7` = 1L
)

EVENT_SET_ID <- "published_50_hallmark_gsva"

# ------------------------------- helpers -------------------------------------
require_columns <- function(x, cols, object_name) {
  miss <- setdiff(cols, names(x))
  if (length(miss) > 0L) {
    stop(
      object_name, " is missing required column(s): ",
      paste(miss, collapse = ", ")
    )
  }
}

locate_existing_file <- function(env_var, candidate_paths, label) {
  env_path <- Sys.getenv(env_var, unset = "")

  candidates <- unique(c(
    if (nzchar(env_path)) path.expand(env_path) else character(),
    candidate_paths
  ))
  candidates <- candidates[nzchar(candidates)]

  hit <- candidates[file.exists(candidates)]

  if (length(hit) == 0L) {
    stop(
      "Could not locate ", label, ".\n",
      "Set environment variable ", env_var, " to the full path."
    )
  }

  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

canonicalize_zero <- function(x, tol = 1e-12) {
  nm <- names(x)
  x <- as.numeric(x)
  names(x) <- nm
  x[is.finite(x) & abs(x) <= tol] <- 0
  x
}

canonicalize_correlation <- function(x, tol = 1e-12) {
  nm <- names(x)
  x <- as.numeric(x)
  names(x) <- nm

  finite <- is.finite(x)
  if (any(x[finite] < -1 - tol | x[finite] > 1 + tol)) {
    stop("Correlation outside [-1,1] beyond numerical tolerance.")
  }

  x[finite & abs(x - 1) <= tol] <- 1
  x[finite & abs(x + 1) <= tol] <- -1
  x[finite & abs(x) <= tol] <- 0

  x
}

mean_pairwise_rms_distance <- function(z_matrix) {
  z_matrix <- as.matrix(z_matrix)

  if (nrow(z_matrix) < 2L) {
    return(NA_real_)
  }

  if (ncol(z_matrix) < 1L) {
    stop("No pathway columns supplied to pairwise-distance calculation.")
  }

  d <- stats::dist(z_matrix, method = "euclidean")
  out <- mean(as.numeric(d)) / sqrt(ncol(z_matrix))

  if (!is.finite(out) || out < -1e-12) {
    stop("Invalid pathway-ITH distance.")
  }

  if (abs(out) <= 1e-12) out <- 0
  out
}

# Exact nonparametric bootstrap functions: shared implementation, sourced
# near the top of this script (functions_bootstrap_exact.R).

# ------------------------------ output paths ---------------------------------
ext_dir <- file.path(
  PATHS$output_dir,
  "results",
  "extensions",
  "rna"
)

ext_qc_dir <- file.path(
  PATHS$output_dir,
  "results",
  "qc",
  "extensions"
)

dir.create(
  ext_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  ext_qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------ locate inputs --------------------------------
mmc4_path <- locate_existing_file(
  "MMC4_XLSX_PATH",
  c(
    file.path(PROJECT_ROOT, "mmc4.xlsx"),
    file.path(dirname(PROJECT_ROOT), "mmc4.xlsx"),
    file.path(script_dir, "mmc4.xlsx")
  ),
  "Supplementary Table mmc4.xlsx"
)

mmc4_md5 <- unname(tools::md5sum(mmc4_path))

if (!identical(mmc4_md5, EXPECTED_MMC4_MD5)) {
  stop(
    "mmc4.xlsx MD5 mismatch.\n",
    "Expected: ", EXPECTED_MMC4_MD5, "\n",
    "Observed: ", mmc4_md5, "\n",
    "File: ", mmc4_path
  )
}

sample_map_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_sample_map.rds"
)

if (!file.exists(sample_map_path)) {
  stop(
    "Missing validated WGS sample map: ",
    sample_map_path
  )
}

# --------------------------- validated DNA crosswalk --------------------------
sample_map <- as.data.table(
  readRDS(sample_map_path)
)

require_columns(
  sample_map,
  c(
    "patient_id",
    "sample_id",
    "n_sectors"
  ),
  "hcc_sample_map"
)

sample_map[, `:=`(
  patient_id = as.character(patient_id),
  sample_id = as.character(sample_id),
  n_sectors = as.integer(n_sectors)
)]

if (nrow(sample_map) != EXPECTED_DNA_SECTORS ||
    uniqueN(sample_map$patient_id) != EXPECTED_DNA_PATIENTS ||
    uniqueN(sample_map$sample_id) != EXPECTED_DNA_SECTORS ||
    anyDuplicated(sample_map$sample_id) > 0L) {
  stop(
    "Validated WGS sample-map regression failed."
  )
}

sheet_names <- readxl::excel_sheets(mmc4_path)

if (!(DNA_SHEET %in% sheet_names) ||
    !(RNA_SHEET %in% sheet_names)) {
  stop(
    "mmc4.xlsx is missing required DNA/RNA merged sheet(s)."
  )
}

cat("Reading mmc4 DNA/RNA supplementary tables ...\n")

dna <- as.data.table(
  readxl::read_excel(
    mmc4_path,
    sheet = DNA_SHEET,
    .name_repair = "minimal"
  )
)

rna <- as.data.table(
  readxl::read_excel(
    mmc4_path,
    sheet = RNA_SHEET,
    .name_repair = "minimal"
  )
)

require_columns(
  dna,
  c(
    "Publication_ID",
    "DNA_lib",
    "Tumor",
    "Sample_type"
  ),
  DNA_SHEET
)

require_columns(
  rna,
  c(
    "Publication_ID",
    "RNA_lib",
    "Tumor",
    "Sample_type"
  ),
  RNA_SHEET
)

dna[, `:=`(
  Publication_ID = as.character(Publication_ID),
  DNA_lib = as.character(DNA_lib),
  Tumor = as.character(Tumor),
  Sample_type = as.character(Sample_type)
)]

rna[, `:=`(
  Publication_ID = as.character(Publication_ID),
  RNA_lib = as.character(RNA_lib),
  Tumor = as.character(Tumor),
  Sample_type = as.character(Sample_type)
)]

if (nrow(dna) != EXPECTED_DNA_SECTORS ||
    uniqueN(dna$Publication_ID) != EXPECTED_DNA_PATIENTS ||
    uniqueN(dna$DNA_lib) != EXPECTED_DNA_SECTORS ||
    anyDuplicated(dna$DNA_lib) > 0L) {
  stop(
    "mmc4 merged_DNA_sample_table regression failed."
  )
}

dna_crosswalk <- merge(
  dna[, .(
    Publication_ID,
    DNA_lib
  )],
  sample_map[, .(
    sample_id,
    patient_id
  )],
  by.x = "DNA_lib",
  by.y = "sample_id",
  all = TRUE,
  sort = FALSE
)

if (nrow(dna_crosswalk) != EXPECTED_DNA_SECTORS ||
    anyNA(dna_crosswalk$Publication_ID) ||
    anyNA(dna_crosswalk$patient_id)) {
  stop(
    "mmc4 Publication_ID -> validated patient crosswalk failed."
  )
}

publication_patient_map <- unique(
  dna_crosswalk[, .(
    Publication_ID,
    patient_id
  )]
)

if (nrow(publication_patient_map) != EXPECTED_DNA_PATIENTS ||
    anyDuplicated(publication_patient_map$Publication_ID) > 0L ||
    anyDuplicated(publication_patient_map$patient_id) > 0L) {
  stop(
    "Publication_ID <-> validated patient_id is not one-to-one."
  )
}

# ------------------------------- RNA QC --------------------------------------
sample_type_values <- sort(unique(
  trimws(rna$Sample_type)
))

if (!identical(
  sample_type_values,
  "Tumor"
)) {
  stop(
    "merged_RNA_sample_table contains unexpected Sample_type values: ",
    paste(sample_type_values, collapse = ", ")
  )
}

if (nrow(rna) != EXPECTED_RNA_SECTORS ||
    uniqueN(rna$Publication_ID) != EXPECTED_RNA_PATIENTS ||
    uniqueN(rna$RNA_lib) != EXPECTED_RNA_SECTORS ||
    anyDuplicated(rna$RNA_lib) > 0L ||
    anyNA(rna$Publication_ID) ||
    anyNA(rna$RNA_lib)) {
  stop(
    "mmc4 merged_RNA_sample_table regression failed: expected ",
    EXPECTED_RNA_PATIENTS, " patients / ",
    EXPECTED_RNA_SECTORS, " tumour RNA libraries."
  )
}

gsva_cols <- grep(
  "^GSVA_HALLMARK_",
  names(rna),
  value = TRUE
)

if (length(gsva_cols) != EXPECTED_GSVA_PATHWAYS ||
    anyDuplicated(gsva_cols) > 0L) {
  stop(
    "Expected exactly ",
    EXPECTED_GSVA_PATHWAYS,
    " unique GSVA_HALLMARK columns; observed ",
    length(gsva_cols), "."
  )
}

# readxl infers a whole Excel column as character if even one numeric cell was
# stored as text. The published workbook contains this formatting issue in a
# small number of Hallmark columns. Coerce only when every non-missing cell is
# losslessly parseable as a finite number; otherwise stop.
gsva_coercion_qc <- list()
coercion_i <- 0L

for (cc in gsva_cols) {
  original <- rna[[cc]]
  original_class <- class(original)[1L]

  if (is.numeric(original)) {
    numeric_values <- as.numeric(original)
    n_non_numeric_after_parse <- 0L
  } else {
    original_chr <- trimws(as.character(original))

    numeric_values <- suppressWarnings(
      as.numeric(original_chr)
    )

    bad_parse <- !is.na(original_chr) &
      nzchar(original_chr) &
      is.na(numeric_values)

    n_non_numeric_after_parse <- sum(bad_parse)

    if (n_non_numeric_after_parse > 0L) {
      stop(
        "GSVA pathway column contains genuinely non-numeric value(s): ",
        cc, ". Examples: ",
        paste(
          head(unique(original_chr[bad_parse]), 10L),
          collapse = ", "
        )
      )
    }

    rna[[cc]] <- numeric_values
  }

  if (anyNA(numeric_values) ||
      any(!is.finite(numeric_values))) {
    stop(
      "GSVA pathway column contains missing/non-finite value(s) after safe numeric coercion: ",
      cc
    )
  }

  coercion_i <- coercion_i + 1L
  gsva_coercion_qc[[coercion_i]] <- data.table(
    pathway_column = cc,
    original_R_class = original_class,
    coerced_to_numeric = !is.numeric(original),
    n_values = length(numeric_values),
    n_missing_after_parse = sum(is.na(numeric_values)),
    n_non_numeric_after_parse = n_non_numeric_after_parse,
    min_value = min(numeric_values),
    max_value = max(numeric_values)
  )
}

gsva_coercion_qc <- rbindlist(
  gsva_coercion_qc,
  use.names = TRUE
)

if (DEBUG) {
  fwrite(
    gsva_coercion_qc,
    file.path(
      ext_qc_dir,
      "21_rna_gsva_numeric_coercion_qc.tsv"
    ),
    sep = "\t"
  )
}

n_coerced_gsva_columns <- gsva_coercion_qc[
  coerced_to_numeric == TRUE,
  .N
]

rna <- merge(
  rna,
  publication_patient_map,
  by = "Publication_ID",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(rna$patient_id) ||
    uniqueN(rna$patient_id) != EXPECTED_RNA_PATIENTS) {
  stop(
    "At least one RNA Publication_ID failed DNA/validated-patient crosswalk."
  )
}

rna_patient_n <- rna[, .(
  n_rna_sectors = .N
), by = .(
  Publication_ID,
  patient_id
)]

sector_dist <- table(
  factor(
    rna_patient_n$n_rna_sectors,
    levels = as.integer(
      names(EXPECTED_RNA_SECTOR_DISTRIBUTION)
    )
  )
)

sector_dist_named <- as.integer(sector_dist)
names(sector_dist_named) <- names(
  EXPECTED_RNA_SECTOR_DISTRIBUTION
)

if (!identical(
  unname(sector_dist_named),
  unname(as.integer(EXPECTED_RNA_SECTOR_DISTRIBUTION))
)) {
  stop(
    "RNA sector-count distribution regression failed.\nExpected: ",
    paste(
      names(EXPECTED_RNA_SECTOR_DISTRIBUTION),
      EXPECTED_RNA_SECTOR_DISTRIBUTION,
      sep = ":",
      collapse = ", "
    ),
    "\nObserved: ",
    paste(
      names(sector_dist_named),
      sector_dist_named,
      sep = ":",
      collapse = ", "
    )
  )
}

fixed_patients <- rna_patient_n[
  n_rna_sectors >= PRIMARY_MIN_N,
  .(
    Publication_ID,
    patient_id,
    n = as.integer(n_rna_sectors)
  )
]

setorder(
  fixed_patients,
  patient_id
)

if (nrow(fixed_patients) != EXPECTED_FIXED_RNA_N_GE5) {
  stop(
    "RNA fixed n>=5 cohort regression failed: expected ",
    EXPECTED_FIXED_RNA_N_GE5,
    ", observed ",
    nrow(fixed_patients), "."
  )
}

if (any(fixed_patients$n <= max(PRIMARY_K))) {
  stop(
    "RNA primary fixed cohort contains a k=n point."
  )
}

# Freeze RNA/sample QC outputs.
if (DEBUG) {
  fwrite(
    rna_patient_n,
    file.path(
      ext_qc_dir,
      "21_rna_patient_sector_counts.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    data.table(
      n_rna_sectors = as.integer(
        names(EXPECTED_RNA_SECTOR_DISTRIBUTION)
      ),
      n_patients = as.integer(
        EXPECTED_RNA_SECTOR_DISTRIBUTION
      )
    ),
    file.path(
      ext_qc_dir,
      "21_rna_sector_count_distribution.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    data.table(
      pathway_column = gsva_cols
    ),
    file.path(
      ext_qc_dir,
      "21_rna_hallmark_pathway_columns.tsv"
    ),
    sep = "\t"
  )
}

# ------------------------- pathway standardization ----------------------------
cat("Standardizing the 50 published Hallmark GSVA pathways ...\n")

gsva_raw <- as.matrix(
  rna[, ..gsva_cols]
)

storage.mode(gsva_raw) <- "double"

pathway_mean <- colMeans(gsva_raw)
pathway_sd <- apply(
  gsva_raw,
  2L,
  stats::sd
)

if (any(!is.finite(pathway_mean)) ||
    any(!is.finite(pathway_sd)) ||
    any(pathway_sd <= 0)) {
  stop(
    "At least one Hallmark pathway has invalid/zero global SD."
  )
}

gsva_z <- sweep(
  gsva_raw,
  2L,
  pathway_mean,
  FUN = "-"
)

gsva_z <- sweep(
  gsva_z,
  2L,
  pathway_sd,
  FUN = "/"
)

if (any(!is.finite(gsva_z))) {
  stop(
    "Non-finite value after GSVA pathway z-standardization."
  )
}

scaling_qc <- data.table(
  pathway_column = gsva_cols,
  global_mean_406_tumour_sectors =
    as.numeric(pathway_mean),
  global_sd_406_tumour_sectors =
    as.numeric(pathway_sd)
)

if (DEBUG) {
  fwrite(
    scaling_qc,
    file.path(
      ext_qc_dir,
      "21_rna_hallmark_global_scaling.tsv"
    ),
    sep = "\t"
  )
}

# Keep row indices frozen after merge/order.
rna[, matrix_row := .I]

# ----------------------- exhaustive subset analysis ---------------------------
cat("Enumerating all RNA-sector subsets for fixed n>=5 cohort ...\n")

subset_rows <- list()
patient_rows <- list()

subset_i <- 0L
patient_i <- 0L

for (pp in seq_len(nrow(fixed_patients))) {
  pid <- fixed_patients$patient_id[[pp]]
  pubid <- fixed_patients$Publication_ID[[pp]]
  n <- fixed_patients$n[[pp]]

  rr <- rna[
    patient_id == pid
  ]

  if (nrow(rr) != n) {
    stop(
      "RNA patient sector-count mismatch for patient ",
      pid, "."
    )
  }

  row_idx <- rr$matrix_row
  raw_patient <- gsva_raw[row_idx, , drop = FALSE]
  z_patient <- gsva_z[row_idx, , drop = FALSE]

  full_centroid_raw <- colMeans(raw_patient)
  full_centroid_z <- colMeans(z_patient)
  full_ith <- mean_pairwise_rms_distance(
    z_patient
  )

  if (!is.finite(full_ith)) {
    stop(
      "Full RNA pathway ITH is not finite for fixed-cohort patient ",
      pid, "."
    )
  }

  if (stats::sd(full_centroid_raw) <= 0) {
    stop(
      "Full RNA centroid has zero across-pathway variance for patient ",
      pid, "."
    )
  }

  for (kval in PRIMARY_K) {
    if (kval >= n) {
      stop(
        "True-downsampling invariant violated for RNA patient ",
        pid, ": k=", kval, ", n=", n
      )
    }

    combos <- utils::combn(
      seq_len(n),
      kval,
      simplify = TRUE
    )

    if (kval == 1L) {
      combos <- matrix(
        combos,
        nrow = 1L
      )
    }

    n_subsets <- ncol(combos)

    metric_rho <- numeric(n_subsets)
    metric_rmse <- numeric(n_subsets)
    metric_ith_error <- rep(
      NA_real_,
      n_subsets
    )

    for (jj in seq_len(n_subsets)) {
      idx <- combos[, jj]

      subset_centroid_raw <- colMeans(
        raw_patient[idx, , drop = FALSE]
      )

      subset_centroid_z <- colMeans(
        z_patient[idx, , drop = FALSE]
      )

      rho <- suppressWarnings(
        stats::cor(
          subset_centroid_raw,
          full_centroid_raw,
          method = "spearman"
        )
      )

      if (!is.finite(rho)) {
        stop(
          "Non-finite landscape Spearman correlation for patient ",
          pid, ", k=", kval, ", subset=", jj, "."
        )
      }

      rmse <- sqrt(
        mean(
          (
            subset_centroid_z -
              full_centroid_z
          )^2
        )
      )

      if (!is.finite(rmse) ||
          rmse < -1e-12) {
        stop(
          "Invalid standardized landscape RMSE for patient ",
          pid, ", k=", kval, "."
        )
      }

      if (abs(rmse) <= 1e-12) {
        rmse <- 0
      }

      ith_error <- NA_real_

      if (kval >= 2L) {
        subset_ith <- mean_pairwise_rms_distance(
          z_patient[idx, , drop = FALSE]
        )

        ith_error <- abs(
          subset_ith -
            full_ith
        )

        if (!is.finite(ith_error) ||
            ith_error < -1e-12) {
          stop(
            "Invalid pathway-ITH error for patient ",
            pid, ", k=", kval, "."
          )
        }

        if (abs(ith_error) <= 1e-12) {
          ith_error <- 0
        }
      }

      metric_rho[[jj]] <- rho
      metric_rmse[[jj]] <- rmse
      metric_ith_error[[jj]] <- ith_error

      subset_i <- subset_i + 1L

      subset_rows[[subset_i]] <- data.table(
        event_set = EVENT_SET_ID,
        Publication_ID = pubid,
        patient_id = pid,
        n = n,
        k = as.integer(kval),
        subset_index = as.integer(jj),
        subset_rna_libraries = paste(
          rr$RNA_lib[idx],
          collapse = ";"
        ),
        landscape_spearman =
          canonicalize_correlation(rho),
        landscape_rmse_z =
          canonicalize_zero(rmse),
        full_pathway_ith_rms_z =
          full_ith,
        pathway_ith_absolute_error =
          if (is.finite(ith_error)) {
            canonicalize_zero(ith_error)
          } else {
            NA_real_
          }
      )
    }

    patient_i <- patient_i + 1L

    patient_rows[[patient_i]] <- data.table(
      event_set = EVENT_SET_ID,
      Publication_ID = pubid,
      patient_id = pid,
      n = n,
      k = as.integer(kval),
      n_subsets = as.integer(n_subsets),
      full_pathway_ith_rms_z =
        full_ith,
      expected_landscape_spearman =
        mean(metric_rho),
      expected_landscape_rmse_z =
        mean(metric_rmse),
      expected_pathway_ith_absolute_error =
        if (kval >= 2L) {
          mean(metric_ith_error)
        } else {
          NA_real_
        }
    )
  }
}

subset_results <- rbindlist(
  subset_rows,
  use.names = TRUE,
  fill = TRUE
)

patient_expected <- rbindlist(
  patient_rows,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  subset_results,
  patient_id,
  k,
  subset_index
)

setorder(
  patient_expected,
  patient_id,
  k
)

patient_expected[
  ,
  expected_landscape_spearman :=
    canonicalize_correlation(
      expected_landscape_spearman
    )
]

patient_expected[
  ,
  expected_landscape_rmse_z :=
    canonicalize_zero(
      expected_landscape_rmse_z
    )
]

patient_expected[
  ,
  expected_pathway_ith_absolute_error :=
    canonicalize_zero(
      expected_pathway_ith_absolute_error
    )
]

if (nrow(patient_expected) !=
    EXPECTED_FIXED_RNA_N_GE5 *
      length(PRIMARY_K) ||
    uniqueN(patient_expected$patient_id) !=
      EXPECTED_FIXED_RNA_N_GE5 ||
    anyDuplicated(
      patient_expected[, .(
        patient_id,
        k
      )]
    ) > 0L) {
  stop(
    "RNA patient-expectation table has unexpected dimensions."
  )
}

if (any(
  patient_expected[
    k == 1L,
    is.finite(
      expected_pathway_ith_absolute_error
    )
  ]
)) {
  stop(
    "RNA pathway-ITH error must be undefined at k=1."
  )
}

# Monotonicity is NOT mathematically guaranteed patient-by-patient for these
# continuous metrics, so do not impose it as a false QC criterion.

if (DEBUG) {
  fwrite(
    subset_results,
    file.path(
      ext_dir,
      "21_rna_subset_metrics_fixed_n_ge5.tsv"
    ),
    sep = "\t"
  )
}

if (DEBUG) {
  fwrite(
    patient_expected,
    file.path(
      ext_dir,
      "21_rna_patient_expected_metrics_fixed_n_ge5.tsv"
    ),
    sep = "\t"
  )
}

# ----------------------------- cohort summaries -------------------------------
metric_specs <- data.table(
  metric = c(
    "landscape_spearman",
    "landscape_rmse_z",
    "pathway_ith_absolute_error"
  ),
  patient_column = c(
    "expected_landscape_spearman",
    "expected_landscape_rmse_z",
    "expected_pathway_ith_absolute_error"
  ),
  direction = c(
    "higher_is_better",
    "lower_is_better",
    "lower_is_better"
  ),
  defined_at_k1 = c(
    TRUE,
    TRUE,
    FALSE
  )
)

summary_rows <- list()
ss <- 0L

for (mm in seq_len(nrow(metric_specs))) {
  metric_id <- metric_specs$metric[[mm]]
  value_col <- metric_specs$patient_column[[mm]]
  direction <- metric_specs$direction[[mm]]

  for (kval in PRIMARY_K) {
    vals <- patient_expected[
      k == kval,
      get(value_col)
    ]

    vals <- as.numeric(
      vals[is.finite(vals)]
    )

    ci <- exact_bootstrap_median_ci(
      vals,
      level = CI_LEVEL,
      min_n = 1L
    )

    ss <- ss + 1L

    summary_rows[[ss]] <- data.table(
      analysis =
        "rna_fixed_n_ge5_true_downsampling",
      event_set = EVENT_SET_ID,
      metric = metric_id,
      direction = direction,
      k = as.integer(kval),
      n_fixed_cohort =
        EXPECTED_FIXED_RNA_N_GE5,
      n_patients_with_defined_endpoint =
        length(vals),
      fraction_fixed_cohort_with_defined_endpoint =
        length(vals) /
          EXPECTED_FIXED_RNA_N_GE5,
      median =
        safe_median(vals),
      ci_lower =
        unname(ci[["lower"]]),
      ci_upper =
        unname(ci[["upper"]]),
      ci_level =
        CI_LEVEL,
      bootstrap_method =
        if (length(vals) > 0L) {
          "exact_nonparametric_percentile"
        } else {
          NA_character_
        }
    )
  }
}

depth_summary <- rbindlist(
  summary_rows,
  use.names = TRUE
)

setorder(
  depth_summary,
  metric,
  k
)

# Explicit regression: no bootstrap label if the endpoint is undefined.
if (nrow(
  depth_summary[
    metric == "pathway_ith_absolute_error" &
      k == 1L &
      (
        n_patients_with_defined_endpoint != 0L |
          !is.na(bootstrap_method)
      )
  ]
) > 0L) {
  stop(
    "Undefined RNA ITH endpoint at k=1 has inconsistent metadata."
  )
}

fwrite(
  depth_summary,
  file.path(
    ext_dir,
    "21_rna_depth_summary_fixed_n_ge5.tsv"
  ),
  sep = "\t"
)

# ---------------------------- full-cohort ITH QC ------------------------------
full_ith_patient <- unique(
  patient_expected[, .(
    Publication_ID,
    patient_id,
    n,
    full_pathway_ith_rms_z
  )]
)

if (nrow(full_ith_patient) != EXPECTED_FIXED_RNA_N_GE5) {
  stop(
    "Unexpected full pathway-ITH patient table size."
  )
}

if (DEBUG) {
  fwrite(
    full_ith_patient,
    file.path(
      ext_dir,
      "21_rna_full_pathway_ith_fixed_n_ge5.tsv"
    ),
    sep = "\t"
  )
}

# ------------------------------- provenance ----------------------------------
provenance <- data.table(
  item = c(
    "mmc4_path",
    "mmc4_md5",
    "rna_sheet",
    "rna_patients",
    "rna_tumour_sectors",
    "hallmark_gsva_pathways",
    "hallmark_columns_safely_coerced_from_character",
    "fixed_n_ge5_rna_patients",
    "fixed_k_min",
    "fixed_k_max",
    "standardization_reference",
    "landscape_primary_similarity",
    "landscape_magnitude_error",
    "pathway_ith_definition",
    "raw_feature_counts_recomputed",
    "monte_carlo_sector_sampling"
  ),
  value = c(
    mmc4_path,
    mmc4_md5,
    RNA_SHEET,
    as.character(
      uniqueN(rna$patient_id)
    ),
    as.character(
      nrow(rna)
    ),
    as.character(
      length(gsva_cols)
    ),
    as.character(
      n_coerced_gsva_columns
    ),
    as.character(
      nrow(fixed_patients)
    ),
    as.character(
      min(PRIMARY_K)
    ),
    as.character(
      max(PRIMARY_K)
    ),
    "per-pathway mean/SD across all 406 published tumour RNA sectors",
    "Spearman(subset pathway centroid, full available-sector pathway centroid)",
    "RMSE in globally z-standardized 50-pathway centroid space",
    "mean pairwise RMS Euclidean distance in globally z-standardized 50-pathway space; endpoint=absolute subset-vs-full error",
    "FALSE",
    "FALSE"
  )
)

if (DEBUG) {
  fwrite(
    provenance,
    file.path(
      ext_qc_dir,
      "21_rna_sampling_provenance.tsv"
    ),
    sep = "\t"
  )
}

# -------------------------- human-readable summary ----------------------------
get_summary <- function(metric_id, kval) {
  x <- depth_summary[
    metric == metric_id &
      k == kval
  ]

  if (nrow(x) != 1L) {
    stop(
      "Expected exactly one RNA summary row for ",
      metric_id, " k=", kval, "."
    )
  }

  x
}

rho_k1 <- get_summary(
  "landscape_spearman",
  1L
)
rho_k4 <- get_summary(
  "landscape_spearman",
  4L
)

rmse_k1 <- get_summary(
  "landscape_rmse_z",
  1L
)
rmse_k4 <- get_summary(
  "landscape_rmse_z",
  4L
)

ith_k2 <- get_summary(
  "pathway_ith_absolute_error",
  2L
)
ith_k4 <- get_summary(
  "pathway_ith_absolute_error",
  4L
)

summary_lines <- c(
  "HCC extensions - Script 21 RNA pathway sampling",
  "=======================================================",
  "",
  paste0(
    "mmc4.xlsx MD5: ",
    mmc4_md5,
    " [PASS]"
  ),
  sprintf(
    "Published RNA supplementary table: %d tumour sectors / %d patients [PASS]",
    nrow(rna),
    uniqueN(rna$patient_id)
  ),
  sprintf(
    "Published Hallmark GSVA columns: %d; %d Excel-text-formatted column(s) safely coerced; all values numeric/complete [PASS]",
    length(gsva_cols),
    n_coerced_gsva_columns
  ),
  paste0(
    "RNA sector-count distribution: ",
    paste(
      names(EXPECTED_RNA_SECTOR_DISTRIBUTION),
      EXPECTED_RNA_SECTOR_DISTRIBUTION,
      sep = ":",
      collapse = ", "
    ),
    " [PASS]"
  ),
  sprintf(
    "Primary RNA fixed n>=5 cohort: %d patients; k=1..4 are all true downsampling",
    nrow(fixed_patients)
  ),
  "",
  "Pathway-landscape preservation:",
  sprintf(
    "  k=1 median Spearman = %.6f [%.6f, %.6f]",
    rho_k1$median,
    rho_k1$ci_lower,
    rho_k1$ci_upper
  ),
  sprintf(
    "  k=4 median Spearman = %.6f [%.6f, %.6f]",
    rho_k4$median,
    rho_k4$ci_lower,
    rho_k4$ci_upper
  ),
  sprintf(
    "  k=1 median standardized centroid RMSE = %.6f [%.6f, %.6f]",
    rmse_k1$median,
    rmse_k1$ci_lower,
    rmse_k1$ci_upper
  ),
  sprintf(
    "  k=4 median standardized centroid RMSE = %.6f [%.6f, %.6f]",
    rmse_k4$median,
    rmse_k4$ci_lower,
    rmse_k4$ci_upper
  ),
  "",
  "Pathway-level ITH recovery:",
  "  k=1: undefined by construction (one sector has no within-patient pairwise heterogeneity)",
  sprintf(
    "  k=2 median absolute ITH error = %.6f [%.6f, %.6f]",
    ith_k2$median,
    ith_k2$ci_lower,
    ith_k2$ci_upper
  ),
  sprintf(
    "  k=4 median absolute ITH error = %.6f [%.6f, %.6f]",
    ith_k4$median,
    ith_k4$ci_lower,
    ith_k4$ci_upper
  ),
  "",
  paste0(
    "Guardrail: RNA uses the 50 Hallmark GSVA scores already published in ",
    "mmc4.xlsx; raw counts were not reprocessed."
  ),
  paste0(
    "No arbitrary RNA recovery threshold or universal optimal k is imposed; ",
    "the analysis reports the empirical learning curves."
  ),
  paste0(
    "Every possible RNA-sector subset was enumerated; patient is the ",
    "statistical unit; cohort CIs use exact-bootstrap."
  ),
  paste0(
    "The supplied sector-level RNA supplement contains 111 tumour-bearing ",
    "patients / 406 tumour sectors; this script follows those actual data."
  ),
  "",
  sprintf(
    "Elapsed time: %.1f seconds",
    proc.time()[["elapsed"]] - t0
  )
)

writeLines(
  summary_lines,
  file.path(
    ext_dir,
    "21_rna_pathway_sampling_summary.txt"
  )
)


cat(
  "\n",
  paste(
    summary_lines,
    collapse = "\n"
  ),
  "\n",
  sep = ""
)
