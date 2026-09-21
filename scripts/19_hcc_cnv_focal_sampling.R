# ==============================================================================
# 19_hcc_cnv_focal_sampling.R
#
# Purpose
#   Sampling-depth analysis of the seven recurrent focal/cytoband CNV events
#   reported in the published 2024 study, using the authoritative per-sector calls
#   supplied in Supplementary Table mmc4.xlsx.
#
# Focal events
#   AMP 11q13.3
#   DEL 10q25.2
#   DEL 13q14.13
#   DEL 17p13.3
#   DEL 1p36.33
#   DEL 21p12
#   DEL 8p23.1
#
# Design
#   - No reconstruction from cnv_gene.tsv.
#   - No new CNV calling threshold.
#   - The WT/AMP/DEL calls in mmc4.xlsx are used exactly as published.
#   - Before sampling, the script reproduces the seven published patient-level
#     trunk-ratio percentages from the full 490-sector table.
#   - Full available-sector reference, not whole-tumour truth.
#   - Primary cohort: fixed n>=5, k=1..4; hence every point is true
#     downsampling (k<n).
#   - Expected recovery is exact from event occupancy via the hypergeometric
#     identity; there is no Monte-Carlo sector sampling.
#   - Patient is the statistical unit.
#   - Cohort summary = median across patients.
#   - 95% CI = exact empirical nonparametric bootstrap of the median using the
#     validated exact-bootstrap algorithm.
#   - Paired comparison with the validated protein-altering SNV/indel recovery
#     is performed on the same patients with a defined focal-CNV endpoint.
#
# Required prerequisites
#   - 00_config.R
#   - validated Script 01/02 intermediates
#   - mmc4.xlsx (Supplementary Table, sheet merged_DNA_sample_table)
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

cat("Script 19: HCC focal CNV sampling\n")
cat("=====================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------- constants -----------------------------------
CI_LEVEL <- 0.95
PRIMARY_K <- 1:4
PRIMARY_MIN_N <- 5L

EXPECTED_N_PATIENTS <- 123L
EXPECTED_N_SECTORS <- 490L
EXPECTED_FIXED_N_GE5 <- 55L

EXPECTED_MMC4_MD5 <- "33c14fc9e52e32cca6ea8c9c771d8985"
MMC4_SHEET <- "merged_DNA_sample_table"

FOCAL_EVENTS <- data.table(
  event_id = c(
    "DEL_10q25.2",
    "AMP_11q13.3",
    "DEL_13q14.13",
    "DEL_17p13.3",
    "DEL_1p36.33",
    "DEL_21p12",
    "DEL_8p23.1"
  ),
  source_column = c(
    "10q25.2",
    "11q13.3",
    "13q14.13",
    "17p13.3",
    "1p36.33",
    "21p12",
    "8p23.1"
  ),
  direction = c(
    "DEL",
    "AMP",
    "DEL",
    "DEL",
    "DEL",
    "DEL",
    "DEL"
  ),
  # Integer percentages displayed in the published 2024 Figure 1B.
  published_trunk_ratio_percent = c(
    25, 30, 25, 0, 60, 24, 36
  )
)

EVENT_SET_ID <- "seven_recurrent_focal_cnv_events"

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

canonicalize_probability <- function(x, tol = 1e-12, label = "probability") {
  nm <- names(x)
  x <- as.numeric(x)
  names(x) <- nm

  finite <- is.finite(x)

  if (any(x[finite] < -tol | x[finite] > 1 + tol)) {
    stop(label, " outside [0,1] beyond numerical tolerance.")
  }

  x[finite & abs(x) <= tol] <- 0
  x[finite & abs(x - 1) <= tol] <- 1
  x
}

canonicalize_zero <- function(x, tol = 1e-12) {
  nm <- names(x)
  x <- as.numeric(x)
  names(x) <- nm
  x[is.finite(x) & abs(x) <= tol] <- 0
  x
}

# Exact nonparametric bootstrap functions: shared implementation, sourced
# near the top of this script (functions_bootstrap_exact.R).

# --------------------- exact expected recovery -------------------------------
expected_recovery_from_occupancies <- function(m, n, k) {
  m <- as.integer(m)

  if (length(m) == 0L) {
    return(c(
      detection = NA_real_,
      classification = NA_real_
    ))
  }

  if (any(m <= 0L | m >= n)) {
    stop(
      "Non-ubiquitous occupancy vector contains invalid m."
    )
  }

  if (k < 1L || k > n) {
    stop("Invalid k in expected-recovery calculation.")
  }

  denom <- choose(n, k)

  p_not_detected <- choose(n - m, k) / denom
  p_apparent_ubiquity <- choose(m, k) / denom

  out <- c(
    detection =
      mean(1 - p_not_detected),
    classification =
      mean(
        1 -
          p_not_detected -
          p_apparent_ubiquity
      )
  )

  canonicalize_probability(
    out,
    label = "expected focal-CNV recovery"
  )
}

# ------------------------------ output paths ---------------------------------
ext_dir <- file.path(
  PATHS$output_dir,
  "results",
  "extensions",
  "cnv"
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

# ------------------------------ locate input ---------------------------------
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

protein_patient_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_unconstrained_patient_summary_protein_altering.rds"
)

required_paths <- c(
  sample_map_path,
  protein_patient_path
)

missing_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_paths) > 0L) {
  stop(
    "Missing validated prerequisite file(s):\n",
    paste0(
      "  ",
      missing_paths,
      collapse = "\n"
    )
  )
}

# --------------------------- validated WGS cohort -----------------------------
sample_map <- as.data.table(
  readRDS(sample_map_path)
)

protein_patient <- as.data.table(
  readRDS(protein_patient_path)
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

require_columns(
  protein_patient,
  c(
    "patient_id",
    "n",
    "k",
    "recall_nonubiquitous_detection_mean",
    "recall_heterogeneity_classification_mean"
  ),
  "validated protein-altering patient summary"
)

sample_map[, `:=`(
  patient_id = as.character(patient_id),
  sample_id = as.character(sample_id),
  n_sectors = as.integer(n_sectors)
)]

protein_patient[
  ,
  patient_id := as.character(patient_id)
]

if (nrow(sample_map) != EXPECTED_N_SECTORS ||
    uniqueN(sample_map$patient_id) != EXPECTED_N_PATIENTS ||
    uniqueN(sample_map$sample_id) != EXPECTED_N_SECTORS ||
    anyDuplicated(sample_map$sample_id) > 0L) {
  stop(
    "Validated WGS sample-map regression failed."
  )
}

fixed_patients <- unique(
  sample_map[
    n_sectors >= PRIMARY_MIN_N,
    .(
      patient_id,
      n = as.integer(n_sectors)
    )
  ]
)

setorder(
  fixed_patients,
  patient_id
)

if (nrow(fixed_patients) != EXPECTED_FIXED_N_GE5) {
  stop(
    "Fixed n>=5 cohort-size regression failed: expected ",
    EXPECTED_FIXED_N_GE5,
    ", observed ",
    nrow(fixed_patients),
    "."
  )
}

if (any(fixed_patients$n <= max(PRIMARY_K))) {
  stop(
    "Primary fixed cohort contains a k=n point."
  )
}

protein_anchor <- stats::median(
  protein_patient[
    n >= 5L &
      k == 4L,
    recall_nonubiquitous_detection_mean
  ],
  na.rm = TRUE
)

if (abs(
  protein_anchor -
    0.843939393939394
) > 1e-12) {
  stop(
    "Validated protein-altering k=4 regression anchor failed."
  )
}

# ------------------------------- read mmc4 -----------------------------------
cat("Reading mmc4.xlsx / merged_DNA_sample_table ...\n")

sheet_names <- readxl::excel_sheets(mmc4_path)

if (!(MMC4_SHEET %in% sheet_names)) {
  stop(
    "Required sheet '", MMC4_SHEET,
    "' not found in mmc4.xlsx."
  )
}

dna <- as.data.table(
  readxl::read_excel(
    mmc4_path,
    sheet = MMC4_SHEET,
    .name_repair = "minimal"
  )
)

required_mmc4_cols <- c(
  "Publication_ID",
  "DNA_lib",
  "Tumor",
  "Sample_type",
  FOCAL_EVENTS$source_column
)

require_columns(
  dna,
  required_mmc4_cols,
  "mmc4.xlsx merged_DNA_sample_table"
)

dna[, `:=`(
  Publication_ID = as.character(Publication_ID),
  DNA_lib = as.character(DNA_lib),
  Tumor = as.character(Tumor),
  Sample_type = as.character(Sample_type)
)]

if (nrow(dna) != EXPECTED_N_SECTORS ||
    uniqueN(dna$Publication_ID) != EXPECTED_N_PATIENTS ||
    uniqueN(dna$DNA_lib) != EXPECTED_N_SECTORS ||
    anyDuplicated(dna$DNA_lib) > 0L) {
  stop(
    "mmc4 merged_DNA_sample_table is not 123 patients / 490 unique DNA libraries."
  )
}

if (anyNA(dna$Publication_ID) ||
    anyNA(dna$DNA_lib) ||
    any(!nzchar(dna$Publication_ID)) ||
    any(!nzchar(dna$DNA_lib))) {
  stop(
    "mmc4 contains missing/empty Publication_ID or DNA_lib."
  )
}

# The sheet should already contain tumour sectors only; freeze this explicitly.
sample_type_values <- sort(unique(
  trimws(dna$Sample_type)
))

if (!identical(
  sample_type_values,
  "Tumor"
)) {
  stop(
    "Unexpected Sample_type domain in merged_DNA_sample_table: ",
    paste(sample_type_values, collapse = ", ")
  )
}

# --------------------------- library crosswalk -------------------------------
crosswalk <- merge(
  dna[, .(
    Publication_ID,
    DNA_lib,
    Tumor
  )],
  sample_map[, .(
    sample_id,
    patient_id,
    n_sectors
  )],
  by.x = "DNA_lib",
  by.y = "sample_id",
  all = TRUE,
  sort = FALSE
)

if (nrow(crosswalk) != EXPECTED_N_SECTORS ||
    anyNA(crosswalk$DNA_lib) ||
    anyNA(crosswalk$Publication_ID) ||
    anyNA(crosswalk$patient_id) ||
    anyNA(crosswalk$n_sectors)) {
  stop(
    "mmc4 DNA_lib -> validated WGS sample crosswalk is not complete."
  )
}

publication_patient_map <- unique(
  crosswalk[, .(
    Publication_ID,
    patient_id
  )]
)

if (nrow(publication_patient_map) != EXPECTED_N_PATIENTS ||
    publication_patient_map[
      ,
      .N,
      by = Publication_ID
    ][
      ,
      max(N)
    ] != 1L ||
    publication_patient_map[
      ,
      .N,
      by = patient_id
    ][
      ,
      max(N)
    ] != 1L) {
  stop(
    "Publication_ID and validated patient_id are not one-to-one."
  )
}

sector_count_regression <- crosswalk[, .(
  n_mmc4 = .N,
  n_validated = unique(n_sectors)[1L]
), by = .(
  Publication_ID,
  patient_id
)]

if (any(
  sector_count_regression$n_mmc4 !=
    sector_count_regression$n_validated
)) {
  stop(
    "mmc4 and validated WGS sector counts disagree for at least one patient."
  )
}

if (DEBUG) {
  fwrite(
    crosswalk,
    file.path(
      ext_qc_dir,
      "19_focal_cnv_mmc4_wgs_crosswalk.tsv"
    ),
    sep = "\t"
  )
}

# ----------------------- event-call domain QC --------------------------------
call_domain_rows <- list()

for (ii in seq_len(nrow(FOCAL_EVENTS))) {
  col <- FOCAL_EVENTS$source_column[[ii]]
  direction <- FOCAL_EVENTS$direction[[ii]]

  vals <- trimws(as.character(dna[[col]]))

  if (anyNA(vals) || any(!nzchar(vals))) {
    stop(
      "Focal CNV column ", col,
      " contains missing/empty calls."
    )
  }

  unexpected <- setdiff(
    sort(unique(vals)),
    c("WT", "AMP", "DEL")
  )

  if (length(unexpected) > 0L) {
    stop(
      "Unexpected call value(s) in focal column ",
      col, ": ",
      paste(unexpected, collapse = ", ")
    )
  }

  tab <- as.data.table(table(
    factor(
      vals,
      levels = c("WT", "AMP", "DEL")
    )
  ))

  setnames(
    tab,
    c("call", "N")
  )

  tab[, `:=`(
    event_id = FOCAL_EVENTS$event_id[[ii]],
    source_column = col,
    target_direction = direction
  )]

  call_domain_rows[[ii]] <- tab
}

call_domain <- rbindlist(
  call_domain_rows,
  use.names = TRUE
)

setcolorder(
  call_domain,
  c(
    "event_id",
    "source_column",
    "target_direction",
    "call",
    "N"
  )
)

if (DEBUG) {
  fwrite(
    call_domain,
    file.path(
      ext_qc_dir,
      "19_focal_cnv_call_domain.tsv"
    ),
    sep = "\t"
  )
}

# -------------------------- long focal call table -----------------------------
long_rows <- vector(
  "list",
  nrow(FOCAL_EVENTS)
)

for (ii in seq_len(nrow(FOCAL_EVENTS))) {
  event_id <- FOCAL_EVENTS$event_id[[ii]]
  col <- FOCAL_EVENTS$source_column[[ii]]
  direction <- FOCAL_EVENTS$direction[[ii]]

  long_rows[[ii]] <- data.table(
    DNA_lib = dna$DNA_lib,
    Publication_ID = dna$Publication_ID,
    event_id = event_id,
    source_column = col,
    direction = direction,
    source_call = trimws(
      as.character(dna[[col]])
    )
  )
}

focal_long <- rbindlist(
  long_rows,
  use.names = TRUE
)

focal_long[, event_positive :=
             source_call == direction]

focal_long <- merge(
  focal_long,
  crosswalk[, .(
    DNA_lib,
    patient_id,
    n_sectors
  )],
  by = "DNA_lib",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(focal_long$patient_id) ||
    anyNA(focal_long$n_sectors)) {
  stop(
    "Focal event table failed WGS crosswalk."
  )
}

if (anyDuplicated(
  focal_long[, .(
    DNA_lib,
    event_id
  )]
) > 0L) {
  stop(
    "Duplicate DNA_lib/event_id in focal call table."
  )
}

if (nrow(focal_long) !=
    EXPECTED_N_SECTORS *
      nrow(FOCAL_EVENTS)) {
  stop(
    "Focal long table is not 490 sectors x 7 events."
  )
}

# ---------------------- full available-sector reference -----------------------
reference <- focal_long[, .(
  n_sectors = unique(n_sectors)[1L],
  full_count = sum(event_positive),
  n_opposite_direction_calls = sum(
    source_call != "WT" &
      source_call != direction
  )
), by = .(
  patient_id,
  Publication_ID,
  event_id,
  source_column,
  direction
)]

reference[, `:=`(
  full_count = as.integer(full_count),
  is_present = full_count > 0L,
  is_ubiquitous =
    full_count > 0L &
      full_count == n_sectors,
  is_nonubiquitous =
    full_count > 0L &
      full_count < n_sectors,
  is_private =
    full_count == 1L
)]

if (nrow(reference) !=
    EXPECTED_N_PATIENTS *
      nrow(FOCAL_EVENTS)) {
  stop(
    "Focal reference is not 123 patients x 7 events."
  )
}

if (any(reference$full_count < 0L) ||
    any(
      reference$full_count >
        reference$n_sectors
    )) {
  stop(
    "Invalid focal-CNV occupancy relative to patient sector count."
  )
}

setorder(
  reference,
  patient_id,
  event_id
)

if (DEBUG) {
  fwrite(
    reference,
    file.path(
      ext_dir,
      "19_focal_cnv_full_reference.tsv"
    ),
    sep = "\t"
  )
}

# -------------------- reproduce published trunk ratios ------------------------
trunk_qc <- reference[, {
  present <- is_present %in% TRUE

  n_present <- sum(present)
  n_trunk <- sum(
    is_ubiquitous %in% TRUE &
      present
  )

  list(
    n_patients_with_event = n_present,
    n_trunk_patients = n_trunk,
    reconstructed_trunk_ratio_percent =
      if (n_present > 0L) {
        100 * n_trunk / n_present
      } else {
        NA_real_
      },
    n_nonubiquitous_patients =
      sum(is_nonubiquitous),
    n_private_patients =
      sum(is_private)
  )
}, by = .(
  event_id,
  source_column,
  direction
)]

trunk_qc <- merge(
  FOCAL_EVENTS,
  trunk_qc,
  by = c(
    "event_id",
    "source_column",
    "direction"
  ),
  all.x = TRUE,
  sort = FALSE
)

trunk_qc[, `:=`(
  rounded_reconstructed_percent =
    round(reconstructed_trunk_ratio_percent),
  published_ratio_match =
    round(reconstructed_trunk_ratio_percent) ==
      published_trunk_ratio_percent
)]

trunk_qc[, focal_event_order := match(
  event_id,
  FOCAL_EVENTS$event_id
)]

setorder(
  trunk_qc,
  focal_event_order
)

trunk_qc[, focal_event_order := NULL]

if (anyNA(trunk_qc$published_ratio_match) ||
    !all(trunk_qc$published_ratio_match)) {
  if (DEBUG) {
    fwrite(
      trunk_qc,
      file.path(
        ext_qc_dir,
        "19_focal_cnv_published_trunk_ratio_qc.tsv"
      ),
      sep = "\t"
    )
  }

  stop(
    "Failed to reproduce all seven published focal-CNV trunk ratios."
  )
}

if (DEBUG) {
  fwrite(
    trunk_qc,
    file.path(
      ext_qc_dir,
      "19_focal_cnv_published_trunk_ratio_qc.tsv"
    ),
    sep = "\t"
  )
}

# ------------------------ event reference summary -----------------------------
event_reference_summary <- reference[, .(
  n_patients_present =
    sum(is_present),
  n_patients_ubiquitous =
    sum(is_ubiquitous),
  n_patients_nonubiquitous =
    sum(is_nonubiquitous),
  n_patients_private =
    sum(is_private),
  total_opposite_direction_sector_calls =
    sum(n_opposite_direction_calls)
), by = .(
  event_id,
  source_column,
  direction
)]

event_reference_summary[, focal_event_order := match(
  event_id,
  FOCAL_EVENTS$event_id
)]

setorder(
  event_reference_summary,
  focal_event_order
)

event_reference_summary[
  ,
  focal_event_order := NULL
]

if (DEBUG) {
  fwrite(
    event_reference_summary,
    file.path(
      ext_dir,
      "19_focal_cnv_event_reference_summary.tsv"
    ),
    sep = "\t"
  )
}

# --------------------- exact sampling expectation -----------------------------
rows <- vector(
  "list",
  nrow(fixed_patients) *
    length(PRIMARY_K)
)

ii <- 0L

for (pp in seq_len(nrow(fixed_patients))) {
  pid <- fixed_patients$patient_id[[pp]]
  n <- fixed_patients$n[[pp]]

  rr <- reference[
    patient_id == pid &
      is_nonubiquitous == TRUE
  ]

  m <- rr$full_count
  n_nonubi <- length(m)

  n_present <- reference[
    patient_id == pid &
      is_present == TRUE,
    .N
  ]

  for (kval in PRIMARY_K) {
    if (kval >= n) {
      stop(
        "True-downsampling invariant violated for patient ",
        pid, ": k=", kval, ", n=", n
      )
    }

    rec <- expected_recovery_from_occupancies(
      m = m,
      n = n,
      k = kval
    )

    ii <- ii + 1L

    rows[[ii]] <- data.table(
      event_set = EVENT_SET_ID,
      patient_id = pid,
      n = n,
      k = as.integer(kval),
      n_reference_events_present =
        as.integer(n_present),
      n_nonubiquitous_reference =
        as.integer(n_nonubi),
      expected_detection_recall =
        unname(rec[["detection"]]),
      expected_classification_recall =
        unname(rec[["classification"]])
    )
  }
}

patient_expected <- rbindlist(
  rows,
  use.names = TRUE
)

patient_expected[
  ,
  expected_detection_recall :=
    canonicalize_probability(
      expected_detection_recall,
      label = "focal expected_detection_recall"
    )
]

patient_expected[
  ,
  expected_classification_recall :=
    canonicalize_probability(
      expected_classification_recall,
      label = "focal expected_classification_recall"
    )
]

if (nrow(patient_expected) !=
    EXPECTED_FIXED_N_GE5 *
      length(PRIMARY_K) ||
    uniqueN(patient_expected$patient_id) !=
      EXPECTED_FIXED_N_GE5 ||
    anyDuplicated(
      patient_expected[, .(
        patient_id,
        k
      )]
    ) > 0L) {
  stop(
    "Focal patient-expectation table has unexpected dimensions."
  )
}

k1_class <- patient_expected[
  k == 1L &
    is.finite(
      expected_classification_recall
    ),
  expected_classification_recall
]

if (any(k1_class != 0)) {
  stop(
    "Focal classification recall at k=1 is not exactly zero."
  )
}

valid_pair <- patient_expected[
  is.finite(expected_detection_recall) &
    is.finite(expected_classification_recall)
]

if (nrow(valid_pair) > 0L &&
    any(
      valid_pair$expected_classification_recall >
        valid_pair$expected_detection_recall +
          1e-12
    )) {
  stop(
    "Focal classification recall exceeds detection recall."
  )
}

setorder(
  patient_expected,
  patient_id,
  k
)

if (DEBUG) {
  fwrite(
    patient_expected,
    file.path(
      ext_dir,
      "19_focal_cnv_patient_expected_recovery_fixed_n_ge5.tsv"
    ),
    sep = "\t"
  )
}

# ----------------------------- cohort summary --------------------------------
metric_map <- c(
  detection =
    "expected_detection_recall",
  heterogeneity_classification =
    "expected_classification_recall"
)

summary_rows <- list()
ss <- 0L

for (metric_id in names(metric_map)) {
  value_col <- metric_map[[metric_id]]

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
        "fixed_n_ge5_true_downsampling",
      event_set = EVENT_SET_ID,
      metric = metric_id,
      k = as.integer(kval),
      n_fixed_cohort =
        EXPECTED_FIXED_N_GE5,
      n_patients_with_defined_endpoint =
        length(vals),
      fraction_fixed_cohort_with_defined_endpoint =
        length(vals) /
          EXPECTED_FIXED_N_GE5,
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

fwrite(
  depth_summary,
  file.path(
    ext_dir,
    "19_focal_cnv_depth_summary_fixed_n_ge5.tsv"
  ),
  sep = "\t"
)

# ----------------------- matched protein comparator ---------------------------
protein_comp <- protein_patient[
  patient_id %in%
    fixed_patients$patient_id &
    n >= PRIMARY_MIN_N &
    k %in% PRIMARY_K,
  .(
    patient_id,
    n = as.integer(n),
    k = as.integer(k),
    protein_detection =
      as.numeric(
        recall_nonubiquitous_detection_mean
      ),
    protein_classification =
      as.numeric(
        recall_heterogeneity_classification_mean
      )
  )
]

setorder(
  protein_comp,
  patient_id,
  k
)

if (nrow(protein_comp) !=
    EXPECTED_FIXED_N_GE5 *
      length(PRIMARY_K) ||
    uniqueN(protein_comp$patient_id) !=
      EXPECTED_FIXED_N_GE5 ||
    anyDuplicated(
      protein_comp[, .(
        patient_id,
        k
      )]
    ) > 0L) {
  stop(
    "Validated protein comparator fixed-cohort regression failed."
  )
}

paired <- merge(
  patient_expected[, .(
    patient_id,
    n,
    k,
    focal_detection =
      expected_detection_recall,
    focal_classification =
      expected_classification_recall
  )],
  protein_comp,
  by = c(
    "patient_id",
    "n",
    "k"
  ),
  all = FALSE
)

paired[
  ,
  diff_detection :=
    canonicalize_zero(
      focal_detection -
        protein_detection
    )
]

paired[
  ,
  diff_classification :=
    canonicalize_zero(
      focal_classification -
        protein_classification
    )
]

paired_summary_rows <- list()
ps <- 0L

for (metric_id in c(
  "detection",
  "heterogeneity_classification"
)) {
  diff_col <- if (
    metric_id == "detection"
  ) {
    "diff_detection"
  } else {
    "diff_classification"
  }

  for (kval in PRIMARY_K) {
    vals <- paired[
      k == kval,
      get(diff_col)
    ]

    vals <- as.numeric(
      vals[is.finite(vals)]
    )

    ci <- exact_bootstrap_median_ci(
      vals,
      level = CI_LEVEL,
      min_n = 1L
    )

    ps <- ps + 1L

    paired_summary_rows[[ps]] <- data.table(
      analysis =
        "paired_fixed_n_ge5_true_downsampling",
      cnv_event_set =
        EVENT_SET_ID,
      comparator_event_set =
        "protein_altering_all_genes",
      metric =
        metric_id,
      k =
        as.integer(kval),
      n_paired_patients =
        length(vals),
      median_focal_minus_protein =
        safe_median(vals),
      ci_lower =
        unname(ci[["lower"]]),
      ci_upper =
        unname(ci[["upper"]]),
      proportion_focal_greater =
        if (length(vals) > 0L) {
          mean(vals > 0)
        } else {
          NA_real_
        },
      proportion_equal =
        if (length(vals) > 0L) {
          mean(abs(vals) <= 1e-12)
        } else {
          NA_real_
        },
      bootstrap_method =
        if (length(vals) > 0L) {
          "exact_nonparametric_percentile"
        } else {
          NA_character_
        }
    )
  }
}

paired_summary <- rbindlist(
  paired_summary_rows,
  use.names = TRUE
)

setorder(
  paired_summary,
  metric,
  k
)

fwrite(
  paired_summary,
  file.path(
    ext_dir,
    "19_focal_cnv_vs_protein_paired_fixed_n_ge5.tsv"
  ),
  sep = "\t"
)

# ----------------------- operational heuristic -------------------------------
det <- depth_summary[
  metric == "detection"
]

hit <- det[
  is.finite(median) &
    median >=
      PARAMS$adequate_recovery_threshold
]

heuristic <- data.table(
  event_set =
    EVENT_SET_ID,
  cohort =
    "fixed_n_ge5_true_downsampling",
  operational_threshold =
    PARAMS$adequate_recovery_threshold,
  first_tested_k_reaching_threshold =
    if (nrow(hit) > 0L) {
      min(hit$k)
    } else {
      NA_integer_
    },
  interpretation =
    paste0(
      "Legacy operational benchmark only; ",
      "not a universal biological optimum."
    )
)

if (DEBUG) {
  fwrite(
    heuristic,
    file.path(
      ext_dir,
      "19_focal_cnv_sampling_operational_heuristic.tsv"
    ),
    sep = "\t"
  )
}

# ------------------------------- provenance ----------------------------------
n_info_fixed <- uniqueN(
  reference[
    patient_id %in%
      fixed_patients$patient_id &
      is_nonubiquitous == TRUE,
    patient_id
  ]
)

provenance <- data.table(
  item = c(
    "mmc4_path",
    "mmc4_md5",
    "mmc4_sheet",
    "mmc4_rows",
    "mmc4_patients",
    "mmc4_DNA_libraries",
    "focal_event_definitions",
    "validated_wgs_patients",
    "validated_wgs_sectors",
    "fixed_n_ge5_patients",
    "fixed_patients_with_nonubiquitous_focal_event",
    "published_trunk_ratio_qc_passed"
  ),
  value = c(
    mmc4_path,
    mmc4_md5,
    MMC4_SHEET,
    as.character(nrow(dna)),
    as.character(
      uniqueN(dna$Publication_ID)
    ),
    as.character(
      uniqueN(dna$DNA_lib)
    ),
    as.character(
      nrow(FOCAL_EVENTS)
    ),
    as.character(
      uniqueN(sample_map$patient_id)
    ),
    as.character(
      nrow(sample_map)
    ),
    as.character(
      nrow(fixed_patients)
    ),
    as.character(
      n_info_fixed
    ),
    "TRUE"
  )
)

if (DEBUG) {
  fwrite(
    provenance,
    file.path(
      ext_qc_dir,
      "19_focal_cnv_sampling_provenance.tsv"
    ),
    sep = "\t"
  )
}

# -------------------------- human-readable summary ----------------------------
k4_detection <- depth_summary[
  metric == "detection" &
    k == 4L
]

k4_classification <- depth_summary[
  metric ==
    "heterogeneity_classification" &
    k == 4L
]

k4_pair <- paired_summary[
  metric == "detection" &
    k == 4L
]

first_k <- heuristic[
  ,
  first_tested_k_reaching_threshold
][[1L]]

heuristic_text <- if (is.na(first_k)) {
  "No tested k=1..4 reaches the legacy 80% median-detection benchmark."
} else {
  paste0(
    "k=", first_k,
    " is the first TESTED depth reaching the legacy 80% ",
    "median-detection benchmark; not a universal biological optimum."
  )
}

trunk_lines <- vapply(
  seq_len(nrow(trunk_qc)),
  function(ii) {
    sprintf(
      "  %s: %.1f%% -> rounded %d%%; published %d%% [PASS]",
      trunk_qc$event_id[[ii]],
      trunk_qc$reconstructed_trunk_ratio_percent[[ii]],
      trunk_qc$rounded_reconstructed_percent[[ii]],
      trunk_qc$published_trunk_ratio_percent[[ii]]
    )
  },
  character(1)
)

summary_lines <- c(
  "HCC extensions - Script 19 focal CNV sampling",
  "=====================================================",
  "",
  paste0(
    "mmc4.xlsx MD5: ",
    mmc4_md5,
    " [PASS]"
  ),
  sprintf(
    "merged_DNA_sample_table: %d patients / %d DNA libraries [PASS]",
    uniqueN(dna$Publication_ID),
    uniqueN(dna$DNA_lib)
  ),
  sprintf(
    "mmc4 DNA_lib -> validated WGS crosswalk: %d/%d sectors [PASS]",
    nrow(crosswalk),
    EXPECTED_N_SECTORS
  ),
  "",
  "Published focal-CNV trunk-ratio regression:",
  trunk_lines,
  "",
  sprintf(
    "Seven recurrent focal CNV events; fixed n>=5 cohort: %d total; %d with >=1 non-ubiquitous focal event",
    nrow(fixed_patients),
    n_info_fixed
  ),
  sprintf(
    "k=4 focal-CNV detection median = %.6f [%.6f, %.6f]",
    k4_detection$median,
    k4_detection$ci_lower,
    k4_detection$ci_upper
  ),
  sprintf(
    "k=4 focal-CNV heterogeneity-classification median = %.6f [%.6f, %.6f]",
    k4_classification$median,
    k4_classification$ci_lower,
    k4_classification$ci_upper
  ),
  sprintf(
    "k=4 paired focal-CNV - protein detection difference = %.6f [%.6f, %.6f] (N=%d)",
    k4_pair$median_focal_minus_protein,
    k4_pair$ci_lower,
    k4_pair$ci_upper,
    k4_pair$n_paired_patients
  ),
  paste0(
    "Operational 80% heuristic: ",
    heuristic_text
  ),
  paste0(
    "Guardrail: this analysis represents the seven recurrent focal CNV ",
    "events explicitly called in the published supplementary table; it is ",
    "not claimed to represent every possible focal CNV in the genome."
  ),
  paste0(
    "Source WT/AMP/DEL calls are used exactly as supplied; ",
    "no CNV threshold is re-estimated."
  ),
  paste0(
    "Sampling expectation is exact/hypergeometric; cohort CI is ",
    "exact-bootstrap; no Monte-Carlo sampling was used."
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
    "19_focal_cnv_sampling_summary.txt"
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
