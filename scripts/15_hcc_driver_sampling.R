# ==============================================================================
# 15_hcc_driver_sampling.R
#
#
# Primary driver-event rule (frozen here)
#   1. Exact SNV/indel events in the supplied 48-driver + DICER1 gene list that
#      are HIGH/MODERATE under the already validated Script 02 annotation rule;
#   2. plus the two canonical TERT promoter hotspots in GRCh38:
#        C228T / c.-124C>T : chr5:1295113 G>A
#        C250T / c.-146C>T : chr5:1295135 G>A
#      recovered from the official source-data snv_indel.tsv.
#
# Important design choices
#   - TERT hotspot calls are subjected to the SAME numeric filters used in the
#     validated WGS pipeline: ALT >= 3, total depth >= 10, VAF >= 0.05.
#   - No other LOW/MODIFIER/intronic/flanking variants in driver genes are used.
#   - The unit of an event remains the exact patient-specific variant.
#   - Full available-sector reference, not whole-tumour truth.
#   - Patient is the statistical unit.
#   - For depth k, expected recall is computed analytically from event occupancy
#     using the exact hypergeometric identities already validated against the
#     exhaustive Script 04 enumeration. Thus there is NO Monte-Carlo sampling.
#   - Primary cohort: fixed n >= 5, k = 1..4 (all points true downsampling).
#   - 95% CIs: exact empirical nonparametric bootstrap of the cohort median,
#     using the exact-bootstrap algorithm.
#
# Required prerequisites
#   - 00_config.R
#   - Script 14 completed successfully
#   - validated Script 01/02 intermediates
#   - official source-data snv_indel.tsv (or a .zip containing it)
#   - cnv_arm_level.tsv supplied previously; used ONLY as a frozen crosswalk
#     from source-data ITH patient IDs to WGS library IDs. No CNV analysis is
#     performed in this script.
#
# Optional environment variables if automatic file discovery fails
#   SNV_INDEL_PATH      = /full/path/snv_indel.tsv[.zip]
#   CNV_ARM_LEVEL_PATH  = /full/path/cnv_arm_level.tsv
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

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

cat("Script 15: HCC driver sampling\n")
cat("================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ constants ------------------------------------
CI_LEVEL <- 0.95
PRIMARY_K <- 1:4
PRIMARY_MIN_N <- 5L
EXPECTED_FIXED_N_GE5 <- 55L

EXPECTED_SNV_INDEL_MD5 <- "45ab1b9aebea41e66cb986dc9474d729"
EXPECTED_CNV_ARM_LEVEL_MD5 <- "6608ba4096df53bfa014877f5e952052"

EXPECTED_TERT_SOURCE_ROWS_RAW <- 137L
EXPECTED_TERT_SOURCE_ROWS_FILTERED <- 133L
EXPECTED_TERT_PATIENTS <- 42L
EXPECTED_TERT_C228T_PATIENTS <- 39L
EXPECTED_TERT_C250T_PATIENTS <- 3L

TERT_HOTSPOTS <- data.table(
  CHROM_SOURCE = "chr5",
  POS_SOURCE = c(1295113L, 1295135L),
  REF_SOURCE = "G",
  ALT_SOURCE = "A",
  hotspot_label = c("TERT_C228T", "TERT_C250T")
)

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

safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  stats::median(x)
}

# Canonicalize probabilities only at machine-precision boundaries. This does
# not alter substantive estimates; it prevents values such as -2.8e-17 from
# propagating into summaries/CIs when the mathematical value is exactly 0.
canonicalize_probability <- function(x, tol = 1e-12, label = "probability") {
  original_names <- names(x)
  x <- as.numeric(x)
  names(x) <- original_names
  finite <- is.finite(x)
  if (any(x[finite] < -tol | x[finite] > 1 + tol)) {
    bad <- x[finite][x[finite] < -tol | x[finite] > 1 + tol]
    stop(
      label, " outside [0,1] beyond numerical tolerance; range = ",
      paste(range(bad), collapse = " to ")
    )
  }
  x[finite & abs(x) <= tol] <- 0
  x[finite & abs(x - 1) <= tol] <- 1
  x
}

canonicalize_zero <- function(x, tol = 1e-12) {
  original_names <- names(x)
  x <- as.numeric(x)
  names(x) <- original_names
  x[is.finite(x) & abs(x) <= tol] <- 0
  x
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

materialize_snv_indel <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (!grepl("\\.zip$", path, ignore.case = TRUE)) return(path)

  listing <- unzip(path, list = TRUE)
  members <- as.character(listing$Name)
  hit <- members[
    basename(members) == "snv_indel.tsv" &
      !grepl("(^|/)__MACOSX(/|$)", members)
  ]
  if (length(hit) != 1L) {
    stop(
      "Expected exactly one non-__MACOSX snv_indel.tsv inside ", path,
      "; observed ", length(hit), "."
    )
  }

  tmp_dir <- tempfile("snv_indel_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  unzip(path, files = hit[[1L]], exdir = tmp_dir, junkpaths = TRUE)
  out <- file.path(tmp_dir, "snv_indel.tsv")
  if (!file.exists(out)) stop("Failed to extract snv_indel.tsv from zip.")
  out
}

# Exact nonparametric bootstrap functions: shared implementation, sourced
# near the top of this script (functions_bootstrap_exact.R).

# Exact expected detection / classification under uniform k-of-n sampling.
expected_recovery_from_occupancies <- function(m, n, k) {
  m <- as.integer(m)
  if (length(m) == 0L) {
    return(c(detection = NA_real_, classification = NA_real_))
  }
  if (any(m <= 0L | m >= n)) {
    stop("Occupancy vector supplied to non-ubiquitous recovery contains invalid m.")
  }
  if (k < 1L || k > n) stop("Invalid k in expected-recovery calculation.")

  denom <- choose(n, k)
  p_not_detected <- choose(n - m, k) / denom
  p_apparent_ubiquity <- choose(m, k) / denom

  out <- c(
    detection = mean(1 - p_not_detected),
    classification = mean(1 - p_not_detected - p_apparent_ubiquity)
  )
  canonicalize_probability(out, label = "expected recovery")
}

# ------------------------------ output paths ---------------------------------
ext_dir <- file.path(PATHS$output_dir, "results", "extensions", "driver")
ext_qc_dir <- file.path(PATHS$output_dir, "results", "qc", "extensions")
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ext_qc_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- prerequisites ---------------------------------
sample_map_path <- file.path(PATHS$intermediate_dir, "hcc_sample_map.rds")
all_filtered_patient_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_unconstrained_patient_summary_all_filtered.rds"
)
protein_patient_path <- file.path(
  PATHS$intermediate_dir,
  "hcc_unconstrained_patient_summary_protein_altering.rds"
)
audit_path <- file.path(PATHS$qc_dir, "01_tumor_sample_map_audit.tsv")
driver_candidate_path <- file.path(
  ext_dir,
  "14_driver_exact_event_candidates_preeligibility.tsv"
)

required_paths <- c(
  sample_map_path, all_filtered_patient_path, protein_patient_path,
  audit_path, driver_candidate_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing validated prerequisite file(s):\n",
    paste0("  ", missing_paths, collapse = "\n")
  )
}

snv_source_in <- locate_existing_file(
  "SNV_INDEL_PATH",
  c(
    file.path(PROJECT_ROOT, "source_data", "snv_indel.tsv"),
    file.path(PROJECT_ROOT, "snv_indel.tsv"),
    file.path(PROJECT_ROOT, "source_data", "snv_indel.tsv.zip"),
    file.path(PROJECT_ROOT, "snv_indel.tsv.zip"),
    file.path(PROJECT_ROOT, "source_data.zip"),
    file.path(dirname(PROJECT_ROOT), "source_data", "snv_indel.tsv"),
    file.path(dirname(PROJECT_ROOT), "snv_indel.tsv"),
    file.path(dirname(PROJECT_ROOT), "snv_indel.tsv.zip"),
    file.path(dirname(PROJECT_ROOT), "source_data.zip"),
    file.path(script_dir, "snv_indel.tsv"),
    file.path(script_dir, "snv_indel.tsv.zip"),
    file.path(script_dir, "source_data.zip")
  ),
  "official source-data snv_indel.tsv (or zip)"
)

cnv_arm_map_path <- locate_existing_file(
  "CNV_ARM_LEVEL_PATH",
  c(
    file.path(PROJECT_ROOT, "cnv_arm_level.tsv"),
    file.path(PROJECT_ROOT, "data", "cnv_arm_level.tsv"),
    file.path(dirname(PROJECT_ROOT), "cnv_arm_level.tsv"),
    file.path(dirname(PROJECT_ROOT), "data", "cnv_arm_level.tsv"),
    file.path(script_dir, "cnv_arm_level.tsv")
  ),
  "cnv_arm_level.tsv crosswalk file"
)

snv_tsv_path <- materialize_snv_indel(snv_source_in)
snv_md5 <- unname(tools::md5sum(snv_tsv_path))
if (!identical(snv_md5, EXPECTED_SNV_INDEL_MD5)) {
  stop(
    "snv_indel.tsv MD5 mismatch.\nExpected: ", EXPECTED_SNV_INDEL_MD5,
    "\nObserved: ", snv_md5,
    "\nFile: ", snv_tsv_path
  )
}

cnv_map_md5 <- unname(tools::md5sum(cnv_arm_map_path))
if (!identical(cnv_map_md5, EXPECTED_CNV_ARM_LEVEL_MD5)) {
  stop(
    "cnv_arm_level.tsv MD5 mismatch.\nExpected: ", EXPECTED_CNV_ARM_LEVEL_MD5,
    "\nObserved: ", cnv_map_md5,
    "\nFile: ", cnv_arm_map_path,
    "\nThis file is used only as the previously validated ITH->WGS library crosswalk."
  )
}

sample_map <- as.data.table(readRDS(sample_map_path))
all_filtered_patient <- as.data.table(readRDS(all_filtered_patient_path))
protein_patient <- as.data.table(readRDS(protein_patient_path))
audit <- fread(audit_path)
driver_candidates <- fread(driver_candidate_path)

require_columns(sample_map, c("patient_id", "sample_id", "n_sectors"), "sample_map")
require_columns(
  all_filtered_patient,
  c(
    "patient_id", "n", "k",
    "recall_nonubiquitous_detection_mean",
    "recall_heterogeneity_classification_mean"
  ),
  "validated all-filtered patient summary"
)
require_columns(
  protein_patient,
  c(
    "patient_id", "n", "k",
    "recall_nonubiquitous_detection_mean",
    "recall_heterogeneity_classification_mean"
  ),
  "validated protein-altering patient summary"
)
require_columns(
  audit,
  c("patient_id", "sample_id", "category", "n_sectors"),
  "01_tumor_sample_map_audit.tsv"
)
require_columns(
  driver_candidates,
  c(
    "event_key", "patient_id", "CHROM", "POS", "REF_ALLELE", "ALT_ALLELE",
    "driver_genes", "contains_DICER1", "any_high_moderate", "full_count",
    "n_sectors", "is_ubiquitous", "is_nonubiquitous", "is_private"
  ),
  "Script 14 driver candidate table"
)

sample_map[, patient_id := as.character(patient_id)]
sample_map[, sample_id := as.character(sample_id)]
all_filtered_patient[, patient_id := as.character(patient_id)]
protein_patient[, patient_id := as.character(patient_id)]
audit[, patient_id := as.character(patient_id)]
audit[, sample_id := as.character(sample_id)]
driver_candidates[, patient_id := as.character(patient_id)]

if (uniqueN(sample_map$patient_id) != 123L || nrow(sample_map) != 490L) {
  stop("Validated sample-map cohort is no longer 123 patients / 490 tumour sectors.")
}
if (anyDuplicated(sample_map$sample_id) > 0L) {
  stop("Duplicate sample_id in validated sample map.")
}
if (nrow(audit) != 490L || anyDuplicated(audit$sample_id) > 0L) {
  stop("Sample-map audit is not one row per 490 validated tumour sectors.")
}

# Validate audit <-> authoritative sample map.
map_regression <- merge(
  sample_map[, .(patient_id, sample_id, n_sectors_map = as.integer(n_sectors))],
  audit[, .(patient_id, sample_id, category, n_sectors_audit = as.integer(n_sectors))],
  by = c("patient_id", "sample_id"),
  all = TRUE
)
if (nrow(map_regression) != 490L || anyNA(map_regression$category) ||
    any(map_regression$n_sectors_map != map_regression$n_sectors_audit)) {
  stop("Validated sample map and Script 01 audit are not identical on cohort membership.")
}

# Primary fixed cohort QC.
fixed_patients <- unique(
  sample_map[n_sectors >= PRIMARY_MIN_N, .(patient_id, n = as.integer(n_sectors))]
)
setorder(fixed_patients, patient_id)
if (nrow(fixed_patients) != EXPECTED_FIXED_N_GE5) {
  stop(
    "fixed n>=5 cohort-size regression failed: expected ", EXPECTED_FIXED_N_GE5,
    ", observed ", nrow(fixed_patients), "."
  )
}
if (any(fixed_patients$n <= max(PRIMARY_K))) {
  stop("Primary fixed cohort contains k=n point; true-downsampling invariant failed.")
}

# Regression anchors for the already validated comparator patient summaries.
# These are central estimates only; Script 15 later recomputes exact CIs on the
# matched driver-eligible patient set where needed.
protein_anchor <- stats::median(
  protein_patient[n >= 5L & k == 4L, recall_nonubiquitous_detection_mean],
  na.rm = TRUE
)
all_filtered_anchor <- stats::median(
  all_filtered_patient[n >= 5L & k == 4L, recall_nonubiquitous_detection_mean],
  na.rm = TRUE
)
if (abs(protein_anchor - 0.843939393939394) > 1e-12 ||
    abs(all_filtered_anchor - 0.840473531382622) > 1e-12) {
  stop(
    "Validated comparator regression anchor failed: unexpected fixed-n>=5 k=4 detection median."
  )
}

# -----------------------------------------------------------------------------
# 1. Freeze crosswalk: source-data ITH patient -> validated WGS library/patient
# -----------------------------------------------------------------------------
cnv_crosswalk_raw <- fread(
  cnv_arm_map_path,
  select = c("PATIENT_ID", "Sample")
)
setnames(cnv_crosswalk_raw, c("PATIENT_ID", "Sample"), c("source_patient_id", "sample_id"))
cnv_crosswalk_raw[, source_patient_id := as.character(source_patient_id)]
cnv_crosswalk_raw[, sample_id := as.character(sample_id)]
cnv_crosswalk_raw <- unique(cnv_crosswalk_raw)

if (nrow(cnv_crosswalk_raw) != 490L ||
    uniqueN(cnv_crosswalk_raw$source_patient_id) != 123L ||
    uniqueN(cnv_crosswalk_raw$sample_id) != 490L ||
    anyDuplicated(cnv_crosswalk_raw$sample_id) > 0L) {
  stop("ITH-to-WGS crosswalk did not resolve to 123 patients / 490 unique sectors.")
}

source_crosswalk <- merge(
  cnv_crosswalk_raw,
  audit[, .(sample_id, patient_id, category, n_sectors = as.integer(n_sectors))],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(source_crosswalk$patient_id) || anyNA(source_crosswalk$category)) {
  stop("At least one source-data crosswalk WGS library is absent from validated Script 01 audit.")
}

patient_crosswalk_qc <- source_crosswalk[, .(
  n_validated_patients = uniqueN(patient_id),
  validated_patient_id = if (uniqueN(patient_id) == 1L) unique(patient_id) else NA_character_,
  n_sectors_crosswalk = .N,
  n_sectors_validated = unique(n_sectors)[1L]
), by = source_patient_id]

if (any(patient_crosswalk_qc$n_validated_patients != 1L) ||
    any(patient_crosswalk_qc$n_sectors_crosswalk != patient_crosswalk_qc$n_sectors_validated)) {
  stop("Source ITH patient -> validated patient crosswalk is not one-to-one / sector-complete.")
}

if (DEBUG) {
  fwrite(
    source_crosswalk,
    file.path(ext_qc_dir, "15_source_ITH_to_validated_WGS_crosswalk.tsv"),
    sep = "\t"
  )
}

# -----------------------------------------------------------------------------
# 2. Recover canonical TERT promoter calls from official source-data MAF table
# -----------------------------------------------------------------------------
tert_cols <- c(
  "PATIENT_ID", "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
  "Variant_Classification", "Variant_Type", "Reference_Allele",
  "Tumor_Seq_Allele2", "VAF", "t_alt_count", "t_ref_count", "SAMPLE_ID"
)

snv_header <- names(fread(snv_tsv_path, nrows = 0L))
miss_snv_cols <- setdiff(tert_cols, snv_header)
if (length(miss_snv_cols) > 0L) {
  stop(
    "snv_indel.tsv lacks required columns: ",
    paste(miss_snv_cols, collapse = ", ")
  )
}

# 138 MB is modest for fread; selecting only required columns keeps memory low.
snv_min <- fread(snv_tsv_path, select = tert_cols)

tert_raw <- snv_min[
  Hugo_Symbol == "TERT" &
    Chromosome == "chr5" &
    Start_Position %in% TERT_HOTSPOTS$POS_SOURCE &
    Reference_Allele == "G" &
    Tumor_Seq_Allele2 == "A"
]
rm(snv_min)
invisible(gc(FALSE))

if (nrow(tert_raw) != EXPECTED_TERT_SOURCE_ROWS_RAW) {
  stop(
    "Canonical TERT source-row regression failed: expected ",
    EXPECTED_TERT_SOURCE_ROWS_RAW, ", observed ", nrow(tert_raw), "."
  )
}
if (uniqueN(tert_raw$PATIENT_ID) != EXPECTED_TERT_PATIENTS) {
  stop("Canonical TERT patient-count regression failed before numeric filtering.")
}
if (any(tert_raw$Variant_Classification != "5'Flank") ||
    any(tert_raw$Variant_Type != "SNP") ||
    any(tert_raw$End_Position != tert_raw$Start_Position)) {
  stop("Unexpected annotation among canonical TERT hotspot rows.")
}

tert_raw[, total_depth := t_alt_count + t_ref_count]
tert_filtered <- tert_raw[
  !is.na(t_alt_count) & !is.na(total_depth) & !is.na(VAF) &
    t_alt_count >= PARAMS$min_alt_count &
    total_depth >= PARAMS$min_total_depth &
    VAF >= PARAMS$min_vaf
]

if (nrow(tert_filtered) != EXPECTED_TERT_SOURCE_ROWS_FILTERED) {
  stop(
    "Canonical TERT filtered-row regression failed: expected ",
    EXPECTED_TERT_SOURCE_ROWS_FILTERED, ", observed ", nrow(tert_filtered), "."
  )
}
if (uniqueN(tert_filtered$PATIENT_ID) != EXPECTED_TERT_PATIENTS) {
  stop("Canonical TERT patient count changed after applying validated numeric filters.")
}

tert_filtered <- merge(
  tert_filtered,
  TERT_HOTSPOTS,
  by.x = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"),
  by.y = c("CHROM_SOURCE", "POS_SOURCE", "REF_SOURCE", "ALT_SOURCE"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(tert_filtered$hotspot_label)) {
  stop("TERT hotspot-label merge failed.")
}

tert_patient_counts <- tert_filtered[, .(n_patients = uniqueN(PATIENT_ID)), by = hotspot_label]
if (tert_patient_counts[hotspot_label == "TERT_C228T", n_patients] != EXPECTED_TERT_C228T_PATIENTS ||
    tert_patient_counts[hotspot_label == "TERT_C250T", n_patients] != EXPECTED_TERT_C250T_PATIENTS) {
  stop("TERT C228T/C250T patient-count regression failed.")
}
tert_hotspots_per_patient <- tert_filtered[, .(n_hotspots = uniqueN(hotspot_label)), by = PATIENT_ID]
if (any(tert_hotspots_per_patient$n_hotspots > 1L)) {
  stop("A source-data patient carries both canonical TERT hotspots; expected mutual exclusivity here.")
}

# Source SAMPLE_ID encodes the tumour-category suffix after PATIENT_ID + "_".
tert_filtered[, PATIENT_ID := as.character(PATIENT_ID)]
tert_filtered[, SAMPLE_ID := as.character(SAMPLE_ID)]
prefix_ok <- vapply(
  seq_len(nrow(tert_filtered)),
  function(ii) startsWith(
    tert_filtered$SAMPLE_ID[[ii]],
    paste0(tert_filtered$PATIENT_ID[[ii]], "_")
  ),
  logical(1)
)
if (!all(prefix_ok)) {
  stop("At least one TERT source SAMPLE_ID does not start with PATIENT_ID + '_'.")
}
tert_filtered[, category := vapply(
  seq_len(.N),
  function(ii) substring(SAMPLE_ID[[ii]], nchar(PATIENT_ID[[ii]]) + 2L),
  character(1)
)]
if (anyNA(tert_filtered$category) || any(!nzchar(tert_filtered$category))) {
  stop("Failed to derive tumour category from at least one TERT source SAMPLE_ID.")
}

# Composite tumour categories (e.g. T3;T4) may correspond to >1 WGS library.
# The official source table then contains one positive row per positive library,
# but collapses the library name. We expand ONLY if the positive-row multiplicity
# equals the number of validated target libraries; otherwise occupancy is
# unidentifiable and the script aborts rather than guessing.
tert_positive_groups <- tert_filtered[, .(
  n_source_positive_rows = .N,
  source_sample_ids = paste(sort(unique(SAMPLE_ID)), collapse = ";"),
  min_alt_count = min(t_alt_count),
  min_total_depth = min(total_depth),
  min_vaf = min(VAF)
), by = .(
  source_patient_id = PATIENT_ID,
  category,
  hotspot_label,
  Start_Position,
  Reference_Allele,
  Tumor_Seq_Allele2
)]

crosswalk_group_sizes <- source_crosswalk[, .(
  n_target_libraries = .N
), by = .(source_patient_id, category)]

tert_mapping_qc <- merge(
  tert_positive_groups,
  crosswalk_group_sizes,
  by = c("source_patient_id", "category"),
  all.x = TRUE,
  sort = FALSE
)
tert_mapping_qc[, mapping_resolved :=
                  !is.na(n_target_libraries) &
                  n_source_positive_rows == n_target_libraries]

if (DEBUG) {
  fwrite(
    tert_mapping_qc,
    file.path(ext_qc_dir, "15_TERT_promoter_source_to_sector_mapping_qc.tsv"),
    sep = "\t"
  )
}

if (any(!tert_mapping_qc$mapping_resolved)) {
  stop(
    "At least one positive TERT patient/category cannot be mapped unambiguously ",
    "to validated WGS sectors. See 15_TERT_promoter_source_to_sector_mapping_qc.tsv."
  )
}

# Because every positive source group has multiplicity equal to its number of
# target WGS libraries, ALL target libraries in that group are positive.
tert_presence <- merge(
  tert_positive_groups[, .(
    source_patient_id, category, hotspot_label,
    POS = as.integer(Start_Position),
    REF_ALLELE = as.character(Reference_Allele),
    ALT_ALLELE = as.character(Tumor_Seq_Allele2)
  )],
  source_crosswalk[, .(
    source_patient_id, category, patient_id, sample_id, n_sectors
  )],
  by = c("source_patient_id", "category"),
  all.x = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)

if (anyNA(tert_presence$sample_id) || anyNA(tert_presence$patient_id)) {
  stop("TERT source->validated sector expansion produced missing mappings.")
}
if (anyDuplicated(tert_presence[, .(sample_id, hotspot_label)]) > 0L) {
  stop("Duplicate validated sample x TERT-hotspot presence after mapping.")
}
if (nrow(tert_presence) != EXPECTED_TERT_SOURCE_ROWS_FILTERED) {
  stop(
    "TERT sector-presence count after mapping should equal filtered source rows (",
    EXPECTED_TERT_SOURCE_ROWS_FILTERED, "); observed ", nrow(tert_presence), "."
  )
}

tert_reference <- tert_presence[, .(
  full_count = uniqueN(sample_id),
  n_sector_values = uniqueN(n_sectors),
  n_sectors = as.integer(n_sectors[[1L]])
), by = .(patient_id, hotspot_label, POS, REF_ALLELE, ALT_ALLELE)]

if (any(tert_reference$n_sector_values != 1L)) {
  stop("Inconsistent validated n_sectors within at least one TERT reference event.")
}
tert_reference[, n_sector_values := NULL]

tert_reference[, `:=`(
  event_key = paste0(
    "TERT_PROMOTER|", patient_id, "|chr5:", POS, ":", REF_ALLELE, ">", ALT_ALLELE
  ),
  driver_gene = "TERT",
  contains_DICER1 = FALSE,
  source = "official_snv_indel_TERT_hotspot",
  is_ubiquitous = full_count == n_sectors,
  is_nonubiquitous = full_count > 0L & full_count < n_sectors,
  is_private = full_count == 1L
)]

if (uniqueN(tert_reference$patient_id) != EXPECTED_TERT_PATIENTS) {
  stop("Mapped TERT reference no longer contains 42 validated patients.")
}
if (any(tert_reference$full_count < 1L | tert_reference$full_count > tert_reference$n_sectors)) {
  stop("Invalid TERT full_count relative to patient sector count.")
}

if (DEBUG) {
  fwrite(
    tert_presence,
    file.path(ext_qc_dir, "15_TERT_promoter_validated_sector_presence.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    tert_reference,
    file.path(ext_qc_dir, "15_TERT_promoter_reference_events.tsv"),
    sep = "\t"
  )
}

# -----------------------------------------------------------------------------
# 3. Freeze final driver event universe
# -----------------------------------------------------------------------------
# Canonical promoter hotspots are supplied authoritatively by source_data and
# are removed from the Script-14 component if they happen to appear there, to
# prevent duplicate events or mixed calling rules.
driver_candidates[, chrom_norm := sub(
  "^chr", "", as.character(CHROM), ignore.case = TRUE
)]
driver_candidates[, is_canonical_tert :=
  grepl("(^|;)TERT(;|$)", as.character(driver_genes)) &
    chrom_norm == "5" &
    as.integer(POS) %in% TERT_HOTSPOTS$POS_SOURCE &
    as.character(REF_ALLELE) == "G" &
    as.character(ALT_ALLELE) == "A"
]

if (driver_candidates[is_canonical_tert == TRUE, .N] != 3L ||
    driver_candidates[is_canonical_tert == TRUE & any_high_moderate == TRUE, .N] != 0L) {
  stop("Script-14 canonical-TERT regression failed: expected 3 non-HIGH/MODERATE C250T events only.")
}

# Strong crosswalk regression: C250T is the canonical promoter hotspot already
# visible in the Script-14 TERT annotation (3 patient-specific events). Its
# source-data reconstruction must reproduce the SAME validated patients and
# occupancies exactly. This independently checks the ITH->WGS sector mapping.
script14_c250 <- driver_candidates[
  is_canonical_tert == TRUE & as.integer(POS) == 1295135L,
  .(
    patient_id = as.character(patient_id),
    script14_full_count = as.integer(full_count),
    script14_n_sectors = as.integer(n_sectors),
    script14_nonubiquitous = as.logical(is_nonubiquitous),
    script14_private = as.logical(is_private)
  )
]
source_c250 <- tert_reference[
  hotspot_label == "TERT_C250T",
  .(
    patient_id = as.character(patient_id),
    source_full_count = as.integer(full_count),
    source_n_sectors = as.integer(n_sectors),
    source_nonubiquitous = as.logical(is_nonubiquitous),
    source_private = as.logical(is_private)
  )
]
tert_c250_regression <- merge(
  script14_c250, source_c250, by = "patient_id", all = TRUE, sort = TRUE
)
tert_c250_regression[, pass :=
  !is.na(script14_full_count) & !is.na(source_full_count) &
    script14_full_count == source_full_count &
    script14_n_sectors == source_n_sectors &
    script14_nonubiquitous == source_nonubiquitous &
    script14_private == source_private
]
if (DEBUG) {
  fwrite(
    tert_c250_regression,
    file.path(ext_qc_dir, "15_TERT_C250T_source_vs_script14_regression.tsv"),
    sep = "\t"
  )
}
if (nrow(tert_c250_regression) != 3L || any(!tert_c250_regression$pass)) {
  stop(
    "TERT C250T source-data reconstruction does not exactly reproduce the 3 ",
    "Script-14 patient occupancies. See 15_TERT_C250T_source_vs_script14_regression.tsv."
  )
}

base_driver <- driver_candidates[
  any_high_moderate == TRUE & is_canonical_tert == FALSE
]

base_driver_reference <- base_driver[, .(
  event_key = as.character(event_key),
  patient_id = as.character(patient_id),
  driver_gene = as.character(driver_genes),
  contains_DICER1 = as.logical(contains_DICER1),
  source = "validated_script02_HIGH_MODERATE_driver_gene",
  full_count = as.integer(full_count),
  n_sectors = as.integer(n_sectors),
  is_ubiquitous = as.logical(is_ubiquitous),
  is_nonubiquitous = as.logical(is_nonubiquitous),
  is_private = as.logical(is_private)
)]

if (anyDuplicated(base_driver_reference$event_key) > 0L) {
  stop("Duplicate event_key in HIGH/MODERATE driver component.")
}

driver49_reference <- rbindlist(
  list(
    base_driver_reference,
    tert_reference[, .(
      event_key, patient_id, driver_gene, contains_DICER1, source,
      full_count, n_sectors, is_ubiquitous, is_nonubiquitous, is_private
    )]
  ),
  use.names = TRUE,
  fill = TRUE
)
if (anyDuplicated(driver49_reference$event_key) > 0L) {
  stop("Duplicate event_key in final 48+DICER1 driver universe.")
}

# Sensitivity matching the filename's original 48-gene set (DICER1 removed).
driver48_reference <- driver49_reference[contains_DICER1 == FALSE]

# Reference integrity.
for (nm in c("driver49_reference", "driver48_reference")) {
  rr <- get(nm)
  if (any(rr$full_count < 1L | rr$full_count > rr$n_sectors)) {
    stop(nm, ": invalid full_count/n_sectors.")
  }
  if (any(rr$is_nonubiquitous != (rr$full_count < rr$n_sectors))) {
    stop(nm, ": non-ubiquitous classification inconsistency.")
  }
  if (any(rr$is_private != (rr$full_count == 1L))) {
    stop(nm, ": private classification inconsistency.")
  }
}

if (DEBUG) {
  fwrite(
    driver49_reference,
    file.path(ext_dir, "15_driver49_final_reference_events.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    driver48_reference,
    file.path(ext_qc_dir, "15_driver48_no_DICER1_reference_events.tsv"),
    sep = "\t"
  )
}

# -----------------------------------------------------------------------------
# 4. Exact patient-level expected recovery curves
# -----------------------------------------------------------------------------
# Only the NEW driver curves are computed here. Comparator curves are taken
# directly from the already validated exhaustive Script 04 / Script 08 patient
# summaries, avoiding any redundant recomputation of the ~5.6M-event reference.
setindexv(driver49_reference, "patient_id")
setindexv(driver48_reference, "patient_id")

make_patient_expectations <- function(ref, event_set_label) {
  ref <- as.data.table(ref)
  require_columns(ref, c("patient_id", "full_count", "n_sectors", "is_nonubiquitous"), event_set_label)

  out <- vector("list", nrow(fixed_patients))
  for (ii in seq_len(nrow(fixed_patients))) {
    pid <- fixed_patients$patient_id[[ii]]
    n <- fixed_patients$n[[ii]]
    rr <- ref[.(pid), on = .(patient_id), nomatch = 0L]
    nonubiq_m <- as.integer(rr[is_nonubiquitous == TRUE, full_count])

    rows <- vector("list", length(PRIMARY_K))
    for (jj in seq_along(PRIMARY_K)) {
      k <- PRIMARY_K[[jj]]
      rec <- expected_recovery_from_occupancies(nonubiq_m, n, k)
      rows[[jj]] <- data.table(
        event_set = event_set_label,
        patient_id = pid,
        n = n,
        k = k,
        n_reference_events = nrow(rr),
        n_nonubiquitous_reference = length(nonubiq_m),
        expected_detection_recall = unname(rec[["detection"]]),
        expected_classification_recall = unname(rec[["classification"]])
      )
    }
    out[[ii]] <- rbindlist(rows)
  }
  rbindlist(out)
}

make_validated_comparator <- function(x, event_set_label) {
  x <- as.data.table(x)
  z <- x[
    patient_id %in% fixed_patients$patient_id &
      n >= PRIMARY_MIN_N &
      k %in% PRIMARY_K,
    .(
      event_set = event_set_label,
      patient_id = as.character(patient_id),
      n = as.integer(n),
      k = as.integer(k),
      n_reference_events = NA_integer_,
      n_nonubiquitous_reference = NA_integer_,
      expected_detection_recall = as.numeric(recall_nonubiquitous_detection_mean),
      expected_classification_recall = as.numeric(recall_heterogeneity_classification_mean)
    )
  ]
  setorder(z, patient_id, k)

  if (nrow(z) != EXPECTED_FIXED_N_GE5 * length(PRIMARY_K) ||
      uniqueN(z$patient_id) != EXPECTED_FIXED_N_GE5 ||
      anyDuplicated(z[, .(patient_id, k)]) > 0L ||
      any(z$k >= z$n)) {
    stop(event_set_label, ": validated comparator fixed-cohort regression failed.")
  }
  z
}

patient_expectations <- rbindlist(
  list(
    make_patient_expectations(driver49_reference, "driver_48_plus_DICER1"),
    make_patient_expectations(driver48_reference, "driver_original_48"),
    make_validated_comparator(protein_patient, "protein_altering_all_genes"),
    make_validated_comparator(all_filtered_patient, "all_filtered_exact_variants")
  ),
  use.names = TRUE,
  fill = TRUE
)
setorder(patient_expectations, event_set, patient_id, k)

# Normalize only floating-point boundary noise in both newly computed driver
# expectations and imported validated comparator expectations.
patient_expectations[, expected_detection_recall := canonicalize_probability(
  expected_detection_recall, label = "expected_detection_recall"
)]
patient_expectations[, expected_classification_recall := canonicalize_probability(
  expected_classification_recall, label = "expected_classification_recall"
)]

# Classification at k=1 must be exactly 0 whenever defined.
if (any(
  patient_expectations[
    k == 1L & is.finite(expected_classification_recall),
    abs(expected_classification_recall)
  ] > 1e-12
)) {
  stop("Heterogeneity-classification recall at k=1 is not zero.")
}
if (any(
  patient_expectations[
    k == 1L & is.finite(expected_classification_recall),
    expected_classification_recall
  ] != 0
)) {
  stop("Numerical canonicalization failed: k=1 classification is not exactly zero.")
}
if (any(
  patient_expectations[
    is.finite(expected_detection_recall) & is.finite(expected_classification_recall),
    expected_classification_recall - expected_detection_recall
  ] > 1e-12
)) {
  stop("Classification recall exceeds detection recall for at least one patient/k.")
}

if (DEBUG) {
  fwrite(
    patient_expectations,
    file.path(ext_dir, "15_driver_patient_expected_recovery_fixed_n_ge5.tsv"),
    sep = "\t"
  )
}

# -----------------------------------------------------------------------------
# 5. Cohort summaries with exact-bootstrap CIs
# -----------------------------------------------------------------------------
metric_map <- c(
  detection = "expected_detection_recall",
  heterogeneity_classification = "expected_classification_recall"
)

summary_rows <- list()
idx <- 0L
for (es in unique(patient_expectations$event_set)) {
  for (metric_id in names(metric_map)) {
    col <- metric_map[[metric_id]]
    for (kval in PRIMARY_K) {
      vals <- patient_expectations[event_set == es & k == kval, get(col)]
      vals <- as.numeric(vals[is.finite(vals)])
      ci <- exact_bootstrap_median_ci(vals, level = CI_LEVEL, min_n = 1L)
      idx <- idx + 1L
      summary_rows[[idx]] <- data.table(
        analysis = "fixed_n_ge5_true_downsampling",
        event_set = es,
        metric = metric_id,
        k = as.integer(kval),
        n_fixed_cohort = EXPECTED_FIXED_N_GE5,
        n_patients_with_defined_endpoint = length(vals),
        fraction_fixed_cohort_with_defined_endpoint = length(vals) / EXPECTED_FIXED_N_GE5,
        median = safe_median(vals),
        ci_lower = unname(ci[["lower"]]),
        ci_upper = unname(ci[["upper"]]),
        ci_level = CI_LEVEL,
        bootstrap_method = if (length(vals) > 0L) "exact_nonparametric_percentile" else NA_character_
      )
    }
  }
}

depth_summary <- rbindlist(summary_rows)
setorder(depth_summary, event_set, metric, k)
fwrite(
  depth_summary,
  file.path(ext_dir, "15_driver_depth_summary_fixed_n_ge5.tsv"),
  sep = "\t"
)

# -----------------------------------------------------------------------------
# 6. Paired comparisons: driver vs current mutation universes
# -----------------------------------------------------------------------------
make_paired_comparison <- function(driver_label, comparator_label, metric_id, col) {
  a <- patient_expectations[event_set == driver_label, .(
    patient_id, k, driver_value = get(col)
  )]
  b <- patient_expectations[event_set == comparator_label, .(
    patient_id, k, comparator_value = get(col)
  )]
  z <- merge(a, b, by = c("patient_id", "k"), all = FALSE)
  z <- z[is.finite(driver_value) & is.finite(comparator_value)]
  z[, difference := canonicalize_zero(driver_value - comparator_value)]

  rows <- vector("list", length(PRIMARY_K))
  for (jj in seq_along(PRIMARY_K)) {
    kval <- PRIMARY_K[[jj]]
    d <- z[k == kval]
    vals <- d$difference
    ci <- exact_bootstrap_median_ci(vals, level = CI_LEVEL, min_n = 1L)
    rows[[jj]] <- data.table(
      analysis = "paired_fixed_n_ge5_true_downsampling",
      driver_event_set = driver_label,
      comparator_event_set = comparator_label,
      metric = metric_id,
      k = kval,
      n_paired_patients = length(vals),
      median_driver_minus_comparator = safe_median(vals),
      ci_lower = unname(ci[["lower"]]),
      ci_upper = unname(ci[["upper"]]),
      proportion_driver_greater = if (length(vals)) mean(vals > 0) else NA_real_,
      proportion_equal = if (length(vals)) mean(abs(vals) <= 1e-12) else NA_real_,
      bootstrap_method = if (length(vals)) "exact_nonparametric_percentile" else NA_character_
    )
  }
  rbindlist(rows)
}

paired_comparisons <- rbindlist(
  lapply(
    names(metric_map),
    function(metric_id) {
      col <- metric_map[[metric_id]]
      rbindlist(list(
        make_paired_comparison(
          "driver_48_plus_DICER1", "protein_altering_all_genes", metric_id, col
        ),
        make_paired_comparison(
          "driver_48_plus_DICER1", "all_filtered_exact_variants", metric_id, col
        ),
        make_paired_comparison(
          "driver_48_plus_DICER1", "driver_original_48", metric_id, col
        )
      ))
    }
  )
)
setorder(paired_comparisons, comparator_event_set, metric, k)
fwrite(
  paired_comparisons,
  file.path(ext_dir, "15_driver_paired_comparisons_fixed_n_ge5.tsv"),
  sep = "\t"
)

# -----------------------------------------------------------------------------
# 7. Minimal heuristic + event-universe QC
# -----------------------------------------------------------------------------
primary_detection <- depth_summary[
  event_set == "driver_48_plus_DICER1" & metric == "detection"
]
setorder(primary_detection, k)

if (nrow(primary_detection) != length(PRIMARY_K) ||
    any(primary_detection$n_patients_with_defined_endpoint < 1L)) {
  stop("Primary driver detection endpoint is undefined or incomplete in the fixed cohort.")
}
if (uniqueN(primary_detection$n_patients_with_defined_endpoint) != 1L) {
  stop("Driver endpoint-defined patient N changes across k in the fixed cohort.")
}

over_threshold <- primary_detection[
  is.finite(median) & median >= PARAMS$adequate_recovery_threshold
]
first_k_ge_80 <- if (nrow(over_threshold) == 0L) NA_integer_ else min(over_threshold$k)

heuristic <- data.table(
  event_set = "driver_48_plus_DICER1",
  cohort = "fixed_n_ge5_true_downsampling",
  threshold = PARAMS$adequate_recovery_threshold,
  first_tested_k_meeting_threshold = first_k_ge_80,
  interpretation = if (is.na(first_k_ge_80)) {
    "No tested k=1..4 reaches the legacy 80% median-detection benchmark."
  } else {
    paste0(
      "k=", first_k_ge_80,
      " is the first TESTED depth reaching the legacy 80% median-detection benchmark; not a universal biological optimum."
    )
  }
)
if (DEBUG) {
  fwrite(
    heuristic,
    file.path(ext_dir, "15_driver_sampling_operational_heuristic.tsv"),
    sep = "\t"
  )
}

universe_qc <- rbindlist(list(
  data.table(
    event_set = "driver_48_plus_DICER1",
    n_exact_patient_events = nrow(driver49_reference),
    n_patients_with_any_event = uniqueN(driver49_reference$patient_id),
    n_nonubiquitous_events = sum(driver49_reference$is_nonubiquitous),
    n_patients_with_nonubiquitous_event = uniqueN(
      driver49_reference[is_nonubiquitous == TRUE, patient_id]
    ),
    n_private_events = sum(driver49_reference$is_private)
  ),
  data.table(
    event_set = "driver_original_48",
    n_exact_patient_events = nrow(driver48_reference),
    n_patients_with_any_event = uniqueN(driver48_reference$patient_id),
    n_nonubiquitous_events = sum(driver48_reference$is_nonubiquitous),
    n_patients_with_nonubiquitous_event = uniqueN(
      driver48_reference[is_nonubiquitous == TRUE, patient_id]
    ),
    n_private_events = sum(driver48_reference$is_private)
  ),
  data.table(
    event_set = "TERT_promoter_hotspots_only",
    n_exact_patient_events = nrow(tert_reference),
    n_patients_with_any_event = uniqueN(tert_reference$patient_id),
    n_nonubiquitous_events = sum(tert_reference$is_nonubiquitous),
    n_patients_with_nonubiquitous_event = uniqueN(
      tert_reference[is_nonubiquitous == TRUE, patient_id]
    ),
    n_private_events = sum(tert_reference$is_private)
  )
))
if (DEBUG) {
  fwrite(
    universe_qc,
    file.path(ext_qc_dir, "15_driver_event_universe_qc.tsv"),
    sep = "\t"
  )
}

# -----------------------------------------------------------------------------
# 8. Human-readable summary and checkpoint
# -----------------------------------------------------------------------------
d4 <- depth_summary[
  event_set == "driver_48_plus_DICER1" & metric == "detection" & k == 4L
]
h4 <- depth_summary[
  event_set == "driver_48_plus_DICER1" &
    metric == "heterogeneity_classification" & k == 4L
]
p4 <- paired_comparisons[
  comparator_event_set == "protein_altering_all_genes" &
    metric == "detection" & k == 4L
]

if (nrow(d4) != 1L || nrow(h4) != 1L || nrow(p4) != 1L ||
    d4$n_patients_with_defined_endpoint < 1L || p4$n_paired_patients < 1L) {
  stop("Final driver summary rows are missing or undefined.")
}

summary_lines <- c(
  "HCC extensions - Script 15 driver sampling",
  "================================================",
  "",
  paste0("snv_indel.tsv MD5: ", snv_md5, " [PASS]"),
  paste0("cnv_arm_level.tsv crosswalk MD5: ", cnv_map_md5, " [PASS]"),
  "Source-data ITH -> validated WGS crosswalk: 123 patients / 490 sectors [PASS]",
  sprintf(
    "Canonical TERT hotspots before/after validated numeric filters: %d / %d rows; %d patients [PASS]",
    nrow(tert_raw), nrow(tert_filtered), uniqueN(tert_filtered$PATIENT_ID)
  ),
  sprintf(
    "TERT C228T/C250T patients: %d / %d [PASS]",
    tert_patient_counts[hotspot_label == "TERT_C228T", n_patients],
    tert_patient_counts[hotspot_label == "TERT_C250T", n_patients]
  ),
  "TERT positive source-category multiplicity -> validated sectors: unambiguous [PASS]",
  "TERT C250T source-data occupancy -> Script-14 occupancy regression: 3/3 [PASS]",
  "",
  sprintf(
    "Final 48+DICER1 driver universe: %d patient-specific exact events; %d non-ubiquitous",
    nrow(driver49_reference), sum(driver49_reference$is_nonubiquitous)
  ),
  sprintf(
    "Patients with >=1 non-ubiquitous driver event (all n): %d",
    uniqueN(driver49_reference[is_nonubiquitous == TRUE, patient_id])
  ),
  sprintf(
    "Primary fixed n>=5 cohort: 55 total; %d with defined non-ubiquitous-driver endpoint",
    d4$n_patients_with_defined_endpoint
  ),
  "",
  sprintf(
    "k=4 driver detection median = %.6f [%.6f, %.6f]",
    d4$median, d4$ci_lower, d4$ci_upper
  ),
  sprintf(
    "k=4 driver heterogeneity-classification median = %.6f [%.6f, %.6f]",
    h4$median, h4$ci_lower, h4$ci_upper
  ),
  sprintf(
    "k=4 paired driver - protein detection difference = %.6f [%.6f, %.6f] (N=%d)",
    p4$median_driver_minus_comparator, p4$ci_lower, p4$ci_upper, p4$n_paired_patients
  ),
  paste0("Operational 80% heuristic: ", heuristic$interpretation),
  "",
  "Interpretation guardrail: the 80% line is a legacy operational benchmark, not a biological optimum.",
  "Sampling expectation is exact/hypergeometric; cohort CI is exact-bootstrap; no Monte-Carlo sampling was used.",
  "",
  sprintf("Elapsed time: %.1f seconds", proc.time()[["elapsed"]] - t0)
)

writeLines(
  summary_lines,
  file.path(ext_qc_dir, "15_driver_sampling_summary.txt")
)



cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
