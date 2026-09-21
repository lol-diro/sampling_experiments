# ==============================================================================
# 01_hcc_metadata_qc.R
#
# Purpose
#   1. Reconstruct the tumour-sector analysis population without silently
#      discarding metadata anomalies.
#   2. Parse the spatial labels conservatively.
#   3. Define a strict primary spatial cohort and a broader order-only
#      sensitivity cohort.
#   4. Audit the clinical-data schema without assuming which columns encode
#      tumour size/stage.
#
# This script does NOT filter mutations and does NOT perform downsampling.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# Locate 00_config.R next to this script when run with Rscript; fall back to the
# current working directory when sourced interactively.
args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
if (length(file_arg) > 0L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  script_dir <- getwd()
}
source(file.path(script_dir, "00_config.R"))

# ------------------------------- helpers -------------------------------------
load_expected_object <- function(path, expected_name) {
  e <- new.env(parent = baseenv())
  loaded <- load(path, envir = e)
  if (!(expected_name %in% loaded)) {
    stop(
      "Expected object '", expected_name, "' not found in ", path,
      ". Objects present: ", paste(loaded, collapse = ", ")
    )
  }
  get(expected_name, envir = e, inherits = FALSE)
}

require_columns <- function(x, cols, object_name) {
  miss <- setdiff(cols, names(x))
  if (length(miss) > 0L) {
    stop(object_name, " is missing required column(s): ", paste(miss, collapse = ", "))
  }
}

extract_integer_tokens <- function(x, token) {
  x <- ifelse(is.na(x), "", as.character(x))
  # Examples handled:
  # T1, T_1, MF1T2, MF2_T1, T4;T5
  pattern <- paste0(token, "\\s*_?\\s*[0-9]+")
  hits <- regmatches(x, gregexpr(pattern, x, perl = TRUE, ignore.case = TRUE))
  lapply(hits, function(z) {
    if (length(z) == 0L) return(integer(0))
    as.integer(gsub("[^0-9]", "", z))
  })
}

collapse_ints <- function(z) {
  if (length(z) == 0L) return(NA_character_)
  paste(z, collapse = ",")
}

single_int_or_na <- function(z) {
  if (length(z) != 1L) return(NA_integer_)
  as.integer(z[[1]])
}

make_primary_reason <- function(any_mf, any_missing_t, any_composite_t,
                                duplicate_position, consecutive_positions) {
  reasons <- character(0)
  if (isTRUE(any_mf)) reasons <- c(reasons, "MF-labelled patient")
  if (isTRUE(any_missing_t)) reasons <- c(reasons, "missing T coordinate")
  if (isTRUE(any_composite_t)) reasons <- c(reasons, "composite/multiple T coordinates")
  if (isTRUE(duplicate_position)) reasons <- c(reasons, "duplicate T coordinate")
  if (!isTRUE(consecutive_positions) &&
      !isTRUE(any_missing_t) && !isTRUE(any_composite_t) && !isTRUE(duplicate_position)) {
    reasons <- c(reasons, "non-consecutive T coordinates")
  }
  if (length(reasons) == 0L) "eligible" else paste(reasons, collapse = "; ")
}

make_order_reason <- function(any_mf, any_missing_t, any_composite_t,
                              duplicate_position) {
  reasons <- character(0)
  if (isTRUE(any_mf)) reasons <- c(reasons, "MF-labelled patient")
  if (isTRUE(any_missing_t)) reasons <- c(reasons, "missing T coordinate")
  if (isTRUE(any_composite_t)) reasons <- c(reasons, "composite/multiple T coordinates")
  if (isTRUE(duplicate_position)) reasons <- c(reasons, "duplicate T coordinate")
  if (length(reasons) == 0L) "eligible" else paste(reasons, collapse = "; ")
}

safe_examples <- function(x, n = 5L) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)
  ux <- unique(as.character(x))
  paste(utils::head(ux, n), collapse = " | ")
}

# ------------------------------- load map ------------------------------------
map_raw <- as.data.table(load_expected_object(INPUT_FILES$sample_map, "wgs_samples_matching"))
require_columns(
  map_raw,
  c("PATIENT_ID", "SAMPLE_ID", "SAMPLE_TYPE", "CATEGORY"),
  "wgs_samples_matching"
)

# Preserve exact source columns as character for deterministic comparisons.
map_raw[, PATIENT_ID := as.character(PATIENT_ID)]
map_raw[, SAMPLE_ID := as.character(SAMPLE_ID)]
map_raw[, SAMPLE_TYPE := as.character(SAMPLE_TYPE)]
map_raw[, CATEGORY := as.character(CATEGORY)]

sample_type_counts <- map_raw[, .N, by = SAMPLE_TYPE][order(-N, SAMPLE_TYPE)]
if (DEBUG) {
  fwrite(sample_type_counts, file.path(PATHS$qc_dir, "01_sample_type_counts.tsv"), sep = "\t")
}

# Exact duplicated rows can be removed, but conflicting duplicate sample IDs are
# not silently resolved.
map_unique <- unique(map_raw[, .(PATIENT_ID, SAMPLE_ID, SAMPLE_TYPE, CATEGORY)])
conflicting_sample_ids <- map_unique[, .N, by = SAMPLE_ID][N > 1L]
if (nrow(conflicting_sample_ids) > 0L) {
  if (DEBUG) {
    fwrite(
      map_unique[SAMPLE_ID %in% conflicting_sample_ids$SAMPLE_ID][order(SAMPLE_ID)],
      file.path(PATHS$qc_dir, "01_conflicting_sample_ids.tsv"), sep = "\t"
    )
  }
  stop(
    "At least one SAMPLE_ID maps to multiple metadata rows. See ",
    file.path(PATHS$qc_dir, "01_conflicting_sample_ids.tsv")
  )
}

# ----------------------- audit raw mutation sample IDs -----------------------
# The legacy analysis used somatic_mutations$PATIENT_ID as sample_id. We verify
# this again but immediately discard the large mutation object after extracting
# unique sample identifiers.
mut_obj <- load_expected_object(INPUT_FILES$mutations, "somatic_mutations")
mut_dt <- as.data.table(mut_obj)
require_columns(mut_dt, "PATIENT_ID", "somatic_mutations")
mut_sample_ids <- unique(as.character(mut_dt$PATIENT_ID))
rm(mut_obj, mut_dt)
invisible(gc())

map_sample_ids <- unique(map_unique$SAMPLE_ID)
mutation_to_map_match_rate <- mean(mut_sample_ids %in% map_sample_ids)

if (!is.finite(mutation_to_map_match_rate) || mutation_to_map_match_rate < 0.95) {
  stop(
    sprintf(
      "Only %.3f of unique somatic_mutations$PATIENT_ID values match metadata SAMPLE_ID. ",
      mutation_to_map_match_rate
    ),
    "Do not proceed until the identifier mapping is resolved."
  )
}

# ------------------------- tumour analysis population ------------------------
tumor_map <- map_unique[SAMPLE_TYPE == "Tumor"]
tumor_map[, has_raw_mutation_record := SAMPLE_ID %chin% mut_sample_ids]

if (DEBUG) {
  fwrite(
    tumor_map[has_raw_mutation_record == FALSE][order(PATIENT_ID, SAMPLE_ID)],
    file.path(PATHS$qc_dir, "01_tumor_samples_without_raw_mutation_record.tsv"),
    sep = "\t"
  )
}

# Keep the legacy analysis population, but only after explicitly reporting any
# tumour samples without a raw mutation record.
analysis_map <- tumor_map[has_raw_mutation_record == TRUE,
  .(
    patient_id = PATIENT_ID,
    sample_id = SAMPLE_ID,
    sample_type = SAMPLE_TYPE,
    category = CATEGORY
  )
]
setorder(analysis_map, patient_id, sample_id)

if (analysis_map[, anyDuplicated(sample_id)] > 0L) {
  stop("Duplicate sample_id remains in the tumour analysis map.")
}

# ------------------------------- parse labels --------------------------------
t_lists <- extract_integer_tokens(analysis_map$category, "T")
mf_lists <- extract_integer_tokens(analysis_map$category, "MF")

analysis_map[, n_t_labels := lengths(t_lists)]
analysis_map[, t_values := vapply(t_lists, collapse_ints, character(1))]
analysis_map[, sector_num := vapply(t_lists, single_int_or_na, integer(1))]
analysis_map[, n_mf_labels := lengths(mf_lists)]
analysis_map[, mf_values := vapply(mf_lists, collapse_ints, character(1))]
analysis_map[, focus_num := vapply(mf_lists, single_int_or_na, integer(1))]
analysis_map[, has_mf := n_mf_labels > 0L]
analysis_map[, category_has_semicolon := grepl(";", category, fixed = TRUE)]

# ----------------------------- patient-level QC ------------------------------
patient_qc <- analysis_map[, {
  n <- .N
  any_mf <- any(has_mf)
  any_missing_t <- any(n_t_labels == 0L)
  any_composite_t <- any(n_t_labels > 1L)
  all_single_t <- all(n_t_labels == 1L)

  positions <- sector_num
  duplicate_position <- if (all_single_t) anyDuplicated(positions) > 0L else FALSE

  ordered_positions <- if (all_single_t && !duplicate_position) sort(positions) else integer(0)
  consecutive_positions <- if (length(ordered_positions) <= 1L && all_single_t && !duplicate_position) {
    TRUE
  } else if (length(ordered_positions) >= 2L) {
    all(diff(ordered_positions) == 1L)
  } else {
    FALSE
  }

  unique_mf <- sort(unique(focus_num[!is.na(focus_num)]))

  # Primary spatial cohort is deliberately conservative:
  #   - no MF label at all
  #   - exactly one T coordinate per sample
  #   - no duplicate coordinate
  #   - consecutive sampled T coordinates
  spatial_primary_eligible <-
    !any_mf && all_single_t && !duplicate_position && consecutive_positions

  # Order-only sensitivity relaxes only the consecutiveness requirement. It is
  # useful later for gapped but otherwise unambiguous single-axis samples.
  spatial_order_eligible <-
    !any_mf && all_single_t && !duplicate_position

  list(
    n_sectors = n,
    any_mf = any_mf,
    n_distinct_mf = length(unique_mf),
    mf_values_patient = if (length(unique_mf) == 0L) NA_character_ else paste(unique_mf, collapse = ","),
    any_missing_t = any_missing_t,
    any_composite_t = any_composite_t,
    duplicate_sector_position = duplicate_position,
    t_positions = if (length(ordered_positions) == 0L) NA_character_ else paste(ordered_positions, collapse = ","),
    consecutive_positions = consecutive_positions,
    spatial_primary_eligible = spatial_primary_eligible,
    spatial_order_eligible = spatial_order_eligible,
    categories = paste(category, collapse = "|")
  )
}, by = patient_id]

patient_qc[, spatial_primary_reason := mapply(
  make_primary_reason,
  any_mf, any_missing_t, any_composite_t,
  duplicate_sector_position, consecutive_positions,
  USE.NAMES = FALSE
)]

patient_qc[, spatial_order_reason := mapply(
  make_order_reason,
  any_mf, any_missing_t, any_composite_t, duplicate_sector_position,
  USE.NAMES = FALSE
)]

setorder(patient_qc, patient_id)

# Join patient-level eligibility back to sample rows for later scripts.
analysis_map <- merge(
  analysis_map,
  patient_qc[, .(
    patient_id, n_sectors,
    spatial_primary_eligible, spatial_order_eligible,
    spatial_primary_reason, spatial_order_reason
  )],
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(analysis_map, patient_id, sector_num, sample_id)

# ------------------------------- summaries -----------------------------------
sector_dist <- patient_qc[, .(n_patients = .N), by = n_sectors][order(n_sectors)]

eligibility_by_n <- patient_qc[, .(
  n_patients = .N,
  n_primary_spatial = sum(spatial_primary_eligible),
  n_order_spatial = sum(spatial_order_eligible)
), by = n_sectors][order(n_sectors)]

problem_rows <- analysis_map[
  n_t_labels != 1L | has_mf == TRUE | spatial_primary_eligible == FALSE
][order(patient_id, sector_num, sample_id)]

# Always written (not gated by DEBUG): Script 15 (Figure 5 driver panel)
# reads this file back as a frozen patient/sector ID crosswalk, not just as
# a diagnostic.
fwrite(analysis_map, file.path(PATHS$qc_dir, "01_tumor_sample_map_audit.tsv"), sep = "\t")
if (DEBUG) {
  fwrite(patient_qc, file.path(PATHS$qc_dir, "01_patient_spatial_qc.tsv"), sep = "\t")
}
if (DEBUG) {
  fwrite(sector_dist, file.path(PATHS$qc_dir, "01_sector_count_distribution.tsv"), sep = "\t")
}
if (DEBUG) {
  fwrite(eligibility_by_n, file.path(PATHS$qc_dir, "01_spatial_eligibility_by_n.tsv"), sep = "\t")
}
if (DEBUG) {
  fwrite(problem_rows, file.path(PATHS$qc_dir, "01_spatial_problem_rows.tsv"), sep = "\t")
}

saveRDS(analysis_map, file.path(PATHS$intermediate_dir, "hcc_sample_map.rds"))
saveRDS(patient_qc, file.path(PATHS$intermediate_dir, "hcc_patient_metadata_qc.rds"))

# --------------------------- clinical schema audit ----------------------------
clinical_obj <- load_expected_object(INPUT_FILES$clinical, "clinical_data")
clinical_dt <- as.data.table(clinical_obj)

clinical_schema <- data.table(
  column = names(clinical_dt),
  class = vapply(clinical_dt, function(x) paste(class(x), collapse = "/"), character(1)),
  n_nonmissing = vapply(clinical_dt, function(x) sum(!is.na(x)), integer(1)),
  n_unique_nonmissing = vapply(clinical_dt, function(x) uniqueN(x[!is.na(x)]), integer(1)),
  examples = vapply(clinical_dt, safe_examples, character(1))
)

# These are only candidate columns identified by name. No biological role is
# assigned automatically; we inspect this output before Script 02.
clinical_schema[, candidate_role := ""]
clinical_schema[grepl("PATIENT|SUBJECT|CASE|ID$", column, ignore.case = TRUE),
                candidate_role := "possible_patient_identifier"]
clinical_schema[grepl("SIZE|DIAM", column, ignore.case = TRUE),
                candidate_role := ifelse(
                  candidate_role == "", "possible_tumour_size_or_diameter",
                  paste(candidate_role, "possible_tumour_size_or_diameter", sep = ";")
                )]
clinical_schema[grepl("STAGE|TNM", column, ignore.case = TRUE),
                candidate_role := ifelse(
                  candidate_role == "", "possible_stage",
                  paste(candidate_role, "possible_stage", sep = ";")
                )]

if (DEBUG) {
  fwrite(clinical_schema, file.path(PATHS$qc_dir, "01_clinical_schema.tsv"), sep = "\t")
}
rm(clinical_obj, clinical_dt)
invisible(gc())

# ----------------------- compare with legacy expectations ---------------------
observed_n_patients <- uniqueN(analysis_map$patient_id)
observed_n_sectors <- nrow(analysis_map)

legacy_total_pass <-
  observed_n_patients == EXPECTED_LEGACY$n_patients &&
  observed_n_sectors == EXPECTED_LEGACY$n_tumor_sectors

obs_dist <- setNames(sector_dist$n_patients, as.character(sector_dist$n_sectors))
exp_dist <- EXPECTED_LEGACY$sector_count_distribution
all_levels <- union(names(exp_dist), names(obs_dist))
obs_full <- setNames(rep(0L, length(all_levels)), all_levels)
exp_full <- obs_full
obs_full[names(obs_dist)] <- as.integer(obs_dist)
exp_full[names(exp_dist)] <- as.integer(exp_dist)
legacy_dist_pass <- identical(as.integer(obs_full), as.integer(exp_full))

mutation_match_summary <- data.table(
  metric = c(
    "unique_mutation_sample_ids",
    "unique_metadata_sample_ids",
    "mutation_to_metadata_match_rate",
    "tumour_samples_in_metadata",
    "tumour_samples_with_raw_mutation_record",
    "tumour_samples_without_raw_mutation_record"
  ),
  value = c(
    length(mut_sample_ids),
    length(map_sample_ids),
    mutation_to_map_match_rate,
    nrow(tumor_map),
    nrow(analysis_map),
    sum(!tumor_map$has_raw_mutation_record)
  )
)
if (DEBUG) {
  fwrite(mutation_match_summary, file.path(PATHS$qc_dir, "01_mutation_sample_match_summary.tsv"), sep = "\t")
}

# ------------------------------- report --------------------------------------
summary_lines <- c(
  "Script 01: HCC metadata QC",
  "=====================================",
  "",
  sprintf("Analysis population: %d patients, %d tumour sectors", observed_n_patients, observed_n_sectors),
  sprintf("Legacy total-count QC: %s", ifelse(legacy_total_pass, "PASS", "FAIL")),
  sprintf("Legacy sector-distribution QC: %s", ifelse(legacy_dist_pass, "PASS", "FAIL")),
  sprintf("Mutation sample-ID -> metadata SAMPLE_ID match rate: %.4f", mutation_to_map_match_rate),
  sprintf("Tumour samples lacking any raw mutation record: %d", sum(!tumor_map$has_raw_mutation_record)),
  "",
  "Sector-count distribution:",
  paste(sprintf("  n=%s: %s patients", sector_dist$n_sectors, sector_dist$n_patients), collapse = "\n"),
  "",
  sprintf(
    "Strict primary spatial cohort: %d/%d patients (%d with n>=2)",
    sum(patient_qc$spatial_primary_eligible), nrow(patient_qc),
    sum(patient_qc$spatial_primary_eligible & patient_qc$n_sectors >= 2L)
  ),
  paste(
    "  Definition: no MF label; exactly one T coordinate per sample;",
    "no duplicate T coordinate; consecutive T coordinates."
  ),
  sprintf(
    "Order-only spatial sensitivity cohort: %d/%d patients (%d with n>=2)",
    sum(patient_qc$spatial_order_eligible), nrow(patient_qc),
    sum(patient_qc$spatial_order_eligible & patient_qc$n_sectors >= 2L)
  ),
  "  This sensitivity allows clean gaps in T numbering but still excludes MF/composite/duplicate coordinates.",
  "",
  sprintf("Patients with any MF label: %d", sum(patient_qc$any_mf)),
  sprintf("Patients with composite/multiple T labels: %d", sum(patient_qc$any_composite_t)),
  sprintf(
    "Otherwise clean no-MF patients with non-consecutive T positions: %d",
    sum(patient_qc$spatial_order_eligible & !patient_qc$consecutive_positions)
  ),
  "",
  "Spatial eligibility by original n:",
  paste(
    sprintf(
      "  n=%s: total=%s, primary=%s, order_sensitivity=%s",
      eligibility_by_n$n_sectors,
      eligibility_by_n$n_patients,
      eligibility_by_n$n_primary_spatial,
      eligibility_by_n$n_order_spatial
    ),
    collapse = "\n"
  ),
  "",
  "Clinical columns were NOT interpreted automatically.",
  "Inspect 01_clinical_schema.tsv before defining tumour-size/stage QC.",
  "",
  "Do not proceed to mutation filtering until these QC outputs are reviewed."
)

writeLines(summary_lines, file.path(PATHS$qc_dir, "01_qc_summary.txt"))


cat(paste(summary_lines, collapse = "\n"), "\n")

if (!legacy_total_pass || !legacy_dist_pass) {
  warning(
    "The reconstructed population does not reproduce the legacy QC targets. ",
    "Review Script 01 outputs before proceeding."
  )
}
