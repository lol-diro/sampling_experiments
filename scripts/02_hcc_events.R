# ==============================================================================
# 02_hcc_events.R
#
# Purpose
#   1. Reproduce the legacy primary SNV/indel filtering exactly.
#   2. Audit mutation annotations needed for the planned protein-altering
#      sensitivity analysis.
#   3. Collapse filtered calls to one presence/absence record per
#      patient x sector x exact variant.
#   4. Construct the full available-sector reference for every patient.
#
# This script does NOT perform downsampling.
# It requires Script 01 to have completed successfully.
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

# ---------------------------- frozen definitions -----------------------------
# Legacy QC target taken from the original analysis report. This is used only
# to verify exact reproduction of the previous primary filtering.
LEGACY_EXPECTED_FILTERED_ROWS <- 10612828L

# Sensitivity definition. The source data use the IMPACT field already present
# in somatic_mutations. We do not infer protein effects from gene names or from
# external databases here.
PROTEIN_ALTERING_IMPACTS <- c("HIGH", "MODERATE")
KNOWN_IMPACT_CATEGORIES <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

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

safe_char <- function(x) {
  if (is.factor(x)) as.character(x) else x
}

fmt_int <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

summarise_reference_set <- function(reference_events, patient_qc, event_set_name) {
  z <- reference_events[, .(
    n_events_reference = .N,
    n_ubiquitous_reference = sum(is_ubiquitous),
    n_nonubiquitous_reference = sum(is_nonubiquitous),
    n_private_reference = sum(is_private),
    mean_occupancy_fraction = mean(occupancy_fraction),
    median_occupancy_fraction = stats::median(occupancy_fraction)
  ), by = patient_id]

  out <- merge(
    patient_qc[, .(
      patient_id, n_sectors,
      spatial_primary_eligible, spatial_order_eligible
    )],
    z,
    by = "patient_id",
    all.x = TRUE,
    sort = FALSE
  )

  count_cols <- c(
    "n_events_reference", "n_ubiquitous_reference",
    "n_nonubiquitous_reference", "n_private_reference"
  )
  for (cc in count_cols) {
    set(out, which(is.na(out[[cc]])), cc, 0L)
  }

  out[, nonubiquitous_fraction_reference := fifelse(
    n_events_reference > 0L,
    n_nonubiquitous_reference / n_events_reference,
    NA_real_
  )]
  out[, event_set := event_set_name]
  setcolorder(out, c(
    "event_set", "patient_id", "n_sectors",
    "spatial_primary_eligible", "spatial_order_eligible",
    "n_events_reference", "n_ubiquitous_reference",
    "n_nonubiquitous_reference", "n_private_reference",
    "nonubiquitous_fraction_reference",
    "mean_occupancy_fraction", "median_occupancy_fraction"
  ))
  setorder(out, patient_id)
  out
}

# -------------------------- require Script 01 outputs -------------------------
sample_map_path <- file.path(PATHS$intermediate_dir, "hcc_sample_map.rds")
patient_qc_path <- file.path(PATHS$intermediate_dir, "hcc_patient_metadata_qc.rds")

if (!file.exists(sample_map_path) || !file.exists(patient_qc_path)) {
  stop(
    "Required Script 01 intermediate file(s) are missing.\n",
    "Run and validate 01_hcc_metadata_qc.R before Script 02.\n",
    "Expected:\n  ", sample_map_path, "\n  ", patient_qc_path
  )
}

sample_map <- as.data.table(readRDS(sample_map_path))
patient_qc <- as.data.table(readRDS(patient_qc_path))

require_columns(
  sample_map,
  c("patient_id", "sample_id", "category", "n_sectors",
    "spatial_primary_eligible", "spatial_order_eligible"),
  "hcc_sample_map"
)
require_columns(
  patient_qc,
  c("patient_id", "n_sectors", "spatial_primary_eligible", "spatial_order_eligible"),
  "hcc_patient_metadata_qc"
)

sample_map[, patient_id := as.character(patient_id)]
sample_map[, sample_id := as.character(sample_id)]
patient_qc[, patient_id := as.character(patient_id)]

if (uniqueN(sample_map$patient_id) != EXPECTED_LEGACY$n_patients ||
    nrow(sample_map) != EXPECTED_LEGACY$n_tumor_sectors) {
  stop(
    "Script 01 intermediate sample map no longer matches the validated cohort: ",
    uniqueN(sample_map$patient_id), " patients / ", nrow(sample_map), " sectors."
  )
}
if (anyDuplicated(sample_map$sample_id) > 0L) {
  stop("Duplicate sample_id found in the validated Script 01 sample map.")
}

# ----------------------------- load mutations --------------------------------
cat("Script 02: HCC event construction\n")
cat("============================================\n\n")
cat("Loading somatic_mutations...\n")

mut <- load_expected_object(INPUT_FILES$mutations, "somatic_mutations")
setDT(mut)

required_mut_cols <- c(
  "PATIENT_ID", "CHROM", "POS", "REF_ALLELE", "ALT_ALLELE",
  "REF_COUNT", "ALT_COUNT", "VAF", "GENE", "IMPACT"
)
require_columns(mut, required_mut_cols, "somatic_mutations")

raw_n_rows <- nrow(mut)

# Record a compact schema before modifying the object.
mutation_schema <- data.table(
  column = required_mut_cols,
  class = vapply(required_mut_cols, function(cc) paste(class(mut[[cc]]), collapse = "/"), character(1)),
  n_missing = vapply(required_mut_cols, function(cc) sum(is.na(mut[[cc]])), integer(1))
)
if (DEBUG) {
  fwrite(mutation_schema, file.path(PATHS$qc_dir, "02_mutation_schema.tsv"), sep = "\t")
}

# Keep only columns required for the sampling analysis. This is intentionally
# done by reference to reduce memory use on the ~10-million-row call table.
drop_cols <- setdiff(names(mut), required_mut_cols)
if (length(drop_cols) > 0L) mut[, (drop_cols) := NULL]

# Standardize identifier/annotation columns without changing numerical filters.
mut[, PATIENT_ID := as.character(PATIENT_ID)]
mut[, CHROM := as.character(CHROM)]
mut[, REF_ALLELE := as.character(REF_ALLELE)]
mut[, ALT_ALLELE := as.character(ALT_ALLELE)]
mut[, GENE := as.character(GENE)]
mut[, IMPACT := as.character(IMPACT)]
setnames(mut, "PATIENT_ID", "sample_id")

# ------------------------ reproduce legacy sample filter ----------------------
validated_sample_ids <- sample_map$sample_id
n_unique_raw_sample_ids <- uniqueN(mut$sample_id)
n_unique_raw_sample_ids_matching <- uniqueN(mut[sample_id %chin% validated_sample_ids, sample_id])

mut <- mut[sample_id %chin% validated_sample_ids]
n_rows_tumour_mapped <- nrow(mut)

lookup_idx <- match(mut$sample_id, sample_map$sample_id)
if (anyNA(lookup_idx)) {
  stop("At least one mutation row retained after sample filtering has no patient mapping.")
}
mut[, patient_id := sample_map$patient_id[lookup_idx]]
rm(lookup_idx)
invisible(gc())

# ----------------------------- numeric filters -------------------------------
# Preserve the legacy definitions exactly:
#   total_depth = REF_COUNT + ALT_COUNT
#   ALT_COUNT >= 3, total_depth >= 10, VAF >= 0.05
# with non-missing ALT_COUNT, total_depth and VAF.
if (!is.numeric(mut$REF_COUNT) || !is.numeric(mut$ALT_COUNT) || !is.numeric(mut$VAF)) {
  stop(
    "REF_COUNT, ALT_COUNT and VAF must be numeric to reproduce the legacy filter. ",
    "Inspect 02_mutation_schema.tsv."
  )
}

mut[, total_depth := REF_COUNT + ALT_COUNT]

numeric_qc <- data.table(
  criterion = c(
    "tumour_mapped_rows",
    "nonmissing_ALT_COUNT_total_depth_VAF",
    paste0("ALT_COUNT_ge_", PARAMS$min_alt_count),
    paste0("total_depth_ge_", PARAMS$min_total_depth),
    paste0("VAF_ge_", PARAMS$min_vaf),
    "all_legacy_numeric_filters"
  ),
  n_rows = c(
    nrow(mut),
    mut[, sum(!is.na(ALT_COUNT) & !is.na(total_depth) & !is.na(VAF))],
    mut[, sum(!is.na(ALT_COUNT) & ALT_COUNT >= PARAMS$min_alt_count)],
    mut[, sum(!is.na(total_depth) & total_depth >= PARAMS$min_total_depth)],
    mut[, sum(!is.na(VAF) & VAF >= PARAMS$min_vaf)],
    mut[, sum(
      !is.na(ALT_COUNT) & !is.na(total_depth) & !is.na(VAF) &
        ALT_COUNT >= PARAMS$min_alt_count &
        total_depth >= PARAMS$min_total_depth &
        VAF >= PARAMS$min_vaf
    )]
  )
)

mut <- mut[
  !is.na(ALT_COUNT) & !is.na(total_depth) & !is.na(VAF) &
    ALT_COUNT >= PARAMS$min_alt_count &
    total_depth >= PARAMS$min_total_depth &
    VAF >= PARAMS$min_vaf
]

n_rows_filtered <- nrow(mut)
legacy_filter_pass <- identical(as.integer(n_rows_filtered), as.integer(LEGACY_EXPECTED_FILTERED_ROWS))

filter_flow <- rbindlist(list(
  data.table(stage = "raw_somatic_mutation_rows", n_rows = raw_n_rows),
  data.table(stage = "rows_in_validated_tumour_samples", n_rows = n_rows_tumour_mapped),
  data.table(stage = "rows_after_legacy_numeric_filters", n_rows = n_rows_filtered)
))
if (DEBUG) {
  fwrite(filter_flow, file.path(PATHS$qc_dir, "02_mutation_filter_flow.tsv"), sep = "\t")
}
if (DEBUG) {
  fwrite(numeric_qc, file.path(PATHS$qc_dir, "02_numeric_filter_qc.tsv"), sep = "\t")
}

if (!legacy_filter_pass) {
  early_lines <- c(
    "Script 02: EARLY STOP",
    "====================================",
    "",
    sprintf("Observed filtered rows: %s", fmt_int(n_rows_filtered)),
    sprintf("Legacy expected rows:   %s", fmt_int(LEGACY_EXPECTED_FILTERED_ROWS)),
    "Legacy filtering QC: FAIL",
    "",
    "Do not construct the reference until this discrepancy is resolved."
  )
  writeLines(early_lines, file.path(PATHS$qc_dir, "02_qc_summary.txt"))
  stop(paste(early_lines, collapse = "\n"))
}

# ------------------------ event-identity integrity QC -------------------------
# The primary unit is the exact variant CHROM:POS:REF:ALT, as in the legacy
# analysis. Missing identity components are not silently converted into strings.
invalid_event_component <-
  is.na(mut$CHROM) | !nzchar(mut$CHROM) |
  is.na(mut$POS) |
  is.na(mut$REF_ALLELE) | !nzchar(mut$REF_ALLELE) |
  is.na(mut$ALT_ALLELE) | !nzchar(mut$ALT_ALLELE)

n_invalid_event_rows <- sum(invalid_event_component)
if (n_invalid_event_rows > 0L) {
  invalid_examples <- utils::head(mut[which(invalid_event_component)], 100L)
  if (DEBUG) {
    fwrite(
      invalid_examples,
      file.path(PATHS$qc_dir, "02_invalid_event_identity_examples.tsv"),
      sep = "\t"
    )
  }
  stop(
    n_invalid_event_rows,
    " filtered row(s) have an incomplete CHROM/POS/REF/ALT identity. ",
    "See 02_invalid_event_identity_examples.tsv."
  )
}
rm(invalid_event_component)

# Non-finite values would be biologically/technically anomalous and were not
# explicitly protected against in the legacy script. Abort rather than silently
# retain them if they exist.
nonfinite_numeric <-
  !is.finite(mut$ALT_COUNT) |
  !is.finite(mut$total_depth) |
  !is.finite(mut$VAF)
if (any(nonfinite_numeric)) {
  if (DEBUG) {
    fwrite(
      utils::head(mut[which(nonfinite_numeric)], 100L),
      file.path(PATHS$qc_dir, "02_nonfinite_numeric_examples.tsv"),
      sep = "\t"
    )
  }
  stop(
    sum(nonfinite_numeric),
    " filtered row(s) contain non-finite ALT_COUNT/total_depth/VAF."
  )
}
rm(nonfinite_numeric)

# ---------------------------- IMPACT audit -----------------------------------
impact_distribution <- mut[, .N, by = IMPACT][order(-N, IMPACT)]
impact_distribution[, fraction := N / sum(N)]
if (DEBUG) {
  fwrite(impact_distribution, file.path(PATHS$qc_dir, "02_impact_distribution.tsv"), sep = "\t")
}

impact_lookup <- unique(mut[, .(IMPACT)])
impact_lookup[, impact_norm := fifelse(
  is.na(IMPACT) | !nzchar(trimws(IMPACT)),
  NA_character_,
  toupper(trimws(IMPACT))
)]

unexpected_impacts <- sort(unique(
  impact_lookup[!is.na(impact_norm) & !(impact_norm %chin% KNOWN_IMPACT_CATEGORIES), impact_norm]
))

if (length(unexpected_impacts) > 0L) {
  if (DEBUG) {
    fwrite(
      impact_lookup[impact_norm %chin% unexpected_impacts],
      file.path(PATHS$qc_dir, "02_unexpected_impact_categories.tsv"),
      sep = "\t"
    )
  }
  stop(
    "Unexpected IMPACT category/categories: ", paste(unexpected_impacts, collapse = ", "),
    ". Protein-altering sensitivity has not been defined for these values."
  )
}

impact_lookup[, protein_altering_row :=
  !is.na(impact_norm) & impact_norm %chin% PROTEIN_ALTERING_IMPACTS]

# Tiny lookup joined by raw IMPACT avoids creating a normalized 10-million-row
# character vector.
mut[impact_lookup, on = .(IMPACT), protein_altering_row := i.protein_altering_row]
if (anyNA(mut$protein_altering_row)) {
  stop("Internal IMPACT mapping error: protein_altering_row contains NA after lookup.")
}

# Count/annotation columns have served their QC purpose. Drop them before the
# large grouping step to reduce peak memory use.
post_qc_drop <- intersect(
  c("REF_COUNT", "ALT_COUNT", "VAF", "GENE", "IMPACT", "total_depth"),
  names(mut)
)
if (length(post_qc_drop) > 0L) mut[, (post_qc_drop) := NULL]
invisible(gc())

# ---------------------- collapse to sample-event presence --------------------
cat("Legacy filter reproduced exactly: ", fmt_int(n_rows_filtered), " rows.\n", sep = "")
cat("Collapsing annotation-level rows to exact variant presence by sector...\n")

# If the source contains multiple annotation rows for one exact variant in one
# sector, they are one biological presence event. Protein-altering status is
# TRUE if any retained annotation row for that sample-event is HIGH/MODERATE.
presence_coords <- mut[, .(
  protein_altering_in_sample = any(protein_altering_row)
), by = .(
  patient_id, sample_id, CHROM, POS, REF_ALLELE, ALT_ALLELE
)]

n_unique_sample_event_presences <- nrow(presence_coords)
n_annotation_rows_collapsed <- n_rows_filtered - n_unique_sample_event_presences

# The filtered mutation table is no longer needed after the exact presence table
# has been constructed.
rm(mut)
invisible(gc())

# Deterministic ordering makes event keys stable for a fixed input dataset.
setorder(
  presence_coords,
  patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE, sample_id
)

# event_key is unique for patient x exact variant. Including patient_id is
# deliberate: all downstream sampling is within patient and this prevents an
# accidental cross-patient event merge.
presence_coords[, event_key := .GRP,
  by = .(patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE)
]

# Build one row per patient-specific exact variant. Event-level protein-altering
# status is intrinsic to the event in this sensitivity analysis: if any retained
# annotation labels it HIGH/MODERATE, all regional presences of that event are
# kept in the protein-altering subset.
event_dictionary <- presence_coords[, .(
  protein_altering = any(protein_altering_in_sample),
  n_samples_with_protein_altering_annotation = sum(protein_altering_in_sample),
  n_samples_present = .N
), by = .(
  event_key, patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE
)]

setcolorder(event_dictionary, c(
  "event_key", "patient_id",
  "CHROM", "POS", "REF_ALLELE", "ALT_ALLELE",
  "protein_altering", "n_samples_present",
  "n_samples_with_protein_altering_annotation"
))
setorder(event_dictionary, event_key)

# Audit whether the HIGH/MODERATE flag was inconsistent across regional calls of
# the same exact event. This does not alter the primary all-variant analysis.
annotation_flag_mixed <- event_dictionary[
  n_samples_with_protein_altering_annotation > 0L &
    n_samples_with_protein_altering_annotation < n_samples_present
]
if (DEBUG) {
  fwrite(
    annotation_flag_mixed,
    file.path(PATHS$qc_dir, "02_events_with_mixed_protein_annotation_across_samples.tsv"),
    sep = "\t"
  )
}

# Compact presence table used by all later HCC sampling scripts.
event_presence <- presence_coords[, .(patient_id, sample_id, event_key)]
setorder(event_presence, patient_id, sample_id, event_key)
rm(presence_coords)
invisible(gc())

duplicate_presence <- event_presence[, .N, by = .(patient_id, sample_id, event_key)][N > 1L]
if (nrow(duplicate_presence) > 0L) {
  stop("Internal error: duplicate patient/sample/event_key remains after presence collapsing.")
}
rm(duplicate_presence)

# -------------------- full available-sector reference ------------------------
cat("Constructing full available-sector patient references...\n")

reference_events <- event_presence[, .(
  full_count = uniqueN(sample_id)
), by = .(patient_id, event_key)]

reference_events <- merge(
  reference_events,
  patient_qc[, .(patient_id, n_sectors)],
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
reference_events <- merge(
  reference_events,
  event_dictionary[, .(event_key, protein_altering)],
  by = "event_key",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(reference_events$n_sectors) || anyNA(reference_events$protein_altering)) {
  stop("Internal error while joining patient/event metadata to reference_events.")
}
if (any(reference_events$full_count < 1L) ||
    any(reference_events$full_count > reference_events$n_sectors)) {
  stop("Reference integrity failure: full_count is outside [1, n_sectors].")
}

reference_events[, occupancy_fraction := full_count / n_sectors]
reference_events[, is_ubiquitous := full_count == n_sectors]
reference_events[, is_nonubiquitous := full_count > 0L & full_count < n_sectors]
reference_events[, is_private := full_count == 1L]
setorder(reference_events, patient_id, event_key)

# Sanity partition: every event is either ubiquitous or non-ubiquitous.
if (any((reference_events$is_ubiquitous + reference_events$is_nonubiquitous) != 1L)) {
  stop("Reference classification integrity failure: ubiquitous/non-ubiquitous is not a partition.")
}
if (any(reference_events$is_private & !reference_events$is_nonubiquitous &
        reference_events$n_sectors > 1L)) {
  stop("Reference classification integrity failure for private events.")
}

# ------------------------ patient/sample summaries ---------------------------
reference_summary_all <- summarise_reference_set(
  reference_events,
  patient_qc,
  "all_filtered"
)
reference_summary_protein <- summarise_reference_set(
  reference_events[protein_altering == TRUE],
  patient_qc,
  "protein_altering"
)
reference_summary <- rbindlist(
  list(reference_summary_all, reference_summary_protein),
  use.names = TRUE,
  fill = TRUE
)
setorder(reference_summary, event_set, patient_id)

# Count detected events per sector for QC, retaining sectors with zero events.
sample_event_counts_all <- event_presence[, .(n_events_all_filtered = .N),
  by = .(patient_id, sample_id)
]
protein_keys <- event_dictionary[protein_altering == TRUE, event_key]
sample_event_counts_protein <- event_presence[event_key %in% protein_keys,
  .(n_events_protein_altering = .N),
  by = .(patient_id, sample_id)
]

sample_event_counts <- merge(
  sample_map[, .(patient_id, sample_id, n_sectors,
                 spatial_primary_eligible, spatial_order_eligible)],
  sample_event_counts_all,
  by = c("patient_id", "sample_id"),
  all.x = TRUE,
  sort = FALSE
)
sample_event_counts <- merge(
  sample_event_counts,
  sample_event_counts_protein,
  by = c("patient_id", "sample_id"),
  all.x = TRUE,
  sort = FALSE
)
sample_event_counts[is.na(n_events_all_filtered), n_events_all_filtered := 0L]
sample_event_counts[is.na(n_events_protein_altering), n_events_protein_altering := 0L]
setorder(sample_event_counts, patient_id, sample_id)

n_zero_event_sectors_all <- sample_event_counts[n_events_all_filtered == 0L, .N]
n_zero_event_sectors_protein <- sample_event_counts[n_events_protein_altering == 0L, .N]
n_zero_event_patients_all <- reference_summary_all[n_events_reference == 0L, .N]
n_zero_event_patients_protein <- reference_summary_protein[n_events_reference == 0L, .N]

# ------------------------------- save outputs --------------------------------
cat("Saving compact event objects...\n")

saveRDS(
  event_presence,
  file.path(PATHS$intermediate_dir, "hcc_event_presence.rds")
)
saveRDS(
  event_dictionary,
  file.path(PATHS$intermediate_dir, "hcc_event_dictionary.rds")
)
saveRDS(
  reference_events,
  file.path(PATHS$intermediate_dir, "hcc_reference_events.rds")
)

# Always written (not gated by DEBUG): Script 12 reads this file
# back as the TRACERx100 comparison reference, not just as a diagnostic.
fwrite(reference_summary, file.path(PATHS$qc_dir, "02_reference_summary_by_patient.tsv"), sep = "\t")
if (DEBUG) {
  fwrite(sample_event_counts, file.path(PATHS$qc_dir, "02_sample_event_counts.tsv"), sep = "\t")
}

# Global event-set summary.
event_set_summary <- rbindlist(list(
  data.table(
    event_set = "all_filtered",
    n_patient_specific_events = nrow(reference_events),
    n_sample_event_presences = nrow(event_presence),
    n_patients_with_events = uniqueN(reference_events$patient_id),
    n_zero_event_patients = n_zero_event_patients_all,
    n_zero_event_sectors = n_zero_event_sectors_all
  ),
  data.table(
    event_set = "protein_altering",
    n_patient_specific_events = reference_events[protein_altering == TRUE, .N],
    n_sample_event_presences = event_presence[event_key %in% protein_keys, .N],
    n_patients_with_events = uniqueN(reference_events[protein_altering == TRUE, patient_id]),
    n_zero_event_patients = n_zero_event_patients_protein,
    n_zero_event_sectors = n_zero_event_sectors_protein
  )
))
if (DEBUG) {
  fwrite(event_set_summary, file.path(PATHS$qc_dir, "02_event_set_summary.tsv"), sep = "\t")
}

# ------------------------------- report --------------------------------------
impact_values_pretty <- impact_distribution[, paste0(
  ifelse(is.na(IMPACT), "<NA>", IMPACT), "=", fmt_int(N)
)]

summary_lines <- c(
  "Script 02: HCC event construction",
  "============================================",
  "",
  sprintf("Raw somatic mutation rows: %s", fmt_int(raw_n_rows)),
  sprintf("Rows mapped to validated tumour sectors: %s", fmt_int(n_rows_tumour_mapped)),
  sprintf("Rows after legacy numeric filters: %s", fmt_int(n_rows_filtered)),
  sprintf("Legacy filtered-row target: %s", fmt_int(LEGACY_EXPECTED_FILTERED_ROWS)),
  sprintf("Legacy filtering QC: %s", ifelse(legacy_filter_pass, "PASS", "FAIL")),
  "",
  sprintf("Unique raw mutation sample IDs: %s", fmt_int(n_unique_raw_sample_ids)),
  sprintf("Unique raw mutation sample IDs matching validated tumour sectors: %s", fmt_int(n_unique_raw_sample_ids_matching)),
  sprintf("Validated tumour sectors: %s", fmt_int(nrow(sample_map))),
  "",
  sprintf("Unique sector x exact-variant presences: %s", fmt_int(n_unique_sample_event_presences)),
  sprintf("Annotation-level filtered rows collapsed as duplicates: %s", fmt_int(n_annotation_rows_collapsed)),
  sprintf("Patient-specific exact variants in reference: %s", fmt_int(nrow(reference_events))),
  "",
  sprintf("Protein-altering definition: IMPACT in {%s}", paste(PROTEIN_ALTERING_IMPACTS, collapse = ", ")),
  sprintf("Patient-specific protein-altering variants: %s", fmt_int(reference_events[protein_altering == TRUE, .N])),
  sprintf("Events with mixed protein-altering flag across regional calls: %s", fmt_int(nrow(annotation_flag_mixed))),
  "",
  sprintf("All-filtered zero-event patients: %s", fmt_int(n_zero_event_patients_all)),
  sprintf("All-filtered zero-event sectors: %s", fmt_int(n_zero_event_sectors_all)),
  sprintf("Protein-altering zero-event patients: %s", fmt_int(n_zero_event_patients_protein)),
  sprintf("Protein-altering zero-event sectors: %s", fmt_int(n_zero_event_sectors_protein)),
  "",
  "Observed IMPACT distribution after legacy numeric filtering:",
  paste0("  ", impact_values_pretty),
  "",
  "Primary reference terminology:",
  "  ubiquitous      = event present in every available tumour sector of the patient",
  "  non-ubiquitous  = event present in >=1 but not all available tumour sectors",
  "  private         = event present in exactly one available tumour sector",
  "",
  "The reference is the FULL AVAILABLE-SECTOR REFERENCE, not whole-tumour ground truth.",
  "No downsampling was performed in Script 02.",
  "Review Script 02 QC outputs before proceeding to sampling-function tests."
)

writeLines(summary_lines, file.path(PATHS$qc_dir, "02_qc_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
