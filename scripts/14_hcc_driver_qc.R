# ==============================================================================
# 14_hcc_driver_qc.R
#
# Purpose
#   QC-only preparation for the collaborator-requested non-ubiquitous driver
#   mutation sampling analysis.
#
#   This script deliberately DOES NOT yet define the final driver-event universe
#   or perform downsampling. It:
#     1. freezes and validates the supplied driver-gene list;
#     2. reproduces Script 02's validated HCC mutation/sample/numeric filters;
#     3. maps filtered annotations in driver genes back to the validated exact
#        event_key universe from Script 02;
#     4. audits IMPACT categories, multi-gene annotations, and TERT specifically;
#     5. writes the exact information needed to choose the final biologically
#        justified driver-event rule without guessing.
#
# Required prerequisites
#   - 00_config.R
#   - Script 01 and Script 02 validated intermediates in the intermediate/ output directory
#   - Driver_48genes_DICER1.csv supplied by collaborators
#
# Driver CSV location
#   Preferred: set environment variable DRIVER_CSV_PATH to the full path.
#   Otherwise the script searches PROJECT_ROOT, V3_DIR, and this script folder.
# ============================================================================== 

suppressPackageStartupMessages(library(data.table))

# ------------------------------ locate files ---------------------------------
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

cat("Script 14: HCC driver QC\n")
cat("===========================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------- constants -----------------------------------
LEGACY_EXPECTED_FILTERED_ROWS <- 10612828L
PROTEIN_ALTERING_IMPACTS <- c("HIGH", "MODERATE")
KNOWN_IMPACT_CATEGORIES <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

EXPECTED_DRIVER_CSV_MD5 <- "e7ba1edf8c34ad2d3532743f6d61ca67"
EXPECTED_DRIVER_ROWS <- 49L
EXPECTED_BASE_DRIVER_ROWS <- 48L
EXPECTED_EXTRA_GENE <- "DICER1"

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
    stop(
      object_name, " is missing required column(s): ",
      paste(miss, collapse = ", ")
    )
  }
}

collapse_sorted_unique <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0L) return(NA_character_)
  paste(sort(x), collapse = ";")
}

find_driver_csv <- function() {
  env_path <- Sys.getenv("DRIVER_CSV_PATH", unset = "")

  candidates <- unique(c(
    if (nzchar(env_path)) path.expand(env_path) else character(),
    file.path(PROJECT_ROOT, "Driver_48genes_DICER1.csv"),
    file.path(PATHS$output_dir, "Driver_48genes_DICER1.csv"),
    file.path(script_dir, "Driver_48genes_DICER1.csv")
  ))

  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop(
      "Could not find Driver_48genes_DICER1.csv.\n",
      "Set DRIVER_CSV_PATH to its full path, e.g.:\n",
      "  Sys.setenv(DRIVER_CSV_PATH='/path/to/Driver_48genes_DICER1.csv')"
    )
  }

  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

# ------------------------------ output paths ---------------------------------
ext_dir <- file.path(PATHS$output_dir, "results", "extensions", "driver")
ext_qc_dir <- file.path(PATHS$output_dir, "results", "qc", "extensions")
dir.create(ext_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ext_qc_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- validated prerequisites ---------------------------
sample_map_path <- file.path(PATHS$intermediate_dir, "hcc_sample_map.rds")
event_dictionary_path <- file.path(PATHS$intermediate_dir, "hcc_event_dictionary.rds")
reference_path <- file.path(PATHS$intermediate_dir, "hcc_reference_events.rds")
presence_path <- file.path(PATHS$intermediate_dir, "hcc_event_presence.rds")

required_paths <- c(
  sample_map_path,
  event_dictionary_path,
  reference_path,
  presence_path,
  INPUT_FILES$mutations
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing validated prerequisite file(s):\n",
    paste0("  ", missing_paths, collapse = "\n")
  )
}

sample_map <- as.data.table(readRDS(sample_map_path))
event_dictionary <- as.data.table(readRDS(event_dictionary_path))
reference_events <- as.data.table(readRDS(reference_path))
event_presence <- as.data.table(readRDS(presence_path))

require_columns(
  sample_map,
  c("patient_id", "sample_id", "n_sectors"),
  "hcc_sample_map"
)
require_columns(
  event_dictionary,
  c(
    "event_key", "patient_id", "CHROM", "POS", "REF_ALLELE", "ALT_ALLELE",
    "protein_altering", "n_samples_present"
  ),
  "hcc_event_dictionary"
)
require_columns(
  reference_events,
  c(
    "event_key", "patient_id", "full_count", "n_sectors",
    "is_ubiquitous", "is_nonubiquitous", "is_private"
  ),
  "hcc_reference_events"
)
require_columns(
  event_presence,
  c("patient_id", "sample_id", "event_key"),
  "hcc_event_presence"
)

sample_map[, patient_id := as.character(patient_id)]
sample_map[, sample_id := as.character(sample_id)]
event_dictionary[, patient_id := as.character(patient_id)]
reference_events[, patient_id := as.character(patient_id)]
event_presence[, patient_id := as.character(patient_id)]
event_presence[, sample_id := as.character(sample_id)]

if (uniqueN(sample_map$patient_id) != EXPECTED_LEGACY$n_patients ||
    nrow(sample_map) != EXPECTED_LEGACY$n_tumor_sectors) {
  stop(
    "Validated sample-map cohort mismatch: observed ",
    uniqueN(sample_map$patient_id), " patients / ", nrow(sample_map),
    " sectors; expected ", EXPECTED_LEGACY$n_patients, " / ",
    EXPECTED_LEGACY$n_tumor_sectors, "."
  )
}
if (anyDuplicated(sample_map$sample_id) > 0L) {
  stop("Duplicate sample_id in validated sample map.")
}
if (anyDuplicated(event_dictionary$event_key) > 0L) {
  stop("Duplicate event_key in validated event dictionary.")
}

# ---------------------------- driver-gene list -------------------------------
driver_csv <- find_driver_csv()
driver_md5 <- unname(tools::md5sum(driver_csv))

if (!identical(driver_md5, EXPECTED_DRIVER_CSV_MD5)) {
  stop(
    "Driver CSV MD5 mismatch.\n",
    "Expected: ", EXPECTED_DRIVER_CSV_MD5, "\n",
    "Observed: ", driver_md5, "\n",
    "File: ", driver_csv, "\n",
    "Do not proceed until the supplied collaborator file is confirmed."
  )
}

driver_raw <- fread(driver_csv, header = TRUE)
if (ncol(driver_raw) != 1L) {
  stop("Driver CSV must contain exactly one gene column; observed ", ncol(driver_raw), ".")
}

driver_genes <- trimws(as.character(driver_raw[[1L]]))
if (anyNA(driver_genes) || any(!nzchar(driver_genes))) {
  stop("Driver CSV contains missing/empty gene symbols.")
}
if (anyDuplicated(driver_genes) > 0L) {
  stop(
    "Driver CSV contains duplicate gene symbol(s): ",
    paste(unique(driver_genes[duplicated(driver_genes)]), collapse = ", ")
  )
}
if (length(driver_genes) != EXPECTED_DRIVER_ROWS) {
  stop(
    "Driver-list size mismatch: expected ", EXPECTED_DRIVER_ROWS,
    " rows (48 original + DICER1), observed ", length(driver_genes), "."
  )
}
if (sum(driver_genes == EXPECTED_EXTRA_GENE) != 1L) {
  stop("Expected exactly one DICER1 row in collaborator driver list.")
}

base_driver_genes <- driver_genes[driver_genes != EXPECTED_EXTRA_GENE]
if (length(base_driver_genes) != EXPECTED_BASE_DRIVER_ROWS) {
  stop("Expected exactly 48 non-DICER1 driver genes.")
}

if (DEBUG) {
  fwrite(
    data.table(
      gene = driver_genes,
      list_membership = ifelse(
        driver_genes == EXPECTED_EXTRA_GENE,
        "DICER1_addition",
        "original_48"
      )
    ),
    file.path(ext_qc_dir, "14_driver_gene_list_frozen.tsv"),
    sep = "\t"
  )
}

# ----------------------------- load mutations --------------------------------
cat("Loading somatic_mutations and reproducing Script 02 filters...\n")

mut <- load_expected_object(INPUT_FILES$mutations, "somatic_mutations")
setDT(mut)

required_mut_cols <- c(
  "PATIENT_ID", "CHROM", "POS", "REF_ALLELE", "ALT_ALLELE",
  "REF_COUNT", "ALT_COUNT", "VAF", "GENE", "IMPACT"
)
require_columns(mut, required_mut_cols, "somatic_mutations")

raw_n_rows <- nrow(mut)

drop_cols <- setdiff(names(mut), required_mut_cols)
if (length(drop_cols) > 0L) mut[, (drop_cols) := NULL]

mut[, PATIENT_ID := as.character(PATIENT_ID)]
mut[, CHROM := as.character(CHROM)]
mut[, REF_ALLELE := as.character(REF_ALLELE)]
mut[, ALT_ALLELE := as.character(ALT_ALLELE)]
mut[, GENE := trimws(as.character(GENE))]
mut[, IMPACT := as.character(IMPACT)]
setnames(mut, "PATIENT_ID", "sample_id")

validated_sample_ids <- sample_map$sample_id
mut <- mut[sample_id %chin% validated_sample_ids]
n_rows_tumour_mapped <- nrow(mut)

lookup_idx <- match(mut$sample_id, sample_map$sample_id)
if (anyNA(lookup_idx)) {
  stop("At least one retained mutation row has no patient mapping.")
}
mut[, patient_id := sample_map$patient_id[lookup_idx]]
rm(lookup_idx)

if (!is.numeric(mut$REF_COUNT) ||
    !is.numeric(mut$ALT_COUNT) ||
    !is.numeric(mut$VAF)) {
  stop("REF_COUNT, ALT_COUNT and VAF must be numeric.")
}

mut[, total_depth := REF_COUNT + ALT_COUNT]
mut <- mut[
  !is.na(ALT_COUNT) & !is.na(total_depth) & !is.na(VAF) &
    ALT_COUNT >= PARAMS$min_alt_count &
    total_depth >= PARAMS$min_total_depth &
    VAF >= PARAMS$min_vaf
]

n_rows_filtered <- nrow(mut)
if (!identical(as.integer(n_rows_filtered), LEGACY_EXPECTED_FILTERED_ROWS)) {
  stop(
    "Legacy mutation-filter regression FAILED: expected ",
    LEGACY_EXPECTED_FILTERED_ROWS, ", observed ", n_rows_filtered, "."
  )
}

invalid_event_component <-
  is.na(mut$CHROM) | !nzchar(mut$CHROM) |
  is.na(mut$POS) |
  is.na(mut$REF_ALLELE) | !nzchar(mut$REF_ALLELE) |
  is.na(mut$ALT_ALLELE) | !nzchar(mut$ALT_ALLELE)

if (any(invalid_event_component)) {
  stop(sum(invalid_event_component), " filtered row(s) have invalid exact-variant identity.")
}

nonfinite_numeric <-
  !is.finite(mut$ALT_COUNT) |
  !is.finite(mut$total_depth) |
  !is.finite(mut$VAF)
if (any(nonfinite_numeric)) {
  stop(sum(nonfinite_numeric), " filtered row(s) contain non-finite numeric values.")
}

mut[, impact_norm := fifelse(
  is.na(IMPACT) | !nzchar(trimws(IMPACT)),
  NA_character_,
  toupper(trimws(IMPACT))
)]

unexpected_impacts <- sort(unique(
  mut[
    !is.na(impact_norm) & !(impact_norm %chin% KNOWN_IMPACT_CATEGORIES),
    impact_norm
  ]
))
if (length(unexpected_impacts) > 0L) {
  stop(
    "Unexpected IMPACT categories after validated filtering: ",
    paste(unexpected_impacts, collapse = ", ")
  )
}

mut[, protein_altering_row :=
      !is.na(impact_norm) & impact_norm %chin% PROTEIN_ALTERING_IMPACTS]

# ----------------------- restrict to supplied driver genes --------------------
# The collaborator list uses uppercase HGNC-style symbols. Abort if the source
# contains case-only matches that would otherwise be silently missed.
case_only_gene_matches <- sort(unique(
  mut[
    !is.na(GENE) &
      toupper(GENE) %chin% driver_genes &
      !(GENE %chin% driver_genes),
    GENE
  ]
))
if (length(case_only_gene_matches) > 0L) {
  stop(
    "Case-only driver-gene symbol mismatch(es) in mutation annotations: ",
    paste(case_only_gene_matches, collapse = ", ")
  )
}

driver_rows <- mut[
  !is.na(GENE) & GENE %chin% driver_genes,
  .(
    patient_id, sample_id, CHROM, POS, REF_ALLELE, ALT_ALLELE,
    GENE, impact_norm, protein_altering_row
  )
]

if (nrow(driver_rows) == 0L) {
  stop("No filtered mutation annotation rows map to the supplied driver genes.")
}

# The ~10.6M-row filtered table is no longer needed.
rm(mut)
invisible(gc(FALSE))

# One row per exact patient-specific variant x driver gene. Keeping gene-level
# mapping separate lets us audit rare overlapping annotations while the eventual
# sampling event remains one exact CHROM:POS:REF:ALT event.
driver_event_gene <- driver_rows[, .(
  n_annotation_rows = .N,
  n_samples_annotated = uniqueN(sample_id),
  impact_categories = collapse_sorted_unique(impact_norm),
  any_high_moderate = any(protein_altering_row),
  all_high_moderate = all(protein_altering_row)
), by = .(
  patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE, GENE
)]

# Collapse directly from annotation rows to one exact event. This avoids any
# ambiguity from concatenated per-gene IMPACT strings when an event overlaps
# more than one supplied driver gene.
driver_event_candidate <- driver_rows[, .(
  driver_genes = collapse_sorted_unique(GENE),
  n_driver_genes = uniqueN(GENE),
  impact_categories = collapse_sorted_unique(impact_norm),
  any_high_moderate = any(protein_altering_row),
  all_high_moderate = all(protein_altering_row),
  n_annotation_rows = .N,
  n_samples_annotated = uniqueN(sample_id),
  contains_TERT = any(GENE == "TERT"),
  contains_DICER1 = any(GENE == "DICER1")
), by = .(
  patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE
)]

# Validated Script 02 coordinates define the authoritative event_key universe.
dict_key <- event_dictionary[, .(
  event_key,
  patient_id,
  CHROM = as.character(CHROM),
  POS,
  REF_ALLELE = as.character(REF_ALLELE),
  ALT_ALLELE = as.character(ALT_ALLELE),
  script02_protein_altering = as.logical(protein_altering),
  script02_n_samples_present = as.integer(n_samples_present)
)]

if (anyDuplicated(
  dict_key[, .(patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE)]
) > 0L) {
  stop("Validated event dictionary has duplicate patient-specific exact coordinates.")
}

driver_event_candidate <- merge(
  driver_event_candidate,
  dict_key,
  by = c("patient_id", "CHROM", "POS", "REF_ALLELE", "ALT_ALLELE"),
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(driver_event_candidate$event_key)) {
  if (DEBUG) {
    fwrite(
      driver_event_candidate[is.na(event_key)],
      file.path(ext_qc_dir, "14_driver_events_unmatched_to_script02.tsv"),
      sep = "\t"
    )
  }
  stop(
    sum(is.na(driver_event_candidate$event_key)),
    " driver exact event(s) did not map to Script 02 event_key universe."
  )
}

# Script 02's protein_altering flag is event-wide across ALL annotations, while
# any_high_moderate here is deliberately DRIVER-GENE-ANNOTATION-specific. They
# should usually agree, but a mismatch can be biologically legitimate for an
# overlapping annotation and is therefore audited rather than silently forced.
impact_flag_mismatch <- driver_event_candidate[
  any_high_moderate != script02_protein_altering
]
if (DEBUG) {
  fwrite(
    impact_flag_mismatch,
    file.path(ext_qc_dir, "14_driver_vs_script02_protein_flag_mismatch.tsv"),
    sep = "\t"
  )
}

# Attach authoritative reference occupancy/class labels.
ref_key <- reference_events[, .(
  event_key,
  full_count,
  n_sectors,
  is_ubiquitous,
  is_nonubiquitous,
  is_private
)]
if (anyDuplicated(ref_key$event_key) > 0L) {
  stop("Duplicate event_key in validated reference_events.")
}

driver_event_candidate <- merge(
  driver_event_candidate,
  ref_key,
  by = "event_key",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(driver_event_candidate$full_count) ||
    anyNA(driver_event_candidate$n_sectors)) {
  stop("Driver candidate failed to map to validated Script 02 reference labels.")
}

# Reconstruct sample-presence count as an independent regression check.
presence_count <- event_presence[
  event_key %in% driver_event_candidate$event_key,
  .(presence_count_from_script02 = uniqueN(sample_id)),
  by = event_key
]
driver_event_candidate <- merge(
  driver_event_candidate,
  presence_count,
  by = "event_key",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(driver_event_candidate$presence_count_from_script02) ||
    any(driver_event_candidate$presence_count_from_script02 != driver_event_candidate$full_count) ||
    any(driver_event_candidate$script02_n_samples_present != driver_event_candidate$full_count)) {
  stop("Driver candidate presence-count regression against Script 02 FAILED.")
}

setorder(driver_event_candidate, patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE)

# ------------------------------- QC outputs ----------------------------------
# 1) IMPACT distribution among supplied driver-gene annotation rows.
driver_impact_qc <- driver_rows[, .N, by = .(GENE, impact_norm)]
setorder(driver_impact_qc, GENE, impact_norm)
if (DEBUG) {
  fwrite(
    driver_impact_qc,
    file.path(ext_qc_dir, "14_driver_filtered_annotation_impact_by_gene.tsv"),
    sep = "\t"
  )
}

# 2) Exact-event candidate table. This is deliberately pre-eligibility.
fwrite(
  driver_event_candidate,
  file.path(ext_dir, "14_driver_exact_event_candidates_preeligibility.tsv"),
  sep = "\t"
)

# 3) Events mapping to >1 supplied driver gene.
multi_gene_events <- driver_event_candidate[n_driver_genes > 1L]
if (DEBUG) {
  fwrite(
    multi_gene_events,
    file.path(ext_qc_dir, "14_driver_multigene_exact_events.tsv"),
    sep = "\t"
  )
}

# 4) TERT audit: all filtered exact events annotated to TERT, without silently
# deciding which non-HIGH/MODERATE events should count as promoter drivers.
tert_events <- driver_event_candidate[contains_TERT == TRUE]
setorder(tert_events, CHROM, POS, REF_ALLELE, ALT_ALLELE, patient_id)
if (DEBUG) {
  fwrite(
    tert_events,
    file.path(ext_qc_dir, "14_TERT_exact_event_audit.tsv"),
    sep = "\t"
  )
}

tert_locus_summary <- tert_events[, .(
  n_patient_specific_events = .N,
  n_patients = uniqueN(patient_id),
  n_nonubiquitous_patient_events = sum(is_nonubiquitous),
  n_private_patient_events = sum(is_private),
  any_high_moderate = any(any_high_moderate),
  impact_categories = collapse_sorted_unique(impact_categories)
), by = .(CHROM, POS, REF_ALLELE, ALT_ALLELE)]
setorder(tert_locus_summary, CHROM, POS, REF_ALLELE, ALT_ALLELE)
if (DEBUG) {
  fwrite(
    tert_locus_summary,
    file.path(ext_qc_dir, "14_TERT_locus_summary.tsv"),
    sep = "\t"
  )
}

# 5) Per-gene reference burden. Because an exact event can theoretically map to
# >1 driver gene, this table is gene-mapping-level and is not summed to derive
# the unique exact-event total.
event_gene_with_key <- merge(
  driver_event_gene,
  driver_event_candidate[, .(
    event_key, patient_id, CHROM, POS, REF_ALLELE, ALT_ALLELE,
    is_ubiquitous, is_nonubiquitous, is_private
  )],
  by = c("patient_id", "CHROM", "POS", "REF_ALLELE", "ALT_ALLELE"),
  all.x = TRUE,
  sort = FALSE
)

per_gene_summary <- event_gene_with_key[, .(
  n_patient_specific_exact_events = uniqueN(event_key),
  n_patients_with_event = uniqueN(patient_id),
  n_nonubiquitous_patient_events = uniqueN(event_key[is_nonubiquitous]),
  n_private_patient_events = uniqueN(event_key[is_private]),
  n_high_moderate_patient_events = uniqueN(event_key[any_high_moderate])
), by = GENE]

all_driver_gene_table <- merge(
  data.table(
    GENE = driver_genes,
    list_membership = ifelse(
      driver_genes == EXPECTED_EXTRA_GENE,
      "DICER1_addition",
      "original_48"
    )
  ),
  per_gene_summary,
  by = "GENE",
  all.x = TRUE,
  sort = FALSE
)

count_cols <- c(
  "n_patient_specific_exact_events",
  "n_patients_with_event",
  "n_nonubiquitous_patient_events",
  "n_private_patient_events",
  "n_high_moderate_patient_events"
)
for (cc in count_cols) {
  set(all_driver_gene_table, which(is.na(all_driver_gene_table[[cc]])), cc, 0L)
}
all_driver_gene_table[, driver_list_order := match(GENE, driver_genes)]
setorder(all_driver_gene_table, driver_list_order)
all_driver_gene_table[, driver_list_order := NULL]
if (DEBUG) {
  fwrite(
    all_driver_gene_table,
    file.path(ext_qc_dir, "14_driver_reference_summary_by_gene.tsv"),
    sep = "\t"
  )
}

# 6) Compact decision table for HIGH/MODERATE vs non-HIGH/MODERATE events.
decision_qc <- driver_event_candidate[, .(
  n_patient_specific_events = .N,
  n_patients = uniqueN(patient_id),
  n_nonubiquitous_patient_events = sum(is_nonubiquitous),
  n_private_patient_events = sum(is_private)
), by = .(
  contains_TERT,
  any_high_moderate,
  impact_categories
)]
setorder(decision_qc, -contains_TERT, -any_high_moderate, impact_categories)
if (DEBUG) {
  fwrite(
    decision_qc,
    file.path(ext_qc_dir, "14_driver_eligibility_decision_qc.tsv"),
    sep = "\t"
  )
}

# 7) Input/filter provenance.
provenance <- data.table(
  item = c(
    "driver_csv_path",
    "driver_csv_md5",
    "raw_somatic_mutation_rows",
    "rows_in_validated_tumour_samples",
    "rows_after_legacy_numeric_filters",
    "supplied_driver_genes_total",
    "supplied_original_driver_genes",
    "supplied_DICER1_rows",
    "driver_gene_annotation_rows_after_filter",
    "unique_patient_specific_driver_exact_event_candidates",
    "candidate_events_high_moderate",
    "candidate_events_non_high_moderate",
    "candidate_events_TERT",
    "candidate_events_DICER1",
    "multigene_candidate_events",
    "driver_vs_script02_protein_flag_mismatches",
    "candidate_events_driver_annotation_missing_in_some_present_sectors"
  ),
  value = c(
    driver_csv,
    driver_md5,
    as.character(raw_n_rows),
    as.character(n_rows_tumour_mapped),
    as.character(n_rows_filtered),
    as.character(length(driver_genes)),
    as.character(length(base_driver_genes)),
    as.character(sum(driver_genes == EXPECTED_EXTRA_GENE)),
    as.character(nrow(driver_rows)),
    as.character(nrow(driver_event_candidate)),
    as.character(driver_event_candidate[any_high_moderate == TRUE, .N]),
    as.character(driver_event_candidate[any_high_moderate == FALSE, .N]),
    as.character(driver_event_candidate[contains_TERT == TRUE, .N]),
    as.character(driver_event_candidate[contains_DICER1 == TRUE, .N]),
    as.character(nrow(multi_gene_events)),
    as.character(nrow(impact_flag_mismatch)),
    as.character(driver_event_candidate[n_samples_annotated < full_count, .N])
  )
)
if (DEBUG) {
  fwrite(
    provenance,
    file.path(ext_qc_dir, "14_driver_qc_provenance.tsv"),
    sep = "\t"
  )
}

# ------------------------- human-readable summary -----------------------------
missing_genes <- all_driver_gene_table[
  n_patient_specific_exact_events == 0L,
  GENE
]

summary_lines <- c(
  "HCC extensions - Script 14 driver QC",
  "===========================================",
  "",
  paste0("Driver CSV: ", driver_csv),
  paste0("Driver CSV MD5: ", driver_md5, " [PASS]"),
  sprintf(
    "Driver list: %d unique genes = %d original + DICER1 [PASS]",
    length(driver_genes), length(base_driver_genes)
  ),
  sprintf(
    "Legacy mutation filtering: %d rows [PASS]",
    n_rows_filtered
  ),
  sprintf(
    "Unique patient-specific exact-event candidates in supplied driver genes: %d",
    nrow(driver_event_candidate)
  ),
  sprintf(
    "  HIGH/MODERATE by validated Script 02 IMPACT rule: %d",
    driver_event_candidate[any_high_moderate == TRUE, .N]
  ),
  sprintf(
    "  non-HIGH/MODERATE: %d",
    driver_event_candidate[any_high_moderate == FALSE, .N]
  ),
  sprintf(
    "  TERT exact-event candidates: %d across %d patients",
    nrow(tert_events), uniqueN(tert_events$patient_id)
  ),
  sprintf(
    "  DICER1 exact-event candidates: %d across %d patients",
    driver_event_candidate[contains_DICER1 == TRUE, .N],
    uniqueN(driver_event_candidate[contains_DICER1 == TRUE, patient_id])
  ),
  sprintf(
    "  exact events annotated to >1 supplied driver gene: %d",
    nrow(multi_gene_events)
  ),
  if (length(missing_genes) == 0L) {
    "Supplied driver genes with zero filtered exact events: none"
  } else {
    paste0(
      "Supplied driver genes with zero filtered exact events (",
      length(missing_genes), "): ",
      paste(missing_genes, collapse = ", ")
    )
  },
  "",
  "Critical QC regressions:",
  "  exact coordinates -> Script 02 event_key: PASS",
  sprintf(
    "  driver-specific HIGH/MODERATE vs Script 02 event-wide flag mismatches: %d [AUDIT]",
    nrow(impact_flag_mismatch)
  ),
  "  exact-event presence counts -> Script 02 reference full_count: PASS",
  "",
  "NO DRIVER ELIGIBILITY RULE HAS BEEN APPLIED YET.",
  "Review 14_TERT_locus_summary.tsv and 14_driver_eligibility_decision_qc.tsv",
  "before defining the final event universe for sampling.",
  "",
  sprintf(
    "Elapsed time: %.1f seconds",
    proc.time()[["elapsed"]] - t0
  )
)

writeLines(
  summary_lines,
  file.path(ext_qc_dir, "14_driver_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
