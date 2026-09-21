# ==============================================================================
# 10_tracerx100_event_reference_qc.R
#
# Purpose
#   Freeze and validate the TRACERx100 primary-tumour mutation event universe
#   BEFORE any downsampling.
#
# Canonical cohort (frozen in Script 09):
#   - 100 patients
#   - 323 primary tumour regions (-R#)
#   - lymph-node (-LN#) and cfDNA samples excluded
#   - primary-only region depth: median 3, range 2-7
#
# Cross-cohort event definition
#   HCC sensitivity analysis used VEP IMPACT in {HIGH, MODERATE}.
#   TRACERx100 does not provide the VEP IMPACT field directly, but its MAF does
#   provide Variant_Classification and Consequence. We therefore define a
#   consequence-matched protein-altering set requiring BOTH an allowed
#   Variant_Classification and a compatible Consequence:
#
#     MODERATE-like
#       Missense_Mutation
#       In_Frame_Ins
#       In_Frame_Del (allowed by definition; absent in this public archive)
#
#     HIGH-like
#       Nonsense_Mutation
#       Splice_Site
#       Frame_Shift_Del
#       Frame_Shift_Ins
#       Translation_Start_Site
#       Nonstop_Mutation
#
#   Splice_Region is deliberately EXCLUDED. In Ensembl VEP,
#   splice_region_variant is LOW impact, whereas splice donor/acceptor variants
#   are HIGH impact.
#
# Unit of event:
#   exact (Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)
#   within patient, matching the exact-variant strategy used in HCC.
#
# Presence:
#   one public MAF row for sample x exact variant = called presence.
#
# RegionSum:
#   Script 09 established that RegionSum is present and internally consistent for
#   all tissue mutation rows. Here we independently expand RegionSum for the
#   selected protein-altering event set and verify that every primary tumour
#   region is represented by an R# count pair. We use the two values only as
#   numerator/denominator fields for QC; no new mutation call is created from
#   RegionSum in this script.
#
# IMPORTANT:
#   - no downsampling is performed here;
#   - no driver-only filter is applied;
#   - zero-denominator patients are retained and explicitly flagged;
#   - downstream recall/classification must be NA, not zero, when the
#     full-reference non-ubiquitous event count is zero.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ------------------------------ locate script --------------------------------
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

cat("Script 10: TRACERx100 event/reference QC\n")
cat("=====================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------ output paths ---------------------------------
tracerx_result_dir <- file.path(PATHS$output_dir, "results", "tracerx100")
dir.create(tracerx_result_dir, recursive = TRUE, showWarnings = FALSE)

sample_map_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_primary_sample_map.rds"
)
primary_mutation_path <- file.path(
  PATHS$intermediate_dir,
  "tracerx100_primary_mutations_raw.rds"
)

for (p in c(sample_map_path, primary_mutation_path)) {
  if (!file.exists(p)) {
    stop("Missing validated Script 09 input:\n", p)
  }
}

# ------------------------------- constants -----------------------------------
PROTEIN_CLASSES <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Splice_Site",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Translation_Start_Site",
  "Nonstop_Mutation",
  "In_Frame_Ins",
  "In_Frame_Del"
)

# The supplied canonical archive has no In_Frame_Del rows, but the category is
# kept in the scientific definition for symmetry with a HIGH/MODERATE-like set.
EXPECTED <- list(
  n_patients = 100L,
  n_primary_regions = 323L,
  n_selected_sample_event_rows = 85185L,
  n_patient_specific_events = 31919L,
  n_nonubiquitous_events = 8946L,
  n_private_events = 5484L,
  n_zero_nonubiquitous_patients = 4L,
  zero_nonubiquitous_patients = c(
    "CRUK0040", "CRUK0059", "CRUK0061", "CRUK0090"
  ),
  n_regionsum_called_pairs = 85185L,
  n_regionsum_absent_pairs = 18547L,
  minimum_regionsum_denominator = 30L,
  maximum_absent_regionsum_numerator = 0L,
  minimum_called_regionsum_numerator = 1L
)

# ------------------------------- helpers -------------------------------------
safe_int <- function(x, context) {
  out <- suppressWarnings(as.integer(x))
  if (anyNA(out)) {
    stop("Non-integer value while parsing ", context, ".")
  }
  out
}

make_event_key <- function(chr, pos, ref, alt) {
  paste(
    as.character(chr),
    as.character(pos),
    as.character(ref),
    as.character(alt),
    sep = ":"
  )
}

expand_regionsum <- function(patient_id, event_key, region_sum) {
  if (is.na(region_sum) || !nzchar(region_sum)) {
    stop(
      "Missing RegionSum for ", patient_id,
      " / ", event_key
    )
  }

  pieces <- strsplit(region_sum, ";", fixed = TRUE)[[1L]]
  pieces <- pieces[nzchar(pieces)]

  if (length(pieces) == 0L) {
    stop("Empty RegionSum after parsing.")
  }

  labels <- sub(":.*$", "", pieces)
  pairs <- sub("^[^:]*:", "", pieces)

  # This script only needs primary tumour labels R#.
  keep <- grepl("^R[0-9]+$", labels)

  labels <- labels[keep]
  pairs <- pairs[keep]

  if (length(labels) == 0L) {
    stop(
      "No primary R# label in RegionSum for ",
      patient_id, " / ", event_key
    )
  }

  if (anyDuplicated(labels) > 0L) {
    stop(
      "Duplicate R# label in RegionSum for ",
      patient_id, " / ", event_key
    )
  }

  split_pair <- strsplit(pairs, "/", fixed = TRUE)

  if (any(lengths(split_pair) != 2L)) {
    stop(
      "Unexpected RegionSum count-pair format for ",
      patient_id, " / ", event_key
    )
  }

  numerator <- vapply(
    split_pair,
    function(z) safe_int(z[[1L]], "RegionSum numerator"),
    integer(1)
  )
  denominator <- vapply(
    split_pair,
    function(z) safe_int(z[[2L]], "RegionSum denominator"),
    integer(1)
  )

  data.table(
    patient_id = patient_id,
    event_key = event_key,
    region_number = safe_int(
      sub("^R", "", labels),
      "RegionSum R# label"
    ),
    regionsum_numerator = numerator,
    regionsum_denominator = denominator
  )
}

# ------------------------------- 1. load --------------------------------------
sample_map <- as.data.table(readRDS(sample_map_path))
mut <- as.data.table(readRDS(primary_mutation_path))

required_sample_cols <- c(
  "patient_id", "sample_id", "region_number"
)
required_mut_cols <- c(
  "patient_id",
  "mutation_sample_id",
  "sample_origin",
  "Chromosome",
  "Start_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Variant_Classification",
  "Consequence",
  "RegionSum"
)

if (!all(required_sample_cols %in% names(sample_map))) {
  stop(
    "Primary sample map lacks required columns: ",
    paste(setdiff(required_sample_cols, names(sample_map)), collapse = ", ")
  )
}
if (!all(required_mut_cols %in% names(mut))) {
  stop(
    "Primary mutation object lacks required columns: ",
    paste(setdiff(required_mut_cols, names(mut)), collapse = ", ")
  )
}

sample_map[
  ,
  `:=`(
    patient_id = as.character(patient_id),
    sample_id = as.character(sample_id),
    region_number = as.integer(region_number)
  )
]

mut[
  ,
  `:=`(
    patient_id = as.character(patient_id),
    mutation_sample_id = as.character(mutation_sample_id),
    Variant_Classification = as.character(Variant_Classification),
    Consequence = as.character(Consequence)
  )
]

if (nrow(sample_map) != EXPECTED$n_primary_regions ||
    uniqueN(sample_map$patient_id) != EXPECTED$n_patients) {
  stop("Script 09 primary sample-map count QC FAILED.")
}

if (any(mut$sample_origin != "primary_tumour_region")) {
  stop(
    "Primary mutation input unexpectedly contains non-primary samples."
  )
}

# Every MAF sample must belong to the canonical primary sample map.
if (!all(unique(mut$mutation_sample_id) %in% sample_map$sample_id)) {
  stop("Primary mutation object contains an unknown sample ID.")
}

# ------------------------ 2. consequence definition ---------------------------
class_counts_all <- mut[
  ,
  .N,
  by = Variant_Classification
][order(-N, Variant_Classification)]

excluded_splice_region_n <- mut[
  Variant_Classification == "Splice_Region",
  .N
]

# Consequence-level semantic mapping.
#
# IMPORTANT:
# Variant_Classification is a MAF-facing label. For cross-cohort comparability
# with HCC VEP IMPACT in {HIGH, MODERATE}, we require BOTH:
#   (i) an allowed MAF Variant_Classification, and
#   (ii) a compatible SO/VEP Consequence.
#
# This avoids silently treating an inconsistent MAF label as protein-altering.
consequence_rules <- list(
  Missense_Mutation = "missense_variant",
  Nonsense_Mutation = "stop_gained",
  Splice_Site = "splice_acceptor_variant|splice_donor_variant",
  Frame_Shift_Del = "frameshift_variant",
  Frame_Shift_Ins = "frameshift_variant",
  Translation_Start_Site = "start_lost",
  Nonstop_Mutation = "stop_lost",
  In_Frame_Ins = "inframe_insertion",
  In_Frame_Del = "inframe_deletion"
)

candidate <- mut[
  Variant_Classification %in% PROTEIN_CLASSES
]

candidate[
  ,
  semantic_match := FALSE
]

for (vc in PROTEIN_CLASSES) {
  pat <- consequence_rules[[vc]]
  candidate[
    Variant_Classification == vc,
    semantic_match := grepl(pat, Consequence)
  ]
}

# Audit all classification/consequence mismatches rather than weakening the
# semantic rule. In the supplied archive this identifies four regional rows of
# one patient-specific OR2L2 event labelled Missense_Mutation even though its
# Consequence is coding_sequence_variant,5_prime_UTR_variant and HGVSp is empty.
semantic_exclusions <- candidate[
  semantic_match == FALSE
]

selected <- candidate[
  semantic_match == TRUE
]

consequence_qc_rows <- list()
cq_i <- 0L

for (vc in PROTEIN_CLASSES) {
  d <- candidate[Variant_Classification == vc]
  pat <- consequence_rules[[vc]]

  if (nrow(d) == 0L) {
    cq_i <- cq_i + 1L
    consequence_qc_rows[[cq_i]] <- data.table(
      Variant_Classification = vc,
      n_candidate_rows = 0L,
      expected_consequence_pattern = pat,
      n_semantic_matches = 0L,
      n_semantic_exclusions = 0L,
      note = "category absent from supplied archive"
    )
    next
  }

  cq_i <- cq_i + 1L
  consequence_qc_rows[[cq_i]] <- data.table(
    Variant_Classification = vc,
    n_candidate_rows = nrow(d),
    expected_consequence_pattern = pat,
    n_semantic_matches = sum(d$semantic_match),
    n_semantic_exclusions = sum(!d$semantic_match),
    note = ""
  )
}

consequence_semantic_qc <- rbindlist(consequence_qc_rows)

# Archive-specific fingerprint of the only semantic mismatch.
# This is deliberately explicit so a future DataHub update cannot silently
# change the protein-altering definition.
if (nrow(semantic_exclusions) != 4L ||
    uniqueN(semantic_exclusions$patient_id) != 1L ||
    unique(semantic_exclusions$patient_id) != "CRUK0055" ||
    uniqueN(
      make_event_key(
        semantic_exclusions$Chromosome,
        semantic_exclusions$Start_Position,
        semantic_exclusions$Reference_Allele,
        semantic_exclusions$Tumor_Seq_Allele2
      )
    ) != 1L ||
    uniqueN(semantic_exclusions$Variant_Classification) != 1L ||
    unique(semantic_exclusions$Variant_Classification) !=
      "Missense_Mutation" ||
    uniqueN(semantic_exclusions$Consequence) != 1L ||
    unique(semantic_exclusions$Consequence) !=
      "coding_sequence_variant,5_prime_UTR_variant") {

  if (DEBUG) {
    fwrite(
      semantic_exclusions,
      file.path(
        PATHS$qc_dir,
        "10_tracerx100_semantic_exclusions_unexpected.tsv"
      ),
      sep = "\t"
    )
  }
  stop(
    "Unexpected Variant_Classification/Consequence semantic mismatch pattern."
  )
}

selected_class_counts <- selected[
  ,
  .N,
  by = Variant_Classification
][order(-N, Variant_Classification)]

# -------------------------- 3. exact event presence ---------------------------
if (anyNA(selected$Chromosome) ||
    anyNA(selected$Start_Position) ||
    anyNA(selected$Reference_Allele) ||
    anyNA(selected$Tumor_Seq_Allele2)) {
  stop("Selected event set contains a missing exact-variant component.")
}

selected[
  ,
  event_key := make_event_key(
    Chromosome,
    Start_Position,
    Reference_Allele,
    Tumor_Seq_Allele2
  )
]

presence <- unique(
  selected[
    ,
    .(
      patient_id,
      sample_id = mutation_sample_id,
      event_key
    )
  ]
)

if (nrow(presence) != nrow(selected)) {
  stop(
    "Selected MAF contains duplicate sample x exact-event rows ",
    "after protein-altering filtering."
  )
}

if (nrow(presence) != EXPECTED$n_selected_sample_event_rows) {
  stop(
    "Selected sample-event fingerprint mismatch: expected ",
    EXPECTED$n_selected_sample_event_rows,
    ", observed ", nrow(presence), "."
  )
}

# Verify patient assignment against the canonical sample map.
presence_check <- merge(
  presence,
  sample_map[, .(sample_id, patient_id_from_map = patient_id)],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(presence_check$patient_id_from_map) ||
    any(presence_check$patient_id !=
        presence_check$patient_id_from_map)) {
  stop("Sample -> patient mapping mismatch in selected event presence.")
}

# ----------------------------- 4. reference -----------------------------------
patient_n <- sample_map[
  ,
  .(n_regions = .N),
  by = patient_id
]

event_reference <- presence[
  ,
  .(
    full_count = uniqueN(sample_id)
  ),
  by = .(patient_id, event_key)
]

event_reference <- merge(
  event_reference,
  patient_n,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(event_reference$n_regions)) {
  stop("Missing primary-region count while building event reference.")
}

if (any(event_reference$full_count < 1L) ||
    any(event_reference$full_count > event_reference$n_regions)) {
  stop("Invalid full_count in event reference.")
}

event_reference[
  ,
  `:=`(
    is_ubiquitous = full_count == n_regions,
    is_nonubiquitous = full_count > 0L & full_count < n_regions,
    is_private = full_count == 1L
  )
]

if (nrow(event_reference) != EXPECTED$n_patient_specific_events) {
  stop(
    "Patient-specific event fingerprint mismatch: expected ",
    EXPECTED$n_patient_specific_events,
    ", observed ", nrow(event_reference), "."
  )
}
if (sum(event_reference$is_nonubiquitous) !=
    EXPECTED$n_nonubiquitous_events) {
  stop("Non-ubiquitous event-count fingerprint mismatch.")
}
if (sum(event_reference$is_private) != EXPECTED$n_private_events) {
  stop("Private-event count fingerprint mismatch.")
}

# Annotation must not vary across regional calls of the same patient event.
annotation_consistency <- selected[
  ,
  .(
    n_variant_classifications = uniqueN(Variant_Classification)
  ),
  by = .(patient_id, event_key)
]

mixed_annotation <- annotation_consistency[
  n_variant_classifications != 1L
]

if (nrow(mixed_annotation) > 0L) {
  stop(
    "At least one patient-specific event has inconsistent ",
    "Variant_Classification across primary regions."
  )
}

# ---------------------- 5. zero-denominator audit -----------------------------
patient_reference_summary <- event_reference[
  ,
  .(
    n_events_reference = .N,
    n_ubiquitous_reference = sum(is_ubiquitous),
    n_nonubiquitous_reference = sum(is_nonubiquitous),
    n_private_reference = sum(is_private)
  ),
  by = patient_id
]

patient_reference_summary <- merge(
  patient_n,
  patient_reference_summary,
  by = "patient_id",
  all.x = TRUE,
  sort = TRUE
)

for (cc in c(
  "n_events_reference",
  "n_ubiquitous_reference",
  "n_nonubiquitous_reference",
  "n_private_reference"
)) {
  set(
    patient_reference_summary,
    i = which(is.na(patient_reference_summary[[cc]])),
    j = cc,
    value = 0L
  )
}

zero_nonubiq <- patient_reference_summary[
  n_nonubiquitous_reference == 0L
]

if (nrow(zero_nonubiq) !=
    EXPECTED$n_zero_nonubiquitous_patients ||
    !identical(
      sort(zero_nonubiq$patient_id),
      sort(EXPECTED$zero_nonubiquitous_patients)
    )) {
  if (DEBUG) {
    fwrite(
      zero_nonubiq,
      file.path(
        PATHS$qc_dir,
        "10_tracerx100_zero_nonubiquitous_unexpected.tsv"
      ),
      sep = "\t"
    )
  }
  stop("Zero-nonubiquitous-patient fingerprint mismatch.")
}

# -------------------------- 6. RegionSum QC -----------------------------------
# One RegionSum string per patient-specific selected event.
event_regionsum <- selected[
  ,
  .(
    n_regionsum_strings = uniqueN(RegionSum),
    RegionSum = RegionSum[[1L]]
  ),
  by = .(patient_id, event_key)
]

if (any(event_regionsum$n_regionsum_strings != 1L)) {
  stop("Inconsistent RegionSum strings within a selected patient event.")
}

regionsum_expanded <- rbindlist(
  lapply(
    seq_len(nrow(event_regionsum)),
    function(ii) {
      expand_regionsum(
        event_regionsum$patient_id[[ii]],
        event_regionsum$event_key[[ii]],
        event_regionsum$RegionSum[[ii]]
      )
    }
  ),
  use.names = TRUE
)

# Every selected patient event must have exactly one RegionSum pair for every
# primary -R# region in that patient.
expected_event_region_grid <- merge(
  event_reference[, .(patient_id, event_key)],
  sample_map[
    ,
    .(
      patient_id,
      sample_id,
      region_number
    )
  ],
  by = "patient_id",
  allow.cartesian = TRUE,
  sort = FALSE
)

event_region_qc <- merge(
  expected_event_region_grid,
  regionsum_expanded,
  by = c("patient_id", "event_key", "region_number"),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(event_region_qc) != nrow(expected_event_region_grid)) {
  stop("RegionSum expansion changed expected event-region row count.")
}

if (anyNA(event_region_qc$regionsum_numerator) ||
    anyNA(event_region_qc$regionsum_denominator)) {
  stop(
    "At least one selected patient-event x primary-region pair ",
    "is absent from RegionSum."
  )
}

# Add called-presence flag using the public MAF rows.
called_key <- copy(presence)
called_key[, called_presence := TRUE]

event_region_qc <- merge(
  event_region_qc,
  called_key,
  by = c("patient_id", "sample_id", "event_key"),
  all.x = TRUE,
  sort = FALSE
)

event_region_qc[
  is.na(called_presence),
  called_presence := FALSE
]

called_qc <- event_region_qc[called_presence == TRUE]
absent_qc <- event_region_qc[called_presence == FALSE]

regionsum_summary <- data.table(
  item = c(
    "called_event_region_pairs",
    "absent_event_region_pairs",
    "minimum_denominator_all_pairs",
    "minimum_numerator_called_pairs",
    "maximum_numerator_absent_pairs"
  ),
  observed = c(
    nrow(called_qc),
    nrow(absent_qc),
    min(event_region_qc$regionsum_denominator),
    min(called_qc$regionsum_numerator),
    max(absent_qc$regionsum_numerator)
  ),
  expected = c(
    EXPECTED$n_regionsum_called_pairs,
    EXPECTED$n_regionsum_absent_pairs,
    EXPECTED$minimum_regionsum_denominator,
    EXPECTED$minimum_called_regionsum_numerator,
    EXPECTED$maximum_absent_regionsum_numerator
  )
)
regionsum_summary[, pass := observed == expected]

if (!all(regionsum_summary$pass)) {
  if (DEBUG) {
    fwrite(
      regionsum_summary,
      file.path(
        PATHS$qc_dir,
        "10_tracerx100_regionsum_selected_event_failure.tsv"
      ),
      sep = "\t"
    )
  }
  stop("Selected-event RegionSum fingerprint QC FAILED.")
}

# ---------------------------- 7. distribution QC ------------------------------
occupancy_counts <- event_reference[
  ,
  .N,
  by = full_count
][order(full_count)]

patient_event_burden_summary <- patient_reference_summary[
  ,
  .(
    n_patients = .N,
    min_events = min(n_events_reference),
    q1_events = as.numeric(
      quantile(n_events_reference, 0.25, names = FALSE)
    ),
    median_events = median(n_events_reference),
    q3_events = as.numeric(
      quantile(n_events_reference, 0.75, names = FALSE)
    ),
    max_events = max(n_events_reference),
    min_nonubiquitous = min(n_nonubiquitous_reference),
    median_nonubiquitous = median(n_nonubiquitous_reference),
    max_nonubiquitous = max(n_nonubiquitous_reference)
  )
]

# ------------------------------- 8. save --------------------------------------
saveRDS(
  presence,
  file.path(
    PATHS$intermediate_dir,
    "tracerx100_event_presence_protein_altering.rds"
  )
)

saveRDS(
  event_reference,
  file.path(
    PATHS$intermediate_dir,
    "tracerx100_reference_events_protein_altering.rds"
  )
)

saveRDS(
  patient_reference_summary,
  file.path(
    PATHS$intermediate_dir,
    "tracerx100_patient_reference_summary_protein_altering.rds"
  )
)

# The complete event x primary-region RegionSum audit is only ~104k rows and is
# useful for transparent downstream callability checks.

if (DEBUG) {
  fwrite(
    selected_class_counts,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_selected_variant_classification_counts.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    semantic_exclusions,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_semantic_excluded_rows.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    consequence_semantic_qc,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_consequence_semantic_qc.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    patient_reference_summary,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_patient_reference_summary.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    zero_nonubiq,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_zero_nonubiquitous_patients.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    occupancy_counts,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_event_occupancy_counts.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    regionsum_summary,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_regionsum_selected_event_qc.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    patient_event_burden_summary,
    file.path(
      PATHS$qc_dir,
      "10_tracerx100_patient_event_burden_summary.tsv"
    ),
    sep = "\t"
  )
}

# ------------------------ 9. human-readable summary ----------------------------
elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  "Script 10: TRACERx100 event/reference QC",
  "=====================================================",
  "",
  "Canonical cohort:",
  sprintf(
    "  %d patients / %d primary tumour regions",
    uniqueN(sample_map$patient_id),
    nrow(sample_map)
  ),
  "",
  "Cross-cohort protein-altering definition:",
  paste0("  included: ", paste(PROTEIN_CLASSES, collapse = ", ")),
  sprintf(
    "  Splice_Region excluded: %d primary-region MAF rows",
    excluded_splice_region_n
  ),
  sprintf(
    "  Variant_Classification candidates excluded by semantic Consequence check: %d rows",
    nrow(semantic_exclusions)
  ),
  "  excluded rows correspond to one CRUK0055 event labelled Missense_Mutation",
  "  but annotated as coding_sequence_variant,5_prime_UTR_variant with no protein consequence",
  "  all retained selected rows satisfy the prespecified Consequence mapping: PASS",
  "",
  "Exact-event reference:",
  sprintf(
    "  sample x event called presences: %d",
    nrow(presence)
  ),
  sprintf(
    "  patient-specific exact variants: %d",
    nrow(event_reference)
  ),
  sprintf(
    "  non-ubiquitous events: %d",
    sum(event_reference$is_nonubiquitous)
  ),
  sprintf(
    "  private events: %d",
    sum(event_reference$is_private)
  ),
  "  mixed regional Variant_Classification for same event: 0",
  "",
  "Zero-denominator audit:",
  sprintf(
    "  patients with zero full-reference non-ubiquitous protein-altering events: %d",
    nrow(zero_nonubiq)
  ),
  paste0(
    "  IDs: ",
    paste(sort(zero_nonubiq$patient_id), collapse = ", ")
  ),
  "  these patients MUST have recall/classification = NA downstream, not 0",
  "",
  "RegionSum selected-event audit:",
  sprintf(
    "  called patient-event x region pairs: %d",
    nrow(called_qc)
  ),
  sprintf(
    "  absent patient-event x region pairs: %d",
    nrow(absent_qc)
  ),
  sprintf(
    "  minimum RegionSum denominator across all pairs: %d",
    min(event_region_qc$regionsum_denominator)
  ),
  sprintf(
    "  minimum RegionSum numerator among called pairs: %d",
    min(called_qc$regionsum_numerator)
  ),
  sprintf(
    "  maximum RegionSum numerator among absent pairs: %d",
    max(absent_qc$regionsum_numerator)
  ),
  "  every selected event has an R# RegionSum pair for every primary region: PASS",
  "",
  "No downsampling was performed.",
  "No driver-only filter was applied.",
  "Proceed to exhaustive TRACERx100 downsampling only if every QC above is PASS.",
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed)
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "10_tracerx100_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
