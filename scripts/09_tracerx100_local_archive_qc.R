# ==============================================================================
# 09_tracerx100_local_archive_qc.R
#
# Purpose
#   Audit the USER-PROVIDED public cBioPortal TRACERx 2017 archive before any
#   downsampling. No network download is performed.
#
# Canonical study for the external validation:
#   Jamal-Hanjani et al., NEJM 2017, TRACERx100
#   cBioPortal study ID: nsclc_tracerx_2017
#
# Key distinction established from the archive + original paper:
#   - 323 primary-tumour region samples: SAMPLE_CLASS == "Tumor" and ID -R#
#   -   4 lymph-node region samples:    SAMPLE_CLASS == "Tumor" and ID -LN#
#   - 120 cfDNA/ctDNA samples:          SAMPLE_CLASS == "cfDNA"
#
# IMPORTANT:
#   cBioPortal SAMPLE_TYPE is "Primary" for all 327 tissue-region samples,
#   including the four -LN# samples. Therefore SAMPLE_TYPE alone must NOT be
#   used to define the primary-tumour cohort.
#
# This script:
#   1. Reads/extracts the local .tar.gz archive.
#   2. Validates study identity and file integrity.
#   3. Reconstructs explicit sample-origin classes.
#   4. Reproduces the published 327 = 323 primary + 4 LN tissue regions.
#   5. Verifies that all 100 patients retain >=2 PRIMARY -R# regions after LN
#      exclusion (primary-only range 2-7, median 3).
#   6. Audits the public MAF at sample and exact-variant level.
#   7. Audits the RegionSum field, which carries per-tissue-region read counts.
#   8. Saves canonical raw objects and a primary-region sample map for Script 10.
#
# NO mutation consequence filter is applied here.
# NO downsampling is performed here.
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

cat("Script 09: TRACERx100 local archive QC\n")
cat("===================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ---------------------------- archive argument --------------------------------
args <- commandArgs(trailingOnly = TRUE)

candidate_archives <- c(
  if (length(args) >= 1L) args[[1L]] else character(),
  file.path(PATHS$input_dir, "nsclc_tracerx_2017.tar.gz"),
  file.path(getwd(), "nsclc_tracerx_2017.tar.gz"),
  file.path(script_dir, "nsclc_tracerx_2017.tar.gz")
)
candidate_archives <- unique(candidate_archives)

existing_archives <- candidate_archives[file.exists(candidate_archives)]

if (length(existing_archives) == 0L) {
  stop(
    "TRACERx100 archive not found.\n",
    "Place nsclc_tracerx_2017.tar.gz in your data folder, next to the\n",
    "wgs_samples_matching.RData / somatic_mutations.RData / clinical_data.RData\n",
    "files (i.e. in:\n  ", PATHS$input_dir, "\n)",
    ", or run e.g.:\n",
    "  Rscript 09_tracerx100_local_archive_qc.R nsclc_tracerx_2017.tar.gz"
  )
}

ARCHIVE <- normalizePath(
  existing_archives[[1L]],
  winslash = "/",
  mustWork = TRUE
)

# ------------------------------ directories ----------------------------------
tracerx_root <- file.path(PATHS$output_dir, "external", "tracerx100_2017")
extract_dir <- file.path(tracerx_root, "study_files")
result_dir <- file.path(PATHS$output_dir, "results", "tracerx100")

dir.create(tracerx_root, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (dir.exists(extract_dir)) {
  unlink(extract_dir, recursive = TRUE, force = TRUE)
}
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- constants -----------------------------------
STUDY_ID <- "nsclc_tracerx_2017"

# Version-specific hard QC for this canonical public TRACERx100 archive.
EXPECTED <- list(
  clinical_samples = 447L,
  patients = 100L,
  tumour_tissue_regions = 327L,
  primary_regions = 323L,
  lymph_node_regions = 4L,
  cfdna_samples = 120L,
  primary_patients = 100L,
  primary_min_regions = 2L,
  primary_median_regions = 3,
  primary_max_regions = 7L,
  all_tissue_min_regions = 2L,
  all_tissue_median_regions = 3,
  all_tissue_max_regions = 8L
)

# These MAF row counts are recorded as reproducibility fingerprints of the
# supplied archive. A mismatch is reported but does not define biology.
EXPECTED_MAF_FINGERPRINT <- list(
  total_rows = 209017L,
  primary_region_rows = 205219L,
  lymph_node_rows = 1064L,
  cfdna_rows = 2734L
)

# ------------------------------- helpers -------------------------------------
read_cbio_table <- function(path) {
  if (!file.exists(path)) {
    stop("Missing cBioPortal file: ", path)
  }

  probe <- readLines(path, n = 50L, warn = FALSE)
  non_comment <- which(
    nzchar(trimws(probe)) &
      !startsWith(trimws(probe), "#")
  )

  if (length(non_comment) == 0L) {
    stop("No tabular header found in: ", path)
  }

  fread(
    path,
    sep = "\t",
    skip = non_comment[[1L]] - 1L,
    quote = "",
    na.strings = c(
      "", "NA", "NaN", "[Not Available]", "Not Available",
      "[Not Applicable]", "[Unknown]"
    ),
    showProgress = TRUE
  )
}

find_unique_file <- function(root, basename_target) {
  hits <- list.files(
    root,
    recursive = TRUE,
    full.names = TRUE
  )
  hits <- hits[basename(hits) == basename_target]

  if (length(hits) != 1L) {
    stop(
      "Expected exactly one ", basename_target,
      " after extraction; observed ", length(hits), "."
    )
  }

  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}

parse_meta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  lines <- lines[grepl(":", lines, fixed = TRUE)]

  # IMPORTANT:
  # data.table() has a formal argument named `key`. Using `key = ...` here is
  # therefore interpreted as a request to set a data.table key rather than as
  # creation of a column named "key". Use unambiguous column names.
  data.table(
    meta_key = trimws(sub(":.*$", "", lines)),
    meta_value = trimws(sub("^[^:]*:", "", lines))
  )
}

safe_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

# ---------------------------- 1. archive audit --------------------------------
tar_listing <- tryCatch(
  utils::untar(ARCHIVE, list = TRUE),
  error = function(e) e
)

if (inherits(tar_listing, "error") || length(tar_listing) == 0L) {
  stop(
    "The supplied file is not a readable tar.gz archive.\n",
    if (inherits(tar_listing, "error")) conditionMessage(tar_listing) else ""
  )
}

required_basenames <- c(
  "meta_study.txt",
  "meta_mutations.txt",
  "data_clinical_sample.txt",
  "data_clinical_patient.txt",
  "data_mutations.txt"
)

for (bb in required_basenames) {
  if (!any(basename(tar_listing) == bb)) {
    stop("Required TRACERx100 file absent from archive: ", bb)
  }
}

utils::untar(ARCHIVE, exdir = extract_dir)

meta_study_path <- find_unique_file(extract_dir, "meta_study.txt")
meta_mutations_path <- find_unique_file(extract_dir, "meta_mutations.txt")
clinical_sample_path <- find_unique_file(extract_dir, "data_clinical_sample.txt")
clinical_patient_path <- find_unique_file(extract_dir, "data_clinical_patient.txt")
mutation_path <- find_unique_file(extract_dir, "data_mutations.txt")

study_meta <- parse_meta(meta_study_path)
mutation_meta <- parse_meta(meta_mutations_path)

study_id_observed <- study_meta[
  meta_key == "cancer_study_identifier",
  meta_value
]

if (length(study_id_observed) != 1L ||
    study_id_observed != STUDY_ID) {
  stop(
    "Study identity QC FAILED. Expected ", STUDY_ID,
    "; observed: ", paste(study_id_observed, collapse = ", ")
  )
}

archive_manifest <- data.table(
  archive_path = ARCHIVE,
  archive_bytes = file.info(ARCHIVE)$size,
  archive_md5 = unname(tools::md5sum(ARCHIVE)),
  study_id = study_id_observed,
  analysis_date = as.character(Sys.Date())
)

# ----------------------------- 2. clinical data -------------------------------
cat("Reading TRACERx100 clinical metadata...\n")
clinical_sample <- read_cbio_table(clinical_sample_path)
clinical_patient <- read_cbio_table(clinical_patient_path)

required_sample_cols <- c(
  "SAMPLE_ID", "PATIENT_ID", "SAMPLE_CLASS", "SAMPLE_TYPE",
  "SAMPLE_COLLECTION_TIMEPOINT"
)

missing_sample_cols <- setdiff(
  required_sample_cols,
  names(clinical_sample)
)

if (length(missing_sample_cols) > 0L) {
  stop(
    "Clinical sample schema missing: ",
    paste(missing_sample_cols, collapse = ", ")
  )
}

clinical_sample[
  ,
  `:=`(
    SAMPLE_ID = as.character(SAMPLE_ID),
    PATIENT_ID = as.character(PATIENT_ID),
    SAMPLE_CLASS = as.character(SAMPLE_CLASS),
    SAMPLE_TYPE = as.character(SAMPLE_TYPE)
  )
]

if (anyDuplicated(clinical_sample$SAMPLE_ID) > 0L) {
  stop("Duplicate SAMPLE_ID in clinical sample table.")
}

# Explicit biological source classification.
clinical_sample[
  ,
  sample_origin := fcase(
    SAMPLE_CLASS == "Tumor" &
      grepl("-R[0-9]+$", SAMPLE_ID),
    "primary_tumour_region",

    SAMPLE_CLASS == "Tumor" &
      grepl("-LN[0-9]+$", SAMPLE_ID),
    "lymph_node_region",

    SAMPLE_CLASS == "cfDNA",
    "cfdna",

    default = "unresolved"
  )
]

origin_counts <- clinical_sample[
  ,
  .(
    n_samples = .N,
    n_patients = uniqueN(PATIENT_ID)
  ),
  by = sample_origin
]
setorder(origin_counts, sample_origin)

if (nrow(clinical_sample) != EXPECTED$clinical_samples) {
  stop(
    "Expected ", EXPECTED$clinical_samples,
    " clinical samples; observed ", nrow(clinical_sample), "."
  )
}
if (uniqueN(clinical_sample$PATIENT_ID) != EXPECTED$patients) {
  stop("Expected 100 clinical patients.")
}
if (clinical_sample[sample_origin == "unresolved", .N] != 0L) {
  stop("At least one sample has unresolved origin.")
}
if (clinical_sample[sample_origin == "primary_tumour_region", .N] !=
    EXPECTED$primary_regions) {
  stop("Primary -R# region count QC FAILED.")
}
if (clinical_sample[sample_origin == "lymph_node_region", .N] !=
    EXPECTED$lymph_node_regions) {
  stop("Lymph-node -LN# region count QC FAILED.")
}
if (clinical_sample[sample_origin == "cfdna", .N] !=
    EXPECTED$cfdna_samples) {
  stop("cfDNA sample-count QC FAILED.")
}
if (clinical_sample[
  sample_origin %in% c("primary_tumour_region", "lymph_node_region"),
  .N
] != EXPECTED$tumour_tissue_regions) {
  stop("Total tissue-region count QC FAILED.")
}

# Critical archive-specific observation: cBioPortal SAMPLE_TYPE alone does not
# distinguish the four lymph-node samples.
ln_samples <- clinical_sample[
  sample_origin == "lymph_node_region",
  .(
    SAMPLE_ID,
    PATIENT_ID,
    SAMPLE_CLASS,
    SAMPLE_TYPE,
    SAMPLE_COLLECTION_TIMEPOINT
  )
]

if (any(ln_samples$SAMPLE_TYPE != "Primary")) {
  stop(
    "Unexpected change: at least one -LN# sample is no longer labelled ",
    "SAMPLE_TYPE='Primary'. Re-audit classification logic."
  )
}

# Primary-only patient depth after LN exclusion.
primary_sample_map <- clinical_sample[
  sample_origin == "primary_tumour_region",
  .(
    patient_id = PATIENT_ID,
    sample_id = SAMPLE_ID,
    region_number = as.integer(
      sub("^.*-R([0-9]+)$", "\\1", SAMPLE_ID)
    )
  )
]

if (anyNA(primary_sample_map$region_number)) {
  stop("Could not parse an -R# region number.")
}
if (anyDuplicated(primary_sample_map$sample_id) > 0L) {
  stop("Duplicate primary sample ID.")
}
if (anyDuplicated(
  primary_sample_map[, .(patient_id, region_number)]
) > 0L) {
  stop("Duplicate primary region number within patient.")
}

primary_depth <- primary_sample_map[
  ,
  .(n_primary_regions = .N),
  by = patient_id
]
setorder(primary_depth, patient_id)

all_tissue_depth <- clinical_sample[
  sample_origin %in% c(
    "primary_tumour_region",
    "lymph_node_region"
  ),
  .(n_tissue_regions = .N),
  by = PATIENT_ID
]

primary_depth_qc <- data.table(
  n_patients = nrow(primary_depth),
  n_primary_regions = nrow(primary_sample_map),
  min_regions = min(primary_depth$n_primary_regions),
  median_regions = median(primary_depth$n_primary_regions),
  max_regions = max(primary_depth$n_primary_regions)
)

all_tissue_depth_qc <- data.table(
  n_patients = nrow(all_tissue_depth),
  n_tissue_regions = sum(all_tissue_depth$n_tissue_regions),
  min_regions = min(all_tissue_depth$n_tissue_regions),
  median_regions = median(all_tissue_depth$n_tissue_regions),
  max_regions = max(all_tissue_depth$n_tissue_regions)
)

if (primary_depth_qc$n_patients != EXPECTED$primary_patients ||
    primary_depth_qc$min_regions != EXPECTED$primary_min_regions ||
    primary_depth_qc$median_regions != EXPECTED$primary_median_regions ||
    primary_depth_qc$max_regions != EXPECTED$primary_max_regions) {
  stop(
    "Primary-only region-depth QC FAILED. Expected 100 patients, ",
    "range 2-7, median 3."
  )
}

if (all_tissue_depth_qc$min_regions != EXPECTED$all_tissue_min_regions ||
    all_tissue_depth_qc$median_regions != EXPECTED$all_tissue_median_regions ||
    all_tissue_depth_qc$max_regions != EXPECTED$all_tissue_max_regions) {
  stop(
    "All-tissue region-depth QC FAILED. Expected range 2-8, median 3."
  )
}

# Every TRACERx100 patient remains analyzable after excluding LN samples.
if (any(primary_depth$n_primary_regions < 2L)) {
  stop("At least one patient has <2 primary -R# regions after LN exclusion.")
}

# ----------------------------- 3. mutation MAF --------------------------------
cat("Reading TRACERx100 mutation MAF...\n")
mutations <- read_cbio_table(mutation_path)

required_mut_cols <- c(
  "Chromosome",
  "Start_Position",
  "Reference_Allele",
  "Tumor_Seq_Allele2",
  "Tumor_Sample_Barcode",
  "Variant_Classification",
  "Variant_Type",
  "RegionSum"
)

missing_mut_cols <- setdiff(required_mut_cols, names(mutations))
if (length(missing_mut_cols) > 0L) {
  stop(
    "Mutation MAF schema missing: ",
    paste(missing_mut_cols, collapse = ", ")
  )
}

mutations[
  ,
  mutation_sample_id := as.character(Tumor_Sample_Barcode)
]

sample_origin_map <- clinical_sample[
  ,
  .(
    mutation_sample_id = SAMPLE_ID,
    patient_id = PATIENT_ID,
    sample_origin
  )
]

mutations <- merge(
  mutations,
  sample_origin_map,
  by = "mutation_sample_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(mutations$sample_origin) ||
    anyNA(mutations$patient_id)) {
  stop("At least one mutation sample ID does not map to clinical metadata.")
}

mutation_sample_ids <- unique(mutations$mutation_sample_id)

if (length(mutation_sample_ids) != EXPECTED$clinical_samples) {
  stop(
    "Expected mutation rows from all 447 clinical sample IDs; observed ",
    length(mutation_sample_ids), " unique mutation sample IDs."
  )
}

# Every primary region must have mutation rows.
primary_without_mutations <- setdiff(
  primary_sample_map$sample_id,
  mutation_sample_ids
)

if (length(primary_without_mutations) > 0L) {
  stop(
    "Primary tumour region(s) with zero mutation rows: ",
    paste(primary_without_mutations, collapse = ", ")
  )
}

mutation_origin_counts <- mutations[
  ,
  .(
    n_mutation_rows = .N,
    n_samples = uniqueN(mutation_sample_id),
    n_patients = uniqueN(patient_id)
  ),
  by = sample_origin
]
setorder(mutation_origin_counts, sample_origin)

# Reproducibility fingerprints: report differences but do not reinterpret data.
maf_fingerprint <- data.table(
  item = c(
    "total_rows",
    "primary_region_rows",
    "lymph_node_rows",
    "cfdna_rows"
  ),
  expected = c(
    EXPECTED_MAF_FINGERPRINT$total_rows,
    EXPECTED_MAF_FINGERPRINT$primary_region_rows,
    EXPECTED_MAF_FINGERPRINT$lymph_node_rows,
    EXPECTED_MAF_FINGERPRINT$cfdna_rows
  ),
  observed = c(
    nrow(mutations),
    mutations[sample_origin == "primary_tumour_region", .N],
    mutations[sample_origin == "lymph_node_region", .N],
    mutations[sample_origin == "cfdna", .N]
  )
)
maf_fingerprint[, match := expected == observed]

# ------------------------- 4. exact-variant audit ------------------------------
if (anyNA(mutations$Chromosome) ||
    anyNA(mutations$Start_Position) ||
    anyNA(mutations$Reference_Allele) ||
    anyNA(mutations$Tumor_Seq_Allele2)) {
  stop("Missing component in exact variant key.")
}

mutations[
  ,
  exact_variant_key := paste(
    Chromosome,
    Start_Position,
    Reference_Allele,
    Tumor_Seq_Allele2,
    sep = ":"
  )
]

sample_variant_duplicates <- mutations[
  ,
  .N,
  by = .(mutation_sample_id, exact_variant_key)
][N > 1L]

if (nrow(sample_variant_duplicates) > 0L) {
  stop("Duplicate sample x exact-variant rows detected in public MAF.")
}

# Primary-region mutation table only; this is the candidate source for Script 10.
primary_mutations <- mutations[
  sample_origin == "primary_tumour_region"
]

primary_patient_event_occupancy <- primary_mutations[
  ,
  .(
    n_primary_regions_called = uniqueN(mutation_sample_id)
  ),
  by = .(patient_id, exact_variant_key)
]

# ---------------------------- 5. RegionSum audit -------------------------------
# The public MAF contains RegionSum for tissue-region mutation rows and leaves it
# missing for cfDNA rows. RegionSum is useful as a callability/depth audit, but
# it is NOT used here to redefine called presence.
region_sum_qc <- data.table(
  sample_origin = c(
    "primary_tumour_region",
    "lymph_node_region",
    "cfdna"
  ),
  n_rows = c(
    mutations[sample_origin == "primary_tumour_region", .N],
    mutations[sample_origin == "lymph_node_region", .N],
    mutations[sample_origin == "cfdna", .N]
  ),
  n_region_sum_nonmissing = c(
    mutations[
      sample_origin == "primary_tumour_region" & !is.na(RegionSum),
      .N
    ],
    mutations[
      sample_origin == "lymph_node_region" & !is.na(RegionSum),
      .N
    ],
    mutations[
      sample_origin == "cfdna" & !is.na(RegionSum),
      .N
    ]
  )
)

if (region_sum_qc[
  sample_origin %in% c("primary_tumour_region", "lymph_node_region"),
  any(n_rows != n_region_sum_nonmissing)
]) {
  stop("At least one tissue-region mutation row lacks RegionSum.")
}

if (region_sum_qc[
  sample_origin == "cfdna",
  n_region_sum_nonmissing
] != 0L) {
  stop("Unexpected non-missing RegionSum in cfDNA mutation rows.")
}

# For a patient-specific tissue exact variant, RegionSum should be the same on
# every row representing that event.
tissue_mutations <- mutations[
  sample_origin %in% c(
    "primary_tumour_region",
    "lymph_node_region"
  )
]

region_sum_consistency <- tissue_mutations[
  ,
  .(n_distinct_region_sum = uniqueN(RegionSum)),
  by = .(patient_id, exact_variant_key)
]

if (region_sum_consistency[n_distinct_region_sum != 1L, .N] > 0L) {
  stop("Inconsistent RegionSum strings across rows of the same patient event.")
}

# Audit RegionSum label count against number of tissue regions for each patient.
tissue_count_map <- all_tissue_depth[
  ,
  .(
    patient_id = PATIENT_ID,
    n_tissue_regions
  )
]

region_sum_unique <- unique(
  tissue_mutations[
    ,
    .(
      patient_id,
      RegionSum
    )
  ]
)

region_sum_unique[
  ,
  n_region_sum_labels := lengths(
    strsplit(RegionSum, ";", fixed = TRUE)
  )
]

region_sum_unique <- merge(
  region_sum_unique,
  tissue_count_map,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

if (region_sum_unique[
  n_region_sum_labels != n_tissue_regions,
  .N
] > 0L) {
  stop(
    "RegionSum label count does not match tissue-region count ",
    "for at least one patient/event."
  )
}

# --------------------------- 6. consequence audit ------------------------------
variant_classification_counts <- mutations[
  ,
  .N,
  by = Variant_Classification
]
setorder(variant_classification_counts, -N, Variant_Classification)

variant_type_counts <- mutations[
  ,
  .N,
  by = Variant_Type
]
setorder(variant_type_counts, -N, Variant_Type)

primary_variant_classification_counts <- primary_mutations[
  ,
  .N,
  by = Variant_Classification
]
setorder(
  primary_variant_classification_counts,
  -N,
  Variant_Classification
)

# ------------------------------ 7. save ---------------------------------------
saveRDS(
  primary_sample_map,
  file.path(
    PATHS$intermediate_dir,
    "tracerx100_primary_sample_map.rds"
  )
)
saveRDS(
  primary_mutations,
  file.path(
    PATHS$intermediate_dir,
    "tracerx100_primary_mutations_raw.rds"
  ),
  compress = TRUE
)

# ------------------------------ 8. QC outputs ---------------------------------
if (DEBUG) {
  fwrite(
    archive_manifest,
    file.path(PATHS$qc_dir, "09_tracerx100_archive_manifest.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    origin_counts,
    file.path(PATHS$qc_dir, "09_tracerx100_sample_origin_counts.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    ln_samples,
    file.path(PATHS$qc_dir, "09_tracerx100_lymph_node_samples.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    primary_sample_map,
    file.path(PATHS$qc_dir, "09_tracerx100_primary_sample_map.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    primary_depth,
    file.path(PATHS$qc_dir, "09_tracerx100_primary_regions_per_patient.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    mutation_origin_counts,
    file.path(PATHS$qc_dir, "09_tracerx100_mutation_origin_counts.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    maf_fingerprint,
    file.path(PATHS$qc_dir, "09_tracerx100_maf_fingerprint.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    region_sum_qc,
    file.path(PATHS$qc_dir, "09_tracerx100_regionsum_qc.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    variant_classification_counts,
    file.path(
      PATHS$qc_dir,
      "09_tracerx100_variant_classification_all_samples.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    primary_variant_classification_counts,
    file.path(
      PATHS$qc_dir,
      "09_tracerx100_variant_classification_primary_regions.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    variant_type_counts,
    file.path(PATHS$qc_dir, "09_tracerx100_variant_type_counts.tsv"),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    primary_patient_event_occupancy[
      ,
      .N,
      by = n_primary_regions_called
    ][order(n_primary_regions_called)],
    file.path(
      PATHS$qc_dir,
      "09_tracerx100_primary_event_occupancy_counts.tsv"
    ),
    sep = "\t"
  )
}
if (DEBUG) {
  fwrite(
    study_meta,
    file.path(PATHS$qc_dir, "09_tracerx100_meta_study.tsv"),
    sep = "\t"
  )
}

# ------------------------ 9. human-readable summary ----------------------------
n_primary_events <- nrow(primary_patient_event_occupancy)
n_primary_multi_region_events <- primary_patient_event_occupancy[
  n_primary_regions_called > 1L,
  .N
]

elapsed <- proc.time()[["elapsed"]] - t0

summary_lines <- c(
  "Script 09: TRACERx100 local archive QC",
  "===================================================",
  "",
  sprintf("Archive: %s", ARCHIVE),
  sprintf("Archive size: %.1f MB", file.info(ARCHIVE)$size / 1024^2),
  sprintf("Archive MD5: %s", unname(tools::md5sum(ARCHIVE))),
  sprintf("Study ID: %s", STUDY_ID),
  "",
  "Clinical/sample structure:",
  sprintf("  patients: %d", uniqueN(clinical_sample$PATIENT_ID)),
  sprintf("  all samples: %d", nrow(clinical_sample)),
  sprintf("  primary tumour -R# regions: %d", nrow(primary_sample_map)),
  sprintf("  lymph-node -LN# regions: %d", nrow(ln_samples)),
  sprintf("  cfDNA/ctDNA samples: %d",
          clinical_sample[sample_origin == "cfdna", .N]),
  "",
  "Primary-only region depth after excluding lymph-node samples:",
  sprintf(
    "  N=%d patients; total=%d regions; min=%d; median=%.1f; max=%d",
    primary_depth_qc$n_patients,
    primary_depth_qc$n_primary_regions,
    primary_depth_qc$min_regions,
    primary_depth_qc$median_regions,
    primary_depth_qc$max_regions
  ),
  "  all 100 patients retain >=2 primary regions: PASS",
  "",
  "All tissue regions (primary + LN), reproducing the published TRACERx100 sampling range:",
  sprintf(
    "  total=%d; min=%d; median=%.1f; max=%d",
    all_tissue_depth_qc$n_tissue_regions,
    all_tissue_depth_qc$min_regions,
    all_tissue_depth_qc$median_regions,
    all_tissue_depth_qc$max_regions
  ),
  "",
  "Critical metadata observation:",
  "  all four -LN# samples have SAMPLE_CLASS='Tumor' and SAMPLE_TYPE='Primary'.",
  "  Therefore primary-vs-LN classification MUST use the explicit sample IDs",
  "  (cross-checked against the original 323-primary + 4-LN publication count),",
  "  not SAMPLE_TYPE alone.",
  "",
  "Mutation MAF:",
  sprintf("  rows: %d", nrow(mutations)),
  sprintf("  unique mutation sample IDs: %d", length(mutation_sample_ids)),
  sprintf(
    "  primary-region mutation rows: %d",
    primary_mutations[, .N]
  ),
  sprintf(
    "  patient-specific exact variants in primary regions: %d",
    n_primary_events
  ),
  sprintf(
    "  exact variants called in >1 primary region: %d (%.1f%%)",
    n_primary_multi_region_events,
    100 * n_primary_multi_region_events / n_primary_events
  ),
  "  duplicate sample x exact-variant rows: 0",
  "  every primary region has mutation rows: PASS",
  "",
  "RegionSum audit:",
  "  present on every tissue-region mutation row: PASS",
  "  absent on every cfDNA mutation row: PASS",
  "  identical across rows of the same patient-specific tissue event: PASS",
  "  number of RegionSum region labels matches each patient's tissue-region count: PASS",
  "",
  "MAF reproducibility fingerprint:",
  paste0(
    "  ",
    maf_fingerprint$item,
    ": observed=",
    maf_fingerprint$observed,
    ", expected=",
    maf_fingerprint$expected,
    " -> ",
    ifelse(maf_fingerprint$match, "MATCH", "DIFFERS"),
    collapse = "\n"
  ),
  "",
  "No mutation consequence filter was applied.",
  "No downsampling was performed.",
  "Proceed to Script 10 only after reviewing the primary-region consequence categories.",
  "",
  sprintf("Elapsed time: %.1f seconds", elapsed)
)

writeLines(
  summary_lines,
  file.path(PATHS$qc_dir, "09_tracerx100_qc_summary.txt")
)


cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
