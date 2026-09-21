# ==============================================================================
# Configuration: portable paths, filter thresholds, and frozen QC targets.
# ==============================================================================

options(stringsAsFactors = FALSE)

# ==============================================================================
# USER CONFIGURATION -- edit the two values below for your machine.
# ------------------------------------------------------------------------------


# INPUT_PATH
#   Absolute path to the directory that directly contains the input data
INPUT_PATH <- "/Users/ciro/Documents/R_Projects/data/planet_sampling/input_data"

# OUTPUT_PATH
#   Where results/intermediate/figures are written.
OUTPUT_PATH <- "/Users/ciro/Documents/R_Projects/data/planet_sampling/output"

# DEBUG
#   FALSE (default): results/ contains only what is needed to reproduce the main 
#   results and figures.
#   TRUE: every table behind every check and every sensitivity/robustness
#   analysis is also written, for debugging or a from-scratch
#   re-verification. 

#   Nothing computed differs between the two cases (TRUE or FALSE) only what
#   is written to disk.
DEBUG <- TRUE

# FIGURE_FORMAT
FIGURE_FORMAT <- "both"
# ==============================================================================

resolve_project_root <- function() {
  configured_root <- if (!is.null(INPUT_PATH) &&
      nzchar(INPUT_PATH)) {
    INPUT_PATH
  } else {
    ""
  }

  env_root <- Sys.getenv("INPUT_DIR", unset = "")

  candidates <- unique(
    c(
      if (nzchar(configured_root)) configured_root else character(),
      if (nzchar(env_root)) env_root else character(),
      getwd(),
      file.path(getwd(), "data")
    )
  )

  required_files <- c(
    "wgs_samples_matching.RData",
    "somatic_mutations.RData",
    "clinical_data.RData"
  )

  for (candidate in candidates) {
    expanded <- path.expand(candidate)

    if (dir.exists(expanded) &&
        all(file.exists(file.path(expanded, required_files)))) {
      return(
        normalizePath(
          expanded,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  }

  stop(
    "Could not identify the project root.\n",
    "Set INPUT_PATH near the top of 00_config.R to the directory that\n",
    "directly contains wgs_samples_matching.RData, somatic_mutations.RData\n",
    "and clinical_data.RData.\n",
    "Example:\n",
    "  INPUT_PATH <- \"/path/to/latest/data\"\n",
    "(Alternatively, set the environment variable INPUT_DIR.)"
  )
}

PROJECT_ROOT <- resolve_project_root()

output_dir_override <- if (!is.null(OUTPUT_PATH) &&
    nzchar(OUTPUT_PATH)) {
  OUTPUT_PATH
} else {
  Sys.getenv("OUTPUT_DIR", unset = "")
}

if (nzchar(output_dir_override)) {
  RESOLVED_OUTPUT_DIR <- path.expand(output_dir_override)
  dir.create(RESOLVED_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  RESOLVED_OUTPUT_DIR <- normalizePath(RESOLVED_OUTPUT_DIR, winslash = "/", mustWork = TRUE)
} else {
  RESOLVED_OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
}

PATHS <- list(
  project_root = PROJECT_ROOT,
  input_dir = PROJECT_ROOT,
  output_dir = RESOLVED_OUTPUT_DIR
)

PATHS$qc_dir <- file.path(PATHS$output_dir, "results", "qc")
PATHS$intermediate_dir <- file.path(PATHS$output_dir, "intermediate")
PATHS$figures_dir <- file.path(PATHS$output_dir, "figures")

for (d in unname(unlist(PATHS[
  c(
    "output_dir",
    "qc_dir",
    "intermediate_dir",
    "figures_dir"
  )
]))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

INPUT_FILES <- list(
  sample_map = file.path(
    PATHS$input_dir,
    "wgs_samples_matching.RData"
  ),
  mutations = file.path(
    PATHS$input_dir,
    "somatic_mutations.RData"
  ),
  clinical = file.path(
    PATHS$input_dir,
    "clinical_data.RData"
  )
)

missing_input <- names(INPUT_FILES)[
  !file.exists(unlist(INPUT_FILES))
]

if (length(missing_input) > 0L) {
  stop(
    "Missing input file(s): ",
    paste(missing_input, collapse = ", "),
    "\nExpected input directory: ",
    PATHS$input_dir
  )
}

# ------------------------------------------------------------------------------
# PARAMS -- verified usage across the pipeline (see docs/PIPELINE_MAP.md).
# Entries marked LOAD-BEARING change the scientific results if edited and are
# frozen to the values reported in the manuscript; do not change them to
# reproduce the published analysis. Entries marked INERT are read into a
# local variable somewhere downstream but never actually change any computed
# result (confirmed by exhaustive search across every script) -- safe to
# leave untouched, documented here only for transparency.
# ------------------------------------------------------------------------------
PARAMS <- list(

  # INERT: captured into a `seed`/`BASE_SEED` provenance column in several
  # output tables (05, 07, 08, 12) for audit purposes only. The pipeline uses
  # exact enumeration and exact-bootstrap distributions throughout -- no
  # function in this codebase draws random numbers -- so this value does not
  # affect any computed result.
  seed = 12345L,

  # LOAD-BEARING: the primary HCC mutation filters (Methods -> "Construction
  # of the HCC somatic-variant reference"). Applied in 02_hcc_events.R.
  min_alt_count = 3L,
  min_total_depth = 10L,
  min_vaf = 0.05,

  # COSMETIC / reporting-only: used in 05 and 12 to annotate where a curve
  # crosses this legacy benchmark. The manuscript (Methods -> "Statistical
  # inference") explicitly states this threshold is not used to define an
  # optimal sector number; changing it only changes an annotation, not any
  # endpoint, recall, or confidence interval.
  adequate_recovery_threshold = 0.80,

  # INERT: not read by any script in this repository. Every output-table
  # column that historically carried a Monte-Carlo replicate count is fixed
  # to NA under the exact-bootstrap method (functions_bootstrap_exact.R).
  # Retained only as documentation of the legacy Monte-Carlo fallback.
  bootstrap_replicates = 100000L,

  # INERT: not read by any script in this repository. The k ranges actually
  # analysed (e.g. k<=4 for the fixed n>=5 cohort, k<=5 for the n>=6
  # exploratory cohort) are literal values inside each script, not driven by
  # this parameter.
  main_k_max = 5L
)

EXPECTED_LEGACY <- list(
  n_patients = 123L,
  n_tumor_sectors = 490L,
  sector_count_distribution = c(
    `1` = 1L,
    `2` = 22L,
    `3` = 24L,
    `4` = 21L,
    `5` = 48L,
    `6` = 4L,
    `7` = 2L,
    `11` = 1L
  )
)

# Frozen raw-input fingerprints from the independent end-to-end audit.
EXPECTED_INPUT_MD5 <- c(
  somatic_mutations.RData =
    "30b6aa77fcf4885741cb703ea2a034cc",
  wgs_samples_matching.RData =
    "b48b33d718bde6cc0d9dae2c6bbac033",
  clinical_data.RData =
    "f26c9f0f5f43560d525f56ac03e10526"
)

EXPECTED_TRACERX_ARCHIVE_MD5 <-
  "c4006c662ece2e125b2f05babcac8345"

set.seed(PARAMS$seed)
