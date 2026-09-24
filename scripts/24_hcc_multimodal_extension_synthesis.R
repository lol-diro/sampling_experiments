# ==============================================================================
# 24_hcc_multimodal_extension_synthesis.R
#
#
# Integrates only validated/frozen modules:
#   1) Protein-altering SNV/indel comparator
#   2) 48-gene + DICER1 + canonical TERT driver mutations
#   3) Six recurrent broad/arm-level CNV drivers from mmc4.xlsx
#   4) Seven recurrent focal/cytoband CNV events from mmc4.xlsx
#   5) Published 50-Hallmark GSVA RNA pathway landscape
#
# Outputs
#   results/extensions/integration/
#     24_genomic_detection_curves.tsv
#     24_genomic_classification_curves.tsv
#     24_genomic_k4_paired_vs_protein.tsv
#     24_rna_pathway_curves.tsv
#     24_multimodal_key_findings.tsv
#     24_multimodal_extension_summary.txt
#
#   figures/
#     Figure5_multimodal_sampling_extensions.pdf
#     Figure5_multimodal_sampling_extensions.png
#
# Figure panels
#   A. Detection of non-ubiquitous genomic events
#   B. Correct heterogeneity classification of genomic events
#   C. Preservation of the RNA Hallmark pathway landscape
#   D. RNA standardized pathway-landscape / pathway-ITH errors
#
# No new inference is performed here. All confidence intervals and paired
# comparisons are read from the already validated exact-bootstrap outputs.
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

cat("Script 24: HCC multimodal synthesis\n")
cat("=======================================================\n\n")

t0 <- proc.time()[["elapsed"]]

# ------------------------------- constants -----------------------------------
PRIMARY_K <- 1:4
OPERATIONAL_THRESHOLD <- 0.80
TOL <- 1e-10

EXPECTED_FIXED_GENOMIC_N <- 55L
EXPECTED_DRIVER_INFORMATIVE_N <- 46L
EXPECTED_BROAD_INFORMATIVE_N <- 31L
EXPECTED_FOCAL_INFORMATIVE_N <- 32L
EXPECTED_FIXED_RNA_N <- 39L

# Regression anchors from validated modules.
ANCHORS <- list(
  protein_detection_k4 = 0.843939393939394,
  driver_detection_k4 = 0.866666666666667,
  broad_detection_k4 = 0.933333333333333,
  focal_detection_k4 = 0.816666666666667,
  rna_spearman_k1 = 0.884024009603842,
  rna_spearman_k4 = 0.989454981992797,
  rna_rmse_k1 = 0.382652186118313,
  rna_rmse_k4 = 0.0956630465295782,
  rna_ith_error_k2 = 0.153060217919537,
  rna_ith_error_k4 = 0.0485518763583384
)

# ------------------------------- helpers -------------------------------------
require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Missing required validated output: ",
      label, "\nExpected: ", path
    )
  }

  normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )
}

# Like require_file(), but returns NULL instead of stopping when the file is
# absent. Used only for the driver-mutation inputs, which are optional in
# this repository (Script 15 needs snv_indel.tsv, which is not always
# available) -- see the "driver_available" branch below.
optional_file <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )
}

require_columns <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))

  if (length(miss) > 0L) {
    stop(
      label, " missing required column(s): ",
      paste(miss, collapse = ", ")
    )
  }
}

assert_close <- function(observed, expected, label, tol = TOL) {
  if (length(observed) != 1L ||
      !is.finite(observed) ||
      abs(observed - expected) > tol) {
    stop(
      label, " regression failed.\nExpected: ",
      format(expected, digits = 17),
      "\nObserved: ",
      format(observed, digits = 17)
    )
  }
}

canonicalize_zero <- function(x, tol = 1e-12) {
  x <- as.numeric(x)
  x[is.finite(x) & abs(x) <= tol] <- 0
  x
}

first_tested_k_reaching <- function(dt, threshold = 0.80) {
  hit <- dt[
    metric == "detection" &
      is.finite(median) &
      median >= threshold,
    k
  ]

  if (length(hit) == 0L) {
    NA_integer_
  } else {
    min(as.integer(hit))
  }
}

metric_row <- function(dt, event_set_id, metric_id, kval) {
  x <- dt[
    event_set == event_set_id &
      metric == metric_id &
      k == kval
  ]

  if (nrow(x) != 1L) {
    stop(
      "Expected exactly one row: event_set=",
      event_set_id,
      ", metric=", metric_id,
      ", k=", kval,
      "; observed ", nrow(x), "."
    )
  }

  x
}

fmt_ci <- function(median, lower, upper, digits = 3L) {
  paste0(
    formatC(median, digits = digits, format = "f"),
    " [",
    formatC(lower, digits = digits, format = "f"),
    ", ",
    formatC(upper, digits = digits, format = "f"),
    "]"
  )
}

# ------------------------------ directories ----------------------------------
driver_dir <- file.path(
  PATHS$output_dir,
  "results",
  "extensions",
  "driver"
)

cnv_dir <- file.path(
  PATHS$output_dir,
  "results",
  "extensions",
  "cnv"
)

rna_dir <- file.path(
  PATHS$output_dir,
  "results",
  "extensions",
  "rna"
)

qc_ext_dir <- file.path(
  PATHS$output_dir,
  "results",
  "qc",
  "extensions"
)

integration_dir <- file.path(
  PATHS$output_dir,
  "results",
  "extensions",
  "integration"
)

figure_dir <- PATHS$figures_dir

dir.create(
  integration_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------- inputs --------------------------------------
# Driver-mutation inputs are optional: Script 15 needs snv_indel.tsv, which is
# not always available. When absent, this script still produces Figure 5 with
# the CNV + RNA panels only (see "driver_available" below), and the
# protein-altering comparator -- currently only available as the copy
# embedded in Script 15's own output -- is left out of those panels too,
# rather than silently mixing it with a different endpoint definition.
driver_depth_path <- optional_file(
  file.path(
    driver_dir,
    "15_driver_depth_summary_fixed_n_ge5.tsv"
  )
)

driver_paired_path <- optional_file(
  file.path(
    driver_dir,
    "15_driver_paired_comparisons_fixed_n_ge5.tsv"
  )
)

driver_available <- !is.null(driver_depth_path) &&
  !is.null(driver_paired_path)

broad_depth_path <- require_file(
  file.path(
    cnv_dir,
    "20_broad_cnv_depth_summary_fixed_n_ge5.tsv"
  ),
  "mmc4 broad-CNV depth summary"
)

broad_paired_path <- require_file(
  file.path(
    cnv_dir,
    "20_broad_cnv_vs_protein_paired_fixed_n_ge5.tsv"
  ),
  "mmc4 broad-CNV paired comparison"
)

focal_depth_path <- require_file(
  file.path(
    cnv_dir,
    "19_focal_cnv_depth_summary_fixed_n_ge5.tsv"
  ),
  "mmc4 focal-CNV depth summary"
)

focal_paired_path <- require_file(
  file.path(
    cnv_dir,
    "19_focal_cnv_vs_protein_paired_fixed_n_ge5.tsv"
  ),
  "mmc4 focal-CNV paired comparison"
)

rna_depth_path <- require_file(
  file.path(
    rna_dir,
    "21_rna_depth_summary_fixed_n_ge5.tsv"
  ),
  "RNA Hallmark pathway depth summary"
)

rna_summary_path <- require_file(
  file.path(
    rna_dir,
    "21_rna_pathway_sampling_summary.txt"
  ),
  "RNA Hallmark pathway summary"
)

# ------------------------------- load ----------------------------------------
if (driver_available) {
  driver_depth <- fread(driver_depth_path)
  driver_paired <- fread(driver_paired_path)

  require_columns(
    driver_depth,
    c(
      "event_set",
      "metric",
      "k",
      "n_fixed_cohort",
      "n_patients_with_defined_endpoint",
      "median",
      "ci_lower",
      "ci_upper",
      "bootstrap_method"
    ),
    "driver_depth"
  )

  require_columns(
    driver_paired,
    c(
      "driver_event_set",
      "comparator_event_set",
      "metric",
      "k",
      "n_paired_patients",
      "median_driver_minus_comparator",
      "ci_lower",
      "ci_upper"
    ),
    "driver_paired"
  )
} else {
  cat(
    "\nNOTE: snv_indel.tsv / cnv_arm_level.tsv / Driver_48genes_DICER1.csv\n",
    "not all found, so Script 15 has not run. Figure 5 will be generated\n",
    "with the CNV + RNA panels only. The protein-altering comparator is\n",
    "currently only available as the copy embedded in Script 15's own\n",
    "output, so it is left out of the affected panels rather than mixed in\n",
    "from a different endpoint definition.\n",
    sep = ""
  )
}

broad_depth <- fread(broad_depth_path)
broad_paired <- fread(broad_paired_path)
focal_depth <- fread(focal_depth_path)
focal_paired <- fread(focal_paired_path)
rna_depth <- fread(rna_depth_path)

require_columns(
  broad_depth,
  c(
    "event_set",
    "metric",
    "k",
    "n_fixed_cohort",
    "n_patients_with_defined_endpoint",
    "median",
    "ci_lower",
    "ci_upper",
    "bootstrap_method"
  ),
  "broad_depth"
)

require_columns(
  focal_depth,
  c(
    "event_set",
    "metric",
    "k",
    "n_fixed_cohort",
    "n_patients_with_defined_endpoint",
    "median",
    "ci_lower",
    "ci_upper",
    "bootstrap_method"
  ),
  "focal_depth"
)

require_columns(
  rna_depth,
  c(
    "event_set",
    "metric",
    "direction",
    "k",
    "n_fixed_cohort",
    "n_patients_with_defined_endpoint",
    "median",
    "ci_lower",
    "ci_upper",
    "bootstrap_method"
  ),
  "rna_depth"
)

# ------------------------------- regression QC -------------------------------
protein_set <- "protein_altering_all_genes"
driver_set <- "driver_48_plus_DICER1"
broad_set <- "six_recurrent_broad_cnv_events"
focal_set <- "seven_recurrent_focal_cnv_events"
rna_set <- "published_50_hallmark_gsva"

# Fixed cohort / informative N.
if (driver_available) {
  for (es in c(
    protein_set,
    driver_set
  )) {
    xx <- driver_depth[
      event_set == es
    ]

    if (uniqueN(xx$n_fixed_cohort) != 1L ||
        unique(xx$n_fixed_cohort) != EXPECTED_FIXED_GENOMIC_N) {
      stop(
        "Genomic fixed-cohort regression failed for ",
        es, "."
      )
    }
  }

  if (unique(
    driver_depth[
      event_set == driver_set &
        metric == "detection",
      n_patients_with_defined_endpoint
    ]
  ) != EXPECTED_DRIVER_INFORMATIVE_N) {
    stop("Driver informative-N regression failed.")
  }
}

if (unique(
  broad_depth[
    metric == "detection",
    n_patients_with_defined_endpoint
  ]
) != EXPECTED_BROAD_INFORMATIVE_N) {
  stop("Broad-CNV informative-N regression failed.")
}

if (unique(
  focal_depth[
    metric == "detection",
    n_patients_with_defined_endpoint
  ]
) != EXPECTED_FOCAL_INFORMATIVE_N) {
  stop("Focal-CNV informative-N regression failed.")
}

if (unique(
  rna_depth$n_fixed_cohort
) != EXPECTED_FIXED_RNA_N) {
  stop("RNA fixed-cohort regression failed.")
}

# Exact-bootstrap provenance should be preserved everywhere it is defined.
bootstrap_defined <- rbindlist(
  c(
    if (driver_available) list(
      driver_depth[
        event_set %in% c(
          protein_set,
          driver_set
        )
      ]
    ),
    list(
      broad_depth,
      focal_depth
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

if (any(
  bootstrap_defined$bootstrap_method !=
    "exact_nonparametric_percentile"
)) {
  stop(
    "A defined genomic summary does not use exact-bootstrap."
  )
}

rna_defined <- rna_depth[
  n_patients_with_defined_endpoint > 0L
]

if (any(
  rna_defined$bootstrap_method !=
    "exact_nonparametric_percentile"
)) {
  stop(
    "A defined RNA summary does not use exact-bootstrap."
  )
}

# Canonicalize numerical zero before integration.
if (driver_available) {
  driver_depth[
    ,
    `:=`(
      median = canonicalize_zero(median),
      ci_lower = canonicalize_zero(ci_lower),
      ci_upper = canonicalize_zero(ci_upper)
    )
  ]
}

broad_depth[
  ,
  `:=`(
    median = canonicalize_zero(median),
    ci_lower = canonicalize_zero(ci_lower),
    ci_upper = canonicalize_zero(ci_upper)
  )
]

focal_depth[
  ,
  `:=`(
    median = canonicalize_zero(median),
    ci_lower = canonicalize_zero(ci_lower),
    ci_upper = canonicalize_zero(ci_upper)
  )
]

# Numeric anchors.
if (driver_available) {
  assert_close(
    metric_row(
      driver_depth,
      protein_set,
      "detection",
      4L
    )$median,
    ANCHORS$protein_detection_k4,
    "Protein-altering k4 detection"
  )

  assert_close(
    metric_row(
      driver_depth,
      driver_set,
      "detection",
      4L
    )$median,
    ANCHORS$driver_detection_k4,
    "Driver k4 detection"
  )
}

assert_close(
  broad_depth[
    metric == "detection" &
      k == 4L,
    median
  ],
  ANCHORS$broad_detection_k4,
  "Broad-CNV k4 detection"
)

assert_close(
  focal_depth[
    metric == "detection" &
      k == 4L,
    median
  ],
  ANCHORS$focal_detection_k4,
  "Focal-CNV k4 detection"
)

assert_close(
  rna_depth[
    metric == "landscape_spearman" &
      k == 1L,
    median
  ],
  ANCHORS$rna_spearman_k1,
  "RNA Spearman k1"
)

assert_close(
  rna_depth[
    metric == "landscape_spearman" &
      k == 4L,
    median
  ],
  ANCHORS$rna_spearman_k4,
  "RNA Spearman k4"
)

assert_close(
  rna_depth[
    metric == "landscape_rmse_z" &
      k == 1L,
    median
  ],
  ANCHORS$rna_rmse_k1,
  "RNA RMSE k1"
)

assert_close(
  rna_depth[
    metric == "landscape_rmse_z" &
      k == 4L,
    median
  ],
  ANCHORS$rna_rmse_k4,
  "RNA RMSE k4"
)

assert_close(
  rna_depth[
    metric == "pathway_ith_absolute_error" &
      k == 2L,
    median
  ],
  ANCHORS$rna_ith_error_k2,
  "RNA pathway ITH error k2"
)

assert_close(
  rna_depth[
    metric == "pathway_ith_absolute_error" &
      k == 4L,
    median
  ],
  ANCHORS$rna_ith_error_k4,
  "RNA pathway ITH error k4"
)

# ------------------------ genomic harmonization -------------------------------
broad <- copy(broad_depth)
focal <- copy(focal_depth)

broad[, modality := "Broad CNV drivers"]
focal[, modality := "Focal CNV events"]

broad[, modality_order := 3L]
focal[, modality_order := 4L]

genomic_parts <- list(
  broad,
  focal
)

if (driver_available) {
  protein <- driver_depth[
    event_set == protein_set
  ]

  driver <- driver_depth[
    event_set == driver_set
  ]

  protein[, modality := "Protein-altering SNV/indel"]
  driver[, modality := "Driver mutations"]

  protein[, modality_order := 1L]
  driver[, modality_order := 2L]

  genomic_parts <- c(
    list(protein, driver),
    genomic_parts
  )
}

genomic_curves <- rbindlist(
  genomic_parts,
  use.names = TRUE,
  fill = TRUE
)

genomic_curves <- genomic_curves[
  k %in% PRIMARY_K &
    metric %in% c(
      "detection",
      "heterogeneity_classification"
    ),
  .(
    modality,
    modality_order,
    event_set,
    metric,
    k = as.integer(k),
    n_fixed_cohort =
      as.integer(n_fixed_cohort),
    n_informative =
      as.integer(
        n_patients_with_defined_endpoint
      ),
    median =
      as.numeric(median),
    ci_lower =
      as.numeric(ci_lower),
    ci_upper =
      as.numeric(ci_upper),
    bootstrap_method
  )
]

setorder(
  genomic_curves,
  metric,
  modality_order,
  k
)

genomic_detection <- genomic_curves[
  metric == "detection"
]

genomic_classification <- genomic_curves[
  metric ==
    "heterogeneity_classification"
]

fwrite(
  genomic_detection,
  file.path(
    integration_dir,
    "24_genomic_detection_curves.tsv"
  ),
  sep = "\t"
)

fwrite(
  genomic_classification,
  file.path(
    integration_dir,
    "24_genomic_classification_curves.tsv"
  ),
  sep = "\t"
)

# ------------------------- paired k4 comparisons ------------------------------
if (driver_available) {
  require_columns(
    driver_paired,
    c(
      "driver_event_set",
      "comparator_event_set",
      "metric",
      "k",
      "n_paired_patients",
      "median_driver_minus_comparator",
      "ci_lower",
      "ci_upper"
    ),
    "driver_paired"
  )
}

require_columns(
  broad_paired,
  c(
    "cnv_event_set",
    "comparator_event_set",
    "metric",
    "k",
    "n_paired_patients",
    "median_broad_minus_protein",
    "ci_lower",
    "ci_upper"
  ),
  "broad_paired"
)

require_columns(
  focal_paired,
  c(
    "cnv_event_set",
    "comparator_event_set",
    "metric",
    "k",
    "n_paired_patients",
    "median_focal_minus_protein",
    "ci_lower",
    "ci_upper"
  ),
  "focal_paired"
)

if (driver_available) {
  driver_k4 <- driver_paired[
    driver_event_set == driver_set &
      comparator_event_set ==
        protein_set &
      metric == "detection" &
      k == 4L
  ]

  if (nrow(driver_k4) != 1L) {
    stop(
      "Paired k4 driver comparison regression failed."
    )
  }
}

broad_k4 <- broad_paired[
  cnv_event_set == broad_set &
    comparator_event_set ==
      protein_set &
    metric == "detection" &
    k == 4L
]

focal_k4 <- focal_paired[
  cnv_event_set == focal_set &
    comparator_event_set ==
      protein_set &
    metric == "detection" &
    k == 4L
]

if (nrow(broad_k4) != 1L ||
    nrow(focal_k4) != 1L) {
  stop(
    "Paired k4 genomic comparison regression failed."
  )
}

paired_k4_parts <- list(
  data.table(
    modality = "Broad CNV drivers",
    comparator =
      "Protein-altering SNV/indel",
    k = 4L,
    n_paired =
      broad_k4$n_paired_patients,
    median_difference =
      broad_k4$
        median_broad_minus_protein,
    ci_lower =
      broad_k4$ci_lower,
    ci_upper =
      broad_k4$ci_upper
  ),
  data.table(
    modality = "Focal CNV events",
    comparator =
      "Protein-altering SNV/indel",
    k = 4L,
    n_paired =
      focal_k4$n_paired_patients,
    median_difference =
      focal_k4$
        median_focal_minus_protein,
    ci_lower =
      focal_k4$ci_lower,
    ci_upper =
      focal_k4$ci_upper
  )
)

if (driver_available) {
  paired_k4_parts <- c(
    list(
      data.table(
        modality = "Driver mutations",
        comparator =
          "Protein-altering SNV/indel",
        k = 4L,
        n_paired =
          driver_k4$n_paired_patients,
        median_difference =
          driver_k4$
            median_driver_minus_comparator,
        ci_lower =
          driver_k4$ci_lower,
        ci_upper =
          driver_k4$ci_upper
      )
    ),
    paired_k4_parts
  )
}

paired_k4 <- rbindlist(
  paired_k4_parts,
  use.names = TRUE
)

paired_k4[
  ,
  ci_excludes_zero :=
    is.finite(ci_lower) &
      is.finite(ci_upper) &
      (
        ci_lower > 0 |
          ci_upper < 0
      )
]

fwrite(
  paired_k4,
  file.path(
    integration_dir,
    "24_genomic_k4_paired_vs_protein.tsv"
  ),
  sep = "\t"
)

# ------------------------------ RNA curves -----------------------------------
rna_curves <- rna_depth[
  event_set == rna_set &
    k %in% PRIMARY_K,
  .(
    modality =
      "Published Hallmark GSVA",
    metric,
    direction,
    k = as.integer(k),
    n_fixed_cohort =
      as.integer(n_fixed_cohort),
    n_informative =
      as.integer(
        n_patients_with_defined_endpoint
      ),
    median =
      as.numeric(median),
    ci_lower =
      as.numeric(ci_lower),
    ci_upper =
      as.numeric(ci_upper),
    bootstrap_method
  )
]

setorder(
  rna_curves,
  metric,
  k
)

fwrite(
  rna_curves,
  file.path(
    integration_dir,
    "24_rna_pathway_curves.tsv"
  ),
  sep = "\t"
)

# ---------------------------- key findings table ------------------------------
get_genomic_key <- function(
  modality_name,
  event_set_id,
  source_dt
) {
  det <- source_dt[
    event_set == event_set_id &
      metric == "detection"
  ]

  cls <- source_dt[
    event_set == event_set_id &
      metric ==
        "heterogeneity_classification"
  ]

  k1_det <- det[k == 1L]
  k4_det <- det[k == 4L]
  k4_cls <- cls[k == 4L]

  if (nrow(k1_det) != 1L ||
      nrow(k4_det) != 1L ||
      nrow(k4_cls) != 1L) {
    stop(
      "Key-finding extraction failed for ",
      modality_name, "."
    )
  }

  data.table(
    domain = "Genomic",
    modality = modality_name,
    fixed_cohort_N =
      as.integer(
        k4_det$n_fixed_cohort
      ),
    informative_N =
      as.integer(
        k4_det$
          n_patients_with_defined_endpoint
      ),
    k1_primary =
      as.numeric(k1_det$median),
    k4_primary =
      as.numeric(k4_det$median),
    k4_primary_ci_lower =
      as.numeric(k4_det$ci_lower),
    k4_primary_ci_upper =
      as.numeric(k4_det$ci_upper),
    k4_secondary =
      as.numeric(k4_cls$median),
    first_tested_k_reaching_80pct_detection =
      first_tested_k_reaching(
        det,
        OPERATIONAL_THRESHOLD
      )
  )
}

key_genomic_parts <- list(
  get_genomic_key(
    "Broad CNV drivers",
    broad_set,
    broad_depth
  ),
  get_genomic_key(
    "Focal CNV events",
    focal_set,
    focal_depth
  )
)

if (driver_available) {
  key_genomic_parts <- c(
    list(
      get_genomic_key(
        "Protein-altering SNV/indel",
        protein_set,
        driver_depth
      ),
      get_genomic_key(
        "Driver mutations",
        driver_set,
        driver_depth
      )
    ),
    key_genomic_parts
  )
}

key_genomic <- rbindlist(
  key_genomic_parts,
  use.names = TRUE
)

key_genomic[
  ,
  `:=`(
    primary_endpoint =
      "Detection of non-ubiquitous reference events",
    secondary_endpoint =
      "Correct heterogeneity classification",
    operational_threshold_note =
      "80% is a legacy operational benchmark, not a universal biological optimum."
  )
]

rna_sp1 <- rna_curves[
  metric == "landscape_spearman" &
    k == 1L
]

rna_sp4 <- rna_curves[
  metric == "landscape_spearman" &
    k == 4L
]

rna_rmse1 <- rna_curves[
  metric == "landscape_rmse_z" &
    k == 1L
]

rna_rmse4 <- rna_curves[
  metric == "landscape_rmse_z" &
    k == 4L
]

rna_ith2 <- rna_curves[
  metric ==
    "pathway_ith_absolute_error" &
    k == 2L
]

rna_ith4 <- rna_curves[
  metric ==
    "pathway_ith_absolute_error" &
    k == 4L
]

key_rna <- data.table(
  domain = "Transcriptomic",
  modality = "Published 50-Hallmark GSVA",
  fixed_cohort_N =
    EXPECTED_FIXED_RNA_N,
  informative_N =
    EXPECTED_FIXED_RNA_N,
  primary_endpoint =
    "Spearman preservation of the 50-pathway centroid",
  secondary_endpoint =
    "Standardized centroid RMSE / pathway-ITH absolute error",
  k1_primary =
    rna_sp1$median,
  k4_primary =
    rna_sp4$median,
  k4_primary_ci_lower =
    rna_sp4$ci_lower,
  k4_primary_ci_upper =
    rna_sp4$ci_upper,
  k4_secondary =
    rna_rmse4$median,
  first_tested_k_reaching_80pct_detection =
    NA_integer_,
  operational_threshold_note =
    paste0(
      "No arbitrary RNA recovery threshold imposed. ",
      "RMSE: ",
      formatC(
        rna_rmse1$median,
        digits = 3,
        format = "f"
      ),
      " -> ",
      formatC(
        rna_rmse4$median,
        digits = 3,
        format = "f"
      ),
      "; pathway-ITH absolute error k2 -> k4: ",
      formatC(
        rna_ith2$median,
        digits = 3,
        format = "f"
      ),
      " -> ",
      formatC(
        rna_ith4$median,
        digits = 3,
        format = "f"
      ),
      "."
    )
)

key_findings <- rbindlist(
  list(
    key_genomic,
    key_rna
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  key_findings,
  file.path(
    integration_dir,
    "24_multimodal_key_findings.tsv"
  ),
  sep = "\t"
)

# (Script 23's reconstruction QC is not part of this package; the exclusion
# of transcriptomic_ITH from sampling is stated in the readable summary below
# as documented fact from the manuscript, not re-verified against a QC
# summary here.)

# ------------------------------- figure ---------------------------------------
if (driver_available) {
  modality_levels <- c(
    "Protein-altering SNV/indel",
    "Driver mutations",
    "Broad CNV drivers",
    "Focal CNV events"
  )
  line_types <- c(1, 2, 3, 4)
  point_symbols <- c(16, 17, 15, 18)
} else {
  modality_levels <- c(
    "Broad CNV drivers",
    "Focal CNV events"
  )
  line_types <- c(3, 4)
  point_symbols <- c(15, 18)
}

draw_errorbar <- function(
  x,
  lower,
  upper,
  width = 0.045
) {
  ok <- is.finite(x) &
    is.finite(lower) &
    is.finite(upper)

  if (!any(ok)) {
    return(invisible(NULL))
  }

  graphics::segments(
    x[ok],
    lower[ok],
    x[ok],
    upper[ok]
  )

  graphics::segments(
    x[ok] - width,
    lower[ok],
    x[ok] + width,
    lower[ok]
  )

  graphics::segments(
    x[ok] - width,
    upper[ok],
    x[ok] + width,
    upper[ok]
  )

  invisible(NULL)
}

plot_genomic_panel <- function(metric_id, panel_title) {
  yy <- genomic_curves[
    metric == metric_id
  ]

  graphics::plot(
    NA,
    xlim = c(0.8, 4.2),
    ylim = c(0, 1),
    xaxs = "i",
    yaxs = "i",
    xaxt = "n",
    xlab = "Number of sampled tumour sectors (k)",
    ylab = if (
      metric_id == "detection"
    ) {
      "Expected recovery"
    } else {
      "Expected correct classification"
    },
    main = panel_title
  )

  graphics::axis(
    1,
    at = PRIMARY_K,
    labels = PRIMARY_K
  )

  graphics::abline(
    h = OPERATIONAL_THRESHOLD,
    lty = 3
  )

  for (ii in seq_along(
    modality_levels
  )) {
    md <- modality_levels[[ii]]

    dd <- yy[
      modality == md
    ]

    setorder(
      dd,
      k
    )

    draw_errorbar(
      dd$k,
      dd$ci_lower,
      dd$ci_upper
    )

    graphics::lines(
      dd$k,
      dd$median,
      lty = line_types[[ii]],
      lwd = 1.7
    )

    graphics::points(
      dd$k,
      dd$median,
      pch = point_symbols[[ii]],
      cex = 0.9
    )
  }

  graphics::legend(
    "bottomright",
    legend = modality_levels,
    lty = line_types,
    pch = point_symbols,
    bty = "n",
    cex = 0.72
  )

  graphics::mtext(
    "Dashed line: legacy 80% operational benchmark",
    side = 3,
    line = 0.15,
    cex = 0.62
  )
}

plot_rna_spearman <- function() {
  dd <- rna_curves[
    metric == "landscape_spearman"
  ]

  setorder(dd, k)

  ylim <- range(
    c(
      dd$ci_lower,
      dd$ci_upper,
      0.80,
      1
    ),
    finite = TRUE
  )

  graphics::plot(
    NA,
    xlim = c(0.8, 4.2),
    ylim = ylim,
    xaxs = "i",
    yaxs = "i",
    xaxt = "n",
    xlab = "Number of sampled tumour sectors (k)",
    ylab = "Spearman correlation",
    main = "C. RNA Hallmark landscape preservation"
  )

  graphics::axis(
    1,
    at = PRIMARY_K,
    labels = PRIMARY_K
  )

  draw_errorbar(
    dd$k,
    dd$ci_lower,
    dd$ci_upper
  )

  graphics::lines(
    dd$k,
    dd$median,
    lty = 1,
    lwd = 1.8
  )

  graphics::points(
    dd$k,
    dd$median,
    pch = 16,
    cex = 0.9
  )

  graphics::abline(
    h = 1,
    lty = 3
  )
}

plot_rna_errors <- function() {
  rmse <- rna_curves[
    metric == "landscape_rmse_z"
  ]

  ith <- rna_curves[
    metric ==
      "pathway_ith_absolute_error" &
      is.finite(median)
  ]

  setorder(rmse, k)
  setorder(ith, k)

  ymax <- max(
    c(
      rmse$ci_upper,
      ith$ci_upper
    ),
    na.rm = TRUE
  )

  graphics::plot(
    NA,
    xlim = c(0.8, 4.2),
    ylim = c(0, ymax * 1.08),
    xaxs = "i",
    yaxs = "i",
    xaxt = "n",
    xlab = "Number of sampled tumour sectors (k)",
    ylab = "Standardized error (lower is better)",
    main = "D. RNA landscape and pathway-ITH error"
  )

  graphics::axis(
    1,
    at = PRIMARY_K,
    labels = PRIMARY_K
  )

  draw_errorbar(
    rmse$k,
    rmse$ci_lower,
    rmse$ci_upper
  )

  graphics::lines(
    rmse$k,
    rmse$median,
    lty = 1,
    lwd = 1.8
  )

  graphics::points(
    rmse$k,
    rmse$median,
    pch = 16,
    cex = 0.9
  )

  draw_errorbar(
    ith$k,
    ith$ci_lower,
    ith$ci_upper
  )

  graphics::lines(
    ith$k,
    ith$median,
    lty = 2,
    lwd = 1.8
  )

  graphics::points(
    ith$k,
    ith$median,
    pch = 17,
    cex = 0.9
  )

  graphics::legend(
    "topright",
    legend = c(
      "Pathway centroid RMSE",
      "Pathway-ITH absolute error"
    ),
    lty = c(1, 2),
    pch = c(16, 17),
    bty = "n",
    cex = 0.75
  )

  graphics::mtext(
    "Pathway-ITH error is undefined at k=1",
    side = 3,
    line = 0.15,
    cex = 0.62
  )
}

draw_figure <- function() {
  old_par <- graphics::par(
    no.readonly = TRUE
  )

  on.exit(
    graphics::par(
      old_par
    ),
    add = TRUE
  )

  graphics::par(
    mfrow = c(2, 2),
    mar = c(4.3, 4.4, 3.1, 1.0),
    oma = c(0.4, 0.4, 1.1, 0.4),
    mgp = c(2.4, 0.75, 0),
    tcl = -0.3,
    las = 1,
    cex.axis = 0.82,
    cex.lab = 0.88,
    cex.main = 0.94
  )

  plot_genomic_panel(
    "detection",
    "A. Detection of non-ubiquitous genomic events"
  )

  plot_genomic_panel(
    "heterogeneity_classification",
    "B. Genomic heterogeneity classification"
  )

  plot_rna_spearman()
  plot_rna_errors()

  graphics::mtext(
    if (driver_available) {
      "HCC: multimodal sampling-depth extensions"
    } else {
      "HCC: multimodal sampling-depth extensions (driver-mutation panel not available)"
    },
    outer = TRUE,
    side = 3,
    line = 0.15,
    font = 2,
    cex = 1.05
  )
}

if (!exists("FIGURE_FORMAT") || is.null(FIGURE_FORMAT) ||
    !FIGURE_FORMAT %in% c("pdf", "png", "both")) {
  FIGURE_FORMAT <- "both"
}

pdf_path <- file.path(
  figure_dir,
  "Figure5_multimodal_sampling_extensions.pdf"
)

png_path <- file.path(
  figure_dir,
  "Figure5_multimodal_sampling_extensions.png"
)

if (FIGURE_FORMAT %in% c("pdf", "both")) {
  grDevices::pdf(
    pdf_path,
    width = 9,
    height = 8.2,
    useDingbats = FALSE
  )

  draw_figure()

  grDevices::dev.off()
}

if (FIGURE_FORMAT %in% c("png", "both")) {
  grDevices::png(
    png_path,
    width = 9,
    height = 8.2,
    units = "in",
    res = 300
  )

  draw_figure()

  grDevices::dev.off()
}

# ----------------------------- readable summary -------------------------------
protein_det <- genomic_detection[
  modality ==
    "Protein-altering SNV/indel"
]

driver_det <- genomic_detection[
  modality ==
    "Driver mutations"
]

broad_det <- genomic_detection[
  modality ==
    "Broad CNV drivers"
]

focal_det <- genomic_detection[
  modality ==
    "Focal CNV events"
]

protein_k4 <- protein_det[k == 4L]
driver_k4_det <- driver_det[k == 4L]
broad_k4_det <- broad_det[k == 4L]
focal_k4_det <- focal_det[k == 4L]

protein_k4_cls <- genomic_classification[
  modality ==
    "Protein-altering SNV/indel" &
    k == 4L
]

driver_k4_cls <- genomic_classification[
  modality ==
    "Driver mutations" &
    k == 4L
]

broad_k4_cls <- genomic_classification[
  modality ==
    "Broad CNV drivers" &
    k == 4L
]

focal_k4_cls <- genomic_classification[
  modality ==
    "Focal CNV events" &
    k == 4L
]

rna_sp_k1 <- rna_curves[
  metric == "landscape_spearman" &
    k == 1L
]

rna_sp_k4 <- rna_curves[
  metric == "landscape_spearman" &
    k == 4L
]

rna_rmse_k1 <- rna_curves[
  metric == "landscape_rmse_z" &
    k == 1L
]

rna_rmse_k4 <- rna_curves[
  metric == "landscape_rmse_z" &
    k == 4L
]

rna_ith_k2 <- rna_curves[
  metric ==
    "pathway_ith_absolute_error" &
    k == 2L
]

rna_ith_k4 <- rna_curves[
  metric ==
    "pathway_ith_absolute_error" &
    k == 4L
]

driver_pair <- paired_k4[
  modality == "Driver mutations"
]

broad_pair <- paired_k4[
  modality == "Broad CNV drivers"
]

focal_pair <- paired_k4[
  modality == "Focal CNV events"
]

summary_lines <- c(
  "HCC collaborator-requested extensions: integrated synthesis",
  "===================================================================",
  "",
  "Validated modules included:",
  if (driver_available) "  - Protein-altering SNV/indel comparator",
  if (driver_available) "  - Driver mutations: 48-gene panel + DICER1 + canonical TERT promoter hotspots",
  "  - Six recurrent broad/arm-level CNV drivers from published mmc4 calls",
  "  - Seven recurrent focal/cytoband CNV events from published mmc4 calls",
  "  - Published 50-Hallmark GSVA RNA pathway landscape",
  if (!driver_available) "  - NOT included: driver mutations / protein-altering comparator (snv_indel.tsv not supplied to Script 15)",
  "",
  "Genomic detection of non-ubiquitous reference events:",
  if (driver_available) paste0(
    "  Protein-altering SNV/indel: k=4 ",
    fmt_ci(
      protein_k4$median,
      protein_k4$ci_lower,
      protein_k4$ci_upper
    ),
    "; first tested k >=80% = ",
    first_tested_k_reaching(
      protein_det,
      OPERATIONAL_THRESHOLD
    ),
    "."
  ),
  if (driver_available) paste0(
    "  Driver mutations: k=4 ",
    fmt_ci(
      driver_k4_det$median,
      driver_k4_det$ci_lower,
      driver_k4_det$ci_upper
    ),
    " (N informative=",
    driver_k4_det$n_informative,
    "); first tested k >=80% = ",
    first_tested_k_reaching(
      driver_det,
      OPERATIONAL_THRESHOLD
    ),
    "."
  ),
  paste0(
    "  Broad CNV drivers: k=4 ",
    fmt_ci(
      broad_k4_det$median,
      broad_k4_det$ci_lower,
      broad_k4_det$ci_upper
    ),
    " (N informative=",
    broad_k4_det$n_informative,
    "); first tested k >=80% = ",
    first_tested_k_reaching(
      broad_det,
      OPERATIONAL_THRESHOLD
    ),
    "."
  ),
  paste0(
    "  Focal CNV events: k=4 ",
    fmt_ci(
      focal_k4_det$median,
      focal_k4_det$ci_lower,
      focal_k4_det$ci_upper
    ),
    " (N informative=",
    focal_k4_det$n_informative,
    "); first tested k >=80% = ",
    first_tested_k_reaching(
      focal_det,
      OPERATIONAL_THRESHOLD
    ),
    "."
  ),
  "",
  "Genomic heterogeneity classification at k=4:",
  if (driver_available) paste0(
    "  Protein-altering SNV/indel: ",
    fmt_ci(
      protein_k4_cls$median,
      protein_k4_cls$ci_lower,
      protein_k4_cls$ci_upper
    )
  ),
  if (driver_available) paste0(
    "  Driver mutations: ",
    fmt_ci(
      driver_k4_cls$median,
      driver_k4_cls$ci_lower,
      driver_k4_cls$ci_upper
    )
  ),
  paste0(
    "  Broad CNV drivers: ",
    fmt_ci(
      broad_k4_cls$median,
      broad_k4_cls$ci_lower,
      broad_k4_cls$ci_upper
    )
  ),
  paste0(
    "  Focal CNV events: ",
    fmt_ci(
      focal_k4_cls$median,
      focal_k4_cls$ci_lower,
      focal_k4_cls$ci_upper
    )
  ),
  "",
  "Paired detection difference versus protein-altering SNV/indel at k=4:",
  if (driver_available) paste0(
    "  Driver mutations: +",
    formatC(
      driver_pair$median_difference,
      digits = 3,
      format = "f"
    ),
    " [",
    formatC(
      driver_pair$ci_lower,
      digits = 3,
      format = "f"
    ),
    ", ",
    formatC(
      driver_pair$ci_upper,
      digits = 3,
      format = "f"
    ),
    "] (N=",
    driver_pair$n_paired,
    "); CI includes zero."
  ),
  paste0(
    "  Broad CNV drivers: +",
    formatC(
      broad_pair$median_difference,
      digits = 3,
      format = "f"
    ),
    " [",
    formatC(
      broad_pair$ci_lower,
      digits = 3,
      format = "f"
    ),
    ", ",
    formatC(
      broad_pair$ci_upper,
      digits = 3,
      format = "f"
    ),
    "] (N=",
    broad_pair$n_paired,
    "); CI excludes zero."
  ),
  paste0(
    "  Focal CNV events: +",
    formatC(
      focal_pair$median_difference,
      digits = 3,
      format = "f"
    ),
    " [",
    formatC(
      focal_pair$ci_lower,
      digits = 3,
      format = "f"
    ),
    ", ",
    formatC(
      focal_pair$ci_upper,
      digits = 3,
      format = "f"
    ),
    "] (N=",
    focal_pair$n_paired,
    "); CI includes zero."
  ),
  if (!driver_available) "",
  if (!driver_available) paste0(
    "NOTE: the driver-mutation panel and its protein-altering comparator ",
    "were not generated because snv_indel.tsv (required by Script 15) was ",
    "not supplied. Figure 5 and the tables above show the CNV + RNA panels only."
  ),
  "",
  "RNA Hallmark pathway landscape (fixed n>=5, N=39):",
  paste0(
    "  Landscape Spearman: k=1 ",
    fmt_ci(
      rna_sp_k1$median,
      rna_sp_k1$ci_lower,
      rna_sp_k1$ci_upper
    ),
    " -> k=4 ",
    fmt_ci(
      rna_sp_k4$median,
      rna_sp_k4$ci_lower,
      rna_sp_k4$ci_upper
    ),
    "."
  ),
  paste0(
    "  Standardized centroid RMSE: k=1 ",
    fmt_ci(
      rna_rmse_k1$median,
      rna_rmse_k1$ci_lower,
      rna_rmse_k1$ci_upper
    ),
    " -> k=4 ",
    fmt_ci(
      rna_rmse_k4$median,
      rna_rmse_k4$ci_lower,
      rna_rmse_k4$ci_upper
    ),
    "."
  ),
  paste0(
    "  Pathway-ITH absolute error: k=2 ",
    fmt_ci(
      rna_ith_k2$median,
      rna_ith_k2$ci_lower,
      rna_ith_k2$ci_upper
    ),
    " -> k=4 ",
    fmt_ci(
      rna_ith_k4$median,
      rna_ith_k4$ci_lower,
      rna_ith_k4$ci_upper
    ),
    "; undefined at k=1."
  ),
  "",
  "Interpretation guardrails:",
  paste0(
    "  - The 80% genomic threshold is retained only as the legacy operational ",
    "benchmark; it is not a biological or universal optimum."
  ),
  if (driver_available) paste0(
    "  - Broad CNV drivers are recovered more readily than protein-altering ",
    "SNV/indels in the paired analysis; focal CNVs and driver mutations do not ",
    "show a clear paired k=4 difference because their 95% CIs include zero."
  ) else paste0(
    "  - Broad CNV drivers are recovered more readily than protein-altering ",
    "SNV/indels in the paired analysis; focal CNVs do not show a clear paired ",
    "k=4 difference because their 95% CI includes zero."
  ),
  paste0(
    "  - RNA pathway preservation improves continuously with sampling depth; ",
    "no arbitrary RNA threshold is imposed."
  ),
  paste0(
    "  - Published quantitative transcriptomic_ITH is excluded from sampling ",
    "because its documented Methods definition could not reproduce the ",
    "supplementary values and no tested explicit alternative was exact."
  ),
  paste0(
    "  - All reference sets are the available sectors in each modality, not ",
    "whole-tumour ground truth."
  ),
  "",
  paste0(
    "Integrated figure: ",
    pdf_path
  ),
  paste0(
    "Integrated PNG: ",
    png_path
  ),
  "",
  sprintf(
    "Elapsed time: %.1f seconds",
    proc.time()[["elapsed"]] - t0
  )
)

summary_path <- file.path(
  integration_dir,
  "24_multimodal_extension_summary.txt"
)

writeLines(
  summary_lines,
  summary_path
)

# ------------------------------- provenance ----------------------------------
path_or_placeholder <- function(x) if (is.null(x)) "not_supplied" else x

provenance <- data.table(
  item = c(
    "driver_depth_source",
    "driver_paired_source",
    "broad_depth_source",
    "broad_paired_source",
    "focal_depth_source",
    "focal_paired_source",
    "rna_depth_source",
    "genomic_fixed_cohort_N",
    "rna_fixed_cohort_N",
    "genomic_operational_threshold",
    "new_inference_performed",
    "new_bootstrap_performed",
    "transcriptomic_ITH_downsampling_included"
  ),
  value = c(
    path_or_placeholder(driver_depth_path),
    path_or_placeholder(driver_paired_path),
    broad_depth_path,
    broad_paired_path,
    focal_depth_path,
    focal_paired_path,
    rna_depth_path,
    as.character(
      EXPECTED_FIXED_GENOMIC_N
    ),
    as.character(
      EXPECTED_FIXED_RNA_N
    ),
    as.character(
      OPERATIONAL_THRESHOLD
    ),
    "FALSE",
    "FALSE",
    "FALSE"
  )
)

fwrite(
  provenance,
  file.path(
    integration_dir,
    "24_multimodal_synthesis_provenance.tsv"
  ),
  sep = "\t"
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
