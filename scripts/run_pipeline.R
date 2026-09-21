# ==============================================================================
# run_pipeline.R
#
# Pure-R orchestrator for the frozen HCC / TRACERx100 sampling pipeline
# (Figures 1-4, Supplementary Figure 1, and -- when their inputs are
# present -- the Figure 5 driver/CNV/RNA extensions). Runs scripts in order
# inside your current R session using source(), stopping immediately if any
# script raises an error (R's normal error propagation -- no special
# handling needed). This does not change any analysis logic in the numbered
# scripts. Runs unattended, start to finish, with no prompts to answer.
#
# HOW TO RUN (R console, RStudio, or `Rscript run_pipeline.R` -- all work,
# since nothing here waits for you to type anything):
#
#   setwd("/full/path/to/sampling_experiments/scripts")
#   source("run_pipeline.R")
#
# Prerequisites:
#   1. INPUT_PATH must be set inside 00_config.R (see README.md).
#   2. Your working directory (getwd()) must be this scripts/ folder, because
#      every numbered script locates its neighbours (00_config.R,
#      functions_sampling.R, functions_bootstrap_exact.R) relative to it.
#   3. For scripts 09-13 (TRACERx100 arm + manuscript figures): place your
#      local copy of nsclc_tracerx_2017.tar.gz directly in INPUT_PATH,
#      alongside the three .RData files.
#   4. For Figure 5 (scripts 14-15, 19-21, 24): place mmc4.xlsx in INPUT_PATH
#      for the CNV/RNA panels; additionally place Driver_48genes_DICER1.csv,
#      snv_indel.tsv[.zip], and cnv_arm_level.tsv there for the driver panel
#      (optional -- Figure 5 is generated with the CNV+RNA panels only if
#      those three are missing). See README.md for details.
#   5. Run verify_inputs.R first (recommended, not enforced here).
# ==============================================================================

run_step <- function(script_name) {
  cat("\n==============================================================\n")
  cat("Running:", script_name, "\n")
  cat("==============================================================\n")
  source(script_name, echo = FALSE)
}

# --- HCC arm ----------------------------------------------------------------
run_step("01_hcc_metadata_qc.R")

cat("\n>>> NOTE: script 01 reports a legacy QC mismatch as a warning rather\n")
cat(">>> than aborting, so the pipeline continues automatically either way.\n")
cat(">>> After it finishes, check qc/01_qc_summary.txt and confirm both\n")
cat(">>> legacy QC lines read PASS before trusting downstream results.\n")

run_step("02_hcc_events.R")
run_step("04_hcc_exhaustive_sampling.R")
run_step("05_hcc_depth_analysis_exact_bootstrap.R")
run_step("06_hcc_spatial_designs.R")
run_step("07_hcc_spatial_analysis_exact_bootstrap.R")
run_step("08_hcc_protein_altering_sensitivity_exact_bootstrap.R")

# --- TRACERx100 arm + cross-cohort validation + manuscript figures ----------
# PATHS is already in the global environment at this point (every numbered
# script sources 00_config.R with source()'s default local = FALSE, so its
# variables land in .GlobalEnv). The archive is expected in the same data
# folder as the three .RData files.
tracerx_archive_path <- file.path(PATHS$input_dir, "nsclc_tracerx_2017.tar.gz")

if (file.exists(tracerx_archive_path)) {
  run_step("09_tracerx100_local_archive_qc.R")
  run_step("10_tracerx100_event_reference_qc.R")
  run_step("11_tracerx100_exhaustive_sampling.R")
  run_step("12_hcc_tracerx100_external_validation_exact_bootstrap.R")
  run_step("13_generate_manuscript_figures.R")
} else {
  cat("\nnsclc_tracerx_2017.tar.gz not found in the data folder (",
      PATHS$input_dir, "): skipping scripts 09-13.\n", sep = "")
  cat("Place the archive there and re-run source(\"run_pipeline.R\") to include it.\n")
  cat("Script 13 (manuscript figures) needs Script 12's output, so it is\n")
  cat("skipped too until the TRACERx100 arm has been run.\n")
}

cat("\nPipeline completed (Figures 1-4 + Suppl. Fig. 1). Outputs are under\n")
cat("the directory printed as PATHS$output_dir by 00_config.R.\n")

# --- Figure 5 extensions: driver mutations, CNV, RNA pathway ---------------
# Each block runs only if its own input files are present in the data
# folder (INPUT_PATH). Script 14-15 (driver) and 19-21 (CNV/RNA) are fully
# independent of each other. Script 24 (Figure 5 itself) runs whenever
# EITHER side is available: with only CNV/RNA, it produces Figure 5 with
# the CNV + RNA panels only (no driver panel, no protein-altering
# comparator in those panels -- see docs/PIPELINE_MAP.md).

driver_inputs_present <- file.exists(file.path(PATHS$input_dir, "Driver_48genes_DICER1.csv")) &&
  (file.exists(file.path(PATHS$input_dir, "snv_indel.tsv")) ||
     file.exists(file.path(PATHS$input_dir, "snv_indel.tsv.zip"))) &&
  file.exists(file.path(PATHS$input_dir, "cnv_arm_level.tsv"))

if (driver_inputs_present) {
  run_step("14_hcc_driver_qc.R")
  run_step("15_hcc_driver_sampling.R")
} else {
  cat("\nDriver_48genes_DICER1.csv / snv_indel.tsv[.zip] / cnv_arm_level.tsv not\n")
  cat("all found in the data folder: skipping scripts 14-15 (driver mutations).\n")
  cat("Figure 5, if generated below, will have no driver panel and no\n")
  cat("protein-altering comparator in the genomic panels.\n")
}

mmc4_present <- file.exists(file.path(PATHS$input_dir, "mmc4.xlsx"))

if (mmc4_present) {
  run_step("19_hcc_cnv_focal_sampling.R")
  run_step("20_hcc_cnv_broad_sampling_mmc4.R")
  run_step("21_hcc_rna_pathway_sampling.R")
} else {
  cat("\nmmc4.xlsx not found in the data folder: skipping scripts 19-21\n")
  cat("(focal CNV, broad CNV, RNA pathway sampling).\n")
}

if (mmc4_present) {
  run_step("24_hcc_multimodal_extension_synthesis.R")
} else {
  cat("\nScript 24 (Figure 5) needs at least Scripts 19, 20 and 21 to have run\n")
  cat("(mmc4.xlsx), so it is skipped until that file is present. The driver\n")
  cat("arm (14-15) alone is not enough to generate Figure 5.\n")
}

cat("\nFigure 5 extension analyses completed (or skipped as noted above).\n")

# --- Consolidated validation log --------------------------------------------
source("generate_validation_summary.R")
generate_validation_summary(PATHS$output_dir)
