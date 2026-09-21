# ==============================================================================
# 13_generate_manuscript_figures.R
#
# Purpose
#   Generate the frozen manuscript figure set directly from validated outputs
#   of Scripts 05, 07, 08 and 12. No analytical quantities are recomputed from
#   raw genomic data in this script.
#
# Main figures
#   Figure 1  Sampling-design schematic requested by collaborators
#             rows: grid-based dispersed / random-unconstrained example
#             columns: k=3 / k=5, n=9 ordered sectors
#
#   Figure 2  HCC sampling-depth learning curve, fixed n>=5 (N=55)
#             A: non-ubiquitous detection recall
#             B: heterogeneity-classification recall
#             all-filtered primary + protein-altering sensitivity
#
#   Figure 3  HCC spatial-design effect, strict fixed n>=5 (N=44)
#             A-B: paired grid-vs-cluster advantage with 95% bootstrap CI
#             C-D: absolute grid / unconstrained / clustered performance
#
#   Figure 4  Cross-cancer external validation, fixed n>=5
#             HCC N=55; TRACERx100 N=15; k=1..4
#             A: detection learning curves
#             B: heterogeneity-classification learning curves
#             C: paired marginal detection gains
#             D: reference non-ubiquitous/private fractions
#
# Supplementary figure
#   Figure S1 HCC ceiling diagnostic: available cohort vs fixed n>=5 detection
#
# Scientific guardrails
#   - Figure 1 random/unconstrained panels are illustrative only; the real unconstrained
#     analysis exhaustively enumerates all k-subsets.
#   - The 0.80 line is the legacy operational threshold, not a biological or
#     universally optimal threshold.
#   - Spatial effects are within-patient paired effects; positive = dispersed
#     grid better than clustered/adjacent.
#   - TRACERx validates learning-curve behavior, not equality or a universal k.
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

cat("Script 13: HCC / TRACERx100 manuscript figures\n")
cat("========================================================\n\n")

figure_dir <- PATHS$figures_dir
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ validated inputs ------------------------------
hcc_depth_all_path <- file.path(
  PATHS$output_dir, "results", "hcc_depth",
  "05_hcc_depth_core_summary_all_filtered.tsv"
)
hcc_depth_pa_path <- file.path(
  PATHS$output_dir, "results", "hcc_sensitivity",
  "08_hcc_protein_altering_depth_summary.tsv"
)
spatial_effect_path <- file.path(
  PATHS$output_dir, "results", "hcc_spatial",
  "07_hcc_spatial_main_effect_summary.tsv"
)
spatial_absolute_path <- file.path(
  PATHS$output_dir, "results", "hcc_spatial",
  "07_hcc_spatial_absolute_design_summary_all_filtered.tsv"
)
cross_curve_path <- file.path(
  PATHS$output_dir, "results", "cross_cohort",
  "12_hcc_tracerx100_fixed_n_ge5_curve_summary.tsv"
)
cross_gain_path <- file.path(
  PATHS$output_dir, "results", "cross_cohort",
  "12_hcc_tracerx100_fixed_n_ge5_paired_gains.tsv"
)
ref_arch_path <- file.path(
  PATHS$output_dir, "results", "cross_cohort",
  "12_hcc_tracerx100_reference_architecture_comparison.tsv"
)

required <- c(
  hcc_depth_all_path,
  hcc_depth_pa_path,
  spatial_effect_path,
  spatial_absolute_path,
  cross_curve_path,
  cross_gain_path,
  ref_arch_path
)
if (any(!file.exists(required))) {
  stop(
    "Missing validated figure input(s):\n",
    paste0("  ", required[!file.exists(required)], collapse = "\n")
  )
}

hcc_all <- fread(hcc_depth_all_path)
hcc_pa <- fread(hcc_depth_pa_path)
sp_eff <- fread(spatial_effect_path)
sp_abs <- fread(spatial_absolute_path)
cross_curve <- fread(cross_curve_path)
cross_gain <- fread(cross_gain_path)
ref_arch <- fread(ref_arch_path)

# ------------------------------ plotting helpers ------------------------------
cols2 <- hcl.colors(2, palette = "Dark 2")
cols3 <- hcl.colors(3, palette = "Dark 2")

# FIGURE_FORMAT (set in 00_config.R): "pdf", "png", or "both".
if (!exists("FIGURE_FORMAT") || is.null(FIGURE_FORMAT) ||
    !FIGURE_FORMAT %in% c("pdf", "png", "both")) {
  FIGURE_FORMAT <- "both"
}

open_devices <- function(stem, width, height) {
  pdf(
    file.path(figure_dir, paste0(stem, ".pdf")),
    width = width,
    height = height,
    useDingbats = FALSE,
    family = "sans"
  )
}

open_png <- function(stem, width, height, res = 300) {
  png(
    file.path(figure_dir, paste0(stem, ".png")),
    width = width,
    height = height,
    units = "in",
    res = res,
    family = "sans"
  )
}

with_two_devices <- function(stem, width, height, plot_fun) {
  if (FIGURE_FORMAT %in% c("pdf", "both")) {
    open_devices(stem, width, height)
    plot_fun()
    dev.off()
  }

  if (FIGURE_FORMAT %in% c("png", "both")) {
    open_png(stem, width, height)
    plot_fun()
    dev.off()
  }
}

error_bars <- function(x, lo, hi, length = 0.04, ...) {
  # Avoid base-R warnings for exactly zero-length confidence intervals
  # (e.g. heterogeneity-classification recall at k=1 is identically 0).
  keep <- is.finite(x) & is.finite(lo) & is.finite(hi) &
    (abs(hi - lo) > sqrt(.Machine$double.eps))

  if (any(keep)) {
    arrows(
      x0 = x[keep], y0 = lo[keep],
      x1 = x[keep], y1 = hi[keep],
      angle = 90,
      code = 3,
      length = length,
      ...
    )
  }
}

panel_letter <- function(letter) {
  mtext(
    letter,
    side = 3,
    adj = -0.12,
    line = 0.2,
    font = 2,
    cex = 1.15
  )
}

# ------------------------------------------------------------------------------
# FIGURE 1. Sampling-design schematic
# ------------------------------------------------------------------------------
plot_sampling_panel <- function(k, design, selected, panel_title) {
  n <- 9L
  x <- seq_len(n)
  y <- rep(0, n)

  plot(
    x, y,
    type = "n",
    xlim = c(0.4, 9.6),
    ylim = c(-0.8, 0.9),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = panel_title,
    cex.main = 1.05
  )

  segments(1, 0, 9, 0, lwd = 2)

  symbols(
    x, y,
    circles = rep(0.14, n),
    inches = FALSE,
    add = TRUE,
    bg = "white",
    fg = "black"
  )

  symbols(
    x[selected], y[selected],
    circles = rep(0.16, length(selected)),
    inches = FALSE,
    add = TRUE,
    bg = cols2[[1L]],
    fg = "black"
  )

  text(
    x,
    rep(-0.34, n),
    labels = paste0("T", x),
    cex = 0.78
  )

  text(
    5,
    0.55,
    labels = paste0(design, "; k = ", k),
    cex = 0.88,
    font = 2
  )
}

figure1_fun <- function() {
  op <- par(
    mfrow = c(2, 2),
    mar = c(1.7, 1.0, 2.1, 0.7),
    oma = c(2.6, 2.8, 1.2, 0.5),
    xpd = NA
  )
  on.exit(par(op), add = TRUE)

  plot_sampling_panel(
    k = 3,
    design = "Grid-based dispersed",
    selected = c(1, 5, 9),
    panel_title = "A"
  )

  plot_sampling_panel(
    k = 5,
    design = "Grid-based dispersed",
    selected = c(1, 3, 5, 7, 9),
    panel_title = "B"
  )

  # Fixed illustrative random/unconstrained selections. They are deliberately not
  # interpreted as analysis results; exhaustive enumeration is used in data.
  plot_sampling_panel(
    k = 3,
    design = "Random/unconstrained example",
    selected = c(2, 4, 8),
    panel_title = "C"
  )

  plot_sampling_panel(
    k = 5,
    design = "Random/unconstrained example",
    selected = c(1, 2, 5, 6, 9),
    panel_title = "D"
  )

  mtext(
    "Ordered HCC-like tumour sectors (schematic; n = 9)",
    side = 1,
    outer = TRUE,
    line = 1.3,
    cex = 0.9
  )
  mtext(
    "Random panels are illustrative examples; quantitative unconstrained analysis exhaustively enumerates all k-subsets",
    side = 3,
    outer = TRUE,
    line = 0.1,
    cex = 0.82
  )
}

with_two_devices(
  "Figure1_sampling_design_schematic",
  width = 8.0,
  height = 5.6,
  plot_fun = figure1_fun
)

# ------------------------------------------------------------------------------
# FIGURE 2. HCC depth: all-filtered primary + protein-altering sensitivity
# ------------------------------------------------------------------------------
get_hcc_depth <- function(metric_id) {
  a <- hcc_all[
    analysis == "fixed_n_ge5" & metric == metric_id
  ][order(k)]
  p <- hcc_pa[
    analysis == "fixed_n_ge5" & metric == metric_id
  ][order(k)]

  if (nrow(a) != 4L || nrow(p) != 4L ||
      !identical(a$k, p$k) ||
      any(a$n_patients != 55L) ||
      any(p$n_patients != 55L)) {
    stop("HCC fixed n>=5 figure-data QC failed for ", metric_id)
  }

  list(all = a, pa = p)
}

plot_hcc_depth_panel <- function(metric_id, ylab, letter, threshold = FALSE) {
  z <- get_hcc_depth(metric_id)
  ks <- z$all$k

  ylim <- c(0, 1)
  plot(
    ks,
    z$all$median,
    type = "n",
    xlim = c(1, 4),
    ylim = ylim,
    xaxt = "n",
    xlab = "Number of sampled sectors (k)",
    ylab = ylab,
    las = 1
  )
  axis(1, at = 1:4)

  if (threshold) {
    abline(h = PARAMS$adequate_recovery_threshold, lty = 3, lwd = 1.2)
    text(
      1.05,
      PARAMS$adequate_recovery_threshold + 0.025,
      labels = "Legacy operational 80% threshold",
      adj = 0,
      cex = 0.72
    )
  }

  error_bars(
    ks - 0.035,
    z$all$bootstrap_ci_lower,
    z$all$bootstrap_ci_upper,
    col = cols2[[1L]],
    lwd = 1.5
  )
  lines(
    ks - 0.035,
    z$all$median,
    type = "b",
    pch = 16,
    lwd = 2,
    col = cols2[[1L]]
  )

  error_bars(
    ks + 0.035,
    z$pa$bootstrap_ci_lower,
    z$pa$bootstrap_ci_upper,
    col = cols2[[2L]],
    lwd = 1.5
  )
  lines(
    ks + 0.035,
    z$pa$median,
    type = "b",
    pch = 17,
    lwd = 2,
    col = cols2[[2L]]
  )

  legend(
    "bottomright",
    legend = c(
      "All filtered SNV/indels",
      "Protein-altering sensitivity"
    ),
    col = cols2,
    pch = c(16, 17),
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.82
  )

  panel_letter(letter)
}

figure2_fun <- function() {
  op <- par(
    mfrow = c(1, 2),
    mar = c(4.3, 4.3, 2.0, 1.0),
    oma = c(0.5, 0.5, 2.2, 0.2)
  )
  on.exit(par(op), add = TRUE)

  plot_hcc_depth_panel(
    "nonubiquitous_detection_recall",
    "Non-ubiquitous detection recall",
    "A",
    threshold = TRUE
  )
  plot_hcc_depth_panel(
    "heterogeneity_classification_recall",
    "Heterogeneity-classification recall",
    "B",
    threshold = FALSE
  )

  mtext(
    "HCC fixed n >= 5 cohort (N = 55; k = 1-4 are all true downsampling)",
    side = 3,
    outer = TRUE,
    line = 0.7,
    cex = 1.0,
    font = 2
  )
}

with_two_devices(
  "Figure2_HCC_sampling_depth",
  width = 9.0,
  height = 4.7,
  plot_fun = figure2_fun
)

# ------------------------------------------------------------------------------
# FIGURE 3. Spatial design
# ------------------------------------------------------------------------------
get_spatial_effect <- function(metric_id) {
  z <- sp_eff[
    analysis == "primary_strict_fixed_n_ge5" &
      contrast == "grid_vs_cluster" &
      metric == metric_id
  ][order(k)]

  if (nrow(z) != 3L ||
      !identical(z$k, 2:4) ||
      any(z$n_patients != 44L)) {
    stop("Spatial paired-effect figure-data QC failed for ", metric_id)
  }

  z
}

plot_spatial_effect_panel <- function(metric_id, ylab, letter) {
  z <- get_spatial_effect(metric_id)

  ylim <- range(
    c(z$bootstrap_ci_lower, z$bootstrap_ci_upper, 0),
    finite = TRUE
  )
  pad <- 0.12 * diff(ylim)
  ylim <- c(ylim[[1L]] - pad, ylim[[2L]] + pad)

  plot(
    z$k,
    z$median_paired_advantage,
    type = "n",
    xlim = c(1.8, 4.2),
    ylim = ylim,
    xaxt = "n",
    xlab = "Number of sampled sectors (k)",
    ylab = ylab,
    las = 1
  )
  axis(1, at = 2:4)
  abline(h = 0, lty = 3, lwd = 1.2)

  error_bars(
    z$k,
    z$bootstrap_ci_lower,
    z$bootstrap_ci_upper,
    col = cols2[[1L]],
    lwd = 1.7
  )
  lines(
    z$k,
    z$median_paired_advantage,
    type = "b",
    pch = 16,
    lwd = 2,
    col = cols2[[1L]]
  )

  text(
    z$k,
    z$bootstrap_ci_upper,
    labels = paste0(round(100 * z$proportion_positive), "% positive"),
    pos = 3,
    cex = 0.68
  )

  panel_letter(letter)
}

get_spatial_absolute <- function(metric_id) {
  z <- sp_abs[
    analysis == "primary_strict_fixed_n_ge5" &
      metric == metric_id &
      design %in% c(
        "grid_dispersed",
        "unconstrained",
        "clustered_adjacent"
      )
  ]

  if (nrow(z) != 9L) {
    stop("Spatial absolute-performance figure-data QC failed for ", metric_id)
  }

  z
}

plot_spatial_absolute_panel <- function(metric_id, ylab, letter) {
  z <- get_spatial_absolute(metric_id)
  design_order <- c(
    "grid_dispersed",
    "unconstrained",
    "clustered_adjacent"
  )
  labels <- c(
    "Grid/dispersed",
    "Unconstrained",
    "Clustered/adjacent"
  )

  plot(
    NA,
    xlim = c(1.8, 4.2),
    ylim = c(0, 1),
    xaxt = "n",
    xlab = "Number of sampled sectors (k)",
    ylab = ylab,
    las = 1
  )
  axis(1, at = 2:4)

  for (ii in seq_along(design_order)) {
    d <- z[design == design_order[[ii]]][order(k)]
    lines(
      d$k,
      d$cohort_median,
      type = "b",
      pch = c(16, 17, 15)[[ii]],
      lty = c(1, 2, 3)[[ii]],
      lwd = 2,
      col = cols3[[ii]]
    )
  }

  legend(
    "bottomright",
    legend = labels,
    pch = c(16, 17, 15),
    lty = c(1, 2, 3),
    lwd = 2,
    col = cols3,
    bty = "n",
    cex = 0.78
  )

  panel_letter(letter)
}

figure3_fun <- function() {
  op <- par(
    mfrow = c(2, 2),
    mar = c(4.0, 4.3, 1.9, 1.0),
    oma = c(0.4, 0.5, 2.2, 0.2)
  )
  on.exit(par(op), add = TRUE)

  plot_spatial_effect_panel(
    "recall_nonubiquitous_detection",
    "Paired grid advantage in detection recall",
    "A"
  )
  plot_spatial_effect_panel(
    "recall_heterogeneity_classification",
    "Paired grid advantage in classification recall",
    "B"
  )
  plot_spatial_absolute_panel(
    "recall_nonubiquitous_detection",
    "Absolute detection recall",
    "C"
  )
  plot_spatial_absolute_panel(
    "recall_heterogeneity_classification",
    "Absolute classification recall",
    "D"
  )

  mtext(
    "HCC spatial design: strict fixed n >= 5 cohort (N = 44)",
    side = 3,
    outer = TRUE,
    line = 0.7,
    cex = 1.0,
    font = 2
  )
}

with_two_devices(
  "Figure3_HCC_spatial_design",
  width = 9.0,
  height = 7.5,
  plot_fun = figure3_fun
)

# ------------------------------------------------------------------------------
# FIGURE 4. HCC vs TRACERx100 external validation
# ------------------------------------------------------------------------------
get_cross_curve <- function(metric_id) {
  z <- cross_curve[
    analysis == "fixed_n_ge5" &
      metric == metric_id &
      cohort %in% c("HCC", "TRACERx100")
  ]

  if (nrow(z) != 8L) {
    stop("Cross-cohort curve figure-data QC failed for ", metric_id)
  }

  z
}

plot_cross_curve_panel <- function(metric_id, ylab, letter, threshold = FALSE) {
  z <- get_cross_curve(metric_id)

  plot(
    NA,
    xlim = c(1, 4),
    ylim = c(0, 1),
    xaxt = "n",
    xlab = "Number of sampled regions (k)",
    ylab = ylab,
    las = 1
  )
  axis(1, at = 1:4)

  if (threshold) {
    abline(h = PARAMS$adequate_recovery_threshold, lty = 3, lwd = 1.2)
  }

  cohorts <- c("HCC", "TRACERx100")
  labels <- c("HCC (N=55)", "TRACERx100 (N=15)")

  for (ii in seq_along(cohorts)) {
    d <- z[cohort == cohorts[[ii]]][order(k)]
    offset <- if (ii == 1L) -0.035 else 0.035

    error_bars(
      d$k + offset,
      d$bootstrap_ci_lower,
      d$bootstrap_ci_upper,
      col = cols2[[ii]],
      lwd = 1.5
    )

    lines(
      d$k + offset,
      d$median,
      type = "b",
      pch = c(16, 17)[[ii]],
      lwd = 2,
      col = cols2[[ii]]
    )
  }

  legend(
    "bottomright",
    legend = labels,
    col = cols2,
    pch = c(16, 17),
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.78
  )

  panel_letter(letter)
}

plot_gain_panel <- function() {
  z <- cross_gain[
    analysis == "fixed_n_ge5" &
      metric == "nonubiquitous_detection_recall"
  ]

  x <- 1:3
  labels_x <- c("1->2", "2->3", "3->4")

  plot(
    NA,
    xlim = c(0.7, 3.3),
    ylim = range(
      c(z$bootstrap_ci_lower, z$bootstrap_ci_upper),
      finite = TRUE
    ) * c(0.9, 1.08),
    xaxt = "n",
    xlab = "Increment in sampled regions",
    ylab = "Median within-patient detection gain",
    las = 1
  )
  axis(1, at = x, labels = labels_x)

  cohorts <- c("HCC", "TRACERx100")
  labels <- c("HCC", "TRACERx100")

  for (ii in seq_along(cohorts)) {
    d <- z[cohort == cohorts[[ii]]][order(k_to)]
    offset <- if (ii == 1L) -0.045 else 0.045

    error_bars(
      x + offset,
      d$bootstrap_ci_lower,
      d$bootstrap_ci_upper,
      col = cols2[[ii]],
      lwd = 1.5
    )
    lines(
      x + offset,
      d$median_patient_gain,
      type = "b",
      pch = c(16, 17)[[ii]],
      lwd = 2,
      col = cols2[[ii]]
    )
  }

  legend(
    "topright",
    legend = labels,
    col = cols2,
    pch = c(16, 17),
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.78
  )

  panel_letter("C")
}

plot_reference_architecture_panel <- function() {
  z <- ref_arch[
    variable %in% c("nonubiquitous_fraction", "private_fraction")
  ]

  if (nrow(z) != 2L) {
    stop("Reference-architecture figure-data QC failed.")
  }

  categories <- c("Non-ubiquitous\nfraction", "Private\nfraction")
  x <- 1:2

  plot(
    NA,
    xlim = c(0.6, 2.4),
    ylim = c(0, 1),
    xaxt = "n",
    xlab = "",
    ylab = "Full-reference event fraction",
    las = 1
  )
  axis(1, at = x, labels = categories)

  hcc_med <- z$hcc_median
  hcc_q1 <- z$hcc_q1
  hcc_q3 <- z$hcc_q3
  trx_med <- z$tracerx_median
  trx_q1 <- z$tracerx_q1
  trx_q3 <- z$tracerx_q3

  error_bars(
    x - 0.06,
    hcc_q1,
    hcc_q3,
    col = cols2[[1L]],
    lwd = 2
  )
  points(
    x - 0.06,
    hcc_med,
    pch = 16,
    cex = 1.1,
    col = cols2[[1L]]
  )

  error_bars(
    x + 0.06,
    trx_q1,
    trx_q3,
    col = cols2[[2L]],
    lwd = 2
  )
  points(
    x + 0.06,
    trx_med,
    pch = 17,
    cex = 1.1,
    col = cols2[[2L]]
  )

  legend(
    "topright",
    legend = c("HCC", "TRACERx100"),
    col = cols2,
    pch = c(16, 17),
    bty = "n",
    cex = 0.78
  )

  panel_letter("D")
}

figure4_fun <- function() {
  op <- par(
    mfrow = c(2, 2),
    mar = c(4.0, 4.3, 1.9, 1.0),
    oma = c(0.4, 0.5, 2.2, 0.2)
  )
  on.exit(par(op), add = TRUE)

  plot_cross_curve_panel(
    "nonubiquitous_detection_recall",
    "Non-ubiquitous detection recall",
    "A",
    threshold = TRUE
  )
  plot_cross_curve_panel(
    "heterogeneity_classification_recall",
    "Heterogeneity-classification recall",
    "B",
    threshold = FALSE
  )
  plot_gain_panel()
  plot_reference_architecture_panel()

  mtext(
    "External validation of sampling-depth behavior in TRACERx100",
    side = 3,
    outer = TRUE,
    line = 0.7,
    cex = 1.0,
    font = 2
  )
}

with_two_devices(
  "Figure4_HCC_TRACERx100_external_validation",
  width = 9.0,
  height = 7.5,
  plot_fun = figure4_fun
)

# ------------------------------------------------------------------------------
# SUPPLEMENTARY FIGURE S1. Ceiling diagnostic
# ------------------------------------------------------------------------------
get_available_detection <- function() {
  hcc_all[
    analysis == "available_cohort" &
      metric == "nonubiquitous_detection_recall" &
      k <= 5L
  ][order(k)]
}

get_fixed_detection <- function() {
  hcc_all[
    analysis == "fixed_n_ge5" &
      metric == "nonubiquitous_detection_recall"
  ][order(k)]
}

figureS1_fun <- function() {
  a <- get_available_detection()
  f <- get_fixed_detection()

  op <- par(
    mar = c(4.3, 4.5, 2.2, 1.0)
  )
  on.exit(par(op), add = TRUE)

  plot(
    a$k,
    a$median,
    type = "n",
    xlim = c(1, 5),
    ylim = c(0, 1.03),
    xaxt = "n",
    xlab = "Number of sampled sectors (k)",
    ylab = "Non-ubiquitous detection recall",
    las = 1,
    main = "HCC ceiling diagnostic"
  )
  axis(1, at = 1:5)
  abline(h = PARAMS$adequate_recovery_threshold, lty = 3)

  error_bars(
    a$k,
    a$bootstrap_ci_lower,
    a$bootstrap_ci_upper,
    col = cols2[[1L]],
    lwd = 1.5
  )
  lines(
    a$k,
    a$median,
    type = "b",
    pch = 16,
    lwd = 2,
    col = cols2[[1L]]
  )

  error_bars(
    f$k,
    f$bootstrap_ci_lower,
    f$bootstrap_ci_upper,
    col = cols2[[2L]],
    lwd = 1.5
  )
  lines(
    f$k,
    f$median,
    type = "b",
    pch = 17,
    lwd = 2,
    col = cols2[[2L]]
  )

  legend(
    "bottomright",
    legend = c(
      "Available cohort (changing N; ceiling-sensitive)",
      "Fixed n>=5 cohort (N=55; k=1-4)"
    ),
    col = cols2,
    pch = c(16, 17),
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.8
  )

  text(
    5,
    a[k == 5L, median],
    labels = "87.3% have n = k at k=5",
    pos = 2,
    cex = 0.75
  )
}

with_two_devices(
  "FigureS1_HCC_available_cohort_ceiling_diagnostic",
  width = 6.2,
  height = 5.2,
  plot_fun = figureS1_fun
)

# ----------------------------- figure manifest --------------------------------
manifest <- data.table(
  figure = c(
    "Figure1",
    "Figure2",
    "Figure3",
    "Figure4",
    "FigureS1"
  ),
  stem = c(
    "Figure1_sampling_design_schematic",
    "Figure2_HCC_sampling_depth",
    "Figure3_HCC_spatial_design",
    "Figure4_HCC_TRACERx100_external_validation",
    "FigureS1_HCC_available_cohort_ceiling_diagnostic"
  ),
  format = FIGURE_FORMAT,
  status = "generated_from_frozen_outputs"
)

fwrite(
  manifest,
  file.path(figure_dir, "13_figure_manifest.tsv"),
  sep = "\t"
)


cat("Generated publication figure set in:\n  ", figure_dir, "\n", sep = "")
cat("Figure 1: sampling-design schematic\n")
cat("Figure 2: HCC sampling depth\n")
cat("Figure 3: HCC spatial design\n")
cat("Figure 4: HCC/TRACERx100 external validation\n")
cat("Figure S1: available-cohort ceiling diagnostic\n")
