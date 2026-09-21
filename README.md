# Multiregion sampling analysis

Reproduces manuscript **Figures 1-4 and Supplementary Figure 1** (the HCC
all-filtered and protein-altering SNV/indel sampling-depth analysis,
spatial designs, and TRACERx100 external validation), plus **Figure 5**
(driver-mutation, CNV, and RNA-pathway sampling extensions) when its own
input files are supplied. Figure 5's driver-mutation panel is itself
optional within Figure 5: if `snv_indel.tsv` isn't available, Figure 5 is
still generated with the CNV + RNA panels only (see §7-§8).

## Step-by-step Instructions

**1. Download this repository.** Either:
- clone it with `git`:
  ```bash
  git clone <REPOSITORY_URL>
  ```
- or, on the repository's page in your browser, click the green **`<> Code`** button → **Download ZIP**, then extract it.

**2. Open R.** After downloading, you have one folder (from `git clone` it's named `sampling_experiments`; from a ZIP download it may be named `sampling_experiments-main` — either way, it's the same folder, containing `README.md` and `scripts/`). Then, either:
- **RStudio:** open the app, then
  - use *File → Open Project...* and pick `scripts/sampling_experiments.Rproj` — this opens RStudio with the working directory already set correctly, so you can skip the `setwd()` line below; or
  - open RStudio without a project and set the working directory yourself:
    ```r
    setwd("/full/path/to/sampling_experiments/scripts")
    ```
- **Terminal:** type `R` and press Enter to start an interactive R session, then set the working directory to the `scripts/` folder:
  ```r
  setwd("/full/path/to/sampling_experiments/scripts")
  ```

**3. Install the R packages this pipeline needs.** Still inside R:
```r
install.packages(c("data.table", "readxl"))
```
`readxl` is only used by the Figure 5 CNV/RNA scripts (19-21); the core pipeline (01-13) never loads it.

**4. Get the input data.** You need one data folder containing all of these together:
- `wgs_samples_matching.RData`
- `somatic_mutations.RData`
- `clinical_data.RData`
- `patient_level_summary.tsv` — the original/legacy pipeline's own output, required by script 04 as an independent regression check; without it, that script (and everything downstream of it) will stop
- `nsclc_tracerx_2017.tar.gz` — only needed for the TRACERx100 arm and the manuscript figures (scripts 09-13); omit it to run just scripts 01-08

Optional, for Figure 5 — each block below runs only if its own files are present:
- `mmc4.xlsx` — needed for scripts 19, 20, 21 (focal CNV, broad CNV, RNA pathway) and for script 24 to generate Figure 5 at all
- `Driver_48genes_DICER1.csv`, `snv_indel.tsv` (or `.zip`), `cnv_arm_level.tsv` — needed **together** for scripts 14-15 (driver mutations); without all three, Figure 5 is still generated (from `mmc4.xlsx` alone) but without the driver panel or the protein-altering comparator in the genomic panels

**5. Configure the pipeline.** Open `00_config.R` (in `scripts/`) and set:

| Parameter | Set it to | Notes |
|---|---|---|
| `INPUT_PATH` | the full path to your data folder from step 4 | **required** |
| `OUTPUT_PATH` | a full path, or leave as `NULL` | optional; `NULL` writes outputs to `<data folder>/output` |
| `DEBUG` | `FALSE` (default) or `TRUE` | optional; `FALSE` writes only what §7/§8 below list (enough to reproduce every figure); `TRUE` additionally writes every internal diagnostic/sensitivity table, for debugging — see §7 |
| `FIGURE_FORMAT` | `"pdf"`, `"png"`, or `"both"` (default) | which file format(s) Scripts 13 and 24 save each figure in |

No other file needs editing to run the pipeline.

**6. Run a script.** In R or RStudio, from the `scripts/` working directory:
```r
# optional but recommended: verify your input files first
source("verify_inputs.R")
verify_inputs(
  processed_data_dir  = "/full/path/to/your_data_folder",
  tracerx_archive_path = "/full/path/to/nsclc_tracerx_2017.tar.gz"   # omit to skip this check
)

# run any single script on its own
source("01_hcc_metadata_qc.R")
```

**7. Run the whole pipeline.** In R, RStudio, or `Rscript`:
```r
source("run_pipeline.R")
```
This sources scripts `01` through `13` in order, unattended, stopping only if a script raises an error. If `nsclc_tracerx_2017.tar.gz` isn't in your data folder, it skips scripts 09-13 automatically. After it finishes, check `qc/01_qc_summary.txt` and confirm both legacy QC lines read `PASS` before trusting downstream results — script 01 reports a mismatch as a warning rather than aborting, so the pipeline does not stop on it.

It then runs the Figure 5 extensions, if their inputs (§4) are present:
- scripts `14`-`15` (driver mutations) and `19`-`21` (CNV/RNA) are fully independent of each other and each print what's missing and are skipped otherwise;
- script `24` (Figure 5 itself) runs whenever `mmc4.xlsx` is available, regardless of whether the driver arm ran — with only CNV/RNA data, Figure 5 is generated with those panels only, and prints a note saying so.

Every check still runs, in full, every time -- nothing about `DEBUG` skips or weakens a single `stop()`/`warning()`. That switch only controls what gets **written to disk**:
- `DEBUG = FALSE` (default): `results/` contains only the tables that are literally read as plotting input by Script 13 (Figures 1-4/S1) or Script 24 (Figure 5) — verified by reading their own `fread()` calls, not guessed from folder names — the exhaustive-subset/patient-summary substrate tables for HCC and TRACERx100 (Scripts 04 and 11), the handful of tables other scripts genuinely re-read as real inputs (not just diagnostics), and the short `*_qc_summary.txt` narratives (enough to see every PASS/FAIL, and enough for `generate_validation_summary.R`). This is a lean run: roughly 25 files without Figure 5, ~35 with it.
- `DEBUG = TRUE`: additionally writes every per-check diagnostic detail table and every sensitivity/robustness table that supports a sentence in the text but was never plotted as its own figure panel (e.g. the legacy 80%-threshold or mean-vs-median tables, or the many driver/CNV/RNA crosswalk and provenance tables). Roughly 200-260 files total. Nothing computed differs between the two settings — only what lands on disk.

Finally, it writes one consolidated `results/VALIDATION_SUMMARY.txt`, collecting every PASS/FAIL line from every script's own summary into a single audit log. You can also regenerate it on its own, any time:
```r
source("generate_validation_summary.R")
generate_validation_summary("/full/path/to/your_data_folder/output")
```

**8. Look at the results.** Everything lands under `<your_data_folder>/output/`:

```
output/
├── results/
│   ├── qc/               ← scripts 01-12: short *_qc_summary.txt narratives always;
│   │                        detailed diagnostic .tsv tables only if DEBUG = TRUE
│   ├── qc/extensions/     ← scripts 14, 15, 19, 20, 21: same rule
│   ├── hcc_depth/         ← scripts 04, 05 (Fig. 2, Suppl. Fig. 1)
│   ├── hcc_spatial/       ← scripts 06, 07 (Fig. 3)
│   ├── hcc_sensitivity/   ← script 08
│   ├── tracerx100/        ← script 11
│   ├── cross_cohort/      ← script 12 (Fig. 4)
│   ├── extensions/
│   │   ├── driver/          ← scripts 14, 15
│   │   ├── cnv/               ← scripts 19, 20
│   │   ├── rna/                 ← script 21
│   │   └── integration/          ← script 24 (Fig. 5)
│   └── VALIDATION_SUMMARY.txt
├── intermediate/          ← .rds hand-off objects between scripts (not for reading)
├── figures/               ← scripts 13, 24: Figures 1-5 + Suppl. Fig. 1, in FIGURE_FORMAT
└── external/              ← script 09's extracted TRACERx100 archive (audit trail)
```

None of these are checked into the repository — they're regenerated every time you run the pipeline. There is no separate "paper package" folder: with the default `DEBUG = FALSE`, `results/` already *is* the paper package.

**With `DEBUG = FALSE` (default), `results/` contains exactly:**
- `qc/01_qc_summary.txt` … `qc/12_hcc_tracerx100_qc_summary.txt`, and (if their inputs are present) `qc/extensions/14_driver_qc_summary.txt`, `15_driver_sampling_summary.txt`, `19_focal_cnv_sampling_summary.txt`, `20_broad_cnv_sampling_mmc4_summary.txt`, `21_rna_pathway_sampling_summary.txt` — one short narrative per script, always kept
- `qc/01_tumor_sample_map_audit.tsv` and `qc/02_reference_summary_by_patient.tsv` — **not** diagnostics-only: re-read by Scripts 15 and 12 respectively as real inputs, so always written regardless of `DEBUG`
- `hcc_depth/04_hcc_exhaustive_subsets_all_filtered.tsv`, `04_hcc_unconstrained_patient_summary_all_filtered.tsv`
- `hcc_depth/05_hcc_depth_core_summary_all_filtered.tsv`
- `hcc_spatial/06_hcc_spatial_selected_subsets_primary_strict.tsv`, `06_hcc_spatial_selected_subsets_order_sensitivity.tsv` — re-read by Script 08
- `hcc_spatial/07_hcc_spatial_main_effect_summary.tsv`, `07_hcc_spatial_absolute_design_summary_all_filtered.tsv`
- `hcc_sensitivity/08_hcc_protein_altering_depth_summary.tsv`
- `cross_cohort/12_hcc_tracerx100_fixed_n_ge5_curve_summary.tsv`, `12_hcc_tracerx100_fixed_n_ge5_paired_gains.tsv`, `12_hcc_tracerx100_reference_architecture_comparison.tsv`
- `tracerx100/11_tracerx100_exhaustive_subsets_protein_altering.tsv`, `11_tracerx100_unconstrained_patient_summary_protein_altering.tsv`
- `VALIDATION_SUMMARY.txt`
- `figures/` — Figure1, Figure2, Figure3, Figure4, FigureS1 (in `FIGURE_FORMAT`), plus `13_figure_manifest.tsv`; Figure5 too if its inputs are present

**Figure 5 only, additionally kept regardless of `DEBUG`** (if `mmc4.xlsx` is present):
- `extensions/driver/14_driver_exact_event_candidates_preeligibility.tsv` — re-read by Script 15 (only if the driver arm ran)
- `extensions/driver/15_driver_depth_summary_fixed_n_ge5.tsv`, `15_driver_paired_comparisons_fixed_n_ge5.tsv` (only if the driver arm ran)
- `extensions/cnv/19_focal_cnv_depth_summary_fixed_n_ge5.tsv`, `19_focal_cnv_vs_protein_paired_fixed_n_ge5.tsv`
- `extensions/cnv/20_broad_cnv_depth_summary_fixed_n_ge5.tsv`, `20_broad_cnv_vs_protein_paired_fixed_n_ge5.tsv`
- `extensions/rna/21_rna_depth_summary_fixed_n_ge5.tsv`
- `extensions/integration/24_genomic_detection_curves.tsv`, `24_genomic_classification_curves.tsv`, `24_genomic_k4_paired_vs_protein.tsv`, `24_rna_pathway_curves.tsv`, `24_multimodal_key_findings.tsv`, `24_multimodal_synthesis_provenance.tsv`, `24_multimodal_extension_summary.txt` — Script 24's own final tables, always kept as the Figure 5 deliverable
- `figures/Figure5_multimodal_sampling_extensions.*` (in `FIGURE_FORMAT`)

Set `DEBUG <- TRUE` in `00_config.R` to additionally get every diagnostic/sensitivity table behind these results (see `docs/PIPELINE_MAP.md` for the full script-by-script breakdown).
