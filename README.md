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

**4. Get the input data.** You need one data folder containing all of these together:
- `wgs_samples_matching.RData`
- `somatic_mutations.RData`
- `clinical_data.RData`
- `patient_level_summary.tsv`
- `nsclc_tracerx_2017.tar.gz`
- `mmc4.xlsx`
- `Driver_48genes_DICER1.csv`, `snv_indel.tsv` (or `.zip`), `cnv_arm_level.tsv`

**5. Configure the pipeline.** Open `00_config.R` (in `scripts/`) and set:

| Parameter | Set it to | Notes |
|---|---|---|
| `INPUT_PATH` | the full path to your data folder from step 4 | **required** |
| `OUTPUT_PATH` | a full path, or leave as `NULL` | optional; `NULL` writes outputs to `<data folder>/output` |
| `DEBUG` | `FALSE` (default) or `TRUE` | optional; if `TRUE` the pipeline writes every internal diagnostic/sensitivity table, for debugging |
| `FIGURE_FORMAT` | `"pdf"`, `"png"`, or `"both"` (default) | the file format(s) of the output figures |


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
