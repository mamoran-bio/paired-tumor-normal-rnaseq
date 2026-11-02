# Paired tumour/normal RNA-seq analysis of TCGA-BRCA

This project quantifies differential gene expression between primary tumours and matched adjacent-normal tissue in 113 patients from TCGA-BRCA. The analysis uses a paired design (`~ patient + condition`), where each patient acts as their own control. The analysis plan was committed to version control (`556a91d`, 2026-09-14) before accessing any count data. Subsequent amendments and the choice of gene filter are documented in the plan. The repository includes the full pipeline, from estimated count curation to publication figures, containerised with a base image pinned by digest.

**Key numbers:**
- **20,620** of **25,836** genes are significantly differentially expressed at FDR < 0.05.
- **15/15** pre-registered biological controls are significant and have the correct direction.
- **97.7%** of the primary gene list is perfectly stable (significant in all 113 leave-one-out refits).
- **99.63%** of the primary list is preserved when evaluated under an alternative gene filter.

## Key Findings

The project consists of the seven analyses defined in the plan, plus one secondary thresholded test. The core results are:

| Analysis | Result | Evidence |
| :--- | :--- | :--- |
| **1. Primary DE** | 20,620 genes significant at FDR < 0.05. | [`primary-paired.csv`](docs/results/primary-paired.csv) |
| **2. Paired vs. Unpaired** | 20,620 genes paired vs. 20,106 unpaired. Jaccard index 0.825. 1,698 genes called only in the unpaired design. | [`unpaired-comparison.csv`](docs/results/unpaired-comparison.csv) |
| **3. Dispersion** | The paired design estimates a lower median dispersion (0.19) than the unpaired design (0.24). | [`08-unpaired-comparison.log`](docs/logs/08-unpaired-comparison.log) |
| **4. Leave-One-Out** | 20,152 of 20,620 genes (97.7%) remain significant across all 113 patient drop-outs, demonstrating that no single patient determines whether these genes are called. | [`leave-one-out.csv`](docs/results/leave-one-out.csv) |
| **5. Subsampling** | The paired design yields slightly more genes at every n (ratio of medians 1.01–1.08), but the difference is small. At n=113: 20,620 vs 20,106. | [`subsampling.csv`](docs/results/subsampling.csv) |
| **6. FDR Sensitivity** | At FDR < 0.01, 19,114 genes pass (compared to 20,620 at FDR < 0.05). | [`06-primary-analysis.log`](docs/logs/06-primary-analysis.log) |
| **7. Filter Sensitivity** | Expanding the universe (`design=` filter) preserves 20,544 of 20,620 genes (99.63%). Of the 9,748 added genes, 5,435 are significant. | [`filter-sensitivity.csv`](docs/results/filter-sensitivity.csv) |
| **Secondary** | Thresholded test (`lfcThreshold = 1`, H0: absolute log2FC ≤ 1) at `padj < 0.01`. 3,555 genes pass; all 3,555 are also in the primary list. | [`secondary-paired-highconfidence.csv`](docs/results/secondary-paired-highconfidence.csv) |

**A note on the paired design:** While the paired design reduces median dispersion (0.19 vs. 0.24), the *number* of significant genes is similar to the unpaired design. A plausible explanation (not tested here) is that the biological contrast is strong (79.8% of the tested genes are already significant), which likely saturates the primary list.

## Data and Design

- **Source:** TCGA-BRCA via `recount3` (GENCODE v26).
- **Sample curation (Amendment 2026-09-16, `2bc08b7`):**
  The rules were applied in the order: 3 → 4 → 2.
  1. Patient ID and sample type parsed from the barcode.
  2. Only primary tumour (`01`) and solid normal (`11`) were kept.
  3. Resequenced runs of the same aliquot were summed, which is equivalent to deeper sequencing (1,256 runs → 1,232 unique aliquots).
  4. Groups of patient × type with more than one aliquot were resolved by keeping the one with the greatest depth (1,232 → 1,215 aliquots). There were 12 such duplicate groups overall; 4 of them were among the 113 paired patients.
  5. Metastatic samples were excluded by the `01`/`11` rule. No patient had a metastatic sample as their only tumour.
  *Flow (from the "Counts transform closed" amendment of 2026-09-23):* 1,256 runs → 1,232 unique aliquots → 1,215 after depth selection → 1,208 after filtering to `01`/`11` → 113 paired patients (226 samples). `counts_raw` final dimensions: 63,856 genes × 226 samples.
- **Filter:** The primary analysis uses `edgeR::filterByExpr()` with `group = condition`, keeping 25,836 genes of the 63,856. This filter choice is a declared deviation from the original plan and is evaluated in Analysis 7.
- **Model:** DESeq2 1.52.0, paired design (`~ patient + condition`), Wald test. `independentFiltering = FALSE` (per the plan's "Independent filtering" section) and `cooksCutoff = FALSE` (per the 2026-10-17 amendment). `minReplicatesForReplace = Inf` was used for comparability. No further filtering occurs inside DESeq2; all 25,836 genes are retained, with 13 non-converged genes flagged in `beta_converged`.

## Pre-registration

The analysis plan (`docs/analysis-plan.md`) was committed on 2026-09-14 (`556a91d`), written before any count data were downloaded. The count matrix was downloaded and inspected by 2026-09-23, but the script that loads counts (`R/02-load-counts.R`) was not committed until 2026-10-05 (`b7a5dd8`).

Key amendments affecting results include:
- `2bc08b7` (2026-09-16): Sample selection rules (metadata only).
- `7fc8995` (2026-10-07): Use of `group = condition` filter and redefinition of Analysis 7. (Written after a timing fit on 63,856 genes, but no differential expression result was inspected).
- `845f682` (2026-10-17): Disabling outlier replacement and Cook's filtering. (Written after timing runs on 2,500, 5,000, and 25,836 genes and a 200-gene smoke test, but recorded before any differential expression result was inspected).
- `628ea87` (2026-10-21): Handling of non-converged genes.
- `690e5c4` (authored 2026-10-25 20:03 UTC−5, i.e., 2026-10-26 01:03 UTC; the plan heading carries the UTC date): Seed derivation for subsampling.
- `7c51a69` (2026-10-28): Amendment specifying the Analysis 7 script, committed before the script itself (`43ec353`).

Predictions were written before the full run for Analysis 1 (primary), the secondary analysis, and Analyses 2, 3, 4, 5, and 7. No prediction was recorded for Analysis 6. These predictions live outside this repository, except for Analysis 7 in PR #21.

## Results files (`docs/results/`)

| File | Description | Rows (excl. header) | Criteria |
| :--- | :--- | :--- | :--- |
| `primary-paired.csv` | Primary differential expression. | 20,620 | `padj < 0.05` |
| `secondary-paired-highconfidence.csv` | Thresholded test results. | 3,555 | `lfcThreshold = 1`, `padj_thresh < 0.01`. |
| `unpaired-comparison.csv` | Summary of paired vs. unpaired metrics. | 7 | Primary, secondary, and sensitivity results. |
| `leave-one-out.csv` | Stability of the primary list. | 20,620 | Fraction of 113 refits with `padj < 0.05`. |
| `subsampling.csv` | Curve of significant genes by cohort size. | 6 | Medians and 10th–90th percentiles per n. |
| `filter-sensitivity.csv` | Comparison of `group=` and `design=` filters. | 14 | Long format (`metric`, `value`). |

**Notes:**
- "Gene list" in these files refers to genes called significant, not the full filtered table.
- In `unpaired-comparison.csv`, "sensitivity" refers to repeating Analysis 2 after excluding the 13 non-converged paired genes, using the original `padj` without recalculating BH.
- `ratio_med` in `subsampling.csv` is the ratio of medians (descriptive), not a pre-registered test.
- `filter-sensitivity.csv` is in long format and is written by `R/13`, unlike the other CSVs which are written by `R/11`.

## Figures

- **`01-volcano-primary.png`**: Volcano plot of raw `log2FoldChange` vs. `-log10(padj)` from `res_primary`. Positive log2FC indicates higher expression in tumour. The two genes with `padj = 0` are drawn at a declared ceiling of 278.7. The lines at |log2FC| = 1 are a visual reference and do not define significance.
- **`02-ma-primary.png`**: MA plot showing expression level (`baseMean`) vs. raw `log2FoldChange` from `res_primary`. Positive log2FC indicates higher expression in tumour.
- **`03-subsampling-curve.png`**: Median significant genes (with 10th–90th percentiles) for paired and unpaired designs as cohort size increases.
- **`04-leave-one-out-stability.png`**: Distribution of stability fractions across the 20,620 primary genes (y-axis in pseudo-log scale). A gene is "stable" if it remains significant (`padj < 0.05`) in all 113 leave-one-patient-out refits (97.7% of the list).

## Interpretation and limitations

- **Conditional vs. marginal effects:** The paired design estimates the *conditional* effect (within the same patient), while the unpaired design estimates the *marginal* effect (cohort average). The Jaccard index (0.825) indicates a high overlap, but the unique gene sets differ (2,212 unique to paired, 1,698 unique to unpaired). Whether the estimand difference explains this is an interpretation and was not tested here.
- **The 79.8% significant rate:** Of the 25,836 tested genes, 20,620 are significant. The hypothesis that this reflects massive cell-type replacement in breast cancer is consistent with the strong down-regulation of adipocyte markers, though no deconvolution was performed to test it. The 15/15 controls rule out gross errors. Only 3,555 genes (13.8% of tested) pass the stricter thresholded test.
- **Filter sensitivity (Analysis 7):** Expanding the gene filter to `design =` retains 9,748 genes not kept by `group =`. Of these, 5,435 are significant. As exploratory results (PR #21), 107 of these significant genes did not converge; among all 9,748 exclusive genes, 153 did not converge. `nonconverged_B` = 173 is the total for the entire fit B, not just the exclusive genes.
- **Outlier protection:** As recorded in the 2026-10-17 amendment, disabling Cook's flagging and outlier replacement removes those protections from the unpaired design. The paired primary fit never had this protection, because its design matrix never creates large enough cells.
- **Count estimates:** The counts from `recount3` are estimated by dividing the `raw_counts` assay (base-level coverage) by the average mapped read length from the STAR QC metadata (`recount_qc.star.average_mapped_length`), then rounding (`compute_read_counts()`). They are not counts from a dedicated counting tool (e.g., featureCounts).
- **Clinical metadata:** ER was positive in 78 and negative in 21 of the 113 patients (1 indeterminate, 13 missing). HER2 was positive in 24, negative in 68, and missing in 21. PAM50 is empty in all 113 samples. No stratification by molecular subtype is possible, and this unmodeled heterogeneity is a limitation.

## Repository layout

```text
.
├── R/
│   ├── 00-helpers.R            # Shared function (check_controls)
│   ├── 01-select-tissue.R      # Curates BRCA metadata without downloading counts
│   ├── 01-count-pairs.R        # Counts paired patients across 33 projects (writes docs/tissue-pair-counts.csv)
│   ├── 02-load-counts.R        # Download and curate counts
│   ├── 03-build-coldata.R      # Prepare sample metadata
│   ├── 04-filter-genes.R       # Apply edgeR::filterByExpr()
│   ├── 05-fit-models.R         # Fit DESeq2 models (writes model-inputs, dds_paired, dds_unpaired)
│   ├── 06-primary-analysis.R   # Extract primary results (writes res_primary)
│   ├── 07-secondary-analysis.R # Extract thresholded results (writes res_secondary)
│   ├── 08-unpaired-comparison.R# Extract unpaired and compare (writes res_unpaired, analysis_2, analysis_3)
│   ├── 09-subsampling.R        # Run one subsampling replicate
│   ├── 09b-aggregate.R         # Aggregate subsampling (writes analysis_5)
│   ├── 10-leave-one-out.R      # Run one LOO replicate
│   ├── 10b-aggregate.R         # Aggregate LOO (writes analysis_4)
│   ├── 11-save-results.R       # Export 5 CSVs and sessionInfo-analysis.txt
│   ├── 12-plots.R              # Generate figures (writes docs/figures/*.png)
│   ├── 13-filter-sensitivity.R # Analysis 7 (writes filter-sensitivity.csv and analysis_7.rds)
│   ├── timing-serial.R         # Benchmarking
│   └── timing-parallel.R       # Benchmarking
├── scripts/
│   ├── run-batch.sh            # Parallelizes R/09 (3 concurrent jobs, --memory=4g)
│   └── run-loo.sh              # Parallelizes R/10 (3 concurrent jobs, --memory=4g)
├── docs/
│   ├── analysis-plan.md        # Pre-registration and amendments
│   ├── tissue-pair-counts.csv  # Output of 01-count-pairs.R
│   ├── figures/                # 4 PNGs
│   ├── results/                # 6 versioned CSVs
│   ├── logs/                   # 24 versioned logs (10 top-level for R/05-R/13; 14 launcher and pilot logs in subsampling/ and loo/)
│   ├── timing/                 # Benchmarking logs (aarch64/ and x86_64/)
│   ├── sessionInfo.txt         # Frozen environment (x86_64 timing run)
│   ├── sessionInfo-arm64.txt   # Frozen environment (arm64 image validation)
│   └── sessionInfo-analysis.txt# Written by R/11 (x86_64)
├── results/                    # Untracked intermediate .rds files (ignored by git)
├── Dockerfile                  # Pins base image by digest; checks package presence
├── LICENSE
├── .gitignore
└── .dockerignore
```

## How to reproduce

1. **Build the Docker image:**
   ```bash
   docker build -t paired-rnaseq:dev .
   ```
   *Note:* The `Dockerfile` pins the Bioconductor 3.23 base image by digest. The package versions depend on `BiocManager::install()` at build time and are not fixed between builds (no `renv.lock` is used). Rebuilding the image later may install newer versions; compare `sessionInfo()` with `docs/sessionInfo-analysis.txt` to verify the environment.

2. **Run the pipeline:**
   Run all commands from the repository root. The scripts are not strictly sequential; `R/02` to `R/04` are sourced by `R/05` (and `R/13`), so they do not need to be run separately. `R/01-count-pairs.R` generates `docs/tissue-pair-counts.csv` and is not required to regenerate the differential expression results. The bash scripts skip fits whose `.rds` already exists (`[SKIP]`), so relaunching them resumes where they left off.

   ```bash
   # 1. Fit models (reuses results/dds_paired.rds if it exists. Delete it to force a re-fit.)
   docker run --rm --user "$(id -u):$(id -g)" -e OPENBLAS_NUM_THREADS=1 \
     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
     Rscript R/05-fit-models.R 2>&1 | tee docs/logs/05-fit-models.log

   # 2. Extract results (R/08 depends on R/06 and R/07)
   docker run --rm --user "$(id -u):$(id -g)" -e OPENBLAS_NUM_THREADS=1 \
     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
     Rscript R/06-primary-analysis.R 2>&1 | tee docs/logs/06-primary-analysis.log
   # (Repeat for R/07 and R/08)

   # 3. Subsampling
   ./scripts/run-batch.sh
   docker run --rm --user "$(id -u):$(id -g)" -e OPENBLAS_NUM_THREADS=1 \
     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
     Rscript R/09b-aggregate.R 2>&1 | tee docs/logs/09b-aggregate.log

   # 4. Leave-one-out (R/10 needs R/06; R/13 needs R/10b)
   ./scripts/run-loo.sh
   docker run --rm --user "$(id -u):$(id -g)" -e OPENBLAS_NUM_THREADS=1 \
     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
     Rscript R/10b-aggregate.R 2>&1 | tee docs/logs/10b-aggregate.log

   # 5. Save results, generate plots, and run Analysis 7
   docker run --rm --user "$(id -u):$(id -g)" -e OPENBLAS_NUM_THREADS=1 \
     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
     Rscript R/11-save-results.R 2>&1 | tee docs/logs/11-save-results.log
   docker run --rm --user "$(id -u):$(id -g)" -e OPENBLAS_NUM_THREADS=1 \
     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
     Rscript R/12-plots.R 2>&1 | tee docs/logs/12-plots.log
   docker run --rm --user "$(id -u):$(id -g)" -e OPENBLAS_NUM_THREADS=1 \
     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
     Rscript R/13-filter-sensitivity.R 2>&1 | tee docs/logs/13-filter-sensitivity.log
   ```

   *Note on file permissions:* The `.sh` scripts do not include `--user`. In this environment, outputs were owned by the user. On a native Linux Docker installation, files may be owned by `root`. If you encounter permission issues, run `sudo chown -R "$(id -u):$(id -g)" results/ docs/logs/` after the runs. Note that `tee` will overwrite the versioned logs.

3. **Intermediate files:**
   The `results/` directory contains untracked `.rds` files. To regenerate them, run the scripts in the order above.

4. **System requirements and timings:**
   - Platform: The full pipeline and all reported results were generated on the x86_64 laptop (i7-8550U, 15.5 GiB).
   - Memory: 12 GB for 3 concurrent jobs (each limited to `--memory=4g`). The peak memory in `R/02` is ~3.4 GiB.
   - Time: `R/05` took 2,534.4 seconds (~42 minutes) for the paired fit. The leave-one-out analysis took ~44 hours of elapsed wall-clock time (including gaps between launches), from 2026-10-26 10:17 to 2026-10-28 06:18 local time (UTC−5). `R/13` took 80.7 minutes.
   - *Note on architecture:* Only `R/02` (as image validation) and the timing scripts (as benchmarks) were run on the Mac (aarch64). No full-universe fit was run there, serial or parallel. Floating-point results are not expected to be bit-identical between architectures.

## Notes and credits

- **Timestamps:** The timestamps printed by R inside the container are in UTC. The git commit timestamps in this repository are in the author's local time (UTC−5).
- **Gene identifiers:** Overlap metrics (in `R/08` and `R/13`) and LOO stability (in `R/10b`) are computed on `gene_id` (Ensembl gene IDs with version suffix, GENCODE v26). `gene_label` disambiguates duplicate symbols and is provided only for readability.
- **License:** This project is licensed under the MIT License. See the `LICENSE` file for details.

### Data source and citations

The results shown here are in whole or part based upon data generated by the TCGA Research Network: https://www.cancer.gov/tcga.

If you use this pipeline or its results, please cite the underlying tools:

- **recount3:** Wilks C, Zheng SC, Chen FY, Charles R, Solomon B, Ling JP, Imada EL, Zhang D, Joseph L, Leek JT, Jaffe AE, Nellore A, Collado-Torres L, Hansen KD, Langmead B (2021). “recount3: summaries and queries for large-scale RNA-seq expression and splicing.” *Genome Biology*. doi:10.1186/s13059-021-02533-6.
- **DESeq2:** Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” *Genome Biology*, 15(12):550. doi:10.1186/s13059-014-0550-8.
- **edgeR:** Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025). “edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” *Nucleic Acids Research*, 53(2), gkaf018. doi:10.1093/nar/gkaf018.
- **apeglm:** Zhu A, Ibrahim JG, Love MI (2018). “Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.” *Bioinformatics*. doi:10.1093/bioinformatics/bty895.
