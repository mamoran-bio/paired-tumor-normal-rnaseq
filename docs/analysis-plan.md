# Analysis plan

This plan was written before any count data were downloaded. The commit date in
the git history is the evidence for that ordering. Amendments are appended at
the end of this document with their date and reason; nothing above is edited
silently.

## Question

Which genes are differentially expressed between primary tumour and matched
adjacent-normal tissue in a single TCGA tissue, and how much does the resulting
gene list depend on modelling the patient pairing?

The gene list is the primary biological output. The paired-versus-unpaired
comparison is reported as a methodological result in its own right.

## Data

Source: TCGA via recount3 — uniformly processed gene-level counts, GENCODE v26.

Inclusion: only patients contributing BOTH a primary tumour sample and a
matched adjacent-normal sample. Patients contributing only one of the two are
excluded, even though that discards most of the cohort.

Adjacent-normal tissue is not healthy tissue. It is sampled near a resected
tumour and may carry early field-effect changes. It is called "adjacent normal"
throughout this repository and never "healthy" or "control".

## Tissue selection rule

The tissue is selected on design grounds only: the TCGA project with the
largest number of patients contributing a matched tumour / adjacent-normal
pair, among projects with at least 20 such pairs. Twenty is set here as a floor
below which the leave-one-patient-out and subsampling analyses below stop being
informative.

Selection uses sample metadata only. No count matrix is downloaded, and no
differential expression is run, before the tissue is fixed and recorded in the
amendment log.

## Design

    ~ patient + condition

Each patient acts as their own control. Tumour and adjacent-normal samples from
one patient share germline background, age, comorbidity and, usually,
processing batch; treating them as independent pushes all of that
between-patient variation into the residual and masks real signal.

This model cannot estimate patient-level covariates — age, sex, stage — because
they are constant within a patient and are absorbed by the patient term. That
is a deliberate cost of the design, not an oversight: the design buys power for
the within-patient contrast and gives up the ability to ask about
between-patient variables.

## Pre-specified filtering

`edgeR::filterByExpr()` with default parameters, using the design matrix above.

Filtering depends only on overall expression level. It never uses the
tumour-normal difference, because filtering on the quantity being tested is
circular and inflates the false positive rate.

## Pre-specified thresholds — primary analysis

FDR: 0.05, Benjamini-Hochberg.

Effect size: no minimum required for the primary list. Log2 fold changes are
reported for every gene and used to rank and interpret results, but are not
used as a selection filter. A post-hoc |log2FC| cutoff applied to a
threshold-zero test invalidates FDR control, because the reported adjusted
p-values test the null "the change is zero", not "the change exceeds the
cutoff".

## Pre-specified secondary analysis

A high-confidence subset is defined in advance as the genes significant at
FDR 0.01 with an effect size tested directly against a threshold, using
`DESeq2::results(dds, lfcThreshold = 1)`.

This subset is reported alongside the primary list, never instead of it. Both
thresholds are fixed before any data are downloaded.

## Planned analyses

1. Primary differential expression under `~ patient + condition`.
2. Paired versus unpaired: the same data fitted under `~ condition`. Report the
   number of genes called at the same FDR under each design, the overlap, and
   the direction of the difference.
3. Dispersion estimates under both designs, to check whether the patient term
   reduces within-condition dispersion as the design predicts.
4. Leave-one-patient-out stability: refit the paired model dropping each
   patient in turn, and report for each called gene the fraction of refits in
   which it remains significant.
5. Subsampling: refit both designs on random subsets of n = 5, 10, 20, ...
   patients, repeated, and report how many genes each design calls as a
   function of n. This locates where the paired design pays off most.
6. Sensitivity to the FDR threshold: how the primary list changes at 0.01
   versus 0.05.
7. Sensitivity to the filtering rule: `filterByExpr()` defaults versus a manual
   rule of at least 10 counts in at least N samples, N = the smaller group size.

## Reporting rules

All analyses listed above are reported whatever their outcome, including
outcomes that weaken the case for the paired design. If analysis 2 shows the
paired design calls fewer genes, or analysis 3 shows no reduction in
dispersion, that is reported as the result and discussed, not dropped.

Gene lists are published as CSV files under `docs/results/`, with the code that
produced them.

## What this analysis cannot establish

### Structural limitations — from the design and the data type

**Cell composition.** Bulk RNA-seq measures a mixture of cell types. A gene can
appear over-expressed in tumour simply because the tumour contains more of the
cell type that normally expresses it — immune infiltrate, stroma, altered
epithelial fraction — with no change in regulation inside any cell. Bulk data
alone cannot separate differential composition from differential regulation;
that needs deconvolution or single-cell data, neither of which is done here.

**Direction of effect.** The design is observational and cross-sectional. An
expression difference between tumour and adjacent normal does not say whether
the change contributes to tumour development, results from it, or both.

**The comparison actually made.** The contrast is tumour versus adjacent
normal, not tumour versus healthy tissue. If field effects extend into the
adjacent-normal sample, genes altered early and uniformly across the organ will
be under-detected, because they are altered on both sides of the comparison.

**Who the result applies to.** Patients with a matched adjacent-normal sample
are those who underwent surgical resection, which skews towards operable
disease and away from advanced or non-surgical presentations. The cohort is not
a random sample of patients with this cancer, and the result does not
automatically transfer to those excluded.

**Batch structure.** TCGA has documented plate and tissue-source-site effects.
The paired design absorbs batch only to the extent that a patient's two samples
were processed together, which is usual but not guaranteed and is not verified
here.

### Tissue-specific limitations

To be completed once the tissue is selected under the rule above, and before
any count data are downloaded. Recorded in the amendment log.

## Amendment log

Amendments are appended here with date and reason. Sections above are not
edited in place.


### 2026-09-16 — Sample selection rules, after inspecting TCGA-BRCA metadata

Metadata for TCGA-BRCA was inspected under the tissue selection rule above:
metadata only, no count data downloaded. These five rules are fixed before
any count matrix is retrieved.

1. There were 10 missing values in submitter_id and sample_type_id within the BRCA dataset, so it was necessary to use the TCGA barcode to identify the patient and sample type, as it contains both pieces of data. Otherwise, patient TCGA-A7-A0DC would have been excluded from the analysis.

2. Looking into the BRCA sample types, I found that "Primary Tumor" and "Primary solid Tumor" are the same because positions 14 and 15 in the barcode detail the sample type (both were "01"). Only samples with "01" (primary tumor) and "11" (solid tissue normal) are included; all other sample types are excluded.

3. The BRCA dataset contains 1,256 rows, but only 1,232 unique aliquots. There were 17 duplicated barcodes: aliquots that were sequenced more than once — 10 twice and 7 three times. These duplicate samples (those that were resequenced) were summed, because they are independent reads from the same library; summing is equivalent to deeper sequencing and preserves the count nature of the data, avoiding giving them double mathematical weight.

4. After summing the resequenced runs, 4 of the 113 paired patients had more than one aliquot of the same sample type; for each of these patients, the aliquot with the greatest sequencing depth was selected for the analysis.

5. Patients with metastatic samples were excluded because this paired analysis is designed for primary tumor vs. adjacent normal tissue, and metastasis represents a distinct biological entity. No patient had a metastatic sample as their only tumor sample, so excluding metastatic samples removes no patient from the analysis.


### 2026-09-18 — Tissue selected

Tissue selected: TCGA-BRCA, 113 paired patients — the largest of the 33 TCGA projects, of which 13 meet the pre-specified floor of 20 pairs. Second is KIRC with 72. Selection ran on metadata only; see docs/tissue-pair-counts.csv and R/01-count-pairs.R.


#### Selection script

- count_pairs() + check_results() functions with five checks.

- Loop over the 33 projects, CSV written sorted by n_paired, descending.

- Five checks output: all ok.

- Verified against the 113 paired patients.

- n_tumor rather than n_tumour


### 2026-09-18 — Tissue-specific limitations

Recorded before any count data were downloaded, using the recount3 metadata
endpoint only.

**Receptor status.** Among the 113 paired patients, ER was Positive in 78 and
Negative in 21 (1 indeterminate, 13 missing); PR was Positive in 69 and
Negative in 30 (1 indeterminate, 13 missing). HER2 was assessed by IHC
(available for 87 patients) and by in situ hybridization (available for 22).
IHC was taken as the primary call where available and non-equivocal, and ISH
was used to resolve equivocal and missing IHC cases following standard clinical
practice; the 2 cases with IHC-positive / ISH-negative results were retained
as positive under this rule. The combined HER2 status is 24 positive, 68
negative, and 21 with no HER2 data. Receptor-defined subtypes (luminal,
HER2-enriched, triple-negative) have distinct expression profiles and are not
controlled in the paired design. PAM50 molecular subtypes are absent from the
TCGA clinical model served by recount3 — the corresponding field is entirely
empty — and are published separately as part of the PanCanAtlas.

**Neoadjuvant treatment and cohort selection.** Of the 113 paired patients,
111 are recorded as not having received neoadjuvant treatment, 1 as having
received it, and 1 has no record. This field cannot distinguish a genuinely
treatment-naive cohort from under-reporting. Additionally, the paired subset is
not a random sample of TCGA-BRCA: patients with banked adjacent-normal tissue
are by construction those who underwent primary surgery, which is the group
least likely to have received neoadjuvant therapy. Downstream interpretation
should assume neither reading.

**Normal tissue composition.** Breast normal tissue composition (epithelial /
stromal / adipose) is not recorded directly, and the patients themselves vary:
55 of 113 are post-menopausal, 28 pre-menopausal, 2 peri-menopausal, 1
indeterminate, and 27 have no menopausal status recorded. Composition changes
with menopausal status, age, and adiposity, all of which differ across the
cohort. The "normal" reference is therefore not a uniform baseline across
patients.

### 2026-09-18 — Implementation decisions before count data

These decisions are recorded before any count matrix is downloaded or any
differential expression is run. Each corresponds to a script under `R/`, in the
order listed below. The scripts are sourced in sequence from a single session,
so objects persist across them.

#### Script layout

- `R/00-setup.R` — load libraries, fix the random seed.
- `R/01-count-pairs.R` — count paired tumor/normal patients across the 33
  TCGA projects (already in the repository).
- `R/02-load-counts.R` — download the BRCA gene-level count matrix from
  recount3 and align it with the 113 paired patients. Produces `counts_raw`
  and `sample_annotation`, consumed by the scripts below.
- `R/03-build-coldata.R` — build the sample-level table and fix factor levels.
- `R/04-filter-genes.R` — apply `edgeR::filterByExpr()` as a pre-filter.
- `R/05-fit-models.R` — fit the paired and unpaired DESeq2 models.
- `R/06-primary-analysis.R` — extract the primary results.
- `R/07-secondary-analysis.R` — thresholded test and shrunken log2 fold
  changes.
- `R/08-unpaired-comparison.R` — unpaired results and the three comparison
  quantities.
- `R/09-subsampling.R` — subsampling across both designs.
- `R/10-leave-one-out.R` — leave-one-patient-out stability.
- `R/11-save-results.R` — write CSVs and capture the R environment.

#### Factor levels and coefficient name

The factor levels for `condition` are set explicitly with `normal` as the
reference, so that the log2 fold change is signed as tumor relative to normal
and the coefficient name is `condition_tumor_vs_normal`. This is done in
`R/03-build-coldata.R`.

The coefficient name is captured once in the same script and reused by every
downstream `results()` call, guarded after model fitting by a single
`stopifnot(coef_name %in% resultsNames(dds))` in `R/05-fit-models.R`. Without
that guard, a name mismatch would silently select no coefficient or the wrong
one. American spelling (`tumor`) is used throughout, consistent with
`n_tumor` in the count-pairs output.

#### Pipeline: DESeq2 as the single estimation engine

The sections above reference `edgeR::filterByExpr()` and
`DESeq2::results()` without declaring which package owns the analysis. This is
now fixed:

- Normalisation and dispersion estimation: DESeq2, via
  `DESeqDataSetFromMatrix()` and `DESeq()`. See `R/05-fit-models.R`.
- Size factors: `DESeq2::estimateSizeFactors()`, median-of-ratios, default
  parameters.
- Dispersion: DESeq2 empirical Bayes, default `fitType = "parametric"`.
- `filterByExpr()` is used only as a pre-filter on the count matrix before
  building the `DESeqDataSet` (`R/04-filter-genes.R`). It is not part of the
  estimation.
- No edgeR model is fitted. edgeR appears only as the source of
  `filterByExpr()`.

Rationale: mixing dispersion models across packages would make the
paired-versus-unpaired comparison (analysis 2) and the dispersion comparison
(analysis 3) depend on two different estimation methods. One engine, one set
of assumptions.

#### Pre-filtering: what `filterByExpr()` receives

`edgeR::filterByExpr()` does not accept a formula. It takes either a `group`
vector or a `design` matrix. The paired design matrix is built once in
`R/04-filter-genes.R` and passed in. The resulting `counts_filtered` matrix is
the input to both the paired and the unpaired `DESeqDataSet`, so the two
designs share an identical gene universe.

#### Independent filtering

DESeq2's independent filtering is disabled in every `results()` call
(`independentFiltering = FALSE`), so that the number of genes entering the FDR
correction is controlled by the pre-filter only and is identical across the
paired and unpaired fits.

Disabling independent filtering reduces power relative to the DESeq2 default,
which would otherwise choose the filter that maximises the number of genes at
the target FDR. That power is traded away in exchange for comparability between
designs: analysis 2 is itself a result of the project, and it must not be
confounded by two different filter statistics. The decision is deliberate and
the cost is accepted.

Both filters depend only on overall expression level. Neither uses the
tumor-normal difference.

#### Primary analysis

The paired model is fit in `R/05-fit-models.R`, and the primary results are
extracted in `R/06-primary-analysis.R`. No shrinkage is applied to the primary
list: the quantity being ranked is the evidence against a zero effect, not the
magnitude of a non-zero one. Log2 fold changes from this table are used for
ranking and interpretation.

#### Secondary analysis

The secondary high-confidence subset is defined in `R/07-secondary-analysis.R`
as genes with `padj < 0.01` from a thresholded test at `lfcThreshold = 1` on
the same fitted model. The null hypothesis of this test is |log2FC| <= 1, and
the `padj` column it returns controls FDR for that null. It is the only source
of significance calls in the secondary analysis.

Shrunken log2 fold changes for ranking and visualisation are computed
separately in the same script with `lfcShrink(type = "apeglm")`. The shrinkage
reuses the already-fitted model; it does not refit. `apeglm` is the DESeq2
default and the recommended estimator for a simple two-group contrast (Zhu,
Ibrahim & Love 2019, Bioinformatics 35:2084). `ashr` is reserved for contrasts
between two non-reference levels and is not used here.

The two tables are not combined into a single filter. Shrunken log2FC is used
to rank and plot genes, and is reported alongside the thresholded `padj`, but
no gene is selected on the basis of a shrunken log2FC. Selecting on a shrunken
estimate while reporting FDR from a different test would break FDR control,
which is the failure mode the primary analysis section already prohibits.

#### Unpaired comparison

Analysis 2 uses a separate fit on the same `counts_filtered` matrix, in
`R/08-unpaired-comparison.R`. The gene universe is not recomputed under
`~ condition`; it is fixed by the paired-design `filterByExpr()`. This keeps
the number of genes entering the FDR correction identical across designs,
which is the point of disabling independent filtering.

#### Subsampling

Analysis 5 is stochastic and is fixed in `R/09-subsampling.R` as:

- `set.seed(20260918)` at the start of the block.
- 200 random subsets per value of n.
- n takes the values 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 113. At
  n = 113 there is a single distinct subset, so the 200 repeats collapse to
  one; this serves as an anchor against the primary result.
- Each subset is drawn once per repeat, without replacement, and the same
  patient subset is passed to both the paired and the unpaired design. The
  paired and unpaired fits therefore see identical patients within a repeat,
  which is required for the comparison to be valid. Sampling is over patient
  IDs, not rows.
- `filterByExpr()` is applied once on the full 113-patient cohort, and the
  resulting gene universe is held constant across all values of n.

At small n, the paired design has few residual degrees of freedom after fitting
patient effects, and dispersion estimation may fail to converge for some genes
or for the whole fit. Failures are recorded; the number of genes called at that
n is reported as NA if the fit fails, and the failure is noted. This is
expected and is part of what the curve shows.

#### Leave-one-patient-out

Analysis 4 is deterministic and is run in `R/10-leave-one-out.R`. Each refit
recalculates size factors and dispersion on the reduced dataset; the gene
universe is held fixed at the primary analysis filter, and "remains
significant" means `padj < 0.05` from the refit.

#### Direction of the difference in analysis 2

"Direction of the difference" is specified in `R/08-unpaired-comparison.R` as
three reported quantities:

1. **Count:** number of genes at FDR < 0.05 under the paired design minus the
   same number under the unpaired design. Sign and magnitude.
2. **Overlap:** Jaccard index between the two gene sets, plus the number of
   genes unique to each design.
3. **Sign concordance:** among genes significant in both designs, the fraction
   whose log2FC has the same sign. A low value would indicate the unpaired
   design is detecting different genes, not just fewer. Reported as NA if no
   genes are shared, which would itself be a reportable result.

These three are reported together. A single number is not sufficient.

#### Expected findings

This section records what the analysis is expected to show, so that a null or
surprising result is recognisable as such rather than as a pipeline failure.

**Magnitude of the primary list.** For a paired tumor versus adjacent-normal
comparison in a solid tumor with a well-powered cohort (n > 100 pairs), the
expected number of genes at FDR < 0.05 is in the range of several thousand
(order 10^3 to 10^4), given whole-transcriptome coverage and typical
tumor-normal expression divergence. A result in the tens would indicate a
pipeline problem, not a biological absence of signal. If the primary list is
very small, that does not mean the biology is absent; it means the analysis
should be checked before interpretation.

**Known positive controls for breast.** Under the assumption that the
113-patient cohort is predominantly ER-positive (consistent with the ER
78+/21- receptor counts above), the following genes are expected to appear as
significantly up-regulated in tumor:

- `ESR1`, `GATA3`, `FOXA1`, `XBP1` — luminal / ER-driven program
- `MKI67`, `TOP2A`, `AURKA`, `BIRC5` — proliferation
- `CCNB1`, `CDC20` — cell cycle

And the following as significantly down-regulated in tumor relative to
adjacent normal, i.e. higher in the adjacent-normal sample:

- `ADIPOQ`, `FABP4`, `PLIN1` — adipose content of normal breast
- `LEP`, `CFD` — adipocyte markers

These are the negative controls. They are chosen because normal breast tissue
is predominantly adipose, so a tumor-versus-normal contrast should show them
lower in tumor. Immune markers are deliberately not listed: breast tumors
often carry more immune infiltrate than the adjacent-normal sample, so the
expected direction is not obvious in advance, and a mis-predicted control is
worse than no control at all.

If a majority of these controls are absent from the primary list, that is a
result to investigate before interpreting the full list. It does not invalidate
the list, but it does mean the list cannot be read biologically until the
absence is understood.

**Paired versus unpaired.** The paired design is expected to call more genes
than the unpaired design at the same FDR, because it removes between-patient
variance from the residual. If the paired design calls fewer, or the same
number, that is the methodological result and is reported as such (as stated
in the Reporting rules above).

#### Reproducibility and code paths

- Dependencies: DESeq2 (estimation and testing), apeglm (shrinkage), edgeR
  (`filterByExpr` only). `ashr` is not required; only `apeglm` is used for
  shrinkage. `renv.lock` captures exact versions if `renv` is used; otherwise
  `sessionInfo()` does.
- R version and package versions are captured with `sessionInfo()` and written
  to `docs/sessionInfo.txt` from `R/11-save-results.R`, immediately before the
  analysis begins. `sessionInfo` is versioned, because it is evidence of the
  environment that produced the results.
- Gene lists are written to `docs/results/` from `R/11-save-results.R`:
  - `docs/results/primary-paired.csv`
  - `docs/results/secondary-paired-highconfidence.csv`
  - `docs/results/unpaired-comparison.csv`
  - `docs/results/subsampling.csv`
  - `docs/results/leave-one-out.csv`

`results/` is excluded by `.gitignore`. The plan body states that gene lists
are published alongside the code, so they must live in a versioned path. Before
the first commit that writes these files, confirm two things:

- that `docs/results/` is not itself matched by a `.gitignore` rule;
- that the directory exists and is tracked.

See the 2026-09-18 sizing amendment for the anchoring fix that resolves this.

#### Consolidation and editorial corrections

The original document contained a duplicated `## Amendment log` header, the
first with an empty table. That was a template artefact. The two headers have
been merged into one, without changing any content.

Two typos were present in the 2026-09-18 tissue-selection entry as originally written: "CVS" for "CSV", and "PAM50 emptu" for "PAM50 empty". Both are corrected in place here. They are typographic, not substantive, and are noted only so that the correction is visible rather than silent.

### 2026-09-18 — Computational sizing of the planned analyses

The subsampling and leave-one-patient-out analyses are sized before
running. The parameters below are provisional and subject to revision
once a single fit has been timed on the real data.

**Cost to be measured.** A paired fit on the full cohort estimates one
coefficient per patient plus the condition effect, on the full gene
universe. Before `R/09-subsampling.R` is run at full size, a single
paired fit will be timed with `system.time()` and the result recorded in
`docs/sessionInfo.txt`. If the projected total exceeds the available
budget, `n_repeats` is reduced further in a follow-up amendment.

**Provisional parameters.** The subsampling curve is computed with:

- `n_repeats = 20` per value of n, instead of the 200 originally
  considered. Twenty replicates give a readable band of variability
  around the median; the curve shape does not need 200.
- n values `c(5, 10, 20, 40, 70, 113)`, instead of twelve values. The
  removed points (30, 50, 60, 80, 90, 100) are in the region where the
  curve is smooth and add no information.
- At n = 113 there is a single distinct subset, so `n_repeats` collapses
  to 1 for that value.

**Provisional cost estimate.** 5 values of n at 20 repeats each = 100
iterations, plus 1 for n = 113 = 101 iterations per design. Two designs
= 202 DESeq2 fits, plus 113 leave-one-patient-out refits = **315 fits
total**. This number is provisional and will be updated after the
timing run, with the measured per-fit time recorded alongside it.

**What is unchanged.** The seed, the sampling over patient IDs, the
shared subset across designs, the fixed gene universe, and the
reporting rules are as described in the previous amendment. Only the
number of replicates and the grid of n values are revised.

**This is an analysis decision, not a technical one.** Fewer replicates
mean a less precise band around the curve, and it is declared here rather
than changed silently after the fact.

#### .gitignore anchoring

`results/` in `.gitignore` matches at any depth, so `docs/results/` was
excluded despite the intent to version it. The pattern is anchored to the
root as `/results/`, which matches only the top-level directory.
Verified with `git check-ignore -v docs/results/primary-paired.csv`
returning nothing.

#### Script numbering

The Script layout section above lists the scripts with a load-counts
script at position 02. The scripts themselves were previously numbered
starting at 02 for build-coldata. The plan is the authoritative
numbering and the scripts are renumbered to match:

| Position | Script | Purpose |
|---|---|---|
| 00 | `R/00-setup.R` | libraries, seed |
| 01 | `R/01-count-pairs.R` | paired patient counts (existing) |
| 02 | `R/02-load-counts.R` | recount3 count matrix, aligned to 113 patients |
| 03 | `R/03-build-coldata.R` | sample table, factor levels, coef name |
| 04 | `R/04-filter-genes.R` | `filterByExpr()` pre-filter |
| 05 | `R/05-fit-models.R` | paired and unpaired fits |
| 06 | `R/06-primary-analysis.R` | primary results |
| 07 | `R/07-secondary-analysis.R` | thresholded test and shrinkage |
| 08 | `R/08-unpaired-comparison.R` | unpaired results, comparison |
| 09 | `R/09-subsampling.R` | subsampling curve |
| 10 | `R/10-leave-one-out.R` | per-patient stability |
| 11 | `R/11-save-results.R` | CSVs and sessionInfo |

Every reference in the sections above points to this numbering.

#### Subsample pre-check

Before `R/09-subsampling.R` runs at full size, a single paired fit is
timed with `system.time()`. The result is recorded in
`docs/sessionInfo.txt`. If the projected total time for the 315 fits
exceeds a reasonable budget, `n_repeats` is reduced further and the
reduction is recorded in a follow-up amendment.

### 2026-09-21 — Counts transform and filtering rule (decisions pending)

Two decisions required before any count matrix is downloaded. They are
recorded here as open, with the measurement that will resolve each.

#### Which counts

recount3 does not serve raw read counts by default. Its base assay
derives from coverage (AUC), and the package exposes at least two
transforms to convert it: `transform_counts()` and `compute_read_counts()`.
DESeq2 assumes count-like values, so the choice matters.

The selected transform will be declared here before `R/02-load-counts.R`
runs. The candidates and what the recount3 vignette says about each will
be consulted first, and the reason for the choice recorded alongside it.
Until then, `R/02-load-counts.R` contains a placeholder marking where the
transform is applied.

#### `filterByExpr()`: design or group

The pre-filter is currently written to pass the full paired design matrix
to `edgeR::filterByExpr()`. That matrix has 115 columns for 226 samples
— near-saturated — and `filterByExpr` infers its minimum group size from
the design. With a matrix this saturated the inferred size may be very
small, which makes the filter permissive and lets through genes expressed
in a handful of samples.

The alternative is to pass the biological factor directly:

```
keep <- edgeR::filterByExpr(counts_raw, group = coldata$condition)

```

Neither is declared correct here. The two counts will be compared on the
real matrix:

```
sum(filterByExpr(counts_raw, design = design_paired))
sum(filterByExpr(counts_raw, group = coldata$condition))

```

If the retained gene counts differ materially, the chosen filter and its
reason are recorded in a follow-up amendment. This comparison belongs
with analysis 7 (sensitivity to the filtering rule) rather than as an
ad-hoc decision.

Until both decisions are made and recorded, `R/04-filter-genes.R` retains
the design-matrix form. The script will be updated alongside the
amendment that closes this entry.


### 2026-09-23 — Deduplication measured on the count matrix

**Superseded by the "Counts transform closed" entry below.** The figures in
this entry were measured on a run that used `transform_counts()`, which that
entry withdraws. They are preserved here as a record of what was measured at
the time; the reproducible figures are in the entry below.

The five rules of 2026-09-16 were applied to the real assay after
`transform_counts()`. Measured counts:

- Barcodes before rule 3: 1256
- Unique aliquots after rule 3 (sum of resequenced runs): 1232
- Groups with more than one aliquot after rule 3: 12 (7 with two, 5
  with three)
- Groups with more than one aliquot among the 113 paired patients: 4
- Aliquots after rule 4 (deepest per patient × type): 1215
- All 12 duplicate groups are sample type "01" (primary tumour)

The "4 of the 113" stated in the 2026-09-16 entry is confirmed. The 8
additional duplicate groups belong to patients who do not enter the
final paired cohort, so they do not affect the analysis. The order
applied was rules 3 → 4 → 2 (filter to "01"/"11") → pairing, so rule 4
also processed duplicates that were subsequently discarded. This is
correct but worth noting: the 1215 figure is larger than the final
cohort by construction.

### 2026-09-23 — Counts transform closed: compute_read_counts, no scaling

The 2026-09-21 entry left the counts transform open, listing three
candidates. It is now closed.

**Decision.** `compute_read_counts(rse)` is applied to the RSE, and its
result is used directly. No library-size scaling is applied anywhere in
`R/02-load-counts.R`.

**Why not `transform_counts()`.** It rescales each column by its AUC —
total coverage, including intergenic — to a target size. The resulting
`colSums` of the gene-level matrix therefore measure the share of
coverage falling in annotated genes, not sequencing depth: measured on
BRCA they ranged from 9.95e6 to 42.5e6 (median 39.4e6), a spread driven
by gene-capture fraction rather than by sequencing depth. This is
incompatible with two of the curation rules:

- **Rule 3** sums the runs of a resequenced aliquot. Each run is
  independently rescaled to the target before summing, so the summed
  aliquot lands at roughly twice the target — double weight, the
  opposite of what the 2026-09-16 entry states the rule does.
- **Rule 4** selects the aliquot with the greatest sequencing depth. On
  a rescaled matrix `colSums` no longer measures sequencing depth, so
  rule 4 would rank aliquots by gene-capture fraction — a different
  quantity from the one the rule names.

On unscaled read counts the library sizes range from 12.1e6 to 176.4e6, a
14.6-fold spread. That is the depth variation rule 4 is written to exploit.

**Why `compute_read_counts()`.** It returns estimated read counts, derived
by dividing the coverage AUC by the average mapped fragment length from
the recount3 QC metadata (`recount_qc.star.average_mapped_length`). They
are not the counts produced by the aligner, but they are on the natural
scale of sequencing depth — the scale the curation rules require and the
scale DESeq2 expects.

**Why not scaling at all.** DESeq2 estimates size factors internally from
the counts it receives. Passing pre-scaled values breaks that correction.

**Measured outputs after the change.**

- RSE dimensions: 63856 genes × 1256 runs
- Read-count assay max (single gene in a single sample): 19,759,117
- Library sizes over all 1256 runs: 12,078,452 – 176,403,775
- Unique aliquots after rule 3: 1232
- Duplicate groups after rule 3: 12 (7 with two, 5 with three)
- Duplicate groups by sample type: 01 = 12
- Aliquots after rule 4: 1215
- Patients with only a metastatic tumour: 0
- Aliquots after rule 2: 1208
- Paired patients: 113
- Duplicate groups among paired patients: 4
- Final samples: 226
- `counts_raw` dimensions: 63856 × 226

All figures above, except the `transform_counts()` quartiles cited under
"Why not", are printed by `R/02-load-counts.R` at run time. The quartiles
are a one-off measurement of the previous, incorrect run, recorded here
as evidence for the decision and not regenerable from the current script.

**On the deduplication counts.** The values 1232 → 1215 → 1208, the 12
duplicate groups, and the 4 among paired patients are identical to those
obtained with the previous `transform_counts()` run. This is expected:
rule 4 keeps exactly one aliquot per patient × type, so the number of kept
aliquots is invariant by construction. What could differ is *which*
aliquot is kept in each of the 12 duplicate groups. That comparison was
not performed. The correction does not depend on the two versions
disagreeing: it is required because `transform_counts()` ranks aliquots by
a quantity (gene-capture fraction) that is not the quantity rule 4 names,
and because DESeq2 must receive unnormalized counts.

### 2026-09-24 — Timing measured

**Closes the 2026-09-18 sizing entry.** A single paired fit on the full
gene universe was timed with `system.time()`, before any optimisations.

**Environment.** 63,856 genes, 226 samples, design `~ patient + condition`
(114 coefficients). BLAS multithreaded, ~7 effective threads.

**Measured.**

- elapsed: 25,524 s (7.09 h)
- user: 94,617 s
- system: 87,996 s
- ratio (user + system) / elapsed: ~7.2, confirming ~7 BLAS threads

The `system` time of 24.4 h indicates sustained memory pressure during the
fit: the host had 4 GB of swap, all of it in use. This is not intrinsic to
the model but is worth recording as a property of the environment.

**Projected totals with these parameters.**

- Leave-one-out (113 fits): 801 h ≈ 33 days
- Subsampling (202 fits): 1,432 h ≈ 60 days
- Total: 93 days. Not viable.

**Two optimisations are required before the analyses run at full size.**

1. `filterByExpr()` on the paired design, reducing the gene universe to
   the filtered set. Expected factor ~3–4×. This is already planned as
   `R/04-filter-genes.R`; the timing above was measured without it.
2. `DESeq(dds, parallel = TRUE)` with `MulticoreParam(workers = 6)`.
   Parallelises `nbinomWaldTest` across genes, which dominates the
   runtime. Expected factor ~5× on that phase.

Combined, the projected per-fit cost is ~25 min, giving ~2 days for
leave-one-out and ~3.5 days for subsampling. The `n_repeats` and
`n_values` parameters currently recorded in the 2026-09-18 sizing entry
assume the optimised cost; they are not viable at the unoptimised cost.

**Pending measurement.** The optimised per-fit time will be measured on a
filtered subset (e.g. 5000 genes with `parallel = TRUE`) before the full
subsampling and leave-one-out analyses are launched. That measurement,
and any further revision to the parameters, will be recorded in a
follow-up entry.

### 2026-10-06 — Computing environment moved to aarch64

**Reason.** The x86_64 machine used until 2026-10-05 is no longer
available. From this date the work runs in a container on an Apple M1
laptop (aarch64, 8 CPUs, 8 GB of RAM, of which 5.3 GiB are allocated to the
Docker virtual machine). No analysis output (gene list or effect estimate)
has been recorded: the counts were loaded and the cohort built on
2026-09-23, and one unoptimised fit was run on 2026-09-24 only to be timed.

**Environment.** `bioconductor/bioconductor_docker:RELEASE_3_23`, pinned in
`Dockerfile` by the digest of the multi-architecture index, so the same
file builds on x86_64 and aarch64. The tag is mutable and no longer
resolves to the digest recorded on 2026-09-14, so the two environments were
compared directly. `docs/sessionInfo.txt` (x86_64, 2026-09-24) and
`docs/sessionInfo-arm64.txt` (aarch64, 2026-10-06) both report R 4.6.1,
OpenBLAS 0.3.26 and LAPACK 3.12.0. Of the 78 packages listed in both
files, 77 have identical versions; `XML` differs (3.99-0.24 and
3.99-0.25). These files list loaded namespaces only, so this establishes
equality for those packages and not for the whole library.

**Data stage validated.** `R/02-load-counts.R` was run in the new
container on 2026-10-06. Every figure it prints matches the 2026-09-23
run:

- RSE: 63,856 genes x 1,256 runs
- read counts: maximum 19,759,117; library size 12,078,452 to 176,403,775
- after rule 3: 1,232 aliquots; 12 duplicate groups (7 with two, 5 with
  three), all of sample type 01
- after rule 4: 1,215 aliquots
- patients with only a metastatic tumour: 0
- after rule 2: 1,208 aliquots
- paired patients: 113; duplicate groups among them: 4
- `counts_raw`: 63,856 genes x 226 samples

Scope. `counts_raw` is not saved to disk, so the matrices were not
compared element by element, and which aliquot was kept in each duplicate
group was not compared. The figures agree; that is evidence, not proof of
identity.

**Consequences.**

1. Floating-point results (dispersion estimates, coefficients, p-values)
   are not expected to be bit-identical across platforms. Integer count
   operations are expected to match, and did. Every result is reported
   together with the platform that produced it.
2. The 2026-09-24 timing was measured on x86_64 with ~7 BLAS threads and
   4 GB of swap. It stands as a measurement of that environment, not of
   this one, and the 93-day projection is not carried over to this
   machine. The optimised timing pending in that entry will be measured on
   aarch64, recorded with the number of workers and the memory available
   to the container.
3. `R/02` completed within the 5.3 GiB of the virtual machine; its peak
   memory was not measured. The number of workers for `parallel = TRUE`
   will be chosen from a measured peak. The `workers = 6` in the
   2026-09-24 entry was set for the previous host.
4. `docs/sessionInfo.txt` is the evidence for the 2026-09-24 timing and is
   frozen; `docs/sessionInfo-arm64.txt` records the environment validated
   on 2026-10-06. Earlier entries say `R/11-save-results.R` writes
   `docs/sessionInfo.txt`. That is amended: `R/11` will write
   `docs/sessionInfo-analysis.txt`, so neither record is overwritten.

**Unchanged.** The design, thresholds, filtering rules, sample selection
rules and the counts decision (2026-09-23) are not affected by this entry.

### 2026-10-06 — Correction to the 2026-09-24 entry

The 2026-09-24 entry is left as written. This entry corrects it.

**What was timed.** `DESeq()` was called once, with default arguments, on
the full matrix (63,856 genes, 226 samples). `parallel = TRUE` was not
used, so the fit was serial in R. The ratio (user + system) / elapsed =
7.2 comes from the linear algebra library: BLAS kept ~7 threads busy. The
timing covers `DESeq()` as a whole; size-factor estimation, dispersion
estimation and the Wald test were not timed separately.

**Withdrawn from that entry.**

1. That `nbinomWaldTest` dominates the runtime. It was not measured; which
   phase dominates is unknown.
2. The description of `parallel = TRUE` as parallelising `nbinomWaldTest`.
   It is incomplete. In the DESeq2 source (`DESeq2:::DESeqParallel`) the
   genes are split among workers twice: for the gene-wise dispersion
   estimates, and for the final dispersion estimates together with the
   Wald test. The dispersion trend is fitted serially in between. Read in
   the development version (1.53.6) and confirmed in the installed
   version (1.52.0).
3. The expected factor of ~5 from six workers. The baseline already kept
   ~7 threads busy, so a factor of ~5 over it is not supported by the
   measurement. The combined projection derived from it (~25 min per fit,
   ~2 days for leave-one-out, ~3.5 days for subsampling) is withdrawn
   with it.
4. That the 24.4 h of `system` time indicate memory pressure. A
   `docker stats` reading taken about 3 h 40 min into that run shows the
   container using 4.287 GiB of a 15.49 GiB limit at 757 % CPU. That does
   not show memory pressure inside the container; swapping at host level,
   caused by other processes, is not excluded. Overhead from
   synchronising BLAS threads is a second candidate. The cause is not
   established.

The factor of ~3-4 expected from `filterByExpr()` is retained as an
expectation, not a measurement. The elapsed, user and system times of
2026-09-24 stand as measured.

**Timing protocol, fixed before it is run.** One gene set is used for all
runs: a seeded random sample of 5,000 genes drawn from those retained by
`filterByExpr()` (the option used is declared with the result). The first
5,000 rows are not used, because they are not a representative sample of
the filtered gene universe. The runs, on aarch64:

- A: serial, default BLAS threads (the configuration of 2026-09-24)
- B: serial, `OPENBLAS_NUM_THREADS=1`
- C: `OPENBLAS_NUM_THREADS=1` with `DESeq(parallel = TRUE)` and
  `MulticoreParam(workers = k)`, for k = 4 and k = 6

In A and B the three phases are timed separately: `estimateSizeFactors()`,
`estimateDispersions()` and `nbinomWaldTest()`. In C only the total is
timed, because `parallel` is an argument of `DESeq()` and not of the
phase functions. For every run, elapsed, user and system time are
recorded (self and child processes), with the peak memory of the
container from `docker stats`. A run that exhausts memory is recorded as
such. The per-fit cost for the filtered gene universe is extrapolated
from the best configuration only, and labelled as an extrapolation.

**Consequence for the sizing.** The `n_repeats` and `n_values` of the
2026-09-18 entry remain provisional. The leave-one-out (113 fits) and
subsampling (202 fits) analyses are not launched until that timing is
recorded. If the measured cost does not fit the budget, the reduction is
made in a follow-up amendment, as the 2026-09-18 entry already requires.

### 2026-10-07 — filterByExpr decision closed: group = condition

The 2026-09-21 entry left open the choice between `design =
design_paired` and `group = coldata$condition` as the filter passed to
`edgeR::filterByExpr()`. It is now closed.

**Decision.** The primary filter is
```
    filterByExpr(counts_raw, group = coldata$condition)
```
**What this changes.** The pre-specified filtering section above names
`filterByExpr()` with default parameters applied to the design matrix.
Moving to `group =` is a deviation from that pre-specification, and this
entry records the reason, as the 2026-09-21 entry requires.

**Timing of the decision.** The choice of `group =` was made before the
two filter outputs were computed. The reason was that `group =`
corresponds to the tumour-versus-normal comparison the analysis makes,
whereas `design =` is calibrated to the patient coefficients, and the
two minimum sample counts (2 and 83) had already been derived from the
design matrix and the edgeR source. The full reasoning in the "Reasons"
section below was written after the filters were run; none of it
depends on the counts of genes kept by each filter. This entry is the
first written record of the decision. There is no earlier record in the
repository; before this entry, the choice existed only in a working
session.

Chronology:

- 2026-09-23: the count matrix was first downloaded and the 113-patient
  cohort was built.
- 2026-09-24: a full `DESeq()` fit under `~ patient + condition` was
  run, for timing purposes only, as recorded in that entry. No
  differential expression results were inspected or saved from that
  fit.
- 2026-10-06: the minimum sample counts (2 for `design =`, 83 for
  `group =`) were derived. The choice of `group =` was made on that
  basis. The two filters were then run on the count matrix, and the
  figures below were obtained, confirming the derivation.

**Measured figures.** All printed by `R/04-filter-genes.R`, reproducible
by running `R/02`, `R/03` and `R/04` in sequence:

- `design_paired`: 226 rows × 114 columns
- leverage: constant at 0.504425 across all 226 samples;
  `1 / max(leverage)` = 1.98
- `keep_design` (minimum ~2 samples): 35,584 genes
- `keep_group` (minimum 83 samples): 25,836 genes
- genes kept by `group` and not by `design`: 0
- crosstab: `design` only = 9,748 · both = 25,836 · neither = 28,272
- minimum CPM threshold: 0.149
- retained after filter: 40.5 % of 63,856 genes

**Reasons.**

`design =` sets the minimum number of samples a gene must be expressed
in to `1 / max(leverage)`. In this design, every patient contributes
exactly two samples, so the design is balanced, and all leverages are
equal at 114 / 226 = 0.504425. The reciprocal gives 1.98, i.e. a gene
passes if it is expressed in 2 of the 226 samples. That minimum is
sized for the patient coefficients, 112 of the 114 in the model, each
estimated from two samples. Those coefficients are not the comparison
of interest, so `design =` calibrates the filter for a part of the
model that is not the comparison.

`group =` uses the size of the smaller group instead. With N = 113
samples per condition and N > 10, edgeR applies the rule
`10 + (113 - 10) * 0.7 = 82.1`, which means a gene must be expressed in
83 or more samples. That minimum corresponds to the tumour-versus-normal
comparison the analysis makes. `group =` ignores the design matrix and
uses only the group sizes, which is the feature of the design that maps
to the contrast.

`group =` does not use the tumour-normal difference. It counts the
number of samples in which a gene is expressed, irrespective of the
condition each sample belongs to, and the minimum depends only on the
two group sizes. A gene expressed in 80 tumours and 3 normals is kept,
as is one expressed in 83 normals and no tumour; a gene expressed in 82
samples of any composition is dropped. The pre-specified rule that the
filter never uses the quantity being tested is respected.

`group =` is a subset of `design =`: every gene kept by `group` is also
kept by `design`, because the same CPM threshold applies to both and
`group` requires a higher minimum sample count. This is a coherence
check between the two calls, not a reason to prefer one over the other.

The filtered matrix, `counts_filtered`, remains the single input to both
the paired and the unpaired `DESeqDataSet`, as the plan requires. Only
its gene universe changes; the sample columns are the same 226.

**Cost.** 9,748 genes that `design =` keeps are dropped by `group =`.
These are genes expressed in fewer than 83 samples in total, regardless
of condition. This is the cost of the stricter threshold. A gene
confined to tumours is kept only if it is expressed in at least 83 of
the 113 tumours, which is 73 % of the tumour samples. A gene present in
half the tumours and in no normal is discarded.

**Correction to the 2026-09-21 entry.** That entry states the paired
design matrix has 115 columns. The measured value is 114:
1 intercept, plus 112 patient indicators (113 patients, one as
reference), plus 1 condition indicator. The 115 figure is not edited in
place; it is corrected here.

**Analysis 7, redefined.** The plan currently specifies analysis 7 as
`filterByExpr()` defaults versus a manual rule of at least 10 counts in
at least N samples. It is redefined as the primary filter
(`group = condition`) versus the paired-design filter
(`design = design_paired`). The comparison requires one additional
paired fit on the 35,584 genes kept by `design =`, because changing the
gene universe changes dispersion estimation and the multiple-testing
correction; it is not a re-use of the fit under the primary filter. The
reported quantities are: how many genes are significant in each fit at
the primary FDR threshold, the overlap between the two lists, and how
many of the 9,748 genes unique to `design =` are significant in that
fit. This is one additional fit, not a re-run of the subsampling or
leave-one-out analyses.

**Consequences.**

- The gene universe for the analysis is 25,836 genes, down from 63,856
  — a factor of 2.47. The 2026-10-06 entry recorded a factor of 3–4 as
  an expectation for the filter; the measured factor replaces it. This
  is a reduction in the number of genes entering each fit, not in
  wall-clock time; the timing protocol measures the time separately.
- The timing protocol (2026-10-06) uses the `group =` filter to define
  the 5,000-gene sample it times on.
- The filter is `group =`, not `design =`, in `R/04-filter-genes.R`. The
  `design =` form is retained in that script for the analysis 7
  comparison.


### 2026-10-15 — Timing results: serial versus parallel DESeq() fits

This entry records the timing measurements. The configuration chosen on
the basis of them is recorded in the entry below.

**Protocol.** All 5,000-gene measurements use the same inputs.

- 5,000 genes randomly sampled from the 25,836 that pass the filter, with
  seed 20260918 and `sort()`, so the 5,000-gene configurations time the
  same set. The 2,500-gene and 25,836-gene runs use different sets.
- Two scripts: `R/timing-serial.R` (three phases timed separately) and
  `R/timing-parallel.R` (a single `system.time()` around
  `DESeq(parallel = TRUE)`).
- Configurations: A (default BLAS threading, ~8 threads), B (single BLAS
  thread), C (`DESeq(parallel = TRUE)` with 4 or 6 workers).
- The image is built on each machine from the same Dockerfile, with the
  base pinned by digest. The Bioconductor image is not used directly.
- Container memory sampled with `docker stats` every ~5–7 s. Container
  peaks reported below are lower bounds. R peaks come from `gc()` and are
  exact, but cover the main process only, not the workers.

#### Mac (aarch64), runs on 10-08 and 10-09

The Docker VM memory was raised from 5.3 GiB (as recorded on 2026-10-06)
to 6.27 GiB before these runs.

| Config | Genes | elapsed | result |
|---|---|---|---|
| A | 5,000 | 1,704 s | — |
| B | 5,000 | 4,098 s | — |
| C, 4 workers (run 1) | 5,000 | — | OOM |
| C, 4 workers (run 2) | 5,000 | — | OOM |
| C, 6 workers | 5,000 | — | OOM |

The `dmesg` line from run 2 is in
`docs/timing/aarch64/dmesg-k4-rep.log`. No full-universe fit was run on
the Mac, serial or parallel.

#### Laptop (x86_64, i7-8550U, 15.5 GiB), runs on 10-10 (4-worker run 2: 10-15)

Docker on this machine runs without a VM and without a per-container
memory limit.

| Config | Genes | elapsed | ratio | R peak | Container peak |
|---|---|---|---|---|---|
| A | 5,000 | 1,495 s | 7.93 | 2,342 MB | 3.4 GiB |
| B | 5,000 | 489 s | 1.00 | 2,344 MB | 3.1 GiB |
| C, 4 workers (run 1, 10-10) | 5,000 | 245.7 s | 3.69 | 1,713 MB | 8.4 GiB |
| C, 4 workers (run 2, 10-15) | 5,000 | 224.9 s | 3.67 | 1,713 MB | not measured |
| C, 6 workers | 5,000 | 215.5 s | 4.88 | 1,714 MB | 9.0 GiB |
| C, 4 workers | 2,500 | 115.3 s | 3.61 | 1,524 MB | 7.3 GiB |
| **B** | **25,836** | **2,551 s** | **1.00** | **2,570 MB** | **3.4 GiB** |

The A and B container peaks are both during the R/02 data load, not
during the fit itself.

The two 4-worker runs at 5,000 genes differ by about 9% in elapsed time
(245.7 s on 10-10, 224.9 s on 10-15): same configuration and machine, on
different days. Run 2 was made to measure run-to-run variability, with
the machine set up as on 10-10. The differences between configurations
(A vs B, 3.1×; B between machines, 8.4×) are far larger than this.

Scaling: from 2,500 to 5,000 genes under C with 4 workers, 2.13× (run 1)
or 1.95× (run 2) against an expected 2×, both within run-to-run
variability. From 5,000 to 25,836 genes under B, 5.22× against an
expected 5.17×.

The full-universe B fit was the primary model (`~ patient + condition`),
run only to time it. The script does not print DESeq2 results and does
not save the fitted object.

#### What remains unexplained

B is 8.4× faster on the laptop than on the Mac (489 s vs 4,098 s at 5,000
genes), while A is roughly comparable (1,495 s vs 1,704 s). The reason a
single BLAS thread behaves so differently across platforms, while the
multithreaded configuration does not, is not established.

#### The 2026-09-24 run

The 2026-09-24 measurement (7.09 h, line 717 of this document) was taken
on this same laptop, before the gene filter was decided and with ~7 BLAS
threads. Its `system` time pattern is consistent with configuration A.
The mechanism is not demonstrated for that specific run.

The plan did not record which image was used for that run. `docker
inspect` on its container shows it was created from
`bioconductor/bioconductor_docker:latest`, image digest
`sha256:979991813693d4256f685a565cd59560d7bca7c4d4d41f66ed5fc38b31338594`
(a single-platform amd64 image). The pinned environment recorded on
2026-10-06 (line 761) applies from that date on.

#### Cost of the full plan (estimate)

The sizing recorded in the 2026-09-18 entry (line 478), declared
provisional by the 2026-10-06 correction until this timing was
recorded, uses `n_repeats = 20` and n values `c(5, 10, 20, 40, 70, 113)`,
with n = 113 collapsing to a single repetition. That gives 101 paired
fits and 101 unpaired fits. Adding 113 leave-one-out fits and one
analysis-7 fit gives 316 fits total.

To convert these into full-fit equivalents, the cost per fit is modelled
as proportional to `samples × coefficients² + coefficients³`, the
standard cost shape for an iteratively reweighted least-squares fit.
Under this model, weighted by each fit's sample count, coefficient count
and gene count:

- Leave-one-out (113 fits, 224 samples each, 113 coefficients,
  25,836 genes): ~110 full-fit equivalents
- Subsampling, paired (101 fits, n = 5 to 113 patients, coefficients from
  6 to 114, 25,836 genes): ~7 full-fit equivalents
- Subsampling, unpaired (101 fits, n = 5 to 113 patients, 2
  coefficients, 25,836 genes): <1 full-fit equivalent
- Analysis 7 (1 fit, 226 samples, 114 coefficients, 35,584 genes): ~1.4
  full-fit equivalents
- **Total: ~120 full-fit equivalents**

This is an estimate from an unvalidated cost model, not a measurement.
It does not account for differences in convergence between fits. It also
ignores per-gene work that does not depend on the number of coefficients,
so small fits are underestimated: each of the 101 unpaired fits still
processes 25,836 genes. The total is approximate.

At 42.5 min per full-fit equivalent (measured, configuration B on the
full universe):

- Serial, one fit at a time: **~85 h**.
- Concurrent, if running several fits at once yields ~2× throughput:
  **~43 h**. The 2× factor is an assumption, not measured. Measuring it
  is the first step of the execution script.

### 2026-10-15 — Computing configuration for the analysis

**Decision.** The analysis runs on the x86_64 laptop. Each fit is a normal
serial `DESeq()` with a single BLAS thread (`OPENBLAS_NUM_THREADS=1`).
Multiple fits run concurrently, each in its own container with its own
memory limit (`docker run --memory`). Each fit writes its result to disk
on completion.

The isolation is a design requirement, not an observed property. Without
per-container limits, a memory kill on Linux can affect any process on the
machine, including another running fit. The limit keeps the failure inside
the fit that caused it.

The results are identical to running the fits one at a time. Each fit is
an independent `DESeq()` call on its own data, with no shared state. How
many fits run at once does not change any result; it will be fixed when
the execution script is written, bounded by available memory (at most 3
with the current machine and dataset).

**Reasons.**

*Speed.* Measured on x86_64 with 5,000 genes and the same random sample:
A (default BLAS threading, ~8 threads) took 1,495 s and B (single BLAS
thread) took 489 s. A is 3.1× slower than B. System time in A was 5,541 s
against 6,310 s user time over the 5,000-gene run, i.e. about half the CPU
spent in the kernel. That is consistent with coordination between BLAS
threads, though the mechanism is not demonstrated.

On the Mac (aarch64) the ordering is reversed: A 1,704 s, B 4,098 s. The
reason for the difference between platforms is not established here.

*Memory.* A single serial fit peaks at ~2.6 GB of R memory on the full
gene universe. The container itself reaches 3.4 GiB during the R/02 data
load, which is the figure that matters for the `--memory` limit. C
(`DESeq(parallel = TRUE)`) measured, on x86_64 at 5,000 genes, 8.4 GiB
with 4 workers and 9.0 GiB with 6 workers (container peaks sampled every
~7 s, so lower bounds), and the memory does not fall in proportion to the
gene count: 7.3 GiB at 2,500 genes. Whether a full-universe C fit would
fit in the laptop's 15.5 GiB is not established. Independent processes
with per-container limits avoid the question.

*Failure isolation.* On 2026-10-09 the Mac lost a C run at 5,000 genes
when one worker died and the whole fit went with it. Independent fits do
not have this failure mode.

**Alternatives not chosen.**

- **A (default BLAS threading).** 3.1× slower than B on x86_64 at 5,000
  genes, with about half the CPU in system time. The mechanism is
  consistent with thread coordination but not measured.
- **B (one fit at a time, single BLAS thread).** Correct and lighter than
  the chosen configuration, but leaves the machine's other cores idle.
- **C (`DESeq(parallel = TRUE)`).** 8.4 GiB with 4 workers and 9.0 GiB
  with 6 workers on x86_64 at 5,000 genes, with no proportional reduction
  at smaller gene counts. Failed on the Mac when a worker died.
- **Gene-block parallelism with new code.** Would require validation for
  a marginal gain over the chosen configuration. Not pursued.
- **Mac as compute host.** Two reasons, in order of weight. First,
  consistency: the analysis compares fits with each other, and mixing
  architectures would blend each platform's floating-point rounding into
  what the comparison measures. Second, throughput: the Mac's best
  configuration projects ~2.4 h per full fit against 42.5 min measured
  on the x86_64 laptop, and its Docker VM of 6.27 GiB allows one full fit
  at a time. The Mac is not excluded because it cannot run fits — A and B
  both ran there, and a serial full fit uses ~2.6 GB (measured on
  x86_64), well within the VM. What fails on the Mac is C. A
  cross-platform verification (one full fit on each machine, compared) is
  left as future work.

**Platform.** This entry changes the 2026-10-06 decision that the Mac
would be the compute host. All fits run on the x86_64 laptop. The Mac
figures come from `docs/timing/aarch64/`, the x86_64 figures from
`docs/timing/x86_64/`. Smoke tests (200 genes) were run before the long
runs and were not saved; no figure in this entry comes from them.

The 2026-09-24 measurement was taken on this same machine under different
conditions: before the gene filter was decided, and with ~7 BLAS threads.
It is superseded by the measurements in this entry.

**Timing of the decision.** The configuration was chosen on 2026-10-10,
after the measurements recorded in the preceding entry. This is the first
written record of the choice.


### 2026-10-16 — Sizing of the subsampling and leave-one-out analyses

**Context.** The 2026-09-18 entry left the subsampling parameters
provisional pending a `system.time()` measurement on a full paired fit.
The 2026-10-06 correction withdrew the factor-5 expectation and blocked
both the leave-one-out and the subsampling analyses until that timing
was recorded. The timing was recorded on 2026-10-15 in `docs/timing/`
and in the entries of that day, not in `docs/sessionInfo.txt` as the
2026-09-18 entry had stated. This entry revises the parameters of both
analyses.

**When this is decided.** Before any differential expression result.
Two full DESeq() fits have been run so far, both for timing only: the
2026-09-24 fit (interactive session) and the 2026-10-10 fit
(`R/timing-serial.R`, full gene universe). Neither printed nor saved
DESeq2 results, and no gene list has been inspected.

**Band.** The 2026-09-18 entry specified `n_repeats` and the grid of n
values but not the band reported around the curve. It is defined here:
the median with a 10th-to-90th percentile band. An 80% central band
excludes the tails, which are sensitive to single extreme draws at the
replicate counts considered here. The 10-90 band gives a readable
summary of the outcomes the analysis can produce without letting
outliers set its edges.

**Repetitions.** `n_repeats` rises from 20 to 50. With 50 replicates,
about 5 values fall beyond each edge of the 10-90 band; with 20, about
2. At 5 values, a single extreme draw moves the edge less than at 2.
The n grid does not change: `c(5, 10, 20, 40, 70, 113)`; at n = 113
there is a single distinct subset, so that value has one replicate.

**Leave-one-out.** The analysis runs on all 113 patients, as specified
in analysis 4 (line 92) and its implementation section (line 374): refit
the paired model dropping each patient in turn, and report for each
called gene the fraction of refits in which it remains significant,
where "remains significant" means `padj < 0.05` from the refit and the
gene universe is fixed at the primary analysis filter. The 113 refits
are not independent: every pair shares 111 patients. What each refit
provides is the only answer for its own patient: whether the primary
list depends on that patient being included. That is the question the
analysis asks, and it requires all 113. Dropping a subset of patients
would answer a different question — whether the list survives a random
subsample — and is not what is reported here.

**Budget.** No total time cap. Runs are organised in batches of up to
14 hours per day, with at least 10 hours of rest between them. If the
measured cost exceeds the estimate, nothing is cut: the analysis takes
more days. This rule does not depend on the measured cost, so no result
can influence it.

**Resulting size.** 616 fits in total:

- paired subsampling: 5 values of n × 50 repeats + 1 at n = 113 = 251
- unpaired subsampling: same grid and repeats = 251
- leave-one-out: 113
- analysis 7 (filter comparison): 1 additional paired fit

Projected to full-fit equivalents from the measured per-fit cost and the
cost model in the 2026-10-15 entry:

- paired subsampling: ~16 full-fit equivalents
- unpaired subsampling: under 1 by the model
- leave-one-out: ~110
- analysis 7: ~1.4
- **total: ~127 full-fit equivalents**

The model underestimates the unpaired fits for the same reason it does
in the 2026-10-15 entry: it ignores per-gene work that does not depend
on the number of coefficients, and each of the 251 unpaired fits still
processes 25,836 genes. The unpaired cost is measured when the
execution starts.

**Duration.** At the measured 42.5 min per full-fit equivalent
(configuration B, full gene universe, measured on 2026-10-10, recorded
2026-10-15):

- serial, one fit at a time: **~90 hours**
- if running several fits concurrently yields ~2× throughput: **~45
  hours**

The 2× factor is an assumption, not measured. This is a projection from
one measured full fit and an unvalidated cost model, not a measurement
of the full analysis.

**Consequence for execution.** The subsampling script and the
leave-one-out script are both written to save each fit on completion,
skip fits whose output file already exists, and not start a new fit
when the daily window is close to closing. Under this rule a planned
stop at the end of a day loses nothing. Only an unplanned interruption
(power cut, machine hang) loses the fits that were in progress.

The 43-minute figure is one full fit in serial. Under concurrency each
fit takes longer in wall-clock terms, so the daily stop must allow for
that. The analysis-7 fit runs on 35,584 genes rather than 25,836, and
is projected at roughly 1 hour on its own.

### 2026-10-17 — Unpaired fit location and outlier handling

**Unpaired fit.** The "Unpaired comparison" section (line 343) places the
unpaired fit in `R/08-unpaired-comparison.R`. This conflicts with the
Script layout table of the 2026-09-18 entry, which assigns both fits to
`R/05-fit-models.R` and is declared authoritative for the numbering. The
table is followed: `R/05` fits both designs on the same data with the
same configuration; `R/08` extracts and compares.

**Outlier handling under DESeq2 defaults.** DESeq2 applies two
different outlier mechanisms by default, and both depend on the cell
sizes of the design matrix (a cell is a group of samples with identical
rows in the design matrix):

1. **Replacement** (`minReplicatesForReplace`, default 7): replaces
   outliers during `DESeq()`. Applies only to samples in cells with ≥ 7
   samples.
2. **Cook's distance flagging** (`cooksCutoff`, default TRUE): sets
   `padj` to NA for genes with outliers during `results()`. Applies only
   to samples in cells with ≥ 3 samples.

In the paired design (`~ patient + condition`), every sample is alone in
its cell, because every patient has their own coefficient. No cell
reaches size 2, let alone 3 or 7. Neither mechanism can ever fire.
In the unpaired design (`~ condition`), all 113 samples of each
condition share a cell. Both mechanisms fire.

This creates two problems for the planned analyses:

- **Analysis 2 (paired vs unpaired).** Comparing the two designs would
  confound the effect of the design with the effect of outlier removal.
  The Cook's flagging violates the plan's rule at lines 300-301 ("the
  number of genes entering the FDR correction is controlled by the
  pre-filter only"): it sets `padj` to NA for genes with outliers,
  affecting the unpaired design but never the paired one. The replacement
  mechanism does not remove genes from the correction but changes the
  estimates of some genes in only one design, confounding the
  comparison.
- **Subsampling curve (unpaired).** At n = 5 (5 samples per cell),
  Cook's flagging applies but replacement does not; from n = 10, both
  apply. The curve would mix two different methods, changing under the
  same line.

**Verification.** A smoke test with 200 genes was run, printing only
diagnostics and no differential expression results:

- Samples in cells ≥ 3: paired 0 · unpaired 226
- Samples in cells ≥ 7: paired 0 · unpaired 226
- Genes replaced during the fit: paired 0 · unpaired 12 (of 200)
- `cooks` assay stored: paired TRUE · unpaired TRUE

**Decision.** To keep the two designs comparable and the subsampling
curve consistent:

- `cooksCutoff = FALSE` in every `results()` call. This disables Cook's
  distance flagging. The rationale is the plan's rule at lines 300-301:
  the number of genes entering the FDR correction must depend only on
  the pre-filter, not on the design.
- `minReplicatesForReplace = Inf` in every `DESeq()` call. This disables
  replacement in both designs. The rationale is comparability between
  designs and across the subsampling curve.

**Cost.** The unpaired design loses DESeq2's protection against
outliers; this is reported as a limitation. The primary analysis (paired
design) is unaffected: it never had that protection, because its design
matrix never creates large enough cells. The decision is recorded here
before any differential expression result is inspected.

### 2026-10-21 — Non-converged genes in the paired fit

**Observation.** The paired fit (`~ patient + condition`) had 13 out of
25,836 genes that did not converge in beta (line 170 of
`docs/logs/05-fit-models.log`). The unpaired fit (`~ condition`) had 0
out of 25,836. For all 13 genes, there is at least one patient with
exactly zero counts in both their tumour and normal samples (ranging
from 1 to 9 patients per gene). Among the 25,823 genes that did
converge, exactly 606 (2.35%) share this same pattern. No gene without
this pattern failed to converge. The 13 non-converged genes are not
low-expression genes: their total read counts range from 13,063 to
1,631,217, and between 14 and 55 of the 226 samples have zero counts for
these genes.

**Test of more iterations.** The 13 genes were refitted alone with
`nbinomWaldTest()`, reusing the size factors and dispersion estimates of
the full fit. The control (`maxit = 100`, the value used in the original
fit) reproduced the failure: 0 of 13 converged. With `maxit = 500` and
`maxit = 2000`, 0 of 13 converged. Run time did not grow with `maxit`
(29, 27.5 and 29 s), which indicates that the iteration limit is not
what stops these fits. For a patient with zero counts in both samples,
the maximum-likelihood estimate of that patient's coefficient is not
finite. This pattern is necessary but not sufficient here, since 606
converged genes share it. Why these 13 fail is not investigated, because
doing so would require inspecting fitted coefficients before the primary
analysis.

**Decision.** All genes, including the 13 non-converged ones, are kept
and flagged. Excluding them would contradict lines 300-302: the number
of genes entering the FDR correction "is controlled by the pre-filter
only and is identical across the paired and unpaired fits". It would
also break the fixed gene universe of the subsampling analysis (lines
365-366) and of the leave-one-out analysis (lines 377-378): different
genes may fail to converge in each refit, so excluding them would make
the universe vary from one fit to another. Excluding them from the
paired design only would also create an artificial asymmetry in
Analysis 2 (lines 381-396): any of them called in the unpaired design
would appear as unique to it. Keeping them is consistent with lines
368-372, which already require convergence failures in the subsampling
to be recorded rather than hidden.

**Cost.** Up to 13 genes with unreliable paired estimates may enter the
called gene lists. They are flagged, and the sensitivity check measures
their impact.

**Rule for all fits.** All genes are retained. Every gene-level output
must carry a boolean column `beta_converged`. For every fit in the
subsampling (`R/09`) in both designs, leave-one-out (`R/10`), and
analysis 7, the number of non-converged genes must be recorded.

**Sensitivity check.** As a secondary view, all Analysis 2 measures
(lines 381-396) are recomputed after excluding any gene flagged in
*either* design (13 genes in the full fits). Each design's `padj` values
(each corrected over 25,836 genes) are used without recalculating
Benjamini-Hochberg. This isolates the effect of the flagged genes.
Recalculating BH over the remaining genes would confound their removal
with a slight change in the correction of all other genes. The number of
flagged genes that survive the primary FDR threshold will also be
reported.

This decision is recorded before inspecting any differential expression
results.
