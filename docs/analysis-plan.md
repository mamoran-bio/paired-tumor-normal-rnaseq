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

Before `R/09-subsampling.R` runs at full size, a single paired fit is timed:
```
system.time({
d <- DESeqDataSetFromMatrix(counts_filtered, coldata, ~ patient + condition)
d <- DESeq(d)
})

```
The result is recorded in `docs/sessionInfo.txt`. If the projected total
time for the 315 fits exceeds a reasonable budget, `n_repeats` is
reduced further and the reduction is recorded in a follow-up amendment.

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
