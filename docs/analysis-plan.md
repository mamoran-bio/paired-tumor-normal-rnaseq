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

Gene lists are published as CSV files under `results/`, with the code that
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

| Date | Change | Reason |
|---|---|---|
| | | |

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
