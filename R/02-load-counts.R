# R/02-load-counts.R
#
# Download the BRCA gene-level count matrix from recount3 and apply the
# five curation rules of the 2026-09-16 amendment.
#
# Rules: (1) patient and type from barcode, (2) keep only "01" and "11",
# (3) sum resequenced runs of the same aliquot, (4) keep the deepest
# aliquot per patient x type, (5) exclude metastatic samples.
#
# Inputs:  none (downloads from recount3)
# Outputs: counts_raw (genes x samples), sample_annotation

library(recount3)

# --- Download RSE ----------------------------------------------------------
proj_info <- subset(available_projects(),
                    project == "BRCA" & file_source == "tcga")
stopifnot(nrow(proj_info) == 1)

rse <- create_rse(proj_info, type = "gene")
stopifnot(is(rse, "RangedSummarizedExperiment"))
message("RSE: ", nrow(rse), " genes x ", ncol(rse), " runs")

# --- Read counts (unnormalized, estimated) ---------------------------------
# compute_read_counts() returns ESTIMATED read counts, derived from the
# coverage data and the recount3 QC metadata (read length). They are not
# the counts produced by the aligner, but they are on the natural scale
# of sequencing depth.
#
# Do NOT use transform_counts(): it rescales each column by its AUC
# (coverage including intergenic) to a target size, so colSums of the
# gene-level matrix measure gene-capture fraction, not depth. That is
# incompatible with two of the curation rules below — rule 3 would give
# resequenced aliquots double weight, and rule 4 would rank by capture
# fraction instead of depth — and it breaks DESeq2's internal size-factor
# estimation, which expects unnormalized counts.

assays(rse)$read_counts <- compute_read_counts(rse)
counts_reads <- assay(rse, "read_counts")
barcode      <- as.data.frame(colData(rse))$tcga.tcga_barcode

# Guards on the barcode. Rule 1 exists precisely because the barcode is
# the only complete identifier; if it were NA or truncated, rules 2-5
# would misbehave silently. substr() on a short string returns "", and
# "" would drop a sample from the analysis without a warning.
stopifnot(!anyNA(barcode))
stopifnot(all(nchar(barcode) >= 15))
stopifnot(length(barcode) == ncol(counts_reads))
stopifnot(all(counts_reads == round(counts_reads)))
message("Read counts: max = ", max(counts_reads),
        " | lib size = ", min(colSums(counts_reads)),
        " - ", max(colSums(counts_reads)))

# --- Rule 3: sum runs of the same aliquot ---------------------------------
counts_aliquot <- t(rowsum(t(counts_reads), group = barcode))
message("After rule 3: ", ncol(counts_aliquot), " unique aliquots")

# --- Rule 4: deepest aliquot per patient x type ---------------------------
bc    <- colnames(counts_aliquot)
key   <- paste(substr(bc, 1, 12), substr(bc, 14, 15), sep = "_")
depth <- colSums(counts_aliquot)

# Duplicate group diagnostics (reproducible evidence for the amendment)
dup <- table(key)[table(key) > 1]
message("Duplicate groups after rule 3: ", length(dup),
        " (", sum(dup == 2), " with two, ", sum(dup == 3), " with three)")

# Positions 14-15 of `key` hold the sample type by construction of key
# (12 + 1 + 2 = 15 characters), not because key is a barcode. The index
# coincides with the barcode's 14-15 by accident of the same layout.
dup_type <- table(substr(names(dup), 14, 15))
message("Duplicate groups by sample type: ",
        paste(names(dup_type), as.integer(dup_type), sep = "=",
              collapse = " / "))

keep_bc <- vapply(
  split(seq_along(bc), key),
  function(i) bc[i[which.max(depth[i])]],
  character(1)
)
keep_bc <- unname(keep_bc)
message("After rule 4: ", length(keep_bc), " aliquots")

# --- Rules 2 and 5 ---------------------------------------------------------
st_pre  <- substr(keep_bc, 14, 15)
pat_pre <- substr(keep_bc, 1, 12)

solo_06 <- setdiff(pat_pre[st_pre == "06"], pat_pre[st_pre == "01"])
message("Patients with only a metastatic tumour: ", length(solo_06))

keep_bc <- keep_bc[st_pre %in% c("01", "11")]
message("After rule 2: ", length(keep_bc), " aliquots")

# --- Pairing ---------------------------------------------------------------
pat <- substr(keep_bc, 1, 12)
st  <- substr(keep_bc, 14, 15)

paired_patients <- intersect(pat[st == "01"], pat[st == "11"])
message("Paired patients: ", length(paired_patients))

dup_paired <- names(dup)[substr(names(dup), 1, 12) %in% paired_patients]
message("Duplicate groups among paired patients: ", length(dup_paired))

final_bc <- keep_bc[pat %in% paired_patients]
stopifnot(length(final_bc) == 2 * length(paired_patients))
message("Final samples: ", length(final_bc))

# --- Assemble outputs ------------------------------------------------------
counts_raw <- counts_aliquot[, final_bc, drop = FALSE]

sample_annotation <- data.frame(
  sample_id = final_bc,
  patient   = substr(final_bc, 1, 12),
  condition = ifelse(substr(final_bc, 14, 15) == "01", "tumor", "normal"),
  stringsAsFactors = FALSE
)

stopifnot(identical(sample_annotation$sample_id, colnames(counts_raw)))

# The paired design requires exactly one tumour and one normal per patient.
# The total-count assert above would still pass if a patient ended up with
# two tumours and no normal; this table catches that case.
stopifnot(all(table(sample_annotation$patient,
                    sample_annotation$condition) == 1))

message("counts_raw: ", nrow(counts_raw), " genes x ", ncol(counts_raw),
        " samples")
