# R/10-leave-one-out.R
#
# Fit the paired DESeq2 model omitting one patient at a time (Leave-One-Out).
# Evaluates how many of the primary significant genes remain significant (padj < 0.05).

suppressPackageStartupMessages(library(DESeq2))

# 1. Loading inputs and primary list
message("Loading inputs...")
inputs <- readRDS("results/model-inputs.rds")
primary_list <- readRDS("results/res_primary.rds")
counts <- inputs$counts_filtered
coldata <- inputs$coldata
coef_name <- inputs$coef_name

# Keep only genes called in the primary analysis (THE FIX)
primary_sig <- primary_list[!is.na(primary_list$padj) & primary_list$padj < 0.05, ]
message("Primary significant genes to evaluate: ", nrow(primary_sig))

# 2. Reading args and pre-flight checks
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Use: Rscript R/10-leave-one-out.R <patient_id>")

omit_patient <- args[1]
message("Omit patient: ", omit_patient)

stopifnot(omit_patient %in% coldata$patient)

# 3. Subset counts and coldata to exclude the patient
keep_idx <- which(coldata$patient != omit_patient)
sub_counts <- counts[, keep_idx, drop = FALSE]
sub_coldata <- droplevels(coldata[keep_idx, , drop = FALSE])

# Sanity checks for the subset
stopifnot(ncol(sub_counts) == 224)             # 113 patients * 2 - 2
stopifnot(nlevels(sub_coldata$patient) == 112) # One patient removed
stopifnot(nlevels(sub_coldata$condition) == 2)

# 4. Fit DESeq2 model
message("Fitting model...")
t_fit <- system.time({
  dds <- DESeqDataSetFromMatrix(sub_counts, sub_coldata, ~ patient + condition)
  dds <- DESeq(dds, minReplicatesForReplace = Inf)
  res <- results(dds, name = coef_name, independentFiltering = FALSE, cooksCutoff = FALSE)
})

n_nonconv <- sum(!mcols(dds)$betaConv, na.rm = TRUE)

# 5. Evaluate stability
sig_genes <- rownames(res)[!is.na(res$padj) & res$padj < 0.05]
still_sig <- rownames(primary_sig) %in% sig_genes

message("Fit time: ", round(t_fit[["elapsed"]], 1), "s | ",
        "Still significant: ", sum(still_sig), " / ", nrow(primary_sig),
        " | Non-converged: ", n_nonconv)

# 6. Save result (Consolidated into a single list)
out_list <- list(
  genes = data.frame(
    omitted_patient = omit_patient,
    gene_id = rownames(primary_sig),
    is_sig = still_sig
  ),
  meta = data.frame(
    omitted_patient = omit_patient,
    fit_time = t_fit[["elapsed"]],
    n_still_sig = sum(still_sig),
    n_nonconv = n_nonconv
  )
)

dir.create("results/loo", showWarnings = FALSE)
saveRDS(out_list, sprintf("results/loo/loo_%s.rds", omit_patient))

message("Done.")