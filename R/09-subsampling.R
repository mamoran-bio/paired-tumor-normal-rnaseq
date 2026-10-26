# R/09-subsampling.R
#
# Fit one subsampling replicate (paired and unpaired) for a given n.

suppressPackageStartupMessages(library(DESeq2))
source("R/00-helpers.R")

# 1. Loading inputs 
inputs <- readRDS("results/model-inputs.rds")
counts_filtered <- inputs$counts_filtered
coldata         <- inputs$coldata
coef_name       <- inputs$coef_name
gene_info       <- inputs$gene_info

# 2. Arguments from bash
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Use: Rscript R/09-subsampling.R <n> <rep>")
n    <- as.integer(args[1])
rep  <- as.integer(args[2])

# 3. FIXED SEED (no collisions)
set.seed(20260918L + n * 1000L + rep)

# 4. Setting random set of paired patients samples
patients <- unique(coldata$patient)
sampled_patients <- sample(patients, size = n, replace = FALSE)
idx <- which(coldata$patient %in% sampled_patients)
sub_counts  <- counts_filtered[, idx, drop = FALSE]
sub_coldata <- droplevels(coldata[idx, , drop = FALSE])

stopifnot(ncol(sub_counts) == 2 * n)
stopifnot(nlevels(sub_coldata$patient) == n)
stopifnot(nlevels(sub_coldata$condition) == 2)

# 5. Paired fit
res_paired <- tryCatch({
  t_fit <- system.time({
    dds_p <- DESeqDataSetFromMatrix(sub_counts, sub_coldata, ~ patient + condition)
    dds_p <- DESeq(dds_p, minReplicatesForReplace = Inf)
    res <- results(dds_p, name = coef_name, independentFiltering = FALSE, cooksCutoff = FALSE)
  })
  
  res_df <- as.data.frame(res)
  res_df$gene_name <- gene_info$gene_name
  
  list(
    n_sig = sum(res_df$padj < 0.05, na.rm = TRUE),
    n_nonconv = sum(!mcols(dds_p)$betaConv, na.rm = TRUE),
    n_ctrl_up = sum(res_df$gene_name %in% controls_up & res_df$padj < 0.05, na.rm = TRUE),
    n_ctrl_dn = sum(res_df$gene_name %in% controls_dn & res_df$padj < 0.05, na.rm = TRUE),
    time = t_fit[["elapsed"]],
    error = NA_character_
  )
}, error = function(e) list(n_sig = NA, n_nonconv = NA, n_ctrl_up = NA, n_ctrl_dn = NA, time = NA, error = conditionMessage(e)))

# 6. Unpaired fit
res_unpaired <- tryCatch({
  t_fit <- system.time({
    dds_u <- DESeqDataSetFromMatrix(sub_counts, sub_coldata, ~ condition)
    dds_u <- DESeq(dds_u, minReplicatesForReplace = Inf)
    res <- results(dds_u, name = coef_name, independentFiltering = FALSE, cooksCutoff = FALSE)
  })
  
  res_df <- as.data.frame(res)
  res_df$gene_name <- gene_info$gene_name
  
  list(
    n_sig = sum(res_df$padj < 0.05, na.rm = TRUE),
    n_nonconv = sum(!mcols(dds_u)$betaConv, na.rm = TRUE),
    n_ctrl_up = sum(res_df$gene_name %in% controls_up & res_df$padj < 0.05, na.rm = TRUE),
    n_ctrl_dn = sum(res_df$gene_name %in% controls_dn & res_df$padj < 0.05, na.rm = TRUE),
    time = t_fit[["elapsed"]],
    error = NA_character_
  )
}, error = function(e) list(n_sig = NA, n_nonconv = NA, n_ctrl_up = NA, n_ctrl_dn = NA, time = NA, error = conditionMessage(e)))

# 7. Save
out_df <- data.frame(
  n = n, rep = rep,
  paired_sig = res_paired$n_sig, paired_nonconv = res_paired$n_nonconv, 
  paired_ctrl_up = res_paired$n_ctrl_up, paired_ctrl_dn = res_paired$n_ctrl_dn, paired_time = res_paired$time, paired_err = res_paired$error,
  unpaired_sig = res_unpaired$n_sig, unpaired_nonconv = res_unpaired$n_nonconv, 
  unpaired_ctrl_up = res_unpaired$n_ctrl_up, unpaired_ctrl_dn = res_unpaired$n_ctrl_dn, unpaired_time = res_unpaired$time, unpaired_err = res_unpaired$error
)

dir.create("results/subsampling", showWarnings = FALSE, recursive = TRUE)
saveRDS(out_df, sprintf("results/subsampling/n%03d_rep%02d.rds", n, rep))
message("Done n=", n, " rep=", rep)