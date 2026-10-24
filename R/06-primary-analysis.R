# R/06-primary-analysis.R
#
# Extract the primary differential expression results from the paired
# model. Merge with gene annotation, flag non-converged genes, and 
# evaluate pre-registered controls.
#
# Inputs:  results/dds_paired.rds
#          results/model-inputs.rds
# Outputs: results/res_primary.rds

library(DESeq2)
source("R/00-helpers.R")

# --- Load inputs -----------------------------------------------------------
message("Loading models and annotations...")
dds_paired <- readRDS("results/dds_paired.rds")
inputs     <- readRDS("results/model-inputs.rds")

gene_info <- inputs$gene_info
coef_name <- inputs$coef_name

# Checkpoint: ensure the model and the annotation are aligned
stopifnot(identical(rownames(dds_paired), gene_info$gene_id))

# --- Extract results -------------------------------------------------------
message("Extracting results...")
res_primary <- results(
  dds_paired,
  name                 = coef_name,   # 'name' is conventional for string inputs
  independentFiltering = FALSE,
  cooksCutoff          = FALSE
)

# Convert to data.frame
res_df <- as.data.frame(res_primary)
res_df$gene_id <- rownames(res_df)

# --- Add convergence flag and count (Amendment 2026-10-21) -----------------
res_df$beta_converged <- mcols(dds_paired)$betaConv
message("Non-converged genes: ", sum(!res_df$beta_converged), " of ", nrow(res_df))

# --- Merge annotation using direct indexing (prevents reordering) ----------
stopifnot(identical(res_df$gene_id, gene_info$gene_id))
res_df$gene_name <- gene_info$gene_name
res_df$gene_type <- gene_info$gene_type

# --- Create gene_label (disambiguate duplicates) ---------------------------
gene_id_base <- sub("\\..*", "", res_df$gene_id)
dup_names <- names(which(table(res_df$gene_name) > 1))
is_dup <- res_df$gene_name %in% dup_names

res_df$gene_label <- ifelse(
  is_dup,
  paste0(res_df$gene_name, ":", gene_id_base),
  res_df$gene_name
)

# --- Summary and FDR sensitivity (Analysis 6) ------------------------------
n_sig_05 <- sum(res_df$padj < 0.05, na.rm = TRUE)
n_sig_01 <- sum(res_df$padj < 0.01, na.rm = TRUE)

message("Primary list (FDR < 0.05): ", n_sig_05, " genes")
message("Sensitivity (FDR < 0.01):  ", n_sig_01, " genes")

# --- Pre-registered controls (plan line 397) ------------------------------
message("--- Positive controls ---")
n_up <- check_controls(res_df, controls_up, +1, pcol = "padj", lcol = "log2FoldChange", alpha = 0.05)
message("--- Negative controls ---")
n_dn <- check_controls(res_df, controls_dn, -1, pcol = "padj", lcol = "log2FoldChange", alpha = 0.05)

# --- Save ------------------------------------------------------------------
saveRDS(res_df, "results/res_primary.rds")
message("Done. Outputs written to results/res_primary.rds")
