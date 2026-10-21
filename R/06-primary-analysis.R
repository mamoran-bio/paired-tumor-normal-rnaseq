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
controls_up <- c("ESR1", "GATA3", "FOXA1", "XBP1",
                 "MKI67", "TOP2A", "AURKA", "BIRC5",
                 "CCNB1", "CDC20")
controls_dn <- c("ADIPOQ", "FABP4", "PLIN1", "LEP", "CFD")

check_controls <- function(genes, expected_sign) {
  found <- 0
  for (g in genes) {
    idx <- which(res_df$gene_name == g)
    if (length(idx) == 0) {
      message(sprintf("Control %-6s | NOT in gene universe", g))
      next
    }
    row <- res_df[idx[1], ]
    sig <- !is.na(row$padj) && row$padj < 0.05
    sign_ok <- sign(row$log2FoldChange) == expected_sign
    message(sprintf("Control %-6s | log2FC %6.2f | padj %9.2e | sig %s | dir %s",
                    g, row$log2FoldChange, row$padj,
                    ifelse(sig, "yes", "no "),
                    ifelse(sign_ok, "ok", "WRONG")))
    if (sig && sign_ok) found <- found + 1
  }
  found
}

n_up <- check_controls(controls_up, +1)
n_dn <- check_controls(controls_dn, -1)
message(sprintf("Controls up:   %d / %d significant with correct direction", n_up, length(controls_up)))
message(sprintf("Controls down: %d / %d significant with correct direction", n_dn, length(controls_dn)))

# --- Save ------------------------------------------------------------------
saveRDS(res_df, "results/res_primary.rds")
message("Done. Outputs written to results/res_primary.rds")
