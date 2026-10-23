# R/07-secondary-analysis.R
#
# Thresholded test (lfcThreshold = 1) + shrunken log2FC on the paired fit.
# Produces one combined table: significance from the thresholded padj,
# shrunken log2FC reported alongside for ranking and plots. No gene is
# selected on the basis of a shrunken estimate (plan line 332).
#
# Inputs:  results/dds_paired.rds, results/model-inputs.rds
# Outputs: results/res_secondary.rds

suppressPackageStartupMessages(library(DESeq2))
source("R/00-helpers.R")

# --- Load inputs -----------------------------------------------------------
message("Loading inputs...")
dds_paired <- readRDS("results/dds_paired.rds")
inputs     <- readRDS("results/model-inputs.rds")
gene_info  <- inputs$gene_info
coef_name  <- inputs$coef_name

stopifnot(identical(rownames(dds_paired), gene_info$gene_id))

# --- 1. Extract models -----------------------------------------------------
message("Running thresholded test (lfcThreshold = 1)...")
res_thresh <- results(
  dds_paired,
  name                 = coef_name,
  lfcThreshold         = 1,
  independentFiltering = FALSE,
  cooksCutoff          = FALSE
)

message("Calculating shrunken log2FC (apeglm)...")
res_shrunk <- lfcShrink(
  dds_paired,
  coef = coef_name, # lfcShrink usa coef, no name
  type = "apeglm"
)

stopifnot(identical(rownames(res_thresh), rownames(res_shrunk)))
stopifnot(identical(rownames(res_thresh), gene_info$gene_id))

# --- 2. Build combined data.frame ------------------------------------------
res_df <- data.frame(
  gene_id        = rownames(res_thresh),
  gene_name      = gene_info$gene_name,
  gene_type      = gene_info$gene_type,
  beta_converged = mcols(dds_paired)$betaConv,
  log2FC_raw     = res_thresh$log2FoldChange,
  lfcSE_raw      = res_thresh$lfcSE,
  stat_thresh    = res_thresh$stat,
  padj_thresh    = res_thresh$padj,
  log2FC_shrunk  = res_shrunk$log2FoldChange,
  lfcSE_shrunk   = res_shrunk$lfcSE,
  stringsAsFactors = FALSE
)

# gene_label (formato idéntico a R/06: gene_name:gene_id_base)
gene_id_base <- sub("\\..*", "", res_df$gene_id)
dup_names <- names(which(table(res_df$gene_name) > 1))
res_df$gene_label <- ifelse(
  res_df$gene_name %in% dup_names,
  paste0(res_df$gene_name, ":", gene_id_base),
  res_df$gene_name
)

# --- 3. Reports ------------------------------------------------------------
message("Non-converged genes: ", sum(!res_df$beta_converged), " of ", nrow(res_df))
message("Secondary list (padj_thresh < 0.01): ",
        sum(res_df$padj_thresh < 0.01, na.rm = TRUE), " genes")

message("--- Positive controls at threshold ---")
check_controls(res_df, controls_up, +1, pcol = "padj_thresh", lcol = "log2FC_raw", alpha = 0.01)

message("--- Negative controls at threshold ---")
check_controls(res_df, controls_dn, -1, pcol = "padj_thresh", lcol = "log2FC_raw", alpha = 0.01)

# --- 4. Save output --------------------------------------------------------
saveRDS(res_df, "results/res_secondary.rds")
message("Done. Outputs written to results/res_secondary.rds")
