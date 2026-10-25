# R/08-unpaired-comparison.R
#
# Analysis 2 (paired vs unpaired) and Analysis 3 (dispersion comparison).
# Evaluates primary (FDR < 0.05) and secondary (lfcThreshold = 1, FDR < 0.01).

suppressPackageStartupMessages(library(DESeq2))
source("R/00-helpers.R")

# 1. Load data
message("Loading data...")
dds_paired   <- readRDS("results/dds_paired.rds")
dds_unpaired <- readRDS("results/dds_unpaired.rds")
inputs       <- readRDS("results/model-inputs.rds")
coef_name    <- inputs$coef_name
gene_info    <- inputs$gene_info

# Load paired results from R/06 and R/07
res_p_pri <- readRDS("results/res_primary.rds")
res_p_sec <- readRDS("results/res_secondary.rds")

# 2. Extract unpaired results 
message("Extracting unpaired results...")

res_u_pri_raw <- results(dds_unpaired, name = coef_name,
                         independentFiltering = FALSE, cooksCutoff = FALSE)
res_u_sec_raw <- results(dds_unpaired, name = coef_name, lfcThreshold = 1, 
                         independentFiltering = FALSE, cooksCutoff = FALSE)

# --- Format Primary Unpaired ---
res_u_pri <- as.data.frame(res_u_pri_raw)
res_u_pri$gene_id        <- rownames(res_u_pri)
res_u_pri$gene_name      <- gene_info$gene_name
res_u_pri$gene_type      <- gene_info$gene_type
res_u_pri$beta_converged <- mcols(dds_unpaired)$betaConv

# --- Format Secondary Unpaired (Decision b: rename for consistency) ---
res_u_sec <- as.data.frame(res_u_sec_raw)
res_u_sec$gene_id        <- rownames(res_u_sec)
res_u_sec$gene_name      <- gene_info$gene_name
res_u_sec$gene_type      <- gene_info$gene_type
res_u_sec$beta_converged <- mcols(dds_unpaired)$betaConv

# Rename columns to match res_p_sec
names(res_u_sec)[names(res_u_sec) == "log2FoldChange"] <- "log2FC_raw"
names(res_u_sec)[names(res_u_sec) == "lfcSE"]          <- "lfcSE_raw"
names(res_u_sec)[names(res_u_sec) == "stat"]           <- "stat_thresh"
names(res_u_sec)[names(res_u_sec) == "padj"]           <- "padj_thresh"

# Gene_label logic (same as R/06 and R/07)
add_gene_label <- function(df) {
  gene_id_base <- sub("\\..*", "", df$gene_id)
  dup_names <- names(which(table(df$gene_name) > 1))
  df$gene_label <- ifelse(
    df$gene_name %in% dup_names,
    paste0(df$gene_name, ":", gene_id_base),
    df$gene_name
  )
  df
}
res_u_pri <- add_gene_label(res_u_pri)
res_u_sec <- add_gene_label(res_u_sec)

saveRDS(res_u_pri, "results/res_unpaired_primary.rds")
saveRDS(res_u_sec, "results/res_unpaired_secondary.rds")

# 3. Function for Analysis 2 metrics (FIXED: using gene_id column)
calc_metrics <- function(df_p, df_u, p_col_p, p_col_u, lfc_col_p, lfc_col_u, alpha) {
  sig_p <- df_p$gene_id[!is.na(df_p[[p_col_p]]) & df_p[[p_col_p]] < alpha]
  sig_u <- df_u$gene_id[!is.na(df_u[[p_col_u]]) & df_u[[p_col_u]] < alpha]
  
  shared <- intersect(sig_p, sig_u)
  unique_p <- setdiff(sig_p, sig_u)
  unique_u <- setdiff(sig_u, sig_p)
  jaccard <- length(shared) / length(union(sig_p, sig_u))
  
  sign_conc <- NA
  if (length(shared) > 0) {
    lfc_p <- df_p[[lfc_col_p]][match(shared, df_p$gene_id)]
    lfc_u <- df_u[[lfc_col_u]][match(shared, df_u$gene_id)]
    sign_conc <- mean(sign(lfc_p) == sign(lfc_u), na.rm = TRUE)
  }
  
  list(
    n_paired = length(sig_p), n_unpaired = length(sig_u),
    jaccard = jaccard, n_shared = length(shared),
    n_unique_p = length(unique_p), n_unique_u = length(unique_u),
    sign_concordance = sign_conc,
    genes_unique_p = unique_p, genes_unique_u = unique_u
  )
}

# Run metrics
m2_primary   <- calc_metrics(res_p_pri, res_u_pri, "padj", "padj", 
                             "log2FoldChange", "log2FoldChange", 0.05)
m2_secondary <- calc_metrics(res_p_sec, res_u_sec, "padj_thresh", "padj_thresh", 
                             "log2FC_raw", "log2FC_raw", 0.01)

# 4. Sensitivity (Exclude flagged genes)
flag_genes <- res_p_pri$gene_id[!res_p_pri$beta_converged]
message("Number of flagged genes to exclude in sensitivity: ", length(flag_genes))

df_p_pri_sens <- res_p_pri[!res_p_pri$gene_id %in% flag_genes, ]
df_u_pri_sens <- res_u_pri[!res_u_pri$gene_id %in% flag_genes, ]
df_p_sec_sens <- res_p_sec[!res_p_sec$gene_id %in% flag_genes, ]
df_u_sec_sens <- res_u_sec[!res_u_sec$gene_id %in% flag_genes, ]

m2_primary_sens   <- calc_metrics(df_p_pri_sens, df_u_pri_sens, "padj", "padj", 
                                  "log2FoldChange", "log2FoldChange", 0.05)
m2_secondary_sens <- calc_metrics(df_p_sec_sens, df_u_sec_sens, "padj_thresh", "padj_thresh", 
                                  "log2FC_raw", "log2FC_raw", 0.01)

# 5. Analysis 3: Dispersion
disp_p <- dispersions(dds_paired)
disp_u <- dispersions(dds_unpaired)
m3_disp <- list(
  paired_median   = median(disp_p, na.rm = TRUE),
  unpaired_median = median(disp_u, na.rm = TRUE)
)

# 6. Unpaired control verification
message("\n--- Controls in unpaired model (primary) ---")
check_controls(res_u_pri, controls_up, +1, alpha = 0.05)
check_controls(res_u_pri, controls_dn, -1, alpha = 0.05)

message("\n--- Controls in unpaired model (secondary) ---")
check_controls(res_u_sec, controls_up, +1, pcol = "padj_thresh", lcol = "log2FC_raw", alpha = 0.01)
check_controls(res_u_sec, controls_dn, -1, pcol = "padj_thresh", lcol = "log2FC_raw", alpha = 0.01)

# ----- Save Results ------------------------------------------------------------------
analysis_2 <- list(
  primary               = m2_primary,
  secondary             = m2_secondary,
  primary_sensitivity   = m2_primary_sens,
  secondary_sensitivity = m2_secondary_sens
)

saveRDS(analysis_2, "results/analysis_2.rds")
saveRDS(m3_disp,   "results/analysis_3.rds")

# ----- Log Summary -------------------------------------------------------------------
message("\n==================================================")
message("ANALYSIS 2 SUMMARY")
message(sprintf("%-20s | %10s | %10s", "Metric", "Primary", "Secondary"))
message(sprintf("%-20s | %10d | %10d", "Paired genes", m2_primary$n_paired, m2_secondary$n_paired))
message(sprintf("%-20s | %10d | %10d", "Unpaired genes", m2_primary$n_unpaired, m2_secondary$n_unpaired))
message(sprintf("%-20s | %10.3f | %10.3f", "Jaccard index", m2_primary$jaccard, m2_secondary$jaccard))
message(sprintf("%-20s | %10.3f | %10.3f", "Sign concordance", m2_primary$sign_concordance, m2_secondary$sign_concordance))
message(sprintf("%-20s | %10d | %10d", "Unique to paired", m2_primary$n_unique_p, m2_secondary$n_unique_p))
message(sprintf("%-20s | %10d | %10d", "Unique to unpaired", m2_primary$n_unique_u, m2_secondary$n_unique_u))
message("--------------------------------------------------")
message("SENSITIVITY (excluding ", length(flag_genes), " flagged genes)")
message(sprintf("%-20s | %10.3f | %10.3f", "Jaccard index", m2_primary_sens$jaccard, m2_secondary_sens$jaccard))
message(sprintf("%-20s | %10.3f | %10.3f", "Sign concordance", m2_primary_sens$sign_concordance, m2_secondary_sens$sign_concordance))
message("==================================================")
message("ANALYSIS 3 SUMMARY (Dispersion)")
message("  Paired median:    ", round(m3_disp$paired_median, 4))
message("  Unpaired median:  ", round(m3_disp$unpaired_median, 4))
message("==================================================")

message("R/08 complete.")