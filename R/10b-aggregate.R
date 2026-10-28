# R/10b-aggregate.R
#
# Aggregate the 113 leave-one-out fits.
# Calculates the stability fraction for each gene in the primary list.
#
# Inputs:  results/loo/loo_*.rds (113 files)
# Outputs: results/analysis_4.rds

message("Loading leave-one-out fits...")
files <- list.files("results/loo", pattern = "^loo_TCGA.*\\.rds$", full.names = TRUE)
stopifnot(length(files) == 113)

# Load all LOO files
loo_list <- lapply(files, readRDS)
message("Loaded ", length(loo_list), " fits.")

# 1. Aggregate gene stability
# Extract the genes data.frames and combine them
genes_list <- lapply(loo_list, function(x) x$genes)
all_genes <- do.call(rbind, genes_list)

# Calculate stability fraction per gene
# (Number of times is_sig is TRUE divided by 113)
stability_df <- aggregate(is_sig ~ gene_id, data = all_genes, FUN = mean)
colnames(stability_df)[2] <- "stability_fraction"

# 2. Aggregate metadata (fit times, non-converged)
meta_list <- lapply(loo_list, function(x) x$meta)
meta_df <- do.call(rbind, meta_list)

# --- Sanity checks ---------------------------------------------------------
stopifnot(nrow(stability_df) == 20620) # All primary genes present
stopifnot(all(stability_df$stability_fraction >= 0 & stability_df$stability_fraction <= 1))
stopifnot(nrow(meta_df) == 113)

# --- Report ----------------------------------------------------------------
message("\n==================================================")
message("LEAVE-ONE-OUT SUMMARY (n=113 refits)")
message("  Median fit time   : ", round(median(meta_df$fit_time) / 60, 1), " min")
message("  Total compute time: ", round(sum(meta_df$fit_time) / 3600, 1), " hours")
message("  Non-converged genes (median per refit): ", median(meta_df$n_nonconv))
message("--------------------------------------------------")
message("Gene Stability:")
message("  Genes with 100% stability (113/113): ", sum(stability_df$stability_fraction == 1))
message("  Genes with >95% stability          : ", sum(stability_df$stability_fraction >= 0.95))
message("  Genes with <50% stability          : ", sum(stability_df$stability_fraction < 0.50))
message("  Median stability across all genes  : ", round(median(stability_df$stability_fraction), 3))
message("==================================================")

# --- Optional: Check the 15 controls ---------------------------------------
inputs <- readRDS("results/model-inputs.rds")
gene_info <- inputs$gene_info
stability_annot <- merge(stability_df, gene_info[, c("gene_id", "gene_name")], by = "gene_id")
controls <- c("ESR1", "GATA3", "FOXA1", "XBP1", "MKI67", "TOP2A", "AURKA", "BIRC5",
              "CCNB1", "CDC20", "ADIPOQ", "FABP4", "PLIN1", "LEP", "CFD")

ctrl_stab <- stability_annot[stability_annot$gene_name %in% controls, ]
message("\nControl stability:")
print(ctrl_stab[order(ctrl_stab$stability_fraction), c("gene_name", "stability_fraction")])

# --- Save ------------------------------------------------------------------
saveRDS(list(genes = stability_df, meta = meta_df), "results/analysis_4.rds")
message("\nSaved to results/analysis_4.rds")