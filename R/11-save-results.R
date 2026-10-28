# R/11-save-results.R
#
# Export the final analysis tables to docs/results/ as CSVs and capture
# the R environment used for the analysis. Writes the five CSVs listed
# in docs/analysis-plan.md lines 455-459, plus the sessionInfo of the
# container in which R/11 ran (line 451).

message("Loading packages to record their versions in sessionInfo...")
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(edgeR)
  library(recount3)
})

message("Saving final results to docs/results/...")
dir.create("docs/results", showWarnings = FALSE, recursive = TRUE)

# Load analysis_2 first to get the expected counts for the sanity checks
analysis_2 <- readRDS("results/analysis_2.rds")

# --- 1. Primary (plan line 455) --------------------------------------------
res_primary <- readRDS("results/res_primary.rds")
primary_list <- res_primary[!is.na(res_primary$padj) & res_primary$padj < 0.05, ]
stopifnot(nrow(primary_list) == analysis_2$primary$n_paired)

write.csv(primary_list,
          "docs/results/primary-paired.csv",
          row.names = FALSE)

# --- 2. Secondary high-confidence (plan line 456) --------------------------
res_secondary <- readRDS("results/res_secondary.rds")
secondary_list <- res_secondary[!is.na(res_secondary$padj_thresh) &
                                res_secondary$padj_thresh < 0.01, ]
stopifnot(nrow(secondary_list) == analysis_2$secondary$n_paired)

write.csv(secondary_list,
          "docs/results/secondary-paired-highconfidence.csv",
          row.names = FALSE)

# --- 3. Unpaired comparison (plan line 457) --------------------------------
# Summary of Analysis 2 metrics (paired vs unpaired) for primary and
# secondary, plus the pre-registered sensitivity (2026-10-21 amendment:
# exclude the 13 non-converged paired genes, use original padj).

row_names <- c("Paired_genes", "Unpaired_genes",
               "Shared_genes",
               "Unique_to_paired", "Unique_to_unpaired",
               "Jaccard", "Sign_concordance")

col_from_list <- function(l) {
  c(l$n_paired, l$n_unpaired, l$n_shared,
    l$n_unique_p, l$n_unique_u,
    l$jaccard, l$sign_concordance)
}

unpaired_comparison <- data.frame(
  metric                = row_names,
  primary               = col_from_list(analysis_2$primary),
  secondary             = col_from_list(analysis_2$secondary),
  primary_sensitivity   = col_from_list(analysis_2$primary_sensitivity),
  secondary_sensitivity = col_from_list(analysis_2$secondary_sensitivity)
)
write.csv(unpaired_comparison,
          "docs/results/unpaired-comparison.csv",
          row.names = FALSE)

# --- 4. Subsampling curve (plan line 458) ----------------------------------
analysis_5 <- readRDS("results/analysis_5.rds")
write.csv(analysis_5$summary,
          "docs/results/subsampling.csv",
          row.names = FALSE)

# --- 5. Leave-one-out stability (plan line 459) ----------------------------
# Join stability_fraction with gene_name for readability.
# Order by stability (least stable first) so the fragile genes are visible.
analysis_4 <- readRDS("results/analysis_4.rds")
inputs    <- readRDS("results/model-inputs.rds")
gene_info <- inputs$gene_info

loo <- merge(analysis_4$genes,
             gene_info[, c("gene_id", "gene_name")],
             by = "gene_id", all.x = TRUE, sort = FALSE)
loo <- loo[, c("gene_id", "gene_name", "stability_fraction")]
loo <- loo[order(loo$stability_fraction, loo$gene_id), ]
write.csv(loo,
          "docs/results/leave-one-out.csv",
          row.names = FALSE)

# --- 6. Session info (plan line 451) ---------------------------------------
message("Saving session info to docs/sessionInfo-analysis.txt...")
writeLines(capture.output(sessionInfo()),
           "docs/sessionInfo-analysis.txt")

message("All CSVs and sessionInfo saved.")