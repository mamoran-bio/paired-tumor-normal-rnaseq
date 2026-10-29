# R/13-filter-sensitivity.R
#
# Fit the paired model on the expanded universe (35,584 genes) using the
# `design =` filter, and compare it with the primary analysis (A).
#
# Inputs:  data via R/02, R/03, R/04; results/res_primary.rds; results/analysis_4.rds
# Outputs: results/analysis_7.rds; docs/results/filter-sensitivity.csv

suppressPackageStartupMessages(library(DESeq2))
source("R/00-helpers.R")

# --- Configuration for smoke test ------------------------------------------
n_genes_smoke <- Sys.getenv("R13_N_GENES", unset = NA)
is_smoke      <- !is.na(n_genes_smoke)
out_dir       <- if (is_smoke) "results/smoke" else "results"

# Ensure directory exists and is writable
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(out_dir), file.access(out_dir, mode = 2) == 0)
if (!is_smoke) stopifnot(file.access("docs/results", mode = 2) == 0)

t0 <- Sys.time()
message("Starting Analysis 7 at ", format(t0, "%Y-%m-%d %H:%M:%S"))

# --- 1. Load inputs (like R/05) --------------------------------------------
message("Generating inputs for Analysis 7 (design = filter)...")
source("R/02-load-counts.R")
source("R/03-build-coldata.R")
source("R/04-filter-genes.R")

# Extract the expanded universe using keep_design from R/04
counts_filtered_B <- counts_raw[keep_design, ]

# --- 2. Pre-flight checks (run in smoke test too) --------------------------
res_df_A <- readRDS("results/res_primary.rds")

stopifnot(
  nrow(counts_filtered_B) == 35584,
  ncol(counts_filtered_B) == 226,
  identical(colnames(counts_filtered_B), rownames(coldata)), # Order check
  nrow(res_df_A) == 25836,
  all(rownames(res_df_A) %in% rownames(counts_filtered_B)),
  "beta_converged" %in% names(res_df_A)
)

unique_universe_B <- setdiff(rownames(counts_filtered_B), rownames(res_df_A))
stopifnot(length(unique_universe_B) == 9748)
stopifnot(!anyNA(res_df_A$beta_converged))

# --- 3. Build gene_info for B ----------------------------------------------
rd <- SummarizedExperiment::rowData(rse)
gene_info_B <- as.data.frame(rd[rownames(counts_filtered_B),
                                c("gene_id", "gene_name", "gene_type")])
stopifnot(identical(gene_info_B$gene_id, rownames(counts_filtered_B)))

# Free memory
rm(rse, counts_reads, counts_aliquot, counts_raw)
invisible(gc())

# --- 4. Smoke test subsetting ----------------------------------------------
if (is_smoke) {
  n_genes <- as.integer(n_genes_smoke)
  message("Smoke test: subsetting B to ~", n_genes, " genes (plus 15 controls)")

  controls_all <- c(controls_up, controls_dn)
  control_genes <- gene_info_B$gene_id[gene_info_B$gene_name %in% controls_all]

  set.seed(20260918)
  sampled_genes <- sort(sample(rownames(counts_filtered_B), n_genes))
  genes_to_keep <- unique(c(sampled_genes, control_genes))

  counts_filtered_B <- counts_filtered_B[genes_to_keep, , drop = FALSE]
  gene_info_B <- gene_info_B[genes_to_keep, , drop = FALSE]
}

# --- 5. Fit model and extract results (like R/05, R/06) --------------------
message("Fitting paired model on expanded universe...")
dds_B <- DESeqDataSetFromMatrix(
  countData = counts_filtered_B,
  colData   = coldata,
  design    = ~ patient + condition
)
dds_B <- DESeq(dds_B, minReplicatesForReplace = Inf)

message("Extracting results...")
res_B <- results(
  dds_B,
  name                 = coef_name,
  independentFiltering = FALSE,
  cooksCutoff          = FALSE
)

# --- Checkpoint & Annotation -----------------------------------------------
stopifnot(identical(rownames(res_B), gene_info_B$gene_id))

res_df_B <- as.data.frame(res_B)
res_df_B$gene_id        <- rownames(res_df_B)
res_df_B$beta_converged <- mcols(dds_B)$betaConv
res_df_B$gene_name      <- gene_info_B$gene_name
res_df_B$gene_type      <- gene_info_B$gene_type

# gene_label: same code as R/06/R/07
gene_id_base <- sub("\\..*", "", res_df_B$gene_id)
dup_names <- names(which(table(res_df_B$gene_name) > 1))
res_df_B$gene_label <- ifelse(
  res_df_B$gene_name %in% dup_names,
  paste0(res_df_B$gene_name, ":", gene_id_base),
  res_df_B$gene_name
)

ckpt <- file.path(out_dir, "analysis_7-checkpoint.rds")
saveRDS(res_df_B, ckpt)

stopifnot(!anyNA(res_df_B$beta_converged),
          identical(names(res_df_B), names(res_df_A)))

# --- 6. Calculate metrics --------------------------------------------------
sig_A <- rownames(res_df_A)[!is.na(res_df_A$padj) & res_df_A$padj < 0.05]
sig_B <- rownames(res_df_B)[!is.na(res_df_B$padj) & res_df_B$padj < 0.05]

stopifnot(length(sig_A) == nrow(readRDS("results/analysis_4.rds")$genes))

shared   <- intersect(sig_A, sig_B)
A_only   <- setdiff(sig_A, sig_B)
B_only   <- setdiff(sig_B, sig_A)
jaccard  <- length(shared) / length(union(sig_A, sig_B))

sig_B_among_unique <- sum(unique_universe_B %in% sig_B)

nonconv_A <- sum(!res_df_A$beta_converged)
nonconv_B <- sum(!res_df_B$beta_converged)

# --- Controls --------------------------------------------------------------
message("--- Controls in A ---")
ctrl_up_A <- check_controls(res_df_A, controls_up, +1, pcol = "padj", lcol = "log2FoldChange", alpha = 0.05)
ctrl_dn_A <- check_controls(res_df_A, controls_dn, -1, pcol = "padj", lcol = "log2FoldChange", alpha = 0.05)

message("--- Controls in B ---")
ctrl_up_B <- check_controls(res_df_B, controls_up, +1, pcol = "padj", lcol = "log2FoldChange", alpha = 0.05)
ctrl_dn_B <- check_controls(res_df_B, controls_dn, -1, pcol = "padj", lcol = "log2FoldChange", alpha = 0.05)

# --- 7. Build output table (14 rows) ---------------------------------------
metrics <- c(
  genes_in_universe_A             = nrow(res_df_A),
  genes_in_universe_B             = nrow(res_df_B),
  genes_unique_to_design_filter   = length(unique_universe_B),
  significant_A                   = length(sig_A),
  significant_B                   = length(sig_B),
  significant_in_both             = length(shared),
  significant_in_A_only           = length(A_only),
  significant_in_B_only           = length(B_only),
  jaccard                         = jaccard,
  significant_B_among_unique_genes = sig_B_among_unique,
  nonconverged_A                  = nonconv_A,
  nonconverged_B                  = nonconv_B,
  controls_passing_A              = ctrl_up_A + ctrl_dn_A,
  controls_passing_B              = ctrl_up_B + ctrl_dn_B
)

out_df <- data.frame(metric = names(metrics), value = unname(metrics))
stopifnot(nrow(out_df) == 14)

# --- 8. Report and save ----------------------------------------------------
message("\n--- Analysis 7 Summary ---")
print(out_df)

saveRDS(list(results = res_df_B, metrics = metrics), file.path(out_dir, "analysis_7.rds"))

if (!is_smoke) {
  write.csv(out_df, "docs/results/filter-sensitivity.csv", row.names = FALSE)
  message("Saved to results/analysis_7.rds and docs/results/filter-sensitivity.csv")
} else {
  message("Smoke test complete. Results saved to ", out_dir)
}

# Cleanup checkpoint
if (file.exists(ckpt)) invisible(file.remove(ckpt))

t1 <- Sys.time()
message("Finished Analysis 7 at ", format(t1, "%Y-%m-%d %H:%M:%S"))
message("Total elapsed time: ", round(difftime(t1, t0, units = "mins"), 1), " minutes")