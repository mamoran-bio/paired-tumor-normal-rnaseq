# R/05-fit-models.R
#
# Fit the paired and unpaired DESeq2 models on the filtered count matrix.
# Produces the two fitted objects that R/06, R/07 and R/08 consume, and
# the input list that R/09 and R/10 consume.
#
# Inputs:  counts_filtered, coldata (via R/02, R/03, R/04)
# Outputs: results/model-inputs.rds
#          results/dds_paired.rds
#          results/dds_unpaired.rds
#          results/sessionInfo.txt
#
# Run with a single BLAS thread (see 2026-10-15 amendment):
#   mkdir -p docs/logs
#   docker run --rm -e OPENBLAS_NUM_THREADS=1 \
#     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
#     Rscript R/05-fit-models.R 2>&1 | tee docs/logs/05-fit-models.log
#
# Smoke test (200 genes, outputs to results/smoke/):
#   docker run --rm -e OPENBLAS_NUM_THREADS=1 -e R05_N_GENES=200 \
#     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
#     Rscript R/05-fit-models.R

source("R/02-load-counts.R")
source("R/03-build-coldata.R")
source("R/04-filter-genes.R")

library(DESeq2)

# --- Configuration ---------------------------------------------------------
n_genes_smoke <- Sys.getenv("R05_N_GENES", unset = NA)
is_smoke      <- !is.na(n_genes_smoke)

out_dir <- if (is_smoke) "results/smoke" else "results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Pre-flight checks -----------------------------------------------------
# Stop before spending 45 minutes if the inputs are not what R/04 produced.
stopifnot(nrow(counts_filtered) == 25836)
stopifnot(ncol(counts_filtered) == 226)
stopifnot(identical(colnames(counts_filtered), rownames(coldata)))
stopifnot(nlevels(coldata$patient) == 113)
stopifnot(identical(levels(coldata$condition), c("normal", "tumor")))

# --- Subset for smoke test if requested ------------------------------------
if (is_smoke) {
  n_genes <- as.integer(n_genes_smoke)
  message("Smoke test: subsetting to ", n_genes, " genes")
  set.seed(20260918)
  genes_to_keep <- sort(sample(rownames(counts_filtered), n_genes))
  counts_filtered <- counts_filtered[genes_to_keep, , drop = FALSE]
}

# --- Save inputs and free memory -------------------------------------------
message("Saving model inputs to ", out_dir, "/model-inputs.rds")
saveRDS(
  list(counts_filtered = counts_filtered,
       coldata = coldata,
       coef_name = coef_name),
  file.path(out_dir, "model-inputs.rds")
)

# Free intermediates from R/02 that are not needed for the fits.
rm(rse, counts_reads, counts_aliquot, counts_raw)
invisible(gc())

# --- Paired fit ------------------------------------------------------------
dds_paired_file <- file.path(out_dir, "dds_paired.rds")

if (file.exists(dds_paired_file) && !is_smoke) {
  message("Loading existing paired fit from ", dds_paired_file)
  dds_paired <- readRDS(dds_paired_file)
  
  # A file left by a run with other inputs must not be reused silently.
  stopifnot(identical(rownames(dds_paired), rownames(counts_filtered)))
  stopifnot(identical(colnames(dds_paired), colnames(counts_filtered)))
  stopifnot(identical(deparse(design(dds_paired)), "~patient + condition"))
} else {
  message("Fitting paired model (~ patient + condition)...")
  dds_paired <- DESeqDataSetFromMatrix(
    countData = counts_filtered,
    colData   = coldata,
    design    = ~ patient + condition
  )
  
  t_paired <- system.time({
    dds_paired <- DESeq(dds_paired)
  })
  
  message("Paired fit elapsed: ", round(t_paired[["elapsed"]], 1), " s")
  
  saveRDS(dds_paired, dds_paired_file)
}

# Guard the coefficient name on both paths (load and fit).
stopifnot(coef_name %in% resultsNames(dds_paired))

# Free memory before the next fit
rm(dds_paired)
invisible(gc())

# --- Unpaired fit ----------------------------------------------------------
message("Fitting unpaired model (~ condition)...")
dds_unpaired <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = coldata,
  design    = ~ condition
)

t_unpaired <- system.time({
  dds_unpaired <- DESeq(dds_unpaired)
})

message("Unpaired fit elapsed: ", round(t_unpaired[["elapsed"]], 1), " s")
stopifnot(coef_name %in% resultsNames(dds_unpaired))

saveRDS(dds_unpaired, file.path(out_dir, "dds_unpaired.rds"))

# --- Session info ----------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

message("Done. Outputs written to ", out_dir)
