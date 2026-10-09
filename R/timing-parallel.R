# R/timing-parallel.R
#
# Time a DESeq() fit with parallel = TRUE on a random sample of genes
# from the filtered universe. Configuration C: 1 BLAS thread, k workers.
#
# Not part of the numbered pipeline. Run via:
#   docker run -it --rm \
#     -e OPENBLAS_NUM_THREADS=1 \
#     -e TIMING_N_GENES=200 \
#     -e TIMING_WORKERS=4 \
#     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
#     Rscript R/timing-parallel.R
#
# Sources: R/02, R/03, R/04
# Uses:    counts_filtered, coldata
# Outputs: printed timing and memory report; sessionInfo()

source("R/02-load-counts.R")
source("R/03-build-coldata.R")
source("R/04-filter-genes.R")

library(DESeq2)
library(BiocParallel)

# --- Configuration ---------------------------------------------------------
n_genes      <- as.integer(Sys.getenv("TIMING_N_GENES", unset = "5000"))
n_workers    <- as.integer(Sys.getenv("TIMING_WORKERS", unset = "4"))
blas_threads <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = "default")

param <- MulticoreParam(workers = n_workers, RNGseed = 20260918)

# --- Measure and release R/02 intermediates --------------------------------
# These sizes are the measurement missing from the R/02 memory decision.
message("R/02 intermediates:")
message("  rse            : ",
        round(as.numeric(object.size(rse)) / 1e6, 1), " MB")
message("  counts_reads   : ",
        round(as.numeric(object.size(counts_reads)) / 1e6, 1), " MB")
message("  counts_aliquot : ",
        round(as.numeric(object.size(counts_aliquot)) / 1e6, 1), " MB")
rm(rse, counts_reads, counts_aliquot)
invisible(gc())

# --- Sample genes ----------------------------------------------------------
set.seed(20260918)
sampled_genes <- sort(sample(rownames(counts_filtered), n_genes))
counts_timing <- counts_filtered[sampled_genes, , drop = FALSE]

# --- Build the DESeqDataSet ------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_timing,
  colData   = coldata,
  design    = ~ patient + condition
)

# --- Reset memory counter --------------------------------------------------
invisible(gc(reset = TRUE))

# --- Time each phase -------------------------------------------------------
message("fit start  UTC ", format(Sys.time(), "%H:%M:%S", tz = "UTC"))

t_fit <- system.time({
  dds <- DESeq(dds, parallel = TRUE, BPPARAM = param)
})

message("fit end    UTC ", format(Sys.time(), "%H:%M:%S", tz = "UTC"))

# --- Peak memory since reset -----------------------------------------------
peak_mb_main <- sum(gc()[, 6])

# --- Report ----------------------------------------------------------------
total_elapsed <- t_fit[["elapsed"]]

cpu_total <- t_fit[["user.self"]] + t_fit[["sys.self"]] +
             t_fit[["user.child"]] + t_fit[["sys.child"]]
ratio <- cpu_total / t_fit[["elapsed"]]

message("")
message("Timing report — configuration C (parallel = TRUE)")
message("  genes         : ", n_genes)
message("  samples       : ", ncol(counts_timing))
message("  workers       : ", n_workers)
message("  BLAS threads  : ", blas_threads)
message("  platform      : ", R.version$platform)
message("  DESeq2        : ", as.character(packageVersion("DESeq2")))
message("")
message(sprintf(
  "%-13s elapsed %8.1f s | user.self %8.1f | sys.self %8.1f | user.child %8.1f | sys.child %8.1f | ratio %.2f",
  "DESeq()",
  t_fit[["elapsed"]],
  t_fit[["user.self"]],
  t_fit[["sys.self"]],
  t_fit[["user.child"]],
  t_fit[["sys.child"]],
  ratio))
message(sprintf("  seconds/gene  : %.4f", total_elapsed / n_genes))
message(sprintf("  peak memory   : %.0f MB (main process only; workers not included)",
                peak_mb_main))

# --- Cleanup and session info ----------------------------------------------
rm(dds, counts_timing)
invisible(gc())

message("")
print(sessionInfo())
