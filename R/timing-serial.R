# R/timing-serial.R
#
# Time the three phases of a serial DESeq() fit on a random sample of
# genes from the filtered universe. Used for configurations A (BLAS
# default) and B (OPENBLAS_NUM_THREADS=1); the configuration is set at
# container start, not inside this script.
#
# Not part of the numbered pipeline. Run via:
#   docker run -it --rm \
#     -e OPENBLAS_NUM_THREADS=1 -e TIMING_N_GENES=200 \
#     -v "$PWD":/home/rstudio/project paired-rnaseq:dev \
#     Rscript R/timing-serial.R
#
# Sources: R/02, R/03, R/04
# Uses:    counts_filtered, coldata
# Outputs: printed timing and memory report; sessionInfo()

source("R/02-load-counts.R")
source("R/03-build-coldata.R")
source("R/04-filter-genes.R")

library(DESeq2)

# --- Configuration ---------------------------------------------------------
n_genes      <- as.integer(Sys.getenv("TIMING_N_GENES", unset = "5000"))
blas_threads <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = "default")

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
t_size_factors <- system.time({
  dds <- estimateSizeFactors(dds)
})

t_disp <- system.time({
  dds <- estimateDispersions(dds)
})

t_wald <- system.time({
  dds <- nbinomWaldTest(dds)
})
# --- Peak memory since reset -----------------------------------------------
peak_mb <- sum(gc()[, 6])

# --- Report ----------------------------------------------------------------
report_phase <- function(label, t) {
  ratio <- (t[["user.self"]] + t[["sys.self"]]) / t[["elapsed"]]
  message(sprintf(
    "%-13s elapsed %8.1f s | user %8.1f | system %8.1f | ratio %.2f",
    label, t[["elapsed"]], t[["user.self"]], t[["sys.self"]], ratio))
}

total_elapsed <- t_size_factors[["elapsed"]] +
                 t_disp[["elapsed"]] +
                 t_wald[["elapsed"]]
t_total <- t_size_factors + t_disp + t_wald

message("")
message("Timing report")
message("  genes         : ", n_genes)
message("  samples       : ", ncol(counts_timing))
message("  BLAS threads  : ", blas_threads)
message("  platform      : ", R.version$platform)
message("  DESeq2        : ", as.character(packageVersion("DESeq2")))
message("")
report_phase("sizeFactors", t_size_factors)
report_phase("dispersions",  t_disp)
report_phase("Wald",         t_wald)
report_phase("total", t_total)
message(sprintf("%-13s elapsed %8.1f s", "total", total_elapsed))
message(sprintf("  seconds/gene  : %.4f", total_elapsed / n_genes))
message(sprintf("  peak memory   : %.0f MB", peak_mb))

# --- Cleanup and session info ----------------------------------------------
rm(dds, counts_timing)
invisible(gc())

message("")
print(sessionInfo())
