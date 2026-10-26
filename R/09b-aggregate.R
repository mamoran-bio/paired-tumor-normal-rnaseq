# R/09b-aggregate.R
#
# Aggregate the 251 subsampling fits into a single analysis_5 object.
# Calculates median and 10-90 % bands for paired and unpaired designs.
#
# Inputs:  results/subsampling/*.rds (251 files)
# Outputs: results/analysis_5.rds (list with raw data and summary table)

message("Loading subsampling fits...")
files <- list.files("results/subsampling", pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(files) == 251)

raw_df <- do.call(rbind, lapply(files, readRDS))
message("Loaded ", nrow(raw_df), " rows across ", length(unique(raw_df$n)), " values of n.")

# Aggregate per n (base R; no new dependencies)
by_n <- split(raw_df, raw_df$n)

summary_df <- do.call(rbind, lapply(by_n, function(x) {
  data.frame(
    n                     = x$n[1],
    reps                  = nrow(x),
    paired_med            = median(x$paired_sig, na.rm = TRUE),
    paired_p10            = quantile(x$paired_sig, 0.10, na.rm = TRUE, names = FALSE),
    paired_p90            = quantile(x$paired_sig, 0.90, na.rm = TRUE, names = FALSE),
    paired_nonconv_med    = median(x$paired_nonconv, na.rm = TRUE),
    paired_ctrl_up_med    = median(x$paired_ctrl_up, na.rm = TRUE),
    paired_ctrl_dn_med    = median(x$paired_ctrl_dn, na.rm = TRUE),
    unpaired_med          = median(x$unpaired_sig, na.rm = TRUE),
    unpaired_p10          = quantile(x$unpaired_sig, 0.10, na.rm = TRUE, names = FALSE),
    unpaired_p90          = quantile(x$unpaired_sig, 0.90, na.rm = TRUE, names = FALSE),
    unpaired_nonconv_med  = median(x$unpaired_nonconv, na.rm = TRUE),
    unpaired_ctrl_up_med  = median(x$unpaired_ctrl_up, na.rm = TRUE),
    unpaired_ctrl_dn_med  = median(x$unpaired_ctrl_dn, na.rm = TRUE)
  )
}))
summary_df$ratio_med <- summary_df$paired_med / summary_df$unpaired_med
summary_df <- summary_df[order(summary_df$n), ]
rownames(summary_df) <- NULL

# Sanity checks
stopifnot(all(summary_df$reps >= 1))
stopifnot(all(summary_df$paired_med >= 0, na.rm = TRUE))
stopifnot(sum(summary_df$reps) == 251)

# Report
message("\n==================================================")
message("SUBSAMPLING CURVE SUMMARY")
message(sprintf("%-4s | %-19s | %-19s | %-5s | %-4s",
                "n", "Paired [p10-p90]", "Unpaired [p10-p90]", "Ratio", "reps"))
message("------------------------------------------------------------")
for (i in seq_len(nrow(summary_df))) {
  s <- summary_df[i, ]
  # CORREGIDO AQUÍ: %5.0f en lugar de %5d
  message(sprintf("%-4d | %5.0f [%5.0f-%5.0f] | %5.0f [%5.0f-%5.0f] | %5.2fx | %4d",
                  s$n,
                  s$paired_med, s$paired_p10, s$paired_p90,
                  s$unpaired_med, s$unpaired_p10, s$unpaired_p90,
                  s$ratio_med, s$reps))
}
message("------------------------------------------------------------")
message("Controls found (median; paired up/dn | unpaired up/dn):")
for (i in seq_len(nrow(summary_df))) {
  s <- summary_df[i, ]
  # CORREGIDO AQUÍ: %.1f en lugar de %2d
  message(sprintf("  n=%-3d | Paired: %.1f/10 up, %.1f/5 dn | Unpaired: %.1f/10 up, %.1f/5 dn",
                  s$n, s$paired_ctrl_up_med, s$paired_ctrl_dn_med,
                  s$unpaired_ctrl_up_med, s$unpaired_ctrl_dn_med))
}
message("Non-converged (median paired): ",
        paste(summary_df$n, summary_df$paired_nonconv_med, sep="=", collapse=", "))
message("==================================================")

# Save
saveRDS(list(raw = raw_df, summary = summary_df), "results/analysis_5.rds")
message("Saved to results/analysis_5.rds")