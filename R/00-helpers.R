# R/00-helpers.R
#
# Shared helper functions and pre-registered breast controls.

controls_up <- c("ESR1", "GATA3", "FOXA1", "XBP1",
                 "MKI67", "TOP2A", "AURKA", "BIRC5",
                 "CCNB1", "CDC20")
controls_dn <- c("ADIPOQ", "FABP4", "PLIN1", "LEP", "CFD")

check_controls <- function(res_df, genes, expected_sign,
                           pcol = "padj", lcol = "log2FoldChange",
                           alpha = 0.05) {
  passed <- 0
  for (g in genes) {
    idx <- which(res_df$gene_name == g)
    if (length(idx) == 0) {
      message(sprintf("Control %-6s | NOT in gene universe", g))
      next
    }
    row <- res_df[idx[1], ]
    pval <- row[[pcol]]
    lfc  <- row[[lcol]]
    sig  <- !is.na(pval) && pval < alpha
    dir_ok <- !is.na(lfc) && sign(lfc) == expected_sign

    message(sprintf("Control %-6s | %s %6.2f | %s %9.2e | sig %s | dir %s",
                    g, lcol, lfc, pcol, pval,
                    ifelse(sig, "yes", "no "),
                    ifelse(dir_ok, "ok", "WRONG")))

    if (sig && dir_ok) passed <- passed + 1
  }
  message(sprintf("Passed: %d / %d", passed, length(genes)))
  invisible(passed)
}
