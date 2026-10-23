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
    
    # AQUÍ ESTÁ EL CAMBIO: Añadimos lcol y pcol al mensaje
    message(sprintf("Control %-6s | %s %6.2f | %s %9.2e | sig %s | dir %s",
                    g, lcol, lfc, pcol, pval,
                    ifelse(sig, "yes", "no "),
                    ifelse(dir_ok, "ok", "WRONG")))
    
    if (sig && dir_ok) passed <- passed + 1
  }
  message(sprintf("Passed: %d / %d", passed, length(genes)))
  invisible(passed)
}
