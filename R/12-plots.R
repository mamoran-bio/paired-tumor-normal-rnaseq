# R/12-plots.R  (styled version, for comparison with your current script)
#
# Generate publication figures from the stored result objects.
# Colours: Okabe-Ito, as fixed in the 2026-10-28 plan amendment.
#
# Inputs:  results/res_primary.rds
#          results/analysis_5.rds
#          results/analysis_4.rds
# Outputs: docs/figures/*.png

suppressPackageStartupMessages(library(ggplot2))

message("Generating figures...")
dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

# --- Style ------------------------------------------------------------------
# Data colours: Okabe-Ito (amendment).
col_sig      <- "#0072B2"  # blue: significant / primary / paired
col_ns       <- "#999999"  # neutral grey: not significant
col_paired   <- "#0072B2"
col_unpaired <- "#D55E00"  # vermillion: unpaired / comparison
# Neutral chrome (text, grid, axes). These are not data colours.
ink_1 <- "#1F2933"; ink_2 <- "#52606D"; ink_3 <- "#7B8794"
grid_col <- "#E4E7EB"; axis_col <- "#CBD2D9"

theme_portfolio <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text                  = element_text(colour = ink_1),
      plot.background       = element_rect(fill = "white", colour = NA),
      plot.title            = element_text(face = "bold", size = rel(1.3),
                                           margin = margin(b = 4)),
      plot.subtitle         = element_text(colour = ink_2, size = rel(0.95),
                                           margin = margin(b = 14)),
      plot.caption          = element_text(colour = ink_3, size = rel(0.8), hjust = 0,
                                           lineheight = 1.25, margin = margin(t = 14)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      plot.margin           = margin(18, 22, 14, 18),
      axis.title            = element_text(colour = ink_2, size = rel(0.9)),
      axis.title.x          = element_text(margin = margin(t = 10)),
      axis.title.y          = element_text(margin = margin(r = 10)),
      axis.text             = element_text(colour = ink_3, size = rel(0.85)),
      axis.line.x           = element_line(colour = axis_col, linewidth = 0.4),
      panel.grid.major      = element_line(colour = grid_col, linewidth = 0.3),
      panel.grid.minor      = element_blank(),
      legend.position       = "bottom",
      legend.title          = element_blank(),
      legend.text           = element_text(colour = ink_2, size = rel(0.9)),
      legend.margin         = margin(t = 4)
    )
}

fmt <- function(x) format(x, big.mark = ",", trim = TRUE)

compact <- function(x) {  # 0, 500, 1K, 2.5K, 20K, 1M
  vapply(x, function(v) {
    if (is.na(v)) NA_character_
    else if (abs(v) >= 1e6) paste0(v / 1e6, "M")
    else if (abs(v) >= 1e3) paste0(v / 1e3, "K")
    else format(v, trim = TRUE)
  }, character(1))
}

save_fig <- function(plot, file) {
  ggsave(file.path("docs/figures", file), plot,
         width = 8, height = 6, dpi = 300, bg = "white")
}

# --- 1. Primary results: data used by figures 1, 2 and 4 --------------------
res_primary <- readRDS("results/res_primary.rds")

n_dropped <- sum(is.na(res_primary$padj) | is.na(res_primary$log2FoldChange))
message("Dropped ", n_dropped, " genes with NA padj or log2FC from Volcano/MA.")

plot_df <- res_primary[!is.na(res_primary$padj) & !is.na(res_primary$log2FoldChange), ]
# padj == 0 means the p-value underflowed (below machine precision); -log10(0) is Inf.
# ggplot2 would keep those points and draw them at the panel edge, i.e. at an
# arbitrary height. We draw them at a stated ceiling instead and say so in the caption.
n_zero <- sum(plot_df$padj == 0)
y_cap  <- max(-log10(plot_df$padj[plot_df$padj > 0])) * 1.05   # a bit above the highest finite value
plot_df$y_val  <- ifelse(plot_df$padj == 0, y_cap, -log10(plot_df$padj))
plot_df$capped <- plot_df$padj == 0
message(n_zero, " genes have padj == 0 (p-value below machine precision); ",
        "drawn as triangles at a ceiling of ", round(y_cap, 1), " in the volcano.")

status_levels <- c("Not significant", "Significant (padj < 0.05)")
plot_df$status <- factor(ifelse(plot_df$padj < 0.05, status_levels[2], status_levels[1]),
                         levels = status_levels)
status_cols <- setNames(c(col_ns, col_sig), status_levels)

n_sig <- sum(plot_df$status == status_levels[2])
n_all <- nrow(plot_df)

# --- 2. Volcano --------------------------------------------------------------
message("Plotting Volcano...")
volcano_df <- plot_df[order(plot_df$status), ]          # grey first, blue on top

p_volcano <- ggplot(volcano_df, aes(log2FoldChange, y_val)) +
  geom_point(data = subset(volcano_df, !capped), aes(colour = status),
             shape = 16, size = 0.7, alpha = 0.35) +
  geom_point(data = subset(volcano_df, capped), aes(colour = status),
             shape = 17, size = 1.8, show.legend = FALSE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, colour = ink_2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = ink_2) +
  scale_colour_manual(values = status_cols) +
  scale_x_continuous(breaks = scales::breaks_width(2)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.06))) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  labs(title = "Tumour vs normal: differential expression",
       subtitle = sprintf("%s of %s genes significant at padj < 0.05 (paired design, 113 patients)",
                          fmt(n_sig), fmt(n_all)),
       x = "log2 fold change (tumour vs normal)", y = "-log10(adjusted p-value)",
       caption = paste0("Positive log2 fold change: higher expression in tumour; negative: higher in normal tissue.\n",
                 "Dashed lines: padj = 0.05 and |log2 fold change| = 1 (visual reference only).",
                        if (n_zero > 0) paste0("\nTriangles: ", n_zero,
                          " genes with padj = 0 (p-value below machine precision), drawn at the top of the axis.")
                        else "")) +
  theme_portfolio()

save_fig(p_volcano, "01-volcano-primary.png")

# --- 3. MA plot --------------------------------------------------------------
message("Plotting MA...")
ma_df <- plot_df[order(plot_df$status, decreasing = TRUE), ]  # grey on top so the null cloud shows

p_ma <- ggplot(ma_df, aes(baseMean, log2FoldChange)) +
  geom_point(aes(colour = status), shape = 16, size = 0.7, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = ink_2) +
  scale_colour_manual(values = status_cols) +
  scale_x_log10(breaks = 10^(0:6), labels = compact) +
  scale_y_continuous(breaks = scales::breaks_width(2)) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  labs(title = "Expression level and effect size",
       subtitle = sprintf("Each point is one gene (n = %s); colour marks padj < 0.05", fmt(n_all)),
       x = "Mean of normalised counts (baseMean, log10 scale)",
       y = "log2 fold change (tumour vs normal)", caption = "Positive log2 fold change: higher expression in tumour; negative: higher in normal tissue.") +
  theme_portfolio()

save_fig(p_ma, "02-ma-primary.png")

# --- 4. Subsampling curve ---------------------------------------------------
message("Plotting Subsampling Curve...")
sub_df <- readRDS("results/analysis_5.rds")$summary
last   <- sub_df[sub_df$n == max(sub_df$n), ]
design_cols <- c("Paired" = col_paired, "Unpaired" = col_unpaired)

p_subsample <- ggplot(sub_df, aes(x = n)) +
  geom_ribbon(aes(ymin = paired_p10,   ymax = paired_p90),   fill = col_paired,   alpha = 0.12) +
  geom_ribbon(aes(ymin = unpaired_p10, ymax = unpaired_p90), fill = col_unpaired, alpha = 0.12) +
  geom_line(aes(y = paired_med,   colour = "Paired",   linetype = "Paired"),   linewidth = 0.8) +
  geom_line(aes(y = unpaired_med, colour = "Unpaired", linetype = "Unpaired"), linewidth = 0.8) +
  geom_point(aes(y = paired_med,   fill = "Paired"),   shape = 21, colour = "white", stroke = 0.8, size = 3) +
  geom_point(aes(y = unpaired_med, fill = "Unpaired"), shape = 21, colour = "white", stroke = 0.8, size = 3) +
  scale_colour_manual(values = design_cols) +
  scale_fill_manual(values = design_cols) +
  scale_linetype_manual(values = c("Paired" = "solid", "Unpaired" = "dashed")) +
  scale_x_continuous(breaks = sub_df$n) +
  scale_y_continuous(labels = compact, expand = expansion(mult = c(0, 0.04))) +
  expand_limits(y = 0) +
  labs(title = "Genes called significant as the cohort grows",
       subtitle = sprintf("At n = %d: %s genes paired vs %s unpaired",
                          last$n, fmt(last$paired_med), fmt(last$unpaired_med)),
       x = "Number of patients (n)", y = "Significant genes (padj < 0.05)",
       colour = "Design", fill = "Design", linetype = "Design",  # same title on all three = one legend
       caption = sprintf(paste0(
         "Lines: median over %d random subsets per n; bands: 10th-90th percentile.\n",
         "At n = %d the full cohort is a single subset, so no band is drawn.\n",
         "Ratio of medians (paired / unpaired): %.2f-%.2f across all n."),
         max(sub_df$reps), last$n, min(sub_df$ratio_med), max(sub_df$ratio_med))) +
  theme_portfolio()

save_fig(p_subsample, "03-subsampling-curve.png")

# --- 5. Leave-one-out stability ---------------------------------------------
message("Plotting LOO Stability...")
loo_df <- readRDS("results/analysis_4.rds")$genes

# The LOO list is the primary list: both counts come from independent stored objects.
stopifnot(nrow(loo_df) == n_sig)

n_loo   <- nrow(loo_df)
n_one   <- sum(loo_df$stability_fraction == 1)
n_below <- n_loo - n_one

# Explicit bin edges. In ggplot2 4.0.3, binwidth = 0.02 with boundary = 1 and
# closed = "left" puts the value 1.0 in the bin [0.98, 1.00], not in its own bin.
# Integer multiples of 0.02 keep the edge at 1.0 exact.
bin_w <- 0.02
lo_i  <- floor(min(loo_df$stability_fraction) / bin_w)

p_loo <- ggplot(loo_df, aes(stability_fraction)) +
  geom_histogram(breaks = (lo_i:51) * bin_w, closed = "left",
                 fill = col_paired, colour = "white", linewidth = 0.25) +
  annotate("text", x = 0.975, y = n_one, hjust = 1, vjust = 1, size = 3.4,
           colour = ink_2, lineheight = 1.15,
           label = sprintf("%s genes (%.1f%%)\nsignificant in all %d refits",
                           fmt(n_one), 100 * n_one / n_loo, 113)) +
  scale_x_continuous(breaks = seq(0.4, 1, by = 0.1)) +
  # pseudo-log keeps the 468 genes below 1.0 visible next to the 20,000 at 1.0
  scale_y_continuous(transform = "pseudo_log",
                     breaks = c(0, 10, 100, 1000, 10000),
                     labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(title = "How stable is the primary gene list?",
       subtitle = sprintf(paste0("Share of the 113 leave-one-patient-out refits in which each of the ",
                                 "%s genes stays significant (padj < 0.05)"), fmt(n_loo)),
       x = "Stability fraction", y = "Number of genes (pseudo-log scale)",
       caption = sprintf(paste0("The rightmost bar is the bin [1.00, 1.02): only genes at exactly 1.0. ",
                                "%s genes (%.1f%%) fall below 1.0."),
                         fmt(n_below), 100 * n_below / n_loo)) +
  theme_portfolio()

# Verify that the bin holding 1.0 contains exactly the genes at 1.0 (and exists).
top <- layer_data(p_loo, 1)
top <- top[top$xmin >= 1 - 1e-9, ]
stopifnot(nrow(top) == 1, top$count == n_one)

save_fig(p_loo, "04-leave-one-out-stability.png")

message("All figures saved under docs/figures/")