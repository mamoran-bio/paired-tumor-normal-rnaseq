# R/04-filter-genes.R
#
# Build the paired design matrix, compare the two forms of
# edgeR::filterByExpr(), and produce `counts_filtered` with the form
# chosen on 2026-10-06. The choice is recorded in the amendment log
# entry of 2026-10-07.
#
# Decision: the filter uses the biological factor (`group = condition`),
# not the paired design matrix.
#
# Rationale in the amendment log entry of 2026-10-07; the summary here is:
#
#   `design =` sets the minimum number of samples a gene must be
#   expressed in to 1 / max(leverage). On this paired design that
#   minimum degenerates to ~2, sized for the patient coefficients
#   rather than for the tumour-versus-normal comparison.
#
#   `group =` sizes the filter to the smaller of the two conditions
#   (113 samples). edgeR's rule for groups this size gives a minimum of
#   10 + (113 - 10) * 0.7 = 82.1, i.e. a gene must be expressed in 83 or
#   more samples. That minimum corresponds to the tumour-versus-normal
#   comparison the analysis makes.
#
#   `group =` does not use the tumour-normal difference: it counts in
#   how many samples a gene is expressed, irrespective of condition, and
#   the minimum depends on group sizes only.
#
#   Cost: a gene expressed in fewer than 83 samples in total is dropped
#   even if all its expression is in tumours. The genes kept by
#   `design =` and not by `group =` are reported below.
#
# Inputs:  counts_raw, coldata (from R/03-build-coldata.R)
# Outputs: design_paired, keep_design, keep_group, counts_filtered

# --- Paired design matrix --------------------------------------------------
design_paired <- model.matrix(~ patient + condition, data = coldata)

# --- Leverage diagnostic ---------------------------------------------------
# Leverage is the diagonal of the hat matrix H = X (X'X)^-1 X'. It
# measures how much each sample weighs in estimating its own
# coefficients. filterByExpr(design =) uses 1 / max(leverage) as the
# minimum number of samples a gene must be expressed in.
#
# hat() is the same computation filterByExpr() performs internally, so
# the diagnostic reproduces that call exactly.
hat_values <- hat(design_paired)

message("design_paired: ", nrow(design_paired), " x ",
        ncol(design_paired))
message("leverage: min ", round(min(hat_values), 6),
        " | max ", round(max(hat_values), 6),
        " | 1 / max = ", round(1 / max(hat_values), 2),
        " samples minimum for filterByExpr(design=)")

# --- Both filters ----------------------------------------------------------
# Both are computed so the comparison below can be reported. Only
# `keep_group` is used to build counts_filtered.
keep_design <- edgeR::filterByExpr(counts_raw, design = design_paired)
keep_group  <- edgeR::filterByExpr(counts_raw, group  = coldata$condition)

message("keep_design: ", sum(keep_design), " genes")
message("keep_group: ", sum(keep_group), " genes")

# --- Coherence check -------------------------------------------------------
# `group` requires a higher minimum sample count than `design`, with the
# same CPM threshold, so every gene kept by `group` must also be kept by
# `design`. If a gene passes `group` and fails `design`, the two calls
# are not doing what they appear to and the numbers below cannot be
# interpreted.
n_group_only <- sum(keep_group & !keep_design)
message("genes kept by group and not by design: ", n_group_only)
stopifnot(n_group_only == 0)

message("crosstab design x group: ",
        "design only = ", sum(keep_design & !keep_group), " | ",
        "both = ", sum(keep_design & keep_group), " | ",
        "neither = ", sum(!keep_design & !keep_group))

# --- Minimum CPM threshold -------------------------------------------------
# filterByExpr's default min.count is 10. Expressing it as CPM over the
# median library size gives the expression threshold a gene must clear.
min_cpm <- 10 / median(colSums(counts_raw)) * 1e6
message("minimum CPM threshold: ", round(min_cpm, 3))

# --- Apply the chosen filter -----------------------------------------------
counts_filtered <- counts_raw[keep_group, , drop = FALSE]

message("Genes before filter: ", nrow(counts_raw),
        " | after: ", nrow(counts_filtered),
        " (", round(100 * nrow(counts_filtered) / nrow(counts_raw), 1),
        "% retained)")
