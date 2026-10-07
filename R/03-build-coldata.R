# R/03-build-coldata.R
#
# Build the sample-level table (`coldata`) and fix the factor levels that
# determine both the sign and the name of the condition coefficient. Also
# defines the coefficient name used by every downstream results() call.
#
# Inputs:  counts_raw, sample_annotation (from R/02-load-counts.R)
# Outputs: coldata, coef_name

# Fail early if the sample table and the matrix disagree on order.
stopifnot(identical(as.character(sample_annotation$sample_id),
                    colnames(counts_raw)))

coldata <- data.frame(
  patient   = factor(sample_annotation$patient),
  condition = factor(sample_annotation$condition,
                     levels = c("normal", "tumor")),
  row.names = sample_annotation$sample_id
)

# Structural checks. These verify the assumptions the model makes, not
# the appearance of the table.
stopifnot(is.factor(coldata$patient))
stopifnot(is.factor(coldata$condition))
stopifnot(identical(levels(coldata$condition), c("normal", "tumor")))

# The paired design requires exactly one tumour and one normal per
# patient. A total-count check would not catch a patient with two
# tumours and another patient with two normals.
stopifnot(all(table(coldata$patient, coldata$condition) == 1))

# Coefficient name fixed here, guarded in R/05 against resultsNames(dds).
coef_name <- "condition_tumor_vs_normal"

message("coldata: ", nrow(coldata), " samples | ",
        nlevels(coldata$patient), " patients | ",
        nlevels(coldata$condition), " conditions")
