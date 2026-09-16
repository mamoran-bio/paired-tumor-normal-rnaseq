# TCGA-BRCA metadata selection from recount3
# Only metadata - count matrices are not downloaded
# Data curation to select patients with paired tissue
# Patient ID and sample type are taken from barcode, not from
# submitter_id / sample_type_id: those have 10 missing values in BRCA,
# one of which hid a complete tumour/normal pair.

library(recount3)

proj <- available_projects()
tcga <- subset(proj, file_source == "tcga")
brca <- subset(tcga, project  == "BRCA")

url <- locate_url(
	project = brca$project,
	project_home = brca$project_home,
	type = "metadata",
	organism = "human"
)

md <- read_metadata(file_retrieve(url))

message("Metadata loaded: ", nrow(md), " rows ", ncol(md), " columns")
