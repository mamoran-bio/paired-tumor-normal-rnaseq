# This function gets all patients with paired tissue (tumor/normal) from a TCGA projects.
# Barcode is used to avoid gaps in submitter_id (first 12 characters) and sample_type_id data (14, 15 characters)
# Any error can by checked in the script's log
# The output is rows with n_rows, n_aliquots, n_patients, n_tumor, n_normal, n_paired

#------------------------------------------------------------------------------------
# Libraries
#------------------------------------------------------------------------------------

library(recount3)

#------------------------------------------------------------------------------------
# Counting paired patients function
#------------------------------------------------------------------------------------

count_pairs <- function(proj_row) {
	proj_name <- proj_row$project
	message("Processing ", proj_name, "...")

	# NA row with the structure of a result you are looking for,
	# These values are the output if there is an error
	na_row <- data.frame(
		project		= proj_name,
		n_rows		= NA_integer_,
		n_aliquots	= NA_integer_,
		n_patients	= NA_integer_,
		n_tumor		= NA_integer_,
		n_normal	= NA_integer_,
		n_paired	= NA_integer_,
		stringsAsFactors = FALSE
		)
	result <- tryCatch({

		# What is extracted from project URL 
		url <- locate_url(
			project		= proj_name,
			project_home	= proj_row$project_home,
			type		= "metadata",
			organism	= "human"
		)

		# Downloading the data
		md <- read_metadata(file_retrieve(url))
		
		# Getting patient and sample type from barcode
		barcode <- md$tcga.tcga_barcode
		patient <- substr(barcode, 1, 12)
		stype	<- substr(barcode, 14, 15)
		
		data.frame(
			project		= proj_name,
			n_rows		= nrow(md),
			n_aliquots	= length(unique(barcode)),
			n_patients	= length(unique(patient)),
			n_tumor		= length(unique(patient[stype == "01"])),
			n_normal	= length(unique(patient[stype == "11"])),
			n_paired	= length(intersect(
							patient[stype == "01"],
							patient[stype == "11"]
						)),
			stringsAsFactors = FALSE
		)
	}, error = function(e) {
		message("FAILED On", proj_name, ": ", conditionMessage(e))
	na_row
	})
	
	result
	
}

#-------------------------------------------------------------------------
# Execution
#-------------------------------------------------------------------------

proj <- available_projects()
tcga <- subset(proj, file_source == "tcga")

results <- do.call(
  rbind,
  lapply(seq_len(nrow(tcga)), function(i) count_pairs(tcga[i, ]))
)


#---------------------------------------------------------------------------
# Check results
#---------------------------------------------------------------------------


check_results <- function(results) {

  n_total <- nrow(results)

  # --- 1. Download failures ---
  failed <- subset(results, is.na(n_paired))
  if (nrow(failed) == 0) {
    message("OK    download: all ", n_total, " projects ran without error")
  } else {
    message("FAIL  download: ", nrow(failed), " of ", n_total,
            " projects failed -> ",
            paste(failed$project, collapse = ", "))
  }

  # --- 2. Internal consistency ---
  bad <- subset(results,
                !is.na(n_paired) & (
                  n_tumor > n_patients |
                  n_normal > n_patients |
                  n_paired > n_tumor |
                  n_paired > n_normal
                ))
  if (nrow(bad) == 0) {
    message("OK    consistency: no row violates the inequalities")
  } else {
    message("FAIL  consistency: ", nrow(bad),
            " inconsistent projects -> ",
            paste(bad$project, collapse = ", "))
  }

  # --- 3. BRCA, the known anchor ---
  brca_pairs <- results$n_paired[results$project == "BRCA"]
  if (length(brca_pairs) == 1 && !is.na(brca_pairs) && brca_pairs == 113) {
    message("OK    BRCA: n_paired = 113, as expected")
  } else {
    message("FAIL  BRCA: n_paired = ", brca_pairs,
            " (expected 113)")
  }

  # --- 4. Column types ---
  cols <- c("n_rows", "n_aliquots", "n_patients",
            "n_tumor", "n_normal", "n_paired")
  bad_types <- cols[!sapply(results[cols], is.integer)]
  if (length(bad_types) == 0) {
    message("OK    types: all 6 count columns are integer")
  } else {
    message("FAIL  types: not integer -> ",
            paste(bad_types, collapse = ", "))
  }

  # --- 5. Duplicated rows ---
  dup <- results$project[duplicated(results$project)]
  if (length(dup) == 0) {
    message("OK    uniqueness: ", n_total, " distinct projects")
  } else {
    message("FAIL  uniqueness: duplicated projects -> ",
            paste(dup, collapse = ", "))
  }

  invisible(NULL)
}

check_results(results)

# Save, sorted by number of paired patients (descending)

results <- results[order(-results$n_paired), ]
write.csv(results, "docs/tissue-pair-counts.csv", row.names = FALSE)
message("Written: docs/tissue-pair-counts.csv")
