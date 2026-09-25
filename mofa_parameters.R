# *****************************************************************************
# User-configurable settings for MOFA_analysis.R.
#
# The variable timepoint specifies the follow-up visit to analyze. The clinical
# column containing visit codes is configured with timepoint_col, and the codes
# and display labels are configured with timepoint_labels.
#
# If use_change_scores is TRUE, omics values are converted to change
# scores by subtracting the value identified by the "baseline" entry in
# timepoint_labels from the target visit.
#
# If use_change_scores is FALSE, omics values from the target visit are used
# directly, without baseline subtraction.
#
# The script supports any combination of input data matrices
# (methylation, proteomics, metabolomics). The file paths for these are
# defined below. To exclude a data type, simply comment out the
# corresponding file path; the script will automatically use only the
# available data matrices. If all omics file paths are omitted, the script
# performs a metadata-only run: it writes the target-visit sample metadata and
# run-provenance files but skips MOFA preprocessing and model training.
# *****************************************************************************

#------------ Configure run by setting the values of the variables below

# Location of directory where to store results and output - modify as desired
# for your local run environment
output_dir <- "D:/R_sessions/PROSPR/MOFA"
# output_dir <- "D:/R_sessions/PROSPR/MOFA/reliable_probes"
# output_dir <- "D:/R_sessions/PROSPR/MOFA/foo"

# Study name used in run-specific documentation.
study_name <- "CALERIE"

# Names of the subject identifier, treatment, and timepoint columns in the
# clinical data file. These settings allow the analysis code to remain
# independent of study-specific clinical column names.
subject_id_col <- "DEID"
treatment_col <- "CR"
timepoint_col <- "fumo"

# Map user-facing timepoint labels to the corresponding codes in the clinical
# data. This vector must be named, its values must match the unique non-missing
# values in timepoint_col, and exactly one entry must be named "baseline".
timepoint_labels <- setNames(c(0, 12, 24, 3, 6, 18, 98, 99), 
		c("baseline", "12mo", "24mo", "3mo", "6mo", "18mo", "BLs1", "BLs2"))

# Map user-facing treatment labels to the corresponding codes in the clinical
# data. This vector must be named and its values must match the unique
# non-missing values in treatment_col. The names are used in sample IDs.
treatment_labels <- setNames(c(0, 1), c("AL", "CR"))

# Specify the target visit using one of the values in timepoint_labels.
timepoint <- 24

# If TRUE, compute (timepoint - baseline) change scores. If FALSE, use the
# source data directly for samples where timepoint_col == timepoint.
use_change_scores <- TRUE

# Set scaling factor controlling how many features to retain per view
# (relative to the smallest view) 
max_factor <- 15

# Set an absolute upper bound on the number of features per view - useful
# primarily when running MOFA with only the methylation data view, to
# avoid performance degradation due to hundreds of thousands of features
max_features <- 15000

# Indicates if the data should be scaled. According to the MOFA paper, scaling
# is not needed. Scaling the data, however, may be helpful for some plots; for
# example, clusters may be more distinguishable in heatmaps when using z-scores.
scale_views <- TRUE

# The maximum number of variational inference iterations allowed during model 
# training in MOFA. The default value is 1000 but we use a higher setting in  
# case it is required for convergence. A small value can be entered here to
# test code changes, so that the MOFA completes quickly.
mofa_max_iter = 4000

#------------------------------------------------------------------------------

# *****************************************************************************
# File locations. Edit according to your execution environment. 
#
# If you want to exclude one or more views from MOFA, simply comment the relevant 
# file name variables.
# *****************************************************************************
# Location of the methylation M-values file.
methy_fn = file.path("E:/Data/PROSPR-FAST/CALERIE/methylation",
		"GRSet_fully_filtered_bmiq_chunk_mvals.rds")
#methy_fn = file.path("E:/Data/PROSPR-FAST/CALERIE/methylation",
#		"CALERIE_GRSet_CPR_processing_reliable_subset_mvals.rds")

# Location of the proteomics data file.
prot_fn = file.path("E:/Data/PROSPR-FAST/CALERIE/proteomics",
		"CALERIE_somascan_processed_unique_human_proteins_matrix.csv")

# Location of the Metabolon data and annotation files.
metab_data_fn = file.path("E:/Data/PROSPR-FAST/CALERIE/metabolomics",
		"CALERIE_metabolon.csv")
metab_annot_fn = file.path("E:/Data/PROSPR-FAST/CALERIE/metabolomics",
		"metabolon_annotations.csv")

# Location of the clinical/phenotype data file. This *MUST* be defined and 
# cannot be commented out.
clin_fn = file.path("E:/Data/PROSPR-FAST/CALERIE/clinical",
		"CALERIE_clinical_w_omics_crosswalks.csv")

# Names of variables specifying file locations, expected to be defined above.
resource_vars <- c("methy_fn", "prot_fn",
		"metab_data_fn", "metab_annot_fn",
		"clin_fn")

# Definitions for the optional omics data views.
#
# omics_input_defs is a named list. Each top-level entry defines one possible
# MOFA data view; the top-level entry name, such as "proteomics", becomes the
# view name used inside the MOFA model and in generated documentation.
#
# Each view definition must contain the following slots:
#
#   moniker:
#     Short label used in result-folder names. For example, the proteomics
#     moniker "pr" appears in folder tags such as "views=pr+mb". In this 
#	  manner, documentation about how MOFA was run is directly embedded in
#	  the results folder.
#
#   file_var:
#     Character string naming the file-path variable for this view. For example,
#     file_var = "prot_fn" tells MOFA_analysis.R to look for an object named
#     prot_fn above. If that object is missing, commented out, empty, or NA, this
#     view is skipped. The value of file_var should also appear in resource_vars.
#
#   barcode_col:
#     Name of the column in the clinical/phenotype table that links each
#     clinical row to the corresponding sample column in this omics matrix. The
#     values in this clinical column must match column names in the matrix
#     returned by reader.
#
#   reader:
#     Function used to read the view-specific input file. It must accept one
#     argument, fn, which is the file path named by file_var, and it must return
#     a numeric feature-by-sample matrix with features in rows and sample
#     barcodes in columns. Definitions for these methods *MUST* be provided in 
#	  the file file_readers.R.
#
#   readme_notes:
#     Named list of text snippets inserted into the run-specific README.txt when
#     this view is included. It must contain one character string named
#     change_scores and one named timepoint_values. MOFA_analysis.R selects
#     change_scores when use_change_scores == TRUE and timepoint_values when
#     use_change_scores == FALSE.
#
# To add a new omics view, define its file-path variable above, add that variable
# name to resource_vars, and add a new named entry here with the same five slots.
omics_input_defs <- list(
		methylation = list(
				moniker = "me",
				file_var = "methy_fn",
				barcode_col = "blood_DNAm_barcode",
				reader = function(fn) readMethylation(fn),
				readme_notes = list(
						change_scores = "Methylation data are reported as M-values.",
						timepoint_values = "Methylation data are reported as M-values."
				)
		),
		proteomics = list(
				moniker = "pr",
				file_var = "prot_fn",
				barcode_col = "proteomics_barcode",
				reader = function(fn) readProteomics(fn),
				readme_notes = list(
						change_scores = paste("Proteomics data, from Somascan, are already",
								"log2 transformed."),
						timepoint_values = paste("Proteomics data, from Somascan, are already",
								"log2 transformed.")
				)
		),
		metabolomics = list(
				moniker = "mb",
				file_var = "metab_data_fn",
				barcode_col = "metabolomics_barcode",
				reader = function(fn) readMetabolomics(fn)[["data"]],
				readme_notes = list(
						change_scores = paste("Metabolomics data in", study_name,
								"are reported in Metabolon median-centered units; accordingly,",
								"we log2 transform them prior to change-score calculation."),
						timepoint_values = paste("Metabolomics data in", study_name,
								"are reported in Metabolon median-centered units; accordingly,",
								"we log2 transform them prior to analysis.")
				)
		)
)
