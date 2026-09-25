# *****************************************************************************
# Open files, prepare data, and run MOFA analysis.
#
# User-configurable analysis settings, input-file locations, and optional omics
# view definitions are maintained in MOFA_parameters.R.
#
# After configuring the analysis settings in MOFA_parameters.R, run this script
# from an R session:
#
#   source("path_to_mofa_code/MOFA_analysis.R", local = TRUE, chdir = TRUE)
#
# The chdir = TRUE argument temporarily sets the working directory to the folder
# containing this script, allowing it to locate the supporting R files. Any
# relative input and output paths are also resolved against that folder.
#
# Alternatively, run the script from the command line:
#
#   Rscript.exe -e "source('path_to_mofa_code/MOFA_analysis.R', local = TRUE, 
#				chdir = TRUE)"
#
# Running from an R session is recommended for capturing the model-training log.
# When MOFA2::run_mofa() uses a separate process for Python execution, its progress
# output may not reach R's standard output stream and therefore may not be
# captured by the script's sink() call. Logging may consequently be incomplete
# when the script is launched with Rscript.exe.
# *****************************************************************************

library(MOFA2)
library(dplyr)
library(tidyr)
library(xlsx)

# Source the file readers that parse the input files
source("file_readers.R", local = TRUE)
# Load user-configurable settings
source("MOFA_parameters.R", local = TRUE)
source("methods.R", local = TRUE)

# Read and validate the clinical data and shared study settings.
pheno_d <- readPhenoData(getResourceFile("clin_fn", required = TRUE))
clinical_context <- validateClinicalConfiguration(pheno_d)
baseline_timepoint <- clinical_context$baseline_timepoint
target_timepoint_label <- clinical_context$target_timepoint_label

validateOmicsInputDefinitions(omics_input_defs)

# Read only omics views whose configured file variables are defined.
loaded_info <- loadConfiguredOmicsViews(omics_input_defs)
loaded_views <- loaded_info$loaded_views
barcode_cols <- loaded_info$barcode_cols

if (length(loaded_views) == 0){
  message("No omics data file names are defined in MOFA_parameters.R. ",
      "Writing sample metadata only; MOFA will not be run.")

  sample_metadata <- buildTargetVisitSampleMetadata(pheno_d)
  sample_metadata <- data.frame(sample_metadata, sample = rownames(sample_metadata))

  res_folder <- file.path(output_dir,
      paste0(format(Sys.time(), "%Y_%b_%d_%H_%M_%S"),
          "_visit=", target_timepoint_label, "_metadata_only"))
  dir.create(res_folder, recursive = TRUE)
  log_folder <- file.path(res_folder, "logs")
  dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)

  sample_metadata_xlsx <- file.path(res_folder, "sample_metadata.xlsx")
  xlsx::write.xlsx(as.data.frame(sample_metadata), sample_metadata_xlsx,
      sheetName = study_name, row.names = FALSE)
  saveRDS(sample_metadata, file.path(res_folder, "sample_metadata.rds"))

  session_file <- file.path(log_folder, "sessionInfo.txt")
  writeLines(capture.output(sessionInfo()), session_file)
  file.copy(
      from = c("MOFA_analysis.R", "MOFA_parameters.R", "file_readers.R", "methods.R"),
      to = log_folder,
      overwrite = TRUE
  )

  readme_file <- file.path(res_folder, "README.txt")
  readme_text <- paste(c(
      paste(study_name, "sample metadata were exported for target visit",
          target_timepoint_label, "."),
      "No omics data file names were defined in MOFA_parameters.R, so",
      "MOFA preprocessing and model training were skipped.",
      "The result files are sample_metadata.xlsx and sample_metadata.rds.",
      "The logs folder contains sessionInfo.txt, which records the R/package",
      "environment, plus copies of MOFA_analysis.R, MOFA_parameters.R,",
      "file_readers.R, and methods.R used for the run."
  ), collapse = " ")
  writeLines(strwrap(readme_text, width = 80), con = readme_file)

  message("Metadata-only run completed. Results written to: ", res_folder)
} else {

# Prepare target-timepoint values or change scores for every loaded view.
data_views <- prepareConfiguredOmicsViews(
  loaded_views = loaded_views,
  pheno_df = pheno_d,
  barcode_cols = barcode_cols,
  timepoint = timepoint,
  use_change_scores = use_change_scores
)
if (length(data_views) == 0L) {
  stop("No omics views have samples at the requested timepoint.", call. = FALSE)
}

# Return the indices of the most variable features in a matrix. If there are
# fewer than two samples, feature variance cannot be estimated, so the first
# max_features rows are returned.
getTopVariableFeatureIndex <- function(view_data, max_features){
	max_features <- min(max_features, nrow(view_data))
	
	if (max_features >= nrow(view_data)){
		return(seq_len(nrow(view_data)))
	}
	if (ncol(view_data) < 2){
		return(seq_len(max_features))
	}
	
	row_vars <- apply(view_data, 1, var, na.rm = TRUE)
	row_vars[is.na(row_vars)] <- -Inf
	order(row_vars, decreasing = TRUE)[seq_len(max_features)]
}

# MOFA recommends that data views have a comparable number of features, to avoid
# one view dominating factor fitting. We therefore reduce large views by keeping
# only the most variable features. A view may have at most max_factor times as
# many features as the smallest view, and never more than max_features.
filterMostVariableFeatures <- function(data_views, max_factor=5,
		max_features = 10000){
	# Get number of features (rows) for each omics view
	feature_counts <- vapply(data_views, nrow, integer(1))
	
	if (any(feature_counts == 0)){
		stop("At least one omics view has zero features after preprocessing.")
	}
	
	min_feat_num <- min(feature_counts)
	
	# Iterate over each omics view (data matrix)
	lapply(data_views, function(view_data){
				# Compute number of features to retain for this view:
				# limited by:
				#   - max_factor * min features across views
				#   - global max_features cap
				#   - total features available in this view
				vf_num <- min(max_factor * min_feat_num, max_features,
						nrow(view_data))
				vf_num <- floor(vf_num)
				
				if (vf_num >= nrow(view_data)){
					return(view_data)
				}
				
				feat_idx <- getTopVariableFeatureIndex(view_data, vf_num)
				view_data[feat_idx, , drop = FALSE]
			})
}

# Fill missing sample columns with NA values.
# Given a features x samples matrix and a desired sample-ID vector, this helper
# adds any absent samples as all-NA columns and returns columns in sids order.
fillEmptyCols <- function(dmat, sids) {
	if (missing(sids) || length(sids) == 0L) {
		stop("'sids' must be a non-empty character vector of sample IDs", call. = FALSE)
	}
	
	rn <- rownames(dmat)
	mat <- as.matrix(dmat)
	if (!is.null(rn)) rownames(mat) <- rn
	
	sids <- as.character(sids)
	if (anyDuplicated(sids)) sids <- unique(sids)
	
	missing_cols <- setdiff(sids, colnames(mat))
	if (length(missing_cols) > 0L) {
		nc_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = length(missing_cols))
		if (!is.null(rownames(mat))) rownames(nc_mat) <- rownames(mat)
		colnames(nc_mat) <- missing_cols
		mat <- cbind(mat, nc_mat)
	}
	
	mat[, sids, drop = FALSE]
}

# Apply feature filtering to each view, enforcing both relative (max_factor)
# and absolute (max_features) limits.
data_views <- filterMostVariableFeatures(data_views, max_factor,
		max_features = max_features)

# Prep data matrices by filling in missing columns. This is a requirement of
# the MOFA package: if a data view is missing data for some samples, these need
# to be added with NA values.
sample_ids <- Reduce(union, lapply(data_views, colnames))
if (length(sample_ids) == 0){
	stop("No samples are available after preprocessing the omics views.")
}

# Prep data matrices by filling in missing columns. This is a requirement of
# the MOFA package: if a data view is missing data for some samples, these
# need to be added with NA values.
data_views <- lapply(data_views, function(view_data){
			view_data <- fillEmptyCols(view_data, sample_ids)
			view_data[, sample_ids, drop = FALSE]
		})

# Create MOFA object, comprising all available data matrices (aka, views).
# All samples are analyzed together without separate MOFA groups.
mdd <- create_mofa(data_views)
views_names(mdd) <- names(data_views)

# Set data options.
data_opts <- get_default_data_options(mdd)
if (scale_views){
	data_opts$center_groups <- TRUE
	data_opts$scale_groups <- TRUE
}
data_opts$use_float32 <- FALSE

# Set model options.
model_opts <- get_default_model_options(mdd)
model_opts$num_factors <- 45 # If we want to manually set the factor num

# Set training options.
train_opts <- get_default_training_options(mdd)
# variance explained threshold for dropping factors. For example, a value of 0.01 
# implies that factors explaining less than 1% of variance (in each view) will be 
# dropped. Here we use the default setiing of -1 (no dropping of factors).
train_opts$drop_factor_threshold <- -1 
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- mofa_max_iter # maximum number of model training iterations

# Prepare the MOFA object and add the metadata data frame, with one row per
# sample.
mo <- prepare_mofa(object = mdd, data_options = data_opts,
		model_options = model_opts, training_options = train_opts)

meta_d <- buildTargetVisitSampleMetadata(pheno_d)

missing_meta <- setdiff(sample_ids, rownames(meta_d))
if (length(missing_meta) > 0){
	stop("Missing phenotypic metadata for samples: ",
			paste(missing_meta, collapse = ", "))
}

meta_d <- meta_d[sample_ids, , drop = FALSE]
meta_d <- data.frame(meta_d, sample = rownames(meta_d))
samples_metadata(mo) <- meta_d

# Train the MOFA model.
#
# Note the parameter setting "use_basilisk = TRUE". This forces a local Python
# environment to be installed in order to download and run the mofapy2 package,
# which performs the MOFA computations.
#
# This setup occurs only once, the first time the code is run on a given
# machine. The initial run may take significantly longer, as it involves
# downloading and configuring several gigabytes of dependencies.
#
# Alternatively, if a suitable Python environment already exists, you can
# install mofapy2 manually and run the model with "use_basilisk = FALSE".

# Tag the output folder with the current date and time, for logging purposes.
# Also save the output of the MOFA method.
score_tag <- ifelse(use_change_scores, "change_scores", "timepoint_values")
views_tag <- paste(sapply(views_names(mo), function(vn){
					omics_input_defs[[vn]]$moniker}), 
		collapse="+")
res_folder <- paste(output_dir,
		paste0(format(Sys.time(), "%Y_%b_%d_%H_%M_%S"),
				"_visit=", target_timepoint_label, "_", score_tag,
				"_views=", views_tag),
		sep = "/")
log_folder <- file.path(res_folder, "logs")
res_file <- file.path(res_folder, "model.hdf5")
log_file <- file.path(log_folder, "mofa_out.txt")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)
sink(file = log_file)
mo.tr <- tryCatch(
		run_mofa(mo, res_file, use_basilisk = TRUE),
		finally = sink()
)

# Export fitted factor values, feature weights, and sample metadata in both
# Excel and RDS formats. Excel workbooks are written with the xlsx package.
factor_values <- get_factors(mo.tr)[[1]]
factor_values_xlsx <- file.path(res_folder, "factor_values.xlsx")
xlsx::write.xlsx(as.data.frame(factor_values), factor_values_xlsx,
		sheetName = study_name, row.names = TRUE)
saveRDS(factor_values, file.path(res_folder, "factor_values.rds"))

view_names <- views_names(mo.tr)
feature_weights <- setNames(lapply(view_names, function(view_name){
			get_weights(mo.tr)[[view_name]]
		}), view_names)
feature_weights_xlsx <- file.path(res_folder, "feature_weights.xlsx")
for (i in seq_along(feature_weights)){
	xlsx::write.xlsx(as.data.frame(feature_weights[[i]]), feature_weights_xlsx,
			sheetName = names(feature_weights)[i], row.names = TRUE,
			append = i > 1)
}
saveRDS(feature_weights, file.path(res_folder, "feature_weights.rds"))

sample_metadata <- samples_metadata(mo.tr)
sample_metadata_xlsx <- file.path(res_folder, "sample_metadata.xlsx")
xlsx::write.xlsx(as.data.frame(sample_metadata), sample_metadata_xlsx,
		sheetName = study_name, row.names = FALSE)
saveRDS(sample_metadata, file.path(res_folder, "sample_metadata.rds"))

# Logging/provenance - write the R session information and copy the analysis
# script, file-reader file, methods file, and parameters file used for the run.
session_file <- file.path(log_folder, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)
file.copy(
		from = c("MOFA_analysis.R", "MOFA_parameters.R",
				"file_readers.R", "methods.R"),
		to = log_folder,
		overwrite = TRUE
)
# Document what was done.
readme_file <- file.path(res_folder, "README.txt")
score_desc <- if (use_change_scores){
			paste0("change scores (", target_timepoint_label, " - baseline)")
		} else {
			paste0("source values from samples where ", timepoint_col, " == ",
					timepoint, " (", target_timepoint_label, ")")
		}

included_views <- paste(names(data_views), collapse = ", ")
# Pull view-specific prose from MOFA_parameters.R for the selected score mode.
data_notes <- vapply(names(data_views), function(view_name){
			omics_input_defs[[view_name]]$readme_notes[[score_tag]]
		}, character(1))

readme_text <- paste(c(
				paste(study_name, "MOFA was run on", score_desc, "."),
				paste("Included omics views:", included_views, "."),
				data_notes,
				"Each included omics data table becomes a separate MOFA view.",
				"All samples were analyzed together without separate MOFA groups.",
				"Result files in this folder include model.hdf5, factor_values.xlsx,",
				"factor_values.rds, feature_weights.xlsx, feature_weights.rds,",
				"sample_metadata.xlsx, and sample_metadata.rds.",
				"The logs folder contains mofa_out.txt (captured model-training output),",
				"sessionInfo.txt (the R/package environment), and copies of",
				"MOFA_analysis.R, MOFA_parameters.R, file_readers.R, and methods.R",
				"used for the run."
		), collapse = " ")

writeLines(strwrap(readme_text, width = 80), con = readme_file)
}
