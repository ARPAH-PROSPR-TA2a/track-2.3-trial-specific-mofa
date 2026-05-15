# *****************************************************************************
# Open files, prepare data, and run MOFA analysis. Three variables control
# how the analysis is configured:
#   timepoint, use_change_scores, use_groups
#
# The variable timepoint specifies the follow-up visit to analyze
# (identified by fumo == timepoint). This determines which samples are
# used as the target condition.
#
# If use_change_scores is TRUE, omics values are converted to change
# scores by subtracting baseline values (fumo == 0) from the target
# visit (fumo == timepoint).
#
# If use_change_scores is FALSE, omics values from the target visit
# (fumo == timepoint) are used directly, without baseline subtraction.
#
# If use_groups is TRUE, treatment groups (caloric restriction or normal
# diet) are modeled as separate MOFA groups. If FALSE, all samples are
# analyzed together without grouping.
#
# If use_residuals is TRUE, then the omics data matrices are residualized
# using as covariates the columns specified in the variable covariates.
# Residualization takes place prior to any downstream computations (e.g.,
# calculation of change scores).
#
# The script supports any combination of input data matrices
# (methylation, proteomics, metabolomics). The file paths for these are
# defined in resources.R. To exclude a data type, simply comment out the
# corresponding file path in resources.R; the script will automatically
# use only the available data matrices.
#
# To ensure that savehistory() captures all executed commands, run this code
# by pasting it into a fresh R GUI session.
# *****************************************************************************

Sys.setenv(R_HISTSIZE = 2000) # to make sure all history is saved at the end.

library(MOFA2)
library(dplyr)
library(tidyr)

########################################
#------------ Configure run by setting the values of the variables below
# Location of source code root directory - modify as needed for your environment
code_dir <- "D:/git/Belsky/track-2.3-trial-specific-mofa/CALERIE"

# Location of directory where to store results and output - modify as desired
# for your local run environment
output_dir <- "D:/R_sessions/PROSPR/MOFA"

# Specify the target visit timepoint, either 12 or 24
timepoint <- 24

# If TRUE, compute (timepoint - baseline) change scores. If FALSE, use the
# source data directly for samples where fumo == timepoint.
use_change_scores <- TRUE

# If TRUE, run MOFA with AL and CR as separate groups. If FALSE, do not pass
# groups to MOFA, so all samples are analyzed together.
use_groups <- TRUE

# If TRUE, each omics matrix is residualized immediately after it is read from
# file and before any downstream processing, such as change-score calculation,
# timepoint extraction, variable-feature filtering, or MOFA object creation.
use_residuals <- FALSE

# Covariates to use when residualizing omics matrices. The var_names column must
# contain column names from the phenotypic data table loaded from calerie_clin_fn.
# The is_factor column indicates whether each variable should be coerced to a
# factor before residualization. Leave this empty when use_residuals is FALSE.
covariates <- data.frame(
		var_names = c("female", "race", "agebl"),
		is_factor = c(TRUE, TRUE, FALSE),
		stringsAsFactors = FALSE
)

# If TRUE, feature-wise residualization is parallelized when the fast matrix-wide
# residualization path cannot be used. The fast path is still preferred whenever
# possible because it is usually faster than launching worker processes.
use_parallel_residualization <- TRUE

# Number of worker processes to use for parallel feature-wise residualization.
# This leaves one CPU core free when more than one core is available. If R cannot
# detect the number of cores, the script falls back to serial residualization.
detected_residualization_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_residualization_cores)){
	detected_residualization_cores <- 1
}
n_residualization_cores <- max(1, detected_residualization_cores - 1)
#------------------------------------------------------------------------------
########################################

# Source files. 
source(file.path(code_dir, "resources.R"))
source(file.path(code_dir, "methods.R"))
# Names of variables specifying file locations, expected to be defined in 
# resources.R
resource_vars <- c("calerie_methy_fn", "calerie_prot_fn",
		"calerie_metab_data_fn", "calerie_metab_annot_fn",
		"calerie_clin_fn")

# Return TRUE only if the variable is defined in resources.R and contains a
# non-empty file name.
resourceFileSpecified <- function(var_name){
	if (!exists(var_name)){
		return(FALSE)
	}	
	value <- get(var_name)
	is.character(value) && length(value) == 1 && !is.na(value) && nzchar(value)
}

# Get a file name defined in resources.R. Required files trigger an error if
# missing; optional files return NULL if missing.
getResourceFile <- function(var_name, required = FALSE){
	if (!resourceFileSpecified(var_name)){
		if (required){
			stop(paste0("Required variable '", var_name,
							"' is not defined in resources.R."))
		}
		return(NULL)
	}
	get(var_name)
}

# Read phenotypic data. This file is required.
pheno_d <- readPhenoData(getResourceFile("calerie_clin_fn", required = TRUE))

# Validate the residualization settings before any omics tables are loaded. This
# catches misspelled metadata variables early, rather than after partial work has
# already been done.
validateResidualizationInputs <- function(pheno_df, covariates){
	# No validation is needed when residualization is disabled.
	if (!use_residuals){
		return(invisible(TRUE))
	}
	
	# The covariates object must be a data frame with exactly the expected fields.
	if (!is.data.frame(covariates) ||
			!all(c("var_names", "is_factor") %in% names(covariates))){
		stop("When use_residuals is TRUE, covariates must be a data frame with columns 'var_names' and 'is_factor'.")
	}
	
	# At least one covariate is required in order to residualize the omics data.
	if (nrow(covariates) == 0){
		stop("use_residuals is TRUE, but no residualization covariates were supplied.")
	}
	
	# Convert user-supplied fields to simple vectors so later checks are reliable.
	covariates$var_names <- as.character(covariates$var_names)
	covariates$is_factor <- as.logical(covariates$is_factor)
	
	# Missing or duplicated covariate names would make the regression formula
	# ambiguous, so treat them as configuration errors.
	if (any(is.na(covariates$var_names) | !nzchar(covariates$var_names))){
		stop("All entries in covariates$var_names must be non-empty strings.")
	}
	if (anyDuplicated(covariates$var_names)){
		stop("Duplicate entries are not allowed in covariates$var_names.")
	}
	if (any(is.na(covariates$is_factor))){
		stop("All entries in covariates$is_factor must be TRUE or FALSE.")
	}
	
	# Each residualization covariate must be present in the phenotypic table.
	missing_covariates <- setdiff(covariates$var_names, names(pheno_df))
	if (length(missing_covariates) > 0){
		stop("Residualization covariates missing from the phenotype table: ",
				paste(missing_covariates, collapse = ", "))
	}
	
	invisible(TRUE)
}

# Escape column names for use in an R formula. This allows residualization to
# work even if phenotype column names contain spaces or other non-standard
# characters.
escapeFormulaName <- function(x){
	paste0("`", gsub("`", "\\`", x), "`")
}

# Residualize a single feature vector against a fixed covariate design matrix.
# This helper is used only when feature-specific missing values prevent the
# faster matrix-wide residualization path. It returns both residuals and any
# warning text so parallel workers do not lose diagnostics.
residualizeOneFeature <- function(feature_idx, data_mat, complete_idx,
		design_mat, view_name = "omics"){
	# Extract the values for this feature among samples with complete covariate
	# information. Feature-specific missing values are handled below.
	y <- data_mat[feature_idx, complete_idx]
	fit_idx <- !is.na(y)
	feature_name <- rownames(data_mat)[feature_idx]
	if (is.null(feature_name) || is.na(feature_name) || !nzchar(feature_name)){
		feature_name <- paste0("feature_", feature_idx)
	}
	
	# Initialize this feature's residual vector. Values remain NA for samples that
	# lack covariates or have missing omics measurements.
	feature_residuals <- rep(NA_real_, length(complete_idx))
	
	# Require enough observations to estimate the linearly independent columns of
	# the design matrix. Sparse features are left as NA and reported to the caller.
	if (sum(fit_idx) <= qr(design_mat[fit_idx, , drop = FALSE])$rank){
		warn_msg <- paste0("Feature '", feature_name, "' in view '", view_name,
				"' has too few non-missing samples for residualization; ",
				"setting residuals to NA.")
		return(list(residuals = feature_residuals, warning = warn_msg))
	}
	
	# lm.fit() is used instead of lm() because the design matrix has already been
	# constructed once and can be reused for every feature.
	fit <- lm.fit(design_mat[fit_idx, , drop = FALSE], y[fit_idx])
	feature_residuals[fit_idx] <- fit$residuals
	
	list(residuals = feature_residuals, warning = NA_character_)
}

# Residualize one omics matrix against the requested phenotype covariates.
# Rows are features and columns are samples. Column names are matched to the
# phenotype table through the view-specific barcode column.
residualizeOmicsMatrix <- function(data_mat, pheno_df, barcode_col, covariates,
		view_name = "omics"){
	# Return the original matrix unchanged when residualization is disabled.
	if (!use_residuals){
		return(data_mat)
	}
	
	# Ensure the phenotype table contains the barcode column needed for this view.
	if (!(barcode_col %in% names(pheno_df))){
		stop("Phenotype table is missing barcode column '", barcode_col,
				"' required for residualizing view '", view_name, "'.")
	}
	
	# Work with a numeric matrix while preserving feature and sample names.
	data_mat <- as.matrix(data_mat)
	storage.mode(data_mat) <- "double"
	
	# Remove omics samples that are not present in the phenotype table before
	# residualization. This prevents unmatched columns from causing a hard error
	# and ensures that all downstream calculations use the same sample set that
	# could be adjusted against the requested phenotype covariates.
	pheno_barcodes <- unique(pheno_df[[barcode_col]])
	matched_samples <- colnames(data_mat) %in% pheno_barcodes
	missing_pheno <- colnames(data_mat)[!matched_samples]
	if (length(missing_pheno) > 0){
		warning("Removing ", length(missing_pheno), " sample(s) from view '",
				view_name, "' because they are not present in phenotype barcode column '",
				barcode_col, "': ", paste(missing_pheno, collapse = ", "))
		data_mat <- data_mat[, matched_samples, drop = FALSE]
	}
	
	# Stop if no omics samples remain after dropping columns without phenotype
	# records; without matched samples, no residualization or downstream MOFA
	# analysis can be performed for this view.
	if (ncol(data_mat) == 0){
		stop("No samples remain in view '", view_name,
				"' after removing samples absent from the phenotype table.")
	}
	
	# Restrict the phenotype table to the retained samples from the omics matrix.
	# Each omics column is then matched back to exactly one phenotype row using the
	# barcode column. Base-R subsetting is used here to avoid
	# non-standard-evaluation issues with a barcode column selected by name.
	pheno_slice <- pheno_df[pheno_df[[barcode_col]] %in% colnames(data_mat),
			, drop = FALSE]
	pheno_slice <- pheno_slice[!duplicated(pheno_slice[[barcode_col]]),
			, drop = FALSE]
	pheno_slice <- pheno_slice[match(colnames(data_mat), pheno_slice[[barcode_col]]),
			, drop = FALSE]
	
	# Coerce user-declared categorical covariates to factors before constructing
	# the regression model matrix.
	factor_vars <- covariates$var_names[covariates$is_factor]
	for (var_name in factor_vars){
		pheno_slice[[var_name]] <- factor(pheno_slice[[var_name]])
	}
	
	# Build a formula containing only the covariates, for example ~ age + sex.
	# The model matrix created from this formula is reused for every feature.
	covariate_formula <- as.formula(paste("~",
					paste(escapeFormulaName(covariates$var_names), collapse = " + ")))
	model_frame <- model.frame(covariate_formula, data = pheno_slice,
			na.action = na.pass)
	
	# Only samples with complete covariate information can be residualized.
	complete_covariates <- complete.cases(model_frame)
	if (sum(complete_covariates) == 0){
		stop("No complete phenotype rows are available for residualizing view '",
				view_name, "'.")
	}
	
	# Drop unused factor levels after removing incomplete covariate rows, then make
	# sure each factor still has at least two observed levels.
	model_frame_complete <- droplevels(model_frame[complete_covariates,
					, drop = FALSE])
	for (var_name in factor_vars){
		if (var_name %in% names(model_frame_complete) &&
				length(levels(model_frame_complete[[var_name]])) < 2){
			stop("Categorical residualization covariate '", var_name,
					"' has fewer than two observed levels in view '", view_name, "'.")
		}
	}
	
	# Construct the design matrix for the complete-covariate samples. R handles
	# numeric covariates as continuous variables and factor covariates through the
	# usual contrast/dummy-variable coding.
	design_mat <- model.matrix(covariate_formula, data = model_frame_complete)
	complete_idx <- which(complete_covariates)
	
	# Initialize the output with NA values. Samples with missing covariates or
	# missing feature values remain NA because no residual can be estimated there.
	residualized_mat <- data_mat
	residualized_mat[,] <- NA_real_
	
	# Extract the complete-covariate part of the omics matrix once. If this matrix
	# has no missing values, all features can be residualized in a single lm.fit()
	# call. This vectorized path is usually faster than any feature-wise loop, even
	# a parallelized one.
	data_complete <- data_mat[, complete_idx, drop = FALSE]
	if (!anyNA(data_complete)){
		# lm.fit() accepts a response matrix with samples in rows and features in
		# columns. Its residual matrix therefore needs to be transposed back to the
		# original feature-by-sample orientation.
		if (nrow(design_mat) <= qr(design_mat)$rank){
			warning("View '", view_name,
					"' has too few complete samples for residualization; setting residuals to NA.")
		} else {
			message("Residualizing ", view_name,
					" with matrix-wide lm.fit() because no omics values are missing.")
			fit <- lm.fit(design_mat, t(data_complete))
			residualized_mat[, complete_idx] <- t(fit$residuals)
		}
	} else {
		# Missing omics values require feature-specific fits, because each feature may
		# use a different subset of samples. These independent fits can be executed in
		# parallel when requested.
		feature_indices <- seq_len(nrow(data_mat))
		
		if (use_parallel_residualization && length(feature_indices) > 1 &&
				n_residualization_cores > 1){
			# Limit the number of workers to the number of features to avoid launching
			# idle processes for small matrices. PSOCK clusters are used because they are
			# portable across Windows, macOS, and Linux.
			n_workers <- min(n_residualization_cores, length(feature_indices))
			message("Residualizing ", view_name, " feature-wise using ",
					n_workers, " parallel worker(s).")
			cl <- parallel::makeCluster(n_workers)
			on.exit(parallel::stopCluster(cl), add = TRUE)
			
			# Export the shared objects once per worker, instead of resending the full
			# omics matrix with every feature. This is important for large omics tables.
			parallel::clusterExport(
					cl,
					varlist = c("data_mat", "complete_idx", "design_mat",
							"view_name", "residualizeOneFeature"),
					envir = environment()
			)
			
			residual_list <- parallel::parLapply(
					cl,
					X = feature_indices,
					fun = function(feature_idx){
						residualizeOneFeature(feature_idx, data_mat, complete_idx,
								design_mat, view_name)
					}
			)
		} else {
			# Serial fallback used when parallelization is disabled, only one core is
			# available, or the matrix has only one feature.
			residual_list <- lapply(
					feature_indices,
					residualizeOneFeature,
					data_mat = data_mat,
					complete_idx = complete_idx,
					design_mat = design_mat,
					view_name = view_name
			)
		}
		
		# Emit warnings collected from serial or parallel feature-wise fits. Returning
		# warning text explicitly avoids losing worker diagnostics from parLapply().
		warn_msgs <- unique(na.omit(vapply(residual_list, `[[`, character(1),
								"warning")))
		for (warn_msg in warn_msgs){
			warning(warn_msg)
		}
		
		# Convert the list of residual vectors back into a feature-by-sample matrix
		# and store it in the complete-covariate sample columns.
		residualized_complete <- do.call(rbind,
				lapply(residual_list, `[[`, "residuals"))
		residualized_mat[, complete_idx] <- residualized_complete
	}
	
	message("Residualized ", view_name, " data using covariates: ",
			paste(covariates$var_names, collapse = ", "))
	residualized_mat
}

# Check the residualization configuration before reading any omics data views.
validateResidualizationInputs(pheno_d, covariates)

# Definitions for the optional omics data views.
omics_input_defs <- list(
		methylation = list(
				moniker = "me",
				file_var = "calerie_methy_fn",
				barcode_col = "blood_DNAm_barcode",
				reader = function(fn) readMethylation(fn)
		),
		proteomics = list(
				moniker = "pr",
				file_var = "calerie_prot_fn",
				barcode_col = "proteomics_barcode",
				reader = function(fn) readProteomics(fn)
		),
		metabolomics = list(
				moniker = "mb",
				file_var = "calerie_metab_data_fn",
				barcode_col = "metabolomics_barcode",
				reader = function(fn) readMetabolomics(fn)[["data"]]
		)
)

# Read only the omics data views whose file names are defined in resources.R.
loaded_views <- list()
barcode_cols <- character(0)
for (view_name in names(omics_input_defs)){
	view_def <- omics_input_defs[[view_name]]
	input_file <- getResourceFile(view_def$file_var)
	
	if (is.null(input_file)){
		next
	}
	
	message("Reading ", view_name, " data from: ", input_file)
	
	# Read the raw omics matrix from file. Residualization, if enabled, is applied
	# immediately after this step so every downstream calculation uses adjusted
	# feature values.
	view_data <- view_def$reader(input_file)
	view_data <- residualizeOmicsMatrix(
			view_data,
			pheno_d,
			view_def$barcode_col,
			covariates,
			view_name = view_name
	)
	
	loaded_views[[view_name]] <- view_data
	barcode_cols[view_name] <- view_def$barcode_col
}

if (length(loaded_views) == 0){
	stop("No omics data file names are defined in resources.R.")
}

# Remove features with zero variance across selected samples. If fewer than two
# samples are present, variance cannot be estimated, so the data are returned
# unchanged.
removeZeroVarianceFeatures <- function(data_mat){
	if (ncol(data_mat) < 2){
		return(data_mat)
	}
	
	row_vars <- apply(data_mat, 1, var, na.rm = TRUE)
	keep <- !is.na(row_vars) & row_vars != 0
	data_mat[keep, , drop = FALSE]
}

# Compute delta (timepoint - baseline) for a single data view
#
# Generate change scores for the target visit against baseline. The helper
# extracts samples with both a baseline measurement and a measurement at the
# requested timepoint, then computes per-feature differences.
#
# Columns can be renamed to the format <DEID>_<AL|CR>_<visit>.
computeDelta <- function(data_mat, pheno_df, barcode_col, timepoint,
		remove_zero_var = TRUE, rename_cols = TRUE){
	tp_label <- paste0("mo_", timepoint)
	
	# Build the pheno slice: keep only rows for visits 0, timepoint and where 
	# the barcode for this view is present in the data matrix's columns.
	ph_slice <- pheno_df %>%
			filter(
					as.character(fumo) %in% c("0", as.character(timepoint)),
					.data[[barcode_col]] %in% colnames(data_mat)
			) %>%
			# normalise visit codes to strings like 'mo_0', 'mo_12' etc.
			mutate(fumo = paste0("mo_", fumo)) %>%
			# collapse duplicates by DEID/CR/fumo keeping the first barcode
			group_by(DEID, CR, fumo) %>%
			summarise(barcode = first(.data[[barcode_col]]), .groups = "drop") %>%
			# pivot to wide so each visit becomes its own column (mo_0, mo_12, mo_24)
			pivot_wider(names_from = fumo, values_from = barcode)
	
	if (!("mo_0" %in% names(ph_slice))){
		ph_slice[["mo_0"]] <- NA_character_
	}
	if (!(tp_label %in% names(ph_slice))){
		ph_slice[[tp_label]] <- NA_character_
	}
	
	# Require baseline and timepoint measurements present
	ph_slice <- ph_slice %>% filter(!is.na(mo_0))
	ph_slice_tp <- ph_slice %>% filter(!is.na(.data[[tp_label]]))
	
	if (nrow(ph_slice_tp) == 0){
		return(list(delta = data_mat[, FALSE, drop = FALSE],
						slice = ph_slice_tp))
	}
	
	# Compute delta (timepoint - baseline)
	delta <- data_mat[, ph_slice_tp[[tp_label]], drop = FALSE] -
			data_mat[, ph_slice_tp[["mo_0"]], drop = FALSE]
	
	# Eemove features with zero variance across selected delta
	# samples (primarily useful for proteomics/metabolomics views)
	if (remove_zero_var){
		delta <- removeZeroVarianceFeatures(delta)
	}
	
	# Rename columns to a standardized format: <DEID>_<AL|CR>_<visit>
	if (rename_cols){
		name_map <- setNames(
				paste(ph_slice_tp$DEID,
						ifelse(ph_slice_tp$CR == 0, "AL", "CR"),
						paste0(timepoint, "m"), sep = "_"),
				ph_slice_tp[[tp_label]]
		)
		colnames(delta) <- name_map[colnames(delta)]
	}
	
	list(delta = delta, slice = ph_slice_tp)
}

# Extract target-visit data for a single data view without subtracting baseline.
extractTimepointData <- function(data_mat, pheno_df, barcode_col, timepoint,
		remove_zero_var = TRUE, rename_cols = TRUE){
	ph_slice_tp <- pheno_df %>%
			filter(
					as.character(fumo) == as.character(timepoint),
					.data[[barcode_col]] %in% colnames(data_mat)
			) %>%
			group_by(DEID, CR) %>%
			summarise(barcode = first(.data[[barcode_col]]), .groups = "drop")
	
	if (nrow(ph_slice_tp) == 0){
		return(list(data = data_mat[, FALSE, drop = FALSE],
						slice = ph_slice_tp))
	}
	
	data_tp <- data_mat[, ph_slice_tp$barcode, drop = FALSE]
	
	if (remove_zero_var){
		data_tp <- removeZeroVarianceFeatures(data_tp)
	}
	
	if (rename_cols){
		name_map <- setNames(
				paste(ph_slice_tp$DEID,
						ifelse(ph_slice_tp$CR == 0, "AL", "CR"),
						paste0(timepoint, "m"), sep = "_"),
				ph_slice_tp$barcode
		)
		colnames(data_tp) <- name_map[colnames(data_tp)]
	}
	
	list(data = data_tp, slice = ph_slice_tp)
}

# Prepare one omics view according to the use_change_scores setting.
prepareOmicsView <- function(data_mat, pheno_df, barcode_col, timepoint,
		use_change_scores){
	if (use_change_scores){
		return(computeDelta(data_mat, pheno_df, barcode_col, timepoint)$delta)
	}
	
	extractTimepointData(data_mat, pheno_df, barcode_col, timepoint)$data
}

# Apply the selected preprocessing to all available omics views.
data_views <- list()
for (view_name in names(loaded_views)){
	data_views[[view_name]] <- prepareOmicsView(
			loaded_views[[view_name]], pheno_d, barcode_cols[[view_name]],
			timepoint, use_change_scores
	)
}

# Drop any view that has no samples at the requested timepoint.
sample_counts <- vapply(data_views, ncol, integer(1))
if (any(sample_counts == 0)){
	warning(
			"Dropping omics views with no samples at the requested timepoint: ",
			paste(names(data_views)[sample_counts == 0], collapse = ", ")
	)
	data_views <- data_views[sample_counts > 0]
}

if (length(data_views) == 0){
	stop("No omics views have samples at the requested timepoint.")
}

# Extract the AL/CR label from sample IDs of the form <DEID>_<AL|CR>_<visit>.
# The treatment label is taken from the penultimate field, so this remains
# valid if DEID itself contains underscores.
getSampleGroups <- function(sample_ids){
	groups <- vapply(strsplit(sample_ids, split = "_", fixed = TRUE),
			function(parts){
				if (length(parts) < 3){
					return(NA_character_)
				}
				parts[length(parts) - 1]
			}, character(1))
	
	if (any(is.na(groups) | !nzchar(groups))){
		stop("Could not infer AL/CR labels from all sample IDs.")
	}
	
	names(groups) <- sample_ids
	groups
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

# Select highly variable features within each group, then take their union.
# If the union exceeds max_features, retain the features with the largest
# within-group variance across all groups. This is how MOFA recommends
# selection of highly variable features be handled when using groups.
getTopVariableFeatureIndexByGroup <- function(view_data, view_groups,
		max_features){
	# Ensure that the group vector matches the number of samples (columns)
	if (length(view_groups) != ncol(view_data)){
		stop("The group vector must contain one label per sample column.")
	}
	
	# Disallow missing or empty group labels
	if (any(is.na(view_groups) | !nzchar(view_groups))){
		stop("Missing group labels are not allowed for grouped HVF selection.")
	}
	
	# Cap max_features to the total number of available features (rows)
	max_features <- min(max_features, nrow(view_data))
	
	# Identify unique group labels (e.g., "AL", "CR")
	group_names <- unique(view_groups)
	
	# For each group:
	#   - subset the data to samples in that group
	#   - compute the top variable features within that subset
	# Combine all selected feature indices across groups
	group_feat_idx <- unlist(lapply(group_names, function(group_name){
						group_cols <- view_groups == group_name
						getTopVariableFeatureIndex(view_data[, group_cols, drop = FALSE],
								max_features)
					}), use.names = FALSE)
	# Take the union of features selected across all groups
	feat_idx <- unique(group_feat_idx)
	
	# If the total number of selected features does not exceed the limit,
	# return them directly
	if (length(feat_idx) <= max_features){
		return(feat_idx)
	}
	
	# For each group, compute the variance of each feature (row-wise)
	group_vars <- lapply(group_names, function(group_name){
				group_data <- view_data[, view_groups == group_name, drop = FALSE]
				if (ncol(group_data) < 2){
					return(rep(NA_real_, nrow(group_data)))
				}
				apply(group_data, 1, var, na.rm = TRUE)
			})
	
	# Combine per-group variance vectors into a matrix:
	# rows = features, columns = groups
	group_var_mat <- do.call(cbind, group_vars)
	group_var_mat[is.na(group_var_mat)] <- -Inf
	
	# For each feature (restricted to those already selected across groups),
	# take the maximum variance observed across all groups
	feature_scores <- apply(group_var_mat[feat_idx, , drop = FALSE], 1, max)
	# Rank features by their maximum group-specific variance and
	# select the top max_features
	feat_idx[order(feature_scores, decreasing = TRUE)[seq_len(max_features)]]
}

# MOFA recommends that data views have a comparable number of features, to avoid
# one view dominating factor fitting. We therefore reduce large views by keeping
# only the most variable features. A view may have at most max_factor times as
# many features as the smallest view, and never more than max_features.
#
# If use_groups is TRUE, highly variable features are selected separately
# within each group and then combined across groups before enforcing the final
# feature cap, per MOFA guidelines.
filterMostVariableFeatures <- function(data_views, max_factor=5, groups = NULL,
		use_groups = FALSE, max_features = 10000){
	# Get number of features (rows) for each omics view
	feature_counts <- vapply(data_views, nrow, integer(1))
	
	if (any(feature_counts == 0)){
		stop("At least one omics view has zero features after preprocessing.")
	}
	if (use_groups && is.null(groups)){
		stop("groups must be provided when use_groups is TRUE.")
	}
	if (use_groups && is.null(names(groups))){
		stop("groups must be a named vector when use_groups is TRUE.")
	}
	
	min_feat_num <- min(feature_counts)
	
	# Iterate over each omics view (data matrix)
	mapply(function(view_data, view_name){
				
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
				
				if (use_groups){
					# Map sample columns to their group labels
					view_groups <- groups[colnames(view_data)]
					# Ensure all samples have valid group assignments
					if (any(is.na(view_groups))){
						stop(paste0("Missing group labels for view '", view_name,
										"'."))
					}
					# Select top variable features using group-aware method
					feat_idx <- getTopVariableFeatureIndexByGroup(view_data,
							view_groups, vf_num)
				} else {
					# Select top variable features across all samples (no grouping)
					feat_idx <- getTopVariableFeatureIndex(view_data, vf_num)
				}
				
				view_data[feat_idx, , drop = FALSE]
			}, data_views, names(data_views), SIMPLIFY = FALSE)
}

# Set scaling factor controlling how many features to retain per view
# (relative to the smallest view) 
max_factor <- 5
# Set an absolute upper bound on the number of features per view - useful
# primarily when running MOFA with only the methylation data view, to
# avoid performance degradation due to hundreds of thousands of features
max_features <- 10000
# Collect the union of all sample IDs across views (column names),
# ensuring that group assignment covers all samples present
sample_ids_for_filter <- Reduce(union, lapply(data_views, colnames))
sample_groups_for_filter <- if (use_groups){
			getSampleGroups(sample_ids_for_filter)
		} else {
			NULL
		}

# Apply feature filtering to each view:
#   - selects most variable features
#   - uses group-aware selection if use_groups == TRUE
#   - enforces both relative (max_factor) and absolute (max_features) limits
data_views <- filterMostVariableFeatures(data_views, max_factor,
		groups = sample_groups_for_filter, use_groups = use_groups,
		max_features = max_features)

# Indicates if the data should be scaled. According to the MOFA paper, scaling
# is not needed. Scaling the data, however, may be helpful for some plots; for
# example, clusters may be more distinguishable in heatmaps when using z-scores.
scale_views <- FALSE

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

# Assign each sample to its AL/CR treatment label. These labels are passed to
# MOFA only when use_groups is TRUE.
grps <- if (use_groups){
			getSampleGroups(sample_ids)
		} else {
			NULL
		}

# Create MOFA object, comprising all available data matrices (aka, views).
if (use_groups){
	mdd <- create_mofa(data_views, groups = grps)
} else {
	mdd <- create_mofa(data_views)
}
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
model_opts$num_factors <- 50 # If we want to manually set the factor num

# Set training options.
train_opts <- get_default_training_options(mdd)
train_opts$drop_factor_threshold <- 0.001
train_opts$convergence_mode <- "slow"
# increase from default value of 1000, just in case more iters are required 
# for convergence
train_opts$maxiter <- 4000

# Prepare the MOFA object and add the metadata data frame, with one row per
# sample. The sample column is always required; the group column is added only
# when running MOFA in multi-group mode.
mo <- prepare_mofa(object = mdd, data_options = data_opts,
		model_options = model_opts, training_options = train_opts)

meta_d <- pheno_d %>%
		filter(as.character(fumo) == as.character(timepoint)) %>%
		mutate(sample_id = paste(DEID, ifelse(CR == 0, "AL", "CR"),
						paste0(timepoint, "m"), sep = "_")) %>%
		distinct(sample_id, .keep_all = TRUE)
rownames(meta_d) <- meta_d$sample_id

missing_meta <- setdiff(sample_ids, rownames(meta_d))
if (length(missing_meta) > 0){
	stop("Missing phenotypic metadata for samples: ",
			paste(missing_meta, collapse = ", "))
}

meta_d <- meta_d[sample_ids, , drop = FALSE]
meta_d <- data.frame(meta_d, sample = rownames(meta_d))
if (use_groups){
	meta_d$group <- grps[rownames(meta_d)]
}
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
group_tag <- ifelse(use_groups, "grouped", "ungrouped")
resid_tag <- ifelse(use_residuals, "_residualized", "")
views_tag <- paste(sapply(views_names(mo), function(vn){
					omics_input_defs[[vn]]$moniker}), 
		collapse="+")
res_folder <- paste(output_dir,
		paste0(format(Sys.time(), "%Y_%b_%d_%H_%M_%S"),
				"_visit=", timepoint, "mo_", score_tag, "_", group_tag,
				resid_tag, "_views=", views_tag),
		sep = "/")
res_file <- paste(res_folder, "model.hdf5", sep = "/")
log_file <- paste(res_folder, "mofa_out.txt", sep = "/")
dir.create(res_folder, recursive = TRUE)
sink(file = log_file)
mo.tr <- run_mofa(mo, res_file, use_basilisk = TRUE)
sink()

# Logging - write out the history of commands executed as well as the session
# info.
session_file <- paste(res_folder, "sessionInfo.txt", sep = "/")
writeLines(capture.output(sessionInfo()), session_file)
hist_file <- paste(res_folder, "MOFA_analysis_calerie.R", sep = "/")
savehistory(hist_file)

# Document what was done.
readme_file <- paste(res_folder, "README.txt", sep = "/")
score_desc <- if (use_change_scores){
			paste0("change scores (", timepoint, "m - baseline)")
		} else {
			paste0("source values from samples where fumo == ", timepoint)
		}

included_views <- paste(names(data_views), collapse = ", ")
data_notes <- character(0)
if (use_residuals){
	# Document the residualization step in the run-specific README so the saved
	# output folder records that MOFA was fit on adjusted omics matrices.
	data_notes <- c(data_notes,
			paste("Before downstream processing, each omics matrix was",
					"residualized using phenotype covariates:",
					paste(covariates$var_names, collapse = ", "), "."),
			paste("Residualization used matrix-wide lm.fit() when possible and",
					"feature-wise fitting when missing omics values required it;",
					"parallel feature-wise residualization was",
					ifelse(use_parallel_residualization, "enabled", "disabled"), "."))
}
if ("methylation" %in% names(data_views)){
	data_notes <- c(data_notes, "Methylation data are reported as M-values.")
}
if ("metabolomics" %in% names(data_views)){
	metab_note <- if (use_change_scores){
				paste("Metabolomics data in CALERIE are reported in Metabolon",
						"median-centered units; accordingly, we log2 transform them",
						"prior to change-score calculation.")
			} else {
				paste("Metabolomics data in CALERIE are reported in Metabolon",
						"median-centered units; accordingly, we log2 transform them",
						"prior to analysis.")
			}
	data_notes <- c(data_notes, metab_note)
}
if ("proteomics" %in% names(data_views)){
	data_notes <- c(data_notes,
			"Proteomics data, from Somascan, are already log2 transformed.")
}

group_desc <- if (use_groups){
			paste("Samples were split into two MOFA groups: subjects receiving",
					"the caloric restriction diet (CR group) and subjects on the",
					"alternative diet (AL group).")
		} else {
			paste("Samples were not split into MOFA groups; all samples were",
					"analyzed together. AL/CR labels remain available in the sample",
					"IDs and phenotype metadata.")
		}

readme_text <- paste(c(
				paste("MOFA was run on", score_desc, "."),
				paste("Included omics views:", included_views, "."),
				data_notes,
				"Each included omics data table becomes a separate MOFA view.",
				group_desc
		), collapse = " ")

sink(file = readme_file)
writeLines(strwrap(readme_text, width = 80))
sink()
