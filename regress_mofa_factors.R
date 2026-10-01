# *****************************************************************************
# Regression of MOFA factors against user-selected metadata variables, adjusting
# for selected sample-level metadata covariates.
#
# This script:
#   1. Loads regression inputs from either:
#        - a trained MOFA model in model.hdf5, or
#        - persisted sample-metadata and factor-value files (.rds or .xlsx).
#   2. Retrieves/reads sample metadata.
#   3. Retrieves/reads predicted sample factor values.
#   4. Converts user-specified metadata columns to categorical variables.
#   5. Fits one regression model per MOFA factor with each Y as the outcome:
#        Y ~ Factor_i + X_covariates
#      using:
#        - logistic regression if Y is categorical, or
#        - linear regression if Y is continuous.
#   6. Fits one linear model per MOFA factor with the factor as the outcome:
#        Factor_i ~ Y + X_covariates
#      and focuses the report on Y as the explanatory variable.
#   7. Repeats both model directions for every Y value and every covariate set
#      listed in regression_parameters.R.
#   8. Returns fitted model objects in a nested list called res.
#   9. Generates a separate PDF report for each model direction and Y/X
#      combination.
#  10. Generates one consolidated Excel workbook containing model-level results,
#      coefficient-level results, and a data dictionary for all fitted models.
#  11. Writes the same three summary tables to an RDS file, as a named list with
#      one data matrix per Excel worksheet.
#
# Edit regression_parameters.R, then execute this script in one of two ways:
#
#   source("C:/my_project/code/regress_mofa_factors.R",
#          local = TRUE, chdir = TRUE)
#
# or, from a command prompt:
#
#   Rscript.exe -e "source('C:/my_project/code/regress_mofa_factors.R',
#                   local = TRUE, chdir = TRUE)"
#
# In both cases, chdir = TRUE temporarily sets the working directory to the
# directory containing this script. This allows the companion parameter file
# to be sourced directly from the same directory.
# *****************************************************************************

# Load MOFA2 functions, including load_model(), samples_metadata(), and get_factors().
library(MOFA2)

# -----------------------------------------------------------------------------
# User-specified inputs
# -----------------------------------------------------------------------------

# Load the user-configurable regression settings from the companion parameter
# file. The supported invocation methods above make this script's directory the
# current working directory before this statement is evaluated.
source("regression_parameters.R", local = TRUE)

# Normalize an optional path setting. NULL and an empty string both mean that
# the corresponding input has not been specified.
normalizeOptionalPath <- function(x, setting_name) {
	if (is.null(x)) {
		return(NULL)
	}
	if (!is.character(x) || length(x) != 1L || is.na(x)) {
		stop(setting_name, " must be NULL or a single non-missing file path.")
	}
	x <- trimws(x)
	if (!nzchar(x)) {
		return(NULL)
	}
	x
}

# Validate the configured metadata column that contains sample identifiers.
normalizeSampleIdColumn <- function(x) {
	if (!is.character(x) || length(x) != 1L || is.na(x)) {
		stop("sample_id_col must be a single non-missing character string.")
	}
	x <- trimws(x)
	if (!nzchar(x)) {
		stop("sample_id_col must be a non-empty character string.")
	}
	x
}

# Validate this setting once when the script starts.
sample_id_col <- normalizeSampleIdColumn(sample_id_col)

# Determine which of the two supported regression input modes is configured.
# Exactly one mode must be selected:
#   1. mofa_res_fname, or
#   2. sample_metadata_fname + factor_values_fname.
resolveRegressionDataSource <- function(mofa_res_fname, sample_metadata_fname,
		factor_values_fname) {
	mofa_res_fname <- normalizeOptionalPath(mofa_res_fname, "mofa_res_fname")
	sample_metadata_fname <- normalizeOptionalPath(sample_metadata_fname,
			"sample_metadata_fname")
	factor_values_fname <- normalizeOptionalPath(factor_values_fname,
			"factor_values_fname")

	model_mode <- !is.null(mofa_res_fname)
	metadata_file_set <- !is.null(sample_metadata_fname)
	factor_file_set <- !is.null(factor_values_fname)
	file_mode <- metadata_file_set || factor_file_set

	if (model_mode && file_mode) {
		stop(
				"Specify exactly one regression input mode: either mofa_res_fname, ",
				"or both sample_metadata_fname and factor_values_fname, but not both."
		)
	}

	if (!model_mode && !file_mode) {
		stop(
				"No regression input source was specified. Set mofa_res_fname, or set ",
				"both sample_metadata_fname and factor_values_fname."
		)
	}

	if (file_mode && !(metadata_file_set && factor_file_set)) {
		missing_setting <- if (metadata_file_set) "factor_values_fname" else
				"sample_metadata_fname"
		stop(
				"File-based regression input requires both sample_metadata_fname and ",
				"factor_values_fname. Missing setting: ", missing_setting, "."
		)
	}

	if (model_mode) {
		if (!file.exists(mofa_res_fname)) {
			stop("Cannot find MOFA result file: ", mofa_res_fname)
		}
		return(list(
				mode = "mofa_model",
				mofa_res_fname = mofa_res_fname,
				sample_metadata_fname = NULL,
				factor_values_fname = NULL,
				output_root = dirname(normalizePath(mofa_res_fname, winslash = "/",
						mustWork = TRUE))
		))
	}

	for (fname in c(sample_metadata_fname, factor_values_fname)) {
		if (!file.exists(fname)) {
			stop("Cannot find regression input file: ", fname)
		}
		ext <- tolower(tools::file_ext(fname))
		if (!(ext %in% c("rds", "xlsx"))) {
			stop(
					"Unsupported regression input file type for '", fname,
					"'. Expected an .rds or .xlsx file."
			)
		}
	}

	list(
			mode = "matrix_files",
			mofa_res_fname = NULL,
			sample_metadata_fname = sample_metadata_fname,
			factor_values_fname = factor_values_fname,
			# The supported analysis and factor-estimation workflows write these files
			# into the same result directory. Use the factor-values directory as the
			# regression output root.
			output_root = dirname(normalizePath(factor_values_fname, winslash = "/",
					mustWork = TRUE))
	)
}

regression_data_source <- resolveRegressionDataSource(
		mofa_res_fname = mofa_res_fname,
		sample_metadata_fname = sample_metadata_fname,
		factor_values_fname = factor_values_fname
)
regression_output_root <- regression_data_source$output_root

# Build a compact, filesystem-safe label for text inserted into output filenames.
# The original variable names used in the regression model are not changed.
makeFilenameLabel <- function(values, empty_label = "none") {
	# Use a clear label when no values are supplied.
	if (is.null(values) || length(values) == 0L) {
		return(empty_label)
	}
	
	# Replace characters that are inconvenient or invalid in filenames.
	safe_values <- gsub("[^A-Za-z0-9_.-]+", "-", values)
	safe_values <- gsub("^-+|-+$", "", safe_values)
	
	# If sanitization produced an empty string, fall back to the requested label.
	safe_values[safe_values == ""] <- empty_label
	
	# Join multiple values with '+', as in covar=female+race.
	paste(safe_values, collapse = "+")
}

# Normalize one covariate-set entry from X. character(0), NULL, and c("") are
# all treated as an unadjusted model; mixed empty and non-empty names are rejected
# because they are usually accidental.
normalizeCovariateSet <- function(covariates, set_label) {
	if (is.null(covariates)) {
		return(character(0))
	}
	if (!is.character(covariates)) {
		stop("Each entry in X must be a character vector of metadata column names. ",
				"Problematic entry: ", set_label, ".")
	}
	if (length(covariates) == 0L) {
		return(character(0))
	}
	if (any(is.na(covariates))) {
		stop("Covariate set ", set_label, " contains NA. Use character(0) ",
				"or c(\"\") for an unadjusted model.")
	}
	covariates <- trimws(covariates)
	blank <- !nzchar(covariates)
	if (all(blank)) {
		return(character(0))
	}
	if (any(blank)) {
		stop("Covariate set ", set_label, " mixes empty and non-empty column ",
				"names. Use character(0) or c(\"\") for an unadjusted model.")
	}
	if (anyDuplicated(covariates)) {
		stop("Covariate set ", set_label, " contains duplicated column names: ",
				paste(unique(covariates[duplicated(covariates)]), collapse = ", "))
	}
	covariates
}

# Validate the batch-level Y and X settings from regression_parameters.R before
# any model fitting starts. runRegression() still validates each scalar Y/X
# combination immediately before fitting.
normalizeRegressionBatchInputs <- function(Y, X) {
	if (!is.character(Y) || length(Y) == 0L || any(is.na(Y)) ||
			any(!nzchar(trimws(Y)))) {
		stop("Y must be a non-empty character vector of metadata column names.")
	}
	Y <- trimws(Y)
	if (anyDuplicated(Y)) {
		stop("Y contains duplicated column names: ",
				paste(unique(Y[duplicated(Y)]), collapse = ", "))
	}
	outcome_labels <- vapply(Y, makeFilenameLabel, character(1),
			empty_label = "outcome")
	if (anyDuplicated(outcome_labels)) {
		stop("Y contains column names that would generate duplicate report ",
				"filenames after sanitization: ",
				paste(unique(outcome_labels[duplicated(outcome_labels)]),
						collapse = ", "))
	}
	
	if (is.null(X)) {
		X <- list(character(0))
	} else if (is.character(X)) {
		# Preserve the previous single-covariate-vector style as a convenience.
		X <- list(X)
	} else if (!is.list(X) || length(X) == 0L) {
		stop("X must be a non-empty list of character vectors. Use ",
				"list(character(0)) for an unadjusted model.")
	}
	
	X <- lapply(seq_along(X), function(i) {
				normalizeCovariateSet(X[[i]], paste0("X[[", i, "]]"))
			})
	
	covariate_labels <- vapply(X, makeFilenameLabel, character(1),
			empty_label = "none")
	if (anyDuplicated(covariate_labels)) {
		stop("X contains covariate sets that would generate duplicate report ",
				"filenames: ", paste(unique(covariate_labels[duplicated(covariate_labels)]),
						collapse = ", "))
	}
	
	list(Y = Y, X = X, outcome_labels = outcome_labels,
			covariate_labels = covariate_labels)
}

makeRegressionAnalysisLabel <- function(outcome, covariates) {
	paste0("Y=", makeFilenameLabel(outcome, empty_label = "outcome"),
			"__covar=", makeFilenameLabel(covariates, empty_label = "none"))
}

# Build the PDF report path for one regression direction. Reports are written
# under separate mode-specific folders, with filenames keyed by Y and covariates
# so batch runs do not overwrite each other.
makeRegressionReportFilename <- function(outcome, covariates, factors_as_Y) {
	file.path(
			regression_output_root, "regression_analysis",
			ifelse(factors_as_Y, "factors_as_Y", "outcome_as_Y"),
			paste0(
					"mofa_factor_regression_report_Y=",
					makeFilenameLabel(outcome, empty_label = "outcome"),
					"_covar=",
					makeFilenameLabel(covariates, empty_label = "none"),
					".pdf"
			)
	)
}

# Normalize the batch configuration once so the run block can iterate over
# scalar Y values and one covariate vector at a time.
regression_inputs <- normalizeRegressionBatchInputs(Y, X)
Y_values <- regression_inputs$Y
X_sets <- regression_inputs$X

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

quoteFormulaName <- function(x) {
	# lm() and glm() formulas need special handling when variable names contain
	# spaces, punctuation, or other non-syntactic characters. Wrapping names in
	# backticks allows formulas such as `clinical response` ~ `F` + `age at baseline`.
	if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
		stop("Formula variable names must be non-empty strings.")
	}
	
	# Backticks inside variable names would break the quoted formula syntax.
	if (grepl("`", x, fixed = TRUE)) {
		stop("Backticks in variable names are not supported: ", x)
	}
	
	# Return the variable name wrapped in backticks for safe use in a formula.
	paste0("`", x, "`")
}

makeRegressionFormula <- function(outcome, covariates, factors_as_Y = FALSE) {
	# Build the model formula used for every MOFA factor. The column F is a
	# temporary column created later for the current MOFA factor.
	if (isTRUE(factors_as_Y)) {
		# Reverse-direction analysis: the current MOFA factor is the dependent
		# variable, while Y and the selected covariates are explanatory variables.
		rhs <- paste(
				vapply(c(outcome, covariates), quoteFormulaName, character(1)),
				collapse = " + "
		)
		
		# Return a formula equivalent to:
		#   F ~ Y + covariate_1 + covariate_2 + ...
		return(stats::as.formula(paste(quoteFormulaName("F"), "~", rhs)))
	}
	
	# In this analysis mode, the metadata variable Y is the dependent variable, and
	# the current MOFA factor plus the selected covariates are explanatory variables.
	rhs <- paste(
			vapply(c("F", covariates), quoteFormulaName, character(1)),
			collapse = " + "
	)
	
	# Return a formula equivalent to:
	#   Y ~ F + covariate_1 + covariate_2 + ...
	stats::as.formula(paste(quoteFormulaName(outcome), "~", rhs))
}

looksLikeFactorName <- function(x) {
	# Detect common MOFA factor labels such as "Factor1", "Factor 1", or "Factor_1".
	# This is used only to infer whether factors are stored in rows or columns.
	grepl("^Factor[ _]*[0-9]+$", x, ignore.case = TRUE)
}

normalizeFactorMatrix <- function(x) {
	# Convert the MOFA factor object into a standard numeric matrix.
	# Expected final layout: rows = samples, columns = MOFA factors.
	mat <- as.matrix(x)
	
	# Determine whether row names or column names appear to represent MOFA factors.
	rownames_are_factors <- !is.null(rownames(mat)) &&
			any(looksLikeFactorName(rownames(mat)))
	colnames_are_factors <- !is.null(colnames(mat)) &&
			any(looksLikeFactorName(colnames(mat)))
	
	# MOFA2 usually returns a sample x factor matrix. This transpose handles the
	# alternative factor x sample layout if encountered.
	if (rownames_are_factors && !colnames_are_factors) {
		mat <- t(mat)
	}
	
	# After normalization, row names should identify samples. These are needed to
	# align factor values with the sample metadata table.
	if (is.null(rownames(mat)) || any(is.na(rownames(mat))) ||
			any(!nzchar(rownames(mat)))) {
		stop("Could not identify sample names in the MOFA factor matrix.")
	}
	
	# If factor column names are missing, create simple names. The returned res list
	# will still be named "Factor 1", "Factor 2", etc. later in the script.
	if (is.null(colnames(mat)) || any(is.na(colnames(mat))) ||
			any(!nzchar(colnames(mat)))) {
		colnames(mat) <- paste("Factor", seq_len(ncol(mat)))
	}

	# Factor values must be numeric. Excel readers can occasionally return a
	# character matrix even when all data cells contain numbers, so safely coerce
	# numeric-looking values while rejecting genuine non-numeric entries.
	if (!is.numeric(mat)) {
		mat_chr <- matrix(as.character(mat), nrow = nrow(mat), ncol = ncol(mat),
				dimnames = dimnames(mat))
		mat_num <- suppressWarnings(matrix(
				as.numeric(mat_chr),
				nrow = nrow(mat_chr),
				ncol = ncol(mat_chr),
				dimnames = dimnames(mat_chr)
		))
		nonempty_original <- !is.na(mat_chr) & nzchar(trimws(mat_chr))
		if (any(is.na(mat_num) & nonempty_original)) {
			stop("The factor values matrix contains non-numeric values.")
		}
		mat <- mat_num
	}
	
	# Return the standardized sample x factor matrix.
	mat
}

getFactorDataFrame <- function(mofa_obj) {
	# Extract predicted factor values from the MOFA object for all factors and all groups.
	# as.data.frame = FALSE preserves matrix/list structure, which is easier to align
	# robustly across one-group and multi-group MOFA models.
	factor_obj <- MOFA2::get_factors(
			mofa_obj,
			factors = "all",
			groups = "all",
			as.data.frame = FALSE
	)
	
	# Multi-group MOFA models may return a list of matrices, one matrix per group.
	if (is.list(factor_obj) && !is.data.frame(factor_obj)) {
		group_names <- names(factor_obj)
		
		# Create fallback group names if the list has no usable names.
		if (is.null(group_names) || any(!nzchar(group_names))) {
			group_names <- paste0("group", seq_along(factor_obj))
		}
		
		# Convert each group's factor matrix into a data frame with sample identifiers.
		parts <- vector("list", length(factor_obj))
		for (i in seq_along(factor_obj)) {
			mat <- normalizeFactorMatrix(factor_obj[[i]])
			
			# Store sample ID and MOFA group next to the factor columns.
			parts[[i]] <- cbind(
					data.frame(
							sample_id = rownames(mat),
							.mofa_group = group_names[i],
							stringsAsFactors = FALSE,
							check.names = FALSE
					),
					as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
			)
		}
		
		# Stack all groups into one data frame.
		factor_d <- do.call(rbind, parts)
		rownames(factor_d) <- NULL
		return(factor_d)
	}
	
	# Single-group MOFA models usually return one matrix rather than a list.
	mat <- normalizeFactorMatrix(factor_obj)
	
	# Convert the sample x factor matrix into a data frame and retain sample names.
	factor_d <- cbind(
			data.frame(
					sample_id = rownames(mat),
					.mofa_group = NA_character_,
					stringsAsFactors = FALSE,
					check.names = FALSE
			),
			as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
	)
	rownames(factor_d) <- NULL
	
	# Return a data frame with columns: sample_id, .mofa_group, Factor 1, Factor 2, ...
	factor_d
}

normalizeSampleMetadata <- function(meta_d, sample_id_col, require_sample_id = FALSE) {
	# Stop if there is no metadata. The requested regression requires Y and X columns.
	if (is.null(meta_d) || nrow(meta_d) == 0L) {
		stop("No sample metadata were found.")
	}

	# Convert to a plain data frame while preserving original column names.
	meta_d <- as.data.frame(meta_d, stringsAsFactors = FALSE, check.names = FALSE)

	# File-based input must contain the configured sample-ID column. For MOFA-object
	# input, create that column from MOFA2's conventional 'sample' column or valid
	# row names if the configured column is not already present.
	if (!(sample_id_col %in% colnames(meta_d))) {
		if (isTRUE(require_sample_id)) {
			stop(
					"The sample metadata file must contain the sample-ID column '",
					sample_id_col, "', as specified by sample_id_col."
			)
		}
		if ("sample" %in% colnames(meta_d)) {
			meta_d[[sample_id_col]] <- meta_d$sample
		} else if (!is.null(rownames(meta_d)) && all(!is.na(rownames(meta_d))) &&
				all(nzchar(rownames(meta_d)))) {
			meta_d[[sample_id_col]] <- rownames(meta_d)
		} else {
			stop(
					"The sample metadata does not contain the configured sample-ID column '",
					sample_id_col, "', and no 'sample' column or valid row names are available."
			)
		}
	}

	# Store sample identifiers as character strings to make matching unambiguous.
	meta_d[[sample_id_col]] <- as.character(meta_d[[sample_id_col]])
	if (any(is.na(meta_d[[sample_id_col]])) || any(!nzchar(meta_d[[sample_id_col]]))) {
		stop(
				"The sample metadata contains missing or empty values in sample-ID column '",
				sample_id_col, "'."
		)
	}

	# Return metadata with the configured sample-ID column guaranteed to exist.
	meta_d
}

getSampleMetadata <- function(mofa_obj, sample_id_col) {
	# Retrieve the sample metadata stored in the MOFA object and normalize its
	# configured sample identifier column for the common regression code path.
	normalizeSampleMetadata(
			MOFA2::samples_metadata(mofa_obj),
			sample_id_col = sample_id_col,
			require_sample_id = FALSE
	)
}

# Read the first worksheet from a regression input workbook. Sample-metadata
# workbooks store every field, including the sample-ID field, as an ordinary
# named worksheet column. Factor-values workbooks instead store sample IDs in
# the first worksheet column because those IDs represent matrix row names.
readRegressionExcelFile <- function(fname, object_label,
		first_column_contains_row_names = FALSE) {
	if (!requireNamespace("xlsx", quietly = TRUE)) {
		stop(
				"Reading ", object_label, " from an .xlsx file requires the R package ",
				"'xlsx'. Install it or use the corresponding .rds export."
		)
	}

	d <- xlsx::read.xlsx(
			fname,
			sheetIndex = 1L,
			header = TRUE,
			stringsAsFactors = FALSE,
			check.names = FALSE
	)

	if (nrow(d) == 0L || ncol(d) == 0L) {
		stop("The ", object_label, " Excel file is empty or has an unexpected layout: ",
				fname)
	}

	# Sample metadata must retain every worksheet column and its header. In
	# particular, do not discard the first metadata column or convert it to row
	# names; it may be the configured sample-ID or categorical column.
	if (!isTRUE(first_column_contains_row_names)) {
		return(d)
	}

	# A factor-values workbook must contain a sample-ID column followed by at
	# least one MOFA factor column.
	if (ncol(d) < 2L) {
		stop("The ", object_label,
				" Excel file must contain sample row identifiers and factor values: ",
				fname)
	}

	row_ids <- as.character(d[[1L]])
	if (any(is.na(row_ids)) || any(!nzchar(row_ids))) {
		stop("The first column of ", object_label,
				" must contain non-empty row identifiers: ", fname)
	}
	if (anyDuplicated(row_ids)) {
		stop("The first column of ", object_label,
				" contains duplicated row identifiers: ", fname)
	}

	d <- d[, -1L, drop = FALSE]
	rownames(d) <- row_ids
	d
}

readRegressionInputFile <- function(fname, object_label,
		first_column_contains_row_names = FALSE) {
	ext <- tolower(tools::file_ext(fname))
	if (ext == "rds") {
		return(readRDS(fname))
	}
	if (ext == "xlsx") {
		return(readRegressionExcelFile(
				fname,
				object_label,
				first_column_contains_row_names = first_column_contains_row_names
		))
	}
	stop("Unsupported file type for ", object_label, ": ", fname)
}

readSampleMetadataFile <- function(fname, sample_id_col) {
	# All worksheet columns are metadata fields and must retain their headers.
	meta_d <- readRegressionInputFile(
			fname,
			"sample metadata",
			first_column_contains_row_names = FALSE
	)
	normalizeSampleMetadata(
			meta_d,
			sample_id_col = sample_id_col,
			require_sample_id = TRUE
	)
}

readFactorValuesFile <- function(fname) {
	# The first worksheet column contains the sample IDs written from matrix row
	# names; restore those IDs as row names before normalizing the factor matrix.
	factor_values <- readRegressionInputFile(
			fname,
			"factor values",
			first_column_contains_row_names = TRUE
	)
	normalizeFactorMatrix(factor_values)
}

getFactorDataFrameFromValues <- function(factor_values) {
	# Convert an already loaded samples x factors object into the same internal
	# layout produced by getFactorDataFrame() for a single-group MOFA model.
	mat <- normalizeFactorMatrix(factor_values)
	factor_d <- cbind(
			data.frame(
					sample_id = rownames(mat),
					.mofa_group = NA_character_,
					stringsAsFactors = FALSE,
					check.names = FALSE
			),
			as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
	)
	rownames(factor_d) <- NULL
	factor_d
}

loadRegressionData <- function(data_source, sample_id_col) {
	if (identical(data_source$mode, "mofa_model")) {
		mo.tr <- MOFA2::load_model(data_source$mofa_res_fname)
		return(list(
				mo.tr = mo.tr,
				sample_metadata = NULL,
				factor_values = NULL
		))
	}

	meta_d <- readSampleMetadataFile(data_source$sample_metadata_fname, sample_id_col)
	factor_values <- readFactorValuesFile(data_source$factor_values_fname)

	# Abort before model fitting if any factor-value sample is absent from the
	# configured sample-ID column in the sample metadata.
	missing_sample_ids <- setdiff(rownames(factor_values), meta_d[[sample_id_col]])
	if (length(missing_sample_ids) > 0L) {
		stop(
				"The factor values matrix contains sample IDs that are not present in ",
				"the '", sample_id_col, "' column of the sample metadata matrix: ",
				paste(missing_sample_ids, collapse = ", ")
		)
	}

	list(
			mo.tr = NULL,
			sample_metadata = meta_d,
			factor_values = factor_values
	)
}

validateRegressionInputs <- function(X, Y, categorical, sample_id_col) {
	# Treat NULL covariates as no covariates.
	if (is.null(X)) {
		X <- character(0)
	}
	
	# Treat NULL categorical variables as none explicitly categorical.
	if (is.null(categorical)) {
		categorical <- character(0)
	}
	
	# Validate X before using it in formulas or metadata subsetting.
	if (!is.character(X)) {
		stop("X must be a character vector of metadata column names.")
	}
	if (any(is.na(X)) || any(!nzchar(X))) {
		stop("X contains missing or empty column names.")
	}
	if (anyDuplicated(X)) {
		stop(
				"X contains duplicated column names: ",
				paste(unique(X[duplicated(X)]), collapse = ", ")
		)
	}
	
	# Validate Y as a single metadata column name.
	if (!is.character(Y) || length(Y) != 1L || is.na(Y) || !nzchar(Y)) {
		stop("Y must be a single non-empty character string naming the dependent variable.")
	}
	
	# Validate categorical before using it to coerce metadata columns.
	if (!is.character(categorical)) {
		stop("categorical must be a character vector of metadata column names.")
	}
	if (any(is.na(categorical)) || any(!nzchar(categorical))) {
		stop("categorical contains missing or empty column names.")
	}
	if (anyDuplicated(categorical)) {
		stop(
				"categorical contains duplicated column names: ",
				paste(unique(categorical[duplicated(categorical)]), collapse = ", ")
		)
	}
	
	# The dependent variable cannot also be included as an adjustment covariate.
	if (Y %in% X) {
		stop("Do not include Y in X; the dependent variable cannot also be a covariate.")
	}
	
	# F is a temporary column used for the current MOFA factor, so reserve that name.
	if (Y == "F" || "F" %in% X || "F" %in% categorical) {
		stop("'F' is reserved for the current MOFA factor and cannot be used as Y, X, or categorical.")
	}
	
	# The configured sample-ID column is used only for alignment; MOFA-object
	# metadata may also retain MOFA2's conventional 'sample' column. Neither
	# identifier should be modeled.
	reserved_sample_cols <- unique(c(sample_id_col, "sample"))
	if (Y %in% reserved_sample_cols || any(X %in% reserved_sample_cols)) {
		stop(
				"The configured sample-ID column ('", sample_id_col,
				"') and the MOFA 'sample' column are reserved for sample alignment and ",
				"cannot be used as Y or as covariates in X."
		)
	}
	
	# Return normalized input vectors so the main function can use them directly.
	list(X = X, Y = Y, categorical = categorical)
}

coerceModelVariables <- function(meta_d, model_vars, categorical) {
	# Convert only variables used in the model. Extra names in categorical are
	# validated elsewhere but do not need to be modified here.
	relevant_categorical <- intersect(categorical, model_vars)
	
	# User-specified categorical variables are explicitly converted to factors.
	# This is especially important for numeric encodings such as 0/1, 1/2/3, or
	# batch numbers that should not be treated as continuous linear predictors.
	for (v in relevant_categorical) {
		meta_d[[v]] <- factor(meta_d[[v]])
	}
	
	# Character and logical model variables are also categorical in R's formula
	# interface. Converting them to factors makes that behavior explicit.
	for (v in model_vars) {
		if (is.character(meta_d[[v]]) || is.logical(meta_d[[v]])) {
			meta_d[[v]] <- factor(meta_d[[v]])
		}
	}
	
	# Return metadata with model variables prepared for lm()/glm().
	meta_d
}

isCategoricalOutcome <- function(meta_d, outcome, categorical) {
	# After coerceModelVariables(), user-specified categorical outcomes and
	# character/logical outcomes are factors. Factor outcomes trigger logistic
	# regression, while numeric outcomes trigger linear regression.
	y <- meta_d[[outcome]]
	
	if (is.factor(y)) {
		return(TRUE)
	}
	
	if (is.numeric(y)) {
		# Numeric variables are treated as continuous unless the user listed the
		# outcome in categorical. Warn about the common binary-numeric case because
		# users often intend logistic regression for 0/1 outcomes.
		vals <- sort(unique(stats::na.omit(y)))
		if (length(vals) == 2L && !(outcome %in% categorical)) {
			warning(
					"Outcome '", outcome,
					"' has exactly two numeric values but is not listed in categorical; ",
					"treating it as continuous and using linear regression. Add '", outcome,
					"' to categorical to use logistic regression."
			)
		}
		return(FALSE)
	}
	
	# Other classes, such as Date or list columns, are not handled by this script.
	stop(
			"Outcome '", outcome, "' has unsupported class: ",
			paste(class(y), collapse = ", "),
			". Use a numeric continuous outcome or a categorical/factor outcome."
	)
}

prepareBinaryOutcome <- function(y, outcome, factor_label) {
	# Binary logistic regression requires a two-level outcome. The outcome should
	# already be a factor, but factor() also handles a logical or character vector
	# if one slips through.
	y <- droplevels(factor(y))
	
	# Binomial glm() with a factor response models the second level relative to the
	# first level. Use relevel() before calling runRegression() if you need a
	# specific reference level.
	if (nlevels(y) != 2L) {
		stop(
				"Model for ", factor_label,
				" cannot be fit because outcome '", outcome,
				"' does not have exactly two non-missing classes after removing missing values."
		)
	}
	
	# Return the checked binary factor outcome.
	y
}

prepareContinuousOutcome <- function(y, outcome, factor_label) {
	# Linear regression requires a numeric dependent variable.
	if (!is.numeric(y)) {
		stop(
				"Model for ", factor_label,
				" cannot be fit because outcome '", outcome,
				"' is not numeric after removing missing values."
		)
	}
	
	# Return the checked numeric outcome.
	y
}

dropUnusedFactorLevels <- function(d) {
	# After complete-case filtering, some factor levels may no longer be present.
	# Dropping unused levels keeps nlevels() checks accurate and avoids stale levels.
	factor_cols <- vapply(d, is.factor, logical(1))
	d[factor_cols] <- lapply(d[factor_cols], droplevels)
	d
}

assertCategoricalCovariatesUsable <- function(d, covariates, factor_label) {
	# lm()/glm() cannot construct contrasts for a factor covariate that has fewer
	# than two observed levels in the complete-case data for the current model.
	for (v in covariates) {
		if (is.factor(d[[v]]) && nlevels(d[[v]]) < 2L) {
			stop(
					"Model for ", factor_label,
					" cannot be fit because categorical covariate '", v,
					"' has fewer than two observed levels after removing missing values."
			)
		}
	}
}

# -----------------------------------------------------------------------------
# Main analysis function
# -----------------------------------------------------------------------------

# Fit one regression model per MOFA factor.
#
# Args:
#   mo.tr: Optional trained MOFA object, usually returned by MOFA2::load_model().
#   sample_metadata: Optional sample metadata data frame/matrix containing the
#     sample-ID column specified by sample_id_col.
#   factor_values: Optional samples x factors matrix with sample IDs as row names.
#     Supply sample_metadata and factor_values together when mo.tr is NULL.
#   sample_id_col: Name of the metadata column containing sample identifiers.
#   X: Character vector of metadata column names to use as covariates.
#   Y: Single character string naming the metadata column to use as the outcome.
#   categorical: Character vector of metadata column names to treat as categorical.
#
# Returns:
#   A named list of fitted model objects. Names are "Factor 1", "Factor 2", etc.
#   Each entry is a glm object for categorical/binary Y in the Y-as-outcome mode,
#   or an lm object for continuous Y or for factors_as_Y mode.
runRegression <- function(mo.tr = NULL, sample_metadata = NULL,
		factor_values = NULL, sample_id_col = "sample_id", X = character(0), Y,
		categorical = character(0), factors_as_Y = FALSE) {
	# Validate and normalize user inputs.
	sample_id_col <- normalizeSampleIdColumn(sample_id_col)
	inputs <- validateRegressionInputs(X, Y, categorical, sample_id_col)
	X <- inputs$X
	Y <- inputs$Y
	categorical <- inputs$categorical
	
	# factors_as_Y must be a single TRUE/FALSE value so model direction is unambiguous.
	if (!is.logical(factors_as_Y) || length(factors_as_Y) != 1L || is.na(factors_as_Y)) {
		stop("factors_as_Y must be a single TRUE/FALSE value.")
	}
	
	# Obtain sample metadata and factor values from exactly one source. Direct calls
	# to runRegression() can therefore use either the historical MOFA-object argument
	# or the new pair of already loaded sample_metadata/factor_values objects.
	using_mofa_obj <- !is.null(mo.tr)
	using_direct_data <- !is.null(sample_metadata) || !is.null(factor_values)
	if (using_mofa_obj && using_direct_data) {
		stop(
				"runRegression() received both mo.tr and direct sample_metadata/factor_values. ",
				"Use exactly one input mode."
		)
	}
	if (!using_mofa_obj && !using_direct_data) {
		stop(
				"runRegression() requires either mo.tr or both sample_metadata and ",
				"factor_values."
		)
	}
	if (using_direct_data && (is.null(sample_metadata) || is.null(factor_values))) {
		stop(
				"Direct runRegression() input requires both sample_metadata and factor_values."
		)
	}

	if (using_mofa_obj) {
		meta_d <- getSampleMetadata(mo.tr, sample_id_col)
		factor_d <- getFactorDataFrame(mo.tr)
	} else {
		meta_d <- normalizeSampleMetadata(
				sample_metadata,
				sample_id_col = sample_id_col,
				require_sample_id = TRUE
		)
		factor_d <- getFactorDataFrameFromValues(factor_values)
	}
	
	# Confirm that metadata contains sample IDs, outcome Y, and all requested covariates.
	required_meta_cols <- c(sample_id_col, Y, X)
	missing_meta_cols <- setdiff(required_meta_cols, colnames(meta_d))
	if (length(missing_meta_cols) > 0L) {
		stop("Missing metadata columns required for the model: ", paste(missing_meta_cols, collapse = ", "))
	}
	
	# Confirm that every user-specified categorical variable exists in the metadata.
	# Some categorical variables may not be used in the current model; that is allowed.
	missing_categorical_cols <- setdiff(categorical, colnames(meta_d))
	if (length(missing_categorical_cols) > 0L) {
		stop("Missing metadata columns listed in categorical: ", paste(missing_categorical_cols, collapse = ", "))
	}
	
	# Prepare model variables. User-specified categorical variables are converted to
	# factors; character and logical model variables are also converted to factors.
	model_vars <- c(Y, X)
	meta_d <- coerceModelVariables(meta_d, model_vars, categorical)
	
	# Decide which type of regression to run. In the Y-as-outcome mode, this decision
	# is based on Y. In factors_as_Y mode, the dependent variable is always the
	# continuous MOFA factor value F, so linear regression is always used.
	if (isTRUE(factors_as_Y)) {
		model_type <- "linear"
	} else {
		outcome_is_categorical <- isCategoricalOutcome(meta_d, Y, categorical)
		model_type <- if (outcome_is_categorical) "logistic" else "linear"
	}
	
	# Sample IDs must be unique so that each factor row maps to exactly one metadata row.
	if (anyDuplicated(meta_d[[sample_id_col]])) {
		stop("Sample identifiers are duplicated in metadata column '", sample_id_col, "'.")
	}
	if (anyDuplicated(factor_d$sample_id)) {
		stop("Sample identifiers are duplicated in the MOFA factor matrix.")
	}
	
	# Match factor rows to metadata rows using sample identifiers. Every sample in
	# the factor-value matrix must be represented in metadata; otherwise abort and
	# list the missing IDs rather than silently dropping samples.
	meta_idx <- match(factor_d$sample_id, meta_d[[sample_id_col]])
	if (anyNA(meta_idx)) {
		missing_sample_ids <- factor_d$sample_id[is.na(meta_idx)]
		stop(
				"The factor values matrix contains sample IDs that are not present in ",
				"the '", sample_id_col, "' column of the sample metadata matrix: ",
				paste(missing_sample_ids, collapse = ", ")
		)
	}
	
	# Reorder metadata to match the factor-value rows. Extra metadata rows are
	# allowed because they are simply not part of the factor matrix being analyzed.
	meta_d <- meta_d[meta_idx, , drop = FALSE]
	
	# Factor columns are all columns except sample identifiers and optional group labels.
	factor_cols <- setdiff(colnames(factor_d), c("sample_id", ".mofa_group"))
	if (length(factor_cols) == 0L) {
		stop("No MOFA factors were found in the factor values data.")
	}
	
	# Construct the model formula once. The column F will be overwritten in each loop
	# with the values of the current MOFA factor.
	fmla <- makeRegressionFormula(Y, X, factors_as_Y = factors_as_Y)
	
	# Pre-allocate the output list and name entries as requested.
	res <- vector("list", length(factor_cols))
	names(res) <- paste("Factor", seq_along(factor_cols))
	
	# Fit one model for each MOFA factor.
	for (i in seq_along(factor_cols)) {
		factor_label <- names(res)[i]
		
		# Start with the metadata columns needed for this model.
		this_d <- meta_d[, c(sample_id_col, Y, X), drop = FALSE]
		
		# Add the current MOFA factor as a column named F. Depending on factors_as_Y,
		# the formula is either Y ~ F + X or F ~ Y + X.
		this_d$F <- factor_d[[factor_cols[i]]]
		
		# Remove rows with missing outcome, missing factor value, or missing covariates.
		model_cols <- c(Y, "F", X)
		complete <- stats::complete.cases(this_d[, model_cols, drop = FALSE])
		this_d <- this_d[complete, , drop = FALSE]
		
		# lm()/glm() need at least two complete samples. In practice, more samples may
		# be required if the model has many covariates or categorical levels.
		if (nrow(this_d) < 2L) {
			stop(
					"Model for ", factor_label,
					" cannot be fit because fewer than two complete samples remain."
			)
		}
		
		# Drop unused levels after complete-case filtering and check categorical predictors.
		this_d <- dropUnusedFactorLevels(this_d)
		
		if (isTRUE(factors_as_Y)) {
			# In factors_as_Y mode, Y and X are all independent variables. Any of them
			# that are categorical must have at least two observed levels.
			assertCategoricalCovariatesUsable(this_d, c(Y, X), factor_label)
			
			# The dependent variable is the current MOFA factor value, which must remain
			# numeric for the linear model F ~ Y + X.
			this_d$F <- prepareContinuousOutcome(this_d$F, "F", factor_label)
			
			# Fit linear regression for this factor and store the lm object.
			res[[i]] <- stats::lm(fmla, data = this_d)
		} else if (model_type == "logistic") {
			# In the Y-as-outcome mode, X variables are independent covariates.
			assertCategoricalCovariatesUsable(this_d, X, factor_label)
			
			# Re-check the outcome after dropping incomplete rows; it must contain exactly
			# two observed classes for binomial logistic regression.
			this_d[[Y]] <- prepareBinaryOutcome(this_d[[Y]], Y, factor_label)
			
			# Fit logistic regression for this factor and store the glm object.
			res[[i]] <- stats::glm(fmla, data = this_d, family = stats::binomial())
		} else {
			# In the Y-as-outcome mode, X variables are independent covariates.
			assertCategoricalCovariatesUsable(this_d, X, factor_label)
			
			# Re-check the outcome after dropping incomplete rows; it must remain numeric
			# for linear regression.
			this_d[[Y]] <- prepareContinuousOutcome(this_d[[Y]], Y, factor_label)
			
			# Fit linear regression for this factor and store the lm object.
			res[[i]] <- stats::lm(fmla, data = this_d)
		}
	}
	
	# Store lightweight metadata on the returned list. The list entries themselves
	# are still the fitted glm/lm objects requested above.
	attr(res, "outcome") <- Y
	attr(res, "covariates") <- X
	attr(res, "categorical") <- categorical
	attr(res, "model_type") <- model_type
	attr(res, "formula") <- fmla
	attr(res, "factors_as_Y") <- factors_as_Y
	attr(res, "association_term") <- if (isTRUE(factors_as_Y)) Y else "F"
	
	# Return the named list of fitted model objects.
	res
}

# -----------------------------------------------------------------------------
# Excel summary helper functions
# -----------------------------------------------------------------------------

# Replace the temporary model variable F in a fitted formula with the displayed
# MOFA factor name. The resulting formula is used as the human-readable model ID.
makeModelId <- function(model, factor_name) {
	fmla <- stats::formula(model)
	parts <- as.character(fmla)
	if (length(parts) != 3L) {
		return(paste(deparse(fmla), collapse = " "))
	}
	lhs <- parts[2L]
	rhs <- parts[3L]
	factor_term <- quoteFormulaName(factor_name)
	lhs <- gsub("(?<![[:alnum:]_.`])F(?![[:alnum:]_.`])", factor_term, lhs,
			perl = TRUE)
	rhs <- gsub("(?<![[:alnum:]_.`])F(?![[:alnum:]_.`])", factor_term, rhs,
			perl = TRUE)
	paste(lhs, "~", rhs)
}

# Return a model-wide p-value against the intercept-only model. For linear models
# this is the usual F-test p-value; for logistic models it is a likelihood-ratio
# chi-square test based on the null and residual deviances.
extractOverallModelPValue <- function(model) {
	model_type <- getModelType(model)
	if (model_type == "linear") {
		fs <- summary(model)$fstatistic
		if (is.null(fs) || anyNA(fs)) return(NA_real_)
		return(stats::pf(fs[[1L]], fs[[2L]], fs[[3L]], lower.tail = FALSE))
	}
	if (model_type == "logistic") {
		df <- model$df.null - model$df.residual
		if (!is.finite(df) || df <= 0L) return(NA_real_)
		return(stats::pchisq(model$null.deviance - model$deviance,
				df = df, lower.tail = FALSE))
	}
	NA_real_
}

# Build the one-row-per-model table used by the Models worksheet.
buildExcelModelRows <- function(res, outcome, covariates, factors_as_Y) {
	summary_d <- summarizeRegressionModels(res)
	association_term <- attr(res, "association_term")
	rows <- lapply(seq_along(res), function(i) {
		model <- res[[i]]
		factor_name <- names(res)[i]
		model_type <- getModelType(model)
		model_id <- makeModelId(model, factor_name)
		mf <- stats::model.frame(model)
		response <- if (isTRUE(factors_as_Y)) factor_name else names(mf)[1L]
		observed <- stats::model.response(mf)
		log_lik <- tryCatch(as.numeric(stats::logLik(model)), error = function(e) NA_real_)
		residual_se <- if (model_type == "linear") unname(summary(model)$sigma) else NA_real_
		events <- non_events <- NA_integer_
		if (model_type == "logistic") {
			y01 <- as.integer(observed) - 1L
			events <- sum(y01 == 1L, na.rm = TRUE)
			non_events <- sum(y01 == 0L, na.rm = TRUE)
		}
		assoc <- extractAssociationCoefficient(model, association_term)
		or <- or_lo <- or_hi <- NA_real_
		if (model_type == "logistic" && is.finite(assoc["estimate"]) &&
				is.finite(assoc["se"])) {
			or <- exp(assoc["estimate"])
			or_lo <- exp(assoc["estimate"] - stats::qnorm(0.975) * assoc["se"])
			or_hi <- exp(assoc["estimate"] + stats::qnorm(0.975) * assoc["se"])
		}
		data.frame(
				model_id = model_id,
				model_direction = ifelse(factors_as_Y, "factor_as_outcome", "metadata_as_outcome"),
				model_type = model_type,
				response = response,
				mofa_factor = factor_name,
				metadata_variable = outcome,
				adjustment_covariates = ifelse(length(covariates) == 0L, "", paste(covariates, collapse = "; ")),
				primary_term = ifelse(factors_as_Y, outcome, factor_name),
				n = stats::nobs(model),
				events = events,
				non_events = non_events,
				primary_estimate = unname(assoc["estimate"]),
				primary_se = unname(assoc["se"]),
				primary_statistic = unname(assoc["statistic"]),
				primary_p_value = unname(assoc["p_value"]),
				primary_ci_lower = unname(assoc["estimate"] - stats::qnorm(0.975) * assoc["se"]),
				primary_ci_upper = unname(assoc["estimate"] + stats::qnorm(0.975) * assoc["se"]),
				odds_ratio = or,
				odds_ratio_ci_lower = or_lo,
				odds_ratio_ci_upper = or_hi,
				model_p_value = extractOverallModelPValue(model),
				log_likelihood = log_lik,
				aic = summary_d$aic[i],
				bic = summary_d$bic[i],
				r_squared = summary_d$r_squared[i],
				adj_r_squared = summary_d$adj_r_squared[i],
				residual_se = residual_se,
				rmse = summary_d$rmse[i],
				mae = summary_d$mae[i],
				mcfadden_pseudo_r2 = summary_d$pseudo_r2[i],
				auc = summary_d$auc[i],
				accuracy_0.5 = summary_d$accuracy[i],
				brier_score = summary_d$brier[i],
				stringsAsFactors = FALSE,
				check.names = FALSE
		)
	})
	do.call(rbind, rows)
}

# Build the coefficient-level table. Each coefficient gets its own row, including
# the intercept and dummy coefficients created for categorical predictors.
buildExcelCoefficientRows <- function(res) {
	rows <- lapply(seq_along(res), function(i) {
		model <- res[[i]]
		factor_name <- names(res)[i]
		model_type <- getModelType(model)
		coef_tab <- summary(model)$coefficients
		z <- stats::qnorm(0.975)
		d <- data.frame(
				model_id = makeModelId(model, factor_name),
				mofa_factor = factor_name,
				term = ifelse(rownames(coef_tab) == "F", factor_name, rownames(coef_tab)),
				estimate = unname(coef_tab[, "Estimate"]),
				se = unname(coef_tab[, "Std. Error"]),
				statistic = unname(coef_tab[, 3L]),
				p_value = unname(coef_tab[, 4L]),
				ci_lower = unname(coef_tab[, "Estimate"] - z * coef_tab[, "Std. Error"]),
				ci_upper = unname(coef_tab[, "Estimate"] + z * coef_tab[, "Std. Error"]),
				odds_ratio = NA_real_,
				odds_ratio_ci_lower = NA_real_,
				odds_ratio_ci_upper = NA_real_,
				stringsAsFactors = FALSE,
				check.names = FALSE
		)
		if (model_type == "logistic") {
			d$odds_ratio <- exp(d$estimate)
			d$odds_ratio_ci_lower <- exp(d$ci_lower)
			d$odds_ratio_ci_upper <- exp(d$ci_upper)
		}
		d
	})
	do.call(rbind, rows)
}

makeRegressionDictionary <- function() {
	data.frame(
			column = c("model_id", "model_direction", "model_type", "response", "mofa_factor",
					"metadata_variable", "adjustment_covariates", "primary_term", "n", "events",
					"non_events", "primary_estimate", "primary_se", "primary_statistic",
					"primary_p_value", "primary_ci_lower", "primary_ci_upper", "odds_ratio",
					"odds_ratio_ci_lower", "odds_ratio_ci_upper", "model_p_value",
					"log_likelihood", "aic", "bic", "r_squared", "adj_r_squared", "residual_se",
					"rmse", "mae", "mcfadden_pseudo_r2", "auc", "accuracy_0.5", "brier_score"),
			description = c(
					"Actual fitted regression formula, with the MOFA factor name substituted for temporary variable F.",
					"Whether metadata or the MOFA factor is the dependent variable.",
					"linear or logistic.", "Dependent variable in the displayed model formula.",
					"MOFA factor used in this model.", "Metadata variable selected through Y.",
					"Adjustment covariates selected through X; semicolon separated.",
					"Explanatory variable of primary interest for this workflow.",
					"Number of observations used to fit the model.",
					"Number of modeled positive-class observations; logistic models only.",
					"Number of modeled reference-class observations; logistic models only.",
					"Coefficient for the primary term; NA for a multi-level categorical primary term.",
					"Standard error of the primary coefficient.", "t statistic for lm or z statistic for glm.",
					"Primary-term p-value; omnibus drop1 test when appropriate.",
					"Lower Wald 95% confidence limit for the primary coefficient.",
					"Upper Wald 95% confidence limit for the primary coefficient.",
					"exp(primary_estimate); logistic models only.",
					"Lower Wald 95% confidence limit for the primary odds ratio.",
					"Upper Wald 95% confidence limit for the primary odds ratio.",
					"Overall model p-value versus an intercept-only model.", "Model log-likelihood.",
					"Akaike information criterion.", "Bayesian information criterion.",
					"Ordinary R-squared; linear models only.", "Adjusted R-squared; linear models only.",
					"Residual standard error; linear models only.", "Root mean squared residual; linear models only.",
					"Mean absolute residual; linear models only.",
					"McFadden pseudo-R-squared; logistic models only.", "ROC area under the curve; logistic models only.",
					"Classification accuracy at probability threshold 0.5; logistic models only.",
					"Brier score; logistic models only."),
			stringsAsFactors = FALSE,
			check.names = FALSE
	)
}

# Order the Models worksheet so nested/expanded versions of the same model are
# adjacent. Within each regression direction and metadata variable, models are
# grouped by MOFA factor and then ordered by increasing number of adjustment
# covariates. Covariate sets of the same size retain their order in X.
sortExcelModelRows <- function(models_d, Y_values, X_sets) {
	if (nrow(models_d) == 0L) {
		return(models_d)
	}

	# Reconstruct the display labels used in adjustment_covariates and associate
	# each configured covariate set with its original position in X.
	covariate_labels <- vapply(X_sets, function(x) {
		if (length(x) == 0L) "" else paste(x, collapse = "; ")
	}, character(1))
	covariate_counts <- lengths(X_sets)
	covariate_set_order <- match(models_d$adjustment_covariates, covariate_labels)
	covariate_count <- covariate_counts[covariate_set_order]

	# Preserve the configured Y order and the natural MOFA-factor order already
	# produced by runRegression(), rather than sorting factor names alphabetically
	# (which would place Factor 10 before Factor 2).
	outcome_order <- match(models_d$metadata_variable, Y_values)
	factor_levels <- unique(models_d$mofa_factor)
	factor_order <- match(models_d$mofa_factor, factor_levels)
	direction_order <- match(models_d$model_direction,
			c("metadata_as_outcome", "factor_as_outcome"))

	ord <- order(
			direction_order,
			outcome_order,
			factor_order,
			covariate_count,
			covariate_set_order,
			na.last = TRUE
	)

	models_d[ord, , drop = FALSE]
}

writeRegressionExcelWorkbook <- function(models_d, coefficients_d, xlsx_file) {
	if (!requireNamespace("xlsx", quietly = TRUE)) {
		stop("Writing the regression Excel summary requires the R package 'xlsx'.")
	}
	out_dir <- dirname(xlsx_file)
	if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
	wb <- xlsx::createWorkbook()
	add_sheet <- function(name, d) {
		sheet <- xlsx::createSheet(wb, sheetName = name)
		xlsx::addDataFrame(d, sheet, row.names = FALSE)
		xlsx::autoSizeColumn(sheet, colIndex = seq_len(ncol(d)))
	}
	add_sheet("Models", models_d)
	add_sheet("Coefficients", coefficients_d)
	add_sheet("Dictionary", makeRegressionDictionary())
	xlsx::saveWorkbook(wb, xlsx_file)
	invisible(xlsx_file)
}

# -----------------------------------------------------------------------------
# PDF report helper functions
# -----------------------------------------------------------------------------

safeLog10P <- function(p) {
	# Convert p-values to -log10(p) for plotting while avoiding Inf when p == 0.
	ifelse(is.na(p), NA_real_, -log10(pmax(p, .Machine$double.xmin)))
}

safeText <- function(x, fallback = "not available") {
	# Convert missing scalar values into readable report text.
	if (length(x) == 0L || is.na(x) || !is.character(x)) {
		return(fallback)
	}
	as.character(x)
}

fmt <- function(x, digits = 3) {
	# Format numbers consistently in tables and narrative text.
	ifelse(is.na(x) | !is.finite(x), NA_character_, formatC(x, digits = digits, format = "fg"))
}

useGgplotReportPlots <- function() {
	# Use ggplot2 for higher-quality report figures when the package is installed.
	# The script intentionally does not require ggplot2, so base R plots remain the
	# fallback on systems where ggplot2 is unavailable.
	requireNamespace("ggplot2", quietly = TRUE)
}

makeFactorPlotLevels <- function(factors) {
	# Preserve the current row order of factors in ggplot2 axes.
	#
	# The same factor labels can appear more than once in long-format plotting
	# data frames, for example after stacking AUC, accuracy, pseudo-R2, and Brier
	# score into one table. The levels argument to factor() must be unique, so we
	# use unique(factors) to retain first-seen order without duplicated levels.
	factor(factors, levels = unique(factors))
}

getBestModelName <- function(summary_d) {
	# Identify the factor with the smallest nominal factor-term p-value.
	# If no p-values are available, use the first factor as a deterministic fallback.
	if (all(is.na(summary_d$factor_p_value))) {
		return(summary_d$factor[1])
	}
	summary_d$factor[which.min(summary_d$factor_p_value)]
}

getBestModel <- function(res, summary_d) {
	# Retrieve the fitted model corresponding to the best factor name.
	best_name <- getBestModelName(summary_d)
	res[[best_name]]
}

getModelType <- function(model) {
	# Determine the fitted model type from the fitted object itself rather than
	# assuming that every entry in res has the same type.
	if (inherits(model, "glm") && identical(model$family$family, "binomial")) {
		return("logistic")
	}
	if (inherits(model, "lm")) {
		return("linear")
	}
	"unknown"
}

aucFromScores <- function(y01, score) {
	# Compute the ROC AUC using the rank-sum/Mann-Whitney formulation. This avoids
	# a dependency on packages such as pROC while handling tied predictions via ranks.
	keep <- !is.na(y01) & !is.na(score)
	y01 <- y01[keep]
	score <- score[keep]
	
	n_pos <- sum(y01 == 1)
	n_neg <- sum(y01 == 0)
	if (n_pos == 0L || n_neg == 0L) {
		return(NA_real_)
	}
	
	r <- rank(score, ties.method = "average")
	(sum(r[y01 == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

extractTermPValue <- function(model, term) {
	# Obtain an omnibus p-value for a model term using drop1(). This is especially
	# useful when the term is categorical and expands to multiple dummy coefficients.
	test_type <- if (inherits(model, "glm")) "Chisq" else "F"
	drop_tab <- try(stats::drop1(model, test = test_type), silent = TRUE)
	if (inherits(drop_tab, "try-error") || !(term %in% rownames(drop_tab))) {
		return(NA_real_)
	}
	
	p_cols <- grep("^Pr\\(", colnames(drop_tab), value = TRUE)
	if (length(p_cols) == 0L) {
		return(NA_real_)
	}
	
	unname(drop_tab[term, p_cols[1]])
}

extractAssociationCoefficient <- function(model, term = "F") {
	# Extract the coefficient, standard error, test statistic, and p-value for the
	# explanatory term of interest. In the Y-as-outcome mode, this is F. In factors_as_Y
	# mode this is Y. If the term expands to multiple coefficients, such as a
	# multi-level factor, the coefficient columns are set to NA and the p-value is
	# taken from the omnibus drop1() test for the term.
	sm <- summary(model)
	coef_tab <- sm$coefficients
	coef_names <- rownames(coef_tab)
	quoted_term <- quoteFormulaName(term)
	
	# Numeric terms usually appear as an exact coefficient row. Binary factors often
	# appear as one row whose name begins with the term name plus the non-reference level.
	exact_matches <- coef_names[coef_names %in% c(term, quoted_term)]
	prefix_matches <- coef_names[startsWith(coef_names, term) | startsWith(coef_names, quoted_term)]
	candidates <- unique(c(exact_matches, prefix_matches))
	
	if (length(candidates) == 1L) {
		row_name <- candidates[1]
		p_value <- extractTermPValue(model, term)
		if (is.na(p_value)) {
			p_value <- unname(coef_tab[row_name, 4])
		}
		
		# For lm(), the statistic column is a t statistic; for glm(), it is a z statistic.
		return(c(
						estimate = unname(coef_tab[row_name, "Estimate"]),
						se = unname(coef_tab[row_name, "Std. Error"]),
						statistic = unname(coef_tab[row_name, 3]),
						p_value = p_value
				))
	}
	
	c(
			estimate = NA_real_,
			se = NA_real_,
			statistic = NA_real_,
			p_value = extractTermPValue(model, term)
	)
}

summarizeRegressionModels <- function(res) {
	# Build one summary row per fitted model. The columns include both common
	# metrics and metrics that apply only to logistic or linear regression.
	association_term <- attr(res, "association_term")
	if (is.null(association_term) || length(association_term) != 1L || is.na(association_term)) {
		association_term <- "F"
	}
	
	rows <- lapply(names(res), function(factor_name) {
				model <- res[[factor_name]]
				model_type <- getModelType(model)
				factor_coef <- extractAssociationCoefficient(model, term = association_term)
				
				# Initialize model-specific metrics as missing values.
				r_squared <- adj_r_squared <- rmse <- mae <- NA_real_
				pseudo_r2 <- auc <- accuracy <- brier <- NA_real_
				
				# Compute linear-regression diagnostics when the model is an lm object.
				if (model_type == "linear") {
					sm <- summary(model)
					observed <- stats::model.response(stats::model.frame(model))
					fitted <- stats::fitted(model)
					residuals <- stats::residuals(model)
					r_squared <- unname(sm$r.squared)
					adj_r_squared <- unname(sm$adj.r.squared)
					rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
					mae <- mean(abs(residuals), na.rm = TRUE)
				}
				
				# Compute logistic-regression diagnostics when the model is binomial glm.
				if (model_type == "logistic") {
					observed <- stats::model.response(stats::model.frame(model))
					y01 <- as.integer(observed) - 1L
					prob <- stats::fitted(model)
					pred_class <- as.integer(prob >= 0.5)
					pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
					auc <- aucFromScores(y01, prob)
					accuracy <- mean(pred_class == y01, na.rm = TRUE)
					brier <- mean((prob - y01)^2, na.rm = TRUE)
				}
				
				data.frame(
						factor = factor_name,
						model_type = model_type,
						n = stats::nobs(model),
						aic = stats::AIC(model),
						bic = stats::BIC(model),
						factor_estimate = factor_coef["estimate"],
						factor_se = factor_coef["se"],
						factor_statistic = factor_coef["statistic"],
						factor_p_value = factor_coef["p_value"],
						factor_ci_lower = factor_coef["estimate"] - stats::qnorm(0.975) * factor_coef["se"],
						factor_ci_upper = factor_coef["estimate"] + stats::qnorm(0.975) * factor_coef["se"],
						r_squared = r_squared,
						adj_r_squared = adj_r_squared,
						rmse = rmse,
						mae = mae,
						pseudo_r2 = pseudo_r2,
						auc = auc,
						accuracy = accuracy,
						brier = brier,
						stringsAsFactors = FALSE,
						check.names = FALSE
				)
			})
	
	# Stack the per-model rows into one data frame.
	summary_d <- do.call(rbind, rows)
	rownames(summary_d) <- NULL
	summary_d
}

plotReportTextPage <- function(title, lines, cex_title = 1.25, cex_text = 0.85) {
	# Draw a simple text page inside the PDF. This is used for the cover page and
	# explanatory pages that describe how to interpret the model diagnostics.
	graphics::plot.new()
	graphics::par(mar = c(1, 1, 1, 1))
	graphics::text(0.02, 0.96, title, adj = c(0, 1), cex = cex_title, font = 2)
	
	if (length(lines) > 0L) {
		y <- 0.90
		line_height <- 0.035
		for (line in lines) {
			graphics::text(0.02, y, line, adj = c(0, 1), cex = cex_text)
			y <- y - line_height
			if (y < 0.04) {
				break
			}
		}
	}
}

plotSummaryTablePage <- function(summary_d, model_type, association_label = "F", rows_per_page = 28L) {
	# Render compact numeric tables. Rows are sorted by the factor-term p-value so
	# the strongest nominal associations are listed first. Additional PDF pages are
	# created automatically when there are more factors than fit on one page.
	if (model_type == "logistic") {
		display_d <- summary_d[, c("factor", "n", "factor_estimate", "factor_p_value",
						"pseudo_r2", "auc", "accuracy", "brier", "aic"), drop = FALSE]
		colnames(display_d) <- c("Factor", "N", paste0(association_label, " log-OR"), paste0(association_label, " p"), "McFadden R2", "AUC", "Accuracy", "Brier", "AIC")
	} else {
		display_d <- summary_d[, c("factor", "n", "factor_estimate", "factor_p_value",
						"r_squared", "adj_r_squared", "rmse", "mae", "aic"), drop = FALSE]
		colnames(display_d) <- c("Factor", "N", paste0(association_label, " beta"), paste0(association_label, " p"), "R2", "Adj R2", "RMSE", "MAE", "AIC")
	}
	
	# Sort by nominal p-value, with missing p-values placed last.
	p_col <- paste0(association_label, " p")
	sort_idx <- order(display_d[[p_col]], na.last = TRUE)
	display_d <- display_d[sort_idx, , drop = FALSE]
	
	# Format numeric columns for readability after sorting.
	for (j in seq_along(display_d)) {
		if (is.numeric(display_d[[j]])) {
			display_d[[j]] <- fmt(display_d[[j]], digits = 3)
		}
	}
	
	# Split the table into as many pages as needed; no factor rows are omitted.
	n_pages <- ceiling(nrow(display_d) / rows_per_page)
	for (page_idx in seq_len(n_pages)) {
		row_start <- (page_idx - 1L) * rows_per_page + 1L
		row_end <- min(page_idx * rows_per_page, nrow(display_d))
		page_d <- display_d[row_start:row_end, , drop = FALSE]
		
		graphics::plot.new()
		graphics::par(mar = c(1, 1, 1, 1))
		graphics::text(
				0.02,
				0.96,
				paste0("Model performance summary (page ", page_idx, " of ", n_pages, ")"),
				adj = c(0, 1),
				cex = 1.20,
				font = 2
		)
		graphics::text(
				0.02,
				0.92,
				paste0("Rows are sorted by the nominal p-value for ", association_label, "; smallest p-value first."),
				adj = c(0, 1),
				cex = 0.8
		)
		
		# Draw the table using monospaced text for reliable alignment in a base PDF.
		table_lines <- capture.output(print(page_d, row.names = FALSE, right = FALSE))
		y <- 0.86
		for (line in table_lines) {
			graphics::text(0.02, y, line, adj = c(0, 1), cex = 0.75, family = "mono")
			y <- y - 0.028
		}
	}
}

plotFactorEffectPage <- function(summary_d, model_type, association_label = "F") {
	# Plot the estimated coefficient with 95% confidence intervals across factors.
	estimate_label <- if (model_type == "logistic") {
				paste0(association_label, " coefficient: log odds ratio")
			} else {
				paste0(association_label, " coefficient: beta")
			}
	
	if (all(is.na(summary_d$factor_estimate)) || all(is.na(summary_d$factor_ci_lower)) ||
			all(is.na(summary_d$factor_ci_upper))) {
		plotReportTextPage(
				title = "Association estimates not shown",
				lines = c(
						paste0("No single coefficient estimate is available for ", association_label, " across the fitted models."),
						"This usually occurs when the explanatory variable of interest is a categorical variable with multiple levels,",
						"which is represented by multiple regression coefficients rather than one scalar beta.",
						"The p-value summary pages still report the omnibus test for the variable."
				)
		)
		return(invisible(NULL))
	}
	
	if (useGgplotReportPlots()) {
		plot_d <- summary_d
		plot_d$factor <- makeFactorPlotLevels(plot_d$factor)
		
		p <- ggplot2::ggplot(
						plot_d,
						ggplot2::aes(x = factor, y = factor_estimate, ymin = factor_ci_lower, ymax = factor_ci_upper)
				) +
				ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
				ggplot2::geom_pointrange() +
				ggplot2::labs(
						title = paste0("Association estimate for ", association_label, " across MOFA factors"),
						subtitle = "Intervals are approximate Wald 95% confidence intervals.",
						x = "MOFA factor",
						y = estimate_label
				) +
				ggplot2::theme_bw() +
				ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
		
		print(p)
		return(invisible(NULL))
	}
	
	x <- seq_len(nrow(summary_d))
	y <- summary_d$factor_estimate
	lo <- summary_d$factor_ci_lower
	hi <- summary_d$factor_ci_upper
	
	graphics::par(mar = c(6, 5, 4, 2))
	graphics::plot(x, y, pch = 19, xaxt = "n", xlab = "MOFA factor", ylab = estimate_label,
			main = paste0("Association estimate for ", association_label, " across MOFA factors"))
	graphics::axis(1, at = x, labels = summary_d$factor, las = 2, cex.axis = 0.75)
	graphics::abline(h = 0, lty = 2)
	graphics::arrows(x0 = x, y0 = lo, x1 = x, y1 = hi, angle = 90, code = 3, length = 0.04)
	graphics::mtext("Intervals are approximate Wald 95% confidence intervals.", side = 3, line = 0.3, cex = 0.75)
}

plotPValuePage <- function(summary_d, association_label = "F") {
	# Plot factor-level p-values on the -log10 scale; the dashed line marks p = 0.05.
	plot_d <- summary_d
	plot_d$neg_log10_p <- safeLog10P(plot_d$factor_p_value)
	
	# For the p-value plot, sort factors by significance to make the ranking clear.
	plot_d <- plot_d[order(plot_d$factor_p_value, na.last = TRUE), , drop = FALSE]
	plot_d$factor <- makeFactorPlotLevels(plot_d$factor)
	
	if (useGgplotReportPlots()) {
		p <- ggplot2::ggplot(plot_d, ggplot2::aes(x = factor, y = neg_log10_p)) +
				ggplot2::geom_col() +
				ggplot2::geom_point() +
				ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
				ggplot2::labs(
						title = paste0("Evidence for association with ", association_label),
						subtitle = "Dashed line marks nominal p = 0.05; no multiple-testing correction is applied here.",
						x = "MOFA factor, sorted by nominal p-value",
						y = "-log10(p-value)"
				) +
				ggplot2::theme_bw() +
				ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
		
		print(p)
		return(invisible(NULL))
	}
	
	y <- plot_d$neg_log10_p
	x <- seq_along(y)
	
	graphics::par(mar = c(6, 5, 4, 2))
	graphics::plot(x, y, type = "h", lwd = 3, xaxt = "n", xlab = "MOFA factor, sorted by nominal p-value",
			ylab = "-log10(p-value)", main = paste0("Evidence for association with ", association_label))
	graphics::points(x, y, pch = 19)
	graphics::axis(1, at = x, labels = plot_d$factor, las = 2, cex.axis = 0.75)
	graphics::abline(h = -log10(0.05), lty = 2)
	graphics::mtext("Dashed line marks nominal p = 0.05; no multiple-testing correction is applied here.",
			side = 3, line = 0.3, cex = 0.75)
}

plotLinearPerformancePage <- function(summary_d) {
	# Plot R-squared and RMSE for linear models across factors.
	if (useGgplotReportPlots()) {
		plot_d <- rbind(
				data.frame(factor = summary_d$factor, metric = "R-squared", value = summary_d$r_squared),
				data.frame(factor = summary_d$factor, metric = "RMSE", value = summary_d$rmse)
		)
		plot_d$factor <- makeFactorPlotLevels(plot_d$factor)
		
		p <- ggplot2::ggplot(plot_d, ggplot2::aes(x = factor, y = value, group = metric)) +
				ggplot2::geom_line() +
				ggplot2::geom_point() +
				ggplot2::facet_wrap(stats::as.formula("~ metric"), scales = "free_y", ncol = 1) +
				ggplot2::labs(
						title = "Linear model performance across MOFA factors",
						x = "MOFA factor",
						y = "Metric value"
				) +
				ggplot2::theme_bw() +
				ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
		
		print(p)
		return(invisible(NULL))
	}
	
	x <- seq_len(nrow(summary_d))
	
	old_par <- graphics::par(no.readonly = TRUE)
	on.exit(graphics::par(old_par), add = TRUE)
	graphics::par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))
	
	graphics::plot(x, summary_d$r_squared, type = "b", pch = 19, xaxt = "n",
			xlab = "MOFA factor", ylab = "R-squared",
			main = "Linear model variance explained")
	graphics::axis(1, at = x, labels = summary_d$factor, las = 2, cex.axis = 0.70)
	
	graphics::plot(x, summary_d$rmse, type = "b", pch = 19, xaxt = "n",
			xlab = "MOFA factor", ylab = "RMSE",
			main = "Linear model prediction error")
	graphics::axis(1, at = x, labels = summary_d$factor, las = 2, cex.axis = 0.70)
}

plotLogisticPerformancePage <- function(summary_d) {
	# Plot AUC, accuracy, pseudo-R2, and Brier score for logistic models.
	if (useGgplotReportPlots()) {
		plot_d <- rbind(
				data.frame(factor = summary_d$factor, metric = "AUC", value = summary_d$auc),
				data.frame(factor = summary_d$factor, metric = "Accuracy", value = summary_d$accuracy),
				data.frame(factor = summary_d$factor, metric = "McFadden pseudo-R2", value = summary_d$pseudo_r2),
				data.frame(factor = summary_d$factor, metric = "Brier score", value = summary_d$brier)
		)
		plot_d$factor <- makeFactorPlotLevels(plot_d$factor)
		
		p <- ggplot2::ggplot(plot_d, ggplot2::aes(x = factor, y = value, group = metric)) +
				ggplot2::geom_line() +
				ggplot2::geom_point() +
				ggplot2::facet_wrap(stats::as.formula("~ metric"), scales = "free_y", ncol = 2) +
				ggplot2::labs(
						title = "Logistic model performance across MOFA factors",
						x = "MOFA factor",
						y = "Metric value"
				) +
				ggplot2::theme_bw() +
				ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
		
		print(p)
		return(invisible(NULL))
	}
	
	x <- seq_len(nrow(summary_d))
	
	old_par <- graphics::par(no.readonly = TRUE)
	on.exit(graphics::par(old_par), add = TRUE)
	graphics::par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
	
	graphics::plot(x, summary_d$auc, type = "b", pch = 19, xaxt = "n",
			xlab = "MOFA factor", ylab = "AUC", main = "Discrimination")
	graphics::axis(1, at = x, labels = summary_d$factor, las = 2, cex.axis = 0.65)
	
	graphics::plot(x, summary_d$accuracy, type = "b", pch = 19, xaxt = "n",
			xlab = "MOFA factor", ylab = "Accuracy", main = "0.5-threshold accuracy")
	graphics::axis(1, at = x, labels = summary_d$factor, las = 2, cex.axis = 0.65)
	
	graphics::plot(x, summary_d$pseudo_r2, type = "b", pch = 19, xaxt = "n",
			xlab = "MOFA factor", ylab = "McFadden pseudo-R2", main = "Relative deviance reduction")
	graphics::axis(1, at = x, labels = summary_d$factor, las = 2, cex.axis = 0.65)
	
	graphics::plot(x, summary_d$brier, type = "b", pch = 19, xaxt = "n",
			xlab = "MOFA factor", ylab = "Brier score", main = "Probability calibration error")
	graphics::axis(1, at = x, labels = summary_d$factor, las = 2, cex.axis = 0.65)
}

plotBestModelDiagnostics <- function(res, summary_d, model_type) {
	# Select the model with the smallest factor p-value for a more detailed plot.
	# Models are retrieved by factor name so this remains correct even if summary
	# tables are sorted differently from the original res list.
	best_name <- getBestModelName(summary_d)
	model <- res[[best_name]]
	
	if (model_type == "linear") {
		# For linear regression, show standard residual diagnostics for the selected model.
		old_par <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(old_par), add = TRUE)
		graphics::par(
				mfrow = c(2, 2),
				mar = c(5.2, 4.6, 2.8, 1.8),  # margins around each individual plot
				oma = c(1, 1, 3, 1),           # outer margins: bottom, left, top, right
				cex = 0.85
		)
		graphics::plot(model,main = "", 
				sub.caption = paste("Diagnostic plots for best model:", best_name))
	} else if (model_type == "logistic") {
		# For logistic regression, show predicted probability distributions by observed class.
		observed <- stats::model.response(stats::model.frame(model))
		prob <- stats::fitted(model)
		
		if (useGgplotReportPlots()) {
			plot_d <- data.frame(
					observed = observed,
					probability = prob,
					stringsAsFactors = FALSE
			)
			
			p <- ggplot2::ggplot(plot_d, ggplot2::aes(x = observed, y = probability)) +
					ggplot2::geom_boxplot() +
					ggplot2::geom_jitter(width = 0.12, height = 0, alpha = 0.45) +
					ggplot2::labs(
							title = paste(best_name, "predicted probabilities"),
							subtitle = "Diagnostic plot for the factor model with the smallest nominal F-term p-value.",
							x = "Observed outcome class",
							y = "Predicted probability"
					) +
					ggplot2::theme_bw()
			
			print(p)
			return(invisible(NULL))
		}
		
		y01 <- as.integer(observed) - 1L
		
		old_par <- graphics::par(no.readonly = TRUE)
		on.exit(graphics::par(old_par), add = TRUE)
		graphics::par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
		
		graphics::boxplot(prob ~ observed, xlab = "Observed outcome class", ylab = "Predicted probability",
				main = paste(best_name, "predicted probabilities"))
		
		graphics::hist(prob[y01 == 0], breaks = 10, xlim = c(0, 1), main = "Predicted probabilities",
				xlab = "Predicted probability")
		graphics::hist(prob[y01 == 1], breaks = 10, add = TRUE)
		graphics::legend("topright", legend = levels(observed), bty = "n", cex = 0.8)
	}
}

createRegressionPdfReport <- function(res, pdf_file, outcome = attr(res, "outcome"),
		covariates = attr(res, "covariates"),
		formula = attr(res, "formula"),
		factors_as_Y = attr(res, "factors_as_Y")) {
	# Create a PDF report summarizing all fitted regression models stored in res.
	# The function uses ggplot2 for higher-quality plots when available, with base R graphics as a fallback.
	if (!is.list(res) || length(res) == 0L) {
		stop("res must be a non-empty list of fitted regression models.")
	}
	
	# Summarize model coefficients and performance metrics into one data frame.
	summary_d <- summarizeRegressionModels(res)
	
	# Stop if the list contains unsupported model types.
	model_types <- unique(summary_d$model_type)
	if (any(model_types == "unknown")) {
		stop("At least one model in res is neither a binomial glm nor an lm object.")
	}
	
	# The current runRegression() implementation uses one model type per run, but
	# the report code checks this explicitly to avoid producing misleading plots.
	if (length(model_types) > 1L) {
		stop("This report expects all models in res to have the same regression type.")
	}
	model_type <- model_types[1]
	factors_as_Y_flag <- isTRUE(factors_as_Y)
	association_label <- if (factors_as_Y_flag) safeText(outcome, fallback = "Y") else "F"
	
	# Create the output directory if it does not already exist.
	out_dir <- dirname(pdf_file)
	if (!dir.exists(out_dir)) {
		dir.create(out_dir, recursive = TRUE)
	}
	
	# Open the PDF graphics device. Each plotting function below adds one page.
	grDevices::pdf(pdf_file, width = 11, height = 8.5, onefile = TRUE)
	on.exit(grDevices::dev.off(), add = TRUE)
	
	# Add the cover page with run-level context and interpretation guidance.
	plotReportTextPage(
			title = "MOFA factor regression report",
			lines = c(
					paste("Metadata variable Y:", safeText(outcome)),
					paste("Regression type:", ifelse(model_type == "logistic", "logistic regression", "linear regression")),
					paste("Model direction:", ifelse(factors_as_Y_flag, "Factor_i ~ Y + X", "Y ~ Factor_i + X")),
					paste("Model formula:", safeText(deparse(formula)[1])),
					paste("Adjustment covariates:", ifelse(length(covariates) == 0L, "none", paste(covariates, collapse = ", "))),
					paste("Number of fitted factor models:", nrow(summary_d)),
					paste("Plotting backend:", ifelse(useGgplotReportPlots(), "ggplot2", "base R graphics")),
					"",
					ifelse(
							factors_as_Y_flag,
							"Each model regresses one MOFA factor on Y plus the selected covariates.",
							"Each model regresses the selected metadata outcome on one MOFA factor plus the selected covariates."
					),
					ifelse(
							factors_as_Y_flag,
							"The report focuses on Y as the explanatory variable of interest, adjusted for X.",
							"The coefficient named F is the association between the current MOFA factor and the outcome, adjusted for X."
					),
					"The report emphasizes effect estimates, nominal p-values, and model-level performance metrics.",
					"P-values shown here are not adjusted for testing multiple MOFA factors."
			)
	)
	
	# Add type-specific explanatory text to make the report self-contained.
	if (factors_as_Y_flag) {
		plotReportTextPage(
				title = "How to interpret the factor-as-outcome linear-regression pages",
				lines = c(
						"Linear regression is used because the dependent variable in each model is a numeric MOFA factor value.",
						paste0("The explanatory variable of interest is Y = ", safeText(outcome), "."),
						"The reported p-value tests the association between Y and the current MOFA factor, adjusted for X.",
						"If Y is continuous or binary categorical, the report can usually show a single coefficient estimate.",
						"If Y is a categorical variable with more than two levels, the p-value is an omnibus test for Y and no single beta is shown.",
						"R-squared is the fraction of MOFA-factor variance explained by the full model.",
						"RMSE and MAE summarize residual error in the fitted MOFA factor values."
				)
		)
	} else if (model_type == "logistic") {
		plotReportTextPage(
				title = "How to interpret the logistic-regression pages",
				lines = c(
						"Logistic regression is used because the selected outcome is categorical/binary.",
						"The F coefficient is on the log-odds scale. exp(F coefficient) is the odds ratio for a one-unit increase in the MOFA factor.",
						"AUC measures class discrimination; 0.5 is no better than random and 1.0 is perfect discrimination.",
						"Accuracy is computed using a simple predicted-probability threshold of 0.5.",
						"McFadden pseudo-R2 compares residual deviance with null deviance; larger values indicate greater improvement over the null model.",
						"The Brier score is the mean squared difference between predicted probabilities and observed 0/1 outcomes; smaller is better.",
						"AIC and BIC are information criteria; smaller values indicate better tradeoffs between fit and model complexity."
				)
		)
	} else {
		plotReportTextPage(
				title = "How to interpret the linear-regression pages",
				lines = c(
						"Linear regression is used because the selected outcome is numeric/continuous.",
						"The F coefficient is the adjusted expected change in the outcome for a one-unit increase in the MOFA factor.",
						"R-squared is the fraction of outcome variance explained by the model.",
						"Adjusted R-squared penalizes model complexity and is often more useful when covariates are included.",
						"RMSE is the root mean squared residual and MAE is the mean absolute residual; smaller values indicate better fit.",
						"AIC and BIC are information criteria; smaller values indicate better tradeoffs between fit and model complexity.",
						"The final diagnostic page shows residual plots for the factor model with the smallest nominal F-term p-value."
				)
		)
	}
	
	# Add a table of model-level performance metrics.
	plotSummaryTablePage(summary_d, model_type, association_label = association_label)
	
	# Add cross-factor plots for effect estimates and p-values.
	plotFactorEffectPage(summary_d, model_type, association_label = association_label)
	plotPValuePage(summary_d, association_label = association_label)
	
	# Add regression-type-specific model performance plots.
	if (model_type == "logistic") {
		plotLogisticPerformancePage(summary_d)
	} else {
		plotLinearPerformancePage(summary_d)
	}
	
	# Add diagnostic plots for the strongest nominal factor association.
	plotBestModelDiagnostics(res, summary_d, model_type)
	
	# Return the summary data invisibly so the caller can inspect or save it.
	invisible(summary_d)
}

# -----------------------------------------------------------------------------
# Run analysis
# -----------------------------------------------------------------------------

# Load the selected data source once. In model mode this loads model.hdf5; in
# file mode it reads sample_metadata.{rds,xlsx} and factor_values.{rds,xlsx}.
regression_data <- loadRegressionData(regression_data_source, sample_id_col)
mo.tr <- regression_data$mo.tr
regression_sample_metadata <- regression_data$sample_metadata
regression_factor_values <- regression_data$factor_values

# Keep fitted-model lists, report-summary tables, and report paths available
# under stable analysis labels of the form Y=<name>__covar=<covariates>.
res <- list(
		outcome_as_Y = list(),
		factors_as_Y = list()
)
report_summary_d <- list(
		outcome_as_Y = list(),
		factors_as_Y = list()
)
regression_report_files <- list(
		outcome_as_Y = list(),
		factors_as_Y = list()
)

# Accumulate rows for one consolidated Excel workbook spanning all Y/X
# combinations and both regression directions.
excel_model_rows <- list()
excel_coefficient_rows <- list()

# Run both regression directions for every requested metadata variable and every
# requested covariate set.
for (outcome_i in Y_values) {
	for (covariates_i in X_sets) {
		analysis_label <- makeRegressionAnalysisLabel(outcome_i, covariates_i)
		
		message("Running outcome-as-Y regression for ", analysis_label)
		res_outcome_as_Y <- runRegression(
				mo.tr = mo.tr,
				sample_metadata = regression_sample_metadata,
				factor_values = regression_factor_values,
				sample_id_col = sample_id_col,
				X = covariates_i,
				Y = outcome_i,
				categorical = categorical,
				factors_as_Y = FALSE
		)
		
		outcome_as_Y_report_fname <- makeRegressionReportFilename(
				outcome = outcome_i,
				covariates = covariates_i,
				factors_as_Y = FALSE
		)
		outcome_as_Y_report_summary_d <- createRegressionPdfReport(
				res = res_outcome_as_Y,
				pdf_file = outcome_as_Y_report_fname,
				outcome = outcome_i,
				covariates = covariates_i,
				formula = attr(res_outcome_as_Y, "formula"),
				factors_as_Y = FALSE
		)
		message("Outcome-as-Y regression PDF report written to: ",
				outcome_as_Y_report_fname)
		
		message("Running factors-as-Y regression for ", analysis_label)
		res_factors_as_Y <- runRegression(
				mo.tr = mo.tr,
				sample_metadata = regression_sample_metadata,
				factor_values = regression_factor_values,
				sample_id_col = sample_id_col,
				X = covariates_i,
				Y = outcome_i,
				categorical = categorical,
				factors_as_Y = TRUE
		)
		
		factors_as_Y_report_fname <- makeRegressionReportFilename(
				outcome = outcome_i,
				covariates = covariates_i,
				factors_as_Y = TRUE
		)
		factors_as_Y_report_summary_d <- createRegressionPdfReport(
				res = res_factors_as_Y,
				pdf_file = factors_as_Y_report_fname,
				outcome = outcome_i,
				covariates = covariates_i,
				formula = attr(res_factors_as_Y, "formula"),
				factors_as_Y = TRUE
		)
		message("Factors-as-Y regression PDF report written to: ",
				factors_as_Y_report_fname)
		
		# Add both directions to the consolidated Excel summary tables.
		excel_model_rows[[paste0(analysis_label, "__outcome_as_Y")]] <-
				buildExcelModelRows(res_outcome_as_Y, outcome_i, covariates_i, FALSE)
		excel_model_rows[[paste0(analysis_label, "__factors_as_Y")]] <-
				buildExcelModelRows(res_factors_as_Y, outcome_i, covariates_i, TRUE)
		excel_coefficient_rows[[paste0(analysis_label, "__outcome_as_Y")]] <-
				buildExcelCoefficientRows(res_outcome_as_Y)
		excel_coefficient_rows[[paste0(analysis_label, "__factors_as_Y")]] <-
				buildExcelCoefficientRows(res_factors_as_Y)
		
		res$outcome_as_Y[[analysis_label]] <- res_outcome_as_Y
		res$factors_as_Y[[analysis_label]] <- res_factors_as_Y
		report_summary_d$outcome_as_Y[[analysis_label]] <- outcome_as_Y_report_summary_d
		report_summary_d$factors_as_Y[[analysis_label]] <- factors_as_Y_report_summary_d
		regression_report_files$outcome_as_Y[[analysis_label]] <-
				outcome_as_Y_report_fname
		regression_report_files$factors_as_Y[[analysis_label]] <-
				factors_as_Y_report_fname
	}
}

# Write one workbook after all requested models have been fitted.
regression_models_d <- do.call(rbind, excel_model_rows)
rownames(regression_models_d) <- NULL

# Arrange the Models worksheet for side-by-side inspection of nested models:
# unadjusted first, followed by models with progressively more covariates.
regression_models_d <- sortExcelModelRows(regression_models_d, Y_values, X_sets)
rownames(regression_models_d) <- NULL

regression_coefficients_d <- do.call(rbind, excel_coefficient_rows)
rownames(regression_coefficients_d) <- NULL
regression_excel_file <- file.path(
		regression_output_root, "regression_analysis", "regression_model_summary.xlsx")
writeRegressionExcelWorkbook(
		models_d = regression_models_d,
		coefficients_d = regression_coefficients_d,
		xlsx_file = regression_excel_file
)
message("Consolidated regression Excel summary written to: ", regression_excel_file)

# Store the same three tables in an RDS file. Each Excel worksheet is represented
# by one named list entry and is converted to a data matrix, as requested. Because
# the tables contain both character and numeric columns, R's matrix representation
# necessarily stores their values using a common (character) type.
regression_rds_data <- list(
		Models = as.matrix(regression_models_d),
		Coefficients = as.matrix(regression_coefficients_d),
		Dictionary = as.matrix(makeRegressionDictionary())
)
regression_rds_file <- file.path(
		regression_output_root, "regression_analysis", "regression_model_summary.rds")
saveRDS(regression_rds_data, regression_rds_file)
message("Consolidated regression RDS summary written to: ", regression_rds_file)

# Example inspection commands; uncomment after running the script if desired.
# names(res$outcome_as_Y)
# attr(res$outcome_as_Y[["Y=deltawtb__covar=none"]], "model_type")
# attr(res$outcome_as_Y[["Y=deltawtb__covar=none"]], "formula")
# summary(res$outcome_as_Y[["Y=deltawtb__covar=none"]][["Factor 1"]])
# lapply(res$factors_as_Y[["Y=deltawtb__covar=none"]], summary)

# If running logistic regression, the following prints the model parameters, in
# terms of odds ratios.
# sapply(res$outcome_as_Y[["Y=CR__covar=none"]], function(model){exp(coef(model))})
