# *****************************************************************************
# Convenience methods used by the MOFA, GLM, and MOFA factor-projection
# workflows. This file is sourced by MOFA_analysis.R, glm_factor_estimate.R,
# and project_mofa_factors.R before the reader functions listed in their
# respective parameter files are evaluated.
#
# The reader methods listed below will likely have to be amended for each study
# to accommodate study-specific differences in how omics and clinical data
# files are formatted.
# *****************************************************************************\

## ---------------------------------------------------------------------------
#' Read methylation data
#'
#' Read a methylation input file for use as an optional MOFA data view.
#'
#' This default reader expects the methylation dataset to be stored as an RDS
#' file. When this function is used in an `omics_input_defs` entry, the saved
#' object should be the numeric feature-by-sample matrix that MOFA_analysis.R
#' will use for the corresponding view: features in rows, sample barcodes in
#' columns, and column names matching the clinical barcode column named by that
#' view's `barcode_col` slot.
#'
#' @param fname Character scalar. Path to the methylation RDS file.
#' @return The object loaded from the RDS file. For direct use as an
#'   `omics_input_defs$reader`, this object is expected to be a numeric
#'   feature-by-sample matrix.
#' @examples
#' # readMethylation("path/to/methylation.rds")
#' @seealso readRDS
readMethylation <- function(fname){
	if (!is.character(fname) || length(fname) != 1) 
		stop("'fname' must be a single file path string")
	if (!file.exists(fname)) 
		stop("Methylation file not found: ", fname)

	# Load the serialized object. The corresponding omics_input_defs entry
	# determines whether this object is passed directly to MOFA_analysis.R or
	# further subset/transformed in a wrapper function.
	res <- tryCatch(
		readRDS(fname),
		error = function(e) stop("Failed to read methylation RDS file: ", e$message)
	)

	return(res)
}

#' Read proteomics data
#'
#' Read a proteomics input file for use as an optional MOFA data view.
#'
#' This default reader expects a CSV file with samples in rows and measured
#' proteins/features in columns. If the first column contains unique,
#' non-numeric sample identifiers, it is used as the source of sample barcodes.
#' The returned object follows the `omics_input_defs$reader` contract used by
#' MOFA_analysis.R: a numeric feature-by-sample matrix with features in rows,
#' sample barcodes in columns, and column names matching the clinical barcode
#' column named by the view's `barcode_col` slot.
#'
#' @param fname Character scalar. Path to the proteomics CSV file.
#' @return A numeric feature-by-sample matrix. The function raises an
#'   informative error if the file cannot be read or no numeric measurement
#'   columns are found.
#' @examples
#' # readProteomics("path/to/proteomics.csv")
#' @seealso read.csv
readProteomics <- function(fname){
	if (!is.character(fname) || length(fname) != 1) 
		stop("'fname' must be a single file path string")
	if (!file.exists(fname)) 
		stop("Proteomics file not found: ", fname)

	dd <- tryCatch(
		read.csv(fname, check.names = FALSE, stringsAsFactors = FALSE),
		error = function(e) stop("Failed to read proteomics CSV: ", e$message)
	)

	# Detect and handle a leading sample identifier column (common layout:
	# first column contains sample IDs and remaining columns are numeric
	# measurements). If the first column is non-numeric and contains unique
	# values, treat it as sample IDs and set rownames accordingly.
	if (ncol(dd) >= 1) {
		first_col <- dd[[1]]
		is_id_col <- is.character(first_col) || is.factor(first_col) || 
				!all(grepl("^-?\\d*\\.?\\d+(e[+-]?\\d+)?$", na.omit(as.character(first_col))))
		if (is_id_col && length(unique(first_col)) == nrow(dd)){
			rownames(dd) <- as.character(first_col)
			dd <- dd[, -1, drop = FALSE]
		}
	}

	# Keep numeric measurement columns, tolerating numeric values stored as text.
	numeric_cols <- vapply(dd, function(col) {
		vals <- na.omit(as.character(col))
		if (length(vals) == 0) return(TRUE)
		all(grepl("^-?\\d*\\.?\\d+(e[+-]?\\d+)?$", vals))
	}, logical(1))

	if (!any(numeric_cols)) stop("No numeric measurement columns found in proteomics file: ", fname)

	mat <- as.matrix(dd[, numeric_cols, drop = FALSE])
	# Coerce to numeric (column-wise); warn if coercion introduces NAs
	mat_num <- apply(mat, 2, function(x) as.numeric(as.character(x)))
	if (any(is.na(mat_num) & !is.na(mat))) warning("NAs introduced while coercing proteomics measurements to numeric for file: ", fname)

	# Transpose into the matrix shape expected by omics_input_defs$reader:
	# rows are features, columns are sample barcodes.
	res <- t(mat_num)
	rownames(res) <- colnames(dd)[numeric_cols]
	if (!is.null(rownames(dd)) && all(rownames(dd) != "")) colnames(res) <- rownames(dd)

	return(res)
}


#' Read metabolomics data
#'
#' Read a metabolomics input file and optional metabolite annotations.
#'
#' This default reader expects a CSV file whose first 24 columns contain
#' per-sample metadata and whose remaining columns contain metabolite
#' measurements. The measurement table is log2 transformed and transposed so
#' that `data` is a numeric feature-by-sample matrix.
#'
#' Unlike `readMethylation()` and `readProteomics()`, this function returns a
#' list because the metabolomics workflow may also need sample metadata and
#' metabolite annotations. In MOFA_parameters.R, the corresponding
#' `omics_input_defs` reader should pass `readMetabolomics(fn)[["data"]]` to
#' MOFA_analysis.R, because the `data` element is the object that satisfies the
#' MOFA view-reader contract: features in rows, sample barcodes in columns, and
#' column names matching the clinical barcode column named by `barcode_col`.
#'
#' The returned list contains:
#' - data: numeric feature-by-sample matrix used as the MOFA data view
#' - meta: data.frame containing the first 24 per-sample metadata columns
#' - annot: optional data.frame of metabolite annotations when `ann_fname` is
#'   provided and readable
#'
#' @param fname Character scalar. Path to the metabolomics CSV file.
#' @param ann_fname Character scalar or NULL. Path to an optional annotation CSV.
#' @return A list with elements `data`, `meta`, and (optionally) `annot`.
#' @details The function expects at least 24 columns (metadata + measurements).
#'   If the file contains fewer columns, an error is raised. Files are read with
#'   latin1 encoding to preserve non-ASCII characters sometimes found in
#'   measurement-unit labels. These requirements are specific to the Metabolon
#'   files data/annotation format.
#' @examples
#' # readMetabolomics("path/to/metabolomics.csv")
#' # readMetabolomics("path/to/metabolomics.csv", ann_fname = "path/to/annot.csv")
#' @seealso read.csv
readMetabolomics <- function(fname, ann_fname = NULL){
	if (!is.character(fname) || length(fname) != 1) stop("'fname' must be a single file path string")
	if (!file.exists(fname)) stop("Metabolomics file not found: ", fname)

	# Specify latin1 encoding because input files may contain non-ASCII
	# measurement-unit characters. Read safely and avoid converting strings
	# to factors.
	dd <- tryCatch(
		read.csv(fname, fileEncoding = "latin1", check.names = FALSE, stringsAsFactors = FALSE),
		error = function(e) stop("Failed to read metabolomics CSV: ", e$message)
	)

	# Separate per-sample metadata from the measurement matrix. The data
	# element returned below is the feature-by-sample matrix used by the
	# metabolomics entry in omics_input_defs.
	meta_cols <- 24L
	if (ncol(dd) < meta_cols) stop(sprintf("Expected at least %d columns (metadata + data); file has %d", meta_cols, ncol(dd)))
	meta <- dd[, seq_len(meta_cols), drop = FALSE]
	data <- dd[, -(seq_len(meta_cols)), drop = FALSE]
	rownames(data) <- meta$CLIENT_IDENTIFIER
	# Metabolon measurements are expected to be median centered. Log2
	# transform them to prepare the matrix for MOFA analysis.
	epsilon = min(data[data != 0])/2
	data = log2(data)
	# Transpose so rows are metabolite features and columns are sample
	# barcodes, matching the omics_input_defs$reader contract.
	res <- list(data = t(data), meta = meta)

	# Read in metabolite annotations, if provided.
	if (!is.null(ann_fname)){
		if (!file.exists(ann_fname)) warning("Annotation file not found: ", ann_fname) else {
			annot <- tryCatch(
				read.csv(ann_fname, check.names = FALSE, stringsAsFactors = FALSE),
				error = function(e){ warning("Failed to read annotation file: ", e$message); return(NULL) }
			)
			res$annot <- annot
		}
	}

	return(res)
}


#' Read phenotypic (clinical) patient data
#'
#' Read and return the phenotypic patient/clinical data from a CSV file.
#' MOFA_analysis.R links this table to each omics view using the clinical
#' barcode-column names configured in each `omics_input_defs` entry.
#'
#' @param fname Character scalar. Path to the phenotypic CSV file.
#' @param select Character vector, optional. If provided, only these column
#'   names will be returned (columns not found will be ignored with a warning).
#' @return A data.frame containing the phenotypic data (possibly subsetted).
#' @examples
#' # readPhenoData("path/to/pheno.csv")
#' # readPhenoData("pheno.csv", select = c("age", "sex", "visit"))
readPhenoData <- function(fname, study = NULL, select = NULL){
	if (!is.character(fname) || length(fname) != 1) stop("'fname' must be a single file path string")
	if (!file.exists(fname)) stop("Phenotypic file not found: ", fname)

	df <- tryCatch(
		read.csv(fname, check.names = FALSE, stringsAsFactors = FALSE),
		error = function(e) stop("Failed to read phenotypic CSV: ", e$message)
	)

	if (!is.null(select)){
		if (!is.character(select)) stop("'select' must be a character vector of column names")
		missing_cols <- setdiff(select, colnames(df))
		if (length(missing_cols) > 0) warning("The following requested columns were not found and will be ignored: ", paste(missing_cols, collapse = ", "))
		keep <- intersect(select, colnames(df))
		df <- df[, keep, drop = FALSE]
	}

	return(df)
}
