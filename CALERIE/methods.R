# *****************************************************************************
# Convenience methods, to use in other places, as needed. This file should be
# sourced within other R scripts, where the methods herein are utilized.
# *****************************************************************************

## ---------------------------------------------------------------------------
#' Read methylation data
#'
#' Read methylation data for the specified study. For the CALERIE study the
#' methylation dataset is expected to be stored as an RDS file and this
#' function will load and return the object saved within that file.
#'
#' @param fname Character scalar. Path to the methylation file (RDS format for
#'   the CALERIE study).
#' @param study Character scalar. Study identifier. Defaults to "CALERIE".
#' @return The object loaded from the RDS file (often a data.frame, matrix or
#'   list) for CALERIE. The function will raise an error if the file does not
#'   exist or cannot be read.
#' @examples
#' # readMethylation("path/to/methylation.rds")
#' @seealso readRDS
readMethylation <- function(fname, study = c("CALERIE")){
	study <- match.arg(study)

	if (!is.character(fname) || length(fname) != 1) stop("'fname' must be a single file path string")
	if (!file.exists(fname)) stop("Methylation file not found: ", fname)

	# For CALERIE study the methylation data are provided as an RDS object.
	# Use tryCatch to provide a clearer error message if reading fails.
	res <- tryCatch(
		readRDS(fname),
		error = function(e) stop("Failed to read methylation RDS file: ", e$message)
	)

	return(res)
}

#' Read proteomics data
#'
#' Read proteomics data for the specified study. For the CALERIE study the
#' proteomics dataset is expected to be a CSV file where rows correspond to
#' samples and columns to measured proteins/features. This function reads the
#' table and returns a transposed numeric matrix (features x samples) so that
#' columns correspond to sample identifiers when available.
#'
#' @param fname Character scalar. Path to the proteomics file (CSV format for
#'   the CALERIE study).
#' @param study Character scalar. Study identifier. Defaults to "CALERIE".
#' @return A numeric matrix with features as rows and samples as columns. The
#'   function throws an informative error if the file cannot be read or no
#'   numeric measurement columns are found.
#' @examples
#' # readProteomics("path/to/proteomics.csv")
#' @seealso read.csv
readProteomics <- function(fname, study = c("CALERIE")){
	study <- match.arg(study)

	if (!is.character(fname) || length(fname) != 1) stop("'fname' must be a single file path string")
	if (!file.exists(fname)) stop("Proteomics file not found: ", fname)

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

	# Identify numeric measurement columns (tolerant to numeric-like strings)
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

	# Transpose so rows = features, cols = samples. After transposition:
	# rownames(res) should be feature names (original column names), and
	# colnames(res) should be sample IDs (original rownames or numeric index)
	res <- t(mat_num)
	rownames(res) <- colnames(dd)[numeric_cols]
	if (!is.null(rownames(dd)) && all(rownames(dd) != "")) colnames(res) <- rownames(dd)

	return(res)
}


#' Read metabolomics data
#'
#' Read metabolomics data for the specified study. For the CALERIE study the
#' input is expected to be a CSV file where the first 24 columns contain
#' per-sample metadata and the remaining columns contain metabolite measurements.
#' Further, measurements are median centered. Accordingly, the data table is 
#' log2-tranformed before being returned.
#'
#' The function returns a list with elements:
#' - data: a data.frame of metabolite measurements (columns after the first 24)
#' - meta: a data.frame containing the first 24 metadata columns
#' - annot: (optional) a data.frame of metabolite annotations when `ann_fname`
#'          is provided and readable
#'
#' @param fname Character scalar. Path to the metabolomics CSV file.
#' @param ann_fname Character scalar or NULL. Path to an optional annotation CSV.
#' @param study Character scalar. Study identifier. Defaults to "CALERIE".
#' @return A list with elements `data`, `meta`, and (optionally) `annot`.
#' @details The function expects at least 24 columns (metadata + data). If the
#'   file contains fewer columns an error is raised. For CALERIE, files are read
#'   with latin1 encoding to preserve non-ASCII characters (e.g. the mu symbol).
#' @examples
#' # readMetabolomics("path/to/metabolomics.csv")
#' # readMetabolomics("path/to/metabolomics.csv", ann_fname = "path/to/annot.csv")
#' @seealso read.csv
readMetabolomics <- function(fname, ann_fname = NULL, study = c("CALERIE")){
	study <- match.arg(study)

	if (!is.character(fname) || length(fname) != 1) stop("'fname' must be a single file path string")
	if (!file.exists(fname)) stop("Metabolomics file not found: ", fname)

	if (study == "CALERIE"){
		# Specify latin1 encoding because files may contain non-ASCII characters
		# (specifically, the 'µ' character, e.g. in 'µl'). Read safely and avoid
		# converting strings to factors.
		dd <- tryCatch(
			read.csv(fname, fileEncoding = "latin1", check.names = FALSE, stringsAsFactors = FALSE),
			error = function(e) stop("Failed to read metabolomics CSV: ", e$message)
		)

		# Separate metadata from data and return both as separate members of a list.
		meta_cols <- 24L
		if (ncol(dd) < meta_cols) stop(sprintf("Expected at least %d columns (metadata + data); file has %d", meta_cols, ncol(dd)))
		meta <- dd[, seq_len(meta_cols), drop = FALSE]
		data <- dd[, -(seq_len(meta_cols)), drop = FALSE]
		rownames(data) <- meta$CLIENT_IDENTIFIER
		# CALERIE Metabolon data is median centered. Log2 transform them, to prep
		# them for MOFA analsis
		epsilon = min(data[data != 0])/2
		data = log2(data)
		# Transpose the data table so that samples ids appear as column names
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
	} else {
		# Generic fallback: attempt to read the CSV and return as 'data'.
		dd <- tryCatch(
			read.csv(fname, check.names = FALSE, stringsAsFactors = FALSE),
			error = function(e) stop("Failed to read CSV: ", e$message)
		)
		res <- list(data = dd)
	}

	return(res)
}


#' Read phenotypic (clinical) patient data
#'
#' Read and return the phenotypic patient/clinical data from a CSV file.
#' For the CALERIE study the file is expected to be a plain CSV; this
#' function reads the file and optionally returns a subset of columns.
#'
#' @param fname Character scalar. Path to the phenotypic CSV file.
#' @param study Character scalar. Study identifier. Defaults to "CALERIE".
#' @param select Character vector, optional. If provided, only these column
#'   names will be returned (columns not found will be ignored with a warning).
#' @return A data.frame containing the phenotypic data (possibly subsetted).
#' @examples
#' # readPhenoData("path/to/pheno.csv")
#' # readPhenoData("pheno.csv", select = c("age", "sex", "visit"))
readPhenoData <- function(fname, study = c("CALERIE"), select = NULL){
	study <- match.arg(study)

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

## ---------------------------------------------------------------------------
#' Fill missing sample columns with NA
#'
#' Given a data matrix (features x samples) and a vector of desired sample
#' identifiers (`sids`), return a matrix containing all `sids` as columns. Any
#' samples missing from the input `dmat` will be added as columns filled with
#' NA. The returned matrix has columns ordered as `sids`.
#'
#' @param dmat A matrix-like object with features in rows and sample IDs in
#'   columns (column names). If `dmat` is a data.frame it will be coerced to a
#'   matrix for the operation; rownames are preserved.
#' @param sids Character vector of sample IDs that should appear as column
#'   names in the returned matrix. Order of columns in the result matches
#'   `sids`.
#' @return A matrix with the same number of rows as `dmat` and columns equal to
#'   `sids`. Columns that were not present in `dmat` are filled with NA.
#' @examples
#' # dmat has columns A and B; request columns A,B,C -> C will be all NA
#' # m <- fillEmptyCols(dmat, c("A","B","C"))
fillEmptyCols <- function(dmat, sids){
	if (missing(sids) || length(sids) == 0) stop("'sids' must be a non-empty character vector of sample IDs")

	# Coerce to matrix but preserve rownames
	rn <- rownames(dmat)
	mat <- as.matrix(dmat)
	if (!is.null(rn)) rownames(mat) <- rn

	# Ensure sids are character and unique
	sids <- as.character(sids)
	if (anyDuplicated(sids)) sids <- unique(sids)

	# Determine which sample IDs are missing and need NA columns
	existing <- colnames(mat)
	missing_cols <- setdiff(sids, existing)

	if (length(missing_cols) > 0){
		# create NA matrix with same rownames and appropriate column names
		nc_mat <- matrix(NA, nrow = nrow(mat), ncol = length(missing_cols))
		if (!is.null(rownames(mat))) rownames(nc_mat) <- rownames(mat)
		colnames(nc_mat) <- missing_cols
		mat <- cbind(mat, nc_mat)
	}

	# Reorder columns to match sids; if some requested sids are still absent
	# (e.g. when sids contained names not representable), this will introduce NA
	# columns for them by creating them explicitly.
	absent_after <- setdiff(sids, colnames(mat))
	if (length(absent_after) > 0){
		# create explicit NA columns for these and bind
		extra <- matrix(NA, nrow = nrow(mat), ncol = length(absent_after))
		if (!is.null(rownames(mat))) rownames(extra) <- rownames(mat)
		colnames(extra) <- absent_after
		mat <- cbind(mat, extra)
	}

	# Select and return columns in the requested order
	res <- mat[, sids, drop = FALSE]
	return(res)
}

