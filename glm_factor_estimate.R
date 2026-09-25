# *****************************************************************************
# Estimate MOFA-derived factor values in a new dataset using a trained GLM model.
#
# User-configurable settings and input paths are maintained in
# glm_factor_estimate_parameters.R. The script reuses file_readers.R and
# methods.R so clinical/timepoint/change-score handling is consistent with
# MOFA_analysis.R.
#
# Run from an R session as:
#
#   source("path_to_code/glm_factor_estimate.R", local = TRUE, chdir = TRUE)
# 
# Alternatively, run the script from the command line:
#
#   Rscript.exe -e "source('path_to_code/glm_factor_estimate.R', local = TRUE, 
#				chdir = TRUE)"
#
# *****************************************************************************

source("file_readers.R", local = TRUE)
source("glm_factor_estimate_parameters.R", local = TRUE)
source("methods.R", local = TRUE)

#' Predict MOFA-derived factor scores in an external dataset
#'
#' This function applies trained view-specific MOFA factor predictors, as returned
#' by trainMofaFactorPredictor(), to new external data views. Each external data
#' view should be a features x samples matrix or data frame, on the same feature
#' scale as the data used to train the corresponding predictor.
#'
#' @param predictor A "MofaFactorPredictor" object returned by
#'   trainMofaFactorPredictor().
#' @param new_data_views A named list of external data views. Each entry must be a
#'   numeric matrix or data frame with features in rows and samples in columns.
#' @param views Character vector of views to use, or "all" to use all trained
#'   views available in predictor and present in new_data_views. Default is "all".
#' @param factors Character vector of factors to predict, or "all" to use all
#'   trained factors available for the selected views. Default is "all".
#' @param combine_views How to combine multiple view-specific predictions for the
#'   same factor. One of "mean", "median", "best_cv", or "none". If "none",
#'   view-specific predictions are returned without combining. Default is "mean".
#' @param metric_for_best_cv Metric used when combine_views = "best_cv". Higher
#'   values are treated as better for Pearson, Spearman, and R-squared; lower
#'   values are treated as better for RMSE and MAE. Default is "pearson".
#' @param lambda_choice Optional lambda choice for prediction. If NULL, the
#'   lambda choice stored in each trained model is used. Otherwise one of
#'   "lambda.min" or "lambda.1se". Default NULL.
#' @param verbose Logical; if TRUE, progress messages are printed. Default TRUE.
#'
#' @return A list containing view-specific predictions and, unless
#'   combine_views = "none", a combined factors x samples prediction matrix.
#'
#' @details
#' Unknown views in new_data_views are ignored. For each selected model, features
#' are matched by row names. Missing trained features are allowed, but prediction
#' can proceed only if at least one trained feature is present. Features missing
#' from the external data are filled with zero after alignment. This is
#' appropriate only if the training data were centered/scaled so that zero
#' represents the training mean for a feature. If that assumption is not valid,
#' external data should be harmonized before calling this function.
#'
#' @examples
#' \dontrun{
#' pred <- predictMofaFactorScores(
#'   predictor = predictor,
#'   new_data_views = list(RNA = new_rna, Proteomics = new_protein),
#'   combine_views = "mean"
#' )
#' pred$combined
#' }
#'
#' @export
predictMofaFactorScores <- function(
  predictor,
  new_data_views,
  views = "all",
  factors = "all",
  combine_views = c("mean", "median", "best_cv", "none"),
  metric_for_best_cv = "pearson",
  lambda_choice = NULL,
  verbose = TRUE
) {
  requireNamespace("glmnet", quietly = TRUE) ||
    stop("The glmnet package is required but is not installed.", call. = FALSE)

  if (!inherits(predictor, "MofaFactorPredictor")) {
    stop("predictor must be an object returned by trainMofaFactorPredictor().", call. = FALSE)
  }

  if (!is.list(new_data_views) || length(new_data_views) == 0L) {
    stop("new_data_views must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(new_data_views)) || any(names(new_data_views) == "")) {
    stop("new_data_views must be a named list of features x samples matrices.", call. = FALSE)
  }

  combine_views <- match.arg(combine_views)
  if (!is.null(lambda_choice)) {
    lambda_choice <- match.arg(lambda_choice, choices = c("lambda.min", "lambda.1se"))
  }

  trained_views <- names(predictor$results)
  trained_views <- trained_views[vapply(predictor$results, length, integer(1)) > 0L]
  available_views <- intersect(trained_views, names(new_data_views))
  selected_views <- .resolve_requested_names(views, available_views, "views")

  view_specific <- list()

  for (view_name in selected_views) {
    if (verbose) {
      message("Predicting from view: ", view_name)
    }

    X_new <- .as_numeric_matrix(
      new_data_views[[view_name]],
      object_name = paste0("new_data_views[['", view_name, "']]")
    )

    # External views may use the same MOFA-style padding convention, where
    # all-NA columns represent samples with no measurements in this view. Drop
    # those sample columns before feature alignment and prediction.
    X_new <- .drop_all_na_sample_columns(
      X_new,
      object_name = paste0("new_data_views[['", view_name, "']]")
    )

    trained_factors <- names(predictor$results[[view_name]])
    selected_factors <- .resolve_requested_names(factors, trained_factors, "factors")

    view_specific[[view_name]] <- list()

    for (factor_name in selected_factors) {
      model_info <- predictor$results[[view_name]][[factor_name]]
      if (is.null(model_info)) {
        next
      }

      model_features <- model_info$features
      common_features <- intersect(model_features, rownames(X_new))

      if (length(common_features) == 0L) {
        warning(sprintf(
          "Skipping view '%s', factor '%s': no trained features found in external data.",
          view_name, factor_name
        ))
        next
      }

      # Align external features to the exact feature order used for training.
      # Missing features are set to zero. This assumes zero is the centered value.
      X_aligned <- matrix(
        0,
        nrow = length(model_features),
        ncol = ncol(X_new),
        dimnames = list(model_features, colnames(X_new))
      )
      X_aligned[common_features, ] <- X_new[common_features, , drop = FALSE]

      if (any(!is.finite(X_aligned))) {
        stop(sprintf(
          "External data for view '%s' contains non-finite values after alignment.",
          view_name
        ), call. = FALSE)
      }

      pred_lambda <- if (is.null(lambda_choice)) model_info$lambda_choice else lambda_choice
      pred_type <- if (model_info$outcome_type == "binary") "response" else "link"

      pred <- stats::predict(
        model_info$final_model,
        newx = t(X_aligned),
        s = pred_lambda,
        type = pred_type
      )

      pred <- as.numeric(pred)
      names(pred) <- colnames(X_new)

      view_specific[[view_name]][[factor_name]] <- pred
    }
  }

  view_matrices <- .prediction_lists_to_matrices(view_specific)

  if (combine_views == "none") {
    return(list(
      view_specific = view_matrices,
      combined = NULL,
      combine_views = combine_views
    ))
  }

  combined <- .combine_view_predictions(
    view_matrices = view_matrices,
    predictor = predictor,
    method = combine_views,
    metric = metric_for_best_cv
  )

  list(
    view_specific = view_matrices,
    combined = combined,
    combine_views = combine_views,
    metric_for_best_cv = metric_for_best_cv
  )
}

# -----------------------------------------------------------------------------
# Internal prediction helper functions
# -----------------------------------------------------------------------------

.resolve_requested_names <- function(requested, available, label) {
  if (length(requested) == 1L && identical(requested, "all")) {
    return(available)
  }

  missing <- setdiff(requested, available)
  if (length(missing) > 0L) {
    warning(sprintf(
      "Ignoring unknown %s: %s",
      label,
      paste(missing, collapse = ", ")
    ))
  }

  selected <- intersect(requested, available)
  if (length(selected) == 0L) {
    stop("No requested ", label, " are available.", call. = FALSE)
  }

  selected
}

.as_numeric_matrix <- function(x, object_name = "object") {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop(object_name, " must be a matrix or data frame.", call. = FALSE)
  }

  x_mat <- as.matrix(x)

  if (!is.numeric(x_mat)) {
    stop(object_name, " must contain only numeric values.", call. = FALSE)
  }

  if (is.null(rownames(x_mat)) || anyNA(rownames(x_mat)) || any(rownames(x_mat) == "")) {
    stop(object_name, " must have non-empty row names.", call. = FALSE)
  }

  if (is.null(colnames(x_mat)) || anyNA(colnames(x_mat)) || any(colnames(x_mat) == "")) {
    stop(object_name, " must have non-empty column names.", call. = FALSE)
  }

  x_mat
}

.drop_all_na_sample_columns <- function(x, object_name = "object") {
  # Input matrices are expected to be features x samples. Columns whose values
  # are all NA represent samples missing entirely from that data view; these are
  # valid in MOFA input matrices but unusable for supervised prediction.
  all_na_sample <- apply(x, 2, function(z) all(is.na(z)))

  if (any(all_na_sample)) {
    message(
      object_name,
      ": dropping ",
      sum(all_na_sample),
      " all-NA sample column(s): ",
      paste(colnames(x)[all_na_sample], collapse = ", ")
    )
    x <- x[, !all_na_sample, drop = FALSE]
  }

  if (ncol(x) == 0L) {
    stop(object_name, " has no sample columns remaining after dropping all-NA columns.", call. = FALSE)
  }

  x
}

.prediction_lists_to_matrices <- function(view_specific) {
  out <- list()

  for (view_name in names(view_specific)) {
    factor_list <- view_specific[[view_name]]
    if (length(factor_list) == 0L) next

    samples <- Reduce(intersect, lapply(factor_list, names))
    factors <- names(factor_list)

    mat <- matrix(
      NA_real_,
      nrow = length(factors),
      ncol = length(samples),
      dimnames = list(factors, samples)
    )

    for (factor_name in factors) {
      mat[factor_name, ] <- factor_list[[factor_name]][samples]
    }

    out[[view_name]] <- mat
  }

  out
}

.combine_view_predictions <- function(view_matrices, predictor, method, metric) {
  if (length(view_matrices) == 0L) {
    stop("No view-specific predictions are available to combine.", call. = FALSE)
  }

  common_factors <- Reduce(intersect, lapply(view_matrices, rownames))
  common_samples <- Reduce(intersect, lapply(view_matrices, colnames))

  if (length(common_factors) == 0L || length(common_samples) == 0L) {
    stop("View-specific predictions do not share common factors and samples.", call. = FALSE)
  }

  if (method %in% c("mean", "median")) {
    arr <- array(
      NA_real_,
      dim = c(length(common_factors), length(common_samples), length(view_matrices)),
      dimnames = list(common_factors, common_samples, names(view_matrices))
    )

    for (view_name in names(view_matrices)) {
      arr[, , view_name] <- view_matrices[[view_name]][common_factors, common_samples, drop = FALSE]
    }

    if (method == "mean") {
      combined <- apply(arr, c(1, 2), mean, na.rm = TRUE)
    } else {
      combined <- apply(arr, c(1, 2), stats::median, na.rm = TRUE)
    }

    return(matrix(
      combined,
      nrow = length(common_factors),
      ncol = length(common_samples),
      dimnames = list(common_factors, common_samples)
    ))
  }

  if (method == "best_cv") {
    combined <- matrix(
      NA_real_,
      nrow = length(common_factors),
      ncol = length(common_samples),
      dimnames = list(common_factors, common_samples)
    )

    lower_is_better <- metric %in% c("rmse", "mae")

    for (factor_name in common_factors) {
      scores <- sapply(names(view_matrices), function(view_name) {
        res <- predictor$results[[view_name]][[factor_name]]
        if (is.null(res)) return(NA_real_)
        ms <- res$internal_cv$metrics_summary
        row <- ms[ms$metric == metric, , drop = FALSE]
        if (nrow(row) == 0L) return(NA_real_)
        row$mean[1]
      })

      if (all(is.na(scores))) {
        next
      }

      best_view <- if (lower_is_better) {
        names(scores)[which.min(scores)]
      } else {
        names(scores)[which.max(scores)]
      }

      combined[factor_name, ] <- view_matrices[[best_view]][factor_name, common_samples]
    }

    return(combined)
  }

  stop("Unknown combine method: ", method, call. = FALSE)
}


# -----------------------------------------------------------------------------
# Execute the configured external-scoring workflow.
# -----------------------------------------------------------------------------

# Validate and load the trained predictor object.
if (!is.character(glm_predictor_fn) || length(glm_predictor_fn) != 1L ||
    is.na(glm_predictor_fn) || !nzchar(glm_predictor_fn)) {
  stop("'glm_predictor_fn' must be a single non-empty character string.", call. = FALSE)
}
if (!file.exists(glm_predictor_fn)) {
  stop("GLM predictor file does not exist: ", glm_predictor_fn, call. = FALSE)
}

predictor <- readRDS(glm_predictor_fn)
if (!inherits(predictor, "MofaFactorPredictor")) {
  stop("The object in glm_predictor_fn is not a MofaFactorPredictor.", call. = FALSE)
}

# Read and validate clinical data using the same shared logic as MOFA_analysis.R.
pheno_d <- readPhenoData(getResourceFile("clin_fn", required = TRUE))
clinical_context <- validateClinicalConfiguration(pheno_d)
baseline_timepoint <- clinical_context$baseline_timepoint
target_timepoint_label <- clinical_context$target_timepoint_label

validateOmicsInputDefinitions(omics_input_defs)
loaded_info <- loadConfiguredOmicsViews(omics_input_defs)
loaded_views <- loaded_info$loaded_views
barcode_cols <- loaded_info$barcode_cols

if (length(loaded_views) == 0L) {
  stop("No omics data file names are defined in glm_factor_estimate_parameters.R.", call. = FALSE)
}

# Select the target samples and, if requested, compute target-minus-baseline
# change scores. Column names are converted to the same standardized sample IDs
# used by MOFA_analysis.R.
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

# Apply optional feature-wise z-score scaling after timepoint/change-score
# preparation. External data must be placed on the same scale as the data used
# to train the predictor.
if (isTRUE(scale_views)) {
  data_views <- lapply(data_views, scaleOmicsViewFeatures)
}

# Estimate factor values. Factors default to all factor models available in each
# selected view; the user-facing parameter file controls the remaining options.
prediction <- predictMofaFactorScores(
  predictor = predictor,
  new_data_views = data_views,
  views = views,
  factors = "all",
  combine_views = combine_views,
  metric_for_best_cv = metric_for_best_cv,
  lambda_choice = lambda_choice,
  verbose = verbose
)

# Determine the output root. When output_dir is NULL, place results beside the
# predictor RDS file.
if (is.null(output_dir)) {
  output_root <- dirname(normalizePath(glm_predictor_fn))
} else {
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("'output_dir' must be NULL or a single non-empty character path.", call. = FALSE)
  }
  output_root <- output_dir
}

score_tag <- ifelse(use_change_scores, "change_scores", "timepoint_values")
actual_views <- intersect(names(data_views), names(predictor$results))
if (!(length(views) == 1L && identical(views, "all"))) {
  actual_views <- intersect(actual_views, views)
}
if (length(actual_views) == 0L) {
  stop("No loaded data views match the trained predictor views requested for scoring.", call. = FALSE)
}

views_tag <- paste(vapply(actual_views, function(vn) {
  if (!is.null(omics_input_defs[[vn]]$moniker)) omics_input_defs[[vn]]$moniker else vn
}, character(1)), collapse = "+")

lambda_tag <- if (is.null(lambda_choice)) "trained" else gsub("\\.", "", lambda_choice)
metric_tag <- if (combine_views == "best_cv") paste0("_metric=", metric_for_best_cv) else ""
res_folder <- file.path(
  output_root,
  paste0(
    "glm_estimate_visit=", target_timepoint_label,
    "_", score_tag,
    "_views=", views_tag,
    "_combine=", combine_views,
    metric_tag,
    "_lambda=", lambda_tag
  )
)
log_folder <- file.path(res_folder, "logs")
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)

requireNamespace("xlsx", quietly = TRUE) ||
  stop("The xlsx package is required to write the workflow Excel outputs.", call. = FALSE)

# Export sample metadata using the same construction as MOFA_analysis.R, but
# restrict it to the samples for which factor values were actually estimated.
# This matches MOFA_analysis.R, where metadata are ultimately subset to the
# samples retained in the prepared omics views.
estimated_sample_ids <- if (identical(combine_views, "none")) {
  view_sample_ids <- lapply(prediction$view_specific, colnames)
  unique(unlist(view_sample_ids, use.names = FALSE))
} else {
  colnames(prediction$combined)
}

if (length(estimated_sample_ids) == 0L) {
  stop("No sample IDs were returned with the estimated factor values.", call. = FALSE)
}

sample_metadata <- buildTargetVisitSampleMetadata(pheno_d)
missing_metadata <- setdiff(estimated_sample_ids, rownames(sample_metadata))
if (length(missing_metadata) > 0L) {
  stop(
    "Sample metadata are missing for samples with estimated factor values: ",
    paste(missing_metadata, collapse = ", "),
    call. = FALSE
  )
}
sample_metadata <- sample_metadata[estimated_sample_ids, , drop = FALSE]
sample_metadata <- data.frame(sample_metadata, sample = rownames(sample_metadata))
xlsx::write.xlsx(
  as.data.frame(sample_metadata),
  file.path(res_folder, "sample_metadata.xlsx"),
  sheetName = study_name,
  row.names = FALSE
)
saveRDS(sample_metadata, file.path(res_folder, "sample_metadata.rds"))

xlsx_file <- file.path(res_folder, "glm_factor_values.xlsx")
rds_file <- file.path(res_folder, "glm_factor_values.rds")

if (combine_views == "none") {
  # With no cross-view combination there is no single factor matrix. Preserve
  # all view-specific results in the same two output files: an RDS named list
  # and an Excel workbook with one worksheet per view.
  view_tables <- lapply(prediction$view_specific, function(mat) t(mat))
  saveRDS(view_tables, rds_file)

  first_sheet <- TRUE
  for (view_name in names(view_tables)) {
    xlsx::write.xlsx(
      as.data.frame(view_tables[[view_name]]),
      xlsx_file,
      sheetName = substr(view_name, 1L, 31L),
      row.names = TRUE,
      append = !first_sheet
    )
    first_sheet <- FALSE
  }
} else {
  if (is.null(prediction$combined)) {
    stop("No combined factor estimates were returned.", call. = FALSE)
  }

  # predictMofaFactorScores() returns factors x samples; the workflow output is
  # samples x factors, with standardized sample IDs as row names.
  factor_values <- t(prediction$combined)
  xlsx::write.xlsx(
    as.data.frame(factor_values),
    xlsx_file,
    sheetName = substr(study_name, 1L, 31L),
    row.names = TRUE
  )
  saveRDS(factor_values, rds_file)
}

# Record session/source provenance alongside the estimates.
session_file <- file.path(log_folder, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)
file.copy(
  from = c(
    "glm_factor_estimate.R", "glm_factor_estimate_parameters.R",
    "file_readers.R", "methods.R"
  ),
  to = log_folder,
  overwrite = TRUE
)

# Generate run-specific README prose using the same view-specific notes as the
# MOFA workflow.
score_desc <- if (use_change_scores) {
  paste0("change scores (", target_timepoint_label, " - baseline)")
} else {
  paste0(
    "source values from samples where ", timepoint_col, " == ", timepoint,
    " (", target_timepoint_label, ")"
  )
}

included_views <- paste(actual_views, collapse = ", ")
data_notes <- vapply(actual_views, function(view_name) {
  omics_input_defs[[view_name]]$readme_notes[[score_tag]]
}, character(1))

combine_desc <- if (combine_views == "none") {
  "View-specific factor estimates were retained separately and were not combined."
} else if (combine_views == "best_cv") {
  paste0(
    "For each factor, predictions from the view with the best mean internal-CV ",
    metric_for_best_cv, " metric were used."
  )
} else {
  paste0("View-specific factor estimates were combined using the ", combine_views, ".")
}

readme_text <- paste(c(
  paste(study_name, "GLM-derived MOFA factor values were estimated from", score_desc, "."),
  paste("Included omics views:", included_views, "."),
  data_notes,
  paste("The trained predictor was read from", normalizePath(glm_predictor_fn), "."),
  paste("Estimated factor values were written to", basename(xlsx_file), "and",
      basename(rds_file), "."),
  "Sample metadata for the configured target visit were written to sample_metadata.xlsx and sample_metadata.rds.",
  combine_desc,
  if (isTRUE(scale_views)) {
    "Prepared external features were z-score scaled across selected samples before prediction."
  } else {
    "No additional feature scaling was applied after timepoint/change-score preparation."
  },
  "The logs folder contains sessionInfo.txt, which records the R/package",
  "environment, plus copies of glm_factor_estimate.R,",
  "glm_factor_estimate_parameters.R, file_readers.R, and methods.R used for",
  "the run."
), collapse = " ")
writeLines(strwrap(readme_text, width = 80), con = file.path(res_folder, "README.txt"))

message("GLM factor estimation completed. Results written to: ", res_folder)
